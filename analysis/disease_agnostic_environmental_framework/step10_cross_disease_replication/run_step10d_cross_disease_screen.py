"""Run the frozen 29-test screen for the frozen Step 10 disease panel.

The disease panel and all eligibility decisions are read from the Step 10C
lock. Each disease is fitted separately with the same survey-weighted model
family and its own BH family of 29 tests. No Step 10 downstream biology is
implemented here by design.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import sys
import types
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
STEP10_DIR = Path(__file__).resolve().parent
DEFAULT_TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
DEFAULT_REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_PANEL = STEP10_DIR / "step10c_randomized_disease_panel.csv"
FDR_DENOMINATOR = 29

CYCLES = [
    "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
    "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018",
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_statistical_model(path: Path, name: str):
    """Load the shared survey model without requiring its plotting dependency."""
    # The shared legacy model imports matplotlib for unrelated plotting
    # helpers. Step 10 invokes only its harmonization and survey estimator;
    # keeping this small import shim local avoids changing the model or
    # installing a plotting stack just to run the screen.
    if "matplotlib" not in sys.modules:
        mpl = types.ModuleType("matplotlib")
        mpl.use = lambda *_args, **_kwargs: None
        pyplot = types.ModuleType("matplotlib.pyplot")
        mpl.pyplot = pyplot
        sys.modules["matplotlib"] = mpl
        sys.modules["matplotlib.pyplot"] = pyplot
    return load_module(path, name)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def normalize_registry_paths(registry: pd.DataFrame) -> pd.DataFrame:
    """Make provenance paths readable when the runner is executed in WSL."""
    out = registry.copy()
    if "local_xpt" not in out.columns:
        return out
    def one(value: object) -> str:
        text = str(value)
        match = re.match(r"^([A-Za-z]):[\\/](.*)$", text)
        if match:
            return "/mnt/" + match.group(1).lower() + "/" + match.group(2).replace("\\", "/")
        return text
    out["local_xpt"] = out["local_xpt"].map(one)
    return out


def derive_core(model, demo: pd.DataFrame, bmx: pd.DataFrame, smq: pd.DataFrame, cycle_index: int) -> pd.DataFrame:
    d = model.derive_demo(demo, cycle_index)
    b = model.derive_bmx(bmx)
    s = model.derive_smoking(smq)
    return d.merge(b, on="SEQN", how="left", validate="one_to_one").merge(
        s, on="SEQN", how="left", validate="one_to_one"
    )


def build_outcome_frame(model, data_dir: Path, variable: str, cycles: list[str]) -> tuple[pd.DataFrame, list[dict]]:
    specs = {str(spec["cycle"]): (idx, spec) for idx, spec in enumerate(model.CYCLES)}
    frames: list[pd.DataFrame] = []
    audit: list[dict] = []
    for cycle in cycles:
        if cycle not in specs:
            raise ValueError(f"Unknown cycle: {cycle}")
        idx, _ = specs[cycle]
        mcq_path = data_dir / f"{cycle}_MCQ.XPT"
        demo_path = data_dir / f"{cycle}_DEMO.XPT"
        bmx_path = data_dir / f"{cycle}_BMX.XPT"
        smq_path = data_dir / f"{cycle}_SMQ.XPT"
        mcq = model.read_xpt(mcq_path)
        if variable not in mcq.columns:
            raise ValueError(f"Outcome variable {variable} missing in {mcq_path.name}")
        core = derive_core(model, model.read_xpt(demo_path), model.read_xpt(bmx_path), model.read_xpt(smq_path), idx)
        frame = mcq[["SEQN", variable]].merge(core, on="SEQN", how="inner", validate="one_to_one")
        value = numeric(frame[variable])
        frame["outcome"] = value.where(value.isin([1, 2])).map({1: 1.0, 2: 0.0})
        frame["cycle"] = cycle
        frame["cycle_index"] = idx
        # model.derive_demo already applies the cycle-specific strata/PSU
        # offsets used by all existing NHANES disease plug-ins.
        frame["adult"] = numeric(frame["age"]).ge(20)
        frames.append(frame)
        audit.append({
            "cycle": cycle,
            "outcome_variable": variable,
            "source_file": mcq_path.name,
            "outcome_source_rows": int(len(mcq)),
            "merged_core_rows": int(len(frame)),
            "adult_rows": int(frame["adult"].sum()),
            "valid_binary_rows": int((frame["adult"] & frame["outcome"].notna()).sum()),
            "case_rows": int((frame["adult"] & frame["outcome"].eq(1)).sum()),
            "control_rows": int((frame["adult"] & frame["outcome"].eq(0)).sum()),
        })
    out = pd.concat(frames, ignore_index=True)
    out = out.loc[out["adult"]].copy()
    out["creatinine_log2"] = np.nan
    return out, audit


def complete_case(population: pd.DataFrame, urine: bool) -> pd.DataFrame:
    required = ["outcome", "axis_log2", "age", "bmi", "pir", "sex", "race", "smoking", "pooled_weight", "psu", "strata"]
    if urine:
        required.append("creatinine_log2")
    complete = population.dropna(subset=required).copy()
    return complete.loc[complete["pooled_weight"].gt(0)].copy()


def fixed_bh(p_values: pd.Series, denominator: int) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    q = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna() & np.isfinite(p)
    if not valid.any():
        return q
    order = p.loc[valid].sort_values().index
    ranked = p.loc[order].to_numpy(float)
    raw = ranked * denominator / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(raw[::-1])[::-1]
    q.loc[order] = np.minimum(adjusted, 1.0)
    return q


def public_fit(fit: dict) -> dict:
    """Expose generic disease-case fields for the cross-disease screen.

    The shared legacy estimator uses ``CRC_N``/``Control_N`` internally. That
    naming is not appropriate once the outcome is a randomly selected disease,
    and it also made successful and non-estimable rows use different schemas.
    Normalize both fields here without changing the estimator or its numbers.
    """
    out = {key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}}
    if "CRC_N" in out:
        out["case_N"] = out.pop("CRC_N")
    if "Control_N" in out:
        out["control_N"] = out.pop("Control_N")
    return out


def run_one_disease(panel_row: pd.Series, tests: pd.DataFrame, registry: pd.DataFrame, model, reader, data_dir: Path, outdir: Path, order: int) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    disease_id = str(panel_row["disease_id"])
    variable = str(panel_row["source_variable(s)"])
    cycles = str(panel_row["cycles_available"]).split(";")
    outcome, outcome_audit = build_outcome_frame(model, data_dir, variable, cycles)
    results: list[dict] = []
    merge_rows: list[dict] = []
    source_rows: list[dict] = []
    for _, test_row in tests.sort_values("test_id").iterrows():
        exposure, source = reader.read_test_exposure(test_row, registry)
        urine = str(test_row["matrix"]).lower() == "urine"
        base = {
            "disease_id": disease_id,
            "disease_name": str(panel_row["disease_name"]),
            "randomization_order": order,
            "test_id": str(test_row["test_id"]),
            "biomarker": str(test_row["biomarker"]),
            "variable": str(test_row["variable"]),
            "matrix": str(test_row["matrix"]),
            "exposure_axis": str(test_row["exposure_axes"]),
            "frozen_cycle_list": str(test_row["cycles"]),
            "frozen_weight": str(test_row["weight"]),
            "source_registry_status": source.get("status"),
            "source_registry_n": int(source.get("n_raw", 0) or 0),
            "source_cycles_read": ";".join(source.get("cycles", [])),
        }
        if exposure.empty:
            results.append({**base, "status": "not_estimable", "reason": source.get("reason", "empty exposure"), "N": 0, "case_N": 0, "control_N": 0, "P": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "analytic_n": 0, "analytic_cases": 0, "analytic_controls": 0})
            continue
        out = outcome.drop(columns=["axis_log2", "pooled_weight"], errors="ignore")
        merged = exposure.merge(out, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
        population = merged.copy()
        population["pooled_weight"] = numeric(population["test_weight"]) / max(len(source.get("cycles", [])), 1)
        population["axis_log2"] = numeric(population["axis_log2"])
        population["creatinine_log2"] = numeric(population.get("creatinine_log2"))
        continuous = ["axis_log2", "age", "bmi", "pir"]
        if urine:
            continuous.append("creatinine_log2")
            # Reuse the same cycle-specific urine creatinine source used by
            # the existing disease plug-ins.
            creat_parts = []
            for cycle in source["cycles"]:
                cfile = data_dir / f"{cycle}_ALB_CR.XPT"
                if not cfile.exists():
                    continue
                c = model.read_xpt(cfile)[["SEQN", "URXUCR"]].copy()
                c["SEQN"] = pd.to_numeric(c["SEQN"], errors="coerce").astype("Int64")
                c["cycle"] = cycle
                c["creatinine_log2"] = np.log2(numeric(c["URXUCR"]).where(numeric(c["URXUCR"]) > 0))
                creat_parts.append(c[["SEQN", "cycle", "creatinine_log2"]])
            if creat_parts:
                population = population.drop(columns=["creatinine_log2"], errors="ignore").merge(pd.concat(creat_parts, ignore_index=True), on=["SEQN", "cycle"], how="left", validate="one_to_one")
        fit = model.fit_survey_logistic(population, continuous, ["sex", "race", "smoking"], exposure_name="axis_log2", levels=model.LEVELS)
        cc = complete_case(population, urine)
        result = {**base, **public_fit(fit), "model_specification": "disease ~ log2(exposure) + age + sex + race + BMI + smoking + PIR" + (" + log2(creatinine)" if urine else ""), "weight_rule": f"analyte-specific laboratory/subsample weight / {len(source.get('cycles', []))} included cycles", "merge_N": int(len(merged)), "analytic_n": int(len(cc)), "analytic_cases": int(cc["outcome"].sum()) if len(cc) else 0, "analytic_controls": int(len(cc) - cc["outcome"].sum()) if len(cc) else 0}
        results.append(result)
        for cycle in source.get("cycles", []):
            e = exposure.loc[exposure["cycle"].eq(cycle)]
            m = merged.loc[merged["cycle"].eq(cycle)]
            p = population.loc[population["cycle"].eq(cycle)]
            pcc = complete_case(p, urine)
            source_row = next((x for x in source.get("source_rows", []) if x.get("cycle") == cycle), {})
            source_rows.append({"disease_id": disease_id, "test_id": str(test_row["test_id"]), "variable": str(test_row["variable"]), "cycle": cycle, **source_row})
            merge_rows.append({"disease_id": disease_id, "test_id": str(test_row["test_id"]), "variable": str(test_row["variable"]), "cycle": cycle, "exposure_rows": int(len(e)), "exposure_nonmissing": int(e["exposure_raw"].notna().sum()), "outcome_rows": int(len(out.loc[out["cycle"].eq(cycle)])), "merge_rows": int(len(m)), "complete_case_rows": int(len(pcc)), "complete_case_cases": int(pcc["outcome"].sum()) if len(pcc) else 0})
    result_df = pd.DataFrame(results).sort_values("test_id").reset_index(drop=True)
    result_df["BH_FDR"] = fixed_bh(result_df["P"], FDR_DENOMINATOR)
    result_df["FDR_supported"] = result_df["BH_FDR"].lt(0.05)
    result_df["fdr_denominator"] = FDR_DENOMINATOR
    finite = result_df["P"].notna() & np.isfinite(pd.to_numeric(result_df["P"], errors="coerce"))
    tech_estimable = result_df["status"].astype(str).isin(["ok", "converged_with_warning"])
    summary = {
        "disease_id": disease_id,
        "disease_name": str(panel_row["disease_name"]),
        "pooled_cases": int(panel_row["case_count_pooled"]),
        "n_cycles": int(panel_row["n_cycles"]),
        "median_analytic_n": float(pd.to_numeric(result_df["analytic_n"], errors="coerce").replace(0, np.nan).median()),
        "median_analytic_cases": float(pd.to_numeric(result_df["analytic_cases"], errors="coerce").replace(0, np.nan).median()),
        "estimable_tests": int(finite.sum()),
        "technical_fit_ok_or_warning": int(tech_estimable.sum()),
        "nominal_positive_tests": int((finite & pd.to_numeric(result_df["P"], errors="coerce").lt(0.05)).sum()),
        "fdr_positive_tests": int(result_df["FDR_supported"].sum()),
        "branch": "Positive" if bool(result_df["FDR_supported"].any()) else "Negative",
        "technical_warning_count": int((result_df["status"].astype(str) != "ok").sum()),
        "replacement_eligible": bool((finite.sum() < 0.5 * FDR_DENOMINATOR) or (result_df["status"].astype(str).eq("fit_failed").any())),
    }
    result_df.to_csv(outdir / f"disease_{order:02d}_results.csv", index=False)
    pd.DataFrame(outcome_audit).to_csv(outdir / f"disease_{order:02d}_outcome_cycle_audit.csv", index=False)
    pd.DataFrame(merge_rows).to_csv(outdir / f"disease_{order:02d}_merge_audit.csv", index=False)
    pd.DataFrame(source_rows).to_csv(outdir / f"disease_{order:02d}_source_manifest.csv", index=False)
    return result_df, pd.DataFrame(outcome_audit), summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", type=Path, default=DEFAULT_TESTS)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--outdir", type=Path, default=STEP10_DIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tests = pd.read_csv(args.tests, dtype=str, keep_default_na=False)
    registry = normalize_registry_paths(pd.read_csv(args.registry, low_memory=False))
    panel = pd.read_csv(args.panel, dtype=str, keep_default_na=False).sort_values("randomization_order")
    if len(tests) != FDR_DENOMINATOR or tests["test_id"].nunique() != FDR_DENOMINATOR:
        raise ValueError("Step 10 requires exactly the frozen 29-test family")
    if panel.empty:
        raise ValueError("Frozen Step 10 panel is empty")

    model = load_statistical_model(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "step10_survey_model")
    reader = load_module(FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py", "step10_exposure_reader")
    model.DATA_DIR = args.data_dir

    summaries = []
    for _, panel_row in panel.iterrows():
        order = int(panel_row["randomization_order"])
        _, _, summary = run_one_disease(panel_row, tests, registry, model, reader, args.data_dir, args.outdir, order)
        summaries.append(summary)
    summary_df = pd.DataFrame(summaries).sort_values("disease_id").reset_index(drop=True)
    summary_df.to_csv(args.outdir / "step10e_cross_disease_replication_summary.csv", index=False)

    generated = datetime.now(timezone.utc).isoformat()
    panel_lock = args.outdir / "STEP10C_RANDOMIZATION_LOCK.json"
    report = [
        "# Step 10E cross-disease replication report",
        "",
        f"Generated (UTC): {generated}",
        "",
        "The frozen 29-test family was screened separately within each disease. BH-FDR denominator was 29 for every disease; branch assignment used only Positive = >=1 BH-FDR<0.05 and Negative = 0 BH-FDR<0.05.",
        "",
        "No LOCO, clustering, GeneCards, CTD-gene convergence, pathways, STRING, transcriptomics, or disease-specific rescue analysis was performed.",
        "",
        "## Replication summary",
        "",
        "| Disease | Pooled cases | Estimable tests | Nominal P<0.05 | FDR<0.05 | Branch | Technical warnings |",
        "|---|---:|---:|---:|---:|---|---:|",
    ]
    for row in summary_df.itertuples(index=False):
        report.append(f"| {row.disease_name} | {row.pooled_cases} | {row.estimable_tests}/29 | {row.nominal_positive_tests} | {row.fdr_positive_tests} | {row.branch} | {row.technical_warning_count} |")
    report += [
        "",
        "## Primary reproducibility metrics",
        "",
        f"- Successfully processed diseases: **{len(summary_df)}/{len(panel)}**.",
        f"- Mean test-family retention: **{summary_df['estimable_tests'].mean() / FDR_DENOMINATOR:.3f}**.",
        f"- Deterministic branch assignment: **{summary_df['branch'].isin(['Positive', 'Negative']).all()}**.",
        "",
        "## Replacement rule audit",
        "",
        "The frozen replacement rule allowed replacement only for a non-reconstructible outcome, <50% technically estimable tests, no compatible survey design, or an invalid outcome definition. No result-driven replacement was performed.",
    ]
    (args.outdir / "STEP10E_CROSS_DISEASE_REPLICATION_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = {
        "lock_type": "STEP10D_CROSS_DISEASE_SCREEN",
        "generated_utc": generated,
        "tests_path": str(args.tests),
        "tests_sha256": sha256(args.tests),
        "registry_path": str(args.registry),
        "registry_sha256": sha256(args.registry),
        "panel_path": str(args.panel),
        "panel_sha256": sha256(args.panel),
        "randomization_lock_sha256": sha256(panel_lock),
        "frozen_test_count": FDR_DENOMINATOR,
        "fdr_denominator_per_disease": FDR_DENOMINATOR,
        "model_specification": "disease ~ log2(exposure) + age + sex + race + BMI + smoking + PIR; urine tests additionally + log2(creatinine)",
        "branch_rule": {"Positive": ">=1 test with BH-FDR < 0.05", "Negative": "0 tests with BH-FDR < 0.05"},
        "replacement_rule": "R1/R2/R3/R4 only; no result-driven replacement",
        "association_results_generated_after_randomization_lock": True,
        "downstream_analysis_performed": False,
        "output_hashes": {name: sha256(args.outdir / name) for name in ["step10e_cross_disease_replication_summary.csv", "STEP10E_CROSS_DISEASE_REPLICATION_REPORT.md"]},
    }
    (args.outdir / "STEP10D_SCREEN_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
