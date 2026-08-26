"""Assay-specific rebuild of the disease-aware NHANES CRC screen.

The provisional Step 5 implementation used the MBzP/phthalate laboratory
sample to build the participant frame for every test.  This rebuild keeps the
frozen, pre-disease 29-test family but reads each test's own laboratory file,
weight and cycle coverage, then merges it to a cycle-matched NHANES outcome
and covariate frame built independently of any exposure laboratory file.

The old provisional outputs are never overwritten by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
DEFAULT_TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
DEFAULT_REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_OUT = FRAMEWORK / "step05_crc_screen_rebuilt"
FDR_DENOMINATOR = 29


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def read_universal_outcome_cycle(model, spec: dict, cycle_index: int, data_dir: Path) -> pd.DataFrame:
    """Build cycle outcome/covariates without any exposure lab file."""
    cycle = str(spec["cycle"])
    paths = {
        "mcq": data_dir / f"{cycle}_MCQ.XPT",
        "demo": data_dir / f"{cycle}_DEMO.XPT",
        "bmx": data_dir / f"{cycle}_BMX.XPT",
        "smq": data_dir / f"{cycle}_SMQ.XPT",
        "alq": data_dir / f"{cycle}_ALQ.XPT",
        "diq": data_dir / f"{cycle}_DIQ.XPT",
        "paq": data_dir / f"{cycle}_PAQ.XPT",
        "creatinine": data_dir / f"{cycle}_ALB_CR.XPT",
    }
    frames = {key: model.read_xpt(path) for key, path in paths.items()}
    merged = model.cancer_flags(frames["mcq"], frames["demo"])
    merged = merged.merge(model.derive_demo(frames["demo"], cycle_index), on="SEQN", how="inner", validate="one_to_one")
    for supplement in (
        model.derive_bmx(frames["bmx"]),
        model.derive_smoking(frames["smq"]),
        model.derive_alcohol(frames["alq"]),
        model.derive_diabetes(frames["diq"]),
        model.derive_activity(frames["paq"]),
        frames["creatinine"][["SEQN", "URXUCR"]].drop_duplicates("SEQN"),
    ):
        merged = merged.merge(supplement, on="SEQN", how="left", validate="one_to_one")
    merged["SEQN"] = numeric(merged["SEQN"]).astype("Int64")
    merged["cycle"] = cycle
    merged["cycle_index"] = cycle_index
    merged["creatinine_log2"] = np.log2(numeric(merged["URXUCR"]).where(numeric(merged["URXUCR"]) > 0))
    merged["years_since_crc"] = merged["age"] - merged["crc_diagnosis_age"]
    return merged


def complete_case(population: pd.DataFrame, urine: bool) -> pd.DataFrame:
    required = ["outcome", "axis_log2", "age", "bmi", "pir", "sex", "race", "smoking", "pooled_weight", "psu", "strata"]
    if urine:
        required.append("creatinine_log2")
    complete = population.dropna(subset=required).copy()
    return complete[complete["pooled_weight"].gt(0)]


def public_fit(fit: dict) -> dict:
    return {key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}}


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


def fit_test(test_row: pd.Series, registry: pd.DataFrame, outcome_frames: dict[str, pd.DataFrame], model):
    exposure, source = load_module(FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py", "step5_reader").read_test_exposure(test_row, registry)
    variable = str(test_row["variable"])
    urine = str(test_row["matrix"]).lower() == "urine"
    base = {
        "test_id": str(test_row["test_id"]),
        "biomarker": str(test_row["biomarker"]),
        "variable": variable,
        "matrix": str(test_row["matrix"]),
        "exposure_axis": str(test_row["exposure_axes"]),
        "frozen_mapping_count": int(test_row["mapping_count"]),
        "frozen_cycle_list": str(test_row["cycles"]),
        "frozen_weight": str(test_row["weight"]),
        "source_registry_status": source.get("status"),
        "source_registry_n": int(source.get("n_raw", 0) or 0),
        "source_cycles_read": ";".join(source.get("cycles", [])),
        "source_rows_used": int(len(source.get("source_rows", []))),
    }
    if exposure.empty:
        result = {**base, "status": "not_estimable", "reason": source.get("reason", "empty exposure"), "N": 0, "CRC_N": 0, "Control_N": 0, "P": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "merge_N": 0, "analytic_n": 0, "analytic_crc_cases": 0, "analytic_controls": 0}
        return result, pd.DataFrame(), source

    outcome = pd.concat([outcome_frames[cycle] for cycle in source["cycles"]], ignore_index=True)
    merged = exposure.merge(outcome, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
    population = model.population_frames(merged)["CRC_vs_cancer_free"].copy()
    population["pooled_weight"] = population["test_weight"] / len(source["cycles"])
    continuous = ["axis_log2", "age", "bmi", "pir"]
    if urine:
        continuous.append("creatinine_log2")
    categorical = ["sex", "race", "smoking"]
    fit = model.fit_survey_logistic(population, continuous, categorical, exposure_name="axis_log2", levels=model.LEVELS)
    cc = complete_case(population, urine)
    result = {
        **base,
        **public_fit(fit),
        "model_specification": "CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR" + (" + log2(creatinine)" if urine else ""),
        "weight_rule": f"test-specific laboratory/subsample weight / {len(source.get('cycles', []))} included cycles",
        "merge_N": int(len(merged)),
        "analytic_n": int(len(cc)),
        "analytic_crc_cases": int(cc["outcome"].sum()) if len(cc) else 0,
        "analytic_controls": int(len(cc) - cc["outcome"].sum()) if len(cc) else 0,
    }
    cycle_rows = []
    for cycle in source["cycles"]:
        e = exposure.loc[exposure["cycle"].eq(cycle)]
        o = outcome.loc[outcome["cycle"].eq(cycle)]
        m = merged.loc[merged["cycle"].eq(cycle)]
        p = population.loc[population["cycle"].eq(cycle)]
        cycle_rows.append({
            "test_id": str(test_row["test_id"]),
            "variable": variable,
            "cycle": cycle,
            "source_data_file": next((x["data_file"] for x in source["source_rows"] if x["cycle"] == cycle), ""),
            "source_weight_variable": next((x["weight_variable"] for x in source["source_rows"] if x["cycle"] == cycle), ""),
            "exposure_rows": int(len(e)),
            "exposure_nonmissing": int(e["exposure_raw"].notna().sum()),
            "universal_outcome_covariate_rows": int(len(o)),
            "merge_rows": int(len(m)),
            "crc_vs_cancer_free_rows": int(len(p)),
            "complete_case_rows": int(len(complete_case(p, urine))),
            "complete_case_crc_cases": int(complete_case(p, urine)["outcome"].sum()) if len(complete_case(p, urine)) else 0,
        })
    return result, pd.DataFrame(cycle_rows), source


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", type=Path, default=DEFAULT_TESTS)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tests = pd.read_csv(args.tests, dtype=str, keep_default_na=False)
    registry = pd.read_csv(args.registry, low_memory=False)
    if len(tests) != FDR_DENOMINATOR or tests["test_id"].nunique() != FDR_DENOMINATOR:
        raise ValueError("The rebuilt Step 5 requires exactly the frozen 29-test set")
    model = load_module(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "step5_rebuilt_model")
    step5_reader = load_module(FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py", "step5_rebuilt_reader")
    model.DATA_DIR = args.data_dir

    # Build one outcome/covariate frame per NHANES cycle, independently of all
    # exposure laboratory files.  The exposure reader then supplies the
    # assay-specific participant subset and survey weight for each test.
    specs = {str(spec["cycle"]): (index, spec) for index, spec in enumerate(model.CYCLES)}
    required_cycles = sorted(set(cycle for value in tests["cycles"] for cycle in str(value).split(";") if cycle))
    outcome_frames = {}
    outcome_qc = []
    for cycle in required_cycles:
        index, spec = specs[cycle]
        frame = read_universal_outcome_cycle(model, spec, index, args.data_dir)
        outcome_frames[cycle] = frame
        pop = model.population_frames(frame)["CRC_vs_cancer_free"]
        outcome_qc.append({
            "cycle": cycle,
            "universal_frame_rows": int(len(frame)),
            "crc_vs_cancer_free_rows": int(len(pop)),
            "crc_cases": int(pop["outcome"].sum()) if len(pop) else 0,
            "cancer_free_controls": int(len(pop) - pop["outcome"].sum()) if len(pop) else 0,
            "n_psu": int(frame["psu"].nunique()),
            "n_strata": int(frame["strata"].nunique()),
        })

    results, merge_audits, source_rows = [], [], []
    for _, test_row in tests.sort_values("test_id").iterrows():
        # Keep the reader module fixed and local to the provisional Step 5
        # code; only its exposure-file function is reused here.
        exposure, source = step5_reader.read_test_exposure(test_row, registry)
        variable = str(test_row["variable"])
        urine = str(test_row["matrix"]).lower() == "urine"
        base = {
            "test_id": str(test_row["test_id"]),
            "biomarker": str(test_row["biomarker"]),
            "variable": variable,
            "matrix": str(test_row["matrix"]),
            "exposure_axis": str(test_row["exposure_axes"]),
            "frozen_mapping_count": int(test_row["mapping_count"]),
            "frozen_cycle_list": str(test_row["cycles"]),
            "frozen_weight": str(test_row["weight"]),
            "source_registry_status": source.get("status"),
            "source_registry_n": int(source.get("n_raw", 0) or 0),
            "source_cycles_read": ";".join(source.get("cycles", [])),
            "source_rows_used": int(len(source.get("source_rows", []))),
        }
        if exposure.empty:
            result = {**base, "status": "not_estimable", "reason": source.get("reason", "empty exposure"), "N": 0, "CRC_N": 0, "Control_N": 0, "P": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "merge_N": 0, "analytic_n": 0, "analytic_crc_cases": 0, "analytic_controls": 0}
            results.append(result)
            continue
        outcome = pd.concat([outcome_frames[cycle] for cycle in source["cycles"]], ignore_index=True)
        merged = exposure.merge(outcome, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
        population = model.population_frames(merged)["CRC_vs_cancer_free"].copy()
        population["pooled_weight"] = population["test_weight"] / len(source["cycles"])
        continuous = ["axis_log2", "age", "bmi", "pir"]
        if urine:
            continuous.append("creatinine_log2")
        fit = model.fit_survey_logistic(population, continuous, ["sex", "race", "smoking"], exposure_name="axis_log2", levels=model.LEVELS)
        cc = complete_case(population, urine)
        result = {
            **base,
            **public_fit(fit),
            "model_specification": "CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR" + (" + log2(creatinine)" if urine else ""),
            "weight_rule": f"test-specific laboratory/subsample weight / {len(source.get('cycles', []))} included cycles",
            "merge_N": int(len(merged)),
            "analytic_n": int(len(cc)),
            "analytic_crc_cases": int(cc["outcome"].sum()) if len(cc) else 0,
            "analytic_controls": int(len(cc) - cc["outcome"].sum()) if len(cc) else 0,
        }
        results.append(result)
        for cycle in source["cycles"]:
            e = exposure.loc[exposure["cycle"].eq(cycle)]
            o = outcome.loc[outcome["cycle"].eq(cycle)]
            m = merged.loc[merged["cycle"].eq(cycle)]
            p = population.loc[population["cycle"].eq(cycle)]
            pcc = complete_case(p, urine)
            merge_audits.append({
                "test_id": str(test_row["test_id"]),
                "variable": variable,
                "cycle": cycle,
                "source_data_file": next((x["data_file"] for x in source["source_rows"] if x["cycle"] == cycle), ""),
                "source_weight_variable": next((x["weight_variable"] for x in source["source_rows"] if x["cycle"] == cycle), ""),
                "exposure_rows": int(len(e)),
                "exposure_nonmissing": int(e["exposure_raw"].notna().sum()),
                "universal_outcome_covariate_rows": int(len(o)),
                "merge_rows": int(len(m)),
                "crc_vs_cancer_free_rows": int(len(p)),
                "complete_case_rows": int(len(pcc)),
                "complete_case_crc_cases": int(pcc["outcome"].sum()) if len(pcc) else 0,
            })
            source_rows.append({
                "test_id": str(test_row["test_id"]),
                "variable": variable,
                "cycle": cycle,
                **next((x for x in source["source_rows"] if x["cycle"] == cycle), {}),
            })

    result_df = pd.DataFrame(results).sort_values("test_id").reset_index(drop=True)
    result_df["BH_FDR"] = fixed_bh(result_df["P"], FDR_DENOMINATOR)
    result_df["FDR_supported"] = result_df["BH_FDR"].lt(0.05)
    result_df["fdr_denominator"] = FDR_DENOMINATOR
    result_df["outcome_stage"] = "CRC outcome-aware; assay-specific rebuild"
    result_df.to_csv(args.outdir / "full_29_test_crc_screen_rebuilt.csv", index=False)
    result_df[["test_id", "biomarker", "variable", "P", "BH_FDR", "FDR_supported", "status", "N", "CRC_N", "OR", "CI_low", "CI_high"]].sort_values(["BH_FDR", "P", "test_id"], na_position="last").to_csv(args.outdir / "crc_bh_fdr_29_tests_rebuilt.csv", index=False)
    result_df[["test_id", "biomarker", "variable", "source_registry_n", "source_cycles_read", "merge_N", "N", "CRC_N", "Control_N", "analytic_n", "analytic_crc_cases", "analytic_controls", "status"]].to_csv(args.outdir / "sample_size_by_test_rebuilt.csv", index=False)
    pd.DataFrame(outcome_qc).to_csv(args.outdir / "assay_specific_outcome_frame_qc.csv", index=False)
    pd.DataFrame(merge_audits).to_csv(args.outdir / "assay_specific_merge_audit.csv", index=False)
    pd.DataFrame(source_rows).to_csv(args.outdir / "assay_specific_source_manifest.csv", index=False)

    finite_p = result_df["P"].notna() & np.isfinite(pd.to_numeric(result_df["P"], errors="coerce"))
    supported = result_df.loc[result_df["FDR_supported"]].copy()
    nominal = result_df.loc[finite_p & result_df["P"].lt(0.05)].copy()
    report = [
        "# Step 5 assay-specific rebuild report",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Why this rebuild exists",
        "",
        "The provisional Step 5 implementation used the MBzP/phthalate laboratory sample as the participant frame for all 29 tests. This rebuild keeps the frozen pre-disease test family but reads each test's own laboratory file, analyte variable, survey weight and cycle coverage. Outcome and covariates are built independently from cycle-matched NHANES MCQ/DEMO/BMX/SMQ/ALQ/DIQ/PAQ/ALB_CR files.",
        "",
        "The provisional Step 5 outputs are retained unchanged for audit and are not overwritten.",
        "",
        "## Screen summary",
        "",
        f"- Frozen tests entered: **{len(tests)}**.",
        f"- Models with finite P values: **{int(finite_p.sum())}/{len(tests)}**.",
        f"- Nominal P<0.05: **{len(nominal)}**.",
        f"- BH-FDR<0.05 using the frozen denominator {FDR_DENOMINATOR}: **{len(supported)}**.",
        "",
        "| Biomarker | N | CRC cases | OR | 95% CI | P | BH-FDR (29) | Status |",
        "|---|---:|---:|---:|---|---:|---:|---|",
    ]
    for _, row in result_df.sort_values(["BH_FDR", "P"], na_position="last").iterrows():
        fmt = lambda x: "NA" if pd.isna(x) else f"{float(x):.6g}"
        report.append(f"| {row['biomarker']} | {int(row.get('N', 0) or 0)} | {int(row.get('CRC_N', 0) or 0)} | {fmt(row.get('OR'))} | {fmt(row.get('CI_low'))}–{fmt(row.get('CI_high'))} | {fmt(row.get('P'))} | {fmt(row.get('BH_FDR'))} | {row.get('status', '')} |")
    report += [
        "",
        "## URXP25 resolution",
        "",
        "URXP25 now uses PAH_H/PAH_I/PAH_J directly. Its previous N=0 was caused by zero overlap with the provisional phthalate-shaped frame, not by absent exposure values. The corrected isolated model is retained in the same rebuild output and is not hand-selected or excluded.",
        "",
        "## Multiple-testing rule",
        "",
        f"BH-FDR is recomputed once across all **{FDR_DENOMINATOR}** frozen tests, including any test that remains not estimable after the assay-specific rebuild. No result is removed before FDR, and no test is added after seeing CRC outcomes.",
        "",
        "## QC and interpretation firewall",
        "",
        "For each test, source rows, cycle coverage, exposure-to-outcome merge counts, analytic N, CRC cases, survey PSU/strata and fit status are written to the accompanying QC tables. This screen remains a cross-sectional prevalent-CRC analysis and does not establish causality or temporality.",
    ]
    (args.outdir / "STEP5_REBUILT_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = {
        "lock_type": "STEP5_ASSAY_SPECIFIC_REBUILD",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_test_count": FDR_DENOMINATOR,
        "fdr_denominator": FDR_DENOMINATOR,
        "test_set": str(args.tests),
        "test_set_sha256": sha256(args.tests),
        "registry": str(args.registry),
        "registry_sha256": sha256(args.registry),
        "model_script": str(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"),
        "model_script_sha256": sha256(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"),
        "outcome_frame_definition": "cycle-matched MCQ/DEMO/BMX/SMQ/ALQ/DIQ/PAQ/ALB_CR; no exposure laboratory file used to construct outcome frame",
        "exposure_definition": "test-specific registry local_xpt + variable + flag_variable + weight_variable; pooled weight divided by the test's included cycle count",
        "model_specification": "CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR; urine tests additionally + log2(creatinine)",
        "old_provisional_step5_preserved": True,
        "results": result_df[["variable", "N", "CRC_N", "OR", "P", "BH_FDR", "FDR_supported", "status"]].to_dict(orient="records"),
        "output_hashes": {name: sha256(args.outdir / name) for name in ["full_29_test_crc_screen_rebuilt.csv", "crc_bh_fdr_29_tests_rebuilt.csv", "sample_size_by_test_rebuilt.csv", "assay_specific_outcome_frame_qc.csv", "assay_specific_merge_audit.csv", "assay_specific_source_manifest.csv", "STEP5_REBUILT_REPORT.md"]},
    }
    (args.outdir / "STEP5_REBUILT_LOCK.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
