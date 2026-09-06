"""Step 5: uniform CRC screen for the frozen pre-disease test set.

The input test set is the immutable Step 4 unique biomarker table. This
script is the first outcome-aware stage: it reconstructs the NHANES CRC
outcome/covariate frame, fits the same survey-weighted Python estimator to
every frozen analyte, and applies BH-FDR with a fixed denominator of 29.
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
DEFAULT_OUT = FRAMEWORK / "step05_crc_screen"
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


def load_harmonized(model, data_dir: Path) -> tuple[pd.DataFrame, list[dict]]:
    model.DATA_DIR = data_dir
    frames = []
    manifest: list[dict] = []
    for index, spec in enumerate(model.CYCLES):
        frame, source = model.load_cycle(spec, index)
        frames.append(frame)
        manifest.extend(source)
    harmonized = pd.concat(frames, ignore_index=True)
    harmonized = harmonized[(harmonized["age"] >= 20) & harmonized["cancer_outcome_available"]].copy()
    if harmonized.empty:
        raise RuntimeError("The rebuilt NHANES CRC frame is empty")
    harmonized["SEQN"] = pd.to_numeric(harmonized["SEQN"], errors="coerce").astype("Int64")
    return harmonized, manifest


def write_outcome_definition(out: Path, harmonized: pd.DataFrame, manifest: list[dict], model_path: Path) -> None:
    lines = [
        "# CRC outcome definition for Step 5",
        "",
        "This document defines the CRC disease plug-in used only after the Step 4 environmental test set was frozen.",
        "",
        "## Primary population",
        "",
        "- Source: NHANES Medical Conditions Questionnaire (`MCQ`).",
        "- Adult restriction: age >=20 years (`RIDAGEYR`).",
        "- CRC case: `MCQ220=1` and at least one cancer-type code `16` (colon) or `31` (rectum) in the cycle-appropriate cancer-type fields.",
        "- Cancer-free control: `MCQ220=2`.",
        "- Participants reporting a known non-CRC cancer are excluded from the primary CRC-versus-cancer-free comparison.",
        "- Unknown/missing cancer history is excluded from the primary analysis.",
        "- Diagnosis age is retained for later reverse-causation sensitivity analyses; it does not define the primary case/control status.",
        "",
        "## Variables retained in the case/control ledger",
        "",
        "`SEQN`, cycle, age, cancer outcome availability/known status, CRC case, colon case, rectal case, cancer-free status, CRC diagnosis age, and years since CRC diagnosis when available.",
        "",
        "## Primary model",
        "",
        "For every frozen unique NHANES biomarker test: survey-weighted logistic regression of prevalent CRC on log2 biomarker concentration, age, sex, race/ethnicity, BMI, smoking, and poverty-income ratio. Urinary biomarkers additionally include log2 urinary creatinine. Each analyte uses its own laboratory/subsample weight, cycle-specific strata and PSU, and pooled weights divided by the number of included cycles.",
        "",
        f"The BH-FDR denominator is fixed at **{FDR_DENOMINATOR} frozen tests**, including tests whose model is not estimable.",
        "",
        f"Rebuilt adult participants with known cancer history in the harmonized frame: **{len(harmonized):,}**.",
        f"NHANES source records hashed in the Step 5 run manifest: **{len(manifest):,}**.",
        f"Model implementation: `{model_path}`.",
        "",
        "This is a prevalent CRC association screen and does not establish temporality or causality.",
    ]
    (out / "CRC_outcome_definition.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_case_control_ledger(harmonized: pd.DataFrame) -> pd.DataFrame:
    wanted = [
        "SEQN", "cycle", "age", "sex", "race", "crc_case", "colon_case", "rectal_case",
        "both_colon_rectal", "cancer_free", "cancer_known", "cancer_outcome_available",
        "crc_diagnosis_age", "years_since_crc",
    ]
    ledger = harmonized[[c for c in wanted if c in harmonized.columns]].copy()
    ledger["participant_role"] = np.select(
        [ledger["crc_case"].eq(True), ledger["cancer_free"].eq(True)],
        ["CRC_case", "cancer_free_control"],
        default="excluded_known_nonCRC_or_unknown",
    )
    ledger["primary_analysis_eligible"] = ledger["participant_role"].isin(["CRC_case", "cancer_free_control"])
    return ledger.sort_values(["cycle", "SEQN"]).reset_index(drop=True)


def read_test_exposure(test_row: pd.Series, registry: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    variable = str(test_row["variable"])
    rows = registry.loc[registry["variable"].eq(variable)].copy()
    rows = rows.drop_duplicates(["cycle", "data_file", "variable"]).sort_values(["cycle", "data_file"])
    rows = rows.loc[rows["weight_variable"].notna() & rows["weight_variable"].ne("")]
    pieces = []
    source_rows = []
    for row in rows.itertuples(index=False):
        path = Path(row.local_xpt)
        if not path.exists():
            continue
        frame = pd.read_sas(path, format="xport", encoding="latin1")
        flag = str(row.flag_variable) if pd.notna(row.flag_variable) else ""
        needed = ["SEQN", variable, str(row.weight_variable)]
        if flag and flag in frame.columns:
            needed.append(flag)
        if not set(needed).issubset(frame.columns):
            continue
        part = frame[needed].copy()
        part["SEQN"] = pd.to_numeric(part["SEQN"], errors="coerce").astype("Int64")
        part["exposure_raw"] = pd.to_numeric(part[variable], errors="coerce")
        part["test_weight"] = pd.to_numeric(part[str(row.weight_variable)], errors="coerce")
        part["above_lod"] = pd.to_numeric(part[flag], errors="coerce").ne(1) if flag and flag in part.columns else part["exposure_raw"].notna()
        multiplier = 2.0 if str(row.cycle) in {"1999-2000", "2001-2002"} and str(row.weight_variable) == "WTSPH4YR" else 1.0
        part["test_weight"] = part["test_weight"] * multiplier
        part["cycle"] = str(row.cycle)
        pieces.append(part[["SEQN", "cycle", "exposure_raw", "test_weight", "above_lod"]])
        source_rows.append(
            {
                "cycle": str(row.cycle),
                "data_file": str(row.data_file),
                "variable": variable,
                "weight_variable": str(row.weight_variable),
                "weight_multiplier": multiplier,
                "local_xpt": str(path),
            }
        )
    if not pieces:
        return pd.DataFrame(), {"status": "not_estimable", "reason": "no readable weighted analyte file", "cycles": [], "source_rows": []}
    exposure = pd.concat(pieces, ignore_index=True)
    exposure = exposure.sort_values(["cycle", "SEQN"]).drop_duplicates(["cycle", "SEQN"], keep="first")
    cycles = sorted(exposure["cycle"].dropna().unique().tolist())
    exposure["axis_log2"] = np.log2(exposure["exposure_raw"].where(exposure["exposure_raw"] > 0))
    exposure["pooled_weight"] = exposure["test_weight"] / len(cycles)
    return exposure, {"status": "ok", "cycles": cycles, "source_rows": source_rows, "n_raw": len(exposure)}


def public_fit(fit: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}}


def fit_one_test(test_row: pd.Series, registry: pd.DataFrame, harmonized: pd.DataFrame, model) -> tuple[dict[str, object], dict[str, object]]:
    exposure, source = read_test_exposure(test_row, registry)
    base = {
        "test_id": test_row["test_id"],
        "biomarker": test_row["biomarker"],
        "variable": test_row["variable"],
        "matrix": test_row["matrix"],
        "frozen_mapping_count": int(test_row["mapping_count"]),
        "frozen_cycle_list": test_row["cycles"],
        "frozen_weight": test_row["weight"],
        "source_registry_status": source.get("status"),
        "source_registry_n": int(source.get("n_raw", 0) or 0),
        "source_cycles_read": ";".join(source.get("cycles", [])),
    }
    urine = "urine" in str(test_row["matrix"]).lower()
    if exposure.empty:
        result = {**base, "status": "not_estimable", "reason": source.get("reason", "empty exposure"), "N": 0, "CRC_N": 0, "Control_N": 0, "P": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan}
        return result, {**base, "diagnostic_class": "not_estimable", "message": result["reason"]}
    frame = exposure.merge(harmonized, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
    population = model.population_frames(frame)["CRC_vs_cancer_free"]
    continuous = ["axis_log2", "age", "bmi", "pir"]
    if urine:
        continuous.append("creatinine_log2")
    categorical = ["sex", "race", "smoking"]
    fit = model.fit_survey_logistic(population, continuous, categorical, exposure_name="axis_log2", levels=model.LEVELS)
    result = {**base, **public_fit(fit), "model_specification": "CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR" + (" + log2(creatinine)" if urine else ""), "weight_rule": f"analyte-specific weight / {len(source.get('cycles', []))} cycles"}
    required = ["outcome", "axis_log2", "age", "bmi", "pir", *categorical, "pooled_weight", "psu", "strata"]
    if urine:
        required.append("creatinine_log2")
    complete = population.dropna(subset=required).copy()
    complete = complete[complete["pooled_weight"].gt(0)]
    result["analytic_n"] = int(len(complete))
    result["analytic_crc_cases"] = int(complete["outcome"].sum()) if len(complete) else 0
    result["analytic_controls"] = int(len(complete) - complete["outcome"].sum()) if len(complete) else 0
    diagnostic_class = "ok" if fit.get("status") == "ok" else str(fit.get("status", "unknown"))
    diagnostic = {**base, "diagnostic_class": diagnostic_class, "message": fit.get("message", ""), "design_df": fit.get("design_df", np.nan), "PSU_N": fit.get("PSU_N", np.nan), "strata_N": fit.get("strata_N", np.nan), "analytic_n": result["analytic_n"], "analytic_crc_cases": result["analytic_crc_cases"]}
    return result, diagnostic


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
        raise ValueError(f"Frozen Step 4 test set must contain exactly {FDR_DENOMINATOR} unique tests")
    model_path = ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"
    model = load_module(model_path, "step05_survey_model")
    harmonized, source_manifest = load_harmonized(model, args.data_dir)
    ledger = build_case_control_ledger(harmonized)
    ledger.to_csv(args.outdir / "CRC_case_control_ledger.csv", index=False)
    write_outcome_definition(args.outdir, harmonized, source_manifest, model_path)

    results = []
    diagnostics = []
    for _, test_row in tests.sort_values("test_id").iterrows():
        result, diagnostic = fit_one_test(test_row, registry, harmonized, model)
        results.append(result)
        diagnostics.append(diagnostic)
    result_df = pd.DataFrame(results).sort_values("test_id").reset_index(drop=True)
    result_df["BH_FDR"] = fixed_bh(result_df["P"], FDR_DENOMINATOR)
    result_df["FDR_supported"] = result_df["BH_FDR"].lt(0.05)
    result_df["fdr_denominator"] = FDR_DENOMINATOR
    result_df["outcome_stage"] = "CRC outcome-aware"
    result_df.to_csv(args.outdir / "full_29_test_crc_screen.csv", index=False)

    fdr = result_df[["test_id", "biomarker", "variable", "P", "BH_FDR", "FDR_supported", "status", "N", "CRC_N", "OR", "CI_low", "CI_high"]].copy()
    fdr = fdr.sort_values(["BH_FDR", "P", "test_id"], na_position="last").reset_index(drop=True)
    fdr["screen_rank"] = np.arange(1, len(fdr) + 1)
    fdr.to_csv(args.outdir / "crc_bh_fdr_29_tests.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(args.outdir / "model_diagnostics.csv", index=False)
    result_df[["test_id", "biomarker", "variable", "source_registry_n", "N", "CRC_N", "Control_N", "analytic_n", "analytic_crc_cases", "analytic_controls", "status", "source_cycles_read"]].to_csv(args.outdir / "sample_size_by_test.csv", index=False)

    finite = result_df["P"].notna() & np.isfinite(pd.to_numeric(result_df["P"], errors="coerce"))
    supported = result_df.loc[result_df["FDR_supported"]].copy()
    nominal = result_df.loc[finite & result_df["P"].lt(0.05)].copy()
    def one_line(variable: str) -> str:
        row = result_df.loc[result_df["variable"].eq(variable)]
        if row.empty:
            return f"- {variable}: not found in frozen test set"
        r = row.iloc[0]
        return f"- {variable}: OR={r.get('OR', np.nan):.6g}; 95% CI {r.get('CI_low', np.nan):.6g}–{r.get('CI_high', np.nan):.6g}; P={r.get('P', np.nan):.6g}; BH-FDR={r.get('BH_FDR', np.nan):.6g}; status={r.get('status', '')}"
    warning_n = int((result_df["status"].astype(str) != "ok").sum())
    report = [
        "# Step 5 CRC screen report",
        "",
        "## Frozen scope",
        "",
        f"- Frozen unique NHANES tests entered: **{len(tests)}**.",
        f"- Models with finite P values: **{int(finite.sum())}/{len(tests)}**.",
        f"- Nominal P<0.05 tests: **{len(nominal)}**.",
        f"- BH-FDR<0.05 tests with fixed denominator {FDR_DENOMINATOR}: **{len(supported)}**.",
        f"- Model warning/non-ok statuses: **{warning_n}**.",
        "",
        "Primary model: `CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR`; urinary biomarkers additionally include `log2(creatinine)`. Each analyte uses its own laboratory/subsample weight and cycle-pooled strata/PSU.",
        "",
        "## Key prespecified biomarkers",
        "",
        one_line("URXCOP"),
        one_line("LBXPFHS"),
        "",
        "## Nominal and FDR-supported signals",
        "",
        f"Nominal P<0.05: {', '.join(nominal['variable'].tolist()) if len(nominal) else 'none'}.",
        f"BH-FDR<0.05: {', '.join(supported['variable'].tolist()) if len(supported) else 'none'}.",
        "",
        "No test was removed or re-ranked before the 29-test BH-FDR calculation. This screen is an association analysis of prevalent CRC and does not establish temporality or causality.",
    ]
    (args.outdir / "STEP5_CRC_SCREEN_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    lock = {
        "lock_type": "CRC_SCREEN_LOCK",
        "lock_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit_of_frozen_step04_test_set": "f98e353",
        "tests_input": str(args.tests),
        "tests_input_sha256": sha256(args.tests),
        "registry_input": str(args.registry),
        "registry_input_sha256": sha256(args.registry),
        "model_script": str(model_path),
        "model_script_sha256": sha256(model_path),
        "n_frozen_tests": int(len(tests)),
        "fdr_denominator": FDR_DENOMINATOR,
        "model_specification": "CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR; urinary tests additionally + log2(creatinine)",
        "survey_weight_definition": "test-specific laboratory/subsample weight; cycle-specific strata and PSU; pooled weight divided by included cycle count",
        "raw_p_values": {row["variable"]: (None if pd.isna(row["P"]) else float(row["P"])) for row in result_df.to_dict("records")},
        "BH_FDR_values": {row["variable"]: (None if pd.isna(row["BH_FDR"]) else float(row["BH_FDR"])) for row in result_df.to_dict("records")},
        "FDR_supported_signals": supported["variable"].tolist(),
        "status_counts": result_df["status"].value_counts(dropna=False).to_dict(),
    }
    (args.outdir / "CRC_SCREEN_LOCK.json").write_text(json.dumps(lock, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    run_manifest = {
        "analysis": "Step 5 uniform CRC screen for frozen disease-agnostic test set",
        "run_timestamp_utc": lock["lock_timestamp_utc"],
        "n_frozen_tests": int(len(tests)),
        "fdr_denominator": FDR_DENOMINATOR,
        "n_finite_p_values": int(finite.sum()),
        "n_nominal_p_lt_0_05": int(len(nominal)),
        "n_fdr_supported": int(len(supported)),
        "outcome_definition": "MCQ220=1 plus colon/rectal code 16/31 vs MCQ220=2; age >=20",
        "source_manifest_n": len(source_manifest),
        "source_manifest": source_manifest,
        "output_files": ["CRC_outcome_definition.md", "CRC_case_control_ledger.csv", "full_29_test_crc_screen.csv", "crc_bh_fdr_29_tests.csv", "model_diagnostics.csv", "sample_size_by_test.csv", "STEP5_CRC_SCREEN_REPORT.md", "CRC_SCREEN_LOCK.json"],
    }
    (args.outdir / "run_manifest.json").write_text(json.dumps(run_manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(json.dumps({"n_tests": len(tests), "finite_p": int(finite.sum()), "nominal": nominal["variable"].tolist(), "fdr_supported": supported["variable"].tolist(), "mcop": result_df.loc[result_df.variable.eq('URXCOP'), ['OR','CI_low','CI_high','P','BH_FDR','status']].to_dict('records'), "pfhs": result_df.loc[result_df.variable.eq('LBXPFHS'), ['OR','CI_low','CI_high','P','BH_FDR','status']].to_dict('records')}, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
