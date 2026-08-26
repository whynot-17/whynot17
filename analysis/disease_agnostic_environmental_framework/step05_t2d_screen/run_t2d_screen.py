"""Assay-specific T2D plug-in screen for the frozen 29-test family.

This is the first disease plug-in after the outcome-blinded environmental
framework.  It deliberately keeps the 29 Step 4 tests unchanged and builds a
T2D outcome/covariate frame independently of every exposure laboratory file.
The primary outcome combines diagnosed diabetes with HbA1c-defined probable
undiagnosed diabetes, while preserving separate QC categories for controls,
indeterminate participants, and likely early-onset insulin-dependent cases.

No CRC, GeneCards, CTD chemical-gene, transcriptomic, or robustness output is
read by this script.
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
DEFAULT_OUT = FRAMEWORK / "step05_t2d_screen"
FDR_DENOMINATOR = 29

GHB_FILES = {
    "1999-2000": ("LAB10.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/LAB10.XPT"),
    "2001-2002": ("L10_B.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/L10_B.XPT"),
    "2003-2004": ("L10_C.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L10_C.XPT"),
    "2005-2006": ("GHB_D.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/GHB_D.XPT"),
    "2007-2008": ("GHB_E.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/GHB_E.XPT"),
    "2009-2010": ("GHB_F.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/GHB_F.XPT"),
    "2011-2012": ("GHB_G.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/GHB_G.XPT"),
    "2013-2014": ("GHB_H.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/GHB_H.XPT"),
    "2015-2016": ("GHB_I.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/GHB_I.XPT"),
    "2017-2018": ("GHB_J.XPT", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/GHB_J.XPT"),
}


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


def clean_code(series: pd.Series, valid: set[int] | None = None) -> pd.Series:
    value = numeric(series)
    if valid is not None:
        value = value.where(value.isin(valid))
    return value


def age_at_diabetes_diagnosis(diq: pd.DataFrame) -> pd.Series:
    candidates = ["DIQ040Q", "DID040Q", "DIQ040G", "DID040"]
    selected = pd.Series(np.nan, index=diq.index, dtype=float)
    for col in candidates:
        if col not in diq.columns:
            continue
        value = numeric(diq[col])
        # 666 is less than one year; 777/888/999 are non-substantive codes.
        value = value.mask(value.eq(666), 0.5)
        value = value.where(value.between(0, 100))
        selected = selected.combine_first(value)
    return selected


def derive_t2d_status(diq: pd.DataFrame, ghb: pd.DataFrame) -> pd.DataFrame:
    if "SEQN" not in diq.columns or "DIQ010" not in diq.columns:
        raise ValueError("DIQ file must contain SEQN and DIQ010")
    diq_code = numeric(diq["DIQ010"])
    insulin = numeric(diq["DIQ050"]) if "DIQ050" in diq.columns else pd.Series(np.nan, index=diq.index)
    diagnosis_age = age_at_diabetes_diagnosis(diq)
    ghb_part = ghb[["SEQN", "LBXGH"]].drop_duplicates("SEQN") if "LBXGH" in ghb.columns else pd.DataFrame({"SEQN": ghb["SEQN"], "LBXGH": np.nan})
    out = pd.DataFrame(
        {
            "SEQN": diq["SEQN"],
            "diq010": diq_code,
            "insulin_current": insulin.where(insulin.isin([1, 2])),
            "diagnosis_age": diagnosis_age,
        }
    ).merge(ghb_part, on="SEQN", how="left", validate="one_to_one")
    out["hba1c_available"] = numeric(out["LBXGH"]).notna()
    out["hba1c_diabetes_threshold"] = numeric(out["LBXGH"]).ge(6.5)
    out["diagnosed_diabetes"] = out["diq010"].eq(1)
    out["likely_type1_excluded"] = out["diagnosed_diabetes"] & out["insulin_current"].eq(1) & out["diagnosis_age"].lt(20)
    # Keep the objective-glycemia extension conservative: only an explicit
    # "No" to DIQ010 can be classified as probable undiagnosed diabetes or a
    # non-diabetic control.  Borderline, unknown, and missing DIQ010 responses
    # remain indeterminate, and a missing HbA1c cannot silently create a
    # control.
    out["probable_undiagnosed_diabetes"] = out["diq010"].eq(2) & out["hba1c_available"] & out["hba1c_diabetes_threshold"]
    out["t2d_case"] = (out["diagnosed_diabetes"] & ~out["likely_type1_excluded"]) | out["probable_undiagnosed_diabetes"]
    out["non_diabetic_control"] = out["diq010"].eq(2) & out["hba1c_available"] & numeric(out["LBXGH"]).lt(6.5)
    out["indeterminate_diabetes"] = ~(out["t2d_case"] | out["non_diabetic_control"] | out["likely_type1_excluded"])
    out["t2d_eligible"] = out["t2d_case"] | out["non_diabetic_control"]
    out["case_source"] = np.select(
        [out["diagnosed_diabetes"] & ~out["likely_type1_excluded"], out["probable_undiagnosed_diabetes"]],
        ["diagnosed_diabetes", "probable_undiagnosed_diabetes"],
        default="",
    )
    out["t2d_status"] = np.select(
        [out["t2d_case"], out["non_diabetic_control"], out["likely_type1_excluded"]],
        ["T2D_case", "non_diabetic_control", "likely_type1_excluded"],
        default="indeterminate",
    )
    return out


def read_t2d_outcome_cycle(model, spec: dict, cycle_index: int, data_dir: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    cycle = str(spec["cycle"])
    paths = {
        "demo": data_dir / f"{cycle}_DEMO.XPT",
        "bmx": data_dir / f"{cycle}_BMX.XPT",
        "smq": data_dir / f"{cycle}_SMQ.XPT",
        "diq": data_dir / f"{cycle}_DIQ.XPT",
        "ghb": data_dir / f"{cycle}_{GHB_FILES[cycle][0]}",
        # Urinary creatinine is a covariate source, not an exposure source.
        # Keep it in the outcome/covariate frame so every urine assay receives
        # the prespecified adjustment without borrowing a column from the
        # assay-specific exposure file.
        "creatinine": data_dir / f"{cycle}_ALB_CR.XPT",
    }
    frames = {key: model.read_xpt(path) for key, path in paths.items()}
    outcome = derive_t2d_status(frames["diq"], frames["ghb"])
    outcome = outcome.merge(model.derive_demo(frames["demo"], cycle_index), on="SEQN", how="inner", validate="one_to_one")
    outcome = outcome.merge(model.derive_bmx(frames["bmx"]), on="SEQN", how="left", validate="one_to_one")
    outcome = outcome.merge(model.derive_smoking(frames["smq"]), on="SEQN", how="left", validate="one_to_one")
    creatinine = frames["creatinine"][["SEQN", "URXUCR"]].drop_duplicates("SEQN")
    outcome = outcome.merge(creatinine, on="SEQN", how="left", validate="one_to_one")
    outcome["SEQN"] = numeric(outcome["SEQN"]).astype("Int64")
    outcome["cycle"] = cycle
    outcome["cycle_index"] = cycle_index
    outcome = outcome.loc[outcome["age"].ge(20)].copy()
    qc = {
        "cycle": cycle,
        "adult_outcome_frame_rows": int(len(outcome)),
        "diq010_available_n": int(outcome["diq010"].notna().sum()),
        "hba1c_available_n": int(outcome["hba1c_available"].sum()),
        "hba1c_threshold_positive_n": int(outcome["hba1c_diabetes_threshold"].sum()),
        "diagnosed_diabetes_n": int(outcome["diagnosed_diabetes"].sum()),
        "probable_undiagnosed_diabetes_n": int(outcome["probable_undiagnosed_diabetes"].sum()),
        "likely_type1_excluded_n": int(outcome["likely_type1_excluded"].sum()),
        "t2d_case_n": int(outcome["t2d_case"].sum()),
        "non_diabetic_control_n": int(outcome["non_diabetic_control"].sum()),
        "indeterminate_n": int(outcome["indeterminate_diabetes"].sum()),
        "t2d_eligible_n": int(outcome["t2d_eligible"].sum()),
        "diq_file": paths["diq"].name,
        "ghb_file": paths["ghb"].name,
        "ghb_variable": "LBXGH",
        "objective_hba1c_rule": "LBXGH >= 6.5%",
        "fasting_glucose_used": False,
    }
    return outcome, qc


def complete_case(population: pd.DataFrame, urine: bool) -> pd.DataFrame:
    required = ["outcome", "axis_log2", "age", "bmi", "pir", "sex", "race", "smoking", "pooled_weight", "psu", "strata"]
    if urine:
        required.append("creatinine_log2")
    complete = population.dropna(subset=required).copy()
    return complete.loc[complete["pooled_weight"].gt(0)].copy()


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


def read_test_exposure(reader, test_row: pd.Series, registry: pd.DataFrame):
    exposure, source = reader.read_test_exposure(test_row, registry)
    return exposure, source


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
    reader_path = FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py"
    model = load_module(model_path, "t2d_survey_model")
    reader = load_module(reader_path, "t2d_exposure_reader")
    model.DATA_DIR = args.data_dir

    specs = {str(spec["cycle"]): (index, spec) for index, spec in enumerate(model.CYCLES)}
    required_cycles = sorted({cycle for value in tests["cycles"] for cycle in str(value).split(";") if cycle})
    missing_ghb = [cycle for cycle in required_cycles if not (args.data_dir / f"{cycle}_{GHB_FILES[cycle][0]}").exists()]
    if missing_ghb:
        raise FileNotFoundError(f"Missing official glycohemoglobin input files: {missing_ghb}")

    outcome_frames: dict[str, pd.DataFrame] = {}
    outcome_qc: list[dict[str, object]] = []
    for cycle in required_cycles:
        index, spec = specs[cycle]
        frame, qc = read_t2d_outcome_cycle(model, spec, index, args.data_dir)
        outcome_frames[cycle] = frame
        outcome_qc.append(qc)
    pd.DataFrame(outcome_qc).sort_values("cycle").to_csv(args.outdir / "t2d_outcome_qc.csv", index=False)

    results: list[dict[str, object]] = []
    merge_audits: list[dict[str, object]] = []
    source_rows: list[dict[str, object]] = []
    for _, test_row in tests.sort_values("test_id").iterrows():
        exposure, source = read_test_exposure(reader, test_row, registry)
        variable = str(test_row["variable"])
        urine = str(test_row["matrix"]).lower() == "urine"
        base = {
            "test_id": str(test_row["test_id"]),
            "biomarker": str(test_row["biomarker"]),
            "variable": variable,
            "exposure_axis": str(test_row["exposure_axes"]),
            "matrix": str(test_row["matrix"]),
            "frozen_mapping_count": int(test_row["mapping_count"]),
            "frozen_cycle_list": str(test_row["cycles"]),
            "frozen_weight": str(test_row["weight"]),
            "source_registry_status": source.get("status"),
            "source_registry_n": int(source.get("n_raw", 0) or 0),
            "source_cycles_read": ";".join(source.get("cycles", [])),
            "source_rows_used": int(len(source.get("source_rows", []))),
        }
        if exposure.empty:
            results.append({
                **base,
                "status": "not_estimable",
                "reason": source.get("reason", "empty exposure"),
                "N": 0,
                "T2D_cases": 0,
                "Control_N": 0,
                "OR_per_doubling": np.nan,
                "CI_low": np.nan,
                "CI_high": np.nan,
                "P": np.nan,
                "analytic_n": 0,
                "analytic_t2d_cases": 0,
                "analytic_controls": 0,
            })
            continue

        outcome = pd.concat([outcome_frames[cycle] for cycle in source["cycles"]], ignore_index=True)
        merged = exposure.merge(outcome, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
        population = merged.loc[merged["t2d_eligible"]].copy()
        population["outcome"] = population["t2d_case"].astype(int)
        continuous = ["axis_log2", "age", "bmi", "pir"]
        if urine:
            population["creatinine_log2"] = np.log2(numeric(population["URXUCR"]).where(numeric(population["URXUCR"]).gt(0)))
            continuous.append("creatinine_log2")
        fit = model.fit_survey_logistic(
            population,
            continuous,
            ["sex", "race", "smoking"],
            outcome="outcome",
            exposure_name="axis_log2",
            levels=model.LEVELS,
        )
        cc = complete_case(population, urine)
        fit_public = public_fit(fit)
        result = {
            **base,
            **fit_public,
            "OR_per_doubling": fit_public.pop("OR", np.nan),
            "model_specification": "T2D ~ log2(exposure) + age + sex + race + BMI + smoking + PIR" + (" + log2(creatinine)" if urine else ""),
            "weight_rule": f"test-specific laboratory/subsample weight / {len(source.get('cycles', []))} included cycles",
            "N": fit_public.pop("N", len(cc)),
            "T2D_cases": fit_public.pop("CRC_N", int(cc["outcome"].sum()) if len(cc) else 0),
            "Control_N": fit_public.pop("Control_N", int(len(cc) - cc["outcome"].sum()) if len(cc) else 0),
            "OR": np.nan,
            "analytic_n": int(len(cc)),
            "analytic_t2d_cases": int(cc["outcome"].sum()) if len(cc) else 0,
            "analytic_controls": int(len(cc) - cc["outcome"].sum()) if len(cc) else 0,
            "merge_N": int(len(merged)),
            "t2d_eligible_merge_N": int(len(population)),
        }
        # Keep a clean public schema while retaining the exact fit status and CI.
        result["OR"] = result["OR_per_doubling"]
        results.append(result)

        for cycle in source["cycles"]:
            e = exposure.loc[exposure["cycle"].eq(cycle)]
            o = outcome.loc[outcome["cycle"].eq(cycle)]
            m = merged.loc[merged["cycle"].eq(cycle)]
            p = population.loc[population["cycle"].eq(cycle)]
            pcc = complete_case(p, urine)
            source_info = next((x for x in source["source_rows"] if x["cycle"] == cycle), {})
            merge_audits.append({
                "test_id": str(test_row["test_id"]),
                "variable": variable,
                "cycle": cycle,
                "source_data_file": source_info.get("data_file", ""),
                "source_weight_variable": source_info.get("weight_variable", ""),
                "exposure_rows": int(len(e)),
                "exposure_nonmissing": int(e["exposure_raw"].notna().sum()),
                "outcome_frame_rows": int(len(o)),
                "merged_rows": int(len(m)),
                "t2d_eligible_rows": int(len(p)),
                "complete_case_rows": int(len(pcc)),
                "complete_case_t2d_cases": int(pcc["outcome"].sum()) if len(pcc) else 0,
            })
            source_rows.append({"test_id": str(test_row["test_id"]), "variable": variable, "cycle": cycle, **source_info})

    result_df = pd.DataFrame(results).sort_values("test_id").reset_index(drop=True)
    result_df["BH_FDR"] = fixed_bh(result_df["P"], FDR_DENOMINATOR)
    result_df["FDR_supported"] = result_df["BH_FDR"].lt(0.05)
    result_df["suggestive_q_lt_0_10"] = result_df["BH_FDR"].lt(0.10)
    result_df["fdr_denominator"] = FDR_DENOMINATOR
    result_df["outcome_definition"] = "Adult T2D case: diagnosed diabetes except likely early-onset insulin-dependent (<20 years at diagnosis), or DIQ010=No with available HbA1c >=6.5%; control: DIQ010=No with available HbA1c <6.5%; borderline, unknown, missing, and other ambiguous responses indeterminate."
    result_df.to_csv(args.outdir / "t2d_primary_29_tests.csv", index=False)
    result_df.sort_values(["P", "BH_FDR", "test_id"], na_position="last").to_csv(args.outdir / "t2d_primary_ranked.csv", index=False)
    pd.DataFrame(merge_audits).to_csv(args.outdir / "t2d_merge_audit.csv", index=False)
    pd.DataFrame(source_rows).to_csv(args.outdir / "t2d_source_manifest.csv", index=False)

    finite = result_df["P"].notna() & np.isfinite(pd.to_numeric(result_df["P"], errors="coerce"))
    nominal = result_df.loc[finite & result_df["P"].lt(0.05)].copy()
    suggestive = result_df.loc[finite & result_df["BH_FDR"].lt(0.10)].copy()
    supported = result_df.loc[finite & result_df["BH_FDR"].lt(0.05)].copy()
    report = [
        "# T2D disease plug-in screen",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Frozen scope",
        "",
        f"- Frozen Step 4 tests entered: **{len(tests)}**.",
        f"- Models with finite P values: **{int(finite.sum())}/{len(tests)}**.",
        f"- Nominal P<0.05: **{len(nominal)}**.",
        f"- Suggestive BH-FDR q<0.10: **{len(suggestive)}**.",
        f"- Primary BH-FDR q<0.05: **{len(supported)}**.",
        "- BH-FDR denominator is fixed at 29; no test was removed or re-ranked before correction.",
        "",
        "## Outcome definition",
        "",
        "Adults aged >=20 years were classified using the Diabetes Questionnaire and official glycohemoglobin files. Diagnosed diabetes was defined by DIQ010=1, with a conservative exclusion for likely early-onset insulin-dependent cases (current insulin use and reported diagnosis age <20). Probable undiagnosed diabetes was restricted to DIQ010=2 with available LBXGH >=6.5%. Controls were restricted to DIQ010=2 with available LBXGH <6.5%; borderline, unknown, missing, and other ambiguous responses remained indeterminate.",
        "",
        "## Primary model",
        "",
        "`T2D ~ log2(exposure) + age + sex + race/ethnicity + BMI + smoking + PIR`; urinary biomarkers additionally include `log2(urinary creatinine)`. Each test uses its own assay-specific laboratory file, cycle coverage and subsample weight. Outcome and covariates are constructed without using any exposure laboratory file.",
        "",
        "## Primary screen result",
        "",
        "| Biomarker | N | T2D cases | OR/doubling | 95% CI | P | q (29) | Status |",
        "| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |",
    ]
    for _, row in result_df.sort_values(["BH_FDR", "P"], na_position="last").iterrows():
        def fmt(value):
            return "NA" if pd.isna(value) else f"{float(value):.6g}"
        report.append(
            f"| {row['biomarker']} | {int(row.get('N', 0) or 0)} | {int(row.get('T2D_cases', 0) or 0)} | {fmt(row.get('OR_per_doubling'))} | {fmt(row.get('CI_low'))}–{fmt(row.get('CI_high'))} | {fmt(row.get('P'))} | {fmt(row.get('BH_FDR'))} | {row.get('status', '')} |"
        )
    report += [
        "",
        "## Interpretation firewall",
        "",
        "This is the first outcome-aware T2D demonstration of the frozen environmental test set. No GeneCards, disease-specific CTD, transcriptomic, literature, robustness, or mechanistic analysis was performed in this stage. The output is a screening result and does not establish causality.",
        "",
        f"Nominal signals: {', '.join(nominal['variable'].tolist()) if len(nominal) else 'none'}.",
        f"Suggestive q<0.10 signals: {', '.join(suggestive['variable'].tolist()) if len(suggestive) else 'none'}.",
        f"q<0.05 signals: {', '.join(supported['variable'].tolist()) if len(supported) else 'none'}.",
    ]
    (args.outdir / "t2d_screen_summary.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    input_files = []
    for cycle in required_cycles:
        ghb_name, ghb_url = GHB_FILES[cycle]
        input_files.append({
            "cycle": cycle,
            "file": ghb_name,
            "url": ghb_url,
            "local_path": str(args.data_dir / f"{cycle}_{ghb_name}"),
            "bytes": (args.data_dir / f"{cycle}_{ghb_name}").stat().st_size,
            "sha256": sha256(args.data_dir / f"{cycle}_{ghb_name}"),
        })
    manifest = {
        "lock_type": "T2D_SCREEN_29_TESTS",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_test_count": FDR_DENOMINATOR,
        "fdr_denominator": FDR_DENOMINATOR,
        "tests_input": str(args.tests),
        "tests_input_sha256": sha256(args.tests),
        "registry_input": str(args.registry),
        "registry_input_sha256": sha256(args.registry),
        "survey_model_script": str(model_path),
        "survey_model_script_sha256": sha256(model_path),
        "exposure_reader_script": str(reader_path),
        "exposure_reader_script_sha256": sha256(reader_path),
        "outcome_definition": "Adults >=20; diagnosed diabetes DIQ010=1 except likely early-onset insulin-dependent cases with DIQ050=1 and diagnosis age <20; probable undiagnosed diabetes if DIQ010=2 and available LBXGH >=6.5%; control if DIQ010=2 and available LBXGH <6.5%; all other ambiguous categories indeterminate.",
        "objective_glycemia": "Official NHANES glycohemoglobin files; LBXGH >=6.5%; fasting glucose not used because it was not part of the local input package.",
        "model_specification": "T2D ~ log2(exposure) + age + sex + race + BMI + smoking + PIR; urinary tests additionally + log2(creatinine)",
        "weight_rule": "test-specific laboratory/subsample weight divided by included assay cycles; cycle-specific PSU and strata retained",
        "n_tests": int(len(result_df)),
        "n_finite_p": int(finite.sum()),
        "n_nominal_p_lt_0_05": int(len(nominal)),
        "n_suggestive_q_lt_0_10": int(len(suggestive)),
        "n_q_lt_0_05": int(len(supported)),
        "source_ghb_files": input_files,
        "output_files": [
            "t2d_primary_29_tests.csv", "t2d_primary_ranked.csv", "t2d_outcome_qc.csv",
            "t2d_merge_audit.csv", "t2d_source_manifest.csv", "t2d_screen_summary.md",
        ],
        "downstream_biology_run": False,
        "robustness_run": False,
        "gene_cards_or_ctd_disease_data_run": False,
    }
    (args.outdir / "t2d_analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "n_tests": len(result_df),
        "finite_p": int(finite.sum()),
        "nominal": nominal["variable"].tolist(),
        "suggestive": suggestive["variable"].tolist(),
        "q_lt_0_05": supported["variable"].tolist(),
    }, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
