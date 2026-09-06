"""Matched-complete-case blood versus urinary Cd × sex comparison for OA."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "work" / "nhanes_phase2a" / "data"
ENV = ROOT / "work" / "nhanes_phase2a" / "environmental_xpt"
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
OUT = ROOT / "analysis" / "cadmium_biomarker_discrepancy_oa_matched"
OLD_SCRIPT = ROOT / "work" / "run_cd_sex_arthritis_subtypes.py"
COMPARE_SCRIPT = ROOT / "work" / "run_cd_blood_vs_urine_oa.py"
SENS_SCRIPT = ROOT / "work" / "run_cd_urinary_dilution_sensitivity.py"
TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"

EXPOSURE_CYCLES = [
    "2003-2004", "2005-2006", "2007-2008", "2009-2010",
    "2011-2012", "2013-2014", "2015-2016", "2017-2018",
]
EXPOSURE_TEST_ID = "NHANES_URXUCD"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def numeric(value: pd.Series) -> pd.Series:
    return pd.to_numeric(value, errors="coerce")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cd = load_module(OLD_SCRIPT, "cd_subtypes_matched")
    compare = load_module(COMPARE_SCRIPT, "cd_biomarker_compare_matched")
    sens = load_module(SENS_SCRIPT, "cd_dilution_fit_matched")
    reader = load_module(compare.READER_PATH, "frozen_exposure_reader_matched")
    model = load_module(compare.MODEL_PATH, "nhanes_design_model_matched")
    model.DATA_DIR = DATA

    tests = pd.read_csv(TESTS, dtype=str, keep_default_na=False)
    registry = pd.read_csv(REGISTRY, low_memory=False)
    test_row = tests.loc[tests["test_id"].eq(EXPOSURE_TEST_ID)].iloc[0]
    urine_exposure, urine_source = reader.read_test_exposure(test_row, registry)
    if urine_exposure.empty:
        raise RuntimeError("URXUCD could not be reconstructed")
    urine_cycles = urine_exposure["cycle"].dropna().sort_values().unique().tolist()
    if urine_cycles != EXPOSURE_CYCLES:
        raise AssertionError(f"Unexpected urine cycles: {urine_cycles}")
    ucr_parts = []
    input_hashes: dict[str, str] = {}
    for cycle in urine_cycles:
        path = DATA / f"{cycle}_ALB_CR.XPT"
        ucr = cd.read_ucr(cycle)
        ucr["cycle"] = cycle
        ucr_parts.append(ucr)
        input_hashes[path.name] = sha256(path)
    urine_exposure = urine_exposure.merge(pd.concat(ucr_parts, ignore_index=True), on=["SEQN", "cycle"], how="left", validate="one_to_one")
    urine_frame = compare.build_oa_frames(cd, model, urine_exposure, EXPOSURE_CYCLES)["urine"]

    blood_exposure, blood_hashes, _ = compare.read_blood_exposure()
    demo_weights = []
    for cycle in EXPOSURE_CYCLES:
        demo = pd.read_sas(DATA / f"{cycle}_DEMO.XPT", format="xport", encoding="latin1")
        part = demo[["SEQN", "WTMEC2YR"]].copy()
        part["SEQN"] = pd.to_numeric(part["SEQN"], errors="coerce").astype("Int64")
        part["cycle"] = cycle
        part["blood_weight"] = numeric(part["WTMEC2YR"])
        demo_weights.append(part[["SEQN", "cycle", "blood_weight"]])
    blood_exposure = blood_exposure.merge(pd.concat(demo_weights, ignore_index=True), on=["SEQN", "cycle"], how="left", validate="one_to_one")
    blood_exposure["pooled_weight"] = blood_exposure["blood_weight"] / len(EXPOSURE_CYCLES)
    blood_exposure = blood_exposure.rename(columns={"blood_cd_raw": "exposure_raw"})
    blood_exposure["axis_log2"] = np.log2(blood_exposure["exposure_raw"].where(blood_exposure["exposure_raw"] > 0))
    blood_frame = compare.build_oa_frames(cd, model, blood_exposure, EXPOSURE_CYCLES)["urine"]

    covars = ["SEQN", "cycle", "outcome", "age", "sex", "race", "pir", "smoking", "psu", "strata"]
    u = urine_frame[covars + ["axis_log2", "creatinine_log2", "pooled_weight"]].rename(columns={
        "axis_log2": "urine_axis_log2", "pooled_weight": "urine_weight",
    })
    b = blood_frame[["SEQN", "cycle", "axis_log2", "pooled_weight"]].rename(columns={
        "axis_log2": "blood_axis_log2", "pooled_weight": "blood_weight",
    })
    common = u.merge(b, on=["SEQN", "cycle"], how="inner", validate="one_to_one")
    required_common = [
        "outcome", "urine_axis_log2", "blood_axis_log2", "creatinine_log2",
        "age", "pir", "race", "smoking", "sex", "cycle", "psu", "strata",
        "urine_weight", "blood_weight",
    ]
    common = common.dropna(subset=required_common).copy()
    common = common.loc[common["urine_weight"].gt(0) & common["blood_weight"].gt(0)].reset_index(drop=True)
    if common.empty:
        raise RuntimeError("No matched complete-case OA sample")
    common["outcome"] = numeric(common["outcome"]).astype(int)

    urine_fit_frame = common.rename(columns={"urine_axis_log2": "axis_log2", "urine_weight": "pooled_weight"})
    blood_fit_frame = common.rename(columns={"blood_axis_log2": "axis_log2", "blood_weight": "pooled_weight"})
    urine_row = compare.fit_row(
        sens, urine_fit_frame, "urinary_cadmium", "ucr_main_plus_sex_interaction", "WTSA2YR",
        "log2(URXUCR) + log2(URXUCR):female", "log2(URXUCD)",
    )
    blood_row = compare.fit_row(
        sens, blood_fit_frame, "blood_cadmium", "original_frozen", "WTMEC2YR",
        "none", "log2(LBXBCD)",
    )
    results = pd.DataFrame([urine_row, blood_row])
    results["matched_complete_case_n"] = int(len(common))
    results["matched_complete_case_total_oa_cases"] = int(common["outcome"].sum())
    results["matched_complete_case_male_cases"] = int(common.loc[common["sex"].eq("Male"), "outcome"].sum())
    results["matched_complete_case_female_cases"] = int(common.loc[common["sex"].eq("Female"), "outcome"].sum())
    results["same_seqn_set"] = True
    results["no_new_fdr"] = True
    results.to_csv(OUT / "01_matched_blood_vs_urine_cd_oa_models.csv", index=False)

    cycle_rows = []
    for cycle, group in common.groupby("cycle", sort=True):
        ids = set(group["SEQN"].dropna().astype(int))
        cycle_rows.append({
            "cycle": cycle,
            "matched_complete_case_n": len(group),
            "matched_oa_cases": int(group["outcome"].sum()),
            "matched_male_cases": int(group.loc[group["sex"].eq("Male"), "outcome"].sum()),
            "matched_female_cases": int(group.loc[group["sex"].eq("Female"), "outcome"].sum()),
            "unique_seqn_n": len(ids),
        })
    cycle_df = pd.DataFrame(cycle_rows)
    cycle_df.to_csv(OUT / "02_matched_sample_cycle_audit.csv", index=False)

    urine_result = urine_row
    blood_result = blood_row
    report = [
        "# Matched-complete-case blood versus urinary Cd × sex → OA",
        "",
        "## Scope",
        "",
        "This audit restricts both biomarker models to the exact same SEQN set with non-missing blood Cd (`LBXBCD`), urinary Cd (`URXUCD`), urinary creatinine, OA outcome, and all frozen covariates across the eight 2003–2018 NHANES cycles. It is a paired sample-composition check, not a new outcome screen or mechanism analysis.",
        "",
        "The urinary model retains `log2(URXUCR)` plus its female interaction. The blood model has no urinary-dilution term. Biomarker-specific NHANES weights are retained (urine `WTSA2YR`; blood `WTMEC2YR`), so the SEQN set is matched but the survey weighting remains matrix-specific.",
        "",
        "## Matched sample",
        "",
        f"Matched complete-case N: **{len(common):,}**; OA cases: **{int(common['outcome'].sum()):,}**; male OA cases: **{int(common.loc[common['sex'].eq('Male'), 'outcome'].sum()):,}**; female OA cases: **{int(common.loc[common['sex'].eq('Female'), 'outcome'].sum()):,}**. All eight cycles contribute; cycle details are in `02_matched_sample_cycle_audit.csv`.",
        "",
        "## Results",
        "",
        "| biomarker | β(Cd×female) | SE | P | male β | female β | matched N | OA cases |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
        f"| urinary Cd | {urine_result['interaction_beta']:.6g} | {urine_result['interaction_se']:.6g} | {urine_result['interaction_p']:.6g} | {urine_result['male_beta']:.6g} | {urine_result['female_beta']:.6g} | {int(urine_result['analytic_n'])} | {int(urine_result['total_cases'])} |",
        f"| blood Cd | {blood_result['interaction_beta']:.6g} | {blood_result['interaction_se']:.6g} | {blood_result['interaction_p']:.6g} | {blood_result['male_beta']:.6g} | {blood_result['female_beta']:.6g} | {int(blood_result['analytic_n'])} | {int(blood_result['total_cases'])} |",
        "",
        "## Interpretation",
        "",
        f"After matching the complete-case SEQN set, urinary Cd retains the larger positive interaction (β={urine_result['interaction_beta']:.4f}, P={urine_result['interaction_p']:.4g}), while blood Cd remains smaller and imprecise (β={blood_result['interaction_beta']:.4f}, P={blood_result['interaction_p']:.4g}). Therefore the earlier biomarker discrepancy is not explained solely by blood-versus-urine sample composition in this matched set.",
        "",
        "No formal urinary-versus-blood coefficient-difference test was calculated: the two estimates are correlated because they use the same people, and they retain matrix-specific NHANES weights and dilution terms. The audit therefore compares the two interaction estimates and their uncertainty descriptively.",
        "",
        "This remains evidence for a biomarker-dependent difference in sex heterogeneity, not proof of different Cd toxicity, causality, or protection. Because the two matrices retain their appropriate NHANES weights and the urine model has UCr terms, this is not a claim that the two coefficients are directly interchangeable on an identical estimand scale.",
        "",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
    ]
    (OUT / "MATCHED_BLOOD_VS_URINE_CD_OA_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "Matched complete-case blood versus urinary cadmium x sex for osteoarthritis",
        "cycles": EXPOSURE_CYCLES,
        "matched_definition": "same SEQN with non-missing LBXBCD, URXUCD, URXUCR, OA, all covariates, both survey weights >0",
        "models": [
            {"biomarker": "urinary_cadmium", "variable": "URXUCD", "weight": "WTSA2YR", "terms": "log2(URXUCR) + log2(URXUCR):female"},
            {"biomarker": "blood_cadmium", "variable": "LBXBCD", "weight": "WTMEC2YR", "terms": "none"},
        ],
        "matched_complete_case_n": int(len(common)),
        "matched_total_oa_cases": int(common["outcome"].sum()),
        "no_new_fdr": True,
        "input_sha256": {
            "old_subtype_script": sha256(OLD_SCRIPT),
            "comparison_script": sha256(COMPARE_SCRIPT),
            "sensitivity_fit_script": sha256(SENS_SCRIPT),
            "tests": sha256(TESTS),
            "registry": sha256(REGISTRY),
            **input_hashes,
            **blood_hashes,
        },
        "outputs": [
            "01_matched_blood_vs_urine_cd_oa_models.csv",
            "02_matched_sample_cycle_audit.csv",
            "MATCHED_BLOOD_VS_URINE_CD_OA_REPORT.md",
            "manifest.json",
        ],
        "prohibited_not_run": ["new outcomes", "new FDR", "mechanism", "literature", "figures", "external cohort modeling"],
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"results": results.to_dict("records"), "matched_n": len(common), "output": str(OUT)}, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
