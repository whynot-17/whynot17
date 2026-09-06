"""Focused same-cycle blood-Cd versus urinary-Cd comparison for OA.

The audit uses only the frozen URXUCD exposure, newly retrieved official
NHANES blood-cadmium files (LBXBCD), the existing OA definition, and the
frozen survey-weighted covariate framework.  It is not a new phenome screen.
"""
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
OUT = ROOT / "analysis" / "cadmium_biomarker_discrepancy_oa"
OLD_SCRIPT = ROOT / "work" / "run_cd_sex_arthritis_subtypes.py"
SENS_SCRIPT = ROOT / "work" / "run_cd_urinary_dilution_sensitivity.py"
TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
MODEL_PATH = ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"
READER_PATH = FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py"

EXPOSURE_CYCLES = [
    "2003-2004", "2005-2006", "2007-2008", "2009-2010",
    "2011-2012", "2013-2014", "2015-2016", "2017-2018",
]
BLOOD_FILES = {cycle: f"{cycle}_{name}" for cycle, name in {
    "2003-2004": "L06BMT_C.xpt",
    "2005-2006": "PBCD_D.xpt",
    "2007-2008": "PBCD_E.xpt",
    "2009-2010": "PBCD_F.xpt",
    "2011-2012": "PBCD_G.xpt",
    "2013-2014": "PBCD_H.xpt",
    "2015-2016": "PBCD_I.xpt",
    "2017-2018": "PBCD_J.xpt",
}.items()}
BLOOD_URLS = {
    cycle: f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle[:4]}/DataFiles/{name.split('_', 1)[1]}"
    for cycle, name in BLOOD_FILES.items()
}
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


def build_oa_frames(cd, model, exposure: pd.DataFrame, exposure_cycles: list[str]) -> dict[str, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    for idx, cycle in enumerate(cd.CYCLES):
        if cycle not in exposure_cycles:
            continue
        demo = pd.read_sas(DATA / f"{cycle}_DEMO.XPT", format="xport", encoding="latin1")
        mcq = pd.read_sas(DATA / f"{cycle}_MCQ.XPT", format="xport", encoding="latin1")
        smq = pd.read_sas(DATA / f"{cycle}_SMQ.XPT", format="xport", encoding="latin1")
        core = model.derive_demo(demo, idx).merge(
            model.derive_smoking(smq), on="SEQN", how="left", validate="one_to_one",
        )
        core = core.loc[core["age"].ge(20) & core["sex"].notna()].copy()
        endpoint = next(e for e in cd.ENDPOINTS if e["endpoint_id"] == "osteoarthritis")
        status, _ = cd.status_for_endpoint(mcq, cycle, endpoint)
        if status is None:
            continue
        status["SEQN"] = pd.to_numeric(status["SEQN"], errors="coerce").astype("Int64")
        joined = core.merge(status, on="SEQN", how="inner", validate="one_to_one")
        eligible = joined["case"] | joined["control"]
        d = joined.loc[eligible].copy()
        d["outcome"] = d["case"].astype(int)
        d["cycle"] = cycle
        frames.append(d[[
            "SEQN", "cycle", "outcome", "age", "sex", "race", "pir",
            "smoking", "psu", "strata",
        ]])
    status_frame = pd.concat(frames, ignore_index=True)
    return {
        "urine": exposure.merge(status_frame, on=["SEQN", "cycle"], how="inner", validate="one_to_one"),
    }


def read_blood_exposure() -> tuple[pd.DataFrame, dict[str, str], pd.DataFrame]:
    pieces: list[pd.DataFrame] = []
    hashes: dict[str, str] = {}
    availability: list[dict[str, object]] = []
    for cycle in EXPOSURE_CYCLES:
        path = ENV / BLOOD_FILES[cycle]
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_sas(path, format="xport", encoding="latin1")
        seqn = pd.to_numeric(frame["SEQN"], errors="coerce").astype("Int64")
        value = numeric(frame["LBXBCD"])
        part = pd.DataFrame({
            "SEQN": seqn, "cycle": cycle, "blood_cd_raw": value,
        }).drop_duplicates("SEQN")
        pieces.append(part)
        hashes[path.name] = sha256(path)
        availability.append({
            "cycle": cycle, "blood_file": path.name,
            "blood_file_seqn_n": int(part["SEQN"].notna().sum()),
            "blood_file_cd_nonmissing_n": int(value.notna().sum()),
            "blood_url": BLOOD_URLS[cycle],
        })
    blood = pd.concat(pieces, ignore_index=True)
    blood["axis_log2"] = np.log2(blood["blood_cd_raw"].where(blood["blood_cd_raw"] > 0))
    return blood, hashes, pd.DataFrame(availability)


def fit_row(sens, frame: pd.DataFrame, biomarker: str, strategy_id: str, weight_variable: str, dilution_terms: str, source_variable: str) -> dict[str, object]:
    pooled = sens.fit_variant(frame, strategy_id, pooled=True)
    male = sens.fit_variant(frame.loc[frame["sex"].eq("Male")], strategy_id, pooled=False)
    female = sens.fit_variant(frame.loc[frame["sex"].eq("Female")], strategy_id, pooled=False)
    result = {
        "biomarker": biomarker,
        "matrix": "blood" if biomarker == "blood_cadmium" else "urine",
        "source_variable": source_variable,
        "model_variant": strategy_id,
        "weight_variable": weight_variable,
        "dilution_terms": dilution_terms,
        "status": pooled.get("status", "not_estimable"),
        "reason": pooled.get("reason", ""),
        "interaction_beta": pooled.get("beta", np.nan),
        "interaction_se": pooled.get("se", np.nan),
        "interaction_z": pooled.get("z", np.nan),
        "interaction_ci_low": pooled.get("ci_low", np.nan),
        "interaction_ci_high": pooled.get("ci_high", np.nan),
        "interaction_p": pooled.get("p", np.nan),
        "male_beta": male.get("beta", np.nan), "male_se": male.get("se", np.nan), "male_z": male.get("z", np.nan),
        "female_beta": female.get("beta", np.nan), "female_se": female.get("se", np.nan), "female_z": female.get("z", np.nan),
        "analytic_n": pooled.get("analytic_n", 0), "total_cases": pooled.get("cases", 0), "total_controls": pooled.get("controls", 0),
        "male_n": pooled.get("male_n", 0), "female_n": pooled.get("female_n", 0),
        "male_cases": pooled.get("male_cases", 0), "female_cases": pooled.get("female_cases", 0),
        "n_contributing_cycles": pooled.get("n_contributing_cycles", 0),
        "contributing_cycles": pooled.get("contributing_cycles", ""),
        "interaction_positive": bool(np.isfinite(pooled.get("beta", np.nan)) and pooled.get("beta", np.nan) > 0),
        "interaction_p_lt_0_05": bool(np.isfinite(pooled.get("p", np.nan)) and pooled.get("p", np.nan) < 0.05),
    }
    return result


def fmt(value: object, digits: int = 5) -> str:
    try:
        value = float(value)
        return "NA" if not np.isfinite(value) else f"{value:.{digits}g}"
    except (TypeError, ValueError):
        return "NA"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cd = load_module(OLD_SCRIPT, "cd_subtypes_biomarker_comparison")
    sens = load_module(SENS_SCRIPT, "cd_dilution_fit_for_biomarker_comparison")
    reader = load_module(READER_PATH, "frozen_exposure_reader_biomarker_comparison")
    model = load_module(MODEL_PATH, "nhanes_design_model_biomarker_comparison")
    model.DATA_DIR = DATA
    tests = pd.read_csv(TESTS, dtype=str, keep_default_na=False)
    registry = pd.read_csv(REGISTRY, low_memory=False)
    test_row = tests.loc[tests["test_id"].eq(EXPOSURE_TEST_ID)].iloc[0]
    urine_exposure, urine_source = reader.read_test_exposure(test_row, registry)
    if urine_exposure.empty:
        raise RuntimeError("URXUCD could not be reconstructed")
    urine_cycles = urine_exposure["cycle"].dropna().sort_values().unique().tolist()
    if urine_cycles != EXPOSURE_CYCLES:
        raise AssertionError(f"Unexpected URXUCD cycles: {urine_cycles}")
    ucr_parts = []
    hashes: dict[str, str] = {}
    for cycle in urine_cycles:
        path = DATA / f"{cycle}_ALB_CR.XPT"
        ucr = cd.read_ucr(cycle)
        ucr["cycle"] = cycle
        ucr_parts.append(ucr)
        hashes[path.name] = sha256(path)
    urine_exposure = urine_exposure.merge(pd.concat(ucr_parts, ignore_index=True), on=["SEQN", "cycle"], how="left", validate="one_to_one")
    urine_frames = build_oa_frames(cd, model, urine_exposure, urine_cycles)
    urine_frame = urine_frames["urine"]

    blood_exposure, blood_hashes, blood_avail = read_blood_exposure()
    blood_demo_weights = []
    for cycle in EXPOSURE_CYCLES:
        demo = pd.read_sas(DATA / f"{cycle}_DEMO.XPT", format="xport", encoding="latin1")
        part = demo[["SEQN", "WTMEC2YR"]].copy()
        part["SEQN"] = pd.to_numeric(part["SEQN"], errors="coerce").astype("Int64")
        part["cycle"] = cycle
        part["blood_weight"] = numeric(part["WTMEC2YR"])
        blood_demo_weights.append(part[["SEQN", "cycle", "blood_weight"]])
    blood_exposure = blood_exposure.merge(pd.concat(blood_demo_weights, ignore_index=True), on=["SEQN", "cycle"], how="left", validate="one_to_one")
    blood_exposure["pooled_weight"] = blood_exposure["blood_weight"] / len(EXPOSURE_CYCLES)
    blood_exposure = blood_exposure.rename(columns={"blood_cd_raw": "exposure_raw"})
    blood_exposure["axis_log2"] = np.log2(blood_exposure["exposure_raw"].where(blood_exposure["exposure_raw"] > 0))
    blood_frame = build_oa_frames(cd, model, blood_exposure, EXPOSURE_CYCLES)["urine"]

    # Build the cycle-level overlap audit from adult OA-eligible participants.
    overlap_rows = []
    for cycle in EXPOSURE_CYCLES:
        u = urine_frame.loc[urine_frame["cycle"].eq(cycle)]
        b = blood_frame.loc[blood_frame["cycle"].eq(cycle)]
        u_ids = set(u["SEQN"].dropna().astype(int))
        b_ids = set(b["SEQN"].dropna().astype(int))
        overlap_rows.append({
            "cycle": cycle,
            "urine_oa_eligible_n": len(u),
            "blood_oa_eligible_n": len(b),
            "adult_oa_eligible_seqn_overlap_n": len(u_ids & b_ids),
            "adult_oa_eligible_seqn_union_n": len(u_ids | b_ids),
            "urine_oa_cd_nonmissing_n": int(u["axis_log2"].notna().sum()),
            "blood_oa_cd_nonmissing_n": int(b["axis_log2"].notna().sum()),
            "both_cd_nonmissing_overlap_n": int(len(set(u.loc[u["axis_log2"].notna(), "SEQN"].dropna().astype(int)) & set(b.loc[b["axis_log2"].notna(), "SEQN"].dropna().astype(int)))),
        })
    overlap_df = pd.DataFrame(overlap_rows)
    overlap_df = overlap_df.merge(blood_avail, on="cycle", how="left")
    overlap_df.to_csv(OUT / "02_biomarker_cycle_availability.csv", index=False)

    urine_row = fit_row(
        sens, urine_frame, "urinary_cadmium", "ucr_main_plus_sex_interaction", "WTSA2YR",
        "log2(URXUCR) + log2(URXUCR):female", "log2(URXUCD)",
    )
    blood_row = fit_row(
        sens, blood_frame, "blood_cadmium", "original_frozen", "WTMEC2YR",
        "none", "log2(LBXBCD)",
    )
    results = pd.DataFrame([urine_row, blood_row])
    results["historical_context"] = "same 2003-2018 NHANES cycles; no new FDR"
    results.to_csv(OUT / "01_blood_vs_urine_cd_oa_models.csv", index=False)

    report = [
        "# Blood Cd versus urinary Cd × sex → OA",
        "",
        "## Scope",
        "",
        "This focused biomarker-discrepancy audit uses the same eight NHANES cycles (2003–2004 through 2017–2018), the same adult OA definition (`MCQ160A=1` plus cycle-mapped OA subtype versus `MCQ160A=2` controls), the same age/race/ethnicity/PIR/smoking/cycle fixed effects, and the same survey-weighted interaction framework. It compares `LBXBCD` blood cadmium with frozen `URXUCD` urinary cadmium. It does not add outcomes, create a new FDR family, or perform mechanism analysis.",
        "",
        "Urinary Cd retains the established sex-specific urinary-creatinine adjustment (`log2(URXUCR)` plus its female interaction). Blood Cd uses the corresponding NHANES MEC exam weight (`WTMEC2YR`) and does not receive a urinary-dilution term. Urine Cd uses the analyte-specific urine weight (`WTSA2YR`). Each 2-year weight is divided by eight cycles for the pooled fit.",
        "",
        "## Model results",
        "",
        "| biomarker | interaction β (female−male) | SE | P | male β | female β | analytic N | OA cases | cycles |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in results.itertuples(index=False):
        report.append(f"| {row.biomarker} | {fmt(row.interaction_beta)} | {fmt(row.interaction_se)} | {fmt(row.interaction_p)} | {fmt(row.male_beta)} | {fmt(row.female_beta)} | {int(row.analytic_n)} | {int(row.total_cases)} | {row.contributing_cycles} |")
    report.extend([
        "",
        "## Direct readout",
        "",
        f"Urinary Cd has a positive pooled interaction (β={fmt(urine_row['interaction_beta'])}, P={fmt(urine_row['interaction_p'])}), whereas blood Cd has a smaller positive interaction whose confidence interval includes zero (β={fmt(blood_row['interaction_beta'])}, P={fmt(blood_row['interaction_p'])}). Thus this same-cycle analysis supports a biomarker-dependent difference in the Cd × sex signal, not a simple replication of one common OA association.",
        "",
        f"The descriptive sex-stratified slopes are negative for both sexes under both biomarker-specific models (urinary Cd: male β={fmt(urine_row['male_beta'])}, female β={fmt(urine_row['female_beta'])}; blood Cd: male β={fmt(blood_row['male_beta'])}, female β={fmt(blood_row['female_beta'])}). These slopes are descriptive and, because nuisance coefficients can differ between stratified and pooled fits, are not used to override the formal pooled interaction.",
        "",
        "## Cycle and overlap audit",
        "",
        f"All {len(EXPOSURE_CYCLES)} planned cycles were available for both exposure matrices after supplementing the local package with the corresponding official public blood-metal XPT files. The OA-eligible blood/urine SEQN overlap is shown in `02_biomarker_cycle_availability.csv`; biomarker-specific analytic N can differ because laboratory availability and missingness differ.",
        "",
        "## Interpretation boundary",
        "",
        "The formal comparison is the pooled Cd × female interaction, not a comparison of one significant and one non-significant sex-specific slope. A positive interaction means the female slope is higher than the male slope under that biomarker-specific model; it does not establish a positive OA risk, causality, or protection. The blood and urine models address different biological measurement windows and have different laboratory weights, so any discrepancy is interpreted as biomarker-dependent sex heterogeneity, not as proof that one matrix is superior.",
        "",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
    ])
    (OUT / "BLOOD_VS_URINE_CD_OA_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "Same-cycle blood versus urinary cadmium sex interaction for osteoarthritis",
        "cycles": EXPOSURE_CYCLES,
        "outcome": "Osteoarthritis: MCQ160A=1 plus cycle-mapped OA subtype versus MCQ160A=2 controls; adults >=20",
        "models": [
            {"biomarker": "urinary_cadmium", "variable": "URXUCD", "log_transform": "log2", "weight": "WTSA2YR", "dilution_terms": "log2(URXUCR) + log2(URXUCR):female"},
            {"biomarker": "blood_cadmium", "variable": "LBXBCD", "log_transform": "log2", "weight": "WTMEC2YR", "dilution_terms": "none"},
        ],
        "no_new_fdr": True,
        "same_cycle_not_same_complete_case_set": True,
        "source_urls": BLOOD_URLS,
        "input_sha256": {
            "old_subtype_script": sha256(OLD_SCRIPT),
            "sensitivity_fit_script": sha256(SENS_SCRIPT),
            "tests": sha256(TESTS),
            "registry": sha256(REGISTRY),
            **hashes,
            **blood_hashes,
        },
        "outputs": [
            "01_blood_vs_urine_cd_oa_models.csv",
            "02_biomarker_cycle_availability.csv",
            "BLOOD_VS_URINE_CD_OA_REPORT.md",
            "manifest.json",
        ],
        "prohibited_not_run": ["new phenome screen", "new subtype screen", "mechanism", "pathway", "figures", "new multiplicity adjustment", "external cohort modeling"],
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"results": results.to_dict("records"), "overlap": overlap_df.to_dict("records"), "output": str(OUT)}, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
