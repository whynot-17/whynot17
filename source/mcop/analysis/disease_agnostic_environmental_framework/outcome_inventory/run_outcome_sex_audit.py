"""Create a blinded, outcome-only NHANES sex-specific feasibility audit.

This program has an intentionally small input surface.  It can read only
DEMO plus the outcome modules listed below; it never opens an exposure result,
an exposure laboratory file, a Step 5 result, or a CRC/T2D association file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_OUT_DIR = Path(__file__).resolve().parent
CYCLES = [
    "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
    "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018",
]
GHB_FILES = {
    "1999-2000": "LAB10.XPT", "2001-2002": "L10_B.XPT", "2003-2004": "L10_C.XPT",
    "2005-2006": "GHB_D.XPT", "2007-2008": "GHB_E.XPT", "2009-2010": "GHB_F.XPT",
    "2011-2012": "GHB_G.XPT", "2013-2014": "GHB_H.XPT", "2015-2016": "GHB_I.XPT",
    "2017-2018": "GHB_J.XPT",
}

# Every candidate and definition is fixed before counts are computed.  The
# unavailable entries remain in this registry so the frozen panel has exactly
# the requested 18 candidates rather than a post-hoc selected subset.
REGISTRY = [
    ("T2D", "Type 2 diabetes", "DIQ; HbA1c", "DIQ010; DIQ050; diagnosis-age field; LBXGH",
     "Age >=20. Case: diagnosed diabetes excluding current-insulin cases diagnosed before 20, or DIQ010=No with HbA1c >=6.5%. Control: DIQ010=No with HbA1c <6.5%.", "t2d"),
    ("obesity", "Obesity", "BMX", "BMXBMI", "Age >=20. Case: BMI >=30 kg/m2. Control: BMI <30 kg/m2.", "bmi"),
    ("hypertension", "Hypertension", "BPQ", "BPQ020", "Age >=20. Case: BPQ020=Yes. Control: BPQ020=No.", "binary:BPQ:BPQ020"),
    ("congestive_heart_failure", "Congestive heart failure", "MCQ", "MCQ160B", "Age >=20. Case: MCQ160B=Yes. Control: MCQ160B=No.", "binary:MCQ:MCQ160B"),
    ("coronary_heart_disease", "Coronary heart disease", "MCQ", "MCQ160C", "Age >=20. Case: MCQ160C=Yes. Control: MCQ160C=No.", "binary:MCQ:MCQ160C"),
    ("myocardial_infarction", "Myocardial infarction", "MCQ", "MCQ160E", "Age >=20. Case: MCQ160E=Yes. Control: MCQ160E=No.", "binary:MCQ:MCQ160E"),
    ("stroke", "Stroke", "MCQ", "MCQ160F", "Age >=20. Case: MCQ160F=Yes. Control: MCQ160F=No.", "binary:MCQ:MCQ160F"),
    ("asthma", "Asthma", "MCQ", "MCQ010", "Age >=20. Case: MCQ010=Yes. Control: MCQ010=No.", "binary:MCQ:MCQ010"),
    ("chronic_bronchitis", "Chronic bronchitis", "MCQ", "MCQ160K", "Age >=20. Case: MCQ160K=Yes. Control: MCQ160K=No.", "binary:MCQ:MCQ160K"),
    ("emphysema", "Emphysema", "MCQ", "MCQ160G", "Age >=20. Case: MCQ160G=Yes. Control: MCQ160G=No.", "binary:MCQ:MCQ160G"),
    ("objective_CKD", "Objective CKD", "serum creatinine laboratory module", "serum creatinine/eGFR", "Age >=20. Prespecified requirement: objective eGFR-based CKD definition. Not substituted with albuminuria alone.", "unavailable:serum creatinine laboratory module absent locally"),
    ("kidney_stones", "Kidney stones", "kidney conditions questionnaire", "KIQ kidney-stone item", "Age >=20. Physician/history kidney-stone item required.", "unavailable:kidney conditions questionnaire absent locally"),
    ("liver_disease", "Liver disease", "MCQ", "MCQ160L", "Age >=20. Case: MCQ160L=Yes. Control: MCQ160L=No.", "binary:MCQ:MCQ160L"),
    ("thyroid_disease", "Thyroid disease", "MCQ", "MCQ160M", "Age >=20. Case: MCQ160M=Yes. Control: MCQ160M=No.", "binary:MCQ:MCQ160M"),
    ("arthritis", "Arthritis", "MCQ", "MCQ160A", "Age >=20. Case: MCQ160A=Yes. Control: MCQ160A=No.", "binary:MCQ:MCQ160A"),
    ("osteoporosis", "Osteoporosis", "osteoporosis questionnaire", "OSQ osteoporosis item", "Age >=20. Osteoporosis questionnaire item required.", "unavailable:osteoporosis questionnaire absent locally"),
    ("clinically_relevant_depressive_symptoms", "Clinically relevant depressive symptoms", "DPQ", "PHQ-9 items DPQ010-DPQ090", "Age >=20. PHQ-9 total >=10 vs <10 with complete PHQ-9.", "unavailable:DPQ depression questionnaire absent locally"),
    ("any_cancer_history", "Any cancer history", "MCQ", "MCQ220", "Age >=20. Case: MCQ220=Yes. Control: MCQ220=No.", "binary:MCQ:MCQ220"),
]


def num(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce")


def read_xpt(path: Path) -> pd.DataFrame:
    return pd.read_sas(path, format="xport", encoding="latin1")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def core_demo(data_dir: Path, cycle: str) -> tuple[pd.DataFrame | None, str]:
    path = data_dir / f"{cycle}_DEMO.XPT"
    if not path.exists():
        return None, f"Missing demographics file: {path.name}"
    demo = read_xpt(path)
    required = {"SEQN", "RIDAGEYR", "RIAGENDR"}
    if not required.issubset(demo.columns):
        return None, f"Demographics lacks required columns: {', '.join(sorted(required - set(demo.columns)))}"
    d = demo[["SEQN", "RIDAGEYR", "RIAGENDR"]].copy()
    d["age"] = num(d["RIDAGEYR"])
    d["sex"] = num(d["RIAGENDR"]).map({1: "male", 2: "female"})
    return d.loc[d["age"].ge(20) & d["sex"].notna(), ["SEQN", "sex"]], ""


def t2d_status(data_dir: Path, cycle: str) -> tuple[pd.DataFrame | None, str, str]:
    diq_path = data_dir / f"{cycle}_DIQ.XPT"
    ghb_path = data_dir / f"{cycle}_{GHB_FILES[cycle]}"
    if not diq_path.exists() or not ghb_path.exists():
        return None, "", "Missing DIQ or HbA1c source file"
    diq, ghb = read_xpt(diq_path), read_xpt(ghb_path)
    if "DIQ010" not in diq or "LBXGH" not in ghb:
        return None, "", "Missing DIQ010 or LBXGH"
    diagnosis_age = pd.Series(np.nan, index=diq.index, dtype=float)
    for col in ["DIQ040Q", "DID040Q", "DIQ040G", "DID040"]:
        if col in diq:
            x = num(diq[col]).mask(num(diq[col]).eq(666), 0.5).where(num(diq[col]).between(0, 100))
            diagnosis_age = diagnosis_age.combine_first(x)
    insulin = num(diq["DIQ050"]) if "DIQ050" in diq else pd.Series(np.nan, index=diq.index)
    out = diq[["SEQN"]].copy()
    out["diq010"] = num(diq["DIQ010"])
    out["type1_excluded"] = out["diq010"].eq(1) & insulin.eq(1) & diagnosis_age.lt(20)
    out = out.merge(ghb[["SEQN", "LBXGH"]].drop_duplicates("SEQN"), on="SEQN", how="left", validate="one_to_one")
    hba1c = num(out["LBXGH"])
    out["case"] = (out["diq010"].eq(1) & ~out["type1_excluded"]) | (out["diq010"].eq(2) & hba1c.ge(6.5))
    out["control"] = out["diq010"].eq(2) & hba1c.lt(6.5)
    return out[["SEQN", "case", "control"]], f"{diq_path.name};{ghb_path.name}", ""


def binary_status(data_dir: Path, cycle: str, module: str, variable: str) -> tuple[pd.DataFrame | None, str, str]:
    path = data_dir / f"{cycle}_{module}.XPT"
    if not path.exists():
        return None, "", f"Missing source file: {path.name}"
    source = read_xpt(path)
    if variable not in source:
        return None, path.name, f"Source variable unavailable in this cycle: {variable}"
    value = num(source[variable])
    return pd.DataFrame({"SEQN": source["SEQN"], "case": value.eq(1), "control": value.eq(2)}), path.name, ""


def bmi_status(data_dir: Path, cycle: str) -> tuple[pd.DataFrame | None, str, str]:
    path = data_dir / f"{cycle}_BMX.XPT"
    if not path.exists():
        return None, "", f"Missing source file: {path.name}"
    source = read_xpt(path)
    if "BMXBMI" not in source:
        return None, path.name, "Source variable unavailable in this cycle: BMXBMI"
    value = num(source["BMXBMI"])
    return pd.DataFrame({"SEQN": source["SEQN"], "case": value.ge(30), "control": value.lt(30)}), path.name, ""


def derive_status(kind: str, data_dir: Path, cycle: str) -> tuple[pd.DataFrame | None, str, str]:
    if kind == "t2d":
        return t2d_status(data_dir, cycle)
    if kind == "bmi":
        return bmi_status(data_dir, cycle)
    if kind.startswith("binary:"):
        _, module, variable = kind.split(":")
        return binary_status(data_dir, cycle, module, variable)
    return None, "", kind.removeprefix("unavailable:")


def empty_row(outcome_id: str, name: str, component: str, variable: str, definition: str, cycle: str, reason: str) -> dict:
    return {"outcome_id": outcome_id, "outcome_name": name, "cycle": cycle, "available": False,
            "male_cases": 0, "female_cases": 0, "male_controls": 0, "female_controls": 0,
            "eligible_n": 0, "male_eligible_n": 0, "female_eligible_n": 0,
            "source_component": component, "source_variables": variable, "operational_definition": definition,
            "eligibility_note": "Primary audit age >=20 years; source module unavailable or variable absent.",
            "source_file(s)": "", "unavailability_reason": reason}


def audit_one(outcome: tuple[str, str, str, str, str, str], data_dir: Path) -> list[dict]:
    outcome_id, name, component, variable, definition, kind = outcome
    rows: list[dict] = []
    for cycle in CYCLES:
        core, core_reason = core_demo(data_dir, cycle)
        status, source_files, reason = derive_status(kind, data_dir, cycle)
        if core is None:
            rows.append(empty_row(outcome_id, name, component, variable, definition, cycle, core_reason))
            continue
        if status is None:
            rows.append(empty_row(outcome_id, name, component, variable, definition, cycle, reason))
            continue
        frame = core.merge(status, on="SEQN", how="inner", validate="one_to_one")
        eligible = frame["case"] | frame["control"]
        row = {"outcome_id": outcome_id, "outcome_name": name, "cycle": cycle, "available": True,
               "male_cases": int((eligible & frame["case"] & frame["sex"].eq("male")).sum()),
               "female_cases": int((eligible & frame["case"] & frame["sex"].eq("female")).sum()),
               "male_controls": int((eligible & frame["control"] & frame["sex"].eq("male")).sum()),
               "female_controls": int((eligible & frame["control"] & frame["sex"].eq("female")).sum()),
               "eligible_n": int(eligible.sum()), "male_eligible_n": int((eligible & frame["sex"].eq("male")).sum()),
               "female_eligible_n": int((eligible & frame["sex"].eq("female")).sum()),
               "source_component": component, "source_variables": variable, "operational_definition": definition,
               "eligibility_note": "Primary audit age >=20 years. Counts require a substantive case/control definition and RIAGENDR.",
               "source_file(s)": f"{cycle}_DEMO.XPT;{source_files}", "unavailability_reason": ""}
        rows.append(row)
    return rows


def screen(male_cases: int, female_cases: int, n_cycles: int) -> tuple[str, str]:
    if n_cycles < 4:
        return "FAIL", "Fewer than 4 independent cycles are definable"
    low = min(male_cases, female_cases)
    if low >= 100:
        return "PASS", "Both pooled sex-specific case counts >=100 and >=4 cycles"
    if low >= 50:
        return "AMBER", "At least one pooled sex-specific case count is 50-99"
    return "FAIL", "At least one pooled sex-specific case count is <50"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cycle_df = pd.DataFrame([row for outcome in REGISTRY for row in audit_one(outcome, args.data_dir)])
    cycle_df = cycle_df.sort_values(["outcome_id", "cycle"]).reset_index(drop=True)
    pooled = []
    for outcome_id, name, component, variable, definition, _ in REGISTRY:
        part = cycle_df.loc[cycle_df["outcome_id"].eq(outcome_id)]
        available = part.loc[part["available"]]
        male_cases, female_cases = int(available["male_cases"].sum()), int(available["female_cases"].sum())
        status, reason = screen(male_cases, female_cases, len(available))
        pooled.append({"outcome_id": outcome_id, "outcome_name": name, "available_cycles": len(available),
                       "cycles_available": ";".join(available["cycle"]), "male_cases": male_cases,
                       "female_cases": female_cases, "male_controls": int(available["male_controls"].sum()),
                       "female_controls": int(available["female_controls"].sum()), "eligible_n": int(available["eligible_n"].sum()),
                       "source_component": component, "source_variables": variable, "operational_definition": definition,
                       "screen_status": status, "screen_reason": reason})
    pooled_df = pd.DataFrame(pooled)
    frozen = pooled_df.copy()
    frozen.insert(0, "frozen_set_version", "1.0")
    frozen.insert(1, "pre_frozen_candidate", True)
    frozen["selection_for_followup"] = frozen["screen_status"].eq("PASS")

    cycle_name, pooled_name, frozen_name = "outcome_cycle_sex_audit.csv", "outcome_pooled_sex_audit.csv", "frozen_outcome_set_v1.csv"
    cycle_df.to_csv(args.outdir / cycle_name, index=False)
    pooled_df.to_csv(args.outdir / pooled_name, index=False)
    frozen.to_csv(args.outdir / frozen_name, index=False)
    qc = {"audit_type": "OUTCOME_ONLY_SEX_SPECIFIC_FEASIBILITY_AUDIT", "generated_utc": datetime.now(timezone.utc).isoformat(),
          "data_dir": str(args.data_dir), "candidate_count": len(REGISTRY), "cycle_count": len(CYCLES),
          "input_allowlist": ["*_DEMO.XPT", "*_MCQ.XPT", "*_BPQ.XPT", "*_DIQ.XPT", "cycle-specific HbA1c XPT", "*_BMX.XPT"],
          "explicitly_not_read": ["exposure laboratory results", "Step5 exposure associations", "CRC/T2D association results", "candidate chemical results"],
          "screen_rule": "PASS: both pooled sex-specific cases >=100 and >=4 cycles; AMBER: either sex 50-99 (with >=4 cycles); FAIL: either sex <50 or <4 cycles.",
          "status_counts": pooled_df["screen_status"].value_counts().to_dict(),
          "output_hashes": {x: sha256(args.outdir / x) for x in [cycle_name, pooled_name, frozen_name]}}
    (args.outdir / "outcome_inventory_qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    lines = ["# Frozen Outcome Set v1.0: sex-specific feasibility audit", "", "This is an outcome-only inventory. The program reads demographics and the small, explicit source-file allowlist in `outcome_inventory_qc.json`; it does not read exposure measurements, any association result, or candidate-chemical result.", "", "## Prespecified screen", "", "PASS requires male and female pooled cases each >=100 and at least 4 independent definable cycles. AMBER applies when either sex has 50-99 cases (and >=4 cycles). FAIL applies when either sex has <50 cases or fewer than 4 cycles are definable.", "", "## Local-source limitations", "", "Objective CKD, kidney stones, osteoporosis, and clinically relevant depressive symptoms are retained in the frozen 18-candidate registry but are marked unavailable where their required raw modules are absent locally. Albuminuria alone is not substituted for eGFR-defined CKD.", "", "## Outputs", "", f"- `{cycle_name}` — one row per candidate outcome and cycle, including sex-specific cases and controls.", f"- `{pooled_name}` — pooled sex-specific counts and prespecified status.", f"- `{frozen_name}` — all 18 pre-frozen candidates with follow-up flag.", "- `outcome_inventory_qc.json` — provenance, input allowlist, prohibition statement, and hashes.", ""]
    (args.outdir / "README.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
