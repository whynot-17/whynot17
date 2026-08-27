"""Build the outcome-blind NHANES disease registry for Step 10.

Only outcome/covariate files are read here. The script deliberately does not
load the 29-test exposure registry or any outcome-association result.
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
DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
OUT_DIR = Path(__file__).resolve().parent
MIN_CYCLES = 3
MIN_CASES = 100

CYCLES = [
    "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
    "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018",
]

# These are stable physician-diagnosis binary items in the NHANES Medical
# Conditions Questionnaire. MCQ160J was retired after 2003-2004 but still
# satisfies the prespecified >=3-cycle technical rule.
CANDIDATES = [
    ("MCQ160A", "arthritis", "Physician-diagnosed arthritis"),
    ("MCQ160B", "congestive_heart_failure", "Physician-diagnosed congestive heart failure"),
    ("MCQ160C", "coronary_heart_disease", "Physician-diagnosed coronary heart disease"),
    ("MCQ160D", "angina", "Physician-diagnosed angina or angina pectoris"),
    ("MCQ160E", "heart_attack", "Physician-diagnosed heart attack"),
    ("MCQ160F", "stroke", "Physician-diagnosed stroke"),
    ("MCQ160G", "emphysema", "Physician-diagnosed emphysema"),
    ("MCQ160J", "hay_fever", "Physician-diagnosed hay fever"),
    ("MCQ160K", "chronic_bronchitis", "Physician-diagnosed chronic bronchitis"),
    ("MCQ160L", "liver_condition", "Physician-diagnosed liver condition"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_xpt(path: Path) -> pd.DataFrame:
    return pd.read_sas(path, format="xport", encoding="latin1")


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def derive_core(demo: pd.DataFrame, bmx: pd.DataFrame, smq: pd.DataFrame) -> pd.DataFrame:
    d = demo.copy()
    d["age"] = numeric(d["RIDAGEYR"])
    d["sex"] = numeric(d["RIAGENDR"]).map({1: "Male", 2: "Female"})
    race_col = "RIDRETH3" if "RIDRETH3" in d.columns else "RIDRETH1"
    race_map = ({1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White",
                 4: "Non-Hispanic Black", 6: "Other/Multi", 7: "Other/Multi"}
                if race_col == "RIDRETH3" else
                {1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White",
                 4: "Non-Hispanic Black", 5: "Other/Multi"})
    d["race"] = numeric(d[race_col]).map(race_map)
    d["pir"] = numeric(d["INDFMPIR"])
    d["psu_raw"] = numeric(d["SDMVPSU"])
    d["strata_raw"] = numeric(d["SDMVSTRA"])
    b = bmx[["SEQN", "BMXBMI"]].copy()
    b["bmi"] = numeric(b["BMXBMI"])
    b = b[["SEQN", "bmi"]]
    s = smq[[c for c in ["SEQN", "SMQ020", "SMQ040"] if c in smq.columns]].copy()
    ever = numeric(s.get("SMQ020"))
    current = numeric(s.get("SMQ040"))
    s["smoking"] = pd.Series(pd.NA, index=s.index, dtype="object")
    s.loc[ever.eq(2), "smoking"] = "Never"
    s.loc[ever.eq(1) & current.isin([1, 2]), "smoking"] = "Current"
    s.loc[ever.eq(1) & current.eq(3), "smoking"] = "Former"
    s = s[["SEQN", "smoking"]]
    return d[["SEQN", "age", "sex", "race", "pir", "psu_raw", "strata_raw"]].merge(
        b, on="SEQN", how="left", validate="one_to_one"
    ).merge(s, on="SEQN", how="left", validate="one_to_one")


def build_candidate(variable: str, disease_slug: str, disease_name: str, data_dir: Path) -> tuple[dict, list[dict]]:
    per_cycle: list[dict] = []
    source_files: list[str] = []
    for cycle_index, cycle in enumerate(CYCLES):
        mcq_path = data_dir / f"{cycle}_MCQ.XPT"
        demo_path = data_dir / f"{cycle}_DEMO.XPT"
        bmx_path = data_dir / f"{cycle}_BMX.XPT"
        smq_path = data_dir / f"{cycle}_SMQ.XPT"
        mcq = read_xpt(mcq_path)
        if variable not in mcq.columns:
            continue
        core = derive_core(read_xpt(demo_path), read_xpt(bmx_path), read_xpt(smq_path))
        frame = mcq[["SEQN", variable]].merge(core, on="SEQN", how="inner", validate="one_to_one")
        age = numeric(frame["age"])
        value = numeric(frame[variable])
        adult = age.ge(20)
        valid = adult & value.isin([1, 2])
        complete = valid & frame[["age", "sex", "race", "bmi", "smoking", "pir", "psu_raw", "strata_raw"]].notna().all(axis=1)
        cases = int((valid & value.eq(1)).sum())
        controls = int((valid & value.eq(2)).sum())
        per_cycle.append({
            "cycle": cycle,
            "cycle_index": cycle_index,
            "source_variable": variable,
            "source_file": mcq_path.name,
            "adult_rows": int(adult.sum()),
            "valid_binary_rows": int(valid.sum()),
            "cases": cases,
            "controls": controls,
            "ambiguous_or_missing": int((adult & ~valid).sum()),
            "core_complete_rows": int(complete.sum()),
            "core_complete_cases": int((complete & value.eq(1)).sum()),
            "core_complete_controls": int((complete & value.eq(2)).sum()),
        })
        source_files.extend([mcq_path.name, demo_path.name, bmx_path.name, smq_path.name])
    cycles = [row["cycle"] for row in per_cycle]
    pooled_cases = sum(row["cases"] for row in per_cycle)
    pooled_controls = sum(row["controls"] for row in per_cycle)
    pooled_adult = sum(row["adult_rows"] for row in per_cycle)
    core_complete = sum(row["core_complete_rows"] for row in per_cycle)
    core_cases = sum(row["core_complete_cases"] for row in per_cycle)
    reasons = []
    if len(cycles) < MIN_CYCLES:
        reasons.append(f"cycle coverage < {MIN_CYCLES}")
    if pooled_cases < MIN_CASES:
        reasons.append(f"pooled cases < {MIN_CASES}")
    if not per_cycle:
        reasons.append("source variable unavailable")
    eligible = not reasons
    record = {
        "disease_id": f"NHANES_{disease_slug.upper()}",
        "disease_name": disease_name,
        "outcome_definition": f"Adult participants (age >=20) with {variable}=1 versus {variable}=2 in NHANES MCQ.",
        "source_component": "MCQ",
        "source_variable(s)": variable,
        "cycles_available": ";".join(cycles),
        "n_cycles": len(cycles),
        "adult_only": True,
        "minimum_age": 20,
        "case_definition": f"{variable}=1 (yes/physician-diagnosed)",
        "control_definition": f"{variable}=2 (no)",
        "ambiguous_definition": f"Adult responses other than {variable}=1 or 2, including missing/refused/don't know, are indeterminate.",
        "case_count_pooled": pooled_cases,
        "control_count_pooled": pooled_controls,
        "missing_count": pooled_adult - pooled_cases - pooled_controls,
        "covariate_compatibility": "age; sex; race/ethnicity; BMI; smoking; PIR; SDMVPSU; SDMVSTRA",
        "core_complete_case_n": core_complete,
        "core_complete_case_count": core_cases,
        "core_complete_control_count": sum(row["core_complete_controls"] for row in per_cycle),
        "outcome_definition_reproducible": True,
        "non_overlap_worked_examples": True,
        "eligible_for_randomization": eligible,
        "exclusion_reason": "; ".join(reasons) if reasons else "",
        "definition_provenance": "NHANES MCQ physician-diagnosis binary item; values 1/2 only",
        "source_files": ";".join(sorted(set(source_files))),
    }
    return record, per_cycle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    registry: list[dict] = []
    cycle_rows: list[dict] = []
    for variable, slug, name in CANDIDATES:
        record, cycles = build_candidate(variable, slug, name, args.data_dir)
        registry.append(record)
        cycle_rows.extend([{**row, "disease_id": record["disease_id"]} for row in cycles])

    registry_df = pd.DataFrame(registry).sort_values("disease_id").reset_index(drop=True)
    pool_df = registry_df.loc[registry_df["eligible_for_randomization"]].copy()
    registry_df.to_csv(args.outdir / "disease_outcome_registry.csv", index=False)
    pool_df.to_csv(args.outdir / "step10b_eligible_disease_pool.csv", index=False)
    pd.DataFrame(cycle_rows).sort_values(["disease_id", "cycle"]).to_csv(args.outdir / "disease_registry_cycle_audit.csv", index=False)

    counts = registry_df["exclusion_reason"].replace("", "eligible").value_counts().to_dict()
    report = [
        "# Step 10B disease-pool audit",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "This registry was built from NHANES outcome/covariate files only. No frozen exposure test, exposure value, association estimate, P value, effect direction, or FDR result was loaded.",
        "",
        f"- Candidate outcome definitions evaluated: **{len(registry_df)}**.",
        f"- Eligible randomization pool: **{len(pool_df)}**.",
        f"- Prespecified minimum cycle coverage: **{MIN_CYCLES}**.",
        f"- Prespecified minimum pooled adult cases: **{MIN_CASES}**.",
        "- T2D and CRC remain excluded worked examples under the non-overlap rule and are not represented in the randomization pool.",
        "",
        "## Candidate registry",
        "",
        "| Disease | Variable | Cycles | Adult cases | Controls | Core-complete N | Eligible |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for _, row in registry_df.iterrows():
        report.append(f"| {row['disease_name']} | {row['source_variable(s)']} | {row['n_cycles']} | {row['case_count_pooled']} | {row['control_count_pooled']} | {row['core_complete_case_n']} | {row['eligible_for_randomization']} |")
    report += ["", "## Exclusion-reason counts", "", json.dumps(counts, indent=2, ensure_ascii=False)]
    (args.outdir / "STEP10B_DISEASE_POOL_AUDIT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = {
        "lock_type": "STEP10A_OUTCOME_BLINDED_DISEASE_REGISTRY",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "data_dir": str(args.data_dir),
        "candidate_definitions": len(registry_df),
        "eligible_pool_n": len(pool_df),
        "eligibility_rules": {"minimum_cycles": MIN_CYCLES, "minimum_pooled_adult_cases": MIN_CASES, "core_covariates": ["age", "sex", "race", "bmi", "smoking", "pir", "SDMVPSU", "SDMVSTRA"]},
        "outcome_information_used_for_exposure_selection": False,
        "association_results_loaded": False,
        "output_hashes": {name: sha256(args.outdir / name) for name in ["disease_outcome_registry.csv", "step10b_eligible_disease_pool.csv", "disease_registry_cycle_audit.csv", "STEP10B_DISEASE_POOL_AUDIT.md"]},
    }
    (args.outdir / "STEP10A_REGISTRY_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
