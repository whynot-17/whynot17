#!/usr/bin/env python
"""Phase 8-R2: audit the GDSC drug universe before external validation.

This script does not re-rank drugs. It reports the opportunity set available
to the locked Phase 8-R screen and separates coverage limitation from a
biological no-hit. Clinical/context labels are derived from the harmonized
GDSC compound table and official PRISM secondary-screen metadata when a
normalized drug-name match exists. Unmatched records remain unresolved.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
DATA = ROOT / "work" / "phase7_convergent_vulnerability" / "data"
ASSOC = OUT / "phase8R_drug_association_by_trajectory.csv"
RANK = OUT / "phase8R_global_empirical_drug_ranking.csv"
COMPOUND = DATA / "screened_compounds_rel_8.4.csv"
PRISM_META = ROOT / "work" / "phase8R_prism" / "raw" / "secondary-screen-dose-response-curve-parameters.csv"


def drug_key(x: object) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


ONCOLOGY_WORDS = re.compile(
    r"cancer|carcinoma|neoplasm|tumou?r|oncolog|leukemia|lymphoma|myeloma|"
    r"melanoma|sarcoma|glioma|blastoma|metast|hematologic malignancy", re.I)
NONONCOLOGY_WORDS = re.compile(
    r"arthritis|infection|malaria|viral|bacterial|fungal|diabetes|hypertension|"
    r"epilepsy|seizure|pain|asthma|inflammation|autoimmune|transplant|gout|"
    r"obesity|contracept|depression|schizophrenia|hiv|addiction|alcohol|"
    r"smoking|thrombosis|ulcer|parkinson|multiple sclerosis|psoriasis|"
    r"dermatitis|erectile|conjunctivitis|rhinitis|ulcerative colitis|crohn|"
    r"allergy|ophthalm|rheumat|pulmonary|gastroenter|cardiology|endocrinology|"
    r"obstetric|urology", re.I)


def phase_status(value: object) -> str:
    text = str(value or "").strip().lower()
    if text == "launched":
        return "approved"
    if text.startswith("phase"):
        return "investigational"
    if text == "preclinical":
        return "preclinical"
    if text == "withdrawn":
        return "withdrawn"
    return "unresolved"


def prism_metadata() -> pd.DataFrame:
    columns = ["prism_name", "drug_key", "prism_phase", "prism_disease_area", "prism_indication"]
    if not PRISM_META.exists():
        return pd.DataFrame(columns=columns)
    parts = []
    for chunk in pd.read_csv(PRISM_META, usecols=["name", "phase", "disease.area", "indication"],
                             chunksize=200_000, low_memory=False):
        parts.append(chunk)
    m = pd.concat(parts, ignore_index=True).drop_duplicates("name")
    m = m.rename(columns={"name": "prism_name", "phase": "prism_phase",
                          "disease.area": "prism_disease_area", "indication": "prism_indication"})
    m["drug_key"] = m["prism_name"].map(drug_key)
    m["prism_name_count_by_key"] = m.groupby("drug_key")["prism_name"].transform("nunique")
    return m.sort_values(["drug_key", "prism_name"]).drop_duplicates("drug_key")


def classify(row: pd.Series) -> pd.Series:
    status = phase_status(row.get("prism_phase"))
    area = str(row.get("prism_disease_area", "") or "")
    indication = str(row.get("prism_indication", "") or "")
    target_text = " ".join(str(row.get(x, "") or "") for x in ["TARGET", "TARGET_PATHWAY"])
    context_text = " ".join([area, indication, target_text])
    oncology = bool(ONCOLOGY_WORDS.search(context_text))
    nononc = bool(NONONCOLOGY_WORDS.search(context_text))
    if status == "approved" and oncology and nononc:
        context, eligible = "approved_mixed_oncology_nononcology", False
    elif status == "approved" and oncology:
        context, eligible = "approved_oncology", False
    elif status == "approved" and nononc:
        context, eligible = "approved_nononcology_high_confidence", True
    elif status == "approved":
        context, eligible = "approved_indication_unresolved", False
    elif status == "investigational" and oncology:
        context, eligible = "investigational_oncology", False
    elif status == "investigational":
        context, eligible = "investigational_nononcology_or_unresolved", False
    elif status == "preclinical":
        context, eligible = "preclinical_or_unresolved", False
    elif status == "withdrawn":
        context, eligible = "withdrawn", False
    else:
        context, eligible = "unresolved", False
    return pd.Series({"clinical_status": status, "clinical_context": context,
                      "approved_nononcology_high_confidence": eligible})


def main() -> None:
    assoc = pd.read_csv(ASSOC)
    assoc = assoc[(assoc["signature_size"] == 250) & (assoc["scoring_method"] == "weighted")].copy()
    coverage = (assoc.groupby(["database", "DRUG_NAME", "drug_key"], as_index=False)
                .agg(n_trajectory_rows=("n", "size"),
                     n_trajectories_n20=("n", lambda x: int((x >= 20).sum())),
                     n_trajectories_n30=("n", lambda x: int((x >= 30).sum())),
                     max_crc_n=("n", "max"), median_crc_n=("n", "median")))
    coverage["any_trajectory_n20"] = coverage["max_crc_n"] >= 20
    coverage["any_trajectory_n30"] = coverage["max_crc_n"] >= 30
    coverage["three_trajectories_n20"] = coverage["n_trajectories_n20"] >= 3
    coverage["three_trajectories_n30"] = coverage["n_trajectories_n30"] >= 3
    coverage["all_four_trajectories_n20"] = coverage["n_trajectories_n20"] >= 4

    rank = pd.read_csv(RANK)
    rank = rank[(rank["signature_size"] == 250) & (rank["scoring_method"] == "weighted")]
    coverage = coverage.merge(rank[["database", "DRUG_NAME", "primary_discovery_gate"]],
                              on=["database", "DRUG_NAME"], how="left")
    compound = pd.read_csv(COMPOUND) if COMPOUND.exists() else pd.DataFrame()
    if len(compound):
        compound = compound.drop_duplicates("DRUG_NAME")
        coverage = coverage.merge(compound[["DRUG_NAME", "SYNONYMS", "TARGET", "TARGET_PATHWAY"]],
                                  on="DRUG_NAME", how="left")
    for c in ["SYNONYMS", "TARGET", "TARGET_PATHWAY"]:
        if c not in coverage:
            coverage[c] = ""
    coverage["drug_key"] = coverage["DRUG_NAME"].map(drug_key)
    coverage = coverage.merge(prism_metadata(), on="drug_key", how="left")
    coverage[["clinical_status", "clinical_context", "approved_nononcology_high_confidence"]] = coverage.apply(classify, axis=1)
    coverage["annotation_source"] = np.where(coverage["prism_name"].notna(),
        "GDSC screened_compounds_rel_8.4 + normalized PRISM secondary metadata",
        "GDSC screened_compounds_rel_8.4; no normalized PRISM metadata match")
    coverage.to_csv(OUT / "phase8R2_gdsc_drug_universe_audit.csv", index=False)

    def summary_frame(sub: pd.DataFrame, label: str) -> dict[str, object]:
        return {"universe_definition": label, "n_database_drug_entries": int(len(sub)),
                "n_unique_drugs": int(sub["DRUG_NAME"].nunique()),
                "n_approved": int((sub.clinical_status == "approved").sum()),
                "n_approved_nononcology_high_confidence": int(sub.approved_nononcology_high_confidence.sum()),
                "n_investigational": int((sub.clinical_status == "investigational").sum()),
                "n_preclinical": int((sub.clinical_status == "preclinical").sum()),
                "n_withdrawn": int((sub.clinical_status == "withdrawn").sum()),
                "n_unresolved": int((sub.clinical_status == "unresolved").sum()),
                "n_approved_oncology": int((sub.clinical_context == "approved_oncology").sum()),
                "n_approved_mixed_oncology_nononcology": int((sub.clinical_context == "approved_mixed_oncology_nononcology").sum()),
                "n_prism_metadata_matched": int(sub.prism_name.notna().sum())}

    summaries = []
    for db, g in coverage.groupby("database"):
        for label, mask in [("all_drug_entries", np.ones(len(g), dtype=bool)),
                            ("any_trajectory_n20", g.any_trajectory_n20),
                            ("any_trajectory_n30", g.any_trajectory_n30),
                            ("three_trajectories_n20", g.three_trajectories_n20),
                            ("three_trajectories_n30", g.three_trajectories_n30),
                            ("all_four_trajectories_n20", g.all_four_trajectories_n20)]:
            row = summary_frame(g.loc[mask], label); row["database"] = db; summaries.append(row)
    union = coverage.sort_values("database").drop_duplicates("drug_key")
    for label, mask in [("GDSC_union_all_drugs", np.ones(len(union), dtype=bool)),
                        ("GDSC_union_any_trajectory_n20", union.any_trajectory_n20),
                        ("GDSC_union_any_trajectory_n30", union.any_trajectory_n30),
                        ("GDSC_union_three_trajectories_n20", union.three_trajectories_n20),
                        ("GDSC_union_three_trajectories_n30", union.three_trajectories_n30)]:
        row = summary_frame(union.loc[mask], label); row["database"] = "GDSC_union"; summaries.append(row)
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "phase8R2_gdsc_coverage_summary.csv", index=False)

    lines = ["# Phase 8-R2：GDSC drug-universe coverage audit", "",
        "This audit does not re-rank drugs. It quantifies the opportunity set available to the locked Phase 8-R phenotype-first screen and separates coverage limitation from biological non-hit.", "",
        "## Definitions", "",
        "- `all_drug_entries`: distinct drug names present after CRC/expression mapping in the Phase 8-R association table.",
        "- `any_trajectory_n20/n30`: at least one frozen trajectory has at least 20/30 CRC models after self-line exclusion.",
        "- `three_trajectories_n20/n30`: at least three frozen trajectories meet the threshold.",
        "- `all_four_trajectories_n20`: all four biological backgrounds have evidence at n>=20.",
        "- `approved_nononcology_high_confidence`: PRISM metadata phase `Launched` plus a non-oncology context term and no oncology context term. Unmatched records are not counted.", "",
        "## Summary", "", "```text", summary.to_string(index=False), "```", "",
        "Clinical labels are audit annotations only and do not influence phenotype ranking. PRISM metadata is used here only to quantify GDSC opportunity-set coverage; biological replication is run independently.",
        "Raw GDSC/PRISM files are not committed. Only derived audit tables, report and manifest are versioned."]
    (OUT / "phase8R2_gdsc_coverage_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest = {"phase": "8-R2", "component": "GDSC_drug_universe_audit",
        "locked_screen": "Phase 8-R top250 weighted; six trajectories; self-line exclusion retained",
        "outputs": ["phase8R2_gdsc_drug_universe_audit.csv", "phase8R2_gdsc_coverage_summary.csv", "phase8R2_gdsc_coverage_audit.md"],
        "annotation_source": "local GDSC screened compound metadata + normalized PRISM secondary-screen metadata",
        "approved_nononcology_rule": "PRISM phase=Launched, non-oncology context term present, oncology context term absent; unresolved matches excluded",
        "raw_data_policy": "raw data remain local and are not committed"}
    (OUT / "phase8R2_gdsc_coverage_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
