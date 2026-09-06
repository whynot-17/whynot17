#!/usr/bin/env python3
"""Association-free KNHANES outcome feasibility census.

This script deliberately does not read exposure-outcome result files. It uses
only the frozen environmental test crosswalk/QC and a predeclared outcome
feasibility ledger. Counts that require authorized KNHANES microdata are kept
explicitly non-calculable rather than inferred from publications.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
FROZEN_TESTS = REPO / "analysis/disease_agnostic_environmental_framework/step04_testset_freeze/unique_biomarker_test_set.csv"
CROSSWALK = ROOT / "knhanes_29_test_crosswalk.csv"
PRIOR_QC = ROOT / "KNHANES_QC_SUMMARY.json"
RAW_CATALOG = ROOT / "knhanes_raw_data_catalog.csv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def main() -> None:
    tests = read_csv(FROZEN_TESTS)
    crosswalk = read_csv(CROSSWALK)
    prior_qc = json.loads(PRIOR_QC.read_text(encoding="utf-8"))

    expected_tests = int(prior_qc["frozen_test_count_expected"])
    if len(tests) != expected_tests or len(crosswalk) != expected_tests:
        raise RuntimeError(f"Frozen test count mismatch: tests={len(tests)}, crosswalk={len(crosswalk)}, expected={expected_tests}")

    exact = sum(row.get("exact_analyte_confirmed", "").strip().lower() == "yes" and row.get("exact_matrix_confirmed", "").strip().lower() == "yes" for row in crosswalk)
    matrix_mismatch = sum(row.get("knhanes_match_class", "") == "family_public_matrix_mismatch" for row in crosswalk)
    if exact != int(prior_qc["exact_public_match_count"]) or matrix_mismatch != int(prior_qc["family_public_matrix_mismatch_count"]):
        raise RuntimeError("Prior KNHANES exposure crosswalk QC does not reconcile")

    # This ledger is outcome-free by construction. Scores are assigned from
    # survey content/design and published variable-definition precedent only.
    common = {
        "same_person_status": "documented_by_survey_structure_or_published_precedent; exact joint extract pending",
        "core_covariates": "age; sex; BMI/anthropometry; smoking; alcohol; physical activity; SES/education/income",
        "core_covariate_score": 2,
        "survey_design_status": "complex probability sample documented; current weight/strata/PSU field names pending codebook",
        "survey_design_score": 2,
        "access_status": "official registration/consent raw-data route confirmed; download not completed",
        "access_score": 2,
        "case_count_now": "not_calculable_without_registered_raw_data",
        "exposure_overlap_exact_29": 0,
        "exposure_overlap_nearest_matrix_mismatch": 2,
        "exposure_join_status": "conditional; prior public audit found no exact frozen 29-test match",
        "association_results_consulted": "false",
        "selection_basis": "feasibility-only; no exposure-outcome association, P value, FDR, or published effect estimate used",
    }

    candidates = [
        {
            "outcome_id": "hypertension",
            "outcome_label": "Hypertension",
            "outcome_class": "clinical disease outcome",
            "definition_components": "measured systolic/diastolic blood pressure; antihypertensive medication; diagnosis if available",
            "definition_reproducibility_score": 2,
            "expected_event_density": "high",
            "expected_event_density_score": 2,
            "same_person_components_score": 2,
            "cycle_window_basis": "BP examination and health-interview content documented across KNHANES cycles",
            "definition_status": "highly reproducible by direct measurement plus medication; exact cycle coding pending",
            "candidate_role": "primary disease candidate",
            "source_evidence": "KDCA_OFFICIAL_EN;HTN_KNHANES_PRECEDENT",
        },
        {
            "outcome_id": "diabetes",
            "outcome_label": "Diabetes / glycemic disease outcome",
            "outcome_class": "clinical disease outcome",
            "definition_components": "fasting glucose and/or HbA1c; physician diagnosis; glucose-lowering medication/insulin; type-1 exclusion if identifiable",
            "definition_reproducibility_score": 1,
            "expected_event_density": "high",
            "expected_event_density_score": 2,
            "same_person_components_score": 2,
            "cycle_window_basis": "glycemic examination and diabetes interview/medication content documented; fasting-subsample restriction expected",
            "definition_status": "feasible but T2D-specific/type-1 coding and current field names require codebook confirmation",
            "candidate_role": "primary disease candidate",
            "source_evidence": "KDCA_OFFICIAL_EN;DIABETES_KNHANES_PRECEDENT",
        },
        {
            "outcome_id": "dyslipidemia",
            "outcome_label": "Dyslipidemia / hypercholesterolemia",
            "outcome_class": "clinical phenotype outcome",
            "definition_components": "fasting total cholesterol/LDL/HDL/triglycerides; lipid-lowering medication; diagnosis if available",
            "definition_reproducibility_score": 1,
            "expected_event_density": "high",
            "expected_event_density_score": 2,
            "same_person_components_score": 2,
            "cycle_window_basis": "fasting lipid examination and medication/health interview content documented; harmonized window required",
            "definition_status": "feasible; threshold/definition changes across years require prespecified harmonization",
            "candidate_role": "primary clinical-outcome candidate",
            "source_evidence": "KDCA_OFFICIAL_EN;DYSLIPIDEMIA_KNHANES_PRECEDENT",
        },
        {
            "outcome_id": "metabolic_syndrome",
            "outcome_label": "Metabolic syndrome",
            "outcome_class": "derived clinical phenotype",
            "definition_components": "prespecified 3-of-5 rule using waist, triglycerides, HDL, BP, and glucose; medication alternatives where justified",
            "definition_reproducibility_score": 1,
            "expected_event_density": "high",
            "expected_event_density_score": 2,
            "same_person_components_score": 2,
            "cycle_window_basis": "requires synchronized anthropometry, BP, fasting lipid, and glycemic components",
            "definition_status": "derivable but jointly complete fasting/examination frame must be confirmed",
            "candidate_role": "reserve derived phenotype",
            "source_evidence": "KDCA_OFFICIAL_EN;METABOLIC_SYNDROME_DEFINITION_REFERENCE",
        },
        {
            "outcome_id": "ckd",
            "outcome_label": "Chronic kidney disease",
            "outcome_class": "clinical disease outcome",
            "definition_components": "eGFR from creatinine with prespecified chronicity/albuminuria rule; diagnosis/medication if available",
            "definition_reproducibility_score": 1,
            "expected_event_density": "moderate",
            "expected_event_density_score": 1,
            "same_person_components_score": 2,
            "cycle_window_basis": "renal laboratory and urine components require cycle-specific codebook confirmation",
            "definition_status": "conditional; creatinine harmonization, albuminuria availability, and chronicity rule pending",
            "candidate_role": "reserve disease candidate",
            "source_evidence": "KDCA_OFFICIAL_EN;KNHANES_CODEBOOK_PENDING",
        },
        {
            "outcome_id": "obesity",
            "outcome_label": "Obesity / adiposity phenotype",
            "outcome_class": "non-disease phenotype contrast",
            "definition_components": "measured height/weight BMI; waist circumference; prespecified adult cutoffs",
            "definition_reproducibility_score": 2,
            "expected_event_density": "very high",
            "expected_event_density_score": 2,
            "same_person_components_score": 2,
            "cycle_window_basis": "anthropometric examination is a core KNHANES component",
            "definition_status": "highly feasible but retained as a phenotype contrast, not a primary disease replacement",
            "candidate_role": "contrast only; not eligible for primary disease freeze",
            "source_evidence": "KDCA_OFFICIAL_EN",
        },
        {
            "outcome_id": "liver_disease",
            "outcome_label": "Liver disease / liver injury phenotype",
            "outcome_class": "conditional clinical phenotype",
            "definition_components": "liver enzymes and/or imaging/diagnosis; prespecified disease subtype rule",
            "definition_reproducibility_score": 1,
            "expected_event_density": "moderate",
            "expected_event_density_score": 1,
            "same_person_components_score": 1,
            "cycle_window_basis": "health examination includes liver-related content, but subtype and complete frame require codebook confirmation",
            "definition_status": "conditional; disease subtype and imaging/diagnosis availability not yet frozen",
            "candidate_role": "reserve conditional outcome",
            "source_evidence": "KDCA_OFFICIAL_EN;KNHANES_CODEBOOK_PENDING",
        },
    ]

    rows: list[dict[str, object]] = []
    for candidate in candidates:
        row = {**common, **candidate}
        row["feasibility_score"] = sum(
            int(row[key])
            for key in (
                "definition_reproducibility_score",
                "expected_event_density_score",
                "same_person_components_score",
                "core_covariate_score",
                "survey_design_score",
                "access_score",
            )
        )
        row["raw_case_count_available_now"] = "no"
        rows.append(row)

    # The disease/clinical-outcome freeze is a priori and excludes the
    # The primary freeze excludes derived multicomponent phenotypes and the
    # phenotype-only contrast even when they are highly feasible.
    eligible = [r for r in rows if r["candidate_role"] in {"primary disease candidate", "primary clinical-outcome candidate"}]
    eligible.sort(key=lambda r: (-int(r["feasibility_score"]), str(r["outcome_id"])))
    for rank, row in enumerate(eligible, 1):
        row["preassociation_freeze_rank"] = rank
        row["preassociation_freeze_status"] = "provisional_primary_freeze" if rank <= 3 else "not_selected_primary_freeze"
    for row in rows:
        row.setdefault("preassociation_freeze_rank", "not_applicable")
        row.setdefault("preassociation_freeze_status", "reserve_or_contrast")

    fields = [
        "outcome_id", "outcome_label", "outcome_class", "candidate_role", "definition_components",
        "definition_status", "definition_reproducibility_score", "expected_event_density", "expected_event_density_score",
        "same_person_status", "same_person_components_score", "core_covariates", "core_covariate_score",
        "survey_design_status", "survey_design_score", "access_status", "access_score", "feasibility_score",
        "cycle_window_basis", "case_count_now", "raw_case_count_available_now", "exposure_overlap_exact_29",
        "exposure_overlap_nearest_matrix_mismatch", "exposure_join_status", "source_evidence", "selection_basis",
        "association_results_consulted", "preassociation_freeze_rank", "preassociation_freeze_status",
    ]
    write_csv(ROOT / "knhanes_outcome_candidate_census.csv", rows, fields)

    frozen = [r for r in eligible if r["preassociation_freeze_status"] == "provisional_primary_freeze"]
    freeze_fields = [
        "freeze_order", "outcome_id", "outcome_label", "outcome_class", "feasibility_score",
        "freeze_reason", "association_results_consulted", "p_values_consulted", "published_effect_estimates_consulted",
        "frozen_before_association", "requires_raw_file_confirmation", "exposure_overlap_used_for_selection",
    ]
    freeze_rows = []
    for i, row in enumerate(frozen, 1):
        freeze_rows.append({
            "freeze_order": i,
            "outcome_id": row["outcome_id"],
            "outcome_label": row["outcome_label"],
            "outcome_class": row["outcome_class"],
            "feasibility_score": row["feasibility_score"],
            "freeze_reason": "Top three eligible disease/clinical-outcome candidates under predeclared feasibility rubric",
            "association_results_consulted": "false",
            "p_values_consulted": "false",
            "published_effect_estimates_consulted": "false",
            "frozen_before_association": "true",
            "requires_raw_file_confirmation": "true",
            "exposure_overlap_used_for_selection": "false",
        })
    write_csv(ROOT / "knhanes_outcome_freeze_matrix.csv", freeze_rows, freeze_fields)

    catalog_rows = read_csv(RAW_CATALOG)
    summary = {
        "generated_utc": now_utc(),
        "status": "preassociation_outcome_feasibility_census_complete_provisional_freeze",
        "scope": "KNHANES independent population replacement feasibility; no exposure-outcome models",
        "candidate_outcome_count": len(rows),
        "provisional_primary_freeze": [r["outcome_id"] for r in frozen],
        "provisional_primary_freeze_count": len(frozen),
        "frozen_test_count": len(tests),
        "frozen_exposure_exact_match_count": exact,
        "frozen_exposure_nearest_matrix_mismatch_count": matrix_mismatch,
        "raw_catalog_record_count": len(catalog_rows),
        "raw_catalog_total_count_from_prior_audit": int(prior_qc["raw_catalog_total_count"]),
        "case_counts_calculable_now": False,
        "association_models_run": False,
        "association_results_consulted": False,
        "p_values_consulted": False,
        "published_effect_estimates_imported": False,
        "selection_by_association": False,
        "microdata_downloaded": False,
        "user_registration_attempted": False,
        "final_replication_readiness": "conditional_pending_authorized_raw_file_and_codebook_confirmation",
        "primary_freeze_rule": "Top three eligible directly defined disease/clinical-outcome candidates by feasibility-only score; derived multicomponent and non-disease phenotypes remain reserve/contrast",
        "important_boundary": "Exact frozen exposure overlap is recorded for every outcome but is not used to rank or select outcomes; prior audit found 0/29 exact public KNHANES matches.",
        "input_sha256": {
            "frozen_test_set": sha256(FROZEN_TESTS),
            "knhanes_crosswalk": sha256(CROSSWALK),
            "prior_knhanes_qc": sha256(PRIOR_QC),
            "knhanes_raw_catalog": sha256(RAW_CATALOG),
        },
    }
    (ROOT / "KNHANES_OUTCOME_CENSUS_QC_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    manifest = {
        "generated_utc": summary["generated_utc"],
        "script": "run_knhanes_outcome_feasibility_census.py",
        "method": "preassociation feasibility-only outcome census",
        "outcome_firewall": {
            "association_models_run": False,
            "association_results_consulted": False,
            "p_values_consulted": False,
            "published_effect_estimates_imported": False,
            "selection_by_association": False,
        },
        "frozen_input": "analysis/disease_agnostic_environmental_framework/step04_testset_freeze/unique_biomarker_test_set.csv",
        "prior_exposure_audit": "KNHANES_QC_SUMMARY.json",
        "candidate_outcomes": [r["outcome_id"] for r in rows],
        "provisional_primary_freeze": [r["outcome_id"] for r in frozen],
        "outputs": [
            "PLAN_OUTCOME_CENSUS.md",
            "run_knhanes_outcome_feasibility_census.py",
            "knhanes_outcome_candidate_census.csv",
            "knhanes_outcome_freeze_matrix.csv",
            "KNHANES_OUTCOME_CENSUS_QC_SUMMARY.json",
            "STEP10B_P3_OUTCOME_CENSUS_REPORT.md",
            "STEP10B_P3_OUTCOME_CENSUS_MANIFEST.json",
        ],
    }
    report = f"""# Step 10B-P3 — KNHANES disease-feasibility census

Generated: `{summary['generated_utc']}`  \
Status: **pre-association outcome feasibility census complete; provisional freeze only**

## Decision

The predeclared feasibility rubric provisionally freezes three eligible directly defined disease/clinical-outcome candidates before any exposure–outcome association is inspected:

1. **Hypertension** — score 12/12;
2. **Diabetes / glycemic disease outcome** — score 11/12;
3. **Dyslipidemia / hypercholesterolemia** — score 11/12.

This is a **data-feasibility freeze**, not a claim that KNHANES is ready for immediate replication. The authorized raw files and codebooks must still confirm exact fields, same-person exposure/outcome overlap, missingness, harmonized cycle windows, and current weight/strata/PSU variables.

## What was and was not used

- Candidate selection used only outcome definition reproducibility, expected event/phenotype density, same-person component availability, core covariate availability, survey-design documentation, and official access feasibility.
- No exposure–outcome model, P value, FDR, effect estimate, or association direction was read.
- The frozen 29-test exposure crosswalk was carried forward only as a common constraint: **{exact}/29 exact public matches and {matrix_mismatch}/29 related blood-matrix mismatches** in the prior audit. This exposure overlap was not used to rank outcomes.
- Current case counts are not calculable without registered raw data and are therefore not imputed from publications.

## Candidate census

| Outcome | Class | Score | Role |
|---|---|---:|---|
""" + "\n".join(f"| {r['outcome_label']} | {r['outcome_class']} | {r['feasibility_score']}/12 | {r['candidate_role']} |" for r in rows) + f"""

Obesity scored highly but is retained only as a non-disease phenotype contrast. Metabolic syndrome remains a reserve derived phenotype even though its component variables are plausible, because the primary freeze excludes multicomponent derived outcomes. CKD and liver disease remain reserve/conditional outcomes because renal/liver definitions, synchronized components, or codebook-level availability require additional confirmation.

## Why the three provisional candidates are feasible

- **Hypertension:** direct blood-pressure examination plus medication/diagnosis information gives a reproducible outcome pathway and high expected event density. KNHANES blood-pressure measurement and hypertension definitions are documented in published analyses.
- **Diabetes:** glycemic laboratory measures plus diagnosis/medication fields are documented, but T2D-specific coding, type-1 exclusion, and fasting-subsample handling remain pending.
- **Dyslipidemia:** fasting lipid measurements and medication information are documented, but the historical threshold/definition changes require a harmonized cycle window before modeling.

The KDCA survey overview describes health examinations and interviews covering obesity, hypertension, diabetes, dyslipidemia, kidney/liver disease, smoking, drinking, physical activity, and socioeconomic/health information. The official raw-data route is registration/consent based; no personal information was entered and no microdata were downloaded in this census. [KDCA survey overview](https://kdca.go.kr/eng/4428/subview.do) · [official raw-data record](https://www.data.go.kr/tcs/dss/selectFileDataDetailView.do?publicDataPk=15076556)

Published KNHANES precedents support the feasibility of measured hypertension definitions, glycemic/diabetes components, and fasting lipid outcomes, but they are not imported as our external estimates: [hypertension precedent](https://pmc.ncbi.nlm.nih.gov/articles/PMC4661365/) · [diabetes/glycemia precedent](https://pmc.ncbi.nlm.nih.gov/articles/PMC3678002/) · [dyslipidemia precedent](https://pmc.ncbi.nlm.nih.gov/articles/PMC12488789/)

## Gate before association analysis

The three frozen outcomes may enter an external population-replacement analysis only after authorized KNHANES raw data/codebook review confirms:

1. exact outcome variables and coding;
2. exposure laboratory file and outcome components in the same persons;
3. analytic N, event counts, missingness, and cycle harmonization;
4. correct survey weight, strata, and PSU variables;
5. a prespecified outcome definition and exclusions.

Until then, the correct status is **provisional feasibility freeze, access/codebook confirmation pending**.
    """
    (ROOT / "STEP10B_P3_OUTCOME_CENSUS_REPORT.md").write_text(report, encoding="utf-8")

    manifest["outputs"] = [
        {
            "path": name,
            "bytes": (ROOT / name).stat().st_size,
            "sha256": sha256(ROOT / name),
        }
        for name in (
            "PLAN_OUTCOME_CENSUS.md",
            "run_knhanes_outcome_feasibility_census.py",
            "knhanes_outcome_candidate_census.csv",
            "knhanes_outcome_freeze_matrix.csv",
            "KNHANES_OUTCOME_CENSUS_QC_SUMMARY.json",
            "STEP10B_P3_OUTCOME_CENSUS_REPORT.md",
        )
    ]
    # The manifest itself is written last and is intentionally omitted from
    # its own inventory to avoid a self-referential hash.
    (ROOT / "STEP10B_P3_OUTCOME_CENSUS_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
