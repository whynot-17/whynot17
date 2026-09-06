#!/usr/bin/env python
"""Outcome-free KoNEHS population-replacement feasibility audit.

This script intentionally stops before any KoNEHS association model.  It
records source-verified analyte evidence for the frozen 29-test family and
keeps content evidence, variable-level confirmation, and data-access status
as separate fields.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


# <repo>/analysis/disease_agnostic_environmental_framework/step10_cross_database_validation/step10b_p_konehs_replacement
ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
FROZEN = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
RUN_UTC = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")

SOURCES = {
    "K4_EXPOSURE_STATUS": {
        "title": "KoNEHS cycle 4 environmental chemical exposure status (2018–2020)",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11597996/",
        "source_type": "peer_reviewed_primary_overview",
    },
    "K4_PHTHALATE_DIABETES": {
        "title": "KoNEHS cycle 4 phthalate metabolites and diabetes analysis",
        "url": "https://www.sciencedirect.com/science/article/abs/pii/S143846392400066X",
        "source_type": "peer_reviewed_feasibility_precedent",
    },
    "K2_PAH_DIABETES": {
        "title": "KoNEHS cycle 2 urinary PAHs and diabetes analysis",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7602886/",
        "source_type": "peer_reviewed_feasibility_precedent",
    },
    "K3_MCOP": {
        "title": "KoNEHS cycle 3 urinary phthalate metabolite measurement",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10011409/",
        "source_type": "peer_reviewed_exposure_precedent",
    },
    "NIER_OFFICIAL": {
        "title": "National Institute of Environmental Research environmental health research division",
        "url": "https://www.nier.go.kr/front/kor/cmm/nier_res04.do?menuNo=12004",
        "source_type": "official_program_page",
    },
    "DATA_ACCESS_NOTE": {
        "title": "KoNEHS PFAS–diabetes paper data-availability statement",
        "url": "https://www.sciencedirect.com/science/article/abs/pii/S143846392400066X",
        "source_type": "data_access_statement",
    },
}


# The mapping is deliberately conservative.  Public evidence is not treated
# as exact variable confirmation.  In particular, a source documenting blood
# lead does not establish that the frozen urine-lead test is available.
EXPOSURE_EVIDENCE = {
    "LBXPFDE": dict(domain="PFAS", analyte="PFDA/PFDeA", match="exact_public", cycles=["K4"], basis="K4_EXPOSURE_STATUS", note="Cycle-4 overview documents PFDeA/PFDA among measured PFAS."),
    "LBXPFHS": dict(domain="PFAS", analyte="PFHxS", match="exact_public", cycles=["K4"], basis="K4_EXPOSURE_STATUS", note="Cycle-4 overview documents PFHxS."),
    "LBXPFNA": dict(domain="PFAS", analyte="PFNA", match="exact_public", cycles=["K4"], basis="K4_EXPOSURE_STATUS", note="Cycle-4 overview documents PFNA."),
    "URXBPH": dict(domain="bisphenol", analyte="bisphenol A (BPA)", match="exact_public", cycles=["K4"], basis="K4_EXPOSURE_STATUS", note="Cycle-4 chemical overview supports BPA measurement; exact variable name remains pending."),
    "URXCOP": dict(domain="phthalate", analyte="MCOP/MCiOP", match="exact_public", cycles=["K3"], basis="K3_MCOP", note="Cycle-3 exposure paper explicitly reports MCOP measurement; this is not a disease replication result."),
    "URXECP": dict(domain="phthalate", analyte="MECPP", match="exact_public", cycles=["K4"], basis="K4_PHTHALATE_DIABETES", note="Cycle-4 phthalate study reports MECPP."),
    "URXMBP": dict(domain="phthalate", analyte="MnBP", match="exact_public", cycles=["K4"], basis="K4_PHTHALATE_DIABETES", note="Cycle-4 phthalate study reports MnBP."),
    "URXMEP": dict(domain="phthalate", analyte="MEP", match="exact_public", cycles=["K4"], basis="K4_PHTHALATE_DIABETES", note="Cycle-4 phthalate study reports MEP."),
    "URXMHH": dict(domain="phthalate", analyte="MEHHP", match="exact_public", cycles=["K4"], basis="K4_PHTHALATE_DIABETES", note="Cycle-4 phthalate study reports MEHHP."),
    "URXMHP": dict(domain="phthalate", analyte="MEHP", match="not_confirmed_public", cycles=[], basis="", note="Not in the audited cycle-4 eight-metabolite list; exact other-cycle evidence not confirmed in this audit."),
    "URXMIB": dict(domain="phthalate", analyte="MiBP", match="not_confirmed_public", cycles=[], basis="", note="Not in the audited cycle-4 eight-metabolite list; exact other-cycle evidence not confirmed in this audit."),
    "URXMOH": dict(domain="phthalate", analyte="MEOHP", match="exact_public", cycles=["K4"], basis="K4_PHTHALATE_DIABETES", note="Cycle-4 phthalate study reports MEOHP."),
    "URXMZP": dict(domain="phthalate", analyte="MBzP", match="exact_public", cycles=["K4"], basis="K4_PHTHALATE_DIABETES", note="Cycle-4 phthalate study reports MBzP."),
    "URXP02": dict(domain="PAH", analyte="2-naphthol / naphthalene metabolite family", match="family_public", cycles=["K2"], basis="K2_PAH_DIABETES", note="Cycle-2 paper supports a grouped naphthalene-family match, not an exact frozen NHANES variable."),
    "URXP04": dict(domain="PAH", analyte="2-hydroxyfluorene / fluorene metabolite family", match="family_public", cycles=["K2"], basis="K2_PAH_DIABETES", note="Cycle-2 paper supports a grouped fluorene-family match, not an exact frozen NHANES variable."),
    "URXP10": dict(domain="PAH", analyte="1-hydroxypyrene / pyrene metabolite family", match="family_public", cycles=["K2"], basis="K2_PAH_DIABETES", note="Cycle-2 paper supports a grouped pyrene-family match, not an exact frozen NHANES variable."),
    "URXP25": dict(domain="PAH", analyte="1-hydroxyphenanthrene / phenanthrene metabolite family", match="family_public", cycles=["K2"], basis="K2_PAH_DIABETES", note="Cycle-2 paper supports a grouped phenanthrene-family match, not an exact frozen NHANES variable."),
    "URXUBA": dict(domain="metal", analyte="barium", match="not_confirmed_public", cycles=[], basis="", note="Exact barium matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUCD": dict(domain="metal", analyte="urine cadmium", match="exact_public", cycles=["K4"], basis="K4_EXPOSURE_STATUS", note="Cycle-4 overview supports cadmium biomonitoring; exact public variable name remains pending."),
    "URXUCO": dict(domain="metal", analyte="cobalt", match="not_confirmed_public", cycles=[], basis="", note="Exact cobalt matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUCS": dict(domain="metal", analyte="cesium", match="not_confirmed_public", cycles=[], basis="", note="Exact cesium matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUMO": dict(domain="metal", analyte="molybdenum", match="not_confirmed_public", cycles=[], basis="", note="Exact molybdenum matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUPB": dict(domain="metal", analyte="lead", match="family_public", cycles=["K4"], basis="K4_EXPOSURE_STATUS", note="Public overview supports lead biomonitoring, but the audited evidence does not confirm the frozen urine-lead matrix/variable."),
    "URXUSB": dict(domain="metal", analyte="antimony", match="not_confirmed_public", cycles=[], basis="", note="Exact antimony matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUSN": dict(domain="metal", analyte="tin", match="not_confirmed_public", cycles=[], basis="", note="Exact tin matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUSR": dict(domain="metal", analyte="silver", match="not_confirmed_public", cycles=[], basis="", note="Exact silver matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUTL": dict(domain="metal", analyte="thallium", match="not_confirmed_public", cycles=[], basis="", note="Exact thallium matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUTU": dict(domain="metal", analyte="tungsten", match="not_confirmed_public", cycles=[], basis="", note="Exact tungsten matrix/cycle evidence was not confirmed in audited public sources."),
    "URXUUR": dict(domain="metal", analyte="uranium", match="not_confirmed_public", cycles=[], basis="", note="Exact uranium matrix/cycle evidence was not confirmed in audited public sources."),
}

KONEHS_CYCLES = {
    "K1": "2009–2011",
    "K2": "2012–2014",
    "K3": "2015–2017",
    "K4": "2018–2020",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fetch_source_snapshot() -> dict[str, object]:
    snapshot: dict[str, object] = {"retrieved_utc": RUN_UTC, "sources": {}}
    for source_id, metadata in SOURCES.items():
        entry = dict(metadata)
        entry["retrieval_status"] = "not_attempted"
        try:
            req = Request(metadata["url"], headers={"User-Agent": "whynot17-KoNEHS-audit/1.0"})
            with urlopen(req, timeout=8) as response:  # nosec B310: frozen public source URLs
                data = response.read(2_000_001)
                entry["http_status"] = getattr(response, "status", None)
                entry["retrieval_status"] = "ok" if len(data) <= 2_000_000 else "truncated"
                entry["content_bytes_captured"] = min(len(data), 2_000_000)
                entry["sha256_captured"] = hashlib.sha256(data[:2_000_000]).hexdigest()
                entry["content_type"] = response.headers.get("Content-Type", "")
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            entry["retrieval_status"] = "failed_best_effort"
            entry["error_type"] = type(exc).__name__
            entry["error"] = str(exc)[:300]
        snapshot["sources"][source_id] = entry
    return snapshot


def build_crosswalk(frozen: list[dict[str, str]]) -> list[dict[str, object]]:
    rows = []
    for item in frozen:
        variable = item["variable"]
        evidence = EXPOSURE_EVIDENCE.get(variable)
        if evidence is None:
            raise KeyError(f"No frozen mapping for {variable}")
        rows.append({
            "test_id": item["test_id"],
            "frozen_biomarker": item["biomarker"],
            "nhanes_variable": variable,
            "nhanes_matrix": item["matrix"],
            "frozen_nhanes_cycles": item["cycles"],
            "domain": evidence["domain"],
            "konehs_analyte_or_family": evidence["analyte"],
            "konehs_match_class": evidence["match"],
            "konehs_cycles_with_public_evidence": ";".join(evidence["cycles"]),
            "evidence_basis": evidence["basis"] or "none_in_audited_public_sources",
            "evidence_sources": evidence["basis"] or "",
            "adult_scope": "adult scope documented in cited studies where applicable; exact eligible extract pending",
            "konehs_exact_variable_confirmed": "not_yet",
            "konehs_exact_matrix_confirmed": "not_yet",
            "konehs_weight_variable_confirmed": "not_yet",
            "konehs_design_variables_confirmed": "not_yet",
            "konehs_same_person_exposure_t2d_confirmed": "not_yet",
            "konehs_individual_level_access": "controlled_or_request_to_NIER",
            "execution_status": "conditional_feasibility" if evidence["match"] != "not_confirmed_public" else "not_confirmed_publicly",
            "notes": evidence["note"],
        })
    return rows


def build_cycle_readiness(crosswalk: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for item in crosswalk:
        public_cycles = set(filter(None, str(item["konehs_cycles_with_public_evidence"]).split(";")))
        for cycle_id, cycle_label in KONEHS_CYCLES.items():
            has = cycle_id in public_cycles
            rows.append({
                "test_id": item["test_id"],
                "nhanes_variable": item["nhanes_variable"],
                "domain": item["domain"],
                "konehs_cycle": cycle_id,
                "konehs_cycle_label": cycle_label,
                "public_exposure_evidence": "yes" if has else "no_or_not_confirmed",
                "match_class": item["konehs_match_class"] if has else "not_confirmed_for_cycle",
                "evidence_basis": item["evidence_basis"] if has else "",
                "exact_variable_confirmed": "not_yet",
                "exact_matrix_confirmed": "not_yet",
                "weight_and_design_confirmed": "not_yet",
                "joint_exposure_t2d_confirmed": "not_yet",
                "access_status": "controlled_or_request_to_NIER",
            })
    return rows


def build_outcome_audit() -> list[dict[str, object]]:
    specs = [
        ("T2D_definition", "diabetes outcome", "conditional", "Published KoNEHS analyses document self-report/medication and/or HbA1c-based diabetes definitions, but exact field coding and type-1 exclusion must be confirmed before calling it T2D."),
        ("physician_diagnosis", "diabetes outcome", "documented_by_precedent", "Cycle-2 PAH analysis used physician diagnosis as one component of diabetes ascertainment; this is feasibility precedent, not our result."),
        ("medication_or_insulin", "diabetes outcome", "documented_by_precedent", "Cycle-2 PAH analysis used oral hypoglycemic medication or insulin; exact current coding needs the data dictionary."),
        ("HbA1c", "diabetes outcome", "documented_by_precedent", "Cycle-4 diabetes-related analysis reports HbA1c/diabetes variables; field names and missingness need confirmation."),
        ("type1_exclusion", "diabetes outcome", "not_confirmed", "Cycle-2 publication explicitly did not distinguish type 1 from type 2; a reproducible T2D restriction remains pending."),
        ("age_sex_demographics", "covariate", "documented", "KoNEHS questionnaires and adult analyses include demographic variables; exact names/eligibility pending."),
        ("BMI_anthropometry", "covariate", "documented", "Physical examination and anthropometric measures are documented in cycle-4 overview."),
        ("smoking_alcohol_activity", "covariate", "documented", "Lifestyle and health-behavior questionnaire domains are documented; exact harmonized fields pending."),
        ("SES", "covariate", "documented", "Socioeconomic and demographic questionnaire domains are documented; exact fields pending."),
        ("urinary_creatinine", "laboratory/QC", "documented", "Urinary creatinine and creatinine adjustment are documented for urinary chemicals."),
        ("survey_weights", "survey design", "methodology_documented_names_pending", "Weights incorporate design, nonresponse, and post-stratification; exact variable names and cycle-specific usage are pending."),
        ("strata_psu", "survey design", "methodology_documented_names_pending", "Two-stage stratified sampling and survey analysis are documented; exact strata/PSU identifiers are pending."),
        ("laboratory_lod_qc", "laboratory/QC", "documented", "LOD and laboratory QC procedures are documented; exact analyte files and flags require data dictionary."),
        ("same_person_joint_extract", "access/executability", "not_confirmed", "Public sources do not provide a directly downloadable individual-level joint exposure–T2D extract for this audit."),
        ("individual_level_access", "access/executability", "controlled", "Published data-availability statement directs requests to the relevant Korean environmental-health authority; not public-download executable."),
        ("published_analysis_precedent", "overall feasibility", "yes_but_not_replication", "Independent papers demonstrate that KoNEHS exposure, diabetes, covariates, and survey analysis can be linked by authorized analysts; their estimates are not imported."),
    ]
    rows = []
    for key, domain, status, note in specs:
        rows.append({
            "audit_item": key,
            "domain": domain,
            "status": status,
            "evidence_sources": "K4_EXPOSURE_STATUS;K4_PHTHALATE_DIABETES;K2_PAH_DIABETES;DATA_ACCESS_NOTE",
            "required_next_action": "Obtain current KoNEHS data dictionary/access approval and validate exact fields in an authorized extract." if status not in {"documented", "yes_but_not_replication"} else "Retain as feasibility evidence; verify exact coding at access stage.",
            "notes": note,
        })
    return rows


def main() -> int:
    if not FROZEN.exists():
        raise FileNotFoundError(FROZEN)
    frozen = read_csv(FROZEN)
    if len(frozen) != 29:
        raise ValueError(f"Expected exactly 29 frozen tests, got {len(frozen)}")
    crosswalk = build_crosswalk(frozen)
    cycle_rows = build_cycle_readiness(crosswalk)
    outcome_rows = build_outcome_audit()
    source_snapshot = fetch_source_snapshot()

    write_csv(OUT / "konehs_29_test_crosswalk.csv", crosswalk)
    write_csv(OUT / "konehs_cycle_readiness.csv", cycle_rows)
    write_csv(OUT / "konehs_outcome_covariate_design_audit.csv", outcome_rows)
    (OUT / "KONEHS_SOURCE_SNAPSHOT.json").write_text(json.dumps(source_snapshot, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    exact = sum(row["konehs_match_class"] == "exact_public" for row in crosswalk)
    family = sum(row["konehs_match_class"] == "family_public" for row in crosswalk)
    public = exact + family
    cycle_public = sum(row["public_exposure_evidence"] == "yes" for row in cycle_rows)
    qc = {
        "generated_utc": RUN_UTC,
        "status": "analysis_feasible_by_precedent_but_access_controlled",
        "association_models_run": False,
        "published_effect_estimates_imported": False,
        "outcome_free_crosswalk": True,
        "frozen_test_count_expected": 29,
        "frozen_test_count_observed": len(frozen),
        "crosswalk_rows": len(crosswalk),
        "cycle_rows": len(cycle_rows),
        "outcome_covariate_design_rows": len(outcome_rows),
        "public_exact_match_count": exact,
        "public_family_match_count": family,
        "public_content_match_count": public,
        "public_content_match_fraction": round(public / 29, 6),
        "cycle_rows_with_public_exposure_evidence": cycle_public,
        "exact_variable_confirmed_count": 0,
        "same_person_joint_exposure_t2d_confirmed_count": 0,
        "weight_and_design_variable_confirmed_count": 0,
        "controlled_data_boundary": True,
        "decision": "high_priority_conditional_not_access_confirmed",
        "limitations": [
            "Public papers/overviews support a feasibility floor, not exact variable-level coverage.",
            "Published diabetes definitions are not automatically equivalent to reproducible T2D coding.",
            "Exact survey weight/stratum/PSU names require the current data dictionary or authorized extract.",
            "No KoNEHS association results were run or adopted in this audit.",
        ],
    }
    (OUT / "KONEHS_QC_SUMMARY.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    report = f"""# Step 10B-P — KoNEHS population-replacement feasibility audit

Generated: `{RUN_UTC}`  \\
Status: **{qc['status']}**  \\
Decision: **{qc['decision']}**

## Executive result

The audit supports KoNEHS as a **high-priority conditional population-replacement candidate**, not as an already accessible replication dataset. The audited public evidence gives a conservative content-level floor of **{public}/29 frozen tests** ({exact} exact analyte-level public matches and {family} family/matrix-qualified matches). This is not exact variable-level confirmation.

The audit confirms **29/29 tests were carried into the crosswalk**, but it confirms **0/29 exact public variable names**, **0/29 same-person exposure–T2D extracts**, and **0/29 exact weight/design-variable sets**. No exposure–T2D association model was run.

## Why KoNEHS is promising

KoNEHS is a national environmental-health biomonitoring program. The cycle-4 exposure overview documents a nationally sampled survey, environmental chemicals/metabolites, blood and urine biospecimens, urinary creatinine, laboratory QC, and design/nonresponse/post-stratification weights analyzed with a multistage survey procedure. Published cycle-2 and cycle-4 papers demonstrate that environmental biomarkers, diabetes-related outcomes, covariates, and survey analysis can be linked by authorized analysts.

The source-verified exposure floor includes three PFAS, BPA, six cycle-4 phthalate metabolites, cycle-3 MCOP exposure evidence, four grouped PAH-family matches from cycle 2, and conservative public metal evidence. The full row-level rationale is in `konehs_29_test_crosswalk.csv`; cycle-by-cycle status is in `konehs_cycle_readiness.csv`.

## T2D and design boundary

Published KoNEHS diabetes analyses establish operational feasibility, but not one universally frozen T2D definition. The cycle-2 PAH precedent used physician diagnosis and diabetes medication/insulin and explicitly did not distinguish type 1 from type 2. Therefore the T2D outcome remains **conditional** until the current data dictionary or authorized extract confirms diagnosis, medication, HbA1c/glucose fields, type-1 handling, missingness, and the exact analytic rule.

The same boundary applies to core covariates, urinary creatinine, survey weights, strata, and PSU identifiers: methodology is documented, but exact variable names and cycle-specific construction are pending.

## Access status

The publicly available sources audited here do not provide a direct, unrestricted individual-level joint exposure–T2D extract. The published data-availability boundary directs requests to the relevant Korean environmental-health authority. Thus KoNEHS is **analysis-feasible by precedent but access-controlled**. It should be promoted to primary population replacement only after an approved data request and exact data-dictionary crosswalk.

## Frozen exclusions

- No KoNEHS association result was computed.
- Published PFAS–diabetes, PAH–diabetes, and exposure papers are feasibility precedents only; their estimates are not this project's replication.
- No candidate was selected because it had a published KoNEHS result.
- The 29-test family was read from the existing frozen Step 4 file and was not altered.

## Gate for promotion

KoNEHS may become a primary epidemiologic replacement only if an authorized extract confirms: exact analyte and matrix; same-person exposure plus diabetes outcome; reproducible T2D coding; age/sex/BMI/smoking/alcohol/SES and other prespecified covariates; survey weights and design variables; and adequate access/permission for reproducible analysis.

## Files

- `konehs_29_test_crosswalk.csv`
- `konehs_cycle_readiness.csv`
- `konehs_outcome_covariate_design_audit.csv`
- `KONEHS_SOURCE_SNAPSHOT.json`
- `KONEHS_QC_SUMMARY.json`
- `STEP10B_P_KONEHS_MANIFEST.json`
"""
    (OUT / "STEP10B_P_KONEHS_AUDIT_REPORT.md").write_text(report, encoding="utf-8")

    manifest_files = [
        "PLAN.md",
        "run_step10b_p_konehs_audit.py",
        "konehs_29_test_crosswalk.csv",
        "konehs_cycle_readiness.csv",
        "konehs_outcome_covariate_design_audit.csv",
        "KONEHS_SOURCE_SNAPSHOT.json",
        "KONEHS_QC_SUMMARY.json",
        "STEP10B_P_KONEHS_AUDIT_REPORT.md",
    ]
    manifest = {
        "generated_utc": RUN_UTC,
        "script": "run_step10b_p_konehs_audit.py",
        "frozen_input": str(FROZEN.relative_to(ROOT)).replace("\\", "/"),
        "frozen_test_count": len(frozen),
        "models_run": False,
        "published_effect_estimates_imported": False,
        "controlled_microdata_downloaded": False,
        "source_ids": list(SOURCES),
        "outputs": [],
    }
    for name in manifest_files:
        path = OUT / name
        data = path.read_bytes()
        manifest["outputs"].append({"path": name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    (OUT / "STEP10B_P_KONEHS_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
