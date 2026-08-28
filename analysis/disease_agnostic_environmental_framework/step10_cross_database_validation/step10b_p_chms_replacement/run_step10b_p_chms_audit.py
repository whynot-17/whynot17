#!/usr/bin/env python3
"""Step 10B-P: CHMS population-replacement feasibility audit.

This is a documentation/data-dictionary audit, not an exposure--T2D analysis.
It deliberately separates public content-level evidence from exact variable-level
and person-level access claims. No disease association result is read.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[4]
TESTSET = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
OUTDIR = Path(__file__).resolve().parent

SOURCE_URLS = {
    "content_summary_cycles_1_6": "https://www.statcan.gc.ca/en/statistical-programs/document/5071_D9_V2",
    "data_dictionaries_and_documents": "https://www.statcan.gc.ca/en/statistical-programs/document/5071_D4_V3",
    "cycle_1_user_guide": "https://www.statcan.gc.ca/en/statistical-programs/document/5071_D2_T1_V1",
    "accessing_chms_information_online": "https://www.statcan.gc.ca/en/statistical-programs/document/5071_D5_V2",
    "cycle_1_6_environmental_lab_release": "https://www150.statcan.gc.ca/n1/daily-quotidien/191113/dq191113a-eng.pdf",
    "cycle_7_environmental_lab_release": "https://www150.statcan.gc.ca/n1/daily-quotidien/260410/dq260410b-eng.htm",
}

CHMS_CYCLES = {
    "C1": {"label": "2007-2009", "start": 2007, "end": 2009},
    "C2": {"label": "2009-2011", "start": 2009, "end": 2011},
    "C3": {"label": "2012-2013", "start": 2012, "end": 2013},
    "C4": {"label": "2014-2015", "start": 2014, "end": 2015},
    "C5": {"label": "2016-2017", "start": 2016, "end": 2017},
    "C6": {"label": "2018-2019", "start": 2018, "end": 2019},
}

# The public cycles 1–6 content table uses the following environment rows.
# "exact" means the named analyte is listed in the stated matrix. "family" is
# reserved for grouped NHANES PAH variables whose individual CHMS analytes are
# listed but are not a one-variable equivalent. All exact variable names and
# weights remain pending the cycle-specific dictionaries.
EXPOSURE_EVIDENCE = {
    "LBXPFDE": {
        "chms_measure": "perfluorodecanoic acid (PFDA)", "match_class": "exact",
        "matrix_status": "exact blood", "cycles": ["C2", "C5", "C6"],
        "age_scope": "C2 12-79; C5-C6 3-79", "domain": "PFAS",
        "note": "Named in CHMS environmental exposure table; age/subsample restrictions apply.",
    },
    "LBXPFHS": {
        "chms_measure": "perfluorohexane sulfonate (PFHxS)", "match_class": "exact",
        "matrix_status": "exact blood", "cycles": ["C1", "C2", "C5", "C6"],
        "age_scope": "C1 20-79; C2 12-79; C5-C6 3-79", "domain": "PFAS",
        "note": "Named in CHMS environmental exposure table; age/subsample restrictions apply.",
    },
    "LBXPFNA": {
        "chms_measure": "perfluorononanoic acid (PFNA)", "match_class": "exact",
        "matrix_status": "exact blood", "cycles": ["C2", "C5", "C6"],
        "age_scope": "C2 12-79; C5-C6 3-79", "domain": "PFAS",
        "note": "Named in CHMS environmental exposure table; age/subsample restrictions apply.",
    },
    "URXBPH": {
        "chms_measure": "bisphenol A (BPA)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C1", "C2", "C3", "C4", "C5", "C6"],
        "age_scope": "C1 6-79; C2-C6 3-79, with subsample footnotes",
        "domain": "bisphenol", "note": "BPA is explicitly listed; analogues are not treated as equivalent.",
    },
    "URXCOP": {
        "chms_measure": "mono-(carboxyisooctyl) phthalate (MCIOP)", "match_class": "exact",
        "matrix_status": "exact urine; naming variant of MCOP", "cycles": ["C5", "C6"],
        "age_scope": "3-79 with environmental subsample footnote", "domain": "phthalate",
        "note": "Public table names MCIOP; exact CHMS variable and harmonization must be verified.",
    },
    "URXECP": {
        "chms_measure": "mono-(2-ethyl-5-carboxypentyl) phthalate (MECPP)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C5", "C6"],
        "age_scope": "3-79 with environmental subsample footnote", "domain": "phthalate",
        "note": "Named in CHMS phthalate-metabolite table; exact variable and weight pending.",
    },
    "URXMBP": {
        "chms_measure": "mono-n-butyl phthalate (MnBP)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C1", "C2", "C5", "C6"],
        "age_scope": "C1 6-49; C2-C6 3-79 with subsample footnotes", "domain": "phthalate",
        "note": "Named in CHMS phthalate-metabolite table; age/subsample restrictions apply.",
    },
    "URXMEP": {
        "chms_measure": "mono ethyl phthalate (MEP)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C1", "C2", "C5", "C6"],
        "age_scope": "C1 6-49; C2-C6 3-79 with subsample footnotes", "domain": "phthalate",
        "note": "Named in CHMS phthalate-metabolite table; age/subsample restrictions apply.",
    },
    "URXMHH": {
        "chms_measure": "mono-(2-ethyl-5-hydroxyhexyl) phthalate (MEHHP)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C1", "C2", "C5", "C6"],
        "age_scope": "C1 6-49; C2-C6 3-79 with subsample footnotes", "domain": "phthalate",
        "note": "Named in CHMS phthalate-metabolite table; age/subsample restrictions apply.",
    },
    "URXMHP": {
        "chms_measure": "mono-2-ethylhexyl phthalate (MEHP)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C1", "C2", "C5", "C6"],
        "age_scope": "C1 6-49; C2-C6 3-79 with subsample footnotes", "domain": "phthalate",
        "note": "Named in CHMS phthalate-metabolite table; age/subsample restrictions apply.",
    },
    "URXMIB": {
        "chms_measure": "mono-iso-butyl phthalate (MiBP)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C2", "C5", "C6"],
        "age_scope": "C2 3-79; C5-C6 3-79 with subsample footnotes", "domain": "phthalate",
        "note": "Named in CHMS phthalate-metabolite table; age/subsample restrictions apply.",
    },
    "URXMOH": {
        "chms_measure": "mono-(2-ethyl-5-oxohexyl) phthalate (MEOHP)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C1", "C2", "C5", "C6"],
        "age_scope": "C1 6-49; C2-C6 3-79 with subsample footnotes", "domain": "phthalate",
        "note": "Named in CHMS phthalate-metabolite table; age/subsample restrictions apply.",
    },
    "URXMZP": {
        "chms_measure": "mono benzyl phthalate (MBzP)", "match_class": "exact",
        "matrix_status": "exact urine", "cycles": ["C1", "C2", "C5", "C6"],
        "age_scope": "C1 6-49; C2-C6 3-79 with subsample footnotes", "domain": "phthalate",
        "note": "Named in CHMS phthalate-metabolite table; age/subsample restrictions apply.",
    },
    "URXP02": {
        "chms_measure": "1-/2-hydroxynaphthalene within PAH panel", "match_class": "family",
        "matrix_status": "urine family match, not one-variable equivalence", "cycles": ["C2", "C3", "C4"],
        "age_scope": "3-79 with environmental subsample footnote", "domain": "PAH",
        "note": "CHMS lists hydroxynaphthalene analytes; NHANES URXP02 is a grouped test family.",
    },
    "URXP04": {
        "chms_measure": "2-/3-/9-hydroxyfluorene within PAH panel", "match_class": "family",
        "matrix_status": "urine family match, not one-variable equivalence", "cycles": ["C2", "C3", "C4"],
        "age_scope": "3-79 with environmental subsample footnote", "domain": "PAH",
        "note": "CHMS lists hydroxyfluorene analytes; NHANES URXP04 is a grouped test family.",
    },
    "URXP10": {
        "chms_measure": "1-hydroxypyrene and 3-hydroxybenzo(a)pyrene within PAH panel", "match_class": "family",
        "matrix_status": "urine family match, not one-variable equivalence", "cycles": ["C2", "C3", "C4"],
        "age_scope": "3-79 with environmental subsample footnote", "domain": "PAH",
        "note": "CHMS lists pyrene-related OH-PAHs; exact equivalence to grouped NHANES URXP10 is pending.",
    },
    "URXP25": {
        "chms_measure": "1-/2-/3-/4-/9-hydroxyphenanthrene within PAH panel", "match_class": "family",
        "matrix_status": "urine family match, not one-variable equivalence", "cycles": ["C2", "C3", "C4"],
        "age_scope": "3-79 with environmental subsample footnote", "domain": "PAH",
        "note": "CHMS lists hydroxyphenanthrenes; exact equivalence to grouped NHANES URXP25 is pending.",
    },
    "URXUBA": {
        "chms_measure": "barium listed in hair, not confirmed in urine in cycles 1-6 summary", "match_class": "none",
        "matrix_status": "matrix mismatch / no confirmed urine match", "cycles": [],
        "age_scope": "CHMS cycle 6 hair barium 20-59", "domain": "metals",
        "note": "Do not substitute hair barium for the frozen urine test without a variable-level source.",
    },
    "URXUCD": {
        "chms_measure": "cadmium", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C1", "C2", "C5", "C6"], "age_scope": "C1 6-79; C2 3-79; C5-C6 3-79",
        "domain": "metals", "note": "Named in CHMS metals table; cycle-specific environmental urine weight pending.",
    },
    "URXUCO": {
        "chms_measure": "cobalt", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C2"], "age_scope": "3-79", "domain": "metals",
        "note": "Named in CHMS metals table for cycle 2 only.",
    },
    "URXUCS": {
        "chms_measure": "cesium", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C2"], "age_scope": "3-79", "domain": "metals",
        "note": "Named in CHMS metals table for cycle 2 only.",
    },
    "URXUMO": {
        "chms_measure": "molybdenum", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C1", "C2"], "age_scope": "C1 6-79; C2 3-79", "domain": "metals",
        "note": "Named in CHMS metals table; exact variable and weight pending.",
    },
    "URXUPB": {
        "chms_measure": "lead", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C1", "C2"], "age_scope": "C1 6-79; C2 3-79", "domain": "metals",
        "note": "Named in CHMS metals table; exact variable and weight pending.",
    },
    "URXUSB": {
        "chms_measure": "antimony", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C1", "C2"], "age_scope": "C1 6-79; C2 3-79", "domain": "metals",
        "note": "Named in CHMS metals table; exact variable and weight pending.",
    },
    "URXUSN": {
        "chms_measure": "tin not confirmed in cycles 1-6 public content summary", "match_class": "none",
        "matrix_status": "no confirmed urine match", "cycles": [], "age_scope": "not confirmed",
        "domain": "metals", "note": "Require cycle-specific dictionary search before treating as CHMS-available.",
    },
    "URXUSR": {
        "chms_measure": "silver", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C2"], "age_scope": "3-79", "domain": "metals",
        "note": "Named in CHMS metals table for cycle 2 only.",
    },
    "URXUTL": {
        "chms_measure": "thallium", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C2"], "age_scope": "3-79", "domain": "metals",
        "note": "Named in CHMS metals table for cycle 2 only.",
    },
    "URXUTU": {
        "chms_measure": "tungsten", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C2"], "age_scope": "3-79", "domain": "metals",
        "note": "Named in CHMS metals table for cycle 2 only.",
    },
    "URXUUR": {
        "chms_measure": "uranium", "match_class": "exact", "matrix_status": "exact urine",
        "cycles": ["C1", "C2"], "age_scope": "C1 6-79; C2 3-79", "domain": "metals",
        "note": "Named in CHMS metals table; exact variable and weight pending.",
    },
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch_source_metadata() -> dict:
    result = {}
    for key, url in SOURCE_URLS.items():
        entry = {"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(), "status": "not_attempted"}
        try:
            req = Request(url, headers={"User-Agent": "whynot17-CHMS-audit/1.0"})
            # Metadata capture is best-effort and bounded. The audit must not
            # hang on a government web endpoint or silently turn missing source
            # metadata into an analyte match.
            with urlopen(req, timeout=5) as response:
                body = response.read(2_000_000)
                entry.update({
                    "status": "ok",
                    "http_status": getattr(response, "status", None),
                    "content_type": response.headers.get("Content-Type"),
                    "bytes_read_for_hash": len(body),
                    "read_cap_bytes": 2_000_000,
                    "sha256": sha256_bytes(body),
                })
        except Exception as exc:  # metadata failure must not fabricate coverage
            entry.update({"status": "fetch_error", "error": f"{type(exc).__name__}: {exc}"})
        result[key] = entry
    return result


def read_testset() -> list[dict]:
    with TESTSET.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 29:
        raise RuntimeError(f"Expected 29 frozen tests, found {len(rows)}")
    variables = {row["variable"] for row in rows}
    missing = variables - set(EXPOSURE_EVIDENCE)
    if missing:
        raise RuntimeError(f"Missing CHMS evidence mapping for frozen variables: {sorted(missing)}")
    return rows


def parse_years(cycles_field: str) -> set[int]:
    years: set[int] = set()
    for item in (cycles_field or "").split(";"):
        item = item.strip()
        if not item:
            continue
        start, end = item.split("-")
        years.update(range(int(start), int(end) + 1))
    return years


def calendar_overlap(test_years: set[int], chms_cycle_ids: list[str]) -> list[str]:
    return [cid for cid in chms_cycle_ids if test_years.intersection(range(CHMS_CYCLES[cid]["start"], CHMS_CYCLES[cid]["end"] + 1))]


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise RuntimeError(f"No rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_crosswalk(test_rows: list[dict]) -> tuple[list[dict], list[dict]]:
    crosswalk = []
    cycle_rows = []
    for row in test_rows:
        variable = row["variable"]
        evidence = EXPOSURE_EVIDENCE[variable]
        content_cycles = evidence["cycles"]
        overlaps = calendar_overlap(parse_years(row["cycles"]), content_cycles)
        matrix = row["matrix"]
        required_weight = "environmental_blood_or_urine_subsample_weight" if matrix in {"serum_or_blood", "blood"} else "environmental_urine_subsample_weight"
        crosswalk.append({
            "test_id": row["test_id"],
            "nhanes_biomarker": row["biomarker"],
            "nhanes_variable": variable,
            "nhanes_matrix": matrix,
            "nhanes_frozen_cycles": row["cycles"],
            "domain": evidence["domain"],
            "chms_public_measure": evidence["chms_measure"],
            "chms_match_class": evidence["match_class"],
            "chms_matrix_status": evidence["matrix_status"],
            "chms_public_content_cycles": ";".join(content_cycles),
            "chms_public_content_cycle_labels": ";".join(CHMS_CYCLES[cid]["label"] for cid in content_cycles),
            "calendar_window_overlap_with_nhanes": ";".join(overlaps),
            "calendar_overlap_is_not_same_cycle": "True",
            "age_scope_from_public_content": evidence["age_scope"],
            "exact_chms_variable_confirmed": "not_yet",
            "exact_chms_weight_confirmed": "not_yet",
            "same_person_joint_exposure_t2d_confirmed": "not_yet",
            "required_weight_class": required_weight,
            "access_status": "controlled_microdata_or_RDC; not public raw microdata",
            "audit_note": evidence["note"],
        })
        for cid, details in CHMS_CYCLES.items():
            listed = cid in content_cycles
            cycle_rows.append({
                "test_id": row["test_id"],
                "nhanes_variable": variable,
                "chms_cycle": cid,
                "chms_cycle_label": details["label"],
                "public_content_status": "listed" if listed else "not_listed_or_not_confirmed",
                "content_match_class": evidence["match_class"] if listed else "none",
                "exact_variable_status": "pending_cycle_dictionary",
                "weight_status": "pending_cycle_dictionary",
                "joint_t2d_covariate_status": "pending_individual_microdata",
                "interpretation": "Content table evidence only; not proof of analyzable person-level overlap." if listed else "No named analyte evidence in the cycles 1-6 public content summary.",
            })
    return crosswalk, cycle_rows


def build_outcome_audit() -> list[dict]:
    rows = [
        ("T2D self-reported diagnosis", "household questionnaire / chronic conditions", "C1-C6", "present_in_public_documentation", "not_locked", "Define diagnosed diabetes and inspect type 1/type 2 fields in cycle dictionaries."),
        ("Glucose", "blood; plasma/serum depending cycle", "C1-C6", "present_in_public_content_summary", "definition_pending", "Confirm fasting/random status, units, derived variables, and eligible age by cycle."),
        ("HbA1c", "blood", "C1-C6", "present_in_public_content_summary", "definition_pending", "Confirm field names, assay metadata, and harmonized threshold rule."),
        ("Fasting insulin", "blood; fasting subgroup", "C1-C6", "present_in_public_content_summary", "secondary_only_until_confirmed", "Subsample restriction; do not require it for the primary definition without a prespecified rule."),
        ("Age", "master/MEC questionnaire", "C1-C6", "present_and_core", "likely_available", "Confirm exact age-at-MEC field and age restrictions after dictionary retrieval."),
        ("Sex/gender", "master/household", "C1-C6", "present_and_core", "likely_available", "Use the CHMS-defined field; do not assume NHANES coding equivalence."),
        ("Anthropometry/BMI", "MEC physical measures", "C1-C6", "present_and_core", "likely_available", "Confirm derived BMI field and valid measurement flags."),
        ("Smoking", "household questionnaire", "C1-C6", "documented_domain", "pending_dictionary", "Confirm current/former/never coding and age eligibility."),
        ("Alcohol", "household questionnaire", "C1-C6", "documented_domain", "pending_dictionary", "Confirm frequency/quantity variables and missingness."),
        ("Physical activity", "household/MEC questionnaire", "C1-C6", "documented_domain", "pending_dictionary", "Confirm harmonized measure across cycles."),
        ("Education/household income/SES", "household questionnaire", "C1-C6", "documented_domain", "pending_dictionary", "Confirm common SES representation and release-specific suppression rules."),
        ("Urine creatinine", "urine laboratory component", "C1-C6", "present_in_public_documentation", "pending_dictionary", "Confirm co-release with each environmental urine file and exact units."),
        ("Complex-survey weights", "full-sample and subsample weights", "C1-C6", "confirmed_methodology_level", "not_ready_for_model", "Use the environmental subsample weight for the exposure file; exact variable and cycle-combination instructions required."),
        ("Variance/design variables", "strata/PSU or bootstrap replicate design", "C1-C6", "confirmed_methodology_level", "not_ready_for_model", "Retrieve cycle-specific design variables and combining/replicate-weight instructions."),
        ("Individual-level access", "CHMS microdata", "C1-C6", "controlled_access", "conditional", "Public dashboards/content tables are insufficient for a joint exposure-T2D regression."),
    ]
    return [{
        "domain": domain,
        "source_or_file_domain": source,
        "cycles": cycles,
        "public_documentation_status": public_status,
        "analysis_readiness": readiness,
        "required_next_action": action,
        "outcome_firewall_status": "not_applicable_to_population_replacement_audit",
    } for domain, source, cycles, public_status, readiness, action in rows]


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    test_rows = read_testset()
    crosswalk, cycle_rows = build_crosswalk(test_rows)
    outcome_rows = build_outcome_audit()

    crosswalk_path = OUTDIR / "chms_29_test_crosswalk.csv"
    cycle_path = OUTDIR / "chms_cycle_readiness.csv"
    outcome_path = OUTDIR / "chms_outcome_covariate_design_audit.csv"
    write_csv(crosswalk_path, crosswalk)
    write_csv(cycle_path, cycle_rows)
    write_csv(outcome_path, outcome_rows)

    source_payload = {
        "audit": "STEP10B_P_CHMS_REPLACEMENT_AUDIT",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "source_scope": "Official Statistics Canada public documentation and content tables; no restricted microdata retrieved.",
        "content_evidence_anchors": {
            "environmental_exposure_table": "Content summary cycles 1-6, Table 3k, environmental exposure section",
            "diabetes_table": "Content summary cycles 1-6, Table 3f, diabetes section",
            "documents_and_dictionaries": "CHMS documents page: separate environment lab blood and urine, environment urine main subsample, and cycle-combination documentation",
            "weight_methodology": "CHMS Cycle 1 user guide: subsample file weights are used when a subsample is linked to the master file",
        },
        "sources": fetch_source_metadata(),
        "interpretation": {
            "exact": "Named analyte appears in the stated matrix in the public cycles 1-6 content summary.",
            "family": "Named component analytes appear, but the frozen NHANES grouped test is not a one-variable equivalent.",
            "none": "No matching analyte/matrix was confirmed from the public cycles 1-6 content summary.",
            "calendar_overlap": "Calendar-window overlap is descriptive and does not establish same-cycle comparability.",
        },
    }
    source_path = OUTDIR / "CHMS_SOURCE_SNAPSHOT.json"
    write_json(source_path, source_payload)

    counts = {}
    for row in crosswalk:
        counts[row["chms_match_class"]] = counts.get(row["chms_match_class"], 0) + 1
    exact_count = counts.get("exact", 0)
    family_count = counts.get("family", 0)
    not_confirmed_count = counts.get("none", 0)
    public_listed_cycle_rows = sum(row["public_content_status"] == "listed" for row in cycle_rows)
    qc = {
        "lock_type": "STEP10B_P_CHMS_POPULATION_REPLACEMENT_AUDIT",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "high_priority_conditional_not_access_confirmed",
        "audit_only": True,
        "no_association_models_run": True,
        "checks": {
            "frozen_test_count": {"observed": len(test_rows), "expected": 29, "pass": len(test_rows) == 29},
            "crosswalk_row_count": {"observed": len(crosswalk), "expected": 29, "pass": len(crosswalk) == 29},
            "cycle_audit_row_count": {"observed": len(cycle_rows), "expected": 174, "pass": len(cycle_rows) == 174},
            "outcome_covariate_design_domains": {"observed": len(outcome_rows), "expected_min": 10, "pass": len(outcome_rows) >= 10},
            "public_content_listed_test_cycle_rows": {"observed": public_listed_cycle_rows, "expected_min": 1, "pass": public_listed_cycle_rows > 0},
        },
        "coverage_summary": {
            "exact_named_analyte_or_naming_variant": exact_count,
            "near_exact_grouped_family": family_count,
            "not_confirmed_in_public_cycles_1_6_summary": not_confirmed_count,
            "exact_or_near_exact_content_level": exact_count + family_count,
            "exact_variable_level_confirmed": 0,
            "same_person_joint_exposure_t2d_confirmed": 0,
        },
        "frozen_gates": {
            "exact_exposure_variables": "pending_cycle_specific_data_dictionaries",
            "t2d_definition": "publicly supported but not yet variable-level locked",
            "core_covariates": "domains appear available; harmonized fields pending",
            "subsample_weights": "methodology confirmed; exact exposure-file weight fields pending",
            "complex_design": "methodology confirmed; exact design/replicate fields pending",
            "individual_access": "controlled_access_required",
        },
        "source_files": {
            "frozen_testset": str(TESTSET.relative_to(ROOT)),
            "crosswalk": crosswalk_path.name,
            "cycle_readiness": cycle_path.name,
            "outcome_audit": outcome_path.name,
            "source_snapshot": source_path.name,
        },
    }
    qc_path = OUTDIR / "CHMS_QC_SUMMARY.json"
    write_json(qc_path, qc)

    report = f"""# Step 10B-P — CHMS population-replacement feasibility audit

Generated: `{qc['generated_utc']}`
Status: **`{qc['status']}`**

## Scope

This audit asks whether CHMS is a viable independent population source for the
frozen 29 human biomarker tests and a T2D analysis. It does **not** inspect any
exposure–T2D association result, change the 29-test family, or promote a
candidate.

## Executive result

The official cycles 1–6 content summary supports **{qc['coverage_summary']['exact_or_near_exact_content_level']}/29**
tests at the content level: **{qc['coverage_summary']['exact_named_analyte_or_naming_variant']}** named analytes or naming
variants and **{qc['coverage_summary']['near_exact_grouped_family']}** grouped PAH-family matches.
**{qc['coverage_summary']['not_confirmed_in_public_cycles_1_6_summary']}** tests are not confirmed in the public
cycles 1–6 summary (`URXUBA` urine barium and `URXUSN` urine tin). These are
not treated as negative biological findings.

The audit confirms **0/29 exact CHMS variable-level mappings** and **0/29
person-level joint exposure–T2D confirmations**, because those require the
cycle-specific data dictionaries and controlled individual-level files. The
correct current disposition is therefore **high-priority conditional**, not
primary replacement ready.

## Exposure crosswalk

The full 29-row crosswalk is in `chms_29_test_crosswalk.csv`; the per-cycle
174-row audit is in `chms_cycle_readiness.csv`. The public content table shows:

- PFDA, PFHxS, PFNA in blood, with cycle- and age-specific restrictions;
- BPA in urine across cycles 1–6;
- nine named phthalate metabolites, including MCIOP/MCiOP as the public naming
  variant relevant to the frozen MCOP mapping;
- four PAH grouped-family matches based on named hydroxylated PAHs;
- ten named metals in urine, mostly concentrated in cycles 1–2, while barium
  is listed in hair rather than urine and tin was not confirmed.

Calendar-window overlap with NHANES is reported descriptively only. It is not
treated as same-cycle equivalence.

## T2D, covariates, and survey design

Public CHMS documentation lists glucose and HbA1c in the diabetes laboratory
domain across cycles 1–6, and also provides a self-reported chronic-condition
questionnaire domain. A reproducible T2D definition is **not yet frozen**:
fasting/random glucose, HbA1c, self-report, type 1 exclusion, age eligibility,
and missingness rules must be resolved from the cycle dictionaries before any
model is specified.

Age, sex/gender, anthropometry/BMI, smoking, alcohol, physical activity,
education/income/SES, and urine creatinine are plausible/common covariate
domains. The CHMS user guide and documents page confirm that full-sample and
subsample weights and instructions for combining cycles exist. The exact
environmental blood/urine weight variables, design variables or replicate
weights, and cycle-specific joins remain pending.

Individual-level microdata are not supplied by the public content dashboard;
access is controlled. Thus no claim is made that an exposure, T2D outcome,
covariate set, and design variables are already jointly available for the same
respondents.

## Gate decision

**Promotion to primary epidemiologic replacement: not yet passed.** The next
action is a controlled data-dictionary/access request for the exact analyte
variables, environmental subsample weights, design/replicate variables, T2D
definition fields, covariates, and respondent-level linkage. If those gates
pass, CHMS becomes the primary external population replacement; otherwise it
remains a high-priority conditional source.

## Official source boundary

Source URLs and retrieval metadata are in `CHMS_SOURCE_SNAPSHOT.json`. The
main content evidence is the Statistics Canada *Content summary for cycles 1
to 6*; the official documents page explicitly lists separate environmental
blood/urine and environment-urine-main-subsample dictionaries and warns that
not every subsample is present in every cycle. No restricted file was obtained.
"""
    report_path = OUTDIR / "STEP10B_P_CHMS_AUDIT_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    output_files = [crosswalk_path, cycle_path, outcome_path, source_path, qc_path, report_path]
    manifest = {
        "manifest_type": "STEP10B_P_CHMS_REPLACEMENT_MANIFEST",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": qc["status"],
        "no_association_models_run": True,
        "inputs": {"frozen_testset": {"path": str(TESTSET.relative_to(ROOT)), "sha256": sha256_bytes(TESTSET.read_bytes())}},
        "outputs": {path.name: {"bytes": path.stat().st_size, "sha256": sha256_bytes(path.read_bytes())} for path in output_files},
        "official_sources": SOURCE_URLS,
        "interpretation_boundary": "Public content-level feasibility only; exact variables, weights, joint person-level availability, and access remain conditional.",
    }
    write_json(OUTDIR / "STEP10B_P_CHMS_MANIFEST.json", manifest)
    print(json.dumps({"status": qc["status"], "coverage_summary": qc["coverage_summary"], "output_dir": str(OUTDIR)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
