#!/usr/bin/env python3
"""Build the Phase 1 T2D exposure-opportunity candidate and mapping gate.

This phase intentionally stops before literature, mechanism, docking, or target
prioritisation.  The only disease-dependent input is the already-frozen T2D
Step 5 result table; all chemical-to-biomarker mappings are read from the
outcome-free Step 2 mapping table.  The script keeps every mapped chemical in
the master table, while separately flagging proxy or internally inconsistent
mapping rows so they cannot be silently promoted to specific exposures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


SCRIPT_VERSION = "phase1_mapping_gate_v1.0"

POSITIVE_REQUIRED = 14
EXPECTED_FDR_DENOMINATOR = "29"

ELEMENT_BY_BIOMARKER = {
    "URXUBA": "barium",
    "URXUMO": "molybdenum",
    "URXUPB": "lead",
    "URXUSN": "tin",
    "URXUSR": "silver",
    "URXUTU": "tungsten",
    "URXUUR": "uranium",
}

OUTPUT_MASTER_FIELDS = [
    "chemical_id",
    "chemical_name",
    "chemical_class",
    "positive_biomarker",
    "exposure_axis",
    "mapping_type",
    "mapping_confidence",
    "mapping_source",
    "mapping_status",
    "mapping_grade",
    "mapping_gate_status",
    "mapping_gate_reason",
    "manual_review_required",
    "NHANES_variable",
    "matrix",
    "NHANES_lab_component",
    "cycle_list",
    "n_cycles_available",
    "weight_variable",
    "n_measured",
    "n_above_lod",
    "pooled_above_lod_pct",
    "mapped",
    "candidate_is_primary",
    "T2D_OR",
    "T2D_P",
    "T2D_FDR",
    "T2D_analytic_N",
    "T2D_case_N",
    "T2D_status",
    "robustness_status",
    "robustness_loco_direction_fraction",
    "robustness_cycle_direction_fraction",
    "robustness_cycle_discordant",
    "robustness_priority_tier",
]

OUTPUT_UNIQUE_FIELDS = [
    "chemical_id",
    "chemical_name",
    "chemical_class",
    "positive_biomarkers",
    "exposure_axes",
    "mapping_types",
    "mapping_grades",
    "mapping_gate_statuses",
    "mapping_gate_disposition",
    "mapping_gate_reasons",
    "manual_review_required",
    "NHANES_variables",
    "matrices",
    "cycle_lists",
    "weight_variables",
    "mapping_row_count",
    "primary_mapping_row_count",
    "direct_mapping_row_count",
    "specific_mapping_row_count",
    "proxy_mapping_row_count",
    "T2D_ORs",
    "T2D_P_values",
    "T2D_FDR_values",
    "T2D_analytic_Ns",
    "T2D_case_Ns",
    "T2D_robustness_statuses",
    "T2D_robustness_tiers",
]

OUTPUT_EXCLUSION_FIELDS = [
    "chemical_id",
    "chemical_name",
    "chemical_class",
    "positive_biomarker",
    "NHANES_variable",
    "mapping_type",
    "mapping_confidence",
    "mapping_source",
    "mapping_grade",
    "mapping_gate_status",
    "mapping_gate_reason",
    "exclusion_scope",
    "recommended_next_action",
]


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def value(row: Mapping[str, str], key: str, default: str = "") -> str:
    return str(row.get(key, default) or default).strip()


def truthy(raw: str) -> bool:
    return value({"x": raw}, "x").lower() in {"true", "1", "yes", "y"}


def split_semicolon(raw: str) -> List[str]:
    return [item.strip() for item in raw.split(";") if item.strip()]


def norm_name(raw: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", raw.lower())


def unique_join(values: Iterable[str]) -> str:
    cleaned = sorted({value({"x": item}, "x") for item in values if value({"x": item}, "x")})
    return ";".join(cleaned)


def positive_tests(step5_rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    positives = [row for row in step5_rows if truthy(value(row, "FDR_supported"))]
    positives.sort(key=lambda row: float(value(row, "BH_FDR", "inf") or "inf"))
    if len(positives) != POSITIVE_REQUIRED:
        raise RuntimeError(f"Expected {POSITIVE_REQUIRED} T2D FDR-positive tests, found {len(positives)}")
    biomarkers = [value(row, "biomarker") for row in positives]
    if len(set(biomarkers)) != len(biomarkers):
        raise RuntimeError("T2D positive test table contains duplicate biomarkers")
    bad_denominators = sorted({value(row, "fdr_denominator") for row in positives if value(row, "fdr_denominator") != EXPECTED_FDR_DENOMINATOR})
    if bad_denominators:
        raise RuntimeError(f"Positive tests do not all use frozen FDR denominator 29: {bad_denominators}")
    return positives


def grade_mapping(row: Mapping[str, str]) -> Dict[str, str]:
    """Assign a transparent A/B/C grade without using literature or outcomes."""

    biomarker = value(row, "human_biomarker")
    name = value(row, "chemical_name")
    lower_type = value(row, "mapping_type").lower()
    lower_name = name.lower()

    if lower_type == "direct serum pfas analyte":
        return {
            "mapping_grade": "A",
            "mapping_gate_status": "pass",
            "mapping_gate_reason": "Directly measured serum PFAS analyte mapping.",
            "manual_review_required": "False",
        }

    if lower_type == "parent-to-validated-urinary-metabolite axis":
        return {
            "mapping_grade": "B",
            "mapping_gate_status": "pass",
            "mapping_gate_reason": "Specific parent-to-validated urinary metabolite axis recorded in the neutral mapping table.",
            "manual_review_required": "False",
        }

    if lower_type == "direct urinary metabolite":
        return {
            "mapping_grade": "B",
            "mapping_gate_status": "pass_conditional",
            "mapping_gate_reason": "Specific urinary metabolite is directly measured; parent-compound inference is not assumed.",
            "manual_review_required": "True",
        }

    if lower_type == "parent pah to validated urinary oh-pah proxy":
        return {
            "mapping_grade": "C",
            "mapping_gate_status": "conditional",
            "mapping_gate_reason": "Validated urinary OH-PAH proxy with limited parent specificity; not eligible as a compound-specific primary exposure without further evidence.",
            "manual_review_required": "True",
        }

    if lower_type == "elemental urinary biomarker for parent metal/species":
        parent = ELEMENT_BY_BIOMARKER.get(biomarker, "")
        normalized_name = norm_name(name)
        normalized_parent = norm_name(parent)
        if parent and normalized_name == normalized_parent:
            return {
                "mapping_grade": "A",
                "mapping_gate_status": "pass",
                "mapping_gate_reason": f"Chemical name is the measured parent element ({parent}) for the elemental urinary analyte.",
                "manual_review_required": "False",
            }
        if parent and normalized_name.startswith(normalized_parent) and re.search(r"-?\d+$", lower_name):
            return {
                "mapping_grade": "C",
                "mapping_gate_status": "conditional",
                "mapping_gate_reason": f"Isotopic/species-specific name under the {parent} elemental assay; urinary assay cannot establish isotope/species specificity.",
                "manual_review_required": "True",
            }
        # Use token boundaries rather than a bare substring check.  For
        # example, ``tin`` occurs inside ``actinium`` and ``dentinol`` but
        # neither name identifies a tin species.
        has_parent_token = bool(parent and re.search(rf"(?<![a-z]){re.escape(parent)}(?![a-z])", lower_name))
        if parent and has_parent_token:
            return {
                "mapping_grade": "C",
                "mapping_gate_status": "conditional",
                "mapping_gate_reason": f"Chemical name contains the measured element ({parent}) but identifies a salt, formulation, complex, or other species not resolved by the elemental assay.",
                "manual_review_required": "True",
            }
        return {
            "mapping_grade": "C",
            "mapping_gate_status": "exclude",
            "mapping_gate_reason": f"Source mapping type is elemental, but chemical name does not identify the measured element ({parent or biomarker}); retain for audit and exclude from specific-exposure shortlist pending mapping review.",
            "manual_review_required": "True",
        }

    return {
        "mapping_grade": "C",
        "mapping_gate_status": "exclude",
        "mapping_gate_reason": "Unrecognized mapping type; retain in master but exclude from the specific-exposure shortlist pending manual review.",
        "manual_review_required": "True",
    }


def make_mapping_rows(
    mappings: Sequence[Dict[str, str]],
    positives_by_biomarker: Mapping[str, Dict[str, str]],
    robust_by_biomarker: Mapping[str, Dict[str, str]],
) -> List[Dict[str, str]]:
    output: List[Dict[str, str]] = []
    for raw in mappings:
        biomarker = value(raw, "human_biomarker")
        if biomarker not in positives_by_biomarker:
            continue
        if not truthy(value(raw, "mapped")):
            continue
        test = positives_by_biomarker[biomarker]
        robust = robust_by_biomarker.get(biomarker, {})
        grade = grade_mapping(raw)
        row = dict(raw)
        row.update(grade)
        row.update(
            {
                "positive_biomarker": biomarker,
                "T2D_OR": value(test, "OR"),
                "T2D_P": value(test, "P"),
                "T2D_FDR": value(test, "BH_FDR"),
                "T2D_analytic_N": value(test, "analytic_n"),
                "T2D_case_N": value(test, "analytic_t2d_cases"),
                "T2D_status": value(test, "status"),
                "robustness_status": value(robust, "status"),
                "robustness_loco_direction_fraction": value(robust, "loco_direction_fraction"),
                "robustness_cycle_direction_fraction": value(robust, "cycle_direction_fraction"),
                "robustness_cycle_discordant": value(robust, "cycle_discordant"),
                "robustness_priority_tier": value(robust, "priority_tier"),
            }
        )
        output.append(row)
    output.sort(key=lambda row: (float(value(row, "T2D_FDR", "inf") or "inf"), value(row, "positive_biomarker"), value(row, "chemical_name")))
    return output


def make_unique_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[value(row, "chemical_id")].append(row)

    output: List[Dict[str, str]] = []
    grade_rank = {"A": 0, "B": 1, "C": 2}
    for chemical_id, group in sorted(grouped.items(), key=lambda item: (min(float(value(r, "T2D_FDR", "inf") or "inf") for r in item[1]), item[0])):
        grades = sorted({value(row, "mapping_grade") for row in group}, key=lambda item: grade_rank.get(item, 9))
        statuses = sorted({value(row, "mapping_gate_status") for row in group})
        reasons = sorted({value(row, "mapping_gate_reason") for row in group})
        has_a_or_b = any(grade in {"A", "B"} for grade in grades)
        has_direct_pass = any(value(row, "mapping_gate_status") == "pass" for row in group)
        has_conditional_only = not has_direct_pass and any(value(row, "mapping_gate_status") == "conditional" for row in group)
        if has_a_or_b and has_direct_pass:
            disposition = "advance_to_literature_audit"
        elif has_a_or_b:
            disposition = "advance_with_parent_specificity_review"
        elif has_conditional_only:
            disposition = "proxy_only_not_primary"
        else:
            disposition = "exclude_pending_mapping_review"

        output.append(
            {
                "chemical_id": chemical_id,
                "chemical_name": value(group[0], "chemical_name"),
                "chemical_class": unique_join(value(row, "chemical_class") for row in group),
                "positive_biomarkers": unique_join(value(row, "positive_biomarker") for row in group),
                "exposure_axes": unique_join(value(row, "exposure_axis") for row in group),
                "mapping_types": unique_join(value(row, "mapping_type") for row in group),
                "mapping_grades": ";".join(grades),
                "mapping_gate_statuses": unique_join(statuses),
                "mapping_gate_disposition": disposition,
                "mapping_gate_reasons": " || ".join(reasons),
                "manual_review_required": "True" if any(value(row, "manual_review_required") == "True" for row in group) else "False",
                "NHANES_variables": unique_join(value(row, "NHANES_variable") for row in group),
                "matrices": unique_join(value(row, "matrix") for row in group),
                "cycle_lists": " || ".join(sorted({value(row, "cycle_list") for row in group})),
                "weight_variables": unique_join(value(row, "weight_variable") for row in group),
                "mapping_row_count": str(len(group)),
                "primary_mapping_row_count": str(sum(value(row, "candidate_is_primary") == "True" for row in group)),
                "direct_mapping_row_count": str(sum(value(row, "mapping_grade") == "A" for row in group)),
                "specific_mapping_row_count": str(sum(value(row, "mapping_grade") == "B" for row in group)),
                "proxy_mapping_row_count": str(sum(value(row, "mapping_grade") == "C" for row in group)),
                "T2D_ORs": unique_join(value(row, "T2D_OR") for row in group),
                "T2D_P_values": unique_join(value(row, "T2D_P") for row in group),
                "T2D_FDR_values": unique_join(value(row, "T2D_FDR") for row in group),
                "T2D_analytic_Ns": unique_join(value(row, "T2D_analytic_N") for row in group),
                "T2D_case_Ns": unique_join(value(row, "T2D_case_N") for row in group),
                "T2D_robustness_statuses": unique_join(value(row, "robustness_status") for row in group),
                "T2D_robustness_tiers": unique_join(value(row, "robustness_priority_tier") for row in group),
            }
        )
    return output


def make_exclusion_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    excluded = [row for row in rows if value(row, "mapping_gate_status") != "pass"]
    output = []
    for row in excluded:
        status = value(row, "mapping_gate_status")
        if status == "conditional":
            scope = "excluded_from_compound_specific_primary_shortlist"
            next_action = "Require independent parent-specific exposure evidence before treating as a compound-specific lead."
        elif status == "pass_conditional":
            scope = "excluded_from_unqualified_parent_compound_inference"
            next_action = "Audit parent-metabolite interpretation before literature prioritisation."
        else:
            scope = "excluded_from_specific_exposure_shortlist"
            next_action = "Manually review chemical-to-analyte identity; do not infer exposure specificity from the current row alone."
        output.append(
            {
                "chemical_id": value(row, "chemical_id"),
                "chemical_name": value(row, "chemical_name"),
                "chemical_class": value(row, "chemical_class"),
                "positive_biomarker": value(row, "positive_biomarker"),
                "NHANES_variable": value(row, "NHANES_variable"),
                "mapping_type": value(row, "mapping_type"),
                "mapping_confidence": value(row, "mapping_confidence"),
                "mapping_source": value(row, "mapping_source"),
                "mapping_grade": value(row, "mapping_grade"),
                "mapping_gate_status": status,
                "mapping_gate_reason": value(row, "mapping_gate_reason"),
                "exclusion_scope": scope,
                "recommended_next_action": next_action,
            }
        )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    repo = args.repo.resolve()

    framework = repo / "analysis" / "disease_agnostic_environmental_framework"
    output_dir = framework / "t2d_exposure_opportunity"
    master_dir = output_dir / "01_candidate_master"
    audit_dir = output_dir / "02_mapping_audit"

    step5_path = framework / "step05_t2d_screen" / "t2d_primary_29_tests.csv"
    mapping_path = framework / "step02_biomarker_mapping" / "chemical_biomarker_mapping.csv"
    robustness_path = framework / "step06_t2d_robustness" / "t2d_robustness_results.csv"

    step5 = read_csv(step5_path)
    mapping = read_csv(mapping_path)
    robustness = read_csv(robustness_path)
    positives = positive_tests(step5)
    positives_by_biomarker = {value(row, "biomarker"): row for row in positives}
    robust_by_biomarker = {value(row, "biomarker"): row for row in robustness}

    missing_robustness = sorted(set(positives_by_biomarker) - set(robust_by_biomarker))
    if missing_robustness:
        raise RuntimeError(f"Missing Step 6 robustness rows for positive biomarkers: {missing_robustness}")

    mapping_rows = make_mapping_rows(mapping, positives_by_biomarker, robust_by_biomarker)
    if not mapping_rows:
        raise RuntimeError("No mapped chemical rows were found for the 14 positive tests")
    mapped_biomarkers = {value(row, "positive_biomarker") for row in mapping_rows}
    missing_mappings = sorted(set(positives_by_biomarker) - mapped_biomarkers)
    if missing_mappings:
        raise RuntimeError(f"Positive biomarkers without mapped chemical rows: {missing_mappings}")

    unique_rows = make_unique_rows(mapping_rows)
    exclusion_rows = make_exclusion_rows(mapping_rows)

    write_csv(master_dir / "all_upstream_chemicals.csv", OUTPUT_MASTER_FIELDS, mapping_rows)
    write_csv(master_dir / "unique_candidate_chemicals.csv", OUTPUT_UNIQUE_FIELDS, unique_rows)
    write_csv(audit_dir / "mapping_specificity_audit.csv", OUTPUT_MASTER_FIELDS, mapping_rows)
    write_csv(audit_dir / "mapping_exclusion_log.csv", OUTPUT_EXCLUSION_FIELDS, exclusion_rows)

    grade_counts = Counter(value(row, "mapping_grade") for row in mapping_rows)
    status_counts = Counter(value(row, "mapping_gate_status") for row in mapping_rows)
    disposition_counts = Counter(value(row, "mapping_gate_disposition") for row in unique_rows)
    per_test_counts = {
        biomarker: sum(value(row, "positive_biomarker") == biomarker for row in mapping_rows)
        for biomarker in sorted(positives_by_biomarker)
    }

    generated_files = [
        master_dir / "all_upstream_chemicals.csv",
        master_dir / "unique_candidate_chemicals.csv",
        audit_dir / "mapping_specificity_audit.csv",
        audit_dir / "mapping_exclusion_log.csv",
    ]
    report = f"""# Phase 1 — T2D exposure-opportunity candidate and mapping gate

## Material Passport

- **Material ID:** `T2D-EXPOSURE-OPPORTUNITY-PHASE1`
- **Status:** `complete_mapping_gate_only`
- **Script:** `{SCRIPT_VERSION}`
- **Generated (UTC):** `{datetime.now(timezone.utc).isoformat()}`
- **Scope:** all upstream mapped chemicals associated with the 14 frozen T2D FDR-positive tests.
- **Not performed:** literature collision audit, mechanism analysis, target nomination, docking, experimental feasibility, and opportunity scoring.

## Frozen input and provenance

- T2D input: `step05_t2d_screen/t2d_primary_29_tests.csv`; positive tests were derived from `FDR_supported=True` and verified to use the frozen 29-test denominator.
- Chemical mapping input: `step02_biomarker_mapping/chemical_biomarker_mapping.csv`; only rows with `mapped=True` and one of the 14 positive biomarkers were included in the candidate master.
- Robustness annotation: `step06_t2d_robustness/t2d_robustness_results.csv`; used only as a downstream annotation, not as a mapping gate.
- No GeneCards, CTD chemical–gene interactions, disease-specific pathway data, or literature counts were used in this phase.

## Counts

- Positive T2D tests: **{len(positives)}**
- Upstream mapped chemical–biomarker rows: **{len(mapping_rows)}**
- Unique upstream chemical IDs: **{len(unique_rows)}**
- Mapping rows by grade: **A={grade_counts.get('A', 0)}, B={grade_counts.get('B', 0)}, C={grade_counts.get('C', 0)}**
- Mapping rows by gate status: **{json.dumps(dict(sorted(status_counts.items())), ensure_ascii=False)}**
- Unique-chemical dispositions: **{json.dumps(dict(sorted(disposition_counts.items())), ensure_ascii=False)}**

## Per-test upstream mapping coverage

| Positive test | Upstream mapped rows |
|---|---:|
""" + "\n".join(f"| `{biomarker}` | {count} |" for biomarker, count in per_test_counts.items()) + f"""

## Gate interpretation

- **Grade A:** direct analyte identity supported by the neutral mapping record (parent elemental analyte or direct serum PFAS analyte).
- **Grade B:** specific urinary metabolite or parent–validated-metabolite relationship; direct urinary metabolites retain a parent-inference caution.
- **Grade C:** family/proxy mapping or an elemental/species name not resolved by the assay. These rows remain in the full master for audit but are not promoted as compound-specific primary candidates.
- A chemical is marked `advance_to_literature_audit` only when it has a direct Grade A mapping with `pass`, or a specific Grade B mapping with the corresponding review flag. C-only mappings are `proxy_only_not_primary`; internally inconsistent elemental rows are excluded pending manual mapping review.

## File outputs

- `01_candidate_master/all_upstream_chemicals.csv`: one row per mapped chemical–positive biomarker relationship.
- `01_candidate_master/unique_candidate_chemicals.csv`: one row per unique chemical ID.
- `02_mapping_audit/mapping_specificity_audit.csv`: row-level grade and gate rationale.
- `02_mapping_audit/mapping_exclusion_log.csv`: proxy/conditional/mismatch rows not eligible for an unqualified compound-specific shortlist.

This is a mapping gate, not a novelty or mechanistic conclusion. A Grade A/B mapping means the exposure-to-analyte link is more actionable; it does not establish a T2D causal association.
"""
    report_path = output_dir / "PHASE1_T2D_EXPOSURE_OPPORTUNITY_REPORT.md"
    report_path.write_text(report, encoding="utf-8")
    generated_files.append(report_path)

    manifest = {
        "material_id": "T2D-EXPOSURE-OPPORTUNITY-PHASE1",
        "script_version": SCRIPT_VERSION,
        "status": "complete_mapping_gate_only",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "positive_test_count": len(positives),
        "positive_biomarkers": [value(row, "biomarker") for row in positives],
        "upstream_mapped_mapping_row_count": len(mapping_rows),
        "unique_chemical_count": len(unique_rows),
        "mapping_grade_counts": dict(sorted(grade_counts.items())),
        "mapping_gate_status_counts": dict(sorted(status_counts.items())),
        "unique_disposition_counts": dict(sorted(disposition_counts.items())),
        "per_test_mapping_row_counts": per_test_counts,
        "inputs": {
            "t2d_primary_29_tests.csv": {"path": str(step5_path.relative_to(repo)), "sha256": sha256_file(step5_path)},
            "chemical_biomarker_mapping.csv": {"path": str(mapping_path.relative_to(repo)), "sha256": sha256_file(mapping_path)},
            "t2d_robustness_results.csv": {"path": str(robustness_path.relative_to(repo)), "sha256": sha256_file(robustness_path)},
        },
        "outputs": {
            str(path.relative_to(repo)): {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in generated_files
        },
        "firewall_note": "Outcome statistics annotate candidates after the frozen T2D screen; no disease information was used to create or map the environmental candidate universe in this phase.",
        "not_performed": ["literature_collision", "mechanism", "target_nomination", "docking", "experimental_feasibility", "opportunity_scoring"],
    }
    (output_dir / "PHASE1_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": "complete_mapping_gate_only",
        "positive_tests": len(positives),
        "mapping_rows": len(mapping_rows),
        "unique_chemicals": len(unique_rows),
        "grade_counts": dict(sorted(grade_counts.items())),
        "gate_status_counts": dict(sorted(status_counts.items())),
        "unique_dispositions": dict(sorted(disposition_counts.items())),
        "output_dir": str(output_dir),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
