"""Audit the multiplicity unit of the outcome-blinded NHANES test set.

This script is deliberately outcome-free.  It reads only the frozen Step 4
test set, the Step 2 outcome-blinded chemical-to-biomarker ledger, and the
available actionability matrix as a supplemental parent-annotation source.
It does not read CRC results, GeneCards, CTD chemical-gene interactions, or
any robustness output.

The primary question is descriptive: does each of the 29 frozen tests
represent one measurable biomarker hypothesis, a parent-chemical hypothesis,
or a multi-chemical proxy?  No multiplicity denominator is changed here.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def find_repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / ".git").exists():
            return parent
    raise FileNotFoundError("Could not locate repository root")


ROOT = find_repo_root()
OUTDIR = Path(__file__).resolve().parent
TESTSET = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
STEP2 = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step02_biomarker_mapping" / "chemical_biomarker_mapping.csv"
LOCK = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step04_testset_freeze" / "PRE_DISEASE_TESTSET_LOCK.json"
SUPPLEMENTAL_PARENT = ROOT / "outputs" / "environmental_crc_267_identity_table_v2.csv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def split_ids(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [x.strip() for x in re.split(r"[;,|\s]+", str(value)) if x.strip()]


def clean_set(values: object) -> list[str]:
    if isinstance(values, pd.Series):
        raw = values.tolist()
    elif isinstance(values, (list, tuple, set)):
        raw = values
    else:
        raw = [values]
    return sorted({str(x).strip() for x in raw if str(x).strip() and str(x).strip().lower() != "nan"})


def unique_join(values: object) -> str:
    return "|".join(clean_set(values))


def bool_string(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def primary_class(variable: str, mapped: pd.DataFrame) -> str:
    if variable.startswith("LBXPF"):
        return "PFAS"
    if variable == "URXBPH":
        return "bisphenols"
    if variable.startswith("URXP"):
        return "PAHs"
    if variable.startswith("URXU"):
        return "metals"
    if variable.startswith("URX"):
        return "phthalates"
    classes = unique_join(mapped.get("chemical_class", []))
    return classes or "unresolved_class"


def mapping_kind(variable: str, mapped: pd.DataFrame, known_parents: list[str]) -> tuple[str, str, str]:
    """Return unit_kind, candidate unit, and a transparent basis."""
    relationships = "|".join(clean_set(mapped.get("mapping_type", [])))
    n_map = len(mapped)
    if known_parents:
        parent_key = "+".join(known_parents)
        return (
            "known_parent_linked",
            f"parent::{parent_key}",
            "parent_compound/metabolite_of recorded in the supplemental identity table",
        )
    if variable == "URXBPH":
        return (
            "family_proxy",
            f"proxy::{variable}",
            f"one bisphenol assay mapped to {n_map} CTD entities; no single parent recorded",
        )
    if variable.startswith("URXP"):
        return (
            "family_proxy",
            f"proxy::{variable}",
            f"urinary OH-PAH proxy mapped to {n_map} CTD entities; parent specificity limited",
        )
    if variable.startswith("URXU"):
        return (
            "elemental_proxy",
            f"elemental_proxy::{variable}",
            f"elemental urinary biomarker mapped to {n_map} metal/species entities",
        )
    if "direct serum PFAS analyte" in relationships:
        return (
            "direct_analyte",
            f"direct::{variable}",
            "direct serum PFAS analyte; the measured chemical is the hypothesis unit",
        )
    if "direct urinary analyte" in relationships and n_map > 1:
        return (
            "family_proxy",
            f"proxy::{variable}",
            f"one urinary analyte mapped to {n_map} CTD entities; no single parent recorded",
        )
    if "direct urinary metabolite" in relationships:
        return (
            "parent_unresolved",
            f"biomarker::{variable}",
            "direct urinary metabolite, but parent_compound was not recorded in the Step 1–4 ledger",
        )
    return (
        "biomarker_level",
        f"biomarker::{variable}",
        "one frozen NHANES test retained as the measurable unit",
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tests = pd.read_csv(TESTSET, dtype=str, keep_default_na=False)
    step2 = pd.read_csv(STEP2, dtype=str, keep_default_na=False)
    supplemental = pd.read_csv(
        SUPPLEMENTAL_PARENT,
        dtype=str,
        keep_default_na=False,
        usecols=["ChemicalID", "parent_compound", "metabolite_of"],
    )

    test_variables = set(tests["variable"])
    mapped = step2.loc[
        step2["mapped"].map(bool_string) & step2["human_biomarker"].isin(test_variables)
    ].copy()
    mapped["chemical_id"] = mapped["chemical_id"].str.strip()
    mapped["human_biomarker"] = mapped["human_biomarker"].str.strip()
    supplemental["ChemicalID"] = supplemental["ChemicalID"].str.strip()
    supplemental_by_id = {key: group for key, group in supplemental.groupby("ChemicalID", sort=False)}

    audit_rows: list[dict[str, object]] = []
    expanded_rows: list[dict[str, object]] = []
    test_to_unit: dict[str, str] = {}
    unit_members: defaultdict[str, list[str]] = defaultdict(list)

    for row in tests.to_dict("records"):
        variable = row["variable"]
        chemical_ids = split_ids(row["chemical_ids"])
        sub = mapped.loc[
            (mapped["human_biomarker"] == variable) & mapped["chemical_id"].isin(chemical_ids)
        ].copy()
        if len(sub) != int(row["mapping_count"]):
            raise AssertionError(
                f"{variable}: Step 4 mapping_count={row['mapping_count']} but Step 2 matched {len(sub)} rows"
            )

        parent_rows = []
        for chemical_id in chemical_ids:
            if chemical_id in supplemental_by_id:
                parent_rows.append(supplemental_by_id[chemical_id])
        parent = pd.concat(parent_rows, ignore_index=True) if parent_rows else pd.DataFrame()
        known_parents: set[str] = set()
        if not parent.empty:
            for col in ["parent_compound", "metabolite_of"]:
                if col in parent:
                    known_parents.update(x for x in parent[col].astype(str).str.strip() if x)

        unit_kind, unit_id, unit_basis = mapping_kind(variable, sub, sorted(known_parents))
        test_to_unit[variable] = unit_id
        unit_members[unit_id].append(variable)

        classes = clean_set(sub["chemical_class"])
        relationships = clean_set(sub["mapping_type"])
        confidence = ";".join(
            f"{key}:{value}" for key, value in sub["mapping_confidence"].value_counts().sort_index().items()
        )
        parent_annotation_n = int(
            parent.get("parent_compound", pd.Series(dtype=str)).astype(str).str.strip().ne("").sum()
            if not parent.empty else 0
        )
        metabolite_annotation_n = int(
            parent.get("metabolite_of", pd.Series(dtype=str)).astype(str).str.strip().ne("").sum()
            if not parent.empty else 0
        )
        proxy_reason = (
            "single-chemical direct analyte"
            if unit_kind == "direct_analyte" and len(chemical_ids) == 1
            else unit_basis
        )
        audit_rows.append(
            {
                "test_id": row["test_id"],
                "biomarker": row["biomarker"],
                "variable": variable,
                "matrix": row["matrix"],
                "cycles": row["cycles"],
                "weight": row["weight"],
                "chemical_mapping_count": len(sub),
                "unique_chemical_id_count": len(set(chemical_ids)),
                "chemical_ids": ";".join(chemical_ids),
                "chemical_names": ";".join(sorted(set(sub["chemical_name"]))),
                "chemical_class": ";".join(classes),
                "primary_class": primary_class(variable, sub),
                "exposure_axis": row["exposure_axes"],
                "mapping_type": ";".join(relationships),
                "mapping_confidence_counts": confidence,
                "recorded_parent_compound": ";".join(sorted(known_parents)),
                "parent_annotation_rows": parent_annotation_n,
                "metabolite_of_annotation_rows": metabolite_annotation_n,
                "parent_annotation_coverage_pct": round(100 * parent_annotation_n / len(sub), 2) if len(sub) else 0.0,
                "unit_kind": unit_kind,
                "candidate_hypothesis_unit": unit_id,
                "candidate_unit_basis": proxy_reason,
                "current_step4_primary_unit": "biomarker_test",
                "current_step4_fdr_family_member": True,
            }
        )
        for chemical_id in chemical_ids:
            detail = sub.loc[sub["chemical_id"] == chemical_id]
            name = detail["chemical_name"].iloc[0] if not detail.empty else ""
            expanded_rows.append(
                {
                    "variable": variable,
                    "test_id": row["test_id"],
                    "chemical_id": chemical_id,
                    "chemical_name": name,
                    "primary_class": primary_class(variable, sub),
                    "step2_mapping_type": ";".join(clean_set(detail["mapping_type"])) if not detail.empty else "",
                    "step2_mapping_confidence": ";".join(clean_set(detail["mapping_confidence"])) if not detail.empty else "",
                    "recorded_parent_compound": ";".join(sorted(known_parents)),
                    "candidate_hypothesis_unit": unit_id,
                }
            )

    audit = pd.DataFrame(audit_rows).sort_values("variable").reset_index(drop=True)
    expanded = pd.DataFrame(expanded_rows).sort_values(["variable", "chemical_id"]).reset_index(drop=True)
    audit.to_csv(OUTDIR / "step4_test_hypothesis_mapping.csv", index=False)
    expanded.to_csv(OUTDIR / "step4_test_chemical_membership.csv", index=False)

    unit_rows: list[dict[str, object]] = []
    for unit_id, members in sorted(unit_members.items()):
        member = audit.loc[audit["variable"].isin(members)]
        ids = set()
        for value in member["chemical_ids"]:
            ids.update(split_ids(value))
        known_parent_values = set(";".join(member["recorded_parent_compound"]).split(";")) - {""}
        known_parent = ";".join(sorted(known_parent_values))
        unit_rows.append(
            {
                "unit_id": unit_id,
                "unit_kind": ";".join(sorted(set(member["unit_kind"]))),
                "unit_label": unit_id.split("::", 1)[-1],
                "member_test_count": len(members),
                "member_biomarkers": ";".join(sorted(members)),
                "parent_compound_if_recorded": known_parent,
                "chemical_mapping_count_sum": int(member["chemical_mapping_count"].sum()),
                "unique_chemical_id_count": len(ids),
                "primary_classes": ";".join(sorted(set(member["primary_class"]))),
                "basis": "; ".join(sorted(set(member["candidate_unit_basis"]))),
                "is_current_primary_fdr_unit": False,
                "status": "secondary_operational_grouping_only",
            }
        )
    units = pd.DataFrame(unit_rows)
    units.to_csv(OUTDIR / "hypothesis_unit_summary.csv", index=False)

    axis_memberships = sum(
        len([label.strip() for label in value.split(";") if label.strip()])
        for value in tests["exposure_axes"]
    )
    axis_labels = sorted({label.strip() for value in tests["exposure_axes"] for label in value.split(";") if label.strip()})
    duplicate_ids = expanded.loc[expanded["chemical_id"].duplicated(False), ["chemical_id", "variable"]]
    duplicate_detail = "; ".join(
        f"{chemical_id}: {', '.join(group['variable'].tolist())}"
        for chemical_id, group in duplicate_ids.groupby("chemical_id", sort=True)
    )
    known_parent_tests = audit.loc[audit["recorded_parent_compound"].ne(""), "variable"].tolist()
    class_counts = audit.groupby("primary_class")["variable"].nunique().sort_index().to_dict()
    kind_counts = audit.groupby("unit_kind")["variable"].nunique().sort_index().to_dict()

    lock_data = json.loads(LOCK.read_text(encoding="utf-8"))
    try:
        git_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        git_commit = "unavailable"
    generated = datetime.now(timezone.utc).isoformat()
    lock = {
        "lock_type": "HYPOTHESIS_UNIT_AUDIT_LOCK",
        "generated_utc": generated,
        "git_commit_at_generation": git_commit,
        "outcome_free": True,
        "sources": {
            "step4_testset": str(TESTSET),
            "step4_testset_sha256": sha256(TESTSET),
            "step4_lock": str(LOCK),
            "step4_lock_sha256": sha256(LOCK),
            "step2_mapping": str(STEP2),
            "step2_mapping_sha256": sha256(STEP2),
            "supplemental_parent_annotation": str(SUPPLEMENTAL_PARENT),
            "supplemental_parent_annotation_sha256": sha256(SUPPLEMENTAL_PARENT),
        },
        "current_frozen_primary": {
            "unit_definition": "unique NHANES biomarker test",
            "test_count": int(len(tests)),
            "planned_fdr_denominator": int(lock_data.get("planned_fdr_denominator", len(tests))),
            "source_statement": "Step 4 collapse was performed only at the unique NHANES test level.",
        },
        "mapping_counts": {
            "actionable_chemical_biomarker_mappings": int(len(mapped)),
            "unique_chemical_ids_across_tests": int(expanded["chemical_id"].nunique()),
            "mapping_rows_across_tests": int(len(expanded)),
            "duplicate_chemical_memberships": duplicate_detail,
        },
        "descriptive_hierarchy": {
            "unique_exposure_axis_labels": int(len(axis_labels)),
            "exposure_axis_labels": axis_labels,
            "axis_memberships": int(axis_memberships),
            "known_parent_annotated_test_count": int(len(known_parent_tests)),
            "known_parent_annotated_tests": known_parent_tests,
            "known_parent_collapsed_operational_unit_count": int(len(units)),
            "class_test_counts": class_counts,
            "unit_kind_test_counts": kind_counts,
        },
        "interpretation_guardrail": "No denominator change is authorized by this audit. Any parent-level or hierarchical multiplicity analysis must be labeled secondary/post hoc unless frozen before CRC outcome access.",
    }
    (OUTDIR / "HYPOTHESIS_UNIT_AUDIT_LOCK.json").write_text(
        json.dumps(lock, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    report = [
        "# Step 4 hypothesis-unit audit",
        "",
        "## Scope",
        "",
        "This audit reads only the outcome-blinded Step 1–4 mapping artifacts. It does not read CRC results, GeneCards, CTD chemical–gene interactions, or any downstream model output. It changes no FDR denominator.",
        "",
        "## Direct findings",
        "",
        f"- Frozen Step 4 biomarker tests: **{len(tests)}**.",
        f"- Actionable chemical–biomarker mapping rows represented by those tests: **{len(mapped)}**.",
        f"- Unique CTD chemical IDs represented across the tests: **{expanded['chemical_id'].nunique()}**.",
        f"- Distinct exposure-axis labels: **{len(axis_labels)}**, with **{axis_memberships}** test-to-axis memberships because some tests carry more than one axis label.",
        f"- Supplemental parent annotations are present for **{len(known_parent_tests)}/{len(tests)}** tests; missing parent metadata is not inferred from metabolite names.",
        f"- Known-parent collapse produces **{len(units)}** operational groups, but this is not a replacement for the frozen 29-test primary family.",
        "",
        "## Current primary unit",
        "",
        "The canonical Step 4 lock explicitly states that collapse was performed only at the unique NHANES test level and sets the planned downstream FDR denominator to 29. Therefore the defensible primary multiplicity family remains **29 measured biomarker tests**. The 29 tests are not 29 independent chemicals: several tests represent multi-chemical proxies, and several phthalate tests can belong to one recorded parent axis.",
        "",
        "## Descriptive grouping",
        "",
        "The operational secondary grouping uses a parent compound only when it is explicitly recorded in the supplemental identity table. Otherwise it retains a biomarker/proxy-level unit rather than guessing a parent.",
        "",
        "| Grouping | Count | Interpretation |",
        "| --- | ---: | --- |",
        f"| Frozen biomarker-test family | {len(tests)} | Current Step 4 primary unit and flat-FDR family |",
        f"| Known-parent operational groups | {len(units)} | Secondary descriptive grouping; incomplete parent annotation |",
        f"| Exposure-axis labels | {len(axis_labels)} | Labels, not a valid denominator because memberships overlap and vary in specificity |",
        "",
        "## Complete 29-test map",
        "",
        "| Biomarker test | Chemical class | Mappings | Recorded parent | Operational unit | Parent annotation |",
        "| --- | --- | ---: | --- | --- | ---: |",
    ]
    for row in audit.itertuples(index=False):
        report.append(
            f"| {row.variable} | {row.primary_class} | {row.chemical_mapping_count} | "
            f"{row.recorded_parent_compound or 'not recorded'} | {row.candidate_hypothesis_unit} | "
            f"{row.parent_annotation_coverage_pct:.1f}% |"
        )
    report += [
        "",
        "## Known parent links",
        "",
        "- `DINP`: URXCOP.",
        "- `DEHP`: URXECP, URXMHH, URXMHP, URXMOH.",
        "- `BBzP`: URXMZP.",
        "- URXMBP, URXMEP, and URXMIB remain parent-unresolved in the current Step 1–4 ledger; no parent was imputed.",
        "- Parent annotation is partial for URXECP and URXMOH (50% of their mapped rows in the supplemental identity table); this is why the 26-unit grouping is explicitly operational rather than a complete parent ontology.",
        "",
        "## Statistical interpretation guardrail",
        "",
        "This audit supports a transparent distinction between biomarker-level and parent-level hypotheses, but it does not rescue a CRC result by changing the denominator after outcome inspection. If a parent-level or hierarchical FDR analysis is later added, it must be reported as a secondary/post hoc reanalysis unless the hierarchy is frozen before CRC outcome access. The flat 29-test result remains the primary reference for the current manuscript state.",
        "",
        "## Files",
        "",
        "- `step4_test_hypothesis_mapping.csv`: one row per frozen biomarker test.",
        "- `step4_test_chemical_membership.csv`: one row per chemical membership in those tests.",
        "- `hypothesis_unit_summary.csv`: operational parent/proxy grouping.",
        "- `HYPOTHESIS_UNIT_AUDIT_LOCK.json`: source hashes and audit lock.",
    ]
    (OUTDIR / "HYPOTHESIS_UNIT_AUDIT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps(lock["descriptive_hierarchy"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
