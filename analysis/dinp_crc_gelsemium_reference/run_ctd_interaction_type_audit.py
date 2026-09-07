#!/usr/bin/env python3
"""Audit the semantic type of CTD evidence for the fresh 41-gene set.

This is a provenance audit, not a new target-discovery analysis.  It reports
whether the CTD records for each DINP--CRC intersection gene describe direct
binding, expression/activity changes, reaction/abundance context, or only
multi-chemical co-treatment.  It keeps the source's original action labels
and flags instead of converting every chemical--gene record into a direct
target claim.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
INTERSECTION = ROOT / "outputs" / "dinp_crc_intersection_41.csv"
CTD = ROOT.parent / "dinp_crc_multi_database_target_convergence" / "outputs" / "source_records" / "ctd_dinp_human_interactions.csv"
OUT_DIR = ROOT / "outputs" / "ctd_interaction_audit_41_genes"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def is_true(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def classify(record_rows: list[dict[str, str]]) -> dict[str, Any]:
    action_labels = sorted({label.strip() for row in record_rows for label in row.get("InteractionActions", "").split("|") if label.strip()})
    pmids = sorted({p.strip() for row in record_rows for p in row.get("PubMedIDs", "").split(";") if p.strip()})
    direct_binding_rows = [row for row in record_rows if "affects^binding" in row.get("InteractionActions", "")]
    expression_rows = [row for row in record_rows if "expression" in row.get("InteractionActions", "").lower()]
    activity_rows = [row for row in record_rows if "activity" in row.get("InteractionActions", "").lower()]
    reaction_rows = [row for row in record_rows if "reaction" in row.get("InteractionActions", "").lower()]
    abundance_rows = [row for row in record_rows if "abundance" in row.get("InteractionActions", "").lower()]
    single_rows = [row for row in record_rows if is_true(row.get("single_chemical_record_flag"))]
    cotreatment_rows = [row for row in record_rows if is_true(row.get("multi_chemical_co_treatment_flag"))]
    single_direct_binding = [row for row in direct_binding_rows if is_true(row.get("single_chemical_record_flag"))]
    single_functional = [row for row in expression_rows + activity_rows if is_true(row.get("single_chemical_record_flag"))]

    if direct_binding_rows:
        interpretation = "CTD contains DINP-specific binding/interaction record(s); activity is present where recorded, but this remains source-level evidence rather than in-vivo proof."
        evidence_class = "single_chemical_binding_or_interaction"
    elif single_functional:
        interpretation = "CTD contains at least one single-chemical functional response record (expression/activity); not necessarily direct binding."
        evidence_class = "single_chemical_functional_response"
    elif record_rows and len(cotreatment_rows) == len(record_rows):
        interpretation = "CTD evidence is co-treatment-only; it should not be described as DINP-specific gene regulation."
        evidence_class = "co_treatment_only"
    elif record_rows:
        interpretation = "CTD contains chemical-gene literature records, but the available records do not establish a direct DINP-specific binding or single-chemical functional response."
        evidence_class = "literature_association_without_direct_single_chemical_support"
    else:
        interpretation = "No CTD record for this intersection gene; its DINP support comes from another exposure source."
        evidence_class = "no_CTD_record"

    return {
        "CTD_record_count": len(record_rows),
        "CTD_unique_PubMed_count": len(pmids),
        "CTD_unique_PubMedIDs": ";".join(pmids),
        "CTD_single_chemical_record_count": len(single_rows),
        "CTD_cotreatment_record_count": len(cotreatment_rows),
        "CTD_single_chemical_support": int(bool(single_rows)),
        "CTD_cotreatment_only": int(bool(record_rows) and len(cotreatment_rows) == len(record_rows)),
        "CTD_direct_binding_record_count": len(direct_binding_rows),
        "CTD_single_chemical_binding_record_count": len(single_direct_binding),
        "CTD_expression_record_count": len(expression_rows),
        "CTD_activity_record_count": len(activity_rows),
        "CTD_reaction_record_count": len(reaction_rows),
        "CTD_abundance_record_count": len(abundance_rows),
        "CTD_action_labels": ";".join(action_labels),
        "CTD_evidence_class": evidence_class,
        "CTD_evidence_interpretation": interpretation,
    }


def main() -> None:
    intersection = read_csv(INTERSECTION)
    ctd_rows = read_csv(CTD)
    genes = {row["gene_symbol"].strip().upper() for row in intersection}
    if len(genes) != 41:
        raise RuntimeError(f"expected 41 intersection genes, found {len(genes)}")

    by_gene: dict[str, list[dict[str, str]]] = {gene: [] for gene in genes}
    for row in ctd_rows:
        gene = (row.get("normalized_gene_symbol") or row.get("GeneSymbol") or "").strip().upper()
        if gene in by_gene:
            by_gene[gene].append(row)

    rows: list[dict[str, Any]] = []
    for base in intersection:
        gene = base["gene_symbol"].strip().upper()
        details = classify(by_gene[gene])
        rows.append(
            {
                "gene_symbol": gene,
                "crc_article_rank": base.get("crc_article_rank", ""),
                "crc_relevance_score": base.get("crc_relevance_score", ""),
                "DINP_source_flags": ";".join(
                    source for source in ("CTD", "ToxCast", "Tox21", "T3DB") if is_true(base.get(source))
                ),
                **details,
            }
        )
    rows.sort(key=lambda row: (row["CTD_evidence_class"], row["gene_symbol"]))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    table_path = OUT_DIR / "ctd_interaction_type_audit.csv"
    write_csv(table_path, rows, fields)

    direct_genes = [row for row in rows if row["CTD_direct_binding_record_count"] > 0]
    single_functional_genes = [row for row in rows if row["CTD_evidence_class"] == "single_chemical_functional_response"]
    cotreatment_only_genes = [row for row in rows if row["CTD_evidence_class"] == "co_treatment_only"]
    no_ctd_genes = [row for row in rows if row["CTD_evidence_class"] == "no_CTD_record"]
    summary_lines = [
        "# CTD interaction-type audit for the fresh 41-gene DINP–CRC intersection",
        "",
        f"Generated (UTC): {now_utc()}",
        "",
        "## Main result",
        "",
        f"The 41-gene intersection contains **{len(direct_genes)}** genes with at least one CTD DINP-specific binding/interaction action, **{len(single_functional_genes)}** additional genes with single-chemical functional-response evidence, **{len(cotreatment_only_genes)}** co-treatment-only genes, and **{len(no_ctd_genes)}** genes without a CTD record.",
        "",
        "A CTD chemical–gene record is not automatically a direct target claim. The original action labels, single-chemical flags, co-treatment flags, and PubMed IDs are retained below so that direct interaction, functional response, and literature association remain separate.",
        "",
        "| Evidence class | Genes |",
        "|---|---:|",
        f"| single-chemical binding or interaction | {len(direct_genes)} |",
        f"| single-chemical functional response without binding action | {len(single_functional_genes)} |",
        f"| co-treatment only | {len(cotreatment_only_genes)} |",
        f"| no CTD record | {len(no_ctd_genes)} |",
        "",
        "## Binding/interaction genes",
        "",
        "| Gene | DINP source flags | CTD binding records | Single-chemical binding records | Actions |",
        "|---|---|---:|---:|---|",
    ]
    for row in sorted(direct_genes, key=lambda row: row["gene_symbol"]):
        summary_lines.append(
            f"| **{row['gene_symbol']}** | {row['DINP_source_flags']} | {row['CTD_direct_binding_record_count']} | {row['CTD_single_chemical_binding_record_count']} | {row['CTD_action_labels']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "Genes with CTD binding/interaction actions are exposure-side mechanistic candidates. They are not all equivalent: direct binding/interaction records are distinct from expression/activity response records and from co-treatment literature. GO/Reactome/KEGG recurrence is a separate pathway-context layer and is not used here to upgrade CTD evidence class.",
            "",
            "## Files",
            "",
            "- `ctd_interaction_type_audit.csv`: one row per fresh intersection gene with CTD evidence decomposition.",
            "- `ctd_interaction_audit_summary.md`: reviewer-facing summary and binding subset.",
        ]
    )
    summary_path = OUT_DIR / "ctd_interaction_audit_summary.md"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest = {
        "generated_utc": now_utc(),
        "script": str(Path(__file__).resolve()),
        "inputs": {
            "intersection": {"path": str(INTERSECTION.resolve()), "sha256": sha256_file(INTERSECTION)},
            "ctd_source": {"path": str(CTD.resolve()), "sha256": sha256_file(CTD)},
        },
        "intersection_gene_count": len(rows),
        "ctd_source_record_count": len(ctd_rows),
        "direct_binding_gene_count": len(direct_genes),
        "single_functional_response_gene_count": len(single_functional_genes),
        "co_treatment_only_gene_count": len(cotreatment_only_genes),
        "no_ctd_gene_count": len(no_ctd_genes),
        "new_statistical_test_performed": False,
        "outputs": {
            "table": {"path": str(table_path.resolve()), "sha256": sha256_file(table_path)},
            "summary": {"path": str(summary_path.resolve()), "sha256": sha256_file(summary_path)},
        },
    }
    manifest_path = OUT_DIR / "ctd_interaction_audit_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print("CTD INTERACTION-TYPE AUDIT: PASS")
    print(f"Intersection genes: {len(rows)}")
    print(f"Direct binding/interaction genes: {len(direct_genes)}")
    print(f"Single-chemical functional-response genes: {len(single_functional_genes)}")
    print(f"Co-treatment-only genes: {len(cotreatment_only_genes)}")
    print(f"No-CTD genes: {len(no_ctd_genes)}")
    print(f"Table: {table_path.resolve()}")
    print(f"Summary: {summary_path.resolve()}")


if __name__ == "__main__":
    main()
