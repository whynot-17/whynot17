#!/usr/bin/env python3
"""Audit DINP evidence and pathway-prioritization components for 41 genes.

The audit separates three evidence layers that should not be collapsed into a
single hidden score:

1. CRC reference relevance from the published 1,893-gene Table S2 ranking;
2. DINP-side source support from CTD, ToxCast/Tox21, and T3DB;
3. GO/Reactome/KEGG recurrence used by the descriptive target ordering.

It is an audit/decomposition layer. It does not re-rank the original
intersection, recalculate enrichment, or correct for annotation density.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"
INTERSECTION = OUTPUTS / "dinp_crc_intersection_41.csv"
EXPOSURE = ROOT.parent / "dinp_crc_multi_database_target_convergence" / "outputs" / "dinp_exposure_gene_matrix_ctd_toxcast_tox21_t3db.csv"
TARGET = OUTPUTS / "target_prioritization_1893" / "go_reactome_kegg_target_prioritization.csv"
OUT_DIR = OUTPUTS / "evidence_audit_41_genes"


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


def integer(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def number(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_flag(value: Any) -> int:
    return 1 if str(value).strip().lower() in {"1", "true", "yes"} else 0


def main() -> None:
    for path in (INTERSECTION, EXPOSURE, TARGET):
        if not path.exists():
            raise FileNotFoundError(path)

    intersection = read_csv(INTERSECTION)
    exposure_rows = {row["gene_symbol"].strip().upper(): row for row in read_csv(EXPOSURE)}
    target_rows = {row["gene_symbol"].strip().upper(): row for row in read_csv(TARGET)}
    if len(intersection) != 41 or len({row["gene_symbol"].strip().upper() for row in intersection}) != 41:
        raise RuntimeError("expected exactly 41 unique intersection genes")
    if len(exposure_rows) != 93:
        raise RuntimeError(f"expected 93 DINP exposure genes, found {len(exposure_rows)}")
    if len(target_rows) != 41:
        raise RuntimeError(f"expected 41 target-prioritization rows, found {len(target_rows)}")

    rows: list[dict[str, Any]] = []
    for source_row in intersection:
        gene = source_row["gene_symbol"].strip().upper()
        exp = exposure_rows.get(gene)
        target = target_rows.get(gene)
        if exp is None or target is None:
            raise RuntimeError(f"missing audit row for {gene}")

        ctd = as_flag(exp.get("CTD"))
        toxcast = as_flag(exp.get("ToxCast"))
        tox21 = as_flag(exp.get("Tox21"))
        t3db = as_flag(exp.get("T3DB"))
        comptox = int(bool(toxcast or tox21))
        source_family_count = ctd + comptox + t3db
        source_count = integer(exp.get("exposure_support_count"))
        if source_count != ctd + toxcast + tox21 + t3db:
            raise RuntimeError(f"source count mismatch for {gene}")

        pathway_term_count = (
            integer(target.get("GO_representative_count"))
            + integer(target.get("Reactome_term_count"))
            + integer(target.get("KEGG_term_count"))
        )
        same_theme_count = integer(target.get("same_theme_two_source_count"))
        cross_resource_count = integer(target.get("source_count"))
        if pathway_term_count != integer(target.get("total_aligned_term_count")):
            raise RuntimeError(f"pathway term count mismatch for {gene}")

        rows.append(
            {
                "gene_symbol": gene,
                "crc_article_rank": integer(source_row.get("crc_article_rank")),
                "crc_relevance_score": number(source_row.get("crc_relevance_score")),
                "crc_rank_percentile_top": f"{(1894 - integer(source_row.get('crc_article_rank'))) / 1893:.6f}",
                "CTD_support": ctd,
                "ToxCast_support": toxcast,
                "Tox21_support": tox21,
                "CompTox_family_support": comptox,
                "T3DB_support": t3db,
                "DINP_raw_source_count": source_count,
                "DINP_source_family_count": source_family_count,
                "DINP_evidence_profile": "; ".join(
                    source
                    for source, present in (("CTD", ctd), ("CompTox", comptox), ("T3DB", t3db))
                    if present
                ),
                "DINP_multi_family_support": int(source_family_count >= 2),
                "GO_representative_count": integer(target.get("GO_representative_count")),
                "Reactome_term_count": integer(target.get("Reactome_term_count")),
                "KEGG_term_count": integer(target.get("KEGG_term_count")),
                "pathway_recurrence_total": pathway_term_count,
                "pathway_source_count": cross_resource_count,
                "same_theme_two_source_count": same_theme_count,
                "same_theme_three_source_count": integer(target.get("same_theme_three_source_count")),
                "same_theme_labels": target.get("same_theme_labels", ""),
                "pathway_theme_breadth": integer(target.get("theme_breadth")),
                "pathway_anchor_class": (
                    "strong_cross-resource_pathway_anchor"
                    if cross_resource_count == 3 and same_theme_count >= 2
                    else "cross-resource_pathway_supported"
                    if cross_resource_count >= 2 and same_theme_count >= 1
                    else "single-resource_or_nonconcordant_pathway_context"
                ),
                "original_biological_priority_rank": integer(target.get("biological_priority_rank")),
            }
        )

    rows.sort(key=lambda row: row["original_biological_priority_rank"])
    max_recurrence = max(row["pathway_recurrence_total"] for row in rows)
    max_breadth = max(row["pathway_theme_breadth"] for row in rows)
    for row in rows:
        row["pathway_recurrence_percentile_within_41"] = f"{sum(other['pathway_recurrence_total'] <= row['pathway_recurrence_total'] for other in rows) / len(rows):.6f}"
        row["pathway_breadth_percentile_within_41"] = f"{sum(other['pathway_theme_breadth'] <= row['pathway_theme_breadth'] for other in rows) / len(rows):.6f}"
        row["annotation_breadth_proxy_note"] = "pathway recurrence only; not a corrected annotation-density estimate"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    matrix_path = OUT_DIR / "dinp_crc_41_gene_evidence_matrix.csv"
    write_csv(matrix_path, rows, fields)

    decomposition_fields = [
        "original_biological_priority_rank",
        "gene_symbol",
        "crc_article_rank",
        "crc_relevance_score",
        "DINP_evidence_profile",
        "DINP_raw_source_count",
        "DINP_source_family_count",
        "DINP_multi_family_support",
        "GO_representative_count",
        "Reactome_term_count",
        "KEGG_term_count",
        "pathway_recurrence_total",
        "pathway_source_count",
        "same_theme_two_source_count",
        "same_theme_three_source_count",
        "same_theme_labels",
        "pathway_theme_breadth",
        "pathway_recurrence_percentile_within_41",
        "pathway_breadth_percentile_within_41",
        "pathway_anchor_class",
        "annotation_breadth_proxy_note",
    ]
    decomposition_path = OUT_DIR / "prioritization_decomposition.csv"
    write_csv(decomposition_path, rows, decomposition_fields)

    source_family_counts = {}
    for row in rows:
        key = row["DINP_source_family_count"]
        source_family_counts[key] = source_family_counts.get(key, 0) + 1
    strong_anchor_count = sum(row["pathway_anchor_class"] == "strong_cross-resource_pathway_anchor" for row in rows)
    multi_dinp_count = sum(row["DINP_multi_family_support"] for row in rows)
    rxra = next(row for row in rows if row["gene_symbol"] == "RXRA")

    summary_lines = [
        "# 41-gene DINP evidence and target-prioritization audit",
        "",
        f"Generated (UTC): {now_utc()}",
        "",
        "## Scope",
        "",
        "This audit decomposes the fresh 41-gene DINP–CRC intersection into CRC reference relevance, DINP source support, and GO/Reactome/KEGG pathway recurrence. It does not change the intersection, recompute enrichment, or introduce a hidden composite score.",
        "",
        f"- Intersection genes: **{len(rows)}**",
        f"- Genes supported by at least two DINP source families (CTD, CompTox, T3DB): **{multi_dinp_count}**",
        f"- Genes classified as strong cross-resource pathway anchors: **{strong_anchor_count}**",
        "",
        "## DINP source-family support",
        "",
        "ToxCast and Tox21 are shown separately and also collapsed into one CompTox family flag for source-family counting. The raw source count is retained and is not interpreted as independent biological truth.",
        "",
        "| DINP source-family count | Genes |",
        "|---:|---:|",
    ]
    for count in sorted(source_family_counts):
        summary_lines.append(f"| {count} | {source_family_counts[count]} |")
    summary_lines.extend(
        [
            "",
            "## Why RXRA ranks first",
            "",
            f"RXRA has CRC article rank **{rxra['crc_article_rank']}**, CRC relevance score **{rxra['crc_relevance_score']:.2f}**, DINP support from **{rxra['DINP_evidence_profile']}** ({rxra['DINP_source_family_count']} source families), and pathway recurrence of **{rxra['pathway_recurrence_total']}** terms across **{rxra['pathway_source_count']}** resources. Its pathway recurrence is a breadth proxy, not an annotation-density correction; this audit does not claim that RXRA is the strongest direct DINP target.",
            "",
            "## Prioritization decomposition",
            "",
            "The original ranking is pathway-context prioritization. It is best read together with the independent DINP evidence profile and CRC relevance columns below.",
            "",
            "| Rank | Gene | CRC rank | CRC score | DINP families | DINP profile | GO | Reactome | KEGG | Pathway anchor |",
            "|---:|---|---:|---:|---:|---|---:|---:|---:|---|",
        ]
    )
    for row in rows[:20]:
        summary_lines.append(
            f"| {row['original_biological_priority_rank']} | **{row['gene_symbol']}** | {row['crc_article_rank']} | {row['crc_relevance_score']:.2f} | {row['DINP_source_family_count']} | {row['DINP_evidence_profile'] or '—'} | {row['GO_representative_count']} | {row['Reactome_term_count']} | {row['KEGG_term_count']} | {row['pathway_anchor_class']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "RXRA, PPARA, PPARG, PPARD, PTGS2, and related genes should be described as cross-resource pathway anchors or prioritized candidates only after this decomposition. A high GO/Reactome/KEGG recurrence can reflect annotation breadth. Direct DINP binding, causal mediation, and in-vivo relevance require independent evidence and are not established here.",
            "",
            "## Files",
            "",
            "- `dinp_crc_41_gene_evidence_matrix.csv`: full source/relevance/pathway decomposition for all 41 genes.",
            "- `prioritization_decomposition.csv`: compact table for reviewer-facing rank decomposition.",
        ]
    )
    summary_path = OUT_DIR / "evidence_audit_summary.md"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest = {
        "generated_utc": now_utc(),
        "script": str(Path(__file__).resolve()),
        "inputs": {
            "intersection": {"path": str(INTERSECTION.resolve()), "sha256": sha256_file(INTERSECTION)},
            "exposure_matrix": {"path": str(EXPOSURE.resolve()), "sha256": sha256_file(EXPOSURE)},
            "target_prioritization": {"path": str(TARGET.resolve()), "sha256": sha256_file(TARGET)},
        },
        "gene_count": len(rows),
        "dinp_source_family_definition": "CTD; CompTox=(ToxCast OR Tox21); T3DB",
        "pathway_recurrence_definition": "GO cleaned representatives + significant Reactome + significant KEGG term membership",
        "annotation_density_correction_performed": False,
        "new_statistical_test_performed": False,
        "outputs": {
            name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
            for name, path in {
                "evidence_matrix": matrix_path,
                "prioritization_decomposition": decomposition_path,
                "summary": summary_path,
            }.items()
        },
    }
    manifest_path = OUT_DIR / "evidence_audit_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print("41-GENE DINP EVIDENCE AUDIT: PASS")
    print(f"Genes audited: {len(rows)}")
    print(f"Multi-source-family DINP support: {multi_dinp_count}")
    print(f"Strong pathway anchors: {strong_anchor_count}")
    print(f"Evidence matrix: {matrix_path.resolve()}")
    print(f"Summary: {summary_path.resolve()}")


if __name__ == "__main__":
    main()
