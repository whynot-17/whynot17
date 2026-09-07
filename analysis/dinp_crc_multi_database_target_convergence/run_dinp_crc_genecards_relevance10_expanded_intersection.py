#!/usr/bin/env python3
"""Intersect the expanded DINP exposure matrix with GeneCards CRC genes.

The exposure side is the source-preserving CTD + ToxCast + Tox21 + T3DB
matrix.  The CRC side is deliberately GeneCards-only and uses the frozen
RelevanceScore >= 10 rule.  Open Targets and other disease databases are not
loaded by this script.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
OUT = HERE / "outputs"
EXPOSURE_MATRIX = OUT / "dinp_exposure_gene_matrix_ctd_toxcast_tox21_t3db.csv"
GENECARDS = OUT / "source_records" / "genecards_crc_top2000_reference.csv"

THRESHOLD = 10.0


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repo_relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if not EXPOSURE_MATRIX.exists():
        raise FileNotFoundError(EXPOSURE_MATRIX)
    if not GENECARDS.exists():
        raise FileNotFoundError(GENECARDS)

    exposure_rows = read_csv(EXPOSURE_MATRIX)
    cards_rows = read_csv(GENECARDS)
    if not exposure_rows:
        raise ValueError("Expanded DINP exposure matrix is empty")

    exposure_by_gene: dict[str, dict[str, str]] = {}
    for row in exposure_rows:
        gene = clean(row.get("gene_symbol")).upper()
        if not gene:
            continue
        if gene in exposure_by_gene:
            raise ValueError(f"Duplicate gene in exposure matrix: {gene}")
        exposure_by_gene[gene] = row

    cards_by_gene: dict[str, dict[str, str]] = {}
    cards_ge_threshold = 0
    for row in cards_rows:
        gene = clean(row.get("gene_symbol") or row.get("GeneSymbol")).upper()
        if not gene:
            continue
        score_raw = clean(row.get("RelevanceScore") or row.get("gene_cards_relevance_score"))
        try:
            score = float(score_raw)
        except ValueError:
            continue
        if score >= THRESHOLD:
            cards_ge_threshold += 1
            if gene in cards_by_gene:
                raise ValueError(f"Duplicate GeneCards gene above threshold: {gene}")
            cards_by_gene[gene] = row

    intersection_genes = sorted(set(exposure_by_gene) & set(cards_by_gene))
    rows: list[dict[str, Any]] = []
    matrix_fields = [field for field in exposure_rows[0] if field != "gene_symbol"]
    output_fields = [
        "gene_symbol",
        "gene_cards_rank",
        "gene_cards_relevance_score",
        "gene_cards_knowledge_score",
        "gene_cards_gene_name",
        "gene_cards_gene_type",
        *matrix_fields,
        "rule",
    ]

    for gene in intersection_genes:
        cards = cards_by_gene[gene]
        exposure = exposure_by_gene[gene]
        rows.append(
            {
                "gene_symbol": gene,
                "gene_cards_rank": clean(cards.get("GeneCards_Rank") or cards.get("Rank")),
                "gene_cards_relevance_score": clean(cards.get("RelevanceScore") or cards.get("gene_cards_relevance_score")),
                "gene_cards_knowledge_score": clean(cards.get("KnowledgeScore") or cards.get("gene_cards_knowledge_score")),
                "gene_cards_gene_name": clean(cards.get("GeneName") or cards.get("Description")),
                "gene_cards_gene_type": clean(cards.get("GeneType") or cards.get("Category")),
                **{field: exposure.get(field, "") for field in matrix_fields},
                "rule": "expanded DINP exposure matrix AND GeneCards CRC RelevanceScore >= 10",
            }
        )

    rows.sort(key=lambda row: (-float(row["gene_cards_relevance_score"]), row["gene_symbol"]))
    output = OUT / "dinp_crc_genecards_relevance10_expanded_intersection.csv"
    write_csv(output, rows, output_fields)

    added_by_expansion = sorted(
        set(intersection_genes)
        & {gene for gene, row in exposure_by_gene.items() if row.get("CTD") == "0"}
    )
    source_files = {
        "expanded_dinp_exposure_matrix": {
            "path": repo_relative(EXPOSURE_MATRIX),
            "sha256": sha256_file(EXPOSURE_MATRIX),
        },
        "genecards_crc_top2000_reference": {
            "path": repo_relative(GENECARDS),
            "sha256": sha256_file(GENECARDS),
        },
    }
    manifest = {
        "generated_at": utc_now(),
        "analysis": "Expanded DINP exposure matrix ∩ GeneCards CRC threshold set",
        "exposure_sources": ["CTD", "ToxCast", "Tox21", "T3DB"],
        "disease_source": "GeneCards only",
        "gene_cards_rule": {
            "field": "RelevanceScore",
            "operator": ">=",
            "value": THRESHOLD,
            "scope": "ordinary CRC top-2000 reference archived locally; not a full-ranking export",
        },
        "counts": {
            "expanded_exposure_genes": len(exposure_by_gene),
            "genecards_input_rows": len(cards_rows),
            "genecards_rows_relevance_ge_10": cards_ge_threshold,
            "expanded_exposure_genecards_intersection": len(rows),
            "new_intersection_genes_from_non_ctd_layers": len(added_by_expansion),
        },
        "new_intersection_genes_from_non_ctd_layers": added_by_expansion,
        "source_files": source_files,
        "open_targets_used": False,
        "interpretation": "This is a GeneCards-threshold intersection, not a causal target set or a merged biological truth score.",
        "output": repo_relative(output),
        "output_sha256": None,
    }
    output_hash = sha256_file(output)
    manifest["output_sha256"] = output_hash
    manifest_path = OUT / "dinp_crc_genecards_relevance10_expanded_intersection_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = [
        "# Expanded DINP–CRC GeneCards threshold intersection",
        "",
        f"Generated: `{manifest['generated_at']}`",
        "",
        "## Frozen rule",
        "",
        "The exposure side is the 93-gene source-preserving DINP matrix from CTD, EPA CompTox/ToxCast/Tox21, and T3DB. The CRC side is GeneCards only, filtered at `RelevanceScore >= 10`. Open Targets is not used.",
        "",
        "## Counts",
        "",
        "| Quantity | Count |",
        "|---|---:|",
        f"| Expanded DINP exposure genes | {len(exposure_by_gene)} |",
        f"| GeneCards input rows | {len(cards_rows)} |",
        f"| GeneCards rows with RelevanceScore ≥10 | {cards_ge_threshold} |",
        f"| Expanded DINP ∩ GeneCards ≥10 | {len(rows)} |",
        f"| Added by non-CTD layers | {len(added_by_expansion)} |",
        "",
        "## Added genes from the expanded exposure side",
        "",
        f"`{', '.join(added_by_expansion) if added_by_expansion else 'None'}`",
        "",
        "The archived GeneCards input is the ordinary CRC top-2000 reference; this output is reproducible within that archived input and is not presented as a full-ranking GeneCards result.",
        "",
        "## Output",
        "",
        "- `dinp_crc_genecards_relevance10_expanded_intersection.csv`",
        "- `dinp_crc_genecards_relevance10_expanded_intersection_manifest.json`",
    ]
    summary_path = OUT / "dinp_crc_genecards_relevance10_expanded_intersection_summary.md"
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("DINP–CRC EXPANDED GENECARDS INTERSECTION: PASS")
    print(f"Expanded DINP exposure genes: {len(exposure_by_gene)}")
    print(f"GeneCards RelevanceScore >= 10 rows: {cards_ge_threshold}")
    print(f"Intersection genes: {len(rows)}")
    print(f"Added by non-CTD layers: {', '.join(added_by_expansion) if added_by_expansion else 'None'}")
    print(f"Output: {output.resolve()}")


if __name__ == "__main__":
    main()
