#!/usr/bin/env python3
"""Filter the DINP CTD gene set against a GeneCards CRC score threshold.

This is an intentionally narrow, source-preserving analysis.  Open Targets,
DisGeNET, CRC outcomes, enrichment statistics, and any downstream ranking are
not loaded or used.  The result is a CTD DINP gene ∩ GeneCards CRC set defined
by the prespecified GeneCards Relevance Score threshold.
"""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
SOURCE = OUT / "source_records"
CTD_MATRIX = OUT / "dinp_exposure_gene_matrix.csv"
GENECARDS = SOURCE / "genecards_crc_top2000_reference.csv"
THRESHOLD = 10.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def first(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = (row.get(name) or "").strip()
        if value:
            return value
    return ""


def main() -> None:
    if not CTD_MATRIX.exists():
        raise FileNotFoundError(f"Missing CTD exposure matrix: {CTD_MATRIX}")
    if not GENECARDS.exists():
        raise FileNotFoundError(f"Missing GeneCards input: {GENECARDS}")

    ctd_rows = read_csv(CTD_MATRIX)
    ctd_rows = [r for r in ctd_rows if (r.get("CTD") or "0").strip() == "1"]
    ctd_by_gene = {
        first(r, "gene_symbol", "GeneSymbol").upper(): r
        for r in ctd_rows
        if first(r, "gene_symbol", "GeneSymbol")
    }
    if len(ctd_by_gene) != 86:
        raise ValueError(f"Expected 86 unique CTD DINP genes; found {len(ctd_by_gene)}")

    cards_rows = read_csv(GENECARDS)
    cards_by_gene: dict[str, dict[str, str]] = {}
    for row in cards_rows:
        symbol = first(row, "GeneSymbol", "gene_symbol", "Symbol").upper()
        if not symbol:
            continue
        score_raw = first(row, "RelevanceScore", "relevance_score")
        try:
            score = float(score_raw)
        except ValueError:
            continue
        if score >= THRESHOLD and symbol not in cards_by_gene:
            cards_by_gene[symbol] = row

    hits: list[dict[str, str]] = []
    for symbol, cards in ctd_by_gene.items():
        if symbol not in cards_by_gene:
            continue
        ctd = cards_by_gene[symbol]
        hits.append(
            {
                "gene_symbol": symbol,
                "gene_cards_rank": first(ctd, "GeneCards_Rank", "Rank"),
                "gene_cards_relevance_score": first(ctd, "RelevanceScore", "relevance_score"),
                "gene_cards_knowledge_score": first(ctd, "KnowledgeScore", "knowledge_score"),
                "gene_cards_gene_name": first(ctd, "GeneName", "Description"),
                "gene_cards_gene_type": first(ctd, "GeneType", "Category"),
                "ctd_raw_row_count": first(cards, "CTD_raw_row_count"),
                "ctd_direction": first(cards, "CTD_direction"),
                "ctd_action": first(cards, "CTD_action"),
                "ctd_single_chemical_evidence": first(cards, "CTD_single_chemical_evidence"),
                "ctd_cotreatment_evidence": first(cards, "CTD_cotreatment_evidence"),
                "rule": "CTD DINP human gene AND GeneCards CRC RelevanceScore >= 10",
            }
        )

    hits.sort(key=lambda r: (-float(r["gene_cards_relevance_score"]), r["gene_symbol"]))
    OUT.mkdir(parents=True, exist_ok=True)
    result_path = OUT / "dinp_crc_genecards_relevance10_intersection.csv"
    with result_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(hits[0]) if hits else ["gene_symbol"])
        writer.writeheader()
        writer.writerows(hits)

    now = datetime.now(timezone.utc).isoformat()
    manifest = {
        "analysis": "DINP CTD ∩ GeneCards CRC threshold intersection",
        "generated_at_utc": now,
        "gene_cards_threshold": {"field": "RelevanceScore", "operator": ">=", "value": THRESHOLD},
        "inputs": {
            "ctd_exposure_matrix": str(CTD_MATRIX),
            "ctd_exposure_matrix_sha256": sha256(CTD_MATRIX),
            "genecards_input": str(GENECARDS),
            "genecards_input_sha256": sha256(GENECARDS),
            "genecards_scope": "archived ordinary CRC search top-2000 reference; not a full ranking export",
        },
        "firewall": {
            "open_targets_loaded": False,
            "disgenet_loaded": False,
            "crc_outcome_data_loaded": False,
            "p_values_or_fdr_used": False,
            "downstream_enrichment_used": False,
        },
        "counts": {
            "ctd_unique_dinp_human_genes": len(ctd_by_gene),
            "genecards_rows": len(cards_rows),
            "genecards_rows_relevance_ge_10": len(cards_by_gene),
            "intersection_genes": len(hits),
        },
        "output": {"path": str(result_path), "sha256": sha256(result_path)},
    }
    (OUT / "dinp_crc_genecards_relevance10_intersection_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# DINP–CRC GeneCards threshold intersection",
        "",
        f"- GeneCards rule: `RelevanceScore >= {THRESHOLD:g}`",
        "- Exposure-side input: CTD human DINP genes only",
        "- Disease-side input: GeneCards only",
        "- Open Targets: not used",
        "- DisGeNET: not used",
        "- Outcome statistics, P values, FDR, enrichment, and downstream ranking: not used",
        "",
        "## Counts",
        "",
        f"- CTD DINP human genes: **{len(ctd_by_gene)}**",
        f"- GeneCards input rows: **{len(cards_rows)}**",
        f"- GeneCards rows with RelevanceScore >= 10: **{len(cards_by_gene)}**",
        f"- CTD ∩ GeneCards (RelevanceScore >= 10): **{len(hits)} genes**",
        "",
        "## Scope note",
        "",
        "The available GeneCards file is the archived ordinary CRC top-2000 reference. "
        "This thresholded result is therefore reproducible within that archived input; "
        "it is not presented as a full-ranking export until a full raw GeneCards file is available.",
        "",
        "## Output",
        "",
        "- `dinp_crc_genecards_relevance10_intersection.csv`",
        "- `dinp_crc_genecards_relevance10_intersection_manifest.json`",
    ]
    (OUT / "dinp_crc_genecards_relevance10_intersection_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest["counts"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
