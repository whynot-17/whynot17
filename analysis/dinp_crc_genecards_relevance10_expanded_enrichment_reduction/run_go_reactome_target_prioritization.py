#!/usr/bin/env python3
"""Prioritize DINP–CRC candidate genes after GO–Reactome alignment.

The output deliberately separates biological convergence from downstream
macrophage/network/structure follow-up.  No new enrichment test is run and no
docking/MD result is used as biological proof.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"
ALIGNED_DIR = OUTPUTS / "go_reactome_aligned"
OUT_DIR = OUTPUTS / "go_reactome_target_prioritization"

ALIGNED_GLOBAL = ALIGNED_DIR / "go_reactome_global_driver_contribution.csv"
ALIGNED_TERMS = ALIGNED_DIR / "go_reactome_aligned_terms.csv"
THEME_ALIGNMENT = ALIGNED_DIR / "go_reactome_theme_alignment.csv"
MACROPHAGE = ROOT.parent / "dinp_crc_81gene_macrophage_driver_decomposition" / "outputs" / "macrophage_driver_candidates.csv"
NETWORK = ROOT.parent / "dinp_crc_81gene_macrophage_network_prioritization" / "outputs" / "network_target_evidence_matrix.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def num(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def integer(value: str | None, default: int = 0) -> int:
    try:
        return int(float(value)) if value not in (None, "") else default
    except ValueError:
        return default


def main() -> None:
    global_rows = read_csv(ALIGNED_GLOBAL)
    aligned_rows = read_csv(ALIGNED_TERMS)
    theme_rows = read_csv(THEME_ALIGNMENT)
    macrophage_rows = read_csv(MACROPHAGE) if MACROPHAGE.exists() else []
    network_rows = read_csv(NETWORK) if NETWORK.exists() else []

    if len(global_rows) != 18:
        raise RuntimeError(f"Expected 18 aligned genes, found {len(global_rows)}")

    theme_gene_sources: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in aligned_rows:
        for gene in (row.get("intersection_genes") or "").split(";"):
            gene = gene.strip()
            if gene:
                theme_gene_sources[row["theme_id"]][row["source"]].add(gene)

    theme_labels = {row["theme_id"]: row["theme_label"] for row in theme_rows}
    theme_order = {row["theme_id"]: integer(row.get("theme_order"), 999) for row in theme_rows}
    shared_theme_ids: set[str] = set()
    for theme_id, sources in theme_gene_sources.items():
        if sources.get("GO:BP") and sources.get("REAC"):
            shared_theme_ids.add(theme_id)

    macrophage_by_gene = {row["gene_symbol"]: row for row in macrophage_rows}
    network_by_gene = {row["gene_symbol"]: row for row in network_rows}

    prioritization_rows: list[dict[str, Any]] = []
    for row in global_rows:
        gene = row["gene_symbol"]
        go_n = integer(row.get("GO_representative_count"))
        reac_n = integer(row.get("Reactome_term_count"))
        go_themes = {
            theme_id
            for theme_id, sources in theme_gene_sources.items()
            if gene in sources.get("GO:BP", set())
        }
        reac_themes = {
            theme_id
            for theme_id, sources in theme_gene_sources.items()
            if gene in sources.get("REAC", set())
        }
        shared_gene_themes = sorted(go_themes & reac_themes, key=lambda value: (theme_order.get(value, 999), value))
        shared_theme_labels = "; ".join(theme_labels[theme_id] for theme_id in shared_gene_themes)
        role_labels: list[str] = []
        if "lipid_eicosanoid" in reac_themes:
            role_labels.append("direct lipid/eicosanoid pathway context")
        if "signaling_receptor" in reac_themes:
            role_labels.append("nuclear-receptor pathway context")
        if "matrix_migration" in reac_themes:
            role_labels.append("matrix-remodeling pathway context")
        if not role_labels:
            role_labels.append("GO-only pathway context")

        macrophage = macrophage_by_gene.get(gene, {})
        network = network_by_gene.get(gene, {})
        both_source = bool(go_n and reac_n)
        shared_theme_count = len(shared_gene_themes)

        # This is a transparent lexicographic evidence ordering, not a fitted
        # composite score.  It privileges same-theme cross-resource support,
        # then Reactome recurrence, then GO recurrence.
        biological_key = (
            -shared_theme_count,
            -int(both_source),
            -reac_n,
            -go_n,
            -len(go_themes | reac_themes),
            gene,
        )
        prioritization_rows.append(
            {
                "gene_symbol": gene,
                "GO_representative_count": go_n,
                "Reactome_term_count": reac_n,
                "aligned_term_count": go_n + reac_n,
                "GO_theme_count": len(go_themes),
                "Reactome_theme_count": len(reac_themes),
                "shared_theme_count": shared_theme_count,
                "shared_theme_labels": shared_theme_labels,
                "source_status": "GO+Reactome_same_theme" if shared_theme_count else ("GO+Reactome" if both_source else "GO_only"),
                "biological_role": "; ".join(role_labels),
                "GO_recurrence_fraction": f"{go_n / 43:.4f}",
                "Reactome_recurrence_fraction": f"{reac_n / 6:.4f}",
                "macrophage_driver_rank": macrophage.get("driver_priority_rank", ""),
                "macrophage_paired_t_BH_FDR": macrophage.get("paired_t_BH_FDR", ""),
                "macrophage_cohen_dz": macrophage.get("cohen_dz_paired", ""),
                "macrophage_detection_fraction": macrophage.get("tumor_cell_detection_fraction", ""),
                "network_priority_rank": network.get("network_priority_rank", ""),
                "candidate_role": network.get("candidate_role", ""),
                "docking_protein_eligible": network.get("docking_protein_eligible", ""),
                "best_pdb_id": network.get("best_pdb_id", ""),
                "pdb_structure_count": network.get("pdb_structure_count", ""),
                "known_measured_target_evidence": network.get("known_measured_target_evidence", ""),
                "_biological_key": biological_key,
            }
        )

    prioritization_rows.sort(key=lambda row: row["_biological_key"])
    for rank, row in enumerate(prioritization_rows, start=1):
        row["biological_priority_rank"] = rank
        row.pop("_biological_key", None)
    priority_by_gene = {row["gene_symbol"]: row for row in prioritization_rows}

    fields = [
        "biological_priority_rank",
        "gene_symbol",
        "GO_representative_count",
        "Reactome_term_count",
        "aligned_term_count",
        "GO_theme_count",
        "Reactome_theme_count",
        "shared_theme_count",
        "shared_theme_labels",
        "source_status",
        "biological_role",
        "GO_recurrence_fraction",
        "Reactome_recurrence_fraction",
        "macrophage_driver_rank",
        "macrophage_paired_t_BH_FDR",
        "macrophage_cohen_dz",
        "macrophage_detection_fraction",
        "network_priority_rank",
        "candidate_role",
        "docking_protein_eligible",
        "best_pdb_id",
        "pdb_structure_count",
        "known_measured_target_evidence",
    ]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUT_DIR / "go_reactome_target_prioritization.csv", prioritization_rows, fields)

    shared_rows = [row for row in prioritization_rows if row["source_status"] == "GO+Reactome_same_theme"]
    shared_fields = fields
    write_csv(OUT_DIR / "go_reactome_cross_resource_target_shortlist.csv", shared_rows, shared_fields)

    # Keep the seven-gene macrophage/STRING structural branch intact even
    # when a candidate is not part of the cleaned 18-gene GO–Reactome input
    # universe (currently PTGER4 and PTGES3).  Such rows are explicitly
    # labelled as outside the aligned biological universe rather than being
    # silently dropped or assigned zero biological evidence.
    structural_candidates = sorted(set(macrophage_by_gene) | set(network_by_gene))
    structural_rows: list[dict[str, Any]] = []
    for gene in structural_candidates:
        if gene in priority_by_gene:
            structural_rows.append(dict(priority_by_gene[gene]))
            continue
        macrophage = macrophage_by_gene.get(gene, {})
        network = network_by_gene.get(gene, {})
        extra = {field: "" for field in fields}
        extra.update(
            {
                "biological_priority_rank": "",
                "gene_symbol": gene,
                "GO_representative_count": 0,
                "Reactome_term_count": 0,
                "aligned_term_count": 0,
                "GO_theme_count": 0,
                "Reactome_theme_count": 0,
                "shared_theme_count": 0,
                "shared_theme_labels": "",
                "source_status": "outside_aligned_18_gene_universe",
                "biological_role": "No GO–Reactome evidence in the aligned 18-gene universe",
                "GO_recurrence_fraction": "0.0000",
                "Reactome_recurrence_fraction": "0.0000",
                "macrophage_driver_rank": macrophage.get("driver_priority_rank", ""),
                "macrophage_paired_t_BH_FDR": macrophage.get("paired_t_BH_FDR", ""),
                "macrophage_cohen_dz": macrophage.get("cohen_dz_paired", ""),
                "macrophage_detection_fraction": macrophage.get("tumor_cell_detection_fraction", ""),
                "network_priority_rank": network.get("network_priority_rank", ""),
                "candidate_role": network.get("candidate_role", ""),
                "docking_protein_eligible": network.get("docking_protein_eligible", ""),
                "best_pdb_id": network.get("best_pdb_id", ""),
                "pdb_structure_count": network.get("pdb_structure_count", ""),
                "known_measured_target_evidence": network.get("known_measured_target_evidence", ""),
            }
        )
        structural_rows.append(extra)
    structural_rows.sort(
        key=lambda row: (
            integer(row.get("network_priority_rank"), 999),
            integer(row.get("macrophage_driver_rank"), 999),
            row["gene_symbol"],
        )
    )
    for rank, row in enumerate(structural_rows, start=1):
        row = dict(row)
        row["structural_followup_rank"] = rank
        # Rebuild in place so the output keeps the same dict instances only
        # conceptually; no source rows are mutated.
        structural_rows[rank - 1] = row
    structural_fields = ["structural_followup_rank"] + fields
    write_csv(OUT_DIR / "go_reactome_structural_followup_prioritization.csv", structural_rows, structural_fields)

    top_biology = prioritization_rows[:8]
    summary = [
        "# GO–Reactome target prioritization",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Boundary",
        "",
        "This is an evidence-ordering layer over the existing 18-gene DINP–CRC intersection. No enrichment P values or FDR values were recalculated. Database recurrence is not causality, and structural/docking availability is not evidence of in-vivo binding.",
        "",
        "Biological priority is ordered lexicographically by: (1) number of same-theme GO–Reactome concordances, (2) presence in both sources, (3) Reactome term recurrence, (4) GO representative-term recurrence, and (5) total theme breadth. This avoids a hidden weighted composite score.",
        "",
        "## Cross-resource biological anchors",
        "",
        f"- 18 total genes; **{len(shared_rows)}** have GO and Reactome support in the same theme.",
        f"- Shared GO–Reactome themes: **{len(shared_theme_ids)}** — lipid/eicosanoid metabolism, nuclear-receptor signaling, and matrix remodeling.",
        "",
        "| Rank | Gene | GO | Reactome | Same-theme concordance | Role |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for row in top_biology:
        summary.append(
            f"| {row['biological_priority_rank']} | **{row['gene_symbol']}** | "
            f"{row['GO_representative_count']} | {row['Reactome_term_count']} | "
            f"{row['shared_theme_count']} | {row['biological_role']} |"
        )
    summary.extend(
        [
            "",
            "### Interpretation",
            "",
            "- **PPARD / PPARG** are the strongest multi-theme cross-resource anchors because they recur in lipid/eicosanoid and nuclear-receptor themes.",
            "- **PTGS2** is the strongest direct lipid/eicosanoid recurrence anchor, with support across all four informative lipid/eicosanoid Reactome terms and extensive GO recurrence.",
            "- **MMP9 / MMP2 / TIMP1** form the matrix-remodeling cross-resource branch; MMP9 additionally has macrophage and STRING-network support.",
            "- **CYP1A2 / CYP3A4 / PTGS1** are direct lipid/eicosanoid pathway-context candidates, but their prioritization is pathway-contextual rather than proof of a DINP-specific target interaction.",
            "",
            "## Structural follow-up branch",
            "",
            "The structural branch is reported separately because the existing macrophage/STRING workflow covers seven candidates and is not the same universe as the 18-gene GO–Reactome alignment.",
            "",
            "| Follow-up rank | Gene | Macrophage driver rank | STRING/network rank | Biological source status | Role |",
            "|---:|---|---:|---:|---|---|",
        ]
    )
    for row in structural_rows:
        summary.append(
            f"| {row['structural_followup_rank']} | **{row['gene_symbol']}** | "
            f"{row['macrophage_driver_rank'] or '—'} | {row['network_priority_rank'] or '—'} | "
            f"{row['source_status']} | {row['candidate_role'] or row['biological_role']} |"
        )
    summary.extend(
        [
            "",
            "## Recommended use",
            "",
            "For the biological mechanism narrative, carry forward the cross-resource anchors first: **PPARD/PPARG–PTGS2**, with the **MMP9/MMP2/TIMP1 matrix branch** as a complementary tumor-microenvironment axis.",
            "",
            "For protein-level follow-up, retain MMP9 and STAT3 as network-bridge candidates and PTGER4/PTGES3 as direct prostaglandin-context candidates. NEAT1 remains a non-protein state node and is not docking-eligible.",
            "",
            "These are prioritization recommendations only. They do not establish that DINP binds any target or that any target mediates the epidemiologic association.",
            "",
            "## Files",
            "",
            "- `go_reactome_target_prioritization.csv`: all 18 genes with aligned evidence and follow-up annotations.",
            "- `go_reactome_cross_resource_target_shortlist.csv`: genes with same-theme GO–Reactome concordance.",
            "- `go_reactome_structural_followup_prioritization.csv`: macrophage/network/structure follow-up branch.",
        ]
    )
    (OUT_DIR / "go_reactome_target_prioritization_summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_aligned_global": str(ALIGNED_GLOBAL),
        "input_aligned_terms": str(ALIGNED_TERMS),
        "input_theme_alignment": str(THEME_ALIGNMENT),
        "candidate_gene_count": len(global_rows),
        "same_theme_cross_resource_gene_count": len(shared_rows),
        "shared_theme_count": len(shared_theme_ids),
        "structural_followup_candidate_count": len(structural_rows),
        "biological_ordering": [
            "same_theme_concordance_count descending",
            "both-source support descending",
            "Reactome recurrence descending",
            "GO recurrence descending",
            "theme breadth descending",
            "gene symbol ascending",
        ],
        "statistical_boundary": "No new enrichment tests; no P/FDR changes; prioritization is descriptive.",
        "outputs": [
            "go_reactome_target_prioritization.csv",
            "go_reactome_cross_resource_target_shortlist.csv",
            "go_reactome_structural_followup_prioritization.csv",
            "go_reactome_target_prioritization_summary.md",
        ],
    }
    (OUT_DIR / "go_reactome_target_prioritization_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print("GO–Reactome target prioritization: PASS")
    print(f"Candidate genes: {len(global_rows)}")
    print(f"Same-theme cross-resource genes: {len(shared_rows)}")
    print(f"Structural follow-up candidates: {len(structural_rows)}")
    print(f"Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
