#!/usr/bin/env python3
"""Align cleaned GO:BP representatives with informative Reactome terms.

This is an additive interpretation layer.  It does not recompute enrichment,
change any P values/FDR values, or treat GO and Reactome terms as independent
tests.  Reactome's root term is explicitly excluded from biological themes.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"
GO_DIR = OUTPUTS / "go_bp_cleaning"
OUT_DIR = OUTPUTS / "go_reactome_aligned"
GO_REPRESENTATIVES = GO_DIR / "go_bp_clean_representatives.csv"
RAW_TERMS = OUTPUTS / "primary_106_terms_with_genes.csv"


THEMES: list[tuple[str, str]] = [
    ("lipid_eicosanoid", "Lipid / fatty-acid / eicosanoid metabolism"),
    ("inflammation_defense", "Inflammation / defense / stimulus response"),
    ("chemical_oxygen_stress", "Chemical / oxygen / abiotic stress response"),
    ("apoptosis_cell_death", "Apoptosis / programmed cell death"),
    ("mirna_regulation", "miRNA transcription / metabolism"),
    ("signaling_receptor", "Signaling / receptor / hormone response"),
    ("matrix_migration", "Matrix remodeling / cell migration"),
    ("development_reproduction", "Development / reproduction / circulation"),
    ("general_metabolism", "General metabolism / biosynthesis / storage"),
]
THEME_LABEL = dict(THEMES)
THEME_ORDER = {theme_id: index + 1 for index, (theme_id, _label) in enumerate(THEMES)}

# Explicit cross-resource alignment.  These are the same semantic assignments
# already present in the prior mixed-source audit, but the root term is not
# carried forward and no new Reactome hierarchy is inferred here.
REACTOME_THEME_MAP = {
    "REAC:R-HSA-556833": "lipid_eicosanoid",  # Metabolism of lipids
    "REAC:R-HSA-8978868": "lipid_eicosanoid",  # Fatty acid metabolism
    "REAC:R-HSA-9018678": "lipid_eicosanoid",  # SPM biosynthesis
    "REAC:R-HSA-9018677": "lipid_eicosanoid",  # DHA-derived SPM biosynthesis
    "REAC:R-HSA-383280": "signaling_receptor",  # Nuclear receptor pathway
    "REAC:R-HSA-1592389": "matrix_migration",  # MMP activation
}
REACTOME_ROOT_ID = "REAC:0000000"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def genes(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item.strip() for item in value.split(";") if item.strip()}


def fmt_genes(values: set[str]) -> str:
    return ";".join(sorted(values))


def fmt_float(value: float | None) -> str:
    return "" if value is None else f"{value:.4f}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_gene_counts(rows: list[dict[str, Any]], source: str) -> dict[str, Counter[str]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        if row["source"] != source:
            continue
        for gene in row["intersection_genes_set"]:
            counts[row["theme_id"]][gene] += 1
    return counts


def driver_string(
    counts: Counter[str],
    theme_term_count: int,
    limit: int = 6,
) -> str:
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:limit]
    return "; ".join(f"{gene}:{count}/{theme_term_count}" for gene, count in ordered)


def main() -> None:
    go_rows_raw = read_csv(GO_REPRESENTATIVES)
    raw_rows = read_csv(RAW_TERMS)

    if len(go_rows_raw) != 43:
        raise RuntimeError(f"Expected 43 cleaned GO representatives, found {len(go_rows_raw)}")

    aligned: list[dict[str, Any]] = []
    for row in go_rows_raw:
        theme_id = row["term_theme_id"]
        if theme_id not in THEME_LABEL:
            raise RuntimeError(f"Unknown GO theme: {theme_id}")
        aligned.append(
            {
                "source": "GO:BP",
                "source_status": "GO_clean_representative",
                "theme_id": theme_id,
                "theme_label": THEME_LABEL[theme_id],
                "term_id": row["representative_term_id"],
                "term_name": row["term_name"],
                "global_fdr": row["global_fdr"],
                "intersection_size": row["intersection_size"],
                "intersection_genes": row["intersection_genes"],
                "intersection_genes_set": genes(row["intersection_genes"]),
                "alignment_basis": "GO semantic cleaning representative",
            }
        )

    reactome_rows = [row for row in raw_rows if row.get("source") == "REAC"]
    reactome_root_rows = [row for row in reactome_rows if row.get("term_id") == REACTOME_ROOT_ID]
    informative_reactome: list[dict[str, Any]] = []
    for row in reactome_rows:
        term_id = row["term_id"]
        if term_id == REACTOME_ROOT_ID:
            continue
        if term_id not in REACTOME_THEME_MAP:
            raise RuntimeError(f"Unmapped informative Reactome term: {term_id}")
        theme_id = REACTOME_THEME_MAP[term_id]
        informative_reactome.append(
            {
                "source": "REAC",
                "source_status": "Reactome_informative_term",
                "theme_id": theme_id,
                "theme_label": THEME_LABEL[theme_id],
                "term_id": term_id,
                "term_name": row["term_name"],
                "global_fdr": row["global_fdr"],
                "intersection_size": row["intersection_size"],
                "intersection_genes": row["intersection_genes"],
                "intersection_genes_set": genes(row["intersection_genes"]),
                "alignment_basis": "pre-specified Reactome-to-GO theme map",
            }
        )
    aligned.extend(informative_reactome)

    if len(reactome_rows) != 7 or len(reactome_root_rows) != 1:
        raise RuntimeError(
            f"Expected 7 Reactome rows including 1 root row; found {len(reactome_rows)} "
            f"and {len(reactome_root_rows)} root rows"
        )
    if len(informative_reactome) != 6:
        raise RuntimeError(f"Expected 6 informative Reactome rows, found {len(informative_reactome)}")

    aligned.sort(
        key=lambda row: (
            THEME_ORDER[row["theme_id"]],
            0 if row["source"] == "GO:BP" else 1,
            float(row["global_fdr"]),
            row["term_id"],
        )
    )

    aligned_fieldnames = [
        "source",
        "source_status",
        "theme_id",
        "theme_label",
        "term_id",
        "term_name",
        "global_fdr",
        "intersection_size",
        "intersection_genes",
        "alignment_basis",
    ]
    aligned_export = [
        {key: row[key] for key in aligned_fieldnames}
        for row in aligned
    ]

    go_counts = source_gene_counts(aligned, "GO:BP")
    reactome_counts = source_gene_counts(aligned, "REAC")
    go_terms = Counter(row["theme_id"] for row in aligned if row["source"] == "GO:BP")
    reactome_terms = Counter(row["theme_id"] for row in aligned if row["source"] == "REAC")

    family_rows: list[dict[str, Any]] = []
    driver_rows: list[dict[str, Any]] = []
    for theme_id, theme_label in THEMES:
        go_theme_rows = [row for row in aligned if row["theme_id"] == theme_id and row["source"] == "GO:BP"]
        reactome_theme_rows = [row for row in aligned if row["theme_id"] == theme_id and row["source"] == "REAC"]
        go_union = set().union(*(row["intersection_genes_set"] for row in go_theme_rows)) if go_theme_rows else set()
        reactome_union = (
            set().union(*(row["intersection_genes_set"] for row in reactome_theme_rows))
            if reactome_theme_rows
            else set()
        )
        shared = go_union & reactome_union
        union = go_union | reactome_union
        jaccard = len(shared) / len(union) if union else None
        family_rows.append(
            {
                "theme_order": THEME_ORDER[theme_id],
                "theme_id": theme_id,
                "theme_label": theme_label,
                "GO_representative_count": go_terms[theme_id],
                "Reactome_informative_count": reactome_terms[theme_id],
                "aligned_term_count": go_terms[theme_id] + reactome_terms[theme_id],
                "GO_union_gene_count": len(go_union),
                "Reactome_union_gene_count": len(reactome_union),
                "shared_gene_count": len(shared),
                "union_gene_count": len(union),
                "GO_Reactome_union_Jaccard": fmt_float(jaccard),
                "GO_core_drivers": driver_string(go_counts[theme_id], go_terms[theme_id]),
                "Reactome_core_drivers": driver_string(reactome_counts[theme_id], reactome_terms[theme_id])
                if reactome_terms[theme_id]
                else "",
            }
        )

        all_genes = sorted(go_counts[theme_id] | reactome_counts[theme_id])
        combined = go_counts[theme_id] + reactome_counts[theme_id]
        max_combined = max((combined[gene] for gene in all_genes), default=0)
        for gene in all_genes:
            go_n = go_counts[theme_id][gene]
            reac_n = reactome_counts[theme_id][gene]
            driver_rows.append(
                {
                    "theme_order": THEME_ORDER[theme_id],
                    "theme_id": theme_id,
                    "theme_label": theme_label,
                    "gene_symbol": gene,
                    "GO_representative_count": go_n,
                    "Reactome_term_count": reac_n,
                    "aligned_term_count": go_n + reac_n,
                    "GO_fraction": fmt_float(go_n / go_terms[theme_id]) if go_terms[theme_id] else "",
                    "Reactome_fraction": fmt_float(reac_n / reactome_terms[theme_id]) if reactome_terms[theme_id] else "",
                    "aligned_fraction": fmt_float((go_n + reac_n) / (go_terms[theme_id] + reactome_terms[theme_id]))
                    if go_terms[theme_id] + reactome_terms[theme_id]
                    else "",
                    "source_recurrence": "both" if go_n and reac_n else ("GO_only" if go_n else "Reactome_only"),
                    "combined_top_driver": (go_n + reac_n) == max_combined and max_combined > 0,
                }
            )

    driver_rows.sort(
        key=lambda row: (
            int(row["theme_order"]),
            -int(row["aligned_term_count"]),
            row["gene_symbol"],
        )
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUT_DIR / "go_reactome_aligned_terms.csv", aligned_export, aligned_fieldnames)
    write_csv(
        OUT_DIR / "go_reactome_theme_alignment.csv",
        family_rows,
        list(family_rows[0].keys()),
    )
    write_csv(
        OUT_DIR / "go_reactome_driver_contribution.csv",
        driver_rows,
        list(driver_rows[0].keys()),
    )

    global_go_counts = Counter()
    global_reac_counts = Counter()
    for row in aligned:
        for gene in row["intersection_genes_set"]:
            if row["source"] == "GO:BP":
                global_go_counts[gene] += 1
            else:
                global_reac_counts[gene] += 1
    global_genes = sorted(global_go_counts | global_reac_counts)
    global_driver_rows = []
    for gene in global_genes:
        go_n = global_go_counts[gene]
        reac_n = global_reac_counts[gene]
        global_driver_rows.append(
            {
                "gene_symbol": gene,
                "GO_representative_count": go_n,
                "Reactome_term_count": reac_n,
                "aligned_term_count": go_n + reac_n,
                "theme_count_GO": sum(1 for row in family_rows if gene in go_counts[row["theme_id"]]),
                "theme_count_Reactome": sum(1 for row in family_rows if gene in reactome_counts[row["theme_id"]]),
                "source_recurrence": "both" if go_n and reac_n else ("GO_only" if go_n else "Reactome_only"),
            }
        )
    global_driver_rows.sort(key=lambda row: (-int(row["aligned_term_count"]), row["gene_symbol"]))
    write_csv(
        OUT_DIR / "go_reactome_global_driver_contribution.csv",
        global_driver_rows,
        list(global_driver_rows[0].keys()),
    )

    summary_lines = [
        "# GO:BP–Reactome aligned theme audit",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Scope and statistical boundary",
        "",
        "This is an additive cross-resource alignment layer. The original enrichment P values and global FDR values are unchanged. GO and Reactome terms are retained as source-specific evidence and are not treated as independent tests.",
        "",
        "The 99 significant GO:BP terms were represented by the existing 43-term GO semantic-cleaning output. Reactome contributed 7 significant rows in the original 106-term table; the Reactome root term was excluded, leaving 6 informative terms.",
        "",
        "## Result",
        "",
        f"- Clean GO representatives: **{len(go_rows_raw)}**",
        f"- Reactome rows in original output: **{len(reactome_rows)}**",
        f"- Reactome root terms excluded: **{len(reactome_root_rows)}**",
        f"- Informative Reactome terms aligned: **{len(informative_reactome)}**",
        f"- Aligned source-specific terms: **{len(aligned)}** (43 GO + 6 Reactome)",
        f"- Shared theme families: **{len(THEMES)}**",
        "",
        "## Theme alignment",
        "",
        "| Theme | GO reps | Reactome terms | GO–Reactome shared genes | Jaccard |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in family_rows:
        summary_lines.append(
            f"| {row['theme_label']} | {row['GO_representative_count']} | "
            f"{row['Reactome_informative_count']} | {row['shared_gene_count']} | "
            f"{row['GO_Reactome_union_Jaccard'] or '—'} |"
        )
    summary_lines.extend(
        [
            "",
            "Reactome support is concentrated in three already-defined themes:",
            "",
            "- **Lipid / fatty-acid / eicosanoid metabolism**: 4 Reactome terms, including lipid metabolism, fatty-acid metabolism, and specialized pro-resolving mediator biosynthesis.",
            "- **Signaling / receptor / hormone response**: 1 Reactome term, Nuclear Receptor transcription pathway.",
            "- **Matrix remodeling / cell migration**: 1 Reactome term, Activation of Matrix Metalloproteinases.",
            "",
            "The other six themes have GO:BP support only in the current significant-term universe; this is not evidence that Reactome contradicts them.",
            "",
            "## Driver interpretation",
            "",
            "Driver contribution is descriptive repeated membership among retained/selected terms. It is not a gene-level test, causal attribution, or proof that a gene drives the enrichment.",
            "",
            "| Theme | GO core drivers | Reactome core drivers |",
            "|---|---|---|",
        ]
    )
    for row in family_rows:
        summary_lines.append(
            f"| {row['theme_label']} | {row['GO_core_drivers'] or '—'} | {row['Reactome_core_drivers'] or '—'} |"
        )
    summary_lines.extend(
        [
            "",
            "## Files",
            "",
            "- `go_reactome_aligned_terms.csv`: term-level alignment with source preserved.",
            "- `go_reactome_theme_alignment.csv`: nine-theme cross-resource summary.",
            "- `go_reactome_driver_contribution.csv`: theme-by-gene GO/Reactome counts and fractions.",
            "- `go_reactome_global_driver_contribution.csv`: global source-specific driver counts.",
            "",
            "## Excluded Reactome root",
            "",
            f"Excluded term ID: `{REACTOME_ROOT_ID}` (`REACTOME root term`). It was not assigned to any biological theme and is not included in driver counts.",
        ]
    )
    (OUT_DIR / "go_reactome_alignment_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_go_representatives": str(GO_REPRESENTATIVES),
        "input_go_representatives_sha256": sha256(GO_REPRESENTATIVES),
        "input_primary_106_terms": str(RAW_TERMS),
        "input_primary_106_terms_sha256": sha256(RAW_TERMS),
        "go_representatives": len(go_rows_raw),
        "reactome_rows": len(reactome_rows),
        "reactome_root_excluded": len(reactome_root_rows),
        "reactome_informative_terms": len(informative_reactome),
        "aligned_terms": len(aligned),
        "theme_count": len(THEMES),
        "reactome_theme_map": REACTOME_THEME_MAP,
        "statistical_boundary": "No enrichment P/FDR recalculation; GO and Reactome remain source-specific.",
        "outputs": [
            "go_reactome_aligned_terms.csv",
            "go_reactome_theme_alignment.csv",
            "go_reactome_driver_contribution.csv",
            "go_reactome_global_driver_contribution.csv",
            "go_reactome_alignment_summary.md",
        ],
    }
    (OUT_DIR / "go_reactome_alignment_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print("GO–Reactome alignment: PASS")
    print(f"GO representatives: {len(go_rows_raw)}")
    print(f"Reactome informative terms: {len(informative_reactome)}")
    print(f"Aligned themes: {len(THEMES)}")
    print(f"Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
