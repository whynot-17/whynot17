"""Independent epithelial replication of the MCOP–CRC PPAR/NR program.

This script analyzes the public GSE144735 processed annotation and raw UMI
matrix from GEO.  The primary comparison is core ``Tumor`` versus matched
``Normal`` epithelial cells from the same patient.  ``Border`` samples are
kept out of the primary contrast and are evaluated only as a sensitivity
analysis.

The inferential unit is patient-level pseudobulk.  Scores and paired tests
reuse the frozen Phase 2F-B implementation from
``mcop_phase2f_singlecell_validation.py``; cells are never treated as
independent observations.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mcop_phase2f_singlecell_validation import (
    ALL_GENES,
    NR_GENES,
    SCORE_COLUMNS,
    paired_donor_contrasts,
    score_pseudobulk,
)


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "work" / "mcop_phase2f_external" / "raw"
OUTPUT = ROOT / "outputs"
ANNOTATION = RAW / "GSE144735_annotation.txt.gz"
COUNTS = RAW / "GSE144735_raw_UMI_count_matrix.txt.gz"
GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_target_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Read only the nine frozen genes from GEO's gene-by-cell matrix."""

    target = set(ALL_GENES)
    found: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
        if len(header) < 2 or header[0] != "Index":
            raise ValueError("Unexpected GSE144735 matrix header; expected Index in column 1.")
        cell_ids = header[1:]
        for line in handle:
            fields = line.rstrip("\n\r").split("\t")
            if not fields or fields[0] not in target:
                continue
            if len(fields) != len(header):
                raise ValueError(f"Matrix row length mismatch for gene {fields[0]}.")
            found[fields[0]] = np.asarray(fields[1:], dtype=np.float64)

    missing = sorted(target - set(found))
    if missing:
        raise ValueError(f"Frozen target genes missing from GSE144735 matrix: {missing}")
    expression = pd.DataFrame({gene: found[gene] for gene in ALL_GENES})
    expression.insert(0, "Index", cell_ids)
    return expression, cell_ids


def build_pseudobulk(annotation: pd.DataFrame, expression: pd.DataFrame, include_border: bool) -> pd.DataFrame:
    annotation = annotation.copy()
    annotation["Index"] = annotation["Index"].astype(str)
    expression["Index"] = expression["Index"].astype(str)
    if set(expression["Index"]) != set(annotation["Index"]):
        missing = sorted(set(expression["Index"]) - set(annotation["Index"]))[:5]
        extra = sorted(set(annotation["Index"]) - set(expression["Index"]))[:5]
        raise ValueError(f"Annotation/matrix cell IDs do not match; missing={missing}, extra={extra}")

    meta = annotation.set_index("Index").loc[expression["Index"]].reset_index()
    values = expression[ALL_GENES].reset_index(drop=True)
    meta[ALL_GENES] = values
    meta["group"] = meta["Class"].map({"Tumor": "tumor", "Normal": "normal"})
    if include_border:
        meta.loc[meta["Class"].eq("Border"), "group"] = "tumor_border_included"
    meta["compartment"] = meta["Cell_type"].map({"Epithelial cells": "epithelial"})
    keep_groups = {"tumor", "normal"}
    if include_border:
        keep_groups.add("tumor_border_included")
    meta = meta[meta["group"].isin(keep_groups) & meta["compartment"].eq("epithelial")].copy()
    if meta.empty:
        raise ValueError("No epithelial Tumor/Normal observations remained after annotation gates.")

    rows: list[dict[str, object]] = []
    grouped = meta.groupby(["Patient", "group", "compartment"], observed=True)
    for (patient, group, compartment), subset in grouped:
        row: dict[str, object] = {
            "dataset_id": "GSE144735",
            "donor_key": f"GSE144735::{patient}",
            "donor_id": str(patient),
            "group": str(group),
            "compartment": str(compartment),
            "n_cells": int(len(subset)),
        }
        row.update({gene: float(subset[gene].sum()) for gene in ALL_GENES})
        rows.append(row)
    return pd.DataFrame(rows)


def write_report(annotation: pd.DataFrame, pseudobulk: pd.DataFrame, scored: pd.DataFrame, paired: pd.DataFrame, include_border: bool) -> None:
    lines = [
        "# MCOP–CRC Phase 2F-B：GSE144735 independent epithelial replication",
        "",
        "## Scope",
        "",
        f"- Source: [GSE144735 GEO record]({GEO_URL})",
        "- Primary contrast: core `Tumor` versus matched `Normal` epithelial cells.",
        "- `Border` samples are excluded from the primary contrast and are sensitivity-only.",
        "- Inferential unit: patient-level pseudobulk; no cell-level P values.",
        f"- Cells in annotation: **{len(annotation):,}**; epithelial cells analyzed: **{int((annotation['Cell_type'] == 'Epithelial cells').sum()):,}**.",
        f"- Patient-level pseudobulk rows: **{len(pseudobulk):,}**.",
        "",
        "## Primary paired result",
        "",
        "| score | paired patients | median tumor-minus-normal delta | Wilcoxon P |",
        "|---|---:|---:|---:|",
    ]
    for _, row in paired[paired.get("group_mode", pd.Series(index=paired.index, dtype=str)).eq("core_tumor_vs_normal")].iterrows():
        delta = "NA" if not np.isfinite(row["median_delta_tumor_minus_normal"]) else f"{row['median_delta_tumor_minus_normal']:.3f}"
        p = "NA" if not np.isfinite(row["p_value"]) else f"{row['p_value']:.3g}"
        lines.append(f"| {row['score']} | {int(row['paired_donors'])} | {delta} | {p} |")

    if include_border and "group_mode" in paired.columns:
        lines += [
            "",
            "## Border-included sensitivity",
            "",
            "This sensitivity adds Border epithelial cells to the tumor group; it is not substituted for the core Tumor-versus-Normal primary analysis.",
            "",
            "| score | paired patients | median tumor-plus-border-minus-normal delta | Wilcoxon P |",
            "|---|---:|---:|---:|",
        ]
        for _, row in paired[paired["group_mode"].eq("tumor_plus_border_vs_normal")].iterrows():
            delta = "NA" if not np.isfinite(row["median_delta_tumor_minus_normal"]) else f"{row['median_delta_tumor_minus_normal']:.3f}"
            p = "NA" if not np.isfinite(row["p_value"]) else f"{row['p_value']:.3g}"
            lines.append(f"| {row['score']} | {int(row['paired_donors'])} | {delta} | {p} |")

    lines += [
        "",
        "## Interpretation",
        "",
        "A concordant negative epithelial PPAR/NR delta would provide independent disease-state replication of the Census/TCGA direction. It does not establish DINP/MCOP causality or exposure mediation.",
        "",
        "## Reproducibility",
        "",
        f"- Annotation SHA256: `{sha256(ANNOTATION)}`",
        f"- Matrix SHA256: `{sha256(COUNTS)}`",
        f"- Frozen NR genes: `{','.join(NR_GENES)}`",
        f"- Run UTC: `{datetime.now(timezone.utc).isoformat()}`",
        "",
        "## Output files",
        "",
        "- `mcop_phase2f_external_gse144735_pseudobulk.csv`",
        "- `mcop_phase2f_external_gse144735_scores.csv`",
        "- `mcop_phase2f_external_gse144735_paired_contrasts.csv`",
        "- `mcop_phase2f_external_gse144735_border_sensitivity_pseudobulk.csv`",
        "- `mcop_phase2f_external_gse144735_border_sensitivity_scores.csv`",
        "- `mcop_phase2f_external_gse144735_report.md`",
    ]
    (OUTPUT / "mcop_phase2f_external_gse144735_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-border", action="store_true", help="Add Border epithelial cells to the tumor group as sensitivity analysis.")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if not ANNOTATION.exists() or not COUNTS.exists():
        raise SystemExit("GSE144735 raw files are missing; download the public GEO supplementary files first.")

    annotation = pd.read_csv(ANNOTATION, sep="\t")
    expression, cell_ids = read_target_matrix(COUNTS)
    pseudobulk = build_pseudobulk(annotation, expression, include_border=False)
    scored = score_pseudobulk(pseudobulk)
    paired_primary = paired_donor_contrasts(scored)
    paired_primary["group_mode"] = "core_tumor_vs_normal"
    paired_frames = [paired_primary]
    if args.include_border:
        pseudobulk_border = build_pseudobulk(annotation, expression, include_border=True)
        scored_border = score_pseudobulk(pseudobulk_border)
        paired_border = paired_donor_contrasts(scored_border)
        paired_border["group_mode"] = "tumor_plus_border_vs_normal"
        paired_frames.append(paired_border)
        pseudobulk_border.to_csv(OUTPUT / "mcop_phase2f_external_gse144735_border_sensitivity_pseudobulk.csv", index=False)
        scored_border.to_csv(OUTPUT / "mcop_phase2f_external_gse144735_border_sensitivity_scores.csv", index=False)
    paired = pd.concat(paired_frames, ignore_index=True)

    pseudobulk.to_csv(OUTPUT / "mcop_phase2f_external_gse144735_pseudobulk.csv", index=False)
    scored.to_csv(OUTPUT / "mcop_phase2f_external_gse144735_scores.csv", index=False)
    paired.to_csv(OUTPUT / "mcop_phase2f_external_gse144735_paired_contrasts.csv", index=False)
    write_report(annotation, pseudobulk, scored, paired, include_border=args.include_border)
    manifest = {
        "analysis": "MCOP-CRC Phase 2F-B independent epithelial replication",
        "dataset": "GSE144735",
        "source_url": GEO_URL,
        "primary_contrast": "Tumor versus Normal epithelial cells; Border excluded",
        "unit_of_analysis": "patient-level pseudobulk",
        "n_annotation_cells": int(len(annotation)),
        "n_matrix_cells": int(len(cell_ids)),
        "n_pseudobulk_rows": int(len(pseudobulk)),
        "n_patients": int(annotation["Patient"].nunique()),
        "nr_genes": NR_GENES,
        "all_genes": ALL_GENES,
        "annotation_sha256": sha256(ANNOTATION),
        "matrix_sha256": sha256(COUNTS),
        "include_border_sensitivity_requested": bool(args.include_border),
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUTPUT / "mcop_phase2f_external_gse144735_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
