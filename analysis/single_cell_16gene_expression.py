#!/usr/bin/env python
"""Descriptive expression summary for the 16-gene model in two CRC scRNA-seq cohorts.

Reads GEO processed raw UMI matrices one gene at a time, keeping only target and
marker genes in memory. Expressions are normalized per cell to 10,000 UMIs and
log1p transformed. This is descriptive; cells are not treated as independent
biological replicates for hypothesis testing.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import html
import json
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "work" / "single_cell_GSE132465_GSE200997"
OUT = ROOT / "results" / "2026-10-11" / "single_cell_16genes_GSE132465_GSE200997"
GENES = [
    "ABHD6", "ACAT1", "BMP2", "CPT2", "ELOVL6", "FABP4", "FHL1", "INHBA",
    "PCK1", "PID1", "RETSAT", "SCD", "SCP2", "TIMP1", "TPM1", "UBE2C",
]

# Canonical CRC microenvironment marker signatures for coarse cell-type inference
# in GSE200997, whose per-cell annotation does not contain cell-type labels.
MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CEACAM5"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "C1R"],
    "Endothelial": ["PECAM1", "CLDN5", "VWF", "ENG", "KDR"],
    "Myeloid": ["LYZ", "MARCO", "CD68", "FCGR3A", "S100A8", "S100A9", "CD14"],
    "T_cell": ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B"],
    "B_cell": ["CD79A", "MS4A1", "CD19", "CD37", "CD74"],
    "Plasma": ["MZB1", "JCHAIN", "IGKC", "XBP1", "IGHG1"],
    "NK": ["NKG7", "GNLY", "PRF1", "GZMB"],
    "Mast": ["KIT", "CPA3", "TPSAB1", "TPSB2"],
}

COHORTS = [
    {
        "name": "GSE132465",
        "matrix": "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
        "annotation": "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
        "delimiter": "\t",
        "cell_key": "Index",
        "sample_key": "Sample",
        "condition_key": "Class",
        "original_cell_type_key": "Cell_type",
    },
    {
        "name": "GSE200997",
        "matrix": "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz",
        "annotation": "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz",
        "delimiter": ",",
        "cell_key": "Unnamed: 0",
        "sample_key": "samples",
        "condition_key": "Condition",
        "original_cell_type_key": None,
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_annotation(path: Path, delimiter: str, cell_key: str) -> pd.DataFrame:
    ann = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False)
    if cell_key not in ann.columns:
        raise ValueError(f"Cell key {cell_key!r} not found in {path.name}: {ann.columns.tolist()}")
    ann = ann.set_index(cell_key, drop=False)
    if ann.index.duplicated().any():
        dup = ann.index[ann.index.duplicated()].tolist()[:5]
        raise ValueError(f"Duplicate annotation cell IDs in {path.name}: {dup}")
    return ann


def parse_matrix(path: Path, delimiter: str, required_genes: set[str]):
    """Stream a genes-by-cells text matrix and return cell libraries and selected rows."""
    delim = delimiter.encode("ascii")
    started = time.time()
    selected: dict[str, np.ndarray] = {}
    duplicate_rows = Counter()
    with gzip.open(path, "rb") as f:
        header = f.readline().rstrip(b"\r\n").split(delim)
        cells = [x.decode("utf-8").strip().strip('"') for x in header[1:]]
        if len(cells) != len(set(cells)):
            raise ValueError(f"Duplicate cell barcodes in matrix header: {path.name}")
        libraries = np.zeros(len(cells), dtype=np.float64)
        n_features = 0
        for line_no, raw in enumerate(f, start=2):
            raw = raw.rstrip(b"\r\n")
            if not raw:
                continue
            pieces = raw.split(delim, 1)
            if len(pieces) != 2:
                raise ValueError(f"Malformed row {line_no} in {path.name}")
            gene = pieces[0].decode("utf-8").strip().strip('"')
            counts = np.fromstring(pieces[1].decode("ascii"), sep=delimiter, dtype=np.float64)
            if counts.size != len(cells):
                raise ValueError(
                    f"Row {line_no} ({gene}) has {counts.size} counts, expected {len(cells)}"
                )
            libraries += counts
            n_features += 1
            if gene in required_genes:
                if gene in selected:
                    selected[gene] += counts
                    duplicate_rows[gene] += 1
                else:
                    selected[gene] = counts.copy()
    if not np.all(np.isfinite(libraries)) or np.any(libraries <= 0):
        raise ValueError(f"Invalid per-cell library size in {path.name}")
    return cells, libraries, selected, n_features, dict(duplicate_rows), time.time() - started


def infer_marker_types(log_expr: dict[str, np.ndarray], cells: int):
    scores: dict[str, np.ndarray] = {}
    detected: dict[str, np.ndarray] = {}
    available: dict[str, list[str]] = {}
    for cell_type, genes in MARKERS.items():
        present = [g for g in genes if g in log_expr]
        available[cell_type] = present
        if present:
            scores[cell_type] = np.mean(np.vstack([log_expr[g] for g in present]), axis=0)
            detected[cell_type] = np.sum(
                np.vstack([log_expr[g] > 0 for g in present]), axis=0
            )
        else:
            scores[cell_type] = np.zeros(cells)
            detected[cell_type] = np.zeros(cells, dtype=int)
    names = list(scores)
    score_mat = np.vstack([scores[k] for k in names])
    order = np.argsort(score_mat, axis=0)
    best_idx = order[-1, :]
    best_score = score_mat[best_idx, np.arange(cells)]
    second_score = score_mat[order[-2, :], np.arange(cells)]
    best_name = np.array(names, dtype=object)[best_idx]
    best_n = np.array([detected[n] for n in names])[best_idx, np.arange(cells)]
    margin = best_score - second_score
    assigned = np.where((best_score > 0) & (best_n >= 2), best_name, "Unassigned")
    confident = (assigned != "Unassigned") & (margin >= 0.05)
    return assigned, best_score, margin, best_n, confident, available


def group_expression(cohort, ann, cells, library, selected, marker_selected, feature_count,
                     duplicate_rows, elapsed):
    scale = 10000.0 / library
    log_expr = {g: np.log1p(v * scale) for g, v in selected.items()}
    missing_targets = [g for g in GENES if g not in log_expr]
    if cohort["name"] == "GSE132465":
        labels = ann[cohort["original_cell_type_key"]].astype(str).to_numpy()
        label_method = "GEO original cell_type annotation"
        confidence = np.full(len(cells), "author_annotated", dtype=object)
        marker_score = np.full(len(cells), np.nan)
        marker_margin = np.full(len(cells), np.nan)
        positive_markers = np.full(len(cells), np.nan)
    else:
        marker_log = {g: np.log1p(v * scale) for g, v in marker_selected.items()}
        labels, marker_score, marker_margin, positive_markers, high_conf, _ = infer_marker_types(
            marker_log, len(cells)
        )
        confidence = np.where(high_conf, "higher_confidence", "lower_confidence")
        label_method = "coarse marker-score inference; top mean log1p(CP10K), at least 2 detected markers"

    condition = ann[cohort["condition_key"]].astype(str).to_numpy()
    sample = ann[cohort["sample_key"]].astype(str).to_numpy()
    if len(labels) != len(cells):
        raise ValueError(f"Annotation/cell count mismatch for {cohort['name']}")

    cell_meta = pd.DataFrame({
        "cell_id": cells,
        "sample": sample,
        "condition": condition,
        "cell_type": labels,
        "cell_type_label_method": label_method,
        "label_confidence": confidence,
        "marker_top_score_mean_log1p_CP10K": marker_score,
        "marker_score_margin": marker_margin,
        "positive_markers_in_top_signature": positive_markers,
        "library_umi": library,
    })
    # Per-cell-type counts and marker-inference audit.
    type_rows = []
    for (cond, ctype), sub in cell_meta.groupby(["condition", "cell_type"], dropna=False):
        type_rows.append({
            "cohort": cohort["name"], "condition": cond, "cell_type": ctype,
            "n_cells": len(sub), "n_samples": sub["sample"].nunique(),
            "n_higher_confidence_cells": int((sub["label_confidence"] == "higher_confidence").sum()),
            "median_library_umi": float(sub["library_umi"].median()),
            "cell_type_label_method": label_method,
        })
    type_counts = pd.DataFrame(type_rows)

    # Cell-level descriptive expression summaries by cohort, condition, and cell type.
    expr_rows = []
    for gene in GENES:
        if gene not in selected:
            continue
        values = log_expr[gene]
        detected_values = selected[gene] > 0
        tmp = cell_meta[["condition", "cell_type"]].copy()
        tmp["expr"] = values
        tmp["detected"] = detected_values
        for (cond, ctype), sub in tmp.groupby(["condition", "cell_type"], dropna=False):
            expr_rows.append({
                "cohort": cohort["name"], "gene": gene, "condition": cond,
                "cell_type": ctype, "n_cells": len(sub),
                "pct_cells_detected": float(sub["detected"].mean() * 100),
                "mean_log1p_CP10K": float(sub["expr"].mean()),
                "median_log1p_CP10K": float(sub["expr"].median()),
                "mean_log1p_CP10K_detected_cells": float(sub.loc[sub["detected"], "expr"].mean())
                    if sub["detected"].any() else 0.0,
                "cell_type_label_method": label_method,
            })
    expr_summary = pd.DataFrame(expr_rows)

    # Sample-level summaries preserve the biological sample as the aggregation unit.
    sample_rows = []
    for gene in GENES:
        if gene not in selected:
            continue
        tmp = cell_meta[["sample", "condition", "cell_type"]].copy()
        tmp["expr"] = log_expr[gene]
        tmp["detected"] = selected[gene] > 0
        for (sid, cond, ctype), sub in tmp.groupby(["sample", "condition", "cell_type"], dropna=False):
            sample_rows.append({
                "cohort": cohort["name"], "sample": sid, "condition": cond,
                "cell_type": ctype, "gene": gene, "n_cells": len(sub),
                "pct_cells_detected": float(sub["detected"].mean() * 100),
                "mean_log1p_CP10K": float(sub["expr"].mean()),
                "mean_log1p_CP10K_detected_cells": float(sub.loc[sub["detected"], "expr"].mean())
                    if sub["detected"].any() else 0.0,
            })
    sample_summary = pd.DataFrame(sample_rows)

    # Marker coverage is explicit so users can see whether a cell type signature
    # was inferred from all canonical markers or only a subset represented in GEO.
    marker_coverage = []
    for cell_type, genes in MARKERS.items():
        marker_coverage.append({
            "cohort": cohort["name"], "putative_cell_type": cell_type,
            "canonical_marker_genes": ",".join(genes),
            "genes_present_in_matrix": ",".join([g for g in genes if g in marker_selected]),
            "n_present": sum(g in marker_selected for g in genes),
            "n_canonical": len(genes),
        })
    coverage = pd.DataFrame(marker_coverage)

    audit = {
        "cohort": cohort["name"], "matrix": cohort["matrix"],
        "annotation": cohort["annotation"], "n_cells": len(cells),
        "n_features_in_matrix": feature_count, "n_target_genes_present": len(GENES) - len(missing_targets),
        "missing_target_genes": missing_targets,
        "n_marker_genes_present": len(marker_selected),
        "duplicate_gene_rows_aggregated": duplicate_rows,
        "cell_id_order_checked_and_aligned_by_key": True,
        "all_annotation_barcodes_match_matrix": True,
        "cell_types": label_method,
        "normalization": "per-cell counts per 10,000 total UMI across all matrix rows; log1p",
        "n_samples": int(cell_meta["sample"].nunique()),
        "conditions": dict(Counter(condition)),
        "median_library_umi": float(np.median(library)),
        "matrix_parse_seconds": round(elapsed, 2),
    }
    return type_counts, expr_summary, sample_summary, coverage, cell_meta, audit


def write_dotplot(expr: pd.DataFrame, cohort: str, path: Path) -> None:
    if cohort == "GSE132465":
        cell_types = ["Epithelial cells", "Stromal cells", "Myeloids", "T cells", "B cells", "Mast cells"]
        labels = {"Epithelial cells": "Epithelial", "Stromal cells": "Stromal", "Myeloids": "Myeloid",
                  "T cells": "T cells", "B cells": "B cells", "Mast cells": "Mast"}
    else:
        cell_types = ["Epithelial", "Fibroblast", "Myeloid", "T_cell", "B_cell", "Plasma", "NK",
                      "Endothelial", "Mast", "Unassigned"]
        labels = {x: x.replace("_", " ") for x in cell_types}
    conditions = ["Normal", "Tumor"]
    groups = [(ct, cond) for ct in cell_types for cond in conditions
              if ((expr["cohort"] == cohort) & (expr["cell_type"] == ct) &
                  (expr["condition"] == cond)).any()]
    if not groups:
        return
    w, left, right = 1180, 115, 35
    top, row_h, bottom = 92, 34, 180
    col_w = min(52, (w - left - right) / len(groups))
    h = top + row_h * len(GENES) + bottom
    lookup = {(r.gene, r.cell_type, r.condition): r for r in
              expr[expr["cohort"] == cohort].itertuples(index=False)}
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="{left}" y="34" font-family="Arial,sans-serif" font-size="22" font-weight="bold">{html.escape(cohort)}: 16-gene single-cell expression</text>',
           f'<text x="{left}" y="58" font-family="Arial,sans-serif" font-size="13" fill="#444">Dot size: cells with detected UMI (%)  |  Color: mean log1p(CP10K), capped at 3</text>']
    for yi, gene in enumerate(GENES):
        y = top + yi * row_h
        if yi % 2 == 0:
            svg.append(f'<rect x="{left-8}" y="{y-18}" width="{w-left-right+16}" height="{row_h}" fill="#f7f9fb"/>')
        svg.append(f'<text x="{left-14}" y="{y+5}" text-anchor="end" font-family="Arial,sans-serif" font-size="13">{html.escape(gene)}</text>')
    for xi, (ct, cond) in enumerate(groups):
        x = left + col_w * (xi + 0.5)
        for yi, gene in enumerate(GENES):
            rec = lookup.get((gene, ct, cond))
            if rec is None:
                continue
            pct = float(rec.pct_cells_detected)
            value = max(0.0, min(3.0, float(rec.mean_log1p_CP10K))) / 3.0
            start = (247, 251, 255)
            end = (8, 48, 107)
            rgb = tuple(round(start[i] + value * (end[i] - start[i])) for i in range(3))
            color = "#%02x%02x%02x" % rgb
            radius = 1.2 + math.sqrt(max(0.0, pct) / 100.0) * 7.4
            y = top + yi * row_h
            svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.2f}" fill="{color}" stroke="#667" stroke-width="0.4"/>')
        label_y = top + row_h * len(GENES) + 18
        short = labels[ct]
        svg.append(f'<text x="{x:.1f}" y="{label_y}" transform="rotate(-55 {x:.1f} {label_y})" text-anchor="end" font-family="Arial,sans-serif" font-size="11">{html.escape(short)}</text>')
        cond_y = h - 19
        svg.append(f'<text x="{x:.1f}" y="{cond_y}" text-anchor="middle" font-family="Arial,sans-serif" font-size="10" fill="#555">{cond[0]}</text>')
    # Compact size legend.
    legend_y = h - 50
    legend_x = left + 3
    svg.append(f'<text x="{legend_x}" y="{legend_y-10}" font-family="Arial,sans-serif" font-size="10" fill="#555">Detection</text>')
    for i, pct in enumerate([10, 30, 60, 90]):
        x = legend_x + 68 + i * 48
        r = 1.2 + math.sqrt(pct / 100.0) * 7.4
        svg.append(f'<circle cx="{x}" cy="{legend_y-13}" r="{r:.2f}" fill="#6b8fb3"/>')
        svg.append(f'<text x="{x}" y="{legend_y+5}" text-anchor="middle" font-family="Arial,sans-serif" font-size="9" fill="#555">{pct}%</text>')
    svg.append('</svg>')
    path.write_text("\n".join(svg), encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    all_types, all_expr, all_samples, all_coverage, all_cells, audits = [], [], [], [], [], []
    files = {}
    for c in COHORTS:
        matrix = INPUT / c["matrix"]
        annotation_path = INPUT / c["annotation"]
        if not matrix.is_file() or not annotation_path.is_file():
            raise FileNotFoundError(f"Missing local input: {matrix} or {annotation_path}")
        ann = read_annotation(annotation_path, c["delimiter"], c["cell_key"])
        # Use only the target plus marker rows when retaining expression vectors.
        required = set(GENES)
        for genes in MARKERS.values():
            required.update(genes)
        cache_path = INPUT / (c["name"] + "_target_marker_library_cache.npz")
        used_cache = cache_path.is_file()
        if used_cache:
            cached = np.load(cache_path, allow_pickle=False)
            cells = [x.strip().strip('"') for x in cached["cells"].astype(str).tolist()]
            library = cached["libraries"]
            n_features = int(cached["n_features"])
            cached_genes = cached["genes"].astype(str).tolist()
            selected = {g: cached["counts"][i] for i, g in enumerate(cached_genes)}
            duplicate_rows = json.loads(str(cached["duplicate_rows"]))
            elapsed = 0.0
            print(f"{c['name']}: using cached streamed matrix vectors", flush=True)
        else:
            cells, library, selected, n_features, duplicate_rows, elapsed = parse_matrix(
                matrix, c["delimiter"], required
            )
            cache_genes = list(selected)
            np.savez_compressed(
                cache_path, cells=np.asarray(cells), libraries=library,
                genes=np.asarray(cache_genes),
                counts=np.vstack([selected[g] for g in cache_genes]),
                n_features=np.asarray(n_features),
                duplicate_rows=np.asarray(json.dumps(duplicate_rows)),
            )
            print(f"{c['name']}: cached {len(cache_genes)} selected gene vectors in work/", flush=True)
        matrix_ids = set(cells)
        ann_ids = set(ann.index.astype(str))
        missing_ann = matrix_ids - ann_ids
        missing_matrix = ann_ids - matrix_ids
        if missing_ann or missing_matrix or len(ann) != len(cells):
            raise ValueError(
                f"Cell ID mismatch {c['name']}: matrix_only={len(missing_ann)}, "
                f"annotation_only={len(missing_matrix)}, matrix={len(cells)}, annotation={len(ann)}"
            )
        # Explicitly align annotation rows to matrix columns by barcode.
        ann = ann.loc[cells]
        target = {g: selected[g] for g in GENES if g in selected}
        markers = {g: selected[g] for g in set().union(*[set(v) for v in MARKERS.values()]) if g in selected}
        outputs = group_expression(c, ann, cells, library, target, markers, n_features,
                                   duplicate_rows, elapsed)
        t, e, s, cov, cm, audit = outputs
        audit["work_vector_cache_used"] = used_cache
        all_types.append(t); all_expr.append(e); all_samples.append(s)
        all_coverage.append(cov); all_cells.append(cm); audits.append(audit)
        files[c["name"]] = {
            "matrix": {"file": c["matrix"], "bytes": matrix.stat().st_size, "sha256": sha256(matrix)},
            "annotation": {"file": c["annotation"], "bytes": annotation_path.stat().st_size,
                           "sha256": sha256(annotation_path)},
        }
        print(f"{c['name']}: {len(cells):,} cells; {len(target)}/{len(GENES)} target genes; "
              f"streamed {n_features:,} rows in {elapsed:.1f}s", flush=True)

    types_df = pd.concat(all_types, ignore_index=True)
    expr_df = pd.concat(all_expr, ignore_index=True)
    samples_df = pd.concat(all_samples, ignore_index=True)
    types_df.to_csv(OUT / "cell_type_counts.tsv", sep="\t", index=False)
    expr_df.to_csv(OUT / "expression_by_cell_type.tsv", sep="\t", index=False)
    samples_df.to_csv(OUT / "expression_by_sample_cell_type.tsv", sep="\t", index=False)
    sample_level = samples_df.groupby(
        ["cohort", "gene", "cell_type", "condition"], as_index=False
    ).agg(
        n_samples=("sample", "nunique"),
        mean_of_sample_means_log1p_CP10K=("mean_log1p_CP10K", "mean"),
        median_of_sample_means_log1p_CP10K=("mean_log1p_CP10K", "median"),
        mean_pct_detected_across_samples=("pct_cells_detected", "mean"),
        median_cells_per_sample=("n_cells", "median"),
    )
    wide = sample_level.pivot(index=["cohort", "gene", "cell_type"], columns="condition")
    delta_rows = []
    for idx, row in wide.iterrows():
        record = dict(zip(["cohort", "gene", "cell_type"], idx))
        for metric in ["n_samples", "mean_of_sample_means_log1p_CP10K",
                       "median_of_sample_means_log1p_CP10K", "mean_pct_detected_across_samples",
                       "median_cells_per_sample"]:
            for cond in ["Normal", "Tumor"]:
                key = (metric, cond)
                record[f"{cond.lower()}_{metric}"] = row[key] if key in row.index else np.nan
        record["tumor_minus_normal_mean_of_sample_means_log1p_CP10K"] = (
            record["tumor_mean_of_sample_means_log1p_CP10K"] - record["normal_mean_of_sample_means_log1p_CP10K"]
            if pd.notna(record["tumor_mean_of_sample_means_log1p_CP10K"]) and
               pd.notna(record["normal_mean_of_sample_means_log1p_CP10K"]) else np.nan
        )
        record["tumor_minus_normal_mean_detection_percentage_points"] = (
            record["tumor_mean_pct_detected_across_samples"] - record["normal_mean_pct_detected_across_samples"]
            if pd.notna(record["tumor_mean_pct_detected_across_samples"]) and
               pd.notna(record["normal_mean_pct_detected_across_samples"]) else np.nan
        )
        delta_rows.append(record)
    pd.DataFrame(delta_rows).to_csv(OUT / "tumor_normal_sample_level_descriptive.tsv", sep="\t", index=False)
    write_dotplot(expr_df, "GSE132465", OUT / "dotplot_GSE132465.svg")
    write_dotplot(expr_df, "GSE200997", OUT / "dotplot_GSE200997.svg")
    pd.concat(all_coverage, ignore_index=True).to_csv(OUT / "GSE200997_marker_coverage.tsv", sep="\t", index=False)
    pd.concat(all_cells, ignore_index=True).to_csv(OUT / "cell_annotation_and_label_scores.tsv.gz", sep="\t", index=False)
    metadata = {
        "analysis_date": "2026-10-11", "genes": GENES, "input_files": files,
        "cohort_audits": audits,
        "GSE200997_inference": {
            "method": "For each cell, assign the highest average log1p(CP10K) among canonical marker panels; require at least two detected markers. Higher-confidence is a score margin of at least 0.05 over the runner-up.",
            "warning": "These are coarse putative labels inferred from marker signatures, not author-provided cell labels. Interpret ambiguous cells cautiously.",
        },
        "statistical_scope": "Descriptive expression summaries only. Cells are nested within samples; no cell-level p-values or differential-expression claims are reported.",
    }
    (OUT / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("Wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
