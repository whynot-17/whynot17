#!/usr/bin/env python3
"""Localize the fresh 41-gene DINP--CRC program in a CRC single-cell object.

The input gene set is the fresh ``query_41_genes.csv`` reconstructed from the
Gelsemium elegans--CRC GeneCards intersection and the DINP multi-source set.
This analysis does not reuse the legacy 81-gene program.

The source object is already processed expression data.  For localization, a
gene-wise z score is computed across eligible primary CRC cells and averaged
to an observed-gene program score.  Inference is donor-level and paired where
the same donor has tumor and normal cells.  A gene absent from the source
feature table is excluded from the score, reported explicitly, and never
treated as zero expression.

Boundary: source-labeled tumor epithelial cells are not called definitively
malignant without an independent CNV/malignancy validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
GENE_FILE = ROOT / "analysis" / "dinp_crc_gelsemium_reference" / "outputs" / "query_41_genes.csv"
DEFAULT_H5AD = Path(
    r"D:\mcop and CRC\cellxgene_census\2025-11-08\16023185-de21-4c0d-a9c8-73abdd52d142.h5ad"
)
DEFAULT_OUT = (
    ROOT
    / "analysis"
    / "dinp_crc_gelsemium_reference"
    / "outputs"
    / "singlecell_41_gene_localization"
)
TUMOR = "colon adenocarcinoma"
NORMAL = "normal"
COMPARTMENTS = ("epithelial", "myeloid", "fibroblast", "endothelial")
CHUNK_SIZE = 10_000


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def resolve_cross_platform_path(path: Path) -> Path:
    """Resolve Windows drive paths when the script is launched inside WSL."""
    raw = str(path).strip()
    if re.match(r"^[A-Za-z]:[\\/]", raw):
        drive = raw[0].lower()
        suffix = raw[2:].replace("\\", "/")
        mounted = Path(f"/mnt/{drive}{suffix}")
        if mounted.exists():
            return mounted.resolve()
    return path.expanduser().resolve()


def dense_matrix(value: object) -> np.ndarray:
    if sparse.issparse(value):
        return value.toarray().astype(np.float32, copy=False)
    if hasattr(value, "toarray"):
        return value.toarray().astype(np.float32, copy=False)
    return np.asarray(value, dtype=np.float32)


def read_selected_csr_columns(
    h5ad_path: Path,
    row_count: int,
    column_count: int,
    target_indices: list[int],
    chunk_size: int,
) -> np.ndarray:
    """Read only target columns from a CSR-backed H5AD in row chunks.

    Some older anndata/scipy combinations cannot slice a backed CSRDataset
    after the object has been opened from a mounted Windows volume.  Direct
    chunked CSR reconstruction is numerically equivalent and avoids loading
    the full 370k x 38k matrix into memory.
    """
    result = np.empty((row_count, len(target_indices)), dtype=np.float32)
    with h5py.File(h5ad_path, "r") as handle:
        x_group = handle["X"]
        encoding = str(x_group.attrs.get("encoding-type", ""))
        if encoding != "csr_matrix":
            raise ValueError(f"Expected CSR H5AD X, found {encoding!r}")
        data = x_group["data"]
        indices = x_group["indices"]
        indptr = x_group["indptr"]
        for start in range(0, row_count, chunk_size):
            stop = min(start + chunk_size, row_count)
            ptr = np.asarray(indptr[start : stop + 1], dtype=np.int64)
            nnz_start, nnz_stop = int(ptr[0]), int(ptr[-1])
            block = sparse.csr_matrix(
                (
                    np.asarray(data[nnz_start:nnz_stop], dtype=np.float32),
                    np.asarray(indices[nnz_start:nnz_stop], dtype=np.int64),
                    ptr - nnz_start,
                ),
                shape=(stop - start, column_count),
            )
            result[start:stop] = block[:, target_indices].toarray().astype(
                np.float32, copy=False
            )
            if stop == row_count or (start // chunk_size) % 10 == 0:
                print(f"  expression rows {stop:,}/{row_count:,}", flush=True)
    return result


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(values)
    if valid.any():
        result.loc[valid] = multipletests(
            values.loc[valid].to_numpy(dtype=float), method="fdr_bh"
        )[1]
    return result


def classify_compartment(cell_type: object) -> str | None:
    value = str(cell_type).lower()
    if re.search(
        r"epithelial|colonocyte|enterocyte|goblet|paneth|enteroendocrine|"
        r"tuft|best4|transit amplifying|stem cell",
        value,
    ):
        return "epithelial"
    if re.search(r"macrophage|monocyte|dendritic|granulocyte|myeloid", value):
        return "myeloid"
    if re.search(r"fibroblast|myofibroblast", value):
        return "fibroblast"
    if "endothelial" in value:
        return "endothelial"
    return None


def classify_subtype(row: pd.Series) -> str | None:
    compartment = row["compartment"]
    disease = row["group"]
    cell_type = str(row["cell_type"]).lower()
    cluster_full = str(row.get("ClusterFull", ""))
    if compartment == "epithelial":
        if disease == TUMOR:
            return (
                "source_tumor_labeled_epithelial"
                if cluster_full.lower().startswith("tumor")
                else "other_tumor_epithelial"
            )
        if disease == NORMAL:
            return "normal_epithelial"
    if compartment == "myeloid":
        if "macrophage" in cell_type or "macrophage" in cluster_full.lower():
            return "macrophage"
        if "monocyte" in cell_type or "monocyte" in cluster_full.lower():
            return "monocyte"
        if "dendritic" in cell_type or re.search(r"\bcM0[3-9]|\bmregDC", cluster_full):
            return "dendritic"
        if "granulocyte" in cell_type or "granulocyte" in cluster_full.lower():
            return "granulocyte"
    return None


def summarize_scores(frame: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    grouped = (
        frame.groupby(group_cols, observed=True)["program_score"]
        .agg(
            n_cells="size",
            mean_program_score="mean",
            median_program_score="median",
            sd_cell_program_score="std",
        )
        .reset_index()
    )
    return grouped


def paired_contrasts(donor_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    specs = [
        ("compartment", "epithelial", "epithelial tumor vs normal"),
        ("compartment", "myeloid", "myeloid tumor vs normal"),
        ("compartment", "fibroblast", "fibroblast tumor vs normal"),
        ("compartment", "endothelial", "endothelial tumor vs normal"),
        ("subtype", "source_tumor_labeled_epithelial", "source tumor-labeled epithelial vs normal epithelial"),
        ("subtype", "other_tumor_epithelial", "other tumor epithelial vs normal epithelial"),
        ("subtype", "macrophage", "macrophage tumor vs normal"),
        ("subtype", "monocyte", "monocyte tumor vs normal"),
        ("subtype", "dendritic", "dendritic tumor vs normal"),
        ("subtype", "granulocyte", "granulocyte tumor vs normal"),
    ]
    summary_rows: list[dict] = []
    delta_rows: list[dict] = []
    for family, target, label in specs:
        if family == "compartment":
            tumor = donor_scores.loc[
                donor_scores["compartment"].eq(target)
                & donor_scores["group"].eq(TUMOR)
            ]
            normal = donor_scores.loc[
                donor_scores["compartment"].eq(target)
                & donor_scores["group"].eq(NORMAL)
            ]
        else:
            normal_target = "normal_epithelial" if "epithelial" in target else target
            tumor = donor_scores.loc[
                donor_scores["subtype"].eq(target)
                & donor_scores["group"].eq(TUMOR)
            ]
            normal = donor_scores.loc[
                donor_scores["subtype"].eq(normal_target)
                & donor_scores["group"].eq(NORMAL)
            ]
        # A donor can contribute multiple fine cell types within one broad
        # compartment.  Collapse those rows to one donor-level mean before
        # pairing, rather than allowing duplicate donor index labels.
        tumor_scores = tumor.groupby("donor_id")["program_score"].mean().rename("tumor")
        normal_scores = normal.groupby("donor_id")["program_score"].mean().rename("normal")
        paired = pd.concat([tumor_scores, normal_scores], axis=1).dropna()
        deltas = (paired["tumor"] - paired["normal"]).to_numpy(dtype=float)
        for donor_id, row in paired.iterrows():
            delta_rows.append(
                {
                    "contrast": label,
                    "family": family,
                    "target": target,
                    "donor_id": str(donor_id),
                    "tumor_score": float(row["tumor"]),
                    "normal_score": float(row["normal"]),
                    "tumor_minus_normal": float(row["tumor"] - row["normal"]),
                }
            )
        n = len(deltas)
        mean_delta = float(np.mean(deltas)) if n else np.nan
        median_delta = float(np.median(deltas)) if n else np.nan
        sd_delta = float(np.std(deltas, ddof=1)) if n > 1 else np.nan
        if n > 1 and np.isfinite(sd_delta) and sd_delta > 0:
            t_stat, t_p = stats.ttest_1samp(deltas, 0.0)
            critical = stats.t.ppf(0.975, n - 1)
            half_width = critical * sd_delta / np.sqrt(n)
            ci_low, ci_high = mean_delta - half_width, mean_delta + half_width
        else:
            t_stat, t_p, ci_low, ci_high = np.nan, np.nan, np.nan, np.nan
        if n >= 5 and np.any(np.abs(deltas) > 0):
            try:
                w_stat, w_p = stats.wilcoxon(deltas, alternative="two-sided", method="auto")
            except ValueError:
                w_stat, w_p = np.nan, np.nan
        else:
            w_stat, w_p = np.nan, np.nan
        summary_rows.append(
            {
                "contrast": label,
                "family": family,
                "target": target,
                "n_paired_donors": n,
                "mean_delta_tumor_minus_normal": mean_delta,
                "median_delta_tumor_minus_normal": median_delta,
                "sd_delta": sd_delta,
                "mean_delta_95ci_low": ci_low,
                "mean_delta_95ci_high": ci_high,
                "paired_t_statistic": float(t_stat) if np.isfinite(t_stat) else np.nan,
                "paired_t_p": float(t_p) if np.isfinite(t_p) else np.nan,
                "paired_wilcoxon_statistic": float(w_stat) if np.isfinite(w_stat) else np.nan,
                "paired_wilcoxon_p": float(w_p) if np.isfinite(w_p) else np.nan,
                "direction": "up" if mean_delta > 0 else "down" if mean_delta < 0 else "flat",
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary["paired_t_BH_FDR"] = bh(summary["paired_t_p"])
    summary["paired_wilcoxon_BH_FDR"] = bh(summary["paired_wilcoxon_p"])
    return summary, pd.DataFrame(delta_rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--genes", type=Path, default=GENE_FILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    parser.add_argument("--require-all-genes", action="store_true")
    args = parser.parse_args()

    h5ad_path = resolve_cross_platform_path(args.h5ad)
    gene_path = resolve_cross_platform_path(args.genes)
    out_dir = resolve_cross_platform_path(args.output_dir)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Source H5AD not found: {h5ad_path}")
    if not gene_path.exists():
        raise FileNotFoundError(f"Gene file not found: {gene_path}")
    out_dir.mkdir(parents=True, exist_ok=True)

    genes = (
        pd.read_csv(gene_path)["gene_symbol"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .drop_duplicates()
        .tolist()
    )
    if len(genes) != 41:
        raise ValueError(f"Expected exactly 41 fresh genes, found {len(genes)}")
    pd.DataFrame({"gene_symbol": genes}).to_csv(out_dir / "input_41_genes.csv", index=False)

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        feature_names = adata.var["feature_name"].astype(str).str.upper()
        duplicated = feature_names[feature_names.isin(genes)].value_counts()
        duplicate_targets = duplicated[duplicated > 1].to_dict()
        present = [g for g in genes if (feature_names == g).sum() == 1]
        missing = sorted(set(genes) - set(present))
        if duplicate_targets:
            raise ValueError(f"Duplicated target features: {duplicate_targets}")
        if args.require_all_genes and missing:
            raise ValueError(f"Missing target genes: {missing}")

        target_indices = [int(np.flatnonzero(feature_names.eq(g).to_numpy())[0]) for g in present]
        obs = adata.obs.copy()
        obs["donor_id"] = obs["donor_id"].astype(str)
        obs["cell_type"] = obs["cell_type"].astype(str)
        obs["group"] = np.where(
            obs["disease"].astype(str).eq(TUMOR),
            TUMOR,
            np.where(obs["disease"].astype(str).eq(NORMAL), NORMAL, "outside_scope"),
        )
        obs["compartment"] = obs["cell_type"].map(classify_compartment)
        obs["subtype"] = obs.apply(classify_subtype, axis=1)
        primary = obs["is_primary_data"].astype(str).str.lower().isin(["true", "1"])
        eligible = (
            primary
            & obs["compartment"].isin(COMPARTMENTS)
            & obs["group"].isin([TUMOR, NORMAL])
            & ~obs["donor_id"].isin(["", "nan", "none"])
        )
        eligible_positions = np.flatnonzero(eligible.to_numpy())
        if len(eligible_positions) < 100:
            raise ValueError(f"Too few eligible cells: {len(eligible_positions)}")

        print(f"Reading {len(present)} target genes from {len(adata):,} cells ...", flush=True)
        values_all = read_selected_csr_columns(
            h5ad_path,
            row_count=int(adata.shape[0]),
            column_count=int(adata.shape[1]),
            target_indices=target_indices,
            chunk_size=args.chunk_size,
        )
        values = values_all[eligible_positions]
        meta = obs.iloc[eligible_positions].reset_index(drop=True)

        gene_mean = values.mean(axis=0, dtype=np.float64)
        gene_sd = values.std(axis=0, ddof=1, dtype=np.float64)
        score_genes = [g for g, sd in zip(present, gene_sd) if np.isfinite(sd) and sd > 0]
        usable_idx = np.array([i for i, sd in enumerate(gene_sd) if np.isfinite(sd) and sd > 0], dtype=int)
        if not len(usable_idx):
            raise ValueError("No non-constant target genes available for scoring")
        z = (values[:, usable_idx] - gene_mean[usable_idx]) / gene_sd[usable_idx]
        program_score = np.mean(z, axis=1, dtype=np.float64)

        meta["program_score"] = program_score
        meta["n_score_genes"] = len(score_genes)
        meta["score_gene_set"] = ";".join(score_genes)

        celltype_summary = summarize_scores(meta, ["compartment", "cell_type", "group"])
        compartment_summary = summarize_scores(meta, ["compartment", "group"])
        subtype_meta = meta.loc[meta["subtype"].notna()].copy()
        subtype_summary = summarize_scores(subtype_meta, ["compartment", "subtype", "group"])

        donor_scores = (
            meta.groupby(["donor_id", "group", "compartment", "subtype"], observed=True)
            .agg(
                program_score=("program_score", "mean"),
                sd_cell_score=("program_score", "std"),
                n_cells=("program_score", "size"),
            )
            .reset_index()
        )
        paired_summary, paired_deltas = paired_contrasts(donor_scores)

        gene_rows: list[dict] = []
        for (compartment, group), idx in meta.groupby(["compartment", "group"], observed=True).groups.items():
            vals = values[np.asarray(list(idx), dtype=int)]
            for j, gene in enumerate(present):
                gene_rows.append(
                    {
                        "compartment": compartment,
                        "group": group,
                        "gene_symbol": gene,
                        "n_cells": int(vals.shape[0]),
                        "mean_expression": float(vals[:, j].mean()),
                        "fraction_detected": float(np.mean(vals[:, j] > 0)),
                    }
                )
        gene_summary = pd.DataFrame(gene_rows)

        presence = pd.DataFrame(
            {
                "gene_symbol": genes,
                "present_in_source": [g in present for g in genes],
                "usable_for_score": [g in score_genes for g in genes],
            }
        )
        presence["source_feature_name"] = presence["gene_symbol"].where(presence["present_in_source"], "")

        metadata_audit = (
            meta.groupby(["compartment", "group"], observed=True)
            .agg(n_cells=("program_score", "size"), n_donors=("donor_id", "nunique"), n_cell_types=("cell_type", "nunique"))
            .reset_index()
        )

        celltype_summary.to_csv(out_dir / "singlecell_41_celltype_summary.csv", index=False)
        compartment_summary.to_csv(out_dir / "singlecell_41_compartment_summary.csv", index=False)
        subtype_summary.to_csv(out_dir / "singlecell_41_subtype_summary.csv", index=False)
        donor_scores.to_csv(out_dir / "singlecell_41_donor_scores.csv", index=False)
        paired_summary.to_csv(out_dir / "singlecell_41_paired_contrasts.csv", index=False)
        paired_deltas.to_csv(out_dir / "singlecell_41_paired_deltas.csv", index=False)
        gene_summary.to_csv(out_dir / "singlecell_41_gene_expression_summary.csv", index=False)
        metadata_audit.to_csv(out_dir / "singlecell_41_metadata_audit.csv", index=False)
        presence.to_csv(out_dir / "singlecell_41_gene_presence.csv", index=False)

        manifest = {
            "analysis": "Fresh 41-gene DINP-CRC single-cell localization",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "source_h5ad": str(h5ad_path),
            "source_h5ad_sha256": sha256_file(h5ad_path),
            "source_shape": [int(adata.shape[0]), int(adata.shape[1])],
            "source_title": str(adata.uns.get("title", "")),
            "source_citation": str(adata.uns.get("citation", "")),
            "gene_file": str(gene_path),
            "gene_file_sha256": sha256_file(gene_path),
            "target_genes_requested": len(genes),
            "target_genes_present": len(present),
            "target_genes_missing": missing,
            "score_genes": score_genes,
            "eligible_cells": int(len(meta)),
            "eligible_donors": int(meta["donor_id"].nunique()),
            "eligible_cell_types": int(meta["cell_type"].nunique()),
            "tumor_cells": int((meta["group"] == TUMOR).sum()),
            "normal_cells": int((meta["group"] == NORMAL).sum()),
            "score_definition": "Mean gene-wise z score across eligible primary colorectal cells; missing genes are excluded, not zero-filled.",
            "inference": "Donor-level mean scores with within-donor tumor-minus-normal contrasts when both groups are available.",
            "compartment_rule": "cell_type regex mapping; epithelial, myeloid, fibroblast, endothelial.",
            "malignancy_boundary": "Tumor-prefixed epithelial labels are source-labeled tumor epithelial, not definitive CNV-validated malignancy.",
            "fdr_family": "BH-FDR across the ten prespecified compartment/subtype paired contrasts.",
        }
        (out_dir / "singlecell_41_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

        lines = [
            "# Fresh 41-gene DINP–CRC single-cell localization",
            "",
            f"- Source object: `{h5ad_path.name}`; shape `{adata.shape[0]:,} cells × {adata.shape[1]:,} features`.",
            f"- Source title: {adata.uns.get('title', '')}",
            f"- Gene input: `query_41_genes.csv` ({len(genes)} genes).",
            f"- Source coverage: **{len(present)}/{len(genes)}** genes; missing: `{', '.join(missing) if missing else 'none'}`.",
            f"- Score genes: **{len(score_genes)}**; the score is the mean gene-wise z score over the observed score genes.",
            f"- Eligible scope: **{len(meta):,} cells**, **{meta['donor_id'].nunique():,} donors**, primary CRC/normal cells with mapped compartments.",
            "",
            "## Compartment localization",
            "",
            "```",
            compartment_summary.to_string(index=False),
            "```",
            "",
            "## Donor-aware paired contrasts",
            "",
            "```",
            paired_summary.to_string(index=False),
            "```",
            "",
            "## Interpretation boundary",
            "",
            "This is a localization analysis of the fresh 41-gene program in a versioned CRC single-cell reference. It does not establish that DINP exposure causes the program, that the program mediates the epidemiologic association, or that source-labeled tumor epithelial cells are definitively malignant without independent CNV/malignancy validation. The absent MIR675 feature is not treated as a negative expression measurement.",
            "",
            "## Reproducibility",
            "",
            "The complete input/output hashes and source metadata are in `singlecell_41_manifest.json`; cell-level expression was read from the source H5AD without re-querying Census or reprocessing the whole object.",
        ]
        (out_dir / "SINGLECELL_41_LOCALIZATION_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        print("\n41-GENE SINGLE-CELL LOCALIZATION: PASS")
        print(f"Observed genes: {len(present)}/{len(genes)}; missing: {missing or 'none'}")
        print(f"Eligible cells: {len(meta):,}; donors: {meta['donor_id'].nunique():,}")
        print(f"Output: {out_dir}")
        return 0
    finally:
        adata.file.close()


if __name__ == "__main__":
    raise SystemExit(main())
