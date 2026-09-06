#!/usr/bin/env python
"""Localize the frozen 81-gene DINP--CRC program across CRC compartments.

This is a targeted, donor-aware localization analysis.  It uses the
non-raw expression matrix in the official source H5AD, z-scores each of the
81 frozen genes across eligible cells, and averages the gene-level z-scores to
obtain a program score.  Inference is performed on donor-level means and,
when available, within-donor tumor-minus-normal contrasts.

The analysis deliberately does not call cells malignant.  The primary
epithelial label is ``tumor-derived epithelial`` and is inherited from the
source disease label plus the frozen cell-type classification rules.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs"
OUT_DIR = Path(__file__).resolve().parent / "outputs"
DEFAULT_H5AD = Path(
    r"D:\mcop and CRC\cellxgene_census\2025-11-08\16023185-de21-4c0d-a9c8-73abdd52d142.h5ad"
)
INPUT_GENES = INPUT_DIR / "dinp_crc_intersection.csv"
CENSUS_VERSION = "2025-11-08"
COMPARTMENTS = ["epithelial", "myeloid", "fibroblast", "endothelial"]
TUMOR_LABEL = "colon adenocarcinoma"
NORMAL_LABEL = "normal"
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


def classify_compartment(cell_type: object) -> str | None:
    value = str(cell_type).lower()
    if re.search(
        r"epithelial|colonocyte|enterocyte|goblet|paneth|enteroendocrine|tuft|"
        r"best4|transit amplifying|stem cell",
        value,
    ):
        return "epithelial"
    if re.search(r"macrophage|monocyte|dendritic|granulocyte|myeloid", value):
        return "myeloid"
    if re.search(r"fibroblast|myofibroblast", value):
        return "fibroblast"
    if re.search(r"endothelial", value):
        return "endothelial"
    return None


def dense_matrix(value: object) -> np.ndarray:
    if sparse.issparse(value):
        return value.toarray().astype(np.float64, copy=False)
    if hasattr(value, "toarray"):
        return value.toarray().astype(np.float64, copy=False)
    return np.asarray(value, dtype=np.float64)


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(values)
    if valid.any():
        result.loc[valid] = multipletests(values.loc[valid].to_numpy(dtype=float), method="fdr_bh")[1]
    return result


def summarize_group(accumulator: dict) -> dict:
    n = int(accumulator["n_cells"])
    score_mean = accumulator["score_sum"] / n
    score_var = max(
        (accumulator["score_sumsq"] - accumulator["score_sum"] ** 2 / n) / max(n - 1, 1),
        0.0,
    )
    record = {
        "n_cells": n,
        "mean_program_score": float(score_mean),
        "sd_cell_program_score": float(np.sqrt(score_var)) if n > 1 else np.nan,
    }
    for gene, value in zip(accumulator["genes"], accumulator["gene_sum"] / n):
        record[f"mean_expr_{gene}"] = float(value)
    for gene, value in zip(accumulator["genes"], accumulator["gene_detected"] / n):
        record[f"fraction_detected_{gene}"] = float(value)
    return record


def build_accumulator(genes: list[str]) -> dict:
    return {
        "genes": genes,
        "n_cells": 0,
        "score_sum": 0.0,
        "score_sumsq": 0.0,
        "gene_sum": np.zeros(len(genes), dtype=np.float64),
        "gene_detected": np.zeros(len(genes), dtype=np.float64),
    }


def paired_contrasts(donor_df: pd.DataFrame, score_column: str = "mean_program_score") -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    delta_rows = []
    for compartment in COMPARTMENTS:
        subset = donor_df.loc[donor_df["compartment"].eq(compartment)].copy()
        pivot = subset.pivot_table(index="donor_id", columns="group", values=score_column, aggfunc="first")
        if TUMOR_LABEL not in pivot.columns or NORMAL_LABEL not in pivot.columns:
            deltas = np.array([], dtype=float)
        else:
            paired = pivot[[TUMOR_LABEL, NORMAL_LABEL]].dropna()
            deltas = (paired[TUMOR_LABEL] - paired[NORMAL_LABEL]).to_numpy(dtype=float)
            for donor_id, row in paired.iterrows():
                delta_rows.append(
                    {
                        "compartment": compartment,
                        "donor_id": str(donor_id),
                        "tumor_score": float(row[TUMOR_LABEL]),
                        "normal_score": float(row[NORMAL_LABEL]),
                        "tumor_minus_normal": float(row[TUMOR_LABEL] - row[NORMAL_LABEL]),
                    }
                )
        n = int(len(deltas))
        mean_delta = float(np.mean(deltas)) if n else np.nan
        median_delta = float(np.median(deltas)) if n else np.nan
        sd_delta = float(np.std(deltas, ddof=1)) if n > 1 else np.nan
        if n > 1 and np.isfinite(sd_delta) and sd_delta > 0:
            t_stat, p_ttest = stats.ttest_1samp(deltas, 0.0)
            t_crit = stats.t.ppf(0.975, n - 1)
            half_width = t_crit * sd_delta / np.sqrt(n)
            ci_low, ci_high = mean_delta - half_width, mean_delta + half_width
        else:
            t_stat, p_ttest, ci_low, ci_high = np.nan, np.nan, np.nan, np.nan
        if n >= 5 and np.any(np.abs(deltas) > 0):
            try:
                w_stat, p_wilcoxon = stats.wilcoxon(deltas, alternative="two-sided", method="auto")
            except ValueError:
                w_stat, p_wilcoxon = np.nan, np.nan
        else:
            w_stat, p_wilcoxon = np.nan, np.nan
        summary_rows.append(
            {
                "compartment": compartment,
                "n_paired_donors": n,
                "mean_delta_tumor_minus_normal": mean_delta,
                "median_delta_tumor_minus_normal": median_delta,
                "sd_delta": sd_delta,
                "mean_delta_95ci_low": ci_low,
                "mean_delta_95ci_high": ci_high,
                "paired_t_statistic": float(t_stat) if np.isfinite(t_stat) else np.nan,
                "paired_t_p": float(p_ttest) if np.isfinite(p_ttest) else np.nan,
                "paired_wilcoxon_statistic": float(w_stat) if np.isfinite(w_stat) else np.nan,
                "paired_wilcoxon_p": float(p_wilcoxon) if np.isfinite(p_wilcoxon) else np.nan,
                "direction": "up" if mean_delta > 0 else "down" if mean_delta < 0 else "flat",
            }
        )
    result = pd.DataFrame(summary_rows)
    result["paired_t_BH_FDR"] = bh(result["paired_t_p"])
    result["paired_wilcoxon_BH_FDR"] = bh(result["paired_wilcoxon_p"])
    return result, pd.DataFrame(delta_rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    args = parser.parse_args()
    h5ad_path = args.h5ad.expanduser().resolve()
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Source H5AD not found: {h5ad_path}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    source = pd.read_csv(INPUT_GENES)
    genes = sorted(set(source["gene_symbol"].dropna().astype(str).str.upper()))
    if len(genes) != 81:
        raise ValueError(f"Expected 81 frozen genes, found {len(genes)}")
    pd.DataFrame({"gene_symbol": genes}).to_csv(OUT_DIR / "input_81_genes.csv", index=False)

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        feature_names = adata.var["feature_name"].astype(str).str.upper()
        missing = sorted(set(genes) - set(feature_names))
        duplicated = feature_names[feature_names.isin(genes)].value_counts()
        if missing or (duplicated > 1).any():
            raise ValueError(f"Gene mapping failed; missing={missing}, duplicated={duplicated[duplicated > 1].to_dict()}")
        gene_indices = [int(np.flatnonzero(feature_names.eq(g).to_numpy())[0]) for g in genes]

        obs = adata.obs.copy()
        obs["donor_id"] = obs["donor_id"].astype(str)
        obs["cell_type"] = obs["cell_type"].astype(str)
        obs["compartment"] = obs["cell_type"].map(classify_compartment)
        obs["group"] = np.where(
            obs["disease"].astype(str).eq(TUMOR_LABEL), TUMOR_LABEL,
            np.where(obs["disease"].astype(str).eq(NORMAL_LABEL), NORMAL_LABEL, "outside_scope"),
        )
        eligible = (
            obs["is_primary_data"].eq(True)
            & obs["compartment"].isin(COMPARTMENTS)
            & obs["group"].isin([TUMOR_LABEL, NORMAL_LABEL])
            & ~obs["donor_id"].isin(["", "nan", "None"])
        )
        obs.loc[eligible, "donor_id"] = obs.loc[eligible, "donor_id"].astype(str)
        n_eligible = int(eligible.sum())
        if n_eligible < 100:
            raise ValueError(f"Too few eligible cells: {n_eligible}")

        # A single targeted column read avoids repeatedly scanning the
        # CSR-backed H5AD for every row chunk.  The selected matrix is only
        # 370k x 81 and remains in memory as float32 (~120 MB).
        target_values = dense_matrix(adata.X[:, gene_indices]).astype(np.float32, copy=False)
        eligible_positions = np.flatnonzero(eligible.to_numpy())
        values = target_values[eligible_positions]
        meta = obs.iloc[eligible_positions].reset_index(drop=True)
        n_cells = int(values.shape[0])
        gene_sum = values.sum(axis=0, dtype=np.float64)
        gene_sumsq = np.square(values, dtype=np.float64).sum(axis=0, dtype=np.float64)
        gene_mean = gene_sum / n_cells
        gene_var = np.maximum((gene_sumsq - gene_sum**2 / n_cells) / max(n_cells - 1, 1), 0.0)
        gene_sd = np.sqrt(gene_var)
        gene_sd[gene_sd == 0] = 1.0

        donor_acc: dict[tuple[str, str, str], dict] = {}
        celltype_acc: dict[tuple[str, str, str], dict] = {}
        compartment_acc: dict[tuple[str, str], dict] = {}
        z = (values - gene_mean) / gene_sd
        scores = z.mean(axis=1)

        def accumulate(grouped, accumulator: dict) -> None:
            for key, row_indices in grouped.groups.items():
                indices = np.asarray(row_indices, dtype=int)
                item = build_accumulator(genes)
                item["n_cells"] = int(len(indices))
                item["score_sum"] = float(scores[indices].sum(dtype=np.float64))
                item["score_sumsq"] = float(np.square(scores[indices]).sum(dtype=np.float64))
                item["gene_sum"] = values[indices].sum(axis=0, dtype=np.float64)
                item["gene_detected"] = (values[indices] != 0).sum(axis=0, dtype=np.float64)
                accumulator[key] = item

        accumulate(meta.groupby(["donor_id", "group", "compartment"], sort=True, observed=True), donor_acc)
        accumulate(meta.groupby(["cell_type", "group", "compartment"], sort=True, observed=True), celltype_acc)
        accumulate(meta.groupby(["group", "compartment"], sort=True, observed=True), compartment_acc)

        # Sensitivity: standardize each gene within compartment before scoring.
        # This removes the possibility that a global cell-composition/baseline
        # difference creates an apparent compartment localization.
        within_z = np.zeros_like(values, dtype=np.float32)
        for compartment in COMPARTMENTS:
            compartment_mask = meta["compartment"].eq(compartment).to_numpy()
            compartment_values = values[compartment_mask]
            local_mean = compartment_values.mean(axis=0)
            local_sd = compartment_values.std(axis=0, ddof=1)
            local_sd[local_sd == 0] = 1.0
            within_z[compartment_mask] = (compartment_values - local_mean) / local_sd
        within_scores = within_z.mean(axis=1)
        within_donor_acc: dict[tuple[str, str, str], dict] = {}
        for key, row_indices in meta.groupby(["donor_id", "group", "compartment"], sort=True, observed=True).groups.items():
            indices = np.asarray(row_indices, dtype=int)
            item = build_accumulator(genes)
            item["n_cells"] = int(len(indices))
            item["score_sum"] = float(within_scores[indices].sum(dtype=np.float64))
            item["score_sumsq"] = float(np.square(within_scores[indices]).sum(dtype=np.float64))
            item["gene_sum"] = values[indices].sum(axis=0, dtype=np.float64)
            item["gene_detected"] = (values[indices] != 0).sum(axis=0, dtype=np.float64)
            within_donor_acc[key] = item
    finally:
        adata.file.close()

    donor_rows = []
    for (donor_id, group, compartment), accumulator in sorted(donor_acc.items()):
        donor_rows.append(
            {
                "donor_id": donor_id,
                "group": group,
                "compartment": compartment,
                **summarize_group(accumulator),
            }
        )
    donor_df = pd.DataFrame(donor_rows)
    donor_df.to_csv(OUT_DIR / "donor_compartment_program_scores.csv", index=False)

    within_donor_rows = []
    for (donor_id, group, compartment), accumulator in sorted(within_donor_acc.items()):
        summary = summarize_group(accumulator)
        within_donor_rows.append(
            {
                "donor_id": donor_id,
                "group": group,
                "compartment": compartment,
                "mean_program_score_within_compartment_z": summary["mean_program_score"],
                "sd_cell_program_score": summary["sd_cell_program_score"],
                "n_cells": summary["n_cells"],
            }
        )
    within_donor_df = pd.DataFrame(within_donor_rows)
    within_donor_df.to_csv(OUT_DIR / "donor_compartment_program_scores_within_compartment_z.csv", index=False)

    celltype_rows = []
    for (cell_type, group, compartment), accumulator in sorted(celltype_acc.items()):
        celltype_rows.append(
            {
                "cell_type": cell_type,
                "group": group,
                "compartment": compartment,
                **summarize_group(accumulator),
            }
        )
    pd.DataFrame(celltype_rows).to_csv(OUT_DIR / "celltype_program_score_summary.csv", index=False)

    compartment_rows = []
    for (group, compartment), accumulator in sorted(compartment_acc.items()):
        compartment_rows.append(
            {
                "group": group,
                "compartment": compartment,
                **summarize_group(accumulator),
            }
        )
    pd.DataFrame(compartment_rows).to_csv(OUT_DIR / "pooled_compartment_program_score_summary.csv", index=False)

    paired_summary, paired_deltas = paired_contrasts(donor_df)
    paired_summary.to_csv(OUT_DIR / "paired_compartment_contrasts.csv", index=False)
    paired_deltas.to_csv(OUT_DIR / "paired_donor_deltas.csv", index=False)
    within_paired_summary, within_paired_deltas = paired_contrasts(
        within_donor_df, score_column="mean_program_score_within_compartment_z"
    )
    within_paired_summary.to_csv(OUT_DIR / "paired_compartment_contrasts_within_compartment_z.csv", index=False)
    within_paired_deltas.to_csv(OUT_DIR / "paired_donor_deltas_within_compartment_z.csv", index=False)

    gene_map = pd.DataFrame({
        "gene_symbol": genes,
        "feature_name": genes,
        "feature_index": gene_indices,
        "present_in_source_h5ad": True,
    })
    gene_map.to_csv(OUT_DIR / "gene_mapping_audit.csv", index=False)

    primary_paired = paired_summary.sort_values("paired_t_BH_FDR", na_position="last").head(1)
    source_counts = meta.groupby(["group", "compartment"], sort=True, observed=True).size().reset_index(name="n_cells")
    donor_counts = donor_df.groupby(["group", "compartment"], as_index=False).agg(n_donors=("donor_id", "nunique"), n_cells=("n_cells", "sum"))
    source_counts.to_csv(OUT_DIR / "eligible_cell_counts.csv", index=False)
    donor_counts.to_csv(OUT_DIR / "eligible_donor_counts.csv", index=False)

    manifest = {
        "analysis": "Frozen 81-gene DINP-CRC program single-cell compartment localization",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_gene_file": str(INPUT_GENES),
        "input_gene_count": len(genes),
        "input_gene_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
        "source_h5ad": str(h5ad_path),
        "source_h5ad_sha256": sha256_file(h5ad_path),
        "source_h5ad_bytes": int(h5ad_path.stat().st_size),
        "census_version": CENSUS_VERSION,
        "source_shape": [int(adata.n_obs), int(adata.n_vars)],
        "expression_matrix": "adata.X (non-raw source matrix); adata.raw was not used for scoring",
        "raw_matrix_available": True,
        "primary_filter": "is_primary_data == True; disease in {colon adenocarcinoma, normal}; donor_id present; four frozen compartment rules",
        "tumor_label": TUMOR_LABEL,
        "normal_label": NORMAL_LABEL,
        "compartment_rules": {
            "epithelial": "epithelial|colonocyte|enterocyte|goblet|paneth|enteroendocrine|tuft|best4|transit amplifying|stem cell",
            "myeloid": "macrophage|monocyte|dendritic|granulocyte|myeloid",
            "fibroblast": "fibroblast|myofibroblast",
            "endothelial": "endothelial",
        },
        "eligible_cell_count": n_eligible,
        "gene_standardization": "upper-cased frozen symbols; exact unique feature_name match",
        "score_definition": "mean of gene-wise z-scores across the 81 frozen genes, standardized across all eligible cells",
        "score_sensitivity": "gene-wise standardization repeated within each compartment; reported separately",
        "unit_of_inference": "donor-level mean score; primary contrast is within-donor tumor minus normal",
        "multiple_testing": "BH-FDR across four compartment paired t-tests; Wilcoxon FDR also reported",
        "malignant_annotation_boundary": "not inferred; epithelial cells are labeled tumor-derived epithelial when disease label is tumor",
        "outputs": [str(p) for p in sorted(OUT_DIR.glob("*.csv"))],
        "headline_paired_result": primary_paired.to_dict("records"),
        "within_compartment_standardization_paired_result": within_paired_summary.to_dict("records"),
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    top = paired_summary.sort_values("paired_t_BH_FDR", na_position="last")
    report_lines = [
        "# Frozen 81-gene DINP–CRC program: single-cell compartment localization",
        "",
        f"Generated: {manifest['run_timestamp_utc']}",
        "",
        "## Analysis boundary",
        "",
        f"- Source: official Census-release source H5AD, pinned to `{CENSUS_VERSION}`; source file is not copied into the repository.",
        f"- Source shape: `{adata.n_obs:,} cells × {adata.n_vars:,} features`; eligible cells across the four compartments: `{n_eligible:,}`.",
        f"- Frozen program: 81 genes from `dinp_crc_intersection.csv`; all 81 had unique exact matches in `feature_name`.",
        "- Expression scale: `adata.X` (non-raw source matrix). Each gene was z-scored across all eligible cells, then averaged into the program score.",
        "- Inference: donor-level means; the primary localization contrast is paired tumor-minus-normal within donor. BH-FDR is across the four compartment paired t-tests.",
        "- Sensitivity: the same 81 genes were standardized separately within each compartment before scoring; this is not pooled with the primary score.",
        "- The analysis does not infer malignant status. The epithelial label is tumor-derived epithelial when the source disease label is colon adenocarcinoma.",
        "",
        "## Paired tumor–normal localization",
        "",
        "| Compartment | Paired donors | Mean Δ (tumor−normal) | 95% CI | t-test P | BH-FDR | Wilcoxon P | Direction |",
        "|---|---:|---:|---|---:|---:|---:|---|",
    ]
    for _, row in top.iterrows():
        ci = f"{row['mean_delta_95ci_low']:.3f} to {row['mean_delta_95ci_high']:.3f}" if pd.notna(row["mean_delta_95ci_low"]) else "NA"
        report_lines.append(
            f"| {row['compartment']} | {int(row['n_paired_donors'])} | {row['mean_delta_tumor_minus_normal']:.3f} | {ci} | {row['paired_t_p']:.3g} | {row['paired_t_BH_FDR']:.3g} | {row['paired_wilcoxon_p']:.3g} | {row['direction']} |"
        )
    report_lines += [
        "",
        "## Within-compartment standardization sensitivity",
        "",
        "| Compartment | Paired donors | Mean Δ (tumor−normal) | t-test P | BH-FDR | Direction |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for _, row in within_paired_summary.sort_values("paired_t_BH_FDR", na_position="last").iterrows():
        report_lines.append(
            f"| {row['compartment']} | {int(row['n_paired_donors'])} | {row['mean_delta_tumor_minus_normal']:.3f} | {row['paired_t_p']:.3g} | {row['paired_t_BH_FDR']:.3g} | {row['direction']} |"
        )
    report_lines += [
        "",
        "## Interpretation boundary",
        "",
        "This is a localization/convergence analysis of the frozen DINP–CRC intersection. A compartment-level expression shift does not establish that DINP exposure causes the shift, nor that the program mediates the epidemiologic association.",
        "",
        "Detailed donor-level scores, paired deltas, cell-type summaries, source hashes, and gene mapping audit are retained in the output directory.",
    ]
    (OUT_DIR / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "source_shape": [int(adata.n_obs), int(adata.n_vars)],
        "input_genes": len(genes),
        "eligible_cells": n_eligible,
        "paired_summary": paired_summary.to_dict("records"),
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
