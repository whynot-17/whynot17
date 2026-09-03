#!/usr/bin/env python
"""Subtype localization of the frozen 81-gene DINP--CRC program.

This analysis drills the previously completed broad-compartment result into
source-labeled tumor epithelial cells and myeloid subtypes.  It keeps the
same 81-gene score, source H5AD, primary-data filter, and donor-level paired
inference used by the broad localization analysis.

Important boundary: ``ClusterFull`` labels beginning with ``Tumor`` are
reported as source-labeled tumor epithelial / malignant-candidate cells.  No
CNV inference or independent malignant-cell validation is performed here, so
the analysis does not claim definitive malignant status.
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
TUMOR_LABEL = "colon adenocarcinoma"
NORMAL_LABEL = "normal"
COMPARTMENTS = ["epithelial", "myeloid", "fibroblast", "endothelial"]
LINEAGES = ["epithelial", "myeloid"]


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


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
        result.loc[valid] = multipletests(
            values.loc[valid].to_numpy(dtype=float), method="fdr_bh"
        )[1]
    return result


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


def classify_subtype(row: pd.Series) -> str | None:
    compartment = row["compartment"]
    disease = row["disease_label"]
    cell_type = str(row["cell_type"])
    tumor_labeled = bool(row["source_tumor_label"])
    if compartment == "epithelial":
        if disease == TUMOR_LABEL:
            return "source_tumor_labeled_epithelial" if tumor_labeled else "other_tumor_epithelial"
        if disease == NORMAL_LABEL:
            return "normal_epithelial"
    if compartment == "myeloid":
        lowered = cell_type.lower()
        if lowered == "macrophage":
            return "macrophage"
        if lowered == "monocyte":
            return "monocyte"
        if "dendritic" in lowered:
            return "dendritic"
        if lowered == "granulocyte":
            return "granulocyte"
    return None


def summarize_donor_scores(meta: pd.DataFrame, scores: np.ndarray, score_name: str) -> pd.DataFrame:
    frame = meta[["donor_id", "group", "disease_label", "compartment", "subtype"]].copy()
    frame[score_name] = scores
    frame["n_cells"] = 1
    summary = (
        frame.loc[frame["subtype"].notna()]
        .groupby(["donor_id", "group", "disease_label", "compartment", "subtype"], observed=True)
        .agg(
            **{
                score_name: (score_name, "mean"),
                "sd_cell_score": (score_name, "std"),
                "n_cells": ("n_cells", "sum"),
            }
        )
        .reset_index()
    )
    return summary


def paired_contrasts(
    donor_df: pd.DataFrame, score_column: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    contrast_specs = [
        ("source_tumor_labeled_epithelial", "epithelial", "source tumor-labeled epithelial vs normal epithelial"),
        ("other_tumor_epithelial", "epithelial", "other tumor epithelial vs normal epithelial"),
        ("macrophage", "myeloid", "macrophage tumor vs normal"),
        ("monocyte", "myeloid", "monocyte tumor vs normal"),
        ("dendritic", "myeloid", "dendritic tumor vs normal"),
        ("granulocyte", "myeloid", "granulocyte tumor vs normal"),
    ]
    summary_rows: list[dict] = []
    delta_rows: list[dict] = []
    for subtype, lineage, label in contrast_specs:
        target = donor_df.loc[donor_df["subtype"].eq(subtype)].copy()
        if lineage == "epithelial":
            normal = donor_df.loc[donor_df["subtype"].eq("normal_epithelial")].copy()
        else:
            normal = donor_df.loc[donor_df["subtype"].eq(subtype)].copy()
        target = target.loc[target["group"].eq("tumor")]
        normal = normal.loc[normal["group"].eq("normal")]
        target_scores = target.set_index("donor_id")[score_column].rename("tumor")
        normal_scores = normal.set_index("donor_id")[score_column].rename("normal")
        paired = pd.concat([target_scores, normal_scores], axis=1).dropna()
        deltas = (paired["tumor"] - paired["normal"]).to_numpy(dtype=float)
        for donor_id, row in paired.iterrows():
            delta_rows.append(
                {
                    "contrast": label,
                    "subtype": subtype,
                    "lineage": lineage,
                    "donor_id": str(donor_id),
                    "tumor_score": float(row["tumor"]),
                    "normal_score": float(row["normal"]),
                    "tumor_minus_normal": float(row["tumor"] - row["normal"]),
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
                "contrast": label,
                "subtype": subtype,
                "lineage": lineage,
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


def build_label_audit(obs: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    epithelial = obs["compartment"].eq("epithelial")
    for disease in [TUMOR_LABEL, NORMAL_LABEL]:
        for source_label in [False, True]:
            rows.append(
                {
                    "audit": "epithelial_disease_by_source_tumor_label",
                    "disease": disease,
                    "source_tumor_label": source_label,
                    "n_cells": int((epithelial & obs["disease_label"].eq(disease) & obs["source_tumor_label"].eq(source_label)).sum()),
                }
            )
    for subtype in [
        "source_tumor_labeled_epithelial",
        "other_tumor_epithelial",
        "normal_epithelial",
        "macrophage",
        "monocyte",
        "dendritic",
        "granulocyte",
    ]:
        subset = obs.loc[obs["subtype"].eq(subtype)]
        for disease, count in subset["disease_label"].value_counts().items():
            rows.append(
                {
                    "audit": "subtype_by_disease",
                    "subtype": subtype,
                    "disease": disease,
                    "n_cells": int(count),
                }
            )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
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
            raise ValueError(
                f"Gene mapping failed; missing={missing}, duplicated={duplicated[duplicated > 1].to_dict()}"
            )
        gene_indices = [int(np.flatnonzero(feature_names.eq(g).to_numpy())[0]) for g in genes]

        obs = adata.obs.copy()
        obs["donor_id"] = obs["donor_id"].astype(str)
        obs["cell_type"] = obs["cell_type"].astype(str)
        obs["disease_label"] = obs["disease"].astype(str)
        obs["compartment"] = obs["cell_type"].map(classify_compartment)
        obs["source_tumor_label"] = obs["ClusterFull"].astype(str).str.startswith("Tumor")
        obs["group"] = np.where(
            obs["disease_label"].eq(TUMOR_LABEL), "tumor",
            np.where(obs["disease_label"].eq(NORMAL_LABEL), "normal", "outside_scope"),
        )
        obs["subtype"] = obs.apply(classify_subtype, axis=1)
        eligible = (
            obs["is_primary_data"].eq(True)
            & obs["compartment"].isin(COMPARTMENTS)
            & obs["group"].isin(["tumor", "normal"])
            & ~obs["donor_id"].isin(["", "nan", "None"])
        )
        eligible_positions = np.flatnonzero(eligible.to_numpy())
        meta = obs.iloc[eligible_positions].reset_index(drop=True)
        if len(meta) < 100:
            raise ValueError(f"Too few eligible cells: {len(meta)}")

        target_values = dense_matrix(adata.X[:, gene_indices]).astype(np.float32, copy=False)
        values = target_values[eligible_positions]
        n_cells = int(values.shape[0])
        gene_mean = values.mean(axis=0, dtype=np.float64)
        gene_sd = values.std(axis=0, ddof=1, dtype=np.float64)
        gene_sd[gene_sd == 0] = 1.0
        global_scores = ((values - gene_mean) / gene_sd).mean(axis=1)

        within_scores = np.zeros_like(values, dtype=np.float32)
        for compartment in COMPARTMENTS:
            mask = meta["compartment"].eq(compartment).to_numpy()
            compartment_values = values[mask]
            local_mean = compartment_values.mean(axis=0)
            local_sd = compartment_values.std(axis=0, ddof=1)
            local_sd[local_sd == 0] = 1.0
            within_scores[mask] = (compartment_values - local_mean) / local_sd
        within_compartment_scores = within_scores.mean(axis=1)

        donor_df = summarize_donor_scores(meta, global_scores, "global_program_score")
        within_donor_df = summarize_donor_scores(
            meta, within_compartment_scores, "within_compartment_program_score"
        )
    finally:
        adata.file.close()

    label_audit = build_label_audit(meta)
    label_audit.to_csv(OUT_DIR / "subtype_label_audit.csv", index=False)

    donor_df.to_csv(OUT_DIR / "donor_subtype_program_scores.csv", index=False)
    within_donor_df.to_csv(
        OUT_DIR / "donor_subtype_program_scores_within_compartment_z.csv", index=False
    )

    paired_summary, paired_deltas = paired_contrasts(donor_df, "global_program_score")
    within_paired_summary, within_paired_deltas = paired_contrasts(
        within_donor_df, "within_compartment_program_score"
    )
    paired_summary.to_csv(OUT_DIR / "paired_subtype_contrasts.csv", index=False)
    paired_deltas.to_csv(OUT_DIR / "paired_subtype_deltas.csv", index=False)
    within_paired_summary.to_csv(
        OUT_DIR / "paired_subtype_contrasts_within_compartment_z.csv", index=False
    )
    within_paired_deltas.to_csv(
        OUT_DIR / "paired_subtype_deltas_within_compartment_z.csv", index=False
    )

    cell_counts = (
        meta.loc[meta["subtype"].notna()]
        .groupby(["subtype", "compartment", "group", "disease_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    cell_counts["n_donors"] = cell_counts.apply(
        lambda row: int(
            meta.loc[
                meta["subtype"].eq(row["subtype"])
                & meta["compartment"].eq(row["compartment"])
                & meta["group"].eq(row["group"]),
                "donor_id",
            ].nunique()
        ),
        axis=1,
    )
    cell_counts.to_csv(OUT_DIR / "subtype_cell_counts.csv", index=False)
    donor_counts = (
        donor_df.groupby(["subtype", "compartment", "group", "disease_label"], observed=True)
        .agg(n_donors=("donor_id", "nunique"), n_cells=("n_cells", "sum"))
        .reset_index()
    )
    donor_counts.to_csv(OUT_DIR / "subtype_donor_counts.csv", index=False)

    gene_map = pd.DataFrame(
        {
            "gene_symbol": genes,
            "feature_name": genes,
            "feature_index": gene_indices,
            "present_in_source_h5ad": True,
        }
    )
    gene_map.to_csv(OUT_DIR / "gene_mapping_audit.csv", index=False)

    source_hash = sha256_file(h5ad_path)
    manifest = {
        "analysis": "Frozen 81-gene DINP-CRC program single-cell subtype localization",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_gene_file": str(INPUT_GENES),
        "input_gene_count": len(genes),
        "input_gene_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
        "source_h5ad": str(h5ad_path),
        "source_h5ad_sha256": source_hash,
        "source_h5ad_bytes": int(h5ad_path.stat().st_size),
        "census_version": CENSUS_VERSION,
        "source_shape": [int(adata.n_obs), int(adata.n_vars)],
        "expression_matrix": "adata.X (non-raw source matrix); adata.raw was not used for scoring",
        "primary_data_filter": "is_primary_data == True; disease in {colon adenocarcinoma, normal}; donor_id present",
        "score_definition": "mean of gene-wise z-scores across all eligible cells, matching broad localization analysis",
        "sensitivity_score_definition": "mean of gene-wise z-scores after standardization within each broad compartment",
        "unit_of_inference": "donor-level subtype mean; paired tumor-minus-normal within donor",
        "subtype_definitions": {
            "source_tumor_labeled_epithelial": "epithelial cell-type rule AND disease=colon adenocarcinoma AND ClusterFull startswith Tumor",
            "other_tumor_epithelial": "epithelial cell-type rule AND disease=colon adenocarcinoma AND ClusterFull does not start with Tumor",
            "normal_epithelial": "epithelial cell-type rule AND disease=normal",
            "macrophage": "cell_type == macrophage",
            "monocyte": "cell_type == monocyte",
            "dendritic": "cell_type contains dendritic",
            "granulocyte": "cell_type == granulocyte",
        },
        "malignant_annotation_boundary": "ClusterFull Tumor prefix is source-labeled tumor epithelial / malignant-candidate only; no CNV or independent malignant validation was performed",
        "multiple_testing": "BH-FDR across six subtype paired t-tests; Wilcoxon FDR reported separately",
        "eligible_cell_count_all_four_compartments": n_cells,
        "eligible_subtype_cell_count": int(meta["subtype"].notna().sum()),
        "outputs": [str(p) for p in sorted(OUT_DIR.glob("*.csv"))],
        "source_hash_verified": True,
    }
    (OUT_DIR / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    report_lines = [
        "# Frozen 81-gene DINP–CRC program: single-cell subtype localization",
        "",
        f"Generated: {manifest['run_timestamp_utc']}",
        "",
        "## Analysis boundary",
        "",
        f"- Source: Census release `{CENSUS_VERSION}` H5AD; source file is not copied into the repository.",
        f"- Source shape: `{manifest['source_shape'][0]:,} cells × {manifest['source_shape'][1]:,} features`; four-compartment eligible cells: `{n_cells:,}`.",
        f"- Frozen program: 81 genes from `dinp_crc_intersection.csv`; all 81 had unique exact matches in `feature_name`.",
        "- Primary score: the same global gene-wise z-score score used in the broad-compartment localization analysis.",
        "- Inference: donor-level means and paired tumor-minus-normal contrasts; BH-FDR is across six subtype paired t-tests.",
        "- Sensitivity: genes standardized within each broad compartment before scoring.",
        "",
        "## Label audit",
        "",
        "`ClusterFull` labels beginning with `Tumor` were used to define a source-labeled tumor epithelial / malignant-candidate subgroup. This is not treated as definitive malignant status because no CNV-based or independent malignant-cell validation was performed.",
        "",
        "## Paired tumor–normal subtype contrasts",
        "",
        "| Contrast | Paired donors | Mean Δ (tumor−normal) | 95% CI | t-test P | BH-FDR | Wilcoxon P | Direction |",
        "|---|---:|---:|---|---:|---:|---:|---|",
    ]
    for _, row in paired_summary.sort_values("paired_t_BH_FDR", na_position="last").iterrows():
        ci = (
            f"{row['mean_delta_95ci_low']:.3f} to {row['mean_delta_95ci_high']:.3f}"
            if pd.notna(row["mean_delta_95ci_low"])
            else "NA"
        )
        report_lines.append(
            f"| {row['contrast']} | {int(row['n_paired_donors'])} | {row['mean_delta_tumor_minus_normal']:.3f} | {ci} | {row['paired_t_p']:.3g} | {row['paired_t_BH_FDR']:.3g} | {row['paired_wilcoxon_p']:.3g} | {row['direction']} |"
        )
    report_lines += [
        "",
        "## Within-compartment standardization sensitivity",
        "",
        "| Contrast | Paired donors | Mean Δ (tumor−normal) | t-test P | BH-FDR | Direction |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for _, row in within_paired_summary.sort_values("paired_t_BH_FDR", na_position="last").iterrows():
        report_lines.append(
            f"| {row['contrast']} | {int(row['n_paired_donors'])} | {row['mean_delta_tumor_minus_normal']:.3f} | {row['paired_t_p']:.3g} | {row['paired_t_BH_FDR']:.3g} | {row['direction']} |"
        )
    report_lines += [
        "",
        "## Interpretation boundary",
        "",
        "This analysis localizes the frozen DINP–CRC program in a CRC single-cell reference. It does not establish that DINP exposure causes the program, that the program mediates the epidemiologic association, or that the source Tumor-prefixed cells are definitively malignant without independent validation.",
        "",
        "Detailed donor-level scores, paired deltas, label audit, cell counts, source hash, and gene mapping audit are retained in the output directory.",
    ]
    (OUT_DIR / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "source_shape": manifest["source_shape"],
                "input_genes": len(genes),
                "eligible_cells_all_four_compartments": n_cells,
                "eligible_subtype_cells": int(meta["subtype"].notna().sum()),
                "paired_summary": paired_summary.to_dict("records"),
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
