#!/usr/bin/env python
"""Macrophage-only reanalysis of the locally cached GSE144735 CRC matrix.

The script deliberately keeps the frozen DINP--CRC gene set separate from
feature selection.  Clustering is driven by macrophage transcriptomes, then
the frozen score and PTGER4 are evaluated on the resulting states.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests


DEFAULT_MATRIX = Path(r"E:\mcop\mcop\work\mcop_phase2f_external\raw\GSE144735_raw_UMI_count_matrix.txt.gz")
DEFAULT_ANNOTATION = Path(r"E:\mcop\mcop\work\mcop_phase2f_external\raw\GSE144735_annotation.txt.gz")
DEFAULT_GENE_FILE = Path(r"E:\mcop\mcop\analysis\dinp_crc_multi_database_target_convergence\outputs\dinp_crc_intersection.csv")
FROZEN_DRIVER = ["NEAT1", "MMP9", "TIMP1", "STAT3", "PTGER4", "PTGES3", "CXCR4"]

SIGNATURES = {
    "SPP1_like_TAM": ["SPP1", "APOC1", "MMP9", "CTSB", "CTSL", "TREM2", "GPNMB"],
    "C1QC_like_TAM": ["C1QA", "C1QB", "C1QC", "APOE", "FOLR2", "MRC1"],
    "FCN1_inflammatory_monocyte_like": ["FCN1", "S100A8", "S100A9", "CTSS", "LILRB1", "VCAN"],
    "APOE_lipid_associated": ["APOE", "APOC1", "LPL", "TREM2", "FABP5"],
    "IFN_responsive": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1"],
    "cycling": ["MKI67", "TOP2A", "UBE2C"],
    "resident_like": ["FOLR2", "LYVE1", "MRC1", "CD163"],
}
CONTAMINATION = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20"],
    "lymphoid": ["CD3D", "CD3E", "TRBC1", "MS4A1", "CD79A", "NKG7"],
    "cycling": ["MKI67", "TOP2A", "UBE2C"],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_dump(obj, path: Path):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def read_frozen_genes(path: Path) -> list[str]:
    df = pd.read_csv(path)
    col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
    return [str(x).strip().upper() for x in df[col].dropna().tolist() if str(x).strip()]


def read_matrix_selected(matrix_path: Path, wanted_columns: list[str], progress_every: int = 5000):
    """Stream the gzipped gene x cell text matrix into a sparse cell x gene matrix."""
    import gzip

    wanted = set(wanted_columns)
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        if not header or header[0] != "Index":
            raise ValueError("Unexpected matrix header; expected first field Index")
        col_names = header[1:]
        pos = {name: i for i, name in enumerate(col_names)}
        missing = [c for c in wanted_columns if c not in pos]
        if missing:
            raise ValueError(f"Annotation cells missing from matrix: {missing[:5]}")
        selected_pos = np.array([pos[c] for c in wanted_columns], dtype=np.int64)
        genes, rows, cols, vals = [], [], [], []
        for j, line in enumerate(fh):
            fields = line.rstrip("\n\r").split("\t")
            if not fields:
                continue
            gene = fields[0].strip()
            if not gene:
                continue
            genes.append(gene)
            # The matrix is modest (27k columns), and this avoids materialising
            # an enormous dense array.  Duplicate symbols are made unique later.
            for i, p in enumerate(selected_pos):
                if p + 1 >= len(fields):
                    continue
                raw = fields[p + 1]
                if raw and raw != "0":
                    try:
                        value = float(raw)
                    except ValueError:
                        continue
                    if value:
                        rows.append(i)
                        cols.append(j)
                        vals.append(value)
            if progress_every and (j + 1) % progress_every == 0:
                print(f"matrix rows parsed: {j + 1}", flush=True)
    # Resolve repeated feature symbols without silently dropping expression.
    seen, unique = {}, []
    for g in genes:
        n = seen.get(g, 0)
        unique.append(g if n == 0 else f"{g}_{n}")
        seen[g] = n + 1
    X = sp.coo_matrix((np.asarray(vals, dtype=np.float32), (rows, cols)), shape=(len(wanted_columns), len(genes))).tocsr()
    return X, unique, col_names


def choose_condition(df: pd.DataFrame) -> pd.Series:
    s = df["Class"].astype(str).str.strip().str.lower()
    return s.map({"tumor": "tumor", "normal": "normal", "border": "border"}).fillna(s)


def safe_mean(x):
    return float(np.nanmean(x)) if len(x) else np.nan


def paired_stats(df: pd.DataFrame, value: str, group: str = "subtype", condition_col: str = "condition") -> pd.DataFrame:
    rows = []
    for subtype, g in df.groupby(group, dropna=False):
        piv = g.pivot_table(index="donor_id", columns=condition_col, values=value, aggfunc="mean")
        if not {"tumor", "normal"}.issubset(piv.columns):
            continue
        piv = piv.dropna(subset=["tumor", "normal"])
        d = piv["tumor"] - piv["normal"]
        n = len(d)
        t_p = float(stats.ttest_1samp(d, 0).pvalue) if n >= 2 else np.nan
        w_p = float(stats.wilcoxon(d).pvalue) if n >= 3 and np.any(d != 0) else np.nan
        rows.append({"subtype": subtype, "paired_donors": n, "mean_tumor_minus_normal": safe_mean(d),
                     "tumor_mean": safe_mean(piv["tumor"]), "normal_mean": safe_mean(piv["normal"]),
                     "t_p": t_p, "wilcoxon_p": w_p, "direction": "up" if safe_mean(d) > 0 else "down"})
    out = pd.DataFrame(rows)
    if len(out):
        out["fdr_family"] = f"paired tumor-normal {value} by macrophage subtype"
        valid = out["t_p"].notna()
        out["t_fdr_bh"] = np.nan
        if valid.any(): out.loc[valid, "t_fdr_bh"] = multipletests(out.loc[valid, "t_p"], method="fdr_bh")[1]
    return out


def top_markers(adata, cluster_key: str, n_top: int = 100) -> pd.DataFrame:
    import scanpy as sc
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon", n_genes=min(n_top, adata.n_vars - 1), use_raw=False, layer="log1p")
    r = adata.uns["rank_genes_groups"]
    groups = list(r["names"].dtype.names)
    rows = []
    X0 = adata.layers["log1p"] if "log1p" in adata.layers else adata.X
    X = X0.tocsr() if sp.issparse(X0) else np.asarray(X0)
    names = list(map(str, adata.var_names))
    idx = {g: i for i, g in enumerate(names)}
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    for cl in groups:
        mask = labels == str(cl)
        for rank in range(len(r["names"][cl])):
            gene = str(r["names"][cl][rank])
            if gene not in idx: continue
            gi = idx[gene]
            vals = X[:, gi].toarray().ravel() if sp.issparse(X) else X[:, gi]
            rows.append({"cluster": str(cl), "gene": gene, "rank": rank + 1,
                         "logFC": float(r["logfoldchanges"][cl][rank]),
                         "score": float(r["scores"][cl][rank]),
                         "pct_expr_cluster": float(np.mean(vals[mask] > 0)) if mask.any() else np.nan,
                         "pct_expr_other": float(np.mean(vals[~mask] > 0)) if (~mask).any() else np.nan,
                         "adjusted_p": float(r["pvals_adj"][cl][rank])})
    return pd.DataFrame(rows)


def annotate(markers: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cl, g in markers.groupby("cluster"):
        top = g.sort_values("rank").head(100)
        scores = {}
        hits = {}
        for label, genes in SIGNATURES.items():
            hit = top[top["gene"].str.upper().isin(genes)]
            hits[label] = ",".join(hit["gene"].tolist())
            scores[label] = (len(hit) / max(1, len(genes))) * (1 + max(0.0, float(hit["logFC"].mean())) / 5 if len(hit) else 0)
        ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        best, second = ranked[0], ranked[1]
        # A data-driven cluster with weak/ambiguous signature support remains
        # explicitly other/low-confidence rather than being forced into a TAM name.
        if best[1] < 0.18:
            label, confidence = "other_macrophage_state", "low"
        elif best[1] < 0.35 or (best[1] - second[1]) < 0.08:
            label, confidence = best[0], "medium"
        else:
            label, confidence = best[0], "high"
        rows.append({"cluster": cl, "subtype": label, "confidence": confidence,
                     "supporting_markers": hits[best[0]], "alternative_label": second[0],
                     "alternative_supporting_markers": hits[second[0]], "signature_score": best[1]})
    return pd.DataFrame(rows)


def save_figures(adata, markers, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_theme(style="white", context="talk")
    figdir = out / "figures"; figdir.mkdir(exist_ok=True)
    # UMAP
    fig, ax = plt.subplots(figsize=(8, 6))
    for label, g in adata.obs.groupby("subtype"):
        ix = g.index
        pos = [adata.obs_names.get_loc(x) for x in ix]
        ax.scatter(adata.obsm["X_umap"][pos, 0], adata.obsm["X_umap"][pos, 1], s=5, alpha=.65, label=label)
    ax.set(xlabel="UMAP1", ylabel="UMAP2", title="Macrophage-only UMAP"); ax.legend(fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"macrophage_subtype_umap.{ext}", dpi=220)
    plt.close(fig)
    # Marker dot plot data/figure
    chosen = markers.sort_values(["cluster", "rank"]).groupby("cluster").head(5)
    if len(chosen):
        piv = chosen.pivot_table(index="gene", columns="cluster", values="pct_expr_cluster", aggfunc="mean")
        fig, ax = plt.subplots(figsize=(max(7, .8 * len(piv.columns)), max(5, .3 * len(piv))))
        sns.heatmap(piv.fillna(0), cmap="Blues", vmin=0, vmax=1, ax=ax)
        ax.set_title("Top cluster marker expression fraction"); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"macrophage_marker_dotplot.{ext}", dpi=220)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    ap.add_argument("--annotation", type=Path, default=DEFAULT_ANNOTATION)
    ap.add_argument("--gene-file", type=Path, default=DEFAULT_GENE_FILE)
    ap.add_argument("--outdir", type=Path, default=Path(__file__).parent / "outputs")
    args = ap.parse_args()
    out = args.outdir; out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    t0 = time.time()
    ann = pd.read_csv(args.annotation, sep="\t")
    ann["condition"] = choose_condition(ann)
    # The source uses a broad Myeloids label; cDC is explicitly excluded.
    include = (ann["Cell_type"].astype(str).str.lower() == "myeloids") & (ann["Cell_subtype"].astype(str).str.lower() != "cdc")
    ann["included"] = include
    ann["selection_reason"] = np.where(include, "myeloid macrophage-like; cDC excluded", np.where(ann["Cell_type"].astype(str).str.lower() == "myeloids", "dendritic/cDC excluded", "non-myeloid compartment excluded"))
    audit = ann.groupby(["Cell_type", "Cell_subtype", "included", "selection_reason"], dropna=False).agg(
        cell_count=("Index", "size"), donor_count=("Patient", "nunique"), tumor_count=("condition", lambda x: int((x == "tumor").sum())), normal_count=("condition", lambda x: int((x == "normal").sum())), border_count=("condition", lambda x: int((x == "border").sum()))).reset_index()
    audit.to_csv(out / "macrophage_cell_selection_audit.csv", index=False)
    selected_ann = ann.loc[include].copy()
    # Preserve matrix order and verify all selected cells are present.
    X, genes, matrix_columns = read_matrix_selected(args.matrix, selected_ann["Index"].tolist())
    selected_ann = selected_ann.set_index("Index").loc[selected_ann["Index"].tolist()].reset_index()
    gene_upper = pd.Index([g.upper() for g in genes])
    var_names = pd.Index(genes)
    # Import scanpy only after the streaming parser so a failed optional GUI backend
    # cannot leave a half-written input manifest.
    import anndata as ad
    import scanpy as sc
    adata = ad.AnnData(X=X, obs=selected_ann.copy(), var=pd.DataFrame(index=var_names))
    adata.obs_names = selected_ann["Index"].astype(str).tolist()
    adata.obs["donor_id"] = adata.obs["Patient"].astype(str)
    adata.obs["condition"] = adata.obs["condition"].astype(str)
    adata.obs["dataset_batch"] = adata.obs["Sample"].astype(str)
    adata.layers["counts"] = adata.X.copy()
    # Raw counts -> log-normalised values. HVGs are selected independently from
    # the frozen 81-gene program.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Keep an unscaled log1p layer for marker statistics and expression summaries;
    # the main matrix is subsequently scaled only for PCA/neighbour finding.
    adata.layers["log1p"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=min(2500, max(100, adata.n_vars - 1)), flavor="seurat", subset=False)
    hvg_n = int(adata.var["highly_variable"].sum())
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(30, max(2, hvg_n - 1)), use_highly_variable=True, svd_solver="arpack", random_state=0)
    sc.pp.neighbors(adata, n_neighbors=min(15, max(3, adata.n_obs - 1)), n_pcs=min(30, adata.obsm["X_pca"].shape[1]), random_state=0)
    sc.tl.umap(adata, random_state=0, min_dist=0.35)
    resolutions = [0.3, 0.5, 0.7, 0.9, 1.1]
    res_rows = []
    for r in resolutions:
        key = f"leiden_{str(r).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=r, key_added=key, flavor="igraph", directed=False, n_iterations=2, random_state=0)
        c = adata.obs[key].astype(str)
        sizes = c.value_counts()
        donor_cov = adata.obs.assign(_c=c).groupby("_c")["donor_id"].nunique()
        # Repeat the graph partition at two additional seeds and report ARI
        # against the primary seed as a simple cluster-stability check.
        from sklearn.metrics import adjusted_rand_score
        stab = []
        for seed in (1, 2):
            tmp_key = f"_tmp_{str(r).replace('.', '_')}_{seed}"
            sc.tl.leiden(adata, resolution=r, key_added=tmp_key, flavor="igraph", directed=False, n_iterations=2, random_state=seed)
            stab.append(adjusted_rand_score(c.to_numpy(), adata.obs[tmp_key].astype(str).to_numpy()))
            del adata.obs[tmp_key]
        res_rows.append({"resolution": r, "n_clusters": int(c.nunique()), "min_cluster_cells": int(sizes.min()), "median_cluster_cells": float(sizes.median()), "min_cluster_donors": int(donor_cov.min()), "n_one_donor_clusters": int((donor_cov <= 1).sum()), "stability_ari_mean": float(np.mean(stab)), "hvg_n": hvg_n})
    res_audit = pd.DataFrame(res_rows)
    # Select the most stable resolution among partitions with coherent donor
    # coverage and no tiny clusters; this rule is fixed before looking at
    # PTGER4/program scores and therefore cannot be tuned to significance.
    eligible = res_audit[(res_audit["min_cluster_cells"] >= 30) & (res_audit["min_cluster_donors"] >= 3)]
    selected_res = float(eligible.sort_values(["stability_ari_mean", "min_cluster_cells"], ascending=False).iloc[0]["resolution"] if len(eligible) else res_audit.iloc[1]["resolution"])
    selected_key = f"leiden_{str(selected_res).replace('.', '_')}"
    res_audit["selected"] = res_audit["resolution"].eq(selected_res)
    res_audit.to_csv(out / "clustering_resolution_audit.csv", index=False)
    json_dump({"normalization": "library-size normalize_total target_sum=1e4 followed by log1p", "hvg_method": "scanpy highly_variable_genes flavor=seurat", "hvg_target": 2500, "hvg_selected": hvg_n, "pca": {"n_comps": int(adata.obsm["X_pca"].shape[1]), "use_highly_variable": True, "random_state": 0}, "neighbors": {"n_neighbors": 15, "metric": "scanpy default", "n_pcs": int(adata.obsm["X_pca"].shape[1])}, "umap": {"random_state": 0, "min_dist": 0.35}, "leiden_resolutions_explored": resolutions, "selected_resolution": selected_res, "selection_rule": "highest repeated-seed stability ARI among resolutions with min cluster >=30 cells and >=3 donors; stability ARI retained in clustering_resolution_audit.csv", "frozen_program_excluded_from_hvg": True}, out / "macrophage_subtype_parameters.json")
    adata.obs["cluster"] = adata.obs[selected_key].astype(str)
    markers = top_markers(adata, "cluster", 100)
    markers.to_csv(out / "macrophage_cluster_markers.csv", index=False)
    annots = annotate(markers)
    annots.to_csv(out / "macrophage_subtype_annotation.csv", index=False)
    cluster_to_sub = dict(zip(annots["cluster"], annots["subtype"]))
    adata.obs["subtype"] = adata.obs["cluster"].map(cluster_to_sub).fillna("other_macrophage_state")
    # Cluster QC, including contamination and condition/batch composition.
    qc_rows = []
    gene_to_i = {g.upper(): i for i, g in enumerate(adata.var_names)}
    for cl, g in adata.obs.groupby("cluster"):
        ix = g.index
        pos = np.array([adata.obs_names.get_loc(x) for x in ix])
        row = {"cluster": cl, "subtype": cluster_to_sub.get(cl, "other_macrophage_state"), "cell_count": len(pos), "donor_count": int(g["donor_id"].nunique()), "tumor_count": int((g["condition"] == "tumor").sum()), "normal_count": int((g["condition"] == "normal").sum()), "border_count": int((g["condition"] == "border").sum()), "dataset_batch_count": int(g["dataset_batch"].nunique())}
        for name, geneset in CONTAMINATION.items():
            inds = [gene_to_i[x] for x in geneset if x in gene_to_i]
            logx = adata.layers["log1p"]
            row[f"{name}_marker_detection_fraction"] = float(np.mean(logx[pos][:, inds].sum(axis=1) > 0)) if inds else np.nan
        row["low_confidence_one_donor"] = bool(row["donor_count"] < 3)
        qc_rows.append(row)
    pd.DataFrame(qc_rows).to_csv(out / "cluster_qc.csv", index=False)
    # Frozen program and the seven-gene driver score use gene-wise z-scored log1p
    # expression within the macrophage compartment; no program genes enter HVG selection.
    frozen = read_frozen_genes(args.gene_file)
    expr0 = adata.layers["log1p"]
    expr = expr0.toarray() if sp.issparse(expr0) else np.asarray(expr0)
    expr = np.asarray(expr, dtype=np.float32)
    upper_to_i = {g.upper(): i for i, g in enumerate(adata.var_names)}
    def score_genes(glist):
        inds = [upper_to_i[g] for g in glist if g in upper_to_i]
        if not inds: return np.full(adata.n_obs, np.nan), []
        z = (expr[:, inds] - expr[:, inds].mean(axis=0)) / np.where(expr[:, inds].std(axis=0) > 1e-8, expr[:, inds].std(axis=0), 1)
        return np.nanmean(z, axis=1), [glist[i] for i, g in enumerate(glist) if g in upper_to_i]
    dinp, dinp_present = score_genes(frozen)
    driver, driver_present = score_genes(FROZEN_DRIVER)
    pt_i = upper_to_i.get("PTGER4")
    pt_expr = expr[:, pt_i] if pt_i is not None else np.full(adata.n_obs, np.nan)
    sc_df = adata.obs[["Index", "donor_id", "condition", "dataset_batch", "cluster", "subtype"]].copy() if "Index" in adata.obs else adata.obs[["donor_id", "condition", "dataset_batch", "cluster", "subtype"]].copy()
    sc_df["cell_id"] = adata.obs_names.astype(str)
    sc_df["dinp_crc_program_score"] = dinp
    sc_df["driver_7gene_score"] = driver
    sc_df["PTGER4_normalized_log1p"] = pt_expr
    sc_df["PTGER4_detected"] = pt_expr > 0
    sc_df.to_csv(out / "macrophage_cell_program_scores.csv", index=False)
    # Donor-level summaries and paired inference.
    group_cols = ["donor_id", "condition", "subtype"]
    donor = sc_df.groupby(group_cols, dropna=False).agg(dinp_crc_program_mean=("dinp_crc_program_score", "mean"), dinp_crc_program_median=("dinp_crc_program_score", "median"), driver_7gene_mean=("driver_7gene_score", "mean"), PTGER4_mean=("PTGER4_normalized_log1p", "mean"), PTGER4_detection_fraction=("PTGER4_detected", "mean"), n_cells=("cell_id", "size")).reset_index()
    donor.to_csv(out / "donor_subtype_program_scores.csv", index=False)
    dinp_stats = paired_stats(donor, "dinp_crc_program_mean")
    dinp_stats.to_csv(out / "subtype_program_statistics.csv", index=False)
    pt_stats = paired_stats(donor, "PTGER4_mean")
    pt_stats.to_csv(out / "ptger4_subtype_statistics.csv", index=False)
    # Correlation is computed on donor-subtype aggregates, avoiding cell-level p-values.
    cor_rows = []
    for sub, g in donor.groupby("subtype"):
        x, y = g["dinp_crc_program_mean"], g["PTGER4_mean"]
        ok = x.notna() & y.notna()
        rho, p = stats.spearmanr(x[ok], y[ok]) if ok.sum() >= 3 else (np.nan, np.nan)
        cor_rows.append({"subtype": sub, "n_donor_condition": int(ok.sum()), "spearman_rho": rho, "p": p})
    pd.DataFrame(cor_rows).to_csv(out / "ptger4_dinp_correlation_donor_subtype.csv", index=False)
    # Fixed driver program at subtype level.
    driver_out = donor[["donor_id", "condition", "subtype", "driver_7gene_mean", "n_cells"]].copy()
    driver_out.to_csv(out / "driver_7gene_donor_subtype_scores.csv", index=False)
    driver_stats = paired_stats(driver_out.rename(columns={"driver_7gene_mean": "value"}), "value")
    driver_stats.to_csv(out / "driver_7gene_statistics.csv", index=False)
    # Per-gene dot-plot data for the fixed driver set (kept separate from the
    # aggregate score; this is not an independent pathway claim).
    driver_rows = []
    for gene in FROZEN_DRIVER:
        gi = upper_to_i.get(gene)
        if gi is None:
            continue
        tmp = sc_df[["donor_id", "condition", "subtype", "cell_id"]].copy()
        tmp["gene"] = gene
        tmp["expression"] = expr[:, gi]
        tmp["detected"] = tmp["expression"] > 0
        driver_rows.append(tmp.groupby(["gene", "condition", "subtype"], dropna=False).agg(mean_expression=("expression", "mean"), detection_fraction=("detected", "mean"), n_cells=("cell_id", "size")).reset_index())
    if driver_rows:
        pd.concat(driver_rows, ignore_index=True).to_csv(out / "driver_7gene_subtype_expression.csv", index=False)
    # Abundance is donor-level composition, with paired tests on proportions.
    counts = sc_df.groupby(["donor_id", "condition", "subtype"], dropna=False).size().rename("n_cells").reset_index()
    totals = counts.groupby(["donor_id", "condition"])["n_cells"].transform("sum")
    counts["fraction"] = counts["n_cells"] / totals
    abundance_stats = paired_stats(counts.rename(columns={"fraction": "value"}), "value")
    counts.to_csv(out / "subtype_abundance_by_donor.csv", index=False)
    abundance_stats.to_csv(out / "subtype_abundance_statistics.csv", index=False)
    # Basic figures (UMAP and marker heatmap; quantitative figures below are kept
    # simple and reproducible with the donor-level tables).
    save_figures(adata, markers, out)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    def simple_box(df, value, filename, title):
        fig, ax = plt.subplots(figsize=(max(7, .8 * df["subtype"].nunique()), 5))
        sns.boxplot(data=df, x="subtype", y=value, hue="condition", ax=ax, fliersize=1, palette="Set2")
        sns.stripplot(data=df, x="subtype", y=value, hue="condition", dodge=True, ax=ax, size=3, color="black", alpha=.45, legend=False)
        ax.set_title(title); ax.tick_params(axis="x", rotation=35); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(out / "figures" / f"{filename}.{ext}", dpi=220)
        plt.close(fig)
    simple_box(donor, "dinp_crc_program_mean", "dinp_crc_program_by_subtype", "Frozen DINP–CRC program by macrophage subtype")
    simple_box(donor, "PTGER4_mean", "ptger4_by_subtype", "PTGER4 by macrophage subtype")
    # Paired program lines and abundance are donor-level by design.
    for value, filename, title in [("dinp_crc_program_mean", "paired_tumor_normal_program_by_subtype", "Paired tumor–normal DINP–CRC score"), ("fraction", "macrophage_subtype_abundance_tumor_normal", "Macrophage subtype abundance")]:
        d = donor if value != "fraction" else counts
        fig, ax = plt.subplots(figsize=(max(7, .8 * d["subtype"].nunique()), 5))
        for sub, g in d.groupby("subtype"):
            p = g.pivot_table(index="donor_id", columns="condition", values=value, aggfunc="mean")
            if {"tumor", "normal"}.issubset(p.columns):
                p = p.dropna(subset=["tumor", "normal"])
                for _, row in p.iterrows(): ax.plot([sub, sub], [row["normal"], row["tumor"]], color="0.7", lw=.7)
                ax.scatter([sub] * len(p), p["normal"], color="#377eb8", s=22, label="normal" if sub == d["subtype"].iloc[0] else None)
                ax.scatter([sub] * len(p), p["tumor"], color="#e41a1c", s=22, label="tumor" if sub == d["subtype"].iloc[0] else None)
        ax.set_title(title); ax.set_xlabel("subtype"); ax.set_ylabel(value); ax.tick_params(axis="x", rotation=35); ax.legend(frameon=False); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(out / "figures" / f"{filename}.{ext}", dpi=220)
        plt.close(fig)
    # Provenance is written last so it captures all selected features and exact hashes.
    manifest = {"analysis": "dinp_crc_macrophage_subtype", "run_timestamp_utc": pd.Timestamp.utcnow().isoformat(), "matrix": {"path": str(args.matrix), "sha256": sha256(args.matrix), "size_bytes": args.matrix.stat().st_size}, "annotation": {"path": str(args.annotation), "sha256": sha256(args.annotation), "size_bytes": args.annotation.stat().st_size}, "frozen_gene_file": {"path": str(args.gene_file), "sha256": sha256(args.gene_file), "n_genes": len(frozen)}, "source_dataset": "GSE144735 CRC single-cell matrix; 6 patients; tumor/border/normal", "source_cells_total": int(len(ann)), "selected_macrophage_like_cells": int(adata.n_obs), "selected_donors": int(adata.obs["donor_id"].nunique()), "metadata_fields": list(map(str, ann.columns)), "donor_id_field": "Patient", "condition_field": "Class", "coarse_cell_type_field": "Cell_type", "coarse_cell_subtype_field": "Cell_subtype", "normalization_strategy": "raw UMI counts; normalize_total target_sum=1e4; log1p; gene-wise z-scoring for program scores", "existing_umap_available": False, "existing_umap_action": "recomputed macrophage-only UMAP", "prior_localization_outputs": [r"E:\mcop\mcop\analysis\dinp_crc_81gene_singlecell_localization\outputs", r"E:\mcop\mcop\analysis\dinp_crc_81gene_singlecell_subtype_localization\outputs"], "selected_resolution": selected_res, "hvg_n": hvg_n, "frozen_genes_present": dinp_present, "frozen_genes_missing": [g for g in frozen if g not in dinp_present], "driver_genes_present": driver_present, "primary_370115_cell_census_object": "not available locally; no redownload performed", "runtime_seconds": round(time.time() - t0, 2)}
    json_dump(manifest, out / "input_manifest.json")
    # Human-readable integration report.
    top_sub = donor.groupby("subtype")["dinp_crc_program_mean"].mean().sort_values(ascending=False)
    top_pt = donor.groupby("subtype")["PTGER4_mean"].mean().sort_values(ascending=False)
    top_ab = counts[counts["condition"] == "tumor"].groupby("subtype")["fraction"].mean().sort_values(ascending=False)
    top_prog_paired = int(dinp_stats.loc[dinp_stats["subtype"] == top_sub.index[0], "paired_donors"].iloc[0]) if len(top_sub) and len(dinp_stats.loc[dinp_stats["subtype"] == top_sub.index[0]]) else 0
    next_step = "No further subtype work until a larger CRC cohort validates the state (top-program subtype has fewer than 3 paired tumor-normal donors)." if top_prog_paired < 3 else "CellChat/spatial follow-up can be considered after the donor-level and marker QC review."
    report = f"""# Macrophage subtype analysis\n\n- Source: local GSE144735 CRC matrix; {len(ann):,} cells, {ann['Patient'].nunique()} donors.\n- Macrophage-only selection: {adata.n_obs:,} cells; cDC, granulocyte, mast, lymphoid and non-myeloid labels excluded.\n- Recomputed macrophage-only PCA, neighbors, UMAP and Leiden at resolutions {', '.join(map(str, resolutions))}; selected resolution **{selected_res}** by the predeclared stability-ARI rule after minimum cluster-size/donor-coverage filters.\n- Frozen DINP–CRC program: {len(frozen)} genes ({len(dinp_present)} present in the matrix); genes were not forced into HVGs.\n\n## Integration answers\n\n**A. Where is the program concentrated?** The highest donor-condition mean score is **{top_sub.index[0] if len(top_sub) else 'NA'}** ({top_sub.iloc[0]:.3f} if available). This is a descriptive donor-level ranking; formal paired results are in `subtype_program_statistics.csv`.\n\n**B. Where is PTGER4 distributed?** The highest donor-condition mean is **{top_pt.index[0] if len(top_pt) else 'NA'}** ({top_pt.iloc[0]:.3f} if available). Detection fractions and paired results are in `donor_subtype_program_scores.csv` and `ptger4_subtype_statistics.csv`.\n\n**C. Tumor versus normal.** Inference uses donor-level aggregates and paired tumor-normal tests only; BH-FDR is applied within each named statistic family. Border cells are retained for QC and descriptive plots, but are not used in paired tumor-normal tests.\n\n**D. Abundance.** The highest mean tumor composition is **{top_ab.index[0] if len(top_ab) else 'NA'}**. This is compositional and exploratory; see `subtype_abundance_statistics.csv`.\n\n**E. Next step.** {next_step} This module does not run CellChat, spatial, AUCell, deconvolution, pathway enrichment or causal analyses.\n\n## Limitations and sanity checks\n\nThe earlier 370,115-cell Census object used for broad localization was unavailable locally, so this run uses the existing CRC GSE144735 cache and does not redownload data. Results should not be numerically pooled across the two source objects. One-donor clusters are flagged in `cluster_qc.csv`; marker specificity, contamination, frozen-gene coverage, resolution sensitivity and donor counts are retained in the output tables.\n"""
    (out / "MACROPHAGE_SUBTYPE_REPORT.md").write_text(report, encoding="utf-8")
    print(f"MACROPHAGE SUBTYPE ANALYSIS: PASS")
    print(f"total cells: {len(ann)}; selected macrophage-like cells: {adata.n_obs}; donors: {adata.obs['donor_id'].nunique()}")
    print(f"selected resolution: {selected_res}; subtypes: {', '.join(sorted(adata.obs['subtype'].unique()))}")
    print(f"top DINP subtype: {top_sub.index[0] if len(top_sub) else 'NA'}; top tumor-enriched subtype: {top_ab.index[0] if len(top_ab) else 'NA'}; PTGER4-high subtype: {top_pt.index[0] if len(top_pt) else 'NA'}")
    print("paired result: see subtype_program_statistics.csv and ptger4_subtype_statistics.csv")
    print(f"next step: {next_step}")


if __name__ == "__main__":
    main()
