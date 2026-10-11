#!/usr/bin/env python
"""Prepare de-novo CRC single-cell annotations and CellChat inputs.

Only reads the original GEO count matrices and their original sample/cell
annotations. It never loads the project's prior derived matrices or results.
The full transcriptome is used to recluster myeloid and stromal candidates;
the epithelial cells are retained with the full CellChat ligand/receptor
panel. SCD and PID1 are excluded from variable-feature clustering and only
examined after cell labels are assigned.
"""

from __future__ import annotations

import gzip
import gc
import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite


ROOT = Path("D:/DINP/whynot17-new-project")
RAW = ROOT / "work/single_cell_GSE132465_GSE200997"
OUT = ROOT / "results/2026-10-11/cellchat_patient_consistency"
DB_GENES = ROOT / "work/cellchat_db_gene_symbols.txt"
OUT.mkdir(parents=True, exist_ok=True)
sc.settings.verbosity = 2
sc.settings.n_jobs = 1
sc.settings.figdir = str(OUT)


MYELOID_MODULES = {
    "Monocyte": ["FCN1", "S100A8", "S100A9", "VCAN", "CCR2", "CD14", "S100A12", "LILRB1", "CTSS"],
    "Macrophage": ["C1QA", "C1QB", "C1QC", "APOE", "APOC1", "CSF1R", "CTSD", "FCGR1A", "MSR1", "TREM2", "SPP1", "LPL", "LGALS3"],
    "DC_control": ["FCER1A", "CD1C", "CLEC10A", "XCR1", "CLEC9A", "IL3RA", "CLEC4C"],
}
CAF_MODULES = {
    "Fibroblast_core": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRA", "COL6A1", "COL14A1", "PI16"],
    "myCAF": ["ACTA2", "TAGLN", "MYL9", "POSTN", "CNN1", "TPM2", "MYLK"],
    "iCAF": ["CXCL12", "IL6", "CXCL14", "CCL2", "CFD", "PTGDS", "PI16", "IL1RL1", "CXCL1", "CXCL8"],
    "ECM_CAF": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "FN1", "TNC", "FAP", "THBS2", "POSTN"],
    "apCAF": ["CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DRB1", "CIITA"],
    "Endothelial_control": ["PECAM1", "VWF", "KDR", "EMCN", "ESAM", "CLDN5"],
    "Pericyte_control": ["RGS5", "CSPG4", "NOTCH3", "PDGFRB", "MCAM", "DES"],
    "Smooth_muscle_control": ["MYH11", "CNN1", "ACTG2", "MYLK", "TAGLN"],
}
BROAD_MODULES = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CDH1", "CEACAM5", "MUC1", "MUC2", "TFF3", "KRT7"],
    "Myeloid": ["PTPRC", "LST1", "TYROBP", "LYZ", "FCER1G", "CTSS", "AIF1", "CD68", "CSF1R", "FCGR3A"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRA", "COL6A1", "COL14A1"],
    "T_NK": ["CD3D", "CD3E", "TRBC1", "TRBC2", "CD247", "LCK", "NKG7", "GNLY", "CD8A"],
    "B_Plasma": ["CD79A", "MS4A1", "CD37", "CD19", "JCHAIN", "MZB1", "IGKC"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ESAM", "CLDN5"],
    "Pericyte": ["RGS5", "CSPG4", "NOTCH3", "PDGFRB", "MCAM", "DES"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
}
LIPID_MODULE_NO_SCD = [
    "ACLY", "ACACA", "ACACB", "FASN", "SREBF1", "SREBF2", "HMGCR", "HMGCS1", "SQLE", "LDLR",
    "ACSS2", "ACSS1", "ACSL1", "ACSL3", "ACSL4", "ACSL5", "FADS1", "FADS2", "ELOVL5", "ELOVL6",
    "CPT1A", "CPT2", "ACOX1", "DGAT1", "DGAT2", "PLIN2", "PLIN3", "GPX4", "LPCAT3", "FABP1",
    "APOE", "LPL", "PPARA", "PPARG", "IDI1", "FDPS", "FDFT1", "DHCR7", "DHCR24", "ACAT2",
]
PID1_SCD = ["PID1", "SCD"]
ALL_MODULES = {**MYELOID_MODULES, **CAF_MODULES, **BROAD_MODULES, "Lipid_program_excluding_SCD": LIPID_MODULE_NO_SCD}


def log(message: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {message}", flush=True)


def load_panel() -> set[str]:
    if not DB_GENES.exists():
        raise FileNotFoundError(f"Run {ROOT / 'analysis/cellchat_export_db_genes.R'} first")
    genes = set(x.strip() for x in DB_GENES.read_text(encoding="utf-8").splitlines() if x.strip())
    genes.update(PID1_SCD)
    for module in ALL_MODULES.values():
        genes.update(module)
    return genes


def make_sparse_gene_by_cell(
    data_parts: list[np.ndarray], index_parts: list[np.ndarray], indptr: list[int],
    n_cell: int, gene_names: list[str],
) -> sparse.csr_matrix:
    if data_parts:
        data = np.concatenate(data_parts).astype(np.int32, copy=False)
        indices = np.concatenate(index_parts).astype(np.int32, copy=False)
    else:
        data = np.array([], dtype=np.int32)
        indices = np.array([], dtype=np.int32)
    mat = sparse.csr_matrix((data, indices, np.asarray(indptr, dtype=np.int64)),
                            shape=(len(gene_names), n_cell), dtype=np.int32)
    mat.sum_duplicates()
    mat.eliminate_zeros()
    return mat


def parse_matrix_header(matrix_path: Path, delimiter: str) -> list[str]:
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        line = fh.readline().rstrip("\r\n")
    if delimiter == "\t":
        fields = line.split("\t")
        return fields[1:]
    fields = next(__import__("csv").reader([line]))
    return fields[1:]


def parse_gene_line(line: str, delimiter: str) -> tuple[str, np.ndarray]:
    line = line.rstrip("\r\n")
    if delimiter == "\t":
        gene, values = line.split("\t", 1)
        return gene, np.fromstring(values, sep="\t", dtype=np.int32)
    comma = line.find(",")
    gene = line[1:comma - 1] if line.startswith('"') else line[:comma]
    values = line[comma + 1:]
    return gene, np.fromstring(values, sep=",", dtype=np.int32)


def module_scores(X: sparse.csr_matrix, genes: list[str], modules: dict[str, list[str]]) -> pd.DataFrame:
    lookup = {g: i for i, g in enumerate(genes)}
    out: dict[str, np.ndarray] = {}
    for name, members in modules.items():
        cols = [lookup[g] for g in members if g in lookup]
        if not cols:
            out[name] = np.zeros(X.shape[0], dtype=np.float32)
        else:
            out[name] = np.asarray(X[:, cols].mean(axis=1)).ravel().astype(np.float32)
    return pd.DataFrame(out)


def add_expr_scores(meta: pd.DataFrame, expr: sparse.csr_matrix, genes: list[str], dataset: str) -> None:
    lookup = {g: i for i, g in enumerate(genes)}
    for target in ["SCD", "PID1"]:
        if target in lookup:
            values = np.asarray(expr[:, lookup[target]].toarray()).ravel()
            meta[f"{target}_log1p_norm"] = values
            meta[f"{target}_detected"] = values > 0
        else:
            meta[f"{target}_log1p_norm"] = 0.0
            meta[f"{target}_detected"] = False
    cols = [lookup[g] for g in LIPID_MODULE_NO_SCD if g in lookup]
    meta["lipid_program_score_excluding_SCD"] = (
        np.asarray(expr[:, cols].mean(axis=1)).ravel() if cols else 0.0
    )
    meta["lipid_program_genes_used"] = ",".join(g for g in LIPID_MODULE_NO_SCD if g in lookup)
    meta["cohort"] = dataset


def summarize_targets(meta: pd.DataFrame, outfile: Path) -> None:
    keys = ["cohort", "sample_id", "celltype"]
    cols = ["SCD_log1p_norm", "SCD_detected", "PID1_log1p_norm", "PID1_detected", "lipid_program_score_excluding_SCD"]
    out = meta.groupby(keys, observed=True).agg(
        n_cells=("celltype", "size"),
        SCD_mean_log1p=("SCD_log1p_norm", "mean"), SCD_pct_detected=("SCD_detected", "mean"),
        PID1_mean_log1p=("PID1_log1p_norm", "mean"), PID1_pct_detected=("PID1_detected", "mean"),
        lipid_score_mean=("lipid_program_score_excluding_SCD", "mean"),
        lipid_score_median=("lipid_program_score_excluding_SCD", "median"),
    ).reset_index()
    out["SCD_pct_detected"] *= 100.0
    out["PID1_pct_detected"] *= 100.0
    out.to_csv(outfile, sep="\t", index=False)


def write_matrix_market(expr_cells_by_genes: sparse.csr_matrix, genes: list[str], meta: pd.DataFrame, stem: str) -> None:
    matrix_path = OUT / f"{stem}_lognorm_genes_by_cells.mtx.gz"
    genes_path = OUT / f"{stem}_features.tsv"
    cells_path = OUT / f"{stem}_cells.tsv"
    genes_path.write_text("\n".join(genes) + "\n", encoding="utf-8")
    cells_path.write_text("\n".join(meta["cell"].astype(str)) + "\n", encoding="utf-8")
    tmp = OUT / f"{stem}_matrix.tmp.mtx"
    mmwrite(str(tmp), expr_cells_by_genes.T.tocoo().astype(np.float32), field="real", precision=6)
    with open(tmp, "rb") as source, gzip.open(matrix_path, "wb", compresslevel=6) as target:
        while True:
            chunk = source.read(8 * 1024 * 1024)
            if not chunk:
                break
            target.write(chunk)
    tmp.unlink(missing_ok=True)
    log(f"Wrote CellChat matrix {matrix_path.name}: {expr_cells_by_genes.shape[0]} cells x {expr_cells_by_genes.shape[1]} genes; nnz={expr_cells_by_genes.nnz}")


def plot_umap(adata: ad.AnnData, stem: str, color: str = "subcluster") -> None:
    if "X_umap" not in adata.obsm:
        return
    coords = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(7, 6))
    labels = adata.obs[color].astype(str)
    cats = sorted(labels.unique())
    cmap = plt.get_cmap("tab10")
    for i, label in enumerate(cats):
        mask = labels.to_numpy() == label
        ax.scatter(coords[mask, 0], coords[mask, 1], s=4, alpha=0.65, label=label, color=cmap(i % 10), linewidths=0)
    ax.set(title=stem.replace("_", " "), xlabel="UMAP1", ylabel="UMAP2")
    ax.legend(frameon=False, markerscale=2, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}_umap.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def recluster_group(
    counts_cells_by_genes: sparse.csr_matrix, genes: list[str], obs: pd.DataFrame,
    axis: str, cohort: str,
) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Recluster one parent lineage on full transcriptome, then annotate Leiden clusters by marker programs."""
    if axis == "Myeloid":
        modules = MYELOID_MODULES
        output_names = ["Monocyte", "Macrophage", "Other_myeloid"]
    else:
        modules = CAF_MODULES
        output_names = ["myCAF", "iCAF", "ECM_CAF", "apCAF", "Stromal_other"]
    tmp = ad.AnnData(X=counts_cells_by_genes.copy(), obs=obs.copy(), var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    tmp.var_names_make_unique()
    sc.pp.filter_genes(tmp, min_cells=3)
    sc.pp.normalize_total(tmp, target_sum=10000)
    sc.pp.log1p(tmp)
    scores = module_scores(tmp.X.tocsr(), tmp.var_names.tolist(), modules)
    for col in scores.columns:
        tmp.obs[f"score_{col}"] = scores[col].to_numpy()
    if "SCD" in tmp.var_names:
        tmp = tmp[:, tmp.var_names != "SCD"].copy()
    if "PID1" in tmp.var_names:
        tmp = tmp[:, tmp.var_names != "PID1"].copy()
    n_top = min(2000, tmp.n_vars)
    sc.pp.highly_variable_genes(tmp, n_top_genes=n_top, flavor="seurat", subset=False)
    keep = tmp.var["highly_variable"].to_numpy()
    keep &= ~tmp.var_names.isin(["SCD", "PID1"])
    if keep.sum() < 20 or tmp.n_obs < 50:
        tmp.obs["leiden"] = "0"
        tmp.obsm["X_umap"] = np.zeros((tmp.n_obs, 2), dtype=np.float32)
    else:
        cluster = tmp[:, keep].copy()
        n_pcs = max(2, min(30, cluster.n_vars - 1, cluster.n_obs - 1))
        sc.pp.scale(cluster, max_value=10)
        sc.tl.pca(cluster, n_comps=n_pcs, svd_solver="randomized", random_state=17)
        sc.pp.neighbors(cluster, n_neighbors=min(15, cluster.n_obs - 1), n_pcs=n_pcs, random_state=17)
        sc.tl.leiden(cluster, resolution=0.5, random_state=17, flavor="igraph", directed=False, n_iterations=2)
        sc.tl.umap(cluster, random_state=17)
        tmp.obs["leiden"] = cluster.obs["leiden"].astype(str).reindex(tmp.obs_names).to_numpy()
        tmp.obsm["X_umap"] = cluster.obsm["X_umap"].copy()
        # Keep the input matrix sparse and unscaled for interpretable marker ranking.
        markers = tmp[:, keep].copy()
        markers.obs["leiden"] = tmp.obs["leiden"].astype(str)
        if markers.obs["leiden"].nunique() > 1:
            sc.tl.rank_genes_groups(markers, groupby="leiden", method="wilcoxon", use_raw=False, pts=True)
            ranked = sc.get.rank_genes_groups_df(markers, group=None)
            ranked = ranked.sort_values(["group", "pvals_adj", "scores"], ascending=[True, True, False])
            ranked = ranked.groupby("group", observed=True).head(15)
            ranked.insert(0, "axis", axis)
            ranked.insert(0, "cohort", cohort)
            markers_out = ranked
        else:
            markers_out = pd.DataFrame()
        del cluster, markers
    score_columns = [f"score_{x}" for x in modules]
    by_cluster = tmp.obs.groupby("leiden", observed=True)[score_columns].mean().reset_index()
    by_cluster.insert(0, "axis", axis)
    by_cluster.insert(0, "cohort", cohort)
    counts_by_cluster = tmp.obs.groupby("leiden", observed=True).size().rename("n_cells").reset_index()
    by_cluster = by_cluster.merge(counts_by_cluster, on="leiden", how="left")

    label_map: dict[str, str] = {}
    if axis == "Myeloid":
        for _, row in by_cluster.iterrows():
            vals = {k: float(row[f"score_{k}"]) for k in MYELOID_MODULES}
            best = max(vals, key=vals.get)
            label_map[str(row["leiden"])] = best if best in ("Monocyte", "Macrophage") else "Other_myeloid"
    else:
        states = ["myCAF", "iCAF", "ECM_CAF", "apCAF"]
        controls = ["Endothelial_control", "Pericyte_control", "Smooth_muscle_control"]
        for _, row in by_cluster.iterrows():
            fibro_core = float(row["score_Fibroblast_core"])
            control = max(float(row[f"score_{x}"]) for x in controls)
            if fibro_core <= control:
                label = "Stromal_other"
            else:
                label = max(states, key=lambda x: float(row[f"score_{x}"]))
            label_map[str(row["leiden"])] = label
    tmp.obs["subcluster"] = tmp.obs["leiden"].astype(str).map(label_map).fillna(output_names[-1])
    assignment = tmp.obs[["subcluster", "leiden"]].copy()
    assignment["axis"] = axis
    assignment.index.name = "cell"
    assignment = assignment.reset_index()
    by_cluster["assigned_subcluster"] = by_cluster["leiden"].astype(str).map(label_map)
    cluster_table = by_cluster.rename(columns={"leiden": "cluster"})
    plot_umap(tmp, f"{cohort}_{axis.lower()}")
    label_series = tmp.obs["subcluster"].copy()
    label_series.index.name = "cell"
    if "markers_out" not in locals():
        markers_out = pd.DataFrame()
    del tmp
    gc.collect()
    return label_series, cluster_table, markers_out, assignment


def scan_gse132465(panel_genes: set[str]) -> None:
    log("GSE132465: reading original annotation and streaming the original count matrix")
    matrix_path = RAW / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
    ann_path = RAW / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t", dtype=str)
    ann = ann[ann["Class"].eq("Tumor") & ann["Cell_type"].isin(["Epithelial cells", "Myeloids", "Stromal cells"])].copy()
    ann = ann.rename(columns={"Index": "cell", "Patient": "sample_id", "Cell_type": "original_type", "Cell_subtype": "original_subtype"})
    source_cells = parse_matrix_header(matrix_path, "\t")
    ann_idx = ann.set_index("cell", drop=False)
    keep_set = set(ann["cell"])
    selected_cells = [c for c in source_cells if c in keep_set]
    if len(selected_cells) != len(ann):
        missing = len(ann) - len(selected_cells)
        raise ValueError(f"GSE132465 matrix/annotation cell mismatch: {missing} annotation cells absent from matrix")
    meta = ann_idx.loc[selected_cells].reset_index(drop=True)
    cluster_local = np.flatnonzero(meta["original_type"].isin(["Myeloids", "Stromal cells"]).to_numpy())
    epi_local = np.flatnonzero(meta["original_type"].eq("Epithelial cells").to_numpy())
    source_index = {c: i for i, c in enumerate(source_cells)}
    selected_source = np.asarray([source_index[c] for c in selected_cells], dtype=np.int32)
    cluster_source = selected_source[cluster_local]
    n_source = len(source_cells)
    log(f"GSE132465 selected tumor cells: {len(selected_cells)}; full-transcriptome reclustering set: {len(cluster_local)}; epithelial: {len(epi_local)}")

    cl_data: list[np.ndarray] = []; cl_idx: list[np.ndarray] = []; cl_ptr = [0]; cl_genes: list[str] = []
    pn_data: list[np.ndarray] = []; pn_idx: list[np.ndarray] = []; pn_ptr = [0]; pn_genes: list[str] = []
    libsize = np.zeros(len(selected_cells), dtype=np.int64)
    t0 = time.time()
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        fh.readline()
        for i, line in enumerate(fh):
            gene, vals = parse_gene_line(line, "\t")
            if vals.size != n_source:
                raise ValueError(f"GSE132465 row {i + 1} has {vals.size} values; expected {n_source}")
            target_vals = vals[selected_source]
            libsize += target_vals
            gv = vals[cluster_source]
            nz = np.flatnonzero(gv)
            if nz.size:
                cl_genes_data = gv[nz]
                cl_data.append(cl_genes_data)
                cl_idx.append(nz.astype(np.int32, copy=False))
                cl_ptr.append(cl_ptr[-1] + nz.size)
            else:
                cl_ptr.append(cl_ptr[-1])
            cl_genes.append(gene)
            if gene in panel_genes:
                pn_genes.append(gene)
                pnz = np.flatnonzero(target_vals)
                if pnz.size:
                    pn_data.append(target_vals[pnz])
                    pn_idx.append(pnz.astype(np.int32, copy=False))
                    pn_ptr.append(pn_ptr[-1] + pnz.size)
                else:
                    pn_ptr.append(pn_ptr[-1])
            if (i + 1) % 5000 == 0:
                log(f"GSE132465 streamed {i + 1:,} genes")
    cl_gxc = make_sparse_gene_by_cell(cl_data, cl_idx, cl_ptr, len(cluster_local), cl_genes)
    panel_gxc = make_sparse_gene_by_cell(pn_data, pn_idx, pn_ptr, len(selected_cells), pn_genes)
    log(f"GSE132465 parsed in {time.time() - t0:.1f}s; full subset nnz={cl_gxc.nnz:,}; panel nnz={panel_gxc.nnz:,}")

    cluster_obs = meta.iloc[cluster_local].set_index("cell").copy()
    cluster_obs["broad_type"] = np.where(cluster_obs["original_type"].eq("Myeloids"), "Myeloid", "Stromal")
    cl_cxg = cl_gxc.T.tocsr()
    del cl_gxc, cl_data, cl_idx, pn_data, pn_idx
    gc.collect()
    all_labels = pd.Series(index=cluster_obs.index, dtype="object")
    cluster_tables: list[pd.DataFrame] = []
    marker_tables: list[pd.DataFrame] = []
    assignment_tables: list[pd.DataFrame] = []
    cluster_genes = cl_genes
    for axis in ["Myeloid", "Stromal"]:
        rows = np.flatnonzero(cluster_obs["broad_type"].eq(axis).to_numpy())
        if not len(rows):
            continue
        labels, table, markers, assignment = recluster_group(cl_cxg[rows].tocsr(), cluster_genes, cluster_obs.iloc[rows], axis, "GSE132465")
        all_labels.loc[labels.index] = labels
        cluster_tables.append(table); marker_tables.append(markers); assignment_tables.append(assignment)
    del cl_cxg
    gc.collect()
    meta = meta.set_index("cell", drop=False)
    meta["celltype"] = "Epithelial"
    meta.loc[all_labels.index, "celltype"] = all_labels
    meta["cohort"] = "GSE132465"
    meta["patient_id"] = meta["sample_id"]
    meta["cellchat_include"] = meta["celltype"].isin(["Epithelial", "Monocyte", "Macrophage", "Other_myeloid", "myCAF", "iCAF", "ECM_CAF", "apCAF"])
    meta_out = meta.reset_index(drop=True)
    meta_out.to_csv(OUT / "GSE132465_cell_annotations.tsv.gz", sep="\t", index=False, compression="gzip")
    if cluster_tables:
        pd.concat(cluster_tables, ignore_index=True).to_csv(OUT / "GSE132465_cluster_program_scores.tsv", sep="\t", index=False)
    if marker_tables:
        pd.concat([x for x in marker_tables if not x.empty], ignore_index=True).to_csv(OUT / "GSE132465_cluster_markers.tsv.gz", sep="\t", index=False, compression="gzip")
    if assignment_tables:
        pd.concat(assignment_tables, ignore_index=True).to_csv(OUT / "GSE132465_cluster_assignments.tsv.gz", sep="\t", index=False, compression="gzip")
    counts = meta_out.groupby(["sample_id", "original_type", "celltype"], observed=True).size().rename("n_cells").reset_index()
    counts.to_csv(OUT / "GSE132465_cell_counts_by_patient.tsv", sep="\t", index=False)
    log("GSE132465 final annotated groups:\n" + meta_out["celltype"].value_counts().to_string())

    panel_cells_by_gene = panel_gxc.T.tocsr()
    # Panel rows represent target genes, but normalization uses complete-cell UMI totals.
    scale = np.divide(10000.0, libsize, out=np.zeros_like(libsize, dtype=np.float32), where=libsize > 0)
    panel_cells_by_gene = panel_cells_by_gene.multiply(scale[:, None]).tocsr().astype(np.float32)
    panel_cells_by_gene.data = np.log1p(panel_cells_by_gene.data)
    panel_cells_by_gene.eliminate_zeros()
    pmeta = meta_out.set_index("cell").loc[selected_cells].reset_index()
    pmeta["cell"] = selected_cells
    add_expr_scores(pmeta, panel_cells_by_gene, pn_genes, "GSE132465")
    pmeta["cellchat_include"] = pmeta["celltype"].isin(["Epithelial", "Monocyte", "Macrophage", "Other_myeloid", "myCAF", "iCAF", "ECM_CAF", "apCAF"])
    include = pmeta["cellchat_include"].to_numpy()
    out_expr = panel_cells_by_gene[include].tocsr()
    out_meta = pmeta.loc[include].reset_index(drop=True)
    write_matrix_market(out_expr, pn_genes, out_meta, "GSE132465")
    summarize_targets(out_meta, OUT / "GSE132465_target_gene_patient_summary.tsv")
    meta_out.to_csv(OUT / "GSE132465_cell_annotations.tsv.gz", sep="\t", index=False, compression="gzip")
    del panel_gxc, panel_cells_by_gene, out_expr
    gc.collect()


def classify_gse200997_broad(expr: sparse.csr_matrix, genes: list[str], meta: pd.DataFrame) -> pd.DataFrame:
    scores = module_scores(expr, genes, BROAD_MODULES)
    lookup = {g: i for i, g in enumerate(genes)}
    counts: dict[str, np.ndarray] = {}
    for name, members in BROAD_MODULES.items():
        cols = [lookup[g] for g in members if g in lookup]
        counts[name] = np.asarray((expr[:, cols] > 0).sum(axis=1)).ravel() if cols else np.zeros(expr.shape[0], dtype=np.int32)
    meta = meta.copy()
    for c in scores:
        meta[f"score_{c}"] = scores[c].to_numpy()
        meta[f"nmark_{c}"] = counts[c]
    # Compete target lineages against the non-target controls too. This prevents
    # T/NK or B cells with a weak myeloid score from being forced into Myeloid.
    all_types = list(BROAD_MODULES)
    all_scores = scores[all_types].to_numpy()
    best_idx = all_scores.argmax(axis=1)
    best_name = np.asarray(all_types)[best_idx]
    best_score = all_scores.max(axis=1)
    second_score = np.sort(all_scores, axis=1)[:, -2]
    min_mark = np.choose(best_idx, [counts[name] for name in all_types])
    allowed = np.isin(best_name, ["Epithelial", "Myeloid", "Fibroblast"])
    # Resolve broad lineage from multi-marker programs; ambiguous and non-target cells remain unassigned.
    assigned = allowed & (best_score >= 0.15) & (min_mark >= 2) & ((best_score - second_score) >= 0.05)
    meta["broad_top_marker_class"] = best_name
    meta["broad_top_marker_score"] = best_score
    meta["broad_type"] = np.where(assigned, best_name, "Other_unassigned")
    return meta


def scan_gse200997_broad() -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    log("GSE200997: inferring broad target lineages from its original tumor/normal metadata and multi-marker expression")
    matrix_path = RAW / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
    ann_path = RAW / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"
    ann = pd.read_csv(ann_path, dtype=str)
    ann = ann[ann["Condition"].eq("Tumor")].copy()
    ann = ann.rename(columns={"Unnamed: 0": "cell"})
    source_cells = parse_matrix_header(matrix_path, ",")
    keep_set = set(ann["cell"])
    selected_cells = [c for c in source_cells if c in keep_set]
    if len(selected_cells) != len(ann):
        raise ValueError(f"GSE200997 matrix/annotation cell mismatch: {len(ann)-len(selected_cells)} annotation cells absent")
    ann_idx = ann.set_index("cell", drop=False)
    meta = ann_idx.loc[selected_cells].reset_index(drop=True)
    source_index = {c: i for i, c in enumerate(source_cells)}
    selected_source = np.asarray([source_index[c] for c in selected_cells], dtype=np.int32)
    source_n = len(source_cells)
    marker_set = set().union(*[set(x) for x in BROAD_MODULES.values()])
    data_parts: list[np.ndarray] = []; index_parts: list[np.ndarray] = []; indptr = [0]; gene_names: list[str] = []
    libsize = np.zeros(len(selected_cells), dtype=np.int64)
    t0 = time.time()
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        fh.readline()
        for i, line in enumerate(fh):
            gene, vals = parse_gene_line(line, ",")
            if vals.size != source_n:
                raise ValueError(f"GSE200997 row {i + 1} has {vals.size} values; expected {source_n}")
            target = vals[selected_source]
            libsize += target
            if gene in marker_set:
                gene_names.append(gene)
                nz = np.flatnonzero(target)
                if nz.size:
                    data_parts.append(target[nz]); index_parts.append(nz.astype(np.int32, copy=False)); indptr.append(indptr[-1] + nz.size)
                else:
                    indptr.append(indptr[-1])
            if (i + 1) % 5000 == 0:
                log(f"GSE200997 broad-lineage pass streamed {i + 1:,} genes")
    gxc = make_sparse_gene_by_cell(data_parts, index_parts, indptr, len(selected_cells), gene_names)
    expr = gxc.T.tocsr().astype(np.float32)
    scale = np.divide(10000.0, libsize, out=np.zeros_like(libsize, dtype=np.float32), where=libsize > 0)
    expr = expr.multiply(scale[:, None]).tocsr()
    expr.data = np.log1p(expr.data)
    meta = classify_gse200997_broad(expr, gene_names, meta)
    meta["cohort"] = "GSE200997"
    meta["sample_id"] = meta["samples"]
    meta["patient_id"] = ""
    meta["original_type"] = "Marker_inferred"
    meta["original_subtype"] = ""
    counts = meta.groupby(["sample_id", "broad_type"], observed=True).size().rename("n_cells").reset_index()
    counts.to_csv(OUT / "GSE200997_broad_inferred_counts.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "GSE200997_broad_marker_annotations.tsv.gz", sep="\t", index=False, compression="gzip")
    log(f"GSE200997 marker inference completed in {time.time()-t0:.1f}s; broad labels:\n" + meta["broad_type"].value_counts().to_string())
    return meta, selected_source, selected_cells, libsize


def scan_gse200997_full_target(meta: pd.DataFrame, selected_cells: list[str], panel_genes: set[str]) -> None:
    matrix_path = RAW / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
    log("GSE200997: streaming full transcriptome for marker-inferred myeloid/fibroblast candidates and CellChat panel for target cells")
    meta = meta.copy()
    keep = meta["broad_type"].isin(["Epithelial", "Myeloid", "Fibroblast"])
    meta = meta.loc[keep].copy().reset_index(drop=True)
    cell_keep = set(meta["cell"])
    cluster_keep = set(meta.loc[meta["broad_type"].isin(["Myeloid", "Fibroblast"]), "cell"])
    source_cells = parse_matrix_header(matrix_path, ",")
    target_cells = [c for c in source_cells if c in cell_keep]
    cluster_cells = [c for c in target_cells if c in cluster_keep]
    meta = meta.set_index("cell", drop=False).loc[target_cells].reset_index(drop=True)
    source_index = {c: i for i, c in enumerate(source_cells)}
    target_source = np.asarray([source_index[c] for c in target_cells], dtype=np.int32)
    cluster_local = np.flatnonzero(meta["broad_type"].isin(["Myeloid", "Fibroblast"]).to_numpy())
    cluster_source = target_source[cluster_local]
    source_n = len(source_cells)
    cl_data: list[np.ndarray] = []; cl_idx: list[np.ndarray] = []; cl_ptr = [0]; cl_genes: list[str] = []
    pn_data: list[np.ndarray] = []; pn_idx: list[np.ndarray] = []; pn_ptr = [0]; pn_genes: list[str] = []
    libsize = np.zeros(len(target_cells), dtype=np.int64)
    t0 = time.time()
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        fh.readline()
        for i, line in enumerate(fh):
            gene, vals = parse_gene_line(line, ",")
            if vals.size != source_n:
                raise ValueError(f"GSE200997 row {i + 1} has {vals.size} values; expected {source_n}")
            target = vals[target_source]
            libsize += target
            cv = vals[cluster_source]
            nz = np.flatnonzero(cv)
            if nz.size:
                cl_data.append(cv[nz]); cl_idx.append(nz.astype(np.int32, copy=False)); cl_ptr.append(cl_ptr[-1] + nz.size)
            else:
                cl_ptr.append(cl_ptr[-1])
            cl_genes.append(gene)
            if gene in panel_genes:
                pn_genes.append(gene)
                pnz = np.flatnonzero(target)
                if pnz.size:
                    pn_data.append(target[pnz]); pn_idx.append(pnz.astype(np.int32, copy=False)); pn_ptr.append(pn_ptr[-1] + pnz.size)
                else:
                    pn_ptr.append(pn_ptr[-1])
            if (i + 1) % 5000 == 0:
                log(f"GSE200997 target-cell pass streamed {i + 1:,} genes")
    cl_gxc = make_sparse_gene_by_cell(cl_data, cl_idx, cl_ptr, len(cluster_local), cl_genes)
    panel_gxc = make_sparse_gene_by_cell(pn_data, pn_idx, pn_ptr, len(target_cells), pn_genes)
    cl_cxg = cl_gxc.T.tocsr()
    del cl_gxc, cl_data, cl_idx, pn_data, pn_idx
    gc.collect()
    labels = pd.Series(index=meta.loc[meta["broad_type"].isin(["Myeloid", "Fibroblast"]), "cell"], dtype="object")
    cluster_obs = meta.iloc[cluster_local].set_index("cell", drop=False)
    cluster_tables: list[pd.DataFrame] = []; marker_tables: list[pd.DataFrame] = []; assignments: list[pd.DataFrame] = []
    for axis, broad in [("Myeloid", "Myeloid"), ("Stromal", "Fibroblast")]:
        rows = np.flatnonzero(cluster_obs["broad_type"].eq(broad).to_numpy())
        if not len(rows):
            continue
        l, t, m, a = recluster_group(cl_cxg[rows].tocsr(), cl_genes, cluster_obs.iloc[rows], axis, "GSE200997")
        labels.loc[l.index] = l
        cluster_tables.append(t); marker_tables.append(m); assignments.append(a)
    del cl_cxg
    gc.collect()
    meta = meta.set_index("cell", drop=False)
    meta["celltype"] = meta["broad_type"]
    meta.loc[labels.index, "celltype"] = labels
    meta["cellchat_include"] = meta["celltype"].isin(["Epithelial", "Monocyte", "Macrophage", "Other_myeloid", "myCAF", "iCAF", "ECM_CAF", "apCAF"])
    meta_out = meta.reset_index(drop=True)
    meta_out.to_csv(OUT / "GSE200997_cell_annotations.tsv.gz", sep="\t", index=False, compression="gzip")
    if cluster_tables:
        pd.concat(cluster_tables, ignore_index=True).to_csv(OUT / "GSE200997_cluster_program_scores.tsv", sep="\t", index=False)
    if marker_tables:
        marker_tables = [x for x in marker_tables if not x.empty]
        if marker_tables:
            pd.concat(marker_tables, ignore_index=True).to_csv(OUT / "GSE200997_cluster_markers.tsv.gz", sep="\t", index=False, compression="gzip")
    if assignments:
        pd.concat(assignments, ignore_index=True).to_csv(OUT / "GSE200997_cluster_assignments.tsv.gz", sep="\t", index=False, compression="gzip")
    meta_out.groupby(["sample_id", "broad_type", "celltype"], observed=True).size().rename("n_cells").reset_index().to_csv(
        OUT / "GSE200997_cell_counts_by_sample.tsv", sep="\t", index=False
    )
    log("GSE200997 final annotated groups:\n" + meta_out["celltype"].value_counts().to_string())

    panel_cells_by_gene = panel_gxc.T.tocsr()
    scale = np.divide(10000.0, libsize, out=np.zeros_like(libsize, dtype=np.float32), where=libsize > 0)
    panel_cells_by_gene = panel_cells_by_gene.multiply(scale[:, None]).tocsr().astype(np.float32)
    panel_cells_by_gene.data = np.log1p(panel_cells_by_gene.data)
    panel_cells_by_gene.eliminate_zeros()
    pmeta = meta_out.set_index("cell").loc[target_cells].reset_index()
    pmeta["cell"] = target_cells
    add_expr_scores(pmeta, panel_cells_by_gene, pn_genes, "GSE200997")
    pmeta["cellchat_include"] = pmeta["celltype"].isin(["Epithelial", "Monocyte", "Macrophage", "Other_myeloid", "myCAF", "iCAF", "ECM_CAF", "apCAF"])
    include = pmeta["cellchat_include"].to_numpy()
    out_expr = panel_cells_by_gene[include].tocsr()
    out_meta = pmeta.loc[include].reset_index(drop=True)
    write_matrix_market(out_expr, pn_genes, out_meta, "GSE200997")
    summarize_targets(out_meta, OUT / "GSE200997_target_gene_sample_summary.tsv")
    del panel_gxc, panel_cells_by_gene, out_expr
    gc.collect()


def main() -> None:
    log(f"Scanpy {sc.__version__}; output: {OUT}")
    panel_genes = load_panel()
    log(f"CellChat + marker panel genes requested: {len(panel_genes)}")
    scan_gse132465(panel_genes)
    gc.collect()
    gse2_meta, _, gse2_cells, _ = scan_gse200997_broad()
    scan_gse200997_full_target(gse2_meta, gse2_cells, panel_genes)
    manifest = {
        "GSE132465_source_annotation": "GEO processed cell annotation; tumor epithelial, myeloid and stromal broad labels retained as candidate compartments",
        "GSE132465_subclustering": "Full raw transcriptome within original tumor myeloid and stromal cells; Scanpy HVG/PCA/neighbors/Leiden; marker-program cluster labels",
        "GSE200997_source_annotation": "GEO processed metadata has tumor/normal, sample and MSI only; target lineages are assigned only when their multi-marker score wins against epithelial, myeloid, fibroblast, T/NK, B/plasma, endothelial, pericyte and mast controls; ambiguous cells remain unassigned",
        "GSE200997_subclustering": "Full raw transcriptome within marker-inferred tumor myeloid and fibroblast candidates; Scanpy HVG/PCA/neighbors/Leiden; marker-program cluster labels",
        "target_genes": ["SCD", "PID1"],
        "target_gene_usage": "Neither SCD nor PID1 used for cell clustering or annotation. SCD is summarized against a lipid-metabolism marker score that excludes SCD; PID1 is summarized after myeloid subtype annotation.",
        "normalization_for_CellChat": "log1p(raw UMI / complete-cell UMI total * 10000), restricted to CellChatDB.human ligand/receptor symbols plus marker genes",
        "minimum_cluster_group_cells_for_CellChat": 10,
        "GSE200997_unit": "Tumor specimen/sample; no patient identifier is supplied in processed metadata",
    }
    (OUT / "preparation_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    log("CellChat matrix preparation finished")


if __name__ == "__main__":
    main()
