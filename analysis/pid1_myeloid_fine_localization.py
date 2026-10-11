#!/usr/bin/env python
"""Independent PID1 fine-localization analysis in CRC tumor myeloid cells.

Inputs are the original GEO processed raw UMI matrices and original GEO cell
annotations only. Outputs are written to a new, task-specific result directory.
"""
from __future__ import annotations

import gzip
import itertools
import json
import math
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scrublet as scr
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr, wilcoxon
from statsmodels.formula.api import ols, glm
from statsmodels.genmod.families import Binomial
from statsmodels.stats.multitest import multipletests


ROOT = Path(r"D:\DINP\whynot17-new-project")
SRC = ROOT / "work" / "single_cell_GSE132465_GSE200997"
RUN_NAME = "PID1_myeloid_finelocalization"
OUT = ROOT / "results" / "2026-10-11" / RUN_NAME
WORK = ROOT / "work" / RUN_NAME
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 2
sc.settings.figdir = str(OUT)
np.random.seed(20261011)
sns.set_theme(style="whitegrid", context="notebook")

G132_MATRIX = SRC / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
G132_META = SRC / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
G200_MATRIX = SRC / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
G200_META = SRC / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"

CORE_MYELOID = ["PTPRC", "LYZ", "LST1", "TYROBP", "FCER1G", "CTSS", "CD68", "FCGR1A", "CSF1R"]
MONO_CLASSICAL = ["CD14", "FCN1", "VCAN", "S100A8", "S100A9", "S100A12", "CCR2", "CTSD"]
MONO_FCGR3A = ["FCGR3A", "MS4A7", "IFITM3", "LILRB1", "CX3CR1", "LST1", "FCER1G"]
MONO_INTERMEDIATE = ["CD14", "FCGR3A", "FCN1", "MS4A7", "LST1", "CTSS"]
MAC_C1QC = ["C1QA", "C1QB", "C1QC", "FOLR2", "VSIG4", "MRC1", "CD163"]
MAC_SPP1 = ["SPP1", "TREM2", "APOC1", "GPNMB", "LPL", "APOE"]
MAC_INFLAM = ["IL1B", "TNF", "CXCL8", "CXCL2", "CXCL3", "CCL3", "PTGS2", "IL1RN"]
MAC_IFN = ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "RSAD2", "STAT1"]
DC_MARKERS = ["FCER1A", "CD1C", "CLEC10A", "CLEC9A", "IRF8", "CLEC4C"]
DC_PDC_MARKERS = ["LILRA4", "GZMB", "IL3RA", "CLEC4C", "TCF4", "IRF7"]
DC_MATURE_MARKERS = ["LAMP3", "CCR7", "FSCN1", "CCL22", "CD83", "BIRC3", "EBI3"]

COMPETING_LINEAGES = {
    "T_NK": ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "CD2", "CD7", "IL7R", "NKG7", "GNLY", "CD4", "CD8A", "CD8B", "CCL5", "KLRB1", "LTB"],
    "B_plasma": ["MS4A1", "CD79A", "CD79B", "CD19", "CD22", "CD37", "BANK1", "JCHAIN", "MZB1", "SDC1"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CEACAM5"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "EMCN"],
    "mast": ["TPSAB1", "TPSB2", "KIT"],
    "neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "MPO"],
    "erythroid": ["HBA1", "HBA2", "HBB", "ALAS2"],
}
STRICT_COMPETING_LINEAGES = {
    "T_NK": ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2"],
    "B_plasma": ["MS4A1", "CD79A", "CD79B", "CD19", "CD22", "BANK1"],
    "epithelial": ["EPCAM", "CEACAM5", "KRT8", "KRT18", "KRT19"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "EMCN"],
    "mast": ["TPSAB1", "TPSB2", "KIT"],
    "neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "MPO"],
    "erythroid": ["HBA1", "HBA2", "HBB", "ALAS2"],
}

# Fixed, PID1-free exploratory marker programs. They are not claimed to be
# standardized pathway signatures; observed genes are recorded in the outputs.
PROGRAMS = {
    "recruitment_chemotaxis": ["CCR2", "CCL2", "CCL7", "CXCR4", "CX3CR1", "SELL", "ITGA4", "ITGB1", "SELPLG", "ICAM1", "ITGAL"],
    "inflammation": ["IL1B", "TNF", "IL6", "CXCL8", "CXCL2", "CXCL3", "CCL3", "CCL4", "NFKBIA", "PTGS2", "IL1RN", "TNFAIP3"],
    "interferon_response": ["ISG15", "IFI6", "IFI27", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "RSAD2", "STAT1", "IRF7", "CXCL10"],
    "glycolysis": ["SLC2A3", "HK2", "PFKP", "ALDOA", "GAPDH", "ENO1", "PKM", "LDHA"],
    "lipid_cholesterol": ["SCD", "FASN", "ACSS2", "HMGCR", "LDLR", "APOE", "LPL", "FABP5", "CPT1A", "CD36", "ABCA1", "ABCG1", "HMOX1"],
    "antigen_presentation": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "CD74", "B2M", "TAP1", "TAP2"],
}
for _name in PROGRAMS:
    PROGRAMS[_name] = [g for g in PROGRAMS[_name] if g.upper() != "PID1"]


def say(message: str) -> None:
    print(pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"), message, flush=True)


def norm_symbol(value: object) -> str:
    return str(value).strip().upper()


def read_metadata(path: Path, sep: str, cohort: str) -> pd.DataFrame:
    meta = pd.read_csv(path, sep=sep, compression="gzip", dtype=str, keep_default_na=False)
    first = meta.columns[0]
    meta = meta.rename(columns={first: "cell_id"})
    if cohort == "GSE132465":
        meta = meta.rename(columns={"Patient": "unit_id", "Sample": "sample_id", "Class": "condition", "Cell_type": "original_broad_type", "Cell_subtype": "original_subtype"})
        meta["unit_type"] = "patient"
        tumor = meta[meta["condition"].str.casefold().eq("tumor")].copy()
        candidate = tumor[tumor["original_broad_type"].str.casefold().str.contains("myeloid", na=False)].copy()
    else:
        meta = meta.rename(columns={"samples": "sample_id", "Condition": "condition", "Location": "location", "MSI_Status": "msi_status", "bulk_prediction": "cms_bulk", "prediction": "cms_cell"})
        meta["unit_id"] = meta["sample_id"]
        meta["unit_type"] = "sample"
        tumor = meta[meta["condition"].str.casefold().eq("tumor")].copy()
        candidate = tumor.copy()
    return meta, tumor, candidate


def read_header_and_features(matrix_path: Path, sep: str):
    header = pd.read_csv(matrix_path, sep=sep, compression="gzip", nrows=0).columns.tolist()
    feature_col = header[0]
    features = pd.read_csv(matrix_path, sep=sep, compression="gzip", usecols=[feature_col], dtype=str)[feature_col].fillna("").astype(str).tolist()
    return header, feature_col, features


def marker_screen_gse200997(matrix_path: Path, tumor_cells: list[str], marker_genes: set[str]):
    """Read only curated marker rows across all tumor cells from the original matrix."""
    sep = ","
    header = pd.read_csv(matrix_path, sep=sep, compression="gzip", nrows=0).columns.tolist()
    cell_to_pos = {cell: i for i, cell in enumerate(header)}
    selected = [cell for cell in header[1:] if cell in set(tumor_cells)]
    positions = np.asarray([cell_to_pos[c] for c in selected], dtype=np.int64)
    values: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt", encoding="utf-8", newline="") as handle:
        _ = handle.readline()
        for line in handle:
            gene, sep_found, _rest = line.partition(sep)
            # GEO's processed matrix quotes the first (gene-symbol) field.
            symbol = norm_symbol(gene.strip().strip('"'))
            if not sep_found or symbol not in marker_genes or symbol in values:
                continue
            fields = line.rstrip("\r\n").split(sep)
            values[symbol] = np.asarray([fields[j] for j in positions], dtype=np.float32)
    return selected, values


def module_sum(marker_values: dict[str, np.ndarray], genes: list[str], n_cells: int) -> np.ndarray:
    present = [norm_symbol(g) for g in genes if norm_symbol(g) in marker_values]
    if not present:
        return np.zeros(n_cells, dtype=np.float32)
    stacked = np.vstack([np.log1p(marker_values[g]) for g in present])
    return stacked.mean(axis=0)


def broad_myeloid_screen(matrix_path: Path, tumor_meta: pd.DataFrame) -> tuple[list[str], pd.DataFrame]:
    all_markers = set(CORE_MYELOID + MONO_CLASSICAL + MONO_FCGR3A + MONO_INTERMEDIATE + MAC_C1QC + MAC_SPP1 + MAC_INFLAM + MAC_IFN + DC_MARKERS)
    for panel in COMPETING_LINEAGES.values():
        all_markers.update(panel)
    cells, values = marker_screen_gse200997(matrix_path, tumor_meta["cell_id"].tolist(), all_markers)
    n = len(cells)
    core_present = [norm_symbol(g) for g in CORE_MYELOID if norm_symbol(g) in values]
    core_detect = np.vstack([values[g] > 0 for g in core_present]).sum(axis=0) if core_present else np.zeros(n)
    lyz = values.get("LYZ", np.zeros(n)) > 0
    orthogonal_support = np.zeros(n, dtype=bool)
    for gene in ["LST1", "TYROBP", "FCER1G", "CTSS", "CD68", "FCGR1A", "CSF1R"]:
        if norm_symbol(gene) in values:
            orthogonal_support |= values[norm_symbol(gene)] > 0
    # PTPRC is shared by all leukocytes; require a myeloid-specific support
    # marker alongside LYZ or at least three distinct core-myeloid markers.
    myeloid_support = (lyz & orthogonal_support) | (core_detect >= 3)
    myeloid_score = module_sum(values, CORE_MYELOID, n)
    comp_scores = {name: module_sum(values, genes, n) for name, genes in COMPETING_LINEAGES.items()}
    comp_hits = {name: np.vstack([(values.get(norm_symbol(g), np.zeros(n)) > 0) for g in genes if norm_symbol(g) in values]).sum(axis=0) if any(norm_symbol(g) in values for g in genes) else np.zeros(n) for name, genes in COMPETING_LINEAGES.items()}
    max_comp = np.max(np.vstack(list(comp_scores.values())), axis=0) if comp_scores else np.zeros(n)
    candidate = (core_detect >= 2) & myeloid_support & (myeloid_score >= max_comp)
    lineage_conflict = np.zeros(n, dtype=bool)
    conflict_names = []
    for name, hits in comp_hits.items():
        lineage_conflict |= hits >= 2
        conflict_names.append(np.where(hits >= 2, name, ""))
    if conflict_names:
        conflict_label = []
        for i in range(n):
            conflict_label.append(";".join(sorted({arr[i] for arr in conflict_names if arr[i]})))
    else:
        conflict_label = [""] * n
    if len(values) < 20:
        raise RuntimeError(f"GSE200997 marker screen parsed only {len(values)} marker rows; check quoted CSV gene names")
    screen = pd.DataFrame({
        "cell_id": cells,
        "myeloid_marker_hits": core_detect,
        "myeloid_module_score": myeloid_score,
        "strongest_competing_score": max_comp,
        "candidate_myeloid": candidate,
        "multi_lineage_marker_conflict": candidate & lineage_conflict,
        "competing_marker_conflict": conflict_label,
    })
    for name, score in comp_scores.items():
        screen[f"score_{name}"] = score
    return screen.loc[screen["candidate_myeloid"], "cell_id"].tolist(), screen


def read_selected_counts(matrix_path: Path, sep: str, cell_ids: list[str], meta: pd.DataFrame, cohort: str) -> ad.AnnData:
    header = pd.read_csv(matrix_path, sep=sep, compression="gzip", nrows=0).columns.tolist()
    feature_col = header[0]
    wanted = set(cell_ids)
    cells = [c for c in header[1:] if c in wanted]
    if not cells:
        raise ValueError(f"No selected cell IDs matched {cohort} raw matrix")
    dtype = {feature_col: str}
    dtype.update({c: np.int32 for c in cells})
    row_blocks = []
    genes: list[str] = []
    reader = pd.read_csv(matrix_path, sep=sep, compression="gzip", usecols=[feature_col] + cells, dtype=dtype, chunksize=192, low_memory=False)
    for block in reader:
        genes.extend(block[feature_col].astype(str).tolist())
        if block.columns[1:].tolist() != cells:
            block = block[[*cells]]
        arr = block[cells].to_numpy(dtype=np.int32, copy=False)
        row_blocks.append(sparse.csr_matrix(arr))
    gene_by_cell = sparse.vstack(row_blocks, format="csr", dtype=np.int32)
    X = gene_by_cell.T.tocsr()
    obs = meta.set_index("cell_id", drop=False).loc[cells].copy()
    obs.index = cells
    var = pd.DataFrame({"gene_symbol": genes, "gene_upper": [norm_symbol(g) for g in genes]}, index=pd.Index(genes, dtype=str))
    result = ad.AnnData(X=X, obs=obs, var=var)
    result.var_names_make_unique(join="_")
    result.obs["cohort"] = cohort
    result.obs["sample_id"] = result.obs["sample_id"].astype(str)
    result.obs["unit_id"] = result.obs["unit_id"].astype(str)
    result.layers["counts"] = result.X.copy()
    result.uns["input_source"] = str(matrix_path)
    return result


def make_gene_lookup(adata: ad.AnnData) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for i, gene in enumerate(adata.var["gene_upper"].tolist()):
        lookup.setdefault(gene, i)
    return lookup


def sum_for_genes(adata: ad.AnnData, genes: list[str], layer: str = "counts") -> np.ndarray:
    lookup = make_gene_lookup(adata)
    idx = [lookup[norm_symbol(g)] for g in genes if norm_symbol(g) in lookup]
    if not idx:
        return np.zeros(adata.n_obs, dtype=np.float32)
    X = adata.layers[layer][:, idx] if layer in adata.layers else adata.X[:, idx]
    return np.asarray(X.sum(axis=1)).ravel()


def marker_conflicts(adata: ad.AnnData) -> tuple[np.ndarray, list[str], dict[str, np.ndarray]]:
    lookup = make_gene_lookup(adata)
    count_hits: dict[str, np.ndarray] = {}
    evidence: dict[str, np.ndarray] = {}
    for name, genes in COMPETING_LINEAGES.items():
        idx = [lookup[norm_symbol(g)] for g in genes if norm_symbol(g) in lookup]
        if idx:
            X = adata.layers["counts"][:, idx]
            count_hits[name] = np.asarray((X > 0).sum(axis=1)).ravel()
            evidence[f"marker_umi_{name}"] = np.asarray(X.sum(axis=1)).ravel()
        else:
            count_hits[name] = np.zeros(adata.n_obs, dtype=int)
            evidence[f"marker_umi_{name}"] = np.zeros(adata.n_obs, dtype=int)
        evidence[f"marker_hits_{name}"] = count_hits[name]
    for name, genes in STRICT_COMPETING_LINEAGES.items():
        idx = [lookup[norm_symbol(g)] for g in genes if norm_symbol(g) in lookup]
        if idx:
            X = adata.layers["counts"][:, idx]
            evidence[f"strict_marker_hits_{name}"] = np.asarray((X > 0).sum(axis=1)).ravel()
            evidence[f"strict_marker_umi_{name}"] = np.asarray(X.sum(axis=1)).ravel()
        else:
            evidence[f"strict_marker_hits_{name}"] = np.zeros(adata.n_obs, dtype=int)
            evidence[f"strict_marker_umi_{name}"] = np.zeros(adata.n_obs, dtype=int)
    core_idx = [lookup[norm_symbol(g)] for g in CORE_MYELOID if norm_symbol(g) in lookup]
    core_hits = np.asarray((adata.layers["counts"][:, core_idx] > 0).sum(axis=1)).ravel() if core_idx else np.zeros(adata.n_obs)
    evidence["myeloid_marker_hits"] = core_hits
    conflict = np.zeros(adata.n_obs, dtype=bool)
    names = []
    for i in range(adata.n_obs):
        hit = [name for name, values in count_hits.items() if values[i] >= 2]
        names.append(";".join(hit))
        conflict[i] = core_hits[i] >= 2 and bool(hit)
    return conflict, names, evidence


def robust_sample_outliers(obs: pd.DataFrame, values: pd.Series, key: str, min_n: int = 30, nmads: float = 4.0) -> np.ndarray:
    out = np.zeros(len(obs), dtype=bool)
    vals = np.asarray(values, dtype=float)
    for _unit, idx in obs.groupby(key, observed=True).indices.items():
        idx = np.asarray(idx, dtype=int)
        if len(idx) < min_n:
            continue
        local = vals[idx]
        med = np.nanmedian(local)
        mad = np.nanmedian(np.abs(local - med))
        if not np.isfinite(mad) or mad <= 0:
            continue
        out[idx] = np.abs(local - med) > nmads * 1.4826 * mad
    return out


def add_qc_and_doublet_calls(adata: ad.AnnData, all_cells_per_sample: dict[str, int], cohort: str) -> pd.DataFrame:
    X = adata.layers["counts"].tocsr()
    total = np.asarray(X.sum(axis=1)).ravel().astype(float)
    ngenes = np.asarray((X > 0).sum(axis=1)).ravel().astype(int)
    mito_idx = np.flatnonzero(adata.var["gene_upper"].str.startswith(("MT-", "MT.")).to_numpy())
    mt_sum = np.asarray(X[:, mito_idx].sum(axis=1)).ravel().astype(float) if len(mito_idx) else np.zeros(adata.n_obs)
    pct_mt = np.divide(mt_sum, total, out=np.zeros_like(total), where=total > 0) * 100
    adata.obs["total_umi"] = total
    adata.obs["n_genes"] = ngenes
    adata.obs["pct_mito"] = pct_mt
    adata.obs["log10_total_umi"] = np.log10(total + 1)
    out_total = robust_sample_outliers(adata.obs, pd.Series(np.log10(total + 1)), "sample_id", nmads=4.0)
    out_genes = robust_sample_outliers(adata.obs, pd.Series(ngenes), "sample_id", nmads=4.0)
    adata.obs["qc_low_library"] = total < 500
    adata.obs["qc_low_gene_count"] = ngenes < 200
    adata.obs["qc_high_mito"] = pct_mt > 20
    adata.obs["qc_sample_outlier"] = out_total | out_genes
    adata.obs["qc_pass"] = ~(adata.obs[["qc_low_library", "qc_low_gene_count", "qc_high_mito", "qc_sample_outlier"]].any(axis=1).to_numpy())
    conflict, names, marker_evidence = marker_conflicts(adata)
    adata.obs["marker_lineage_conflict"] = conflict
    adata.obs["competing_lineages_detected"] = names
    for col, values in marker_evidence.items():
        adata.obs[col] = values
    adata.obs["scrublet_score"] = np.nan
    adata.obs["scrublet_predicted_doublet"] = False
    adata.obs["scrublet_run"] = False
    rows = []
    for sample, idx in adata.obs.groupby("sample_id", observed=True).indices.items():
        idx = np.asarray(idx, dtype=int)
        qidx = idx[adata.obs.iloc[idx]["qc_pass"].to_numpy()]
        n = len(qidx)
        threshold = np.nan
        n_pred = 0
        run_ok = False
        status = "not_run_lt80_qc_pass_myeloid"
        if n >= 80:
            sample_total = all_cells_per_sample.get(str(sample), n)
            expected_rate = min(0.10, max(0.01, 0.008 * (sample_total / 1000.0)))
            n_pcs = min(20, n - 2, adata.n_vars - 2)
            try:
                detector = scr.Scrublet(X[qidx], expected_doublet_rate=expected_rate, sim_doublet_ratio=2.0, random_state=20261011)
                scores, calls = detector.scrub_doublets(
                    use_approx_neighbors=False,
                    n_prin_comps=n_pcs,
                    min_gene_variability_pctl=85,
                    verbose=False,
                )
                adata.obs.iloc[qidx, adata.obs.columns.get_loc("scrublet_score")] = scores
                adata.obs.iloc[qidx, adata.obs.columns.get_loc("scrublet_predicted_doublet")] = calls
                adata.obs.iloc[qidx, adata.obs.columns.get_loc("scrublet_run")] = True
                threshold = float(detector.threshold_) if detector.threshold_ is not None else np.nan
                n_pred = int(np.sum(calls))
                run_ok = True
                status = "scrublet_completed"
            except Exception as exc:  # preserve unit-level diagnostics and continue
                status = f"scrublet_failed:{type(exc).__name__}:{str(exc)[:160]}"
        rows.append({
            "cohort": cohort, "sample_id": str(sample), "n_candidate_cells": len(idx), "n_qc_pass": n,
            "n_cells_in_original_sample": sample_total if n >= 80 else all_cells_per_sample.get(str(sample), np.nan),
            "expected_doublet_rate_for_scrublet": expected_rate if n >= 80 else np.nan,
            "scrublet_threshold": threshold, "n_scrublet_calls": n_pred, "scrublet_fraction": n_pred / n if n else np.nan,
            "scrublet_completed": run_ok, "status": status,
        })
    # Multi-lineage marker co-detection is retained as an auditable flag: low-level
    # epithelial/stromal RNA can occur in tumor myeloid cells and is not by itself a
    # doublet call. Scrublet calls and sample-level QC determine automatic exclusion.
    adata.obs["exclude_doublet"] = adata.obs["scrublet_predicted_doublet"].astype(bool)
    adata.obs["analysis_keep"] = adata.obs["qc_pass"].astype(bool) & ~adata.obs["exclude_doublet"].astype(bool)
    return pd.DataFrame(rows)


def add_log_norm_and_cluster(adata: ad.AnnData, cohort: str) -> ad.AnnData:
    # The full matrix is retained in layers['counts']; only log1p CPM is used for embeddings.
    work = adata[adata.obs["analysis_keep"].to_numpy()].copy()
    work.X = work.layers["counts"].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    work.layers["lognorm"] = work.X.copy()
    work.raw = work
    sc.pp.filter_genes(work, min_cells=3)
    n_top = min(2000, max(500, work.n_vars - 1))
    try:
        sc.pp.highly_variable_genes(work, n_top_genes=n_top, flavor="seurat", batch_key="sample_id")
    except Exception as exc:
        say(f"{cohort}: batch-aware HVG fallback because {type(exc).__name__}: {exc}")
        sc.pp.highly_variable_genes(work, n_top_genes=n_top, flavor="seurat")
    hvg = work.var["highly_variable"].to_numpy()
    emb = work[:, hvg].copy()
    sc.pp.scale(emb, max_value=10)
    n_pcs = min(30, emb.n_obs - 2, emb.n_vars - 2)
    sc.tl.pca(emb, n_comps=n_pcs, svd_solver="arpack", random_state=20261011)
    sc.pp.neighbors(emb, n_neighbors=min(20, max(5, int(math.sqrt(emb.n_obs)))), n_pcs=n_pcs, random_state=20261011)
    sc.tl.leiden(emb, resolution=0.6, key_added="leiden", random_state=20261011, flavor="igraph")
    sc.tl.umap(emb, min_dist=0.35, random_state=20261011)
    work.obs["leiden"] = emb.obs["leiden"].astype(str).values
    work.obsm["X_pca"] = emb.obsm["X_pca"].copy()
    work.obsm["X_umap"] = emb.obsm["X_umap"].copy()
    work.obsp["connectivities"] = emb.obsp["connectivities"].copy()
    work.obsp["distances"] = emb.obsp["distances"].copy()
    work.uns["neighbors"] = emb.uns["neighbors"].copy()
    work.uns["pid1_hvg_n"] = int(hvg.sum())
    work.uns["pid1_hvg_method"] = "Scanpy seurat HVG, batch_key=sample_id when available"
    return work


def gene_zscore_matrix(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    lookup = make_gene_lookup(adata)
    genes_present = [norm_symbol(g) for g in genes if norm_symbol(g) in lookup and norm_symbol(g) != "PID1"]
    if not genes_present:
        return np.zeros(adata.n_obs, dtype=np.float32), []
    idx = [lookup[g] for g in genes_present]
    mat = adata.layers["lognorm"][:, idx]
    arr = mat.toarray() if sparse.issparse(mat) else np.asarray(mat)
    mu = np.mean(arr, axis=0)
    sd = np.std(arr, axis=0, ddof=0)
    sd[sd < 1e-6] = 1.0
    score = ((arr - mu) / sd).mean(axis=1)
    return score.astype(np.float32), genes_present


def annotate_fine_states(adata: ad.AnnData, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_sets = {
        "mono_classical_score": MONO_CLASSICAL,
        "mono_fcgr3a_score": MONO_FCGR3A,
        "mono_intermediate_score": MONO_INTERMEDIATE,
        "mac_c1qc_score": MAC_C1QC,
        "mac_spp1_score": MAC_SPP1,
        "mac_inflammatory_score": MAC_INFLAM,
        "mac_ifn_score": MAC_IFN,
        "dc_score": DC_MARKERS,
        "dc_pdc_score": DC_PDC_MARKERS,
        "dc_mature_score": DC_MATURE_MARKERS,
    }
    measured = []
    for score_name, genes in score_sets.items():
        score, present = gene_zscore_matrix(adata, genes)
        adata.obs[score_name] = score
        measured.append({"cohort": cohort, "program": score_name, "genes_requested": ";".join(genes), "genes_measured": ";".join(present), "n_genes_measured": len(present), "PID1_excluded": True})
    cluster_scores = adata.obs.groupby("leiden", observed=True)[list(score_sets)].mean()
    cluster_counts = adata.obs["leiden"].value_counts().to_dict()
    competing_stats = {}
    for lineage in COMPETING_LINEAGES:
        hit_col = f"marker_hits_{lineage}"
        umi_col = f"marker_umi_{lineage}"
        if hit_col in adata.obs and umi_col in adata.obs:
            competing_stats[lineage] = adata.obs.groupby("leiden", observed=True)[[hit_col, umi_col]].mean()
    rows = []
    for cluster, scores in cluster_scores.iterrows():
        mono_a = float(scores["mono_classical_score"])
        mono_b = float(scores["mono_fcgr3a_score"])
        mono_i = float(scores["mono_intermediate_score"])
        mono_score = max(mono_a, mono_b, mono_i)
        mac_options = {
            "Macrophage_C1QC_FOLR2": float(scores["mac_c1qc_score"]),
            "Macrophage_SPP1_TREM2": float(scores["mac_spp1_score"]),
            "Macrophage_inflammatory": float(scores["mac_inflammatory_score"]),
            "Macrophage_IFN": float(scores["mac_ifn_score"]),
        }
        mac_label = max(mac_options, key=mac_options.get)
        mac_score = mac_options[mac_label]
        dc_score = float(scores["dc_score"])
        # Flag clusters whose cell-level counts strongly support a competing
        # lineage; low-level marker co-detection remains an audit flag only.
        dominant_competitor = ""
        dominant_hits = 0.0
        dominant_umi = 0.0
        for lineage in STRICT_COMPETING_LINEAGES:
            hit_col = f"strict_marker_hits_{lineage}"
            umi_col = f"strict_marker_umi_{lineage}"
            if hit_col not in adata.obs or umi_col not in adata.obs:
                continue
            strict_stats = adata.obs.groupby("leiden", observed=True)[[hit_col, umi_col]].mean()
            hits = float(strict_stats.loc[cluster, hit_col])
            umi = float(strict_stats.loc[cluster, umi_col])
            min_hits, min_umi = (3.0, 8.0) if lineage in {"epithelial", "fibroblast"} else (2.0, 5.0)
            if hits >= min_hits and umi >= min_umi and (hits, umi) > (dominant_hits, dominant_umi):
                dominant_competitor, dominant_hits, dominant_umi = lineage, hits, umi
        dc_options = {
            "DC_like": float(scores["dc_score"]),
            "DC_pDC_like": float(scores["dc_pdc_score"]),
            "DC_mature_migratory": float(scores["dc_mature_score"]),
        }
        dc_label = max(dc_options, key=dc_options.get)
        dc_score = dc_options[dc_label]
        mono_inflammatory = mono_a >= 0.40 and float(scores["mac_c1qc_score"]) < 0 and float(scores["mac_inflammatory_score"]) >= 0.35
        if dominant_competitor:
            label = f"Non_myeloid_{dominant_competitor}_contamination"
        elif dc_score > max(mono_score, mac_score) + 0.10 and dc_score > 0.15:
            label = dc_label
        elif mono_inflammatory:
            label = "Monocyte_CD14_classical_inflammatory"
        elif mono_score >= mac_score:
            if mono_i >= max(mono_a, mono_b) - 0.15 and min(mono_a, mono_b) > 0.05:
                label = "Monocyte_intermediate_CD14_FCGR3A"
            elif mono_b > mono_a:
                label = "Monocyte_FCGR3A"
            else:
                label = "Monocyte_CD14_classical"
        elif mac_score > 0.05:
            if mac_label == "Macrophage_C1QC_FOLR2":
                label = "Macrophage_C1QC"
            elif mac_label == "Macrophage_SPP1_TREM2" and float(scores["mac_spp1_score"])-float(scores["mac_c1qc_score"]) < 0.20:
                label = "Macrophage_C1QC"
            else:
                label = mac_label
        else:
            label = "Myeloid_unresolved"
        rows.append({"cohort": cohort, "leiden": str(cluster), "n_cells": int(cluster_counts.get(cluster, 0)), "annotation": label,
                     "mono_classical_score": mono_a, "mono_fcgr3a_score": mono_b, "mono_intermediate_score": mono_i,
                     "mac_c1qc_score": float(scores["mac_c1qc_score"]), "mac_spp1_score": float(scores["mac_spp1_score"]),
                     "mac_inflammatory_score": float(scores["mac_inflammatory_score"]), "mac_ifn_score": float(scores["mac_ifn_score"]),
                     "dc_score": dc_options["DC_like"], "dc_pdc_score": dc_options["DC_pDC_like"], "dc_mature_score": dc_options["DC_mature_migratory"],
                     "dominant_nonmyeloid_lineage": dominant_competitor, "mean_competing_marker_hits": dominant_hits,
                     "mean_competing_marker_umi": dominant_umi,
                     "cluster_rule": "PID1-free state scores plus cluster-level competing-lineage marker audit; strong non-myeloid identity flagged"})
    cluster_annotation = pd.DataFrame(rows)
    mapping = dict(zip(cluster_annotation["leiden"], cluster_annotation["annotation"]))
    adata.obs["pid1_subtype"] = adata.obs["leiden"].astype(str).map(mapping).fillna("Myeloid_unresolved").astype("category")
    adata.obs["pid1_cellclass"] = np.where(adata.obs["pid1_subtype"].astype(str).str.startswith("Monocyte"), "Monocyte",
                                  np.where(adata.obs["pid1_subtype"].astype(str).str.startswith("Macrophage"), "Macrophage",
                                  np.where(adata.obs["pid1_subtype"].astype(str).str.startswith("Non_myeloid_"), "Non_myeloid", adata.obs["pid1_subtype"].astype(str))))
    adata.obs["cellchat_eligible"] = ~(
        adata.obs["pid1_cellclass"].astype(str).eq("Non_myeloid")
        | adata.obs["pid1_subtype"].astype(str).eq("Myeloid_unresolved")
    )
    cluster_annotation["n_units"] = cluster_annotation["leiden"].map(adata.obs.groupby("leiden", observed=True)["unit_id"].nunique())
    return cluster_annotation, pd.DataFrame(measured)


def summarize_unit_expression(adata: ad.AnnData, cohort: str) -> pd.DataFrame:
    lookup = make_gene_lookup(adata)
    pid1_idx = lookup.get("PID1")
    if pid1_idx is None:
        raise KeyError(f"PID1 missing from {cohort} features")
    pid1 = adata.layers["lognorm"][:, pid1_idx]
    pid1 = np.asarray(pid1.toarray()).ravel() if sparse.issparse(pid1) else np.asarray(pid1).ravel()
    adata.obs["PID1_log1p"] = pid1
    adata.obs["PID1_detected"] = pid1 > 0
    rows = []
    frame = adata.obs.copy()
    for (unit, sample, state), idx in frame.groupby(["unit_id", "sample_id", "pid1_subtype"], observed=True).indices.items():
        ix = np.asarray(idx, dtype=int)
        vals = pid1[ix]
        rows.append({"cohort": cohort, "unit_type": str(frame.iloc[ix[0]]["unit_type"]), "unit_id": str(unit), "sample_id": str(sample),
                     "pid1_subtype": str(state), "n_cells": len(ix), "PID1_mean_log1p": float(np.mean(vals)),
                     "PID1_median_log1p": float(np.median(vals)), "PID1_detected_pct": float(np.mean(vals > 0) * 100)})
    return pd.DataFrame(rows)


def paired_pid1_tests(summary: pd.DataFrame, cohort: str) -> pd.DataFrame:
    rows = []
    all_states = set(summary["pid1_subtype"].astype(str))
    for threshold in [10, 20, 50]:
        # Aggregate all fine states within the broad class, weighting by the number of cells.
        # Each sample/patient remains one inferential unit.
        def collapse(prefix: str) -> pd.DataFrame:
            subset=summary[summary.pid1_subtype.astype(str).str.startswith(prefix)]
            rows=[]
            for (unit_id,sample_id),x in subset.groupby(["unit_id","sample_id"],observed=True):
                rows.append({"unit_id":unit_id,"sample_id":sample_id,f"n_cells_{prefix.lower()}":int(x.n_cells.sum()),
                             f"mean_{prefix.lower()}":float(np.average(x.PID1_mean_log1p,weights=x.n_cells))})
            return pd.DataFrame(rows)
        mono = collapse("Monocyte").rename(columns={"n_cells_monocyte":"n_cells_mono","mean_monocyte":"mean_mono"})
        mac = collapse("Macrophage").rename(columns={"n_cells_macrophage":"n_cells_mac","mean_macrophage":"mean_mac"})
        pair=mono.merge(mac,on=["unit_id","sample_id"],how="inner")
        ncol_mono,ncol_mac,mean_mono,mean_mac="n_cells_mono","n_cells_mac","mean_mono","mean_mac"
        paired=pair[(pair[ncol_mono]>=threshold)&(pair[ncol_mac]>=threshold)].copy()
        delta=paired[mean_mono]-paired[mean_mac]
        p=float(wilcoxon(delta, alternative="two-sided", zero_method="wilcox").pvalue) if len(delta)>=2 and np.any(delta!=0) else np.nan
        rows.append({"cohort":cohort,"comparison":"all reannotated monocytes vs all reannotated macrophages","min_cells_per_class_per_unit":threshold,
                     "n_paired_units":len(paired),"median_monocyte_mean_log1p":float(paired[mean_mono].median()) if len(paired) else np.nan,
                     "median_macrophage_mean_log1p":float(paired[mean_mac].median()) if len(paired) else np.nan,
                     "median_paired_delta_mono_minus_mac":float(delta.median()) if len(delta) else np.nan,
                     "n_delta_positive":int((delta>0).sum()),"n_delta_negative":int((delta<0).sum()),"wilcoxon_p":p})
        # Repeated-unit pairwise tests for the fine monocyte states and macrophage states.
        states = sorted([s for s in all_states if s.startswith("Monocyte") or s.startswith("Macrophage")])
        for s1,s2 in itertools.combinations(states,2):
            a=summary[summary.pid1_subtype.astype(str).eq(s1)].set_index(["unit_id","sample_id"])
            b=summary[summary.pid1_subtype.astype(str).eq(s2)].set_index(["unit_id","sample_id"])
            both=a[["n_cells","PID1_mean_log1p"]].join(b[["n_cells","PID1_mean_log1p"]],lsuffix="_a",rsuffix="_b",how="inner")
            both=both[(both.n_cells_a>=threshold)&(both.n_cells_b>=threshold)]
            d=both.PID1_mean_log1p_a-both.PID1_mean_log1p_b
            pp=float(wilcoxon(d,alternative="two-sided",zero_method="wilcox").pvalue) if len(d)>=2 and np.any(d!=0) else np.nan
            rows.append({"cohort":cohort,"comparison":f"{s1} vs {s2}","min_cells_per_class_per_unit":threshold,"n_paired_units":len(both),
                         "median_monocyte_mean_log1p":float(both.PID1_mean_log1p_a.median()) if len(both) else np.nan,
                         "median_macrophage_mean_log1p":float(both.PID1_mean_log1p_b.median()) if len(both) else np.nan,
                         "median_paired_delta_mono_minus_mac":float(d.median()) if len(d) else np.nan,
                         "n_delta_positive":int((d>0).sum()),"n_delta_negative":int((d<0).sum()),"wilcoxon_p":pp})
    result=pd.DataFrame(rows)
    result["p_BH"] = np.nan
    for threshold, idx in result.groupby("min_cells_per_class_per_unit").groups.items():
        valid=[i for i in idx if pd.notna(result.loc[i,"wilcoxon_p"])]
        if valid:
            result.loc[valid,"p_BH"]=multipletests(result.loc[valid,"wilcoxon_p"],method="fdr_bh")[1]
    return result


def add_program_scores_and_tests(adata: ad.AnnData, cohort: str) -> tuple[pd.DataFrame,pd.DataFrame]:
    mon=adata.obs[adata.obs["pid1_cellclass"].astype(str).eq("Monocyte")].copy()
    if len(mon)<100:
        return pd.DataFrame(),pd.DataFrame()
    idx=mon.index
    sub=adata[idx].copy()
    measured_rows=[]
    program_cols=[]
    for program, genes in PROGRAMS.items():
        score, present=gene_zscore_matrix(sub, genes)
        col=f"score_{program}"
        sub.obs[col]=score
        program_cols.append(col)
        measured_rows.append({"cohort":cohort,"program":program,"genes_requested":";".join(genes),"genes_measured":";".join(present),"n_genes_measured":len(present),"PID1_excluded":True})
    # Regression uses the cell as measurement, sample fixed effects, and cluster-robust SE by unit.
    # Outcome is normalized log1p PID1; technical adjustments are log library size and mitochondrial fraction.
    regression=[]
    for col in program_cols:
        n_genes=int(next((x["n_genes_measured"] for x in measured_rows if f"score_{x['program']}"==col),0))
        if n_genes<3:
            regression.append({"cohort":cohort,"program":col.removeprefix("score_"),"model":"cell-level OLS; sample fixed effects; cluster-robust by unit","n_cells":len(sub),"n_units":sub.obs.unit_id.nunique(),"beta_per_1SD_score":np.nan,"se_cluster":np.nan,"p_value":np.nan,"status":"fewer than 3 measured program genes"})
            continue
        tmp=sub.obs[["PID1_log1p","log10_total_umi","pct_mito","sample_id","unit_id",col]].dropna().copy()
        tmp=tmp.rename(columns={col:"program_score"})
        try:
            fit=ols("PID1_log1p ~ program_score + log10_total_umi + pct_mito + C(sample_id)",data=tmp).fit(cov_type="cluster",cov_kwds={"groups":tmp["unit_id"],"use_correction":True})
            regression.append({"cohort":cohort,"program":col.removeprefix("score_"),"model":"cell-level OLS; sample fixed effects; cluster-robust by unit","n_cells":int(fit.nobs),"n_units":tmp.unit_id.nunique(),"beta_per_1SD_score":float(fit.params.get("program_score",np.nan)),"se_cluster":float(fit.bse.get("program_score",np.nan)),"p_value":float(fit.pvalues.get("program_score",np.nan)),"status":"ok"})
        except Exception as exc:
            regression.append({"cohort":cohort,"program":col.removeprefix("score_"),"model":"cell-level OLS; sample fixed effects; cluster-robust by unit","n_cells":len(tmp),"n_units":tmp.unit_id.nunique(),"beta_per_1SD_score":np.nan,"se_cluster":np.nan,"p_value":np.nan,"status":f"failed:{type(exc).__name__}:{str(exc)[:120]}"})
    reg=pd.DataFrame(regression)
    reg["p_BH"] = multipletests(reg.p_value.fillna(1).values,method="fdr_bh")[1] if len(reg) else []
    # Unit-level correlations are a distinct, non-cell-pseudoreplicated consistency check.
    unit_rows=[]
    scores_obs=sub.obs.copy()
    for col in program_cols:
        unit=scores_obs.groupby(["unit_id","sample_id"],observed=True).agg(PID1_mean=("PID1_log1p","mean"),program_mean=(col,"mean"),n_cells=(col,"size"),median_log10_umi=("log10_total_umi","median")).reset_index()
        unit=unit[unit.n_cells>=10]
        rho,p=spearmanr(unit.PID1_mean,unit.program_mean) if len(unit)>=5 else (np.nan,np.nan)
        unit_rows.append({"cohort":cohort,"program":col.removeprefix("score_"),"n_units":len(unit),"spearman_rho":rho,"p_value":p,"median_unit_cells":float(unit.n_cells.median()) if len(unit) else np.nan})
    units=pd.DataFrame(unit_rows)
    units["p_BH"]=multipletests(units.p_value.fillna(1).values,method="fdr_bh")[1] if len(units) else []
    # Store scores on the original AnnData in cell order for downstream annotation/CellChat.
    for col in program_cols:
        adata.obs.loc[sub.obs.index,col]=sub.obs[col].values
    return reg, pd.concat([pd.DataFrame(measured_rows).assign(result_type="gene_coverage"), units.assign(result_type="unit_correlation")],ignore_index=True,sort=False)


def plot_umap(adata: ad.AnnData, cohort: str) -> None:
    lookup=make_gene_lookup(adata)
    pididx=lookup.get("PID1")
    pid=adata.layers["lognorm"][:,pididx]
    adata.obs["PID1_log1p"]=np.asarray(pid.toarray()).ravel() if sparse.issparse(pid) else np.asarray(pid).ravel()
    fig,axes=plt.subplots(1,3,figsize=(15,4.8))
    coords=adata.obsm["X_umap"]
    sca=axes[0].scatter(coords[:,0],coords[:,1],c=adata.obs["PID1_log1p"],s=5,cmap="viridis",rasterized=True)
    axes[0].set_title(f"{cohort}: PID1 expression")
    fig.colorbar(sca,ax=axes[0],fraction=.046,pad=.04,label="mean log1p normalized PID1")
    palette=sns.color_palette("tab10",n_colors=adata.obs["pid1_subtype"].nunique())
    sns.scatterplot(x=coords[:,0],y=coords[:,1],hue=adata.obs["pid1_subtype"].astype(str),s=5,linewidth=0,ax=axes[1],palette=palette,legend=False)
    axes[1].set_title("Marker-based myeloid states")
    sca2=axes[2].scatter(coords[:,0],coords[:,1],c=np.log10(adata.obs["total_umi"].to_numpy()+1),s=5,cmap="magma",rasterized=True)
    axes[2].set_title("Library size (log10 UMI)")
    fig.colorbar(sca2,ax=axes[2],fraction=.046,pad=.04,label="log10(total UMI + 1)")
    for ax in axes: ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel("UMAP1");ax.set_ylabel("UMAP2")
    fig.tight_layout();fig.savefig(OUT/f"{cohort}_PID1_myeloid_UMAP.png",dpi=220,bbox_inches="tight");plt.close(fig)


def plot_marker_dot(adata: ad.AnnData, cohort: str) -> None:
    genes=["CD14","FCN1","VCAN","S100A8","S100A9","S100A12","FCGR3A","MS4A7","LST1","C1QA","C1QC","FOLR2","VSIG4","SPP1","TREM2","APOE","IL1B","HLA-DRA","CD1C","CLEC10A","LILRA4","IL3RA","GZMB","LAMP3","CCR7","FSCN1","PID1"]
    lookup=make_gene_lookup(adata); genes=[g for g in genes if norm_symbol(g) in lookup]
    states=sorted(adata.obs["pid1_subtype"].astype(str).unique())
    mean=np.zeros((len(states),len(genes))); pct=np.zeros_like(mean)
    for i,state in enumerate(states):
        ix=np.flatnonzero(adata.obs["pid1_subtype"].astype(str).to_numpy()==state)
        for j,g in enumerate(genes):
            v=adata.layers["lognorm"][ix,lookup[norm_symbol(g)]]
            v=np.asarray(v.toarray()).ravel() if sparse.issparse(v) else np.asarray(v).ravel()
            mean[i,j]=float(v.mean());pct[i,j]=float(np.mean(v>0)*100)
    fig,ax=plt.subplots(figsize=(max(12,.52*len(genes)),max(4,.58*len(states))))
    sca=ax.scatter(np.tile(np.arange(len(genes)),len(states)),np.repeat(np.arange(len(states)),len(genes)),s=12+1.4*pct.ravel(),c=mean.ravel(),cmap="viridis",edgecolors="none")
    ax.set_xticks(np.arange(len(genes)),genes,rotation=60,ha="right");ax.set_yticks(np.arange(len(states)),states);ax.invert_yaxis()
    ax.set_xlabel("Marker / target gene");ax.set_ylabel("Marker-defined state");ax.set_title(f"{cohort}: myeloid marker and PID1 dot plot")
    cbar=fig.colorbar(sca,ax=ax,pad=.02);cbar.set_label("Mean log1p normalized expression")
    ax.scatter([],[],s=18,c="#777",label="10% detected");ax.scatter([],[],s=82,c="#777",label="50% detected");ax.legend(frameon=False,title="Cell detection",loc="upper left",bbox_to_anchor=(1.14,1))
    fig.tight_layout();fig.savefig(OUT/f"{cohort}_PID1_myeloid_DotPlot.png",dpi=220,bbox_inches="tight");plt.close(fig)


def plot_violin_detection(adata: ad.AnnData, summary: pd.DataFrame, cohort: str) -> None:
    order=sorted(summary.pid1_subtype.unique())
    obs=adata.obs[["pid1_subtype","PID1_log1p"]].copy();obs["pid1_subtype"]=obs.pid1_subtype.astype(str)
    fig,ax=plt.subplots(figsize=(max(10,.9*len(order)),5))
    sns.violinplot(data=obs,x="pid1_subtype",y="PID1_log1p",order=order,inner="quartile",cut=0,scale="width",ax=ax,color="#7aa6c2")
    ax.set_xlabel("");ax.set_ylabel("PID1 mean log1p normalized expression");ax.set_title(f"{cohort}: cell-level PID1 distribution")
    ax.tick_params(axis="x",rotation=55);fig.tight_layout();fig.savefig(OUT/f"{cohort}_PID1_myeloid_ViolinPlot.png",dpi=220,bbox_inches="tight");plt.close(fig)
    # Replicate-level detection: each patient/sample contributes one observation per state.
    fig,ax=plt.subplots(figsize=(max(10,.9*len(order)),5))
    replicate=summary[summary.n_cells>=10].copy()
    sns.stripplot(data=replicate,x="pid1_subtype",y="PID1_detected_pct",order=order,jitter=.22,size=4,alpha=.65,color="#2c6e91",ax=ax)
    sns.pointplot(data=replicate,x="pid1_subtype",y="PID1_detected_pct",order=order,errorbar="se",join=False,markers="D",color="#d1495b",ax=ax)
    ax.set_xlabel("");ax.set_ylabel("PID1-positive cells per patient/sample (%)");ax.set_title(f"{cohort}: replicate-level PID1 detection rate (units with ≥10 cells)")
    ax.tick_params(axis="x",rotation=55);fig.tight_layout();fig.savefig(OUT/f"{cohort}_PID1_detection_rate.png",dpi=220,bbox_inches="tight");plt.close(fig)


def pseudotime_exploration(adata: ad.AnnData, cohort: str) -> tuple[pd.DataFrame,pd.DataFrame]:
    mask=adata.obs["pid1_cellclass"].astype(str).isin(["Monocyte","Macrophage"]).to_numpy()
    sub=adata[mask].copy()
    states=sub.obs["pid1_cellclass"].astype(str)
    per_state=sub.obs.groupby("pid1_subtype",observed=True).size()
    n_mono=int(states.eq("Monocyte").sum());n_mac=int(states.eq("Macrophage").sum())
    eligible_units=int(sub.obs["unit_id"].nunique())
    status="not_supported_small_or_single_state"
    detail=f"n_monocytes={n_mono}; n_macrophages={n_mac}; units={eligible_units}"
    result=[]; plot_path=""
    if n_mono>=200 and n_mac>=200 and eligible_units>=8 and sub.obs["pid1_subtype"].nunique()>=3:
        conn=sub.obsp.get("connectivities")
        if conn is None:
            status="not_supported_no_neighbor_graph"
        else:
            ncomp,labels=connected_components(conn,directed=False)
            largest=max(np.bincount(labels)) / len(labels)
            detail+=f"; graph_components={ncomp}; largest_component_fraction={largest:.3f}"
            if largest<0.90:
                status="not_supported_disconnected_graph"
            else:
                root_candidates=np.flatnonzero(sub.obs["pid1_subtype"].astype(str).eq("Monocyte_CD14_classical").to_numpy())
                if len(root_candidates)==0:
                    root_candidates=np.flatnonzero(states.eq("Monocyte").to_numpy())
                if len(root_candidates)<100:
                    status="not_supported_no_clear_monocyte_root"
                else:
                    centroid=np.median(sub.obsm["X_pca"][root_candidates,:10],axis=0)
                    root=root_candidates[np.argmin(np.sum((sub.obsm["X_pca"][root_candidates,:10]-centroid)**2,axis=1))]
                    try:
                        sc.tl.diffmap(sub,n_comps=min(15,sub.n_obs-2),random_state=20261011)
                        sub.uns["iroot"]=int(root)
                        sc.tl.dpt(sub,n_dcs=min(10,sub.obsm["X_diffmap"].shape[1]-1))
                        pt=np.asarray(sub.obs["dpt_pseudotime"],dtype=float)
                        if not np.isfinite(pt).all() or np.nanstd(pt)<0.05:
                            status="not_supported_no_pseudotime_gradient"
                        else:
                            status="exploratory_DPT_computed"
                            tmp=sub.obs.copy();tmp["DPT_pseudotime"]=pt
                            f=ols("PID1_log1p ~ DPT_pseudotime + log10_total_umi + pct_mito + C(sample_id)",data=tmp).fit(cov_type="cluster",cov_kwds={"groups":tmp["unit_id"],"use_correction":True})
                            result.append({"cohort":cohort,"analysis":"cell-level OLS with sample fixed effects, cluster-robust by unit","status":status,"n_cells":len(sub),"n_units":sub.obs.unit_id.nunique(),"largest_component_fraction":largest,"PID1_beta_per_DPT_unit":float(f.params.get("DPT_pseudotime",np.nan)),"se_cluster":float(f.bse.get("DPT_pseudotime",np.nan)),"p_value":float(f.pvalues.get("DPT_pseudotime",np.nan)),"caveat":"DPT is an expression ordering, not lineage tracing"})
                            sub.obs["DPT_pseudotime"]=pt
                            fig,axes=plt.subplots(1,2,figsize=(10,4.2))
                            c0=axes[0].scatter(sub.obsm["X_umap"][:,0],sub.obsm["X_umap"][:,1],c=pt,s=6,cmap="viridis",rasterized=True);fig.colorbar(c0,ax=axes[0],label="DPT pseudotime");axes[0].set_title(f"{cohort}: exploratory DPT")
                            c1=axes[1].scatter(sub.obs["DPT_pseudotime"],sub.obs["PID1_log1p"],c=sub.obs["log10_total_umi"],s=6,alpha=.45,cmap="magma",rasterized=True);fig.colorbar(c1,ax=axes[1],label="log10 UMI");axes[1].set_xlabel("DPT pseudotime");axes[1].set_ylabel("PID1 log1p normalized");axes[1].set_title("Not lineage tracing")
                            for ax in axes:ax.set_xlabel(ax.get_xlabel() or "UMAP1");ax.set_xticks([]) if ax is axes[0] else None;ax.set_yticks([]) if ax is axes[0] else None
                            fig.tight_layout();plot_path=OUT/f"{cohort}_PID1_DPT_exploration.png";fig.savefig(plot_path,dpi=220,bbox_inches="tight");plt.close(fig)
                    except Exception as exc:
                        status=f"DPT_failed:{type(exc).__name__}:{str(exc)[:150]}"
    if not result:
        result=[{"cohort":cohort,"analysis":"conditional DPT exploration","status":status,"n_cells":len(sub),"n_units":eligible_units,"largest_component_fraction":np.nan,"PID1_beta_per_DPT_unit":np.nan,"se_cluster":np.nan,"p_value":np.nan,"caveat":detail+"; if computed, DPT is not lineage tracing"}]
    return pd.DataFrame(result),pd.DataFrame([{"cohort":cohort,"pseudotime_status":status,"n_monocytes":n_mono,"n_macrophages":n_mac,"n_units":eligible_units,"status_detail":detail,"plot":str(plot_path)}])


def export_annotations(adata: ad.AnnData, cohort: str) -> None:
    keep_cols=["cohort","cell_id","unit_id","unit_type","sample_id","condition","pid1_subtype","pid1_cellclass","cellchat_eligible","leiden","total_umi","n_genes","pct_mito","qc_pass","scrublet_score","scrublet_predicted_doublet","marker_lineage_conflict","competing_lineages_detected","PID1_log1p","PID1_detected"]
    if cohort=="GSE200997": keep_cols += [x for x in ["location","msi_status"] if x in adata.obs]
    for col in keep_cols:
        if col not in adata.obs: adata.obs[col]=np.nan
    adata.obs[keep_cols].to_csv(OUT/f"{cohort}_PID1_CellChat_myeloid_annotations.tsv.gz",sep="\t",index=False,compression="gzip")


def cohort_run(cohort: str, matrix: Path, metadata_path: Path, sep: str) -> dict:
    say(f"{cohort}: read original cell annotation")
    meta, tumor_meta, candidate_meta = read_metadata(metadata_path,sep,cohort)
    all_cells_per_sample=meta.groupby("sample_id",observed=True).size().to_dict()
    if cohort=="GSE132465":
        candidate_ids=candidate_meta["cell_id"].tolist()
        screen=None
    else:
        say(f"{cohort}: broad myeloid marker screen across tumor cells only")
        candidate_ids,screen=broad_myeloid_screen(matrix,tumor_meta)
        screen.to_csv(OUT/f"{cohort}_broad_myeloid_marker_screen.tsv.gz",sep="\t",index=False,compression="gzip")
        candidate_meta=candidate_meta[candidate_meta.cell_id.isin(candidate_ids)].copy()
    say(f"{cohort}: {len(candidate_ids)} candidate myeloid cells; read raw counts for candidates")
    adata=read_selected_counts(matrix,sep,candidate_ids,candidate_meta,cohort)
    say(f"{cohort}: raw matrix loaded {adata.n_obs} cells x {adata.n_vars} genes")
    if screen is not None:
        screen_idx = screen.set_index("cell_id")
        aligned = screen_idx.reindex(adata.obs_names)
        for col in ["myeloid_marker_hits", "myeloid_module_score", "strongest_competing_score", "multi_lineage_marker_conflict", "competing_marker_conflict"]:
            if col in aligned:
                adata.obs[f"screen_{col}"] = aligned[col].to_numpy()
    scrub_summary=add_qc_and_doublet_calls(adata,all_cells_per_sample,cohort)
    scrub_summary.to_csv(OUT/f"{cohort}_QC_doublet_by_sample.tsv",sep="\t",index=False)
    audit_cols=["cell_id","unit_id","unit_type","sample_id","condition","original_broad_type","original_subtype",
                "total_umi","n_genes","pct_mito","log10_total_umi","qc_low_library","qc_low_gene_count","qc_high_mito","qc_sample_outlier","qc_pass",
                "scrublet_run","scrublet_score","scrublet_predicted_doublet","marker_lineage_conflict","competing_lineages_detected","exclude_doublet","analysis_keep"]
    audit_cols += [c for c in adata.obs.columns if c.startswith("marker_hits_") or c.startswith("marker_umi_") or c.startswith("strict_marker_") or c.startswith("screen_")]
    for col in audit_cols:
        if col not in adata.obs:
            adata.obs[col] = np.nan
    adata.obs[audit_cols].to_csv(OUT/f"{cohort}_cell_QC_doublet_marker_audit.tsv.gz",sep="\t",index=False,compression="gzip")
    qc_cols=["total_umi","n_genes","pct_mito","qc_low_library","qc_low_gene_count","qc_high_mito","qc_sample_outlier","qc_pass","marker_lineage_conflict","scrublet_predicted_doublet","analysis_keep"]
    qcs=adata.obs.groupby("sample_id",observed=True)[qc_cols].agg({x:("mean" if x in ["qc_low_library","qc_low_gene_count","qc_high_mito","qc_sample_outlier","qc_pass","marker_lineage_conflict","scrublet_predicted_doublet","analysis_keep"] else "median") for x in qc_cols}).reset_index()
    qcs.insert(0,"cohort",cohort);qcs.to_csv(OUT/f"{cohort}_QC_summary_by_sample.tsv",sep="\t",index=False)
    flow=adata.obs.groupby("sample_id",observed=True).agg(
        n_candidate_cells=("cell_id","size"),n_qc_pass=("qc_pass","sum"),n_scrublet_doublets=("scrublet_predicted_doublet","sum"),
        n_marker_conflict_flags=("marker_lineage_conflict","sum"),n_retained=("analysis_keep","sum")).reset_index()
    flow.insert(0,"cohort",cohort)
    for col in ["n_qc_pass","n_scrublet_doublets","n_marker_conflict_flags","n_retained"]:
        flow[col]=flow[col].astype(int)
    flow.to_csv(OUT/f"{cohort}_QC_flow_by_sample.tsv",sep="\t",index=False)
    sc_keep=adata[adata.obs["analysis_keep"].to_numpy()].copy()
    if sc_keep.n_obs<100:
        raise RuntimeError(f"{cohort}: fewer than 100 cells remained after QC/doublet review")
    say(f"{cohort}: {sc_keep.n_obs} retained after QC and doublet review; reclustering")
    clustered=add_log_norm_and_cluster(sc_keep,cohort)
    cluster_tab,marker_cov=annotate_fine_states(clustered,cohort)
    cluster_tab.to_csv(OUT/f"{cohort}_cluster_marker_score_annotation.tsv",sep="\t",index=False)
    marker_cov.to_csv(OUT/f"{cohort}_cluster_program_gene_coverage.tsv",sep="\t",index=False)
    # Cluster marker testing is descriptive and used for annotation audit, not a patient-level inferential test.
    sc.tl.rank_genes_groups(clustered,groupby="leiden",method="wilcoxon",use_raw=False,pts=True,n_genes=min(100,clustered.n_vars))
    markers=sc.get.rank_genes_groups_df(clustered,group=None)
    markers.insert(0,"cohort",cohort)
    markers.to_csv(OUT/f"{cohort}_cluster_rank_markers.tsv.gz",sep="\t",index=False,compression="gzip")
    summary=summarize_unit_expression(clustered,cohort)
    summary.to_csv(OUT/f"{cohort}_patient_sample_PID1_by_subtype.tsv",sep="\t",index=False)
    tests=paired_pid1_tests(summary,cohort)
    tests.to_csv(OUT/f"{cohort}_PID1_unit_level_sensitivity.tsv",sep="\t",index=False)
    reg,program_meta=add_program_scores_and_tests(clustered,cohort)
    reg.to_csv(OUT/f"{cohort}_PID1_monocyte_program_adjusted_models.tsv",sep="\t",index=False)
    program_meta.to_csv(OUT/f"{cohort}_PID1_program_coverage_and_unit_correlations.tsv",sep="\t",index=False)
    ptdf,ptstatus=pseudotime_exploration(clustered,cohort)
    ptdf.to_csv(OUT/f"{cohort}_PID1_pseudotime_association.tsv",sep="\t",index=False)
    ptstatus.to_csv(OUT/f"{cohort}_pseudotime_eligibility.tsv",sep="\t",index=False)
    # Cohort-level state counts; report both distinct units and absolute cells.
    counttab=clustered.obs.groupby("pid1_subtype",observed=True).agg(n_cells=("cell_id","size"),n_units=("unit_id","nunique")).reset_index()
    counttab.insert(0,"cohort",cohort)
    counttab.to_csv(OUT/f"{cohort}_subtype_cell_unit_counts.tsv",sep="\t",index=False)
    plot_umap(clustered,cohort);plot_marker_dot(clustered,cohort);plot_violin_detection(clustered,summary,cohort)
    export_annotations(clustered,cohort)
    # Keep a compact in-memory output object for cross-cohort-independent summary.
    info={"cohort":cohort,"n_annotation_cells":len(meta),"n_tumor_cells":len(tumor_meta),"n_candidate_mye":len(candidate_ids),
          "n_qc_pass":int(adata.obs.qc_pass.sum()),"n_retained":clustered.n_obs,"n_clusters":clustered.obs.leiden.nunique(),
          "n_units":clustered.obs.unit_id.nunique(),"n_scrublet_samples":int(scrub_summary.scrublet_completed.sum()),
          "n_scrublet_failed_or_skipped":int((~scrub_summary.scrublet_completed).sum()),"pseudotime_status":ptstatus.iloc[0].pseudotime_status,
          "hvg_n":int(clustered.uns.get("pid1_hvg_n",np.nan))}
    clustered.write_h5ad(WORK/f"{cohort}_retained.h5ad")
    return {"info":info,"summary":summary,"tests":tests,"counts":counttab,"regression":reg,"program_meta":program_meta,"pseudotime":ptdf,"screen":screen,"cluster_markers":markers}


def main() -> None:
    say(f"Starting; scanpy={sc.__version__}; output={OUT}")
    results=[]
    results.append(cohort_run("GSE132465",G132_MATRIX,G132_META,"\t"))
    results.append(cohort_run("GSE200997",G200_MATRIX,G200_META,","))
    pd.DataFrame([r["info"] for r in results]).to_csv(OUT/"cohort_qc_summary.tsv",sep="\t",index=False)
    pd.concat([r["summary"] for r in results],ignore_index=True).to_csv(OUT/"PID1_patient_sample_subtype_summary_all_cohorts.tsv",sep="\t",index=False)
    pd.concat([r["tests"] for r in results],ignore_index=True).to_csv(OUT/"PID1_unit_level_threshold_sensitivity_all_cohorts.tsv",sep="\t",index=False)
    pd.concat([r["counts"] for r in results],ignore_index=True).to_csv(OUT/"PID1_subtype_cell_unit_counts_all_cohorts.tsv",sep="\t",index=False)
    pd.concat([r["regression"] for r in results],ignore_index=True).to_csv(OUT/"PID1_monocyte_program_adjusted_models_all_cohorts.tsv",sep="\t",index=False)
    pd.concat([r["program_meta"] for r in results],ignore_index=True).to_csv(OUT/"PID1_program_coverage_and_unit_correlations_all_cohorts.tsv",sep="\t",index=False)
    pd.concat([r["pseudotime"] for r in results],ignore_index=True).to_csv(OUT/"PID1_pseudotime_association_all_cohorts.tsv",sep="\t",index=False)
    pd.concat([r["cluster_markers"] for r in results],ignore_index=True).to_csv(OUT/"PID1_cluster_rank_markers_all_cohorts.tsv.gz",sep="\t",index=False,compression="gzip")
    manifest={"run_name":RUN_NAME,"date":"2026-10-11","inputs":"Original GEO processed raw UMI matrices and original cell annotations only; no prior analysis results or caches read.",
              "outputs":str(OUT),"scanpy":sc.__version__,"pid1_excluded_from_clustering_and_programs":True,
              "gse200997_candidate_screen":"Tumor-only screen on original raw UMI marker rows; require >=2 core myeloid genes, LYZ plus another myeloid-specific support marker or >=3 core myeloid genes, and a myeloid log1p-marker module score >= every competing lineage module score. PTPRC alone is not sufficient.",
              "qc":"UMI>=500, genes>=200, mitochondrial<=20%; sample-specific >4 MAD outlier filter; Scrublet within sample when >=80 QC-passing myeloid candidates; cross-lineage marker conflict screen.",
              "doublets":"Scrublet 0.2.3 on QC-passing myeloid candidates, exact sklearn-neighbor mode (use_approx_neighbors=False), sample-adaptive expected rate 0.8% per 1,000 captured cells clipped to 1-10%; multi-lineage marker co-detection is exported as an audit flag; clusters dominated by strong competing-lineage markers are annotated as contamination and excluded from monocyte/macrophage contrasts.",
              "cluster":"2000 seurat HVGs, sample-batch-aware when possible; PCA, neighbor graph, Leiden resolution 0.6; each cohort independently clustered. Annotations use PID1-free monocyte/macrophage/DC marker programs, strong competing-lineage clusters are flagged, and cluster marker tables are provided for audit.",
              "gene_programs":PROGRAMS,"paper_audit":{"Nature_Cancer":"https://www.nature.com/articles/s43018-026-01189-0","PubMed":"https://pubmed.ncbi.nlm.nih.gov/42303776/"}}
    (OUT/"analysis_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    say(f"Completed both cohorts; outputs={OUT}")


if __name__=="__main__":
    main()
