#!/usr/bin/env python
"""Map SCD and PID1 to annotated CRC single-cell subpopulations.

Only the raw GEO count matrices and their matching cell annotations are read.
For GSE132465 the author-provided epithelial/stromal subtypes are used. For
GSE200997 broad lineages are inferred from marker expression and fibroblast
states are scored from predefined CAF marker modules. Cells are summarized by
sample before inferential tests; individual cells are not treated as replicates.
"""
from __future__ import annotations

import gzip
import itertools
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "work" / "single_cell_GSE132465_GSE200997"
OUT = ROOT / "results" / "2026-10-11" / "SCD_PID1_cell_subclusters"
OUT.mkdir(parents=True, exist_ok=True)

COHORTS = {
    "GSE132465": {
        "matrix": "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
        "annotation": "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
        "delimiter": "\t",
        "cell_key": "Index",
        "sample_key": "Sample",
        "condition_key": "Class",
        "annotation_sep": "\t",
    },
    "GSE200997": {
        "matrix": "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz",
        "annotation": "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz",
        "delimiter": ",",
        "cell_key": "Unnamed: 0",
        "sample_key": "samples",
        "condition_key": "Condition",
        "annotation_sep": ",",
    },
}

LINEAGE_MARKERS = {
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

# A priori lipid-metabolism panel. SCD is reported separately and excluded from
# the module score so that the SCD result does not tautologically drive it.
LIPID_PANEL = [
    "ACACA", "ACACB", "ACLY", "ACSS2", "FASN", "SREBF1", "SREBF2",
    "FADS1", "FADS2", "FADS3", "SCD5", "ELOVL1", "ELOVL2", "ELOVL3",
    "ELOVL5", "ELOVL6", "ELOVL7", "ACSL1", "ACSL3", "ACSL4", "ACSL5",
    "ACSL6", "CPT1A", "CPT2", "ACADVL", "ACOX1", "HADHA", "HADHB",
    "CD36", "FABP1", "FABP2", "FABP5", "DGAT1", "DGAT2", "GPAT3",
    "GPAT4", "HMGCR", "FDFT1", "SQLE", "INSIG1", "INSIG2",
]

# Standard CAF marker modules. Scores describe transcriptional programs; they
# are not presented as experimentally validated cell labels.
CAF_MODULES = {
    "myCAF": ["ACTA2", "TAGLN", "MYL9", "TPM2", "POSTN", "CNN1", "ACTG2"],
    "iCAF": ["IL6", "CXCL12", "CXCL14", "CCL2", "CFD", "PTGDS", "PI16", "IL1RL1"],
    "ECM_CAF": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "FN1", "TNC", "FAP"],
    "apCAF": ["CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DRB1", "CIITA"],
}

ALL_GENES = set(["SCD", "PID1"] + LIPID_PANEL)
for vals in LINEAGE_MARKERS.values():
    ALL_GENES.update(vals)
for vals in CAF_MODULES.values():
    ALL_GENES.update(vals)


def read_annotation(cohort: str, cfg: dict) -> pd.DataFrame:
    path = INPUT / cfg["annotation"]
    ann = pd.read_csv(path, sep=cfg["annotation_sep"], dtype=str, keep_default_na=False)
    key = cfg["cell_key"]
    if key not in ann.columns:
        raise ValueError(f"Missing cell key {key!r} in {path.name}: {ann.columns.tolist()}")
    ann = ann.set_index(key, drop=False)
    if ann.index.duplicated().any():
        raise ValueError(f"Duplicate cell IDs in {path.name}")
    return ann


def stream_selected_matrix(cohort: str, cfg: dict):
    """Read the GEO genes-by-cells matrix, retaining only panel genes."""
    path = INPUT / cfg["matrix"]
    cache_path = OUT / f"_temporary_panel_cache_{cohort}.npz"
    if cache_path.exists():
        with np.load(cache_path, allow_pickle=False) as cached:
            cells = cached["cells"].astype(str).tolist()
            libraries = cached["libraries"]
            features = int(cached["features"][0])
            elapsed = float(cached["elapsed"][0])
            duplicates = json.loads(str(cached["duplicates"][0]))
            counts = {name[6:]: cached[name] for name in cached.files if name.startswith("count__")}
        print(f"{cohort}: reused current-run panel cache", flush=True)
        return cells, libraries, counts, features, duplicates, elapsed
    delim = cfg["delimiter"].encode("ascii")
    t0 = time.time()
    selected: dict[str, np.ndarray] = {}
    duplicate_rows: dict[str, int] = {}
    with gzip.open(path, "rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(delim)
        cells = [x.decode("utf-8").strip().strip('"') for x in header[1:]]
        if len(cells) != len(set(cells)):
            raise ValueError(f"Duplicate cell IDs in {path.name}")
        libraries = np.zeros(len(cells), dtype=np.float64)
        features = 0
        for row_no, raw in enumerate(handle, start=2):
            raw = raw.rstrip(b"\r\n")
            if not raw:
                continue
            pieces = raw.split(delim, 1)
            if len(pieces) != 2:
                raise ValueError(f"Malformed matrix row {row_no} in {path.name}")
            gene = pieces[0].decode("utf-8").strip().strip('"')
            vals = np.fromstring(pieces[1].decode("ascii"), sep=cfg["delimiter"], dtype=np.float64)
            if vals.size != len(cells):
                raise ValueError(f"Row {row_no} ({gene}) has {vals.size}, expected {len(cells)}")
            libraries += vals
            features += 1
            if gene in ALL_GENES:
                if gene in selected:
                    selected[gene] += vals
                    duplicate_rows[gene] = duplicate_rows.get(gene, 0) + 1
                else:
                    selected[gene] = vals.copy()
            if features % 2500 == 0:
                print(f"{cohort}: scanned {features} genes; elapsed {time.time()-t0:.1f}s", flush=True)
    if np.any(libraries <= 0) or not np.all(np.isfinite(libraries)):
        raise ValueError(f"Invalid library sizes in {cohort}")
    np.savez_compressed(
        cache_path,
        cells=np.asarray(cells, dtype=str),
        libraries=libraries,
        features=np.asarray([features], dtype=np.int64),
        elapsed=np.asarray([time.time() - t0], dtype=np.float64),
        duplicates=np.asarray([json.dumps(duplicate_rows)], dtype=str),
        **{f"count__{gene}": vals.astype(np.float32) for gene, vals in selected.items()},
    )
    return cells, libraries, selected, features, duplicate_rows, time.time() - t0


def normalized_log_expression(libraries: np.ndarray, counts: dict[str, np.ndarray]):
    scale = 10000.0 / libraries
    return {gene: np.log1p(vals * scale).astype(np.float32) for gene, vals in counts.items()}


def mean_score(expr: dict[str, np.ndarray], genes: list[str], n_cells: int):
    present = [g for g in genes if g in expr]
    if not present:
        return np.zeros(n_cells, dtype=np.float32), present
    return np.vstack([expr[g] for g in present]).mean(axis=0), present


def gene_zscores(expr: dict[str, np.ndarray], genes: list[str], mask: np.ndarray):
    present = [g for g in genes if g in expr]
    z = {}
    for g in present:
        vals = expr[g][mask].astype(np.float64)
        sd = vals.std(ddof=0)
        z[g] = ((expr[g] - vals.mean()) / sd).astype(np.float32) if sd > 0 else np.zeros_like(expr[g])
    return z, present


def module_zscore(zgenes: dict[str, np.ndarray], genes: list[str], n_cells: int):
    present = [g for g in genes if g in zgenes]
    if not present:
        return np.zeros(n_cells, dtype=np.float32), present
    return np.vstack([zgenes[g] for g in present]).mean(axis=0), present


def infer_broad_lineages(expr: dict[str, np.ndarray], n_cells: int):
    scores = {}
    detected = {}
    for name, genes in LINEAGE_MARKERS.items():
        present = [g for g in genes if g in expr]
        if present:
            mat = np.vstack([expr[g] for g in present])
            scores[name] = mat.mean(axis=0)
            detected[name] = (mat > 0).sum(axis=0)
        else:
            scores[name] = np.zeros(n_cells)
            detected[name] = np.zeros(n_cells, dtype=int)
    names = list(scores)
    score_mat = np.vstack([scores[k] for k in names])
    order = np.argsort(score_mat, axis=0)
    best = order[-1, :]
    best_score = score_mat[best, np.arange(n_cells)]
    second_score = score_mat[order[-2, :], np.arange(n_cells)]
    best_name = np.array(names, dtype=object)[best]
    best_n = np.array([detected[n] for n in names])[best, np.arange(n_cells)]
    assigned = np.where((best_score > 0) & (best_n >= 2), best_name, "Unassigned")
    confident = (assigned != "Unassigned") & ((best_score - second_score) >= 0.05)
    return assigned, best_score, best_score - second_score, confident, scores


def summarize_cells(cohort: str, cfg: dict, ann: pd.DataFrame, cells: list[str],
                    libraries: np.ndarray, counts: dict[str, np.ndarray], features: int,
                    duplicates: dict[str, int], elapsed: float):
    # Match the GEO annotation to the matrix column order with explicit checks.
    missing = sorted(set(cells) - set(ann.index))
    extra = sorted(set(ann.index) - set(cells))
    if missing or extra:
        raise ValueError(f"Cell ID mismatch in {cohort}: missing annotation={len(missing)}, extra={len(extra)}")
    ann = ann.loc[cells].copy()
    n = len(cells)
    expr = normalized_log_expression(libraries, counts)
    condition = ann[cfg["condition_key"]].to_numpy()
    sample = ann[cfg["sample_key"]].to_numpy()
    result = pd.DataFrame({
        "cohort": cohort,
        "cell_id": cells,
        "sample": sample,
        "condition": condition,
        "library_size_UMI": libraries,
    })

    # Scores are returned for inspection and used only to identify broad lineages.
    broad, broad_score, broad_margin, broad_confident, broad_scores = infer_broad_lineages(expr, n)
    result["inferred_broad_lineage"] = broad
    result["inferred_lineage_score"] = broad_score
    result["inferred_lineage_margin"] = broad_margin
    result["inferred_lineage_confident"] = broad_confident

    if cohort == "GSE132465":
        cell_type = ann["Cell_type"].to_numpy()
        source_subtype = ann["Cell_subtype"].to_numpy()
        epi_mask = (condition == "Tumor") & (cell_type == "Epithelial cells") & np.isin(source_subtype, ["CMS1", "CMS2", "CMS3", "CMS4"])
        fibro_labels = ["Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3"]
        fibro_mask = (condition == "Tumor") & (cell_type == "Stromal cells") & np.isin(source_subtype, fibro_labels)
        result["source_cell_type"] = cell_type
        result["source_subtype"] = source_subtype
        result["compartment"] = np.where(epi_mask, "Tumor_epithelial", np.where(fibro_mask, "Tumor_fibroblast", "Other"))
        result["subtype"] = np.where(epi_mask, source_subtype, np.where(fibro_mask, source_subtype, ""))
        result["subtype_source"] = np.where(epi_mask | fibro_mask, "GEO original Cell_type/Cell_subtype", "")
    else:
        # The GSE200997 per-cell prediction field is a CMS label, not a broad
        # cell lineage. Lineages are inferred from raw marker expression.
        pred = ann["prediction"].to_numpy() if "prediction" in ann.columns else np.full(n, "NA")
        epi_mask = (condition == "Tumor") & (broad == "Epithelial") & broad_confident
        fibro_mask = (condition == "Tumor") & (broad == "Fibroblast") & broad_confident
        result["source_cell_type"] = broad
        result["source_subtype"] = pred
        result["compartment"] = np.where(epi_mask, "Tumor_epithelial", np.where(fibro_mask, "Tumor_fibroblast", "Other"))
        result["subtype"] = ""
        result.loc[epi_mask, "subtype"] = np.where(np.isin(pred[epi_mask], ["CMS1", "CMS2", "CMS3", "CMS4"]), pred[epi_mask], "CMS_unassigned")
        result.loc[epi_mask, "subtype_source"] = "GEO per-cell CMS prediction + marker-inferred epithelial lineage"
        result.loc[fibro_mask, "subtype_source"] = "Marker-inferred fibroblast lineage; CAF program assigned by module score"
        result.loc[fibro_mask, "subtype"] = "CAF_state_pending"

    # Lipid program: within-cohort gene-wise z scores among tumor epithelial cells;
    # SCD itself is excluded from this module score.
    lipid_z, lipid_present = gene_zscores(expr, LIPID_PANEL, epi_mask)
    lipid_module, lipid_used = module_zscore(lipid_z, LIPID_PANEL, n)
    result["lipid_module_z_excluding_SCD"] = lipid_module
    result["SCD_log1p_10K"] = expr.get("SCD", np.zeros(n, dtype=np.float32))
    result["SCD_detected"] = (result["SCD_log1p_10K"].to_numpy() > 0)

    # CAF states are assigned only in the marker-inferred GSE200997 fibroblasts.
    caf_scores = {}
    caf_genes_used = {}
    if cohort == "GSE200997":
        fibro_z, _ = gene_zscores(expr, sorted(set(sum(CAF_MODULES.values(), []))), fibro_mask)
        for state, genes in CAF_MODULES.items():
            score, used = module_zscore(fibro_z, genes, n)
            caf_scores[state] = score
            caf_genes_used[state] = used
            result[f"CAF_{state}_z"] = score
        if caf_scores:
            names = list(caf_scores)
            score_mat = np.vstack([caf_scores[k] for k in names])
            assignment = np.array(names, dtype=object)[np.argmax(score_mat, axis=0)]
            result.loc[fibro_mask, "subtype"] = assignment[fibro_mask]
    else:
        # Carry CAF marker scores in the annotated stromal subtypes as supportive
        # descriptors; the author labels remain the primary subtype definition.
        fib_z, _ = gene_zscores(expr, sorted(set(sum(CAF_MODULES.values(), []))), fibro_mask)
        for state, genes in CAF_MODULES.items():
            score, used = module_zscore(fib_z, genes, n)
            caf_scores[state] = score
            caf_genes_used[state] = used
            result[f"CAF_{state}_z"] = score

    result["PID1_log1p_10K"] = expr.get("PID1", np.zeros(n, dtype=np.float32))
    result["PID1_detected"] = (result["PID1_log1p_10K"].to_numpy() > 0)

    target = result[result["compartment"].isin(["Tumor_epithelial", "Tumor_fibroblast"])].copy()
    # Retain the broad lineage evidence for GSE200997 and original labels for GSE132465.
    if cohort == "GSE200997":
        target["CMS_prediction"] = ann.loc[target["cell_id"], "prediction"].to_numpy()
    else:
        target["CMS_prediction"] = ""
    target["n_lipid_panel_genes_available"] = len(lipid_present)
    target["n_CAF_module_genes_available"] = sum(len(x) for x in caf_genes_used.values())

    metadata = {
        "cohort": cohort,
        "matrix_file": cfg["matrix"],
        "annotation_file": cfg["annotation"],
        "cells": n,
        "matrix_features": features,
        "target_genes_found": sorted(counts),
        "duplicate_target_gene_rows": duplicates,
        "elapsed_seconds": round(elapsed, 2),
        "lipid_panel_found": lipid_present,
        "lipid_panel_used": lipid_used,
        "caf_module_genes_found": caf_genes_used,
    }
    return target, result, metadata


def paired_subtype_tests(sample_means: pd.DataFrame):
    rows = []
    metrics = {
        "Tumor_epithelial": ["SCD_log1p_10K", "lipid_module_z_excluding_SCD"],
        "Tumor_fibroblast": ["PID1_log1p_10K"] + [f"CAF_{x}_z" for x in CAF_MODULES],
    }
    for (cohort, compartment), group in sample_means.groupby(["cohort", "compartment"], sort=True):
        metric_list = metrics.get(compartment, [])
        subtypes = sorted(group["subtype"].unique())
        for metric in metric_list:
            if metric not in group.columns:
                continue
            for a, b in itertools.combinations(subtypes, 2):
                ga = group[(group.subtype == a) & (group.n_cells >= 10)].set_index("sample")[metric]
                gb = group[(group.subtype == b) & (group.n_cells >= 10)].set_index("sample")[metric]
                common = ga.index.intersection(gb.index)
                if len(common) < 3:
                    rows.append({"cohort": cohort, "compartment": compartment, "metric": metric,
                                 "subtype_A": a, "subtype_B": b, "n_paired_samples": len(common),
                                 "median_paired_difference_A_minus_B": np.nan, "p_value": np.nan,
                                 "test_status": "<3 paired samples with >=10 cells per subtype"})
                    continue
                diff = ga.loc[common].to_numpy() - gb.loc[common].to_numpy()
                try:
                    p = float(wilcoxon(diff, zero_method="wilcox", alternative="two-sided").pvalue)
                    status = "paired Wilcoxon on sample means"
                except ValueError:
                    p = 1.0
                    status = "all paired differences zero"
                rows.append({"cohort": cohort, "compartment": compartment, "metric": metric,
                             "subtype_A": a, "subtype_B": b, "n_paired_samples": len(common),
                             "median_paired_difference_A_minus_B": float(np.median(diff)),
                             "p_value": p, "test_status": status})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["BH_FDR_within_cohort_compartment_metric"] = np.nan
        valid = out.p_value.notna()
        for _, idx in out.loc[valid].groupby(["cohort", "compartment", "metric"]).groups.items():
            vals = out.loc[idx, "p_value"].to_numpy(float)
            order = np.argsort(vals)
            q = np.empty(len(vals), dtype=float)
            ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
            ranked = np.minimum.accumulate(ranked[::-1])[::-1]
            q[order] = np.minimum(ranked, 1.0)
            out.loc[idx, "BH_FDR_within_cohort_compartment_metric"] = q
    return out


def main():
    all_cells = []
    metadata = []
    for cohort, cfg in COHORTS.items():
        print(f"Reading {cohort} annotations and raw matrix", flush=True)
        ann = read_annotation(cohort, cfg)
        cells, libs, counts, features, duplicates, elapsed = stream_selected_matrix(cohort, cfg)
        selected, full, meta = summarize_cells(cohort, cfg, ann, cells, libs, counts, features, duplicates, elapsed)
        all_cells.append(selected)
        metadata.append(meta)
        print(f"{cohort}: {len(selected)} tumor epithelial/fibroblast cells retained", flush=True)

    cells_df = pd.concat(all_cells, ignore_index=True)
    cells_df.to_csv(OUT / "cell_subpopulation_expression.tsv.gz", sep="\t", index=False, compression="gzip")

    group_cols = ["cohort", "compartment", "subtype", "subtype_source"]
    summary_rows = []
    sample_rows = []
    for key, g in cells_df.groupby(group_cols, dropna=False, sort=True):
        cohort, compartment, subtype, source = key
        rec = {"cohort": cohort, "compartment": compartment, "subtype": subtype,
               "subtype_source": source, "n_cells": len(g), "n_samples": g["sample"].nunique()}
        if compartment == "Tumor_epithelial":
            rec.update({"SCD_mean_log1p_10K": g.SCD_log1p_10K.mean(),
                        "SCD_pct_detected": 100 * g.SCD_detected.mean(),
                        "lipid_module_mean_z_excluding_SCD": g.lipid_module_z_excluding_SCD.mean(),
                        "lipid_module_median_z_excluding_SCD": g.lipid_module_z_excluding_SCD.median()})
        else:
            rec.update({"PID1_mean_log1p_10K": g.PID1_log1p_10K.mean(),
                        "PID1_pct_detected": 100 * g.PID1_detected.mean()})
            for state in CAF_MODULES:
                col = f"CAF_{state}_z"
                if col in g.columns:
                    rec[f"{col}_mean"] = g[col].mean()
        summary_rows.append(rec)
        for sample_name, sg in g.groupby("sample", sort=True):
            srec = {"cohort": cohort, "sample": sample_name, "compartment": compartment,
                    "subtype": subtype, "subtype_source": source, "n_cells": len(sg)}
            for metric in ["SCD_log1p_10K", "SCD_detected", "lipid_module_z_excluding_SCD",
                           "PID1_log1p_10K", "PID1_detected"] + [f"CAF_{x}_z" for x in CAF_MODULES]:
                if metric in sg.columns:
                    srec[metric] = sg[metric].mean()
            sample_rows.append(srec)

    summary = pd.DataFrame(summary_rows)
    samples = pd.DataFrame(sample_rows)
    summary.to_csv(OUT / "subpopulation_expression_summary.tsv", sep="\t", index=False)
    samples.to_csv(OUT / "sample_subpopulation_means.tsv", sep="\t", index=False)
    tests = paired_subtype_tests(samples)
    tests.to_csv(OUT / "paired_subtype_tests.tsv", sep="\t", index=False)
    with (OUT / "analysis_metadata.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    lines = [
        "# SCD / PID1 单细胞亚群定位分析",
        "",
        "## 输入与标签来源",
        "",
        "- 只读取本项目 `work/single_cell_GSE132465_GSE200997/` 中的两套 GEO 原始 UMI 矩阵和配套注释；未读取历史结果文件。",
        "- 表达量按每细胞 UMI 总量归一到 10,000 后 log1p。",
        "- GSE132465 使用 GEO 原始 `Cell_type` / `Cell_subtype` 注释。上皮以 CMS1–4 分层；成纤维细胞限于 stromal 中的 Myofibroblasts、Stromal 1、2、3。原注释的内皮、周细胞和平滑肌不当作成纤维细胞。",
        "- GSE200997 没有细胞类型字段。按上皮/成纤维等经典 marker 均值推断宽泛谱系，仅保留 marker 判定可信的肿瘤细胞；上皮 CMS 使用 GEO 的逐细胞 `prediction` 字段。成纤维细胞按 myCAF、iCAF、ECM_CAF、apCAF 预设 marker 模块中最高的基因标准分数赋予描述性状态。该状态是模块推断，不等同于作者验证的细胞标签或独立聚类。",
        "- 脂质代谢分数为预先列出的脂质合成、去饱和、延长、转运、氧化与储存基因的跨细胞标准分数均值，并排除 SCD 本身。",
        "",
        "## 统计单位",
        "",
        "细胞级结果用于描述；亚群间检验先按样本汇总每个亚群的均值，再对同一肿瘤样本内的亚群做配对 Wilcoxon 检验。每个亚群/样本至少 10 个细胞，且至少 3 个配对样本才检验；P 值在同一队列、区室、指标的亚群两两比较内做 BH 校正。",
        "",
        "## 结果文件",
        "",
        "- `subpopulation_expression_summary.tsv`：每亚群细胞数、覆盖样本数、SCD/PID1 表达和模块分数。",
        "- `sample_subpopulation_means.tsv`：用于复核独立样本层面的均值。",
        "- `paired_subtype_tests.tsv`：样本配对检验。",
        "- `cell_subpopulation_expression.tsv.gz`：纳入分析的细胞及其亚群、目标表达和模块分数。",
        "- `analysis_metadata.json`：原始输入、矩阵行数、面板基因覆盖情况。",
        "",
        "## 解读限制",
        "",
        "GSE132465 的原始亚群标签与 GSE200997 的推断标签来源不同，因此跨队列一致性仅作为方向性支持。GSE200997 的 CAF 状态标签是预设 marker 模块评分，不应写成无监督聚类结果。低细胞数亚群的表达与检验需谨慎解释。",
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for cohort in COHORTS:
        (OUT / f"_temporary_panel_cache_{cohort}.npz").unlink(missing_ok=True)
    print(f"Wrote outputs to {OUT}", flush=True)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
