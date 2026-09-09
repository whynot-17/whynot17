from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import mannwhitneyu, wilcoxon


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs" / "DINP_CRC_multicohort_scRNA_CEBPB_CD36"
TARGETS = ["CEBPB", "CD36"]
MIN_CELLS_PER_SAMPLE_CELLTYPE = 50
MIN_CROSS_CELLS_PER_SAMPLE_COMPARTMENT = 10
MIN_QC_GENES = 200
MIN_QC_UMI = 500

GSE132465_REPO_OUTPUT = ROOT / "analysis" / "dinp_crc_scRNA_GSE132465_CEBPB_CD36"
GSE132465_LOCAL_OUTPUT = ROOT / "outputs" / "DINP_CRC_scRNA_GSE132465_CEBPB_CD36"
GSE132465_CELL_LEVEL = (
    GSE132465_REPO_OUTPUT / "GSE132465_CEBPB_CD36_cell_level_expression.csv"
    if (GSE132465_REPO_OUTPUT / "GSE132465_CEBPB_CD36_cell_level_expression.csv").exists()
    else GSE132465_LOCAL_OUTPUT / "GSE132465_CEBPB_CD36_cell_level_expression.csv"
)
GSE132465_LOCALIZATION = ROOT / "outputs" / "DINP_CRC_scRNA_GSE132465_CEBPB_CD36" / "GSE132465_CEBPB_CD36_cell_type_localization.csv"
GSE132465_ANNOTATION = ROOT / "work" / "scRNA_GSE132465" / "cell_annotation.txt.gz"
GSE132465_SOFT = ROOT / "work" / "scRNA_GSE132465" / "GSE132465_family.soft.gz"

GSE200997_DIR = ROOT / "work" / "scRNA_GSE200997"
GSE200997_ANNOTATION = GSE200997_DIR / "cell_annotation.csv.gz"
GSE200997_COUNTS = GSE200997_DIR / "raw_UMI_count_matrix.csv.gz"

GSE188711_DIR = ROOT / "work" / "scRNA_GSE188711" / "raw"

MARKER_PANELS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CEACAM5", "MUC1"],
    "T_NK": ["CD3D", "CD3E", "TRAC", "CD247", "CD8A", "NKG7", "GNLY", "TRBC1"],
    "B_Plasma": ["CD79A", "MS4A1", "CD37", "CD74", "CD79B", "CD19", "JCHAIN", "MZB1", "SDC1"],
    "Myeloid": ["LYZ", "LST1", "TYROBP", "FCER1G", "CTSS", "CD68", "C1QA", "C1QB", "C1QC"],
    "Neutrophil": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "S100A12", "LILRB3", "CXCR2"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "THY1", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "ENG", "CLDN5", "RAMP2"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "CPA3", "GATA2", "MS4A2"],
}
MACROPHAGE_MARKERS = ["C1QA", "C1QB", "C1QC", "APOE", "CD68", "MSR1", "MRC1", "VSIG4"]
MONOCYTE_MARKERS = ["LYZ", "FCN1", "S100A8", "S100A9", "FCGR3A", "CTSD", "LILRB1", "CTSS"]
ALL_MARKERS = sorted(set(gene for genes in MARKER_PANELS.values() for gene in genes) | set(MACROPHAGE_MARKERS) | set(MONOCYTE_MARKERS) | set(TARGETS))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_path(path: Path) -> str:
    """Return portable repository-relative paths in the audit manifest."""
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def bh_adjust(values: pd.Series) -> pd.Series:
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    result = np.full(array.shape, np.nan, dtype=float)
    valid = np.flatnonzero(np.isfinite(array))
    if len(valid) == 0:
        return pd.Series(result, index=values.index)
    order = valid[np.argsort(array[valid])]
    adjusted = array[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result[order] = np.minimum(adjusted, 1.0)
    return pd.Series(result, index=values.index)


def safe_wilcoxon(left: np.ndarray, right: np.ndarray) -> tuple[float, float]:
    differences = right - left
    if len(differences) < 3:
        return np.nan, np.nan
    if np.allclose(differences, 0):
        return 0.0, 1.0
    test = wilcoxon(right, left, zero_method="wilcox", alternative="two-sided", method="auto")
    return float(test.statistic), float(test.pvalue)


def safe_mannwhitney(left: np.ndarray, right: np.ndarray) -> tuple[float, float]:
    if len(left) < 2 or len(right) < 2:
        return np.nan, np.nan
    test = mannwhitneyu(left, right, alternative="two-sided", method="auto")
    return float(test.statistic), float(test.pvalue)


def add_target_metrics(cells: pd.DataFrame, counts: dict[str, np.ndarray], total_umi: np.ndarray, n_genes: np.ndarray) -> pd.DataFrame:
    output = cells.copy()
    output["total_umi"] = total_umi.astype(np.int64)
    output["n_genes"] = n_genes.astype(np.int32)
    output["qc_pass"] = (output["total_umi"] >= MIN_QC_UMI) & (output["n_genes"] >= MIN_QC_GENES)
    denominator = output["total_umi"].to_numpy(dtype=float)
    for gene in TARGETS:
        values = counts[gene].astype(np.int32)
        output[f"{gene}_umi"] = values
        output[f"{gene}_detected"] = values > 0
        output[f"{gene}_cpm"] = np.divide(values * 1_000_000.0, denominator, out=np.zeros(len(output), dtype=float), where=denominator > 0)
        output[f"{gene}_log1p_cpm"] = np.log1p(output[f"{gene}_cpm"])
    return output


def infer_marker_cell_types(cells: pd.DataFrame, marker_counts: dict[str, np.ndarray]) -> pd.DataFrame:
    output = cells.copy()
    denominator = output["total_umi"].to_numpy(dtype=float)
    scores = {}
    for panel, genes in MARKER_PANELS.items():
        available = [gene for gene in genes if gene in marker_counts]
        if not available:
            scores[panel] = np.zeros(len(output), dtype=float)
            continue
        panel_values = []
        for gene in available:
            cpm = np.divide(marker_counts[gene].astype(float) * 1_000_000.0, denominator, out=np.zeros(len(output), dtype=float), where=denominator > 0)
            panel_values.append(np.log1p(cpm))
        scores[panel] = np.vstack(panel_values).mean(axis=0)
        output[f"marker_score_{panel}"] = scores[panel]
    score_frame = pd.DataFrame(scores, index=output.index)
    output["marker_top_score"] = score_frame.max(axis=1)
    output["cell_type"] = score_frame.idxmax(axis=1)
    output.loc[output["marker_top_score"].le(0), "cell_type"] = "Unclassified"
    output["marker_score_margin"] = score_frame.apply(lambda row: row.nlargest(2).iloc[0] - row.nlargest(2).iloc[1] if len(row) >= 2 else row.iloc[0], axis=1)
    output["myeloid_subtype"] = ""
    myeloid_mask = output["cell_type"].eq("Myeloid")
    macro_available = [gene for gene in MACROPHAGE_MARKERS if gene in marker_counts]
    mono_available = [gene for gene in MONOCYTE_MARKERS if gene in marker_counts]
    if macro_available and mono_available:
        macro_values = []
        mono_values = []
        for gene in macro_available:
            cpm = np.divide(marker_counts[gene].astype(float) * 1_000_000.0, denominator, out=np.zeros(len(output), dtype=float), where=denominator > 0)
            macro_values.append(np.log1p(cpm))
        for gene in mono_available:
            cpm = np.divide(marker_counts[gene].astype(float) * 1_000_000.0, denominator, out=np.zeros(len(output), dtype=float), where=denominator > 0)
            mono_values.append(np.log1p(cpm))
        macro_score = np.vstack(macro_values).mean(axis=0)
        mono_score = np.vstack(mono_values).mean(axis=0)
        output.loc[myeloid_mask & (macro_score >= mono_score), "myeloid_subtype"] = "Macrophage-like"
        output.loc[myeloid_mask & (mono_score > macro_score), "myeloid_subtype"] = "Monocyte-like"
        output["marker_score_Macrophage_like"] = macro_score
        output["marker_score_Monocyte_like"] = mono_score
    output["epithelial_context"] = ""
    output.loc[output["cell_type"].eq("Epithelial") & output["Class"].eq("Tumor"), "epithelial_context"] = "Tumor epithelial proxy"
    output.loc[output["cell_type"].eq("Epithelial") & output["Class"].eq("Normal"), "epithelial_context"] = "Normal epithelial"
    return output


def extract_dense_matrix(path: Path, delimiter: bytes = b",") -> tuple[list[str], dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray, np.ndarray]:
    targets: dict[str, np.ndarray] = {}
    markers: dict[str, np.ndarray] = {}
    with gzip.open(path, "rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(delimiter)
        cell_ids = [value.decode("utf-8").strip('"') for value in header[1:]]
        total_umi = np.zeros(len(cell_ids), dtype=np.int64)
        n_genes = np.zeros(len(cell_ids), dtype=np.int32)
        wanted = set(ALL_MARKERS)
        for raw_line in handle:
            line = raw_line.rstrip(b"\r\n")
            gene_bytes, payload = line.split(delimiter, 1)
            gene = gene_bytes.decode("utf-8").strip('"').upper()
            values = np.fromstring(payload, dtype=np.int64, sep=delimiter.decode("utf-8"))
            if len(values) != len(cell_ids):
                raise AssertionError(f"Unexpected values length for {gene}: {len(values)}")
            total_umi += values
            n_genes += values > 0
            if gene in wanted:
                if gene in TARGETS:
                    targets[gene] = values.astype(np.int32)
                markers[gene] = values.astype(np.int32)
    missing_targets = set(TARGETS) - set(targets)
    if missing_targets:
        raise AssertionError(f"Missing target rows in {path.name}: {sorted(missing_targets)}")
    return cell_ids, targets, markers, total_umi, n_genes


def load_gse132465() -> tuple[pd.DataFrame, dict]:
    cells = pd.read_csv(GSE132465_CELL_LEVEL, dtype={"cell_id": str, "Patient": str, "Class": str, "Sample": str, "Cell_type": str, "Cell_subtype": str}).fillna("")
    cells = cells.rename(columns={"Cell_type": "cell_type", "Cell_subtype": "cell_subtype"})
    cells["cohort"] = "GSE132465"
    cells["Side"] = ""
    region_by_sample = {}
    current = {}
    with gzip.open(GSE132465_SOFT, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("!Sample_title = "):
                if current:
                    region_by_sample[current["sample"]] = current.get("region", "")
                current = {"sample": line.split("=", 1)[1].strip()}
            elif current and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split("=", 1)[1].strip()
                if ":" in value:
                    key, val = value.split(":", 1)
                    current[key.strip()] = val.strip()
        if current:
            region_by_sample[current["sample"]] = current.get("region", "")
    left_regions = {"rectum", "sigmoid", "rectosigmoid"}
    right_regions = {"ascending", "hepatic flexure"}
    cells["region"] = cells["Sample"].map(region_by_sample).fillna("")
    cells.loc[cells["region"].isin(left_regions), "Side"] = "Left"
    cells.loc[cells["region"].isin(right_regions), "Side"] = "Right"
    patient_side = cells.loc[cells["Side"].ne(""), ["Patient", "Side"]].drop_duplicates("Patient").set_index("Patient")["Side"]
    cells.loc[cells["Side"].eq(""), "Side"] = cells.loc[cells["Side"].eq(""), "Patient"].map(patient_side).fillna("").to_numpy()
    cells["cell_subtype"] = cells["cell_subtype"].astype(str)
    cells["cell_type_source"] = "GEO curated"
    cells["cross_compartment"] = cells["cell_type"].map({
        "Epithelial cells": "Epithelial",
        "Myeloids": "Myeloid",
        "Stromal cells": "Stromal",
        "T cells": "T_NK",
        "B cells": "B_Plasma",
        "Mast cells": "Mast",
    }).fillna("Unclassified")
    cells["qc_pass"] = True
    manifest = {
        "cohort": "GSE132465",
        "annotation_source": repo_path(GSE132465_ANNOTATION),
        "annotation_sha256": sha256(GSE132465_ANNOTATION),
        "source_cell_level_output": repo_path(GSE132465_CELL_LEVEL),
        "cell_type_source": "GEO curated annotation from prior single-cell analysis",
        "side_mapping": "rectum/sigmoid/rectosigmoid=Left; ascending/hepatic flexure=Right; normal samples inherit matched patient side",
    }
    return cells, manifest


def load_gse200997() -> tuple[pd.DataFrame, dict]:
    annotation = pd.read_csv(GSE200997_ANNOTATION, dtype=str).fillna("").rename(columns={"Unnamed: 0": "cell_id"})
    cell_ids, targets, markers, total_umi, n_genes = extract_dense_matrix(GSE200997_COUNTS, delimiter=b",")
    if cell_ids != annotation["cell_id"].tolist():
        raise AssertionError("GSE200997 annotation cell order does not match count matrix")
    cells = annotation[["cell_id", "samples", "Condition", "Location", "MSI_Status", "bulk_prediction", "prediction"]].copy()
    cells = cells.rename(columns={"samples": "Sample", "Condition": "Class", "Location": "Side"})
    cells["cohort"] = "GSE200997"
    cells["Patient"] = cells["Sample"].str.replace(r"^[TB]_", "", regex=True)
    cells["cell_subtype"] = cells["prediction"].replace("", np.nan).fillna("")
    cells = add_target_metrics(cells, targets, total_umi, n_genes)
    cells = infer_marker_cell_types(cells, markers)
    cells["cell_type_source"] = "Marker-inferred from raw UMI"
    cells["cross_compartment"] = cells["cell_type"].map({
        "Epithelial": "Epithelial",
        "Myeloid": "Myeloid",
        "Neutrophil": "Myeloid",
        "Fibroblast": "Stromal",
        "Endothelial": "Stromal",
        "T_NK": "T_NK",
        "B_Plasma": "B_Plasma",
        "Mast": "Mast",
    }).fillna("Unclassified")
    manifest = {
        "cohort": "GSE200997",
        "annotation_source": repo_path(GSE200997_ANNOTATION),
        "annotation_sha256": sha256(GSE200997_ANNOTATION),
        "count_source": repo_path(GSE200997_COUNTS),
        "count_sha256": sha256(GSE200997_COUNTS),
        "cell_type_source": "Marker-inferred; GEO file contains sample/condition/location/CMS metadata but no cell-type labels",
        "marker_panels": MARKER_PANELS,
        "qc_threshold": {"min_umi": MIN_QC_UMI, "min_genes": MIN_QC_GENES},
    }
    return cells, manifest


def read_gse188711_sample(sample_id: str, sample_name: str, side: str) -> pd.DataFrame:
    matrix_path = next(GSE188711_DIR.glob(f"{sample_id}_matrix_*.mtx.gz"))
    feature_path = next(GSE188711_DIR.glob(f"{sample_id}_features_*.tsv.gz"))
    barcode_path = next(GSE188711_DIR.glob(f"{sample_id}_barcodes_*.tsv.gz"))
    feature_table = pd.read_csv(feature_path, sep="\t", header=None, compression="gzip", dtype=str)
    gene_names = feature_table.iloc[:, 1].str.strip().str.upper().tolist()
    barcodes = pd.read_csv(barcode_path, sep="\t", header=None, compression="gzip", dtype=str).iloc[:, 0].tolist()
    matrix = mmread(matrix_path).tocsr()
    if matrix.shape != (len(gene_names), len(barcodes)):
        raise AssertionError(f"Unexpected matrix shape for {sample_id}: {matrix.shape}")
    gene_to_row: dict[str, int] = {}
    for index, gene in enumerate(gene_names):
        gene_to_row.setdefault(gene, index)
    total_umi = np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64)
    n_genes = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(np.int32)
    targets = {gene: np.asarray(matrix[gene_to_row[gene], :].toarray()).ravel().astype(np.int32) for gene in TARGETS if gene in gene_to_row}
    markers = {gene: np.asarray(matrix[gene_to_row[gene], :].toarray()).ravel().astype(np.int32) for gene in ALL_MARKERS if gene in gene_to_row}
    if set(targets) != set(TARGETS):
        raise AssertionError(f"GSE188711 missing targets in {sample_id}: {set(TARGETS) - set(targets)}")
    cells = pd.DataFrame({
        "cell_id": [f"{sample_name}_{barcode}" for barcode in barcodes],
        "Sample": sample_name,
        "Patient": sample_name,
        "Class": "Tumor",
        "Side": side,
        "cohort": "GSE188711",
        "cell_subtype": "",
    })
    cells = add_target_metrics(cells, targets, total_umi, n_genes)
    cells = infer_marker_cell_types(cells, markers)
    cells["cell_type_source"] = "Marker-inferred from 10X UMI"
    cells["cross_compartment"] = cells["cell_type"].map({
        "Epithelial": "Epithelial",
        "Myeloid": "Myeloid",
        "Neutrophil": "Myeloid",
        "Fibroblast": "Stromal",
        "Endothelial": "Stromal",
        "T_NK": "T_NK",
        "B_Plasma": "B_Plasma",
        "Mast": "Mast",
    }).fillna("Unclassified")
    return cells


def load_gse188711() -> tuple[pd.DataFrame, dict]:
    sample_map = {
        "GSM5688706": ("WGC", "Left"),
        "GSM5688707": ("JCA", "Left"),
        "GSM5688708": ("LS-CRC3", "Left"),
        "GSM5688709": ("RS-CRC1", "Right"),
        "GSM5688710": ("R_CRC3", "Right"),
        "GSM5688711": ("R_CRC4", "Right"),
    }
    frames = [read_gse188711_sample(sample_id, sample_name, side) for sample_id, (sample_name, side) in sample_map.items()]
    cells = pd.concat(frames, ignore_index=True)
    manifest = {
        "cohort": "GSE188711",
        "source_tar": repo_path(ROOT / "work" / "scRNA_GSE188711" / "GSE188711_RAW.tar"),
        "source_tar_sha256": sha256(ROOT / "work" / "scRNA_GSE188711" / "GSE188711_RAW.tar"),
        "cell_type_source": "Marker-inferred from raw 10X UMI; GEO series has no cell-type annotation file",
        "sample_side_map": {sample_id: {"sample": sample, "side": side} for sample_id, (sample, side) in sample_map.items()},
        "marker_panels": MARKER_PANELS,
        "qc_threshold": {"min_umi": MIN_QC_UMI, "min_genes": MIN_QC_GENES},
    }
    return cells, manifest


def target_long(cells: pd.DataFrame) -> pd.DataFrame:
    id_columns = ["cohort", "cell_id", "Patient", "Sample", "Class", "Side", "cell_type", "cell_subtype", "myeloid_subtype", "epithelial_context", "cell_type_source", "cross_compartment", "total_umi", "n_genes", "qc_pass"]
    rows = []
    for gene in TARGETS:
        frame = cells[id_columns].copy()
        frame["gene_symbol"] = gene
        frame["umi"] = cells[f"{gene}_umi"]
        frame["detected"] = cells[f"{gene}_detected"]
        frame["cpm"] = cells[f"{gene}_cpm"]
        frame["log1p_cpm"] = cells[f"{gene}_log1p_cpm"]
        rows.append(frame)
    return pd.concat(rows, ignore_index=True)


def summarize_by_group(cells: pd.DataFrame, group_cols: list[str], min_cells: int = MIN_CELLS_PER_SAMPLE_CELLTYPE) -> pd.DataFrame:
    usable = cells[cells["qc_pass"]].copy()
    rows = []
    for gene in TARGETS:
        grouped = usable.groupby(group_cols, dropna=False, sort=True)
        summary = grouped.agg(
            n_cells=("cell_id", "size"),
            total_library_umi=("total_umi", "sum"),
            total_umi=(f"{gene}_umi", "sum"),
            detected_cells=(f"{gene}_detected", "sum"),
            mean_cpm_per_cell=(f"{gene}_cpm", "mean"),
            mean_log1p_cpm=(f"{gene}_log1p_cpm", "mean"),
        ).reset_index()
        summary["gene_symbol"] = gene
        summary["pseudobulk_cpm"] = summary["total_umi"] * 1_000_000.0 / summary["total_library_umi"].replace(0, np.nan)
        summary["log1p_pseudobulk_cpm"] = np.log1p(summary["pseudobulk_cpm"])
        summary["detection_rate"] = summary["detected_cells"] / summary["n_cells"]
        summary["passes_min_cells"] = summary["n_cells"].ge(min_cells)
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def paired_validation(pseudobulk: pd.DataFrame, cohort: str, group_column: str = "cell_type") -> pd.DataFrame:
    subset = pseudobulk[(pseudobulk["cohort"] == cohort) & pseudobulk["passes_min_cells"]].copy()
    rows = []
    for gene in TARGETS:
        for group in sorted(subset[group_column].dropna().unique()):
            current = subset[(subset["gene_symbol"] == gene) & (subset[group_column] == group)]
            tumor = current[current["Class"] == "Tumor"].set_index("Patient")
            normal = current[current["Class"] == "Normal"].set_index("Patient")
            patients = sorted(set(tumor.index) & set(normal.index))
            tumor_values = tumor.loc[patients, "log1p_pseudobulk_cpm"].to_numpy(float) if patients else np.array([])
            normal_values = normal.loc[patients, "log1p_pseudobulk_cpm"].to_numpy(float) if patients else np.array([])
            statistic, p_value = safe_wilcoxon(normal_values, tumor_values)
            delta = tumor_values - normal_values
            rows.append({
                "cohort": cohort,
                "gene_symbol": gene,
                "cell_group": group,
                "matched_patient_n": len(patients),
                "matched_patients": ";".join(patients),
                "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm": float(np.median(delta)) if len(delta) else np.nan,
                "mean_delta_tumor_minus_normal_log1p_pseudobulk_cpm": float(np.mean(delta)) if len(delta) else np.nan,
                "tumor_gt_normal_n": int(np.sum(delta > 0)) if len(delta) else 0,
                "tumor_eq_normal_n": int(np.sum(delta == 0)) if len(delta) else 0,
                "tumor_lt_normal_n": int(np.sum(delta < 0)) if len(delta) else 0,
                "wilcoxon_statistic": statistic,
                "wilcoxon_p_value": p_value,
            })
    output = pd.DataFrame(rows)
    output["wilcoxon_fdr_bh"] = bh_adjust(output["wilcoxon_p_value"])
    output["direction"] = np.where(output["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"].gt(0), "Tumor higher", np.where(output["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"].lt(0), "Normal higher", "No difference"))
    return output


def sidedness_validation(pseudobulk: pd.DataFrame, cohort: str, group_column: str = "cell_type") -> pd.DataFrame:
    subset = pseudobulk[(pseudobulk["cohort"] == cohort) & (pseudobulk["Class"] == "Tumor") & pseudobulk["passes_min_cells"] & pseudobulk["Side"].isin(["Left", "Right"])].copy()
    rows = []
    for gene in TARGETS:
        for group in sorted(subset[group_column].dropna().unique()):
            current = subset[(subset["gene_symbol"] == gene) & (subset[group_column] == group)]
            left = current[current["Side"] == "Left"].set_index("Sample")["log1p_pseudobulk_cpm"].to_numpy(float)
            right = current[current["Side"] == "Right"].set_index("Sample")["log1p_pseudobulk_cpm"].to_numpy(float)
            statistic, p_value = safe_mannwhitney(left, right)
            rows.append({
                "cohort": cohort,
                "gene_symbol": gene,
                "cell_group": group,
                "left_sample_n": len(left),
                "right_sample_n": len(right),
                "median_right_minus_left_log1p_pseudobulk_cpm": float(np.median(right) - np.median(left)) if len(left) and len(right) else np.nan,
                "mannwhitney_statistic": statistic,
                "mannwhitney_p_value": p_value,
            })
    output = pd.DataFrame(rows)
    output["mannwhitney_fdr_bh"] = bh_adjust(output["mannwhitney_p_value"])
    output["direction"] = np.where(output["median_right_minus_left_log1p_pseudobulk_cpm"].gt(0), "Right higher", np.where(output["median_right_minus_left_log1p_pseudobulk_cpm"].lt(0), "Left higher", "No difference"))
    return output


def cross_compartment_paired_validation(cells: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate broad compartments with a documented low-cell threshold.

    GSE200997 has sparse normal myeloid cells after marker inference. The
    primary cell-type analysis keeps a 50-cell/sample threshold; this separate
    broad-compartment analysis uses 10 cells/sample and is labeled exploratory.
    """
    cross_pb = summarize_by_group(
        cells,
        ["cohort", "Patient", "Sample", "Class", "Side", "cross_compartment"],
        min_cells=MIN_CROSS_CELLS_PER_SAMPLE_COMPARTMENT,
    )
    rows = []
    for cohort in ["GSE200997", "GSE132465"]:
        validated = paired_validation(cross_pb, cohort, group_column="cross_compartment")
        validated["cross_compartment"] = validated["cell_group"]
        rows.append(validated)
    return pd.concat(rows, ignore_index=True), cross_pb


def cross_cohort_replication(effects: pd.DataFrame) -> pd.DataFrame:
    output_rows = []
    for (gene, compartment), group in effects.groupby(["gene_symbol", "cross_compartment"], sort=True):
        by_cohort = group.set_index("cohort")
        first = by_cohort.loc["GSE200997"] if "GSE200997" in by_cohort.index else None
        second = by_cohort.loc["GSE132465"] if "GSE132465" in by_cohort.index else None
        d1 = first["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"] if first is not None else np.nan
        d2 = second["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"] if second is not None else np.nan
        direction_agree = bool(np.isfinite(d1) and np.isfinite(d2) and ((d1 > 0 and d2 > 0) or (d1 < 0 and d2 < 0) or (d1 == 0 and d2 == 0)))
        output_rows.append({
            "gene_symbol": gene,
            "cross_compartment": compartment,
            "gse200997_matched_n": first["matched_patient_n"] if first is not None else 0,
            "gse200997_median_delta": d1,
            "gse200997_p": first["wilcoxon_p_value"] if first is not None else np.nan,
            "gse132465_matched_n": second["matched_patient_n"] if second is not None else 0,
            "gse132465_median_delta": d2,
            "gse132465_p": second["wilcoxon_p_value"] if second is not None else np.nan,
            "same_direction": direction_agree,
            "replication_status": "Direction concordant" if direction_agree else "Not concordant/insufficient",
        })
    return pd.DataFrame(output_rows)


def plot_localization(localization: pd.DataFrame) -> None:
    subset = localization.copy()
    subset["panel"] = subset.apply(
        lambda row: f"{row['cohort']} | {row['Class']}" + (f" | {row['Side']}" if row["Side"] else ""),
        axis=1,
    )
    panels = sorted(subset["panel"].unique())
    cell_types = sorted(subset["cell_type"].unique())
    fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex=True, constrained_layout=True)
    for ax, gene in zip(axes, TARGETS):
        current = subset[subset["gene_symbol"].eq(gene)].copy()
        current["x"] = current["cell_type"].map({value: index for index, value in enumerate(cell_types)})
        current["y"] = current["panel"].map({value: index for index, value in enumerate(panels)})
        scatter = ax.scatter(current["x"], current["y"], s=20 + 350 * current["detection_rate"], c=current["mean_log1p_cpm"], cmap="viridis", edgecolors="black", linewidths=0.25)
        ax.set_title(gene)
        ax.set_yticks(range(len(panels)))
        ax.set_yticklabels(panels, fontsize=8)
        ax.grid(axis="x", alpha=0.15)
        fig.colorbar(scatter, ax=ax, label="Mean log1p CPM")
    axes[-1].set_xticks(range(len(cell_types)))
    axes[-1].set_xticklabels(cell_types, rotation=35, ha="right")
    axes[-1].set_xlabel("Cell type (GSE200997/GSE188711 marker-inferred; GSE132465 curated)")
    fig.suptitle("CEBPB/CD36 localization across CRC single-cell cohorts\nDot size = detection rate")
    fig.savefig(OUTPUT_DIR / "multicohort_target_localization.png", dpi=220)
    plt.close(fig)


def plot_sidedness(sidedness: pd.DataFrame) -> None:
    subset = sidedness.copy()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True, constrained_layout=True)
    for ax, gene in zip(axes, TARGETS):
        current = subset[subset["gene_symbol"].eq(gene)].copy()
        current["x"] = np.arange(len(current))
        colors = np.where(current["direction"].eq("Right higher"), "#dc2626", np.where(current["direction"].eq("Left higher"), "#2563eb", "#777777"))
        ax.scatter(current["x"], current["median_right_minus_left_log1p_pseudobulk_cpm"], c=colors, s=45)
        for i, (_, row) in enumerate(current.iterrows()):
            ax.text(i, row["median_right_minus_left_log1p_pseudobulk_cpm"], f"{row['cohort']}\n{row['cell_group']}", fontsize=7, rotation=50, ha="left", va="bottom")
        ax.axhline(0, color="black", linewidth=0.7)
        ax.set_title(gene)
        ax.set_ylabel("Median Right − Left\nlog1p pseudobulk CPM")
        ax.set_xticks([])
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Left-vs-right CRC sidedness sensitivity")
    fig.savefig(OUTPUT_DIR / "multicohort_left_right_sidedness.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    gse132465, manifest_132465 = load_gse132465()
    gse200997, manifest_200997 = load_gse200997()
    gse188711, manifest_188711 = load_gse188711()
    cells = pd.concat([gse132465, gse200997, gse188711], ignore_index=True, sort=False).fillna("")
    cells["Cell_type"] = cells["cell_type"]

    localization = summarize_by_group(cells, ["cohort", "Class", "Side", "cell_type", "cell_type_source"])
    pseudobulk = summarize_by_group(cells, ["cohort", "Patient", "Sample", "Class", "Side", "cell_type", "cell_type_source"])
    paired_200997 = paired_validation(pseudobulk, "GSE200997")
    paired_132465 = paired_validation(pseudobulk, "GSE132465")
    paired = pd.concat([paired_200997, paired_132465], ignore_index=True)
    side_200997 = sidedness_validation(pseudobulk, "GSE200997")
    side_132465 = sidedness_validation(pseudobulk, "GSE132465")
    side_188711 = sidedness_validation(pseudobulk, "GSE188711")
    sidedness = pd.concat([side_200997, side_132465, side_188711], ignore_index=True)
    cross_effects, cross_pseudobulk = cross_compartment_paired_validation(cells)
    replication = cross_cohort_replication(cross_effects)

    long = target_long(cells)
    long.to_csv(OUTPUT_DIR / "DINP_CRC_multicohort_CEBPB_CD36_cell_level_long.csv", index=False)
    localization.to_csv(OUTPUT_DIR / "DINP_CRC_multicohort_CEBPB_CD36_localization.csv", index=False)
    pseudobulk.to_csv(OUTPUT_DIR / "DINP_CRC_multicohort_CEBPB_CD36_patient_sample_pseudobulk.csv", index=False)
    paired.to_csv(OUTPUT_DIR / "DINP_CRC_GSE200997_GSE132465_paired_tumor_normal_validation.csv", index=False)
    cross_effects.to_csv(OUTPUT_DIR / "DINP_CRC_GSE200997_GSE132465_cross_compartment_paired_validation.csv", index=False)
    cross_pseudobulk.to_csv(OUTPUT_DIR / "DINP_CRC_GSE200997_GSE132465_cross_compartment_patient_sample_pseudobulk.csv", index=False)
    sidedness.to_csv(OUTPUT_DIR / "DINP_CRC_GSE200997_GSE132465_GSE188711_left_right_sensitivity.csv", index=False)
    replication.to_csv(OUTPUT_DIR / "DINP_CRC_GSE200997_GSE132465_cross_cohort_replication.csv", index=False)
    plot_localization(localization)
    plot_sidedness(sidedness)

    tumor_200997 = int(((cells["cohort"] == "GSE200997") & (cells["Class"] == "Tumor") & cells["qc_pass"]).sum())
    normal_200997 = int(((cells["cohort"] == "GSE200997") & (cells["Class"] == "Normal") & cells["qc_pass"]).sum())
    tumor_132465 = int(((cells["cohort"] == "GSE132465") & (cells["Class"] == "Tumor") & cells["qc_pass"]).sum())
    normal_132465 = int(((cells["cohort"] == "GSE132465") & (cells["Class"] == "Normal") & cells["qc_pass"]).sum())
    tumor_188711 = int(((cells["cohort"] == "GSE188711") & cells["qc_pass"]).sum())
    raw_188711 = int((cells["cohort"] == "GSE188711").sum())
    failed_188711 = raw_188711 - tumor_188711
    manifest = {
        "targets": TARGETS,
        "qc_threshold_for_raw_matrix_cohorts": {"min_umi": MIN_QC_UMI, "min_genes": MIN_QC_GENES},
        "pseudobulk": {
            "grouping": "cohort x patient x sample x class x side x cell_type",
            "normalization": "aggregate target UMI / aggregate library UMI x 1e6 = pseudobulk CPM",
            "minimum_cells_per_sample_celltype": MIN_CELLS_PER_SAMPLE_CELLTYPE,
            "minimum_cells_per_sample_cross_compartment_exploratory": MIN_CROSS_CELLS_PER_SAMPLE_COMPARTMENT,
            "tumor_normal_test": "paired Wilcoxon across matched patients, BH across tests within cohort",
            "left_right_test": "sample-level Mann-Whitney U among tumor samples, BH within cohort",
        },
        "cohorts": {
            "GSE132465": manifest_132465,
            "GSE200997": manifest_200997,
            "GSE188711": manifest_188711,
        },
        "cell_counts_after_qc": {
            "GSE200997_tumor": tumor_200997,
            "GSE200997_normal": normal_200997,
            "GSE132465_tumor": tumor_132465,
            "GSE132465_normal": normal_132465,
            "GSE188711_tumor": tumor_188711,
            "combined": int(cells["qc_pass"].sum()),
        },
        "cell_type_source_warning": "GSE200997 and GSE188711 labels are marker-inferred; GSE132465 labels are GEO-curated. Tumor epithelial is a sample-context proxy, not CNV-confirmed malignant status.",
        "qc": {
            "combined_cells": int(len(cells)),
            "combined_qc_pass_cells": int(cells["qc_pass"].sum()),
            "long_rows": int(len(long)),
            "localization_rows": int(len(localization)),
            "pseudobulk_rows": int(len(pseudobulk)),
            "paired_rows": int(len(paired)),
            "cross_compartment_paired_rows": int(len(cross_effects)),
            "cross_compartment_pseudobulk_rows": int(len(cross_pseudobulk)),
            "sidedness_rows": int(len(sidedness)),
            "replication_rows": int(len(replication)),
            "target_total_umi": {gene: int(cells[f"{gene}_umi"].sum()) for gene in TARGETS},
        },
    }
    (OUTPUT_DIR / "DINP_CRC_multicohort_CEBPB_CD36_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    report = [
        "# DINP–CRC CEBPB/CD36 multicohort single-cell analysis",
        "",
        "## Frozen cohort design",
        "",
        "- **GSE200997**: discovery cohort for cell-type localization, tumor–normal pseudobulk, and left–right sensitivity.",
        "- **GSE132465**: independent replication cohort for cell-type localization and matched tumor–normal pseudobulk; tumor region is also used for a secondary left–right sensitivity check.",
        "- **GSE188711**: sidedness-only supplementary cohort; no normal samples, so it is not used for tumor–normal validation.",
        "",
        "## Important annotation audit",
        "",
        "GSE132465 uses its GEO-curated cell-type labels. The downloaded GSE200997 annotation file contains sample, condition, location, MSI, and CMS metadata but no cell-type labels; GSE188711 provides raw 10X MTX/TSV files without a cell-type annotation table. Their major cell compartments were therefore inferred with fixed canonical marker-score panels and are explicitly labeled `Marker-inferred` in the outputs. Tumor epithelial cells are treated as a tumor-sample epithelial proxy, not as CNV-confirmed malignant cells.",
        "",
        "## Statistical design",
        "",
        f"Raw UMI matrices were filtered at ≥{MIN_QC_UMI} total UMI and ≥{MIN_QC_GENES} detected genes for GSE200997/GSE188711; GSE132465 is the already-filtered GEO processed cell table used in the previous analysis. Pseudobulk was aggregated at patient × sample × class × cell type and normalized as aggregate target UMI / aggregate library UMI × 1e6. Tumor–normal comparisons use paired Wilcoxon tests across matched patients; left–right comparisons use sample-level Mann–Whitney U tests. A minimum of {MIN_CELLS_PER_SAMPLE_CELLTYPE} cells per sample × cell type is required, and BH correction is applied within each analysis family.",
        "",
        "## Cohort sizes after QC",
        "",
        f"- GSE200997: {tumor_200997:,} tumor cells and {normal_200997:,} normal cells retained.",
        f"- GSE132465: {tumor_132465:,} tumor cells and {normal_132465:,} normal cells retained.",
        f"- GSE188711: {raw_188711:,} raw barcodes; {tumor_188711:,} cells retained after QC ({failed_188711:,} failed) for sidedness analysis. The GEO raw supplementary files contain more barcodes than the source paper's reported 27,927 high-quality cells; this re-analysis does not claim to reproduce the paper's original filtering exactly.",
        "",
        "## Key outputs",
        "",
        "- `DINP_CRC_multicohort_CEBPB_CD36_localization.csv`: localization by cohort, class, and cell type.",
        "- `DINP_CRC_multicohort_CEBPB_CD36_patient_sample_pseudobulk.csv`: patient/sample-level pseudobulk table.",
        "- `DINP_CRC_GSE200997_GSE132465_paired_tumor_normal_validation.csv`: primary tumor–normal validation.",
        f"- `DINP_CRC_GSE200997_GSE132465_cross_compartment_paired_validation.csv`: exploratory broad-compartment tumor–normal validation using ≥{MIN_CROSS_CELLS_PER_SAMPLE_COMPARTMENT} cells per sample × compartment to retain sparse GSE200997 myeloid context.",
        "- `DINP_CRC_GSE200997_GSE132465_cross_compartment_patient_sample_pseudobulk.csv`: broad-compartment pseudobulk source table for the exploratory validation.",
        "- `DINP_CRC_GSE200997_GSE132465_GSE188711_left_right_sensitivity.csv`: secondary sidedness analysis.",
        "- `DINP_CRC_GSE200997_GSE132465_cross_cohort_replication.csv`: direction-concordance audit between the two tumor–normal cohorts.",
        "",
        "## Provisional interpretation",
        "",
        "Use the cross-cohort replication table, not cell-level p-values, to decide whether CEBPB/CD36 have reproducible tumor–normal effects within the same broad compartment. The primary cell-type validation requires ≥50 cells per sample × cell type; the separate broad-compartment myeloid analysis uses ≥10 cells and is exploratory because marker-inferred normal myeloid cells are sparse in GSE200997. Myeloid localization and epithelial tumor-sample localization remain descriptive unless the direction is reproduced at patient level in both cohorts.",
    ]
    (OUTPUT_DIR / "DINP_CRC_multicohort_CEBPB_CD36_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"combined_cells={len(cells)} qc_pass={int(cells['qc_pass'].sum())}")
    print("cohort counts:")
    print(cells.groupby(["cohort", "Class", "cell_type_source"]).size().to_string())
    print("cross-cohort replication:")
    print(replication.to_string(index=False))
    print("sidedness head:")
    print(sidedness.head(30).to_string(index=False))


if __name__ == "__main__":
    main()
