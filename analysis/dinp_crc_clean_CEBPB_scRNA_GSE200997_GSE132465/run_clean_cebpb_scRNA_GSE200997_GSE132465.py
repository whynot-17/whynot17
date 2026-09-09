from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs" / "DINP_CRC_clean_CEBPB_scRNA_GSE200997_GSE132465"

TARGET = "CEBPB"
MIN_QC_GENES = 200
MIN_QC_UMI = 500
MIN_CELLS_PER_SAMPLE_CELLTYPE = 50

GSE132465_DIR = ROOT / "work" / "scRNA_GSE132465"
GSE132465_ANNOTATION = GSE132465_DIR / "cell_annotation.txt.gz"
GSE132465_COUNTS = GSE132465_DIR / "raw_UMI_count_matrix.txt.gz"

GSE200997_DIR = ROOT / "work" / "scRNA_GSE200997"
GSE200997_ANNOTATION = GSE200997_DIR / "cell_annotation.csv.gz"
GSE200997_COUNTS = GSE200997_DIR / "raw_UMI_count_matrix.csv.gz"

REFERENCE_LABELS = [
    "Epithelial cells",
    "Myeloids",
    "Stromal cells",
    "T cells",
    "B cells",
    "Mast cells",
]

# Canonical broad-compartment markers. They are used for an independent audit of
# the GSE200997 label transfer; no marker-max assignment is used as the final rule.
MARKER_PANELS = {
    "Epithelial cells": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CEACAM5", "MUC1"],
    "T cells": ["CD3D", "CD3E", "TRAC", "CD247", "CD8A", "NKG7", "GNLY", "TRBC1"],
    "B cells": ["CD79A", "MS4A1", "CD37", "CD74", "CD79B", "CD19", "JCHAIN", "MZB1", "SDC1"],
    "Myeloids": ["LYZ", "LST1", "TYROBP", "FCER1G", "CTSS", "CD68", "C1QA", "C1QB", "C1QC"],
    "Stromal cells": [
        "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "THY1", "PDGFRA",
        "PECAM1", "VWF", "EMCN", "KDR", "ENG", "CLDN5", "RAMP2",
    ],
    "Mast cells": ["TPSAB1", "TPSB2", "KIT", "CPA3", "GATA2", "MS4A2"],
}
# Keep the target out of the annotation feature set to prevent circularity:
# cell-type mapping must not use the same CEBPB expression that is validated.
ALL_MARKERS = sorted({gene for genes in MARKER_PANELS.values() for gene in genes})
COUNT_GENES = set(ALL_MARKERS) | {TARGET}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_path(path: Path) -> str:
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


def safe_wilcoxon(normal: np.ndarray, tumor: np.ndarray) -> tuple[float, float]:
    delta = tumor - normal
    if len(delta) < 3:
        return np.nan, np.nan
    if np.allclose(delta, 0):
        return 0.0, 1.0
    result = wilcoxon(tumor, normal, zero_method="wilcox", alternative="two-sided", method="auto")
    return float(result.statistic), float(result.pvalue)


def extract_dense_matrix(
    path: Path,
    delimiter: bytes,
    wanted_genes: set[str],
) -> tuple[list[str], dict[str, np.ndarray], np.ndarray, np.ndarray, list[str]]:
    """Stream a GEO dense gene x cell matrix, retaining only target/marker rows."""
    gene_counts: dict[str, np.ndarray] = {}
    with gzip.open(path, "rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(delimiter)
        cell_ids = [value.decode("utf-8").strip().strip('"') for value in header[1:]]
        total_umi = np.zeros(len(cell_ids), dtype=np.int64)
        n_genes = np.zeros(len(cell_ids), dtype=np.int32)
        all_genes: list[str] = []
        for raw_line in handle:
            line = raw_line.rstrip(b"\r\n")
            if not line:
                continue
            gene_bytes, payload = line.split(delimiter, 1)
            gene = gene_bytes.decode("utf-8").strip().strip('"').upper()
            values = np.fromstring(payload, dtype=np.int64, sep=delimiter.decode("utf-8"))
            if len(values) != len(cell_ids):
                raise AssertionError(f"Unexpected values length for {gene}: {len(values)} != {len(cell_ids)}")
            total_umi += values
            n_genes += values > 0
            if gene in wanted_genes:
                # Sum duplicate feature rows if a GEO matrix contains them.
                gene_counts[gene] = gene_counts.get(gene, np.zeros(len(cell_ids), dtype=np.int64)) + values
            all_genes.append(gene)
    missing = wanted_genes - set(gene_counts)
    if missing:
        raise AssertionError(f"Missing requested genes in {path.name}: {sorted(missing)}")
    return cell_ids, gene_counts, total_umi, n_genes, all_genes


def expression_from_counts(gene_counts: dict[str, np.ndarray], total_umi: np.ndarray, genes: list[str]) -> np.ndarray:
    denominator = total_umi.astype(float)
    rows = []
    for gene in genes:
        cpm = np.divide(
            gene_counts[gene].astype(float) * 1_000_000.0,
            denominator,
            out=np.zeros(len(total_umi), dtype=float),
            where=denominator > 0,
        )
        rows.append(np.log1p(cpm))
    return np.vstack(rows).T


def add_target_metrics(cells: pd.DataFrame, gene_counts: dict[str, np.ndarray], total_umi: np.ndarray, n_genes: np.ndarray) -> pd.DataFrame:
    output = cells.copy()
    output["total_umi"] = total_umi.astype(np.int64)
    output["n_genes"] = n_genes.astype(np.int32)
    output["qc_pass"] = output["total_umi"].ge(MIN_QC_UMI) & output["n_genes"].ge(MIN_QC_GENES)
    denominator = output["total_umi"].to_numpy(dtype=float)
    values = gene_counts[TARGET].astype(np.int64)
    output[f"{TARGET}_umi"] = values
    output[f"{TARGET}_detected"] = values > 0
    output[f"{TARGET}_cpm"] = np.divide(
        values * 1_000_000.0,
        denominator,
        out=np.zeros(len(output), dtype=float),
        where=denominator > 0,
    )
    output[f"{TARGET}_log1p_cpm"] = np.log1p(output[f"{TARGET}_cpm"])
    return output


def load_gse132465() -> tuple[pd.DataFrame, dict, dict[str, np.ndarray], np.ndarray]:
    annotation = pd.read_csv(GSE132465_ANNOTATION, sep="\t", dtype=str).fillna("")
    for column in annotation.columns:
        annotation[column] = annotation[column].astype(str).str.strip()
    annotation = annotation.rename(columns={"Index": "cell_id", "Cell_type": "cell_type", "Cell_subtype": "cell_subtype"})
    cell_ids, counts, total_umi, n_genes, all_genes = extract_dense_matrix(GSE132465_COUNTS, b"\t", COUNT_GENES)
    if cell_ids != annotation["cell_id"].tolist():
        raise AssertionError("GSE132465 annotation cell order does not match the raw count matrix")
    cells = annotation[["cell_id", "Patient", "Class", "Sample", "cell_type", "cell_subtype"]].copy()
    cells["cohort"] = "GSE132465"
    cells["cell_type_source"] = "GEO curated annotation"
    cells["annotation_status"] = "Curated"
    cells["analysis_cell_type"] = cells["cell_type"]
    cells = add_target_metrics(cells, counts, total_umi, n_genes)
    manifest = {
        "cohort": "GSE132465",
        "annotation_source": repo_path(GSE132465_ANNOTATION),
        "annotation_sha256": sha256(GSE132465_ANNOTATION),
        "count_source": repo_path(GSE132465_COUNTS),
        "count_sha256": sha256(GSE132465_COUNTS),
        "cell_type_source": "GEO curated annotation read directly from the source annotation file",
        "raw_gene_rows": len(all_genes),
        "duplicate_gene_rows": len(all_genes) - len(set(all_genes)),
        "qc_threshold": {"min_umi": MIN_QC_UMI, "min_genes": MIN_QC_GENES},
    }
    return cells, manifest, counts, total_umi


def cosine_similarity_matrix(query: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    query_norm = np.linalg.norm(query, axis=1, keepdims=True)
    centroid_norm = np.linalg.norm(centroids, axis=1, keepdims=True).T
    numerator = query @ centroids.T
    denominator = query_norm @ centroid_norm
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def leave_one_out_median_centroids(group_data: np.ndarray) -> np.ndarray:
    """Return an exact leave-one-out median centroid for every row in a group.

    The production classifier uses median centroids, so calibration must remove
    each reference cell and recompute the median rather than switching to a mean.
    This implementation uses sorted values and vectorized rank arithmetic so the
    result is identical to np.median(group_data without row_i, axis=0).
    """
    n_rows, n_features = group_data.shape
    if n_rows < 2:
        median = np.median(group_data, axis=0)
        return np.repeat(median[None, :], n_rows, axis=0)
    order = np.argsort(group_data, axis=0, kind="mergesort")
    sorted_values = np.take_along_axis(group_data, order, axis=0)
    ranks = np.empty_like(order)
    columns = np.arange(n_features)
    ranks[order, columns] = np.arange(n_rows)[:, None]
    flat_sorted = sorted_values.ravel(order="F")
    offsets = np.arange(n_features) * n_rows
    remaining = n_rows - 1
    if remaining % 2 == 1:
        middle = remaining // 2
        indices = middle + (ranks <= middle)
        return flat_sorted[indices + offsets]
    lower = remaining // 2 - 1
    upper = remaining // 2
    lower_indices = lower + (ranks <= lower)
    upper_indices = upper + (ranks <= upper)
    return 0.5 * (
        flat_sorted[lower_indices + offsets] + flat_sorted[upper_indices + offsets]
    )


def fit_reference_mapping(
    ref_cells: pd.DataFrame,
    ref_counts: dict[str, np.ndarray],
    ref_total_umi: np.ndarray,
) -> tuple[dict, np.ndarray, np.ndarray, list[str]]:
    """Fit a transparent nearest-centroid marker reference using GSE132465 labels."""
    features = [gene for gene in ALL_MARKERS if gene in ref_counts]
    ref_expr = expression_from_counts(ref_counts, ref_total_umi, features)
    ref_labels = ref_cells["cell_type"].to_numpy(dtype=str)
    if sorted(set(ref_labels)) != sorted(REFERENCE_LABELS):
        raise AssertionError(f"Unexpected GSE132465 cell types: {sorted(set(ref_labels))}")
    ref_mean = ref_expr.mean(axis=0)
    ref_std = ref_expr.std(axis=0)
    ref_std[ref_std == 0] = 1.0
    ref_z = (ref_expr - ref_mean) / ref_std
    centroids = np.vstack([np.median(ref_z[ref_labels == label], axis=0) for label in REFERENCE_LABELS])

    # Leave-one-out calibration uses the same median-centroid estimator as the
    # production classifier. Mean-centroid calibration is deliberately avoided.
    loo_centroids = np.repeat(centroids[None, :, :], len(ref_cells), axis=0)
    for index, label in enumerate(REFERENCE_LABELS):
        mask = ref_labels == label
        rows = np.flatnonzero(mask)
        loo_centroids[rows, index, :] = leave_one_out_median_centroids(ref_z[mask])
    loo_scores = np.empty((len(ref_cells), len(REFERENCE_LABELS)), dtype=float)
    for start in range(0, len(ref_cells), 10000):
        stop = min(start + 10000, len(ref_cells))
        for row in range(start, stop):
            loo_scores[row, :] = cosine_similarity_matrix(ref_z[row : row + 1], loo_centroids[row])[0]
    loo_order = np.argsort(-loo_scores, axis=1)
    loo_top1 = loo_scores[np.arange(len(ref_cells)), loo_order[:, 0]]
    loo_top2 = loo_scores[np.arange(len(ref_cells)), loo_order[:, 1]]
    loo_margin = loo_top1 - loo_top2
    sim_threshold_5 = float(max(0.05, np.nanquantile(loo_top1, 0.05)))
    margin_threshold_5 = float(max(0.02, np.nanquantile(loo_margin, 0.05)))
    sim_threshold_10 = float(max(0.05, np.nanquantile(loo_top1, 0.10)))
    margin_threshold_10 = float(max(0.02, np.nanquantile(loo_margin, 0.10)))
    label_array = np.array(REFERENCE_LABELS, dtype=object)
    model = {
        "features": features,
        "labels": REFERENCE_LABELS,
        "mean": ref_mean,
        "std": ref_std,
        "centroids": centroids,
        "sim_threshold_5pct_reference": sim_threshold_5,
        "margin_threshold_5pct_reference": margin_threshold_5,
        "sim_threshold_10pct_reference": sim_threshold_10,
        "margin_threshold_10pct_reference": margin_threshold_10,
        "reference_cells": int(len(ref_cells)),
        "reference_loo_top1_median": float(np.median(loo_top1)),
        "reference_loo_margin_median": float(np.median(loo_margin)),
        "reference_loo_label_accuracy": float(np.mean(label_array[np.argmax(loo_scores, axis=1)] == ref_labels)),
    }
    return model, ref_expr, centroids, features


def apply_reference_mapping(
    cells: pd.DataFrame,
    counts: dict[str, np.ndarray],
    total_umi: np.ndarray,
    model: dict,
) -> pd.DataFrame:
    output = cells.copy()
    features = model["features"]
    expr = expression_from_counts(counts, total_umi, features)
    query_z = (expr - model["mean"]) / model["std"]
    scores = cosine_similarity_matrix(query_z, model["centroids"])
    order = np.argsort(-scores, axis=1)
    top1_index = order[:, 0]
    top2_index = order[:, 1]
    top1 = scores[np.arange(len(output)), top1_index]
    top2 = scores[np.arange(len(output)), top2_index]
    margin = top1 - top2

    panel_scores: dict[str, np.ndarray] = {}
    for panel, genes in MARKER_PANELS.items():
        available = [gene for gene in genes if gene in counts]
        panel_expr = expression_from_counts(counts, total_umi, available)
        panel_scores[panel] = panel_expr.mean(axis=1)
        output[f"marker_score_{panel}"] = panel_scores[panel]
    panel_score_frame = pd.DataFrame(panel_scores)
    marker_panel_count = (panel_score_frame > 0.25).sum(axis=1)
    top1_label = np.array(model["labels"], dtype=object)[top1_index]
    top2_label = np.array(model["labels"], dtype=object)[top2_index]
    top1_pass = top1 >= model["sim_threshold_5pct_reference"]
    margin_pass = margin >= model["margin_threshold_5pct_reference"]
    is_ambiguous = top1_pass & ~margin_pass & (top2 >= model["sim_threshold_5pct_reference"])
    is_low_confidence = ~top1_pass | (~margin_pass & ~is_ambiguous)
    status = np.where(is_low_confidence, "Low_confidence", np.where(is_ambiguous, "Doublet_like_or_ambiguous", "Confident"))

    output["top1_reference_label"] = top1_label
    output["top2_reference_label"] = top2_label
    output["top1_reference_similarity"] = top1
    output["top2_reference_similarity"] = top2
    output["reference_score_margin"] = margin
    output["marker_panel_count_gt_0_25"] = marker_panel_count.astype(np.int16)
    output["annotation_status"] = status
    output["cell_type"] = np.where(status == "Confident", top1_label, "Unclassified")
    output["analysis_cell_type"] = np.where(status == "Confident", top1_label, "")
    output["cell_type_source"] = "Reference-mapped to GSE132465 curated labels"
    return output


def apply_reference_threshold_variant(
    mapped_cells: pd.DataFrame,
    similarity_threshold: float,
    margin_threshold: float,
    require_margin: bool,
) -> pd.DataFrame:
    """Rebuild the analysis inclusion mask without changing stored raw scores."""
    output = mapped_cells.copy()
    top1_pass = output["top1_reference_similarity"].ge(similarity_threshold)
    margin_pass = output["reference_score_margin"].ge(margin_threshold)
    confident = top1_pass & (margin_pass if require_margin else True)
    output["analysis_cell_type"] = np.where(
        confident,
        output["top1_reference_label"],
        "",
    )
    return output


def load_gse200997(model: dict) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    annotation = pd.read_csv(GSE200997_ANNOTATION, dtype=str).fillna("")
    for column in annotation.columns:
        annotation[column] = annotation[column].astype(str).str.strip()
    annotation = annotation.rename(columns={"Unnamed: 0": "cell_id"})
    # Deliberately select no Location/side field: this clean module has no sidedness analysis.
    cell_ids, counts, total_umi, n_genes, all_genes = extract_dense_matrix(GSE200997_COUNTS, b",", COUNT_GENES)
    if cell_ids != annotation["cell_id"].tolist():
        raise AssertionError("GSE200997 annotation cell order does not match the raw count matrix")
    cells = annotation[["cell_id", "samples", "Condition"]].copy()
    cells = cells.rename(columns={"samples": "Sample", "Condition": "Class"})
    cells["cohort"] = "GSE200997"
    cells["Patient"] = cells["Sample"].str.replace(r"^[TB]_", "", regex=True)
    sample_map = cells[["Patient", "Class", "Sample"]].drop_duplicates()
    duplicate_patient_class = sample_map.groupby(["Patient", "Class"], dropna=False)["Sample"].nunique()
    if duplicate_patient_class.gt(1).any():
        bad = duplicate_patient_class[duplicate_patient_class.gt(1)].to_dict()
        raise AssertionError(f"GSE200997 has >1 sample for Patient + Class: {bad}")
    cells["cell_subtype"] = ""
    cells = add_target_metrics(cells, counts, total_umi, n_genes)
    cells = apply_reference_mapping(cells, counts, total_umi, model)
    audit = cells[
        [
            "cell_id", "Patient", "Sample", "Class", "cohort", "qc_pass", "total_umi", "n_genes",
            "top1_reference_label", "top2_reference_label", "top1_reference_similarity",
            "top2_reference_similarity", "reference_score_margin", "marker_panel_count_gt_0_25",
            "annotation_status", "cell_type", "analysis_cell_type",
        ]
    ].copy()
    manifest = {
        "cohort": "GSE200997",
        "annotation_source": repo_path(GSE200997_ANNOTATION),
        "annotation_sha256": sha256(GSE200997_ANNOTATION),
        "count_source": repo_path(GSE200997_COUNTS),
        "count_sha256": sha256(GSE200997_COUNTS),
        "cell_type_source": "Reference-mapped to GSE132465 curated labels; source file used only sample and condition metadata",
        "raw_gene_rows": len(all_genes),
        "duplicate_gene_rows": len(all_genes) - len(set(all_genes)),
        "qc_threshold": {"min_umi": MIN_QC_UMI, "min_genes": MIN_QC_GENES},
        "reference_mapping": {
            "features": model["features"],
            "sim_threshold_5pct": model["sim_threshold_5pct_reference"],
            "margin_threshold_5pct": model["margin_threshold_5pct_reference"],
            "sim_threshold_10pct": model["sim_threshold_10pct_reference"],
            "margin_threshold_10pct": model["margin_threshold_10pct_reference"],
            "reference_loo_top1_median": model["reference_loo_top1_median"],
            "reference_loo_margin_median": model["reference_loo_margin_median"],
            "reference_loo_label_accuracy": model["reference_loo_label_accuracy"],
        },
        "excluded_metadata_fields": ["Location", "Side", "left_right", "tumor_region"],
    }
    return cells, manifest, audit


def summarize_localization(cells: pd.DataFrame) -> pd.DataFrame:
    usable = cells[cells["qc_pass"] & cells["analysis_cell_type"].ne("")].copy()
    result = (
        usable.groupby(["cohort", "Class", "analysis_cell_type", "cell_type_source"], dropna=False, sort=True)
        .agg(
            n_cells=("cell_id", "size"),
            n_patients=("Patient", "nunique"),
            n_samples=("Sample", "nunique"),
            mean_log1p_cpm=(f"{TARGET}_log1p_cpm", "mean"),
            median_log1p_cpm=(f"{TARGET}_log1p_cpm", "median"),
            detection_rate=(f"{TARGET}_detected", "mean"),
            total_umi=(f"{TARGET}_umi", "sum"),
        )
        .reset_index()
        .rename(columns={"analysis_cell_type": "cell_type"})
    )
    return result


def summarize_pseudobulk(cells: pd.DataFrame) -> pd.DataFrame:
    usable = cells[cells["qc_pass"] & cells["analysis_cell_type"].ne("")].copy()
    group_cols = ["cohort", "Patient", "Sample", "Class", "analysis_cell_type", "cell_type_source"]
    result = (
        usable.groupby(group_cols, dropna=False, sort=True)
        .agg(
            n_cells=("cell_id", "size"),
            total_library_umi=("total_umi", "sum"),
            total_umi=(f"{TARGET}_umi", "sum"),
            detected_cells=(f"{TARGET}_detected", "sum"),
            mean_cpm_per_cell=(f"{TARGET}_cpm", "mean"),
            mean_log1p_cpm=(f"{TARGET}_log1p_cpm", "mean"),
        )
        .reset_index()
        .rename(columns={"analysis_cell_type": "cell_type"})
    )
    result["gene_symbol"] = TARGET
    result["pseudobulk_cpm"] = result["total_umi"] * 1_000_000.0 / result["total_library_umi"].replace(0, np.nan)
    result["log1p_pseudobulk_cpm"] = np.log1p(result["pseudobulk_cpm"])
    result["detection_rate"] = result["detected_cells"] / result["n_cells"]
    result["passes_min_cells"] = result["n_cells"].ge(MIN_CELLS_PER_SAMPLE_CELLTYPE)
    return result


def patient_epithelial_cell_count_audit(cells: pd.DataFrame) -> pd.DataFrame:
    """Audit QC and epithelial inclusion counts at patient/sample/condition level."""
    rows = []
    for (cohort, patient, sample, condition), group in cells.groupby(
        ["cohort", "Patient", "Sample", "Class"], dropna=False, sort=True
    ):
        qc = group[group["qc_pass"]]
        rows.append({
            "cohort": cohort,
            "Patient": patient,
            "Sample": sample,
            "Class": condition,
            "qc_pass_total_cells": int(len(qc)),
            "confident_epithelial_cells": int(qc["analysis_cell_type"].eq("Epithelial cells").sum()),
            "low_confidence_cells": int(qc["annotation_status"].eq("Low_confidence").sum()),
            "ambiguous_cells": int(qc["annotation_status"].eq("Doublet_like_or_ambiguous").sum()),
        })
    output = pd.DataFrame(rows)
    output["epithelial_pass_50"] = output["confident_epithelial_cells"].ge(MIN_CELLS_PER_SAMPLE_CELLTYPE)
    return output


def gse200997_epithelial_inclusion_audit(count_audit: pd.DataFrame) -> pd.DataFrame:
    """Explain exactly why a GSE200997 patient enters or fails paired analysis."""
    source = count_audit[count_audit["cohort"].eq("GSE200997")].copy()
    rows = []
    for patient in sorted(source["Patient"].unique()):
        tumor = source[(source["Patient"] == patient) & source["Class"].eq("Tumor")]
        normal = source[(source["Patient"] == patient) & source["Class"].eq("Normal")]
        tumor_count = int(tumor.iloc[0]["confident_epithelial_cells"]) if not tumor.empty else 0
        normal_count = int(normal.iloc[0]["confident_epithelial_cells"]) if not normal.empty else 0
        tumor_pass = bool(not tumor.empty and tumor_count >= MIN_CELLS_PER_SAMPLE_CELLTYPE)
        normal_pass = bool(not normal.empty and normal_count >= MIN_CELLS_PER_SAMPLE_CELLTYPE)
        reasons = []
        if tumor.empty:
            reasons.append("tumor sample missing")
        elif not tumor_pass:
            reasons.append("tumor epithelial <50")
        if normal.empty:
            reasons.append("normal sample missing")
        elif not normal_pass:
            reasons.append("normal epithelial <50")
        included = tumor_pass and normal_pass
        rows.append({
            "Patient": patient,
            "Tumor Sample": tumor.iloc[0]["Sample"] if not tumor.empty else "",
            "Normal Sample": normal.iloc[0]["Sample"] if not normal.empty else "",
            "Tumor epithelial cells": tumor_count,
            "Normal epithelial cells": normal_count,
            "Tumor pass50": tumor_pass,
            "Normal pass50": normal_pass,
            "paired included/excluded": "included" if included else "excluded",
            "exclusion reason": "; ".join(reasons) if reasons else "",
        })
    return pd.DataFrame(rows)


def gse200997_annotation_status_summary(cells: pd.DataFrame) -> pd.DataFrame:
    """Condition-stratified annotation failure and epithelial retention audit."""
    source = cells[cells["cohort"].eq("GSE200997") & cells["qc_pass"]].copy()
    rows = []
    for condition in ["Tumor", "Normal"]:
        current = source[source["Class"].eq(condition)]
        total = len(current)
        rows.append({
            "cohort": "GSE200997",
            "Class": condition,
            "qc_pass_cells": int(total),
            "confident_pct": float(current["annotation_status"].eq("Confident").mean() * 100) if total else np.nan,
            "low_confidence_pct": float(current["annotation_status"].eq("Low_confidence").mean() * 100) if total else np.nan,
            "ambiguous_pct": float(current["annotation_status"].eq("Doublet_like_or_ambiguous").mean() * 100) if total else np.nan,
            "epithelial_confident_pct": float(current["analysis_cell_type"].eq("Epithelial cells").mean() * 100) if total else np.nan,
            "confident_cells": int(current["annotation_status"].eq("Confident").sum()),
            "low_confidence_cells": int(current["annotation_status"].eq("Low_confidence").sum()),
            "ambiguous_cells": int(current["annotation_status"].eq("Doublet_like_or_ambiguous").sum()),
            "confident_epithelial_cells": int(current["analysis_cell_type"].eq("Epithelial cells").sum()),
        })
    return pd.DataFrame(rows)


def gse200997_label_composition(cells: pd.DataFrame) -> pd.DataFrame:
    """Export top1 and final confident label composition by condition."""
    source = cells[cells["cohort"].eq("GSE200997") & cells["qc_pass"]].copy()
    rows = []
    for condition in ["Tumor", "Normal"]:
        current = source[source["Class"].eq(condition)]
        top1 = current["top1_reference_label"].value_counts(dropna=False)
        for label, count in top1.items():
            rows.append({
                "cohort": "GSE200997", "Class": condition,
                "mapping_scope": "Top1 reference label among all QC cells",
                "cell_type": label, "n_cells": int(count),
                "pct_within_condition": float(count / len(current) * 100) if len(current) else np.nan,
            })
        confident = current[current["analysis_cell_type"].ne("")]
        final = confident["analysis_cell_type"].value_counts(dropna=False)
        for label, count in final.items():
            rows.append({
                "cohort": "GSE200997", "Class": condition,
                "mapping_scope": "Confident final labels",
                "cell_type": label, "n_cells": int(count),
                "pct_within_condition": float(count / len(current) * 100) if len(current) else np.nan,
            })
    return pd.DataFrame(rows)


def assert_unique_patient_class(pseudobulk: pd.DataFrame) -> None:
    """Prevent silent set_index behavior when a patient has duplicate condition rows."""
    usable = pseudobulk[pseudobulk["passes_min_cells"]].copy()
    duplicates = (
        usable.groupby(["cohort", "cell_type", "Patient", "Class"], dropna=False)
        .size()
        .reset_index(name="n_rows")
    )
    duplicates = duplicates[duplicates["n_rows"].gt(1)]
    if not duplicates.empty:
        raise AssertionError(
            "Duplicate Patient + Class records entering paired analysis:\n"
            + duplicates.to_string(index=False)
        )


def paired_effect(pseudobulk: pd.DataFrame, cohort: str, cell_type: str) -> dict:
    current = pseudobulk[
        pseudobulk["passes_min_cells"] &
        pseudobulk["cohort"].eq(cohort) &
        pseudobulk["cell_type"].eq(cell_type)
    ].copy()
    duplicate_rows = current.groupby(["Patient", "Class"], dropna=False).size()
    if duplicate_rows.gt(1).any():
        raise AssertionError(
            f"Duplicate Patient + Class records entering paired analysis for {cohort}/{cell_type}: "
            f"{duplicate_rows[duplicate_rows.gt(1)].to_dict()}"
        )
    tumor = current[current["Class"].eq("Tumor")].set_index("Patient")
    normal = current[current["Class"].eq("Normal")].set_index("Patient")
    patients = sorted(set(tumor.index) & set(normal.index))
    tumor_values = tumor.loc[patients, "log1p_pseudobulk_cpm"].to_numpy(float) if patients else np.array([])
    normal_values = normal.loc[patients, "log1p_pseudobulk_cpm"].to_numpy(float) if patients else np.array([])
    statistic, p_value = safe_wilcoxon(normal_values, tumor_values)
    delta = tumor_values - normal_values
    median_delta = float(np.median(delta)) if len(delta) else np.nan
    return {
        "cohort": cohort,
        "gene_symbol": TARGET,
        "cell_type": cell_type,
        "matched_patient_n": len(patients),
        "matched_patients": ";".join(patients),
        "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm": median_delta,
        "mean_delta_tumor_minus_normal_log1p_pseudobulk_cpm": float(np.mean(delta)) if len(delta) else np.nan,
        "tumor_gt_normal_n": int(np.sum(delta > 0)) if len(delta) else 0,
        "tumor_eq_normal_n": int(np.sum(delta == 0)) if len(delta) else 0,
        "tumor_lt_normal_n": int(np.sum(delta < 0)) if len(delta) else 0,
        "wilcoxon_statistic": statistic,
        "wilcoxon_p_value": p_value,
        "minimum_cell_rule": MIN_CELLS_PER_SAMPLE_CELLTYPE,
        "direction": "Tumor higher" if median_delta > 0 else ("Normal higher" if median_delta < 0 else "No difference"),
    }


def paired_validation(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    assert_unique_patient_class(pseudobulk)
    usable = pseudobulk[pseudobulk["passes_min_cells"]].copy()
    rows = []
    for cohort in ["GSE200997", "GSE132465"]:
        for cell_type in REFERENCE_LABELS:
            rows.append(paired_effect(pseudobulk, cohort, cell_type))
    output = pd.DataFrame(rows)
    output["fdr_bh_secondary_celltype_family"] = output.groupby("cohort")["wilcoxon_p_value"].transform(bh_adjust)
    output["direction"] = np.where(
        output["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"].gt(0),
        "Tumor higher",
        np.where(output["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"].lt(0), "Normal higher", "No difference"),
    )
    output["primary_epithelial_hypothesis"] = output["cell_type"].eq("Epithelial cells")
    primary = output[output["primary_epithelial_hypothesis"]].copy()
    primary_fdr = bh_adjust(primary["wilcoxon_p_value"])
    output["fdr_bh_primary_epithelial_across_cohorts"] = np.nan
    output.loc[primary.index, "fdr_bh_primary_epithelial_across_cohorts"] = primary_fdr.to_numpy()
    return output


def threshold_sensitivity_analysis(gse200997: pd.DataFrame, model: dict) -> pd.DataFrame:
    variants = [
        {
            "threshold_variant": "strict",
            "similarity_threshold": model["sim_threshold_5pct_reference"],
            "margin_threshold": model["margin_threshold_5pct_reference"],
            "require_margin": True,
            "rule": "reference 5th percentile similarity + 5th percentile margin",
        },
        {
            "threshold_variant": "relaxed",
            "similarity_threshold": model["sim_threshold_5pct_reference"],
            "margin_threshold": model["margin_threshold_5pct_reference"],
            "require_margin": False,
            "rule": "reference 5th percentile similarity only; margin not required",
        },
        {
            "threshold_variant": "very_strict",
            "similarity_threshold": model["sim_threshold_10pct_reference"],
            "margin_threshold": model["margin_threshold_10pct_reference"],
            "require_margin": True,
            "rule": "reference 10th percentile similarity + 10th percentile margin",
        },
    ]
    rows = []
    for variant in variants:
        current = apply_reference_threshold_variant(
            gse200997,
            variant["similarity_threshold"],
            variant["margin_threshold"],
            variant["require_margin"],
        )
        pseudobulk = summarize_pseudobulk(current)
        effect = paired_effect(pseudobulk, "GSE200997", "Epithelial cells")
        epithelial_cells = current[
            current["qc_pass"] & current["analysis_cell_type"].eq("Epithelial cells")
        ]
        rows.append({
            "threshold_variant": variant["threshold_variant"],
            "rule": variant["rule"],
            "similarity_threshold": variant["similarity_threshold"],
            "margin_threshold": variant["margin_threshold"],
            "require_margin": variant["require_margin"],
            "confident_epithelial_cells_total": int(len(epithelial_cells)),
            "confident_epithelial_cells_tumor": int((epithelial_cells["Class"] == "Tumor").sum()),
            "confident_epithelial_cells_normal": int((epithelial_cells["Class"] == "Normal").sum()),
            "matched_patient_n": effect["matched_patient_n"],
            "matched_patients": effect["matched_patients"],
            "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm": effect["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"],
            "wilcoxon_p_value": effect["wilcoxon_p_value"],
            "direction": effect["direction"],
        })
    return pd.DataFrame(rows)


def marker_gate_epithelial_sensitivity(gse200997: pd.DataFrame) -> pd.DataFrame:
    """A non-transfer, conservative epithelial gate for sensitivity auditing."""
    epithelial_score = gse200997["marker_score_Epithelial cells"]
    non_epithelial_columns = [
        "marker_score_T cells", "marker_score_B cells", "marker_score_Myeloids",
        "marker_score_Stromal cells", "marker_score_Mast cells",
    ]
    non_epithelial_score = gse200997[non_epithelial_columns].max(axis=1)
    gate = epithelial_score.ge(0.5) & epithelial_score.ge(non_epithelial_score + 0.1)
    current = gse200997.copy()
    current["analysis_cell_type"] = np.where(gate, "Epithelial cells", "")
    pseudobulk = summarize_pseudobulk(current)
    effect = paired_effect(pseudobulk, "GSE200997", "Epithelial cells")
    gated = current[current["qc_pass"] & gate]
    return pd.DataFrame([{
        "method": "canonical_marker_gate",
        "markers": "EPCAM/KRT8/KRT18/KRT19/KRT20/CEACAM5/MUC1",
        "epithelial_score_threshold": 0.5,
        "non_epithelial_exclusion_margin": 0.1,
        "cebpb_excluded_from_gate": True,
        "gated_epithelial_cells_total": int(len(gated)),
        "gated_epithelial_cells_tumor": int((gated["Class"] == "Tumor").sum()),
        "gated_epithelial_cells_normal": int((gated["Class"] == "Normal").sum()),
        "matched_patient_n": effect["matched_patient_n"],
        "matched_patients": effect["matched_patients"],
        "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm": effect["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"],
        "wilcoxon_p_value": effect["wilcoxon_p_value"],
        "direction": effect["direction"],
    }])


def cross_cohort_replication(validation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cell_type in REFERENCE_LABELS:
        current = validation[validation["cell_type"] == cell_type].set_index("cohort")
        row = {"gene_symbol": TARGET, "cell_type": cell_type}
        for cohort in ["GSE200997", "GSE132465"]:
            if cohort in current.index:
                record = current.loc[cohort]
                prefix = cohort.lower()
                row[f"{prefix}_matched_n"] = int(record["matched_patient_n"])
                row[f"{prefix}_median_delta"] = record["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"]
                row[f"{prefix}_p"] = record["wilcoxon_p_value"]
                row[f"{prefix}_direction"] = record["direction"]
            else:
                prefix = cohort.lower()
                row[f"{prefix}_matched_n"] = 0
                row[f"{prefix}_median_delta"] = np.nan
                row[f"{prefix}_p"] = np.nan
                row[f"{prefix}_direction"] = "Insufficient"
        d1 = row["gse200997_median_delta"]
        d2 = row["gse132465_median_delta"]
        agree = bool(np.isfinite(d1) and np.isfinite(d2) and ((d1 > 0 and d2 > 0) or (d1 < 0 and d2 < 0) or (d1 == 0 and d2 == 0)))
        row["same_direction"] = agree
        p1 = row["gse200997_p"]
        p2 = row["gse132465_p"]
        if not agree:
            row["replication_status"] = "Not concordant/insufficient"
        elif np.isfinite(p1) and np.isfinite(p2) and p1 < 0.05 and p2 < 0.05:
            row["replication_status"] = "Direction and nominal evidence concordant"
        else:
            row["replication_status"] = "Direction concordant; statistical replication not established"
        rows.append(row)
    return pd.DataFrame(rows)


def plot_localization(localization: pd.DataFrame) -> None:
    plot_data = localization.copy()
    plot_data["panel"] = plot_data["cohort"] + " | " + plot_data["Class"]
    panel_order = sorted(plot_data["panel"].unique())
    type_order = REFERENCE_LABELS
    x_map = {value: index for index, value in enumerate(type_order)}
    y_map = {value: index for index, value in enumerate(panel_order)}
    fig, ax = plt.subplots(figsize=(12, 6.5), constrained_layout=True)
    current = plot_data.copy()
    current["x"] = current["cell_type"].map(x_map)
    current["y"] = current["panel"].map(y_map)
    scatter = ax.scatter(
        current["x"], current["y"],
        s=30 + 420 * current["detection_rate"].fillna(0),
        c=current["mean_log1p_cpm"], cmap="viridis", edgecolors="black", linewidths=0.25,
    )
    ax.set_xticks(range(len(type_order)))
    ax.set_xticklabels(type_order, rotation=30, ha="right")
    ax.set_yticks(range(len(panel_order)))
    ax.set_yticklabels(panel_order)
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Cohort | condition")
    ax.set_title("CEBPB localization after clean cohort rerun\nDot size = detection rate; color = mean log1p CPM")
    ax.grid(axis="x", alpha=0.15)
    fig.colorbar(scatter, ax=ax, label="Mean log1p CPM")
    fig.savefig(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_localization.png", dpi=240)
    plt.close(fig)


def plot_epithelial_validation(pseudobulk: pd.DataFrame) -> None:
    current = pseudobulk[
        pseudobulk["passes_min_cells"] &
        pseudobulk["cell_type"].eq("Epithelial cells") &
        pseudobulk["Class"].isin(["Tumor", "Normal"])
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True, constrained_layout=True)
    for ax, cohort in zip(axes, ["GSE200997", "GSE132465"]):
        cohort_data = current[current["cohort"] == cohort]
        tumor = cohort_data[cohort_data["Class"] == "Tumor"].set_index("Patient")
        normal = cohort_data[cohort_data["Class"] == "Normal"].set_index("Patient")
        patients = sorted(set(tumor.index) & set(normal.index))
        for index, patient in enumerate(patients):
            ax.plot(
                [0, 1],
                [normal.loc[patient, "log1p_pseudobulk_cpm"], tumor.loc[patient, "log1p_pseudobulk_cpm"]],
                color="#6b7280", alpha=0.45, linewidth=0.8,
            )
        if patients:
            ax.scatter(np.zeros(len(patients)), normal.loc[patients, "log1p_pseudobulk_cpm"], color="#2563eb", s=22, label="Normal")
            ax.scatter(np.ones(len(patients)), tumor.loc[patients, "log1p_pseudobulk_cpm"], color="#dc2626", s=22, label="Tumor")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Normal", "Tumor"])
        ax.set_title(f"{cohort} (matched n={len(patients)})")
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("CEBPB log1p pseudobulk CPM")
    axes[-1].legend(frameon=False)
    fig.suptitle("CEBPB in curated/reference-mapped epithelial cells")
    fig.savefig(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_epithelial_tumor_normal.png", dpi=240)
    plt.close(fig)


def fmt_report_number(value: object, digits: int = 4) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}g}"


def write_report(
    cells: pd.DataFrame,
    localization: pd.DataFrame,
    pseudobulk: pd.DataFrame,
    validation: pd.DataFrame,
    replication: pd.DataFrame,
    annotation_audit: pd.DataFrame,
    epithelial_count_audit: pd.DataFrame,
    inclusion_audit: pd.DataFrame,
    status_summary: pd.DataFrame,
    label_composition: pd.DataFrame,
    threshold_sensitivity: pd.DataFrame,
    marker_gate_sensitivity: pd.DataFrame,
    manifest: dict,
) -> None:
    primary = validation[validation["primary_epithelial_hypothesis"]].copy()
    report = [
        "# DINP–CRC clean CEBPB single-cell rerun",
        "",
        "## Frozen design",
        "",
        "This rerun is CEBPB-only and uses GSE200997 plus GSE132465. GSE188711 was excluded. No sidedness/left-right field is read, derived, grouped, or exported.",
        "",
        "- GSE132465: raw UMI matrix reconstructed independently from `work/scRNA_GSE132465/raw_UMI_count_matrix.txt.gz`; curated `cell_type` labels were read from the source annotation file.",
        "- GSE200997: raw UMI matrix reconstructed independently from `work/scRNA_GSE200997/raw_UMI_count_matrix.csv.gz`; its source annotation was used only for cell ID, sample, and Tumor/Normal condition.",
        "- GSE200997 cell types: nearest-centroid mapping to GSE132465 curated broad labels using canonical marker expression, with reference-calibrated similarity and margin thresholds. Low-confidence and ambiguous cells are retained in the audit table but excluded from the primary pseudobulk analysis.",
        "- Primary validation: matched patient-level epithelial pseudobulk, Tumor versus Normal, paired Wilcoxon, minimum 50 cells per sample × cell type.",
        "- Secondary localization/compartment family: BH-FDR across the six prespecified cell types within each cohort. Primary epithelial family: BH across the two epithelial cohort tests only.",
        "",
        "## QC and annotation audit",
        "",
        f"- Total cells loaded: {len(cells):,}; QC-passing cells: {int(cells['qc_pass'].sum()):,}.",
        f"- QC thresholds: total UMI ≥ {MIN_QC_UMI}; detected genes ≥ {MIN_QC_GENES}.",
        f"- GSE200997 cells retained as confident reference-mapped labels: {int((cells['cohort'].eq('GSE200997') & cells['annotation_status'].eq('Confident')).sum()):,}; low-confidence: {int((cells['cohort'].eq('GSE200997') & cells['annotation_status'].eq('Low_confidence')).sum()):,}; ambiguous/doublet-like: {int((cells['cohort'].eq('GSE200997') & cells['annotation_status'].eq('Doublet_like_or_ambiguous')).sum()):,}.",
        f"- Reference LOO median-centroid calibration: accuracy={manifest['reference_mapping']['reference_loo_label_accuracy']:.3f}; top1 median={manifest['reference_mapping']['reference_loo_top1_median']:.4f}; margin median={manifest['reference_mapping']['reference_loo_margin_median']:.4f}.",
        f"- Strict thresholds: similarity={manifest['reference_mapping']['sim_threshold_5pct']:.4f}; margin={manifest['reference_mapping']['margin_threshold_5pct']:.4f}. Very-strict thresholds: similarity={manifest['reference_mapping']['sim_threshold_10pct']:.4f}; margin={manifest['reference_mapping']['margin_threshold_10pct']:.4f}.",
        "- GSE200997 epithelial labels are therefore reference-mapped epithelial cells, not CNV-confirmed malignant cells; the Tumor sample context is not treated as a malignancy annotation.",
        "",
        "## GSE200997 condition-dependent annotation audit",
        "",
    ]
    for _, row in status_summary.iterrows():
        report.append(
            f"- {row['Class']}: QC cells={int(row['qc_pass_cells'])}; Confident={row['confident_pct']:.2f}%; Low-confidence={row['low_confidence_pct']:.2f}%; Ambiguous={row['ambiguous_pct']:.2f}%; Confident epithelial={row['epithelial_confident_pct']:.2f}%.",
        )
    report.extend([
        "",
        "The full label composition is exported separately for both top1 reference labels among all QC cells and confident final labels. These condition-stratified percentages are the audit for potential condition-dependent filtering.",
        "",
        "## GSE200997 epithelial inclusion audit",
        "",
        "The prior contaminated module reported 7 matched epithelial patients. In this clean strict reference-mapped analysis, only patients passing the 50-cell rule in both conditions enter the paired test; the table below records every exclusion reason.",
    ])
    for _, row in inclusion_audit.iterrows():
        reason = row["exclusion reason"] or "both conditions pass 50 cells"
        report.append(
            f"- {row['Patient']}: Tumor={int(row['Tumor epithelial cells'])}, Normal={int(row['Normal epithelial cells'])}; Tumor pass50={str(row['Tumor pass50']).lower()}, Normal pass50={str(row['Normal pass50']).lower()}; {row['paired included/excluded']}; {reason}.",
        )
    report.extend([
        "",
        "## Threshold sensitivity",
        "",
    ])
    for _, row in threshold_sensitivity.iterrows():
        report.append(
            f"- {row['threshold_variant']}: matched n={int(row['matched_patient_n'])}; median Δ={fmt_report_number(row['median_delta_tumor_minus_normal_log1p_pseudobulk_cpm'])}; P={fmt_report_number(row['wilcoxon_p_value'])}; direction={row['direction']}; included epithelial cells={int(row['confident_epithelial_cells_total'])}.",
        )
    marker_row = marker_gate_sensitivity.iloc[0]
    report.extend([
        "",
        "## Non-transfer epithelial gate sensitivity",
        "",
        f"A separate canonical-marker gate excluding CEBPB from classification retained {int(marker_row['gated_epithelial_cells_total'])} cells and produced matched n={int(marker_row['matched_patient_n'])}, median Δ={fmt_report_number(marker_row['median_delta_tumor_minus_normal_log1p_pseudobulk_cpm'])}, P={fmt_report_number(marker_row['wilcoxon_p_value'])}, direction={marker_row['direction']}.",
        "",
        "## Interpretation",
        "",
        "GSE132465 showed a tumor-high epithelial CEBPB signal. GSE200997 had the same nominal direction, but did not provide statistically supportive replication under strict reference-mapped epithelial analysis (matched n=4; P=0.625); therefore, cross-cohort epithelial replication was not established.",
        "",
        "## Primary epithelial validation",
        "",
    ])
    for _, row in primary.iterrows():
        report.append(
            f"- {row['cohort']}: matched n={int(row['matched_patient_n'])}; median Δ(Tumor−Normal)={fmt_report_number(row['median_delta_tumor_minus_normal_log1p_pseudobulk_cpm'])}; Wilcoxon P={fmt_report_number(row['wilcoxon_p_value'])}; primary epithelial BH-FDR={fmt_report_number(row['fdr_bh_primary_epithelial_across_cohorts'])}; {row['direction']}.",
        )
    epithelial_rep = replication[replication["cell_type"].eq("Epithelial cells")]
    if not epithelial_rep.empty:
        report.extend(["", "## Cross-cohort epithelial replication", ""])
        report.append(
            f"- Direction concordance: {bool(epithelial_rep.iloc[0]['same_direction'])}; status: {epithelial_rep.iloc[0]['replication_status']}."
        )
    report.extend([
        "",
        "## Files",
        "",
        "- `DINP_CRC_clean_CEBPB_cell_level.csv`: cell-level CEBPB metrics and annotation/QC status.",
        "- `DINP_CRC_clean_GSE200997_annotation_audit.csv`: complete GSE200997 reference-mapping audit, including top-two labels, similarities, margins, and confidence status.",
        "- `DINP_CRC_clean_CEBPB_patient_epithelial_cell_count_audit.csv`: patient × sample × condition QC and epithelial cell-count audit for both cohorts.",
        "- `DINP_CRC_clean_GSE200997_epithelial_inclusion_audit.csv`: explicit GSE200997 paired inclusion/exclusion reasons explaining matched n.",
        "- `DINP_CRC_clean_GSE200997_annotation_status_by_condition.csv`: Tumor/Normal condition-stratified confidence and epithelial retention audit.",
        "- `DINP_CRC_clean_GSE200997_label_composition_by_condition.csv`: top1 and confident final label composition by condition.",
        "- `DINP_CRC_clean_GSE200997_epithelial_threshold_sensitivity.csv`: strict/relaxed/very-strict reference threshold sensitivity.",
        "- `DINP_CRC_clean_GSE200997_marker_gate_epithelial_sensitivity.csv`: non-transfer canonical epithelial-gate sensitivity.",
        "- `DINP_CRC_clean_CEBPB_localization.csv`: cohort/condition/cell-type localization summary.",
        "- `DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv`: patient/sample-level pseudobulk.",
        "- `DINP_CRC_clean_CEBPB_tumor_normal_validation_all_celltypes.csv`: secondary six-cell-type paired validation with cohort-specific BH-FDR.",
        "- `DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv`: prespecified epithelial validation with BH-FDR across the two cohorts.",
        "- `DINP_CRC_clean_CEBPB_cross_cohort_replication.csv`: cross-cohort direction/replication audit.",
        "- `DINP_CRC_clean_CEBPB_manifest.json`: source hashes, thresholds, and exact analysis settings.",
        "",
        "The previous three-cohort CEBPB/CD36 module remains frozen and is not used as input here.",
    ])
    (OUTPUT_DIR / "DINP_CRC_clean_CEBPB_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    gse132465, manifest_132465, counts_132465, total_umi_132465 = load_gse132465()
    model, _, _, _ = fit_reference_mapping(gse132465, counts_132465, total_umi_132465)
    gse200997, manifest_200997, annotation_audit = load_gse200997(model)
    cells = pd.concat([gse132465, gse200997], ignore_index=True, sort=False).fillna("")

    epithelial_count_audit = patient_epithelial_cell_count_audit(cells)
    inclusion_audit = gse200997_epithelial_inclusion_audit(epithelial_count_audit)
    status_summary = gse200997_annotation_status_summary(cells)
    label_composition = gse200997_label_composition(cells)
    threshold_sensitivity = threshold_sensitivity_analysis(gse200997, model)
    marker_gate_sensitivity = marker_gate_epithelial_sensitivity(gse200997)
    localization = summarize_localization(cells)
    pseudobulk = summarize_pseudobulk(cells)
    validation = paired_validation(pseudobulk)
    primary = validation[validation["primary_epithelial_hypothesis"]].copy()
    replication = cross_cohort_replication(validation)

    cell_columns = [
        "cohort", "cell_id", "Patient", "Sample", "Class", "cell_type", "analysis_cell_type", "cell_subtype",
        "cell_type_source", "annotation_status", "total_umi", "n_genes", "qc_pass",
        f"{TARGET}_umi", f"{TARGET}_detected", f"{TARGET}_cpm", f"{TARGET}_log1p_cpm",
    ]
    extra_columns = [column for column in cells.columns if column.startswith("marker_score_") or column in {
        "top1_reference_label", "top2_reference_label", "top1_reference_similarity", "top2_reference_similarity",
        "reference_score_margin", "marker_panel_count_gt_0_25",
    }]
    cells[cell_columns + extra_columns].to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_cell_level.csv", index=False)
    annotation_audit.to_csv(OUTPUT_DIR / "DINP_CRC_clean_GSE200997_annotation_audit.csv", index=False)
    epithelial_count_audit.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_patient_epithelial_cell_count_audit.csv", index=False)
    inclusion_audit.to_csv(OUTPUT_DIR / "DINP_CRC_clean_GSE200997_epithelial_inclusion_audit.csv", index=False)
    status_summary.to_csv(OUTPUT_DIR / "DINP_CRC_clean_GSE200997_annotation_status_by_condition.csv", index=False)
    label_composition.to_csv(OUTPUT_DIR / "DINP_CRC_clean_GSE200997_label_composition_by_condition.csv", index=False)
    threshold_sensitivity.to_csv(OUTPUT_DIR / "DINP_CRC_clean_GSE200997_epithelial_threshold_sensitivity.csv", index=False)
    marker_gate_sensitivity.to_csv(OUTPUT_DIR / "DINP_CRC_clean_GSE200997_marker_gate_epithelial_sensitivity.csv", index=False)
    localization.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_localization.csv", index=False)
    pseudobulk.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv", index=False)
    validation.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_tumor_normal_validation_all_celltypes.csv", index=False)
    primary.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv", index=False)
    replication.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_cross_cohort_replication.csv", index=False)
    plot_localization(localization)
    plot_epithelial_validation(pseudobulk)

    manifest = {
        "analysis": "clean_CEBPB_only_GSE200997_GSE132465",
        "targets": [TARGET],
        "excluded_cohorts": ["GSE188711"],
        "excluded_metadata_fields": ["Location", "Side", "left_right", "tumor_region"],
        "qc_threshold": {"min_umi": MIN_QC_UMI, "min_genes": MIN_QC_GENES},
        "pseudobulk": {
            "normalization": "aggregate CEBPB UMI / aggregate library UMI × 1e6 = pseudobulk CPM",
            "grouping": "cohort × patient × sample × class × cell_type",
            "minimum_cells_per_sample_celltype": MIN_CELLS_PER_SAMPLE_CELLTYPE,
            "tumor_normal_test": "paired Wilcoxon across matched patients",
            "secondary_fdr_family": "six prespecified cell types within each cohort",
            "primary_fdr_family": "epithelial tests across the two cohorts",
        },
        "reference_mapping": {
            "reference_cohort": "GSE132465",
            "labels": REFERENCE_LABELS,
            "features": model["features"],
            "similarity": "cosine similarity after gene-wise standardization fit on GSE132465 marker expression",
            "confidence_rule": "Confident if top1 similarity >= reference 5th percentile and top1-top2 margin >= reference 5th percentile; otherwise low-confidence or ambiguous if top2 also clears similarity threshold",
            "calibration_estimator": "exact leave-one-out median centroid; production classifier also uses median centroid",
            "sim_threshold_5pct": model["sim_threshold_5pct_reference"],
            "margin_threshold_5pct": model["margin_threshold_5pct_reference"],
            "sim_threshold_10pct": model["sim_threshold_10pct_reference"],
            "margin_threshold_10pct": model["margin_threshold_10pct_reference"],
            "reference_loo_top1_median": model["reference_loo_top1_median"],
            "reference_loo_margin_median": model["reference_loo_margin_median"],
            "reference_loo_label_accuracy": model["reference_loo_label_accuracy"],
        },
        "cohorts": {"GSE132465": manifest_132465, "GSE200997": manifest_200997},
        "cell_counts": {
            "combined_loaded": int(len(cells)),
            "combined_qc_pass": int(cells["qc_pass"].sum()),
            "GSE132465_qc_pass": int((cells["cohort"].eq("GSE132465") & cells["qc_pass"]).sum()),
            "GSE200997_qc_pass": int((cells["cohort"].eq("GSE200997") & cells["qc_pass"]).sum()),
            "GSE200997_annotation_status": cells.loc[cells["cohort"].eq("GSE200997"), "annotation_status"].value_counts().to_dict(),
            "GSE200997_annotation_status_by_condition": status_summary.to_dict("records"),
            "GSE200997_strict_epithelial_included_patients": inclusion_audit.loc[inclusion_audit["paired included/excluded"].eq("included"), "Patient"].tolist(),
        },
        "output_rows": {
            "cell_level": int(len(cells)),
            "annotation_audit": int(len(annotation_audit)),
            "patient_epithelial_cell_count_audit": int(len(epithelial_count_audit)),
            "gse200997_epithelial_inclusion_audit": int(len(inclusion_audit)),
            "gse200997_annotation_status_by_condition": int(len(status_summary)),
            "gse200997_label_composition_by_condition": int(len(label_composition)),
            "gse200997_threshold_sensitivity": int(len(threshold_sensitivity)),
            "gse200997_marker_gate_sensitivity": int(len(marker_gate_sensitivity)),
            "localization": int(len(localization)),
            "pseudobulk": int(len(pseudobulk)),
            "validation": int(len(validation)),
            "primary_epithelial": int(len(primary)),
            "replication": int(len(replication)),
        },
        "source_note": "The old CEBPB/CD36 multicohort output was not read as input; both cohorts were reconstructed from raw matrices and source annotation files.",
    }
    (OUTPUT_DIR / "DINP_CRC_clean_CEBPB_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    write_report(
        cells, localization, pseudobulk, validation, replication, annotation_audit,
        epithelial_count_audit, inclusion_audit, status_summary, label_composition,
        threshold_sensitivity, marker_gate_sensitivity, manifest,
    )
    print(json.dumps({
        "output_dir": str(OUTPUT_DIR),
        "cells": len(cells),
        "qc_pass": int(cells["qc_pass"].sum()),
        "gse200997_status": cells.loc[cells["cohort"].eq("GSE200997"), "annotation_status"].value_counts().to_dict(),
        "gse200997_status_by_condition": status_summary.to_dict("records"),
        "gse200997_inclusion_audit": inclusion_audit.to_dict("records"),
        "threshold_sensitivity": threshold_sensitivity.to_dict("records"),
        "marker_gate_sensitivity": marker_gate_sensitivity.to_dict("records"),
        "primary_epithelial": primary[["cohort", "matched_patient_n", "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm", "wilcoxon_p_value", "fdr_bh_primary_epithelial_across_cohorts", "direction"]].to_dict("records"),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
