from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "work"))
import run_clean_cebpb_scRNA_GSE200997_GSE132465 as base  # noqa: E402


TARGET = base.TARGET
REFERENCE_LABELS = base.REFERENCE_LABELS
MIN_CELLS_PER_SAMPLE_CELLTYPE = base.MIN_CELLS_PER_SAMPLE_CELLTYPE
OUTPUT_DIR = ROOT / "outputs" / "DINP_CRC_clean_CEBPB_scRNA_GSE200997_GSE132465_GSE294300"
GSE294300_DIR = ROOT / "work" / "scRNA_GSE294300"
GSE294300_CELL_BATCH = GSE294300_DIR / "GSE294300_cell_batch.tsv.gz"
GSE294300_RAW_DIR = Path(r"F:\GSE294300_10x")

# The accession-to-batch mapping is taken from the official GSE294300 sample
# records.  The batch strings intentionally preserve the FTP filenames because
# capitalization/underscore differences are present in the deposited files.
GSE294300_SAMPLES = {
    **{f"nor_{p}": 8901920 + i for i, p in enumerate(["01", "02", "03"])},
    "Nor_04": 8901923,
    "Nor6": 8901924,
    "Nor7": 8901925,
    "Nor8": 8901926,
    **{f"nor_{p}": 8901927 + i for i, p in enumerate(["09", "10", "11", "15", "16", "17", "18", "19", "20", "21", "23"])},
    **{f"tum_{p}": 8901938 + i for i, p in enumerate(["01", "02", "03"])},
    "Tum_04": 8901941,
    "Tum6": 8901942,
    "Tum7": 8901943,
    "Tum8": 8901944,
    **{f"tum_{p}": 8901945 + i for i, p in enumerate(["09", "10", "11", "15", "16", "17", "18", "19", "20", "21", "23"])},
}


def _sample_files(batch: str) -> tuple[int, Path, Path, Path]:
    if batch not in GSE294300_SAMPLES:
        raise AssertionError(f"No official GSM mapping configured for batch {batch}")
    gsm = GSE294300_SAMPLES[batch]
    prefix = f"GSM{gsm}_{batch}"
    paths = tuple(GSE294300_RAW_DIR / f"{prefix}_{suffix}" for suffix in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"))
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing GSE294300 raw file(s): " + ", ".join(missing))
    return gsm, paths[0], paths[1], paths[2]


def _read_10x_sample(
    batch: str,
    batch_cells: pd.DataFrame,
    model: dict,
) -> tuple[pd.DataFrame, dict]:
    gsm, barcode_path, feature_path, matrix_path = _sample_files(batch)
    raw_barcodes = pd.read_csv(barcode_path, header=None, names=["cell_id"], dtype=str)["cell_id"].str.strip().tolist()
    features = pd.read_csv(feature_path, sep="\t", header=None, dtype=str)
    if features.shape[1] < 2:
        raise AssertionError(f"Unexpected feature format in {feature_path}")
    gene_symbols = features.iloc[:, 1].fillna(features.iloc[:, 0]).astype(str).str.strip().str.upper().tolist()

    with gzip.open(matrix_path, "rt") as handle:
        matrix = mmread(handle).tocsr()
    if matrix.shape[0] != len(gene_symbols) or matrix.shape[1] != len(raw_barcodes):
        raise AssertionError(
            f"10x dimensions do not match metadata for {batch}: matrix={matrix.shape}, "
            f"features={len(gene_symbols)}, barcodes={len(raw_barcodes)}"
        )
    if batch_cells["Cell"].duplicated().any():
        raise AssertionError(f"Duplicate cell IDs in cell_batch metadata for {batch}")
    batch_ids = batch_cells["Cell"].astype(str).tolist()
    # The deposited tumor cell-batch table uses a '-2' barcode suffix while
    # the raw 10x barcode files use '-1'.  Match on the 10x barcode core, but
    # retain the deposited cell ID in all exported tables.
    barcode_key = lambda value: re.sub(r"-\d+$", "", str(value))
    raw_keys = [barcode_key(value) for value in raw_barcodes]
    batch_keys = [barcode_key(value) for value in batch_ids]
    if len(set(raw_keys)) != len(raw_keys) or len(set(batch_keys)) != len(batch_keys):
        raise AssertionError(f"Barcode cores are not unique for {batch}")
    if not set(batch_keys).issubset(set(raw_keys)):
        missing_in_raw = sorted(set(batch_keys) - set(raw_keys))[:5]
        raise AssertionError(f"cell_batch IDs missing from raw matrix for {batch}: {missing_in_raw}")
    raw_index = {cell_id: i for i, cell_id in enumerate(raw_keys)}
    order = np.asarray([raw_index[cell_id] for cell_id in batch_keys], dtype=np.int64)

    matrix_csc = matrix.tocsc()
    total_umi = np.asarray(matrix_csc.sum(axis=0)).ravel().astype(np.int64)[order]
    n_genes = np.asarray(matrix_csc.getnnz(axis=0)).ravel().astype(np.int32)[order]
    gene_rows: dict[str, list[int]] = {}
    for index, gene in enumerate(gene_symbols):
        if gene in base.COUNT_GENES:
            gene_rows.setdefault(gene, []).append(index)
    missing = base.COUNT_GENES - set(gene_rows)
    if missing:
        raise AssertionError(f"Missing requested genes in {feature_path.name}: {sorted(missing)}")
    gene_counts = {
        gene: np.asarray(matrix[rows, :].sum(axis=0)).ravel().astype(np.int64)[order]
        for gene, rows in gene_rows.items()
    }

    norm_batch = batch.lower().replace("_", "")
    condition = "Normal" if norm_batch.startswith("nor") else "Tumor" if norm_batch.startswith("tum") else ""
    if not condition:
        raise AssertionError(f"Cannot infer Tumor/Normal from GSE294300 batch {batch}")
    digits = "".join(ch for ch in batch if ch.isdigit())
    patient = f"patient_{int(digits):02d}"
    cells = pd.DataFrame({
        "cell_id": batch_ids,
        "Patient": patient,
        "Sample": f"GSE294300_{batch}",
        "Class": condition,
        "cell_subtype": "",
        "cohort": "GSE294300",
    })
    cells = base.add_target_metrics(cells, gene_counts, total_umi, n_genes)
    cells = base.apply_reference_mapping(cells, gene_counts, total_umi, model)
    audit = cells[[
        "cell_id", "Patient", "Sample", "Class", "cohort", "qc_pass", "total_umi", "n_genes",
        "top1_reference_label", "top2_reference_label", "top1_reference_similarity",
        "top2_reference_similarity", "reference_score_margin", "marker_panel_count_gt_0_25",
        "annotation_status", "cell_type", "analysis_cell_type",
    ]].copy()
    source_meta = {
        "batch": batch,
        "gsm": f"GSM{gsm}",
        "patient": patient,
        "class": condition,
        "cells": int(len(cells)),
        "raw_barcodes": int(len(raw_barcodes)),
        "selected_cell_batch_barcodes": int(len(batch_ids)),
        "barcode_matching": "barcode core with terminal '-N' suffix removed; tumor metadata uses -2 while raw 10x files use -1",
        "raw_matrix_subset_rule": "retain only barcodes listed in GSE294300_cell_batch.tsv.gz",
        "barcodes": str(barcode_path),
        "features": str(feature_path),
        "matrix": str(matrix_path),
        "matrix_size_bytes": int(matrix_path.stat().st_size),
    }
    return cells, {"audit": audit, "source": source_meta}


def load_gse294300(model: dict) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    batch_table = pd.read_csv(GSE294300_CELL_BATCH, sep="\t", dtype=str).fillna("")
    batch_table.columns = [str(c).strip() for c in batch_table.columns]
    if not {"Cell", "batch"}.issubset(batch_table.columns):
        raise AssertionError(f"Unexpected GSE294300 cell-batch columns: {batch_table.columns.tolist()}")
    if batch_table["Cell"].duplicated().any():
        raise AssertionError("GSE294300 cell_batch.tsv.gz contains duplicate cell IDs")
    if set(batch_table["batch"].unique()) != set(GSE294300_SAMPLES):
        raise AssertionError(
            "GSE294300 batch metadata does not match configured sample map: "
            f"metadata_only={sorted(set(batch_table['batch']) - set(GSE294300_SAMPLES))}; "
            f"map_only={sorted(set(GSE294300_SAMPLES) - set(batch_table['batch']))}"
        )

    parts = []
    audits = []
    sources = []
    for batch in batch_table["batch"].drop_duplicates().tolist():
        batch_cells = batch_table[batch_table["batch"].eq(batch)][["Cell", "batch"]].copy()
        current, details = _read_10x_sample(batch, batch_cells, model)
        parts.append(current)
        audits.append(details["audit"])
        sources.append(details["source"])
        print(f"Loaded GSE294300 {batch}: {len(current):,} cells", flush=True)
    cells = pd.concat(parts, ignore_index=True, sort=False).fillna("")
    annotation_audit = pd.concat(audits, ignore_index=True, sort=False)
    manifest = {
        "cohort": "GSE294300",
        "cell_batch_source": base.repo_path(GSE294300_CELL_BATCH),
        "cell_batch_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294300/suppl/GSE294300_cell_batch.tsv.gz",
        "cell_batch_sha256": base.sha256(GSE294300_CELL_BATCH),
        "cell_batch_rows": int(len(batch_table)),
        "raw_matrix_cell_selection": "The raw 10x matrices contain unfiltered barcodes; analysis cells are restricted to the exact barcode subset in cell_batch.tsv.gz.",
        "sample_count": int(len(sources)),
        "patient_count": int(batch_table["batch"].str.extract(r"(\d+)$")[0].nunique()),
        "sample_sources": sources,
        "cell_type_source": "Reference-mapped to GSE132465 curated broad labels; GSE294300 deposited metadata supplied batch only",
        "raw_count_source": "NCBI GEO per-sample 10x matrix/features/barcodes supplements",
        "raw_count_url_template": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8901nnn/GSM{gsm}/suppl/GSM{gsm}_{batch}_{kind}",
        "qc_threshold": {"min_umi": base.MIN_QC_UMI, "min_genes": base.MIN_QC_GENES},
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
    }
    return cells, manifest, annotation_audit


def summarize_localization(cells: pd.DataFrame) -> pd.DataFrame:
    usable = cells[cells["qc_pass"] & cells["analysis_cell_type"].ne("")].copy()
    return (
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
    rows = []
    for (cohort, patient, sample, condition), group in cells.groupby(["cohort", "Patient", "Sample", "Class"], dropna=False, sort=True):
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
    result = pd.DataFrame(rows)
    result["epithelial_pass_50"] = result["confident_epithelial_cells"].ge(MIN_CELLS_PER_SAMPLE_CELLTYPE)
    return result


def epithelial_inclusion_audit(count_audit: pd.DataFrame, cohort: str) -> pd.DataFrame:
    source = count_audit[count_audit["cohort"].eq(cohort)].copy()
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
        rows.append({
            "cohort": cohort,
            "Patient": patient,
            "Tumor Sample": tumor.iloc[0]["Sample"] if not tumor.empty else "",
            "Normal Sample": normal.iloc[0]["Sample"] if not normal.empty else "",
            "Tumor epithelial cells": tumor_count,
            "Normal epithelial cells": normal_count,
            "Tumor pass50": tumor_pass,
            "Normal pass50": normal_pass,
            "paired included/excluded": "included" if tumor_pass and normal_pass else "excluded",
            "exclusion reason": "; ".join(reasons) if reasons else "",
        })
    return pd.DataFrame(rows)


def annotation_status_summary(cells: pd.DataFrame, cohorts: list[str]) -> pd.DataFrame:
    rows = []
    for cohort in cohorts:
        source = cells[cells["cohort"].eq(cohort) & cells["qc_pass"]]
        for condition in ["Tumor", "Normal"]:
            current = source[source["Class"].eq(condition)]
            total = len(current)
            rows.append({
                "cohort": cohort,
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


def label_composition(cells: pd.DataFrame, cohorts: list[str]) -> pd.DataFrame:
    rows = []
    for cohort in cohorts:
        source = cells[cells["cohort"].eq(cohort) & cells["qc_pass"]]
        for condition in ["Tumor", "Normal"]:
            current = source[source["Class"].eq(condition)]
            for scope, values in [
                ("Top1 reference label among all QC cells", current["top1_reference_label"]),
                ("Confident final labels", current.loc[current["analysis_cell_type"].ne(""), "analysis_cell_type"]),
            ]:
                for label, count in values.value_counts(dropna=False).items():
                    rows.append({
                        "cohort": cohort, "Class": condition, "mapping_scope": scope,
                        "cell_type": label, "n_cells": int(count),
                        "pct_within_condition": float(count / len(current) * 100) if len(current) else np.nan,
                    })
    return pd.DataFrame(rows)


def assert_unique_patient_class(pseudobulk: pd.DataFrame) -> None:
    usable = pseudobulk[pseudobulk["passes_min_cells"]]
    duplicates = usable.groupby(["cohort", "cell_type", "Patient", "Class"], dropna=False).size().reset_index(name="n_rows")
    duplicates = duplicates[duplicates["n_rows"].gt(1)]
    if not duplicates.empty:
        raise AssertionError("Duplicate Patient + Class records entering paired analysis:\n" + duplicates.to_string(index=False))


def paired_effect(pseudobulk: pd.DataFrame, cohort: str, cell_type: str) -> dict:
    current = pseudobulk[
        pseudobulk["passes_min_cells"] & pseudobulk["cohort"].eq(cohort) & pseudobulk["cell_type"].eq(cell_type)
    ].copy()
    duplicate_rows = current.groupby(["Patient", "Class"], dropna=False).size()
    if duplicate_rows.gt(1).any():
        raise AssertionError(f"Duplicate Patient + Class records for {cohort}/{cell_type}: {duplicate_rows[duplicate_rows.gt(1)].to_dict()}")
    tumor = current[current["Class"].eq("Tumor")].set_index("Patient")
    normal = current[current["Class"].eq("Normal")].set_index("Patient")
    patients = sorted(set(tumor.index) & set(normal.index))
    tumor_values = tumor.loc[patients, "log1p_pseudobulk_cpm"].to_numpy(float) if patients else np.array([])
    normal_values = normal.loc[patients, "log1p_pseudobulk_cpm"].to_numpy(float) if patients else np.array([])
    statistic, p_value = base.safe_wilcoxon(normal_values, tumor_values)
    delta = tumor_values - normal_values
    median_delta = float(np.median(delta)) if len(delta) else np.nan
    return {
        "cohort": cohort, "gene_symbol": TARGET, "cell_type": cell_type,
        "matched_patient_n": len(patients), "matched_patients": ";".join(patients),
        "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm": median_delta,
        "mean_delta_tumor_minus_normal_log1p_pseudobulk_cpm": float(np.mean(delta)) if len(delta) else np.nan,
        "tumor_gt_normal_n": int(np.sum(delta > 0)) if len(delta) else 0,
        "tumor_eq_normal_n": int(np.sum(delta == 0)) if len(delta) else 0,
        "tumor_lt_normal_n": int(np.sum(delta < 0)) if len(delta) else 0,
        "wilcoxon_statistic": statistic, "wilcoxon_p_value": p_value,
        "minimum_cell_rule": MIN_CELLS_PER_SAMPLE_CELLTYPE,
        "direction": "Tumor higher" if median_delta > 0 else ("Normal higher" if median_delta < 0 else "No difference"),
    }


def paired_validation(pseudobulk: pd.DataFrame, cohorts: list[str]) -> pd.DataFrame:
    assert_unique_patient_class(pseudobulk)
    rows = [paired_effect(pseudobulk, cohort, cell_type) for cohort in cohorts for cell_type in REFERENCE_LABELS]
    output = pd.DataFrame(rows)
    output["fdr_bh_secondary_celltype_family"] = output.groupby("cohort")["wilcoxon_p_value"].transform(base.bh_adjust)
    output["primary_epithelial_hypothesis"] = output["cell_type"].eq("Epithelial cells")
    primary = output[output["primary_epithelial_hypothesis"]]
    output["fdr_bh_primary_epithelial_across_cohorts"] = np.nan
    output.loc[primary.index, "fdr_bh_primary_epithelial_across_cohorts"] = base.bh_adjust(primary["wilcoxon_p_value"]).to_numpy()
    return output


def threshold_sensitivity_analysis(mapped_cells: pd.DataFrame, model: dict, cohort: str) -> pd.DataFrame:
    variants = [
        ("strict", model["sim_threshold_5pct_reference"], model["margin_threshold_5pct_reference"], True, "reference 5th percentile similarity + 5th percentile margin"),
        ("relaxed", model["sim_threshold_5pct_reference"], model["margin_threshold_5pct_reference"], False, "reference 5th percentile similarity only; margin not required"),
        ("very_strict", model["sim_threshold_10pct_reference"], model["margin_threshold_10pct_reference"], True, "reference 10th percentile similarity + 10th percentile margin"),
    ]
    rows = []
    for name, sim, margin, require_margin, rule in variants:
        current = base.apply_reference_threshold_variant(mapped_cells, sim, margin, require_margin)
        effect = paired_effect(summarize_pseudobulk(current), cohort, "Epithelial cells")
        epi = current[current["qc_pass"] & current["analysis_cell_type"].eq("Epithelial cells")]
        rows.append({
            "cohort": cohort, "threshold_variant": name, "rule": rule,
            "similarity_threshold": sim, "margin_threshold": margin, "require_margin": require_margin,
            "confident_epithelial_cells_total": int(len(epi)),
            "confident_epithelial_cells_tumor": int((epi["Class"] == "Tumor").sum()),
            "confident_epithelial_cells_normal": int((epi["Class"] == "Normal").sum()),
            "matched_patient_n": effect["matched_patient_n"], "matched_patients": effect["matched_patients"],
            "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm": effect["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"],
            "wilcoxon_p_value": effect["wilcoxon_p_value"], "direction": effect["direction"],
        })
    return pd.DataFrame(rows)


def marker_gate_sensitivity(mapped_cells: pd.DataFrame, cohort: str) -> pd.DataFrame:
    epithelial_score = mapped_cells["marker_score_Epithelial cells"]
    non_epi = mapped_cells[["marker_score_T cells", "marker_score_B cells", "marker_score_Myeloids", "marker_score_Stromal cells", "marker_score_Mast cells"]].max(axis=1)
    gate = epithelial_score.ge(0.5) & epithelial_score.ge(non_epi + 0.1)
    current = mapped_cells.copy()
    current["analysis_cell_type"] = np.where(gate, "Epithelial cells", "")
    effect = paired_effect(summarize_pseudobulk(current), cohort, "Epithelial cells")
    epi = current[current["qc_pass"] & gate]
    return pd.DataFrame([{
        "cohort": cohort, "method": "canonical_marker_gate",
        "markers": "EPCAM/KRT8/KRT18/KRT19/KRT20/CEACAM5/MUC1",
        "epithelial_score_threshold": 0.5, "non_epithelial_exclusion_margin": 0.1,
        "cebpb_excluded_from_gate": True, "gated_epithelial_cells_total": int(len(epi)),
        "gated_epithelial_cells_tumor": int((epi["Class"] == "Tumor").sum()),
        "gated_epithelial_cells_normal": int((epi["Class"] == "Normal").sum()),
        "matched_patient_n": effect["matched_patient_n"], "matched_patients": effect["matched_patients"],
        "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm": effect["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"],
        "wilcoxon_p_value": effect["wilcoxon_p_value"], "direction": effect["direction"],
    }])


def cross_cohort_replication(validation: pd.DataFrame, cohorts: list[str]) -> pd.DataFrame:
    rows = []
    for cell_type in REFERENCE_LABELS:
        current = validation[validation["cell_type"].eq(cell_type)].set_index("cohort")
        row = {"gene_symbol": TARGET, "cell_type": cell_type}
        deltas = []
        p_values = []
        for cohort in cohorts:
            prefix = cohort.lower()
            if cohort in current.index:
                record = current.loc[cohort]
                row[f"{prefix}_matched_n"] = int(record["matched_patient_n"])
                row[f"{prefix}_median_delta"] = record["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"]
                row[f"{prefix}_p"] = record["wilcoxon_p_value"]
                row[f"{prefix}_direction"] = record["direction"]
                deltas.append(record["median_delta_tumor_minus_normal_log1p_pseudobulk_cpm"])
                p_values.append(record["wilcoxon_p_value"])
            else:
                row[f"{prefix}_matched_n"] = 0
                row[f"{prefix}_median_delta"] = np.nan
                row[f"{prefix}_p"] = np.nan
                row[f"{prefix}_direction"] = "Insufficient"
        finite = [d for d in deltas if np.isfinite(d)]
        agree = bool(finite and all((d > 0) == (finite[0] > 0) for d in finite))
        row["same_direction"] = agree
        row["replication_status"] = (
            "Direction and nominal evidence concordant" if agree and len(p_values) == len(cohorts) and all(np.isfinite(p) and p < 0.05 for p in p_values)
            else "Direction concordant; statistical replication not established" if agree
            else "Not concordant/insufficient"
        )
        rows.append(row)
    return pd.DataFrame(rows)


def plot_epithelial_validation(pseudobulk: pd.DataFrame) -> None:
    current = pseudobulk[pseudobulk["passes_min_cells"] & pseudobulk["cell_type"].eq("Epithelial cells") & pseudobulk["Class"].isin(["Tumor", "Normal"])]
    cohorts = ["GSE200997", "GSE132465", "GSE294300"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True, constrained_layout=True)
    for ax, cohort in zip(axes, cohorts):
        data = current[current["cohort"].eq(cohort)]
        tumor = data[data["Class"].eq("Tumor")].set_index("Patient")
        normal = data[data["Class"].eq("Normal")].set_index("Patient")
        patients = sorted(set(tumor.index) & set(normal.index))
        for patient in patients:
            ax.plot([0, 1], [normal.loc[patient, "log1p_pseudobulk_cpm"], tumor.loc[patient, "log1p_pseudobulk_cpm"]], color="#6b7280", alpha=0.4, linewidth=0.8)
        if patients:
            ax.scatter(np.zeros(len(patients)), normal.loc[patients, "log1p_pseudobulk_cpm"], color="#2563eb", s=22, label="Normal")
            ax.scatter(np.ones(len(patients)), tumor.loc[patients, "log1p_pseudobulk_cpm"], color="#dc2626", s=22, label="Tumor")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Normal", "Tumor"]); ax.set_title(f"{cohort}\nmatched n={len(patients)}"); ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("CEBPB log1p pseudobulk CPM")
    axes[-1].legend(frameon=False)
    fig.suptitle("CEBPB in epithelial cells: three-cohort clean rerun")
    fig.savefig(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_epithelial_tumor_normal.png", dpi=240)
    plt.close(fig)


def _fmt(value: object) -> str:
    return "NA" if value is None or pd.isna(value) else f"{float(value):.4g}"


def write_report(cells: pd.DataFrame, validation: pd.DataFrame, primary: pd.DataFrame, replication: pd.DataFrame, count_audit: pd.DataFrame, inclusion: dict[str, pd.DataFrame], status: pd.DataFrame, threshold: pd.DataFrame, manifest: dict) -> None:
    report = [
        "# DINP–CRC clean CEBPB single-cell rerun with GSE294300",
        "",
        "## Frozen design",
        "",
        "CEBPB-only clean rerun using GSE200997, GSE132465, and GSE294300. GSE188711 was removed because it contained only three patients. No Left/Right/sidedness metadata, CD36, or old multicohort outputs were used.",
        "",
        "- GSE132465 uses source curated broad cell-type labels and raw UMI counts.",
        "- GSE200997 and GSE294300 use raw UMI counts and transparent nearest-centroid mapping to the GSE132465 curated labels; CEBPB is excluded from annotation features.",
        "- Primary test: matched patient-level epithelial pseudobulk, Tumor versus Normal, paired Wilcoxon, minimum 50 confident epithelial cells per sample.",
        "- Primary epithelial BH-FDR is recalculated across all three cohort tests; secondary cell-type FDR is calculated within cohort across six prespecified compartments.",
        "",
        "## GSE294300 source and audit",
        "",
        f"- Loaded cells: {int((cells['cohort'] == 'GSE294300').sum()):,}; QC-passing: {int(((cells['cohort'] == 'GSE294300') & cells['qc_pass']).sum()):,}.",
        "- GSE294300 is a paired 18-patient, 36-sample tumor/adjacent-normal 10x study; the deposited cell-batch barcodes were matched to raw 10x barcode cores after terminal suffix normalization, and raw matrices were restricted to that deposited cell subset. The official cell-batch URL and per-sample source URL template are recorded in the manifest.",
        "- Cell types in GSE294300 are reference-mapped, not deposited curated labels; the mapping audit and condition-stratified label composition are exported separately.",
        "",
        "## Annotation-status audit",
        "",
    ]
    for _, row in status.iterrows():
        report.append(f"- {row['cohort']} {row['Class']}: QC={int(row['qc_pass_cells'])}; Confident={row['confident_pct']:.2f}%; Low-confidence={row['low_confidence_pct']:.2f}%; Ambiguous={row['ambiguous_pct']:.2f}%; confident epithelial={row['epithelial_confident_pct']:.2f}%.")
    report.extend(["", "## Patient-level epithelial inclusion", ""])
    for cohort, table in inclusion.items():
        included = table[table["paired included/excluded"].eq("included")]
        report.append(f"- {cohort}: {len(included)} matched patient pairs pass the 50-cell rule.")
    report.extend(["", "## Threshold sensitivity", ""])
    for _, row in threshold.iterrows():
        report.append(f"- {row['cohort']} {row['threshold_variant']}: matched n={int(row['matched_patient_n'])}; median Δ={_fmt(row['median_delta_tumor_minus_normal_log1p_pseudobulk_cpm'])}; P={_fmt(row['wilcoxon_p_value'])}; {row['direction']}.")
    report.extend(["", "## Primary epithelial validation", ""])
    for _, row in primary.iterrows():
        report.append(f"- {row['cohort']}: matched n={int(row['matched_patient_n'])}; median Δ={_fmt(row['median_delta_tumor_minus_normal_log1p_pseudobulk_cpm'])}; P={_fmt(row['wilcoxon_p_value'])}; BH-FDR={_fmt(row['fdr_bh_primary_epithelial_across_cohorts'])}; {row['direction']}.")
    epi_rep = replication[replication["cell_type"].eq("Epithelial cells")]
    if not epi_rep.empty:
        report.extend(["", "## Cross-cohort epithelial audit", "", f"- Direction concordance across the three cohorts: {bool(epi_rep.iloc[0]['same_direction'])}; status: {epi_rep.iloc[0]['replication_status']}."])
    report.extend([
        "",
        "## Interpretation",
        "",
        "Adding GSE294300 increases the paired cohort evidence base without changing the frozen annotation or 50-cell rules. The result should be interpreted as cross-cohort epithelial replication only if direction is concordant and the patient-level evidence supports it after the three-cohort primary BH correction.",
        "",
        "## Key files",
        "",
        "- `DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv`: primary three-cohort epithelial validation.",
        "- `DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv`: patient/sample/cell-type pseudobulk.",
        "- `DINP_CRC_clean_GSE294300_annotation_audit.csv`: all GSE294300 mapped cells and confidence metrics.",
        "- `DINP_CRC_clean_GSE294300_epithelial_inclusion_audit.csv`: explicit paired inclusion/exclusion audit.",
        "- `DINP_CRC_clean_GSE294300_annotation_status_by_condition.csv`: Tumor/Normal annotation-status audit.",
        "- `DINP_CRC_clean_CEBPB_manifest.json`: source and parameter manifest.",
    ])
    (OUTPUT_DIR / "DINP_CRC_clean_CEBPB_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    gse132465, m132465, c132465, u132465 = base.load_gse132465()
    model, _, _, _ = base.fit_reference_mapping(gse132465, c132465, u132465)
    gse200997, m200997, a200997 = base.load_gse200997(model)
    gse294300, m294300, a294300 = load_gse294300(model)
    cohorts = ["GSE200997", "GSE132465", "GSE294300"]
    cells = pd.concat([gse132465, gse200997, gse294300], ignore_index=True, sort=False).fillna("")
    count_audit = patient_epithelial_cell_count_audit(cells)
    inclusion = {cohort: epithelial_inclusion_audit(count_audit, cohort) for cohort in cohorts}
    mapped_cohorts = ["GSE200997", "GSE294300"]
    status = annotation_status_summary(cells, mapped_cohorts)
    labels = label_composition(cells, mapped_cohorts)
    thresholds = pd.concat([threshold_sensitivity_analysis(gse200997, model, "GSE200997"), threshold_sensitivity_analysis(gse294300, model, "GSE294300")], ignore_index=True)
    marker_sensitivity = pd.concat([marker_gate_sensitivity(gse200997, "GSE200997"), marker_gate_sensitivity(gse294300, "GSE294300")], ignore_index=True)
    localization = summarize_localization(cells)
    pseudobulk = summarize_pseudobulk(cells)
    validation = paired_validation(pseudobulk, cohorts)
    primary = validation[validation["primary_epithelial_hypothesis"]].copy()
    replication = cross_cohort_replication(validation, cohorts)

    cell_columns = ["cohort", "cell_id", "Patient", "Sample", "Class", "cell_type", "analysis_cell_type", "cell_subtype", "cell_type_source", "annotation_status", "total_umi", "n_genes", "qc_pass", f"{TARGET}_umi", f"{TARGET}_detected", f"{TARGET}_cpm", f"{TARGET}_log1p_cpm"]
    extra = [c for c in cells.columns if c.startswith("marker_score_") or c in {"top1_reference_label", "top2_reference_label", "top1_reference_similarity", "top2_reference_similarity", "reference_score_margin", "marker_panel_count_gt_0_25"}]
    cells[cell_columns + extra].to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_cell_level.csv", index=False)
    a200997.to_csv(OUTPUT_DIR / "DINP_CRC_clean_GSE200997_annotation_audit.csv", index=False)
    a294300.to_csv(OUTPUT_DIR / "DINP_CRC_clean_GSE294300_annotation_audit.csv", index=False)
    count_audit.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_patient_epithelial_cell_count_audit.csv", index=False)
    for cohort, table in inclusion.items():
        table.to_csv(OUTPUT_DIR / f"DINP_CRC_clean_{cohort}_epithelial_inclusion_audit.csv", index=False)
    status.to_csv(OUTPUT_DIR / "DINP_CRC_clean_reference_mapped_annotation_status_by_condition.csv", index=False)
    labels.to_csv(OUTPUT_DIR / "DINP_CRC_clean_reference_mapped_label_composition_by_condition.csv", index=False)
    thresholds.to_csv(OUTPUT_DIR / "DINP_CRC_clean_reference_mapped_epithelial_threshold_sensitivity.csv", index=False)
    marker_sensitivity.to_csv(OUTPUT_DIR / "DINP_CRC_clean_reference_mapped_marker_gate_epithelial_sensitivity.csv", index=False)
    localization.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_localization.csv", index=False)
    pseudobulk.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv", index=False)
    validation.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_tumor_normal_validation_all_celltypes.csv", index=False)
    primary.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv", index=False)
    replication.to_csv(OUTPUT_DIR / "DINP_CRC_clean_CEBPB_cross_cohort_replication.csv", index=False)
    plot_epithelial_validation(pseudobulk)

    manifest = {
        "analysis": "clean_CEBPB_only_GSE200997_GSE132465_GSE294300",
        "targets": [TARGET],
        "included_cohorts": cohorts,
        "excluded_cohorts": ["GSE188711"],
        "excluded_metadata_fields": ["Location", "Side", "left_right", "tumor_region"],
        "qc_threshold": {"min_umi": base.MIN_QC_UMI, "min_genes": base.MIN_QC_GENES},
        "pseudobulk": {"normalization": "aggregate CEBPB UMI / aggregate library UMI × 1e6 = pseudobulk CPM", "grouping": "cohort × patient × sample × class × cell_type", "minimum_cells_per_sample_celltype": MIN_CELLS_PER_SAMPLE_CELLTYPE, "tumor_normal_test": "paired Wilcoxon across matched patients", "secondary_fdr_family": "six prespecified cell types within each cohort", "primary_fdr_family": "epithelial tests across all three cohorts"},
        "reference_mapping": {"reference_cohort": "GSE132465", "labels": REFERENCE_LABELS, "features": model["features"], "similarity": "cosine similarity after gene-wise standardization fit on GSE132465 marker expression", "confidence_rule": "top1 similarity >= reference 5th percentile and top1-top2 margin >= reference 5th percentile", "calibration_estimator": "exact leave-one-out median centroid", "sim_threshold_5pct": model["sim_threshold_5pct_reference"], "margin_threshold_5pct": model["margin_threshold_5pct_reference"], "sim_threshold_10pct": model["sim_threshold_10pct_reference"], "margin_threshold_10pct": model["margin_threshold_10pct_reference"], "reference_loo_top1_median": model["reference_loo_top1_median"], "reference_loo_margin_median": model["reference_loo_margin_median"], "reference_loo_label_accuracy": model["reference_loo_label_accuracy"]},
        "cohorts": {"GSE132465": m132465, "GSE200997": m200997, "GSE294300": m294300},
        "cell_counts": {"combined_loaded": int(len(cells)), "combined_qc_pass": int(cells["qc_pass"].sum()), "per_cohort_loaded": cells.groupby("cohort").size().to_dict(), "per_cohort_qc_pass": cells.groupby("cohort")["qc_pass"].sum().astype(int).to_dict(), "primary_included_patients": {cohort: inclusion[cohort].loc[inclusion[cohort]["paired included/excluded"].eq("included"), "Patient"].tolist() for cohort in cohorts}},
        "source_note": "GSE188711 was excluded; no old CEBPB/CD36 or sidedness output was read as input. GSE294300 raw counts were read from official NCBI GEO per-sample supplements.",
    }
    (OUTPUT_DIR / "DINP_CRC_clean_CEBPB_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    write_report(cells, validation, primary, replication, count_audit, inclusion, status, thresholds, manifest)
    print(json.dumps({"output_dir": str(OUTPUT_DIR), "cells": len(cells), "qc_pass": int(cells["qc_pass"].sum()), "primary_epithelial": primary[["cohort", "matched_patient_n", "median_delta_tumor_minus_normal_log1p_pseudobulk_cpm", "wilcoxon_p_value", "fdr_bh_primary_epithelial_across_cohorts", "direction"]].to_dict("records")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
