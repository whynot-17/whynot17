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
DATA_DIR = ROOT / "work" / "scRNA_GSE132465"
OUTPUT_DIR = ROOT / "outputs" / "DINP_CRC_scRNA_GSE132465_CEBPB_CD36"
ANNOTATION_PATH = DATA_DIR / "cell_annotation.txt.gz"
COUNTS_PATH = DATA_DIR / "raw_UMI_count_matrix.txt.gz"
TARGETS = ["CEBPB", "CD36"]
MIN_CELLS_PER_SAMPLE_CELLTYPE = 50
DATASET = "GSE132465"
ANNOTATION_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
COUNTS_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_annotation() -> pd.DataFrame:
    ann = pd.read_csv(ANNOTATION_PATH, sep="\t", dtype=str).fillna("")
    required = {"Index", "Patient", "Class", "Sample", "Cell_type", "Cell_subtype"}
    missing = required - set(ann.columns)
    if missing:
        raise AssertionError(f"Annotation columns missing: {sorted(missing)}")
    ann = ann.rename(columns={"Index": "cell_id"})
    ann["cell_id"] = ann["cell_id"].str.strip()
    if ann["cell_id"].duplicated().any():
        raise AssertionError("Annotation contains duplicated cell IDs")
    return ann


def extract_target_rows(cell_ids: list[str]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    found: dict[str, np.ndarray] = {}
    total_umi = np.zeros(len(cell_ids), dtype=np.int64)
    with gzip.open(COUNTS_PATH, "rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(b"\t")
        matrix_cell_ids = [value.decode("utf-8") for value in header[1:]]
        if matrix_cell_ids != cell_ids:
            if set(matrix_cell_ids) != set(cell_ids):
                raise AssertionError("Annotation cell IDs do not match raw matrix columns")
            raise AssertionError("Raw matrix and annotation cell order differ")
        for raw_line in handle:
            line = raw_line.rstrip(b"\r\n")
            gene_bytes, payload = line.split(b"\t", 1)
            values = np.fromstring(payload, dtype=np.int64, sep="\t")
            if len(values) != len(cell_ids):
                raise AssertionError(f"Unexpected column count for {gene_bytes.decode('utf-8')}")
            total_umi += values
            gene = gene_bytes.decode("utf-8").strip().upper()
            if gene in TARGETS:
                found[gene] = values.astype(np.int32)
    missing = set(TARGETS) - set(found)
    if missing:
        raise AssertionError(f"Targets absent from count matrix: {sorted(missing)}")
    return found, total_umi


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    output = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return pd.Series(output, index=p_values.index)
    valid_indices = np.flatnonzero(valid)
    order = valid_indices[np.argsort(values[valid_indices])]
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output[order] = np.minimum(adjusted, 1.0)
    return pd.Series(output, index=p_values.index)


def make_cell_level(ann: pd.DataFrame, counts: dict[str, np.ndarray]) -> pd.DataFrame:
    output = ann[["cell_id", "Patient", "Class", "Sample", "Cell_type", "Cell_subtype"]].copy()
    output["total_umi"] = ann.attrs["total_umi"]
    for gene in TARGETS:
        output[f"{gene}_umi"] = counts[gene]
        output[f"{gene}_detected"] = counts[gene] > 0
        output[f"{gene}_log1p_umi"] = np.log1p(counts[gene].astype(float))
        output[f"{gene}_cpm"] = np.divide(
            counts[gene].astype(float) * 1_000_000.0,
            output["total_umi"].to_numpy(dtype=float),
            out=np.zeros(len(output), dtype=float),
            where=output["total_umi"].to_numpy(dtype=float) > 0,
        )
        output[f"{gene}_log1p_cpm"] = np.log1p(output[f"{gene}_cpm"])
    return output


def summarize_localization(cells: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for gene in TARGETS:
        umi_col = f"{gene}_umi"
        det_col = f"{gene}_detected"
        log_col = f"{gene}_log1p_umi"
        cpm_col = f"{gene}_cpm"
        log_cpm_col = f"{gene}_log1p_cpm"
        grouped = cells.groupby(group_cols, dropna=False, sort=True)
        summary = grouped.agg(
            n_cells=("cell_id", "size"),
            total_library_umi=("total_umi", "sum"),
            total_umi=(umi_col, "sum"),
            detected_cells=(det_col, "sum"),
            mean_umi_per_cell=(umi_col, "mean"),
            median_umi_per_cell=(umi_col, "median"),
            mean_log1p_umi=(log_col, "mean"),
            mean_cpm_per_cell=(cpm_col, "mean"),
            median_cpm_per_cell=(cpm_col, "median"),
            mean_log1p_cpm=(log_cpm_col, "mean"),
        ).reset_index()
        summary["gene_symbol"] = gene
        summary["detection_rate"] = summary["detected_cells"] / summary["n_cells"]
        summary["pseudobulk_cpm"] = summary["total_umi"] * 1_000_000.0 / summary["total_library_umi"].replace(0, np.nan)
        summary["log1p_pseudobulk_cpm"] = np.log1p(summary["pseudobulk_cpm"])
        class_cols = [column for column in ["Class"] if column in group_cols]
        if class_cols:
            class_total = summary.groupby(class_cols, dropna=False)["total_umi"].transform("sum")
            summary["target_umi_fraction_within_class"] = summary["total_umi"] / class_total.replace(0, np.nan)
        else:
            summary["target_umi_fraction_within_class"] = summary["total_umi"] / summary["total_umi"].sum()
        rows.append(summary)
    output = pd.concat(rows, ignore_index=True)
    return output.sort_values([*group_cols, "gene_symbol"]).reset_index(drop=True)


def make_sample_pseudobulk(cells: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["Patient", "Class", "Sample", "Cell_type"]
    rows = []
    for gene in TARGETS:
        umi_col = f"{gene}_umi"
        det_col = f"{gene}_detected"
        grouped = cells.groupby(group_cols, dropna=False, sort=True)
        summary = grouped.agg(
            n_cells=("cell_id", "size"),
            total_library_umi=("total_umi", "sum"),
            total_umi=(umi_col, "sum"),
            detected_cells=(det_col, "sum"),
        ).reset_index()
        summary["gene_symbol"] = gene
        summary["mean_umi_per_cell"] = summary["total_umi"] / summary["n_cells"]
        summary["log1p_mean_umi_per_cell"] = np.log1p(summary["mean_umi_per_cell"])
        summary["pseudobulk_cpm"] = summary["total_umi"] * 1_000_000.0 / summary["total_library_umi"].replace(0, np.nan)
        summary["log1p_pseudobulk_cpm"] = np.log1p(summary["pseudobulk_cpm"])
        summary["detection_rate"] = summary["detected_cells"] / summary["n_cells"]
        summary["passes_min_cells"] = summary["n_cells"].ge(MIN_CELLS_PER_SAMPLE_CELLTYPE)
        rows.append(summary)
    return pd.concat(rows, ignore_index=True).sort_values(group_cols + ["gene_symbol"]).reset_index(drop=True)


def make_paired_validation(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cell_types = sorted(pseudobulk["Cell_type"].unique())
    for gene in TARGETS:
        for cell_type in cell_types:
            subset = pseudobulk[
                pseudobulk["gene_symbol"].eq(gene)
                & pseudobulk["Cell_type"].eq(cell_type)
                & pseudobulk["passes_min_cells"]
            ].copy()
            tumor = subset[subset["Class"].eq("Tumor")].set_index("Patient")
            normal = subset[subset["Class"].eq("Normal")].set_index("Patient")
            matched_patients = sorted(set(tumor.index) & set(normal.index))
            tumor_values = tumor.loc[matched_patients, "log1p_pseudobulk_cpm"].to_numpy(dtype=float) if matched_patients else np.array([])
            normal_values = normal.loc[matched_patients, "log1p_pseudobulk_cpm"].to_numpy(dtype=float) if matched_patients else np.array([])
            differences = tumor_values - normal_values
            if len(differences) < 3:
                statistic, p_value = np.nan, np.nan
            elif np.allclose(differences, 0):
                statistic, p_value = 0.0, 1.0
            else:
                test = wilcoxon(tumor_values, normal_values, zero_method="wilcox", alternative="two-sided", method="auto")
                statistic, p_value = float(test.statistic), float(test.pvalue)
            rows.append(
                {
                    "gene_symbol": gene,
                    "Cell_type": cell_type,
                    "matched_patient_n": len(matched_patients),
                    "matched_patients": ";".join(matched_patients),
                    "median_tumor_log1p_pseudobulk_cpm": float(np.median(tumor_values)) if len(tumor_values) else np.nan,
                    "median_normal_log1p_pseudobulk_cpm": float(np.median(normal_values)) if len(normal_values) else np.nan,
                    "median_paired_delta_tumor_minus_normal": float(np.median(differences)) if len(differences) else np.nan,
                    "tumor_gt_normal_n": int(np.sum(differences > 0)) if len(differences) else 0,
                    "tumor_eq_normal_n": int(np.sum(differences == 0)) if len(differences) else 0,
                    "tumor_lt_normal_n": int(np.sum(differences < 0)) if len(differences) else 0,
                    "wilcoxon_statistic": statistic,
                    "wilcoxon_p_value": p_value,
                    "passes_min_cells_filter": True,
                }
            )
    output = pd.DataFrame(rows)
    output["wilcoxon_fdr_bh"] = bh_adjust(output["wilcoxon_p_value"])
    output["direction"] = np.where(
        output["median_paired_delta_tumor_minus_normal"].gt(0),
        "Tumor higher",
        np.where(output["median_paired_delta_tumor_minus_normal"].lt(0), "Normal higher", "No difference"),
    )
    return output


def plot_localization(localization: pd.DataFrame) -> None:
    cell_types = sorted(localization["Cell_type"].unique())
    class_order = ["Tumor", "Normal"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True, constrained_layout=True)
    for ax, gene in zip(axes, TARGETS):
        subset = localization[localization["gene_symbol"].eq(gene)].copy()
        subset["x"] = subset["Cell_type"].map({value: index for index, value in enumerate(cell_types)})
        subset["y"] = subset["Class"].map({value: index for index, value in enumerate(class_order)})
        size = 25 + 450 * subset["detection_rate"].fillna(0)
        scatter = ax.scatter(
            subset["x"],
            subset["y"],
            s=size,
            c=subset["mean_log1p_umi"],
            cmap="viridis",
            edgecolors="black",
            linewidths=0.35,
        )
        ax.set_title(gene)
        ax.set_xticks(range(len(cell_types)))
        ax.set_xticklabels(cell_types, rotation=35, ha="right")
        ax.set_yticks(range(len(class_order)))
        ax.set_yticklabels(class_order)
        ax.set_xlabel("Curated cell type")
        ax.grid(axis="x", alpha=0.15)
        fig.colorbar(scatter, ax=ax, label="Mean log1p CPM")
    axes[0].set_ylabel("Sample class")
    fig.suptitle("GSE132465 CRC single-cell target localization\nDot size = detection rate")
    fig.savefig(OUTPUT_DIR / "scRNA_target_localization_dotplot.png", dpi=220)
    plt.close(fig)


def plot_paired_validation(pseudobulk: pd.DataFrame, paired: pd.DataFrame) -> None:
    cell_types = sorted(pseudobulk["Cell_type"].unique())
    fig, axes = plt.subplots(2, len(cell_types), figsize=(3.0 * len(cell_types), 7), squeeze=False, constrained_layout=True)
    for row_index, gene in enumerate(TARGETS):
        for col_index, cell_type in enumerate(cell_types):
            ax = axes[row_index, col_index]
            subset = pseudobulk[
                pseudobulk["gene_symbol"].eq(gene)
                & pseudobulk["Cell_type"].eq(cell_type)
                & pseudobulk["passes_min_cells"]
            ]
            tumor = subset[subset["Class"].eq("Tumor")].set_index("Patient")
            normal = subset[subset["Class"].eq("Normal")].set_index("Patient")
            matched = sorted(set(tumor.index) & set(normal.index))
            for patient in matched:
                ax.plot(
                    [0, 1],
                    [normal.loc[patient, "log1p_pseudobulk_cpm"], tumor.loc[patient, "log1p_pseudobulk_cpm"]],
                    color="#777777",
                    alpha=0.45,
                    linewidth=0.8,
                )
            if matched:
                ax.scatter(np.zeros(len(matched)), normal.loc[matched, "log1p_pseudobulk_cpm"], color="#3b82f6", s=16, zorder=3, label="Normal")
                ax.scatter(np.ones(len(matched)), tumor.loc[matched, "log1p_pseudobulk_cpm"], color="#dc2626", s=16, zorder=3, label="Tumor")
            result = paired[(paired["gene_symbol"].eq(gene)) & (paired["Cell_type"].eq(cell_type))].iloc[0]
            p_text = "NA" if pd.isna(result["wilcoxon_p_value"]) else f"p={result['wilcoxon_p_value']:.3g}"
            ax.set_title(f"{cell_type}\n{p_text}, n={int(result['matched_patient_n'])}", fontsize=9)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Normal", "Tumor"], rotation=45, ha="right")
            ax.grid(axis="y", alpha=0.2)
            if col_index == 0:
                ax.set_ylabel(f"{gene}\nlog1p(pseudobulk CPM)")
    fig.suptitle("Patient-level paired pseudobulk validation (SMC01–SMC10 where eligible)")
    fig.savefig(OUTPUT_DIR / "scRNA_patient_pseudobulk_paired_validation.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if not ANNOTATION_PATH.exists() or not COUNTS_PATH.exists():
        raise FileNotFoundError("Download GSE132465 annotation and raw UMI files before running")
    ann = load_annotation()
    cell_ids = ann["cell_id"].tolist()
    counts, total_umi = extract_target_rows(cell_ids)
    ann.attrs["total_umi"] = total_umi
    cells = make_cell_level(ann, counts)
    localization = summarize_localization(cells, ["Class", "Cell_type"])
    subtype_localization = summarize_localization(cells, ["Class", "Cell_type", "Cell_subtype"])
    pseudobulk = make_sample_pseudobulk(cells)
    paired = make_paired_validation(pseudobulk)

    cells.to_csv(OUTPUT_DIR / "GSE132465_CEBPB_CD36_cell_level_expression.csv", index=False)
    localization.to_csv(OUTPUT_DIR / "GSE132465_CEBPB_CD36_cell_type_localization.csv", index=False)
    subtype_localization.to_csv(OUTPUT_DIR / "GSE132465_CEBPB_CD36_cell_subtype_localization.csv", index=False)
    pseudobulk.to_csv(OUTPUT_DIR / "GSE132465_CEBPB_CD36_patient_celltype_pseudobulk.csv", index=False)
    paired.to_csv(OUTPUT_DIR / "GSE132465_CEBPB_CD36_paired_pseudobulk_validation.csv", index=False)
    plot_localization(localization)
    plot_paired_validation(pseudobulk, paired)

    tumor_cells = int((ann["Class"] == "Tumor").sum())
    normal_cells = int((ann["Class"] == "Normal").sum())
    matched_patients = sorted(set(ann.loc[ann["Class"].eq("Tumor"), "Patient"]) & set(ann.loc[ann["Class"].eq("Normal"), "Patient"]))
    top_rows = []
    tumor_loc = localization[localization["Class"].eq("Tumor")].copy()
    for gene in TARGETS:
        gene_rows = tumor_loc[tumor_loc["gene_symbol"].eq(gene)].sort_values(["detection_rate", "mean_log1p_umi", "total_umi"], ascending=False)
        top_rows.append({
            "gene_symbol": gene,
            "top_tumor_cell_type_by_detection": str(gene_rows.iloc[0]["Cell_type"]),
            "top_tumor_detection_rate": float(gene_rows.iloc[0]["detection_rate"]),
            "top_tumor_cell_type_by_mean_log1p_cpm": str(gene_rows.sort_values("mean_log1p_cpm", ascending=False).iloc[0]["Cell_type"]),
            "top_tumor_mean_log1p_cpm": float(gene_rows["mean_log1p_cpm"].max()),
        })
    summary = {
        "dataset": DATASET,
        "annotation_source": ANNOTATION_URL,
        "count_source": COUNTS_URL,
        "annotation_sha256": sha256(ANNOTATION_PATH),
        "count_matrix_sha256": sha256(COUNTS_PATH),
        "cell_n": int(len(ann)),
        "tumor_cell_n": tumor_cells,
        "normal_cell_n": normal_cells,
        "patient_n": int(ann["Patient"].nunique()),
        "tumor_patient_n": int(ann.loc[ann["Class"].eq("Tumor"), "Patient"].nunique()),
        "normal_patient_n": int(ann.loc[ann["Class"].eq("Normal"), "Patient"].nunique()),
        "matched_patient_n": len(matched_patients),
        "matched_patients": matched_patients,
        "cell_types": sorted(ann["Cell_type"].unique()),
        "target_genes": TARGETS,
        "target_rows_found": sorted(counts),
        "minimum_cells_per_sample_celltype_for_pseudobulk": MIN_CELLS_PER_SAMPLE_CELLTYPE,
        "pseudobulk_unit": "sum raw target UMI per patient x class x cell type, normalized by aggregate library UMI to pseudobulk CPM",
        "paired_test": "two-sided Wilcoxon signed-rank test on log1p(pseudobulk CPM) across matched patients; BH correction across 12 gene x cell-type tests",
        "top_tumor_localization": top_rows,
        "qc": {
            "annotation_cell_ids_unique": bool(ann["cell_id"].is_unique),
            "cell_level_rows": int(len(cells)),
            "pseudobulk_rows": int(len(pseudobulk)),
            "paired_rows": int(len(paired)),
        "total_library_umi_sum": int(cells["total_umi"].sum()),
        "targets_nonzero_umi": {gene: int(cells[f"{gene}_umi"].sum()) for gene in TARGETS},
            "targets_detected_cells": {gene: int(cells[f"{gene}_detected"].sum()) for gene in TARGETS},
        },
    }
    (OUTPUT_DIR / "GSE132465_CEBPB_CD36_analysis_manifest.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    report = [
        "# DINP–CRC CEBPB/CD36 single-cell localization and patient-level pseudobulk validation",
        "",
        f"Dataset: **{DATASET}**. The public GEO series contains {len(ann):,} annotated cells from {ann['Patient'].nunique()} patients, including {tumor_cells:,} tumor cells from {ann.loc[ann['Class'].eq('Tumor'), 'Patient'].nunique()} patients and {normal_cells:,} normal-mucosa cells from {ann.loc[ann['Class'].eq('Normal'), 'Patient'].nunique()} patients. Matched tumor/normal patients: {len(matched_patients)} ({', '.join(matched_patients)}).",
        "",
        "## Scope and design",
        "",
        "- Localization uses the GEO-provided curated cell-type and cell-subtype annotations; no cell type was inferred from target expression alone.",
        "- CEBPB and CD36 raw UMI rows were extracted from the public matrix after exact cell-ID concordance checking.",
        f"- Pseudobulk aggregation is at patient × class × cell type, with a minimum of {MIN_CELLS_PER_SAMPLE_CELLTYPE} cells per sample × cell type. Raw target UMI sums are normalized by aggregate library UMI to pseudobulk CPM; mean UMI per cell and detection rate are also retained.",
        "- Validation uses only matched patients and a two-sided Wilcoxon signed-rank test on log1p(pseudobulk CPM), with Benjamini–Hochberg correction across the 12 gene × cell-type tests.",
        "- This is an expression-localization/validation analysis; it does not establish direct DINP binding or causal CRC activity.",
        "",
        "## Key localization result",
        "",
        "| Gene | Highest tumor detection-rate cell type | Detection rate | Highest tumor mean-expression cell type | Mean log1p CPM |",
        "|---|---|---:|---|---:|",
    ]
    for item in top_rows:
        report.append(f"| {item['gene_symbol']} | {item['top_tumor_cell_type_by_detection']} | {item['top_tumor_detection_rate']:.4f} | {item['top_tumor_cell_type_by_mean_log1p_cpm']} | {item['top_tumor_mean_log1p_cpm']:.4f} |")
    report.extend([
        "",
        "## Patient-level paired validation",
        "",
        "See `GSE132465_CEBPB_CD36_paired_pseudobulk_validation.csv` for every gene × cell-type test and `scRNA_patient_pseudobulk_paired_validation.png` for the matched-patient trajectories.",
        "",
        "| Gene | Cell type | Matched patients | Median tumor − normal log1p pseudobulk CPM | Wilcoxon P | BH FDR | Direction |",
        "|---|---|---:|---:|---:|---:|---|",
    ])
    for _, row in paired.sort_values(["gene_symbol", "wilcoxon_fdr_bh", "Cell_type"], na_position="last").iterrows():
        p = "NA" if pd.isna(row["wilcoxon_p_value"]) else f"{row['wilcoxon_p_value']:.4g}"
        fdr = "NA" if pd.isna(row["wilcoxon_fdr_bh"]) else f"{row['wilcoxon_fdr_bh']:.4g}"
        delta = "NA" if pd.isna(row["median_paired_delta_tumor_minus_normal"]) else f"{row['median_paired_delta_tumor_minus_normal']:.4f}"
        report.append(f"| {row['gene_symbol']} | {row['Cell_type']} | {int(row['matched_patient_n'])} | {delta} | {p} | {fdr} | {row['direction']} |")
    report.extend([
        "",
        "## Provenance and QC",
        "",
        f"- Annotation download: `{ANNOTATION_URL}`",
        f"- Raw UMI download: `{COUNTS_URL}`",
        f"- Target rows found: {', '.join(sorted(counts))}; cell-level rows: {len(cells):,}; patient/cell-type pseudobulk rows: {len(pseudobulk):,}.",
        f"- SHA256 values and machine-readable parameters are in `GSE132465_CEBPB_CD36_analysis_manifest.json`.",
        "- The large public raw matrix is not copied into the repository; the download URLs and SHA256 hashes make the analysis auditable.",
    ])
    (OUTPUT_DIR / "GSE132465_CEBPB_CD36_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"dataset={DATASET} cells={len(ann)} patients={ann['Patient'].nunique()} matched_patients={len(matched_patients)}")
    print(f"targets={','.join(TARGETS)} target_nonzero_umi=" + ",".join(f"{gene}:{int(cells[f'{gene}_umi'].sum())}" for gene in TARGETS))
    print("top tumor localization:")
    print(pd.DataFrame(top_rows).to_string(index=False))
    print("paired validation:")
    print(paired[["gene_symbol", "Cell_type", "matched_patient_n", "median_paired_delta_tumor_minus_normal", "wilcoxon_p_value", "wilcoxon_fdr_bh", "direction"]].to_string(index=False))


if __name__ == "__main__":
    main()
