"""Minimal CRC spatial test: PUFA/AA pressure, neuronal signal and ferroptosis stress.

This is an exploratory, one-section validation using the P2CRC Visium HD sample
from GSE280315 (stage III-B). It intentionally avoids a full spatial model.

Primary outputs:
  * per-bin scores and inferred region labels;
  * region summaries;
  * descriptive Spearman correlations;
  * two requested OLS models with HC3 standard errors;
  * a coordinate map and a reproducibility report.

Important limitation: this is one tissue section. Spot/bin-level p-values are
not patient-level evidence and should not be interpreted as an independent-
sample replication or a stage comparison.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import statsmodels.formula.api as smf
from scipy.ndimage import binary_dilation
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "work" / "data" / "gse280315"
DEFAULT_OUT = ROOT / "outputs"
H5_PATH = DATA_DIR / "GSM8594568_P2CRC_filtered_feature_bc_matrix.h5"
META_PATH = DATA_DIR / "GSM8594568_P2CRC_Metadata.parquet.gz"
POSITION_PATH = DATA_DIR / "GSM8594568_P2CRC_tissue_positions.parquet.gz"
ANNOTATION_PATH = DATA_DIR / "8um_squares_annotation.csv"

GENE_SETS = {
    "pufa_incorporation": ["ACSL4", "LPCAT3", "AGPAT3"],
    "aa_routing": ["PLA2G4A", "PTGS2", "PTGES", "ALOX5", "ALOX15"],
    "neuronal": ["TUBB3", "ELAVL3", "ELAVL4", "SNAP25", "UCHL1", "SYP", "PHOX2B"],
    "ferroptosis_stress": ["TFRC", "ACSL4", "LPCAT3", "ALOX15", "POR", "CYB5R1"],
    "ferroptosis_defense": ["SLC7A11", "GPX4", "AIFM2", "GCH1"],
}
FERROPTOSIS_STRESS_NONOVERLAP = ["TFRC", "POR", "CYB5R1"]


def decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def read_gz_parquet(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    with gzip.open(path, "rb") as handle:
        table = pq.read_table(pa.BufferReader(handle.read()), columns=columns)
    return table.to_pandas()


def read_targeted_h5(path: Path, gene_sets: dict[str, list[str]], chunk_size: int = 20_000) -> tuple[np.ndarray, np.ndarray, list[str], dict[str, list[str]]]:
    requested = list(dict.fromkeys(gene for genes in gene_sets.values() for gene in genes))
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        features = [decode(x) for x in matrix["features"]["name"][:]]
        feature_lookup: dict[str, int] = {}
        for index, name in enumerate(features):
            feature_lookup.setdefault(name.upper(), index)
        available = {gene: gene in feature_lookup for gene in requested}
        target_indices = [feature_lookup.get(gene.upper(), -1) for gene in requested]
        n_genes, n_spots = [int(x) for x in matrix["shape"][:]]
        del n_genes
        values = np.zeros((n_spots, len(requested)), dtype=np.float32)
        totals = np.zeros(n_spots, dtype=np.float64)
        indptr = matrix["indptr"][:]
        data_ds = matrix["data"]
        indices_ds = matrix["indices"]

        for start in range(0, n_spots, chunk_size):
            stop = min(start + chunk_size, n_spots)
            pointers = indptr[start : stop + 1]
            lo, hi = int(pointers[0]), int(pointers[-1])
            column_counts = np.diff(pointers)
            local_columns = np.repeat(np.arange(stop - start, dtype=np.int64), column_counts)
            chunk_data = np.asarray(data_ds[lo:hi], dtype=np.float32)
            chunk_indices = np.asarray(indices_ds[lo:hi], dtype=np.int64)
            totals[start:stop] = np.bincount(local_columns, weights=chunk_data, minlength=stop - start)
            for gene_position, feature_index in enumerate(target_indices):
                if feature_index < 0:
                    continue
                selected = chunk_indices == feature_index
                if np.any(selected):
                    values[start:stop, gene_position][local_columns[selected]] = chunk_data[selected]
        barcodes = np.asarray([decode(x) for x in matrix["barcodes"][:]], dtype=object)

    present_by_set = {
        name: [gene for gene in genes if available.get(gene, False)]
        for name, genes in gene_sets.items()
    }
    return values, totals, barcodes.tolist(), present_by_set


def zscore(values: np.ndarray) -> np.ndarray:
    mean = np.nanmean(values, axis=0)
    sd = np.nanstd(values, axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return np.clip((values - mean) / sd, -5.0, 5.0)


def score_module(log_cpm: np.ndarray, genes: list[str], requested: list[str]) -> np.ndarray:
    positions = [requested.index(gene) for gene in genes if gene in requested]
    if not positions:
        return np.full(log_cpm.shape[0], np.nan)
    return np.nanmean(zscore(log_cpm[:, positions]), axis=1)


def infer_regions(meta: pd.DataFrame, labels: pd.Series) -> pd.Series:
    """Map pathologist labels and infer a conservative neoplasm edge proxy."""
    rows = meta["array_row"].to_numpy(int)
    cols = meta["array_col"].to_numpy(int)
    shape = (int(rows.max()) + 1, int(cols.max()) + 1)
    label_grid = np.full(shape, "unannotated", dtype=object)
    label_grid[rows, cols] = labels.to_numpy(object)
    neoplasm = label_grid == "Neoplasm"
    non_neoplasm = ~neoplasm & (label_grid != "unannotated")
    edge = binary_dilation(non_neoplasm, structure=np.ones((3, 3), dtype=bool)) & neoplasm

    result = np.full(len(meta), "unannotated", dtype=object)
    raw = labels.to_numpy(str)
    result[raw == "Neoplasm"] = "tumor_core"
    result[(raw == "Neoplasm") & edge[rows, cols]] = "invasive_front_proxy"
    result[raw == "Non-neoplastic Epithelium"] = "adjacent_non_tumor"
    result[raw == "Connective Tissue"] = "connective_tissue"
    result[raw == "Smooth Muscle"] = "muscularis_smooth_muscle"
    result[np.isin(raw, ["Vessel", "Veins"])] = "vascular"
    result[raw == "Outside"] = "outside"
    return pd.Series(result, index=meta.index, name="tumor_region")


def correlation_row(name: str, subset_name: str, x: pd.Series, y: pd.Series, subset: pd.Series) -> dict[str, object]:
    keep = subset & x.notna() & y.notna()
    if keep.sum() < 4 or x[keep].nunique() < 2 or y[keep].nunique() < 2:
        return {"analysis": name, "subset": subset_name, "n": int(keep.sum()), "spearman_rho": np.nan, "p_value": np.nan}
    rho, p_value = spearmanr(x[keep], y[keep])
    return {"analysis": name, "subset": subset_name, "n": int(keep.sum()), "spearman_rho": float(rho), "p_value": float(p_value)}


def run_ols(data: pd.DataFrame, formula: str, model_name: str) -> pd.DataFrame:
    keep = data[[term for term in ["neuronal_score", "pufa_aa_pressure_score", "ferroptosis_stress_score", "ferroptosis_stress_nonoverlap_score", "neural_containing", "tumor_region"] if term in data]].notna().all(axis=1)
    subset = data.loc[keep].copy()
    if subset.shape[0] < 20 or subset["tumor_region"].nunique() < 2:
        return pd.DataFrame([{"model": model_name, "term": "MODEL_NOT_ESTIMABLE", "estimate": np.nan, "std_error": np.nan, "p_value_hc3": np.nan, "n": subset.shape[0], "r_squared": np.nan}])
    fit = smf.ols(formula, data=subset).fit(cov_type="HC3")
    rows = []
    for term, estimate in fit.params.items():
        rows.append({
            "model": model_name,
            "term": term,
            "estimate": float(estimate),
            "std_error": float(fit.bse[term]),
            "p_value_hc3": float(fit.pvalues[term]),
            "n": int(subset.shape[0]),
            "r_squared": float(fit.rsquared),
        })
    return pd.DataFrame(rows)


def make_map(data: pd.DataFrame, path: Path) -> None:
    rng = np.random.default_rng(20260902)
    max_points = 100_000
    if len(data) > max_points:
        selected = rng.choice(len(data), size=max_points, replace=False)
        plot_data = data.iloc[selected]
    else:
        plot_data = data
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    panels = [
        ("pufa_aa_pressure_score", "PUFA/AA pressure"),
        ("neuronal_score", "Neuronal score"),
        ("ferroptosis_stress_score", "Ferroptosis stress"),
    ]
    for axis, (column, title) in zip(axes[0], panels, strict=True):
        points = axis.scatter(plot_data["array_col"], -plot_data["array_row"], c=plot_data[column], s=1, cmap="viridis", alpha=0.65, linewidths=0)
        axis.set_title(title)
        axis.set_aspect("equal")
        axis.set_xlabel("array col")
        axis.set_ylabel("array row")
        fig.colorbar(points, ax=axis, fraction=0.046, pad=0.04)
    region_colors = {
        "tumor_core": "#d73027",
        "invasive_front_proxy": "#fc8d59",
        "adjacent_non_tumor": "#91cf60",
        "connective_tissue": "#4575b4",
        "muscularis_smooth_muscle": "#984ea3",
        "vascular": "#4daf4a",
    }
    for region, color in region_colors.items():
        subset = plot_data[plot_data["tumor_region"] == region]
        axes[1, 0].scatter(subset["array_col"], -subset["array_row"], s=1, color=color, label=region, alpha=0.65, linewidths=0)
    axes[1, 0].set_title("Pathologist region / edge proxy")
    axes[1, 0].set_aspect("equal")
    axes[1, 0].legend(markerscale=6, fontsize=7, loc="best")
    axes[1, 1].scatter(plot_data["pufa_aa_pressure_score"], plot_data["neuronal_score"], c=plot_data["neural_containing"], s=1, alpha=0.25, cmap="coolwarm", linewidths=0)
    axes[1, 1].set(xlabel="PUFA/AA pressure", ylabel="neuronal score", title="Pressure vs neuronal signal")
    axes[1, 2].scatter(plot_data["pufa_aa_pressure_score"], plot_data["ferroptosis_stress_score"], c=plot_data["neural_containing"], s=1, alpha=0.25, cmap="coolwarm", linewidths=0)
    axes[1, 2].set(xlabel="PUFA/AA pressure", ylabel="ferroptosis stress", title="Pressure vs ferroptosis stress")
    fig.suptitle("GSE280315 P2CRC: minimal spatial validation", fontsize=14)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_report(out: Path, audit: dict[str, object], regions: pd.DataFrame, correlations: pd.DataFrame, models: pd.DataFrame, present: dict[str, list[str]]) -> None:
    region_lines = ["| region | bins | PUFA/AA pressure | neuronal | ferroptosis stress | non-overlap stress | defense |", "|---|---:|---:|---:|---:|---:|---:|"]
    for row in regions.itertuples(index=False):
        region_lines.append(f"| {row.tumor_region} | {row.n_bins} | {row.pufa_aa_pressure_mean:.3f} | {row.neuronal_mean:.3f} | {row.ferroptosis_stress_mean:.3f} | {row.ferroptosis_stress_nonoverlap_mean:.3f} | {row.ferroptosis_defense_mean:.3f} |")
    corr_lines = ["| analysis | subset | n | Spearman rho | p |", "|---|---|---:|---:|---:|"]
    for row in correlations.itertuples(index=False):
        corr_lines.append(f"| {row.analysis} | {row.subset} | {row.n} | {row.spearman_rho:.4g} | {row.p_value:.4g} |")
    model_lines = ["| model | term | estimate | HC3 p | n | R² |", "|---|---|---:|---:|---:|---:|"]
    for row in models.itertuples(index=False):
        model_lines.append(f"| {row.model} | {row.term} | {row.estimate:.4g} | {row.p_value_hc3:.4g} | {row.n} | {row.r_squared:.4g} |")
    lines = [
        "# Minimal CRC spatial PUFA/AA–neuronal–ferroptosis validation",
        "",
        "## Frozen scope",
        "",
        "- Dataset: GSE280315 P2CRC Visium HD, one stage III-B CRC section.",
        "- Scores are within-section z-scored means of log1p(CPM+1) expression; they are transcriptional proxies, not lipid concentrations or measured flux.",
        "- Neoplasm bins touching any non-neoplasm pathologist-labeled bin in the 8-neighbor grid are labeled `invasive_front_proxy`; other Neoplasm bins are `tumor_core`.",
        "- `neural_containing` means at least two neuronal genes detected in the bin.",
        "- The fixed ferroptosis-stress score overlaps with PUFA pressure through ACSL4/LPCAT3/ALOX15; `ferroptosis_stress_nonoverlap_score` (TFRC/POR/CYB5R1) is reported as a sensitivity analysis.",
        "- Requested models: `neuronal_score ~ pufa_aa_pressure_score + tumor_region` and `ferroptosis_stress_score ~ pufa_aa_pressure_score * neural_containing + tumor_region`, with HC3 standard errors.",
        "",
        "## Data audit",
        "",
        f"- H5 bins: **{audit['n_bins_h5_raw']:,}**; analyzed non-zero-UMI bins: **{audit['n_bins_analyzed']:,}**; annotation overlap: **{audit['annotation_overlap']:,}**; coordinate overlap: **{audit['coordinate_overlap']:,}**.",
        f"- Stage: **{audit['stage']}**; early-versus-late stage comparison: **not estimable in this one-section run**.",
        f"- Neural-containing bins (≥2 neuronal genes): **{audit['n_neural_containing']:,}** ({audit['neural_containing_fraction']:.2%}).",
        "",
        "## Genes used",
        "",
        *[f"- {name}: {', '.join(present[name]) or 'none detected'}" for name in present],
        "",
        "## Region summary",
        "",
        *region_lines,
        "",
        "## Descriptive correlations",
        "",
        *corr_lines,
        "",
        "## Requested exploratory models",
        "",
        *model_lines,
        "",
        "## Interpretation guardrails",
        "",
        "This is a one-section spatial screen. Bin-level p-values are descriptive and can be strongly anti-conservative because neighboring bins are not independent biological replicates. A positive association would motivate multi-patient spatial replication and/or direct neuronal and lipid-peroxidation measurements; it does not establish ENS depletion, AA flux, ferroptosis, causality or stage progression.",
        "",
    ]
    (out / "crc_spatial_pufa_neuronal_ferroptosis_minimal_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    required = [H5_PATH, META_PATH, POSITION_PATH, ANNOTATION_PATH]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required local files: " + "; ".join(missing))

    all_gene_sets = {**GENE_SETS, "ferroptosis_stress_nonoverlap": FERROPTOSIS_STRESS_NONOVERLAP}
    requested = list(dict.fromkeys(gene for genes in all_gene_sets.values() for gene in genes))
    raw, totals, barcodes, present = read_targeted_h5(H5_PATH, all_gene_sets)
    meta = read_gz_parquet(META_PATH, ["barcode", "tissue", "X", "Y", "DeconvolutionClass", "DeconvolutionLabel1", "DeconvolutionLabel2", "Periphery", "UnsupervisedL1", "UnsupervisedL2"])
    positions = read_gz_parquet(POSITION_PATH, ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    annotation = pd.read_csv(ANNOTATION_PATH, sep="\t", header=None, names=["barcode", "pathologist_label"])
    annotation["barcode"] = annotation["barcode"].astype(str)

    h5_barcodes = pd.Index(barcodes, name="barcode")
    meta = meta.set_index("barcode").reindex(h5_barcodes).reset_index()
    positions = positions.set_index("barcode").reindex(h5_barcodes).reset_index()
    labels = annotation.set_index("barcode")["pathologist_label"].reindex(h5_barcodes).fillna("unannotated")
    if meta["barcode"].isna().any() or positions["barcode"].isna().any():
        raise ValueError("Barcode reindexing produced missing rows")
    if labels.eq("unannotated").any():
        raise ValueError("Some H5 barcodes lack pathologist annotation")

    annotation_overlap = int(len(set(barcodes) & set(annotation["barcode"])))
    coordinate_overlap = int(len(set(barcodes) & set(positions["barcode"])))
    valid = np.isfinite(totals) & (totals > 0)
    raw = raw[valid]
    totals = totals[valid]
    barcodes = np.asarray(barcodes, dtype=object)[valid].tolist()
    meta = meta.loc[valid].reset_index(drop=True)
    positions = positions.loc[valid].reset_index(drop=True)
    labels = labels.iloc[np.flatnonzero(valid)].reset_index(drop=True)
    log_cpm = np.log1p(raw / totals[:, None] * 10_000.0)
    data = pd.DataFrame({"barcode": barcodes, "array_row": positions["array_row"].to_numpy(int), "array_col": positions["array_col"].to_numpy(int), "pxl_row": positions["pxl_row_in_fullres"].to_numpy(float), "pxl_col": positions["pxl_col_in_fullres"].to_numpy(float), "total_umi": totals, "pathologist_label": labels.to_numpy(str)})
    for name, genes in GENE_SETS.items():
        data[f"{name}_score"] = score_module(log_cpm, present[name], requested)
    data["ferroptosis_stress_nonoverlap_score"] = score_module(log_cpm, present["ferroptosis_stress_nonoverlap"], requested)
    data["pufa_aa_pressure_score"] = data[["pufa_incorporation_score", "aa_routing_score"]].mean(axis=1)
    data["ferroptosis_net_score"] = data["ferroptosis_stress_score"] - data["ferroptosis_defense_score"]
    data["tumor_region"] = infer_regions(positions, labels).to_numpy()
    data["neural_detected_count"] = (raw[:, [requested.index(gene) for gene in present["neuronal"]]] > 0).sum(axis=1)
    data["neural_containing"] = (data["neural_detected_count"] >= 2).astype(int)
    data["stage"] = "III-B"

    region_summary = data.groupby("tumor_region", sort=True).agg(
        n_bins=("barcode", "size"),
        pufa_aa_pressure_mean=("pufa_aa_pressure_score", "mean"),
        neuronal_mean=("neuronal_score", "mean"),
        ferroptosis_stress_mean=("ferroptosis_stress_score", "mean"),
        ferroptosis_stress_nonoverlap_mean=("ferroptosis_stress_nonoverlap_score", "mean"),
        ferroptosis_defense_mean=("ferroptosis_defense_score", "mean"),
    ).reset_index()
    all_bins = pd.Series(True, index=data.index)
    neural_bins = data["neural_containing"].eq(1)
    correlations = pd.DataFrame([
        correlation_row("PUFA/AA pressure vs neuronal score", "all tissue bins", data["pufa_aa_pressure_score"], data["neuronal_score"], all_bins),
        correlation_row("PUFA/AA pressure vs neuronal score", "neural-containing bins", data["pufa_aa_pressure_score"], data["neuronal_score"], neural_bins),
        correlation_row("PUFA/AA pressure vs ferroptosis stress", "neural-containing bins", data["pufa_aa_pressure_score"], data["ferroptosis_stress_score"], neural_bins),
        correlation_row("PUFA/AA pressure vs ferroptosis stress (non-overlap sensitivity)", "neural-containing bins", data["pufa_aa_pressure_score"], data["ferroptosis_stress_nonoverlap_score"], neural_bins),
    ])
    models = pd.concat([
        run_ols(data, "neuronal_score ~ pufa_aa_pressure_score + C(tumor_region)", "neuronal_score_model"),
        run_ols(data, "ferroptosis_stress_score ~ pufa_aa_pressure_score * neural_containing + C(tumor_region)", "ferroptosis_stress_interaction_model"),
        run_ols(data, "ferroptosis_stress_nonoverlap_score ~ pufa_aa_pressure_score * neural_containing + C(tumor_region)", "ferroptosis_stress_nonoverlap_interaction_model"),
    ], ignore_index=True)
    audit = {
        "dataset": "GSE280315",
        "sample": "P2CRC",
        "stage": "III-B",
        "technology": "Visium HD 8um bins",
        "n_bins_h5_raw": int(len(h5_barcodes)),
        "n_bins_analyzed": int(len(data)),
        "annotation_overlap": annotation_overlap,
        "coordinate_overlap": coordinate_overlap,
        "n_neural_containing": int(data["neural_containing"].sum()),
        "neural_containing_fraction": float(data["neural_containing"].mean()),
        "region_counts": data["tumor_region"].value_counts().to_dict(),
        "stage_comparison": "not_estimable; one CRC section, stage III-B",
        "h5_source": str(H5_PATH),
        "metadata_source": str(META_PATH),
        "position_source": str(POSITION_PATH),
        "annotation_source": str(ANNOTATION_PATH),
        "gene_sets": all_gene_sets,
        "present_genes": present,
        "normalization": "targeted raw UMI extraction; log1p(CPM per 10,000); within-section gene z-score; arithmetic mean module scores",
        "region_definition": "pathologist Neoplasm bins touching any non-Neoplasm 8-neighbor bin are invasive_front_proxy; remaining Neoplasm bins are tumor_core",
        "neural_definition": "at least two neuronal genes detected",
    }
    data.to_csv(args.out_dir / "crc_spatial_minimal_P2CRC_bin_scores.csv", index=False)
    region_summary.to_csv(args.out_dir / "crc_spatial_minimal_P2CRC_region_summary.csv", index=False)
    correlations.to_csv(args.out_dir / "crc_spatial_minimal_P2CRC_correlations.csv", index=False)
    models.to_csv(args.out_dir / "crc_spatial_minimal_P2CRC_model_terms.csv", index=False)
    (args.out_dir / "crc_spatial_minimal_P2CRC_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    make_map(data, args.out_dir / "crc_spatial_minimal_P2CRC_map.png")
    write_report(args.out_dir, audit, region_summary, correlations, models, present)
    print(json.dumps({"n_bins": len(data), "regions": audit["region_counts"], "neural_containing": audit["n_neural_containing"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
