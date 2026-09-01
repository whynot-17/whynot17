"""Compartment source map for the CRC antioxidant-buffering signal.

This extends the previous GCH1 source-map workflow to SLC7A11, GCH1 and a
five-gene antioxidant-buffering score in GSE200997 and GSE132465.  The primary
statistical unit is the patient within a tumor compartment; individual cells
are never treated as independent replicates.

Compartments:
    malignant_epithelial, myeloid, T_B and stromal.

GSE200997 uses the existing marker-defined malignant/immune/stromal calls and
raw-UMI lineage marker scores to split immune-like cells into myeloid versus
T/B. GSE132465 uses the official tumor Cell_type annotation.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind

from gch1_malignant_epithelial_pseudobulk import (
    GSE132465_ANNOT,
    GSE132465_MATRIX,
    GSE132465_SOFT,
    GSE200997_ANNOT,
    GSE200997_MATRIX,
    read_matrix_targets,
    parse_side_map,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
GSE200997_DIR = ROOT / "work" / "gse200997"
GSE200997_SELECTION = GSE200997_DIR / "tumor_marker_scores.csv"

BUFFER_GENES = ["SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH"]
CORE_COMPARTMENTS = ["malignant_epithelial", "myeloid", "T_B", "stromal"]
MYELOID_MARKERS = ["LST1", "TYROBP", "FCER1G", "CTSS", "LILRB1", "AIF1", "CSTA", "LGALS3", "S100A8", "S100A9", "FCGR3A"]
TB_MARKERS = ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD79A", "MS4A1", "CD37", "CD74", "HLA-DRA", "CD52", "LTB", "IL7R"]
GENES = sorted(set(BUFFER_GENES + MYELOID_MARKERS + TB_MARKERS))
METRICS = {
    "SLC7A11": "SLC7A11_log2_CPM_plus1",
    "GCH1": "GCH1_log2_CPM_plus1",
    "antioxidant_buffering_score": "antioxidant_buffering_score",
}


def zscore(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def bh(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out.tolist()
    idx = np.where(ok)[0]
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[idx] = restored
    return out.tolist()


def cell_marker_score(counts: dict[str, np.ndarray], totals: np.ndarray, genes: list[str], mask: np.ndarray) -> tuple[np.ndarray, list[str]]:
    available = [gene for gene in genes if gene in counts]
    if not available:
        return np.full(len(totals), np.nan), []
    scores = []
    denominator = np.maximum(totals.astype(float), 1.0)
    for gene in available:
        log_cpm = np.log2(counts[gene].astype(float) / denominator * 1e6 + 1.0)
        reference = log_cpm[mask]
        sd = reference.std(ddof=1)
        if np.isfinite(sd) and sd > 0:
            scores.append((log_cpm - reference.mean()) / sd)
    if not scores:
        return np.full(len(totals), np.nan), available
    return np.vstack(scores).mean(axis=0), available


def aggregate_compartment(
    ann: pd.DataFrame,
    counts: dict[str, np.ndarray],
    totals: np.ndarray,
    selected: np.ndarray,
    cohort: str,
    compartment: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    rows = []
    work = ann.copy()
    work["_selected"] = selected
    for (patient, sample, side), group in work[work["_selected"]].groupby(["patient", "sample", "side"], sort=True):
        idx = group.index.to_numpy()
        total_umi = int(totals[idx].sum())
        row = {
            "cohort": cohort,
            "patient": str(patient),
            "sample": str(sample),
            "side": str(side),
            "compartment": compartment,
            "n_cells": int(len(idx)),
            "selected_cells_total_UMI": total_umi,
        }
        for gene in BUFFER_GENES:
            umi = int(counts[gene][idx].sum())
            cpm = umi / total_umi * 1e6 if total_umi else np.nan
            row[f"{gene}_UMI"] = umi
            row[f"{gene}_CPM"] = cpm
            row[f"{gene}_log2_CPM_plus1"] = np.log2(cpm + 1) if pd.notna(cpm) else np.nan
            row[f"{gene}_detection_fraction"] = float(np.mean(counts[gene][idx] > 0)) if len(idx) else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    if out.empty:
        return out, {"cohort": cohort, "compartment": compartment, "n_patients": 0, "n_cells": 0}
    out = out.sort_values(["side", "patient"]).reset_index(drop=True)
    for gene in BUFFER_GENES:
        out[f"{gene}_z"] = zscore(out[f"{gene}_log2_CPM_plus1"])
    zcols = [f"{gene}_z" for gene in BUFFER_GENES]
    out["antioxidant_available_genes"] = out[zcols].notna().sum(axis=1)
    out["antioxidant_buffering_score"] = out[zcols].mean(axis=1, skipna=True)
    out.loc[out["antioxidant_available_genes"] < 4, "antioxidant_buffering_score"] = np.nan
    return out, {
        "cohort": cohort,
        "compartment": compartment,
        "n_patients": int(out["patient"].nunique()),
        "n_cells": int(out["n_cells"].sum()),
        "side_counts": out["side"].value_counts().to_dict(),
    }


def load_200997() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    ann = pd.read_csv(GSE200997_ANNOT).rename(columns={"Unnamed: 0": "cell_id"})
    ann["cell_id"] = ann["cell_id"].astype(str)
    ann["sample"] = ann["samples"].astype(str)
    ann["side"] = ann["Location"].astype(str).str.title()
    ann["condition"] = ann["Condition"].astype(str).str.title()
    ann["patient"] = ann["sample"].str.extract(r"T_cac(\d+)", expand=False).map(
        lambda x: f"Patient {int(x)}" if pd.notna(x) else np.nan
    )
    selection = pd.read_csv(GSE200997_SELECTION)
    selection["cell_id"] = selection["cell_id"].astype(str)
    if selection["cell_id"].tolist() != ann["cell_id"].tolist():
        selection = selection.set_index("cell_id").reindex(ann["cell_id"]).reset_index()
    ann["malignant_epithelial"] = selection["malignant_epithelial"].fillna(False).astype(bool).to_numpy()
    immune_like = (
        ~ann["malignant_epithelial"]
        & (selection["immune_score"].to_numpy(float) >= selection["epithelial_score"].to_numpy(float))
        & (selection["immune_score"].to_numpy(float) >= selection["stromal_score"].to_numpy(float))
        & (selection["immune_score"].to_numpy(float) > 0)
    )
    stromal_like = (
        ~ann["malignant_epithelial"]
        & ~immune_like
        & (selection["stromal_score"].to_numpy(float) > selection["epithelial_score"].to_numpy(float))
        & (selection["stromal_score"].to_numpy(float) > selection["immune_score"].to_numpy(float))
        & (selection["stromal_score"].to_numpy(float) > 0)
    )
    tumor = ann["condition"].eq("Tumor").to_numpy()
    cell_ids, counts, totals, n_rows = read_matrix_targets(GSE200997_MATRIX, b",", len(ann), set(GENES))
    if cell_ids != ann["cell_id"].tolist():
        raise ValueError("GSE200997 matrix cell order does not match annotation")
    reference = tumor & ~ann["malignant_epithelial"].to_numpy()
    myeloid_score, myeloid_available = cell_marker_score(counts, totals, MYELOID_MARKERS, reference)
    tb_score, tb_available = cell_marker_score(counts, totals, TB_MARKERS, reference)
    myeloid = immune_like & tumor & np.isfinite(myeloid_score) & (myeloid_score >= tb_score) & (myeloid_score > 0)
    tb = immune_like & tumor & np.isfinite(tb_score) & (tb_score > myeloid_score) & (tb_score > 0)
    malignant = ann["malignant_epithelial"].to_numpy() & tumor
    stromal = stromal_like & tumor
    lineage = np.full(len(ann), "unresolved", dtype=object)
    lineage[malignant] = "malignant_epithelial"
    lineage[myeloid] = "myeloid"
    lineage[tb] = "T_B"
    lineage[stromal] = "stromal"
    ann["lineage"] = lineage
    ann["myeloid_marker_score"] = myeloid_score
    ann["T_B_marker_score"] = tb_score
    rows, compartment_info = [], {}
    for compartment in CORE_COMPARTMENTS:
        selected = ann["lineage"].eq(compartment).to_numpy()
        pb, info = aggregate_compartment(ann, counts, totals, selected, "GSE200997", compartment)
        if not pb.empty:
            rows.append(pb)
        compartment_info[compartment] = info
    pb = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    audit = ann[["cell_id", "sample", "patient", "side", "condition", "lineage", "myeloid_marker_score", "T_B_marker_score"]].copy()
    info = {
        "cohort": "GSE200997",
        "annotation_cells": int(len(ann)),
        "matrix_gene_rows": int(n_rows),
        "compartment_info": compartment_info,
        "lineage_cell_counts": ann["lineage"].value_counts(dropna=False).to_dict(),
        "marker_genes_available": {"myeloid": myeloid_available, "T_B": tb_available},
        "selection": "Existing marker-defined malignant epithelial/immune/stromal calls; immune-like cells split by dominant raw-UMI myeloid versus T/B marker score.",
    }
    return pb, audit, info


def load_132465() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    ann = pd.read_csv(GSE132465_ANNOT, sep="\t").rename(columns={"Index": "cell_id"})
    ann["cell_id"] = ann["cell_id"].astype(str)
    ann["sample"] = ann["Sample"].astype(str)
    ann["patient"] = ann["Patient"].astype(str)
    side_map = parse_side_map(GSE132465_SOFT)
    ann["side"] = ann["sample"].map(side_map.set_index("sample")["side"])
    tumor = ann["Class"].eq("Tumor")
    type_map = {
        "Epithelial cells": "malignant_epithelial",
        "Myeloids": "myeloid",
        "T cells": "T_B",
        "B cells": "T_B",
        "Stromal cells": "stromal",
    }
    ann["lineage"] = ann["Cell_type"].map(type_map).fillna("unresolved")
    cell_ids, counts, totals, n_rows = read_matrix_targets(GSE132465_MATRIX, b"\t", len(ann), set(GENES))
    if cell_ids != ann["cell_id"].tolist():
        raise ValueError("GSE132465 matrix cell order does not match annotation")
    rows, compartment_info = [], {}
    for compartment in CORE_COMPARTMENTS:
        selected = (tumor & ann["lineage"].eq(compartment)).to_numpy()
        pb, info = aggregate_compartment(ann, counts, totals, selected, "GSE132465", compartment)
        if not pb.empty:
            rows.append(pb)
        compartment_info[compartment] = info
    pb = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    audit = ann[["cell_id", "sample", "patient", "side", "Class", "Cell_type", "Cell_subtype", "lineage"]].copy()
    info = {
        "cohort": "GSE132465",
        "annotation_cells": int(len(ann)),
        "matrix_gene_rows": int(n_rows),
        "compartment_info": compartment_info,
        "lineage_cell_counts": ann["lineage"].value_counts(dropna=False).to_dict(),
        "selection": "Official GEO annotation: Class=Tumor; Epithelial cells, Myeloids, T cells/B cells and Stromal cells mapped to the four compartments.",
        "side_mapping": "Existing GSE132465 GEO region map: cecum/ascending/hepatic flexure/transverse=Right; splenic flexure/descending/sigmoid/rectosigmoid/rectum=Left.",
    }
    return pb, audit, info


def sidedness_rows(pb: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cohort, compartment), data in pb.groupby(["cohort", "compartment"], sort=True):
        for metric, column in METRICS.items():
            right = data.loc[data["side"].eq("Right"), column].dropna().to_numpy(float)
            left = data.loc[data["side"].eq("Left"), column].dropna().to_numpy(float)
            row = {
                "cohort": cohort,
                "compartment": compartment,
                "metric": metric,
                "n_right": int(len(right)),
                "n_left": int(len(left)),
                "right_mean": float(right.mean()) if len(right) else np.nan,
                "left_mean": float(left.mean()) if len(left) else np.nan,
                "right_minus_left": float(right.mean() - left.mean()) if len(right) and len(left) else np.nan,
                "welch_p": float(ttest_ind(right, left, equal_var=False).pvalue) if len(right) >= 2 and len(left) >= 2 else np.nan,
                "mannwhitney_p": float(mannwhitneyu(right, left, alternative="two-sided").pvalue) if len(right) >= 2 and len(left) >= 2 else np.nan,
            }
            rows.append(row)
    out = pd.DataFrame(rows)
    if len(out):
        out["welch_fdr_within_cohort"] = out.groupby("cohort")["welch_p"].transform(lambda x: bh(x.tolist()))
        out["mannwhitney_fdr_within_cohort"] = out.groupby("cohort")["mannwhitney_p"].transform(lambda x: bh(x.tolist()))
    return out


def write_figure(sidedness: pd.DataFrame) -> None:
    cohorts = ["GSE200997", "GSE132465"]
    metrics = list(METRICS)
    fig, axes = plt.subplots(1, len(metrics), figsize=(15, 4.8), constrained_layout=True)
    for ax, metric in zip(axes, metrics, strict=False):
        table = sidedness[sidedness["metric"].eq(metric)].pivot(index="compartment", columns="cohort", values="right_minus_left").reindex(index=CORE_COMPARTMENTS, columns=cohorts)
        im = ax.imshow(table.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(len(cohorts)), cohorts, rotation=25, ha="right")
        ax.set_yticks(range(len(CORE_COMPARTMENTS)), CORE_COMPARTMENTS)
        ax.set_title(metric.replace("_", "\n"))
        for i, compartment in enumerate(CORE_COMPARTMENTS):
            for j, cohort in enumerate(cohorts):
                value = table.loc[compartment, cohort]
                if pd.notna(value):
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="right − left")
    fig.suptitle("CRC antioxidant-buffering source map: patient-level pseudobulk", fontsize=13)
    fig.savefig(OUT / "sc_antioxidant_compartment_source_map.png", dpi=220)
    plt.close(fig)


def write_report(sidedness: pd.DataFrame, info: dict[str, object]) -> None:
    lines = [
        "# CRC antioxidant-buffering compartment source map",
        "",
        "## Question",
        "",
        "Does the right-sided SLC7A11/GCH1/antioxidant-buffering signal localize to malignant epithelial cells, myeloid cells, T/B cells or stromal cells?",
        "",
        "## Design",
        "",
        "- Raw UMI counts were aggregated to patient × sample × compartment pseudobulk.",
        "- The statistical unit is the patient; cells were not treated as independent observations.",
        "- SLC7A11 and GCH1 are log2(CPM + 1) pseudobulk values.",
        "- Antioxidant buffering is the mean within-cohort/compartment z-score of SLC7A11, GPX4, AIFM2, GCH1 and DHODH, requiring at least 4 of 5 genes.",
        "- GSE200997 immune-like cells were split using dominant raw-UMI myeloid versus T/B marker scores; GSE132465 used official cell-type labels.",
        "",
        "## Right-versus-left results",
        "",
        "| Cohort | Compartment | Metric | Right n | Left n | Right−left | Welch P | MW P |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in sidedness.iterrows():
        delta = "NA" if pd.isna(row.right_minus_left) else f"{row.right_minus_left:.3f}"
        wp = "NA" if pd.isna(row.welch_p) else f"{row.welch_p:.3g}"
        mp = "NA" if pd.isna(row.mannwhitney_p) else f"{row.mannwhitney_p:.3g}"
        lines.append(f"| {row.cohort} | {row.compartment} | {row.metric} | {row.n_right} | {row.n_left} | {delta} | {wp} | {mp} |")
    lines += [
        "",
        "## Interpretation guardrails",
        "",
        "A replicated malignant-epithelial right-minus-left signal supports a tumor-cell-associated state. A signal confined to myeloid, T/B or stromal compartments supports a microenvironmental source. These data do not establish AA concentration, AA flux, BH4 abundance, ferroptosis resistance or functional dependency.",
        "",
        "## Provenance",
        "",
        f"- GSE200997: {json.dumps(info['GSE200997'], ensure_ascii=False, default=str)}",
        f"- GSE132465: {json.dumps(info['GSE132465'], ensure_ascii=False, default=str)}",
        "",
        "## Files",
        "",
        "- `sc_antioxidant_compartment_patient_pseudobulk.csv`",
        "- `sc_antioxidant_compartment_cell_lineage_audit.csv`",
        "- `sc_antioxidant_compartment_sidedness.csv`",
        "- `sc_antioxidant_compartment_cell_counts.csv`",
        "- `sc_antioxidant_compartment_source_map.png`",
    ]
    (OUT / "sc_antioxidant_compartment_source_map_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pb200, audit200, info200 = load_200997()
    pb132, audit132, info132 = load_132465()
    pb = pd.concat([pb200, pb132], ignore_index=True)
    audit200["cohort"] = "GSE200997"
    audit132["cohort"] = "GSE132465"
    audit = pd.concat([audit200, audit132], ignore_index=True)
    sidedness = sidedness_rows(pb)
    cell_counts = audit.groupby(["cohort", "lineage"], dropna=False).size().reset_index(name="n_cells")
    patient_counts = pb.groupby(["cohort", "compartment"], sort=True).agg(n_patient_rows=("patient", "size"), n_patients=("patient", "nunique"), n_right=("side", lambda x: int((x == "Right").sum())), n_left=("side", lambda x: int((x == "Left").sum()))).reset_index()
    cell_counts.to_csv(OUT / "sc_antioxidant_compartment_cell_counts.csv", index=False)
    audit.to_csv(OUT / "sc_antioxidant_compartment_cell_lineage_audit.csv", index=False)
    pb.to_csv(OUT / "sc_antioxidant_compartment_patient_pseudobulk.csv", index=False)
    sidedness.to_csv(OUT / "sc_antioxidant_compartment_sidedness.csv", index=False)
    patient_counts.to_csv(OUT / "sc_antioxidant_compartment_patient_counts.csv", index=False)
    write_figure(sidedness)
    info = {"GSE200997": info200, "GSE132465": info132}
    manifest = {
        "analysis": "SLC7A11/GCH1/antioxidant-buffering compartment source map",
        "access_date": "2026-09-02",
        "cohorts": info,
        "compartments": CORE_COMPARTMENTS,
        "buffer_genes": BUFFER_GENES,
        "statistical_unit": "patient within compartment",
        "outputs": {
            "pseudobulk": str(OUT / "sc_antioxidant_compartment_patient_pseudobulk.csv"),
            "cell_lineage_audit": str(OUT / "sc_antioxidant_compartment_cell_lineage_audit.csv"),
            "sidedness": str(OUT / "sc_antioxidant_compartment_sidedness.csv"),
            "cell_counts": str(OUT / "sc_antioxidant_compartment_cell_counts.csv"),
            "patient_counts": str(OUT / "sc_antioxidant_compartment_patient_counts.csv"),
            "figure": str(OUT / "sc_antioxidant_compartment_source_map.png"),
            "report": str(OUT / "sc_antioxidant_compartment_source_map_report.md"),
        },
    }
    (OUT / "sc_antioxidant_compartment_source_map_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    write_report(sidedness, info)
    print(json.dumps({"sidedness": sidedness.to_dict("records"), "cell_counts": cell_counts.to_dict("records"), "outputs": manifest["outputs"]}, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
