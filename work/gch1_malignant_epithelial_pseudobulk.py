"""Patient-level single-cell validation of the GCH1 AA-routing association.

The analysis is intentionally restricted to GCH1.  It uses raw UMI matrices
from GSE200997 and GSE132465, selects tumor epithelial cells, aggregates by
patient, and tests:

1. right-versus-left GCH1 expression;
2. GCH1 coupling to a PLA2G4A/PTGS2/PTGES AA-routing proxy within right-sided
   CRC;
3. a descriptive right-sided proxy-high versus proxy-low comparison.

Cells are never treated as independent replicates.  The AA-routing proxy is a
transcriptomic proxy, not measured AA concentration.
"""

from __future__ import annotations

import csv
import gzip
import itertools
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, ttest_ind

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"

GSE200997_DIR = ROOT / "work" / "gse200997"
GSE132465_DIR = ROOT / "work" / "gse132465"
GSE200997_ANNOT = GSE200997_DIR / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"
GSE200997_MATRIX = GSE200997_DIR / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
GSE200997_SELECTION = GSE200997_DIR / "tumor_marker_scores.csv"
GSE132465_ANNOT = GSE132465_DIR / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
GSE132465_MATRIX = GSE132465_DIR / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
GSE132465_SOFT = GSE132465_DIR / "GSE132465_family.soft.gz"

GENES = ["GCH1", "PLA2G4A", "PTGS2", "PTGES"]
AA_PROXY = ["PLA2G4A", "PTGS2", "PTGES"]


def zscore(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / sd


def parse_side_map(path: Path) -> pd.DataFrame:
    text = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    rows = []
    for block in text.split("^SAMPLE = ")[1:]:
        title_match = re.search(r"!Sample_title = ([^\r\n]+)", block)
        region_match = re.search(r"!Sample_characteristics_ch1 = region: ([^\r\n]+)", block)
        if not title_match:
            continue
        sample = title_match.group(1).strip()
        if not sample.endswith("-T"):
            continue
        region = region_match.group(1).strip().lower() if region_match else None
        if region in {"cecum", "ascending", "hepatic flexure", "transverse"}:
            side = "Right"
        elif region in {"splenic flexure", "descending", "sigmoid", "rectosigmoid", "rectum"}:
            side = "Left"
        else:
            side = None
        rows.append({"sample": sample, "region": region, "side": side})
    out = pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)
    if out.empty or out["side"].isna().any():
        raise ValueError(f"Unmapped GSE132465 tumor regions: {out.to_dict('records')}")
    return out


def read_matrix_targets(path: Path, delimiter: bytes, n_cells: int, targets: set[str]) -> tuple[list[str], dict[str, np.ndarray], np.ndarray, int]:
    found: dict[str, np.ndarray] = {}
    totals = np.zeros(n_cells, dtype=np.int64)
    n_rows = 0
    with gzip.open(path, "rb") as handle:
        header = next(handle).decode("utf-8", errors="replace").rstrip("\r\n")
        if delimiter == b",":
            fields = next(csv.reader([header]))
        else:
            fields = header.split("\t")
        cell_ids = [str(x).strip().strip('"') for x in fields[1:]]
        if len(cell_ids) != n_cells:
            raise ValueError(f"{path.name}: matrix has {len(cell_ids)} cells; annotation has {n_cells}")
        sep_text = delimiter.decode("ascii")
        for raw_line in handle:
            if not raw_line.strip():
                continue
            sep = raw_line.find(delimiter)
            if sep < 0:
                continue
            gene = raw_line[:sep].decode("utf-8", errors="replace").strip().strip('"')
            vals = np.fromstring(raw_line[sep + 1 :], sep=sep_text, dtype=np.int32)
            if len(vals) != n_cells:
                raise ValueError(f"{path.name}: {gene} has {len(vals)} values; expected {n_cells}")
            totals += vals
            if gene in targets:
                found[gene] = vals.astype(np.int32, copy=True)
            n_rows += 1
    missing = sorted(targets.difference(found))
    if missing:
        raise ValueError(f"{path.name}: missing target genes {missing}")
    return cell_ids, found, totals, n_rows


def aggregate_patient_pseudobulk(
    ann: pd.DataFrame,
    counts: dict[str, np.ndarray],
    totals: np.ndarray,
    selected: np.ndarray,
    cohort: str,
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
            "n_malignant_epithelial_cells": int(len(idx)),
            "selected_cells_total_UMI": total_umi,
        }
        for gene in GENES:
            umi = int(counts[gene][idx].sum())
            cpm = umi / total_umi * 1e6 if total_umi else np.nan
            row[f"{gene}_UMI"] = umi
            row[f"{gene}_CPM"] = cpm
            row[f"{gene}_log2_CPM_plus1"] = np.log2(cpm + 1) if pd.notna(cpm) else np.nan
            row[f"{gene}_detection_fraction"] = float(np.mean(counts[gene][idx] > 0)) if len(idx) else np.nan
        rows.append(row)
    out = pd.DataFrame(rows).sort_values(["side", "patient"]).reset_index(drop=True)
    log_cols = [f"{g}_log2_CPM_plus1" for g in GENES]
    for gene in GENES:
        out[f"{gene}_z"] = zscore(out[f"{gene}_log2_CPM_plus1"])
    out["AA_routing_proxy_core"] = out[[f"{g}_z" for g in AA_PROXY]].mean(axis=1)
    out["GCH1_z"] = out["GCH1_z"]
    median = out["AA_routing_proxy_core"].median()
    out["AA_routing_proxy_core_high"] = out["AA_routing_proxy_core"] >= median
    info = {
        "cohort": cohort,
        "n_patients": int(len(out)),
        "side_counts": out["side"].value_counts().to_dict(),
        "n_cells": int(out["n_malignant_epithelial_cells"].sum()),
        "target_genes": GENES,
    }
    return out, info


def load_gse200997() -> tuple[pd.DataFrame, dict[str, object]]:
    ann = pd.read_csv(GSE200997_ANNOT).rename(columns={"Unnamed: 0": "cell_id"})
    ann["cell_id"] = ann["cell_id"].astype(str)
    ann["sample"] = ann["samples"].astype(str)
    ann["side"] = ann["Location"].astype(str).str.title()
    ann["condition"] = ann["Condition"].astype(str).str.title()
    ann["patient"] = ann["sample"].str.extract(r"T_cac(\d+)", expand=False).map(
        lambda x: f"Patient {int(x)}" if pd.notna(x) else np.nan
    )
    selection = pd.read_csv(GSE200997_SELECTION, usecols=["cell_id", "malignant_epithelial"])
    selection["cell_id"] = selection["cell_id"].astype(str)
    if selection["cell_id"].tolist() != ann["cell_id"].tolist():
        selection = selection.set_index("cell_id").reindex(ann["cell_id"]).reset_index()
    selected = selection["malignant_epithelial"].fillna(False).astype(bool).to_numpy()
    selected &= ann["condition"].eq("Tumor").to_numpy()
    cell_ids, counts, totals, n_rows = read_matrix_targets(GSE200997_MATRIX, b",", len(ann), set(GENES))
    if cell_ids != ann["cell_id"].tolist():
        raise ValueError("GSE200997 matrix cell order does not match annotation")
    pb, info = aggregate_patient_pseudobulk(ann, counts, totals, selected, "GSE200997")
    if info["side_counts"] != {"Left": 8, "Right": 8}:
        raise ValueError(f"GSE200997 expected 8/8 patients; got {info['side_counts']}")
    info.update({
        "annotation_cells": int(len(ann)),
        "selected_cells": int(selected.sum()),
        "matrix_gene_rows": int(n_rows),
        "selection_rule": "Previously generated marker-defined malignant_epithelial flag, restricted to tumor cells; see gse200997 tumor_marker_scores.csv.",
    })
    return pb, info


def load_gse132465() -> tuple[pd.DataFrame, dict[str, object]]:
    ann = pd.read_csv(GSE132465_ANNOT, sep="\t").rename(columns={"Index": "cell_id"})
    ann["cell_id"] = ann["cell_id"].astype(str)
    ann["sample"] = ann["Sample"].astype(str)
    ann["patient"] = ann["Patient"].astype(str)
    side_map = parse_side_map(GSE132465_SOFT)
    ann["side"] = ann["sample"].map(side_map.set_index("sample")["side"])
    selected = (ann["Class"].eq("Tumor") & ann["Cell_type"].eq("Epithelial cells")).to_numpy()
    cell_ids, counts, totals, n_rows = read_matrix_targets(GSE132465_MATRIX, b"\t", len(ann), set(GENES))
    if cell_ids != ann["cell_id"].tolist():
        raise ValueError("GSE132465 matrix cell order does not match annotation")
    pb, info = aggregate_patient_pseudobulk(ann, counts, totals, selected, "GSE132465")
    if info["side_counts"] != {"Left": 12, "Right": 11}:
        raise ValueError(f"GSE132465 expected 12/11 patients; got {info['side_counts']}")
    info.update({
        "annotation_cells": int(len(ann)),
        "selected_cells": int(selected.sum()),
        "matrix_gene_rows": int(n_rows),
        "selection_rule": "Official annotation: Class=Tumor and Cell_type=Epithelial cells; malignancy inferred from tumor sample origin.",
        "side_mapping": "cecum/ascending/hepatic flexure/transverse=Right; splenic flexure/descending/sigmoid/rectosigmoid/rectum=Left.",
    })
    return pb, info


def compare_sidedness(pb: pd.DataFrame) -> dict[str, object]:
    right = pb.loc[pb["side"].eq("Right"), "GCH1_log2_CPM_plus1"].dropna().to_numpy(float)
    left = pb.loc[pb["side"].eq("Left"), "GCH1_log2_CPM_plus1"].dropna().to_numpy(float)
    welch = ttest_ind(right, left, equal_var=False)
    mw = mannwhitneyu(right, left, alternative="two-sided")
    return {
        "cohort": pb["cohort"].iloc[0],
        "n_right": int(len(right)),
        "n_left": int(len(left)),
        "right_mean_log2_CPM_plus1": float(right.mean()),
        "left_mean_log2_CPM_plus1": float(left.mean()),
        "right_minus_left_mean_log2_CPM_plus1": float(right.mean() - left.mean()),
        "welch_p": float(welch.pvalue),
        "mannwhitney_p": float(mw.pvalue),
    }


def coupling(pb: pd.DataFrame) -> dict[str, object]:
    d = pb[pb["side"].eq("Right")][["AA_routing_proxy_core", "GCH1_z"]].dropna()
    rho, p = spearmanr(d["AA_routing_proxy_core"], d["GCH1_z"])
    return {
        "cohort": pb["cohort"].iloc[0],
        "scope": "Right",
        "n_patients": int(len(d)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
    }


def high_low(pb: pd.DataFrame) -> dict[str, object]:
    d = pb[pb["side"].eq("Right")]
    high = d.loc[d["AA_routing_proxy_core_high"], "GCH1_log2_CPM_plus1"].dropna().to_numpy(float)
    low = d.loc[~d["AA_routing_proxy_core_high"], "GCH1_log2_CPM_plus1"].dropna().to_numpy(float)
    welch = ttest_ind(high, low, equal_var=False)
    mw = mannwhitneyu(high, low, alternative="two-sided")
    return {
        "cohort": pb["cohort"].iloc[0],
        "scope": "Right",
        "n_high": int(len(high)),
        "n_low": int(len(low)),
        "high_minus_low_mean_log2_CPM_plus1": float(high.mean() - low.mean()),
        "welch_p": float(welch.pvalue),
        "mannwhitney_p": float(mw.pvalue),
    }


def write_figure(pb_all: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), constrained_layout=True)
    colors = {"Left": "#4C78A8", "Right": "#E45756"}
    offsets = {"GSE200997": -0.18, "GSE132465": 0.18}
    for cohort_i, cohort in enumerate(["GSE200997", "GSE132465"]):
        d = pb_all[pb_all.cohort.eq(cohort)]
        for side_i, side in enumerate(["Left", "Right"]):
            vals = d.loc[d.side.eq(side), "GCH1_log2_CPM_plus1"].to_numpy(float)
            x = side_i + offsets[cohort]
            axes[0].boxplot(
                vals, positions=[x], widths=0.25, patch_artist=True,
                boxprops={"facecolor": colors[side], "alpha": 0.2, "edgecolor": colors[side]},
                medianprops={"color": colors[side], "linewidth": 1.8},
                whiskerprops={"color": colors[side]}, capprops={"color": colors[side]},
                flierprops={"marker": "", "markersize": 0},
            )
            jitter = np.linspace(-0.07, 0.07, len(vals)) if len(vals) > 1 else np.array([0.0])
            axes[0].scatter(np.full(len(vals), x) + jitter, vals, s=28, color=colors[side], edgecolor="white", linewidth=0.5, zorder=3, label=cohort if side == "Left" else None)
    axes[0].set_xticks([0, 1], ["Left", "Right"])
    axes[0].set_ylabel("GCH1 log2(CPM + 1)\nmalignant epithelial pseudobulk")
    axes[0].set_title("Patient-level GCH1 expression")
    axes[0].legend(frameon=False, title="Dataset")
    axes[0].spines[["top", "right"]].set_visible(False)

    for cohort, marker in [("GSE200997", "o"), ("GSE132465", "s")]:
        d = pb_all[(pb_all.cohort.eq(cohort)) & (pb_all.side.eq("Right"))]
        axes[1].scatter(d["AA_routing_proxy_core"], d["GCH1_z"], s=42, marker=marker, color="#D55E00" if cohort == "GSE200997" else "#0072B2", alpha=0.8, label=cohort)
    axes[1].axhline(0, color="0.75", lw=1)
    axes[1].axvline(0, color="0.75", lw=1)
    axes[1].set_xlabel("AA-routing proxy\n(PLA2G4A/PTGS2/PTGES mean z)")
    axes[1].set_ylabel("GCH1 within-cohort z")
    axes[1].set_title("Coupling within right-sided CRC")
    axes[1].legend(frameon=False)
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.suptitle("Malignant epithelial single-cell validation of GCH1", fontsize=14)
    fig.savefig(OUT / "gch1_malignant_epithelial_pseudobulk.png", dpi=220)
    plt.close(fig)


def write_report(sidedness: pd.DataFrame, correlations: pd.DataFrame, splits: pd.DataFrame, infos: dict[str, object]) -> None:
    lines = [
        "# GCH1 in malignant epithelial single-cell pseudobulk",
        "",
        "## Question",
        "",
        "Does the right-sided CRC GCH1 signal remain detectable in malignant/tumor epithelial cells, and does it couple to the AA-routing transcriptomic proxy within right-sided patients?",
        "",
        "## Design",
        "",
        "- GSE200997: marker-defined malignant epithelial-like cells from tumor samples.",
        "- GSE132465: official `Class=Tumor` and `Cell_type=Epithelial cells` annotation; malignancy inferred from tumor sample origin.",
        "- Raw UMI counts were summed within each patient and epithelial compartment, then normalized as log2(CPM + 1).",
        "- Statistical unit: patient; individual cells were not treated as replicates.",
        "- AA-routing proxy: within-cohort mean z-score of PLA2G4A, PTGS2 and PTGES. It is not measured AA concentration.",
        "",
        "## Right-versus-left GCH1",
        "",
        "| Cohort | Right n | Left n | Right−left mean log2(CPM+1) | Welch P | Mann–Whitney P |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in sidedness.iterrows():
        lines.append(f"| {row.cohort} | {row.n_right} | {row.n_left} | {row.right_minus_left_mean_log2_CPM_plus1:.3f} | {row.welch_p:.3g} | {row.mannwhitney_p:.3g} |")
    lines += [
        "",
        "## AA-routing coupling within right-sided patients",
        "",
        "| Cohort | n | Spearman ρ | P |",
        "|---|---:|---:|---:|",
    ]
    for _, row in correlations.iterrows():
        lines.append(f"| {row.cohort} | {row.n_patients} | {row.spearman_rho:.3f} | {row.spearman_p:.3g} |")
    lines += [
        "",
        "## Right-sided proxy-high versus proxy-low",
        "",
        "| Cohort | High n | Low n | High−low GCH1 mean | Welch P | Mann–Whitney P |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in splits.iterrows():
        lines.append(f"| {row.cohort} | {row.n_high} | {row.n_low} | {row.high_minus_low_mean_log2_CPM_plus1:.3f} | {row.welch_p:.3g} | {row.mannwhitney_p:.3g} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "Persistence of a right-sided GCH1 signal in both epithelial pseudobulk datasets would support a tumor-cell-associated state rather than a purely bulk-composition explanation. Positive within-right coupling would support, but not prove, an AA-associated GCH1/BH4 adaptive program.",
        "Neither result demonstrates tissue AA enrichment, BH4 abundance, ferroptosis resistance, causality or GCH1 dependency. Those require lipid/metabolite measurements and functional perturbation.",
        "",
        "## Provenance",
        "",
        f"- GSE200997: {infos['GSE200997']}",
        f"- GSE132465: {infos['GSE132465']}",
        "- Sidedness for GSE132465 was assigned from GEO sample region metadata.",
        "",
        "## Files",
        "",
        "- `gch1_sc_patient_pseudobulk.csv`",
        "- `gch1_sc_sidedness.csv`",
        "- `gch1_sc_coupling.csv`",
        "- `gch1_sc_high_low.csv`",
        "- `gch1_malignant_epithelial_pseudobulk.png`",
    ]
    (OUT / "gch1_malignant_epithelial_pseudobulk_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pb200, info200 = load_gse200997()
    pb132, info132 = load_gse132465()
    pb_all = pd.concat([pb200, pb132], ignore_index=True)
    pb_all.to_csv(OUT / "gch1_sc_patient_pseudobulk.csv", index=False)
    sidedness = pd.DataFrame([compare_sidedness(pb200), compare_sidedness(pb132)])
    correlations = pd.DataFrame([coupling(pb200), coupling(pb132)])
    splits = pd.DataFrame([high_low(pb200), high_low(pb132)])
    sidedness.to_csv(OUT / "gch1_sc_sidedness.csv", index=False)
    correlations.to_csv(OUT / "gch1_sc_coupling.csv", index=False)
    splits.to_csv(OUT / "gch1_sc_high_low.csv", index=False)
    write_figure(pb_all)
    infos = {"GSE200997": info200, "GSE132465": info132}
    manifest = {
        "analysis": "GCH1 malignant epithelial patient-level pseudobulk validation",
        "access_date": "2026-09-01",
        "gene": "GCH1",
        "aa_proxy": AA_PROXY,
        "aa_proxy_definition": "mean within-cohort z-score of PLA2G4A, PTGS2 and PTGES",
        "statistical_unit": "patient",
        "cohorts": infos,
        "outputs": {
            "pseudobulk": str(OUT / "gch1_sc_patient_pseudobulk.csv"),
            "sidedness": str(OUT / "gch1_sc_sidedness.csv"),
            "coupling": str(OUT / "gch1_sc_coupling.csv"),
            "high_low": str(OUT / "gch1_sc_high_low.csv"),
            "figure": str(OUT / "gch1_malignant_epithelial_pseudobulk.png"),
            "report": str(OUT / "gch1_malignant_epithelial_pseudobulk_report.md"),
        },
    }
    (OUT / "gch1_malignant_epithelial_pseudobulk_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    write_report(sidedness, correlations, splits, infos)
    print(json.dumps({
        "sidedness": sidedness.to_dict("records"),
        "right_coupling": correlations.to_dict("records"),
        "right_high_low": splits.to_dict("records"),
        "outputs": manifest["outputs"],
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
