"""Map the cellular source of the bulk GCH1 signal in sided CRC.

This analysis uses raw UMI counts from GSE200997 and GSE132465.  It creates
patient-level pseudobulk for broad tumor compartments and tests GCH1 by
sidedness and against the PLA2G4A/PTGS2/PTGES AA-routing proxy within
right-sided patients.  Cells are not treated as independent replicates.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, ttest_ind

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
GSE200997_DIR = ROOT / "work" / "gse200997"
GSE132465_DIR = ROOT / "work" / "gse132465"

sys.path.insert(0, str(ROOT / "work"))
from gch1_malignant_epithelial_pseudobulk import (  # noqa: E402
    GSE132465_ANNOT,
    GSE132465_MATRIX,
    GSE132465_SOFT,
    GSE200997_ANNOT,
    GSE200997_MATRIX,
    GENES,
    aggregate_patient_pseudobulk,
    parse_side_map,
    read_matrix_targets,
)

GSE200997_SELECTION = GSE200997_DIR / "tumor_marker_scores.csv"


def bh(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q.tolist()
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    restored = np.empty_like(adj)
    restored[order] = np.minimum(adj, 1.0)
    q[ok] = restored
    return q.tolist()


def load_200997() -> tuple[pd.DataFrame, dict[str, object]]:
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
    ann["malignant_epithelial_like"] = selection["malignant_epithelial"].fillna(False).astype(bool).to_numpy()
    ann["immune_like"] = (
        ~ann["malignant_epithelial_like"]
        & (selection["immune_score"].to_numpy(float) >= selection["epithelial_score"].to_numpy(float))
        & (selection["immune_score"].to_numpy(float) >= selection["stromal_score"].to_numpy(float))
        & (selection["immune_score"].to_numpy(float) > 0)
    )
    ann["stromal_like"] = (
        ~ann["malignant_epithelial_like"]
        & ~ann["immune_like"]
        & (selection["stromal_score"].to_numpy(float) > selection["epithelial_score"].to_numpy(float))
        & (selection["stromal_score"].to_numpy(float) > selection["immune_score"].to_numpy(float))
        & (selection["stromal_score"].to_numpy(float) > 0)
    )
    cell_ids, counts, totals, n_rows = read_matrix_targets(GSE200997_MATRIX, b",", len(ann), set(GENES))
    if cell_ids != ann["cell_id"].tolist():
        raise ValueError("GSE200997 matrix cell order does not match annotation")
    rows = []
    compartments = ["malignant_epithelial_like", "immune_like", "stromal_like"]
    for compartment in compartments:
        selected = ann["condition"].eq("Tumor").to_numpy() & ann[compartment].to_numpy()
        if selected.sum() == 0:
            continue
        pb, _ = aggregate_patient_pseudobulk(ann, counts, totals, selected, "GSE200997")
        pb["compartment"] = compartment
        rows.append(pb)
    result = pd.concat(rows, ignore_index=True)
    info = {
        "cohort": "GSE200997",
        "annotation_cells": int(len(ann)),
        "matrix_gene_rows": int(n_rows),
        "compartment_cell_counts": {c: int((ann["condition"].eq("Tumor") & ann[c]).sum()) for c in compartments},
        "selection": "malignant_epithelial_like from existing marker rule; immune_like/stromal_like assigned by dominant lineage marker score.",
    }
    return result, info


def load_132465() -> tuple[pd.DataFrame, dict[str, object]]:
    ann = pd.read_csv(GSE132465_ANNOT, sep="\t").rename(columns={"Index": "cell_id"})
    ann["cell_id"] = ann["cell_id"].astype(str)
    ann["sample"] = ann["Sample"].astype(str)
    ann["patient"] = ann["Patient"].astype(str)
    side_map = parse_side_map(GSE132465_SOFT)
    ann["side"] = ann["sample"].map(side_map.set_index("sample")["side"])
    cell_ids, counts, totals, n_rows = read_matrix_targets(GSE132465_MATRIX, b"\t", len(ann), set(GENES))
    if cell_ids != ann["cell_id"].tolist():
        raise ValueError("GSE132465 matrix cell order does not match annotation")
    tumor = ann["Class"].eq("Tumor")
    mapping = {
        "tumor_epithelial": ann["Cell_type"].eq("Epithelial cells"),
        "T_cells": ann["Cell_type"].eq("T cells"),
        "B_cells": ann["Cell_type"].eq("B cells"),
        "myeloid": ann["Cell_type"].eq("Myeloids"),
        "stromal": ann["Cell_type"].eq("Stromal cells"),
        "mast": ann["Cell_type"].eq("Mast cells"),
    }
    rows = []
    for compartment, mask in mapping.items():
        selected = (tumor & mask).to_numpy()
        if selected.sum() == 0:
            continue
        pb, _ = aggregate_patient_pseudobulk(ann, counts, totals, selected, "GSE132465")
        pb["compartment"] = compartment
        rows.append(pb)
    result = pd.concat(rows, ignore_index=True)
    info = {
        "cohort": "GSE132465",
        "annotation_cells": int(len(ann)),
        "matrix_gene_rows": int(n_rows),
        "compartment_cell_counts": {c: int((tumor & m).sum()) for c, m in mapping.items()},
        "selection": "Official GEO annotation: Class=Tumor and Cell_type-specific compartment.",
    }
    return result, info


def sidedness_rows(pb: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for (cohort, compartment), d in pb.groupby(["cohort", "compartment"], sort=True):
        right = d.loc[d.side.eq("Right"), "GCH1_log2_CPM_plus1"].dropna().to_numpy(float)
        left = d.loc[d.side.eq("Left"), "GCH1_log2_CPM_plus1"].dropna().to_numpy(float)
        if len(right) < 2 or len(left) < 2:
            continue
        rows.append({
            "cohort": cohort,
            "compartment": compartment,
            "n_right": int(len(right)),
            "n_left": int(len(left)),
            "right_mean": float(right.mean()),
            "left_mean": float(left.mean()),
            "right_minus_left": float(right.mean() - left.mean()),
            "welch_p": float(ttest_ind(right, left, equal_var=False).pvalue),
            "mannwhitney_p": float(mannwhitneyu(right, left, alternative="two-sided").pvalue),
        })
    out = pd.DataFrame(rows)
    if len(out):
        out["welch_fdr_within_cohort"] = out.groupby("cohort")["welch_p"].transform(lambda x: bh(x.tolist()))
        out["mannwhitney_fdr_within_cohort"] = out.groupby("cohort")["mannwhitney_p"].transform(lambda x: bh(x.tolist()))
    return out


def coupling_rows(pb: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for (cohort, compartment), d in pb.groupby(["cohort", "compartment"], sort=True):
        pair = d.loc[d.side.eq("Right"), ["AA_routing_proxy_core", "GCH1_z"]].dropna()
        if len(pair) < 4 or pair["AA_routing_proxy_core"].nunique() < 2 or pair["GCH1_z"].nunique() < 2:
            continue
        rho, p = spearmanr(pair["AA_routing_proxy_core"], pair["GCH1_z"])
        rows.append({
            "cohort": cohort,
            "compartment": compartment,
            "n_right": int(len(pair)),
            "spearman_rho": float(rho),
            "spearman_p": float(p),
        })
    out = pd.DataFrame(rows)
    if len(out):
        out["spearman_fdr_within_cohort"] = out.groupby("cohort")["spearman_p"].transform(lambda x: bh(x.tolist()))
    return out


def write_figure(sidedness: pd.DataFrame, coupling: pd.DataFrame) -> None:
    cohorts = ["GSE200997", "GSE132465"]
    all_compartments = sorted(set(sidedness["compartment"]).union(coupling["compartment"]))
    fig, axes = plt.subplots(1, 2, figsize=(14, max(5.5, 0.55 * len(all_compartments) + 1.5)), constrained_layout=True)
    side_pivot = sidedness.pivot(index="compartment", columns="cohort", values="right_minus_left").reindex(index=all_compartments, columns=cohorts)
    im = axes[0].imshow(side_pivot.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    axes[0].set_xticks(range(len(cohorts)), cohorts, rotation=25, ha="right")
    axes[0].set_yticks(range(len(all_compartments)), all_compartments)
    axes[0].set_title("GCH1 right − left")
    axes[0].set_xlabel("Cohort")
    axes[0].set_ylabel("Tumor compartment")
    for i, comp in enumerate(all_compartments):
        for j, cohort in enumerate(cohorts):
            val = side_pivot.loc[comp, cohort]
            if pd.notna(val):
                axes[0].text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="mean log2(CPM+1) difference")

    corr_pivot = coupling.pivot(index="compartment", columns="cohort", values="spearman_rho").reindex(index=all_compartments, columns=cohorts)
    im2 = axes[1].imshow(corr_pivot.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    axes[1].set_xticks(range(len(cohorts)), cohorts, rotation=25, ha="right")
    axes[1].set_yticks(range(len(all_compartments)), all_compartments)
    axes[1].set_title("GCH1–AA proxy coupling in right-sided CRC")
    axes[1].set_xlabel("Cohort")
    for i, comp in enumerate(all_compartments):
        for j, cohort in enumerate(cohorts):
            val = corr_pivot.loc[comp, cohort]
            if pd.notna(val):
                axes[1].text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04, label="Spearman ρ")
    fig.suptitle("Cellular source map of the GCH1 bulk signal", fontsize=14)
    fig.savefig(OUT / "gch1_compartment_source_map.png", dpi=220)
    plt.close(fig)


def write_report(sidedness: pd.DataFrame, coupling: pd.DataFrame, info: dict[str, object]) -> None:
    lines = [
        "# GCH1 compartment source map",
        "",
        "## Question",
        "",
        "Does the bulk right-sided GCH1 signal localize to malignant/tumor epithelial cells, or is it mainly contributed by immune/stromal compartments?",
        "",
        "## Design",
        "",
        "- Patient-level pseudobulk from raw UMI counts; cells were not treated as independent replicates.",
        "- GSE200997: marker-defined malignant epithelial-like, immune-like and stromal-like compartments.",
        "- GSE132465: official tumor epithelial, T-cell, B-cell, myeloid, stromal and mast-cell compartments.",
        "- AA-routing proxy: mean within-compartment z-score of PLA2G4A, PTGS2 and PTGES.",
        "",
        "## Right-versus-left GCH1 by compartment",
        "",
        "| Cohort | Compartment | Right n | Left n | Right−left mean | Welch P | MW P |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in sidedness.iterrows():
        lines.append(f"| {row.cohort} | {row.compartment} | {row.n_right} | {row.n_left} | {row.right_minus_left:.3f} | {row.welch_p:.3g} | {row.mannwhitney_p:.3g} |")
    lines += [
        "",
        "## Coupling within right-sided patients",
        "",
        "| Cohort | Compartment | Right n | Spearman ρ | P | FDR |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for _, row in coupling.iterrows():
        lines.append(f"| {row.cohort} | {row.compartment} | {row.n_right} | {row.spearman_rho:.3f} | {row.spearman_p:.3g} | {row.spearman_fdr_within_cohort:.3g} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "A consistent epithelial signal would support a tumor-cell-associated GCH1 state. A signal confined to immune or stromal compartments would instead argue that the bulk result is not a malignant epithelial-intrinsic program. These analyses do not establish GCH1 dependency, BH4 abundance, AA flux or causality.",
        "",
        "## Provenance",
        "",
        f"- GSE200997: {info['GSE200997']}",
        f"- GSE132465: {info['GSE132465']}",
        "- GSE132465 sidedness was assigned from GEO sample region metadata.",
        "",
        "## Files",
        "",
        "- `gch1_compartment_source_map_patient_pseudobulk.csv`",
        "- `gch1_compartment_source_map_sidedness.csv`",
        "- `gch1_compartment_source_map_coupling.csv`",
        "- `gch1_compartment_source_map.png`",
    ]
    (OUT / "gch1_compartment_source_map_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pb200, info200 = load_200997()
    pb132, info132 = load_132465()
    pb = pd.concat([pb200, pb132], ignore_index=True)
    sidedness = sidedness_rows(pb)
    coupling = coupling_rows(pb)
    pb.to_csv(OUT / "gch1_compartment_source_map_patient_pseudobulk.csv", index=False)
    sidedness.to_csv(OUT / "gch1_compartment_source_map_sidedness.csv", index=False)
    coupling.to_csv(OUT / "gch1_compartment_source_map_coupling.csv", index=False)
    write_figure(sidedness, coupling)
    info = {"GSE200997": info200, "GSE132465": info132}
    manifest = {
        "analysis": "GCH1 compartment source map",
        "access_date": "2026-09-01",
        "gene": "GCH1",
        "aa_proxy": ["PLA2G4A", "PTGS2", "PTGES"],
        "statistical_unit": "patient",
        "cohorts": info,
        "outputs": {
            "pseudobulk": str(OUT / "gch1_compartment_source_map_patient_pseudobulk.csv"),
            "sidedness": str(OUT / "gch1_compartment_source_map_sidedness.csv"),
            "coupling": str(OUT / "gch1_compartment_source_map_coupling.csv"),
            "figure": str(OUT / "gch1_compartment_source_map.png"),
            "report": str(OUT / "gch1_compartment_source_map_report.md"),
        },
    }
    (OUT / "gch1_compartment_source_map_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    write_report(sidedness, coupling, info)
    print(json.dumps({
        "compartment_cell_counts": {k: v["compartment_cell_counts"] for k, v in info.items()},
        "sidedness": sidedness.to_dict("records"),
        "coupling": coupling.to_dict("records"),
        "outputs": manifest["outputs"],
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
