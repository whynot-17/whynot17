"""Single-cell pseudobulk validation of SLC7A11, MBOAT2 and MBOAT1.

The analysis is restricted to tumor malignant/tumor epithelial cells in
GSE200997 and GSE132465.  It tests right-versus-left expression and coupling
to the PLA2G4A/PTGS2/PTGES AA-routing proxy.  ALOX5/ALOX5AP are added only as a
sensitivity proxy.  Patients, not cells, are the statistical units.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, ttest_ind

from gch1_malignant_epithelial_pseudobulk import (
    GSE132465_ANNOT,
    GSE132465_MATRIX,
    GSE132465_SOFT,
    GSE200997_ANNOT,
    GSE200997_MATRIX,
    aggregate_patient_pseudobulk,
    parse_side_map,
    read_matrix_targets,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
GSE200997_DIR = ROOT / "work" / "gse200997"
GSE132465_DIR = ROOT / "work" / "gse132465"
GSE200997_SELECTION = GSE200997_DIR / "tumor_marker_scores.csv"

CANDIDATES = ["SLC7A11", "MBOAT2", "MBOAT1"]
AA_PROXY_CORE = ["PLA2G4A", "PTGS2", "PTGES"]
AA_PROXY_EXPANDED = AA_PROXY_CORE + ["ALOX5", "ALOX5AP"]
GENES = list(dict.fromkeys(CANDIDATES + AA_PROXY_EXPANDED))


def zscore(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def make_pseudobulk(ann: pd.DataFrame, counts: dict[str, np.ndarray], totals: np.ndarray, selected: np.ndarray, cohort: str) -> pd.DataFrame:
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
            "n_epithelial_cells": int(len(idx)),
            "selected_cells_total_UMI": total_umi,
        }
        for gene in GENES:
            umi = int(counts[gene][idx].sum())
            cpm = umi / total_umi * 1e6 if total_umi else np.nan
            row[f"{gene}_UMI"] = umi
            row[f"{gene}_log2_CPM_plus1"] = np.log2(cpm + 1) if pd.notna(cpm) else np.nan
        rows.append(row)
    out = pd.DataFrame(rows).sort_values(["side", "patient"]).reset_index(drop=True)
    for gene in GENES:
        out[f"{gene}_z"] = zscore(out[f"{gene}_log2_CPM_plus1"])
    out["AA_proxy_core"] = out[[f"{g}_z" for g in AA_PROXY_CORE]].mean(axis=1)
    out["AA_proxy_expanded"] = out[[f"{g}_z" for g in AA_PROXY_EXPANDED]].mean(axis=1)
    return out


def load_200997() -> tuple[pd.DataFrame, dict[str, object]]:
    ann = pd.read_csv(GSE200997_ANNOT).rename(columns={"Unnamed: 0": "cell_id"})
    ann["cell_id"] = ann["cell_id"].astype(str)
    ann["sample"] = ann["samples"].astype(str)
    ann["side"] = ann["Location"].astype(str).str.title()
    ann["condition"] = ann["Condition"].astype(str).str.title()
    ann["patient"] = ann["sample"].str.extract(r"T_cac(\d+)", expand=False).map(lambda x: f"Patient {int(x)}" if pd.notna(x) else np.nan)
    selection = pd.read_csv(GSE200997_SELECTION)
    selection["cell_id"] = selection["cell_id"].astype(str)
    if selection["cell_id"].tolist() != ann["cell_id"].tolist():
        selection = selection.set_index("cell_id").reindex(ann["cell_id"]).reset_index()
    selected = selection["malignant_epithelial"].fillna(False).astype(bool).to_numpy() & ann["condition"].eq("Tumor").to_numpy()
    cell_ids, counts, totals, n_rows = read_matrix_targets(GSE200997_MATRIX, b",", len(ann), set(GENES))
    if cell_ids != ann["cell_id"].tolist():
        raise ValueError("GSE200997 matrix cell order does not match annotation")
    pb = make_pseudobulk(ann, counts, totals, selected, "GSE200997")
    info = {
        "cohort": "GSE200997",
        "annotation_cells": int(len(ann)),
        "selected_cells": int(selected.sum()),
        "matrix_gene_rows": int(n_rows),
        "side_counts": pb["side"].value_counts().to_dict(),
        "selection_rule": "Existing marker-defined malignant_epithelial flag restricted to tumor cells.",
    }
    return pb, info


def load_132465() -> tuple[pd.DataFrame, dict[str, object]]:
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
    pb = make_pseudobulk(ann, counts, totals, selected, "GSE132465")
    info = {
        "cohort": "GSE132465",
        "annotation_cells": int(len(ann)),
        "selected_cells": int(selected.sum()),
        "matrix_gene_rows": int(n_rows),
        "side_counts": pb["side"].value_counts().to_dict(),
        "selection_rule": "Official annotation: Class=Tumor and Cell_type=Epithelial cells.",
    }
    return pb, info


def sidedness(pb: pd.DataFrame, proxy: str) -> pd.DataFrame:
    rows = []
    for cohort, d in pb.groupby("cohort", sort=True):
        right = d[d.side.eq("Right")]
        left = d[d.side.eq("Left")]
        for gene in CANDIDATES:
            r = right[f"{gene}_z"].dropna().to_numpy(float)
            l = left[f"{gene}_z"].dropna().to_numpy(float)
            rows.append({
                "cohort": cohort,
                "proxy": proxy,
                "gene": gene,
                "n_right": int(len(r)),
                "n_left": int(len(l)),
                "right_minus_left_mean_z": float(r.mean() - l.mean()),
                "welch_p": float(ttest_ind(r, l, equal_var=False).pvalue),
                "mannwhitney_p": float(mannwhitneyu(r, l, alternative="two-sided").pvalue),
            })
    out = pd.DataFrame(rows)
    out["welch_fdr_within_cohort"] = out.groupby("cohort")["welch_p"].transform(lambda x: _bh(x.tolist()))
    out["mannwhitney_fdr_within_cohort"] = out.groupby("cohort")["mannwhitney_p"].transform(lambda x: _bh(x.tolist()))
    return out


def _bh(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out.tolist()
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adj = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty_like(adj)
    restored[order] = np.minimum(adj, 1.0)
    out[ok] = restored
    return out.tolist()


def coupling(pb: pd.DataFrame, proxy: str) -> pd.DataFrame:
    rows = []
    for cohort, d in pb.groupby("cohort", sort=True):
        right = d[d.side.eq("Right")]
        for gene in CANDIDATES:
            pair = right[[proxy, f"{gene}_z"]].dropna()
            rho, p = spearmanr(pair[proxy], pair[f"{gene}_z"])
            rows.append({
                "cohort": cohort,
                "proxy": proxy,
                "gene": gene,
                "n_right": int(len(pair)),
                "spearman_rho": float(rho),
                "spearman_p": float(p),
            })
    out = pd.DataFrame(rows)
    out["spearman_fdr_within_cohort"] = out.groupby("cohort")["spearman_p"].transform(lambda x: _bh(x.tolist()))
    return out


def write_figure(sidedness_df: pd.DataFrame, coupling_df: pd.DataFrame) -> None:
    cohorts = ["GSE200997", "GSE132465"]
    genes = CANDIDATES
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    s = sidedness_df[sidedness_df.proxy.eq("AA_proxy_core")].pivot(index="gene", columns="cohort", values="right_minus_left_mean_z").reindex(index=genes, columns=cohorts)
    im = axes[0].imshow(s.to_numpy(dtype=float), cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    axes[0].set_xticks(range(len(cohorts)), cohorts, rotation=25, ha="right")
    axes[0].set_yticks(range(len(genes)), genes)
    axes[0].set_title("Epithelial expression: right − left")
    axes[0].set_xlabel("Cohort")
    for i, gene in enumerate(genes):
        for j, cohort in enumerate(cohorts):
            val = s.loc[gene, cohort]
            if pd.notna(val):
                axes[0].text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="mean within-cohort z")

    c = coupling_df[(coupling_df.proxy.eq("AA_proxy_core"))].pivot(index="gene", columns="cohort", values="spearman_rho").reindex(index=genes, columns=cohorts)
    im2 = axes[1].imshow(c.to_numpy(dtype=float), cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    axes[1].set_xticks(range(len(cohorts)), cohorts, rotation=25, ha="right")
    axes[1].set_yticks(range(len(genes)), genes)
    axes[1].set_title("AA-routing coupling within right-sided CRC")
    axes[1].set_xlabel("Cohort")
    for i, gene in enumerate(genes):
        for j, cohort in enumerate(cohorts):
            val = c.loc[gene, cohort]
            if pd.notna(val):
                axes[1].text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=9)
    fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04, label="Spearman ρ")
    fig.suptitle("Malignant/tumor epithelial validation of membrane-defense candidates", fontsize=13)
    fig.savefig(OUT / "sc_membrane_defense_candidate_pseudobulk.png", dpi=220)
    plt.close(fig)


def write_report(sidedness_df: pd.DataFrame, coupling_df: pd.DataFrame, info: dict[str, object]) -> None:
    lines = [
        "# SLC7A11, MBOAT2 and MBOAT1: epithelial single-cell pseudobulk",
        "",
        "## Scope",
        "",
        "This bounded analysis tests only SLC7A11, MBOAT2 and MBOAT1 in malignant/tumor epithelial patient-level pseudobulk from GSE200997 and GSE132465.",
        "",
        "- Right-versus-left expression was tested in each cohort.",
        "- AA-routing coupling was tested within right-sided patients.",
        "- The primary proxy is PLA2G4A/PTGS2/PTGES; the sensitivity proxy adds ALOX5/ALOX5AP.",
        "- Cells were not treated as independent observations.",
        "- The proxy is not measured AA concentration and no dependency/causality is inferred.",
        "",
        "## Primary right-versus-left results",
        "",
        "| Cohort | Gene | Right n | Left n | Right−left mean z | Welch P | FDR |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in sidedness_df[sidedness_df.proxy.eq("AA_proxy_core")].iterrows():
        lines.append(f"| {row.cohort} | {row.gene} | {row.n_right} | {row.n_left} | {row.right_minus_left_mean_z:.3f} | {row.welch_p:.3g} | {row.welch_fdr_within_cohort:.3g} |")
    lines += [
        "",
        "## Primary coupling within right-sided patients",
        "",
        "| Cohort | Gene | Right n | Spearman ρ | P | FDR |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for _, row in coupling_df[coupling_df.proxy.eq("AA_proxy_core")].iterrows():
        lines.append(f"| {row.cohort} | {row.gene} | {row.n_right} | {row.spearman_rho:.3f} | {row.spearman_p:.3g} | {row.spearman_fdr_within_cohort:.3g} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "A candidate is prioritized only when its direction is reasonably consistent across both single-cell cohorts. A positive coupling in one cohort without replication is treated as exploratory.",
        "These results test epithelial transcriptomic state, not membrane lipid composition, AA flux, ferroptosis resistance or functional dependency.",
        "",
        "## Provenance",
        "",
        f"- GSE200997: {info['GSE200997']}",
        f"- GSE132465: {info['GSE132465']}",
        "- Pseudobulk normalization: summed raw UMI counts per patient and compartment, followed by log2(CPM + 1).",
        "",
        "## Files",
        "",
        "- `sc_membrane_defense_candidate_patient_pseudobulk.csv`",
        "- `sc_membrane_defense_candidate_sidedness.csv`",
        "- `sc_membrane_defense_candidate_coupling.csv`",
        "- `sc_membrane_defense_candidate_pseudobulk.png`",
    ]
    (OUT / "sc_membrane_defense_candidate_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pb200, info200 = load_200997()
    pb132, info132 = load_132465()
    pb = pd.concat([pb200, pb132], ignore_index=True)
    sidedness_tables = [sidedness(pb200, "AA_proxy_core"), sidedness(pb200, "AA_proxy_expanded"), sidedness(pb132, "AA_proxy_core"), sidedness(pb132, "AA_proxy_expanded")]
    coupling_tables = []
    for data in [pb200, pb132]:
        coupling_tables.extend([coupling(data, "AA_proxy_core"), coupling(data, "AA_proxy_expanded")])
    sidedness_df = pd.concat(sidedness_tables, ignore_index=True)
    coupling_df = pd.concat(coupling_tables, ignore_index=True)
    pb.to_csv(OUT / "sc_membrane_defense_candidate_patient_pseudobulk.csv", index=False)
    sidedness_df.to_csv(OUT / "sc_membrane_defense_candidate_sidedness.csv", index=False)
    coupling_df.to_csv(OUT / "sc_membrane_defense_candidate_coupling.csv", index=False)
    write_figure(sidedness_df, coupling_df)
    info = {"GSE200997": info200, "GSE132465": info132}
    manifest = {
        "analysis": "single-cell malignant/tumor epithelial validation of SLC7A11, MBOAT2 and MBOAT1",
        "access_date": "2026-09-01",
        "candidates": CANDIDATES,
        "aa_proxy_core": AA_PROXY_CORE,
        "aa_proxy_expanded": AA_PROXY_EXPANDED,
        "statistical_unit": "patient",
        "cohorts": info,
        "outputs": {
            "pseudobulk": str(OUT / "sc_membrane_defense_candidate_patient_pseudobulk.csv"),
            "sidedness": str(OUT / "sc_membrane_defense_candidate_sidedness.csv"),
            "coupling": str(OUT / "sc_membrane_defense_candidate_coupling.csv"),
            "figure": str(OUT / "sc_membrane_defense_candidate_pseudobulk.png"),
            "report": str(OUT / "sc_membrane_defense_candidate_report.md"),
        },
    }
    (OUT / "sc_membrane_defense_candidate_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    write_report(sidedness_df, coupling_df, info)
    print(json.dumps({
        "primary_sidedness": sidedness_df[sidedness_df.proxy.eq("AA_proxy_core")].to_dict("records"),
        "primary_coupling": coupling_df[coupling_df.proxy.eq("AA_proxy_core")].to_dict("records"),
        "outputs": manifest["outputs"],
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
