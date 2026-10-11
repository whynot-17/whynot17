#!/usr/bin/env python
"""Summarize post-annotation SCD/PID1 localization without using them to cluster."""

from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("D:/DINP/whynot17-new-project")
OUT = ROOT / "results/2026-10-11/cellchat_patient_consistency"


def summarize(path: Path, cohort: str, unit_name: str) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t")
    rows = []
    for celltype, g in d.groupby("celltype", observed=True):
        rows.append({
            "cohort": cohort,
            "celltype": celltype,
            "n_units_with_cells": g["sample_id"].nunique(),
            "n_cells_total": int(g["n_cells"].sum()),
            "median_unit_mean_SCD_log1p": g["SCD_mean_log1p"].median(),
            "mean_unit_pct_SCD_detected": g["SCD_pct_detected"].mean(),
            "median_unit_mean_PID1_log1p": g["PID1_mean_log1p"].median(),
            "mean_unit_pct_PID1_detected": g["PID1_pct_detected"].mean(),
            "median_unit_lipid_score_no_SCD": g["lipid_score_median"].median(),
        })
    return pd.DataFrame(rows)


def main() -> None:
    p132 = OUT / "GSE132465_target_gene_patient_summary.tsv"
    p200 = OUT / "GSE200997_target_gene_sample_summary.tsv"
    a = pd.read_csv(p132, sep="\t")
    b = pd.read_csv(p200, sep="\t")
    # These prep artifacts use 0..1 for detection fractions in the first pass.
    # Normalize to percentages in the report irrespective of source version.
    for df in (a, b):
        for col in ("SCD_pct_detected", "PID1_pct_detected"):
            if df[col].max(skipna=True) <= 1.0:
                df[col] *= 100.0
    a.to_csv(p132, sep="\t", index=False)
    b.to_csv(p200, sep="\t", index=False)
    summary = pd.concat([
        summarize(p132, "GSE132465", "patient"),
        summarize(p200, "GSE200997", "tumor_sample"),
    ], ignore_index=True)
    summary.to_csv(OUT / "target_gene_state_summary.tsv", sep="\t", index=False)

    tests = []
    for cohort, df, unit_type in [("GSE132465", a, "patient"), ("GSE200997", b, "tumor_sample")]:
        epi = df[df["celltype"].eq("Epithelial")].dropna(subset=["SCD_mean_log1p", "lipid_score_mean"])
        if len(epi) >= 3:
            rho, p = spearmanr(epi["SCD_mean_log1p"], epi["lipid_score_mean"])
            loo = [spearmanr(epi.drop(index=i)["SCD_mean_log1p"], epi.drop(index=i)["lipid_score_mean"]).statistic for i in epi.index]
            tests.append({"cohort": cohort, "comparison": "Epithelial SCD mean vs lipid-program score excluding SCD",
                          "unit_type": unit_type, "n_units": len(epi), "statistic": rho, "p_value": p,
                          "median_difference": np.nan, "leave_one_out_stat_min": min(loo), "leave_one_out_stat_max": max(loo),
                          "note": "Spearman correlation of per-unit mean log1p expression and custom non-SCD lipid marker score"})
        my = df[df["celltype"].isin(["Monocyte", "Macrophage"])].pivot_table(
            index="sample_id", columns="celltype", values="PID1_mean_log1p", aggfunc="first")
        my = my.dropna(subset=["Monocyte", "Macrophage"])
        if len(my) >= 3:
            result = wilcoxon(my["Monocyte"], my["Macrophage"], alternative="two-sided", zero_method="wilcox")
            tests.append({"cohort": cohort, "comparison": "PID1 mean: Monocyte minus Macrophage",
                          "unit_type": unit_type, "n_units": len(my), "statistic": result.statistic,
                          "p_value": result.pvalue, "median_difference": float((my["Monocyte"] - my["Macrophage"]).median()),
                          "leave_one_out_stat_min": np.nan, "leave_one_out_stat_max": np.nan,
                          "note": "Paired patient/sample-level test; cell-level measurements are summarized within unit first"})
    test_df = pd.DataFrame(tests)
    if len(test_df):
        test_df["p_BH"] = test_df.groupby("cohort", observed=True)["p_value"].transform(lambda x: __import__("statsmodels.stats.multitest", fromlist=["multipletests"]).multipletests(x, method="fdr_bh")[1])
    test_df.to_csv(OUT / "target_gene_patient_level_tests.tsv", sep="\t", index=False)

    broad_path = OUT / "GSE200997_broad_marker_annotations.tsv.gz"
    if broad_path.exists():
        broad = pd.read_csv(broad_path, sep="\t")
        competition = broad.groupby(["sample_id", "broad_top_marker_class"], observed=True).size().rename("n_cells").reset_index()
        totals = broad.groupby("sample_id", observed=True).size().rename("n_tumor_cells").reset_index()
        competition = competition.merge(totals, on="sample_id", how="left")
        competition["pct_sample_cells"] = 100.0 * competition["n_cells"] / competition["n_tumor_cells"]
        competition.to_csv(OUT / "GSE200997_broad_marker_competition_counts.tsv", sep="\t", index=False)

    # Patient-level descriptive scatter: SCD was excluded from the score by construction.
    epi = a[a["celltype"].eq("Epithelial")]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(epi["lipid_score_mean"], epi["SCD_mean_log1p"], color="#4477AA", s=35, alpha=.8)
    if len(epi) >= 3:
        rho, p = spearmanr(epi["lipid_score_mean"], epi["SCD_mean_log1p"])
        ax.text(.03, .97, f"Spearman ρ={rho:.2f}, p={p:.3g}; n={len(epi)} patients", transform=ax.transAxes, va="top")
    ax.set_xlabel("Mean lipid marker score (SCD excluded)")
    ax.set_ylabel("Mean SCD log1p-normalized expression")
    ax.set_title("GSE132465 tumor epithelium, patient-level")
    fig.tight_layout()
    fig.savefig(OUT / "GSE132465_SCD_vs_lipid_score.png", dpi=220)
    plt.close(fig)

    my = a[a["celltype"].isin(["Monocyte", "Macrophage"])].pivot_table(index="sample_id", columns="celltype", values="PID1_mean_log1p", aggfunc="first").dropna()
    fig, ax = plt.subplots(figsize=(6, 5))
    for _, row in my.iterrows():
        ax.plot([0, 1], [row["Monocyte"], row["Macrophage"]], color="#777777", alpha=.35, linewidth=.8)
    if len(my):
        ax.scatter(np.zeros(len(my)), my["Monocyte"], color="#CC6677", s=22, label="Monocyte")
        ax.scatter(np.ones(len(my)), my["Macrophage"], color="#4477AA", s=22, label="Macrophage")
    ax.set_xticks([0, 1], ["Monocyte", "Macrophage"])
    ax.set_ylabel("Patient mean PID1 log1p-normalized expression")
    ax.set_title("GSE132465; post-annotation PID1 localization")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "GSE132465_PID1_myeloid_subtypes.png", dpi=220)
    plt.close(fig)
    print(summary.to_string(index=False))
    print("\nPatient/sample-level tests:\n", test_df.to_string(index=False))


if __name__ == "__main__":
    main()
