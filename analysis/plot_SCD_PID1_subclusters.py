#!/usr/bin/env python
"""Plot the current SCD/PID1 cell-subpopulation summary tables."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "2026-10-11" / "SCD_PID1_cell_subclusters"
SUMMARY = OUT / "subpopulation_expression_summary.tsv"
FIGURE = OUT / "SCD_PID1_subpopulation_summary.png"

df = pd.read_csv(SUMMARY, sep="\t")
colors = {"GSE132465": "#2468a2", "GSE200997": "#d05c3b"}
fig, axes = plt.subplots(1, 3, figsize=(16, 7))
fig.subplots_adjust(left=0.13, right=0.98, bottom=0.13, top=0.83, wspace=0.58)

epi = df[df.compartment == "Tumor_epithelial"].copy()
for ax, metric, title, xlabel, cmap in [
    (axes[0], "SCD_mean_log1p_10K", "A. SCD in tumor epithelial subtypes", "Mean log1p normalized SCD", None),
    (axes[1], "lipid_module_mean_z_excluding_SCD", "B. Lipid-metabolism module", "Mean gene-wise z score (SCD excluded)", None),
]:
    for cohort, g in epi.groupby("cohort", sort=True):
        y = [f"{r.subtype}  |  {cohort}" for r in g.itertuples()]
        if metric == "SCD_mean_log1p_10K":
            sizes = 35 + g.SCD_pct_detected.fillna(0).to_numpy() * 2.1
        else:
            sizes = 85 + g.SCD_pct_detected.fillna(0).to_numpy() * 1.2
        ax.scatter(g[metric], y, s=sizes, color=colors[cohort], alpha=0.82,
                   edgecolor="white", linewidth=0.8, label=cohort)
    if metric == "lipid_module_mean_z_excluding_SCD":
        ax.axvline(0, color="#777777", linestyle="--", linewidth=1)
    ax.set_title(title, loc="left", fontsize=11, weight="bold")
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", alpha=0.2)
    ax.tick_params(axis="y", labelsize=8)

fib = df[df.compartment == "Tumor_fibroblast"].copy()
ax = axes[2]
for cohort, g in fib.groupby("cohort", sort=True):
    y = [f"{r.subtype}  |  {cohort}" for r in g.itertuples()]
    sizes = 35 + g.PID1_pct_detected.fillna(0).to_numpy() * 5.0
    ax.scatter(g.PID1_mean_log1p_10K, y, s=sizes, color=colors[cohort], alpha=0.82,
               edgecolor="white", linewidth=0.8, label=cohort)
ax.set_title("C. PID1 in tumor fibroblast states", loc="left", fontsize=11, weight="bold")
ax.set_xlabel("Mean log1p normalized PID1")
ax.grid(axis="x", alpha=0.2)
ax.tick_params(axis="y", labelsize=8)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=2, frameon=False)
fig.suptitle("SCD and PID1 by CRC single-cell subpopulation", y=1.02, fontsize=15, weight="bold")
fig.text(0.5, 0.035, "Bubble size represents target-gene detection rate.", ha="center", fontsize=9, color="#555555")
fig.savefig(FIGURE, dpi=220, bbox_inches="tight")
print(FIGURE)
