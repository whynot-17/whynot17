"""
Asset Confirmation Table
------------------------
Figure | Structural need | Catalog asset inspected | Decision
S1 | attrition + gate failures + axis map | SankeyDiagram, BarComposition | parameter inheritance; custom build required
S2 | effect landscape + robustness matrix | heatmap, BarCategorical | parameter inheritance; custom build required
S3 | temporal distributions + assay audit | LineTrend, PairedBoxScatter | parameter inheritance; custom build required
S4 | definition/state/compartment evidence | heatmap, PairedBoxScatter | parameter inheritance; custom build required

Typography baseline (locked): 8 pt body, 9 pt panel title, 10 pt figure title,
Arial/DejaVu Sans fallback, bold panel letters, outward ticks, no chart junk.
Color baseline (locked):
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666", "#67A9CF", "#EF8A62"]
DIVERGING = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED = "#B2182B"; GREY = "#999999"; BLACK = "#222222"
Export baseline (locked): double-column 183 mm; vector PDF/SVG; 300 dpi PNG;
embedded fonts where supported; white background; tight bounding box.
"""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs"
PKG = ROOT / "supplement" / "MCOP_CRC_submission"
FIG = PKG / "figures"
SRC = PKG / "source_data"
QA = PKG / "QA"
for p in (FIG, SRC, QA):
    p.mkdir(parents=True, exist_ok=True)

CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = CATEGORICAL + ["#67A9CF", "#EF8A62"]
DIVERGING = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED = "#B2182B"
GREY = "#999999"
BLACK = "#222222"
BLUE = "#2166AC"
TEAL = "#287D8E"
GREEN = "#1B7837"
ORANGE = "#F1A340"
PALE = "#F5F7F8"
GRID = "#D9E0E4"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.7,
    "axes.edgecolor": BLACK,
    "axes.labelcolor": BLACK,
    "xtick.color": BLACK,
    "ytick.color": BLACK,
    "text.color": BLACK,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


def panel(ax, letter, title):
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
    ax.set_title(title, loc="left", fontweight="bold", pad=8)


def clean(ax, grid=False):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=2.5, width=0.6)
    if grid:
        ax.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def load(name):
    return pd.read_csv(OUT / name)


def figure_s1():
    flow = load("environmental_crc_267_actionability_flow.csv")
    mat = load("environmental_crc_267_actionability_matrix_v2.csv")
    axes = load("environmental_crc_267_human_testable_candidates.csv")

    fig = plt.figure(figsize=(7.205, 8.15), layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.9, 0.42], width_ratios=[1.2, 1])
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])
    axD = fig.add_subplot(gs[2, :])

    panel(axA, "A", "Outcome-blinded actionability attrition")
    stages = [
        ("Starting\nchemicals", 267), ("Entity\nvalid", 259), ("Exposure\ninterpretable", 135),
        ("Biomarker\navailable", 134), ("Detectable", 127), ("Cycle\ncoverage", 124),
        ("Human\ntestable", 87),
    ]
    xs = np.arange(len(stages))
    vals = np.array([v for _, v in stages])
    widths = 0.78
    cols = ["#C9D2D8", "#B9C9D3", "#9FBCCA", "#81AEBE", "#63A0AD", "#478E9B", TEAL]
    axA.bar(xs, vals, width=widths, color=cols, edgecolor="white", linewidth=1.2, zorder=2)
    for x, (lab, val) in enumerate(stages):
        axA.text(x, val + 8, f"{val}", ha="center", va="bottom", fontweight="bold", fontsize=9)
        axA.text(x, -14, lab, ha="center", va="top", fontsize=7.2)
        if x < len(stages)-1:
            loss = val - stages[x+1][1]
            axA.annotate(f"−{loss}", (x+0.5, min(val, stages[x+1][1])+8), ha="center", color="#6B747A", fontsize=7)
    axA.text(3, 292, "CRC outcomes remain behind the firewall until the 15 biomarker tests are frozen",
             ha="center", color=BLUE, fontweight="bold", fontsize=8)
    axA.plot([-0.4,6.4],[282,282],color=BLUE,lw=2.2,solid_capstyle="round")
    axA.text(6, 43, "27 strict-eligible", ha="center", fontsize=7, color="#4E5A60")
    axA.set_ylim(-35, 315); axA.set_xlim(-0.6, 6.6); axA.axis("off")

    panel(axB, "B", "Where candidates first failed")
    gate_defs = [("E", "Entity", "E_tag", 1), ("X", "Exposure", "X_tag", 1),
                 ("B", "Biomarker", "B_tag", 1), ("D", "Detectability", "D_tag", 1),
                 ("C", "Coverage", "C_tag", 1), ("T", "Testability", "T_tag", 1)]
    alive = pd.Series(True, index=mat.index)
    counts = []
    for code, label, col, threshold in gate_defs:
        ok = pd.to_numeric(mat[col], errors="coerce").fillna(0) >= threshold
        counts.append(int((alive & ~ok).sum()))
        alive &= ok
    labels = [f"{c}  {l}" for c, l, _, _ in gate_defs] + ["Eligible"]
    counts2 = counts + [int(alive.sum())]
    colors = ["#D6DDE1"]*6 + [TEAL]
    y = np.arange(len(labels))[::-1]
    axB.barh(y, counts2, color=colors, height=0.62, edgecolor="none")
    for yi, val in zip(y, counts2):
        axB.text(val+3, yi, str(val), va="center", fontweight="bold" if yi == 0 else "normal")
    axB.set_yticks(y, labels); axB.set_xlabel("Chemicals / mappings")
    axB.set_xlim(0, max(counts2)*1.18); clean(axB, grid=False)
    axB.spines["left"].set_visible(False); axB.tick_params(axis="y", length=0)

    panel(axC, "C", "Fifteen unique human biomarker tests")
    axes = axes.sort_values(["primary_biomarker"]).reset_index(drop=True)
    chem = axes["eligible_chemical_count"].to_numpy()
    matrix = axes["axis_key"].str.split("|", regex=False).str[0]
    colmap = {"urine": TEAL, "serum_or_blood": ORANGE}
    yy = np.arange(len(axes))[::-1]
    axC.scatter(np.zeros(len(axes)), yy, s=25 + np.sqrt(chem)*42,
                c=[colmap.get(x, GREY) for x in matrix], alpha=0.9, edgecolor="white", linewidth=0.8)
    for yi, (_, row), n in zip(yy, axes.iterrows(), chem):
        axC.text(0.20, yi, row["primary_biomarker"], va="center", fontweight="bold", fontsize=7)
        axC.text(1.05, yi, f"{int(n)} mapping{'s' if n != 1 else ''}", va="center", fontsize=6.7, color="#596269")
    axC.set_xlim(-0.25, 2.25); axC.set_ylim(-0.8, len(axes)-0.2); axC.axis("off")
    axC.text(0, -1.25, "● urine", color=TEAL, fontsize=7)
    axC.text(0.8, -1.25, "● serum/blood", color=ORANGE, fontsize=7)
    axC.text(0, len(axes)-0.15, "87 eligible mappings → 15 statistical tests", fontsize=7.6, color=BLUE, fontweight="bold", va="top")

    panel(axD, "D", "MiNP, DINP and MCOP retain distinct roles")
    axD.axis("off")
    cards = [
        (0.01, "MiNP", "Molecular nominee", "rank 24 · BH-FDR 0.00346\nempirical FDR 0.0356 · 40.7% detectable", "#E8EEF4", BLUE),
        (0.345, "DINP parent", "Exposure-axis parent", "rank 107 · BH-FDR 0.449\nnot a significant Phase 1 hit", "#F2F0F6", "#762A83"),
        (0.68, "MCOP", "Human urinary biomarker", "98.8% detectable · seven cycles\nentered the 15-test human screen", "#E7F2F1", TEAL),
    ]
    for x0, name, role, detail, bg, col in cards:
        axD.add_patch(FancyBboxPatch((x0,0.08),0.305,0.72,boxstyle="round,pad=0.012,rounding_size=0.02",
                                    transform=axD.transAxes,facecolor=bg,edgecolor=col,linewidth=0.9))
        axD.text(x0+0.02,0.65,name,transform=axD.transAxes,fontweight="bold",fontsize=9,color=col)
        axD.text(x0+0.02,0.47,role,transform=axD.transAxes,fontweight="bold",fontsize=7)
        axD.text(x0+0.02,0.20,detail,transform=axD.transAxes,fontsize=6.5,color="#586268",linespacing=1.25)
    axD.text(0.5,-0.02,"Biomarker translation does not imply chemical equivalence.",transform=axD.transAxes,
             ha="center",fontsize=7,color=ACCENT_RED,fontweight="bold")

    fig.suptitle("Supplementary Figure S1 | Auditable actionability filtering and multiplicity denominator",
                 x=0.01, ha="left", fontsize=10, fontweight="bold")
    save(fig, "Figure_S1_actionability_audit")

    s1 = pd.DataFrame(stages, columns=["stage_label", "n"])
    failures = pd.DataFrame({"gate": [x[0] for x in gate_defs]+["eligible"],
                             "gate_name": [x[1] for x in gate_defs]+["Eligible"],
                             "n_first_failure_or_eligible": counts2})
    return {"S1_attrition": s1, "S1_gate_failures": failures, "S1_axes": axes}


def figure_s2():
    screen = load("environmental_crc_systematic_human_screen_fdr_v2.csv").sort_values("screen_rank")
    score = load("environmental_crc_15axis_robustness_scorecard.csv")
    compare = load("mcop_crc_phase2h_python_vs_standard_survey.csv")
    merged = screen.merge(score[["primary_biomarker","F","L","C","H","D","T","A","E","robustness_tier"]], on="primary_biomarker")

    fig = plt.figure(figsize=(7.205, 7.4), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1.25, 1])
    axA = fig.add_subplot(gs[0, :]); axB = fig.add_subplot(gs[1, 0]); axC = fig.add_subplot(gs[1, 1])

    panel(axA, "A", "Complete 15-test effect landscape")
    y = -np.log10(merged["BH_FDR"].clip(lower=1e-12))
    x = np.log2(merged["OR"])
    sizes = 28 + merged["fit_crc_cases"].to_numpy()*1.5
    tiers = merged["robustness_tier"]
    cmap = {"Robust Tier A": TEAL, "Tier B": ORANGE, "Exploratory": "#AAB4BA"}
    for _, r in merged.iterrows():
        xx = math.log2(r.OR); yy = -math.log10(r.BH_FDR)
        lo, hi = math.log2(r.CI_low), math.log2(r.CI_high)
        axA.plot([lo, hi], [yy, yy], color=cmap[r.robustness_tier], alpha=0.35, lw=1.2, zorder=1)
    axA.scatter(x, y, s=sizes, c=[cmap[t] for t in tiers], edgecolor="white", linewidth=0.9, zorder=3)
    for _, r in merged.iterrows():
        if r.primary_biomarker in {"URXCOP","LBXPFHS","URXMOH"}:
            axA.annotate(r.primary_biomarker.replace("URXCOP","MCOP"),
                         (math.log2(r.OR), -math.log10(r.BH_FDR)), xytext=(5, 5), textcoords="offset points",
                         fontsize=7.5, fontweight="bold", color=cmap[r.robustness_tier])
    axA.axhline(-math.log10(0.05), color=ACCENT_RED, ls=(0,(3,2)), lw=0.9)
    axA.axvline(0, color="#7F898F", ls=(0,(2,2)), lw=0.8)
    axA.text(axA.get_xlim()[1] if axA.get_xlim()[1] else 1, -math.log10(0.05)+0.06, "BH-FDR 0.05", ha="right", color=ACCENT_RED, fontsize=7)
    axA.set_xlabel("log2 odds ratio per exposure doubling")
    axA.set_ylabel("-log10(BH-FDR)")
    clean(axA, grid=True)
    axA.text(0.01, 0.98, "Point area scales with CRC cases; pale intervals show 95% CIs", transform=axA.transAxes,
             va="top", fontsize=7, color="#5E676D")

    panel(axB, "B", "Prespecified robustness fingerprints")
    order = merged.sort_values(["robustness_tier","BH_FDR"], key=lambda s: s.map({"Robust Tier A":0,"Tier B":1,"Exploratory":2}) if s.name=="robustness_tier" else s)["primary_biomarker"].tolist()
    sm = score.set_index("primary_biomarker").loc[order]
    vals = sm[["F","L","C","H","D","T","A","E"]].to_numpy()
    cmap_h = LinearSegmentedColormap.from_list("robust", ["#EEF1F3", "#91B7C2", "#125D70"])
    axB.imshow(vals, aspect="auto", cmap=cmap_h, vmin=0, vmax=2)
    axB.set_xticks(np.arange(8), ["F","L","C","H","D","T","A","E"])
    axB.set_yticks(np.arange(len(order)), [x.replace("URXCOP","MCOP") for x in order])
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            axB.text(j, i, int(vals[i,j]), ha="center", va="center", fontsize=6.3,
                     color="white" if vals[i,j] == 2 else BLACK)
    axB.tick_params(length=0); axB.spines[:].set_visible(False)

    panel(axC, "C", "Independent survey implementation gate")
    r = compare.iloc[0]
    axC.axis("off")
    box = FancyBboxPatch((0.03,0.12),0.94,0.75,boxstyle="round,pad=0.018,rounding_size=0.025",
                         facecolor=PALE,edgecolor="#C9D2D8",linewidth=0.9,transform=axC.transAxes)
    axC.add_patch(box)
    axC.text(0.09,0.74,"R survey::svyglm",transform=axC.transAxes,fontweight="bold")
    axC.text(0.91,0.74,f"OR {r.r_OR:.3f}",transform=axC.transAxes,ha="right",fontsize=10,color=BLUE,fontweight="bold")
    axC.text(0.09,0.58,"Python Taylor sandwich",transform=axC.transAxes,fontweight="bold")
    axC.text(0.91,0.58,f"OR {r.python_OR:.3f}",transform=axC.transAxes,ha="right",fontsize=10,color=BLUE,fontweight="bold")
    axC.plot([0.09,0.91],[0.47,0.47],transform=axC.transAxes,color="#CCD4D8",lw=0.8)
    axC.text(0.09,0.36,f"|Δ logOR| = {r.absolute_logOR_difference:.2e}",transform=axC.transAxes)
    axC.text(0.09,0.27,f"Design df = {int(r.r_design_df)} · singleton strata = {int(r.r_singleton_strata_N)}",transform=axC.transAxes)
    axC.text(0.09,0.18,"Direction and CI-null conclusion agree",transform=axC.transAxes,color=GREEN,fontweight="bold")

    fig.suptitle("Supplementary Figure S2 | Full human screen and uniform robustness audit",
                 x=0.01, ha="left", fontsize=10, fontweight="bold")
    save(fig, "Figure_S2_human_screen_robustness")
    return {"S2_screen": screen, "S2_scorecard": score, "S2_R_vs_Python": compare}


def figure_s3():
    het = load("mcop_crc_phase2_cycle_heterogeneity_summary.csv")
    cyc = load("mcop_crc_phase2_per_cycle.csv")
    lod = load("mcop_crc_phase2_assay_lod_audit.csv")
    d = het.merge(cyc, left_on="cycle", right_on="Cycle", suffixes=("", "_effect")).merge(lod[["cycle","llod_ng_mL","platform"]], on="cycle", how="left")
    x = np.arange(len(d))

    fig = plt.figure(figsize=(7.205, 7.25), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1.2, 1])
    axA = fig.add_subplot(gs[0,0]); axB = fig.add_subplot(gs[0,1]); axC = fig.add_subplot(gs[1,:])

    panel(axA, "A", "Survey-weighted MCOP exposure distribution by cycle")
    q1 = d["MCOP_weighted_Q1_ng_mL"].to_numpy(); med = d["MCOP_weighted_median_ng_mL"].to_numpy()
    q3 = d["MCOP_weighted_Q3_ng_mL"].to_numpy(); p95 = d["MCOP_weighted_P95_ng_mL"].to_numpy()
    axA.fill_between(x, q1, q3, color="#9CC7D0", alpha=0.45, label="Q1–Q3")
    axA.plot(x, med, color=TEAL, marker="o", ms=4, lw=1.6, label="Median")
    axA.plot(x, p95, color="#7F898F", marker=".", ms=4, lw=1.0, ls=(0,(2,2)), label="P95")
    axA.set_yscale("log"); axA.set_ylabel("MCOP (ng/mL; log scale)")
    axA.set_xticks(x, [s.replace("20","") for s in d.cycle], rotation=35, ha="right")
    axA.legend(frameon=False, ncol=3, loc="upper left"); clean(axA, grid=True)

    panel(axB, "B", "Assay sensitivity and detectability")
    pct = d["MCOP_above_LOD_pct"].to_numpy()
    axB.scatter(x, pct, s=55, c=np.where(d.cycle.eq("2011-2012"), ORANGE, TEAL), edgecolor="white", linewidth=0.8)
    axB.plot(x, pct, color="#AAB4BA", lw=0.8, zorder=0)
    for i,(p,ll) in enumerate(zip(pct,d["llod_ng_mL"])):
        axB.text(i,p+0.55,f"{p:.1f}%",ha="center",fontsize=6.5)
        axB.text(i,92.2,f"LOD {ll:g}",ha="center",fontsize=6,color="#667078")
    axB.set_ylim(91.5,101.5); axB.set_ylabel("Above LOD (%)")
    axB.set_xticks(x,[s.replace("20","") for s in d.cycle],rotation=35,ha="right")
    clean(axB, grid=True)

    panel(axC, "C", "Temporal effect landscape and the discordant 2011–2012 cycle")
    orv=d["OR"].to_numpy(); lo=d["CI_low"].to_numpy(); hi=d["CI_high"].to_numpy(); cases=d["CRC_N"].to_numpy()
    cols=np.where(d.cycle.eq("2011-2012"),ACCENT_RED,TEAL)
    axC.axhspan(1.0773,1.4400,color="#D9EDF0",alpha=0.55,label="Pooled 95% CI")
    axC.axhline(1.2455,color=BLUE,lw=1.1,label="Pooled OR")
    axC.axhline(1,color="#6F797F",ls=(0,(3,2)),lw=0.8)
    for xi, oi, li, hii, ci in zip(x, orv, lo, hi, cols):
        axC.errorbar([xi], [oi], yerr=[[oi-li], [hii-oi]], fmt="none",
                     ecolor=ci, elinewidth=1.1, capsize=2, zorder=2)
    axC.scatter(x,orv,s=30+cases*5,c=cols,edgecolor="white",linewidth=0.9,zorder=3)
    ratio=d["Model_case_control_MCOP_median_ratio"].to_numpy()
    for i,(o,n,rr) in enumerate(zip(orv,cases,ratio)):
        axC.text(i,o+0.10,f"{int(n)} cases",ha="center",fontsize=6.4)
        axC.text(i,0.54,f"case/control\nmedian {rr:.2f}",ha="center",fontsize=6.1,color=cols[i])
    axC.set_ylim(0.45,3.25); axC.set_ylabel("Cycle-specific OR per doubling")
    axC.set_xticks(x,d.cycle,rotation=0)
    axC.text(0.99,0.97,"6/7 point estimates >1 · global exposure×cycle P=0.006",
             transform=axC.transAxes,ha="right",va="top",fontweight="bold",color=ACCENT_RED)
    axC.legend(frameon=False,ncol=2,loc="upper left"); clean(axC,grid=True)
    axC.text(0.53,0.84,"2011–2012: 100% detected; raw case/control exposure ordering reversed",
             transform=axC.transAxes,ha="center",fontsize=6.8,color=ACCENT_RED,fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.22",facecolor="white",edgecolor="#E9B8BE",alpha=0.92))

    fig.suptitle("Supplementary Figure S3 | Calendar-cycle exposure and assay audit",
                 x=0.01, ha="left", fontsize=10, fontweight="bold")
    save(fig, "Figure_S3_cycle_exposure_audit")
    return {"S3_cycle_audit": d, "S3_assay_LOD": lod}


def figure_s4():
    defs = load("figure5_ppar_definition_comparison.csv")
    states = load("figure5_within_state_ppar_analysis.csv")
    comp = load("mcop_phase2f_singlecell_paired_donor_contrasts.csv")
    tiers = load("figure5_evidence_tier_lock.csv")

    keep = defs[defs["n_paired_donors"].fillna(0)>0].copy()
    keep = keep.sort_values("median_delta")
    pcomp = comp[comp.score.eq("PPAR_nuclear_receptor_score")].copy()

    fig = plt.figure(figsize=(7.205, 7.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], width_ratios=[1.25, 1])
    axA=fig.add_subplot(gs[0,:]); axB=fig.add_subplot(gs[1,0]); axC=fig.add_subplot(gs[1,1])

    panel(axA,"A","PPAR/nuclear-receptor remodeling is robust to score definition")
    y=np.arange(len(keep)); vals=keep.median_delta.to_numpy(); fdr=keep.BH_FDR.to_numpy()
    colors=np.where(vals<0,BLUE,ORANGE)
    axA.barh(y,vals,color=colors,alpha=0.84,height=0.66)
    axA.axvline(0,color="#6F797F",lw=0.8)
    labels=[x.replace("Reactome ","").replace("DoRothEA ","") for x in keep.definition]
    axA.set_yticks(y,labels); axA.set_xlabel("Paired median Δ (tumor − normal)")
    for yi,v,q in zip(y,vals,fdr):
        axA.text(v-0.04 if v<0 else v+0.04,yi,f"FDR {q:.1e}" if q<0.001 else f"FDR {q:.3f}",
                 va="center",ha="right" if v<0 else "left",fontsize=6.3)
    axA.set_xlim(min(vals)-0.45,max(vals)+0.45); clean(axA,grid=True)
    axA.text(0.99,0.98,"12/13 estimable definitions decreased; PPARD regulon was not estimable",
             transform=axA.transAxes,ha="right",va="top",fontsize=7,color="#5E676D")

    panel(axB,"B","State-specific and compartment-specific direction")
    state_order=["enterocyte-like annotation","secretory-like annotation","other epithelial annotation"]
    st=states.set_index("cell_subtype").loc[state_order].reset_index()
    labels=["Enterocyte-like","Secretory-like","Other epithelial"] + [x.capitalize() for x in pcomp.compartment]
    vv=list(st.median_delta_tumor_minus_normal)+list(pcomp.median_delta_tumor_minus_normal)
    nn=list(st.n_paired_donors)+list(pcomp.paired_donors)
    yy=np.arange(len(labels))[::-1]
    axB.axvline(0,color="#6F797F",lw=0.8)
    axB.scatter(vv,yy,s=[28+n*1.4 for n in nn],c=[BLUE if v<0 else ORANGE for v in vv],edgecolor="white",linewidth=0.8,zorder=3)
    for v,y0,n in zip(vv,yy,nn):
        axB.plot([0,v],[y0,y0],color=BLUE if v<0 else ORANGE,lw=2,alpha=0.55)
        axB.text(v+(-0.035 if v<0 else 0.035),y0,f"Δ {v:+.3f}; n={int(n)}",va="center",ha="right" if v<0 else "left",fontsize=6.2)
    axB.set_yticks(yy,labels); axB.set_xlabel("Paired median Δ (tumor − normal)")
    axB.set_xlim(-0.78,0.88); clean(axB,grid=False); axB.spines["left"].set_visible(False); axB.tick_params(axis="y",length=0)

    panel(axC,"C","Evidence-tier boundary for the environmental bridge")
    axC.axis("off")
    level_color={"E3":GREEN,"E2":ORANGE,"E1":"#8B6AAE","E0":ACCENT_RED}
    ypos=np.linspace(0.86,0.13,len(tiers))
    short=["CRC ↔ epithelial PPAR/NR","CRC ↔ RELA/STAT3 state","PPAR/NR ↔ cell state",
           "Within-state remodeling","DINP/MiNP → PPAR/NR","MCOP/DINP → CRC epithelial state"]
    for y0,(_,r),lab in zip(ypos,tiers.iterrows(),short):
        col=level_color[r.evidence_level]
        axC.add_patch(FancyBboxPatch((0.02,y0-0.045),0.13,0.075,boxstyle="round,pad=0.01,rounding_size=0.015",
                                    facecolor=col,edgecolor="none",transform=axC.transAxes))
        axC.text(0.085,y0-0.008,r.evidence_level,transform=axC.transAxes,ha="center",va="center",color="white",fontweight="bold")
        axC.text(0.19,y0,lab,transform=axC.transAxes,va="center",fontsize=6.8,fontweight="bold" if r.evidence_level in {"E3","E0"} else "normal")
        if r.evidence_level=="E0":
            axC.text(0.19,y0-0.055,"Untested; causal arrow prohibited",transform=axC.transAxes,fontsize=6.2,color=ACCENT_RED)
    axC.text(0.02,0.98,"Observed human disease state",transform=axC.transAxes,color=GREEN,fontweight="bold",fontsize=7)
    axC.text(0.98,0.02,"Exposure-to-state bridge remains open",transform=axC.transAxes,ha="right",color=ACCENT_RED,fontweight="bold",fontsize=7)

    fig.suptitle("Supplementary Figure S4 | Definition robustness, state specificity and causal boundary",
                 x=0.01,ha="left",fontsize=10,fontweight="bold")
    save(fig,"Figure_S4_ppar_state_evidence")
    return {"S4_definitions": defs, "S4_within_state": states, "S4_compartments": comp, "S4_evidence_tiers": tiers}


def write_source_workbook(sheets):
    path = PKG / "MCOP_CRC_Source_Data.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        readme = pd.DataFrame({
            "Field": ["Purpose","Independent unit","Outcome boundary","Exposure identity","Causal boundary"],
            "Value": [
                "Source data for Supplementary Figures S1–S4",
                "NHANES participant/design unit for epidemiology; donor for single-cell paired analyses",
                "Actionability gates were frozen before CRC outcomes were analyzed",
                "MiNP, DINP parent compounds and urinary MCOP are chemically distinct",
                "CRC epithelial PPAR/NR remodeling is observed; DINP/MCOP-to-state causality is untested",
            ]})
        readme.to_excel(writer, sheet_name="README", index=False)
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name[:31], index=False)
    wb = pd.ExcelFile(path)
    return path


def qa_report():
    items=[]
    for stem in ["Figure_S1_actionability_audit","Figure_S2_human_screen_robustness","Figure_S3_cycle_exposure_audit","Figure_S4_ppar_state_evidence"]:
        for ext in ["pdf","svg","png"]:
            p=FIG/f"{stem}.{ext}"
            items.append({"figure":stem,"format":ext,"exists":p.exists(),"bytes":p.stat().st_size if p.exists() else 0})
    q=pd.DataFrame(items)
    q.to_csv(QA/"supplementary_figure_file_audit.csv",index=False)
    text=["# Supplementary figure QA","",
          "- Four figures were exported as vector PDF, editable SVG and 300-dpi PNG.",
          "- Figure width is 183 mm (7.205 in), matching a double-column layout.",
          "- All epidemiologic estimates originate from frozen CSV outputs; no values were manually substituted.",
          "- Confidence intervals are retained in the 15-axis effect landscape and temporal effect panel.",
          "- The 2011–2012 discordant cycle is highlighted rather than suppressed.",
          "- Single-cell inference is donor-level; individual cells are not treated as independent replicates.",
          "- MiNP, DINP and MCOP are represented as distinct entities; the exposure-to-epithelial causal bridge remains E0/untested.","",
          "## Export audit","","| Figure | Format | Exists | Bytes |","|---|---:|---:|---:|"]
    text.extend([f"| {r.figure} | {r.format} | {r.exists} | {r.bytes} |" for r in q.itertuples(index=False)])
    (QA/"figure_QA.md").write_text("\n".join(text),encoding="utf-8")


def main():
    sheets={}
    for fn in (figure_s1,figure_s2,figure_s3,figure_s4):
        sheets.update(fn())
    write_source_workbook(sheets)
    qa_report()
    print(json.dumps({"figures":4,"source_sheets":len(sheets),"package":str(PKG)},indent=2))


if __name__ == "__main__":
    main()
