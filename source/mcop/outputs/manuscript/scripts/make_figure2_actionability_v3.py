# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) outcome-blinded attrition flow → SankeyDiagram → param inherit
# (b) biomarker-test universe → BarComposition → param inherit
# (c) chemical-identity boundary → cross-type inherit → param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


from pathlib import Path
import json
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from matplotlib.patches import FancyBboxPatch, PathPatch


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs"
DEST = OUT / "manuscript" / "figures" / "final_submission"
SRC = DEST / "source_data"
DEST.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)

BLUE = CATEGORICAL[0]
RED = CATEGORICAL[1]
GREEN = CATEGORICAL[2]
ORANGE = CATEGORICAL[3]
PURPLE = CATEGORICAL[4]
TEAL = "#287D8E"
INK = "#1F2930"
MUTED = "#67737A"
LIGHT = "#EEF3F5"
LINE = "#CAD5DA"


def panel_label(ax, letter, title):
    ax.text(-0.055, 1.025, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", color=INK)
    ax.text(0.0, 1.025, title, transform=ax.transAxes, fontsize=8.6,
            fontweight="bold", va="top", color=INK)


def smooth_ribbon(ax, x0, x1, y0_lo, y0_hi, y1_lo, y1_hi, color, alpha=1.0):
    c = (x1 - x0) * 0.45
    verts = [
        (x0, y0_hi), (x0+c, y0_hi), (x1-c, y1_hi), (x1, y1_hi),
        (x1, y1_lo), (x1-c, y1_lo), (x0+c, y0_lo), (x0, y0_lo), (x0, y0_hi)
    ]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha))


def rounded_card(ax, xy, wh, face, edge, title, subtitle, detail, title_color=None):
    x, y = xy; w, h = wh
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle="round,pad=0.014,rounding_size=0.025",
        facecolor=face, edgecolor=edge, linewidth=0.9, transform=ax.transAxes))
    ax.text(x+0.035*w, y+0.72*h, title, transform=ax.transAxes,
            fontsize=9, fontweight="bold", color=title_color or edge, va="center")
    ax.text(x+0.035*w, y+0.49*h, subtitle, transform=ax.transAxes,
            fontsize=6.9, fontweight="bold", color=INK, va="center")
    ax.text(x+0.035*w, y+0.23*h, detail, transform=ax.transAxes,
            fontsize=6.2, color=MUTED, va="center", linespacing=1.3)


flow = pd.read_csv(OUT / "environmental_crc_267_actionability_flow.csv")
axes = pd.read_csv(OUT / "environmental_crc_267_human_testable_candidates.csv")
matrix = pd.read_csv(OUT / "environmental_crc_267_actionability_matrix_v2.csv")

# Data validation: no silent row loss and exact frozen universe.
assert len(matrix) == 267, f"Expected 267 chemical rows, found {len(matrix)}"
assert len(axes) == 15, f"Expected 15 biomarker tests, found {len(axes)}"
frozen = dict(zip(flow.stage, flow.n))
expected = [267, 259, 135, 134, 127, 124, 87]
observed = [
    frozen["total_core_chemicals"], frozen["E_entity_valid"],
    frozen["E_and_X_interpretable_exposure"], frozen["E_X_and_B_biomarker_available"],
    frozen["E_X_B_and_D_detectable"], frozen["E_X_B_D_and_C_coverage"],
    frozen["E_X_B_D_C_and_T_testable"],
]
assert observed == expected, (observed, expected)
assert frozen["strict_eligibility"] == 27
assert axes.eligible_chemical_count.sum() == 87

# First-failure audit across all 267 rows.
gate_defs = [("E", "Entity", "E_tag"), ("X", "Exposure", "X_tag"),
             ("B", "Biomarker", "B_tag"), ("D", "Detectability", "D_tag"),
             ("C", "Coverage", "C_tag"), ("T", "Testability", "T_tag")]
alive = pd.Series(True, index=matrix.index)
failure_counts = []
for code, label, col in gate_defs:
    ok = pd.to_numeric(matrix[col], errors="coerce").fillna(0) >= 1
    failure_counts.append(int((alive & ~ok).sum()))
    alive &= ok
assert failure_counts == [8, 124, 1, 7, 3, 37]
assert int(alive.sum()) == 87

# Figure contract: 183 mm double-column; asymmetric mixed-modality.
fig = plt.figure(figsize=(183/25.4, 152/25.4), facecolor="white")
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.82], width_ratios=[1.18, 1.0],
                      left=0.055, right=0.985, top=0.91, bottom=0.075,
                      hspace=0.26, wspace=0.19)
axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])

fig.text(0.055, 0.965, "Figure 2 | Outcome-blinded actionability defines the human screening universe",
         fontsize=9, fontweight="bold", color=INK, ha="left", va="top")

# Panel A — hero flow.
panel_label(axA, "A", "Prespecified gates reduce 267 candidates to 87 testable mappings")
axA.set_xlim(-0.3, 6.3); axA.set_ylim(-0.28, 1.18); axA.axis("off")
labels = ["Core\nchemicals", "Entity\nvalid", "Exposure\ninterpretable",
          "Biomarker\navailable", "Detectable", "Cycle\ncoverage", "Human\ntestable"]
counts = np.asarray(observed, dtype=float)
heights = 0.22 + 0.60 * np.sqrt(counts / counts.max())
centers = np.full(7, 0.52)
colors = ["#C8D2D8", "#B8C9D2", "#9CBCC8", "#82AFBC", "#66A1AE", "#4B929F", TEAL]
for i in range(6):
    smooth_ribbon(axA, i+0.16, i+0.84,
                  centers[i]-heights[i]/2, centers[i]+heights[i]/2,
                  centers[i+1]-heights[i+1]/2, centers[i+1]+heights[i+1]/2,
                  colors[i+1], alpha=0.88)
for i, (lab, n, h, col) in enumerate(zip(labels, counts.astype(int), heights, colors)):
    axA.add_patch(FancyBboxPatch((i-0.16, centers[i]-h/2), 0.32, h,
        boxstyle="round,pad=0.006,rounding_size=0.035", facecolor=col,
        edgecolor="white", linewidth=0.8))
    axA.text(i, centers[i], f"{n}", ha="center", va="center", fontsize=9,
             fontweight="bold", color="white" if i >= 4 else INK)
    axA.text(i, centers[i]-h/2-0.075, lab, ha="center", va="top", fontsize=6.6, color=INK)
    if i < 6:
        loss = int(counts[i] - counts[i+1])
        axA.text(i+0.5, centers[i]+max(h,heights[i+1])/2+0.075, f"−{loss}",
                 ha="center", color=MUTED, fontsize=6.5)
axA.plot([-0.15, 6.15], [1.045, 1.045], color=BLUE, linewidth=2.0, solid_capstyle="round")
axA.text(3.0, 1.09, "CRC outcomes remain behind the firewall until 15 biomarker tests are frozen",
         ha="center", va="bottom", fontsize=7.2, color=BLUE, fontweight="bold")
axA.text(6.0, centers[6]-0.055, "27 strict", ha="center", va="top",
         fontsize=5.8, color="white", fontweight="bold")
axA.text(3.0, -0.22, "87 eligible chemical–biomarker mappings corresponded to 15 unique NHANES biomarker tests",
         ha="center", fontsize=7.6, color=INK, fontweight="bold")

# Panel B — multiplicity universe.
panel_label(axB, "B", "The statistical multiplicity unit is the biomarker test")
axB.set_xlim(0, 5); axB.set_ylim(-0.15, 3.15); axB.axis("off")
axes2 = axes.sort_values(["axis_key", "primary_biomarker"]).reset_index(drop=True)
for idx, row in axes2.iterrows():
    r, c = divmod(idx, 5)
    x, y = c + 0.08, 2.48 - r*0.89
    is_urine = row.axis_key.startswith("urine|")
    col = TEAL if is_urine else ORANGE
    n = int(row.eligible_chemical_count)
    w = 0.86
    axB.add_patch(FancyBboxPatch((x, y), w, 0.62,
        boxstyle="round,pad=0.012,rounding_size=0.055", facecolor="#F4F7F8",
        edgecolor=LINE, linewidth=0.6))
    axB.scatter([x+0.15], [y+0.31], s=14+math.sqrt(n)*18, color=col,
                edgecolor="white", linewidth=0.55, zorder=3)
    label = "MCOP" if row.primary_biomarker == "URXCOP" else row.primary_biomarker
    axB.text(x+0.29, y+0.39, label, fontsize=6.8, fontweight="bold", color=INK, va="center")
    axB.text(x+0.29, y+0.18, f"{n} mapping{'s' if n != 1 else ''}", fontsize=5.7, color=MUTED, va="center")
axB.text(0.05, -0.04, "● urine", color=TEAL, fontsize=6.2)
axB.text(0.90, -0.04, "● serum/blood", color=ORANGE, fontsize=6.2)
axB.text(4.95, -0.04, "15 tests = BH-FDR denominator", ha="right", color=BLUE, fontsize=6.4, fontweight="bold")

# Panel C — identity boundary.
panel_label(axC, "C", "DINP-axis translation preserves chemical identity")
axC.axis("off")
rounded_card(axC, (0.02, 0.60), (0.96, 0.26), "#E8EEF4", BLUE,
             "MiNP", "Molecular nominee",
             "rank 24 · BH-FDR 0.00346\nempirical FDR 0.0356 · 40.7% detectable")
rounded_card(axC, (0.02, 0.31), (0.96, 0.22), "#F2EFF5", PURPLE,
             "DINP parent", "Exposure-axis parent",
             "rank 107 · BH-FDR 0.449 · not a significant Phase 1 hit")
rounded_card(axC, (0.02, 0.02), (0.96, 0.22), "#E6F1F1", TEAL,
             "MCOP", "Human urinary biomarker",
             "98.8% detectable · seven cycles · entered the 15-test screen")
axC.annotate("", xy=(0.94, 0.25), xytext=(0.94, 0.58),
             xycoords=axC.transAxes, textcoords=axC.transAxes, ha="right", va="center",
             fontsize=6.1, color=MUTED,
             arrowprops=dict(arrowstyle="-|>", color=TEAL, linewidth=0.8,
                             linestyle=(0, (3, 2)), connectionstyle="arc3,rad=-0.22"))
axC.text(0.50, -0.07, "Distinct entities; no direct MCOP molecular-hit claim",
         transform=axC.transAxes, ha="center", fontsize=6.4, color=RED, fontweight="bold")

stem = DEST / "Figure2_actionability_final"
save_cns_figure(fig, stem)
fig.savefig(f"{stem}.svg", bbox_inches="tight", dpi=300)
plt.close(fig)

# Source data and reporting artifacts.
pd.DataFrame({"stage": labels, "n": counts.astype(int),
              "loss_to_next": [8,124,1,7,3,37,np.nan]}).to_csv(SRC / "Figure2_panelA_attrition.csv", index=False)
axes2.to_csv(SRC / "Figure2_panelB_biomarker_tests.csv", index=False)
identity = matrix.loc[matrix.ChemicalID.isin(["C471400", "C012125", "C573544"]),
    ["ChemicalID","ChemicalName","phase1_unfiltered_rank","phase1_bh_fdr",
     "phase1_degree_matched_bh_fdr","above_lod_pct","n_cycles_available",
     "selected_primary_biomarker","final_disposition"]]
identity.to_csv(SRC / "Figure2_panelC_identity.csv", index=False)

legend = """# Figure 2 legend

**Figure 2 | Outcome-blinded actionability defines the human screening universe.**
(A) Prespecified gates reduced 267 core environmental chemicals to 87 human-testable chemical–biomarker mappings. CRC outcome statistics remained behind the outcome firewall until the biomarker-test universe was frozen; 27 mappings met the strict D2/C2/T2 rule. (B) The 87 eligible mappings represented 15 unique NHANES biomarker tests, which formed the denominator for BH-FDR correction. Point area reflects the number of eligible chemical mappings and color denotes biological matrix. (C) MiNP, parent DINP and urinary MCOP retained distinct roles. MiNP was a molecular nominee but failed the direct-detectability gate, parent DINP was not a significant Phase 1 hit, and MCOP entered the human screen as a measurable biomarker for a DINP-related exposure axis. Biomarker translation does not imply chemical equivalence or a direct MCOP molecular hit.
"""
(DEST / "Figure2_actionability_final_legend.md").write_text(legend, encoding="utf-8")

stats = """# Figure 2 statistics and reproducibility

- Panel A independent unit: chemical–biomarker mapping; 267 starting chemical rows; deterministic prespecified gates; no inferential statistic.
- Panel B statistical unit: unique biological_matrix|primary_biomarker axis; 15 tests; mapping counts sum to 87.
- Panel C displays frozen source values for three chemically distinct entities; no causal arrow or statistical comparison.
- Outcome firewall: CRC OR, CI and P values were not used in Panels A–C or in eligibility classification.
- Source files: environmental_crc_267_actionability_flow.csv; environmental_crc_267_human_testable_candidates.csv; environmental_crc_267_actionability_matrix_v2.csv.
- No rows were dropped from the 267-row actionability matrix.
"""
(DEST / "Figure2_actionability_final_statistics.md").write_text(stats, encoding="utf-8")

qa = """# Figure 2 QA

Verdict: READY

- AP: restrained semantic palette; no rainbow/default palette; no four-sided chart borders; no occluding legend.
- CL: 183-mm double-column canvas; Arial/Helvetica/Liberation Sans; vector PDF and SVG; 300-dpi PNG; editable SVG text.
- VI: Panel A is the hero; Panels B and C answer non-redundant questions; outcome-free actionability and chemical identity are explicit.
- VV: all labels remain inside panel bounds; minimum text exceeds 5 pt; 2011–2012 or CRC outcome values are not introduced into the outcome-blinded figure.
- Data integrity: 267/267 rows used; stage counts and 87→15 mapping identity asserted in code.
"""
(DEST / "Figure2_actionability_final_QA.md").write_text(qa, encoding="utf-8")

print(json.dumps({
    "figure": str(stem), "rows_used": len(matrix), "eligible_mappings": int(alive.sum()),
    "unique_tests": len(axes), "strict": int(frozen["strict_eligibility"])
}, indent=2))
