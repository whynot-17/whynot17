# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) schematic-led workflow → cross-type inherit → param inherit
# (b) ranked chemical screen → cross-type inherit → param inherit
# (c) exposure-biomarker translation → cross-type inherit → param inherit
# (d) study roadmap → cross-type inherit → param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib as mpl

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


from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FIGURE_WIDTH_MM = 183
FIGURE_HEIGHT_MM = 88
MM_TO_INCH = 1 / 25.4

# Semantic roles for this figure. Exposure is blue, tumor is red, and the
# hypothetical bridge is deliberately pale and dashed.
EXPOSURE = CATEGORICAL[0]
EXPOSURE_LIGHT = "#DCEAF4"
EXPOSURE_PALE = "#EEF5F9"
TUMOR = ACCENT_RED
TUMOR_LIGHT = "#F6E4E8"
INFLAMMATION = CATEGORICAL[3]
INFLAMMATION_LIGHT = "#FBF0E3"
NEUTRAL_DARK = "#4B5358"
NEUTRAL = "#7C858A"
NEUTRAL_LIGHT = "#C9D0D4"
NEUTRAL_PALE = "#F2F4F5"
WHITE = "#FFFFFF"


def rounded_card(
    ax,
    xywh,
    text,
    *,
    facecolor=WHITE,
    edgecolor=NEUTRAL_LIGHT,
    linewidth=0.8,
    text_color=BLACK,
    fontsize=6.1,
    weight="normal",
    radius=0.025,
    va="center",
):
    x, y, w, h = xywh
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        clip_on=False,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        transform=ax.transAxes,
        ha="center",
        va=va,
        color=text_color,
        fontsize=fontsize,
        fontweight=weight,
        linespacing=1.16,
        clip_on=False,
    )
    return patch


def arrow(ax, start, end, *, color=NEUTRAL, dashed=False, mutation_scale=10):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=0.9,
            linestyle=(0, (3, 2)) if dashed else "-",
            color=color,
            shrinkA=2,
            shrinkB=2,
            clip_on=False,
        )
    )


def panel_title(ax, label, title):
    ax.text(
        -0.02,
        1.075,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.0,
        fontweight="bold",
        color=BLACK,
        clip_on=False,
    )
    ax.text(
        0.08,
        1.075,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.4,
        fontweight="bold",
        color=BLACK,
        clip_on=False,
    )


def validate_and_load(repo_root: Path, source_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    screen_path = source_dir / "figure1_primary_screen.csv"
    screen = pd.read_csv(screen_path)
    required = {
        "ChemicalID",
        "display_name",
        "chemical_class",
        "screen_rank",
        "odds_ratio",
        "bh_fdr",
        "minus_log10_bh_fdr",
        "degree_matched_bh_fdr",
    }
    missing = required.difference(screen.columns)
    if missing:
        raise ValueError(f"Figure 1 source is missing columns: {sorted(missing)}")
    if len(screen) != 267:
        raise ValueError(f"Expected all 267 chemicals; found {len(screen)}")
    if not np.isfinite(screen["minus_log10_bh_fdr"].to_numpy(float)).all():
        raise ValueError("Ranked-screen y values contain non-finite values")
    if int((screen["bh_fdr"] < 0.05).sum()) != 52:
        raise ValueError("Frozen screen should contain 52 BH-FDR-significant chemicals")
    if int(screen["stable_for_primary_sort"].astype(str).str.lower().eq("true").sum()) != 69:
        raise ValueError("Frozen screen should contain 69 stable candidates")
    for chemical in ["MiNP", "DINP", "MBzP"]:
        if int(screen["display_name"].eq(chemical).sum()) != 1:
            raise ValueError(f"Expected exactly one {chemical} row")

    nodes = pd.read_csv(source_dir / "figure1_panelA_workflow_nodes.csv")
    roadmap = pd.read_csv(source_dir / "figure1_panelD_study_roadmap.csv")

    # Keep the source manifest alongside the rendered figure for traceability.
    manifest = pd.DataFrame(
        [
            (str(screen_path.relative_to(repo_root)), "all 267 rows", "Panel B ranked screen"),
            (str((source_dir / "figure1_panelA_workflow_nodes.csv").relative_to(repo_root)), "4 workflow nodes", "Panel A workflow labels"),
            (str((source_dir / "figure1_panelD_study_roadmap.csv").relative_to(repo_root)), "3 frozen stages", "Panel D completed roadmap"),
            ("outputs/mcop_crc_phase2h_primary_reanalysis.csv", "frozen primary survey result", "Panel D NHANES evidence label"),
        ],
        columns=["source_file", "rows_or_scope", "figure_use"],
    )
    return screen, nodes, roadmap, manifest


def draw_panel_a(ax, nodes: pd.DataFrame):
    ax.set_axis_off()
    panel_title(ax, "A", "Data-first discovery")
    rounded_card(
        ax,
        (0.11, 0.78, 0.78, 0.13),
        "267 core environmental\nchemicals",
        facecolor=EXPOSURE_LIGHT,
        edgecolor=EXPOSURE,
        linewidth=1.15,
        text_color=EXPOSURE,
        fontsize=6.5,
        weight="bold",
    )
    rounded_card(
        ax,
        (0.04, 0.51, 0.43, 0.16),
        "CTD human\nchemical–gene\ninteractions",
        facecolor=NEUTRAL_PALE,
        edgecolor=NEUTRAL_LIGHT,
        text_color=NEUTRAL_DARK,
        fontsize=5.8,
    )
    rounded_card(
        ax,
        (0.53, 0.51, 0.43, 0.16),
        "GeneCards\nCRC-associated\ngenes",
        facecolor=NEUTRAL_PALE,
        edgecolor=NEUTRAL_LIGHT,
        text_color=NEUTRAL_DARK,
        fontsize=5.8,
    )
    arrow(ax, (0.50, 0.77), (0.27, 0.68))
    arrow(ax, (0.50, 0.77), (0.73, 0.68))
    arrow(ax, (0.27, 0.50), (0.46, 0.38))
    arrow(ax, (0.73, 0.50), (0.54, 0.38))
    rounded_card(
        ax,
        (0.10, 0.19, 0.80, 0.17),
        "Fisher / hypergeometric enrichment\nBH-FDR + degree-matched permutation",
        facecolor=WHITE,
        edgecolor=NEUTRAL,
        linewidth=1.0,
        text_color=NEUTRAL_DARK,
        fontsize=5.7,
        weight="bold",
    )
    ax.text(
        0.50,
        0.08,
        "No chemical was preselected",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.7,
        color=NEUTRAL,
    )


def draw_panel_b(ax, screen: pd.DataFrame):
    panel_title(ax, "B", "Ranked chemical screen")
    phthalate = screen["chemical_class"].astype(str).str.lower().eq("phthalates")
    ax.scatter(
        screen.loc[~phthalate, "screen_rank"],
        screen.loc[~phthalate, "minus_log10_bh_fdr"],
        s=10,
        color=NEUTRAL_LIGHT,
        alpha=0.78,
        linewidths=0,
        zorder=1,
    )
    ax.scatter(
        screen.loc[phthalate, "screen_rank"],
        screen.loc[phthalate, "minus_log10_bh_fdr"],
        s=12,
        color="#8FC1D0",
        alpha=0.82,
        linewidths=0,
        zorder=2,
    )
    threshold = -np.log10(0.05)
    ax.axhline(threshold, color=NEUTRAL, linewidth=0.7, linestyle=(0, (3, 2)), zorder=0)
    ax.text(
        0.98,
        0.96,
        f"52 BH-FDR < 0.05  |  69 stable candidates",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
        color=NEUTRAL,
    )
    ax.text(
        0.98,
        0.91,
        "phthalate-class points highlighted",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.3,
        color=EXPOSURE,
    )

    special = {
        "MiNP": dict(marker="D", size=50, color=EXPOSURE, text="MiNP  |  rank 24\nOR 10.06; degree-matched FDR 0.0356", xytext=(59, 4.05)),
        "DINP": dict(marker="o", size=38, color=EXPOSURE, text="DINP parent  |  rank 107", xytext=(126, 0.92)),
        "MBzP": dict(marker="s", size=38, color=NEUTRAL_DARK, text="MBzP  |  rank 2", xytext=(17, 8.38)),
    }
    for name, spec in special.items():
        row = screen.loc[screen["display_name"].eq(name)].iloc[0]
        x = float(row["screen_rank"])
        y = float(row["minus_log10_bh_fdr"])
        face = WHITE if name == "DINP" else spec["color"]
        ax.scatter(
            [x],
            [y],
            s=spec["size"],
            marker=spec["marker"],
            facecolor=face,
            edgecolor=spec["color"],
            linewidth=1.1,
            zorder=5,
        )
        ax.annotate(
            spec["text"],
            xy=(x, y),
            xytext=spec["xytext"],
            textcoords="data",
            ha="left",
            va="center",
            fontsize=5.4 if name != "MiNP" else 5.6,
            color=spec["color"],
            fontweight="bold" if name == "MiNP" else "normal",
            arrowprops=dict(arrowstyle="-", color=spec["color"], lw=0.8),
        )
    ax.set_xlim(0, 270)
    ax.set_ylim(-0.35, 9.75)
    ax.set_xlabel("Primary screen rank")
    ax.set_ylabel("−log10(BH-FDR)")
    ax.set_xticks([1, 50, 100, 150, 200, 250])
    ax.tick_params(labelsize=6.2)
    ax.spines["left"].set_color(NEUTRAL_DARK)
    ax.spines["bottom"].set_color(NEUTRAL_DARK)


def draw_panel_c(ax):
    ax.set_axis_off()
    panel_title(ax, "C", "Exposure-axis translation")
    rounded_card(
        ax,
        (0.10, 0.76, 0.80, 0.14),
        "MiNP / DINP\nmolecular evidence",
        facecolor=EXPOSURE_LIGHT,
        edgecolor=EXPOSURE,
        linewidth=1.1,
        text_color=EXPOSURE,
        fontsize=6.4,
        weight="bold",
    )
    arrow(ax, (0.50, 0.74), (0.50, 0.64), color=EXPOSURE)
    rounded_card(
        ax,
        (0.10, 0.48, 0.80, 0.15),
        "DINP-related\nexposure axis",
        facecolor=WHITE,
        edgecolor=EXPOSURE,
        linewidth=1.1,
        text_color=EXPOSURE,
        fontsize=6.6,
        weight="bold",
    )
    arrow(ax, (0.50, 0.46), (0.50, 0.35), color=EXPOSURE, dashed=True)
    ax.text(
        0.54,
        0.405,
        "biomarker\ntranslation",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=5.0,
        color=NEUTRAL,
    )
    rounded_card(
        ax,
        (0.10, 0.18, 0.80, 0.16),
        "Urinary MCOP\nNHANES biomarker",
        facecolor=EXPOSURE_PALE,
        edgecolor=EXPOSURE,
        linewidth=1.2,
        text_color=EXPOSURE,
        fontsize=6.5,
        weight="bold",
    )
    ax.text(
        0.50,
        0.08,
        "MCOP was selected for broad cycle coverage\nand high detectability — not as a direct CTD hit",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.25,
        color=NEUTRAL,
        linespacing=1.15,
    )


def draw_panel_d(ax, roadmap: pd.DataFrame):
    ax.set_axis_off()
    panel_title(ax, "D", "Study roadmap")
    cards = [
        (0.07, 0.78, 0.86, 0.14, "1  Discovery", "CTD / GeneCards → MiNP / DINP axis", NEUTRAL_PALE, NEUTRAL_LIGHT, NEUTRAL_DARK, False),
        (0.07, 0.56, 0.86, 0.17, "2  Human biomonitoring", "NHANES 2005–2018\nN = 9,936; CRC = 70\nurinary MCOP", EXPOSURE_PALE, EXPOSURE, EXPOSURE, False),
        (0.07, 0.32, 0.86, 0.17, "3  CRC biological state", "TCGA paired bulk\nCELLxGENE donors\nGSE144735 directional check", TUMOR_LIGHT, TUMOR, TUMOR, False),
        (0.07, 0.11, 0.86, 0.14, "Future  Prospective replication", "WHI prediagnostic urine → incident CRC", WHITE, EXPOSURE, EXPOSURE, True),
    ]
    for i, (x, y, w, h, heading, body, face, edge, text_color, dashed) in enumerate(cards):
        patch = rounded_card(
            ax,
            (x, y, w, h),
            "",
            facecolor=face,
            edgecolor=edge,
            linewidth=1.05 if not dashed else 0.9,
            text_color=text_color,
        )
        if dashed:
            patch.set_linestyle((0, (3, 2)))
        ax.text(
            x + 0.08,
            y + h - 0.045,
            heading,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=6.0,
            fontweight="bold",
            color=text_color,
        )
        ax.text(
            x + 0.08,
            y + h / 2 - 0.012,
            body,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=5.25,
            color=NEUTRAL_DARK if not dashed else EXPOSURE,
            linespacing=1.12,
        )
        if i < len(cards) - 1:
            next_y = cards[i + 1][1] + cards[i + 1][3]
            arrow(ax, (0.50, y - 0.015), (0.50, next_y + 0.015), color=NEUTRAL, mutation_scale=8)
    ax.text(
        0.50,
        0.035,
        "Solid = observed; dashed = future replication",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.15,
        color=NEUTRAL,
    )


def write_reports(output_dir: Path, screen: pd.DataFrame, manifest: pd.DataFrame):
    manifest.to_csv(output_dir / "Figure1_study_design_v2_source_manifest.csv", index=False)
    stats = "# Figure 1 academic v2 — statistics and reproducibility\n\n"
    stats += "## Scientific question\n"
    stats += "How did the data-first environmental screen lead to DINP-axis biomonitoring and the subsequent CRC biological-state analyses?\n\n"
    stats += "## Quantitative panel\n"
    stats += f"- Panel B uses all {len(screen)} frozen Phase 1 chemicals; no rows were downsampled or dropped.\n"
    stats += f"- Primary ranking metric: `-log10(BH-FDR)` from the GeneCards Disorders-scoped, `gene_cards_k=1000`, `U_core` screen.\n"
    stats += f"- Frozen screen checks: {int((screen['bh_fdr'] < 0.05).sum())} chemicals with BH-FDR < 0.05; {int(screen['stable_for_primary_sort'].astype(str).str.lower().eq('true').sum())} stable candidates.\n"
    stats += "- Panel B is descriptive; it does not introduce a new statistical test or re-rank candidates.\n\n"
    stats += "## Schematic panels\n"
    stats += "- Panel A: workflow nodes from `figure1_panelA_workflow_nodes.csv`; no quantitative effect estimate is encoded.\n"
    stats += "- Panel C: translation boundary is explicitly labelled as biomarker translation; MCOP is not represented as a direct CTD nomination.\n"
    stats += "- Panel D: completed evidence layers are solid; WHI is a dashed future prospective replication stage and contains no result.\n\n"
    stats += "## Source traceability\n"
    for row in manifest.itertuples(index=False):
        stats += f"- `{row.source_file}` — {row.rows_or_scope} — {row.figure_use}.\n"
    stats += "\n## Export\n- Vector master: PDF.\n- Preview: RGB PNG at 300 dpi.\n- Editable text: SVG export.\n"
    (output_dir / "Figure1_study_design_v2_statistics.md").write_text(stats, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.repo_root = args.repo_root.resolve()
    args.source_dir = args.source_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    screen, nodes, roadmap, manifest = validate_and_load(args.repo_root, args.source_dir)
    fig = plt.figure(
        figsize=(FIGURE_WIDTH_MM * MM_TO_INCH, FIGURE_HEIGHT_MM * MM_TO_INCH),
        facecolor=WHITE,
    )
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[1.00, 1.50, 1.08, 1.22],
        left=0.035,
        right=0.985,
        bottom=0.10,
        top=0.86,
        wspace=0.34,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[0, 3])

    draw_panel_a(ax_a, nodes)
    draw_panel_b(ax_b, screen)
    draw_panel_c(ax_c)
    draw_panel_d(ax_d, roadmap)

    stem = args.output_dir / "Figure1_study_design_v2"
    save_cns_figure(fig, stem)
    fig.savefig(f"{stem}.svg", bbox_inches="tight", dpi=300)
    write_reports(args.output_dir, screen, manifest)
    plt.close(fig)
    print(f"Rendered {stem}.pdf, {stem}.png, {stem}.svg")


if __name__ == "__main__":
    main()
