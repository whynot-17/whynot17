from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import FuncFormatter, FixedLocator, NullLocator
from matplotlib.transforms import Bbox


FIGURE_WIDTH_IN = 190 / 25.4

COLORS = {
    "mcop": "#1F6F78",
    "mcop_light": "#79AEB3",
    "mcop_pale": "#DDEBED",
    "tumor": "#8E2C43",
    "tumor_light": "#CF8A99",
    "inflammation": "#C46543",
    "normal": "#A7ADB2",
    "normal_light": "#D6D9DC",
    "neutral": "#6E7478",
    "neutral_dark": "#303437",
    "neutral_pale": "#F2F3F4",
    "warning": "#A26A21",
    "white": "#FFFFFF",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7.2,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.75,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column in out:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def format_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.2e}"
    if value < 0.01:
        return f"{value:.4f}"
    return f"{value:.3f}"


def short_cycle(value: str) -> str:
    left, right = value.split("-")
    return f"{left}\u2013{right[-2:]}"


def set_or_axis(ax: plt.Axes, limits: tuple[float, float], ticks: list[float]) -> None:
    if limits[0] <= 0 or limits[1] <= 0 or any(tick <= 0 for tick in ticks):
        raise ValueError("Log-scale odds-ratio axes require strictly positive limits and ticks")
    ax.set_xscale("log")
    ax.set_xlim(*limits)
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.axvline(1.0, color=COLORS["neutral"], lw=0.8, ls=(0, (3, 2)), zorder=0)
    ax.grid(False)


def add_panel_label(fig: plt.Figure, ax: plt.Axes, label: str, dx: float = 0.032) -> None:
    pos = ax.get_position()
    fig.text(
        max(0.005, pos.x0 - dx),
        min(0.975, pos.y1 + 0.012),
        label,
        fontsize=9.0,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=COLORS["neutral_dark"],
    )


def export_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.svg")
    fig.savefig(output_dir / f"{stem}.pdf")
    fig.savefig(output_dir / f"{stem}.png", dpi=600)
    plt.close(fig)


def export_panel(
    fig: plt.Figure,
    axes: plt.Axes | list[plt.Axes] | tuple[plt.Axes, ...],
    output_dir: Path,
    stem: str,
) -> None:
    panel_axes = [axes] if isinstance(axes, plt.Axes) else list(axes)
    hidden_axes = [axis for axis in fig.axes if axis not in panel_axes and axis.get_visible()]
    try:
        for axis in hidden_axes:
            axis.set_visible(False)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        bbox_px = Bbox.union([axis.get_tightbbox(renderer) for axis in panel_axes])
        bbox_in = bbox_px.transformed(fig.dpi_scale_trans.inverted())
        padded = Bbox.from_extents(
            max(0.0, bbox_in.x0 - 0.08),
            max(0.0, bbox_in.y0 - 0.08),
            min(fig.get_figwidth(), bbox_in.x1 + 0.08),
            min(fig.get_figheight(), bbox_in.y1 + 0.24),
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{stem}.svg", bbox_inches=padded)
        fig.savefig(output_dir / f"{stem}.pdf", bbox_inches=padded)
        fig.savefig(output_dir / f"{stem}.png", dpi=600, bbox_inches=padded)
    finally:
        for axis in hidden_axes:
            axis.set_visible(True)


def rounded_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    facecolor: str,
    edgecolor: str,
    linewidth: float = 0.9,
    radius: float = 0.025,
    linestyle: str | tuple = "-",
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        linestyle=linestyle,
        clip_on=False,
    )
    ax.add_patch(patch)
    return patch


def axes_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = COLORS["neutral"],
    linewidth: float = 1.0,
    linestyle: str | tuple = "-",
    mutation_scale: float = 8.0,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=linewidth,
            linestyle=linestyle,
            color=color,
            shrinkA=2.5,
            shrinkB=2.5,
            clip_on=False,
        )
    )


def forest_rows(
    ax: plt.Axes,
    rows: pd.DataFrame,
    y_positions: np.ndarray,
    colors: list[str],
    markers: list[str] | None = None,
    sizes: list[float] | None = None,
    open_markers: list[bool] | None = None,
) -> None:
    markers = markers or ["o"] * len(rows)
    sizes = sizes or [4.4] * len(rows)
    open_markers = open_markers or [False] * len(rows)
    for (_, row), y, color, marker, size, is_open in zip(
        rows.iterrows(), y_positions, colors, markers, sizes, open_markers
    ):
        face = COLORS["white"] if is_open else color
        ax.errorbar(
            row["OR"],
            y,
            xerr=[[row["OR"] - row["CI_low"]], [row["CI_high"] - row["OR"]]],
            fmt=marker,
            ms=size,
            mfc=face,
            mec=color,
            mew=0.9,
            ecolor=color,
            elinewidth=1.05,
            capsize=1.8,
            capthick=0.8,
            zorder=3,
        )


def paired_plot(
    ax: plt.Axes,
    normal: np.ndarray,
    tumor: np.ndarray,
    ylabel: str | None,
    title: str,
    point_size: float = 12,
    title_size: float = 8.0,
) -> None:
    normal = np.asarray(normal, dtype=float)
    tumor = np.asarray(tumor, dtype=float)
    if len(normal) != len(tumor):
        raise ValueError("Paired arrays must have equal length")
    for left, right in zip(normal, tumor):
        ax.plot([0, 1], [left, right], color=COLORS["normal_light"], lw=0.65, alpha=0.72, zorder=1)
    ax.scatter(
        np.zeros_like(normal),
        normal,
        s=point_size,
        facecolor=COLORS["normal"],
        edgecolor=COLORS["white"],
        linewidth=0.35,
        zorder=3,
    )
    ax.scatter(
        np.ones_like(tumor),
        tumor,
        s=point_size,
        facecolor=COLORS["tumor"],
        edgecolor=COLORS["white"],
        linewidth=0.35,
        zorder=3,
    )
    for x, values, color in [(0, normal, COLORS["neutral_dark"]), (1, tumor, COLORS["tumor"] )]:
        median = float(np.median(values))
        ax.plot([x - 0.16, x + 0.16], [median, median], color=color, lw=1.8, zorder=4)
    combined = np.r_[normal, tumor]
    spread = float(np.nanmax(combined) - np.nanmin(combined))
    pad = 0.10 * spread if spread > 0 else 0.5
    ax.set_ylim(float(np.nanmin(combined) - pad), float(np.nanmax(combined) + pad))
    ax.set_xlim(-0.28, 1.28)
    ax.set_xticks([0, 1], ["Normal", "Tumor"])
    ax.set_ylabel(ylabel or "")
    ax.set_title(
        title,
        loc="left",
        pad=5.0,
        fontweight="bold",
        fontsize=title_size,
    )
    ax.grid(False)


def prepare_figure_1_data(repo_root: Path, source_dir: Path) -> pd.DataFrame:
    raw = numeric(
        pd.read_csv(repo_root / "outputs" / "environmental_toxicology_crc_phase1_ranked_core.csv"),
        [
            "gene_cards_k",
            "odds_ratio",
            "bh_fdr",
            "fisher_p",
            "crc_overlap",
            "n_ctd_human_genes",
        ],
    )
    mask = (
        raw["scope"].eq("GeneCards_Disorders")
        & raw["gene_cards_k"].eq(1000)
        & raw["background"].eq("U_core")
    )
    screen = raw.loc[mask].copy()
    if len(screen) != 267:
        raise ValueError(f"Expected 267 chemicals in the frozen primary screen, found {len(screen)}")
    if int((screen["bh_fdr"] < 0.05).sum()) != 52:
        raise ValueError("Frozen primary screen must contain 52 BH-FDR-significant chemicals")
    if int(screen["stable_for_primary_sort"].astype(str).str.lower().eq("true").sum()) != 69:
        raise ValueError("Frozen primary screen must contain 69 stable candidates")

    permutation = numeric(
        pd.read_csv(
            repo_root
            / "outputs"
            / "environmental_toxicology_crc_phase1_degree_matched_permutation.csv"
        ),
        ["degree_matched_empirical_p", "degree_matched_bh_fdr"],
    )
    screen = screen.merge(
        permutation[
            ["ChemicalID", "degree_matched_empirical_p", "degree_matched_bh_fdr"]
        ],
        on="ChemicalID",
        how="left",
        validate="one_to_one",
    )
    screen = screen.sort_values(["bh_fdr", "fisher_p", "ChemicalID"]).reset_index(drop=True)
    screen["screen_rank"] = np.arange(1, len(screen) + 1)
    screen["minus_log10_bh_fdr"] = -np.log10(screen["bh_fdr"].clip(lower=1e-300))
    screen["display_name"] = screen["ChemicalName"]
    screen.loc[screen["ChemicalID"].eq("C471400"), "display_name"] = "MiNP"
    screen.loc[screen["ChemicalID"].eq("C012125"), "display_name"] = "DINP"
    screen.loc[screen["ChemicalID"].eq("C103325"), "display_name"] = "MBzP"

    minp = screen.loc[screen["ChemicalID"].eq("C471400")]
    if len(minp) != 1:
        raise ValueError("Frozen primary screen must contain exactly one MiNP row")
    minp_row = minp.iloc[0]
    if not np.isclose(float(minp_row["odds_ratio"]), 10.06442831215971):
        raise ValueError("MiNP odds ratio differs from the manuscript lock")
    if not np.isclose(float(minp_row["degree_matched_bh_fdr"]), 0.035564435564435566):
        raise ValueError("MiNP degree-matched FDR differs from the manuscript lock")

    source_dir.mkdir(parents=True, exist_ok=True)
    screen[
        [
            "ChemicalID",
            "ChemicalName",
            "display_name",
            "chemical_class",
            "screen_rank",
            "odds_ratio",
            "bh_fdr",
            "minus_log10_bh_fdr",
            "crc_overlap",
            "n_ctd_human_genes",
            "stable_for_primary_sort",
            "degree_matched_empirical_p",
            "degree_matched_bh_fdr",
        ]
    ].to_csv(source_dir / "figure1_primary_screen.csv", index=False)

    pd.DataFrame(
        [
            ("chemical_universe", "267 environmental chemicals", "input"),
            ("ctd", "CTD human chemical-gene interactions", "evidence"),
            ("genecards", "GeneCards CRC genes", "evidence"),
            ("audit", "Enrichment + degree-matched permutation", "analysis"),
        ],
        columns=["node_id", "label", "node_type"],
    ).to_csv(source_dir / "figure1_panelA_workflow_nodes.csv", index=False)
    pd.DataFrame(
        [
            ("MiNP", "DINP-related exposure axis", "molecular nomination"),
            ("DINP-related exposure axis", "MCOP", "biomarker translation"),
        ],
        columns=["source", "target", "relation"],
    ).to_csv(source_dir / "figure1_panelC_translation_links.csv", index=False)
    pd.DataFrame(
        [
            (1, "Discovery", "CTD/GeneCards; MiNP/DINP axis"),
            (2, "Human biomonitoring", "NHANES 2005-2018; N=9,936; CRC=70; urinary MCOP"),
            (3, "Biological interpretation", "TCGA; CELLxGENE; GSE144735"),
        ],
        columns=["stage", "stage_label", "content"],
    ).to_csv(source_dir / "figure1_panelD_study_roadmap.csv", index=False)
    return screen


def build_figure_1(
    repo_root: Path,
    source_dir: Path,
    output_dir: Path,
    panel_dir: Path,
) -> None:
    screen = prepare_figure_1_data(repo_root, source_dir)
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, 3.75))
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[1.00, 1.42, 1.10, 1.28],
        wspace=0.34,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[0, 3])
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.14, top=0.90)

    # A — discovery universe and audited workflow
    ax_a.set_axis_off()
    ax_a.set_title("Data-first discovery", loc="left", pad=5, fontweight="bold", fontsize=7.0)
    rounded_box(
        ax_a,
        0.11,
        0.78,
        0.78,
        0.12,
        facecolor=COLORS["mcop_pale"],
        edgecolor=COLORS["mcop"],
    )
    ax_a.text(
        0.50,
        0.84,
        "267 environmental\nchemicals",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=6.3,
        fontweight="bold",
        color=COLORS["mcop"],
    )
    rounded_box(
        ax_a,
        0.05,
        0.50,
        0.42,
        0.16,
        facecolor=COLORS["neutral_pale"],
        edgecolor=COLORS["normal"],
    )
    rounded_box(
        ax_a,
        0.53,
        0.50,
        0.42,
        0.16,
        facecolor=COLORS["neutral_pale"],
        edgecolor=COLORS["normal"],
    )
    ax_a.text(
        0.26,
        0.58,
        "CTD human\nchemical-gene\ninteractions",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=5.8,
        color=COLORS["neutral_dark"],
    )
    ax_a.text(
        0.74,
        0.58,
        "GeneCards\nCRC-associated\ngenes",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=5.8,
        color=COLORS["neutral_dark"],
    )
    axes_arrow(ax_a, (0.50, 0.78), (0.27, 0.67))
    axes_arrow(ax_a, (0.50, 0.78), (0.73, 0.67))
    axes_arrow(ax_a, (0.27, 0.49), (0.46, 0.36))
    axes_arrow(ax_a, (0.73, 0.49), (0.54, 0.36))
    rounded_box(
        ax_a,
        0.12,
        0.18,
        0.76,
        0.17,
        facecolor=COLORS["white"],
        edgecolor=COLORS["neutral"],
    )
    ax_a.text(
        0.50,
        0.265,
        "Enrichment + BH-FDR\nDegree-matched\npermutation",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=5.5,
        fontweight="bold",
        color=COLORS["neutral_dark"],
    )
    ax_a.text(
        0.50,
        0.08,
        "No chemical was preselected",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=5.5,
        color=COLORS["neutral"],
    )

    # B — all 267 chemicals, with transparent contextual anchors
    phthalate_mask = screen["chemical_class"].eq("phthalates")
    ax_b.scatter(
        screen.loc[~phthalate_mask, "screen_rank"],
        screen.loc[~phthalate_mask, "minus_log10_bh_fdr"],
        s=8,
        facecolor=COLORS["normal_light"],
        edgecolor="none",
        alpha=0.82,
        zorder=1,
    )
    ax_b.scatter(
        screen.loc[phthalate_mask, "screen_rank"],
        screen.loc[phthalate_mask, "minus_log10_bh_fdr"],
        s=10,
        facecolor=COLORS["mcop_light"],
        edgecolor="none",
        alpha=0.78,
        zorder=2,
    )
    ax_b.axhline(
        -np.log10(0.05),
        color=COLORS["neutral"],
        lw=0.7,
        ls=(0, (3, 2)),
        zorder=0,
    )
    special = {
        "C471400": ("MiNP  rank 24", COLORS["mcop"], "D", (13, 11)),
        "C012125": ("DINP parent", COLORS["mcop"], "o", (8, 8)),
        "C103325": ("MBzP", COLORS["neutral_dark"], "s", (7, -13)),
    }
    for chemical_id, (label, color, marker, offset) in special.items():
        row = screen.loc[screen["ChemicalID"].eq(chemical_id)].iloc[0]
        ax_b.scatter(
            row["screen_rank"],
            row["minus_log10_bh_fdr"],
            s=32 if chemical_id == "C471400" else 24,
            marker=marker,
            facecolor=COLORS["white"] if chemical_id == "C012125" else color,
            edgecolor=color,
            linewidth=0.9,
            zorder=4,
        )
        ax_b.annotate(
            label,
            (row["screen_rank"], row["minus_log10_bh_fdr"]),
            xytext=offset,
            textcoords="offset points",
            ha="left",
            va="bottom" if offset[1] >= 0 else "top",
            fontsize=5.5,
            fontweight="bold" if chemical_id == "C471400" else "normal",
            color=color,
            arrowprops={"arrowstyle": "-", "lw": 0.55, "color": color},
            zorder=5,
        )
    ax_b.set_xlim(-5, 276)
    ymax = float(screen["minus_log10_bh_fdr"].max())
    ax_b.set_ylim(-0.35, ymax + 0.65)
    ax_b.set_xticks([1, 50, 100, 150, 200, 250])
    ax_b.set_xlabel("Primary screen rank")
    ax_b.set_ylabel("-log10(BH-FDR)")
    ax_b.set_title("Ranked chemical screen", loc="left", pad=5, fontweight="bold", fontsize=7.0)
    ax_b.text(
        0.98,
        0.97,
        "52 FDR < 0.05\n69 stable candidates",
        transform=ax_b.transAxes,
        ha="right",
        va="top",
        fontsize=5.5,
        color=COLORS["neutral"],
    )
    ax_b.text(
        0.98,
        0.74,
        "MiNP: OR 10.06; FDR 0.00346\nDegree-matched FDR 0.0356",
        transform=ax_b.transAxes,
        ha="right",
        va="top",
        fontsize=5.3,
        color=COLORS["mcop"],
    )

    # C — translation boundary
    ax_c.set_axis_off()
    ax_c.set_title("MiNP/DINP to MCOP", loc="left", pad=5, fontweight="bold", fontsize=7.0)
    nodes = [
        (0.15, 0.72, 0.70, 0.16, "MiNP molecular evidence", COLORS["mcop_pale"], COLORS["mcop"]),
        (0.15, 0.45, 0.70, 0.16, "DINP-related\nexposure axis", COLORS["white"], COLORS["mcop"]),
        (0.15, 0.18, 0.70, 0.16, "Urinary MCOP\nNHANES biomarker", COLORS["mcop_pale"], COLORS["mcop"]),
    ]
    for x, y, width, height, label, face, edge in nodes:
        rounded_box(ax_c, x, y, width, height, facecolor=face, edgecolor=edge)
        ax_c.text(
            x + width / 2,
            y + height / 2,
            label,
            transform=ax_c.transAxes,
            ha="center",
            va="center",
            fontsize=6.2,
            fontweight="bold",
            color=COLORS["mcop"],
        )
    axes_arrow(ax_c, (0.50, 0.71), (0.50, 0.62), color=COLORS["mcop"])
    axes_arrow(ax_c, (0.50, 0.44), (0.50, 0.35), color=COLORS["mcop"])
    ax_c.text(
        0.50,
        0.08,
        "Broad cycle availability\nand high detectability",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=5.4,
        color=COLORS["neutral"],
    )

    # D — manuscript roadmap
    ax_d.set_axis_off()
    ax_d.set_title("Study roadmap", loc="left", pad=5, fontweight="bold", fontsize=7.0)
    roadmap = [
        (0.08, 0.70, "1  Discovery", "CTD/GeneCards\nMiNP/DINP axis", COLORS["neutral_pale"], COLORS["normal"]),
        (0.08, 0.41, "2  Human biomonitoring", "NHANES 2005-2018\nN = 9,936; CRC = 70\nurinary MCOP", COLORS["mcop_pale"], COLORS["mcop"]),
        (0.08, 0.10, "3  Biological interpretation", "TCGA paired bulk\nCELLxGENE donors\nGSE144735", "#F5E8EC", COLORS["tumor"]),
    ]
    for x, y, header, body, face, edge in roadmap:
        rounded_box(ax_d, x, y, 0.84, 0.19, facecolor=face, edgecolor=edge)
        ax_d.text(
            x + 0.04,
            y + 0.145,
            header,
            transform=ax_d.transAxes,
            ha="left",
            va="center",
            fontsize=5.9,
            fontweight="bold",
            color=edge,
        )
        ax_d.text(
            x + 0.04,
            y + 0.073,
            body,
            transform=ax_d.transAxes,
            ha="left",
            va="center",
            fontsize=5.4,
            color=COLORS["neutral_dark"],
        )
    axes_arrow(ax_d, (0.50, 0.69), (0.50, 0.61), color=COLORS["neutral"])
    axes_arrow(ax_d, (0.50, 0.40), (0.50, 0.31), color=COLORS["neutral"])
    ax_d.text(
        0.50,
        0.02,
        "Candidate bridge - not mechanism proven",
        transform=ax_d.transAxes,
        ha="center",
        va="bottom",
        fontsize=5.1,
        color=COLORS["neutral"],
    )

    for label, axis in zip("ABCD", [ax_a, ax_b, ax_c, ax_d]):
        add_panel_label(fig, axis, label, dx=0.027)
    export_panel(fig, ax_a, panel_dir, "Figure1_panelA_discovery_universe_v1")
    export_panel(fig, ax_b, panel_dir, "Figure1_panelB_ranked_screen_v1")
    export_panel(fig, ax_c, panel_dir, "Figure1_panelC_biomarker_translation_v1")
    export_panel(fig, ax_d, panel_dir, "Figure1_panelD_study_roadmap_v1")
    export_figure(fig, output_dir, "Figure1_study_design_v1")


def build_figure_2(source_dir: Path, output_dir: Path, panel_dir: Path) -> None:
    compare = numeric(
        pd.read_csv(source_dir / "figure2_primary_python_vs_r.csv"),
        [
            "python_OR",
            "python_CI_low",
            "python_CI_high",
            "python_P",
            "r_OR",
            "r_CI_low",
            "r_CI_high",
            "r_P_standard",
            "r_N",
            "r_CRC_N",
            "relative_logOR_change_pct",
        ],
    ).iloc[0]
    loco = numeric(
        pd.read_csv(source_dir / "figure2_loco.csv"),
        ["OR", "CI_low", "CI_high", "P", "N", "CRC_N"],
    )
    per_cycle = numeric(
        pd.read_csv(source_dir / "figure2_per_cycle.csv"),
        ["OR", "CI_low", "CI_high", "P", "N", "CRC_N"],
    )
    interaction = numeric(
        pd.read_csv(source_dir / "figure2_cycle_interaction.csv"),
        ["P_F", "N", "CRC_N"],
    ).iloc[0]

    assert len(loco) == 7 and len(per_cycle) == 7
    assert int(compare["r_N"]) == 9936 and int(compare["r_CRC_N"]) == 70
    assert np.all(loco["CI_low"] > 1.0)

    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, 5.35))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.86, 1.70],
        width_ratios=[1.30, 1.00],
        hspace=0.42,
        wspace=0.45,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    fig.subplots_adjust(left=0.125, right=0.975, bottom=0.095, top=0.94)

    # A — hero estimate
    set_or_axis(ax_a, (0.95, 1.52), [1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
    ax_a.errorbar(
        compare["r_OR"],
        0.34,
        xerr=[
            [compare["r_OR"] - compare["r_CI_low"]],
            [compare["r_CI_high"] - compare["r_OR"]],
        ],
        fmt="D",
        ms=6.2,
        mfc=COLORS["mcop"],
        mec=COLORS["mcop"],
        ecolor=COLORS["mcop"],
        elinewidth=1.6,
        capsize=2.6,
        zorder=3,
    )
    ax_a.set_ylim(0, 1)
    ax_a.set_yticks([])
    ax_a.set_xlabel("Odds ratio per doubling of MCOP")
    ax_a.set_title("Primary complex-survey estimate", loc="left", pad=5, fontweight="bold")
    ax_a.text(
        0.02,
        0.83,
        f"N = {int(compare['r_N']):,}  |  CRC cases = {int(compare['r_CRC_N'])}",
        transform=ax_a.transAxes,
        color=COLORS["neutral"],
        va="top",
    )
    ax_a.text(
        0.02,
        0.64,
        f"OR {compare['r_OR']:.3f}  (95% CI {compare['r_CI_low']:.3f}\u2013{compare['r_CI_high']:.3f})",
        transform=ax_a.transAxes,
        fontsize=8.3,
        fontweight="bold",
        color=COLORS["mcop"],
        va="top",
    )
    ax_a.text(
        0.98,
        0.64,
        f"P = {format_p(compare['r_P_standard'])}",
        transform=ax_a.transAxes,
        ha="right",
        va="top",
        color=COLORS["neutral_dark"],
    )

    # B — implementation check
    impl = pd.DataFrame(
        {
            "label": ["R survey::svyglm", "Python Taylor-sandwich"],
            "OR": [compare["r_OR"], compare["python_OR"]],
            "CI_low": [compare["r_CI_low"], compare["python_CI_low"]],
            "CI_high": [compare["r_CI_high"], compare["python_CI_high"]],
        }
    )
    yb = np.array([1.0, 0.0])
    forest_rows(
        ax_b,
        impl,
        yb,
        [COLORS["neutral_dark"], COLORS["mcop"]],
        markers=["o", "s"],
        sizes=[4.4, 4.4],
    )
    set_or_axis(ax_b, (1.04, 1.49), [1.1, 1.2, 1.3, 1.4])
    ax_b.set_ylim(-0.65, 1.55)
    ax_b.set_yticks(yb, impl["label"])
    ax_b.set_xlabel("Odds ratio per doubling")
    ax_b.set_title("Independent implementation check", loc="left", pad=5, fontweight="bold")
    ax_b.text(
        0.98,
        0.04,
        f"Relative logOR difference = {compare['relative_logOR_change_pct']:.1e}%",
        transform=ax_b.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.0,
        color=COLORS["neutral"],
    )

    # C — LOCO forest
    pooled = pd.DataFrame(
        {
            "Dropped_cycle": ["Pooled"],
            "OR": [compare["r_OR"]],
            "CI_low": [compare["r_CI_low"]],
            "CI_high": [compare["r_CI_high"]],
        }
    )
    loco_plot = pd.concat([pooled, loco[["Dropped_cycle", "OR", "CI_low", "CI_high"]]], ignore_index=True)
    loco_plot["display"] = ["Pooled"] + [f"Drop {short_cycle(v)}" for v in loco["Dropped_cycle"]]
    yc = np.arange(len(loco_plot))[::-1]
    forest_rows(
        ax_c,
        loco_plot,
        yc,
        [COLORS["mcop"]] + [COLORS["mcop_light"]] * 7,
        markers=["D"] + ["o"] * 7,
        sizes=[5.2] + [3.9] * 7,
    )
    set_or_axis(ax_c, (0.96, 1.72), [1.0, 1.2, 1.4, 1.6])
    ax_c.set_ylim(-0.7, len(loco_plot) - 0.3)
    ax_c.set_yticks(yc, loco_plot["display"])
    ax_c.get_yticklabels()[0].set_fontweight("bold")
    ax_c.set_xlabel("Odds ratio per doubling of MCOP")
    ax_c.set_title("Leave-one-cycle-out stability", loc="left", pad=5, fontweight="bold")
    ax_c.text(
        0.98,
        0.02,
        "All seven 95% CIs exclude 1",
        transform=ax_c.transAxes,
        ha="right",
        va="bottom",
        color=COLORS["mcop"],
        fontweight="bold",
    )

    # D — cycle-specific transparency panel
    per_cycle = per_cycle.sort_values("Cycle").reset_index(drop=True)
    yd = np.arange(len(per_cycle))[::-1]
    d_colors = [
        COLORS["tumor"] if cycle == "2011-2012" else COLORS["neutral"]
        for cycle in per_cycle["Cycle"]
    ]
    open_markers = per_cycle["status"].eq("converged_with_warning").tolist()
    forest_rows(
        ax_d,
        per_cycle,
        yd,
        d_colors,
        sizes=[3.8] * len(per_cycle),
        open_markers=open_markers,
    )
    set_or_axis(ax_d, (0.52, 3.35), [0.6, 1.0, 2.0, 3.0])
    labels = [
        f"{short_cycle(cycle)}  ({int(n)} cases)"
        for cycle, n in zip(per_cycle["Cycle"], per_cycle["CRC_N"])
    ]
    ax_d.set_ylim(-1.15, len(per_cycle) - 0.25)
    ax_d.set_yticks(yd, labels)
    for tick, cycle in zip(ax_d.get_yticklabels(), per_cycle["Cycle"]):
        if cycle == "2011-2012":
            tick.set_color(COLORS["tumor"])
            tick.set_fontweight("bold")
    ax_d.set_xlabel("Cycle-specific odds ratio")
    ax_d.set_title("Per-cycle estimates", loc="left", pad=5, fontweight="bold")
    ax_d.text(
        0.98,
        0.97,
        f"Global MCOP-by-cycle interaction P = {interaction['P_F']:.4f}",
        transform=ax_d.transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
        color=COLORS["tumor"],
        fontweight="bold",
    )
    ax_d.text(
        0.98,
        0.02,
        "Open point: convergence warning",
        transform=ax_d.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.3,
        color=COLORS["neutral"],
    )

    for label, ax in zip("ABCD", [ax_a, ax_b, ax_c, ax_d]):
        add_panel_label(fig, ax, label)
    export_panel(fig, ax_a, panel_dir, "Figure2_panelA_primary_estimate_v1")
    export_panel(fig, ax_b, panel_dir, "Figure2_panelB_implementation_check_v1")
    export_panel(fig, ax_c, panel_dir, "Figure2_panelC_loco_v1")
    export_panel(fig, ax_d, panel_dir, "Figure2_panelD_per_cycle_v1")
    export_figure(fig, output_dir, "Figure2_nhanes_primary_v1")


def build_figure_3(source_dir: Path, output_dir: Path, panel_dir: Path) -> None:
    compare = numeric(
        pd.read_csv(source_dir / "figure2_primary_python_vs_r.csv"),
        ["r_OR", "r_CI_low", "r_CI_high", "r_P_standard", "r_N", "r_CRC_N"],
    ).iloc[0]
    sensitivity = numeric(
        pd.read_csv(source_dir / "figure3_sensitivity.csv"),
        ["OR", "CI_low", "CI_high", "P", "N", "CRC_N"],
    )
    coexposure = numeric(
        pd.read_csv(source_dir / "figure3_coexposure.csv"),
        ["OR", "CI_low", "CI_high", "P", "N", "CRC_N"],
    )
    rcs = numeric(
        pd.read_csv(source_dir / "figure3_rcs_curve_with_ci.csv"),
        ["mcop_log2", "mcop_ng_ml", "or_vs_median", "ci_low", "ci_high"],
    )
    rcs_meta = json.loads(
        (source_dir / "figure3_rcs_curve_with_ci_metadata.json").read_text(encoding="utf-8")
    )
    quartiles = numeric(
        pd.read_csv(source_dir / "figure3_weighted_quartiles.csv"),
        ["OR", "CI_low", "CI_high", "P", "P_trend"],
    )
    quartiles = quartiles.loc[
        (quartiles["method"] == "survey_weighted_cutpoints")
        & quartiles["Quartile"].notna()
        & quartiles["Quartile"].ne("")
    ].copy()

    requested = [
        "Age_ge_40_7_cycle",
        "Exclude_diagnosis_lt_1y",
        "Exclude_diagnosis_lt_2y",
        "Exclude_diagnosis_lt_5y",
        "Exclude_top_1pct",
        "Exclude_top_2.5pct",
        "Creatinine_normalized_MCOP",
    ]
    sens = sensitivity.loc[sensitivity["Analysis"].isin(requested)].copy()
    sens = sens.set_index("Analysis").loc[requested].reset_index()
    sens["display"] = [
        "Age >=40 y",
        "Exclude diagnosis <1 y",
        "Exclude diagnosis <2 y",
        "Exclude diagnosis <5 y",
        "Exclude top 1%",
        "Exclude top 2.5%",
        "Creatinine-normalized",
    ]
    primary = pd.DataFrame(
        {
            "Analysis": ["Primary"],
            "OR": [compare["r_OR"]],
            "CI_low": [compare["r_CI_low"]],
            "CI_high": [compare["r_CI_high"]],
            "N": [compare["r_N"]],
            "CRC_N": [compare["r_CRC_N"]],
            "display": ["Primary model"],
        }
    )
    sens_plot = pd.concat([primary, sens], ignore_index=True)
    assert len(sens_plot) == 8

    coex = coexposure.loc[coexposure["Model_role"].eq("coexposure_adjusted")].copy()
    coex = coex.set_index("Secondary_exposure").loc[
        ["MEHHP", "MEOHP", "MECPP", "MBzP", "PhthalateBurden_excl_MCOP"]
    ].reset_index()
    coex["display"] = ["+ MEHHP", "+ MEOHP", "+ MECPP", "+ MBzP", "+ phthalate burden"]
    coex_primary = pd.DataFrame(
        {
            "display": ["Primary model"],
            "OR": [compare["r_OR"]],
            "CI_low": [compare["r_CI_low"]],
            "CI_high": [compare["r_CI_high"]],
        }
    )
    coex_plot = pd.concat([coex_primary, coex[["display", "OR", "CI_low", "CI_high"]]], ignore_index=True)

    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, 5.75))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.45, 1.08],
        width_ratios=[1.33, 1.00],
        hspace=0.46,
        wspace=0.46,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    fig.subplots_adjust(left=0.205, right=0.975, bottom=0.105, top=0.94)

    # A — grouped sensitivity forest
    ya = np.array([10.0, 8.2, 6.6, 5.6, 4.6, 2.8, 1.8, 0.2])
    ax_a.axhspan(7.65, 8.75, color=COLORS["neutral_pale"], zorder=-2)
    ax_a.axhspan(4.05, 7.15, color=COLORS["mcop_pale"], alpha=0.45, zorder=-2)
    ax_a.axhspan(1.25, 3.35, color=COLORS["neutral_pale"], zorder=-2)
    forest_rows(
        ax_a,
        sens_plot,
        ya,
        [COLORS["mcop"]] + [COLORS["mcop_light"]] * 7,
        markers=["D"] + ["o"] * 7,
        sizes=[5.0] + [3.9] * 7,
    )
    set_or_axis(ax_a, (0.95, 1.62), [1.0, 1.2, 1.4, 1.6])
    labels = [
        f"{label}  ({int(n):,}/{int(cases)})"
        for label, n, cases in zip(sens_plot["display"], sens_plot["N"], sens_plot["CRC_N"])
    ]
    ax_a.set_ylim(-0.75, 10.85)
    ax_a.set_yticks(ya, labels)
    ax_a.tick_params(axis="y", labelsize=5.8, pad=2)
    ax_a.get_yticklabels()[0].set_fontweight("bold")
    ax_a.set_xlabel("Odds ratio per doubling of MCOP")
    ax_a.set_title("Sensitivity analyses (N / CRC cases)", loc="left", pad=5, fontweight="bold")

    # B — co-exposure forest
    yb = np.arange(len(coex_plot))[::-1]
    forest_rows(
        ax_b,
        coex_plot,
        yb,
        [COLORS["neutral_dark"]] + [COLORS["mcop"]] * 5,
        markers=["D"] + ["o"] * 5,
        sizes=[4.8] + [4.0] * 5,
    )
    set_or_axis(ax_b, (0.96, 1.55), [1.0, 1.2, 1.4])
    ax_b.set_ylim(-0.7, len(coex_plot) - 0.25)
    ax_b.set_yticks(yb, coex_plot["display"])
    ax_b.get_yticklabels()[0].set_fontweight("bold")
    ax_b.set_xlabel("MCOP odds ratio per doubling")
    ax_b.set_title("Pairwise co-exposure adjustment", loc="left", pad=5, fontweight="bold")
    ax_b.text(
        0.98,
        0.03,
        "All adjusted 95% CIs exclude 1",
        transform=ax_b.transAxes,
        ha="right",
        va="bottom",
        color=COLORS["mcop"],
        fontweight="bold",
    )

    # C — survey-weighted RCS with reconstructed CI
    ax_c.fill_between(
        rcs["mcop_ng_ml"].to_numpy(float),
        rcs["ci_low"].to_numpy(float),
        rcs["ci_high"].to_numpy(float),
        color=COLORS["mcop_light"],
        alpha=0.22,
        linewidth=0,
        zorder=1,
    )
    ax_c.plot(
        rcs["mcop_ng_ml"],
        rcs["or_vs_median"],
        color=COLORS["mcop"],
        lw=1.7,
        zorder=2,
    )
    ax_c.axhline(1.0, color=COLORS["neutral"], lw=0.8, ls=(0, (3, 2)))
    if not (
        (rcs["mcop_ng_ml"] > 0).all()
        and (rcs["ci_low"] > 0).all()
        and (rcs["ci_high"] > 0).all()
        and (rcs["or_vs_median"] > 0).all()
    ):
        raise ValueError("RCS log-scale coordinates and intervals must be strictly positive")
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlim(float(rcs["mcop_ng_ml"].min()), float(rcs["mcop_ng_ml"].max()))
    ax_c.set_ylim(0.10, 5.3)
    ax_c.xaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 100]))
    ax_c.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax_c.xaxis.set_minor_locator(NullLocator())
    ax_c.yaxis.set_major_locator(FixedLocator([0.1, 0.25, 0.5, 1, 2, 4]))
    ax_c.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax_c.yaxis.set_minor_locator(NullLocator())
    ax_c.set_xlabel("Urinary MCOP (ng/mL; log scale)")
    ax_c.set_ylabel("Adjusted OR vs survey-weighted median")
    ax_c.set_title("Survey-weighted restricted cubic spline", loc="left", pad=5, fontweight="bold")
    ax_c.text(
        0.03,
        0.96,
        f"P overall = {rcs_meta['overall_P_F']:.4f}\nP nonlinear = {rcs_meta['nonlinear_P_F']:.3f}",
        transform=ax_c.transAxes,
        ha="left",
        va="top",
        color=COLORS["neutral_dark"],
    )
    ax_c.text(
        0.98,
        0.04,
        "Displayed over weighted 5th\u201395th percentiles",
        transform=ax_c.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.8,
        color=COLORS["neutral"],
    )

    # D — secondary quartile analysis
    quartiles = quartiles.sort_values("Quartile").reset_index(drop=True)
    qx = np.arange(1, 5)
    q_or = quartiles["OR"].to_numpy(float)
    q_low = quartiles["CI_low"].to_numpy(float)
    q_high = quartiles["CI_high"].to_numpy(float)
    ax_d.axhline(1.0, color=COLORS["neutral"], lw=0.8, ls=(0, (3, 2)))
    for x, estimate, low, high in zip(qx, q_or, q_low, q_high):
        if x == 1:
            ax_d.scatter(x, 1.0, marker="s", s=20, color=COLORS["normal"], zorder=3)
        else:
            ax_d.errorbar(
                x,
                estimate,
                yerr=[[estimate - low], [high - estimate]],
                fmt="o",
                ms=4.2,
                color=COLORS["mcop_light"],
                ecolor=COLORS["mcop_light"],
                elinewidth=1.0,
                capsize=2.0,
                zorder=3,
            )
    if not (np.all(q_or > 0) and np.all(q_low > 0) and np.all(q_high > 0)):
        raise ValueError("Quartile odds ratios and intervals must be strictly positive")
    ax_d.set_yscale("log")
    ax_d.set_ylim(0.28, 4.8)
    ax_d.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4]))
    ax_d.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax_d.yaxis.set_minor_locator(NullLocator())
    ax_d.set_xticks(qx, ["Q1", "Q2", "Q3", "Q4"])
    ax_d.set_ylabel("Odds ratio vs Q1")
    ax_d.set_xlabel("Survey-weighted MCOP quartile")
    ax_d.set_title("Secondary categorical analysis", loc="left", pad=5, fontweight="bold")
    p_trend_values = quartiles.loc[quartiles["P_trend"].notna(), "P_trend"]
    if p_trend_values.empty:
        raise ValueError("Quartile trend P value is missing")
    p_trend = float(p_trend_values.iloc[0])
    ax_d.text(
        0.98,
        0.96,
        f"P trend = {p_trend:.3f}\n70 CRC cases across four groups",
        transform=ax_d.transAxes,
        ha="right",
        va="top",
        color=COLORS["neutral"],
    )

    for label, ax in zip("ABCD", [ax_a, ax_b, ax_c, ax_d]):
        add_panel_label(fig, ax, label)
    export_panel(fig, ax_a, panel_dir, "Figure3_panelA_sensitivity_v1")
    export_panel(fig, ax_b, panel_dir, "Figure3_panelB_coexposure_v1")
    export_panel(fig, ax_c, panel_dir, "Figure3_panelC_rcs_v1")
    export_panel(fig, ax_d, panel_dir, "Figure3_panelD_quartiles_v1")
    export_figure(fig, output_dir, "Figure3_robustness_v1")


def paired_from_long(
    frame: pd.DataFrame,
    index_columns: list[str],
    score: str,
) -> pd.DataFrame:
    wide = frame.pivot_table(
        index=index_columns,
        columns="group",
        values=score,
        aggfunc="first",
    )
    if not {"normal", "tumor"}.issubset(wide.columns):
        raise ValueError(f"Missing paired groups for {score}")
    paired = wide[["normal", "tumor"]]
    complete_mask = paired.notna().all(axis=1)
    return paired.loc[complete_mask].reset_index()


def build_figure_4(source_dir: Path, output_dir: Path, panel_dir: Path) -> None:
    bulk_scores = numeric(
        pd.read_csv(source_dir / "figure4_bulk_scores.csv"),
        ["PPAR_nuclear_receptor_score", "inflammatory_RELA_STAT3_score"],
    )
    bulk_manifest = pd.read_csv(source_dir / "figure4_bulk_sample_manifest.csv")
    tcga_summary = numeric(
        pd.read_csv(source_dir / "figure4_tcga_paired_summary.csv"),
        ["paired_n", "median_delta_tumor_minus_normal", "p_value"],
    )
    census_scores = numeric(
        pd.read_csv(source_dir / "figure4_census_donor_scores.csv"),
        ["PPAR_nuclear_receptor_score", "RELA_STAT3_score"],
    )
    census_summary = numeric(
        pd.read_csv(source_dir / "figure4_census_paired_summary.csv"),
        ["paired_donors", "median_delta_tumor_minus_normal", "p_value"],
    )
    gse_scores = numeric(
        pd.read_csv(source_dir / "figure4_gse144735_scores.csv"),
        ["PPAR_nuclear_receptor_score", "RELA_STAT3_score"],
    )
    gse_summary = numeric(
        pd.read_csv(source_dir / "figure4_gse144735_paired_summary.csv"),
        ["paired_donors", "median_delta_tumor_minus_normal", "p_value"],
    )

    # TCGA matched pairs
    bulk = bulk_scores.loc[
        bulk_scores["contrast"].eq("TCGA_primary_vs_TCGA_solid_normal")
    ].merge(
        bulk_manifest[["sample_id", "patient_id"]],
        on="sample_id",
        how="left",
        validate="many_to_one",
    )
    bulk["group"] = bulk["group"].map(
        {
            "TCGA_CRC_primary_tumor": "tumor",
            "TCGA_CRC_solid_normal": "normal",
        }
    )
    bulk_eligible_mask = bulk["patient_id"].notna() & bulk["group"].notna()
    bulk_eligible = bulk.loc[bulk_eligible_mask].copy()
    tcga_pairs = paired_from_long(
        bulk_eligible,
        ["patient_id"],
        "PPAR_nuclear_receptor_score",
    )
    assert len(tcga_pairs) == 32

    # Census pairs and compartment deltas
    census_ppar = paired_from_long(
        census_scores,
        ["dataset_id", "donor_key", "compartment"],
        "PPAR_nuclear_receptor_score",
    )
    census_rela = paired_from_long(
        census_scores,
        ["dataset_id", "donor_key", "compartment"],
        "RELA_STAT3_score",
    )
    epi_ppar = census_ppar.loc[census_ppar["compartment"].eq("epithelial")].copy()
    epi_rela = census_rela.loc[census_rela["compartment"].eq("epithelial")].copy()
    assert len(epi_ppar) == 36 and len(epi_rela) == 36
    expected_pairs = {"epithelial": 36, "endothelial": 33, "fibroblast": 31, "myeloid": 35}
    observed_pairs = census_ppar.groupby("compartment").size().to_dict()
    assert observed_pairs == expected_pairs

    # GSE144735 matched pairs
    gse_pairs = paired_from_long(
        gse_scores,
        ["dataset_id", "donor_key"],
        "PPAR_nuclear_receptor_score",
    )
    assert len(gse_pairs) == 6

    tcga_stat = tcga_summary.loc[
        tcga_summary["score"].eq("PPAR_nuclear_receptor_score")
    ].iloc[0]
    census_ppar_stat = census_summary.loc[
        census_summary["compartment"].eq("epithelial")
        & census_summary["score"].eq("PPAR_nuclear_receptor_score")
    ].iloc[0]
    census_rela_stat = census_summary.loc[
        census_summary["compartment"].eq("epithelial")
        & census_summary["score"].eq("RELA_STAT3_score")
    ].iloc[0]
    gse_stat = gse_summary.loc[
        gse_summary["group_mode"].eq("core_tumor_vs_normal")
        & gse_summary["score"].eq("PPAR_nuclear_receptor_score")
    ].iloc[0]

    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, 5.90))
    outer = fig.add_gridspec(
        2,
        5,
        height_ratios=[1.02, 1.08],
        hspace=0.48,
        wspace=0.42,
    )
    ax_a = fig.add_subplot(outer[0, 0:2])
    b_grid = outer[0, 2:5].subgridspec(1, 2, wspace=0.32)
    ax_b1 = fig.add_subplot(b_grid[0, 0])
    ax_b2 = fig.add_subplot(b_grid[0, 1])
    c_grid = outer[1, 0:3].subgridspec(1, 2, width_ratios=[1.0, 0.46], wspace=0.05)
    ax_c = fig.add_subplot(c_grid[0, 0])
    ax_c_stats = fig.add_subplot(c_grid[0, 1])
    ax_d = fig.add_subplot(outer[1, 3:5])
    fig.subplots_adjust(left=0.145, right=0.975, bottom=0.105, top=0.875)
    paired_plot(
        ax_a,
        tcga_pairs["normal"].to_numpy(float),
        tcga_pairs["tumor"].to_numpy(float),
        "Standardized PPAR/NR score",
        (
            "TCGA matched tumor\u2013normal"
            f"\nn = {int(tcga_stat['paired_n'])}; median delta = {tcga_stat['median_delta_tumor_minus_normal']:.3f}; "
            f"P = {format_p(tcga_stat['p_value'])}"
        ),
        point_size=12,
        title_size=7.0,
    )

    paired_plot(
        ax_b1,
        epi_ppar["normal"].to_numpy(float),
        epi_ppar["tumor"].to_numpy(float),
        "Standardized module score",
        (
            "PPAR/NR"
            f"\nmedian delta {census_ppar_stat['median_delta_tumor_minus_normal']:.3f}; "
            f"P {format_p(census_ppar_stat['p_value'])}"
        ),
        point_size=10,
        title_size=6.2,
    )
    paired_plot(
        ax_b2,
        epi_rela["normal"].to_numpy(float),
        epi_rela["tumor"].to_numpy(float),
        None,
        (
            "RELA/STAT3"
            f"\nmedian delta +{census_rela_stat['median_delta_tumor_minus_normal']:.3f}; "
            f"P {format_p(census_rela_stat['p_value'])}"
        ),
        point_size=10,
        title_size=6.2,
    )

    # C — compartment-specific paired deltas
    compartments = ["epithelial", "endothelial", "fibroblast", "myeloid"]
    positions = np.arange(len(compartments))[::-1]
    comp_colors = {
        "epithelial": COLORS["tumor"],
        "endothelial": COLORS["normal"],
        "fibroblast": COLORS["normal"],
        "myeloid": COLORS["mcop"],
    }
    all_deltas: list[np.ndarray] = []
    for compartment, pos in zip(compartments, positions):
        subset = census_ppar.loc[census_ppar["compartment"].eq(compartment)].copy()
        delta = (subset["tumor"] - subset["normal"]).to_numpy(float)
        all_deltas.append(delta)
        color = comp_colors[compartment]
        bp = ax_c.boxplot(
            [delta],
            positions=[pos],
            orientation="horizontal",
            widths=0.46,
            patch_artist=True,
            showfliers=False,
            manage_ticks=False,
            boxprops={"facecolor": mpl.colors.to_rgba(color, 0.16), "edgecolor": color, "linewidth": 0.9},
            medianprops={"color": color, "linewidth": 1.6},
            whiskerprops={"color": color, "linewidth": 0.8},
            capprops={"color": color, "linewidth": 0.8},
        )
        _ = bp
        # Deterministic sunflower-phase jitter prevents overplotting without altering values.
        jitter = 0.065 * np.sin(np.arange(len(delta), dtype=float) * 2.3999632297)
        ax_c.scatter(
            delta,
            np.full(len(delta), pos) + jitter,
            s=9.5,
            facecolor=mpl.colors.to_rgba(color, 0.78),
            edgecolor=COLORS["white"],
            linewidth=0.25,
            zorder=3,
        )
    ax_c.axvline(0, color=COLORS["neutral"], lw=0.8, ls=(0, (3, 2)), zorder=0)
    flat = np.concatenate(all_deltas)
    span = float(flat.max() - flat.min())
    ax_c.set_xlim(float(flat.min() - 0.05 * span), float(flat.max() + 0.08 * span))
    labels = [f"{name.capitalize()}  (n = {expected_pairs[name]})" for name in compartments]
    ax_c.set_yticks(positions, labels)
    ax_c.tick_params(axis="y", labelsize=5.8, pad=2)
    ax_c.set_xlabel("Tumor minus normal PPAR/NR score")
    ax_c.set_title("Compartment-specific paired donor deltas", loc="left", pad=5, fontweight="bold")
    summary_lookup = census_summary.loc[
        census_summary["score"].eq("PPAR_nuclear_receptor_score")
    ].set_index("compartment")
    ax_c_stats.set_xlim(0, 1)
    ax_c_stats.set_ylim(-0.5, 3.5)
    ax_c_stats.set_title("Median delta / P", loc="left", pad=5, fontsize=6.2, color=COLORS["neutral"])
    ax_c_stats.axis("off")
    for compartment, pos in zip(compartments, positions):
        row = summary_lookup.loc[compartment]
        ax_c_stats.text(
            0.02,
            pos,
            f"{row['median_delta_tumor_minus_normal']:+.3f}\nP {format_p(row['p_value'])}",
            ha="left",
            va="center",
            fontsize=5.8,
            color=comp_colors[compartment],
        )

    paired_plot(
        ax_d,
        gse_pairs["normal"].to_numpy(float),
        gse_pairs["tumor"].to_numpy(float),
        "Standardized PPAR/NR score",
        (
            "GSE144735 matched epithelium"
            f"\nn = {int(gse_stat['paired_donors'])}; median delta = {gse_stat['median_delta_tumor_minus_normal']:.3f}; "
            f"P = {format_p(gse_stat['p_value'])}"
        ),
        point_size=15,
        title_size=7.0,
    )
    ax_d.text(
        0.98,
        0.04,
        "Directionally concordant; underpowered",
        transform=ax_d.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.0,
        color=COLORS["neutral"],
    )

    add_panel_label(fig, ax_a, "A")
    add_panel_label(fig, ax_b1, "B")
    add_panel_label(fig, ax_c, "C")
    add_panel_label(fig, ax_d, "D")
    b1_pos = ax_b1.get_position()
    fig.text(
        b1_pos.x0,
        b1_pos.y1 + 0.064,
        "Census matched epithelium (n = 36)",
        ha="left",
        va="bottom",
        fontsize=7.4,
        fontweight="bold",
        color=COLORS["neutral_dark"],
    )
    export_panel(fig, ax_a, panel_dir, "Figure4_panelA_tcga_pairs_v1")
    export_panel(fig, [ax_b1, ax_b2], panel_dir, "Figure4_panelB_census_epithelium_v1")
    export_panel(fig, [ax_c, ax_c_stats], panel_dir, "Figure4_panelC_compartments_v1")
    export_panel(fig, ax_d, panel_dir, "Figure4_panelD_gse144735_v1")
    export_figure(fig, output_dir, "Figure4_ppar_singlecell_v1")


def prepare_figure_5_data(source_dir: Path) -> None:
    pd.DataFrame(
        [
            ("A1", "267 environmental chemicals", "discovery", "observed"),
            ("A2", "CTD x GeneCards CRC enrichment", "discovery", "observed"),
            ("A3", "MiNP/DINP exposure axis nominated", "discovery", "observed"),
            ("B1", "Urinary MCOP in NHANES 2005-2018", "human biomonitoring", "observed"),
            ("B2", "Prevalent CRC: OR 1.246 per doubling", "human biomonitoring", "observed association"),
            ("C1", "CRC tumor epithelium", "biological state", "observed"),
            ("C2", "PPAR/NR down; RELA/STAT3 up", "biological state", "observed disease-state remodeling"),
        ],
        columns=["node_id", "label", "layer", "evidence_status"],
    ).to_csv(source_dir / "figure5_evidence_nodes.csv", index=False)
    pd.DataFrame(
        [
            ("A1", "A2", "screening workflow", "solid", "observed"),
            ("A2", "A3", "audited nomination", "solid", "observed"),
            ("B1", "B2", "MCOP-CRC association", "solid bidirectional", "observed association"),
            ("C1", "C2", "CRC-state remodeling", "solid bidirectional", "observed disease-state relationship"),
            ("A3", "C2", "candidate mechanistic bridge", "dashed", "hypothetical"),
        ],
        columns=["source", "target", "relation", "connector", "evidence_status"],
    ).to_csv(source_dir / "figure5_evidence_links.csv", index=False)


def build_figure_5(source_dir: Path, output_dir: Path, panel_dir: Path) -> None:
    prepare_figure_5_data(source_dir)
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, 3.62))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.08, 1.0], wspace=0.24)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.20, top=0.88)
    for axis in [ax_a, ax_b, ax_c]:
        axis.set_axis_off()

    # A — discovery layer
    rounded_box(
        ax_a,
        0.02,
        0.02,
        0.96,
        0.94,
        facecolor=COLORS["white"],
        edgecolor=COLORS["normal_light"],
        linewidth=0.9,
    )
    rounded_box(
        ax_a,
        0.08,
        0.82,
        0.84,
        0.11,
        facecolor=COLORS["neutral_pale"],
        edgecolor=COLORS["normal"],
    )
    ax_a.text(
        0.50,
        0.875,
        "DATA-DRIVEN DISCOVERY",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=6.3,
        fontweight="bold",
        color=COLORS["neutral_dark"],
    )
    rounded_box(
        ax_a,
        0.16,
        0.62,
        0.68,
        0.12,
        facecolor=COLORS["neutral_pale"],
        edgecolor=COLORS["normal"],
    )
    ax_a.text(0.50, 0.68, "267 environmental chemicals", transform=ax_a.transAxes, ha="center", va="center", fontsize=6.2)
    axes_arrow(ax_a, (0.50, 0.61), (0.50, 0.53))
    rounded_box(
        ax_a,
        0.16,
        0.40,
        0.68,
        0.12,
        facecolor=COLORS["neutral_pale"],
        edgecolor=COLORS["normal"],
    )
    ax_a.text(0.50, 0.46, "CTD x GeneCards CRC enrichment", transform=ax_a.transAxes, ha="center", va="center", fontsize=6.0)
    axes_arrow(ax_a, (0.50, 0.39), (0.50, 0.31), color=COLORS["mcop"])
    rounded_box(
        ax_a,
        0.10,
        0.13,
        0.80,
        0.17,
        facecolor=COLORS["mcop_pale"],
        edgecolor=COLORS["mcop"],
        linewidth=1.2,
    )
    ax_a.text(
        0.50,
        0.215,
        "MiNP/DINP exposure axis\nnominated",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=7.0,
        fontweight="bold",
        color=COLORS["mcop"],
    )
    ax_a.text(
        0.50,
        0.065,
        "FDR + degree-matched permutation",
        transform=ax_a.transAxes,
        ha="center",
        va="center",
        fontsize=5.3,
        color=COLORS["neutral"],
    )

    # B — human association layer (visual center)
    rounded_box(
        ax_b,
        0.02,
        0.02,
        0.96,
        0.94,
        facecolor=COLORS["white"],
        edgecolor=COLORS["mcop_light"],
        linewidth=1.1,
    )
    rounded_box(
        ax_b,
        0.08,
        0.82,
        0.84,
        0.11,
        facecolor=COLORS["mcop_pale"],
        edgecolor=COLORS["mcop"],
    )
    ax_b.text(
        0.50,
        0.875,
        "HUMAN BIOMONITORING",
        transform=ax_b.transAxes,
        ha="center",
        va="center",
        fontsize=6.3,
        fontweight="bold",
        color=COLORS["mcop"],
    )
    ax_b.text(
        0.50,
        0.715,
        "Urinary MCOP  |  NHANES 2005-2018",
        transform=ax_b.transAxes,
        ha="center",
        va="center",
        fontsize=6.4,
        fontweight="bold",
        color=COLORS["neutral_dark"],
    )
    rounded_box(
        ax_b,
        0.10,
        0.40,
        0.80,
        0.24,
        facecolor=COLORS["mcop_pale"],
        edgecolor=COLORS["mcop"],
        linewidth=1.3,
    )
    ax_b.text(
        0.50,
        0.555,
        "OR 1.246 per doubling",
        transform=ax_b.transAxes,
        ha="center",
        va="center",
        fontsize=8.2,
        fontweight="bold",
        color=COLORS["mcop"],
    )
    ax_b.text(
        0.50,
        0.465,
        "95% CI 1.077-1.440  |  P = 0.0034",
        transform=ax_b.transAxes,
        ha="center",
        va="center",
        fontsize=6.0,
        color=COLORS["neutral_dark"],
    )
    ax_b.text(
        0.50,
        0.33,
        "9,936 adults  |  70 CRC cases",
        transform=ax_b.transAxes,
        ha="center",
        va="center",
        fontsize=6.1,
        color=COLORS["neutral_dark"],
    )
    rounded_box(ax_b, 0.10, 0.14, 0.36, 0.10, facecolor=COLORS["neutral_pale"], edgecolor=COLORS["normal"])
    rounded_box(ax_b, 0.54, 0.14, 0.36, 0.10, facecolor=COLORS["neutral_pale"], edgecolor=COLORS["normal"])
    ax_b.text(0.28, 0.19, "R/Python agree", transform=ax_b.transAxes, ha="center", va="center", fontsize=5.5, fontweight="bold")
    ax_b.text(0.72, 0.19, "7/7 LOCO CIs > 1", transform=ax_b.transAxes, ha="center", va="center", fontsize=5.5, fontweight="bold")
    ax_b.text(
        0.50,
        0.065,
        "Observed association - not causal proof",
        transform=ax_b.transAxes,
        ha="center",
        va="center",
        fontsize=5.3,
        color=COLORS["neutral"],
    )

    # C — CRC biological-state layer
    rounded_box(
        ax_c,
        0.02,
        0.02,
        0.96,
        0.94,
        facecolor=COLORS["white"],
        edgecolor=COLORS["tumor_light"],
        linewidth=1.0,
    )
    rounded_box(
        ax_c,
        0.08,
        0.82,
        0.84,
        0.11,
        facecolor="#F5E8EC",
        edgecolor=COLORS["tumor"],
    )
    ax_c.text(
        0.50,
        0.875,
        "CRC BIOLOGICAL STATE",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=6.3,
        fontweight="bold",
        color=COLORS["tumor"],
    )
    ax_c.text(
        0.50,
        0.72,
        "Matched tumor-derived epithelium",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=6.2,
        fontweight="bold",
        color=COLORS["neutral_dark"],
    )
    rounded_box(ax_c, 0.10, 0.47, 0.80, 0.16, facecolor="#F5E8EC", edgecolor=COLORS["tumor"])
    ax_c.text(
        0.50,
        0.55,
        "PPAR/NR program  ↓",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=6.6,
        fontweight="bold",
        color=COLORS["tumor"],
    )
    rounded_box(ax_c, 0.10, 0.27, 0.80, 0.13, facecolor="#F8EEE8", edgecolor=COLORS["inflammation"])
    ax_c.text(
        0.50,
        0.335,
        "RELA/STAT3 program  ↑",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=6.4,
        fontweight="bold",
        color=COLORS["inflammation"],
    )
    ax_c.text(
        0.50,
        0.19,
        "Myeloid PPAR/NR: opposite direction",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=5.6,
        color=COLORS["mcop"],
        fontweight="bold",
    )
    ax_c.text(
        0.50,
        0.085,
        "Paired TCGA + CELLxGENE\nGSE144735 directional only",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=5.3,
        color=COLORS["neutral"],
    )

    # The only cross-layer causal-looking connector is explicitly hypothetical.
    a_pos = ax_a.get_position()
    c_pos = ax_c.get_position()
    bridge_y = 0.105
    fig.add_artist(
        FancyArrowPatch(
            (a_pos.x0 + 0.55 * a_pos.width, bridge_y),
            (c_pos.x0 + 0.50 * c_pos.width, bridge_y),
            transform=fig.transFigure,
            arrowstyle="-|>",
            mutation_scale=8.0,
            linewidth=1.0,
            linestyle=(0, (4, 3)),
            color=COLORS["mcop_light"],
            connectionstyle="arc3,rad=0",
            clip_on=False,
        )
    )
    fig.text(
        0.51,
        bridge_y + 0.018,
        "candidate mechanistic bridge",
        ha="center",
        va="bottom",
        fontsize=5.8,
        fontweight="bold",
        color=COLORS["mcop"],
    )
    fig.text(
        0.51,
        0.035,
        "Unproven: temporality, direct DINP-to-epithelium perturbation, and mediation",
        ha="center",
        va="bottom",
        fontsize=5.4,
        color=COLORS["neutral"],
    )

    for label, axis in zip("ABC", [ax_a, ax_b, ax_c]):
        add_panel_label(fig, axis, label, dx=0.027)
    export_panel(fig, ax_a, panel_dir, "Figure5_panelA_discovery_layer_v1")
    export_panel(fig, ax_b, panel_dir, "Figure5_panelB_human_biomonitoring_v1")
    export_panel(fig, ax_c, panel_dir, "Figure5_panelC_crc_state_v1")
    export_figure(fig, output_dir, "Figure5_integrated_model_v1")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    panel_dir = args.output_dir / "panels"
    configure_style()
    build_figure_1(args.repo_root, args.source_dir, args.output_dir, panel_dir)
    build_figure_2(args.source_dir, args.output_dir, panel_dir)
    build_figure_3(args.source_dir, args.output_dir, panel_dir)
    build_figure_4(args.source_dir, args.output_dir, panel_dir)
    build_figure_5(args.source_dir, args.output_dir, panel_dir)
    print("Exported Figure 1-5 plus individual panels as SVG, PDF and 600 dpi PNG.")


if __name__ == "__main__":
    main()
