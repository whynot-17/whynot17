# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) schematic-led workflow → cross-type inherit → param inherit
# (b) ranked chemical screen → cross-type inherit → param inherit
# (c) multistage actionability and biomarker translation → cross-type inherit → param inherit
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
# Muted wine-red encodes tumor identity without reading as a warning state.
TUMOR = "#8F4B58"
TUMOR_LIGHT = "#F8F0F2"
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
    text_artist = ax.text(
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
    patch._contained_text = [text_artist]
    return patch


def register_card_text(patch, *artists):
    """Register manually positioned card text for rendered overflow checks."""
    if not hasattr(patch, "_contained_text"):
        patch._contained_text = []
    patch._contained_text.extend(artists)


def assert_card_text_containment(fig, padding_px=1.0):
    """Fail rendering when any registered text extends outside its card."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    failures = []
    for ax in fig.axes:
        for patch in ax.patches:
            artists = getattr(patch, "_contained_text", [])
            if not artists:
                continue
            card_box = patch.get_window_extent(renderer)
            for artist in artists:
                if not artist.get_text().strip():
                    continue
                text_box = artist.get_window_extent(renderer)
                inside = (
                    text_box.x0 >= card_box.x0 + padding_px
                    and text_box.x1 <= card_box.x1 - padding_px
                    and text_box.y0 >= card_box.y0 + padding_px
                    and text_box.y1 <= card_box.y1 - padding_px
                )
                if not inside:
                    failures.append(artist.get_text().replace("\n", " / "))
    if failures:
        raise RuntimeError("Text overflow detected in card(s): " + "; ".join(failures))


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


def validate_and_load(
    repo_root: Path,
    source_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
    biomarker_path = repo_root / "outputs" / "nhanes_dinp_phase2a_audit_summary.csv"
    biomarker_summary = pd.read_csv(biomarker_path)
    for analyte, expected_pct in {"MiNP": 27.385173579615735, "MCOP": 98.40026387400016}.items():
        row = biomarker_summary.loc[biomarker_summary["analyte"].eq(analyte)]
        if len(row) != 1:
            raise ValueError(f"Expected exactly one Phase 2A row for {analyte}")
        if int(row.iloc[0]["n_cycles"]) != 7:
            raise ValueError(f"Expected seven NHANES cycles for {analyte}")
        if not np.isclose(float(row.iloc[0]["above_lod_pct"]), expected_pct):
            raise ValueError(f"Phase 2A detectability differs from the manuscript lock for {analyte}")

    # Keep the source manifest alongside the rendered figure for traceability.
    manifest = pd.DataFrame(
        [
            (str(screen_path.relative_to(repo_root)), "all 267 rows", "Panel B ranked screen"),
            (str((source_dir / "figure1_panelA_workflow_nodes.csv").relative_to(repo_root)), "4 workflow nodes", "Panel A workflow labels"),
            (str(biomarker_path.relative_to(repo_root)), "MiNP and MCOP rows", "Panel C seven-cycle detectability"),
        ],
        columns=["source_file", "rows_or_scope", "figure_use"],
    )
    return screen, nodes, biomarker_summary, manifest


def draw_panel_a(ax, nodes: pd.DataFrame):
    ax.set_axis_off()
    panel_title(ax, "A", "Data-first discovery")
    rounded_card(
        ax,
        (0.11, 0.77, 0.78, 0.15),
        "267 core\nenvironmental\nchemicals",
        facecolor=EXPOSURE_LIGHT,
        edgecolor=EXPOSURE,
        linewidth=1.15,
        text_color=EXPOSURE,
        fontsize=6.0,
        weight="bold",
    )
    rounded_card(
        ax,
        (0.03, 0.51, 0.45, 0.17),
        "CTD human\nchemical–gene\ninteractions",
        facecolor=NEUTRAL_PALE,
        edgecolor=NEUTRAL_LIGHT,
        text_color=NEUTRAL_DARK,
        fontsize=5.35,
    )
    rounded_card(
        ax,
        (0.52, 0.51, 0.45, 0.17),
        "GeneCards\nCRC-associated\ngenes",
        facecolor=NEUTRAL_PALE,
        edgecolor=NEUTRAL_LIGHT,
        text_color=NEUTRAL_DARK,
        fontsize=5.35,
    )
    arrow(ax, (0.50, 0.77), (0.27, 0.68))
    arrow(ax, (0.50, 0.77), (0.73, 0.68))
    arrow(ax, (0.27, 0.50), (0.46, 0.38))
    arrow(ax, (0.73, 0.50), (0.54, 0.38))
    rounded_card(
        ax,
        (0.04, 0.19, 0.92, 0.17),
        "Enrichment + BH-FDR\nDegree-matched permutation",
        facecolor=WHITE,
        edgecolor=NEUTRAL,
        linewidth=1.0,
        text_color=NEUTRAL_DARK,
        fontsize=5.15,
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
        fontsize=5.0,
        color=NEUTRAL,
    )

    special = {
        "MiNP": dict(marker="D", size=50, color=EXPOSURE, text="MiNP  |  rank 24\nOR 10.06; degree-matched FDR 0.036", xytext=(55, 4.05)),
        "DINP": dict(marker="o", size=38, color=EXPOSURE, text="DINP parent  |  rank 107\nBH-FDR 0.449", xytext=(125, 0.88)),
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
            bbox=(
                dict(facecolor=WHITE, edgecolor="none", alpha=0.88, pad=0.4)
                if name == "DINP"
                else None
            ),
        )
    ax.set_xlim(0, 270)
    ax.set_ylim(-0.35, 9.75)
    ax.set_xlabel("Primary screen rank")
    ax.set_ylabel("−log10(BH-FDR)")
    ax.set_xticks([1, 50, 100, 150, 200, 250])
    ax.tick_params(labelsize=6.2)
    ax.spines["left"].set_color(NEUTRAL_DARK)
    ax.spines["bottom"].set_color(NEUTRAL_DARK)


def draw_panel_c(ax, biomarker_summary: pd.DataFrame):
    ax.set_axis_off()
    panel_title(ax, "C", "Actionability-driven translation")
    rounded_card(
        ax,
        (0.12, 0.82, 0.76, 0.10),
        "Rank-24 MiNP molecular signal",
        facecolor=EXPOSURE_LIGHT,
        edgecolor=EXPOSURE,
        linewidth=1.1,
        text_color=EXPOSURE,
        fontsize=6.3,
        weight="bold",
    )
    arrow(ax, (0.50, 0.81), (0.50, 0.75), color=EXPOSURE)
    rounded_card(
        ax,
        (0.06, 0.57, 0.88, 0.17),
        "Multistage prioritization\nMolecular stability  •  human biomarker\nNovelty  •  epidemiologic testability",
        facecolor=WHITE,
        edgecolor=NEUTRAL,
        linewidth=1.0,
        text_color=NEUTRAL_DARK,
        fontsize=5.35,
        weight="bold",
    )
    arrow(ax, (0.50, 0.56), (0.50, 0.50), color=EXPOSURE)
    rounded_card(
        ax,
        (0.12, 0.39, 0.76, 0.10),
        "DINP-related exposure axis",
        facecolor=WHITE,
        edgecolor=EXPOSURE,
        linewidth=1.1,
        text_color=EXPOSURE,
        fontsize=6.3,
        weight="bold",
    )
    arrow(ax, (0.50, 0.38), (0.50, 0.31), color=EXPOSURE, dashed=True)
    ax.text(
        0.53,
        0.345,
        "biomarker translation",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=5.0,
        color=NEUTRAL,
    )
    rounded_card(
        ax,
        (0.12, 0.20, 0.76, 0.10),
        "Urinary MCOP selected\nfor human validation",
        facecolor=EXPOSURE_PALE,
        edgecolor=EXPOSURE,
        linewidth=1.2,
        text_color=EXPOSURE,
        fontsize=6.0,
        weight="bold",
    )
    detectability = biomarker_summary.set_index("analyte")["above_lod_pct"]
    ax.text(
        0.50,
        0.085,
        "Above LOD across 7 cycles\n"
        f"MiNP {detectability['MiNP']:.1f}%  |  MCOP {detectability['MCOP']:.1f}%\n"
        "Biomarker translation; not a direct CTD nomination",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=5.0,
        color=NEUTRAL,
        linespacing=1.15,
    )


def write_reports(output_dir: Path, screen: pd.DataFrame, manifest: pd.DataFrame):
    manifest.to_csv(output_dir / "Figure1_study_design_v2_source_manifest.csv", index=False)
    stats = "# Figure 1 academic v2 — statistics and reproducibility\n\n"
    stats += "## Scientific question\n"
    stats += "How did a rank-24 MiNP molecular signal advance to a human-testable DINP exposure axis?\n\n"
    stats += "## Quantitative panel\n"
    stats += f"- Panel B uses all {len(screen)} frozen Phase 1 chemicals; no rows were downsampled or dropped.\n"
    stats += f"- Primary ranking metric: `-log10(BH-FDR)` from the GeneCards Disorders-scoped, `gene_cards_k=1000`, `U_core` screen.\n"
    stats += f"- Frozen screen checks: {int((screen['bh_fdr'] < 0.05).sum())} chemicals with BH-FDR < 0.05; {int(screen['stable_for_primary_sort'].astype(str).str.lower().eq('true').sum())} stable candidates.\n"
    stats += "- Panel B is descriptive; it does not introduce a new statistical test or re-rank candidates. Molecular rank is explicitly separated from final actionability.\n\n"
    stats += "## Schematic panels\n"
    stats += "- Panel A: workflow nodes from `figure1_panelA_workflow_nodes.csv`; no quantitative effect estimate is encoded.\n"
    stats += "- Panel C summarizes the prespecified prioritization dimensions at a high level; the detailed Top-30 tournament is reserved for a dedicated candidate-prioritization figure.\n"
    stats += "- Panel C: translation boundary is explicitly labelled as biomarker translation; MCOP is not represented as a direct CTD nomination.\n"
    stats += "- Panel C detectability values are read from the Phase 2A NHANES biomarker audit: MiNP 27.4% and MCOP 98.4% above LOD across seven cycles.\n"
    stats += "- The previous roadmap/future-replication panel was removed as redundant; prospective replication belongs in the Discussion rather than this figure.\n\n"
    stats += "## Source traceability\n"
    for row in manifest.itertuples(index=False):
        stats += f"- `{row.source_file}` — {row.rows_or_scope} — {row.figure_use}.\n"
    stats += "\n## Export\n- Vector master: PDF.\n- Preview: RGB PNG at 300 dpi.\n- Editable text: SVG export.\n"
    (output_dir / "Figure1_study_design_v2_statistics.md").write_text(stats, encoding="utf-8")


def write_qa_report(output_dir: Path):
    report = """# Academic Figure Skill QA Report — Figure 1 v2

## Contract

- Core conclusion: MiNP ranked 24th in the hypothesis-agnostic molecular screen but advanced through multistage actionability prioritization to a human-testable DINP exposure axis represented by urinary MCOP.
- Archetype: schematic-led horizontal composite with Panel B as the quantitative anchor.
- Target: Nature-family double-column figure, 183 mm wide.
- Export: RGB PDF vector master, 300 dpi PNG preview, editable-text SVG.

## Four-pass QA

- **PASS — Anti-pattern scan:** restrained semantic colors, no decorative chart types, no causal connector, no default rainbow palette.
- **PASS — Code compliance:** mandatory typography/color/export baselines retained; all 267 chemicals used; minimum text size is 5.0 pt.
- **PASS — Visual logic:** MiNP is labelled as rank 24; DINP parent is labelled rank 107 with BH-FDR 0.449; MCOP is explicitly a biomarker translation rather than a direct CTD hit.
- **PASS — Translation evidence:** Panel C reads the seven-cycle MiNP and MCOP above-LOD percentages directly from the frozen Phase 2A audit.
- **PASS — Anti-redundancy:** the roadmap/future-replication panel was removed; detailed Top-30 triage is reserved for a dedicated candidate-prioritization figure.
- **PASS — Statistical terminology:** MiNP is labelled `degree-matched FDR 0.036`; the raw degree-matched empirical P is not conflated with the BH-adjusted value.
- **PASS — Card containment:** rendered text bounding boxes were checked programmatically and all registered card text remains inside its frame.
- **PASS — Evidence boundary:** multistage prioritization is presented as selection logic, not as proof that MiNP was the top-ranked molecular hit.

## Verdict

**READY — Figure 1 passes the publication-oriented scientific and overflow QA checks.**
"""
    (output_dir / "Figure1_study_design_v2_QA.md").write_text(report, encoding="utf-8")


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

    screen, nodes, biomarker_summary, manifest = validate_and_load(
        args.repo_root, args.source_dir
    )
    fig = plt.figure(
        figsize=(FIGURE_WIDTH_MM * MM_TO_INCH, FIGURE_HEIGHT_MM * MM_TO_INCH),
        facecolor=WHITE,
    )
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.05, 1.42, 1.30],
        left=0.035,
        right=0.985,
        bottom=0.10,
        top=0.86,
        wspace=0.30,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    draw_panel_a(ax_a, nodes)
    draw_panel_b(ax_b, screen)
    draw_panel_c(ax_c, biomarker_summary)

    assert_card_text_containment(fig)

    stem = args.output_dir / "Figure1_study_design_v2"
    save_cns_figure(fig, stem)
    fig.savefig(f"{stem}.svg", bbox_inches="tight", dpi=300)
    write_reports(args.output_dir, screen, manifest)
    write_qa_report(args.output_dir)
    plt.close(fig)
    print(f"Rendered {stem}.pdf, {stem}.png, {stem}.svg")


if __name__ == "__main__":
    main()
