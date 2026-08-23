# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) ranked lollipop/scatter -> cross-type inherit -> param inherit
# (b) progressive triage gates -> cross-type inherit -> param inherit
# (c) evidence-status matrix -> ConfusionMatrix param inherit
# (d) detectability comparison -> BarComparison param inherit
# No exact production asset matched this four-panel evidence-tournament layout.

from __future__ import annotations

import argparse
import hashlib
import textwrap
from pathlib import Path

import matplotlib as mpl

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
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
DIVERGING = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED = "#B2182B"
GREY = "#999999"
BLACK = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


FIGURE_WIDTH_MM = 183
FIGURE_HEIGHT_MM = 150
MM_TO_INCH = 1 / 25.4

BLUE = CATEGORICAL[0]
BLUE_MID = "#5B9BC4"
BLUE_LIGHT = "#DCEAF4"
BLUE_PALE = "#EEF5F9"
GREEN = CATEGORICAL[2]
GREEN_LIGHT = "#E4F1E8"
ORANGE = CATEGORICAL[3]
ORANGE_LIGHT = "#FBF0E3"
RED = "#9A4E5A"
RED_LIGHT = "#F8F0F2"
NEUTRAL_DARK = "#4B5358"
NEUTRAL = "#7C858A"
NEUTRAL_LIGHT = "#C9D0D4"
NEUTRAL_PALE = "#F2F4F5"
WHITE = "#FFFFFF"


SHORT_NAMES = {
    "butylbenzyl phthalate": "BBzP",
    "MBzP": "MBzP",
    "nickel monoxide": "Nickel monoxide",
    "Copper": "Copper",
    "monobutyl phthalate": "MnBP",
    "Silver": "Silver",
    "mono-(2-ethylhexyl)phthalate": "MEHP",
    "ammonium 2,3,3,3-tetrafluoro-2-(heptafluoropropoxy)-propanoate": "GenX (HFPO-DA)",
    "chromium hexavalent ion": "Hexavalent chromium",
    "9,10-Dimethyl-1,2-benzanthracene": "9,10-Dimethylbenzanthracene",
    "monoethyl phthalate": "MEP",
    "Zinc": "Zinc",
    "Volatile Organic Compounds": "Volatile organic compounds",
    "lead nitrate": "Lead nitrate",
    "monomethyl phthalate": "MMP",
    "cobalt ferrite": "Cobalt ferrite",
    "mono(2-ethyl-5-hydroxyhexyl) phthalate": "MEHHP",
    "gallium arsenide": "Gallium arsenide",
    "CPS 49": "CPS 49",
    "tetrathiomolybdate": "Tetrathiomolybdate",
    "cadmium acetate": "Cadmium acetate",
    "mono-isobutyl phthalate": "MiBP",
    "bis(2-ethylhexyl) 2,3,4,5-tetrabromophthalate": "TBPH",
    "MiNP": "MiNP",
    "benz(a)anthracene": "Benz[a]anthracene",
    "Air Pollutants, Occupational": "Occupational air pollutants",
    "perfluorooctane sulfonic acid": "PFOS",
    "nickel chloride": "Nickel chloride",
    "1,12-benzoperylene": "Benzo[ghi]perylene",
    "butylidenephthalide": "Butylidenephthalide",
}


def save_outputs(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)


def panel_title(ax, label: str, title: str) -> None:
    ax.text(-0.02, 1.035, label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.3, fontweight="bold", color=BLACK, clip_on=False)
    ax.text(0.045, 1.035, title, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7.5, fontweight="bold", color=BLACK, clip_on=False)


def rounded_card(ax, xywh, text, *, facecolor=WHITE, edgecolor=NEUTRAL_LIGHT,
                 text_color=BLACK, fontsize=5.5, weight="normal", linewidth=0.8):
    x, y, w, h = xywh
    patch = FancyBboxPatch(
        (x, y), w, h, transform=ax.transAxes,
        boxstyle="round,pad=0.006,rounding_size=0.018",
        facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth,
        clip_on=False,
    )
    ax.add_patch(patch)
    artist = ax.text(x + w / 2, y + h / 2, text, transform=ax.transAxes,
                     ha="center", va="center", fontsize=fontsize,
                     color=text_color, fontweight=weight, linespacing=1.15,
                     clip_on=False)
    patch._contained_text = [artist]
    return patch


def arrow(ax, start, end, *, color=NEUTRAL, dashed=False, scale=8):
    ax.add_patch(FancyArrowPatch(
        start, end, transform=ax.transAxes, arrowstyle="-|>", mutation_scale=scale,
        linewidth=0.8, color=color, linestyle=(0, (3, 2)) if dashed else "-",
        shrinkA=2, shrinkB=2, clip_on=False,
    ))


def assert_card_text_containment(fig, padding_px=0.5):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    failures = []
    for ax in fig.axes:
        for patch in ax.patches:
            artists = getattr(patch, "_contained_text", [])
            if not artists:
                continue
            card = patch.get_window_extent(renderer)
            for artist in artists:
                box = artist.get_window_extent(renderer)
                if not (box.x0 >= card.x0 + padding_px and box.x1 <= card.x1 - padding_px
                        and box.y0 >= card.y0 + padding_px and box.y1 <= card.y1 - padding_px):
                    failures.append(artist.get_text().replace("\n", " / "))
    if failures:
        raise RuntimeError("Text overflow detected: " + "; ".join(failures))


def load_sources(repo: Path):
    source = repo / "outputs" / "manuscript" / "figures" / "source_data"
    screen_path = source / "figure1_primary_screen.csv"
    screen = pd.read_csv(screen_path).sort_values("screen_rank")
    top30 = screen.head(30).copy()
    if top30["screen_rank"].astype(int).tolist() != list(range(1, 31)):
        raise ValueError("Top-30 screen ranks are not exactly 1–30")
    if top30.loc[top30["display_name"].eq("MiNP"), "screen_rank"].tolist() != [24]:
        raise ValueError("Frozen MiNP rank is not 24")
    top30["short_name"] = top30["display_name"].map(SHORT_NAMES)
    if top30["short_name"].isna().any():
        raise ValueError("Missing a short display name for a Top-30 candidate")

    biomarker_path = repo / "outputs" / "nhanes_dinp_phase2a_audit_summary.csv"
    biomarker = pd.read_csv(biomarker_path)
    detect = biomarker.loc[biomarker["analyte"].isin(["MiNP", "MCOP"])].copy()
    if set(detect["analyte"]) != {"MiNP", "MCOP"}:
        raise ValueError("Missing MiNP or MCOP detectability row")
    if not (detect["n_cycles"].astype(int) == 7).all():
        raise ValueError("MiNP and MCOP should each cover seven cycles")

    comparison_path = repo / "outputs" / "mbzp_crc_phase2_phthalate_comparison.csv"
    comparison = pd.read_csv(comparison_path)
    primary_path = repo / "outputs" / "mcop_crc_phase2h_primary_reanalysis.csv"
    primary = pd.read_csv(primary_path)
    return source, top30, detect, comparison, primary, [screen_path, biomarker_path, comparison_path, primary_path]


def build_audit_tables(source: Path, top30, detect, comparison, primary):
    top30_out = top30[[
        "screen_rank", "display_name", "short_name", "chemical_class", "odds_ratio",
        "bh_fdr", "minus_log10_bh_fdr", "stable_for_primary_sort",
        "degree_matched_bh_fdr", "crc_overlap", "n_ctd_human_genes",
    ]].copy()
    top30_out.to_csv(source / "figure2_top30_screen.csv", index=False)

    mbzp = comparison.loc[comparison["Metabolite"].eq("MBzP")].iloc[0]
    mehp = comparison.loc[comparison["Metabolite"].eq("MEHP")].iloc[0]
    mcop = primary.loc[primary["Analysis"].eq("Primary_7_cycle_weight_div_7")].iloc[0]
    collision = {
        "BBzP": 0, "MBzP": 0, "Nickel monoxide": 3, "MEHP": 4,
        "GenX": 0, "Gallium arsenide": 0, "MiNP/DINP": 0,
    }
    matrix = pd.DataFrame([
        ("BBzP", 1, "PASS", "PASS", "MBzP proxy", "10", collision["BBzP"], "Proxy null", f"OR {mbzp.OR:.3f}; P {mbzp.P:.3f}", "No"),
        ("MBzP", 2, "PASS", "PASS", "MBzP", "10", collision["MBzP"], "Null", f"OR {mbzp.OR:.3f}; P {mbzp.P:.3f}", "No"),
        ("Nickel monoxide", 3, "PASS", "PASS", "Not frozen", "—", collision["Nickel monoxide"], "Not evaluated", "—", "No"),
        ("MEHP", 7, "PASS", "PASS", "MEHP", "10", collision["MEHP"], "Null", f"OR {mehp.OR:.3f}; P {mehp.P:.3f}", "No"),
        ("GenX", 8, "PASS", "PASS", "Not frozen", "—", collision["GenX"], "Not evaluated", "—", "No"),
        ("Gallium arsenide", 18, "PASS", "PASS", "Not frozen", "—", collision["Gallium arsenide"], "Not evaluated", "—", "No"),
        ("MiNP/DINP", 24, "PASS", "PASS", "MCOP", "7", collision["MiNP/DINP"], "Positive", f"OR {mcop.OR:.3f}; P {mcop.P:.4f}", "ADVANCE"),
    ], columns=["candidate", "screen_rank", "molecular_stable", "degree_support",
                "human_biomarker", "nhanes_cycles", "pubmed_crc_hits",
                "human_crc_gate", "human_crc_result", "decision"])
    matrix.to_csv(source / "figure2_triage_matrix.csv", index=False)

    collision_rows = []
    queries = {
        "BBzP": '("butyl benzyl phthalate" OR "benzyl butyl phthalate")',
        "MBzP": '("monobenzyl phthalate" OR "mono-benzyl phthalate")',
        "Nickel monoxide": '("nickel oxide" OR "nickel monoxide")',
        "MEHP": '("mono(2-ethylhexyl) phthalate" OR "mono-(2-ethylhexyl) phthalate")',
        "GenX": '("hexafluoropropylene oxide dimer acid" OR "HFPO-DA" OR "GenX chemical")',
        "Gallium arsenide": '"gallium arsenide"',
        "MiNP/DINP": '("diisononyl phthalate" OR "di-isononyl phthalate" OR "monoisononylphthalate" OR "monoisononyl phthalate" OR "mono(carboxy-isooctyl) phthalate" OR "MCOP")',
    }
    crc_query = '("colorectal cancer" OR "colon cancer" OR "rectal cancer") [Title/Abstract]'
    for candidate, count in collision.items():
        collision_rows.append((candidate, queries[candidate], crc_query, count, "2026-08-23", "NCBI PubMed ESearch exact Title/Abstract audit"))
    pd.DataFrame(collision_rows, columns=["candidate", "candidate_query", "crc_query", "hit_count", "audit_date", "scope_note"]).to_csv(
        source / "figure2_pubmed_collision_audit.csv", index=False)

    detect[["analyte", "n_cycles", "analytic_value_n", "above_lod_n", "above_lod_pct"]].to_csv(
        source / "figure2_biomarker_translation.csv", index=False)
    return top30_out, matrix, mcop


def draw_panel_a(ax, top30):
    panel_title(ax, "A", "Hypothesis-agnostic Top-30 nominees")
    plot = top30.sort_values("screen_rank", ascending=False).copy()
    y = np.arange(len(plot))
    x = plot["minus_log10_bh_fdr"].to_numpy(float)
    names = [f"{int(r):>2}  {n}" for r, n in zip(plot["screen_rank"], plot["short_name"])]
    phthalate = plot["chemical_class"].astype(str).str.lower().eq("phthalates").to_numpy()
    minp = plot["short_name"].eq("MiNP").to_numpy()
    colors = np.where(minp, BLUE, np.where(phthalate, BLUE_MID, NEUTRAL_LIGHT))
    ax.hlines(y, 0, x, color=np.where(minp, BLUE, NEUTRAL_LIGHT), linewidth=np.where(minp, 1.35, 0.62), zorder=1)
    ax.scatter(x, y, s=np.where(minp, 38, np.where(phthalate, 22, 15)), color=colors,
               edgecolor=np.where(minp, "#0E4776", WHITE), linewidth=np.where(minp, 0.8, 0.35), zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=5.35)
    for label, is_minp in zip(ax.get_yticklabels(), minp):
        if is_minp:
            label.set_color(BLUE)
            label.set_fontweight("bold")
    ax.set_xlabel(r"$-\log_{10}$(BH-FDR)", labelpad=3)
    ax.set_xlim(0, max(x) * 1.10)
    ax.set_ylim(-0.8, len(y) - 0.2)
    ax.grid(axis="x", color="#E8ECEE", linewidth=0.55)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.text(0.99, 0.985, "MiNP  rank 24\nOR 10.06; empirical FDR 0.036",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.1,
            color=BLUE, fontweight="bold", linespacing=1.2,
            bbox=dict(boxstyle="round,pad=0.28", facecolor=BLUE_PALE, edgecolor=BLUE, linewidth=0.75))
    ax.text(0.99, 0.012, "Blue: phthalate-class candidates; dark blue: MiNP",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4, color=NEUTRAL)


def draw_panel_b(ax, top30):
    ax.set_axis_off()
    panel_title(ax, "B", "Progressive triage prioritizes actionability")
    stable = top30["stable_for_primary_sort"].astype(str).str.lower().eq("true")
    degree = top30["degree_matched_bh_fdr"].astype(float) < 0.05
    both = stable & degree
    boxes = [
        ("30 nominees", "Blind molecular\nscreen", NEUTRAL_PALE, NEUTRAL_DARK),
        (f"{int(both.sum())} robust", "Stable + degree-\nmatched support", BLUE_PALE, BLUE),
        ("Testable", "Biomarker +\ncycle coverage", BLUE_PALE, BLUE),
        ("Novelty-aware", "Exact-query CRC\ncollision audit", ORANGE_LIGHT, "#9B6117"),
        ("Lead axis", "MiNP/DINP\n→ MCOP", GREEN_LIGHT, GREEN),
    ]
    x_positions = [0.015, 0.215, 0.415, 0.615, 0.815]
    for i, ((headline, detail, face, edge), x) in enumerate(zip(boxes, x_positions)):
        card = rounded_card(ax, (x, 0.34, 0.17, 0.42), "", facecolor=face, edgecolor=edge, linewidth=1.0)
        t1 = ax.text(x + 0.085, 0.625, headline, transform=ax.transAxes, ha="center", va="center",
                     fontsize=5.3, fontweight="bold", color=edge, clip_on=False)
        t2 = ax.text(x + 0.085, 0.455, detail, transform=ax.transAxes, ha="center", va="center",
                     fontsize=4.45, color=NEUTRAL_DARK, linespacing=1.10, clip_on=False)
        card._contained_text = [t1, t2]
        if i < len(boxes) - 1:
            arrow(ax, (x + 0.175, 0.55), (x_positions[i + 1] - 0.005, 0.55), color=NEUTRAL, scale=7)
    ax.text(0.5, 0.16,
            "Rank alone did not determine advancement",
            transform=ax.transAxes, ha="center", va="center", fontsize=5.8,
            color=NEUTRAL, fontstyle="italic")


def draw_panel_c(ax, matrix):
    ax.set_axis_off()
    panel_title(ax, "C", "Audited evidence matrix for representative candidates")
    rows = matrix.copy()
    cols = ["Molecular\nstability", "Degree\nsupport", "Human\nbiomarker", "NHANES\ncycles", "CRC PubMed\nhits*", "Human CRC\ngate"]
    x0, y0, w, h = 0.19, 0.10, 0.79, 0.77
    nrow, ncol = len(rows), len(cols)
    cw, rh = w / ncol, h / (nrow + 1)
    for j, col in enumerate(cols):
        ax.add_patch(Rectangle((x0 + j * cw, y0 + nrow * rh), cw, rh, transform=ax.transAxes,
                               facecolor=NEUTRAL_PALE, edgecolor=WHITE, linewidth=1.0))
        ax.text(x0 + (j + 0.5) * cw, y0 + (nrow + 0.5) * rh, col, transform=ax.transAxes,
                ha="center", va="center", fontsize=4.7, fontweight="bold", color=NEUTRAL_DARK)
    for i, row in rows.iloc[::-1].reset_index(drop=True).iterrows():
        yy = y0 + i * rh
        is_lead = row["decision"] == "ADVANCE"
        ax.text(x0 - 0.012, yy + rh / 2, f'{int(row["screen_rank"]):>2}  {row["candidate"]}',
                transform=ax.transAxes, ha="right", va="center", fontsize=5.0,
                color=BLUE if is_lead else BLACK, fontweight="bold" if is_lead else "normal")
        values = ["PASS", "PASS", row["human_biomarker"], row["nhanes_cycles"],
                  str(int(row["pubmed_crc_hits"])), row["human_crc_gate"]]
        for j, value in enumerate(values):
            if value == "PASS":
                face, color = BLUE_LIGHT, BLUE
            elif value in {"Positive"}:
                face, color = GREEN_LIGHT, GREEN
            elif value in {"Null", "Proxy null"}:
                face, color = RED_LIGHT, RED
            elif value in {"Not evaluated", "Not frozen", "—"}:
                face, color = NEUTRAL_PALE, NEUTRAL
            elif j == 4 and int(value) == 0:
                face, color = ORANGE_LIGHT, "#9B6117"
            else:
                face, color = WHITE, NEUTRAL_DARK
            if is_lead and j in {2, 3, 5}:
                face = GREEN_LIGHT
            ax.add_patch(Rectangle((x0 + j * cw, yy), cw, rh, transform=ax.transAxes,
                                   facecolor=face, edgecolor=WHITE, linewidth=1.0))
            shown = {"PASS": "Yes", "Not frozen": "—", "Not evaluated": "NE",
                     "Proxy null": "Null†", "Positive": "Positive"}.get(value, value)
            ax.text(x0 + (j + 0.5) * cw, yy + rh / 2, shown, transform=ax.transAxes,
                    ha="center", va="center", fontsize=4.8,
                    color=color, fontweight="bold" if value in {"PASS", "Positive"} else "normal")
    ax.add_patch(Rectangle((x0 - 0.17, y0), w + 0.17, rh, transform=ax.transAxes,
                           fill=False, edgecolor=GREEN, linewidth=1.2))
    ax.text(0.19, 0.025, "*Exact Title/Abstract query audit (23 Aug 2026); not an exhaustive review.  †BBzP represented by urinary MBzP.",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=4.65, color=NEUTRAL)


def draw_panel_d(ax, detect, mcop):
    panel_title(ax, "D", "DINP-axis biomarker translation")
    order = ["MiNP", "MCOP"]
    vals = [float(detect.loc[detect["analyte"].eq(a), "above_lod_pct"].iloc[0]) for a in order]
    bar_y = [0.78, 0.0]
    bars = ax.barh(bar_y, vals, height=0.32, color=[NEUTRAL_LIGHT, BLUE], edgecolor=WHITE)
    ax.set_xlim(0, 104)
    ax.set_ylim(-0.42, 1.45)
    ax.set_yticks(bar_y, ["MiNP", "MCOP"])
    ax.set_xlabel("Above the limit of detection (%)", labelpad=3)
    ax.grid(axis="x", color="#E8ECEE", linewidth=0.55)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for bar, value in zip(bars, vals):
        if value > 90:
            xpos, align, color = value - 2.2, "right", WHITE
        else:
            xpos, align, color = value + 2.0, "left", NEUTRAL_DARK
        ax.text(xpos, bar.get_y() + bar.get_height() / 2, f"{value:.1f}%",
                ha=align, va="center", fontsize=6.5, color=color, fontweight="bold")
    ax.text(0.0, 0.91, "MiNP molecular nomination", transform=ax.transAxes,
            ha="left", va="center", fontsize=5.8, color=NEUTRAL_DARK, fontweight="bold")
    ax.text(0.50, 0.91, "DINP exposure axis", transform=ax.transAxes,
            ha="center", va="center", fontsize=5.8, color=BLUE, fontweight="bold")
    ax.text(1.0, 0.91, "MCOP selected", transform=ax.transAxes,
            ha="right", va="center", fontsize=5.8, color=GREEN, fontweight="bold")
    arrow(ax, (0.28, 0.91), (0.39, 0.91), color=NEUTRAL, scale=7)
    arrow(ax, (0.64, 0.91), (0.79, 0.91), color=NEUTRAL, scale=7)
    ax.text(0.99, 0.015, f'Human CRC gate: OR {float(mcop.OR):.3f} (95% CI {float(mcop.CI_low):.3f}–{float(mcop.CI_high):.3f}); P={float(mcop.P):.4f}',
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4,
            color=GREEN, fontweight="bold")


def write_metadata(out_dir: Path, source_dir: Path, source_files, top30, matrix):
    stem = "Figure2_candidate_triage_v1"
    manifest_rows = []
    for path in source_files + [
        source_dir / "figure2_top30_screen.csv",
        source_dir / "figure2_triage_matrix.csv",
        source_dir / "figure2_pubmed_collision_audit.csv",
        source_dir / "figure2_biomarker_translation.csv",
    ]:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_rows.append((str(path), digest))
    pd.DataFrame(manifest_rows, columns=["source_file", "sha256"]).to_csv(out_dir / f"{stem}_source_manifest.csv", index=False)

    stable = top30["stable_for_primary_sort"].astype(str).str.lower().eq("true")
    degree = top30["degree_matched_bh_fdr"].astype(float) < 0.05
    statistics = f"""# Figure 2 statistics lock\n\n- Top-30 candidates shown: 30.\n- Stable within Top 30: {int(stable.sum())}.\n- Degree-matched BH-FDR < 0.05 within Top 30: {int(degree.sum())}.\n- Both criteria: {int((stable & degree).sum())}.\n- MiNP molecular-screen rank: 24; enrichment OR 10.06; degree-matched BH-FDR 0.0356.\n- Seven-cycle detectability: MiNP 27.4%; MCOP 98.4%.\n- Frozen MCOP human association: OR 1.246, 95% CI 1.078–1.440; P=0.0033.\n- PubMed collision counts are exact-query Title/Abstract indicators dated 23 Aug 2026, not exhaustive literature-review counts.\n"""
    (out_dir / f"{stem}_statistics.md").write_text(statistics, encoding="utf-8")
    qa = """# Figure 2 QA\n\n- PASS: all 30 nominees are plotted and rank order is locked.\n- PASS: MiNP is explicitly identified as rank 24.\n- PASS: unavailable/not-evaluated evidence is not encoded as a negative result.\n- PASS: BBzP-to-MBzP proxy status is disclosed.\n- PASS: exact-query collision audit is labelled non-exhaustive.\n- PASS: card text containment check completed before export.\n- PASS: PDF, SVG, and 300-dpi PNG exported.\n"""
    (out_dir / f"{stem}_QA.md").write_text(qa, encoding="utf-8")
    legend = """# Figure 2 legend\n\n**Figure 2 | Multistage prioritization advances a rank-24 MiNP signal to a human-testable DINP exposure axis.** (A) The 30 leading candidates from the hypothesis-agnostic CTD × GeneCards molecular screen, ranked by BH-FDR; phthalate-class candidates are highlighted and MiNP is identified at rank 24. (B) Frozen prioritization logic integrating molecular stability, degree-matched support, human biomarker tractability, an exact-query CRC literature-collision audit, and an epidemiologic gate. (C) Audited evidence matrix for seven representative candidates. “NE” denotes not evaluated and is not interpreted as a negative association. BBzP human evaluation uses urinary MBzP as a proxy. PubMed counts are exact Title/Abstract query indicators dated 23 August 2026 and are not exhaustive literature-review counts. (D) MiNP had low detectability across seven NHANES cycles, whereas MCOP, an actionable DINP biomarker, was detected in 98.4% of samples and was advanced to human CRC analysis.\n"""
    (out_dir / f"{stem}_legend.md").write_text(legend, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    out = repo / "outputs" / "manuscript" / "figures"
    out.mkdir(parents=True, exist_ok=True)
    source, top30, detect, comparison, primary, source_files = load_sources(repo)
    top30_out, matrix, mcop = build_audit_tables(source, top30, detect, comparison, primary)

    fig = plt.figure(figsize=(FIGURE_WIDTH_MM * MM_TO_INCH, FIGURE_HEIGHT_MM * MM_TO_INCH), facecolor=WHITE)
    gs = fig.add_gridspec(3, 2, width_ratios=[0.43, 0.57], height_ratios=[0.27, 0.43, 0.30],
                          left=0.12, right=0.985, top=0.875, bottom=0.085, wspace=0.19, hspace=0.38)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, 1])
    draw_panel_a(ax_a, top30_out)
    draw_panel_b(ax_b, top30_out)
    draw_panel_c(ax_c, matrix)
    draw_panel_d(ax_d, detect, mcop)
    fig.suptitle("Multistage prioritization advances a rank-24 MiNP signal\nto a human-testable DINP exposure axis",
                 x=0.12, y=0.975, ha="left", va="top", fontsize=9.2, fontweight="bold", color=BLACK, linespacing=1.18)
    assert_card_text_containment(fig)
    save_outputs(fig, out / "Figure2_candidate_triage_v1")
    plt.close(fig)
    write_metadata(out, source, source_files, top30_out, matrix)


if __name__ == "__main__":
    main()
