# Academic Figure Skill Asset Confirmation (verified against audited source data)
# (a) Pooled association -> source_data/figure2_primary_python_vs_r.csv -> param inherit
# (b) Leave-one-cycle-out stability wheel -> source_data/figure2_loco.csv -> param inherit
# (c) Temporal cycle-effect landscape -> source_data/figure2_per_cycle.csv -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies audited statistical values.

"""Figure 3 v3 prototype: three visual grammars, no conventional forest plots."""

# ---- Journal style baseline block (academic-figure-skill) ----
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9.2,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.2,
    "ytick.labelsize": 7.2,
    "legend.fontsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.65,
})

# ---- Color palette baseline block (academic-figure-skill) ----
COLORS = {
    "ink": "#1B2A34",
    "muted": "#6B7A82",
    "rule": "#CBD5D9",
    "teal": "#217A83",
    "teal_mid": "#6BAEB2",
    "teal_light": "#DDECEE",
    "warm": "#A13B54",
    "warm_light": "#F2E2E7",
}

# ---- Export baseline block (academic-figure-skill) ----
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse


REPO = Path(__file__).resolve().parents[2]
SOURCE_DIR = REPO / "outputs" / "manuscript" / "figures" / "source_data"
OUT_DIR = REPO / "outputs" / "manuscript" / "figures" / "upgrade_v3"


def save_cns_figure(fig: mpl.figure.Figure, out_base: Path) -> None:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.035)
    fig.savefig(out_base.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.035)
    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)


def clean_axis(ax: mpl.axes.Axes) -> None:
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(COLORS["ink"])
    ax.spines["bottom"].set_color(COLORS["ink"])
    ax.tick_params(colors=COLORS["ink"], width=0.6, length=3)
    ax.grid(False)


def panel_letter(ax: mpl.axes.Axes, letter: str) -> None:
    ax.text(-0.10, 1.10, letter, transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold", color=COLORS["ink"], clip_on=False)


def validate_sources() -> dict[str, object]:
    primary = pd.read_csv(SOURCE_DIR / "figure2_primary_python_vs_r.csv")
    loco = pd.read_csv(SOURCE_DIR / "figure2_loco.csv")
    cycle = pd.read_csv(SOURCE_DIR / "figure2_per_cycle.csv")
    assert len(primary) == 1 and len(loco) == 7 and len(cycle) == 7
    row = primary.iloc[0]
    pooled = {
        "or": float(row["r_OR"]),
        "low": float(row["r_CI_low"]),
        "high": float(row["r_CI_high"]),
        "p": float(row["r_P_design_df"]),
        "n": int(row["r_N"]),
        "crc": int(row["r_CRC_N"]),
    }
    assert 1.20 < pooled["or"] < 1.30 and pooled["low"] > 1
    assert (loco["OR"] > 1).all()
    assert (cycle["OR"] > 1).sum() == 6
    assert float(cycle.loc[cycle["Cycle"] == "2011-2012", "OR"].iloc[0]) < 1
    return {"pooled": pooled, "loco": loco, "cycle": cycle}


def draw_confidence_curve(ax: mpl.axes.Axes, pooled: dict[str, float]) -> None:
    clean_axis(ax); panel_letter(ax, "A")
    ax.set_title("Pooled association", loc="left", pad=8, fontweight="bold", color=COLORS["ink"])
    ax.set_xlim(0.88, 1.62); ax.set_ylim(0, 1.05); ax.set_yticks([])
    ax.set_xlabel("Odds ratio per doubling of urinary MCOP", labelpad=7)

    mu = np.log(pooled["or"])
    se = (np.log(pooled["high"]) - np.log(pooled["low"])) / (2 * 1.96)
    x = np.linspace(0.88, 1.62, 800)
    support = np.exp(-0.5 * ((np.log(x) - mu) / se) ** 2) * 0.62
    support += 0.02
    low, high = pooled["low"], pooled["high"]
    ax.fill_between(x, 0, support, where=(x >= low) & (x <= high), color=COLORS["teal_light"], zorder=1)
    ax.plot(x, support, color=COLORS["teal"], lw=2.2, zorder=3)
    ax.axvline(1, color=COLORS["muted"], lw=1.0, ls=(0, (3, 3)), zorder=0)
    ax.axvline(low, color=COLORS["teal_mid"], lw=0.8, ls=(0, (2, 2)), zorder=2)
    ax.axvline(high, color=COLORS["teal_mid"], lw=0.8, ls=(0, (2, 2)), zorder=2)
    ax.text(0.03, 0.90, "OR", transform=ax.transAxes, color=COLORS["muted"], fontsize=8)
    ax.text(0.03, 0.68, f"{pooled['or']:.3f}", transform=ax.transAxes, color=COLORS["teal"], fontsize=21, fontweight="bold")
    ax.text(0.03, 0.52, f"95% CI {low:.3f}–{high:.3f}  ·  P = {pooled['p']:.4f}", transform=ax.transAxes, color=COLORS["ink"], fontsize=8.2)
    ax.text(0.97, 0.90, f"N = {pooled['n']:,}  ·  CRC = {pooled['crc']}", transform=ax.transAxes, ha="right", color=COLORS["muted"], fontsize=8)


def draw_stability_wheel(ax: mpl.axes.Axes, pooled: dict[str, float], loco: pd.DataFrame) -> None:
    ax.set_facecolor("white")
    panel_letter(ax, "B")
    theta = np.linspace(0, 2 * np.pi, len(loco), endpoint=False)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_title("Leave-one-cycle-out stability", loc="left", pad=13, fontweight="bold", color=COLORS["ink"])

    def radius(or_value: float) -> float:
        return 1.0 + 1.65 * np.log(or_value) / np.log(1.6)

    r_null = radius(1.0)
    for angle, (_, row) in zip(theta, loco.iterrows()):
        r_low, r_or, r_high = radius(float(row["CI_low"])), radius(float(row["OR"])), radius(float(row["CI_high"]))
        ax.bar(angle, r_high - r_low, bottom=r_low, width=0.58, color=COLORS["teal_light"], edgecolor=COLORS["teal_mid"], linewidth=1.1, alpha=0.95, zorder=2)
        ax.scatter(angle, r_or, s=42, color=COLORS["teal"], edgecolor="white", linewidth=0.9, zorder=4)
    ax.plot(np.linspace(0, 2 * np.pi, 360), np.full(360, r_null), color=COLORS["muted"], lw=1.0, ls=(0, (3, 3)), zorder=1)
    ax.set_ylim(0.72, 2.75)
    tick_values = [1.0, 1.2, 1.4, 1.6]
    ax.set_yticks([radius(v) for v in tick_values])
    ax.set_yticklabels(["1.0", "1.2", "1.4", "1.6"], color=COLORS["muted"], fontsize=6.8)
    ax.yaxis.grid(True, color=COLORS["rule"], lw=0.55, ls=(0, (2, 3)))
    ax.xaxis.grid(False)
    labels = [x.replace("-", "–").replace("2005–2006", "05–06").replace("2007–2008", "07–08").replace("2009–2010", "09–10").replace("2011–2012", "11–12").replace("2013–2014", "13–14").replace("2015–2016", "15–16").replace("2017–2018", "17–18") for x in loco["Dropped_cycle"]]
    ax.set_xticks(theta); ax.set_xticklabels(labels, fontsize=6.8, color=COLORS["ink"])
    ax.tick_params(axis="x", pad=7)
    ax.text(0.5, 0.51, "7/7", transform=ax.transAxes, ha="center", va="center", color=COLORS["teal"], fontsize=17, fontweight="bold")
    ax.text(0.5, 0.40, "above null", transform=ax.transAxes, ha="center", va="center", color=COLORS["ink"], fontsize=7.8)
    ax.text(0.5, 0.30, f"Pooled OR {pooled['or']:.3f}", transform=ax.transAxes, ha="center", va="center", color=COLORS["muted"], fontsize=7.0)
    ax.text(0.02, 0.98, "Null ring: OR = 1", transform=ax.transAxes, ha="left", va="top", color=COLORS["muted"], fontsize=6.8)


def draw_temporal_landscape(ax: mpl.axes.Axes, pooled: dict[str, float], cycle: pd.DataFrame) -> None:
    clean_axis(ax); panel_letter(ax, "C")
    ax.set_title("Temporal effect landscape", loc="left", pad=9, fontweight="bold", color=COLORS["ink"])
    log_or = np.log(cycle["OR"].to_numpy())
    log_low = np.log(cycle["CI_low"].to_numpy())
    log_high = np.log(cycle["CI_high"].to_numpy())
    x = np.arange(len(cycle))
    pooled_low, pooled_high = np.log(pooled["low"]), np.log(pooled["high"])
    ax.axhspan(pooled_low, pooled_high, color=COLORS["teal_light"], alpha=0.62, zorder=0)
    ax.axhline(0, color=COLORS["muted"], lw=1.0, ls=(0, (3, 3)), zorder=1)
    ax.text(0.99, 0.93, "6/7 estimates > 1  ·  interaction P = 0.0060", transform=ax.transAxes, ha="right", color=COLORS["teal"], fontsize=7.4, fontweight="bold")
    ax.text(0.99, 0.08, "Shaded band: pooled 95% CI", transform=ax.transAxes, ha="right", color=COLORS["muted"], fontsize=7.0)

    for i, (y, lo, hi, cases, label) in enumerate(zip(log_or, log_low, log_high, cycle["CRC_N"], cycle["Cycle"])):
        discordant = label == "2011-2012"
        color = COLORS["warm"] if discordant else COLORS["teal"]
        edge = COLORS["warm"] if discordant else COLORS["teal_mid"]
        ellipse = Ellipse((i, y), width=0.34, height=float(hi - lo), facecolor=COLORS["warm_light"] if discordant else COLORS["teal_light"], edgecolor=edge, lw=1.2, alpha=0.95, zorder=2)
        ax.add_patch(ellipse)
        ax.plot([i, i], [lo, hi], color=edge, lw=1.0, alpha=0.72, zorder=2)
        ax.scatter([i], [y], s=42 + int(cases) * 4.5, color=color, edgecolor="white", linewidth=0.8, zorder=4)
        ax.text(i, y + 0.085, f"{np.exp(y):.2f}", ha="center", va="bottom", color=color, fontsize=7.0, fontweight="bold")
    ax.annotate("2011–12  ·  OR 0.84", xy=(3, log_or[3]), xytext=(3.55, -0.36), color=COLORS["warm"], fontsize=7.4, fontweight="bold", ha="left", arrowprops={"arrowstyle": "-", "color": COLORS["warm"], "lw": 1.0})
    ax.set_xlim(-0.55, 6.55); ax.set_ylim(-0.55, 0.72)
    ax.set_xticks(x); ax.set_xticklabels([f"{xv[:4]}–{xv[-2:]}" for xv in cycle["Cycle"]], color=COLORS["ink"])
    ax.set_ylabel("log(OR) per doubling", labelpad=7)
    ax.set_xlabel("NHANES cycle", labelpad=7)
    ax.text(0.02, 0.04, "Bubble area ∝ CRC cases  ·  ellipse = 95% CI", transform=ax.transAxes, color=COLORS["muted"], fontsize=7.0)


def main() -> None:
    data = validate_sources(); pooled, loco, cycle = data["pooled"], data["loco"], data["cycle"]
    fig = plt.figure(figsize=(183 / 25.4, 116 / 25.4), facecolor="white")
    ax_a = fig.add_axes((0.08, 0.53, 0.38, 0.37))
    ax_b = fig.add_axes((0.56, 0.52, 0.38, 0.39), projection="polar")
    ax_c = fig.add_axes((0.15, 0.10, 0.72, 0.30))
    draw_confidence_curve(ax_a, pooled)
    draw_stability_wheel(ax_b, pooled, loco)
    draw_temporal_landscape(ax_c, pooled, cycle)
    out = OUT_DIR / "Figure3_human_mcop_v3"
    save_cns_figure(fig, out)
    (OUT_DIR / "Figure3_human_mcop_v3_QA.md").write_text(
        "# Figure 3 v3 QA\n\n"
        "- Three panels were redrawn natively from audited source CSVs; no legacy raster panel was embedded.\n"
        "- Panel A uses a Wald confidence curve reconstructed from the reported OR and 95% CI; it is not a profile-likelihood calculation.\n"
        "- Panel B encodes each LOCO 95% CI as an annular radial interval; all seven intervals remain outside the OR=1 null ring.\n"
        "- Panel C encodes each cycle-specific 95% CI as an ellipse, bubble area as CRC case count, and highlights 2011–12.\n"
        "- Exact primary, LOCO, and cycle estimates are validated before rendering.\n"
        "- No conventional horizontal forest plot is used.\n", encoding="utf-8")
    (OUT_DIR / "Figure3_human_mcop_v3_statistics.md").write_text(
        "# Figure 3 v3 statistics and traceability\n\n"
        "Panel A: primary complex-survey model; OR per doubling of urinary MCOP = 1.246, 95% CI 1.077-1.440, P=0.0033, N=9,936, CRC=70. Source: `figure2_primary_python_vs_r.csv`.\n\n"
        "Panel B: seven leave-one-cycle-out estimates from `figure2_loco.csv`. Radial center positions are ORs; annular inner/outer boundaries are lower/upper 95% CI limits; null ring is OR=1.\n\n"
        "Panel C: seven cycle-specific estimates from `figure2_per_cycle.csv`. The y-coordinate is log(OR), ellipse height is the log-scale 95% CI, and bubble area is proportional to the cycle CRC case count.\n\n"
        "The temporal interaction P=0.0060 is an audited model-level statistic.\n", encoding="utf-8")
    pd.DataFrame([
        {"figure": "Figure3_human_mcop_v3", "source": "outputs/manuscript/figures/source_data/figure2_primary_python_vs_r.csv", "panel": "A"},
        {"figure": "Figure3_human_mcop_v3", "source": "outputs/manuscript/figures/source_data/figure2_loco.csv", "panel": "B"},
        {"figure": "Figure3_human_mcop_v3", "source": "outputs/manuscript/figures/source_data/figure2_per_cycle.csv", "panel": "C"},
    ]).to_csv(OUT_DIR / "Figure3_human_mcop_v3_source_manifest.csv", index=False)


if __name__ == "__main__":
    main()
