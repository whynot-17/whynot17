"""SLC7A11 virtual knockout redo in manually curated COAD models.

This script is intentionally analysis-only: it does not download data and it
does not run until invoked explicitly by the user.

Question
--------
Does SLC7A11 CRISPR/Chronos gene effect differ between manually curated
right- and left-sided colorectal adenocarcinoma (COAD) models?  Secondary
analyses ask whether transcriptional PUFA-pressure and AA-routing proxies are
associated with SLC7A11 dependency differently by side.

Interpretation of the dependency score
---------------------------------------
DepMap Chronos gene effect is more negative when knockout causes a larger
fitness loss.  Therefore, a negative Right-minus-Left contrast is compatible
with stronger SLC7A11 dependency in right-sided models.

Inputs are expected to be local files:
    work/depmap_Model.csv
    work/data/CRISPRGeneEffect_26Q1.csv
    work/depmap_expression_24Q4.csv

Example invocation (run manually after inspection):
    python work/depmap_slc7a11_virtual_ko_redo.py

The manual sidedness map is deliberately embedded below so the script remains
portable.  Only explicit anatomical primary-site labels are used in the
primary analysis; published but less certain assignments are sensitivity-only.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, spearmanr, t


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL_FILE = ROOT / "work" / "depmap_Model.csv"
DEFAULT_EFFECT_FILE = ROOT / "work" / "data" / "CRISPRGeneEffect_26Q1.csv"
DEFAULT_EXPRESSION_FILE = ROOT / "work" / "depmap_expression_24Q4.csv"
DEFAULT_OUT_DIR = ROOT / "outputs"

GENE = "SLC7A11"
PRESSURE_GENES = ("ACSL4", "LPCAT3", "ALOX5", "ALOX12", "ALOX15")
AA_ROUTING_GENES = ("PLA2G4A", "PTGS2", "ALOX5", "ALOX5AP")
PUFA_MIN_AVAILABLE = 4
AA_ROUTING_MIN_AVAILABLE = 3
SEED = 20260901


# The map follows the existing project curation rules:
# right = cecum/ascending/hepatic flexure/transverse/ileocecal valve;
# left = descending/sigmoid/rectosigmoid.  Generic "colon" is not enough.
# Medium-confidence entries are retained for sensitivity analysis only.
SIDE_CURATION: dict[str, dict[str, str]] = {
    "LS513": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Ileocecal valve",
        "source": "https://www.cellosaurus.org/CVCL_1386",
        "note": "Ileocecal valve classified as proximal/right-sided.",
    },
    "COLO-320": {
        "side": "Left",
        "confidence": "high",
        "primary_site": "Colon, sigmoid",
        "source": "https://www.cellosaurus.org/CVCL_1989",
        "note": "Sigmoid classified as distal/left-sided.",
    },
    "CL-11": {
        "side": "Left",
        "confidence": "high",
        "primary_site": "Left colon",
        "source": "https://www.dsmz.de/collection/catalogue/details/culture/ACC-467",
        "note": "DSMZ states primary colorectal cancer of the left colon.",
    },
    "LS1034": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Cecum",
        "source": "https://www.cellosaurus.org/CVCL_1382",
        "note": "Cecum classified as proximal/right-sided.",
    },
    "NCI-H747": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Cecum; sampled from common duct-node metastasis",
        "source": "https://www.atcc.org/products/ccl-252",
        "note": "Primary tumor is cecal; metastatic collection site does not change primary sidedness.",
    },
    "NCI-H716": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Cecum; cells harvested from ascites",
        "source": "https://pubmed.ncbi.nlm.nih.gov/1359704/",
        "note": "Original characterization describes a poorly differentiated adenocarcinoma of the caecum.",
    },
    "CL-40": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Right colon",
        "source": "https://www.dsmz.de/collection/catalogue/details/culture/ACC-535",
        "note": "DSMZ states primary colorectal cancer of the right colon.",
    },
    "SNU-C1": {
        "side": "Left",
        "confidence": "high",
        "primary_site": "Colon descendens; cultured from peritoneum",
        "source": "https://aacrjournals.org/clincancerres/article/5/3/643/199260/Thymidylate-Synthase-Level-as-the-Main-Predictive",
        "note": "Primary site is descending colon; peritoneal culture site is not used for sidedness.",
    },
    "SNU-C4": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Colon transversum",
        "source": "https://aacrjournals.org/clincancerres/article/5/3/643/199260/Thymidylate-Synthase-Level-as-the-Main-Predictive",
        "note": "Transverse colon classified as proximal/right-sided for this analysis.",
    },
    "SNU-C5": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Cecum",
        "source": "https://www.cellosaurus.org/CVCL_5112",
        "note": "Cecum classified as proximal/right-sided.",
    },
    "HCT 116": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Colon ascendens",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3816225/",
        "note": "Ascending colon classified as proximal/right-sided.",
    },
    "LS411N": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Cecum",
        "source": "https://www.atcc.org/products/crl-2159",
        "note": "Primary tumor biopsy from cecal carcinoma.",
    },
    "C10": {
        "side": "Left",
        "confidence": "high",
        "primary_site": "Colon, descending",
        "source": "https://www.cellosaurus.org/CVCL_5245",
        "note": "Descending colon classified as distal/left-sided.",
    },
    "C75": {
        "side": "Left",
        "confidence": "high",
        "primary_site": "Colon, sigmoid",
        "source": "https://www.cellosaurus.org/CVCL_5248",
        "note": "Sigmoid classified as distal/left-sided.",
    },
    "C84": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Cecum",
        "source": "https://www.cellosaurus.org/CVCL_5250",
        "note": "Cecum classified as proximal/right-sided.",
    },
    "SNU-1544": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Colon, ascending",
        "source": "https://www.cellosaurus.org/CVCL_5027",
        "note": "Ascending colon classified as proximal/right-sided.",
    },
    "SNU-1235": {
        "side": "Right",
        "confidence": "high",
        "primary_site": "Colon, ascending",
        "source": "https://www.cellosaurus.org/CVCL_5018",
        "note": "Ascending colon classified as proximal/right-sided.",
    },
    "HT-29": {
        "side": "Right",
        "confidence": "medium",
        "primary_site": "Colon, segment not explicit in Cellosaurus",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4104680/",
        "note": "Published sidedness assignment; retained only for sensitivity analysis.",
    },
    "DLD-1": {
        "side": "Left",
        "confidence": "medium",
        "primary_site": "Colon, segment not explicit in Cellosaurus",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4104680/",
        "note": "Published sidedness assignment; retained only for sensitivity analysis.",
    },
    "SW 620": {
        "side": "Left",
        "confidence": "medium",
        "primary_site": "Colon primary of paired SW480/SW620 patient; sampled from lymph node",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4104680/",
        "note": "Primary site is used rather than the metastatic collection site; sensitivity-only label.",
    },
    "HT-55": {
        "side": "Rectum",
        "confidence": "high",
        "primary_site": "Rectum",
        "source": "https://www.cellosaurus.org/CVCL_1294",
        "note": "Excluded from right-vs-left colon comparison.",
    },
    "SW 626": {
        "side": "Exclude",
        "confidence": "high",
        "primary_site": "Ovary metastasis; problematic/misclassified line",
        "source": "https://www.cellosaurus.org/CVCL_1725",
        "note": "Excluded because the record is flagged as problematic and records an ovarian metastatic site.",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-file", type=Path, default=DEFAULT_MODEL_FILE)
    parser.add_argument("--effect-file", type=Path, default=DEFAULT_EFFECT_FILE)
    parser.add_argument("--expression-file", type=Path, default=DEFAULT_EXPRESSION_FILE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Required local input is missing: {path}\n"
            "This script does not download data automatically."
        )


def normalize_name(value: object) -> str:
    return "".join(ch.lower() for ch in str(value) if ch.isalnum())


def curation_lookup() -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for name, annotation in SIDE_CURATION.items():
        lookup[normalize_name(name)] = annotation
    return lookup


def find_gene_column(path: Path, gene: str) -> tuple[int, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle))
    candidates = [
        (i, name)
        for i, name in enumerate(header)
        if name == gene or name.startswith(f"{gene} (")
    ]
    if not candidates:
        raise ValueError(f"Could not find gene {gene} in {path}")
    return candidates[0]


def load_gene_column(path: Path, gene: str, output_name: str) -> tuple[pd.DataFrame, str]:
    gene_index, column = find_gene_column(path, gene)
    frame = pd.read_csv(path, usecols=[0, gene_index], low_memory=False)
    frame = frame.rename(columns={frame.columns[0]: "ModelID", frame.columns[1]: output_name})
    frame["ModelID"] = frame["ModelID"].astype(str)
    frame[output_name] = pd.to_numeric(frame[output_name], errors="coerce")
    return frame, column


def zscore(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    sd = numeric.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index)
    return (numeric - numeric.mean()) / sd


def add_manual_sidedness(models: pd.DataFrame) -> pd.DataFrame:
    lookup = curation_lookup()
    annotations = []
    for _, row in models.iterrows():
        names = [row.get("CellLineName", ""), row.get("StrippedCellLineName", "")]
        annotation = next(
            (lookup[normalize_name(name)] for name in names if normalize_name(name) in lookup),
            None,
        )
        annotation = annotation or {}
        annotations.append(
            {
                "side": annotation.get("side", "Unknown"),
                "side_confidence": annotation.get("confidence", "none"),
                "curated_primary_site": annotation.get("primary_site", "Not explicitly assigned"),
                "side_evidence_source": annotation.get("source", ""),
                "curation_note": annotation.get(
                    "note",
                    "Generic colon/large intestine or metastatic collection site is insufficient for sidedness.",
                ),
            }
        )
    return pd.concat([models.reset_index(drop=True), pd.DataFrame(annotations)], axis=1)


def build_expression_state_score(
    expression_file: Path,
    coad_ids: pd.Series,
    genes: tuple[str, ...],
    score_name: str,
    complete_name: str,
    min_available: int,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Build a within-COAD z-score state with an explicit completeness gate."""
    score = pd.DataFrame({"ModelID": coad_ids.astype(str)})
    columns: dict[str, str] = {}
    for gene in genes:
        values, column = load_gene_column(expression_file, gene, f"{gene}_expression")
        columns[gene] = column
        score = score.merge(values, on="ModelID", how="left")
        score[f"{gene}_z"] = zscore(score[f"{gene}_expression"])
    z_columns = [f"{gene}_z" for gene in genes]
    raw_score = score[z_columns].mean(axis=1, skipna=True)
    score[complete_name] = score[z_columns].notna().sum(axis=1)
    score[f"{score_name}_raw"] = raw_score
    score[score_name] = raw_score.where(score[complete_name] >= min_available)
    return score, columns


def exact_directional_p(values: np.ndarray, groups: np.ndarray, seed: int) -> tuple[float, str]:
    """Test whether the Right mean is more negative than the Left mean."""
    right = values[groups == "Right"]
    left = values[groups == "Left"]
    if len(right) == 0 or len(left) == 0:
        return float("nan"), "not_estimable"
    observed = float(right.mean() - left.mean())
    n = len(values)
    n_right = len(right)
    total = math.comb(n, n_right)
    if total <= 250_000:
        count = 0
        for indices in itertools.combinations(range(n), n_right):
            mask = np.zeros(n, dtype=bool)
            mask[list(indices)] = True
            difference = float(values[mask].mean() - values[~mask].mean())
            count += difference <= observed + 1e-12
        return float((count + 1) / (total + 1)), f"exact_{total}_partitions"
    rng = np.random.default_rng(seed)
    count = 0
    n_permutations = 100_000
    for _ in range(n_permutations):
        mask = np.zeros(n, dtype=bool)
        mask[rng.choice(n, size=n_right, replace=False)] = True
        difference = float(values[mask].mean() - values[~mask].mean())
        count += difference <= observed
    return float((count + 1) / (n_permutations + 1)), f"monte_carlo_{n_permutations}"


def cliffs_delta(right: np.ndarray, left: np.ndarray) -> float:
    """Standard Cliff's delta: P(Right > Left) - P(Right < Left)."""
    if len(right) == 0 or len(left) == 0:
        return float("nan")
    pairwise = right[:, None] - left[None, :]
    return float((np.sum(pairwise > 0) - np.sum(pairwise < 0)) / pairwise.size)


def finite_or_none(value: object) -> float | int | None:
    if value is None:
        return None
    try:
        converted = float(value)
    except (TypeError, ValueError):
        return None
    return converted if np.isfinite(converted) else None


def format_p(value: object) -> str:
    """Compact plot annotation for a p-value or missing value."""
    if value is None:
        return "NA"
    try:
        converted = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{converted:.3g}" if np.isfinite(converted) else "NA"


def side_stats(panel: pd.DataFrame, confidence: set[str], seed: int) -> dict[str, object]:
    data = panel[
        panel["side"].isin(["Right", "Left"])
        & panel["side_confidence"].isin(confidence)
    ].dropna(subset=["SLC7A11_gene_effect"])
    right = data.loc[data["side"].eq("Right"), "SLC7A11_gene_effect"].to_numpy(float)
    left = data.loc[data["side"].eq("Left"), "SLC7A11_gene_effect"].to_numpy(float)
    permutation_p, permutation_method = exact_directional_p(
        data["SLC7A11_gene_effect"].to_numpy(float), data["side"].to_numpy(), seed
    )
    mw_p = (
        float(mannwhitneyu(right, left, alternative="less", method="auto").pvalue)
        if len(right) and len(left)
        else float("nan")
    )
    return {
        "n_total": int(len(data)),
        "n_right": int(len(right)),
        "n_left": int(len(left)),
        "right_mean_gene_effect": finite_or_none(right.mean() if len(right) else np.nan),
        "left_mean_gene_effect": finite_or_none(left.mean() if len(left) else np.nan),
        "right_median_gene_effect": finite_or_none(np.median(right) if len(right) else np.nan),
        "left_median_gene_effect": finite_or_none(np.median(left) if len(left) else np.nan),
        "right_minus_left_gene_effect": finite_or_none(
            (right.mean() - left.mean()) if len(right) and len(left) else np.nan
        ),
        "right_more_negative_exact_p": finite_or_none(permutation_p),
        "exact_test_method": permutation_method,
        "right_more_negative_mann_whitney_p": finite_or_none(mw_p),
        "cliffs_delta_right_vs_left": finite_or_none(cliffs_delta(right, left)),
        "right_fraction_effect_le_minus_0_5": finite_or_none(
            np.mean(right <= -0.5) if len(right) else np.nan
        ),
        "left_fraction_effect_le_minus_0_5": finite_or_none(
            np.mean(left <= -0.5) if len(left) else np.nan
        ),
    }


def fit_state_interaction(
    panel: pd.DataFrame,
    confidence: set[str],
    state_column: str,
    state_name: str,
) -> dict[str, object]:
    """Fit KO effect ~ side + state + side:state with HC3 contrasts."""
    data = panel[
        panel["side"].isin(["Right", "Left"])
        & panel["side_confidence"].isin(confidence)
    ].dropna(subset=[state_column, "SLC7A11_gene_effect"]).copy()
    if len(data) < 8 or data["side"].nunique() < 2:
        return {"state": state_name, "n": int(len(data)), "estimable": False}
    data["right_indicator"] = data["side"].eq("Right").astype(float)
    data["state_centered"] = data[state_column] - data[state_column].mean()
    data["side_x_state"] = data["right_indicator"] * data["state_centered"]
    design = sm.add_constant(data[["right_indicator", "state_centered", "side_x_state"]])
    try:
        fit = sm.OLS(data["SLC7A11_gene_effect"], design).fit(cov_type="HC3")
    except (np.linalg.LinAlgError, ValueError):
        return {"state": state_name, "n": int(len(data)), "estimable": False}
    if fit.df_resid <= 0:
        return {"state": state_name, "n": int(len(data)), "estimable": False}

    params = fit.params.to_numpy(float)
    covariance = np.asarray(fit.cov_params(), dtype=float)
    right_contrast = np.array([0.0, 0.0, 1.0, 1.0])
    right_variance = float(right_contrast @ covariance @ right_contrast)
    right_se = math.sqrt(right_variance) if right_variance > 0 else float("nan")
    right_slope = float(right_contrast @ params)
    right_t = right_slope / right_se if np.isfinite(right_se) and right_se > 0 else float("nan")
    right_p = float(2 * t.sf(abs(right_t), fit.df_resid)) if np.isfinite(right_t) else float("nan")
    return {
        "state": state_name,
        "n": int(len(data)),
        "n_right": int(data["right_indicator"].sum()),
        "n_left": int((1 - data["right_indicator"]).sum()),
        "estimable": True,
        "side_main_effect_right_minus_left": finite_or_none(fit.params.get("right_indicator")),
        "side_main_effect_p_hc3": finite_or_none(fit.pvalues.get("right_indicator")),
        "left_specific_slope": finite_or_none(fit.params.get("state_centered")),
        "left_specific_slope_p_linear_contrast": finite_or_none(fit.pvalues.get("state_centered")),
        "right_specific_slope": finite_or_none(right_slope),
        "right_specific_slope_p_linear_contrast": finite_or_none(right_p),
        "right_minus_left_slope_beta_interaction": finite_or_none(fit.params.get("side_x_state")),
        "right_minus_left_slope_p_interaction": finite_or_none(fit.pvalues.get("side_x_state")),
        "r_squared": finite_or_none(fit.rsquared),
    }


def state_cells(
    panel: pd.DataFrame,
    confidence: set[str],
    state_column: str,
    state_name: str,
) -> dict[str, object]:
    data = panel[
        panel["side"].isin(["Right", "Left"])
        & panel["side_confidence"].isin(confidence)
    ].dropna(subset=[state_column, "SLC7A11_gene_effect"]).copy()
    available = panel[state_column].dropna()
    if available.empty:
        return {"state": state_name, "global_cutoff": None, "n_labeled_models": 0, "cells": {}}
    cutoff = float(available.median())
    data["state_group"] = np.where(data[state_column] >= cutoff, "High", "Low")
    cells: dict[str, object] = {}
    for side in ["Right", "Left"]:
        for group in ["High", "Low"]:
            values = data.loc[
                data["side"].eq(side) & data["state_group"].eq(group),
                "SLC7A11_gene_effect",
            ]
            cells[f"{side}_{group}"] = {
                "n": int(len(values)),
                "mean_gene_effect": finite_or_none(values.mean()),
                "median_gene_effect": finite_or_none(values.median()),
            }
    return {"state": state_name, "global_cutoff": cutoff, "n_labeled_models": int(len(data)), "cells": cells}


def make_figure(panel: pd.DataFrame, summary: dict[str, object], out_dir: Path) -> None:
    primary = panel[
        panel["side"].isin(["Right", "Left"])
        & panel["side_confidence"].eq("high")
    ].dropna(subset=["SLC7A11_gene_effect"])
    colors = {"Right": "#d95f5f", "Left": "#4c78a8"}
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6), dpi=220)

    for index, side in enumerate(["Right", "Left"]):
        values = primary.loc[primary["side"].eq(side), "SLC7A11_gene_effect"].to_numpy(float)
        axes[0].scatter(
            np.full(len(values), index + 1),
            values,
            color=colors[side],
            edgecolor="white",
            linewidth=0.5,
            s=62,
            zorder=3,
            label=f"{side} (n={len(values)})",
        )
        if len(values):
            axes[0].hlines(np.mean(values), index + 0.78, index + 1.22, color=colors[side], linewidth=3)
    axes[0].set_xticks([1, 2], ["Right", "Left"])
    axes[0].set_ylabel("SLC7A11 Chronos gene effect\nmore negative = stronger dependency")
    primary_p = summary["primary_side_comparison"].get("right_more_negative_exact_p")
    axes[0].set_title(
        f"Manual sidedness (high confidence)\nexact one-sided p={format_p(primary_p)}",
        fontsize=10,
    )
    axes[0].axhline(-0.5, color="gray", linestyle="--", linewidth=0.8)
    axes[0].legend(frameon=False, fontsize=8)

    def plot_state(axis: plt.Axes, state_column: str, x_label: str, summary_key: str) -> None:
        data = panel[
            panel["side"].isin(["Right", "Left"])
            & panel["side_confidence"].isin(["high", "medium"])
        ].dropna(subset=[state_column, "SLC7A11_gene_effect"])
        for side in ["Right", "Left"]:
            side_data = data[data["side"].eq(side)]
            axis.scatter(
                side_data[state_column],
                side_data["SLC7A11_gene_effect"],
                color=colors[side],
                edgecolor="white",
                linewidth=0.5,
                s=58,
                label=side,
            )
            if len(side_data) >= 3:
                beta = np.polyfit(side_data[state_column], side_data["SLC7A11_gene_effect"], 1)
                xx = np.linspace(side_data[state_column].min(), side_data[state_column].max(), 50)
                axis.plot(xx, beta[0] * xx + beta[1], color=colors[side], linewidth=1.8)
        interaction_p = summary[summary_key].get("right_minus_left_slope_p_interaction")
        right_p = summary[summary_key].get("right_specific_slope_p_linear_contrast")
        axis.set_xlabel(x_label)
        axis.set_ylabel("SLC7A11 Chronos gene effect")
        axis.set_title(
            f"{x_label.split(chr(10))[0]} × sidedness\n"
            f"Δslope p={format_p(interaction_p)}; Right-slope p={format_p(right_p)}",
            fontsize=10,
        )
        axis.axhline(-0.5, color="gray", linestyle="--", linewidth=0.8)
        axis.legend(frameon=False)

    plot_state(
        axes[1],
        "PUFA_pressure_score",
        "PUFA-pressure proxy\n≥4/5 genes available",
        "primary_pressure_interaction",
    )
    plot_state(
        axes[2],
        "AA_routing_score",
        "AA-routing proxy\n≥3/4 genes available",
        "primary_aa_routing_interaction",
    )
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(out_dir / "depmap_slc7a11_virtual_ko_redo.png", bbox_inches="tight")
    plt.close(fig)


def write_report(summary: dict[str, object], out_dir: Path) -> None:
    primary = summary["primary_side_comparison"]
    sensitivity = summary["sensitivity_side_comparison"]
    pressure_interaction = summary["primary_pressure_interaction"]
    routing_interaction = summary["primary_aa_routing_interaction"]
    lines = [
        "# SLC7A11 virtual knockout redo",
        "",
        "## Question",
        "",
        "Does SLC7A11 CRISPR/Chronos gene effect differ between manually curated right- and left-sided COAD models? More-negative gene effect means a larger fitness loss after SLC7A11 knockout.",
        "",
        "## Design",
        "",
        "- Primary unit: one DepMap COAD cell model.",
        "- Primary sidedness analysis: explicit high-confidence anatomical primary-site labels only.",
        "- Sensitivity analysis: high-confidence plus three medium-confidence published assignments.",
        "- Rectal and problematic/misclassified ovarian-site records are not included in the Right-versus-Left comparison.",
        "- Secondary models test `SLC7A11_gene_effect ~ side + state + side:state` for both PUFA-pressure and AA-routing states; all state-by-side terms are exploratory and do not establish causality.",
        "",
        "## Primary result",
        "",
        f"- Right: n={primary['n_right']}; mean gene effect={primary['right_mean_gene_effect']}",
        f"- Left: n={primary['n_left']}; mean gene effect={primary['left_mean_gene_effect']}",
        f"- Right-minus-Left mean difference={primary['right_minus_left_gene_effect']}",
        f"- Exact one-sided permutation p (Right more negative)={primary['right_more_negative_exact_p']} ({primary['exact_test_method']})",
        f"- One-sided Mann–Whitney p={primary['right_more_negative_mann_whitney_p']}",
        f"- Cliff's delta (Right versus Left)={primary['cliffs_delta_right_vs_left']}",
        "",
        "## Sensitivity result",
        "",
        f"- Right: n={sensitivity['n_right']}; Left: n={sensitivity['n_left']}",
        f"- Right-minus-Left mean difference={sensitivity['right_minus_left_gene_effect']}",
        f"- Exact one-sided permutation p={sensitivity['right_more_negative_exact_p']} ({sensitivity['exact_test_method']})",
        "",
        "## Pressure interaction",
        "",
        f"- Matched labeled models={pressure_interaction.get('n')}; Right-specific slope={pressure_interaction.get('right_specific_slope')}; linear-contrast P={pressure_interaction.get('right_specific_slope_p_linear_contrast')}",
        f"- Right-minus-Left slope interaction beta={pressure_interaction.get('right_minus_left_slope_beta_interaction')}; interaction P={pressure_interaction.get('right_minus_left_slope_p_interaction')}",
        "- The PUFA-pressure proxy is the mean within-COAD z-score of the supplied DepMap expression values for ACSL4, LPCAT3, ALOX5, ALOX12 and ALOX15, retained only when at least 4/5 genes are available. It is not a measurement of AA concentration or flux.",
        "",
        "## AA-routing interaction",
        "",
        f"- Matched labeled models={routing_interaction.get('n')}; Right-specific slope={routing_interaction.get('right_specific_slope')}; linear-contrast P={routing_interaction.get('right_specific_slope_p_linear_contrast')}",
        f"- Right-minus-Left slope interaction beta={routing_interaction.get('right_minus_left_slope_beta_interaction')}; interaction P={routing_interaction.get('right_minus_left_slope_p_interaction')}",
        "- The AA-routing proxy is the mean within-COAD z-score of PLA2G4A, PTGS2, ALOX5 and ALOX5AP, retained when at least 3/4 genes are available. It is a transcriptional proxy, not a direct eicosanoid or AA-flux measurement.",
        "",
        "## Interpretation guardrails",
        "",
        "A negative Right-minus-Left effect estimate is compatible with stronger SLC7A11 knockout sensitivity in right-sided models. A null result does not prove equal biology because the curated cell-line sample is small and sidedness labels are incomplete. A positive result or a pressure interaction remains an association, not proof that AA caused the dependency.",
        "",
        "## Provenance",
        "",
        json.dumps(summary["provenance"], ensure_ascii=False, indent=2),
    ]
    (out_dir / "depmap_slc7a11_virtual_ko_redo_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    for path in [args.model_file, args.effect_file, args.expression_file]:
        require_file(path)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    models = pd.read_csv(args.model_file, low_memory=False)
    required_model_columns = {"ModelID", "CellLineName", "StrippedCellLineName", "OncotreeCode"}
    missing = required_model_columns.difference(models.columns)
    if missing:
        raise ValueError(f"Model file is missing required columns: {sorted(missing)}")
    coad = models[models["OncotreeCode"].eq("COAD")].copy()
    coad = coad[
        [
            "ModelID",
            "CellLineName",
            "StrippedCellLineName",
            "OncotreeCode",
            "Sex",
            "PrimaryOrMetastasis",
            "SampleCollectionSite",
        ]
    ]
    coad["ModelID"] = coad["ModelID"].astype(str)
    coad = add_manual_sidedness(coad)

    effects, effect_column = load_gene_column(args.effect_file, GENE, "SLC7A11_gene_effect")
    panel = coad.merge(effects, on="ModelID", how="inner")
    pressure, pressure_columns = build_expression_state_score(
        args.expression_file,
        panel["ModelID"],
        PRESSURE_GENES,
        "PUFA_pressure_score",
        "pressure_complete_genes",
        PUFA_MIN_AVAILABLE,
    )
    panel = panel.merge(pressure, on="ModelID", how="left")
    routing, routing_columns = build_expression_state_score(
        args.expression_file,
        panel["ModelID"],
        AA_ROUTING_GENES,
        "AA_routing_score",
        "aa_routing_complete_genes",
        AA_ROUTING_MIN_AVAILABLE,
    )
    panel = panel.merge(routing, on="ModelID", how="left")
    panel["pressure_group_global_median"] = np.where(
        panel["PUFA_pressure_score"].notna(),
        np.where(
            panel["PUFA_pressure_score"] >= panel["PUFA_pressure_score"].median(),
            "High",
            "Low",
        ),
        "Unavailable",
    )

    primary = side_stats(panel, {"high"}, args.seed)
    sensitivity = side_stats(panel, {"high", "medium"}, args.seed)
    primary_interaction = fit_state_interaction(
        panel, {"high"}, "PUFA_pressure_score", "PUFA-pressure"
    )
    sensitivity_interaction = fit_state_interaction(
        panel, {"high", "medium"}, "PUFA_pressure_score", "PUFA-pressure"
    )
    primary_aa_routing_interaction = fit_state_interaction(
        panel, {"high"}, "AA_routing_score", "AA-routing"
    )
    sensitivity_aa_routing_interaction = fit_state_interaction(
        panel, {"high", "medium"}, "AA_routing_score", "AA-routing"
    )
    summary: dict[str, object] = {
        "target_gene": GENE,
        "analysis_set": "DepMap OncotreeCode=COAD",
        "primary_side_comparison": primary,
        "sensitivity_side_comparison": sensitivity,
        "primary_pressure_interaction": primary_interaction,
        "sensitivity_pressure_interaction": sensitivity_interaction,
        "primary_aa_routing_interaction": primary_aa_routing_interaction,
        "sensitivity_aa_routing_interaction": sensitivity_aa_routing_interaction,
        "pressure_cells_primary": state_cells(panel, {"high"}, "PUFA_pressure_score", "PUFA-pressure"),
        "pressure_cells_sensitivity": state_cells(
            panel, {"high", "medium"}, "PUFA_pressure_score", "PUFA-pressure"
        ),
        "aa_routing_cells_primary": state_cells(
            panel, {"high"}, "AA_routing_score", "AA-routing"
        ),
        "aa_routing_cells_sensitivity": state_cells(
            panel, {"high", "medium"}, "AA_routing_score", "AA-routing"
        ),
        "n_coad_models": int(len(coad)),
        "n_matched_dependency_models": int(len(panel)),
        "side_counts": panel["side"].value_counts(dropna=False).to_dict(),
        "confidence_counts": panel["side_confidence"].value_counts(dropna=False).to_dict(),
        "provenance": {
            "model_file": str(args.model_file),
            "dependency_file": str(args.effect_file),
            "expression_file": str(args.expression_file),
            "dependency_column": effect_column,
            "pressure_expression_columns": pressure_columns,
            "aa_routing_expression_columns": routing_columns,
            "pufa_min_available_genes": PUFA_MIN_AVAILABLE,
            "aa_routing_min_available_genes": AA_ROUTING_MIN_AVAILABLE,
            "dependency_release": "DepMap Public 26Q1 CRISPRGeneEffect / Chronos gene effect",
            "expression_release": "DepMap 24Q4 expression",
            "manual_curation_entries": len(SIDE_CURATION),
            "seed": args.seed,
            "no_external_download": True,
        },
    }

    panel.to_csv(args.out_dir / "depmap_slc7a11_virtual_ko_redo_panel.csv", index=False)
    panel[
        [
            "ModelID",
            "CellLineName",
            "side",
            "side_confidence",
            "curated_primary_site",
            "SLC7A11_gene_effect",
            "PUFA_pressure_score",
            "pressure_complete_genes",
            "AA_routing_score",
            "aa_routing_complete_genes",
            "pressure_group_global_median",
        ]
    ].to_csv(args.out_dir / "depmap_slc7a11_virtual_ko_redo_sidedness.csv", index=False)
    (args.out_dir / "depmap_slc7a11_virtual_ko_redo_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (args.out_dir / "depmap_slc7a11_virtual_ko_redo_manifest.json").write_text(
        json.dumps(summary["provenance"], ensure_ascii=False, indent=2), encoding="utf-8"
    )
    make_figure(panel, summary, args.out_dir)
    write_report(summary, args.out_dir)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
