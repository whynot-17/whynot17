from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
DATA = WORK / "data"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

MODEL = WORK / "depmap_Model.csv"
EXPRESSION = WORK / "depmap_expression_24Q4.csv"
GENE_EFFECT = DATA / "CRISPRGeneEffect_26Q1.csv"
GENE = "DHODH"


def find_column(path: Path, prefix: str) -> tuple[int, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle))
    candidates = [(i, name) for i, name in enumerate(header) if name == prefix or name.startswith(f"{prefix} (")]
    if not candidates:
        raise ValueError(f"Could not find {prefix} in {path}")
    return candidates[0]


def load_gene_column(path: Path, gene: str, output_name: str) -> pd.DataFrame:
    gene_index, column = find_column(path, gene)
    frame = pd.read_csv(path, usecols=[0, gene_index])
    frame = frame.rename(columns={frame.columns[0]: "ModelID", frame.columns[1]: output_name})
    frame["ModelID"] = frame["ModelID"].astype(str)
    frame[output_name] = pd.to_numeric(frame[output_name], errors="coerce")
    return frame


def fit_association(panel: pd.DataFrame) -> dict[str, object]:
    d = panel.dropna(subset=["DHODH_expression", "DHODH_gene_effect"]).copy()
    x = d["DHODH_expression"].to_numpy(float)
    y = d["DHODH_gene_effect"].to_numpy(float)
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_rho, spearman_p = spearmanr(x, y)
    X = sm.add_constant(d[["DHODH_expression"]])
    fit = sm.OLS(d["DHODH_gene_effect"], X).fit(cov_type="HC3")
    return {
        "n": int(len(d)),
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "spearman_rho": float(spearman_rho),
        "spearman_p": float(spearman_p),
        "ols_beta_expression": float(fit.params["DHODH_expression"]),
        "ols_se_hc3": float(fit.bse["DHODH_expression"]),
        "ols_p_hc3": float(fit.pvalues["DHODH_expression"]),
        "ols_r2": float(fit.rsquared),
        "effect_mean": float(y.mean()),
        "effect_median": float(np.median(y)),
        "n_effect_le_minus_0_5": int((y <= -0.5).sum()),
        "n_effect_le_minus_1": int((y <= -1.0).sum()),
    }


def make_figure(panel: pd.DataFrame) -> None:
    d = panel.dropna(subset=["DHODH_expression", "DHODH_gene_effect"]).copy()
    colors = {"Male": "#3b6fb6", "Female": "#c45a72"}
    fig, ax = plt.subplots(figsize=(6.6, 5.0), dpi=220)
    for sex in ["Male", "Female"]:
        z = d[d.Sex.eq(sex)]
        ax.scatter(z.DHODH_expression, z.DHODH_gene_effect, s=48, alpha=0.85, color=colors[sex], label=f"{sex} (n={len(z)})", edgecolor="white", linewidth=0.4)
    if len(d) >= 3:
        b1, b0 = np.polyfit(d.DHODH_expression.to_numpy(), d.DHODH_gene_effect.to_numpy(), 1)
        xx = np.linspace(d.DHODH_expression.min(), d.DHODH_expression.max(), 100)
        ax.plot(xx, b0 + b1 * xx, color="#333333", linewidth=1.6, label="OLS trend")
    ax.axhline(-0.5, color="#888888", linestyle="--", linewidth=0.8)
    ax.set_xlabel("DHODH expression (DepMap 24Q4 log2 TPM+1)")
    ax.set_ylabel("DHODH CRISPR gene effect (26Q1)\nmore negative = stronger fitness loss")
    ax.set_title("CRC: DHODH expression versus DHODH dependency")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "depmap_dhodh_expression_dependency.png", bbox_inches="tight")
    plt.close(fig)


def write_report(panel: pd.DataFrame, summary: dict[str, object], metadata: dict[str, object]) -> None:
    sex_counts = panel.Sex.value_counts().to_dict()
    lines = [
        "# DHODH expression versus dependency in CRC",
        "",
        "## Question",
        "",
        "Do CRC models with higher DHODH expression show a more negative DHODH CRISPR/Chronos gene effect? This is the first, constitutive dependency check and does not claim an acquired OXA-R dependency.",
        "",
        "## Design",
        "",
        "- CRC restriction: DepMap `OncotreePrimaryDisease = Colorectal Adenocarcinoma` with known biological sex.",
        "- Expression: DepMap 24Q4 expression matrix, DHODH column.",
        "- Dependency: DepMap 26Q1 `CRISPRGeneEffect`/Chronos gene effect, DHODH column.",
        "- Statistical unit: one DepMap model; more-negative gene effect means greater fitness loss after knockout.",
        "- Prespecified association: Spearman/Pearson correlation and HC3-robust linear regression; no genome-wide multiple-testing correction because only one gene pair was tested.",
        "",
        "## Matched panel",
        "",
        f"n={summary['n']} CRC models; sex counts={sex_counts}.",
        "",
        "## Result",
        "",
        f"Spearman rho={summary['spearman_rho']:.4g}, P={summary['spearman_p']:.4g}; Pearson r={summary['pearson_r']:.4g}, P={summary['pearson_p']:.4g}.",
        f"HC3-robust OLS beta for expression={summary['ols_beta_expression']:.4g}, SE={summary['ols_se_hc3']:.4g}, P={summary['ols_p_hc3']:.4g}, R²={summary['ols_r2']:.4g}.",
        f"Mean gene effect={summary['effect_mean']:.4g}; models with effect ≤−0.5: {summary['n_effect_le_minus_0_5']}; ≤−1: {summary['n_effect_le_minus_1']}.",
        "",
        "## Interpretation",
        "",
        "A negative expression–effect association would be compatible with higher DHODH expression marking greater constitutive DHODH dependency. A null result would mean that DHODH upregulation in OXA-R cannot be justified as a general CRC DHODH dependency from this analysis alone.",
        "",
        "Even a positive association is not evidence that OXA resistance caused the dependency. The next step, only if warranted, is to project the OXA-R DHODH-centered transition state into the CRC panel and test whether that context strengthens the association.",
        "",
        "## Provenance",
        "",
        json.dumps(metadata, ensure_ascii=False, indent=2),
    ]
    (OUT / "depmap_dhodh_expression_dependency_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    model = pd.read_csv(MODEL, low_memory=False)
    model = model[
        model["OncotreePrimaryDisease"].eq("Colorectal Adenocarcinoma")
        & model["Sex"].isin(["Male", "Female"])
    ][["ModelID", "CellLineName", "Sex", "OncotreeCode"]].copy()
    model["ModelID"] = model["ModelID"].astype(str)
    expression = load_gene_column(EXPRESSION, GENE, "DHODH_expression")
    effect = load_gene_column(GENE_EFFECT, GENE, "DHODH_gene_effect")
    panel = model.merge(expression, on="ModelID", how="inner").merge(effect, on="ModelID", how="inner")
    summary = fit_association(panel)
    metadata = {
        "gene": GENE,
        "model_file": str(MODEL),
        "expression_file": str(EXPRESSION),
        "dependency_file": str(GENE_EFFECT),
        "expression_column": find_column(EXPRESSION, GENE)[1],
        "dependency_column": find_column(GENE_EFFECT, GENE)[1],
        "disease_filter": "OncotreePrimaryDisease == Colorectal Adenocarcinoma",
        "sex_filter": ["Male", "Female"],
        "n_model_metadata": int(len(model)),
        "n_expression": int(len(expression)),
        "n_dependency": int(len(effect)),
        "n_matched": int(len(panel)),
        "sex_counts": panel.Sex.value_counts().to_dict(),
        "dependency_release": "DepMap 26Q1 CRISPRGeneEffect / Chronos gene effect",
        "expression_release": "DepMap 24Q4 expression",
    }
    panel.to_csv(OUT / "depmap_dhodh_expression_dependency_panel.csv", index=False)
    (OUT / "depmap_dhodh_expression_dependency_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "depmap_dhodh_expression_dependency_manifest.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    make_figure(panel)
    write_report(panel, summary, metadata)
    print(json.dumps({"summary": summary, "metadata": metadata}, indent=2))
    print(panel[["ModelID", "CellLineName", "Sex", "DHODH_expression", "DHODH_gene_effect"]].sort_values("DHODH_gene_effect").to_string(index=False))


if __name__ == "__main__":
    main()
