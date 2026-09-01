"""Three-layer test of ferroptosis sidedness in colorectal cancer.

Goal
----
Answer the bounded first question before mechanism work:
    Does right-sided CRC (RCRC) show a stronger ferroptotic state/liability than
    left-sided CRC (LCRC)?

Layers
------
1. Transcriptomic ferroptosis state in GSE39582 and TCGA-COAD:
   driver score, defense score, and net propensity = driver - defense.
2. Lipid-peroxidation liability in the same cohorts:
   PUFA incorporation, peroxide generation, antioxidant buffering, and a net
   lipid-peroxidation liability score.
3. Functional dependency in manually curated DepMap COAD models:
   Chronos effects for canonical ferroptosis-defense genes and a composite
   defense-dependency index.

Important guardrail
-------------------
These scores are not direct measurements of ferroptotic cell death. They are
used to establish direction and consistency across independent evidence layers.
A direct sidedness claim still requires functional assays such as lipid ROS and
ferroptosis-rescue experiments.
"""

from __future__ import annotations

import io
import json
import math
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, ttest_ind

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
DATA = WORK / "data"
OUT = ROOT / "outputs"

sys.path.insert(0, str(WORK))
from phase1_bulk_analysis import parse_gpl_probe_map, read_gse_targeted, tcga_case_metadata  # noqa: E402
from depmap_slc7a11_virtual_ko_redo import (  # noqa: E402
    add_manual_sidedness,
    load_gene_column,
)

GSE_MATRIX = DATA / "GSE39582_series_matrix.clean.txt.gz"
GPL570_ANNOT = DATA / "GPL570.annot.gz"
TCGA_CASES = DATA / "tcga_coad_cases.json"
TCGA_SEED_EXPRESSION = DATA / "tcga_coad_gene_expression_uqfpkm.tsv"
TCGA_TARGET_EXPRESSION = DATA / "tcga_ferroptosis_three_layer_expression_uqfpkm.tsv"
DEPMAP_MODEL = WORK / "depmap_Model.csv"
DEPMAP_EFFECT = DATA / "CRISPRGeneEffect_26Q1.csv"

# Direction-aware, deliberately compact canonical sets. Mixed-role KEGG/FerrDb
# totals are not used as a single score because drivers and suppressors cancel.
FERROPTOSIS_DRIVER = [
    "ACSL4", "LPCAT3", "TFRC", "NCOA4", "SAT1", "POR", "CYB5R1", "ALOX12", "ALOX15"
]
FERROPTOSIS_DEFENSE = [
    "SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH", "NFE2L2", "GCLC", "GCLM", "FTH1"
]
PUFA_INCORPORATION = ["ACSL4", "LPCAT3", "AGPAT3"]
PEROXIDE_GENERATION = ["POR", "CYB5R1", "ALOX12", "ALOX15"]
ANTIOXIDANT_BUFFERING = ["SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH"]
DEPMAP_DEFENSE_GENES = ["SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH"]

ALL_TRANSCRIPT_GENES = list(dict.fromkeys(
    FERROPTOSIS_DRIVER
    + FERROPTOSIS_DEFENSE
    + PUFA_INCORPORATION
    + PEROXIDE_GENERATION
    + ANTIOXIDANT_BUFFERING
))

# Stable Ensembl gene IDs used by GDC gene_expression/values.
TCGA_GENE_IDS = {
    "ACSL4": "ENSG00000068366",
    "LPCAT3": "ENSG00000111684",
    "TFRC": "ENSG00000072274",
    "NCOA4": "ENSG00000138293",
    "SAT1": "ENSG00000130066",
    "POR": "ENSG00000127948",
    "CYB5R1": "ENSG00000159348",
    "ALOX12": "ENSG00000108839",
    "ALOX15": "ENSG00000161905",
    "SLC7A11": "ENSG00000151012",
    "GPX4": "ENSG00000167468",
    "AIFM2": "ENSG00000042286",
    "GCH1": "ENSG00000131979",
    "DHODH": "ENSG00000102967",
    "NFE2L2": "ENSG00000116044",
    "GCLC": "ENSG00000001084",
    "GCLM": "ENSG00000023909",
    "FTH1": "ENSG00000167996",
    "AGPAT3": "ENSG00000160216",
}


def require(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required local input: {path}")


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / sd


def bh_adjust(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out.tolist()
    idx = np.where(ok)[0]
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty_like(adj)
    tmp[order] = np.minimum(adj, 1.0)
    out[idx] = tmp
    return out.tolist()


def score_mean_z(expr: pd.DataFrame, genes: list[str], min_fraction: float = 0.70) -> tuple[pd.Series, pd.Series]:
    available = [g for g in genes if g in expr.columns]
    if not available:
        return pd.Series(np.nan, index=expr.index), pd.Series(0, index=expr.index)
    z = pd.DataFrame({g: zscore(expr[g]) for g in available}, index=expr.index)
    n = z.notna().sum(axis=1)
    minimum = max(1, math.ceil(len(genes) * min_fraction))
    score = z.mean(axis=1, skipna=True).where(n >= minimum)
    return score, n


def load_gse() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    require(GSE_MATRIX)
    require(GPL570_ANNOT)
    probe_map, probe_counts = parse_gpl_probe_map(GPL570_ANNOT, ALL_TRANSCRIPT_GENES)
    expr, meta, info = read_gse_targeted(GSE_MATRIX, probe_map, ALL_TRANSCRIPT_GENES)
    meta = meta.copy()
    meta["side"] = meta["side"].astype(str).str.lower()
    keep = meta["side"].isin(["left", "right"]) & ~meta["is_normal_like"].fillna(False)
    expr = expr.loc[keep].apply(pd.to_numeric, errors="coerce")
    meta = meta.loc[keep].copy()
    return expr, meta, {**info, "probe_counts": probe_counts, "n": int(len(expr))}


def tcga_case_ids() -> list[str]:
    require(TCGA_SEED_EXPRESSION)
    header = pd.read_csv(TCGA_SEED_EXPRESSION, sep="\t", nrows=0)
    return [str(x) for x in header.columns[1:]]


def load_tcga() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    require(TCGA_CASES)
    case_ids = tcga_case_ids()
    missing_ids = sorted(set(ALL_TRANSCRIPT_GENES) - set(TCGA_GENE_IDS))
    if missing_ids:
        raise ValueError(f"Missing Ensembl mappings: {missing_ids}")

    if TCGA_TARGET_EXPRESSION.exists() and TCGA_TARGET_EXPRESSION.stat().st_size > 1000:
        text = TCGA_TARGET_EXPRESSION.read_text(encoding="utf-8")
        source = "cached targeted GDC expression"
    else:
        payload = {
            "case_ids": case_ids,
            "gene_ids": [TCGA_GENE_IDS[g] for g in ALL_TRANSCRIPT_GENES],
            "tsv_units": "uqfpkm",
            "format": "tsv",
        }
        request = urllib.request.Request(
            "https://api.gdc.cancer.gov/gene_expression/values",
            data=json.dumps(payload).encode("utf-8"),
            headers={"Accept": "text/tab-separated-values", "Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(request, timeout=180) as response:
            text = response.read().decode("utf-8")
        TCGA_TARGET_EXPRESSION.write_text(text, encoding="utf-8")
        source = "GDC gene_expression/values endpoint"

    matrix = pd.read_csv(io.StringIO(text), sep="\t")
    gene_col = matrix.columns[0]
    id_to_symbol = {v: k for k, v in TCGA_GENE_IDS.items()}
    matrix["symbol"] = matrix[gene_col].astype(str).map(id_to_symbol)
    matrix = matrix[matrix["symbol"].notna()].set_index("symbol")
    expr = matrix.drop(columns=[gene_col], errors="ignore").T.apply(pd.to_numeric, errors="coerce")
    expr = np.log2(expr.clip(lower=0) + 1.0)
    expr.index = expr.index.astype(str)

    meta = tcga_case_metadata(TCGA_CASES)
    common = expr.index.intersection(meta.index)
    expr = expr.loc[common]
    meta = meta.loc[common].copy()
    meta["side"] = meta["side"].astype(str).str.lower()
    keep = meta["side"].isin(["left", "right"])
    expr = expr.loc[keep]
    meta = meta.loc[keep]
    return expr, meta, {"source": source, "n": int(len(expr)), "genes": list(expr.columns)}


def build_scores(expr: pd.DataFrame, meta: pd.DataFrame, cohort: str) -> pd.DataFrame:
    out = meta.copy()
    out.index = expr.index
    driver, driver_n = score_mean_z(expr, FERROPTOSIS_DRIVER)
    defense, defense_n = score_mean_z(expr, FERROPTOSIS_DEFENSE)
    pufa, pufa_n = score_mean_z(expr, PUFA_INCORPORATION)
    peroxide, peroxide_n = score_mean_z(expr, PEROXIDE_GENERATION)
    buffer, buffer_n = score_mean_z(expr, ANTIOXIDANT_BUFFERING)

    out["ferroptosis_driver_score"] = driver
    out["ferroptosis_defense_score"] = defense
    out["ferroptosis_net_propensity"] = driver - defense
    out["pufa_incorporation_score"] = pufa
    out["peroxide_generation_score"] = peroxide
    out["antioxidant_buffering_score"] = buffer
    out["lipid_peroxidation_liability"] = ((pufa + peroxide) / 2.0) - buffer
    out["driver_genes_available"] = driver_n
    out["defense_genes_available"] = defense_n
    out["pufa_genes_available"] = pufa_n
    out["peroxide_genes_available"] = peroxide_n
    out["buffer_genes_available"] = buffer_n
    out["cohort"] = cohort
    out["sample_id"] = out.index.astype(str)
    return out


def standardized_mean_diff(right: np.ndarray, left: np.ndarray) -> float:
    if len(right) < 2 or len(left) < 2:
        return float("nan")
    vr = np.var(right, ddof=1)
    vl = np.var(left, ddof=1)
    pooled = math.sqrt(((len(right) - 1) * vr + (len(left) - 1) * vl) / (len(right) + len(left) - 2))
    if pooled == 0 or not np.isfinite(pooled):
        return float("nan")
    return float((np.mean(right) - np.mean(left)) / pooled)


def sidedness_stats(scores: pd.DataFrame, columns: list[str], layer: str) -> pd.DataFrame:
    rows = []
    for col in columns:
        r = pd.to_numeric(scores.loc[scores["side"].eq("right"), col], errors="coerce").dropna().to_numpy()
        l = pd.to_numeric(scores.loc[scores["side"].eq("left"), col], errors="coerce").dropna().to_numpy()
        welch = float(ttest_ind(r, l, equal_var=False).pvalue) if len(r) >= 2 and len(l) >= 2 else np.nan
        mw = float(mannwhitneyu(r, l, alternative="two-sided").pvalue) if len(r) and len(l) else np.nan
        rows.append({
            "cohort": scores["cohort"].iloc[0],
            "layer": layer,
            "metric": col,
            "n_right": len(r),
            "n_left": len(l),
            "right_mean": float(np.mean(r)) if len(r) else np.nan,
            "left_mean": float(np.mean(l)) if len(l) else np.nan,
            "right_minus_left": float(np.mean(r) - np.mean(l)) if len(r) and len(l) else np.nan,
            "standardized_mean_difference": standardized_mean_diff(r, l),
            "welch_p": welch,
            "mannwhitney_p": mw,
        })
    out = pd.DataFrame(rows)
    out["welch_fdr_within_cohort_layer"] = bh_adjust(out["welch_p"].tolist())
    return out


def stage_numeric(meta: pd.DataFrame) -> pd.Series:
    candidates = ["stage", "char_tnm_stage", "ajcc_pathologic_stage"]
    for col in candidates:
        if col in meta.columns:
            x = meta[col].astype(str).str.upper()
            vals = x.str.extract(r"([1-4]|IV|III|II|I)", expand=False)
            mapping = {"I": 1, "II": 2, "III": 3, "IV": 4}
            return pd.to_numeric(vals.replace(mapping), errors="coerce")
    return pd.Series(np.nan, index=meta.index)


def mmr_indicator(meta: pd.DataFrame) -> pd.Series:
    candidates = ["char_mmr_status", "mmr_status", "msi_status"]
    for col in candidates:
        if col in meta.columns:
            s = meta[col].astype(str).str.lower()
            return pd.Series(np.where(s.str.contains("dmmr|msi-h|msih|high"), 1.0,
                                      np.where(s.str.contains("pmmr|mss|stable"), 0.0, np.nan)), index=meta.index)
    return pd.Series(np.nan, index=meta.index)


def adjusted_side_models(scores: pd.DataFrame, columns: list[str], layer: str) -> pd.DataFrame:
    rows = []
    base = pd.DataFrame(index=scores.index)
    base["right"] = scores["side"].eq("right").astype(float)
    base["stage"] = stage_numeric(scores)
    base["mmr_high"] = mmr_indicator(scores)
    for col in columns:
        frame = base.copy()
        frame["y"] = pd.to_numeric(scores[col], errors="coerce")
        covars = ["right"]
        if frame["stage"].notna().sum() >= max(20, int(0.5 * len(frame))):
            covars.append("stage")
        if frame["mmr_high"].notna().sum() >= max(20, int(0.5 * len(frame))):
            covars.append("mmr_high")
        data = frame[["y"] + covars].dropna()
        if len(data) < 20 or data["right"].nunique() < 2:
            rows.append({"cohort": scores["cohort"].iloc[0], "layer": layer, "metric": col, "n": len(data), "estimable": False})
            continue
        design = sm.add_constant(data[covars])
        fit = sm.OLS(data["y"], design).fit(cov_type="HC3")
        rows.append({
            "cohort": scores["cohort"].iloc[0],
            "layer": layer,
            "metric": col,
            "n": int(len(data)),
            "estimable": True,
            "covariates": "+".join(covars),
            "right_beta": float(fit.params["right"]),
            "right_p_hc3": float(fit.pvalues["right"]),
            "r_squared": float(fit.rsquared),
        })
    return pd.DataFrame(rows)


def depmap_layer() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    require(DEPMAP_MODEL)
    require(DEPMAP_EFFECT)
    models = pd.read_csv(DEPMAP_MODEL, low_memory=False)
    coad = models[models["OncotreeCode"].eq("COAD")].copy()
    needed = ["ModelID", "CellLineName", "StrippedCellLineName", "OncotreeCode", "Sex", "PrimaryOrMetastasis", "SampleCollectionSite"]
    coad = coad[needed]
    coad["ModelID"] = coad["ModelID"].astype(str)
    coad = add_manual_sidedness(coad)

    panel = coad.copy()
    source_columns = {}
    for gene in DEPMAP_DEFENSE_GENES:
        values, source = load_gene_column(DEPMAP_EFFECT, gene, f"{gene}_gene_effect")
        source_columns[gene] = source
        panel = panel.merge(values, on="ModelID", how="left")

    # Standardize each dependency within all COAD models. More negative Chronos
    # means stronger KO fitness loss; therefore the index is negated so that a
    # larger value means stronger collective defense dependency.
    zcols = []
    for gene in DEPMAP_DEFENSE_GENES:
        zcol = f"{gene}_effect_z"
        panel[zcol] = zscore(panel[f"{gene}_gene_effect"])
        zcols.append(zcol)
    panel["ferroptosis_defense_dependency_index"] = -panel[zcols].mean(axis=1, skipna=True)
    panel["dependency_genes_available"] = panel[zcols].notna().sum(axis=1)
    panel.loc[panel["dependency_genes_available"] < 4, "ferroptosis_defense_dependency_index"] = np.nan

    rows = []
    metrics = [f"{g}_gene_effect" for g in DEPMAP_DEFENSE_GENES] + ["ferroptosis_defense_dependency_index"]
    for confidence_name, confidence in [("high", {"high"}), ("high_plus_medium", {"high", "medium"})]:
        d = panel[panel["side"].isin(["Right", "Left"]) & panel["side_confidence"].isin(confidence)]
        for metric in metrics:
            r = pd.to_numeric(d.loc[d["side"].eq("Right"), metric], errors="coerce").dropna().to_numpy()
            l = pd.to_numeric(d.loc[d["side"].eq("Left"), metric], errors="coerce").dropna().to_numpy()
            p = float(mannwhitneyu(r, l, alternative="two-sided").pvalue) if len(r) and len(l) else np.nan
            rows.append({
                "confidence_set": confidence_name,
                "metric": metric,
                "n_right": len(r),
                "n_left": len(l),
                "right_mean": float(np.mean(r)) if len(r) else np.nan,
                "left_mean": float(np.mean(l)) if len(l) else np.nan,
                "right_minus_left": float(np.mean(r) - np.mean(l)) if len(r) and len(l) else np.nan,
                "mannwhitney_p_two_sided": p,
                "interpretation": (
                    "For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. "
                    "For composite index, positive Right-minus-Left = stronger right-sided dependency."
                ),
            })
    stats = pd.DataFrame(rows)
    for conf in stats["confidence_set"].unique():
        mask = stats["confidence_set"].eq(conf)
        stats.loc[mask, "fdr_within_confidence_set"] = bh_adjust(stats.loc[mask, "mannwhitney_p_two_sided"].tolist())
    info = {"n_coad": len(coad), "dependency_columns": source_columns}
    return panel, stats, info


def markdown_table(frame: pd.DataFrame) -> str:
    """Render a compact Markdown table without requiring the tabulate package."""
    if frame.empty:
        return "_No rows._"

    def format_cell(value: object) -> str:
        if value is None:
            return ""
        if isinstance(value, (float, np.floating)):
            if not np.isfinite(value):
                return ""
            return f"{value:.4g}"
        text = str(value)
        if text.lower() in {"nan", "nat", "none"}:
            return ""
        return text.replace("|", "\\|").replace("\n", " ")

    headers = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(format_cell(value) for value in row) + " |")
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    gse_expr, gse_meta, gse_info = load_gse()
    tcga_expr, tcga_meta, tcga_info = load_tcga()
    tcga_missing_genes = sorted(set(ALL_TRANSCRIPT_GENES) - set(tcga_expr.columns))
    gse = build_scores(gse_expr, gse_meta, "GSE39582")
    tcga = build_scores(tcga_expr, tcga_meta, "TCGA-COAD")
    samples = pd.concat([gse, tcga], axis=0, ignore_index=False, sort=False)

    layer1_metrics = ["ferroptosis_driver_score", "ferroptosis_defense_score", "ferroptosis_net_propensity"]
    layer2_metrics = ["pufa_incorporation_score", "peroxide_generation_score", "antioxidant_buffering_score", "lipid_peroxidation_liability"]

    stats = []
    adjusted = []
    for cohort in [gse, tcga]:
        stats.append(sidedness_stats(cohort, layer1_metrics, "transcript_ferroptosis"))
        stats.append(sidedness_stats(cohort, layer2_metrics, "lipid_peroxidation"))
        adjusted.append(adjusted_side_models(cohort, layer1_metrics, "transcript_ferroptosis"))
        adjusted.append(adjusted_side_models(cohort, layer2_metrics, "lipid_peroxidation"))
    stats = pd.concat(stats, ignore_index=True)
    adjusted = pd.concat(adjusted, ignore_index=True)

    dep_panel, dep_stats, dep_info = depmap_layer()

    samples.to_csv(OUT / "crc_sidedness_ferroptosis_three_layer_sample_scores.csv", index=False)
    stats.to_csv(OUT / "crc_sidedness_ferroptosis_three_layer_sidedness.csv", index=False)
    adjusted.to_csv(OUT / "crc_sidedness_ferroptosis_three_layer_adjusted.csv", index=False)
    dep_panel.to_csv(OUT / "crc_sidedness_ferroptosis_three_layer_depmap_panel.csv", index=False)
    dep_stats.to_csv(OUT / "crc_sidedness_ferroptosis_three_layer_depmap_stats.csv", index=False)

    manifest = {
        "question": "Which side of CRC shows the stronger ferroptosis-related state/liability?",
        "layer1": {
            "driver_genes": FERROPTOSIS_DRIVER,
            "defense_genes": FERROPTOSIS_DEFENSE,
            "net_definition": "driver mean-z minus defense mean-z",
        },
        "layer2": {
            "pufa_incorporation": PUFA_INCORPORATION,
            "peroxide_generation": PEROXIDE_GENERATION,
            "antioxidant_buffering": ANTIOXIDANT_BUFFERING,
            "liability_definition": "mean(PUFA incorporation, peroxide generation) minus antioxidant buffering",
        },
        "layer3": {
            "depmap_defense_genes": DEPMAP_DEFENSE_GENES,
            "chronos_direction": "more negative = larger fitness loss after knockout",
            "composite_direction": "larger composite index = stronger ferroptosis-defense dependency",
        },
        "cohorts": {
            "GSE39582": gse_info,
            "TCGA-COAD": {**tcga_info, "missing_target_genes": tcga_missing_genes},
            "DepMap": dep_info,
        },
        "guardrail": "No transcriptomic or dependency score is a direct ferroptotic-death measurement.",
    }
    (OUT / "crc_sidedness_ferroptosis_three_layer_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")

    lines = [
        "# CRC sidedness: three-layer ferroptosis screen",
        "",
        "## Question",
        "Which side shows the stronger ferroptosis-related state/liability before any SLC7A11 mechanism is assumed?",
        "",
        "## Data availability note",
        (
            "The TCGA targeted-expression cache lacks the following requested genes: "
            f"{', '.join(tcga_missing_genes) if tcga_missing_genes else 'none'}. "
            "Accordingly, TCGA metrics whose minimum gene-availability rule is not met are not estimable and must not be interpreted as null sidedness results."
        ),
        "",
        "## Layer 1 — transcriptomic ferroptosis state",
        "Positive Right−Left for `ferroptosis_net_propensity` supports a more ferroptosis-prone right-sided transcriptional state; defense and driver components are reported separately.",
        "",
        markdown_table(stats[stats.layer.eq("transcript_ferroptosis")]),
        "",
        "## Layer 2 — lipid-peroxidation liability",
        "Positive Right−Left for `lipid_peroxidation_liability` supports greater right-sided PUFA/peroxide liability after subtracting antioxidant buffering.",
        "",
        markdown_table(stats[stats.layer.eq("lipid_peroxidation")]),
        "",
        "## Adjusted sensitivity models",
        "HC3 OLS uses sidedness plus stage and/or MMR/MSI only when those covariates are sufficiently available in the cohort.",
        "",
        markdown_table(adjusted),
        "",
        "## Layer 3 — DepMap functional dependency",
        "For individual defense genes, a negative Right−Left Chronos effect means stronger right-sided dependency. For the composite index, positive Right−Left means stronger right-sided collective dependency.",
        "",
        markdown_table(dep_stats),
        "",
        "## Decision rule",
        "Do not label either side as having 'stronger ferroptosis' from one score. A sidedness conclusion requires directional agreement across the net transcriptomic propensity, lipid-peroxidation liability, and functional dependency layers. Discordant layers are biologically informative and should be reported as such.",
        "",
        "## Guardrail",
        "These analyses establish ferroptosis-related state/liability, not direct cell death. Direct C11-BODIPY/lipid-ROS plus Fer-1 or Lip-1 rescue remains the functional endpoint.",
    ]
    (OUT / "crc_sidedness_ferroptosis_three_layer_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "gse_n": len(gse),
        "tcga_n": len(tcga),
        "depmap_n": dep_info["n_coad"],
        "outputs": [
            "crc_sidedness_ferroptosis_three_layer_report.md",
            "crc_sidedness_ferroptosis_three_layer_sidedness.csv",
            "crc_sidedness_ferroptosis_three_layer_adjusted.csv",
            "crc_sidedness_ferroptosis_three_layer_depmap_stats.csv",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
