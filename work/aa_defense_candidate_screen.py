"""Screen candidate ferroptosis-defense axes against an AA-routing proxy in CRC.

This is a targeted transcriptomic screen in GSE39582 and TCGA-COAD.  It tests
the user's candidate genes (AIFM2/FSP1, GCH1, SCD, MBOAT1/2, GPX4, SLC7A11,
DHODH) in three bounded ways:

1. right-versus-left expression;
2. AA-routing-proxy high versus low expression;
3. continuous AA-routing coupling, especially within right-sided CRC.

There is no patient-level AA lipidomics in these inputs.  The AA-high label is
therefore a transcriptomic proxy and is never interpreted as measured AA.
No dependency or causal inference is performed here.
"""

from __future__ import annotations

import io
import json
import sys
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, ttest_ind

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
DATA = ROOT / "work" / "data"

sys.path.insert(0, str(ROOT / "work"))
from phase1_bulk_analysis import (  # noqa: E402
    parse_gpl_probe_map,
    read_gse_targeted,
    tcga_case_metadata,
)

GSE_MATRIX = DATA / "GSE39582_series_matrix.clean.txt.gz"
GPL570_ANNOT = DATA / "GPL570.annot.gz"
TCGA_CASES = DATA / "tcga_coad_cases.json"
TCGA_SEED_EXPRESSION = DATA / "tcga_coad_gene_expression_uqfpkm.tsv"
TCGA_EXPRESSION = DATA / "tcga_aa_defense_candidate_expression_uqfpkm.tsv"

CANDIDATES = ["AIFM2", "GCH1", "SCD", "MBOAT1", "MBOAT2", "GPX4", "SLC7A11", "DHODH"]
AA_PROXY_CORE = ["PLA2G4A", "PTGS2", "PTGES"]
AA_PROXY_EXPANDED = AA_PROXY_CORE + ["ALOX5", "ALOX5AP"]
ALL_TARGETS = list(dict.fromkeys(CANDIDATES + AA_PROXY_EXPANDED))

# Stable Ensembl IDs for the GDC gene_expression/values endpoint.
TCGA_GENE_IDS = {
    "AIFM2": "ENSG00000042286",
    "GCH1": "ENSG00000131979",
    "SCD": "ENSG00000099194",
    "MBOAT1": "ENSG00000172197",
    "MBOAT2": "ENSG00000143797",
    "GPX4": "ENSG00000167468",
    "SLC7A11": "ENSG00000151012",
    "DHODH": "ENSG00000102967",
    "PLA2G4A": "ENSG00000116711",
    "PTGS2": "ENSG00000073756",
    "PTGES": "ENSG00000148344",
    "ALOX5": "ENSG00000012779",
    "ALOX5AP": "ENSG00000132965",
}


def zscore(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / sd


def bh_adjust(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q.tolist()
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    q[ok] = restored
    return q.tolist()


def load_gse() -> tuple[pd.DataFrame, dict[str, object]]:
    probe_map, probe_counts = parse_gpl_probe_map(GPL570_ANNOT, ALL_TARGETS)
    expr, meta, info = read_gse_targeted(GSE_MATRIX, probe_map, ALL_TARGETS)
    meta = meta.copy()
    meta["side"] = meta["side"].astype(str).str.lower()
    keep = meta["side"].isin(["left", "right"]) & ~meta["is_normal_like"].fillna(False)
    expr = expr.loc[keep]
    meta = meta.loc[keep].copy()
    info = {**info, "probe_counts": probe_counts, "n_side_known_primary": int(len(expr))}
    return expr, meta, info


def tcga_case_ids() -> list[str]:
    header = pd.read_csv(TCGA_SEED_EXPRESSION, sep="\t", nrows=0)
    return [str(x) for x in header.columns[1:]]


def load_tcga_expression() -> tuple[pd.DataFrame, dict[str, object]]:
    case_ids = tcga_case_ids()
    if TCGA_EXPRESSION.exists() and TCGA_EXPRESSION.stat().st_size > 1000:
        text = TCGA_EXPRESSION.read_text(encoding="utf-8")
        source = "cached GDC gene_expression/values response"
    else:
        payload = {
            "case_ids": case_ids,
            "gene_ids": [TCGA_GENE_IDS[g] for g in ALL_TARGETS],
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
        TCGA_EXPRESSION.write_text(text, encoding="utf-8")
        source = "GDC gene_expression/values endpoint"

    matrix = pd.read_csv(io.StringIO(text), sep="\t")
    gene_col = matrix.columns[0]
    id_to_symbol = {v: k for k, v in TCGA_GENE_IDS.items()}
    matrix[gene_col] = matrix[gene_col].astype(str)
    matrix["symbol"] = matrix[gene_col].map(id_to_symbol)
    matrix = matrix[matrix["symbol"].notna()].copy().set_index("symbol")
    expr = matrix.drop(columns=[gene_col], errors="ignore").T.apply(pd.to_numeric, errors="coerce")
    expr.index.name = "case_id"
    info = {
        "path": str(TCGA_EXPRESSION),
        "source": source,
        "n_cases": int(len(expr)),
        "n_genes": int(len(expr.columns)),
        "missing_genes": sorted(set(ALL_TARGETS) - set(expr.columns)),
    }
    return expr, info


def add_proxy_scores(expr: pd.DataFrame, log2_input: bool, cohort: str, meta: pd.DataFrame) -> pd.DataFrame:
    x = expr.apply(pd.to_numeric, errors="coerce").copy()
    if not log2_input:
        x = np.log2(np.clip(x, a_min=0, a_max=None) + 1.0)
    z = pd.DataFrame({gene: zscore(x[gene]) for gene in x.columns}, index=x.index)
    out = meta.copy()
    out.index = expr.index
    for gene in CANDIDATES:
        out[gene] = z[gene] if gene in z else np.nan
    core = [g for g in AA_PROXY_CORE if g in z]
    expanded = [g for g in AA_PROXY_EXPANDED if g in z]
    out["AA_routing_proxy_core"] = z[core].mean(axis=1) if core else np.nan
    out["AA_routing_proxy_expanded"] = z[expanded].mean(axis=1) if expanded else np.nan
    out["cohort"] = cohort
    out["sample_id"] = out.index.astype(str)
    for proxy in ["AA_routing_proxy_core", "AA_routing_proxy_expanded"]:
        median = out[proxy].median()
        out[f"{proxy}_high"] = out[proxy] >= median
    return out


def compare_sidedness(data: pd.DataFrame, proxy_name: str) -> pd.DataFrame:
    rows = []
    left = data[data.side.eq("left")]
    right = data[data.side.eq("right")]
    for gene in CANDIDATES:
        l = pd.to_numeric(left[gene], errors="coerce").dropna()
        r = pd.to_numeric(right[gene], errors="coerce").dropna()
        if len(l) < 2 or len(r) < 2:
            p = np.nan
            mw_p = np.nan
        else:
            p = float(ttest_ind(r, l, equal_var=False).pvalue)
            mw_p = float(mannwhitneyu(r, l, alternative="two-sided").pvalue)
        rows.append({
            "cohort": data.cohort.iloc[0],
            "comparison": "Right_minus_Left",
            "proxy": proxy_name,
            "gene": gene,
            "n_left": int(len(l)),
            "n_right": int(len(r)),
            "right_minus_left_mean_z": float(r.mean() - l.mean()) if len(l) and len(r) else np.nan,
            "welch_p": p,
            "mannwhitney_p": mw_p,
        })
    out = pd.DataFrame(rows)
    out["welch_fdr_within_cohort"] = bh_adjust(out["welch_p"].tolist())
    out["mannwhitney_fdr_within_cohort"] = bh_adjust(out["mannwhitney_p"].tolist())
    return out


def compare_proxy_split(data: pd.DataFrame, proxy_name: str, scope: str) -> pd.DataFrame:
    d = data if scope == "all_side_known" else data[data.side.eq(scope)]
    high_col = f"{proxy_name}_high"
    rows = []
    for gene in CANDIDATES:
        high = pd.to_numeric(d.loc[d[high_col], gene], errors="coerce").dropna()
        low = pd.to_numeric(d.loc[~d[high_col], gene], errors="coerce").dropna()
        p = float(ttest_ind(high, low, equal_var=False).pvalue) if len(high) >= 2 and len(low) >= 2 else np.nan
        rows.append({
            "cohort": data.cohort.iloc[0],
            "scope": scope,
            "comparison": "proxy_high_minus_low",
            "proxy": proxy_name,
            "gene": gene,
            "n_high": int(len(high)),
            "n_low": int(len(low)),
            "high_minus_low_mean_z": float(high.mean() - low.mean()) if len(high) and len(low) else np.nan,
            "welch_p": p,
        })
    out = pd.DataFrame(rows)
    out["welch_fdr_within_cohort_scope"] = bh_adjust(out["welch_p"].tolist())
    return out


def correlation_table(data: pd.DataFrame, proxy_name: str, scope: str) -> pd.DataFrame:
    d = data if scope == "all_side_known" else data[data.side.eq(scope)]
    rows = []
    for gene in CANDIDATES:
        pair = d[[proxy_name, gene]].dropna()
        if len(pair) >= 4 and pair[proxy_name].nunique() > 1 and pair[gene].nunique() > 1:
            rho, p = spearmanr(pair[proxy_name], pair[gene])
        else:
            rho, p = np.nan, np.nan
        rows.append({
            "cohort": data.cohort.iloc[0],
            "scope": scope,
            "proxy": proxy_name,
            "gene": gene,
            "n": int(len(pair)),
            "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
            "spearman_p": float(p) if np.isfinite(p) else np.nan,
        })
    out = pd.DataFrame(rows)
    out["spearman_fdr_within_cohort_scope"] = bh_adjust(out["spearman_p"].tolist())
    return out


def write_figure(sidedness: pd.DataFrame, correlations: pd.DataFrame) -> None:
    primary_side = sidedness[sidedness.proxy.eq("AA_routing_proxy_core")]
    primary_corr = correlations[
        correlations.proxy.eq("AA_routing_proxy_core") & correlations.scope.eq("right")
    ]
    cohorts = ["GSE39582", "TCGA-COAD"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), constrained_layout=True)

    pivot = primary_side.pivot(index="gene", columns="cohort", values="right_minus_left_mean_z").reindex(CANDIDATES)
    im = axes[0].imshow(pivot[cohorts].to_numpy(dtype=float), cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    axes[0].set_xticks(range(len(cohorts)), cohorts, rotation=25, ha="right")
    axes[0].set_yticks(range(len(CANDIDATES)), CANDIDATES)
    axes[0].set_title("Right − left expression")
    axes[0].set_xlabel("Cohort")
    axes[0].set_ylabel("Candidate axis")
    for i in range(len(CANDIDATES)):
        for j in range(len(cohorts)):
            value = pivot.loc[CANDIDATES[i], cohorts[j]]
            if pd.notna(value):
                axes[0].text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="mean within-cohort z")

    corr_pivot = primary_corr.pivot(index="gene", columns="cohort", values="spearman_rho").reindex(CANDIDATES)
    corr_pivot = corr_pivot.reindex(columns=cohorts)
    im2 = axes[1].imshow(corr_pivot.to_numpy(dtype=float), cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    axes[1].set_xticks(range(len(cohorts)), cohorts, rotation=25, ha="right")
    axes[1].set_yticks(range(len(CANDIDATES)), CANDIDATES)
    axes[1].set_xlabel("Cohort")
    axes[1].set_title("AA-routing coupling within right-sided CRC")
    for i in range(len(CANDIDATES)):
        for j in range(len(cohorts)):
            value = corr_pivot.loc[CANDIDATES[i], cohorts[j]]
            if pd.notna(value):
                axes[1].text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04, label="Spearman ρ")
    fig.suptitle("CRC AA-routing proxy and candidate ferroptosis-defense axes", fontsize=14)
    fig.savefig(OUT / "aa_defense_candidate_screen.png", dpi=220)
    plt.close(fig)


def write_report(sidedness: pd.DataFrame, splits: pd.DataFrame, correlations: pd.DataFrame, manifests: dict[str, object]) -> None:
    lines = [
        "# AA-routing proxy versus candidate ferroptosis-defense axes",
        "",
        "## Scope",
        "",
        "This is a targeted transcriptomic screen of AIFM2/FSP1, GCH1, SCD, MBOAT1, MBOAT2, GPX4, SLC7A11 and DHODH in GSE39582 and TCGA-COAD.",
        "It tests sidedness, proxy-high/low differences and continuous coupling. It does not measure tissue AA, lipid flux, ferroptotic death or genetic dependency.",
        "",
        "## AA-high definition",
        "",
        "The primary AA-routing proxy is the within-cohort mean z-score of PLA2G4A, PTGS2 and PTGES. The sensitivity proxy adds ALOX5 and ALOX5AP.",
        "Because patient-level AA lipidomics are unavailable, AA-high/low means high/low transcriptomic routing proxy, not measured AA concentration.",
        "",
        "## Primary right-sided results",
        "",
    ]
    primary = sidedness[sidedness.proxy.eq("AA_routing_proxy_core")].copy()
    lines.append("| Cohort | Candidate | Left n | Right n | Right−left mean z | Welch P | FDR |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for _, row in primary.iterrows():
        lines.append(f"| {row.cohort} | {row.gene} | {row.n_left} | {row.n_right} | {row.right_minus_left_mean_z:.3f} | {row.welch_p:.3g} | {row.welch_fdr_within_cohort:.3g} |")
    lines += [
        "",
        "## Primary coupling within right-sided CRC",
        "",
        "| Cohort | Candidate | n | Spearman ρ | P | FDR |",
        "|---|---|---:|---:|---:|---:|",
    ]
    right_corr = correlations[(correlations.proxy.eq("AA_routing_proxy_core")) & correlations.scope.eq("right")]
    for _, row in right_corr.iterrows():
        lines.append(f"| {row.cohort} | {row.gene} | {row.n} | {row.spearman_rho:.3f} | {row.spearman_p:.3g} | {row.spearman_fdr_within_cohort_scope:.3g} |")
    lines += [
        "",
        "## Primary AA-routing proxy high versus low",
        "",
        "The split analysis is descriptive and complements the continuous correlation; it should not replace the continuous test.",
        "",
        "| Cohort | Scope | Candidate | High n | Low n | High−low mean z | Welch P | FDR |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    primary_split = splits[splits.proxy.eq("AA_routing_proxy_core")]
    for _, row in primary_split.iterrows():
        lines.append(f"| {row.cohort} | {row.scope} | {row.gene} | {row.n_high} | {row.n_low} | {row.high_minus_low_mean_z:.3f} | {row.welch_p:.3g} | {row.welch_fdr_within_cohort_scope:.3g} |")
    lines += [
        "",
        "## Interpretation rule",
        "",
        "A candidate becomes interesting only if its direction is reasonably reproducible across cohorts and its continuous coupling is positive within right-sided tumors. A positive RNA association is still not proof of AA-induced defense or dependency.",
        "The next step, if a candidate survives this screen, is a context-specific functional test rather than another expression-only expansion.",
        "",
        "## Sensitivity analysis",
        "",
        "The expanded AA-routing proxy results are saved in the machine-readable tables. Adding ALOX5/ALOX5AP is interpreted cautiously because these genes can be strongly contributed by myeloid cells in bulk CRC.",
        "",
        "## Provenance",
        "",
        f"- GSE39582: {manifests['gse']}",
        f"- TCGA-COAD: {manifests['tcga']}",
        "- Expression values were converted to within-cohort z-scores before score construction.",
        "- Statistical unit: bulk tumor sample/case.",
        "- FDR is controlled separately within each cohort/comparison family across the eight pre-specified candidates.",
        "",
        "## Files",
        "",
        "- `aa_defense_candidate_sample_scores.csv`",
        "- `aa_defense_candidate_sidedness.csv`",
        "- `aa_defense_candidate_proxy_splits.csv`",
        "- `aa_defense_candidate_coupling.csv`",
        "- `aa_defense_candidate_screen.png`",
    ]
    (OUT / "aa_defense_candidate_screen_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gse_expr, gse_meta, gse_info = load_gse()
    tcga_expr, tcga_info = load_tcga_expression()
    tcga_meta = tcga_case_metadata(TCGA_CASES)
    common = tcga_expr.index.intersection(tcga_meta.index)
    tcga_expr = tcga_expr.loc[common]
    tcga_meta = tcga_meta.loc[common].copy()
    tcga_keep = tcga_meta["side"].isin(["left", "right"])
    tcga_expr = tcga_expr.loc[tcga_keep]
    tcga_meta = tcga_meta.loc[tcga_keep]

    gse_scores = add_proxy_scores(gse_expr, True, "GSE39582", gse_meta)
    tcga_scores = add_proxy_scores(tcga_expr, False, "TCGA-COAD", tcga_meta)
    combined = pd.concat([gse_scores, tcga_scores], axis=0, ignore_index=False, sort=False)
    combined.to_csv(OUT / "aa_defense_candidate_sample_scores.csv", index=False)

    sidedness = pd.concat(
        [compare_sidedness(gse_scores, "AA_routing_proxy_core"), compare_sidedness(gse_scores, "AA_routing_proxy_expanded"), compare_sidedness(tcga_scores, "AA_routing_proxy_core"), compare_sidedness(tcga_scores, "AA_routing_proxy_expanded")],
        ignore_index=True,
    )
    splits = []
    correlations = []
    for data in [gse_scores, tcga_scores]:
        for proxy in ["AA_routing_proxy_core", "AA_routing_proxy_expanded"]:
            for scope in ["all_side_known", "right"]:
                splits.append(compare_proxy_split(data, proxy, scope))
                correlations.append(correlation_table(data, proxy, scope))
    splits = pd.concat(splits, ignore_index=True)
    correlations = pd.concat(correlations, ignore_index=True)
    sidedness.to_csv(OUT / "aa_defense_candidate_sidedness.csv", index=False)
    splits.to_csv(OUT / "aa_defense_candidate_proxy_splits.csv", index=False)
    correlations.to_csv(OUT / "aa_defense_candidate_coupling.csv", index=False)
    write_figure(sidedness, correlations)

    manifest = {
        "analysis": "AA-routing proxy versus candidate ferroptosis-defense axes",
        "access_date": "2026-09-01",
        "candidates": CANDIDATES,
        "aa_proxy_core": AA_PROXY_CORE,
        "aa_proxy_expanded": AA_PROXY_EXPANDED,
        "aa_proxy_definition": "mean within-cohort z-score",
        "aa_high_definition": "within-cohort median split of transcriptomic proxy; not measured AA",
        "gse": {**gse_info, "n_analyzed": int(len(gse_scores)), "side_counts": gse_scores.side.value_counts().to_dict()},
        "tcga": {**tcga_info, "n_analyzed": int(len(tcga_scores)), "side_counts": tcga_scores.side.value_counts().to_dict()},
        "outputs": {
            "sample_scores": str(OUT / "aa_defense_candidate_sample_scores.csv"),
            "sidedness": str(OUT / "aa_defense_candidate_sidedness.csv"),
            "proxy_splits": str(OUT / "aa_defense_candidate_proxy_splits.csv"),
            "coupling": str(OUT / "aa_defense_candidate_coupling.csv"),
            "figure": str(OUT / "aa_defense_candidate_screen.png"),
            "report": str(OUT / "aa_defense_candidate_screen_report.md"),
        },
    }
    (OUT / "aa_defense_candidate_screen_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    write_report(sidedness, splits, correlations, manifest)

    primary_right_corr = correlations[
        correlations.proxy.eq("AA_routing_proxy_core") & correlations.scope.eq("right")
    ][["cohort", "gene", "n", "spearman_rho", "spearman_p", "spearman_fdr_within_cohort_scope"]]
    print(json.dumps({
        "gse": {"n": len(gse_scores), "sides": gse_scores.side.value_counts().to_dict()},
        "tcga": {"n": len(tcga_scores), "sides": tcga_scores.side.value_counts().to_dict()},
        "primary_right_coupling": primary_right_corr.to_dict("records"),
        "outputs": manifest["outputs"],
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
