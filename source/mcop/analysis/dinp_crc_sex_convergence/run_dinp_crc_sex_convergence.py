#!/usr/bin/env python3
"""Compare the frozen DINP–CRC intersection in male versus female CRC tumors.

This is a tumor-state expression analysis, not an MCOP epidemiologic
interaction test. The exposure-side gene evidence remains upstream and is not
reconstructed from expression data.
"""
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
import xenaPython as xena
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
UPSTREAM = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs"
OUT = HERE / "outputs"
XENA_HUB = "https://toil.xenahubs.net"
EXPRESSION_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"
CRC_TISSUES = {"Colon Adenocarcinoma": "COAD", "Rectum Adenocarcinoma": "READ"}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").fillna(1.0).to_numpy(float)
    order = np.argsort(p)
    adjusted = np.ones(len(p), dtype=float)
    running = 1.0
    n = max(len(p), 1)
    for rank in range(len(p) - 1, -1, -1):
        idx = order[rank]
        running = min(running, p[idx] * n / (rank + 1))
        adjusted[idx] = min(running, 1.0)
    return pd.Series(adjusted, index=values.index)


def decode_codes(field: str) -> list[str]:
    raw = xena.field_codes(XENA_HUB, PHENOTYPE_DATASET, [field])[0]["code"]
    return str(raw).split("\t")


def decode_value(value: object, codes: list[str]) -> str:
    if value is None or str(value).lower() in {"nan", "none"}:
        return ""
    try:
        return codes[int(float(value))]
    except (ValueError, TypeError, IndexError):
        return str(value)


def load_frozen_genes() -> tuple[list[str], str]:
    path = UPSTREAM / "dinp_crc_intersection.csv"
    frame = pd.read_csv(path, dtype=str).fillna("")
    genes = sorted(set(frame["gene_symbol"].str.upper().str.strip()) - {""})
    if len(genes) != 81:
        raise ValueError(f"Expected the frozen 81-gene intersection, found {len(genes)}")
    return genes, sha256_file(path)


def load_sample_metadata() -> pd.DataFrame:
    samples = xena.dataset_samples(XENA_HUB, EXPRESSION_DATASET, None)
    fields = ["_gender", "_study", "primary disease or tissue", "_sample_type"]
    values = [xena.dataset_probe_values(XENA_HUB, PHENOTYPE_DATASET, samples, [field])[1][0] for field in fields]
    code_map = {field: decode_codes(field) for field in fields}
    metadata = pd.DataFrame({"sample_id": samples})
    for field, field_values in zip(fields, values):
        metadata[field] = [decode_value(value, code_map[field]) for value in field_values]
    metadata["tumor_type"] = metadata["primary disease or tissue"].map(CRC_TISSUES).fillna("")
    metadata["group"] = "outside_scope"
    in_scope = (
        metadata["_study"].eq("TCGA")
        & metadata["tumor_type"].ne("")
        & metadata["_sample_type"].eq("Primary Tumor")
        & metadata["_gender"].isin(["Male", "Female"])
    )
    metadata.loc[in_scope, "group"] = "TCGA_CRC_primary_tumor"
    return metadata


def load_expression(sample_ids: list[str], genes: list[str]) -> pd.DataFrame:
    values = xena.dataset_gene_probe_avg(XENA_HUB, EXPRESSION_DATASET, sample_ids, genes)
    expression = pd.DataFrame(index=sample_ids)
    expression.index.name = "sample_id"
    returned: set[str] = set()
    for record in values:
        gene = str(record.get("gene", "")).upper().strip()
        scores = record.get("scores", [[]])
        if gene and scores and scores[0] and len(scores[0]) == len(sample_ids):
            expression[gene] = pd.to_numeric(scores[0], errors="coerce")
            returned.add(gene)
    missing = sorted(set(genes) - returned)
    if missing:
        raise ValueError(f"Xena did not return frozen genes: {missing}")
    return expression.reindex(columns=genes)


def fit_adjusted_ols(data: pd.DataFrame, gene: str) -> dict[str, Any]:
    subset = data[["sex_male", "read_indicator", gene]].dropna().copy()
    if subset["sex_male"].nunique() < 2:
        return {"adjusted_n": len(subset), "adjusted_beta_male_minus_female": np.nan, "adjusted_se": np.nan, "adjusted_ci_low": np.nan, "adjusted_ci_high": np.nan, "adjusted_p_value": 1.0, "adjusted_status": "not_estimable"}
    X = sm.add_constant(subset[["sex_male", "read_indicator"]].astype(float), has_constant="add")
    fit = sm.OLS(subset[gene].astype(float), X).fit(cov_type="HC3")
    beta = float(fit.params["sex_male"])
    se = float(fit.bse["sex_male"])
    ci = fit.conf_int().loc["sex_male"]
    return {
        "adjusted_n": len(subset),
        "adjusted_beta_male_minus_female": beta,
        "adjusted_se": se,
        "adjusted_ci_low": float(ci.iloc[0]),
        "adjusted_ci_high": float(ci.iloc[1]),
        "adjusted_p_value": float(fit.pvalues["sex_male"]),
        "adjusted_status": "ok",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "source_records").mkdir(parents=True, exist_ok=True)
    genes, upstream_hash = load_frozen_genes()
    metadata = load_sample_metadata()
    included = metadata.loc[metadata["group"].eq("TCGA_CRC_primary_tumor"), "sample_id"].tolist()
    expression = load_expression(included, genes)
    cohort = metadata[metadata["sample_id"].isin(expression.index)].copy()
    cohort["sex_male"] = (cohort["_gender"] == "Male").astype(int)
    cohort["read_indicator"] = (cohort["tumor_type"] == "READ").astype(int)
    merged = cohort.set_index("sample_id").join(expression, how="inner")
    merged.reset_index().to_csv(OUT / "tcga_crc_81_gene_expression.csv", index=False)

    rows: list[dict[str, Any]] = []
    for gene in genes:
        male = merged.loc[merged["_gender"].eq("Male"), gene].dropna().to_numpy(float)
        female = merged.loc[merged["_gender"].eq("Female"), gene].dropna().to_numpy(float)
        statistic, mw_p = mannwhitneyu(male, female, alternative="two-sided")
        adjusted = fit_adjusted_ols(merged, gene)
        rows.append({
            "gene_symbol": gene,
            "male_n": len(male),
            "female_n": len(female),
            "male_median": float(np.median(male)),
            "female_median": float(np.median(female)),
            "median_delta_male_minus_female": float(np.median(male) - np.median(female)),
            "mann_whitney_U": float(statistic),
            "mann_whitney_p": float(mw_p),
            **adjusted,
        })
    result = pd.DataFrame(rows)
    result["BH_FDR_mann_whitney_81_genes"] = bh(result["mann_whitney_p"])
    result["BH_FDR_adjusted_OLS_81_genes"] = bh(result["adjusted_p_value"])
    result["direction"] = np.where(result["adjusted_beta_male_minus_female"] > 0, "higher_in_male", "lower_in_male")
    result = result.sort_values(["BH_FDR_adjusted_OLS_81_genes", "adjusted_p_value", "gene_symbol"]).reset_index(drop=True)
    result.insert(0, "rank_by_adjusted_OLS_FDR", np.arange(1, len(result) + 1))
    result.to_csv(OUT / "tcga_crc_sex_gene_results.csv", index=False)
    result.head(20).to_csv(OUT / "tcga_crc_sex_top_genes.csv", index=False)

    sample_audit = cohort[["sample_id", "_gender", "tumor_type", "primary disease or tissue", "_sample_type", "_study", "group"]].copy()
    sample_audit.to_csv(OUT / "tcga_crc_sex_sample_audit.csv", index=False)
    sample_audit.groupby(["_gender", "tumor_type"], dropna=False).size().reset_index(name="n").to_csv(OUT / "tcga_crc_sex_group_counts.csv", index=False)

    counts = {
        "TCGA_CRC_primary_tumor": int(len(cohort)),
        "male": int((cohort["_gender"] == "Male").sum()),
        "female": int((cohort["_gender"] == "Female").sum()),
        "COAD": int((cohort["tumor_type"] == "COAD").sum()),
        "READ": int((cohort["tumor_type"] == "READ").sum()),
        "genes_tested": len(genes),
        "nominal_p_lt_0_05": int((result["adjusted_p_value"] < 0.05).sum()),
        "adjusted_OLS_FDR_lt_0_05": int((result["BH_FDR_adjusted_OLS_81_genes"] < 0.05).sum()),
        "MW_FDR_lt_0_05": int((result["BH_FDR_mann_whitney_81_genes"] < 0.05).sum()),
    }
    manifest = {
        "analysis": "DINP–CRC sex-stratified molecular convergence",
        "generated_at": utc(),
        "upstream_input": "accessible 81-gene CTD × (GeneCards/Open Targets) intersection",
        "upstream_intersection_sha256": upstream_hash,
        "expression_source": {"hub": XENA_HUB, "dataset": EXPRESSION_DATASET},
        "phenotype_source": {"dataset": PHENOTYPE_DATASET, "fields": ["_gender", "_study", "primary disease or tissue", "_sample_type"]},
        "gene_family_size": len(genes),
        "genes": genes,
        "counts": counts,
        "primary_model": "male versus female expression difference with HC3-robust OLS adjustment for COAD versus READ; BH-FDR across the frozen 81-gene family",
        "secondary_model": "two-sided Mann–Whitney U without tissue adjustment; BH-FDR across the frozen 81-gene family",
        "interpretation_boundary": "tumor molecular-state sex comparison; not a sex-specific MCOP epidemiologic interaction or causal exposure analysis",
        "python": platform.python_version(),
        "xenaPython": version("xenaPython"),
    }
    (OUT / "sex_convergence_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    group_table = sample_audit.groupby(["_gender", "tumor_type"], dropna=False).size().reset_index(name="n")
    top = result.head(10)
    report = [
        "# DINP–CRC male-versus-female molecular convergence",
        "",
        f"Generated: `{manifest['generated_at']}`",
        "",
        "## Frozen interpretation",
        "",
        "This analysis compares the frozen 81-gene accessible DINP–CRC intersection between male and female TCGA CRC primary tumors. It is a tumor molecular-state comparison, not evidence that MCOP has a male-specific epidemiologic effect.",
        "",
        "## Cohort",
        "",
        f"- TCGA CRC primary tumors analyzed: **{counts['TCGA_CRC_primary_tumor']}** ({counts['male']} male, {counts['female']} female).",
        f"- Tissue composition: **{counts['COAD']} COAD** and **{counts['READ']} READ**.",
        "",
        "| Sex | Tumor type | N |",
        "|---|---|---:|",
    ]
    for _, row in group_table.iterrows():
        report.append(f"| {row['_gender']} | {row['tumor_type']} | {int(row['n'])} |")
    report += [
        "",
        "## Main result",
        "",
        f"- Frozen gene family: **{counts['genes_tested']} genes**.",
        f"- Adjusted OLS nominal P<0.05: **{counts['nominal_p_lt_0_05']} genes**; adjusted-OLS BH-FDR<0.05: **{counts['adjusted_OLS_FDR_lt_0_05']} genes**.",
        f"- Mann–Whitney BH-FDR<0.05: **{counts['MW_FDR_lt_0_05']} genes**.",
        "- The adjusted OLS effect is male minus female, with COAD/READ included as a tissue covariate and HC3-robust standard errors.",
        "",
        "### Top adjusted-OLS results",
        "",
        "| Rank | Gene | Male−female beta | 95% CI | P | BH-FDR |",
        "|---:|---|---:|---|---:|---:|",
    ]
    for _, row in top.iterrows():
        report.append(f"| {int(row['rank_by_adjusted_OLS_FDR'])} | `{row['gene_symbol']}` | {row['adjusted_beta_male_minus_female']:.4g} | {row['adjusted_ci_low']:.4g} to {row['adjusted_ci_high']:.4g} | {row['adjusted_p_value']:.4g} | {row['BH_FDR_adjusted_OLS_81_genes']:.4g} |")
    report += [
        "",
        "## Boundaries",
        "",
        "- Gene-level tests are descriptive molecular-state analyses; they do not prove sex-specific exposure susceptibility.",
        "- The 81-gene family was inherited unchanged from the upstream three-source convergence output.",
        "- No NHANES MCOP×sex interaction statistic is reinterpreted here.",
        "- Expression values use the Xena-delivered scale; no new exposure or disease outcome was used to redefine the 81-gene family.",
        "",
        "## Files",
        "",
        "- `tcga_crc_sex_gene_results.csv`: all 81 gene-level results.",
        "- `tcga_crc_sex_top_genes.csv`: top 20 adjusted-OLS rows.",
        "- `tcga_crc_81_gene_expression.csv`: expression and phenotype data used for the cohort.",
        "- `tcga_crc_sex_sample_audit.csv` and `tcga_crc_sex_group_counts.csv`: sample QC.",
        "- `sex_convergence_manifest.json`: frozen inputs and provenance.",
        "",
        "## Reproducibility",
        "",
        "- Script: `analysis/dinp_crc_sex_convergence/run_dinp_crc_sex_convergence.py`",
        f"- Xena hub: `{XENA_HUB}`",
        f"- Expression dataset: `{EXPRESSION_DATASET}`",
        f"- Phenotype dataset: `{PHENOTYPE_DATASET}`",
        f"- Upstream 81-gene intersection SHA-256: `{upstream_hash}`",
        "",
        "**Status: completed as a sex-stratified CRC molecular-state audit; no epidemiologic sex-specific MCOP claim is made.**",
    ]
    (OUT / "summary.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
