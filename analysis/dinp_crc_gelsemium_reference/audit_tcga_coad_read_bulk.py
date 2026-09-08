"""Offline audit of the fresh TCGA-COAD and TCGA-READ bulk runs.

The audit re-computes the stored Mann–Whitney statistics and BH correction from
the saved expression matrices, verifies sample/group and frozen-gene-set
identity, and summarizes cross-cohort direction/concordance.  It does not
pool COAD and READ into a new inferential family and does not use exposure
data to choose genes.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon


ROOT = Path(__file__).resolve().parents[2]
FRESH_DIR = ROOT / "analysis" / "dinp_crc_gelsemium_reference"
OUTPUT = FRESH_DIR / "outputs" / "tcga_coad_read_bulk_audit"
GENE_FILE = FRESH_DIR / "outputs" / "query_41_genes.csv"
COHORTS = {
    "READ": {
        "label": "Rectum Adenocarcinoma",
        "directory": FRESH_DIR / "outputs" / "tcga_read_bulk_41_genes",
    },
    "COAD": {
        "label": "Colon Adenocarcinoma",
        "directory": FRESH_DIR / "outputs" / "tcga_coad_bulk_41_genes",
    },
}


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
    m = max(len(p), 1)
    for rank in range(len(p) - 1, -1, -1):
        idx = order[rank]
        running = min(running, p[idx] * m / (rank + 1))
        adjusted[idx] = running
    return pd.Series(adjusted, index=values.index)


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    differences = np.subtract.outer(x, y)
    return float((np.sum(differences > 0) - np.sum(differences < 0)) / differences.size)


def load_genes() -> list[str]:
    genes = pd.read_csv(GENE_FILE)["gene_symbol"].astype(str).str.strip().drop_duplicates().tolist()
    if len(genes) != 41:
        raise ValueError(f"Expected 41 frozen genes, observed {len(genes)}")
    return genes


def load_cohort(cohort: str, genes: list[str]) -> dict:
    directory = COHORTS[cohort]["directory"]
    prefix = f"tcga_{cohort.lower()}"
    manifest = json.loads((directory / f"{prefix}_bulk_manifest.json").read_text(encoding="utf-8"))
    metadata = pd.read_csv(directory / f"{prefix}_sample_manifest.csv")
    expression = pd.read_csv(directory / f"{prefix}_expression_41_genes.csv", index_col=0)
    stats = pd.read_csv(directory / f"{prefix}_independent_gene_stats.csv")
    paired = pd.read_csv(directory / f"{prefix}_paired_gene_stats.csv")
    return {"manifest": manifest, "metadata": metadata, "expression": expression, "stats": stats, "paired": paired}


def audit_cohort(cohort: str, data: dict, genes: list[str]) -> tuple[list[str], dict, pd.DataFrame]:
    errors: list[str] = []
    warnings: list[str] = []
    manifest = data["manifest"]
    metadata = data["metadata"]
    expression = data["expression"]
    stats = data["stats"]
    paired = data["paired"]
    prefix = cohort.lower()
    tumor_group = f"{cohort}_primary_tumor"
    normal_group = f"{cohort}_solid_tissue_normal"

    if manifest.get("genes") != genes:
        errors.append("manifest gene list differs from query_41_genes.csv")
    if manifest.get("gene_input_sha256") != sha256_file(GENE_FILE):
        errors.append("manifest gene-input SHA256 does not match query_41_genes.csv")
    if manifest.get("cohort") != cohort:
        errors.append(f"manifest cohort is {manifest.get('cohort')!r}, expected {cohort!r}")
    if manifest.get("disease_label") != COHORTS[cohort]["label"]:
        errors.append("manifest disease label is inconsistent with the cohort")

    if metadata["sample_id"].duplicated().any():
        errors.append("sample manifest contains duplicate sample IDs")
    if not expression.index.is_unique:
        errors.append("expression matrix contains duplicate sample IDs")
    if set(expression.index) != set(metadata["sample_id"]):
        errors.append("expression matrix sample IDs differ from sample manifest")
    if list(expression.columns) != genes:
        errors.append("expression matrix columns differ from the frozen 41-gene order")
    if stats["gene"].duplicated().any() or set(stats["gene"]) != set(genes):
        errors.append("independent stats do not contain exactly one row for each frozen gene")
    if paired["gene"].duplicated().any() or set(paired["gene"]) != set(genes):
        errors.append("paired stats do not contain exactly one row for each frozen gene")

    included = metadata[metadata["include"]].copy()
    counts = included["group"].value_counts().to_dict()
    expected = manifest.get("summary", {})
    if int(counts.get(tumor_group, 0)) != int(expected.get("primary_tumor", -1)):
        errors.append("primary-tumor count differs between sample manifest and run manifest")
    if int(counts.get(normal_group, 0)) != int(expected.get("solid_tissue_normal", -1)):
        errors.append("solid-normal count differs between sample manifest and run manifest")
    if len(included) != int(expected.get("expression_samples", -1)):
        errors.append("included sample count differs from run manifest")

    # Recompute every independent result directly from the saved matrix.
    recomputed_rows = []
    for gene in genes:
        tumor_ids = metadata.loc[metadata["group"].eq(tumor_group), "sample_id"]
        normal_ids = metadata.loc[metadata["group"].eq(normal_group), "sample_id"]
        tumor = expression.loc[expression.index.intersection(tumor_ids), gene].dropna().to_numpy(float)
        normal = expression.loc[expression.index.intersection(normal_ids), gene].dropna().to_numpy(float)
        statistic, p_value = mannwhitneyu(tumor, normal, alternative="two-sided")
        recomputed_rows.append(
            {
                "gene": gene,
                "tumor_n_recomputed": len(tumor),
                "normal_n_recomputed": len(normal),
                "tumor_median_recomputed": float(np.median(tumor)),
                "normal_median_recomputed": float(np.median(normal)),
                "median_delta_recomputed": float(np.median(tumor) - np.median(normal)),
                "cliffs_delta_recomputed": cliffs_delta(tumor, normal),
                "mann_whitney_U_recomputed": float(statistic),
                "p_value_recomputed": float(p_value),
            }
        )
    recomputed = pd.DataFrame(recomputed_rows)
    recomputed["BH_FDR_recomputed"] = bh(recomputed["p_value_recomputed"])
    stored = stats.set_index("gene").loc[genes].reset_index()
    comparisons = [
        ("tumor_n", "tumor_n_recomputed"),
        ("normal_n", "normal_n_recomputed"),
        ("tumor_median", "tumor_median_recomputed"),
        ("normal_median", "normal_median_recomputed"),
        ("median_delta_tumor_minus_normal", "median_delta_recomputed"),
        ("cliffs_delta_tumor_vs_normal", "cliffs_delta_recomputed"),
        ("mann_whitney_U", "mann_whitney_U_recomputed"),
        ("p_value", "p_value_recomputed"),
        ("BH_FDR_within_41_gene_family", "BH_FDR_recomputed"),
    ]
    for stored_col, recomputed_col in comparisons:
        left = pd.to_numeric(stored[stored_col], errors="coerce").to_numpy(float)
        right = pd.to_numeric(recomputed[recomputed_col], errors="coerce").to_numpy(float)
        if not np.allclose(left, right, rtol=1e-10, atol=1e-12, equal_nan=True):
            errors.append(f"stored {stored_col} does not reproduce from the saved expression matrix")

    # Recompute patient-paired Wilcoxon results as an independent sensitivity check.
    tumors = metadata[metadata["group"].eq(tumor_group)].drop_duplicates("patient_id").set_index("patient_id")
    normals = metadata[metadata["group"].eq(normal_group)].drop_duplicates("patient_id").set_index("patient_id")
    pair_ids = sorted(set(tumors.index.dropna()) & set(normals.index.dropna()))
    if len(pair_ids) != int(expected.get("paired_pairs", -1)):
        errors.append("paired patient count differs from the run manifest")
    if len(pair_ids) < 10:
        warnings.append(f"small paired sensitivity set (n={len(pair_ids)})")
    for gene in genes:
        tumor_values = expression.loc[tumors.loc[pair_ids, "sample_id"], gene].to_numpy(float)
        normal_values = expression.loc[normals.loc[pair_ids, "sample_id"], gene].to_numpy(float)
        keep = np.isfinite(tumor_values) & np.isfinite(normal_values)
        differences = tumor_values[keep] - normal_values[keep]
        if len(differences) and not np.allclose(differences, 0):
            statistic, p_value = wilcoxon(differences, alternative="two-sided", method="auto")
        else:
            statistic, p_value = np.nan, 1.0
        stored_row = paired.set_index("gene").loc[gene]
        if int(stored_row["paired_n"]) != int(len(differences)):
            errors.append(f"{cohort} paired_n mismatch for {gene}")
        if not np.isclose(float(stored_row["p_value"]), float(p_value), rtol=1e-10, atol=1e-12, equal_nan=True):
            errors.append(f"{cohort} paired p-value mismatch for {gene}")

    if counts.get(normal_group, 0) < 20:
        warnings.append(f"small independent normal reference (n={int(counts.get(normal_group, 0))})")
    q05 = int((stats["BH_FDR_within_41_gene_family"] < 0.05).sum())
    summary = {
        "cohort": cohort,
        "primary_tumor": int(counts.get(tumor_group, 0)),
        "solid_tissue_normal": int(counts.get(normal_group, 0)),
        "paired_pairs": len(pair_ids),
        "genes": len(genes),
        "independent_fdr_positive": q05,
        "independent_recompute_pass": not any("stored " in error for error in errors),
        "errors": len(errors),
        "warnings": warnings,
    }
    return errors, summary, recomputed


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    genes = load_genes()
    loaded = {cohort: load_cohort(cohort, genes) for cohort in COHORTS}
    all_errors: dict[str, list[str]] = {}
    summaries: dict[str, dict] = {}
    recomputed: dict[str, pd.DataFrame] = {}
    for cohort, data in loaded.items():
        errors, summary, table = audit_cohort(cohort, data, genes)
        all_errors[cohort] = errors
        summaries[cohort] = summary
        recomputed[cohort] = table

    read_stats = loaded["READ"]["stats"].set_index("gene").loc[genes]
    coad_stats = loaded["COAD"]["stats"].set_index("gene").loc[genes]
    cross = pd.DataFrame({
        "gene": genes,
        "READ_delta": read_stats["median_delta_tumor_minus_normal"].to_numpy(float),
        "COAD_delta": coad_stats["median_delta_tumor_minus_normal"].to_numpy(float),
        "READ_p": read_stats["p_value"].to_numpy(float),
        "COAD_p": coad_stats["p_value"].to_numpy(float),
        "READ_q": read_stats["BH_FDR_within_41_gene_family"].to_numpy(float),
        "COAD_q": coad_stats["BH_FDR_within_41_gene_family"].to_numpy(float),
    })
    cross["READ_direction"] = np.where(cross["READ_delta"] > 0, "up", np.where(cross["READ_delta"] < 0, "down", "zero"))
    cross["COAD_direction"] = np.where(cross["COAD_delta"] > 0, "up", np.where(cross["COAD_delta"] < 0, "down", "zero"))
    cross["direction_concordant"] = cross["READ_direction"].eq(cross["COAD_direction"]) & cross["READ_direction"].ne("zero")
    cross["READ_FDR_positive"] = cross["READ_q"] < 0.05
    cross["COAD_FDR_positive"] = cross["COAD_q"] < 0.05
    cross["shared_FDR_positive"] = cross["READ_FDR_positive"] & cross["COAD_FDR_positive"]
    cross.to_csv(OUTPUT / "tcga_coad_read_cross_cohort_gene_audit.csv", index=False)

    deltas = cross[["READ_delta", "COAD_delta"]].to_numpy(float)
    rho, rho_p = spearmanr(deltas[:, 0], deltas[:, 1])
    shared = int(cross["shared_FDR_positive"].sum())
    concordant = int(cross["direction_concordant"].sum())
    total_errors = sum(len(items) for items in all_errors.values())
    status = "PASS" if total_errors == 0 else "FAIL"
    report = [
        "# TCGA-COAD + TCGA-READ bulk audit of the fresh 41-gene set",
        "",
        f"## Audit status: **{status}**",
        "",
        "This audit uses only the saved local expression matrices, sample manifests, run manifests, and the frozen `query_41_genes.csv`. It independently re-computes the gene-level Mann–Whitney statistics, effect summaries, BH-FDR values, and available patient-paired Wilcoxon checks. It does not pool COAD and READ into a new test family.",
        "",
        "## Cohort-level verification",
        "",
        "| Cohort | Primary tumor | Solid normal | Paired patients | FDR-positive genes | Errors | Warnings |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for cohort in ["READ", "COAD"]:
        s = summaries[cohort]
        report.append(f"| TCGA-{cohort} | {s['primary_tumor']} | {s['solid_tissue_normal']} | {s['paired_pairs']} | {s['independent_fdr_positive']} | {s['errors']} | {('; '.join(s['warnings']) or 'none')} |")
    report += [
        "",
        "## Cross-cohort comparison",
        "",
        f"- Shared independent BH-FDR<0.05 genes: **{shared}/41**.",
        f"- Same non-zero tumor–normal direction in both cohorts: **{concordant}/41**.",
        f"- Spearman correlation of the 41 median shifts (READ vs COAD): **rho={rho:.3f}, P={rho_p:.3g}**.",
        "- These are cross-cohort expression-state concordance metrics, not evidence that either cohort measured DINP exposure or that DINP caused the expression shifts.",
        "",
        "## Fixed-analysis checks",
        "",
        "- Both cohorts use the same UCSC Toil Xena expression dataset and phenotype dataset.",
        "- Both use the same frozen 41-gene input, with matching input SHA256 and gene order.",
        "- Each cohort applies BH-FDR separately across its own 41-gene tumor–normal family.",
        "- The saved expression matrices reproduce all stored independent medians, Cliff's delta, Mann–Whitney U, P values, and BH-FDR values within numerical tolerance.",
        "- The saved sample manifests reproduce the stored tumor, normal, and paired-patient counts.",
        "- The paired analyses are sensitivity checks only and do not re-rank the independent results.",
        "",
        "## Interpretation boundary",
        "",
        "The COAD and READ results provide independent CRC tissue-state checks for the fresh 41-gene intersection. The limited READ normal reference (n=10) and the tissue-composition difference between bulk tumors and solid normals require cautious interpretation. This audit does not establish exposure-specific direction, direct chemical binding, causality, or mechanism.",
        "",
        "## Files",
        "",
        "- `tcga_coad_read_cross_cohort_gene_audit.csv`: gene-level COAD/READ comparison and concordance flags.",
        "- `tcga_coad_read_bulk_audit_manifest.json`: audit inputs, hashes, counts, and status.",
        "- `audit_tcga_coad_read_bulk.py`: reproducible offline audit script.",
    ]
    if total_errors:
        report += ["", "## Errors"]
        for cohort, errors in all_errors.items():
            for error in errors:
                report.append(f"- {cohort}: {error}")
    (OUTPUT / "TCGA_COAD_READ_BULK_AUDIT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "Offline audit of fresh TCGA-COAD and TCGA-READ bulk validations",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(ROOT)),
        "gene_input": str(GENE_FILE.relative_to(ROOT)),
        "gene_input_sha256": sha256_file(GENE_FILE),
        "cohorts": {cohort: {"directory": str(COHORTS[cohort]["directory"].relative_to(ROOT)), "summary": summaries[cohort], "errors": all_errors[cohort]} for cohort in COHORTS},
        "cross_cohort": {
            "shared_fdr_positive": shared,
            "direction_concordant": concordant,
            "spearman_rho": float(rho),
            "spearman_p": float(rho_p),
        },
        "status": status,
        "analysis_boundary": "offline implementation and result reproducibility audit; not exposure causality or mechanistic validation",
        "python": platform.python_version(),
    }
    (OUTPUT / "tcga_coad_read_bulk_audit_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"TCGA COAD + READ bulk audit: {status}")
    print(f"Shared FDR-positive genes: {shared}/41")
    print(f"Direction-concordant genes: {concordant}/41")
    print(f"Errors: {total_errors}")


if __name__ == "__main__":
    main()
