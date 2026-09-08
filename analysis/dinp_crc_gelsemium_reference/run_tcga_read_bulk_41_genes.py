"""TCGA-READ bulk expression validation for the fresh 41-gene DINP–CRC set.

This runner is deliberately separate from the legacy 81-gene/9-gene analyses.
It uses the public UCSC Toil Xena TCGA/GTEx RSEM gene-TPM hub, restricts the
sample universe to TCGA Rectum Adenocarcinoma (READ), and compares READ primary
tumours with READ solid-tissue normals when available.  It is a disease-state
expression check, not evidence that DINP caused the observed expression shift.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

try:
    import xenaPython as xena
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "xenaPython is required. Install with: pip install xenaPython==1.0.14"
    ) from exc


ROOT = Path(__file__).resolve().parents[2]
FRESH_DIR = ROOT / "analysis" / "dinp_crc_gelsemium_reference"
OUTPUT = FRESH_DIR / "outputs" / "tcga_read_bulk_41_genes"
GENE_FILE = FRESH_DIR / "outputs" / "query_41_genes.csv"
EVIDENCE_FILE = FRESH_DIR / "outputs" / "evidence_audit_41_genes" / "dinp_crc_41_gene_evidence_matrix.csv"

XENA_HUB = "https://toil.xenahubs.net"
EXPRESSION_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"
DISEASE_LABEL = "Rectum Adenocarcinoma"
COHORT = "READ"


def bh(values: pd.Series) -> pd.Series:
    """Benjamini–Hochberg adjustment, retaining the input index."""
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


def decode_codes(field: str) -> list[str]:
    response = xena.field_codes(XENA_HUB, PHENOTYPE_DATASET, [field])[0]["code"]
    return str(response).split("\t")


def decode_value(value: object, codes: list[str]) -> str:
    if value is None or str(value).lower() in {"nan", "none"}:
        return "NaN"
    try:
        return codes[int(float(value))]
    except (ValueError, TypeError, IndexError):
        return str(value)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def tumor_group() -> str:
    return f"{COHORT}_primary_tumor"


def normal_group() -> str:
    return f"{COHORT}_solid_tissue_normal"


def load_genes() -> list[str]:
    if not GENE_FILE.exists():
        raise FileNotFoundError(f"Fresh 41-gene input is missing: {GENE_FILE}")
    genes = (
        pd.read_csv(GENE_FILE)["gene_symbol"]
        .astype(str)
        .str.strip()
        .replace("", np.nan)
        .dropna()
        .drop_duplicates()
        .tolist()
    )
    if len(genes) != 41:
        raise ValueError(f"Expected exactly 41 fresh genes; observed {len(genes)}")
    return genes


def load_sample_manifest() -> pd.DataFrame:
    samples = xena.dataset_samples(XENA_HUB, EXPRESSION_DATASET, None)
    fields = ["_study", "primary disease or tissue", "_sample_type"]
    field_values = {
        field: xena.dataset_probe_values(XENA_HUB, PHENOTYPE_DATASET, samples, [field])[1][0]
        for field in fields
    }
    code_maps = {field: decode_codes(field) for field in fields}
    metadata = pd.DataFrame({"sample_id": samples})
    for field in fields:
        metadata[field] = [decode_value(value, code_maps[field]) for value in field_values[field]]

    metadata["study"] = metadata["_study"]
    metadata["disease"] = metadata["primary disease or tissue"]
    metadata["sample_type"] = metadata["_sample_type"]
    metadata["group"] = "outside_scope"
    read_mask = metadata["study"].eq("TCGA") & metadata["disease"].eq(DISEASE_LABEL)
    metadata.loc[read_mask & metadata["sample_type"].eq("Primary Tumor"), "group"] = tumor_group()
    metadata.loc[read_mask & metadata["sample_type"].eq("Solid Tissue Normal"), "group"] = normal_group()
    metadata["include"] = metadata["group"].ne("outside_scope")
    metadata["dataset"] = EXPRESSION_DATASET
    metadata["patient_id"] = metadata["sample_id"].where(metadata["study"].eq("TCGA"), np.nan)
    metadata["patient_id"] = metadata["patient_id"].astype("string").str.split("-").str[:3].str.join("-")
    return metadata


def load_expression(sample_ids: list[str], genes: list[str]) -> pd.DataFrame:
    records = xena.dataset_gene_probe_avg(XENA_HUB, EXPRESSION_DATASET, sample_ids, genes)
    expression = pd.DataFrame(index=sample_ids)
    for record in records:
        gene = record["gene"]
        scores = record.get("scores", [[]])
        values = scores[0] if scores and scores[0] else []
        expression[gene] = values if len(values) == len(sample_ids) else np.nan
    expression.index.name = "sample_id"
    return expression.apply(pd.to_numeric, errors="coerce")


def cliff_delta(x: np.ndarray, y: np.ndarray) -> float:
    """P(X>Y)-P(X<Y), with ties contributing zero."""
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    comparisons = np.subtract.outer(x, y)
    return float((np.sum(comparisons > 0) - np.sum(comparisons < 0)) / comparisons.size)


def independent_stats(expression: pd.DataFrame, metadata: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    tumor_ids = metadata.loc[metadata["group"].eq(tumor_group()), "sample_id"]
    normal_ids = metadata.loc[metadata["group"].eq(normal_group()), "sample_id"]
    rows = []
    for gene in genes:
        tumor = expression.loc[expression.index.intersection(tumor_ids), gene].dropna().to_numpy(float)
        normal = expression.loc[expression.index.intersection(normal_ids), gene].dropna().to_numpy(float)
        if len(tumor) and len(normal):
            statistic, p_value = mannwhitneyu(tumor, normal, alternative="two-sided")
        else:
            statistic, p_value = np.nan, 1.0
        rows.append(
            {
                "comparison": f"{COHORT}_primary_tumor_vs_{COHORT}_solid_tissue_normal",
                "gene": gene,
                "tumor_n": int(len(tumor)),
                "normal_n": int(len(normal)),
                "tumor_median": float(np.median(tumor)) if len(tumor) else np.nan,
                "normal_median": float(np.median(normal)) if len(normal) else np.nan,
                "median_delta_tumor_minus_normal": float(np.median(tumor) - np.median(normal)) if len(tumor) and len(normal) else np.nan,
                "cliffs_delta_tumor_vs_normal": cliff_delta(tumor, normal),
                "mann_whitney_U": float(statistic) if np.isfinite(statistic) else np.nan,
                "p_value": float(p_value),
            }
        )
    result = pd.DataFrame(rows)
    result["BH_FDR_within_41_gene_family"] = bh(result["p_value"])
    return result


def paired_stats(expression: pd.DataFrame, metadata: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    scoped = metadata[metadata["group"].isin([tumor_group(), normal_group()])].copy()
    tumors = scoped[scoped["group"].eq(tumor_group())].drop_duplicates("patient_id")
    normals = scoped[scoped["group"].eq(normal_group())].drop_duplicates("patient_id")
    paired_ids = sorted(set(tumors["patient_id"].dropna()) & set(normals["patient_id"].dropna()))
    rows = []
    if not paired_ids:
        return pd.DataFrame(columns=["comparison", "gene", "paired_n", "median_delta_tumor_minus_normal", "p_value", "BH_FDR_within_41_gene_family"])
    tumor_map = tumors.set_index("patient_id")["sample_id"]
    normal_map = normals.set_index("patient_id")["sample_id"]
    for gene in genes:
        tumor_values = expression.loc[[tumor_map[p] for p in paired_ids], gene].to_numpy(float)
        normal_values = expression.loc[[normal_map[p] for p in paired_ids], gene].to_numpy(float)
        keep = np.isfinite(tumor_values) & np.isfinite(normal_values)
        differences = tumor_values[keep] - normal_values[keep]
        if len(differences) and not np.allclose(differences, 0):
            statistic, p_value = wilcoxon(differences, alternative="two-sided", method="auto")
        else:
            statistic, p_value = np.nan, 1.0
        rows.append(
            {
                "comparison": f"{COHORT}_patient_paired_primary_tumor_vs_solid_normal",
                "gene": gene,
                "paired_n": int(len(differences)),
                "tumor_median": float(np.median(tumor_values[keep])) if keep.any() else np.nan,
                "normal_median": float(np.median(normal_values[keep])) if keep.any() else np.nan,
                "median_delta_tumor_minus_normal": float(np.median(differences)) if len(differences) else np.nan,
                "wilcoxon_statistic": float(statistic) if np.isfinite(statistic) else np.nan,
                "p_value": float(p_value),
            }
        )
    result = pd.DataFrame(rows)
    result["BH_FDR_within_41_gene_family"] = bh(result["p_value"])
    return result


def write_report(metadata: pd.DataFrame, independent: pd.DataFrame, paired: pd.DataFrame, genes: list[str], expression: pd.DataFrame) -> dict:
    counts = metadata.loc[metadata["include"], "group"].value_counts().to_dict()
    significant = independent[independent["BH_FDR_within_41_gene_family"] < 0.05]
    direction = int((independent["median_delta_tumor_minus_normal"] > 0).sum())
    report = [
        f"# TCGA-{COHORT} bulk validation of the fresh 41-gene DINP–CRC intersection",
        "",
        "## Status",
        "",
        "This is a new bulk expression check using only the fresh 41 genes reconstructed from the published Gelsemium elegans–CRC GeneCards supplement and the frozen DINP multi-source matrix. It does not reuse the legacy 81-gene or 9-gene analyses.",
        "",
        f"- {COHORT} primary tumour samples: **{int(counts.get(tumor_group(), 0))}**",
        f"- {COHORT} solid-tissue normal samples: **{int(counts.get(normal_group(), 0))}**",
        f"- Genes queried: **{len(genes)}**",
        f"- Independent tumor–normal tests with BH-FDR < 0.05: **{len(significant)} / {len(genes)}**",
        f"- Positive median shifts (tumor > normal): **{direction} / {len(genes)}**",
        f"- Patient-matched pairs available: **{len(set(metadata.loc[metadata['group'].eq(tumor_group()), 'patient_id'].dropna()) & set(metadata.loc[metadata['group'].eq(normal_group()), 'patient_id'].dropna()))}**",
        "",
        "## Interpretation boundary",
        "",
        f"The analysis tests whether the 41-gene DINP–CRC intersection is expressed differently in TCGA {COHORT} tumour tissue than in {COHORT} solid-tissue normal tissue. It is a disease-state and tissue-context validation only; it does not establish DINP exposure, direct target binding, causality, or temporal direction.",
        "",
        "The primary inferential family is the 41 queried genes, with BH-FDR applied across all 41 independent Mann–Whitney tests. Patient-paired results, when available, are reported separately and are not used to rescue or re-rank the independent analysis.",
        "",
        "## Data source",
        "",
        f"- UCSC Toil Xena hub: `{XENA_HUB}`",
        f"- Expression dataset: `{EXPRESSION_DATASET}`",
        f"- Phenotype dataset: `{PHENOTYPE_DATASET}`",
        f"- Disease field: `primary disease or tissue == {DISEASE_LABEL}`",
        f"- Groups: `_sample_type == Primary Tumor` versus `_sample_type == Solid Tissue Normal`",
        "- Expression values are analyzed on the Xena-delivered gene-TPM scale; no cross-platform normalization is applied.",
        "",
        "## Outputs",
        "",
        f"- `tcga_{COHORT.lower()}_sample_manifest.csv`: complete Xena phenotype mapping and inclusion flag.",
        f"- `tcga_{COHORT.lower()}_expression_41_genes.csv`: downloaded expression matrix for the scoped {COHORT} samples.",
        f"- `tcga_{COHORT.lower()}_independent_gene_stats.csv`: independent tumor–normal comparison with 41-gene BH-FDR.",
        f"- `tcga_{COHORT.lower()}_paired_gene_stats.csv`: available patient-matched sensitivity comparison.",
        f"- `tcga_{COHORT.lower()}_bulk_manifest.json`: source, input hash and run metadata.",
        "",
        "## Caveat",
        "",
        f"{COHORT} solid-tissue normal availability is limited in TCGA. If the normal reference is small, the estimate is treated as a precision-limited tissue comparison rather than a definitive universal CRC direction. The other TCGA CRC subtype, GTEx colon, and single-cell analyses are separate contexts and are not silently pooled into this {COHORT}-only run.",
    ]
    (OUTPUT / f"TCGA_{COHORT}_BULK_41_GENE_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return {
        "primary_tumor": int(counts.get(tumor_group(), 0)),
        "solid_tissue_normal": int(counts.get(normal_group(), 0)),
        "queried_genes": len(genes),
        "independent_fdr_positive": int(len(significant)),
        "positive_median_shifts": direction,
        "paired_pairs": int(len(set(metadata.loc[metadata["group"].eq(tumor_group()), "patient_id"].dropna()) & set(metadata.loc[metadata["group"].eq(normal_group()), "patient_id"].dropna()))),
        "expression_samples": int(expression.shape[0]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run fresh 41-gene TCGA READ or COAD bulk validation")
    parser.add_argument("--cohort", choices=["READ", "COAD"], default="READ")
    args = parser.parse_args()
    global COHORT, DISEASE_LABEL, OUTPUT
    COHORT = args.cohort
    DISEASE_LABEL = "Rectum Adenocarcinoma" if COHORT == "READ" else "Colon Adenocarcinoma"
    OUTPUT = FRESH_DIR / "outputs" / f"tcga_{COHORT.lower()}_bulk_41_genes"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    genes = load_genes()
    metadata = load_sample_manifest()
    scoped_ids = metadata.loc[metadata["include"], "sample_id"].tolist()
    if not scoped_ids:
        raise RuntimeError(f"No TCGA-{COHORT} samples were selected from the Xena phenotype table")
    expression = load_expression(scoped_ids, genes)
    metadata = metadata[metadata["sample_id"].isin(expression.index)].copy()
    expression = expression.reindex(metadata["sample_id"])
    independent = independent_stats(expression, metadata, genes)
    paired = paired_stats(expression, metadata, genes)

    # Attach fresh evidence labels if the audit table is available. This is
    # descriptive metadata only and never changes the expression statistics.
    if EVIDENCE_FILE.exists():
        evidence = pd.read_csv(EVIDENCE_FILE)
        evidence_columns = [
            "gene_symbol",
            "DINP_evidence_profile",
            "DINP_source_family_count",
            "DINP_multi_family_support",
            "pathway_anchor_class",
            "original_biological_priority_rank",
        ]
        columns = [c for c in evidence_columns if c in evidence.columns]
        if "gene_symbol" in columns:
            evidence_subset = evidence[columns].rename(columns={"gene_symbol": "gene"})
            independent = independent.merge(evidence_subset, on="gene", how="left")
            paired = paired.merge(evidence_subset, on="gene", how="left")

    prefix = f"tcga_{COHORT.lower()}"
    metadata.to_csv(OUTPUT / f"{prefix}_sample_manifest.csv", index=False)
    expression.to_csv(OUTPUT / f"{prefix}_expression_41_genes.csv")
    independent.to_csv(OUTPUT / f"{prefix}_independent_gene_stats.csv", index=False)
    paired.to_csv(OUTPUT / f"{prefix}_paired_gene_stats.csv", index=False)
    summary = write_report(metadata, independent, paired, genes, expression)

    manifest = {
        "analysis": f"TCGA-{COHORT} fresh 41-gene bulk expression validation",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(ROOT)),
        "xena_hub": XENA_HUB,
        "expression_dataset": EXPRESSION_DATASET,
        "phenotype_dataset": PHENOTYPE_DATASET,
        "cohort": COHORT,
        "disease_label": DISEASE_LABEL,
        "groups": {
            "tumor": f"TCGA + {DISEASE_LABEL} + Primary Tumor",
            "normal": f"TCGA + {DISEASE_LABEL} + Solid Tissue Normal",
        },
        "gene_input": str(GENE_FILE.relative_to(ROOT)),
        "gene_input_sha256": sha256_file(GENE_FILE),
        "gene_set_hash": sha256_text("|".join(genes)),
        "genes": genes,
        "summary": summary,
        "analysis_boundary": "bulk tumor-normal disease-state validation; not exposure causality or direct target confirmation",
        "python": platform.python_version(),
    }
    (OUTPUT / f"tcga_{COHORT.lower()}_bulk_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
