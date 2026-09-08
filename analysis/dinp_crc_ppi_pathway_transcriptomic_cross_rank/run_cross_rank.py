from __future__ import annotations

from pathlib import Path
import json
import math

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT


def rank_score_desc(values: pd.Series) -> pd.Series:
    """Convert a descending rank to a stable percentile in (0, 1]."""
    numeric = pd.to_numeric(values, errors="coerce")
    n = int(numeric.notna().sum())
    if n <= 0:
        return pd.Series(0.0, index=values.index)
    ranks = numeric.rank(method="average", ascending=False, na_option="bottom")
    score = 1.0 - ((ranks - 1.0) / float(n))
    return score.fillna(0.0).clip(lower=0.0, upper=1.0)


def require_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def main() -> None:
    overlap_path = OUTPUTS / "DINP_CRC_overlap.csv"
    hub_path = OUTPUTS / "DINP_CRC_STRING_hub_ranking.csv"
    membership_path = OUTPUTS / "DINP_CRC_pathway_membership.csv"
    validation_path = OUTPUTS / "DINP_CRC_TCGA_GEO_cross_dataset_summary.csv"

    overlap = pd.read_csv(overlap_path, dtype=str).fillna("")
    hub = pd.read_csv(hub_path, dtype=str).fillna("")
    membership = pd.read_csv(membership_path, dtype=str).fillna("")
    validation = pd.read_csv(validation_path, dtype=str).fillna("")

    require_columns(overlap, ["gene_symbol"], "CRC/DINP overlap")
    require_columns(
        hub,
        [
            "gene_symbol",
            "hub_rank",
            "mapping_status",
            "in_network",
            "degree",
            "weighted_degree",
            "betweenness_centrality",
            "closeness_centrality",
            "eigenvector_centrality",
            "pagerank",
        ],
        "STRING hub ranking",
    )
    require_columns(
        membership,
        ["gene_symbol", "analysis_category", "term", "description", "fdr_numeric"],
        "pathway membership",
    )
    require_columns(
        validation,
        [
            "gene_symbol",
            "independent_datasets_measured",
            "independent_datasets_significant_fdr_lt_0_05",
            "independent_datasets_tumor_high",
            "independent_datasets_normal_high",
            "direction_concordance_fraction",
            "consensus_direction",
            "significant_datasets",
            "significant_directions",
            "sensitivity_contrasts_significant_fdr_lt_0_05",
        ],
        "TCGA/GEO cross-dataset summary",
    )

    genes = (
        overlap[["gene_symbol"]]
        .assign(gene_symbol=lambda frame: frame["gene_symbol"].str.strip().str.upper())
        .query("gene_symbol != ''")
        .drop_duplicates("gene_symbol")
    )

    # PPI metrics: the hub ranking already integrates the six centrality measures.
    # The score is a percentile among mapped nodes; unmapped genes are explicitly 0.
    hub = hub.copy()
    hub["gene_symbol"] = hub["gene_symbol"].str.strip().str.upper()
    for column in [
        "hub_rank",
        "degree",
        "weighted_degree",
        "betweenness_centrality",
        "closeness_centrality",
        "eigenvector_centrality",
        "pagerank",
    ]:
        hub[column] = pd.to_numeric(hub[column], errors="coerce")
    hub = hub.sort_values(["hub_rank", "gene_symbol"], na_position="last").drop_duplicates(
        "gene_symbol", keep="first"
    )
    mapped = hub["mapping_status"].eq("mapped") & hub["in_network"].eq("yes")
    hub["ppi_hub_score"] = 0.0
    # Lower STRING hub_rank is stronger evidence, so invert it before applying
    # the helper whose convention is "larger value = stronger".
    hub.loc[mapped, "ppi_hub_score"] = rank_score_desc(-hub.loc[mapped, "hub_rank"])

    # Pathway membership is based on significant, nonredundant representative terms.
    # Deduplication here makes the aggregation auditable even if an upstream file changes.
    membership["gene_symbol"] = membership["gene_symbol"].str.strip().str.upper()
    membership = membership[membership["gene_symbol"].ne("")].copy()
    membership = membership.drop_duplicates(["gene_symbol", "analysis_category", "term"])
    pathway_counts = (
        membership.groupby("gene_symbol", as_index=False)
        .agg(
            representative_term_count=("term", "nunique"),
            representative_description_count=("description", "nunique"),
            pathway_min_fdr=("fdr_numeric", lambda values: pd.to_numeric(values, errors="coerce").min()),
        )
    )
    category_counts = pd.crosstab(membership["gene_symbol"], membership["analysis_category"])
    category_counts = category_counts.rename(
        columns={
            "GO_BP": "GO_BP_term_count",
            "GO_CC": "GO_CC_term_count",
            "GO_MF": "GO_MF_term_count",
            "KEGG": "KEGG_term_count",
        }
    ).reset_index()
    pathway_counts = pathway_counts.merge(category_counts, on="gene_symbol", how="left")
    for column in ["GO_BP_term_count", "GO_CC_term_count", "GO_MF_term_count", "KEGG_term_count"]:
        if column not in pathway_counts.columns:
            pathway_counts[column] = 0
    pathway_counts["pathway_min_fdr"] = pd.to_numeric(pathway_counts["pathway_min_fdr"], errors="coerce")

    # Transcriptomic validation uses only the four independent primary units:
    # TCGA-COAD, GSE74602, GSE10950, and GSE156355. The two TCGA/GTEx-derived
    # contrasts are retained as sensitivity context but excluded from this score.
    validation["gene_symbol"] = validation["gene_symbol"].str.strip().str.upper()
    for column in [
        "independent_datasets_measured",
        "independent_datasets_significant_fdr_lt_0_05",
        "independent_datasets_tumor_high",
        "independent_datasets_normal_high",
        "direction_concordance_fraction",
        "sensitivity_contrasts_significant_fdr_lt_0_05",
    ]:
        validation[column] = pd.to_numeric(validation[column], errors="coerce")
    validation = validation.drop_duplicates("gene_symbol", keep="first")
    validation["primary_dataset_replication_fraction"] = (
        validation["independent_datasets_significant_fdr_lt_0_05"] / 4.0
    ).clip(lower=0.0, upper=1.0)
    validation["transcriptomic_component_score"] = (
        validation["primary_dataset_replication_fraction"]
        * validation["direction_concordance_fraction"].fillna(0.0).clip(lower=0.0, upper=1.0)
    ).fillna(0.0)

    result = genes.merge(
        overlap[
            [
                "gene_symbol",
                "protein_name",
                "source_databases",
                "DINP_evidence_grade",
                "DINP_evidence_partition",
                "DINP_database_source",
                "DINP_evidence_type",
                "DINP_independent_source_count",
            ]
        ].assign(gene_symbol=lambda frame: frame["gene_symbol"].str.strip().str.upper()),
        on="gene_symbol",
        how="left",
    ).merge(
        hub[
            [
                "gene_symbol",
                "mapping_status",
                "in_network",
                "hub_tier",
                "hub_rank",
                "degree",
                "weighted_degree",
                "betweenness_centrality",
                "closeness_centrality",
                "eigenvector_centrality",
                "pagerank",
                "ppi_hub_score",
            ]
        ],
        on="gene_symbol",
        how="left",
    )
    result = result.merge(
        pathway_counts[
            [
                "gene_symbol",
                "representative_term_count",
                "representative_description_count",
                "GO_BP_term_count",
                "GO_CC_term_count",
                "GO_MF_term_count",
                "KEGG_term_count",
                "pathway_min_fdr",
            ]
        ],
        on="gene_symbol",
        how="left",
    )
    result = result.merge(
        validation[
            [
                "gene_symbol",
                "independent_datasets_measured",
                "independent_datasets_significant_fdr_lt_0_05",
                "independent_datasets_tumor_high",
                "independent_datasets_normal_high",
                "direction_concordance_fraction",
                "consensus_direction",
                "significant_datasets",
                "significant_directions",
                "sensitivity_contrasts_significant_fdr_lt_0_05",
                "primary_dataset_replication_fraction",
                "transcriptomic_component_score",
            ]
        ],
        on="gene_symbol",
        how="left",
    )

    numeric_columns = [
        "hub_rank",
        "degree",
        "weighted_degree",
        "betweenness_centrality",
        "closeness_centrality",
        "eigenvector_centrality",
        "pagerank",
        "ppi_hub_score",
        "representative_term_count",
        "representative_description_count",
        "GO_BP_term_count",
        "GO_CC_term_count",
        "GO_MF_term_count",
        "KEGG_term_count",
        "pathway_min_fdr",
        "independent_datasets_measured",
        "independent_datasets_significant_fdr_lt_0_05",
        "independent_datasets_tumor_high",
        "independent_datasets_normal_high",
        "direction_concordance_fraction",
        "sensitivity_contrasts_significant_fdr_lt_0_05",
        "primary_dataset_replication_fraction",
        "transcriptomic_component_score",
    ]
    for column in numeric_columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    result["representative_term_count"] = result["representative_term_count"].fillna(0).astype(int)
    for column in ["GO_BP_term_count", "GO_CC_term_count", "GO_MF_term_count", "KEGG_term_count"]:
        result[column] = result[column].fillna(0).astype(int)
    result["independent_datasets_measured"] = result["independent_datasets_measured"].fillna(0).astype(int)
    result["independent_datasets_significant_fdr_lt_0_05"] = result[
        "independent_datasets_significant_fdr_lt_0_05"
    ].fillna(0).astype(int)
    result["independent_datasets_tumor_high"] = result["independent_datasets_tumor_high"].fillna(0).astype(int)
    result["independent_datasets_normal_high"] = result["independent_datasets_normal_high"].fillna(0).astype(int)
    result["sensitivity_contrasts_significant_fdr_lt_0_05"] = result[
        "sensitivity_contrasts_significant_fdr_lt_0_05"
    ].fillna(0).astype(int)
    result["ppi_hub_score"] = result["ppi_hub_score"].fillna(0.0)
    result["transcriptomic_component_score"] = result["transcriptomic_component_score"].fillna(0.0)
    result["primary_dataset_replication_fraction"] = result["primary_dataset_replication_fraction"].fillna(0.0)
    result["direction_concordance_fraction"] = result["direction_concordance_fraction"].fillna(0.0)

    # Compute pathway percentile after merging onto the complete 97-gene universe.
    # Genes with no representative term are then set to an explicit zero.
    result["pathway_component_score"] = rank_score_desc(
        np.log1p(result["representative_term_count"])
    )
    result.loc[result["representative_term_count"].eq(0), "pathway_component_score"] = 0.0

    result["cross_rank_score"] = np.cbrt(
        result["ppi_hub_score"]
        * result["pathway_component_score"]
        * result["transcriptomic_component_score"]
    )
    result["ppi_component_rank"] = result["ppi_hub_score"].rank(method="min", ascending=False).astype(int)
    result["pathway_component_rank"] = result["pathway_component_score"].rank(method="min", ascending=False).astype(int)
    result["transcriptomic_component_rank"] = result["transcriptomic_component_score"].rank(
        method="min", ascending=False
    ).astype(int)

    result["ppi_pathway_candidate"] = (
        result["hub_rank"].le(20) & result["representative_term_count"].gt(0)
    ).fillna(False)
    result["cross_support_flag"] = (
        result["ppi_pathway_candidate"]
        & result["independent_datasets_significant_fdr_lt_0_05"].ge(2)
    )
    result["candidate_tier"] = "Tier 4: other"
    result.loc[result["independent_datasets_significant_fdr_lt_0_05"].ge(2) & result["representative_term_count"].gt(0), "candidate_tier"] = "Tier 2: replicated + pathway"
    result.loc[result["ppi_pathway_candidate"] & result["independent_datasets_significant_fdr_lt_0_05"].lt(2), "candidate_tier"] = "Tier 3: PPI/pathway only"
    result.loc[result["cross_support_flag"], "candidate_tier"] = "Tier 1: cross-supported"
    result.loc[
        result["independent_datasets_significant_fdr_lt_0_05"].ge(2)
        & result["representative_term_count"].eq(0),
        "candidate_tier",
    ] = "Tier 2b: replicated, no representative pathway"

    result = result.sort_values(
        [
            "cross_rank_score",
            "transcriptomic_component_score",
            "ppi_hub_score",
            "pathway_component_score",
            "gene_symbol",
        ],
        ascending=[False, False, False, False, True],
    ).reset_index(drop=True)
    result.insert(0, "cross_rank", np.arange(1, len(result) + 1))

    output_columns = [
        "cross_rank",
        "gene_symbol",
        "protein_name",
        "source_databases",
        "DINP_evidence_grade",
        "DINP_evidence_partition",
        "DINP_database_source",
        "DINP_evidence_type",
        "DINP_independent_source_count",
        "candidate_tier",
        "cross_support_flag",
        "ppi_pathway_candidate",
        "cross_rank_score",
        "ppi_hub_score",
        "pathway_component_score",
        "transcriptomic_component_score",
        "ppi_component_rank",
        "pathway_component_rank",
        "transcriptomic_component_rank",
        "hub_rank",
        "hub_tier",
        "mapping_status",
        "in_network",
        "degree",
        "weighted_degree",
        "betweenness_centrality",
        "closeness_centrality",
        "eigenvector_centrality",
        "pagerank",
        "representative_term_count",
        "representative_description_count",
        "GO_BP_term_count",
        "GO_CC_term_count",
        "GO_MF_term_count",
        "KEGG_term_count",
        "pathway_min_fdr",
        "independent_datasets_measured",
        "independent_datasets_significant_fdr_lt_0_05",
        "independent_datasets_tumor_high",
        "independent_datasets_normal_high",
        "direction_concordance_fraction",
        "consensus_direction",
        "significant_datasets",
        "significant_directions",
        "sensitivity_contrasts_significant_fdr_lt_0_05",
        "primary_dataset_replication_fraction",
    ]
    result[output_columns].to_csv(OUTPUTS / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv", index=False)
    result[output_columns].head(20).to_csv(
        OUTPUTS / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank_top20.csv", index=False
    )

    tier_counts = result["candidate_tier"].value_counts().to_dict()
    summary_lines = [
        "# DINP–CRC PPI/pathway/transcriptomic cross-ranking",
        "",
        "## Scope",
        "",
        "This prioritization cross-ranks the 97 DINP–CRC overlapping genes using three evidence dimensions:",
        "",
        "- PPI hub score: percentile of the existing STRING hub ranking among mapped in-network genes.",
        "- Pathway score: percentile of `log1p(representative_term_count)` across the 97 genes, using significant nonredundant GO/KEGG representative-term membership.",
        "- Transcriptomic score: `(independent primary datasets significant at FDR < 0.05 / 4) × direction concordance fraction`.",
        "",
        "The four independent primary datasets are TCGA-COAD, GSE74602, GSE10950, and GSE156355. TCGA paired and TCGA-vs-GTEx contrasts remain sensitivity fields and are excluded from the primary score.",
        "",
        "## Composite score",
        "",
        "`cross_rank_score = (ppi_hub_score × pathway_component_score × transcriptomic_component_score)^(1/3)`.",
        "",
        "This is a transparent prioritization score, not a probability of causality and not a meta-analytic p-value. A missing STRING mapping or no significant representative pathway membership contributes a zero component; the raw metrics are retained for audit.",
        "",
        "## Summary",
        "",
        f"- Unique input genes: {len(result)}",
        f"- PPI-mapped/in-network genes: {int(result['in_network'].eq('yes').sum())}",
        f"- Genes in the top-20 PPI hub set with at least one representative pathway term: {int(result['ppi_pathway_candidate'].sum())}",
        f"- Tier 1 cross-supported genes (top-20 PPI + pathway membership + ≥2/4 primary datasets significant): {int(result['cross_support_flag'].sum())}",
        f"- Primary transcriptomic replication ≥2/4: {int(result['independent_datasets_significant_fdr_lt_0_05'].ge(2).sum())}",
        f"- Primary transcriptomic replication ≥3/4: {int(result['independent_datasets_significant_fdr_lt_0_05'].ge(3).sum())}",
        "",
        "Candidate-tier counts:",
        "",
    ]
    for tier in sorted(tier_counts):
        summary_lines.append(f"- {tier}: {tier_counts[tier]}")
    summary_lines.extend(
        [
            "",
            "## Top 20",
            "",
            "| Rank | Gene | Tier | Score | PPI hub rank | Degree | Representative terms | Significant primary datasets | Consensus direction |",
            "|---:|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for _, row in result.head(20).iterrows():
        hub_rank = "" if pd.isna(row["hub_rank"]) else str(int(row["hub_rank"]))
        summary_lines.append(
            f"| {int(row['cross_rank'])} | {row['gene_symbol']} | {row['candidate_tier']} | {row['cross_rank_score']:.4f} | {hub_rank} | {int(row['degree']) if not pd.isna(row['degree']) else ''} | {int(row['representative_term_count'])} | {int(row['independent_datasets_significant_fdr_lt_0_05'])}/4 | {row['consensus_direction']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Interpretation guardrails",
            "",
            "- The rank rewards convergence across network centrality, enriched pathway membership, and reproducible tumor-versus-normal expression differences.",
            "- Transcriptomic validation is disease-state association; it does not establish that DINP caused the expression change.",
            "- PPI degree and pathway-term counts are database/network-dependent and should not be interpreted as independent biological replicates.",
            "- Direction is retained as consensus context; mixed-direction genes should be reviewed before MR instrument selection.",
        ]
    )
    (OUTPUTS / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank_summary.md").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )

    log = {
        "inputs": {
            "overlap": str(overlap_path),
            "ppi_hub_ranking": str(hub_path),
            "pathway_membership": str(membership_path),
            "transcriptomic_summary": str(validation_path),
        },
        "input_gene_count": int(len(result)),
        "pathway_membership_rows_after_term_deduplication": int(len(membership)),
        "pathway_membership_gene_count": int(membership["gene_symbol"].nunique()),
        "independent_primary_dataset_count": 4,
        "independent_primary_datasets": ["TCGA-COAD", "GSE74602", "GSE10950", "GSE156355"],
        "sensitivity_contrasts_excluded_from_primary_score": ["TCGA-COAD_paired", "TCGA-COAD_vs_GTEx-colon"],
        "score_definition": {
            "ppi_hub_score": "descending percentile among mapped and in-network STRING genes",
            "pathway_component_score": "descending percentile of log1p(representative_term_count) across the 97 genes; zero if no representative term",
            "transcriptomic_component_score": "(independent primary datasets significant at FDR<0.05 / 4) * direction concordance fraction",
            "cross_rank_score": "cube root of product of the three component scores",
        },
        "qc": {
            "duplicate_gene_symbols_in_output": int(result["gene_symbol"].duplicated().sum()),
            "cross_rank_unique": int(result["cross_rank"].nunique()),
            "rank_score_range": {
                "ppi": [float(result["ppi_hub_score"].min()), float(result["ppi_hub_score"].max())],
                "pathway": [float(result["pathway_component_score"].min()), float(result["pathway_component_score"].max())],
                "transcriptomic": [float(result["transcriptomic_component_score"].min()), float(result["transcriptomic_component_score"].max())],
                "cross": [float(result["cross_rank_score"].min()), float(result["cross_rank_score"].max())],
            },
        },
    }
    (OUTPUTS / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank_log.md").write_text(
        "# DINP–CRC cross-ranking audit log\n\n```json\n"
        + json.dumps(log, indent=2, ensure_ascii=False)
        + "\n```\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
