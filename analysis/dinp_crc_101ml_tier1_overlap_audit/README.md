# DINP–CRC stable-ML × Tier 1 overlap and enrichment audit

This module compares the 17 default stable-ML genes with the 16 pre-existing Tier 1 PPI/pathway/transcriptomic cross-supported genes using the frozen 97-gene DINP–CRC overlap universe.

## Primary test

The 2×2 table is evaluated with Fisher's exact test. The audit reports two-sided, overlap-enrichment (`greater`), and overlap-deficit (`less`) p-values, expected overlap, observed/expected fold enrichment, odds ratio, Jaccard index, and the explicit overlap gene list. The 97-gene universe is used for every comparison.

Primary result:

- Stable-ML set: 17 genes
- Tier 1 set: 16 genes
- Observed overlap: 1 gene (`CEBPB`)
- Expected overlap: 2.804 genes
- Observed/expected: 0.357
- Odds ratio: 0.271
- Fisher two-sided p: 0.290
- Fisher overlap-greater p: 0.966
- Fisher overlap-less p: 0.176

## Pairwise audit

In addition to stable-ML × Tier 1, the pairwise table audits overlap with:

- PPI top-20 hub + pathway membership
- transcriptomic replication in ≥2/4 independent primary datasets
- top 20 genes from the PPI/pathway/transcriptomic cross-ranking

Two-sided p-values across the 10 pairwise comparisons are Benjamini–Hochberg adjusted.

## Files

- `DINP_CRC_101ML_Tier1_fisher_exact.csv`: primary stable-ML × Tier 1 Fisher result.
- `DINP_CRC_101ML_overlap_pairwise_audit.csv`: all pairwise overlap/enrichment comparisons.
- `DINP_CRC_101ML_Tier1_overlap_membership.csv`: 97-row gene-level membership audit.
- `DINP_CRC_101ML_Tier1_overlap_summary.md`: human-readable report.
- `DINP_CRC_101ML_Tier1_overlap_log.md`: input definitions and QC metadata.
- `run_101ml_tier1_overlap_audit.py`: reproducible script.

This is a set-overlap audit, not a causal inference, model-performance test, or MR-validity test.
