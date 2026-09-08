# DINP–CRC PPI/pathway/transcriptomic cross-ranking

This module cross-ranks the 97 DINP–CRC overlapping genes using three evidence dimensions:

1. STRING PPI hub centrality (`DINP_CRC_STRING_hub_ranking.csv`).
2. Significant nonredundant GO/KEGG representative-term membership (`DINP_CRC_pathway_membership.csv`).
3. Replicated tumor-versus-normal expression validation across four independent primary datasets: TCGA-COAD, GSE74602, GSE10950, and GSE156355 (`DINP_CRC_TCGA_GEO_cross_dataset_summary.csv`).

The two TCGA-derived sensitivity contrasts (TCGA paired and TCGA-vs-GTEx colon) are retained in the ranked table but excluded from the primary transcriptomic score.

## Score

`cross_rank_score = (ppi_hub_score × pathway_component_score × transcriptomic_component_score)^(1/3)`

- `ppi_hub_score`: descending percentile among mapped, in-network STRING genes; lower STRING `hub_rank` is stronger.
- `pathway_component_score`: descending percentile of `log1p(representative_term_count)` across all 97 genes; zero when no representative term is present.
- `transcriptomic_component_score`: `(number of significant primary datasets / 4) × direction concordance fraction`.

This is a transparent prioritization score, not a causal probability or meta-analytic p-value. PPI and pathway metrics are network/database-dependent and transcriptomic validation is disease-state association rather than proof that DINP caused the expression change.

## Key result

- 97 unique genes ranked; no duplicate gene symbols.
- 16 genes are Tier 1 cross-supported: PPI top-20 hub, pathway membership, and significant in at least 2/4 independent primary transcriptomic datasets.
- 69/97 genes are significant in at least 2/4 independent primary transcriptomic datasets.

## Files

- `DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv`: complete ranked table with raw metrics, component scores, evidence fields, and tiers.
- `DINP_CRC_PPI_pathway_transcriptomic_cross_rank_top20.csv`: top 20 ranked genes.
- `DINP_CRC_PPI_pathway_transcriptomic_cross_rank_summary.md`: methods, summary, and top-20 table.
- `DINP_CRC_PPI_pathway_transcriptomic_cross_rank_log.md`: input paths, score definitions, and QC log.
- `run_cross_rank.py`: reproducible ranking script.

Run from this module directory with Python and pandas:

```bash
python run_cross_rank.py
```
