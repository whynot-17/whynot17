# DINP–CRC PPI/pathway/transcriptomic cross-ranking

## Scope

This prioritization cross-ranks the 97 DINP–CRC overlapping genes using three evidence dimensions:

- PPI hub score: percentile of the existing STRING hub ranking among mapped in-network genes.
- Pathway score: percentile of `log1p(representative_term_count)` across the 97 genes, using significant nonredundant GO/KEGG representative-term membership.
- Transcriptomic score: `(independent primary datasets significant at FDR < 0.05 / 4) × direction concordance fraction`.

The four independent primary datasets are TCGA-COAD, GSE74602, GSE10950, and GSE156355. TCGA paired and TCGA-vs-GTEx contrasts remain sensitivity fields and are excluded from the primary score.

## Composite score

`cross_rank_score = (ppi_hub_score × pathway_component_score × transcriptomic_component_score)^(1/3)`.

This is a transparent prioritization score, not a probability of causality and not a meta-analytic p-value. A missing STRING mapping or no significant representative pathway membership contributes a zero component; the raw metrics are retained for audit.

## Summary

- Unique input genes: 97
- PPI-mapped/in-network genes: 96
- Genes in the top-20 PPI hub set with at least one representative pathway term: 20
- Tier 1 cross-supported genes (top-20 PPI + pathway membership + ≥2/4 primary datasets significant): 16
- Primary transcriptomic replication ≥2/4: 69
- Primary transcriptomic replication ≥3/4: 46

Candidate-tier counts:

- Tier 1: cross-supported: 16
- Tier 2: replicated + pathway: 51
- Tier 2b: replicated, no representative pathway: 2
- Tier 3: PPI/pathway only: 4
- Tier 4: other: 24

## Top 20

| Rank | Gene | Tier | Score | PPI hub rank | Degree | Representative terms | Significant primary datasets | Consensus direction |
|---:|---|---|---:|---:|---:|---:|---:|---|
| 1 | IL1B | Tier 1: cross-supported | 0.8737 | 6 | 29 | 26 | 3/4 | tumor_high |
| 2 | IL1A | Tier 1: cross-supported | 0.8282 | 15 | 16 | 21 | 3/4 | tumor_high |
| 3 | CD36 | Tier 2: replicated + pathway | 0.7780 | 23 | 13 | 17 | 3/4 | normal_high |
| 4 | STAT3 | Tier 1: cross-supported | 0.7702 | 5 | 29 | 28 | 2/4 | normal_high |
| 5 | ESR1 | Tier 1: cross-supported | 0.7594 | 14 | 17 | 14 | 3/4 | normal_high |
| 6 | RXRA | Tier 2: replicated + pathway | 0.7314 | 30 | 10 | 15 | 3/4 | normal_high |
| 7 | IL6 | Tier 1: cross-supported | 0.7149 | 1 | 38 | 33 | 2/4 | tumor_high |
| 8 | NOS3 | Tier 2: replicated + pathway | 0.7115 | 45 | 7 | 21 | 3/4 | tumor_high |
| 9 | TP53 | Tier 1: cross-supported | 0.7099 | 3 | 32 | 33 | 2/4 | tumor_high |
| 10 | PPARA | Tier 1: cross-supported | 0.7079 | 16 | 16 | 15 | 3/4 | normal_high |
| 11 | CEBPB | Tier 1: cross-supported | 0.7055 | 17 | 15 | 12 | 3/4 | tumor_high |
| 12 | MMP9 | Tier 1: cross-supported | 0.7023 | 8 | 21 | 15 | 2/4 | tumor_high |
| 13 | IL18 | Tier 1: cross-supported | 0.6966 | 20 | 13 | 12 | 3/4 | normal_high |
| 14 | RELA | Tier 1: cross-supported | 0.6816 | 12 | 17 | 28 | 2/4 | tumor_high |
| 15 | PPARG | Tier 1: cross-supported | 0.6714 | 7 | 23 | 20 | 2/4 | normal_high |
| 16 | IFNG | Tier 1: cross-supported | 0.6664 | 9 | 20 | 20 | 2/4 | tumor_high |
| 17 | LEPR | Tier 2: replicated + pathway | 0.6580 | 43 | 7 | 14 | 3/4 | normal_high |
| 18 | FAS | Tier 2: replicated + pathway | 0.6560 | 54 | 4 | 18 | 3/4 | normal_high |
| 19 | SLC2A1 | Tier 2: replicated + pathway | 0.6542 | 53 | 4 | 17 | 3/4 | tumor_high |
| 20 | ADIPOQ | Tier 1: cross-supported | 0.6525 | 18 | 15 | 14 | 2/4 | normal_high |

## Interpretation guardrails

- The rank rewards convergence across network centrality, enriched pathway membership, and reproducible tumor-versus-normal expression differences.
- Transcriptomic validation is disease-state association; it does not establish that DINP caused the expression change.
- PPI degree and pathway-term counts are database/network-dependent and should not be interpreted as independent biological replicates.
- Direction is retained as consensus context; mixed-direction genes should be reviewed before MR instrument selection.
