# TCGA-COAD + READ WGCNA audit

## Scope

This is a sensitivity audit of the full-transcriptome WGCNA. The canonical network remains the 8,000-gene signed-bicor analysis. The 15,000-gene and unsigned-Pearson runs are sensitivity analyses; none changes the primary result.

## Audit 1 — candidate eligibility and module placement

- Fresh target set: 41 genes.
- Canonical 8,000-gene network input: 16/41.
- 15,000-gene sensitivity network input: 26/41.
- Unsigned-Pearson 8,000-gene network input: 16/41.

The candidate-level CSV records full-transcriptome MAD rank, inclusion status, module, kME, and projection-only module for every target. Projection-only assignments are not counted as module membership or Fisher enrichment.

### Priority genes

| Gene | Full MAD rank | 8k module | 8k kME | 15k module | 15k kME | Unsigned module | Unsigned kME |
|---|---:|---|---:|---|---:|---|---:|
| RXRA | 20182 | — | — | — | — | — | — |
| PPARA | 19522 | — | — | — | — | — | — |
| PPARD | 19269 | — | — | — | — | — | — |
| PPARG | 15865 | — | — | — | — | — | — |
| PTGS2 | 4741 | turquoise | 0.435 | brown | 0.498 | turquoise | 0.438 |
| NR1I2 | 14953 | — | — | grey | — | — | — |
| CYP2C9 | 3204 | grey | — | grey | — | yellow | 0.378 |
| CYP3A4 | 6741 | grey | — | grey | — | grey | — |

## Audit 2 — 41-gene pairwise co-expression

- Samples: 380 combined COAD/READ primary tumors.
- Pairwise correlations: 820 unique gene pairs.
- Pairs with |r| >= 0.70: 5.
- Pairs with |r| >= 0.80: 1.

The square correlation matrix and complete pair table are provided as audit outputs. This analysis is descriptive and does not establish chemical causality.

## Audit 3 — WGCNA parameter sensitivity

| Analysis | Network | Correlation | Genes | Targets in network | Non-grey modules | Best non-grey overlap | Best raw P | Best module BH-FDR |
|---|---|---|---:|---:|---:|---:|---:|---:|
| signed_bicor_8000 | signed | bicor | 8000 | 16 | 3 | turquoise (9/16) | 0.0642 | 0.257 |
| signed_bicor_15000 | signed | bicor | 15000 | 26 | 3 | brown (7/26) | 0.038 | 0.152 |
| unsigned_pearson_8000 | unsigned | pearson | 8000 | 16 | 4 | yellow (2/16) | 0.107 | 0.29 |

## Interpretation

The 15,000-gene network increased target inclusion, but retained only three non-grey modules and did not produce module-level BH-FDR<0.05 enrichment. The unsigned-Pearson network produced four non-grey modules, but likewise had no module-level BH-FDR<0.05 enrichment. Therefore the canonical negative module-enrichment result is not explained solely by the 8,000-gene filter or by signed-bicor choice.

The audit does not prove absence of co-expression biology: the 41-gene pairwise matrix and projection results show descriptive relationships, and bulk modules can reflect tissue composition. It supports reporting WGCNA as non-confirmatory for a unified DINP–CRC module under these cohorts and parameterizations.
