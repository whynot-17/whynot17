# DINP–CRC 101ML × Tier 1 overlap and enrichment audit

## Primary comparison

The primary test compares the 17 default stable-ML genes with the 16 pre-existing Tier 1 cross-supported genes inside the frozen 97-gene DINP–CRC overlap universe.

| Quantity | Value |
|---|---:|
| Background universe | 97 |
| Stable-ML genes | 17 |
| Tier 1 genes | 16 |
| Observed overlap | 2 |
| Expected overlap under random overlap | 2.804 |
| Observed / expected | 0.713 |
| Odds ratio | 0.6286 |
| Fisher exact p, two-sided | 0.729793 |
| Fisher exact p, overlap greater | 0.823763 |
| Fisher exact p, overlap less | 0.434864 |
| BH q, two-sided pairwise audit | 0.729793 |
| Jaccard index | 0.0645 |
| Overlap genes | CEBPB;IL1A |

## 2×2 table

Rows are stable-ML status; columns are Tier 1 status.

| | Tier 1 yes | Tier 1 no |
|---|---:|---:|
| Stable ML yes | 2 | 15 |
| Stable ML no | 14 | 66 |

## Pairwise overlap audit

The full pairwise table tests the stable-ML and Tier 1 sets against PPI/pathway candidates, transcriptomic replication ≥2/3, and the top-20 cross-ranked set, always using the same 97-gene background. BH q-values adjust the two-sided Fisher tests across all pairwise comparisons.

| Set A | Set B | A | B | Overlap | Expected | Fold enrichment | Odds ratio | Fisher p | BH q |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 17 stable-ML genes | 16 Tier 1 cross-supported genes | 17 | 16 | 2 | 2.804 | 0.713 | 0.629 | 0.7298 | 0.7298 |
| 16 Tier 1 cross-supported genes | PPI top-20 hub + pathway membership | 16 | 20 | 16 | 3.299 | 4.850 | inf | 6.109e-15 | 6.109e-14 |
| 16 Tier 1 cross-supported genes | Top 20 cross-ranked genes | 16 | 20 | 13 | 3.299 | 3.941 | 45.810 | 7.295e-09 | 3.647e-08 |
| PPI top-20 hub + pathway membership | Top 20 cross-ranked genes | 20 | 20 | 13 | 4.124 | 3.153 | 18.571 | 7.194e-07 | 2.398e-06 |
| Transcriptomic replicated in ≥2/3 primary datasets | Top 20 cross-ranked genes | 69 | 20 | 20 | 14.227 | 1.406 | inf | 0.0005922 | 0.001481 |
| 16 Tier 1 cross-supported genes | Transcriptomic replicated in ≥2/3 primary datasets | 16 | 69 | 16 | 11.381 | 1.406 | inf | 0.004736 | 0.009472 |
| 17 stable-ML genes | Transcriptomic replicated in ≥2/3 primary datasets | 17 | 69 | 16 | 12.093 | 1.323 | 8.151 | 0.02028 | 0.0338 |
| 17 stable-ML genes | Top 20 cross-ranked genes | 17 | 20 | 5 | 3.505 | 1.426 | 1.806 | 0.3336 | 0.4766 |
| PPI top-20 hub + pathway membership | Transcriptomic replicated in ≥2/3 primary datasets | 20 | 69 | 16 | 14.227 | 1.125 | 1.811 | 0.4135 | 0.5169 |
| 17 stable-ML genes | PPI top-20 hub + pathway membership | 17 | 20 | 2 | 3.505 | 0.571 | 0.459 | 0.511 | 0.5678 |

## Gene-level membership

The 97-row membership table preserves stable-ML status, Tier 1 status, PPI/pathway candidate status, transcriptomic replication status, and both ranking positions for every background gene.

## Interpretation

The observed stable-ML × Tier 1 overlap is interpreted against the 97-gene background only. Fisher's exact test evaluates set overlap, not biological causality, model performance, or MR validity. The two sets were constructed independently, so limited overlap is informative rather than a failure by itself.
