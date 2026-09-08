# DINP–CRC 101ML × Tier 1 overlap and enrichment audit

## Primary comparison

The primary test compares the 17 default stable-ML genes with the 16 pre-existing Tier 1 cross-supported genes inside the frozen 97-gene DINP–CRC overlap universe.

| Quantity | Value |
|---|---:|
| Background universe | 97 |
| Stable-ML genes | 17 |
| Tier 1 genes | 16 |
| Observed overlap | 1 |
| Expected overlap under random overlap | 2.804 |
| Observed / expected | 0.357 |
| Odds ratio | 0.2708 |
| Fisher exact p, two-sided | 0.290365 |
| Fisher exact p, overlap greater | 0.966008 |
| Fisher exact p, overlap less | 0.176237 |
| BH q, two-sided pairwise audit | 0.362956 |
| Jaccard index | 0.0312 |
| Overlap genes | CEBPB |

## 2×2 table

Rows are stable-ML status; columns are Tier 1 status.

| | Tier 1 yes | Tier 1 no |
|---|---:|---:|
| Stable ML yes | 1 | 16 |
| Stable ML no | 15 | 65 |

## Pairwise overlap audit

The full pairwise table tests the stable-ML and Tier 1 sets against PPI/pathway candidates, transcriptomic replication ≥2/4, and the top-20 cross-ranked set, always using the same 97-gene background. BH q-values adjust the two-sided Fisher tests across all pairwise comparisons.

| Set A | Set B | A | B | Overlap | Expected | Fold enrichment | Odds ratio | Fisher p | BH q |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 17 stable-ML genes | 16 Tier 1 cross-supported genes | 17 | 16 | 1 | 2.804 | 0.357 | 0.271 | 0.2904 | 0.363 |
| 16 Tier 1 cross-supported genes | PPI top-20 hub + pathway membership | 16 | 20 | 16 | 3.299 | 4.850 | inf | 6.109e-15 | 6.109e-14 |
| 16 Tier 1 cross-supported genes | Top 20 cross-ranked genes | 16 | 20 | 14 | 3.299 | 4.244 | 87.500 | 1.445e-10 | 7.226e-10 |
| PPI top-20 hub + pathway membership | Top 20 cross-ranked genes | 20 | 20 | 14 | 4.124 | 3.395 | 27.611 | 3.489e-08 | 1.163e-07 |
| Transcriptomic replicated in ≥2/4 primary datasets | Top 20 cross-ranked genes | 69 | 20 | 20 | 14.227 | 1.406 | inf | 0.0005922 | 0.001481 |
| 17 stable-ML genes | Transcriptomic replicated in ≥2/4 primary datasets | 17 | 69 | 17 | 12.093 | 1.406 | inf | 0.002359 | 0.004719 |
| 16 Tier 1 cross-supported genes | Transcriptomic replicated in ≥2/4 primary datasets | 16 | 69 | 16 | 11.381 | 1.406 | inf | 0.004736 | 0.007894 |
| 17 stable-ML genes | PPI top-20 hub + pathway membership | 17 | 20 | 1 | 3.505 | 0.285 | 0.201 | 0.1828 | 0.2612 |
| PPI top-20 hub + pathway membership | Transcriptomic replicated in ≥2/4 primary datasets | 20 | 69 | 16 | 14.227 | 1.125 | 1.811 | 0.4135 | 0.4595 |
| 17 stable-ML genes | Top 20 cross-ranked genes | 17 | 20 | 3 | 3.505 | 0.856 | 0.794 | 1 | 1 |

## Gene-level membership

The 97-row membership table preserves stable-ML status, Tier 1 status, PPI/pathway candidate status, transcriptomic replication status, and both ranking positions for every background gene.

## Interpretation

The observed stable-ML × Tier 1 overlap is interpreted against the 97-gene background only. Fisher's exact test evaluates set overlap, not biological causality, model performance, or MR validity. The two sets were constructed independently, so limited overlap is informative rather than a failure by itself.
