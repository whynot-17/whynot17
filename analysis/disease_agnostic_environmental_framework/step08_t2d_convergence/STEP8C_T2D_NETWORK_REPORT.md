# Step 8C — STRING network convergence

- Status: **complete_network_convergence**
- STRING species: Homo sapiens (9606)
- Network: functional associations, combined score >= 700, no added interactors
- Input: frozen Step 7 overlap genes per Tier A axis
- Randomization background: frozen union of all 11 Step 7 cluster genes
- Permutations: 1,000 degree-stratified samples per axis
- Network-priority score: normalized degree, betweenness, eigenvector, within-module connectivity, and pathway recurrence; descriptive ranking only

## Axis summary

| Axis | Input genes | Mapped IDs | Network nodes | Edges | Components | Louvain modules | Observed/expected | Empirical P |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| cluster_11 | 245 | 240 | 138 | 190 | 13 | 19 | 1.515 | 0.000999 |
| cluster_5 | 65 | 62 | 51 | 193 | 2 | 4 | 4.567 | 0.000999 |
| cluster_6 | 919 | 801 | 549 | 2536 | 13 | 24 | 1.612 | 0.000999 |
| cluster_8 | 3016 | 2917 | 2224 | 13364 | 36 | 50 | 2.023 | 0.000999 |

## Interpretation boundary

Network-prioritized genes are not causal targets.  The network-priority score is a descriptive ranking aid that averages normalized degree, betweenness, eigenvector centrality, within-module connectivity, and pathway recurrence; module membership and annotations remain separate audit fields.  STRING functional associations may include evidence beyond direct physical binding, so results are described as high-confidence functional network convergence rather than definitive direct PPI.

All axes were analyzed separately.  No transcriptomic data, T2D expression data, or flagship selection entered this stage.
