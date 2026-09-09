# DINP–CRC CytoHubba 8-algorithm PPI sensitivity analysis

## Prespecified sensitivity rule

The original CRC report explicitly states that genes ranked in the top 25 in at least 7 of 8 CytoHubba algorithms were considered hub genes. Its text does not enumerate the eight algorithms and elsewhere refers to all 12 CytoHubba methods; therefore this analysis records an operationally fixed canonical eight-method panel rather than claiming an exact undocumented panel replication.

Panel: MCC, MNC, EPC, Degree, Closeness, Betweenness, Radiality, Stress.

A gene is a sensitivity consensus hub only when it is in the Top 25 for at least 7/8 algorithms. The threshold and Top25 cutoff were not tuned to the results.

## Network

- STRING network edge threshold: combined score ≥0.700.
- Edges: 385; mapped/in-network nodes: 78; largest connected component: 76 nodes.
- Unresolved/self-loop edges skipped: 0.
- Degree, MNC, MCC, shortest-path centralities, and Stress use the unweighted topology; STRING combined score is retained only as the network-construction filter/edge metadata.
- EPC: 1000 reproducible Monte Carlo edge-percolation runs, edge-removal threshold 0.5, seed 20260909.

## Results

- Consensus hubs (≥7/8 Top25): 18
- Original degree Top20: 20
- Consensus hubs retained from original degree Top20: 16
- Original Tier 1 genes retained as consensus hubs: 12/16
- Stable-ML genes that are consensus hubs: 2/17
- CEBPB consensus status: yes
- CD36 consensus status: yes

## Consensus hubs

| Consensus rank | Gene | Algorithms in Top25 | Count | Original degree rank | Tier 1 | Stable ML |
|---:|---|---|---:|---:|---|---|
| 1 | AKT1 | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 4 | no | no |
| 2 | ESR1 | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 14 | yes | no |
| 3 | IL1B | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 6 | yes | no |
| 4 | IL6 | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 1 | yes | no |
| 5 | MMP9 | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 8 | yes | no |
| 6 | NFE2L2 | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 24 | no | no |
| 7 | PPARG | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 7 | yes | no |
| 8 | RELA | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 12 | yes | no |
| 9 | SIRT1 | MCC;MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 8 | 10 | no | no |
| 10 | ADIPOQ | MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 7 | 18 | yes | no |
| 11 | CD36 | MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 7 | 23 | no | yes |
| 12 | CEBPB | MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 7 | 17 | yes | yes |
| 13 | IFNG | MCC;MNC;EPC;Degree;Closeness;Radiality;Stress | 7 | 9 | yes | no |
| 14 | IL4 | MCC;MNC;EPC;Degree;Closeness;Radiality;Stress | 7 | 11 | no | no |
| 15 | PPARA | MNC;EPC;Degree;Closeness;Betweenness;Radiality;Stress | 7 | 16 | yes | no |
| 16 | STAT3 | MCC;MNC;Degree;Closeness;Betweenness;Radiality;Stress | 7 | 5 | yes | no |
| 17 | TNF | MCC;MNC;Degree;Closeness;Betweenness;Radiality;Stress | 7 | 2 | no | no |
| 18 | TP53 | MCC;MNC;Degree;Closeness;Betweenness;Radiality;Stress | 7 | 3 | yes | no |

## Audit interpretation

The degree Top20 and 8-algorithm consensus are separate hub definitions. The sensitivity analysis does not overwrite the original degree-based PPI result or the original Tier 1 label.

The individual algorithm rankings and the gene-level comparison table should be reviewed before using consensus status as a downstream prioritization criterion. Consensus hub status is network-topology evidence, not direct DINP binding or causal evidence.
