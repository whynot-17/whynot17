# DINP–CRC PPI hub and pathway membership audit log

- Generated locally: 2026-09-08T21:08:12.745686+08:00
- Input overlap: `DINP_CRC_overlap.csv` (97 genes).
- PPI network: `DINP_CRC_STRING_network_edges.csv`, STRING functional edges at combined score >=700.
- Hub ranking is continuous and deterministic; it is sorted by degree, weighted degree, betweenness, eigenvector centrality, then gene symbol.
- `top20` is a reporting convenience only; it is not a statistical cutoff.
- Pathway membership uses significant (`FDR < 0.05`) nonredundant STRING representative terms from `DINP_CRC_GO_KEGG_representatives.csv`.
- Membership is reported as gene-to-term evidence; it does not imply DINP-specific causal regulation.

## Count reconciliation

- Input genes represented in hub table: 97
- Mapped STRING proteins ranked: 96
- Input genes not mapped to STRING: 1
- Nodes with >=1 high-confidence edge: 78
- Network edges used: 385
- Significant representative terms: 110
- Gene–term membership rows: 1209
- Top-20 hub–term membership rows: 449

## Output files

- `DINP_CRC_STRING_hub_ranking.csv`: all mapped STRING proteins with degree and centrality metrics.
- `DINP_CRC_STRING_hubs_top20.csv`: top 20 ranked hubs for focused review.
- `DINP_CRC_pathway_membership.csv`: all input-gene memberships in significant representative GO/KEGG terms.
- `DINP_CRC_hub_pathway_membership.csv`: the same membership table restricted to top-20 hubs.
- `DINP_CRC_hub_pathway_summary.csv`: one row per top-20 hub with term counts and member terms.
