# DINP target mining and DINP–CRC PPI/pathway analysis

Frozen analysis outputs generated on 2026-09-08 from the DINP target-mining and CRC-overlap workflow.

## Scope

- DINP target master table with database and literature evidence partitions.
- CRC disease-gene union using GeneCards relevance score ≥10, OMIM and TTD.
- DINP–CRC overlap: 97 genes.
- STRING v12.0 functional network at combined score ≥700: 96 of 97 genes mapped, 385 edges, and 78 nodes with at least one edge.
- STRING GO/KEGG enrichment: 1,463 returned terms; 1,462 with STRING FDR <0.05.
- Redundancy reduction and deterministic ranking: 110 significant representative terms.
- PPI hub ranking by degree, weighted degree, betweenness, eigenvector centrality and PageRank.
- Gene-to-term pathway membership for the significant nonredundant representative terms.

## Important interpretation limits

STRING edges are functional associations and should not be described as direct biochemical binding. Pathway membership is an annotation/enrichment result and does not establish DINP-specific causal regulation. GPX1 is retained in the hub table as an input gene but is explicitly marked as unmapped by STRING and is not assigned a hub rank.

Raw downloaded data and intermediate API payloads are intentionally not included. Source URLs, filters, evidence grading and count reconciliation are recorded in the audit logs under `outputs/`.

## Key outputs

- `DINP_target_master.csv` and `DINP_target_evidence_summary.csv`
- `CRC_gene_set_GeneCards10_OMIM_TTD.csv` and `DINP_CRC_overlap.csv`
- `DINP_CRC_STRING_mapping*.csv`, `DINP_CRC_STRING_network_edges.csv`, `DINP_CRC_STRING_node_degree.csv`, and `DINP_CRC_STRING_ppi_enrichment.csv`
- `DINP_CRC_GO_KEGG_enrichment.csv` plus the GO/KEGG category tables
- `DINP_CRC_GO_KEGG_reduced_ranked.csv` and representative-term tables
- `DINP_CRC_STRING_hub_ranking.csv` and `DINP_CRC_STRING_hubs_top20.csv`
- `DINP_CRC_pathway_membership.csv`, `DINP_CRC_hub_pathway_membership.csv`, and `DINP_CRC_hub_pathway_summary.csv`

The hub/membership script is `run_ppi_hub_pathway_membership.py`. It consumes the versioned CSV outputs in this module and does not download new data.
