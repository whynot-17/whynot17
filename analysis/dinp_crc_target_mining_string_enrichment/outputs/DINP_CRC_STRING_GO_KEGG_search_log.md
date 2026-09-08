# DINP–CRC STRING PPI and GO/KEGG enrichment audit log

- Generated locally: `2026-09-08T20:43:53.546163+08:00`
- Input: `DINP_CRC_overlap.csv`; SHA-256 `205832dc85eea7b6b6837a4b8cb3c1d632cd3761ece2c1972ee08cdcf3d1fbe7`
- Input genes: **97**; species: Homo sapiens; NCBI taxon ID: `9606`
- Scope: STRING identifier mapping, thresholded functional PPI network, STRING PPI enrichment, and STRING GO/KEGG functional enrichment only.
- No new nodes were added to the PPI network (`add_nodes=0`).

## Parameters

- STRING API base: `https://version-12-0.string-db.org/api`
- STRING version response: `[{"string_version": "12.0", "stable_address": "https://version-12-0.string-db.org"}]`
- Network type: `functional` (STRING functional associations; all evidence channels contribute to the combined score).
- Network minimum combined score: `700/1000` (high-confidence threshold; equivalent to >0.7 in the prior DEP–CRC workflow).
- Mapping was performed first; subsequent calls used the returned STRING protein IDs.
- GO/KEGG enrichment used STRING's default human background for the submitted STRING protein set. A custom assay-detectable background was not available for this overlap list and was not invented.

## Retrieval and count reconciliation

| Item | Count/status |
|---|---:|
| Submitted DINP–CRC genes | 97 |
| STRING mapping rows returned | 96 |
| Submitted genes mapped to STRING | 96 |
| Submitted genes not mapped by STRING v12.0 | 1 (GPX1) |
| Network edges at score >= 700 | 385 |
| Network nodes with >=1 edge | 78 |
| GO/KEGG rows returned after category filter | 1463 |
| GO/KEGG rows with STRING FDR <0.05 | 1462 |
| GO-BP / GO-MF / GO-CC / KEGG rows | 1228 / 65 / 41 / 129 |
| Significant GO-BP / GO-MF / GO-CC / KEGG rows | 1227 / 65 / 41 / 129 |

The PPI enrichment statistics are in `DINP_CRC_STRING_ppi_enrichment.csv`; the complete GO/KEGG table retains STRING p-values, FDRs, mapped-gene counts, background-gene counts, and input gene lists.

GPX1 was retained in the input and mapping audit but was not accepted by STRING v12.0's public identifier API (empty mapping response). PPI/enrichment API calls therefore use the 96 successfully mapped STRING proteins; this is a STRING coverage limitation, not an exclusion from the upstream DINP–CRC overlap file.

## PPI enrichment response

- `number_of_nodes`: `96`
- `number_of_edges`: `385`
- `average_node_degree`: `8.02`
- `local_clustering_coefficient`: `0.501`
- `expected_number_of_edges`: `93`
- `p_value`: `0.0`

## Files

- `DINP_CRC_STRING_mapping.csv`: every mapping row returned by STRING.
- `DINP_CRC_STRING_mapping_selected.csv`: one selected STRING ID for each of the 97 input genes.
- `DINP_CRC_STRING_network_edges.csv`: high-confidence functional STRING edges among input proteins only.
- `DINP_CRC_STRING_node_degree.csv`: degree for all 97 input nodes, including zero-degree nodes.
- `DINP_CRC_STRING_ppi_enrichment.csv`: STRING network-vs-background interaction enrichment statistics.
- `DINP_CRC_GO_KEGG_enrichment.csv`: combined GO and KEGG enrichment output.
- `DINP_CRC_GO_enrichment.csv` and `DINP_CRC_KEGG_enrichment.csv`: category-specific copies.
- Raw JSON responses are retained under `work/` for auditability.

## Interpretation boundary

STRING edges are functional associations unless a physical network is explicitly requested; they should not be described as direct biochemical binding. Enrichment is over-representation reported by STRING, and statistical significance is not evidence of DINP causality. Review source/evidence partitions in `DINP_CRC_overlap.csv` before using individual genes for MR or mechanistic claims.
