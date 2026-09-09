# DINP–CRC CytoHubba 8-algorithm hub sensitivity

This directory contains a prespecified sensitivity analysis of the frozen 97-gene DINP–CRC overlap set. It uses the same STRING PPI edge table as the primary analysis, filters edges at combined score ≥0.700, and applies eight fixed CytoHubba-style network-ranking algorithms:

`MCC`, `MNC`, `EPC`, `Degree`, `Closeness`, `Betweenness`, `Radiality`, and `Stress`.

For each algorithm, the top 25 genes are retained. A consensus hub is defined as a gene appearing in the top 25 for at least 7 of 8 algorithms. This threshold and cutoff were fixed before inspecting the resulting gene list.

## Important reproducibility note

The cited CRC study reports the rule “top 25 in at least 7 of 8 CytoHubba algorithms,” but the accessible text does not enumerate those eight algorithms and another methods sentence refers to all 12 CytoHubba methods. Therefore this directory transparently records an operational canonical eight-method panel; it should not be described as an exact reconstruction of an undocumented panel. See the [CRC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7463304/) and the original [CytoHubba paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4290687/).

The implementation uses unweighted graph topology after the STRING score filter. EPC uses 1,000 reproducible Monte Carlo edge-percolation runs with removal threshold 0.5 and seed `20260909`. All parameters and QC checks are recorded in `DINP_CRC_STRING_cytohubba_8_sensitivity_log.md`.

## Results

- 97 input overlap genes; 78 mapped/in-network nodes; 385 STRING edges; largest component 76 nodes.
- 18 consensus hubs meet the ≥7/8 rule.
- CEBPB: consensus hub, 7/8 algorithms; original degree Top20 and Tier 1; stable-ML gene.
- CD36: consensus hub, 7/8 algorithms; stable-ML gene; not original degree Top20 or Tier 1.
- 12/16 original Tier 1 genes are retained by the consensus rule.
- 2/17 stable-ML genes are consensus hubs: CEBPB and CD36.

## Files

- `DINP_CRC_STRING_cytohubba_8_consensus_hubs.csv`: consensus hubs and annotations.
- `DINP_CRC_STRING_cytohubba_8_algorithm_rankings.csv`: all 97 genes with ranks/scores for every algorithm.
- `DINP_CRC_STRING_cytohubba_8_top25_by_algorithm.csv`: long-format Top25 membership for each algorithm.
- `DINP_CRC_STRING_cytohubba_8_sensitivity_comparison.csv`: CEBPB, CD36, Tier 1, and stable-ML comparison table.
- `DINP_CRC_STRING_cytohubba_8_overlap_audit.csv`: set sizes and explicit gene intersections.
- `DINP_CRC_STRING_cytohubba_8_sensitivity_summary.md`: human-readable summary.
- `DINP_CRC_STRING_cytohubba_8_sensitivity_log.md`: machine-readable audit log and QC.
- `run_8_algorithm_hub_sensitivity.py`: reproducible analysis script.

The three input tables copied here document the exact frozen background and STRING network used by the script. This analysis does not overwrite the original degree-based Top20, Tier 1, or stable-ML definitions.
