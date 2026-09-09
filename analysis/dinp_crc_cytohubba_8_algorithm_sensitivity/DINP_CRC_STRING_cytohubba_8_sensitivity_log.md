# DINP–CRC CytoHubba 8-algorithm sensitivity audit log

```json
{
  "network_input": "DINP_CRC_STRING_network_edges.csv filtered to score >= 0.700",
  "network_edge_count": 385,
  "network_node_count": 78,
  "largest_connected_component": 76,
  "input_background_n": 97,
  "algorithm_panel": [
    "MCC",
    "MNC",
    "EPC",
    "Degree",
    "Closeness",
    "Betweenness",
    "Radiality",
    "Stress"
  ],
  "top_k": 25,
  "consensus_min_algorithms": 7,
  "epc": {
    "runs": 1000,
    "edge_removal_threshold": 0.5,
    "seed": 20260909
  },
  "unresolved_edges_skipped": [],
  "consensus_hub_n": 18,
  "original_degree_top20_n": 20,
  "tier1_n": 16,
  "stable_ml_n": 17,
  "cebpb_consensus_hub": true,
  "cd36_consensus_hub": true,
  "tier1_consensus_overlap": [
    "ADIPOQ",
    "CEBPB",
    "ESR1",
    "IFNG",
    "IL1B",
    "IL6",
    "MMP9",
    "PPARA",
    "PPARG",
    "RELA",
    "STAT3",
    "TP53"
  ],
  "stable_ml_consensus_overlap": [
    "CD36",
    "CEBPB"
  ],
  "qc": {
    "ranking_rows": 97,
    "unique_gene_symbols": 97,
    "long_ranking_rows": 624,
    "long_rows_expected": 624,
    "one_top25_set_per_algorithm": {
      "MCC": 25,
      "MNC": 25,
      "EPC": 25,
      "Degree": 25,
      "Closeness": 25,
      "Betweenness": 25,
      "Radiality": 25,
      "Stress": 25
    }
  }
}
```
