# DINP–CRC 101ML × Tier 1 overlap audit log

```json
{
  "background_definition": "97 unique DINP-CRC overlap gene symbols from DINP_CRC_overlap.csv",
  "set_definitions": {
    "stable_ml": "stable_ml_flag == True in DINP_CRC_101ML_gene_stability.csv",
    "tier1": "cross_support_flag == True in DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv",
    "ppi_pathway_candidate": "ppi_pathway_candidate == True",
    "transcriptomic_replicated_ge2": "independent_datasets_significant_fdr_lt_0_05 >= 2",
    "crossrank_top20": "cross_rank <= 20"
  },
  "stable_ml_n": 17,
  "tier1_n": 16,
  "primary_overlap_n": 2,
  "primary_overlap_genes": [
    "CEBPB",
    "IL1A"
  ],
  "fisher_table": [
    [
      2,
      15
    ],
    [
      14,
      66
    ]
  ],
  "multiple_testing": "Benjamini-Hochberg adjustment across the 10 two-sided pairwise Fisher tests",
  "qc": {
    "background_n": 97,
    "membership_rows": 97,
    "membership_gene_symbols_unique": true,
    "pairwise_comparisons": 10
  }
}
```
