# DINP–CRC cross-ranking audit log

```json
{
  "inputs": {
    "overlap": "C:\\Users\\21634\\Documents\\Codex\\crossrank-push\\analysis\\dinp_crc_ppi_pathway_transcriptomic_cross_rank\\DINP_CRC_overlap.csv",
    "ppi_hub_ranking": "C:\\Users\\21634\\Documents\\Codex\\crossrank-push\\analysis\\dinp_crc_ppi_pathway_transcriptomic_cross_rank\\DINP_CRC_STRING_hub_ranking.csv",
    "pathway_membership": "C:\\Users\\21634\\Documents\\Codex\\crossrank-push\\analysis\\dinp_crc_ppi_pathway_transcriptomic_cross_rank\\DINP_CRC_pathway_membership.csv",
    "transcriptomic_summary": "C:\\Users\\21634\\Documents\\Codex\\crossrank-push\\analysis\\dinp_crc_ppi_pathway_transcriptomic_cross_rank\\DINP_CRC_TCGA_GEO_cross_dataset_summary.csv"
  },
  "input_gene_count": 97,
  "pathway_membership_rows_after_term_deduplication": 1209,
  "pathway_membership_gene_count": 93,
  "independent_primary_dataset_count": 4,
  "independent_primary_datasets": [
    "TCGA-COAD",
    "GSE74602",
    "GSE10950",
    "GSE156355"
  ],
  "sensitivity_contrasts_excluded_from_primary_score": [
    "TCGA-COAD_paired",
    "TCGA-COAD_vs_GTEx-colon"
  ],
  "score_definition": {
    "ppi_hub_score": "descending percentile among mapped and in-network STRING genes",
    "pathway_component_score": "descending percentile of log1p(representative_term_count) across the 97 genes; zero if no representative term",
    "transcriptomic_component_score": "(independent primary datasets significant at FDR<0.05 / 4) * direction concordance fraction",
    "cross_rank_score": "cube root of product of the three component scores"
  },
  "qc": {
    "duplicate_gene_symbols_in_output": 0,
    "cross_rank_unique": 97,
    "rank_score_range": {
      "ppi": [
        0.0,
        1.0
      ],
      "pathway": [
        0.0,
        1.0
      ],
      "transcriptomic": [
        0.0,
        0.75
      ],
      "cross": [
        0.0,
        0.8737094376549739
      ]
    }
  }
}
```
