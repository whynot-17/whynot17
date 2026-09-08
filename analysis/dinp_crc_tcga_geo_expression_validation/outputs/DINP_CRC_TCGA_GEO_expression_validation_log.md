# DINP–CRC TCGA/GEO expression validation audit log

- Run date: 2026-09-09T00:09:15.965129+08:00
- Input: `DINP_CRC_overlap.csv`; genes: 97
- No raw TCGA/GEO matrix is committed; only target-gene values, sample metadata and derived statistics are retained in outputs.
- GEO values were aggregated to HGNC symbols by the median across probes mapping to the same target symbol within each sample.
- TCGA values were queried by HGNC symbol from the UCSC Toil/Xena RSEM TPM dataset.
- FDR correction is per contrast across the target list; it is not pooled across datasets.

## Dataset manifests

[
  {
    "dataset_id": "TCGA-COAD+GTEx-colon",
    "platform": "Toil/Xena TcgaTargetGtex_rsem_gene_tpm",
    "source_url": "https://toil.xenahubs.net",
    "expression_dataset": "TcgaTargetGtex_rsem_gene_tpm",
    "phenotype_dataset": "TcgaTargetGTEX_phenotype.txt",
    "tcga_tumor_count": 288,
    "tcga_solid_normal_count": 41,
    "gtex_colon_normal_count": 308,
    "target_genes_with_any_value": 95
  },
  {
    "dataset_id": "GSE74602",
    "platform": "GPL6104",
    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74602",
    "description": "30 paired colorectal tumor and normal samples; Illumina humanRef-8 v2.0",
    "sample_count": 60,
    "tumor_count": 30,
    "normal_count": 30,
    "paired_count": 30,
    "target_genes_with_any_value": 92,
    "target_genes_with_complete_group": 92
  },
  {
    "dataset_id": "GSE10950",
    "platform": "GPL6104",
    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10950",
    "description": "24 paired colon tumor and normal mucosa samples; Illumina humanRef-8 v2.0",
    "sample_count": 48,
    "tumor_count": 24,
    "normal_count": 24,
    "paired_count": 24,
    "target_genes_with_any_value": 92,
    "target_genes_with_complete_group": 92
  },
  {
    "dataset_id": "GSE156355",
    "platform": "GPL21185",
    "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156355",
    "description": "6 paired colorectal tumor and adjacent normal tissues; Agilent human microarray",
    "sample_count": 12,
    "tumor_count": 6,
    "normal_count": 6,
    "paired_count": 6,
    "target_genes_with_any_value": 94,
    "target_genes_with_complete_group": 94
  }
]

## Output integrity

- `DINP_CRC_TCGA_GEO_expression_validation.csv` SHA-256: `0f637842ea58d552733a08ef83de94166a85a116f5d8f6bf00a0fa1eb5424478`
- `DINP_CRC_TCGA_GEO_cross_dataset_summary.csv` SHA-256: `ee8a2207a1feec5545771bab4e13cc302364c5b05662506969b0fa28831d5a74`
- `DINP_CRC_TCGA_GEO_sample_manifest.csv` SHA-256: `b6e6b0eccef4e7298056280afea620a95182319bbe83dacce3c10ec6943bc364`
- `DINP_CRC_TCGA_GEO_target_expression_long.csv` SHA-256: `450445f2bc56bb2002d26fce7cf88afbcb15e50776c3f8900706a1aa7fda7935`
- `DINP_CRC_TCGA_GEO_expression_validation_manifest.json` SHA-256: `407d5bcb3aa105e12cabf4e177f1e9a4ed77e278b5d9c8ca9ed2139e1dc74520`
