# DINP–CRC TCGA/GEO transcriptomic validation

This module evaluates tumor-versus-normal expression for the frozen 97-gene `DINP_CRC_overlap.csv` list. It is a CRC disease-state expression validation layer, not evidence that DINP causes the observed expression changes.

## Design

- TCGA-COAD: 288 primary tumor and 41 solid-tissue normal samples from the UCSC Toil/Xena `TcgaTargetGtex_rsem_gene_tpm` matrix.
- TCGA matched sensitivity: one primary-tumor and one solid-normal sample per matched patient where available; two-sided paired Wilcoxon tests.
- GTEx sensitivity: 308 transverse/sigmoid colon normal samples compared with TCGA-COAD primary tumor.
- GEO primary validation: GSE74602 (30 pairs), GSE10950 (24 pairs), and GSE156355 (6 pairs).
- GEO probe values were aggregated to HGNC symbols by the within-sample median across probes mapping to the same target.
- All-sample comparisons use two-sided Mann–Whitney U tests; paired comparisons use two-sided Wilcoxon signed-rank tests.
- BH FDR is calculated separately within each contrast across the 97 targets.

## Current result summary

- TCGA-COAD: 77 of 95 measured targets have FDR <0.05.
- GSE10950: 66 of 92 measured targets have FDR <0.05.
- GSE74602: 64 of 92 measured targets have FDR <0.05.
- GSE156355: 0 of 94 measured targets have FDR <0.05; this small six-pair cohort is retained as a negative/sensitivity result.
- Across the four independent primary datasets, 69 targets are significant in at least two datasets and 46 in at least three datasets.

The cross-dataset summary is sign-consistency based; expression scales are not pooled into a causal meta-analysis. Normal-reference choice and tissue composition can affect the direction and significance of tumor-versus-normal contrasts.

## Outputs

- `outputs/DINP_CRC_overlap.csv`: frozen 97-gene input list used by this module.
- `outputs/DINP_CRC_TCGA_GEO_expression_validation.csv`: gene-level statistics for all contrasts.
- `outputs/DINP_CRC_TCGA_GEO_cross_dataset_summary.csv`: independent-dataset significance and direction summary.
- `outputs/DINP_CRC_TCGA_GEO_sample_manifest.csv`: sample groups and pair identifiers.
- `outputs/DINP_CRC_TCGA_GEO_target_expression_long.csv`: target-gene expression values retained for audit.
- `outputs/DINP_CRC_TCGA_GEO_expression_validation_report.md`: human-readable report.
- `outputs/DINP_CRC_TCGA_GEO_expression_validation_log.md`: source, filtering and SHA-256 audit log.
- `outputs/DINP_CRC_TCGA_GEO_expression_validation_manifest.json`: machine-readable run manifest.

Raw GEO family SOFT files are not committed. To rerun the script, place the downloaded `GSE74602_family.soft.gz`, `GSE10950_family.soft.gz`, and `GSE156355_family.soft.gz` files under `work/geo_validation/`, then run:

```powershell
python analysis/dinp_crc_tcga_geo_expression_validation/run_tcga_geo_expression_validation.py
```

The script queries the public UCSC Toil/Xena hub for the TCGA/GTEx target-gene values.
