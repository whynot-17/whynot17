# MCOP–CRC Phase 2F-B：GSE144735 independent epithelial replication

## Scope

- Source: [GSE144735 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735)
- Primary contrast: core `Tumor` versus matched `Normal` epithelial cells.
- `Border` samples are excluded from the primary contrast and are sensitivity-only.
- Inferential unit: patient-level pseudobulk; no cell-level P values.
- Cells in annotation: **27,414**; epithelial cells analyzed: **6,168**.
- Patient-level pseudobulk rows: **12**.

## Primary paired result

| score | paired patients | median tumor-minus-normal delta | Wilcoxon P |
|---|---:|---:|---:|
| PPAR_nuclear_receptor_score | 6 | -0.312 | 0.688 |
| RELA_STAT3_score | 6 | 0.372 | 0.438 |
| DINP_axis_9_gene_score | 6 | -0.112 | 1 |

## Border-included sensitivity

This sensitivity adds Border epithelial cells to the tumor group; it is not substituted for the core Tumor-versus-Normal primary analysis.

| score | paired patients | median tumor-plus-border-minus-normal delta | Wilcoxon P |
|---|---:|---:|---:|
| PPAR_nuclear_receptor_score | 6 | -0.326 | 0.688 |
| RELA_STAT3_score | 6 | 0.552 | 0.312 |
| DINP_axis_9_gene_score | 6 | -0.110 | 1 |

## Interpretation

A concordant negative epithelial PPAR/NR delta would provide independent disease-state replication of the Census/TCGA direction. It does not establish DINP/MCOP causality or exposure mediation.

## Reproducibility

- Annotation SHA256: `4eb9da239bef95b9003c0f6d0eeb41c6598ecc794cd25f57ad06d9bad81e881b`
- Matrix SHA256: `046cbbd6c7501a4ddbca6d09fc0bb923e62d4d9421d5781efdeacb2c0a5b3b31`
- Frozen NR genes: `PPARA,PPARD,PPARG,NR1I2,NR1I3,NR1H2,NR1H3`
- Run UTC: `2026-08-23T01:51:11.116916+00:00`

## Output files

- `mcop_phase2f_external_gse144735_pseudobulk.csv`
- `mcop_phase2f_external_gse144735_scores.csv`
- `mcop_phase2f_external_gse144735_paired_contrasts.csv`
- `mcop_phase2f_external_gse144735_border_sensitivity_pseudobulk.csv`
- `mcop_phase2f_external_gse144735_border_sensitivity_scores.csv`
- `mcop_phase2f_external_gse144735_report.md`
