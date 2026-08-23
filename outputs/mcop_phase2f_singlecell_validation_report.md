# MCOP–CRC Phase 2F-B：CELLxGENE Census 单细胞验证

## 当前判定：**PARTIAL**

This is a targeted partial run: the paired Census dataset was processed across the available compartments (endothelial, epithelial, fibroblast, myeloid). It cannot trigger the frozen multi-dataset GREEN/YELLOW/RED mechanism gate; the independent dataset replication remains a separate analysis.

The targeted dataset result is informative but cannot trigger the frozen multi-dataset mechanism gate.

本轮的疾病细胞标签只能支持 **tumor-derived epithelial**；本脚本没有把它写成 malignant epithelial，也没有做 CNV 推断或使用未经核验的 malignant 标签。

## 冻结规则

- Census release: `2025-11-08`；organism: `homo_sapiens`。
- 所有 metadata/expression 查询均包含 `is_primary_data == True`。
- 统计单位是 donor-level pseudobulk；没有把 cell 当作独立样本做 P 值。
- 主 score：`PPAR_nuclear_receptor_score`；基因：`PPARA, PPARD, PPARG, NR1I2, NR1I3, NR1H2, NR1H3`。
- Secondary：`RELA_STAT3_score`、`DINP_axis_9_gene_score`。
- Myeloid、fibroblast、endothelial 仅用于 compartment localization。

## Census 标签审计

- Tumor disease labels: `colon adenocarcinoma; colorectal carcinoma || metastatic malignant neoplasm`
- Normal disease labels: `normal`
- Tissue-general candidates: `large intestine; colon`
- Relevant observations after primary/tissue gate: **1,145,214**
- Datasets: **23**；donors with usable IDs: **289**；cell types: **184**

| group | cells | datasets | donor IDs |
|---|---:|---:|---:|
| normal | 856,427 | 22 | 257 |
| tumor | 288,787 | 2 | 68 |

## Primary donor-level result

The pooled result below is a donor-level descriptive comparison. Dataset-level effects and leave-one-dataset-out results are the required stability checks; pooled P values are not cell-level P values.

| compartment | score | tumor donors | normal donors | median delta | P |
|---|---|---:|---:|---:|---:|
| epithelial | PPAR_nuclear_receptor_score | 62 | 36 | -0.385 | 7.4e-12 |
| epithelial | RELA_STAT3_score | 62 | 36 | 1.379 | 5.87e-13 |
| epithelial | DINP_axis_9_gene_score | 62 | 36 | 0.021 | 0.909 |

## Dataset direction check

A dataset is eligible for this check only when it has at least two tumor and two normal donors in the same compartment.

| compartment | score | eligible datasets | positive | negative | positive fraction |
|---|---|---:|---:|---:|---:|
| endothelial | DINP_axis_9_gene_score | 1 | 1 | 0 | 1.00 |
| endothelial | PPAR_nuclear_receptor_score | 1 | 1 | 0 | 1.00 |
| endothelial | RELA_STAT3_score | 1 | 1 | 0 | 1.00 |
| epithelial | DINP_axis_9_gene_score | 1 | 0 | 1 | 0.00 |
| epithelial | PPAR_nuclear_receptor_score | 1 | 0 | 1 | 0.00 |
| epithelial | RELA_STAT3_score | 1 | 1 | 0 | 1.00 |
| fibroblast | DINP_axis_9_gene_score | 1 | 0 | 1 | 0.00 |
| fibroblast | PPAR_nuclear_receptor_score | 1 | 0 | 1 | 0.00 |
| fibroblast | RELA_STAT3_score | 1 | 1 | 0 | 1.00 |
| myeloid | DINP_axis_9_gene_score | 1 | 1 | 0 | 1.00 |
| myeloid | PPAR_nuclear_receptor_score | 1 | 1 | 0 | 1.00 |
| myeloid | RELA_STAT3_score | 1 | 1 | 0 | 1.00 |

## Paired-donor check

When the same donor ID contains both tumor and normal observations, the paired donor result is reported separately. This is not substituted for the multi-dataset gate.

| dataset | compartment | score | paired donors | median delta | P |
|---|---|---|---:|---:|---:|
| 16023185-de21-4c0d-a9c8-73abdd52d142 | endothelial | PPAR_nuclear_receptor_score | 33 | -0.019 | 0.609 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | endothelial | RELA_STAT3_score | 33 | 0.061 | 0.000952 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | endothelial | DINP_axis_9_gene_score | 33 | -0.000 | 0.358 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | epithelial | PPAR_nuclear_receptor_score | 36 | -0.419 | 4.29e-07 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | epithelial | RELA_STAT3_score | 36 | 1.167 | 1.08e-07 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | epithelial | DINP_axis_9_gene_score | 36 | 0.011 | 0.636 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | fibroblast | PPAR_nuclear_receptor_score | 31 | 0.172 | 0.504 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | fibroblast | RELA_STAT3_score | 31 | 0.051 | 0.347 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | fibroblast | DINP_axis_9_gene_score | 31 | 0.184 | 0.433 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | myeloid | PPAR_nuclear_receptor_score | 35 | 0.610 | 7.97e-09 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | myeloid | RELA_STAT3_score | 35 | -0.029 | 0.728 |
| 16023185-de21-4c0d-a9c8-73abdd52d142 | myeloid | DINP_axis_9_gene_score | 35 | 0.471 | 1.75e-10 |

## Leave-one-dataset-out

`mcop_phase2f_singlecell_leave_one_dataset_out.csv` contains the full table. For the primary epithelial PPAR/NR score, the report uses only leave-one-out rows with at least two donors per group; this prevents a formally computed but uninformative result from being called stable.

## Interpretation boundaries

- A stable epithelial result would support a CRC-associated epithelial program, not exposure causality or DINP-to-CRC mediation.
- A stable myeloid/fibroblast/endothelial result with unstable epithelial direction would redirect the mechanism toward microenvironmental PPAR/NR remodeling.
- Tissue composition, dataset-specific annotation, treatment history and dissociation effects remain possible explanations. The analysis is not a substitute for prediagnostic tissue or perturbation validation.
- Targeted normalization uses the nine frozen genes because the standardized Census observation metadata do not provide a universal total-UMI field. This is explicitly a targeted score and should not be overinterpreted as full transcriptome normalization.

## Output files

- `mcop_phase2f_singlecell_label_discovery.json`
- `mcop_phase2f_singlecell_dataset_donor_cell_audit.csv`
- `mcop_phase2f_singlecell_donor_audit.csv`
- `mcop_phase2f_singlecell_donor_pseudobulk.csv`
- `mcop_phase2f_singlecell_donor_scores.csv`
- `mcop_phase2f_singlecell_dataset_effects.csv`
- `mcop_phase2f_singlecell_pooled_contrasts.csv`
- `mcop_phase2f_singlecell_leave_one_dataset_out.csv`
- `mcop_phase2f_singlecell_direction_summary.csv`
- `mcop_phase2f_singlecell_paired_donor_contrasts.csv`

## Reproducibility

- Run UTC: `2026-08-23T01:50:37.063914+00:00`
- Python: `3.12.3`
- Script: `work/scripts/mcop_phase2f_singlecell_validation.py`
- NR genes: `PPARA,PPARD,PPARG,NR1I2,NR1I3,NR1H2,NR1H3`
- All genes: `PPARA,PPARD,PPARG,NR1I2,NR1I3,NR1H2,NR1H3,RELA,STAT3`

Primary epithelial PPAR/NR eligible datasets=1; positive=0; leave-one-dataset-out stable=False.

## Scope note

This is a targeted/partial run (explicit dataset, compartment, or pseudobulk reuse restriction). Do not apply the frozen GREEN/YELLOW/RED mechanism gate until the full multi-dataset, four-compartment run is complete.
