# MCOP–CRC Phase 2F-A：CRC PPAR/nuclear-receptor disease-state validation

## 判定

本轮只验证 Phase 2E 冻结的 DINP-axis molecular program，不把 MCOP 当作 CTD-specific 分子发现物，也不进行泛 GO/PPI 扩张。

- PPAR/nuclear-receptor score 在 TCGA 内部对照中下降：tumor median=0.0031，normal median=0.671，delta=-0.668，P=1.24e-16，BH-FDR=3.72e-16；病人级配对分析同方向，delta=-0.533，P=7.4e-05。
- 但换成 GTEx colon normal 后 PPAR/NR score 方向变为上升：delta=0.162，P=0.000336，BH-FDR=0.000504。9 个基因中只有 5/9 个基因在两个 normal reference 下方向一致。
- 因此 Phase 2F-A 当前判定是 **reference-dependent / provisional，不升级为机制验证通过**。TCGA 内部配对结果支持疾病状态相关性，但 TCGA–GTEx 对照不支持一个稳定的统一方向。
- 样本量：TCGA primary tumor=380；TCGA CRC solid normal=51；GTEx transverse/sigmoid colon normal=308。
- 这一步可以回答“该 program 是否在 CRC tumor state 中呈现疾病相关表达状态”，但不能回答“MCOP 是否导致该 program 改变”。

## 1. 数据与对照

| Contrast | Tumor | Normal | 用途 |
|---|---:|---:|---|
| TCGA primary CRC | 380 | 51 | primary tumor vs within-TCGA solid normal |
| TCGA vs GTEx | 380 | 308 | external normal-tissue reference |

Expression source: UCSC Toil Xena hub `TcgaTargetGtex_rsem_gene_tpm`; phenotype source: `TcgaTargetGTEX_phenotype.txt`. Values are analyzed on the Xena-delivered scale; no cross-platform re-normalization is applied beyond within-contrast z-scoring for pathway scores.

## 2. Frozen gene sets

- PPAR/nuclear-receptor core: `PPARA; PPARD; PPARG; NR1I2; NR1I3; NR1H2; NR1H3`
- inflammatory complement: `RELA; STAT3`
- combined DINP-axis score: all `PPARA; PPARD; PPARG; NR1I2; NR1I3; NR1H2; NR1H3; RELA; STAT3`

## 3. How to read the outputs

`mcop_phase2f_bulk_gene_stats.csv` reports gene-level medians, tumor-minus-normal shifts, Mann–Whitney P values and BH-FDR separately for the TCGA-normal and GTEx-normal contrasts.

`mcop_phase2f_bulk_pathway_score_stats.csv` reports sample-level mean z-scores for the PPAR/NR, RELA+STAT3 and 9-gene programs. The score is a disease-state readout, not a causal mediation score.

`mcop_phase2f_tcga_paired_gene_stats.csv` and `mcop_phase2f_tcga_paired_pathway_score_stats.csv` use available patient-matched TCGA primary tumor/solid-normal pairs; these are the most internally comparable public checks in this run.

## 4. Single-cell status

CELLxGENE Census was selected for Phase 2F-B because its versioned metadata and primary-data filter are appropriate for cross-dataset single-cell queries. The current Windows runtime could not install `cellxgene-census` because the required TileDB-SOMA package had no compatible wheel and fell back to a local CMake build without a compiler. Therefore **no single-cell result is claimed in this commit**. The next single-cell run must pin the Census release and include `is_primary_data == True` before comparing malignant epithelial, myeloid and stromal compartments.

## 5. Interpretation boundary

The TCGA-internal paired shift supports CRC disease-state relevance, but the reference discordance prevents a clean mechanism upgrade. The next step should resolve normal-reference and tissue-composition effects before any TCGA mechanistic narrative is written. This result still does not prove exposure-specific direction or replace an exposure-linked tissue/perturbation experiment.

## Files

- `mcop_phase2f_bulk_gene_stats.csv`
- `mcop_phase2f_bulk_pathway_score_stats.csv`
- `mcop_phase2f_bulk_pathway_scores_by_sample.csv`
- `mcop_phase2f_bulk_sample_manifest.csv`
- `mcop_phase2f_bulk_coad_read_descriptive.csv`
- `mcop_phase2f_tcga_paired_gene_stats.csv`
- `mcop_phase2f_tcga_paired_pathway_score_stats.csv`
- `mcop_phase2f_bulk_direction_concordance.csv`
- `mcop_phase2f_figure_bulk_expression.png`

## Reproducibility

- Run UTC: `2026-08-22T13:36:03.868344+00:00`
- Script: `work\scripts\mcop_phase2f_crc_expression_validation.py`
- Xena hub: `https://toil.xenahubs.net`
- Expression dataset: `TcgaTargetGtex_rsem_gene_tpm`
- Phenotype dataset: `TcgaTargetGTEX_phenotype.txt`
- Gene-set hash: `12e6ae5dae94fb0bf1a1271fdfd0a2f8185153185e985bcdb2218f6b65349942`
- Python: `3.12.13`

**Phase 2F-A status: bulk tumor-normal analysis completed; result is reference-dependent/provisional. Single-cell validation remains a separately labeled Phase 2F-B task.**
