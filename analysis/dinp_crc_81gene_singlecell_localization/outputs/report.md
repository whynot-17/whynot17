# Frozen 81-gene DINP–CRC program: single-cell compartment localization

Generated: 2026-09-03T02:15:42.014672+00:00

## Analysis boundary

- Source: official Census-release source H5AD, pinned to `2025-11-08`; source file is not copied into the repository.
- Source shape: `370,115 cells × 38,361 features`; eligible cells across the four compartments: `223,160`.
- Frozen program: 81 genes from `dinp_crc_intersection.csv`; all 81 had unique exact matches in `feature_name`.
- Expression scale: `adata.X` (non-raw source matrix). Each gene was z-scored across all eligible cells, then averaged into the program score.
- Inference: donor-level means; the primary localization contrast is paired tumor-minus-normal within donor. BH-FDR is across the four compartment paired t-tests.
- Sensitivity: the same 81 genes were standardized separately within each compartment before scoring; this is not pooled with the primary score.
- The analysis does not infer malignant status. The epithelial label is tumor-derived epithelial when the source disease label is colon adenocarcinoma.

## Paired tumor–normal localization

| Compartment | Paired donors | Mean Δ (tumor−normal) | 95% CI | t-test P | BH-FDR | Wilcoxon P | Direction |
|---|---:|---:|---|---:|---:|---:|---|
| myeloid | 36 | 0.076 | 0.045 to 0.107 | 1.51e-05 | 6.06e-05 | 1.15e-06 | up |
| epithelial | 36 | 0.034 | 0.017 to 0.051 | 0.000196 | 0.000391 | 0.000137 | up |
| fibroblast | 32 | -0.132 | -0.210 to -0.054 | 0.0016 | 0.00213 | 0.00275 | down |
| endothelial | 34 | 0.016 | -0.016 to 0.048 | 0.312 | 0.312 | 0.196 | up |

## Within-compartment standardization sensitivity

| Compartment | Paired donors | Mean Δ (tumor−normal) | t-test P | BH-FDR | Direction |
|---|---:|---:|---:|---:|---|
| epithelial | 36 | 0.041 | 7.19e-05 | 0.000288 | up |
| myeloid | 36 | 0.052 | 0.000413 | 0.000827 | up |
| endothelial | 34 | 0.011 | 0.407 | 0.542 | up |
| fibroblast | 32 | 0.004 | 0.858 | 0.858 | up |

## Interpretation boundary

This is a localization/convergence analysis of the frozen DINP–CRC intersection. A compartment-level expression shift does not establish that DINP exposure causes the shift, nor that the program mediates the epidemiologic association.

Detailed donor-level scores, paired deltas, cell-type summaries, source hashes, and gene mapping audit are retained in the output directory.
