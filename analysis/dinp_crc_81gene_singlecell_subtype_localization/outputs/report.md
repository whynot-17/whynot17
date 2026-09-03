# Frozen 81-gene DINP–CRC program: single-cell subtype localization

Generated: 2026-09-03T07:12:32.443450+00:00

## Analysis boundary

- Source: Census release `2025-11-08` H5AD; source file is not copied into the repository.
- Source shape: `370,115 cells × 38,361 features`; four-compartment eligible cells: `223,160`.
- Frozen program: 81 genes from `dinp_crc_intersection.csv`; all 81 had unique exact matches in `feature_name`.
- Primary score: the same global gene-wise z-score score used in the broad-compartment localization analysis.
- Inference: donor-level means and paired tumor-minus-normal contrasts; BH-FDR is across six subtype paired t-tests.
- Sensitivity: genes standardized within each broad compartment before scoring.

## Label audit

`ClusterFull` labels beginning with `Tumor` were used to define a source-labeled tumor epithelial / malignant-candidate subgroup. This is not treated as definitive malignant status because no CNV-based or independent malignant-cell validation was performed.

## Paired tumor–normal subtype contrasts

| Contrast | Paired donors | Mean Δ (tumor−normal) | 95% CI | t-test P | BH-FDR | Wilcoxon P | Direction |
|---|---:|---:|---|---:|---:|---:|---|
| macrophage tumor vs normal | 35 | 0.105 | 0.074 to 0.135 | 5.67e-08 | 3.4e-07 | 3.73e-08 | up |
| source tumor-labeled epithelial vs normal epithelial | 36 | 0.035 | 0.018 to 0.052 | 0.000159 | 0.000478 | 0.000148 | up |
| monocyte tumor vs normal | 33 | 0.044 | 0.021 to 0.067 | 0.000421 | 0.000842 | 0.000316 | up |
| dendritic tumor vs normal | 34 | 0.047 | -0.000 to 0.093 | 0.0516 | 0.0775 | 0.00168 | up |
| granulocyte tumor vs normal | 3 | 0.071 | -0.120 to 0.263 | 0.251 | 0.302 | nan | up |
| other tumor epithelial vs normal epithelial | 36 | -0.009 | -0.033 to 0.015 | 0.461 | 0.461 | 0.592 | down |

## Within-compartment standardization sensitivity

| Contrast | Paired donors | Mean Δ (tumor−normal) | t-test P | BH-FDR | Direction |
|---|---:|---:|---:|---:|---|
| macrophage tumor vs normal | 35 | 0.079 | 8.05e-07 | 4.83e-06 | up |
| source tumor-labeled epithelial vs normal epithelial | 36 | 0.042 | 6.16e-05 | 0.000185 | up |
| monocyte tumor vs normal | 33 | 0.028 | 0.00411 | 0.00823 | up |
| granulocyte tumor vs normal | 3 | 0.060 | 0.16 | 0.241 | up |
| dendritic tumor vs normal | 34 | 0.026 | 0.258 | 0.31 | up |
| other tumor epithelial vs normal epithelial | 36 | -0.005 | 0.724 | 0.724 | down |

## Interpretation boundary

This analysis localizes the frozen DINP–CRC program in a CRC single-cell reference. It does not establish that DINP exposure causes the program, that the program mediates the epidemiologic association, or that the source Tumor-prefixed cells are definitively malignant without independent validation.

Detailed donor-level scores, paired deltas, label audit, cell counts, source hash, and gene mapping audit are retained in the output directory.
