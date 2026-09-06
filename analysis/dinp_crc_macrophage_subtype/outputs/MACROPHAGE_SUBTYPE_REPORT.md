# Macrophage subtype analysis

- Source: local GSE144735 CRC matrix; 27,414 cells, 6 donors.
- Macrophage-only selection: 2,514 cells; cDC, granulocyte, mast, lymphoid and non-myeloid labels excluded.
- Recomputed macrophage-only PCA, neighbors, UMAP and Leiden at resolutions 0.3, 0.5, 0.7, 0.9, 1.1; selected resolution **0.5** by the predeclared stability-ARI rule after minimum cluster-size/donor-coverage filters.
- Frozen DINP–CRC program: 81 genes (79 present in the matrix); genes were not forced into HVGs.

## Integration answers

**A. Where is the program concentrated?** The highest donor-condition mean score is **SPP1_like_TAM** (0.030 if available). This is a descriptive donor-level ranking; formal paired results are in `subtype_program_statistics.csv`.

**B. Where is PTGER4 distributed?** The highest donor-condition mean is **C1QC_like_TAM** (0.189 if available). Detection fractions and paired results are in `donor_subtype_program_scores.csv` and `ptger4_subtype_statistics.csv`.

**C. Tumor versus normal.** Inference uses donor-level aggregates and paired tumor-normal tests only; BH-FDR is applied within each named statistic family. Border cells are retained for QC and descriptive plots, but are not used in paired tumor-normal tests.

**D. Abundance.** The highest mean tumor composition is **C1QC_like_TAM**. This is compositional and exploratory; see `subtype_abundance_statistics.csv`.

**E. Next step.** No further subtype work until a larger CRC cohort validates the state (top-program subtype has fewer than 3 paired tumor-normal donors). This module does not run CellChat, spatial, AUCell, deconvolution, pathway enrichment or causal analyses.

## Limitations and sanity checks

The earlier 370,115-cell Census object used for broad localization was unavailable locally, so this run uses the existing CRC GSE144735 cache and does not redownload data. Results should not be numerically pooled across the two source objects. One-donor clusters are flagged in `cluster_qc.csv`; marker specificity, contamination, frozen-gene coverage, resolution sensitivity and donor counts are retained in the output tables.
