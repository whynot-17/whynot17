# DINP–CRC 101-model machine-learning discovery

This module implements the frozen discovery → validation design:

`97 DINP–CRC overlap genes → 101 ML model/algorithm combinations → stable ML genes → external GEO validation → independent comparison with the pre-existing PPI/pathway cross-ranking`

## Design

- Training/discovery cohort: TCGA-COAD tumor versus normal, with stratified 5-fold cross-validation.
- External validation cohorts: GSE10950 and GSE74602.
- Input universe: all 97 overlap genes. Existing Tier 1 labels and cross-ranking scores were not used during feature selection or model fitting.
- Expression features used for fitting: 95. ACP3 and CCN2 were retained in the 97-gene audit universe but had no usable expression values in the shared matrix and were not imputed as features.
- Before model fitting, the 95 measured genes are rank/percentile-normalized within each sample to the 0–1 interval, reducing dependence on absolute RNA-seq versus microarray scales.
- Missing expression cells among measured features were median-imputed inside each training fold by the pipeline.

## 101-model catalog

The catalog contains 11 feature-selection configurations × 9 classifiers, plus two full-feature baselines:

- Selectors: all features; ANOVA F-test top 10/20/40/60/80/95; mutual-information top 10/20/40/60.
- Classifiers: logistic L2, logistic L1, linear SVM, RBF SVM, random forest, extra-trees, histogram gradient boosting, gradient boosting, and distance-weighted kNN.
- Baselines: shrinkage LDA and Gaussian naive Bayes.

Scaling, imputation, and supervised feature selection are fitted within each training fold. The final model for each combination is then refit on all rank-normalized TCGA-COAD samples before GSE prediction.

## Stability and external validation

Gene stability is the support frequency across the 10 independent selective selector configurations: six ANOVA F-test selectors and four mutual-information selectors. The nine classifiers attached to a given selector do not create nine independent gene-selection votes.

The 101 model configurations are retained for classification-performance robustness. External-validation-qualified models are selective models with rank-normalized ROC-AUC ≥0.75 in both GSE10950 and GSE74602.

`ml_priority_score` is the geometric mean of:

1. support frequency across the 10 independent selectors;
2. inclusion frequency across qualified model configurations;
3. a 0–1 rescaling of the mean external ROC-AUC of qualified multigene models containing the gene.

The stable-ML flag requires support in at least 8/10 independent selectors. Qualified-model inclusion and AUC remain secondary performance-robustness metrics and are not described as a single-gene AUC.

## Key result

- 101 models completed.
- Stable-ML gene count is determined by the corrected ≥8/10 selector-support rule.
- The overlap with the pre-existing Tier 1 PPI/pathway/transcriptomic cross-supported set is reported from the corrected output.

The Tier 1 comparison is post hoc and independent; it is not part of the ML objective.

## Files

- `DINP_CRC_101ML_model_catalog.csv`: exact 101-model definitions.
- `DINP_CRC_101ML_model_results.csv`: internal CV and external GEO metrics for every model.
- `DINP_CRC_101ML_gene_stability.csv`: full 97-gene stability, external-validation, ML-priority, and Tier 1 comparison table.
- `DINP_CRC_101ML_stable_genes_top20.csv`: default stable-ML subset.
- `DINP_CRC_101ML_summary.md`: methods, thresholds, and top-20 table.
- `DINP_CRC_101ML_log.md`: reproducibility and QC metadata.
- `run_101_ml.py`: reproducible runner.
- `DINP_CRC_TCGA_GEO_target_expression_long.csv`, `DINP_CRC_TCGA_GEO_sample_manifest.csv`: frozen expression input and sample labels.
- `DINP_CRC_overlap.csv`: frozen 97-gene input universe.
- `DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv`: independent comparison table.

Run from this directory:

```bash
python run_101_ml.py
```

This is tumor-versus-normal expression classification. It does not demonstrate DINP causality or establish MR instruments.
