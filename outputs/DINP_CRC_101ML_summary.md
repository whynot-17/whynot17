# DINP–CRC 101-model machine-learning discovery

## Frozen design

- Input universe: all 97 DINP–CRC overlap genes; Tier 1/cross-ranking labels were not used to select features or fit models.
- Training/discovery cohort: TCGA-COAD (tumor vs normal), stratified 5-fold cross-validation.
- Primary external validation: GSE10950 and GSE74602. GSE156355 was not used in this discovery/validation run.
- Models: 101 total = 11 feature-selection configurations × 9 classifiers + 2 full-feature baselines; 10 selective selector configurations are treated as the independent gene-stability units.
- Feature representation: within_sample_rank_percentile_0_to_1; this removes dependence on absolute RNA-seq versus microarray expression scales before fitting and testing.
- Feature selection, imputation, and scaling were fitted inside each training fold; the final model was refit on all rank-normalized TCGA-COAD samples before external validation.

## 101-model catalog

Feature selectors: all features; ANOVA F-test top 10/20/40/60/80/95; mutual-information top 10/20/40/60.
Classifiers: logistic L2, logistic L1, linear SVM, RBF SVM, random forest, extra-trees, histogram gradient boosting, gradient boosting, and distance-weighted kNN. Additional baselines: shrinkage LDA and Gaussian naive Bayes.

## Stability and external-validation ranking

Gene stability is defined as the fraction of the 10 independent selective selector configurations (six ANOVA and four mutual-information selectors) that retain a gene on the full TCGA training set. The 101 model configurations are used for classification-performance robustness. External-validation-qualified models are selective models with ROC-AUC ≥0.75 in both GSE10950 and GSE74602. A gene-level qualified-model AUC is the mean AUC of qualified multigene models containing that gene; it is not a single-gene AUC. The ML priority score integrates selector support with qualified-model inclusion and qualified-model AUC; it does not use Tier 1 labels.

- Frozen input genes: 97
- Measured expression features used for fitting: 95
- Input genes without expression values: 2 (ACP3, CCN2)
- Models completed: 101
- Independent selective selector configurations used for gene stability: 10
- External-validation-qualified selective models (ROC-AUC ≥0.75 in both GEO cohorts): 68
- Stable ML genes (support in ≥80% of the 10 independent selectors, i.e. ≥8/10): 17
- Stable ML genes overlapping pre-existing Tier 1: 2

## Top 20 ML-priority genes

| ML rank | Gene | Stable flag | Selector support | Qualified-model inclusion | Qualified mean AUC of models containing gene | Existing cross-rank | Tier 1 overlap | Direction |
|---:|---|---|---:|---:|---:|---:|---|---|
| 1 | ANGPT2 | yes | 10/10 | 1.000 | 0.983 | 57 | no | tumor_high |
| 2 | CNR1 | yes | 10/10 | 1.000 | 0.983 | 62 | no | normal_high |
| 3 | EPAS1 | yes | 10/10 | 1.000 | 0.983 | 36 | no | normal_high |
| 4 | INHBA | yes | 10/10 | 1.000 | 0.983 | 71 | no | normal_high |
| 5 | TIMP1 | yes | 10/10 | 1.000 | 0.983 | 38 | no | tumor_high |
| 6 | UBE2C | yes | 10/10 | 1.000 | 0.983 | 69 | no | tumor_high |
| 7 | VEGFA | yes | 10/10 | 1.000 | 0.983 | 91 | no | tumor_high |
| 8 | CD36 | yes | 10/10 | 1.000 | 0.983 | 3 | no | normal_high |
| 9 | BMP2 | yes | 10/10 | 1.000 | 0.983 | 58 | no | normal_high |
| 10 | CAT | yes | 9/10 | 0.897 | 0.982 | 72 | no | normal_high |
| 11 | MKI67 | yes | 8/10 | 0.794 | 0.984 | 77 | no | tumor_high |
| 12 | ACLY | yes | 8/10 | 0.794 | 0.984 | 39 | no | tumor_high |
| 13 | IL1A | yes | 8/10 | 0.794 | 0.984 | 2 | yes | tumor_high |
| 14 | CEBPB | yes | 8/10 | 0.794 | 0.984 | 9 | yes | tumor_high |
| 15 | RXRA | yes | 8/10 | 0.794 | 0.984 | 7 | no | normal_high |
| 16 | FAS | yes | 8/10 | 0.794 | 0.984 | 16 | no | normal_high |
| 17 | FASN | yes | 8/10 | 0.794 | 0.984 | 22 | no | tumor_high |
| 18 | AGTR1 | no | 7/10 | 0.691 | 0.983 | 41 | no | normal_high |
| 19 | NR3C1 | no | 7/10 | 0.691 | 0.983 | 49 | no | normal_high |
| 20 | PPARD | no | 7/10 | 0.691 | 0.983 | 43 | no | normal_high |

## Interpretation

This stage is discovery plus external classification validation for the tumor/normal expression phenotype. It does not establish DINP causality, and the high separability of tumor versus normal should not be confused with exposure-response evidence.

Tier 1/cross-rank overlap is an independent post hoc concordance check. It is not part of the ML feature-selection objective and should be reviewed before any MR instrument decision.
The qualified-model inclusion frequency is a model-configuration metric and should not be described as an independent selector-support probability.
