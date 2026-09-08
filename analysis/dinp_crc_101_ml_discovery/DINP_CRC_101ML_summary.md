# DINP–CRC 101-model machine-learning discovery

## Frozen design

- Input universe: all 97 DINP–CRC overlap genes; Tier 1/cross-ranking labels were not used to select features or fit models.
- Training/discovery cohort: TCGA-COAD (tumor vs normal), stratified 5-fold cross-validation.
- Primary external validation: GSE10950 and GSE74602. GSE156355 was not used in this discovery/validation run.
- Models: 101 total = 11 feature-selection configurations × 9 classifiers + 2 full-feature baselines.
- Feature selection and scaling were fitted inside each training fold; the final model was refit on all TCGA-COAD samples before external validation.

## 101-model catalog

Feature selectors: all features; ANOVA F-test top 10/20/40/60/80/95; mutual-information top 10/20/40/60.
Classifiers: logistic L2, logistic L1, linear SVM, RBF SVM, random forest, extra-trees, histogram gradient boosting, gradient boosting, and distance-weighted kNN. Additional baselines: shrinkage LDA and Gaussian naive Bayes.

## Stability and external-validation ranking

`model_selection_frequency_selective_90` is the fraction of the 90 non-all-feature models whose full-TCGA fitted selector retained a gene. The CV-fold frequency is reported separately. External-validation-qualified models are selective models with ROC-AUC ≥0.75 in both GSE10950 and GSE74602. The ML priority score is the geometric mean of all-selective-model frequency, qualified-model frequency, and qualified-model external ROC-AUC score; it does not use Tier 1 labels.

- Frozen input genes: 97
- Measured expression features used for fitting: 95
- Input genes without expression values: 2 (ACP3, CCN2)
- Models completed: 101
- External-validation-qualified selective models (ROC-AUC ≥0.75 in both GEO cohorts): 18
- Default stable-ML flag (selection frequency ≥80% in all selective models and ≥80% in qualified models): 17
- Stable ML genes overlapping pre-existing Tier 1: 1

## Top 20 ML-priority genes

| ML rank | Gene | Stable flag | Selective frequency | Qualified frequency | Qualified mean external AUC | Existing cross-rank | Tier 1 overlap | Direction |
|---:|---|---|---:|---:|---:|---:|---|---|
| 1 | CD36 | yes | 1.000 | 1.000 | 0.881 | 3 | no | normal_high |
| 2 | CNR1 | yes | 1.000 | 1.000 | 0.881 | 59 | no | normal_high |
| 3 | INHBA | yes | 1.000 | 1.000 | 0.881 | 73 | no | mixed_or_flat |
| 4 | TIMP1 | yes | 1.000 | 1.000 | 0.881 | 34 | no | tumor_high |
| 5 | UBE2C | yes | 1.000 | 1.000 | 0.881 | 69 | no | tumor_high |
| 6 | VEGFA | yes | 1.000 | 1.000 | 0.881 | 91 | no | tumor_high |
| 7 | BMP2 | yes | 1.000 | 1.000 | 0.881 | 57 | no | normal_high |
| 8 | FASN | yes | 0.900 | 0.944 | 0.884 | 23 | no | tumor_high |
| 9 | ACACA | yes | 0.900 | 0.944 | 0.884 | 21 | no | tumor_high |
| 10 | FABP1 | yes | 0.900 | 0.944 | 0.884 | 46 | no | normal_high |
| 11 | NR3C1 | yes | 0.900 | 0.944 | 0.884 | 48 | no | normal_high |
| 12 | RXRA | yes | 0.900 | 0.944 | 0.884 | 6 | no | normal_high |
| 13 | MKI67 | yes | 0.800 | 0.944 | 0.884 | 77 | no | tumor_high |
| 14 | CEBPB | yes | 0.800 | 0.889 | 0.887 | 11 | yes | tumor_high |
| 15 | ANGPT2 | yes | 0.800 | 0.889 | 0.887 | 55 | no | tumor_high |
| 16 | EPAS1 | yes | 0.800 | 0.889 | 0.887 | 30 | no | normal_high |
| 17 | AGTR1 | yes | 0.800 | 0.889 | 0.887 | 41 | no | normal_high |
| 18 | ADIPOQ | no | 0.700 | 0.889 | 0.887 | 20 | yes | normal_high |
| 19 | CAT | no | 0.700 | 0.889 | 0.887 | 70 | no | normal_high |
| 20 | FAS | no | 0.700 | 0.889 | 0.887 | 18 | no | normal_high |

## Interpretation

This stage is discovery plus external classification validation for the tumor/normal expression phenotype. It does not establish DINP causality, and the high separability of tumor versus normal should not be confused with exposure-response evidence.

Tier 1/cross-rank overlap is an independent post hoc concordance check. It is not part of the ML feature-selection objective and should be reviewed before any MR instrument decision.
