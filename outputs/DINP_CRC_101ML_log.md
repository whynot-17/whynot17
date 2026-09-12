# DINP–CRC 101-model ML audit log

```json
{
  "seed": 20260909,
  "python": "3.13.9",
  "scikit_learn": "1.9.0",
  "input_gene_count": 97,
  "measured_feature_count": 95,
  "unmeasured_input_genes": [
    "ACP3",
    "CCN2"
  ],
  "training_dataset": "TCGA-COAD",
  "external_validation_datasets": [
    "GSE10950",
    "GSE74602"
  ],
  "excluded_from_this_run": [
    "GSE156355",
    "TCGA paired sensitivity",
    "TCGA-vs-GTEx sensitivity"
  ],
  "feature_normalization": "within_sample_rank_percentile_0_to_1",
  "training_sample_counts": {
    "normal": 41,
    "tumor": 288
  },
  "external_sample_counts": {
    "GSE10950": {
      "normal": 24,
      "tumor": 24
    },
    "GSE74602": {
      "normal": 30,
      "tumor": 30
    }
  },
  "model_count": 101,
  "selective_model_count": 90,
  "independent_selective_selector_count": 10,
  "external_validation_qualified_model_count": 68,
  "cv": {
    "method": "StratifiedKFold",
    "n_splits": 5,
    "shuffle": true,
    "random_state": 20260909
  },
  "feature_selection_leakage_control": "rank normalization is computed within each sample; selectors, imputation and StandardScaler are fitted within each training fold",
  "external_validation_qualified_definition": "selective model with rank-normalized GSE10950 ROC-AUC >= 0.75 and rank-normalized GSE74602 ROC-AUC >= 0.75",
  "stable_ml_definition": "full-TCGA support >= 0.80 in the 10 independent selective selector configurations (>=8/10); qualified-model inclusion is reported separately as performance robustness",
  "qualified_model_gene_metric_definition": "mean external metric across qualified multigene model configurations containing the gene; not a single-gene AUC",
  "tier1_used_in_model_fit": false,
  "qc": {
    "model_ids_unique": true,
    "gene_symbols_unique": true,
    "output_gene_rows": 97,
    "external_auc_range": [
      0.5,
      1.0
    ],
    "stable_ml_count": 17,
    "stable_ml_tier1_overlap_count": 2,
    "selector_support_set_consistent_across_classifiers": true,
    "selector_full_support_counts_range": [
      1,
      10
    ]
  }
}
```
