# DINP–CRC 101-model ML audit log

```json
{
  "seed": 20260909,
  "python": "3.13.9",
  "scikit_learn": "1.7.2",
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
  "external_validation_qualified_model_count": 18,
  "cv": {
    "method": "StratifiedKFold",
    "n_splits": 5,
    "shuffle": true,
    "random_state": 20260909
  },
  "feature_selection_leakage_control": "selectors and StandardScaler fitted within each training fold",
  "external_validation_qualified_definition": "selective model with GSE10950 ROC-AUC >= 0.75 and GSE74602 ROC-AUC >= 0.75",
  "stable_ml_definition": "full-TCGA selection frequency >= 0.80 in all 90 selective models and >= 0.80 in external-validation-qualified selective models",
  "tier1_used_in_model_fit": false,
  "qc": {
    "model_ids_unique": true,
    "gene_symbols_unique": true,
    "output_gene_rows": 97,
    "external_auc_range": [
      0.45215277777777774,
      0.9457465277777779
    ],
    "stable_ml_count": 17,
    "stable_ml_tier1_overlap_count": 1
  }
}
```
