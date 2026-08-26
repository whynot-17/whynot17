# Step 8D — T2D transcriptomic module directionality

- Status: **complete_module_directionality**
- Frozen Step 8C modules tested: **97**
- Primary unit: biological sample within each GEO series; no tissue pooling
- Module score: mean within-dataset gene-wise z-scores; minimum 3 mapped genes
- Test: T2D minus comparator Welch t test; BH-FDR within dataset

## Dataset audit

| Accession | Tissue | T2D | Control | Expression genes | Modules tested |
|---|---|---:|---:|---:|---:|
| GSE23343 | liver | 10 | 7 | 22189 | 46 |
| GSE21340 | skeletal_muscle | 5 | 10 | 5601 | 39 |
| GSE71416 | adipose | 14 | 6 | 22189 | 46 |
| GSE25724 | pancreatic_islet | 6 | 7 | 13299 | 42 |

## Tissue-stratified axis summary

This table audits module signs within each dataset and Tier A axis. It is not a pooled test and is not used to select a flagship axis.

| Accession | Tissue | Axis | Tested modules | Positive | Negative | Median delta | q<0.05 |
|---|---|---|---:|---:|---:|---:|---:|
| GSE21340 | skeletal_muscle | cluster_11 | 10 | 1 | 9 | -0.127 | 0 |
| GSE21340 | skeletal_muscle | cluster_5 | 3 | 1 | 2 | -0.036 | 0 |
| GSE21340 | skeletal_muscle | cluster_6 | 12 | 7 | 5 | 0.004 | 0 |
| GSE21340 | skeletal_muscle | cluster_8 | 14 | 3 | 11 | -0.067 | 0 |
| GSE23343 | liver | cluster_11 | 11 | 9 | 2 | 0.257 | 0 |
| GSE23343 | liver | cluster_5 | 3 | 2 | 1 | 0.176 | 0 |
| GSE23343 | liver | cluster_6 | 12 | 10 | 2 | 0.152 | 0 |
| GSE23343 | liver | cluster_8 | 20 | 18 | 2 | 0.146 | 0 |
| GSE25724 | pancreatic_islet | cluster_11 | 10 | 0 | 10 | -0.484 | 5 |
| GSE25724 | pancreatic_islet | cluster_5 | 3 | 3 | 0 | 0.566 | 2 |
| GSE25724 | pancreatic_islet | cluster_6 | 12 | 8 | 4 | 0.165 | 4 |
| GSE25724 | pancreatic_islet | cluster_8 | 17 | 8 | 9 | -0.032 | 9 |
| GSE71416 | adipose | cluster_11 | 11 | 1 | 10 | -0.283 | 0 |
| GSE71416 | adipose | cluster_5 | 3 | 1 | 2 | -0.149 | 1 |
| GSE71416 | adipose | cluster_6 | 12 | 1 | 11 | -0.202 | 0 |
| GSE71416 | adipose | cluster_8 | 20 | 10 | 10 | -0.004 | 2 |

## Directionality boundary

Cross-dataset synthesis is descriptive and retains tissue context. A positive or negative module score indicates a relative transcriptomic shift in the public series; it does not establish pathway activation, exposure causality, mediation, or a T2D-specific mechanism. No flagship axis is selected in Step 8D.

## Reproducibility

The series matrix and platform annotation URLs, sample-label rules, and SHA-256 checksums are recorded in `STEP8D_MANIFEST.json`. Raw GEO files are local-only and excluded from version control.
