# CRC sidedness stable DEG discovery

## Question

Which genes show reproducible right-versus-left differences across bulk CRC cohorts and persist in tumor epithelial single-cell pseudobulk?

## Frozen design

- Bulk: TCGA-COAD, GSE39582, GSE103479 and GSE41258.
- Single cell: GSE200997, GSE132465, GSE144735 and GSE188711.
- Bulk arrays/RNA-seq were analyzed within cohort; no cross-platform expression matrix was concatenated.
- Single-cell raw UMI counts were summed within patient and tumor-epithelial compartment before testing. No cell-level p-values were used.
- Differential expression is Welch right-versus-left testing on the native normalized/log expression scale or patient pseudobulk log2(CPM+1).
- Meta-analysis uses Hedges-g standardized effects with a DerSimonian–Laird random-effects model.
- Strict stable right-high/left-high: at least 3 contributing cohorts, all effects in the same direction, random-effects BH-FDR < 0.05.

## Counts and strict lists

- Bulk cohort-level tables: 4; strict right-high genes: **1278**; strict left-high genes: **1198**.
- Single-cell cohort-level tables: 4; strict right-high genes: **0**; strict left-high genes: **0**.
- Strict bulk–single-cell overlap: right-high **0**, left-high **0**.

The strict lists are intentionally conservative. The full meta tables retain effect direction, heterogeneity, cohort count and FDR for genes that are biologically interesting but underpowered.

## Bulk cohort audit

| Cohort | n right | n left | genes | source |
|---|---:|---:|---:|---|
| GSE39582 | 232 | 351 | 21750 | C:\Users\21634\Documents\Codex\2026-08-31\1-janssen-kp-et-al-extrinsic\work\data\GSE39582_series_matrix.clean.txt.gz |
| GSE103479 | 56 | 59 | 25382 | C:\Users\21634\Documents\Codex\2026-08-31\1-janssen-kp-et-al-extrinsic\work\data\GSE103479_log2_RMA_annotated.txt.gz |
| GSE41258 | 105 | 132 | 13100 | C:\Users\21634\Documents\Codex\2026-08-31\1-janssen-kp-et-al-extrinsic\work\data\GSE41258_series_matrix.txt.gz |
| TCGA-COAD | 171 | 99 | 20530 | C:\Users\21634\Documents\Codex\2026-08-31\1-janssen-kp-et-al-extrinsic\work\data\TCGA_COAD_HiSeqV2.gz |

## Single-cell cohort audit

| Cohort | n right | n left | selected cells | selection | warning |
|---|---:|---:|---:|---|---|
| GSE200997 | 8 | 8 | 6681 | existing marker-defined malignant epithelial cells |  |
| GSE132465 | 11 | 9 | 17469 | official GEO tumor epithelial annotation; malignancy inferred from tumor origin |  |
| GSE144735 | 2 | 4 | 2212 | official GEO tumor epithelial annotation; malignancy inferred from tumor origin |  |
| GSE188711 | 3 | 3 | 2370 | existing primary putative malignant epithelial Leiden cluster 4 | 3 left and 3 right patients; directional sensitivity only |

## Interpretation guardrails

This is a phenotype-first sidedness screen. It does not by itself establish causality, cell-state mechanism, AA/PUFA flux, ferroptosis, treatment response or therapeutic dependency. GSE132465/GSE144735 tumor epithelial cells are annotation-based tumor-epithelial proxies; GSE188711 is a low-powered 3-versus-3 sensitivity cohort.

## Outputs

- `crc_sidedness_bulk_de_<cohort>.csv`: all gene-level within-cohort statistics.
- `crc_sidedness_bulk_random_effects_meta.csv`: bulk random-effects meta-analysis.
- `crc_sidedness_sc_de_<cohort>.csv`: all gene-level single-cell pseudobulk statistics.
- `crc_sidedness_sc_random_effects_meta.csv`: single-cell random-effects meta-analysis.
- `crc_sidedness_stable_deg_bulk_sc_intersection.csv`: all genes with bulk/sc flags.
- `crc_sidedness_stable_deg_strict_intersection.csv`: strict bulk–single-cell overlap only.
- `crc_sidedness_bulk_side_audit.csv` and `crc_sidedness_sc_patient_audit.csv`: sample/patient sidedness and selection audit.
