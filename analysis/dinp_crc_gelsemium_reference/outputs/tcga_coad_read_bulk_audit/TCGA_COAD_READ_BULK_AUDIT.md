# TCGA-COAD + TCGA-READ bulk audit of the fresh 41-gene set

## Audit status: **PASS**

This audit uses only the saved local expression matrices, sample manifests, run manifests, and the frozen `query_41_genes.csv`. It independently re-computes the gene-level Mann–Whitney statistics, effect summaries, BH-FDR values, and available patient-paired Wilcoxon checks. It does not pool COAD and READ into a new test family.

## Cohort-level verification

| Cohort | Primary tumor | Solid normal | Paired patients | FDR-positive genes | Errors | Warnings |
|---|---:|---:|---:|---:|---:|---|
| TCGA-READ | 92 | 10 | 6 | 18 | 0 | small paired sensitivity set (n=6); small independent normal reference (n=10) |
| TCGA-COAD | 288 | 41 | 26 | 31 | 0 | none |

## Cross-cohort comparison

- Shared independent BH-FDR<0.05 genes: **17/41**.
- Same non-zero tumor–normal direction in both cohorts: **33/41**.
- Spearman correlation of the 41 median shifts (READ vs COAD): **rho=0.787, P=1.06e-09**.
- These are cross-cohort expression-state concordance metrics, not evidence that either cohort measured DINP exposure or that DINP caused the expression shifts.

## Fixed-analysis checks

- Both cohorts use the same UCSC Toil Xena expression dataset and phenotype dataset.
- Both use the same frozen 41-gene input, with matching input SHA256 and gene order.
- Each cohort applies BH-FDR separately across its own 41-gene tumor–normal family.
- The saved expression matrices reproduce all stored independent medians, Cliff's delta, Mann–Whitney U, P values, and BH-FDR values within numerical tolerance.
- The saved sample manifests reproduce the stored tumor, normal, and paired-patient counts.
- The paired analyses are sensitivity checks only and do not re-rank the independent results.

## Interpretation boundary

The COAD and READ results provide independent CRC tissue-state checks for the fresh 41-gene intersection. The limited READ normal reference (n=10) and the tissue-composition difference between bulk tumors and solid normals require cautious interpretation. This audit does not establish exposure-specific direction, direct chemical binding, causality, or mechanism.

## Files

- `tcga_coad_read_cross_cohort_gene_audit.csv`: gene-level COAD/READ comparison and concordance flags.
- `tcga_coad_read_bulk_audit_manifest.json`: audit inputs, hashes, counts, and status.
- `audit_tcga_coad_read_bulk.py`: reproducible offline audit script.
