# DINP–CRC CEBPB/CD36 single-cell localization and patient-level pseudobulk validation

Dataset: **GSE132465**. The public GEO series contains 63,689 annotated cells from 23 patients, including 47,285 tumor cells from 23 patients and 16,404 normal-mucosa cells from 10 patients. Matched tumor/normal patients: 10 (SMC01, SMC02, SMC03, SMC04, SMC05, SMC06, SMC07, SMC08, SMC09, SMC10).

## Scope and design

- Localization uses the GEO-provided curated cell-type and cell-subtype annotations; no cell type was inferred from target expression alone.
- CEBPB and CD36 raw UMI rows were extracted from the public matrix after exact cell-ID concordance checking.
- Pseudobulk aggregation is at patient × class × cell type, with a minimum of 50 cells per sample × cell type. Raw target UMI sums are normalized by aggregate library UMI to pseudobulk CPM; mean UMI per cell and detection rate are also retained.
- Validation uses only matched patients and a two-sided Wilcoxon signed-rank test on log1p(pseudobulk CPM), with Benjamini–Hochberg correction across the 12 gene × cell-type tests.
- This is an expression-localization/validation analysis; it does not establish direct DINP binding or causal CRC activity.

## Key localization result

| Gene | Highest tumor detection-rate cell type | Detection rate | Highest tumor mean-expression cell type | Mean log1p CPM |
|---|---|---:|---|---:|
| CEBPB | Myeloids | 0.7419 | Myeloids | 4.6473 |
| CD36 | Myeloids | 0.1809 | Myeloids | 0.9760 |

## Patient-level paired validation

See `GSE132465_CEBPB_CD36_paired_pseudobulk_validation.csv` for every gene × cell-type test and `scRNA_patient_pseudobulk_paired_validation.png` for the matched-patient trajectories.

| Gene | Cell type | Matched patients | Median tumor − normal log1p pseudobulk CPM | Wilcoxon P | BH FDR | Direction |
|---|---|---:|---:|---:|---:|---|
| CD36 | B cells | 9 | -0.4584 | 0.003906 | 0.03906 | Normal higher |
| CD36 | Stromal cells | 7 | -2.6776 | 0.01562 | 0.07812 | Normal higher |
| CD36 | T cells | 10 | -0.1752 | 0.04688 | 0.1172 | Normal higher |
| CD36 | Epithelial cells | 7 | -0.5555 | 0.0625 | 0.125 | Normal higher |
| CD36 | Myeloids | 3 | -0.6239 | 0.25 | 0.3125 | Normal higher |
| CD36 | Mast cells | 0 | NA | NA | NA | No difference |
| CEBPB | Epithelial cells | 7 | 0.6070 | 0.04688 | 0.1172 | Tumor higher |
| CEBPB | B cells | 9 | -0.2955 | 0.09766 | 0.1562 | Normal higher |
| CEBPB | Stromal cells | 7 | -0.2661 | 0.1094 | 0.1562 | Normal higher |
| CEBPB | T cells | 10 | 0.0265 | 0.4922 | 0.5469 | Tumor higher |
| CEBPB | Myeloids | 3 | -0.0647 | 1 | 1 | Normal higher |
| CEBPB | Mast cells | 0 | NA | NA | NA | No difference |

## Provenance and QC

- Annotation download: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz`
- Raw UMI download: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz`
- Target rows found: CD36, CEBPB; cell-level rows: 63,689; patient/cell-type pseudobulk rows: 354.
- SHA256 values and machine-readable parameters are in `GSE132465_CEBPB_CD36_analysis_manifest.json`.
- The large public raw matrix is not copied into the repository; the download URLs and SHA256 hashes make the analysis auditable.
