# DINP–CRC CEBPB/CD36 multicohort single-cell analysis

## Frozen cohort design

- **GSE200997**: discovery cohort for cell-type localization, tumor–normal pseudobulk, and left–right sensitivity.
- **GSE132465**: independent replication cohort for cell-type localization and matched tumor–normal pseudobulk; tumor region is also used for a secondary left–right sensitivity check.
- **GSE188711**: sidedness-only supplementary cohort; no normal samples, so it is not used for tumor–normal validation.

## Important annotation audit

GSE132465 uses its GEO-curated cell-type labels. The downloaded GSE200997 annotation file contains sample, condition, location, MSI, and CMS metadata but no cell-type labels; GSE188711 provides raw 10X MTX/TSV files without a cell-type annotation table. Their major cell compartments were therefore inferred with fixed canonical marker-score panels and are explicitly labeled `Marker-inferred` in the outputs. Tumor epithelial cells are treated as a tumor-sample epithelial proxy, not as CNV-confirmed malignant cells.

## Statistical design

Raw UMI matrices were filtered at ≥500 total UMI and ≥200 detected genes for GSE200997/GSE188711; GSE132465 is the already-filtered GEO processed cell table used in the previous analysis. Pseudobulk was aggregated at patient × sample × class × cell type and normalized as aggregate target UMI / aggregate library UMI × 1e6. Tumor–normal comparisons use paired Wilcoxon tests across matched patients; left–right comparisons use sample-level Mann–Whitney U tests. A minimum of 50 cells per sample × cell type is required, and BH correction is applied within each analysis family.

## Cohort sizes after QC

- GSE200997: 31,573 tumor cells and 18,273 normal cells retained.
- GSE132465: 47,285 tumor cells and 16,404 normal cells retained.
- GSE188711: 35,693 raw barcodes; 33,482 cells retained after QC (2,211 failed) for sidedness analysis. The GEO raw supplementary files contain more barcodes than the source paper's reported 27,927 high-quality cells; this re-analysis does not claim to reproduce the paper's original filtering exactly.

## Key outputs

- `DINP_CRC_multicohort_CEBPB_CD36_localization.csv`: localization by cohort, class, and cell type.
- `DINP_CRC_multicohort_CEBPB_CD36_patient_sample_pseudobulk.csv`: patient/sample-level pseudobulk table.
- `DINP_CRC_GSE200997_GSE132465_paired_tumor_normal_validation.csv`: primary tumor–normal validation.
- `DINP_CRC_GSE200997_GSE132465_cross_compartment_paired_validation.csv`: exploratory broad-compartment tumor–normal validation using ≥10 cells per sample × compartment to retain sparse GSE200997 myeloid context.
- `DINP_CRC_GSE200997_GSE132465_cross_compartment_patient_sample_pseudobulk.csv`: broad-compartment pseudobulk source table for the exploratory validation.
- `DINP_CRC_GSE200997_GSE132465_GSE188711_left_right_sensitivity.csv`: secondary sidedness analysis.
- `DINP_CRC_GSE200997_GSE132465_cross_cohort_replication.csv`: direction-concordance audit between the two tumor–normal cohorts.

## Provisional interpretation

Use the cross-cohort replication table, not cell-level p-values, to decide whether CEBPB/CD36 have reproducible tumor–normal effects within the same broad compartment. The primary cell-type validation requires ≥50 cells per sample × cell type; the separate broad-compartment myeloid analysis uses ≥10 cells and is exploratory because marker-inferred normal myeloid cells are sparse in GSE200997. Myeloid localization and epithelial tumor-sample localization remain descriptive unless the direction is reproduced at patient level in both cohorts.
