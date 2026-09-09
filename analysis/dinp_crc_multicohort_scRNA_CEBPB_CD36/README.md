# DINP–CRC CEBPB/CD36 multicohort single-cell validation

This directory contains the frozen three-cohort analysis for the two prioritized candidates, **CEBPB** and **CD36**:

- **GSE200997**: discovery; cell-type localization, matched tumor–normal validation, and secondary left/right sensitivity.
- **GSE132465**: independent replication; cell-type localization and matched tumor–normal validation.
- **GSE188711**: sidedness-only sensitivity; no normal samples and only three left/right samples per side.

The analysis uses patient/sample-level pseudobulk. Cell-level observations are not treated as independent biological replicates. The primary cell-type pseudobulk threshold is 50 cells per sample × cell type. Because GSE200997 marker-inferred normal myeloid cells are sparse, a separate broad-compartment analysis uses a minimum of 10 cells per sample × compartment and is explicitly exploratory.

## Annotation audit

GSE132465 uses the previously downloaded GEO-curated cell labels. The public GSE200997 annotation file used here contains sample, condition, location, MSI and CMS metadata but no cell-type labels; GSE188711 provides raw 10X MTX/TSV files without a cell-type annotation table. Major compartments in those two cohorts were therefore inferred from fixed canonical marker panels and are labeled `Marker-inferred` in the tables. “Tumor epithelial” is a tumor-sample epithelial proxy, not CNV-confirmed malignant epithelium.

The GSE188711 raw supplementary files contain 35,693 barcodes, whereas the source publication reports 27,927 high-quality cells. This re-analysis applies the stated UMI/gene QC and does not claim to reproduce the original filtering exactly; GSE188711 is used only for sidedness sensitivity.

## Re-running

The raw GEO files are intentionally not committed. Download/prepare the source files under the repository `work/` paths expected by the script, and make the prior GSE132465 cell-level output available under `outputs/DINP_CRC_scRNA_GSE132465_CEBPB_CD36/`. Then run:

```powershell
& C:\Users\21634\anaconda3\python.exe analysis/dinp_crc_multicohort_scRNA_CEBPB_CD36/run_multicohort_crc_scRNA_cebpb_cd36.py
```

Source records:

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188711

## Key files

- `DINP_CRC_multicohort_CEBPB_CD36_localization.csv`: target localization by cohort/class/side/cell type.
- `DINP_CRC_multicohort_CEBPB_CD36_patient_sample_pseudobulk.csv`: primary patient/sample-level pseudobulk table.
- `DINP_CRC_GSE200997_GSE132465_paired_tumor_normal_validation.csv`: primary paired tumor–normal tests.
- `DINP_CRC_GSE200997_GSE132465_cross_compartment_paired_validation.csv`: exploratory broad-compartment paired tests.
- `DINP_CRC_GSE200997_GSE132465_cross_cohort_replication.csv`: cross-cohort direction-concordance audit.
- `DINP_CRC_GSE200997_GSE132465_GSE188711_left_right_sensitivity.csv`: secondary sidedness analysis.
- `DINP_CRC_multicohort_CEBPB_CD36_report.md`: methods, caveats, cohort counts and interpretation guardrails.
- `DINP_CRC_multicohort_CEBPB_CD36_manifest.json`: checksums and machine-readable provenance.
