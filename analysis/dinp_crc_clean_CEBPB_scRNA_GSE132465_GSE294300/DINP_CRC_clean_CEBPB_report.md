# DINP–CRC clean CEBPB epithelial validation

## Analysis design

This final result package contains the two-cohort CEBPB epithelial validation based on GSE132465 and GSE294300. Cell-type annotation uses the GSE132465 curated broad-label reference; CEBPB is excluded from annotation features.

- Patient/sample-level pseudobulk from raw UMI counts.
- Tumor versus Normal paired Wilcoxon test.
- Minimum 50 confident epithelial cells per sample.
- BH-FDR across the two prespecified cohort-level epithelial tests.

## Primary results

- GSE132465: matched n=7; median Δ(Tumor−Normal)=0.607; P=0.046875; BH-FDR=0.046875; Tumor higher.
- GSE294300: matched n=10; median Δ(Tumor−Normal)=0.6655; P=0.0371094; BH-FDR=0.046875; Tumor higher.

## Summary

CEBPB shows a concordant tumor-higher epithelial signal in both retained cohorts. Both cohort-level tests remain below BH-FDR 0.05 in this two-cohort result family.

## Files

- `DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv`
- `DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv`
- `DINP_CRC_clean_CEBPB_tumor_normal_validation_all_celltypes.csv`
- `DINP_CRC_clean_CEBPB_cross_cohort_replication.csv`
- `DINP_CRC_clean_CEBPB_manifest.json`
- `DINP_CRC_clean_CEBPB_epithelial_tumor_normal.png
