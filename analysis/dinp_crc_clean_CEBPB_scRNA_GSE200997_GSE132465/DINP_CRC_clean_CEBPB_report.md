# DINP–CRC clean CEBPB single-cell rerun

## Frozen design

This rerun is CEBPB-only and uses GSE200997 plus GSE132465. GSE188711 was excluded. No sidedness/left-right field is read, derived, grouped, or exported.

- GSE132465: raw UMI matrix reconstructed independently from `work/scRNA_GSE132465/raw_UMI_count_matrix.txt.gz`; curated `cell_type` labels were read from the source annotation file.
- GSE200997: raw UMI matrix reconstructed independently from `work/scRNA_GSE200997/raw_UMI_count_matrix.csv.gz`; its source annotation was used only for cell ID, sample, and Tumor/Normal condition.
- GSE200997 cell types: nearest-centroid mapping to GSE132465 curated broad labels using canonical marker expression, with reference-calibrated similarity and margin thresholds. Low-confidence and ambiguous cells are retained in the audit table but excluded from the primary pseudobulk analysis.
- Primary validation: matched patient-level epithelial pseudobulk, Tumor versus Normal, paired Wilcoxon, minimum 50 cells per sample × cell type.
- Secondary localization/compartment family: BH-FDR across the six prespecified cell types within each cohort. Primary epithelial family: BH across the two epithelial cohort tests only.

## QC and annotation audit

- Total cells loaded: 113,548; QC-passing cells: 113,535.
- QC thresholds: total UMI ≥ 500; detected genes ≥ 200.
- GSE200997 cells retained as confident reference-mapped labels: 39,909; low-confidence: 9,615; ambiguous/doublet-like: 335.
- Reference leave-one-out label accuracy: 0.971; similarity threshold: 0.3130; margin threshold: 0.1338.
- GSE200997 epithelial labels are therefore reference-mapped epithelial cells, not CNV-confirmed malignant cells; the Tumor sample context is not treated as a malignancy annotation.

## Primary epithelial validation

- GSE200997: matched n=3; median Δ(Tumor−Normal)=-0.0161; Wilcoxon p=1.0; primary epithelial BH-FDR=1.0; Normal higher.
- GSE132465: matched n=7; median Δ(Tumor−Normal)=0.6070; Wilcoxon p=0.046875; primary epithelial BH-FDR=0.09375; Tumor higher.

## Cross-cohort epithelial replication

- Direction concordance: False; status: Not concordant/insufficient.

## Files

- `DINP_CRC_clean_CEBPB_cell_level.csv`: cell-level CEBPB metrics and annotation/QC status.
- `DINP_CRC_clean_GSE200997_annotation_audit.csv`: complete GSE200997 reference-mapping audit, including top-two labels, similarities, margins, and confidence status.
- `DINP_CRC_clean_CEBPB_localization.csv`: cohort/condition/cell-type localization summary.
- `DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv`: patient/sample-level pseudobulk.
- `DINP_CRC_clean_CEBPB_tumor_normal_validation_all_celltypes.csv`: secondary six-cell-type paired validation with cohort-specific BH-FDR.
- `DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv`: prespecified epithelial validation with BH-FDR across the two cohorts.
- `DINP_CRC_clean_CEBPB_cross_cohort_replication.csv`: cross-cohort direction/replication audit.
- `DINP_CRC_clean_CEBPB_manifest.json`: source hashes, thresholds, and exact analysis settings.

The previous three-cohort CEBPB/CD36 module remains frozen and is not used as input here.
