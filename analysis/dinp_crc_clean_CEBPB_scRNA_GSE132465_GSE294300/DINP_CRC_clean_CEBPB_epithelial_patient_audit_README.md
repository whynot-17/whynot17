# Patient-level CEBPB epithelial audit

- Source frozen table: F:/Codex/archived-runs/2026-09-08/codex-dinp-target-mining-crc-mr/analysis/dinp_crc_clean_CEBPB_scRNA_GSE132465_GSE294300/DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv
- Filter: gene_symbol=CEBPB; cell_type=Epithelial cells; cohorts=GSE132465 and GSE294300.
- CPM and log1p CPM values are copied from the frozen patient/sample pseudobulk table.
- delta_pseudobulk_cpm and delta_log1p_tumor_minus_normal are direct Tumor minus Normal differences of frozen pseudobulk columns.
- Complete Tumor/Normal patient pairs in audit: 28.
- Primary-included pairs passing the frozen minimum-cell rule: 17.
- No p-values, FDRs, or inferential statistics were recomputed.

## Primary pairs by cohort
- GSE132465: 7.
- GSE294300: 10.
