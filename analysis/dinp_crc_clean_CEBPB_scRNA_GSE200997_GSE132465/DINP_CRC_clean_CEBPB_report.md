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
- GSE200997 cells retained as confident reference-mapped labels: 46,013; low-confidence: 2,620; ambiguous/doublet-like: 1,226.
- Reference LOO median-centroid calibration: accuracy=0.942; top1 median=0.5602; margin median=0.5263.
- Strict thresholds: similarity=0.1750; margin=0.0633. Very-strict thresholds: similarity=0.2577; margin=0.1376.
- GSE200997 epithelial labels are therefore reference-mapped epithelial cells, not CNV-confirmed malignant cells; the Tumor sample context is not treated as a malignancy annotation.

## GSE200997 condition-dependent annotation audit

- Tumor: QC cells=31573; Confident=92.37%; Low-confidence=5.35%; Ambiguous=2.29%; Confident epithelial=20.46%.
- Normal: QC cells=18273; Confident=92.17%; Low-confidence=5.08%; Ambiguous=2.74%; Confident epithelial=12.55%.

The full label composition is exported separately for both top1 reference labels among all QC cells and confident final labels. These condition-stratified percentages are the audit for potential condition-dependent filtering.

## GSE200997 epithelial inclusion audit

The prior contaminated module reported 7 matched epithelial patients. In this clean strict reference-mapped analysis, only patients passing the 50-cell rule in both conditions enter the paired test; the table below records every exclusion reason.
- cac1: Tumor=1211, Normal=0; Tumor pass50=true, Normal pass50=false; excluded; normal sample missing.
- cac10: Tumor=46, Normal=304; Tumor pass50=false, Normal pass50=true; excluded; tumor epithelial <50.
- cac11: Tumor=752, Normal=379; Tumor pass50=true, Normal pass50=true; included; both conditions pass 50 cells.
- cac12: Tumor=273, Normal=0; Tumor pass50=true, Normal pass50=false; excluded; normal sample missing.
- cac13: Tumor=22, Normal=0; Tumor pass50=false, Normal pass50=false; excluded; tumor epithelial <50; normal sample missing.
- cac14: Tumor=915, Normal=127; Tumor pass50=true, Normal pass50=true; included; both conditions pass 50 cells.
- cac15: Tumor=51, Normal=392; Tumor pass50=true, Normal pass50=true; included; both conditions pass 50 cells.
- cac16: Tumor=76, Normal=0; Tumor pass50=true, Normal pass50=false; excluded; normal sample missing.
- cac2: Tumor=985, Normal=0; Tumor pass50=true, Normal pass50=false; excluded; normal sample missing.
- cac3: Tumor=671, Normal=0; Tumor pass50=true, Normal pass50=false; excluded; normal sample missing.
- cac4: Tumor=47, Normal=363; Tumor pass50=false, Normal pass50=true; excluded; tumor epithelial <50.
- cac5: Tumor=2, Normal=0; Tumor pass50=false, Normal pass50=false; excluded; tumor epithelial <50; normal sample missing.
- cac6: Tumor=49, Normal=155; Tumor pass50=false, Normal pass50=true; excluded; tumor epithelial <50.
- cac7: Tumor=889, Normal=573; Tumor pass50=true, Normal pass50=true; included; both conditions pass 50 cells.
- cac8: Tumor=289, Normal=0; Tumor pass50=true, Normal pass50=false; excluded; normal sample missing.
- cac9: Tumor=181, Normal=0; Tumor pass50=true, Normal pass50=false; excluded; normal sample missing.

## Threshold sensitivity

- strict: matched n=4; median Δ=0.1975; P=0.625; direction=Tumor higher; included epithelial cells=8752.
- relaxed: matched n=4; median Δ=0.2047; P=0.625; direction=Tumor higher; included epithelial cells=8839.
- very_strict: matched n=3; median Δ=-0.008762; P=1; direction=Normal higher; included epithelial cells=8147.

## Non-transfer epithelial gate sensitivity

A separate canonical-marker gate excluding CEBPB from classification retained 9894 cells and produced matched n=7, median Δ=0.3116, P=0.2969, direction=Tumor higher.

## Interpretation

GSE132465 showed a tumor-high epithelial CEBPB signal. GSE200997 had the same nominal direction, but did not provide statistically supportive replication under strict reference-mapped epithelial analysis (matched n=4; P=0.625); therefore, cross-cohort epithelial replication was not established.

## Primary epithelial validation

- GSE200997: matched n=4; median Δ(Tumor−Normal)=0.1975; Wilcoxon P=0.625; primary epithelial BH-FDR=0.625; Tumor higher.
- GSE132465: matched n=7; median Δ(Tumor−Normal)=0.607; Wilcoxon P=0.04688; primary epithelial BH-FDR=0.09375; Tumor higher.

## Cross-cohort epithelial replication

- Direction concordance: True; status: Direction concordant; statistical replication not established.

## Files

- `DINP_CRC_clean_CEBPB_cell_level.csv`: cell-level CEBPB metrics and annotation/QC status.
- `DINP_CRC_clean_GSE200997_annotation_audit.csv`: complete GSE200997 reference-mapping audit, including top-two labels, similarities, margins, and confidence status.
- `DINP_CRC_clean_CEBPB_patient_epithelial_cell_count_audit.csv`: patient × sample × condition QC and epithelial cell-count audit for both cohorts.
- `DINP_CRC_clean_GSE200997_epithelial_inclusion_audit.csv`: explicit GSE200997 paired inclusion/exclusion reasons explaining matched n.
- `DINP_CRC_clean_GSE200997_annotation_status_by_condition.csv`: Tumor/Normal condition-stratified confidence and epithelial retention audit.
- `DINP_CRC_clean_GSE200997_label_composition_by_condition.csv`: top1 and confident final label composition by condition.
- `DINP_CRC_clean_GSE200997_epithelial_threshold_sensitivity.csv`: strict/relaxed/very-strict reference threshold sensitivity.
- `DINP_CRC_clean_GSE200997_marker_gate_epithelial_sensitivity.csv`: non-transfer canonical epithelial-gate sensitivity.
- `DINP_CRC_clean_CEBPB_localization.csv`: cohort/condition/cell-type localization summary.
- `DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv`: patient/sample-level pseudobulk.
- `DINP_CRC_clean_CEBPB_tumor_normal_validation_all_celltypes.csv`: secondary six-cell-type paired validation with cohort-specific BH-FDR.
- `DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv`: prespecified epithelial validation with BH-FDR across the two cohorts.
- `DINP_CRC_clean_CEBPB_cross_cohort_replication.csv`: cross-cohort direction/replication audit.
- `DINP_CRC_clean_CEBPB_manifest.json`: source hashes, thresholds, and exact analysis settings.

The previous three-cohort CEBPB/CD36 module remains frozen and is not used as input here.
