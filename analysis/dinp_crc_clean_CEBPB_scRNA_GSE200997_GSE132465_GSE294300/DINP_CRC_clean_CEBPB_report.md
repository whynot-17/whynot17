# DINP–CRC clean CEBPB single-cell rerun with GSE294300

## Frozen design

CEBPB-only clean rerun using GSE200997, GSE132465, and GSE294300. GSE188711 was removed because it contained only three patients. No Left/Right/sidedness metadata, CD36, or old multicohort outputs were used.

- GSE132465 uses source curated broad cell-type labels and raw UMI counts.
- GSE200997 and GSE294300 use raw UMI counts and transparent nearest-centroid mapping to the GSE132465 curated labels; CEBPB is excluded from annotation features.
- Primary test: matched patient-level epithelial pseudobulk, Tumor versus Normal, paired Wilcoxon, minimum 50 confident epithelial cells per sample.
- Primary epithelial BH-FDR is recalculated across all three cohort tests; secondary cell-type FDR is calculated within cohort across six prespecified compartments.

## GSE294300 source and audit

- Loaded cells: 214,854; QC-passing: 214,854.
- GSE294300 is a paired 18-patient, 36-sample tumor/adjacent-normal 10x study; the deposited cell-batch barcodes were matched to raw 10x barcode cores after terminal suffix normalization, and raw matrices were restricted to that deposited cell subset. The official cell-batch URL and per-sample source URL template are recorded in the manifest.
- Cell types in GSE294300 are reference-mapped, not deposited curated labels; the mapping audit and condition-stratified label composition are exported separately.

## Annotation-status audit

- GSE200997 Tumor: QC=31573; Confident=92.37%; Low-confidence=5.35%; Ambiguous=2.29%; confident epithelial=20.46%.
- GSE200997 Normal: QC=18273; Confident=92.17%; Low-confidence=5.08%; Ambiguous=2.74%; confident epithelial=12.55%.
- GSE294300 Tumor: QC=106715; Confident=91.87%; Low-confidence=5.79%; Ambiguous=2.33%; confident epithelial=9.57%.
- GSE294300 Normal: QC=108139; Confident=93.93%; Low-confidence=4.10%; Ambiguous=1.97%; confident epithelial=5.81%.

## Patient-level epithelial inclusion

- GSE200997: 4 matched patient pairs pass the 50-cell rule.
- GSE132465: 7 matched patient pairs pass the 50-cell rule.
- GSE294300: 10 matched patient pairs pass the 50-cell rule.

## Threshold sensitivity

- GSE200997 strict: matched n=4; median Δ=0.1975; P=0.625; Tumor higher.
- GSE200997 relaxed: matched n=4; median Δ=0.2047; P=0.625; Tumor higher.
- GSE200997 very_strict: matched n=3; median Δ=-0.008762; P=1; Normal higher.
- GSE294300 strict: matched n=10; median Δ=0.6655; P=0.03711; Tumor higher.
- GSE294300 relaxed: matched n=10; median Δ=0.6336; P=0.03711; Tumor higher.
- GSE294300 very_strict: matched n=10; median Δ=0.7002; P=0.03711; Tumor higher.

## Primary epithelial validation

- GSE200997: matched n=4; median Δ=0.1975; P=0.625; BH-FDR=0.625; Tumor higher.
- GSE132465: matched n=7; median Δ=0.607; P=0.04688; BH-FDR=0.07031; Tumor higher.
- GSE294300: matched n=10; median Δ=0.6655; P=0.03711; BH-FDR=0.07031; Tumor higher.

## Cross-cohort epithelial audit

- Direction concordance across the three cohorts: True; status: Direction concordant; statistical replication not established.

## Interpretation

Adding GSE294300 increases the paired cohort evidence base without changing the frozen annotation or 50-cell rules. The result should be interpreted as cross-cohort epithelial replication only if direction is concordant and the patient-level evidence supports it after the three-cohort primary BH correction.

## Key files

- `DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv`: primary three-cohort epithelial validation.
- `DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv`: patient/sample/cell-type pseudobulk.
- `DINP_CRC_clean_GSE294300_annotation_audit.csv`: all GSE294300 mapped cells and confidence metrics.
- `DINP_CRC_clean_GSE294300_epithelial_inclusion_audit.csv`: explicit paired inclusion/exclusion audit.
- `DINP_CRC_clean_GSE294300_annotation_status_by_condition.csv`: Tumor/Normal annotation-status audit.
- `DINP_CRC_clean_CEBPB_manifest.json`: source and parameter manifest.
