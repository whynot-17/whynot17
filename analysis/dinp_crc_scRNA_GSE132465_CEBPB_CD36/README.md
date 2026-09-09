# CEBPB/CD36 CRC single-cell localization and patient-level pseudobulk validation

This analysis evaluates the two current candidates, **CEBPB** and **CD36**, in the public CRC single-cell series **GSE132465**. GEO reports 63,689 cells from 23 CRC patients, including 23 primary tumor samples and 10 matched normal mucosa samples. The analysis uses the GEO-provided curated cell-type annotations and does not infer a cell type from CEBPB/CD36 expression alone.

## Analysis design

- Extract CEBPB and CD36 from the public raw UMI matrix after exact cell-ID/order concordance checking against the GEO annotation table.
- Report target detection rate, raw UMI, CPM-normalized expression, and target-UMI fraction by sample class, curated cell type, and cell subtype.
- Aggregate counts as patient × tumor/normal class × cell type pseudobulk.
- Normalize aggregate target UMI by aggregate library UMI to pseudobulk CPM.
- Use the 10 matched tumor/normal patients and a two-sided paired Wilcoxon test on `log1p(pseudobulk CPM)`; BH correction is applied across the 12 gene × cell-type tests.
- Require at least 50 cells per sample × cell type for pseudobulk eligibility. Cell-level rows, cell counts, raw target UMI, and library UMI remain in the output for audit.

## Main findings

- Both CEBPB and CD36 localize most strongly to **Myeloids** in tumor samples by detection rate and mean CPM.
- CEBPB is nominally higher in tumor epithelial-cell pseudobulk, but the BH-adjusted result is not significant (FDR 0.117).
- CD36 is higher in normal B-cell pseudobulk after BH correction (FDR 0.039); stromal cells also trend normal-higher (FDR 0.078).
- These are expression-localization and validation results, not evidence of direct DINP binding or causal CRC activity.

## Files

- `GSE132465_CEBPB_CD36_report.md`: human-readable report.
- `GSE132465_CEBPB_CD36_analysis_manifest.json`: data hashes, parameters, sample counts, and QC.
- `GSE132465_CEBPB_CD36_cell_level_expression.csv`: target expression and library-size metrics for all 63,689 cells.
- `GSE132465_CEBPB_CD36_cell_type_localization.csv`: class × cell-type localization summary.
- `GSE132465_CEBPB_CD36_cell_subtype_localization.csv`: class × cell-type × subtype summary.
- `GSE132465_CEBPB_CD36_patient_celltype_pseudobulk.csv`: patient-level pseudobulk table.
- `GSE132465_CEBPB_CD36_paired_pseudobulk_validation.csv`: matched-patient tests and FDR.
- `scRNA_target_localization_dotplot.png`: target localization plot.
- `scRNA_patient_pseudobulk_paired_validation.png`: paired patient-level validation plot.
- `run_crc_singlecell_cebpb_cd36.py`: reproducible script.

The large raw GEO matrix is intentionally not stored in Git. Download the two source files listed in the manifest to `work/scRNA_GSE132465/`, then run the script from the repository root.
