# GSE178341 primary validation: PGE2 → PTGER4 axis

This module applies the frozen PGE2 synthesis source → PTGER4 receiver
robustness protocol to GSE178341 as the primary validation cohort. The prior
GSE144735 discovery and robustness modules remain unchanged.

GSE178341 contributes 370,115 cells from 62 donors. The public processed GEO
H5 matrix is used because raw FASTQ access is controlled through dbGaP. The
wrapper extracts only the PGE2 target genes and the fixed SPP1/C1QC/FCN1/
resident marker genes, joins the official c295 `clTopLevel`, `clMidwayPr`, and
`cl295v11SubFull` annotations, and writes the small target-gene cache to
`E:\mcop`. No H5 or other raw data are copied into Git.

The primary unit is donor × analysis_group pseudobulk with ≥20 cells and
threshold sensitivity at ≥10 and ≥30. Scores A/B/D exclude PTGES3; Score C is
the original six-gene sensitivity. Receiver evidence uses PTGER4 mean,
detection, expressing-cell mean, pseudobulk expression, and donor rank.

The official c295 annotation has one macrophage-like class (`cM02`). To make
the discovery labels comparable without inventing official subclusters, the
wrapper applies a transparent marker-defined transfer rule inside official
Macro/Mono cells. A cell is assigned SPP1-like, C1QC-like, FCN1-like, or
resident-like only when its gene-wise z-scored marker score is ≥0.15 and its
margin over the second score is ≥0.05; ambiguous cells remain
`other_myeloid`. The validation report keeps this distinction explicit.

## Run

```powershell
E:\chatgpt\sc_env\Scripts\python.exe `
  analysis/dinp_crc_pge2_ptger4_axis_gse178341_validation/run_pge2_ptger4_gse178341_validation.py
```

The report is `outputs/PGE2_PTGER4_GSE178341_VALIDATION_REPORT.md`. Tables and
PNG/PDF/SVG figures are under `outputs/`. The large input files and generated
cell-score cache remain on `E:\mcop`.

The results are validation evidence at the expression level. They do not
demonstrate extracellular PGE2 production, receptor binding, spatial
proximity, or functional communication.
