# CNP0004138 Qin et al. 2023: paired scRNA and spatial reanalysis

This repository contains an exploratory reuse analysis of the public processed single-cell RNA-seq and spatial data for Qin et al., *Cell Reports Medicine* (2023). The Chinese report documents cohort accounting, patient-level paired tests, candidate epithelial/macrophage/CAF programs, spatial summaries, and interpretation limits.

## Main findings

- Across 27 paired patients, epithelial-like cell fraction did not show a detectable systematic change after chemotherapy. Several epithelial program scores decreased; these marker-defined states are not CNV-confirmed malignant cells.
- CAF-like fraction increased. Candidate CAF programs changed in different directions, including increases in PI16/SLIT2/ADAM28 and decreases in BMP4/FAP/MMP1.
- Patient-level post-treatment co-abundance linked epithelial EMT and CAF_FAP scores. This is a between-patient association and does not establish direct cell interaction.
- Spatial summaries are descriptive for two CR and two NR sections. They do not support patient-level response-group inference.
- The public barcode metadata and files supplied for this analysis do not include the XFZ-to-CR/PR/NR map. NR-specific scRNA tests could not be run.

## Reproduce

The processed expression matrices are not included because of their size. Obtain the study data from [CNSA project CNP0004138](https://db.cngb.org/search/project/CNP0004138/) and place them under raw/scRNA/ and raw/spatial/ in the repository root, retaining the filenames expected by src/analyze_reuse.py.

Install Python 3.12 or later and dependencies from requirements.txt, then run:

    python -m pip install -r requirements.txt
    python src/analyze_reuse.py
    python src/summarize_outputs.py

The scripts write derived files to results/, figures/, and report/. The analysis uses patient-stage samples as inputs to paired tests; individual cells are not treated as independent replicates.

## Repository contents

- analysis_report_zh.md: full Chinese report.
- results/: aggregate paired tests, co-abundance tables, spatial permutation summaries, and cell-count audit. Patient-by-sample and cell-by-cell tables are excluded.
- figures/: the three analysis figures.
- src/: analysis scripts with repository-relative paths.
- data_manifest.md: study/data provenance and cohort accounting.

## Interpretation limits

Cell states are exploratory marker-defined programs; the analysis does not reproduce the authors' Seurat clustering, inferCNV, NMF, or image deconvolution. The processed matrix retained 378,344 barcodes after the stated basic QC thresholds, whereas the article reports 333,285 high-quality cells. Spatial bins are approximately 25 μm and the four sections provide local descriptive evidence only. See the report before interpreting any result.
