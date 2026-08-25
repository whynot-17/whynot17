# MCOP–CRC submission supplement

This directory is the submission-facing Supplementary Information package for the latest MCOP–CRC manuscript.

## Deliverables

- `MCOP_CRC_Supplementary_Information.docx` — editable Supplementary Methods, Results notes, figure legends, table index and consistency checklist.
- `MCOP_CRC_Supplementary_Information.pdf` — A4 review PDF rendered from the Word content.
- `MCOP_CRC_Supplement_Tables.xlsx` — authoritative Tables S1–S9.
- `MCOP_CRC_Source_Data.xlsx` — panel-level source data for Supplementary Figures S1–S4.
- `figures/` — PDF, SVG and 300-dpi PNG exports.
- `scripts/` — deterministic figure, Word and audit builders.
- `QA/` — export audit, numeric consistency audit and SHA-256 file manifest.

## Scientific boundary

- The 267→87→15 actionability sequence is outcome-blinded.
- MCOP is a urinary biomarker for a DINP-related exposure axis; it is not a significant direct CTD molecular hit.
- MiNP, DINP parent and MCOP remain chemically distinct.
- The NHANES result is a cross-sectional association with prevalent CRC.
- CRC epithelial PPAR/nuclear-receptor remodeling is an independent disease-state observation; the exposure-to-state causal bridge is untested.
- Individual single cells are not treated as independent replicates; donor-level inference is primary.

## Build

Run `scripts/make_supplementary_figures.py`, then `scripts/build_supplement_docx.js`, and finally `scripts/final_submission_audit.py`. Large raw Census caches and H5AD/XPT files are intentionally excluded from version control.

