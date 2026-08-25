# MCOP–CRC Supplement Reproducibility Package

This directory contains the pre-submission revalidation package for the frozen MCOP–CRC analysis.

## Contents

- `scripts/` — rerun, audit, and table-building scripts.
- `results/` — rerun outputs organized by analysis phase.
- `supplementary/` — Tables S1–S9 as CSV and one Excel workbook.
- `logs/` — audit results, input/script inventory, hashes, and environment metadata.
- `reproducibility_report.md` — final audit report and reporting boundaries.

## Frozen analysis boundaries

- NHANES MCOP primary analysis uses the seven phthalate cycles from 2005–06 through 2017–18, analyte-specific phthalate subsample weights divided by 7, complex-survey inference, and the frozen covariate model.
- The 267-chemical actionability screen is outcome-blinded for eligibility; CRC odds ratios are used only after the human-testable axes are defined.
- MCOP is treated as the biomarker for a DINP-related exposure axis. The molecular screen nominated the DINP/MiNP axis; it is not represented as a direct CTD mechanism for MCOP.
- Phase2G final local expression results use the verified official H5AD source. The live Census staged cache is retained as an execution log and is not represented as complete when its donor run is incomplete.
- NHANES results remain cross-sectional associations. No causal or mediation claim is introduced by this package.

## Re-running

Use the isolated Windows environment described in `logs/analysis_environment.txt` and run:

```powershell
python scripts/build_supplement_tables.py
python scripts/final_reproducibility_audit.py
```

The raw NHANES XPT data, the external H5AD, and large Census/raw-expression caches are intentionally not duplicated in this package or committed to GitHub. Their paths, sizes, and hashes are recorded in `logs/data_inventory.csv`.
