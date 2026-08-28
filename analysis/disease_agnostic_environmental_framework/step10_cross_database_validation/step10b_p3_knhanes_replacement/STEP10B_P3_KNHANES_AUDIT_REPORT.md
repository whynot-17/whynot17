# Step 10B-P3 — KNHANES executable population-replacement audit

Generated: `2026-08-28T12:29:08Z`  \
Status: **catalog_access_confirmed_exact_29_test_crosswalk_not_confirmed**  \
Decision: **not_currently_actionable_for_frozen_29_exact_replication**

## Executive result

KNHANES has a real official raw-data access route, but the present audit does **not** establish an exact match to any of the 29 frozen biomarker tests. The crosswalk contains **29/29 tests**, with **0/29 exact public matches**, **2/29 related-family blood-matrix mismatches**, and **27/29 not confirmed in audited public sources**.

This is not a failed data-access test. The official KDCA portal and Korea public-data record expose a raw-data catalogue and a registration/consent workflow. The metadata API returned a catalogue with `105` records (the API-reported total was `105` when available). However, this audit did not enter personal information or download a raw file, so file-level access is recorded as `public_registration_path_confirmed_file_download_not_completed`.

## Exact crosswalk result

The audited KNHANES literature documents blood lead, blood cadmium and blood mercury, and urinary total arsenic in specific 2008–2009 analyses. These are valuable environmental-health precedents, but the frozen panel contains serum PFAS and urine phthalate/PAH/metal tests. Blood lead/cadmium were therefore retained as explicit **matrix mismatches**, not upgraded to exact matches. KoNEHS phthalate/PFAS/MCOP evidence was not transferred to KNHANES.

No exact KNHANES evidence was confirmed for frozen serum PFDA/PFDeA, PFHxS, PFNA; urinary BPA, phthalate metabolites, PAH families; or the frozen urinary metals. The row-level evidence and notes are in `knhanes_29_test_crosswalk.csv`; annual metadata status is in `knhanes_year_readiness.csv`.

## Outcome, covariate and design feasibility

The official KDCA survey overview lists diabetes among chronic-disease/health-examination content and includes smoking, drinking, physical activity, obesity and other health domains. Published KNHANES studies demonstrate use of fasting glucose/glycemia, diabetes ascertainment, covariates and complex probability-sample methods. Exact current field names, type-1 handling, missingness, weights, strata and PSU variables remain codebook-level items to verify after download.

Thus the outcome/design side is **feasible by precedent**, but it cannot rescue a missing exact exposure matrix. The official raw catalogue also lists an environmental/indoor-air-related module in 2021; that catalogue metadata does not demonstrate coverage of the frozen 29 tests.

## What this means for population replacement

KNHANES is preferable to a purely hypothetical source because its official raw-data route is demonstrable and registration-based access is available. Nevertheless, under the frozen 29-test rule it is currently **not an exact biomarker population-replacement dataset**. It can be promoted only if the authorized download/codebook reveals additional environmental laboratory files containing one or more frozen analytes with the required matrix and a same-person diabetes frame.

## Frozen exclusions

- No association model was run.
- No published PFAS, arsenic, metal or diabetes effect estimate was adopted as our replication.
- No candidate was selected because it had a prior Korean paper.
- No KoNEHS measurement was treated as KNHANES evidence.
- No raw microdata were downloaded.

## Files

- `knhanes_29_test_crosswalk.csv`
- `knhanes_year_readiness.csv`
- `knhanes_outcome_covariate_design_audit.csv`
- `knhanes_raw_data_catalog.csv`
- `KNHANES_ACCESS_PROBE.json`
- `KNHANES_SOURCE_SNAPSHOT.json`
- `KNHANES_QC_SUMMARY.json`
- `STEP10B_P3_KNHANES_MANIFEST.json`
