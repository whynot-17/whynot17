# Step 10B-P3 — KNHANES executable population-replacement audit

## Objective

Audit whether the Korea National Health and Nutrition Examination Survey (KNHANES) can provide an independently executable population-replacement analysis for the **29 frozen human biomarker tests**. This stage freezes the data-access and analyte crosswalk only; it does not run exposure–diabetes association models.

## Frozen boundary

The 29-test family is read directly from `step04_testset_freeze/unique_biomarker_test_set.csv`. No T2D result, published Korean association estimate, GeneCards result, or candidate rank is used to alter the test family or upgrade a match.

## What is audited

For every frozen test, record:

1. exact analyte and specimen-matrix evidence in KNHANES;
2. relevant year/cycle evidence;
3. whether the official raw-data catalogue exposes a downloadable database record;
4. diabetes, fasting glucose, HbA1c, medication and insulin availability;
5. core covariates and survey design/weight fields;
6. the difference between registration-based raw-data access and a completed file download.

## Evidence classes

- `exact_public`: exact analyte and matrix documented in an audited KNHANES source.
- `family_public_matrix_mismatch`: related chemical family is documented, but the matrix or analyte is not the frozen test.
- `not_confirmed_public`: no exact or defensible family/matrix match was confirmed in this audit.

The audit deliberately does not treat KoNEHS measurements as KNHANES measurements. It also does not treat the existence of a Korean paper as proof that the same frozen biomarker is available in the public KNHANES file.

## Access boundary

The KDCA KNHANES portal exposes an official raw-data catalogue and a registration/consent workflow requiring user information. The catalogue/API availability is recorded as an access-path test. A raw file is **not** downloaded in this audit and no personal information is entered. File-level access remains `registration_required_and_not_completed_in_audit` until the user completes the official registration/download step.

## Promotion rule

KNHANES can be promoted to a primary population replacement only after an authorized download confirms exact analyte/matrix availability, same-person biomarker and diabetes fields, prespecified T2D coding, covariates, weight/strata/PSU variables, and a reproducible analysis frame.

## Outputs

- `knhanes_29_test_crosswalk.csv`
- `knhanes_year_readiness.csv`
- `knhanes_outcome_covariate_design_audit.csv`
- `knhanes_raw_data_catalog.csv`
- `KNHANES_ACCESS_PROBE.json`
- `KNHANES_SOURCE_SNAPSHOT.json`
- `KNHANES_QC_SUMMARY.json`
- `STEP10B_P3_KNHANES_AUDIT_REPORT.md`
- `STEP10B_P3_KNHANES_MANIFEST.json`
