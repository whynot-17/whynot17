# Step 10B-P — KoNEHS population-replacement feasibility audit

Generated: `2026-08-28T11:51:55Z`  \
Status: **analysis_feasible_by_precedent_but_access_controlled**  \
Decision: **high_priority_conditional_not_access_confirmed**

## Executive result

The audit supports KoNEHS as a **high-priority conditional population-replacement candidate**, not as an already accessible replication dataset. The audited public evidence gives a conservative content-level floor of **17/29 frozen tests** (12 exact analyte-level public matches and 5 family/matrix-qualified matches). This is not exact variable-level confirmation.

The audit confirms **29/29 tests were carried into the crosswalk**, but it confirms **0/29 exact public variable names**, **0/29 same-person exposure–T2D extracts**, and **0/29 exact weight/design-variable sets**. No exposure–T2D association model was run.

## Why KoNEHS is promising

KoNEHS is a national environmental-health biomonitoring program. The cycle-4 exposure overview documents a nationally sampled survey, environmental chemicals/metabolites, blood and urine biospecimens, urinary creatinine, laboratory QC, and design/nonresponse/post-stratification weights analyzed with a multistage survey procedure. Published cycle-2 and cycle-4 papers demonstrate that environmental biomarkers, diabetes-related outcomes, covariates, and survey analysis can be linked by authorized analysts.

The source-verified exposure floor includes three PFAS, BPA, six cycle-4 phthalate metabolites, cycle-3 MCOP exposure evidence, four grouped PAH-family matches from cycle 2, and conservative public metal evidence. The full row-level rationale is in `konehs_29_test_crosswalk.csv`; cycle-by-cycle status is in `konehs_cycle_readiness.csv`.

## T2D and design boundary

Published KoNEHS diabetes analyses establish operational feasibility, but not one universally frozen T2D definition. The cycle-2 PAH precedent used physician diagnosis and diabetes medication/insulin and explicitly did not distinguish type 1 from type 2. Therefore the T2D outcome remains **conditional** until the current data dictionary or authorized extract confirms diagnosis, medication, HbA1c/glucose fields, type-1 handling, missingness, and the exact analytic rule.

The same boundary applies to core covariates, urinary creatinine, survey weights, strata, and PSU identifiers: methodology is documented, but exact variable names and cycle-specific construction are pending.

## Access status

The publicly available sources audited here do not provide a direct, unrestricted individual-level joint exposure–T2D extract. The published data-availability boundary directs requests to the relevant Korean environmental-health authority. Thus KoNEHS is **analysis-feasible by precedent but access-controlled**. It should be promoted to primary population replacement only after an approved data request and exact data-dictionary crosswalk.

## Frozen exclusions

- No KoNEHS association result was computed.
- Published PFAS–diabetes, PAH–diabetes, and exposure papers are feasibility precedents only; their estimates are not this project's replication.
- No candidate was selected because it had a published KoNEHS result.
- The 29-test family was read from the existing frozen Step 4 file and was not altered.

## Gate for promotion

KoNEHS may become a primary epidemiologic replacement only if an authorized extract confirms: exact analyte and matrix; same-person exposure plus diabetes outcome; reproducible T2D coding; age/sex/BMI/smoking/alcohol/SES and other prespecified covariates; survey weights and design variables; and adequate access/permission for reproducible analysis.

## Files

- `konehs_29_test_crosswalk.csv`
- `konehs_cycle_readiness.csv`
- `konehs_outcome_covariate_design_audit.csv`
- `KONEHS_SOURCE_SNAPSHOT.json`
- `KONEHS_QC_SUMMARY.json`
- `STEP10B_P_KONEHS_MANIFEST.json`
