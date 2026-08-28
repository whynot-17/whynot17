# Step 10B-P — KoNEHS population-replacement feasibility audit

## Objective

Audit whether the Korean National Environmental Health Survey (KoNEHS) can serve as an independent population-replacement source for the **29 frozen human biomarker tests**. This stage is a feasibility and data-dictionary audit only. It does not run exposure–diabetes models and does not import effect estimates from published KoNEHS studies.

## Frozen boundary

The 29-test family is read from `step04_testset_freeze/unique_biomarker_test_set.csv` and is not reselected. No T2D result, GeneCards result, CTD disease-gene result, or published KoNEHS association estimate is used to alter the crosswalk.

## Audit domains

For each frozen test, record separately:

1. exact analyte or family evidence in an audited KoNEHS cycle;
2. whether the exact public variable name and matrix are available;
3. whether an individual-level exposure–diabetes joint extract is publicly available;
4. whether survey-weight and design-variable names are confirmed;
5. whether the operational diabetes definition is reproducible as T2D rather than broad diabetes.

KoNEHS cycles are represented as K1 (2009–2011), K2 (2012–2014), K3 (2015–2017), and K4 (2018–2020). A public-paper or official-overview match is a **source-verified feasibility floor**, not a claim that every test is measured in every cycle.

## Decision rule

- `analysis-feasible_by_precedent_but_access_controlled`: a published or official source documents the relevant exposure/outcome/design ingredients, but exact variable crosswalk and individual-data access are not confirmed.
- `high-priority_conditional`: content-level overlap is substantial and a data request is justified.
- Promotion to `primary epidemiologic replacement` requires all of: exact analyte and matrix; same-person exposure plus diabetes outcome; reproducible outcome coding; core covariates; survey design/weights; and legally available individual-level data.

## Explicit exclusions

- Do not treat published PFAS–diabetes or PAH–diabetes effect estimates as this project's replication.
- Do not use PFHxS, PAH, MCOP, or any other candidate's publication history to select a result after seeing the crosswalk.
- Do not download controlled KoNEHS microdata.
- Do not infer an exact urine metal variable from a source that only documents a blood measurement.

## Outputs

- `konehs_29_test_crosswalk.csv`
- `konehs_cycle_readiness.csv`
- `konehs_outcome_covariate_design_audit.csv`
- `KONEHS_SOURCE_SNAPSHOT.json`
- `KONEHS_QC_SUMMARY.json`
- `STEP10B_P_KONEHS_AUDIT_REPORT.md`
- `STEP10B_P_KONEHS_MANIFEST.json`
