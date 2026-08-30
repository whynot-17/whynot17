# Sex-divergence statistical analysis plan v1.0

## Status and scope lock

This plan is locked before any of the 29 by 14 exposure--outcome estimates, sex-stratified estimates, or exposure-by-sex interaction tests are fitted, read, or interpreted. It uses the 29-test frozen exposure set and the 14 PASS outcomes in Frozen Outcome Set v1.0. The four pre-frozen but unavailable outcomes are not substituted, imputed, or added to the primary analysis.

The primary analytic grid contains 406 prespecified exposure--outcome pairs (29 exposures x 14 outcomes). The same 406 pairs define the formal interaction family.

## Frozen inputs

| Input | Locked value |
|---|---|
| Exposure set | 29 tests in `step04_testset_freeze/unique_biomarker_test_set.csv` (SHA-256 `6b284b23d69f74991fc28978c85c1e39d10f096e2fd4d42290046df8224c3a9e`) |
| Outcome set | 14 `selection_for_followup=True` rows in `outcome_inventory/frozen_outcome_set_v1.csv` (SHA-256 `3d2ac6f399ee02e9954300e1bbbc20d62a1afb265ca531a13399322194c89786`) |
| Organ systems | 7 systems and within-system weights in `frozen_outcome_organ_system_ontology_v1.csv` |
| Target population | Participants age >=20 years, with sex from DEMO.RIAGENDR and the outcome-specific eligibility/control definition already frozen in Outcome Set v1.0 |

## Primary pairwise estimands

For every exposure--outcome pair, fit the pre-specified survey-weighted logistic model separately in males and females. The exposure coefficient and its design-based standard error are reported as `beta_male`, `se_male`, `beta_female`, and `se_female`. Exposure is transformed, censored/handled at its detection limit, and weighted only according to its existing frozen test-level specification; no outcome result may change those choices.

The formal sex interaction is estimated in the corresponding pooled-sex survey-weighted logistic model:

`logit Pr(outcome=1) = intercept + exposure + female + exposure:female + covariates + cycle fixed effects`.

The coefficient of `exposure:female`, its design-based standard error, Wald statistic, two-sided P value, and confidence interval are the sole primary test of sex heterogeneity. The reported signed contrast is `beta_female - beta_male`; it must agree in sign with the interaction parameter under the common exposure coding. Separate male/female fits are descriptive estimates, not a substitute for the interaction test.

All models use the appropriate frozen NHANES examination/laboratory weight for that test, the matching strata and PSU variables, and the test's available cycles. In pooled analyses, each two-year weight is divided by the number of contributing two-year cycles for that test. The outcome is defined from outcome sources independently of the test laboratory file. The analytic sample is the intersection of a defined binary outcome, measured exposure, positive appropriate survey weight, known sex, and complete pre-specified covariates.

## Covariates and missing data

The primary covariate set is age (linear and quadratic terms), race/ethnicity, poverty-income ratio, and smoking category, plus cycle fixed effects. BMI is deliberately not included in the universal primary adjustment set because it defines obesity and can be downstream of environmental exposures; it is reserved for a clearly labeled sensitivity analysis if run later. No data-driven covariate selection, outcome-specific covariate changes, imputation, or association-informed eligibility changes are permitted in the primary analysis. Primary models use complete cases for this fixed set and report their sex-specific analytic N.

## Multiplicity and estimability

The 406 interaction P values form one pre-specified family. Benjamini--Hochberg FDR correction is applied with a fixed denominator of 406. If a prespecified pair is not estimable after the locked sample rules, it remains in the denominator and is assigned a conservative P value of 1 for the FDR calculation; its non-estimability reason is reported. There is no selective FDR family based on nominal sex-stratified, main-effect, or organ-system findings.

Main-effect and sex-stratified P values are descriptive unless a separately declared secondary family is created after the primary interaction analysis has been reported. No exposure, outcome, system, or sex is removed based on observed effect size or P value.

## Disease-level and organ-system divergence

Let `d(e,j) = beta_female(e,j) - beta_male(e,j)` for exposure `e` and outcome `j` on the shared log-odds scale.

Disease-level divergence is reported without outcome-count reweighting as:

`DLD(e) = (1/14) sum_j |d(e,j)|`.

For each organ system `s`, the system divergence is the weighted mean of its member outcomes:

`SD_s(e) = sum_(j in s) w_j |d(e,j)|`, where the pre-frozen weights `w_j` sum to 1 inside each system.

The overall organ-system divergence gives every system equal influence:

`OSD(e) = (1/7) sum_s SD_s(e)`.

The signed system landing contrast is also reported:

`SL_s(e) = sum_(j in s) w_j d(e,j)`.

`SL_s(e) > 0` denotes stronger positive exposure association in females under the frozen coding; `SL_s(e) < 0` denotes stronger positive association in males. These divergence summaries are descriptive rankings and do not create additional confirmatory P values. Their interpretation must retain the outcome-level estimates and the 406-test interaction FDR results.

## Reporting firewall

Before results are opened, the analysis must write an analysis manifest containing input hashes, 29 x 14 grid, sex-specific analytic N, outcome availability, model convergence/estimability, and the fixed FDR denominator. It must not read CRC/T2D association output, prior Step 5 exposure associations, or candidate-chemical results. Any future change to this plan requires a new versioned plan and a separate post-lock exploratory label.
