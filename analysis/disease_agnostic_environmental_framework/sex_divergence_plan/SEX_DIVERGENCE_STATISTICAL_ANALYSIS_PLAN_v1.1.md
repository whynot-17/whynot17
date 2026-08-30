# Sex-divergence statistical analysis plan v1.1

## Pre-association amendment

This is a pre-association amendment to v1.0. At the time of amendment, no 29 by 14 exposure--outcome association, sex-stratified association, or exposure-by-sex interaction result has been fitted, read, or interpreted. All v1.0 rules remain in force except the formal divergence estimand defined below.

The amendment corrects an estimand mismatch: separate male and female models allow nuisance-covariate coefficients to differ by sex, whereas the primary pooled interaction model constrains them to the shared primary parameterization. Therefore `beta_female_stratified - beta_male_stratified` is not algebraically required to equal the pooled interaction coefficient and is not required to have the same sign.

## Formal pairwise estimand

For each of the fixed 406 exposure--outcome pairs, fit the frozen pooled survey-weighted logistic model:

`logit Pr(outcome=1) = intercept + exposure + female + exposure:female + covariates + cycle fixed effects`.

Define the formal sex-divergence estimand as:

`d(e,j) = beta_exposure:female_pooled(e,j)`.

Under the frozen pooled-model parameterization, `d(e,j)` is the female exposure slope minus the male exposure slope. Its design-based standard error, confidence interval, Wald statistic, and two-sided P value are the sole formal evidence for heterogeneity for that pair.

The male and female survey-weighted stratified exposure coefficients and standard errors remain required outputs for interpretability. They are descriptive estimates, not the formal divergence estimand, and are not algebraically required to equal or share the sign of `d(e,j)` because nuisance-covariate coefficients may differ between stratified fits.

## Propagated divergence definitions

All formal magnitude summaries use the pooled interaction coefficient above, never a difference of separate-fit coefficients:

`DLD(e) = (1/14) sum_j |d(e,j)|`.

For organ system `s`, with pre-frozen within-system weights `w_j` that sum to 1:

`SD_s(e) = sum_(j in s) w_j |d(e,j)|`.

The equal-system overall organ-system divergence is:

`OSD(e) = (1/7) sum_s SD_s(e)`.

The signed system landing contrast is:

`SL_s(e) = sum_(j in s) w_j d(e,j)`.

Positive `d(e,j)` or `SL_s(e)` denotes a stronger positive exposure association in females under the frozen coding; negative values denote a stronger positive association in males. DLD, SD, OSD, and SL are effect-size magnitude summaries and descriptive rankings. They are not significance scores, do not include P values or FDR values, and must be presented alongside the fixed-family interaction FDR evidence.

## Rules retained from v1.0

The 29-test exposure set, 14-outcome set, 7-system ontology, sample rules, covariate set, survey design, outcome definitions, complete-case primary analysis, and Benjamini--Hochberg interaction family remain unchanged. The fixed FDR denominator remains 406; non-estimable pairs receive P=1 and remain in that denominator. No exposure, outcome, system, threshold, or model rule was selected from association results.
