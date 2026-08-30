# Primary Robustness Analysis Plan v1.0

## Status and scope

This plan is frozen after the locked primary package was produced and before any robustness model is fitted, read, or interpreted. It is a uniform post-primary stress test, not a candidate-selection procedure. Every applicable analysis is run across the complete pre-specified 29 x 14 = 406 pair grid. No CTD analysis, literature search, mechanism analysis, candidate prioritization, or new outcome/exposure selection is in scope.

Primary inference remains the existing pooled exposure-by-female interaction with its fixed 406-family BH correction. Robustness outputs are secondary diagnostic evidence and do not replace, refit, or revise the primary BH family.

## Fixed common model

All robustness models retain the locked primary population, outcome definition, exposure log2 coding, test-specific NHANES weight, cycle-offset strata/PSU, age (linear and quadratic), race/ethnicity, PIR, smoking, cycle fixed effects, pooled exposure-by-female interaction estimand, and Taylor-style stratified PSU sandwich variance. Male/female stratified estimates remain descriptive only.

## Uniform sensitivity analyses

### 1. Urinary-creatinine adjustment

For every pair with a urine-matrix exposure, refit the fixed pooled primary model after adding `log2(URXUCR)` and `female x log2(URXUCR)`. The latter permits sex-specific urine-dilution adjustment while retaining the same exposure-by-female estimand. This sensitivity is not applicable to non-urine exposures. It is performed for all applicable urine pairs, not selected primary findings.

### 2. Leave-one-contributing-cycle-out (LOCO)

For every pair, refit the fixed primary model after omitting each contributing cycle in turn. Record the interaction estimate, standard error, P value, analytic N, cases, and fit status for every leave-one-cycle-out refit. The pair-level LOCO diagnostic records the minimum and maximum interaction coefficient, whether any successful refit reverses the primary sign, and the number of successful refits.

### 3. Cycle heterogeneity

For every pair with at least three contributing cycles, fit a pooled model that adds exposure-by-female-by-cycle interaction deviations. Test the joint null that all cycle-specific deviations are zero using a design-based Wald test with the same degrees-of-freedom convention. Fewer than three contributing cycles are recorded as not applicable; this result is a heterogeneity diagnostic, not a new discovery family.

### 4. Upper-tail and influence stability

For every pair, run two pre-defined exposure perturbations, each within exposure cycle and before merging to outcomes:

- Winsorize log2 exposure at its contributing-cycle 1st and 99th percentiles.
- Exclude participants in the contributing-cycle upper 1% of log2 exposure.

Both use the unchanged fixed pooled model and are recorded for every pair.

### 5. Assay/LOD stability

For every pair whose frozen source data include a detection-limit flag, refit after restricting to observations flagged above the detection limit. For tests without a usable flag, record `not_applicable`; no proxy LOD definition is created. The output records the retained sample size and case count.

## Reporting and interpretation firewall

The robustness package must contain a row for every pair and each applicable analysis, including failures and non-applicability reasons. It must report coefficients and design-based uncertainty but must not run a new FDR correction, composite robustness score, candidate rank, or biological interpretation. Any later synthesis must display the locked primary FDR result separately from these diagnostics.
