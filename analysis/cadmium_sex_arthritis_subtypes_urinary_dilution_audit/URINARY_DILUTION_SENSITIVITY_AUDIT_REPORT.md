# Cd × Sex urinary-dilution sensitivity audit

## Scope

This bounded audit compares only the frozen URXUCD exposure for Any arthritis, osteoarthritis (OA), and gout. It uses existing NHANES local data, the established survey-weighted model implementation, and the same age, race/ethnicity, PIR, smoking, cycle fixed effects, weights, PSU/strata variance, and complete-case conventions. No new endpoint family, FDR, mechanism, literature, or figure analysis was performed.

The four estimable strategies are descriptive sensitivity comparisons; no new multiplicity adjustment was calculated. The historical six-test subtype q values are retained only as references.

## Specific-gravity availability

All 94 relevant local XPT files were scanned for specific-gravity variable names. No specific-gravity variable was identified (hits=0). The specific-gravity rows are therefore recorded as `not_applicable`; no proxy was substituted. Files with parser errors were retained in the manifest and did not contain usable SG metadata.

## Interaction comparison

| endpoint_id | strategy_id | status | interaction_beta | interaction_se | interaction_p | male_beta | female_beta | analytic_n | total_cases |
|---|---|---|---|---|---|---|---|---|---|
| any_arthritis | cd_cr_ratio | ok | 0.18659 | 0.050756 | 0.00035108 | -0.16539 | -0.057732 | 12719 | 3422 |
| any_arthritis | original_frozen | ok | 0.14128 | 0.037591 | 0.00026202 | -0.045315 | 0.044198 | 12723 | 3423 |
| any_arthritis | specific_gravity | not_applicable |  |  |  |  |  | 0 | 0 |
| any_arthritis | ucr_main_only | ok | 0.13624 | 0.037988 | 0.00048071 | -0.15916 | -0.051394 | 12719 | 3422 |
| any_arthritis | ucr_main_plus_sex_interaction | ok | 0.18947 | 0.049409 | 0.00019891 | -0.15916 | -0.051394 | 12719 | 3422 |
| gout | cd_cr_ratio | ok | 0.30476 | 0.12851 | 0.019755 | -0.097312 | -0.041168 | 9829 | 491 |
| gout | original_frozen | ok | 0.21447 | 0.12541 | 0.090546 | -0.050626 | 0.034088 | 9833 | 491 |
| gout | specific_gravity | not_applicable |  |  |  |  |  | 0 | 0 |
| gout | ucr_main_only | ok | 0.21219 | 0.12672 | 0.097365 | -0.098073 | -0.035182 | 9829 | 491 |
| gout | ucr_main_plus_sex_interaction | ok | 0.31182 | 0.1339 | 0.02202 | -0.098073 | -0.035182 | 9829 | 491 |
| osteoarthritis | cd_cr_ratio | ok | 0.16829 | 0.076204 | 0.029051 | -0.24653 | -0.097014 | 10606 | 1309 |
| osteoarthritis | original_frozen | ok | 0.15298 | 0.054509 | 0.0058183 | -0.12459 | 0.019999 | 10610 | 1310 |
| osteoarthritis | specific_gravity | not_applicable |  |  |  |  |  | 0 | 0 |
| osteoarthritis | ucr_main_only | ok | 0.14799 | 0.054826 | 0.0079206 | -0.24978 | -0.091685 | 10606 | 1309 |
| osteoarthritis | ucr_main_plus_sex_interaction | ok | 0.17474 | 0.07473 | 0.020975 | -0.24978 | -0.091685 | 10606 | 1309 |

## Pattern summary

| endpoint_id | n_estimable_strategies | all_estimable_interactions_positive | all_estimable_nominal_p_lt_0_05 | min_interaction_beta | max_interaction_beta | strategies_with_nominal_p_lt_0_05 |
|---|---|---|---|---|---|---|
| any_arthritis | 4 | True | True | 0.13624 | 0.18947 | cd_cr_ratio;original_frozen;ucr_main_only;ucr_main_plus_sex_interaction |
| osteoarthritis | 4 | True | True | 0.14799 | 0.17474 | cd_cr_ratio;original_frozen;ucr_main_only;ucr_main_plus_sex_interaction |
| gout | 4 | True | False | 0.21219 | 0.31182 | cd_cr_ratio;ucr_main_plus_sex_interaction |

## Interpretation

The formal quantity compared across strategies is the pooled Cd × female interaction. Positive interaction means the female slope is higher than the male slope under that model specification; it does not by itself establish positive female susceptibility, causality, or protection. Male and female slopes are descriptive and are not interpreted by the one-side-significance rule.

The earlier subtype report now labels Any arthritis, OA, and gout as `FDR-supported sex heterogeneity, but not a positive female-susceptibility pattern`, because their fixed six-endpoint q values are below 0.05 while both sex-specific Cd slopes are negative. This audit does not convert those inverse slopes into protective claims.

## Data and provenance

Run timestamp (UTC): 2026-08-30T14:29:55.060336+00:00
Input exposure cycles: 2003-2004;2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018
Output CSV: 01_urinary_dilution_strategy_comparison.csv
