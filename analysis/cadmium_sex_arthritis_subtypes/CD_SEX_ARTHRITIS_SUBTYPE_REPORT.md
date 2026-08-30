# Cd × Sex × Arthritis Subtype Analysis

## Scope and stop rule

This focused analysis uses the frozen `URXUCD` urinary cadmium exposure and six prespecified outcomes: any arthritis (reference), osteoarthritis, rheumatoid arthritis, psoriatic arthritis, other arthritis, and gout. It stops after variable audit, sample counts, survey-weighted interaction models, six-test BH-FDR, cycle consistency, and prespecified LOCO diagnostics. No CTD, gene, enrichment, PPI, GTEx, single-cell, mediation, hormone, literature, or figure analysis was run.

## Outcome variable mapping

`MCQ160A` is used for any arthritis and as the case gate for arthritis subtypes. The cycle-specific subtype variable is `MCQ190` for 1999–2008, `MCQ191` for 2009–2010, and `MCQ195` for 2011–2018. PsA is not available in the pre-2009 mapping and is not treated as a negative in those cycles. `MCQ160N` gout is available from 2007–2008 onward. CDC codebook URLs are retained in `01_arthritis_subtype_variable_audit.csv`.

## Model

The main model is a survey-weighted logistic regression with `log2(URXUCD)`, female, their interaction, age (centered linear and quadratic), race/ethnicity, PIR, smoking, cycle fixed effects, and urinary creatinine adjustment as `log2(URXUCR)` plus its female interaction. The formal sex difference is the pooled `URXUCD × female` coefficient. Sex-stratified coefficients are descriptive. Weights, cycle-offset strata/PSU identifiers, complete-case handling, and Taylor-style stratified PSU sandwich variance follow the existing NHANES pipeline. The six interaction P values use a fixed BH denominator of 6; the historical any-arthritis fixed-406 q is retained only as a reference.

## Analytic sample and interaction results

| endpoint_id | outcome_name | endpoint_class | interaction_beta | interaction_se | interaction_p | subtype_family_bh_q | interaction_fdr_supported | male_beta | female_beta | analytic_n | total_cases | male_cases | female_cases | interpretation_class | directional_pattern |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| any_arthritis | Any arthritis | reference | 0.1895 | 0.04941 | 0.0001989 | 0.001193 | True | -0.1592 | -0.05139 | 12719 | 3422 | 1451 | 1971 | FDR-supported sex heterogeneity, but not a positive female-susceptibility pattern | female_enhanced_relative_slope_both_inverse |
| gout | Gout | primary_related_endpoint | 0.3118 | 0.1339 | 0.02202 | 0.04404 | True | -0.09807 | -0.03518 | 9829 | 491 | 350 | 141 | FDR-supported sex heterogeneity, but not a positive female-susceptibility pattern | female_enhanced_relative_slope_both_inverse |
| osteoarthritis | Osteoarthritis (OA) | primary_subtype | 0.1747 | 0.07473 | 0.02098 | 0.04404 | True | -0.2498 | -0.09168 | 10606 | 1309 | 494 | 815 | FDR-supported sex heterogeneity, but not a positive female-susceptibility pattern | female_enhanced_relative_slope_both_inverse |
| other_arthritis | Other arthritis | primary_subtype | 0.09737 | 0.1443 | 0.5012 | 0.5012 | False | -0.04168 | -0.01766 | 9692 | 395 | 174 | 221 | No convincing sex heterogeneity | no_fdr_supported_directional_difference |
| psoriatic_arthritis | Psoriatic arthritis (PsA) | exploratory | 0.8549 | 0.6277 | 0.1771 | 0.2126 | False | -0.177 | 1.084 | 5989 | 22 | 15 | 7 | Unstable / non-estimable | no_fdr_supported_directional_difference |
| rheumatoid_arthritis | Rheumatoid arthritis (RA) | primary_subtype | 0.1749 | 0.1046 | 0.09722 | 0.1458 | False | -0.1057 | -0.07756 | 9950 | 653 | 287 | 366 | No convincing sex heterogeneity | no_fdr_supported_directional_difference |

Interpretation class counts: {'FDR-supported sex heterogeneity, but not a positive female-susceptibility pattern': 3, 'No convincing sex heterogeneity': 2, 'Unstable / non-estimable': 1}.

Among the arthritis subtypes, OA is the only subtype whose interaction remains below the fixed six-endpoint BH threshold; RA and Other arthritis do not, and PsA is underpowered. The reference any-arthritis interaction and the related gout interaction also pass the six-test threshold. For these FDR-supported positive interactions, both sex-stratified Cd slopes are negative and the female slope is less negative than the male slope, so the result is a female-enhanced relative slope rather than a positive female-susceptibility association under the prespecified rule.

Case-count cautions are explicit in `02_subtype_sample_counts.csv`; PsA is exploratory and is not eligible for strong interpretation if either sex has fewer than 50 cases or total cases are fewer than 100. No subtype is called protective. Directional wording is limited to susceptible, attenuated/weaker, or no convincing heterogeneity.

## Cycle and LOCO diagnostics

`05_cycle_consistency.csv` reports each contributing cycle's interaction estimate and whether its direction matches the pooled estimate. `06_loco_results.csv` reports each leave-one-cycle-out refit when the prespecified minimum of 100 total cases and 50 cases in each sex was met; otherwise it records `LOCO not performed because of inadequate subtype case counts`.

## Interpretation boundary

The analysis localizes the existing Cd × sex arthritis signal to the prespecified outcome definitions only. It does not establish mechanism, causality, temporality, or clinical protection/resilience.

Run timestamp (UTC): 2026-08-30T14:26:32.136262+00:00
