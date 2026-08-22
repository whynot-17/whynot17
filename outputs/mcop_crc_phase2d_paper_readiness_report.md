# MCOP–CRC Phase 2D：paper-readiness audit

## Material Passport

- Origin Skill: academic-research-suite + experiment-agent
- Origin Mode: run / validation audit
- Origin Date: 2026-08-22
- Verification Status: ANALYZED — local NHANES rerun completed; no WHI data used
- Scope: female/male effect, sex interaction, diagnosis timing, co-exposure specificity, paper figures

## Frozen scope

This Phase 2D run does not screen another chemical and does not run TCGA/PPI/GO. All regression models reuse the established NHANES survey-logistic implementation, cancer-free-control definition, phthalate weights, pooled strata/PSU identifiers and Model 2 covariates: age, BMI, PIR, urinary creatinine, sex, race and smoking unless sex is constant within a sex-specific analysis.

## Gate results

- Female point estimate >1: **PASS**; OR=1.137, 95% CI 0.9462-1.367, P=0.1684, N=5065, CRC=33. This is directionally positive but not a close point-estimate replication of the pooled OR≈1.25.
- Male point estimate: OR=1.352, 95% CI 1.091-1.676, P=0.006381, N=4871, CRC=37.
- MCOP × sex interaction: OR ratio male/female=1.134, 95% CI 0.8845-1.454, P=0.3183.
- Interpretation: the female estimate remains positive but is weaker than the pooled estimate and its CI includes the null; the formal interaction does not support a statistically clear sex difference.
- Recent-diagnosis exclusions retain positive MCOP direction: **PASS**; ORs=1.241; 1.263; 1.265.
- MCOP co-exposure models retain positive MCOP direction: **PASS**; MCOP ORs=1.222; 1.219; 1.202; 1.247; 1.235.

These are association audits. A positive result does not establish that MCOP or DINP caused CRC; the NHANES design still permits reverse causation and survivor bias.

## Diagnosis timing / reverse-causation audit

The case-level timing export contains exam age, CRC diagnosis age and years since diagnosis. The model rows exclude cases with known diagnosis-to-exam interval below 1, 2 or 5 years; cases with missing diagnosis age are retained and counted explicitly as missing timing.

| Analysis | N | CRC cases | Excluded known-timing CRC | OR | 95% CI | P |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Primary_all_cases | 9936 | 70 | 0 | 1.246 | 1.078-1.44 | 0.003311 |
| Exclude_diagnosis_lt_1y | 9928 | 62 | 13 | 1.241 | 1.062-1.45 | 0.006912 |
| Exclude_diagnosis_lt_2y | 9921 | 55 | 23 | 1.263 | 1.071-1.489 | 0.005858 |
| Exclude_diagnosis_lt_5y | 9910 | 44 | 43 | 1.265 | 1.059-1.51 | 0.01013 |

The ≤5-year versus >5-year case groups are descriptive only and are not treated as independent etiologic tests.

## Phthalate specificity

Each co-exposure model is restricted to cycles with MCOP and that co-exposure available. The burden is the mean z-score of nine non-MCOP urinary phthalate metabolites, requiring at least two measured metabolites. Same-complete-case MCOP-only comparator rows are included in the CSV.

| Model | Cycles | N | CRC | MCOP OR | 95% CI | MCOP P | Co-exposure OR | Co-exposure P | MCOP BH-FDR |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MCOP_plus_MEHHP | 2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 | 9936 | 70 | 1.222 | 1.057-1.413 | 0.007071 | 1.184 | 0.1057 | 0.009799 |
| MCOP_plus_MEOHP | 2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 | 9936 | 70 | 1.219 | 1.055-1.409 | 0.007839 | 1.208 | 0.05468 | 0.009799 |
| MCOP_plus_MECPP | 2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 | 9936 | 70 | 1.202 | 1.035-1.395 | 0.0163 | 1.229 | 0.08236 | 0.0163 |
| MCOP_plus_MBzP | 2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 | 9936 | 70 | 1.247 | 1.082-1.439 | 0.002696 | 0.9658 | 0.7546 | 0.009799 |
| MCOP_plus_PhthalateBurden_excl_MCOP | 2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 | 9936 | 70 | 1.235 | 1.069-1.428 | 0.004608 | 1.283 | 0.475 | 0.009799 |

## Figures

- Figure 1: workflow from CTD × GeneCards to NHANES MCOP audit and future WHI replication.
- Figure 2: primary model and stability/reverse-causation sensitivities.
- Figure 3: seven cycle-specific estimates and pooled estimate.
- Figure 4: restricted cubic spline.
- Figure 5: female/male estimates with interaction audit.
- Figure 6: MCOP specificity under co-exposure adjustment.

WHI is shown only as a future prospective replication line; no WHI data were accessed or analyzed in this run.

## Statistical fallacy scan (11/11)

1. Simpson's paradox — sex-specific estimates are compared with pooled and interaction results; direction is reported rather than hidden.
2. Ecological fallacy — not applicable; the analysis unit is the individual NHANES participant.
3. Berkson's paradox — caution: the analytic population is conditioned on observed cancer outcome and phthalate subsample availability.
4. Collider bias — caution: conditioning on cancer ascertainment/survival and complete covariates may induce selection effects.
5. Base-rate neglect — not applicable to the reported OR models; unweighted case counts are reported.
6. Regression to the mean — not applicable to this cross-sectional exposure/outcome comparison.
7. Survivorship bias — RED_FLAG/major limitation: CRC is prevalent at examination, so survivors are the observed cases.
8. Look-elsewhere effect — caution: sex, timing and co-exposure models are a targeted audit set and are all reported.
9. Garden of forking paths — caution: the audit was frozen before execution; prior MCOP analyses and new outputs remain distinguishable.
10. Correlation ≠ causation — caution: use association language only.
11. Reverse causality — caution: timing exclusions are an audit, not proof of temporal causality.

## Decision

The Phase 2D gate is **provisionally passed** under the pre-specified direction-only criteria. This does not upgrade NHANES to prospective evidence; it determines whether DINP-axis molecular validation is justified as the next phase.

## Output files

- mcop_crc_phase2_female_male.csv
- mcop_crc_phase2_sex_interaction.csv
- mcop_crc_phase2_diagnosis_timing.csv
- mcop_crc_phase2_reverse_causation_sensitivity.csv
- mcop_crc_phase2_coexposure_models.csv
- mcop_crc_phase2d_figure1_workflow.png through mcop_crc_phase2d_figure6_coexposure.png