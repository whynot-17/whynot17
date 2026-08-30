# Matched-complete-case blood versus urinary Cd × sex → OA

## Scope

This audit restricts both biomarker models to the exact same SEQN set with non-missing blood Cd (`LBXBCD`), urinary Cd (`URXUCD`), urinary creatinine, OA outcome, and all frozen covariates across the eight 2003–2018 NHANES cycles. It is a paired sample-composition check, not a new outcome screen or mechanism analysis.

The urinary model retains `log2(URXUCR)` plus its female interaction. The blood model has no urinary-dilution term. Biomarker-specific NHANES weights are retained (urine `WTSA2YR`; blood `WTMEC2YR`), so the SEQN set is matched but the survey weighting remains matrix-specific.

## Matched sample

Matched complete-case N: **10,218**; OA cases: **1,268**; male OA cases: **481**; female OA cases: **787**. All eight cycles contribute; cycle details are in `02_matched_sample_cycle_audit.csv`.

## Results

| biomarker | β(Cd×female) | SE | P | male β | female β | matched N | OA cases |
|---|---:|---:|---:|---:|---:|---:|---:|
| urinary Cd | 0.189207 | 0.0745098 | 0.0123421 | -0.2589 | -0.096151 | 10218 | 1268 |
| blood Cd | 0.00164181 | 0.0730547 | 0.982106 | -0.138278 | -0.245342 | 10218 | 1268 |

## Interpretation

After matching the complete-case SEQN set, urinary Cd retains the larger positive interaction (β=0.1892, P=0.01234), while blood Cd remains smaller and imprecise (β=0.0016, P=0.9821). Therefore the earlier biomarker discrepancy is not explained solely by blood-versus-urine sample composition in this matched set.

No formal urinary-versus-blood coefficient-difference test was calculated: the two estimates are correlated because they use the same people, and they retain matrix-specific NHANES weights and dilution terms. The audit therefore compares the two interaction estimates and their uncertainty descriptively.

This remains evidence for a biomarker-dependent difference in sex heterogeneity, not proof of different Cd toxicity, causality, or protection. Because the two matrices retain their appropriate NHANES weights and the urine model has UCr terms, this is not a claim that the two coefficients are directly interchangeable on an identical estimand scale.

Run timestamp (UTC): 2026-08-30T14:49:58.494773+00:00
