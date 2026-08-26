# Step 5 assay-specific rebuild report

Generated (UTC): 2026-08-26T10:40:37.995735+00:00

## Why this rebuild exists

The provisional Step 5 implementation used the MBzP/phthalate laboratory sample as the participant frame for all 29 tests. This rebuild keeps the frozen pre-disease test family but reads each test's own laboratory file, analyte variable, survey weight and cycle coverage. Outcome and covariates are built independently from cycle-matched NHANES MCQ/DEMO/BMX/SMQ/ALQ/DIQ/PAQ/ALB_CR files.

The provisional Step 5 outputs are retained unchanged for audit and are not overwritten.

## Screen summary

- Frozen tests entered: **29**.
- Models with finite P values: **29/29**.
- Nominal P<0.05: **5**.
- BH-FDR<0.05 using the frozen denominator 29: **0**.

| Biomarker | N | CRC cases | OR | 95% CI | P | BH-FDR (29) | Status |
|---|---:|---:|---:|---|---:|---:|---|
| URXCOP | 9936 | 70 | 1.24551 | 1.07753–1.43968 | 0.00331136 | 0.0960294 | ok |
| LBXPFHS | 10773 | 78 | 0.764697 | 0.612228–0.955137 | 0.0184578 | 0.25267 | ok |
| URXMOH | 12603 | 97 | 1.20368 | 1.0226–1.41683 | 0.0261383 | 0.25267 | ok |
| URXUSB | 13780 | 101 | 0.806672 | 0.654872–0.993658 | 0.0434879 | 0.268528 | ok |
| URXMHH | 12603 | 97 | 1.19155 | 1.00292–1.41565 | 0.0462979 | 0.268528 | ok |
| URXECP | 11224 | 85 | 1.21093 | 0.987355–1.48513 | 0.0658442 | 0.318247 | ok |
| URXUBA | 13748 | 100 | 0.847143 | 0.696643–1.03016 | 0.0958772 | 0.397206 | ok |
| LBXPFDE | 10773 | 78 | 0.770918 | 0.523677–1.13489 | 0.185421 | 0.645498 | ok |
| URXUTL | 13830 | 100 | 1.34118 | 0.848647–2.11956 | 0.20702 | 0.645498 | ok |
| URXMEP | 13775 | 104 | 0.932921 | 0.834064–1.0435 | 0.222586 | 0.645498 | ok |
| URXMHP | 13779 | 104 | 1.09293 | 0.926513–1.28924 | 0.289567 | 0.683533 | ok |
| URXP10 | 12538 | 100 | 1.14323 | 0.889594–1.46919 | 0.293221 | 0.683533 | ok |
| URXUTU | 13789 | 101 | 1.07779 | 0.923956–1.25724 | 0.338067 | 0.683533 | ok |
| URXUSR | 4398 | 31 | 0.826052 | 0.545513–1.25086 | 0.358913 | 0.683533 | ok |
| URXUUR | 11448 | 85 | 0.893766 | 0.696185–1.14742 | 0.375298 | 0.683533 | ok |
| LBXPFNA | 10773 | 78 | 0.88625 | 0.672566–1.16782 | 0.387998 | 0.683533 | ok |
| URXMBP | 13779 | 104 | 0.895411 | 0.690063–1.16186 | 0.403443 | 0.683533 | ok |
| URXUMO | 13775 | 100 | 0.920431 | 0.750242–1.12923 | 0.424262 | 0.683533 | ok |
| URXUSN | 5708 | 40 | 1.10437 | 0.806532–1.51221 | 0.530073 | 0.809059 | ok |
| URXMZP | 13779 | 104 | 0.954508 | 0.800546–1.13808 | 0.601792 | 0.829201 | ok |
| URXP25 | 4275 | 30 | 1.14991 | 0.626286–2.11134 | 0.645585 | 0.829201 | ok |
| URXBPH | 7025 | 51 | 1.03349 | 0.887349–1.20371 | 0.668284 | 0.829201 | converged_with_warning |
| URXMIB | 12603 | 97 | 0.965853 | 0.821817–1.13513 | 0.671228 | 0.829201 | ok |
| URXUPB | 13865 | 101 | 0.948276 | 0.73173–1.22891 | 0.686235 | 0.829201 | ok |
| URXP04 | 12563 | 100 | 1.03619 | 0.810043–1.32547 | 0.775715 | 0.899829 | ok |
| URXUCS | 13865 | 101 | 0.956441 | 0.630492–1.4509 | 0.833061 | 0.929183 | ok |
| URXUCD | 11482 | 87 | 1.01315 | 0.814404–1.26039 | 0.905961 | 0.973069 | ok |
| URXUCO | 13865 | 101 | 0.995562 | 0.78324–1.26544 | 0.970825 | 0.992261 | ok |
| URXP02 | 11225 | 92 | 0.998974 | 0.810516–1.23125 | 0.992261 | 0.992261 | ok |

## URXP25 resolution

URXP25 now uses PAH_H/PAH_I/PAH_J directly. Its previous N=0 was caused by zero overlap with the provisional phthalate-shaped frame, not by absent exposure values. The corrected isolated model is retained in the same rebuild output and is not hand-selected or excluded.

## Multiple-testing rule

BH-FDR is recomputed once across all **29** frozen tests, including any test that remains not estimable after the assay-specific rebuild. No result is removed before FDR, and no test is added after seeing CRC outcomes.

## QC and interpretation firewall

For each test, source rows, cycle coverage, exposure-to-outcome merge counts, analytic N, CRC cases, survey PSU/strata and fit status are written to the accompanying QC tables. This screen remains a cross-sectional prevalent-CRC analysis and does not establish causality or temporality.
