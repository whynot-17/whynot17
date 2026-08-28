# Step 10C — CCHS public population-replacement feasibility audit

## Decision

`CCHS 2022 PUMF` passes the **publicly downloadable, individual-level microdata** gate. The archive was downloaded directly from Statistics Canada and verified locally by SHA-256. No registration, application, or controlled-access step was used for this package.

The feasibility-only primary freeze is:

> **Current smoking status (`SMKDVSTY`) → self-reported high blood pressure (`CCC_80`) in adults**

This is a source-native exposure demonstration, not an exact replication of the 29 NHANES biomarker family. The exposure is defined as current daily/occasional smoking (codes 01/02) versus not currently smoking (03–06). Hypertension is `CCC_80` yes/no (1/2). Adults are `DHHGAGE` groups 2–5 (18 years and older).

## Package and design audit

- Direct package: `https://www150.statcan.gc.ca/n1/pub/82m0013x/2024001/2022_TXT.zip`
- Local archive: `D:\whynot17\public_sources\cchs_2022\2022_TXT.zip`
- SHA-256: `a116a01fb35cc3204a8b14e36dda70cf6722e75380f88b6d822e492765ac8c41`
- Master PUMF rows: `67,079`; fixed-width record length: 385; malformed-width rows: 0
- Bootstrap rows: `67,079`; fixed-width record length: 8029; malformed-width rows: 0
- Master and bootstrap row counts match: **True**
- The official bootstrap layout card contains `BSW1–BSW1000`.
- Person-level survey weight: `WTS_M`
- Variance plan: use the package bootstrap replicate weights; the replicate scaling/multiplier must be confirmed from the CCHS 2022 user guide before modeling.

## Variable feasibility

The public PUMF label/input cards contain exact fields for age group, sex, adult BMI classification, education, income, smoking status, high blood pressure, diabetes, alcohol, physical activity, `WTS_M`, and the bootstrap file. The main public label card did **not** contain an exact second-hand-smoke variable; the official questionnaire concept must not be treated as an available public-PUMF exposure without a mapped field.

The candidate-combination table reports only variable validity and descriptive event support. It contains no association estimate and was not used to optimize a result.

## Firewall and scope

`association_models_run = false`. No exposure–outcome regression, P value, odds ratio, or FDR was used to select or freeze the combination. The next step, if authorized, is a separate preregistered CCHS association runner using the frozen combination and bootstrap survey variance.

The large raw archive and raw microdata are deliberately excluded from version control; this directory contains only the small audit artifacts and provenance needed to reproduce the package check.
