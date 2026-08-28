# Step 10C — CCHS public population-replacement audit

## Objective

Identify a population-replacement demonstration that uses only an individual-level health dataset that can be downloaded directly, without registration, application, RDC, or controlled-access approval.

## Frozen scope

- Use the CCHS 2022 TXT PUMF package downloaded from the official Statistics Canada PUMF page.
- Do not force exact replication of the 29 NHANES environmental biomarker tests.
- Treat CCHS source-native questionnaire/behavioral exposures as eligible for this database-replacement stress test.
- Do not inspect exposure–outcome association results during this audit.
- Freeze one exposure–outcome combination using only variable availability, valid coverage, outcome support, covariate availability, and survey design availability.

## Primary feasibility candidate

`SMKDVSTY` current smoking status → `CCC_80` self-reported high blood pressure among adults (`DHHGAGE` groups 2–5).

Core covariates: sex, adult BMI classification, education, household income. Survey weight: `WTS_M`; variance: package bootstrap replicates `BSW1–BSW1000`, with scaling to be confirmed before modeling.

## Explicit non-goals

- No regression, P value, odds ratio, confidence interval, or FDR.
- No candidate selection by literature, prior association, or novelty.
- No inference that the questionnaire's second-hand-smoke concept is present in the public PUMF when no exact public variable is labeled/mapped.
- No commitment of the 193 MB raw archive or raw microdata to Git.

## Reproducible command

```powershell
C:\Users\21634\anaconda3\python.exe `
  analysis/disease_agnostic_environmental_framework/step10_cross_database_validation/step10c_public_population_replacement_cchs/run_cchs_public_dictionary_audit.py
```
