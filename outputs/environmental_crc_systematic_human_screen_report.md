# Systematic human screen after 267-chemical actionability filtering

## Scope

This run uses only candidates that passed the prespecified, CRC-outcome-blinded
permissive actionability rule. The actionability matrix was frozen before this
human association was fitted. No human OR, CI, P value, LOCO result, or
cycle-specific effect entered candidate eligibility.

## Current result

MCOP was the only currently eligible direct axis and yielded OR=1.246 (95% CI 1.078–1.440), P=0.003311; across-axis BH-FDR is 0.003311 because one axis was tested.

The current local biomonitoring audit supports **one** direct human-testable
axis: MCOP (ChemicalID C573544). Therefore the BH-FDR reported here is a
one-test correction, not evidence that all 267 chemicals have already received
an equivalent NHANES analysis.

## Model

- Seven NHANES cycles: 2005–06 through 2017–18.
- Primary comparison: CRC versus cancer-free controls.
- Exposure: log2(MCOP), per doubling.
- Covariates: age, sex, race, BMI, smoking, PIR, and log2 urinary creatinine.
- Survey design: cycle-specific phthalate subsample weight divided by 7,
  cycle-specific strata and PSU.
- The frame is rebuilt from the harmonized source and fit through the validated
  Python Taylor-style complex-survey implementation.

## Identity firewall

MCOP is analyzed as the direct candidate/urinary analyte. It is not relabeled as
MiNP, and MiNP's molecular nomination is not used to retroactively select MCOP.
The paper's direct-discovery statement is therefore:

`267 core chemicals → outcome-blinded actionability → MCOP retained → MCOP human screen`

## Limitation

The full 267-chemical biomonitoring/actionability queue is not yet complete:
265 chemicals remain manual-review unknowns. This output supports the direct
MCOP provenance audit and its downstream human test, but a claim that MCOP was
the unique winner of a fully completed 267-axis epidemiologic screen must wait
until the remaining candidates are annotated and run under the same model.

Generated: 2026-08-23T15:20:14.594724+00:00
