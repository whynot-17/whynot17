# CRC outcome definition for Step 5

This document defines the CRC disease plug-in used only after the Step 4 environmental test set was frozen.

## Primary population

- Source: NHANES Medical Conditions Questionnaire (`MCQ`).
- Adult restriction: age >=20 years (`RIDAGEYR`).
- CRC case: `MCQ220=1` and at least one cancer-type code `16` (colon) or `31` (rectum) in the cycle-appropriate cancer-type fields.
- Cancer-free control: `MCQ220=2`.
- Participants reporting a known non-CRC cancer are excluded from the primary CRC-versus-cancer-free comparison.
- Unknown/missing cancer history is excluded from the primary analysis.
- Diagnosis age is retained for later reverse-causation sensitivity analyses; it does not define the primary case/control status.

## Variables retained in the case/control ledger

`SEQN`, cycle, age, cancer outcome availability/known status, CRC case, colon case, rectal case, cancer-free status, CRC diagnosis age, and years since CRC diagnosis when available.

## Primary model

For every frozen unique NHANES biomarker test: survey-weighted logistic regression of prevalent CRC on log2 biomarker concentration, age, sex, race/ethnicity, BMI, smoking, and poverty-income ratio. Urinary biomarkers additionally include log2 urinary creatinine. Each analyte uses its own laboratory/subsample weight, cycle-specific strata and PSU, and pooled weights divided by the number of included cycles.

The BH-FDR denominator is fixed at **29 frozen tests**, including tests whose model is not estimable.

Rebuilt adult participants with known cancer history in the harmonized frame: **17,382**.
NHANES source records hashed in the Step 5 run manifest: **90**.
Model implementation: `C:\Users\21634\Documents\Codex\2026-08-22\non\work\whynot17\work\scripts\mbzp_crc_phase2b.py`.

This is a prevalent CRC association screen and does not establish temporality or causality.
