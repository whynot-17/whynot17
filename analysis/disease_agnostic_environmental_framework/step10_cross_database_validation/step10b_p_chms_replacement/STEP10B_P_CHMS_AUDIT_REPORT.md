# Step 10B-P — CHMS population-replacement feasibility audit

Generated: `2026-08-28T10:03:21.801639+00:00`
Status: **`high_priority_conditional_not_access_confirmed`**

## Scope

This audit asks whether CHMS is a viable independent population source for the
frozen 29 human biomarker tests and a T2D analysis. It does **not** inspect any
exposure–T2D association result, change the 29-test family, or promote a
candidate.

## Executive result

The official cycles 1–6 content summary supports **27/29**
tests at the content level: **23** named analytes or naming
variants and **4** grouped PAH-family matches.
**2** tests are not confirmed in the public
cycles 1–6 summary (`URXUBA` urine barium and `URXUSN` urine tin). These are
not treated as negative biological findings.

The audit confirms **0/29 exact CHMS variable-level mappings** and **0/29
person-level joint exposure–T2D confirmations**, because those require the
cycle-specific data dictionaries and controlled individual-level files. The
correct current disposition is therefore **high-priority conditional**, not
primary replacement ready.

## Exposure crosswalk

The full 29-row crosswalk is in `chms_29_test_crosswalk.csv`; the per-cycle
174-row audit is in `chms_cycle_readiness.csv`. The public content table shows:

- PFDA, PFHxS, PFNA in blood, with cycle- and age-specific restrictions;
- BPA in urine across cycles 1–6;
- nine named phthalate metabolites, including MCIOP/MCiOP as the public naming
  variant relevant to the frozen MCOP mapping;
- four PAH grouped-family matches based on named hydroxylated PAHs;
- ten named metals in urine, mostly concentrated in cycles 1–2, while barium
  is listed in hair rather than urine and tin was not confirmed.

Calendar-window overlap with NHANES is reported descriptively only. It is not
treated as same-cycle equivalence.

## T2D, covariates, and survey design

Public CHMS documentation lists glucose and HbA1c in the diabetes laboratory
domain across cycles 1–6, and also provides a self-reported chronic-condition
questionnaire domain. A reproducible T2D definition is **not yet frozen**:
fasting/random glucose, HbA1c, self-report, type 1 exclusion, age eligibility,
and missingness rules must be resolved from the cycle dictionaries before any
model is specified.

Age, sex/gender, anthropometry/BMI, smoking, alcohol, physical activity,
education/income/SES, and urine creatinine are plausible/common covariate
domains. The CHMS user guide and documents page confirm that full-sample and
subsample weights and instructions for combining cycles exist. The exact
environmental blood/urine weight variables, design variables or replicate
weights, and cycle-specific joins remain pending.

Individual-level microdata are not supplied by the public content dashboard;
access is controlled. Thus no claim is made that an exposure, T2D outcome,
covariate set, and design variables are already jointly available for the same
respondents.

## Gate decision

**Promotion to primary epidemiologic replacement: not yet passed.** The next
action is a controlled data-dictionary/access request for the exact analyte
variables, environmental subsample weights, design/replicate variables, T2D
definition fields, covariates, and respondent-level linkage. If those gates
pass, CHMS becomes the primary external population replacement; otherwise it
remains a high-priority conditional source.

## Official source boundary

Source URLs and retrieval metadata are in `CHMS_SOURCE_SNAPSHOT.json`. The
main content evidence is the Statistics Canada *Content summary for cycles 1
to 6*; the official documents page explicitly lists separate environmental
blood/urine and environment-urine-main-subsample dictionaries and warns that
not every subsample is present in every cycle. No restricted file was obtained.
