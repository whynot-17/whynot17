# Step 10-R — Expanded outcome-pool replication lock

## Purpose

This is an additive patch to the original Step 10 randomized disease
replication. It broadens the outcome sampling space beyond the MCQ160
physician-diagnosis items while preserving the same frozen 29-test exposure
family and the same survey-weighted model family.

## Outcome pool locked before exposure results

The candidate inventory is defined from the available NHANES questionnaire
modules and official binary response coding only:

- the original ten MCQ160 physician-diagnosis outcomes;
- MCQ010, physician-diagnosed asthma;
- BPQ020, ever told of hypertension/high blood pressure;
- BPQ080, ever told of high blood cholesterol.

All candidates use adults aged >=20, response 1 versus 2, at least three
cycles, and at least 100 pooled adult cases. T2D and CRC are not sampled as
replication outcomes. The BPQ source files were obtained from the official
NHANES public data releases and are analyzed with the same MEC-compatible
survey design variables used by the frozen exposure tests.

## Randomization lock

- Pool audit is completed before the randomization file is generated.
- Fixed seed: `20260828`.
- Requested randomized panel: five outcomes.
- Sampling rule: `random.Random(seed).sample(sorted(eligible disease_id), k)`.
- Exposure values, exposure association estimates, P values, FDR values, and
  downstream biology are not read during pool construction or randomization.

## Analysis lock

Each selected disease is screened separately against all 29 tests. The model
is the existing survey-weighted logistic implementation:

`disease ~ log2(exposure) + age + sex + race + BMI + smoking + PIR`

Urine tests additionally include `log2(urinary creatinine)`. The BH-FDR
denominator is fixed at 29 per disease. No replacement is allowed based on an
association result; the only replacement conditions are a non-reconstructible
outcome, <50% technically estimable tests, absent compatible survey design, or
an invalid outcome definition.

## Interpretation boundary

Step 10-R tests whether an outcome-firewalled environmental test family can
be applied across a broader, questionnaire-defined NHANES disease pool. It
does not establish disease causality, temporal ordering, or independence
among outcomes.
