# Step 10 — Randomized cross-disease replication

This directory implements the prespecified cross-disease replication of the
frozen, outcome-independent 29-test NHANES environmental hypothesis family.

## Frozen rules

- T2D and CRC are excluded as worked examples; they are not eligible randomization targets.
- Disease eligibility uses only reproducibility, at least three NHANES cycles,
  pooled adult case support of at least 100, and compatibility with the fixed
  core covariates (age, sex, race/ethnicity, BMI, smoking, and PIR).
- The disease panel is sampled with `RANDOM_SEED = 20260827` from the sorted
  eligible pool, before any exposure association results are read.
- Each selected disease receives the same 29 frozen tests and the same survey
  logistic model family. Each disease has its own BH-FDR family with denominator
  29.
- Positive branch: at least one test has BH-FDR < 0.05. Negative branch: zero
  tests have BH-FDR < 0.05.
- Technical replacement is allowed only for a non-reconstructible outcome,
  fewer than 50% technically estimable tests, no compatible survey design, or
  an outcome-definition error. A visually unattractive result is never a
  replacement reason.
- The replication panel stops after the 29-test screen, BH-FDR, and branch
  assignment. No disease-specific GeneCards, CTD, pathway, network,
  transcriptomic, or rescue analysis is part of Step 10.

## Firewall statement

The registry and randomization scripts do not import, open, or inspect any
29-test exposure association result. The randomization lock records the frozen
pool hash, seed, timestamp, and code hash before the disease screening runner
is allowed to execute.
