# Step 10B-P — CHMS population-replacement feasibility audit

## Objective

Audit the Canadian Health Measures Survey (CHMS) as a possible independent
population source for the frozen 29 human biomarker tests and a future T2D
analysis. This stage is feasibility-only: it must not inspect or compute any
exposure–T2D association, alter the 29-test family, or select a candidate.

## Frozen questions

1. How many of the 29 tests have an exact or near-exact analyte/matrix match in
   official CHMS public content documentation?
2. Can the T2D outcome be reproduced or closely approximated from documented
   self-report and laboratory domains, subject to cycle-specific dictionary
   confirmation?
3. Are exposure, outcome, core covariates, subsample weights, and complex
   survey design variables plausibly available in the same cycle and linked at
   the respondent level?

## Evidence hierarchy

The audit distinguishes four claims:

- `content-level`: a named analyte/domain appears in an official public content
  table;
- `variable-level`: an exact CHMS variable and coding are confirmed in a
  cycle-specific dictionary;
- `joint person-level`: the exposure, T2D outcome, covariates, and design
  fields are observed for linked respondents;
- `accessible`: the required individual-level data can be obtained through the
  relevant Statistics Canada access route.

The first level cannot be silently upgraded to the latter three.

## Sources and current boundary

Primary sources are the Statistics Canada CHMS content summary for cycles 1–6,
the CHMS documents/data-dictionary index, and the CHMS Cycle 1 user guide.
Public content pages are captured as URL/status/hash metadata only; no large or
restricted microdata are downloaded.

## Outputs

- `chms_29_test_crosswalk.csv` — one row per frozen test;
- `chms_cycle_readiness.csv` — one row per test × CHMS cycle;
- `chms_outcome_covariate_design_audit.csv` — T2D/covariate/design/access domains;
- `CHMS_SOURCE_SNAPSHOT.json` — official source retrieval metadata;
- `CHMS_QC_SUMMARY.json` and `STEP10B_P_CHMS_MANIFEST.json` — audit locks;
- `STEP10B_P_CHMS_AUDIT_REPORT.md` — interpretation and gate decision.

## Promotion gate

CHMS is promoted to primary population replacement only after exact analyte
variables, T2D definition fields, same-person covariate/design availability,
environmental subsample weights, cycle-combination instructions, and controlled
individual-level access are verified. Until then, its status is
`high-priority conditional — access/harmonization pending`.
