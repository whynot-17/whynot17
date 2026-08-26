# Step 9B — CRC negative-branch failure attribution

## Scope

This analysis diagnoses the frozen 29-test CRC screen after the assay-specific
rebuild. It does not re-fit any model, alter the 29-test BH-FDR family, or
retroactively select candidates. It uses the recorded primary estimate and
standard error, the Step 9A merge/readiness fields, and the Step 6 robustness
artifacts. Step 6 was frozen to the two CRC FDR-supported signals, so stability
is **not** imputed for the other 27 tests.

## Frozen attribution rules

The labels are intentionally multi-valued.

- **Signal-limited:** observed OR is within 0.90–1.10, a descriptive near-null rule.
- **Power-limited:** the observed OR is outside 0.90–1.10 and the approximate two-sided alpha=0.05, 80%-power MDE derived from the current SE exceeds OR=1.20 in either direction. Near-null tests still retain their continuous MDE, but are not called power-limited merely because their MDE is wide.
- **Stability-limited:** among tests actually audited by the frozen Step 6 scope, the recorded L/C score is 0 or the H heterogeneity tag is 0. Missing Step 6 coverage is reported as not assessed.
- **Technical-warning:** the recorded primary status is not `ok`.
- **Design-limited-interpretation:** all CRC tests carry the same limitation: prevalent cross-sectional outcome, no prediagnostic biospecimen, and no stage/site/treatment/recurrence/follow-up fields. This label limits interpretation and temporality; it is not claimed to explain a P value.
- **Event support:** E2 = >=60 analytic CRC cases, E1 = 30–59, E0 = <30, following the frozen Step 6 event-support rubric.

## Summary

| Metric | Count |
|---|---:|
| Frozen CRC tests | 29 |
| Nominal P<0.05 | 5 |
| BH-FDR<0.05 | 0 |
| Observed near-null / signal-limited | 13 |
| All tests with 20% MDE above reference | 28 |
| Power-limited among non-near-null tests | 16 |
| Primary technical warnings | 1 |
| Stability assessed under frozen Step 6 scope | 2 |
| Stability-limited among assessed tests | 1 |
| E2 (>=60 cases) | 25 |
| E1 (30–59 cases) | 4 |
| E0 (<30 cases) | 0 |

## Interpretation

The CRC screen has complete model-level estimability for all 29 tests, but the
event-support profile is much weaker than T2D. The attribution table separates
three different statements that should not be collapsed: (i) an observed effect
can be near the null, (ii) the current SE may be too large to detect a prespecified
20% effect reliably, and (iii) even a statistically precise cross-sectional
association would still have limited temporality and outcome granularity.

The MDE column is therefore a sensitivity descriptor. It does not prove that the
CRC negative screen is caused by low power, and it does not imply that adding
cases would guarantee a discovery. Conversely, a test that is not flagged
power-limited is not thereby proven biologically null.

Step 6 stability results are not generalized to all 29 tests: only the two
FDR-supported CRC signals were in that frozen audit. This prevents the negative
branch from manufacturing stability claims by absence of analysis.

## Machine-readable outputs

- `step9b_crc_failure_attribution_29_tests.csv`
- `step9b_crc_failure_attribution_summary.csv`
- `STEP9B_MANIFEST.json`

The complete per-test table retains the observed OR, SE, CI, analytic cases,
retention, approximate MDE, event-support class, available stability fields,
and all applicable attribution labels.
