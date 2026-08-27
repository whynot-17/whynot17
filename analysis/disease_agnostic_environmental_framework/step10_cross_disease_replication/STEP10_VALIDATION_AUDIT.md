# Step 10 validation audit

## Status

The valid Step 10 run is complete. The first WSL execution that produced an
all-non-estimable output was invalidated because Windows registry paths had not
yet been normalized for WSL. It was not used for interpretation. The current
outputs were regenerated after that technical correction, followed by this
independent output audit.

## Frozen design checks

- Eligible outcome pool: **10** reproducible binary NHANES outcomes.
- Random panel: **5** outcomes selected with `random.Random(20260827).sample`
  after sorting the eligible disease IDs.
- Randomization occurred before exposure values and association results were
  loaded; the randomization lock records both as `false`.
- Selected panel: angina, chronic bronchitis, congestive heart failure, heart
  attack, and liver condition.
- Frozen exposure family: **29 unique tests**.
- BH-FDR denominator: **29 within each disease**, including all 29 tests.
- Outcome exclusions: T2D and CRC were excluded as worked examples and were
  not eligible randomization targets.
- No downstream disease-specific GeneCards, CTD, pathway, network,
  transcriptomic, or rescue analysis was performed.

## Independent output checks

| Check | Result |
|---|---:|
| Disease result files | 5/5 |
| Tests per disease | 29/29 |
| Unique test IDs per disease | 29/29 |
| Finite P values | 145/145 |
| Technical status `ok` | 145/145 |
| Technical warnings | 0 |
| Replacement-eligible diseases | 0/5 |
| FDR denominator values | 29 for every test |
| Randomization reproduced from seed | True |
| Generic case-field schema (`case_N`, `control_N`) | True |

## Branch results

| Disease | Pooled cases | Nominal P<0.05 | BH-FDR<0.05 | Branch |
|---|---:|---:|---:|---|
| Angina | 1,624 | 6 | 0 | Negative |
| Chronic bronchitis | 3,213 | 5 | 3 | Positive |
| Congestive heart failure | 1,906 | 10 | 6 | Positive |
| Heart attack | 2,465 | 14 | 10 | Positive |
| Liver condition | 2,065 | 4 | 2 | Positive |

Branch assignment was applied mechanically after the disease-specific
29-test screen: Positive means at least one BH-FDR < 0.05; Negative means zero.
No outcome or disease was replaced because of a visually unattractive,
non-significant, or otherwise unfavorable association result.

## Interpretation boundary

This is a randomized cross-disease replication of the framework's two-output
screening behavior, not evidence that any selected environmental biomarker
causes one of these diseases. The positive/negative labels describe the frozen
screen's output branch only. The analysis does not support disease-specific
mechanistic claims, causal claims, or generalization beyond the selected
NHANES outcomes.
