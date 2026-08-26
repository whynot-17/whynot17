# Step 9C — Counterfactual data-improvement sensitivity analysis

## Scope and firewall

This is a planning sensitivity analysis for the CRC negative branch. It does
not re-fit the CRC models, create new P values or BH-FDR values, change the
29-test family, or claim that an enlarged study would produce an association.
The event and retention sections only rescale the recorded Step 9B estimates
under explicit assumptions.

## 1. Event expansion

For each test, cases were scaled by ×1, ×2, ×4, and ×8. The observed log(OR)
was held constant and `SE_new = SE_current / sqrt(case multiplier)`. The output
reports expected case count, CI width, and the approximate 80%-power MDE. The
reference effect is OR=1.20 per exposure doubling; it is a planning reference,
not a clinical threshold.

| Scenario | Tests with approximate MDE <= OR 1.20 |
|---|---:|
| Current | 1 / 29 |
| Cases ×2 | 9 / 29 |
| Cases ×4 | 20 / 29 |
| Cases ×8 | 24 / 29 |
| Not meeting the MDE reference even at ×8 | 5 / 29 |

This simulation answers detectability under the stated scaling assumption. It
does not identify which tests are biologically real, and it does not simulate
selection, confounding, exposure correlation, survey-design changes, or FDR.
Tests whose observed OR is near the null remain near-null in the scenario even
if their precision improves.

## 2. Analytic-retention improvement

For targets of 70% and 85% retention from the exposure-merged frame, expected
analytic N was set to `merge N × target retention` (capped at merge N), expected
case fraction was held at the observed analytic case fraction, and SE was
scaled by the square root of the current/expected case count. This is a
non-differential missingness approximation. It is not a reanalysis of any
participant-level data.

| Scenario | Tests with approximate MDE <= OR 1.20 |
|---|---:|
| Retention 70% | 2 / 29 |
| Retention 85% | 7 / 29 |

The retention table preserves each biomarker's current merge N, current case
fraction, expected analytic N, expected cases, and counterfactual MDE. It does
not assume that better retention repairs outcome timing or phenotype definition.

## 3. Structural outcome improvements

Structural changes cannot be assigned a new P value from the present NHANES
data. Their value is inferential rather than purely statistical:

- prediagnostic biospecimens primarily improve temporal ordering and protection against reverse causation;
- diagnosis date enables lag/timing analyses;
- stage and site improve phenotype resolution;
- treatment/follow-up/recurrence support incidence, prognosis, and longitudinal outcome analyses;
- repeated exposure measurements reduce exposure misclassification and identify persistent exposure.

The full structural matrix and priority matrix are machine-readable outputs.

## Overall interpretation

The CRC negative branch is therefore actionable but not self-exonerating. Event
expansion and better retention can improve detectability under transparent
assumptions, while prediagnostic and longitudinal outcome data address a
different limitation—temporality and interpretation. No row in this report
means that a future study is expected to be positive; it means that a specified
data addition would improve the ability to answer the CRC question.

## Outputs

- `step9c_event_expansion_simulation.csv`
- `step9c_retention_improvement_simulation.csv`
- `step9c_structural_outcome_improvements.csv`
- `step9c_data_improvement_priority_matrix.csv`
- `step9c_counterfactual_summary.csv`
- `STEP9C_MANIFEST.json`
