# Step 9A — Question-specific data-readiness audit

## Scope and status

This is a descriptive audit of the same frozen 29-test environmental biomarker
panel in the T2D and assay-specific rebuilt CRC outcome screens. It is the first
layer of the CRC negative branch. It does **not** re-fit models, change the
29-test multiplicity family, perform failure attribution, or simulate larger
studies. The purpose is to document whether the two disease questions are being
asked under comparable data conditions.

The analysis deliberately preserves a multidimensional readiness profile rather
than collapsing it into a single score.

## Executive comparison

| Metric | T2D | CRC |
|---|---:|---:|
| NHANES cycles in outcome QC | 10 | 10 |
| Outcome frame rows | 55,081 | 96,811 |
| Eligible case/control rows | 49,181 | 50,275 |
| Outcome cases in pooled outcome QC | 7,772 | 420 |
| Outcome case fraction among eligible rows | 15.8% | 0.8% |
| Frozen tests | 29 | 29 |
| Estimable tests | 29 | 29 |
| Fit warnings retained in status | 0 | 1 |
| Analytic N across tests (min / median / max) | 4,406 / 12,994 / 14,343 | 4,275 / 12,603 / 13,865 |
| Analytic cases across tests (min / median / max) | 712 / 1,940 / 2,120 | 30 / 97 / 104 |
| Median analytic case fraction | 14.8% | 0.7% |
| Frozen cycle coverage (min / median / max) | 3 / 9.0 / 10 | 3 / 9.0 / 10 |
| Cycles with complete cases (min / median / max) | 3 / 9.0 / 10 | 3 / 9.0 / 10 |
| Median analytic retention from merged exposure frame | 81.7% | 49.5% |
| Approx. median loss after merge to analytic complete cases | 18.1% | 50.3% |
| Nominal P<0.05 (context only) | 14/29 | 5/29 |
| BH-FDR<0.05 (context only) | 14/29 | 0/29 |

The main readiness contrast is not simply “significant versus null.” CRC has a
much lower case density and fewer outcome events per biomarker-specific analytic
sample, while its primary outcome is prevalent and cross-sectional. T2D has a
larger outcome case pool and all 29 tests are estimable, but it is also a
same-visit cross-sectional classification rather than prospective incidence.

The CRC rebuilt output contains one `converged_with_warning` fit. It remains
estimable because it has a finite model P value and a non-zero analytic sample;
the warning is retained as a technical-quality flag, not silently converted to
`ok`.

## CRC readiness profile

1. **Case density and precision:** the assay-specific rebuilt outcome QC contains 420 CRC cases across the pooled CRC-versus-cancer-free frame, with a median of 97 complete-case CRC events per biomarker model. The median analytic CRC case fraction is 0.7%.
2. **Assay-specific analytic overlap:** the rebuilt screen uses the correct laboratory file and subsample weight for each assay family. Even after that correction, biomarker-specific merge and complete-case retention vary across tests; this is a measurable source of effective sample-size loss rather than a generic “NHANES N.”
3. **Cycle coverage:** the frozen panel has CRC assay coverage ranging from 3 to 10 cycles. The per-test table records the exact cycles and the number of cycles retaining complete cases.
4. **Outcome granularity:** CRC is defined as prevalent CRC versus cancer-free control. Diagnosis age is retained for 119/123 CRC cases in the separate outcome ledger, but diagnosis date, stage, detailed site, treatment, recurrence/progression, and prospective follow-up are unavailable in this screen. The diagnosis-age ledger is a separate provenance artifact and is not used to replace the assay-specific rebuilt outcome totals.
5. **Temporality:** the urine/serum measurement is contemporaneous with the survey and is not a prediagnostic biospecimen. Therefore the screen can identify an association under the frozen outcome definition, but cannot establish that exposure preceded CRC.

## T2D comparator profile

T2D has 7,772 pooled eligible outcome cases and a median of 1,940 cases per biomarker model. All 29 tests are estimable under their assay-specific frames. Its case definition uses diagnosed diabetes plus an objective HbA1c rule for probable undiagnosed disease, with indeterminate categories excluded from the primary case/control comparison. T2D remains cross-sectional and does not provide prospective exposure-to-onset timing.

## What this audit does and does not establish

This audit establishes that the CRC and T2D questions have materially different
readiness profiles: CRC is particularly constrained by case density, effective
case/sample overlap, and outcome temporality/granularity. It does **not** show
that any CRC association is absent, that a larger sample would guarantee an
association, or that the T2D findings are causal. Those questions belong to the
next negative-branch layer (failure attribution) and to external prospective
validation, respectively.

## Reproducibility and QC

- CRC results are from the assay-specific rebuilt 29-test screen, not the superseded phthalate-shaped frame.
- T2D results are from the frozen assay-specific 29-test screen.
- No CRC or T2D model was re-fit by this script.
- The 29-test BH-FDR results are included only as context; no threshold was used to define readiness.
- Per-test audit totals are compared with the recorded model-level analytic N and case counts. Any non-zero audit discrepancy is retained in the output for inspection rather than silently corrected.
- The outcome-granularity table explicitly records the timing and missing follow-up fields.
- The rebuilt assay-specific outcome QC and the legacy CRC case/control ledger have different row/case totals; they are not pooled. The rebuilt QC supplies primary readiness counts, while the legacy ledger is used only to document diagnosis-age availability.

Output files:

- `step9a_readiness_overview.csv`
- `step9a_per_test_readiness.csv`
- `step9a_outcome_granularity.csv`
- `STEP9A_MANIFEST.json`
