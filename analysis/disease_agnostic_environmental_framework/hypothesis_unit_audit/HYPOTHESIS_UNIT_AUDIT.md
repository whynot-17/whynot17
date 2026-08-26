# Step 4 hypothesis-unit audit

## Scope

This audit reads only the outcome-blinded Step 1–4 mapping artifacts. It does not read CRC results, GeneCards, CTD chemical–gene interactions, or any downstream model output. It changes no FDR denominator.

## Direct findings

- Frozen Step 4 biomarker tests: **29**.
- Actionable chemical–biomarker mapping rows represented by those tests: **411**.
- Unique CTD chemical IDs represented across the tests: **409**.
- Distinct exposure-axis labels: **7**, with **31** test-to-axis memberships because some tests carry more than one axis label.
- Supplemental parent annotations are present for **6/29** tests; missing parent metadata is not inferred from metabolite names.
- Known-parent collapse produces **26** operational groups, but this is not a replacement for the frozen 29-test primary family.

## Current primary unit

The canonical Step 4 lock explicitly states that collapse was performed only at the unique NHANES test level and sets the planned downstream FDR denominator to 29. Therefore the defensible primary multiplicity family remains **29 measured biomarker tests**. The 29 tests are not 29 independent chemicals: several tests represent multi-chemical proxies, and several phthalate tests can belong to one recorded parent axis.

## Descriptive grouping

The operational secondary grouping uses a parent compound only when it is explicitly recorded in the supplemental identity table. Otherwise it retains a biomarker/proxy-level unit rather than guessing a parent.

| Grouping | Count | Interpretation |
| --- | ---: | --- |
| Frozen biomarker-test family | 29 | Current Step 4 primary unit and flat-FDR family |
| Known-parent operational groups | 26 | Secondary descriptive grouping; incomplete parent annotation |
| Exposure-axis labels | 7 | Labels, not a valid denominator because memberships overlap and vary in specificity |

## Complete 29-test map

| Biomarker test | Chemical class | Mappings | Recorded parent | Operational unit | Parent annotation |
| --- | --- | ---: | --- | --- | ---: |
| LBXPFDE | PFAS | 1 | not recorded | direct::LBXPFDE | 0.0% |
| LBXPFHS | PFAS | 1 | not recorded | direct::LBXPFHS | 0.0% |
| LBXPFNA | PFAS | 1 | not recorded | direct::LBXPFNA | 0.0% |
| URXBPH | bisphenols | 34 | not recorded | proxy::URXBPH | 0.0% |
| URXCOP | phthalates | 3 | DINP | parent::DINP | 100.0% |
| URXECP | phthalates | 2 | DEHP | parent::DEHP | 50.0% |
| URXMBP | phthalates | 1 | not recorded | biomarker::URXMBP | 0.0% |
| URXMEP | phthalates | 1 | not recorded | biomarker::URXMEP | 0.0% |
| URXMHH | phthalates | 1 | DEHP | parent::DEHP | 100.0% |
| URXMHP | phthalates | 1 | DEHP | parent::DEHP | 100.0% |
| URXMIB | phthalates | 1 | not recorded | biomarker::URXMIB | 0.0% |
| URXMOH | phthalates | 2 | DEHP | parent::DEHP | 50.0% |
| URXMZP | phthalates | 2 | BBzP | parent::BBzP | 100.0% |
| URXP02 | PAHs | 12 | not recorded | proxy::URXP02 | 0.0% |
| URXP04 | PAHs | 5 | not recorded | proxy::URXP04 | 0.0% |
| URXP10 | PAHs | 95 | not recorded | proxy::URXP10 | 0.0% |
| URXP25 | PAHs | 24 | not recorded | proxy::URXP25 | 0.0% |
| URXUBA | metals | 8 | not recorded | elemental_proxy::URXUBA | 0.0% |
| URXUCD | metals | 7 | not recorded | elemental_proxy::URXUCD | 0.0% |
| URXUCO | metals | 50 | not recorded | elemental_proxy::URXUCO | 0.0% |
| URXUCS | metals | 22 | not recorded | elemental_proxy::URXUCS | 0.0% |
| URXUMO | metals | 18 | not recorded | elemental_proxy::URXUMO | 0.0% |
| URXUPB | metals | 29 | not recorded | elemental_proxy::URXUPB | 0.0% |
| URXUSB | metals | 22 | not recorded | elemental_proxy::URXUSB | 0.0% |
| URXUSN | metals | 18 | not recorded | elemental_proxy::URXUSN | 0.0% |
| URXUSR | metals | 16 | not recorded | elemental_proxy::URXUSR | 0.0% |
| URXUTL | metals | 9 | not recorded | elemental_proxy::URXUTL | 0.0% |
| URXUTU | metals | 12 | not recorded | elemental_proxy::URXUTU | 0.0% |
| URXUUR | metals | 13 | not recorded | elemental_proxy::URXUUR | 0.0% |

## Known parent links

- `DINP`: URXCOP.
- `DEHP`: URXECP, URXMHH, URXMHP, URXMOH.
- `BBzP`: URXMZP.
- URXMBP, URXMEP, and URXMIB remain parent-unresolved in the current Step 1–4 ledger; no parent was imputed.
- Parent annotation is partial for URXECP and URXMOH (50% of their mapped rows in the supplemental identity table); this is why the 26-unit grouping is explicitly operational rather than a complete parent ontology.

## Statistical interpretation guardrail

This audit supports a transparent distinction between biomarker-level and parent-level hypotheses, but it does not rescue a CRC result by changing the denominator after outcome inspection. If a parent-level or hierarchical FDR analysis is later added, it must be reported as a secondary/post hoc reanalysis unless the hierarchy is frozen before CRC outcome access. The flat 29-test result remains the primary reference for the current manuscript state.

## Files

- `step4_test_hypothesis_mapping.csv`: one row per frozen biomarker test.
- `step4_test_chemical_membership.csv`: one row per chemical membership in those tests.
- `hypothesis_unit_summary.csv`: operational parent/proxy grouping.
- `HYPOTHESIS_UNIT_AUDIT_LOCK.json`: source hashes and audit lock.
