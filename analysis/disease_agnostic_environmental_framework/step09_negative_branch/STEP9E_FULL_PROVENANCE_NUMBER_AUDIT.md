# Step 9E — Full provenance / number audit

Generated (UTC): 2026-08-26T19:25:34.578109+00:00

## Audit decision

**Canonical count audit: PASS (40 assertions; 0 mismatches).**

Every headline number requested for the current disease-agnostic framework was
recomputed from the source tables rather than copied from narrative reports.
The historical conflict ledger is intentionally separate: a legacy number may
be preserved for provenance while being explicitly marked superseded or
reconciled under a different unit/frame.

## Canonical chain

| Stage | Recomputed canonical count | Meaning |
|---|---:|---|
| Step 1 | 2,042 | disease-agnostic environmental chemical entities |
| Step 3 | 411 | actionable chemical–biomarker mapping rows |
| Step 4 | 409 | unique chemical IDs represented across those mappings |
| Step 4 | 29 | frozen NHANES biomarker tests |
| Step 5 T2D | 29 / 29 | estimable tests |
| Step 5 T2D | 14 | BH-FDR < 0.05 tests |
| Step 6 T2D | 13 / 14 | robust FDR candidates |
| Step 6 T2D | 11 | exposure clusters |
| Step 7 T2D | 5 | GeneCards-enriched clusters |
| Step 7 T2D | 4 | Tier A clusters |
| Step 8 T2D | 1,647 | globally significant pathway terms |
| Step 8 T2D | 321 | redundancy-reduced pathway modules |
| Step 8C T2D | 97 | STRING/Louvain network modules |
| Step 8E T2D | 4 | integrated Tier A axes |

## CRC negative branch

| Quantity | Recomputed value | Correct interpretation |
|---|---:|---|
| Assay-specific pooled CRC cases | 420 | primary Step 9A readiness frame |
| Median analytic CRC cases/test | 97 | median across 29 test-specific complete-case models |
| CRC tests, estimable | 29 / 29 | 28 clean fits + 1 retained convergence warning |
| CRC nominal P < 0.05 | 5 | nominal only |
| CRC BH-FDR < 0.05 | 0 | canonical rebuilt screen |
| CRC signal-limited | 13 / 29 | near-null descriptive class |
| CRC power-limited | 16 / 29 | non-near-null but imprecise for OR=1.20 reference |

The 13 and 16 failure-attribution classes are mutually exclusive under the
locked rules; they are not a robustness count and must not be confused with
the T2D 13/14 robust subset.

## Historical conflicts resolved

- **267 vs 2,042:** 267 is the earlier CRC-specific matrix; 2,042 is the
  current disease-agnostic upstream universe.
- **411 vs 409:** mapping-row memberships versus unique chemical IDs; the
  difference is the repeated membership of D004051 across three tests.
- **420 vs 123:** assay-specific rebuilt CRC outcome QC versus the legacy
  diagnosis-age/case-control ledger; different frames, not additive counts.
- **420 vs 97:** pooled outcome-QC events versus the median test-specific
  complete-case event count.
- **T2D 13 vs CRC 13/16:** different branches and different definitions.
- **Old CRC 2 FDR positives vs current 0:** the old phthalate-shaped frame is
  superseded by the assay-specific rebuild; current CRC screening is 5 nominal
  and 0 BH-FDR-positive.
- **111 vs 20,554 GeneCards genes:** the 111-row exact scoped query is a
  deprecated preflight artifact; the only formal Step 7 input is 20,554 genes.

## Reproducibility outputs

- `step9e_full_provenance_number_audit.csv`: assertion-by-assertion source,
  rule, recomputed value and status.
- `step9e_historical_number_conflict_ledger.csv`: explicit reconciliation of
  legacy/provisional values.
- `step9e_source_file_snapshot.csv`: source hashes, row counts and purpose.
- `run_step09e_full_provenance_number_audit.py`: rerunnable audit script.

Raw expression matrices, Census caches and other large artifacts are not
included in this audit package.
