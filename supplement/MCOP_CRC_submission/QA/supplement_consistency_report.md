# MCOP–CRC Supplement consistency and submission audit

Numeric/package checks: **23 PASS, 0 FAIL**.

## Frozen numeric consistency

| Check | Status | Observed | Expected |
|---|---:|---:|---:|
| Starting chemical universe | PASS | 267 | 267 |
| Human-testable mappings | PASS | 87 | 87 |
| Strict-eligible mappings | PASS | 27 | 27 |
| Unique biomarker tests | PASS | 15 | 15 |
| Unified screen tests | PASS | 15 | 15 |
| BH-FDR-supported tests | PASS | 2 | 2 |
| Robust Tier A axes | PASS | 1 | 1 |
| Primary complete-case N | PASS | 9936 | 9936 |
| Primary CRC events | PASS | 70 | 70 |
| Primary OR | PASS | 1.2455068 | 1.2455068 ±1e-6 |
| Primary design-df P | PASS | 0.0033113 | 0.0033113 ±1e-6 |
| R/Python logOR agreement | PASS | 1.071e-13 | <1e-10 |
| Frozen PPAR/NR paired donors | PASS | 36 | 36 |
| Frozen PPAR/NR median delta | PASS | -0.418601 | -0.418601 ±1e-6 |
| Frozen PPAR/NR FDR | PASS | 9.299838e-07 | 9.299838e-7 ±1e-12 |
| Supplement avoids causal language | PASS | no prohibited direct-causal phrase | no prohibited direct-causal phrase |
| Supplement locks chemical identity | PASS | identity lock present | identity lock present |
| Supplementary table sheets | PASS | 9 | 9 |
| Figure source-data sheets | PASS | 13 | >=12 |
| Figure S1 export triad | PASS | PDF/SVG/PNG | PDF/SVG/PNG |
| Figure S2 export triad | PASS | PDF/SVG/PNG | PDF/SVG/PNG |
| Figure S3 export triad | PASS | PDF/SVG/PNG | PDF/SVG/PNG |
| Figure S4 export triad | PASS | PDF/SVG/PNG | PDF/SVG/PNG |

## Main-manuscript integration items

- **Supplementary Table X placeholders in current manuscript:** 3. Replace with S8 (gate fields), S1/S8 (attrition), and S4 (robustness rubric).
- **Main Figure 2 asset:** requires final integration check. Repository still contains pre-lock candidate-triage Figure2 files; use the outcome-blinded 267→87→15 architecture before submission.

## Interpretation lock

- NHANES is cross-sectional and estimates prevalent CRC association, not incident risk.
- MiNP, DINP parent and MCOP are chemically distinct; MCOP is the urinary biomarker for a DINP-related axis.
- PPAR/NR is an independently observed CRC epithelial disease-state program; the DINP/MCOP-to-state bridge remains untested.
- LOCO analyses are overlapping pooled re-estimates, not seven independent replications.
- The incomplete 35/36 live Census staged cache remains an INFO item; final transcriptomic inference uses the validated official H5AD.