# MBzP–CRC Phase 2B NHANES 人群验证 + 18 基因机制桥接

## 首页六问

1. **Model 2 MBzP OR per doubling:** 0.953725; 95% CI 0.801141–1.13537; P=0.592117.
2. **Q4 vs Q1:** OR=0.963856; 95% CI 0.466448–1.99169; P-trend=0.699774.
3. **Age ≥40:** OR=0.982558; 95% CI 0.827937–1.16605; P=0.839383.
4. **LOCO OR range:** 0.864728–0.989579; direction consistent: YES (all LOCO point estimates are below 1 when YES).
5. **MBzP among phthalates:** OR rank=7.0; P rank=7.0; FDR rank=7.0.
6. **18 overlap genes exported:** **YES** (18 genes).

## Scope and definitions

- Primary population: CRC (colon or rectal cancer type code) versus participants reporting no cancer history (`MCQ220=2`).
- Sensitivity population: CRC versus participants with known cancer outcome who are not CRC, including other cancer histories.
- Exposure: `log2(URXMZP)`; OR is per doubling of urinary MBzP.
- Primary model: age, sex, race, BMI, smoking, PIR, and `log2(URXUCR)`.
- All estimates use CDC-compatible pooled phthalate subsample weights: `WTSPH4YR×2/10` for 1999–2002 and cycle-specific 2-year weights divided by 10 thereafter, with pooled strata and PSU identifiers and Taylor-style stratified PSU sandwich variance. Counts are unweighted.

## Model 0–3

See `mbzp_crc_phase2_main_models.csv`. Model 3 is sensitivity-only when its complete-case CRC count is below 80.

## Sensitivity and specificity

See `mbzp_crc_phase2_sensitivity.csv`, `mbzp_crc_phase2_leave_one_cycle_out.csv`, `mbzp_crc_phase2_phthalate_comparison.csv`, `mbzp_crc_phase2_phthalate_correlations.csv`, and `mbzp_crc_phase2_phthalate_burden.csv`.

## Timing and spline

Diagnosis timing is descriptive and based on available cancer diagnosis-age fields. Restricted cubic spline uses 5th/35th/65th/95th percentile knots only when the Model 2 primary complete-case CRC count is at least 80; otherwise the output records the prespecified non-run status.

## Fixed molecular bridge

The CTD human-interacting MBzP genes intersected with the Phase 1 GeneCards Disorders CRC set to produce exactly 18 genes. The full table retains gene-specific interaction rows and unique PMID counts; no GO, KEGG, PPI, docking, WGCNA, machine learning, or hub-gene fishing was performed.

Run timestamp (UTC): `2026-08-22T08:01:29.263570+00:00`
