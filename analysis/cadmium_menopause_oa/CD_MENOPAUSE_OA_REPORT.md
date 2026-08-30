# Urinary Cd × OA: pre/post-menopause focused audit

## Scope

This is a bounded female-only analysis of urinary Cd (`URXUCD`) and the existing OA endpoint. It uses the same eight NHANES cycles, survey-weighted framework, age/race/ethnicity/PIR/smoking/cycle covariates, and sex-specific urinary-dilution logic adapted within women. No other outcome, biomarker, FDR family, mechanism, or figure analysis was performed.

Menopause classification is cycle-mapped from RHQ: `RHQ031=1` (at least one menstrual period in the past 12 months) defines pre-menopause; `RHQ031=2` plus menopause/hysterectomy reason defines post-menopause. Other reasons, refused/don't know, and missing values are excluded. The exact mapping and counts are in `01_menopause_variable_audit.csv`.

## Sample

The complete-case comparison contains 4,644 women and 751 OA cases across 8 cycles. Pre-menopause: N=2,558, OA cases=81; post-menopause: N=2,086, OA cases=670.

## Results

| estimate | β | SE | P | N | OA cases |
|---|---:|---:|---:|---:|---:|
| female pre-menopause Cd slope | -0.0362115 | 0.160325 | 0.82168 | 2558 | 81 |
| female post-menopause Cd slope | -0.0737834 | 0.0626357 | 0.24108 | 2086 | 670 |
| formal Cd × post-menopause interaction | -0.0294087 | 0.141673 | 0.835896 | 4644 | 751 |

The formal interaction, rather than a comparison of separate significance tests, determines whether the female Cd–OA slope differs between menopause groups. A positive interaction means a higher post-menopause slope under this model; it does not establish positive OA susceptibility, causality, or protection.

Under this prespecified interaction criterion, there is no evidence that menopause status modifies the urinary Cd–OA slope (Cd × post-menopause β=-0.0294, P=0.836). The pre-menopause stratum has only 81 OA cases, so its separate slope is imprecise; the conservative conclusion is that this focused analysis does not support a menopause-stratified difference in the current NHANES sample.

Run timestamp (UTC): 2026-08-30T15:01:20.004295+00:00
