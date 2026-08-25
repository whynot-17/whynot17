# MCOP-CRC Phase 2：DINP 轴最终人群验证

## 判定

本轮只分析 MCOP（NHANES URXCOP），不再分析 MONP 或 MiNP。分析使用既有 20 岁以上 NHANES CRC harmonized frame、cancer-free controls、pooled phthalate weights、cycle-pooled strata/PSU 和已验证的 survey-logistic sandwich 实现。

- Primary continuous OR per MCOP doubling: **1.24551**
- 95% CI: **1.07753-1.43968**; P=0.00331136
- Primary N=9936, CRC cases=70
- Age >=40 OR: **1.22144**; 95% CI 1.04822-1.42329; P=0.0108405
- Q4 vs Q1 OR: **1.23392**; 95% CI 0.60745-2.50649; P-trend=0.0584221
- LOCO OR range: **1.19676-1.33488**; direction consistent: YES
- Quartile pattern is not strictly monotonic; Q4 is not a stronger contrast than Q3.

MCOP exposure availability and model outputs are reported separately below. This is a validation result, not evidence that MCOP is a DINP-specific causal exposure measure.

## Data audit

- MCOP variable: URXCOP; LOD comment variable: URDCOPLC
- Available cycles: 2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018
- Exposure + cancer outcome: 11086
- CRC cases: 76
- MCOP above LOD: 10907
- Official codebook: [CDC NHANES PHTHTE_J](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm)

## Prespecified analyses

1. Continuous log2(MCOP) - full primary cancer-free-control population.
2. Quartiles - Q1 reference, Q4 vs Q1, and linear trend.
3. Age >=40 - same covariate set.
4. Cancer-free controls - primary population definition.
5. Leave-one-cycle-out - one cycle removed at a time.

## Files

- mcop_crc_phase2_audit.csv
- mcop_crc_phase2_main_models.csv
- mcop_crc_phase2_quartiles.csv
- mcop_crc_phase2_leave_one_cycle_out.csv

Run timestamp (UTC): 2026-08-25T15:28:32.542537+00:00