# MCOP–CRC Phase 2H：复杂抽样设计与标准 survey 软件复核

## Scope

本轮只重建 MCOP 覆盖的七个 NHANES 2-year cycle（2005–06 至 2017–18）。MBzP 的十-cycle 历史分析没有修改。病例定义、cancer-free controls、协变量处理和 MCOP 主暴露均沿用 Phase 2 冻结版本。

- Pooled weight rule: **cycle-specific phthalate subsample weight / 7**.
- Complete-case primary N=9936, CRC cases=70.
- Weight-source audit: **PASS**.
- Near-zero SAS missing-weight sentinels are recorded separately; primary complete-case modeling uses no such row.
- Primary complete-case singleton strata: **0**.

## Python `/10` versus `/7` check

- `/7`: OR=1.2455068, 95% CI 1.0775254–1.4396756, P=0.0033113585.
- Legacy `/10` on the same seven-cycle participants: OR=1.2455068, 95% CI 1.0775254–1.4396756, P=0.0033113585.
- Because the Python fitter normalizes all weights by their sample mean before fitting, `/10` versus `/7` is a common multiplicative rescaling and should not change beta, SE, OR or P.

## Independent R `survey::svyglm` gate

- R `svyglm`: OR=1.2455068, 95% CI 1.0773085–1.4399655, standard P=0.0033957984; design-df P=0.0033113352.
- Python: beta=0.2195425, SE=0.073096728; R: beta=0.2195425, SE=0.073096671.
- Relative logOR change: **4.879e-11%**; direction same: **True**; CI null conclusion same: **True**.
- R design df=109; R model residual df=97; Python design df=109.

## Frozen sensitivity audit

The sensitivity CSV contains LOCO, age ≥40, sex-specific effects and formal interaction, diagnosis-age exclusions (<1/<2/<5 years), top-tail exclusions, creatinine normalization, and pairwise co-exposure models. No mechanistic analysis was run in Phase 2H.

Weighted quartiles use the phthalate subsample weights only for cutpoint construction; the unweighted cutpoint analysis is retained unchanged. The RCS sensitivity reports both unweighted and survey-weighted 5th/35th/65th/95th-percentile knot sets.

## Decision

Phase 2H decision: **GREEN**. The gate is GREEN only when R `svyglm` remains positive, the relative logOR change is ≤10%, and the CI conclusion agrees with the Python implementation. This remains a cross-sectional association audit, not causal evidence.

## CDC design references

- [NHANES weighting tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx)
- [NHANES sample design and analysis](https://wwwn.cdc.gov/nchs/nhanes/tutorials/sampledesign.aspx)

## Outputs

- `mcop_crc_phase2h_weight_sources.csv`
- `mcop_crc_phase2h_design_units.csv`
- `mcop_crc_phase2h_python_vs_standard_survey.csv`
- `mcop_crc_phase2h_primary_reanalysis.csv`
- `mcop_crc_phase2h_sensitivity_reanalysis.csv`
- `mcop_crc_phase2h_weighted_quantile_sensitivity.csv`
- `mcop_crc_phase2h_survey_design_audit.md`