# MCOP-CRC Phase 2B：稳定性审计

本轮不更换候选，只审计 MCOP 连续模型的稳健性：逐 cycle、global MCOP×cycle interaction、RCS、极端值排除和 creatinine-normalized exposure。

## Key findings

- Single-cycle effects above 1: 6/7 (2005-2006; 2007-2008; 2009-2010; 2013-2014; 2015-2016; 2017-2018). Below 1: 1/7 (2011-2012).
- The 2011-2012 single-cycle estimate is below 1, but its confidence interval is wide and includes substantial positive effects.
- Global MCOP×cycle interaction: F=3.21924, df=6,109, design-based P=0.00598088; chi-square P=0.00366274.
- RCS nonlinear test: design-based P=0.358283; chi-square P=0.354815.
- The RCS point curve does not show a stable Q3 peak; it rises toward the extreme upper tail, where uncertainty is large.
- Excluding top 1%: OR=1.27079, 95% CI 1.08307-1.49104, P=0.00364689.
- Excluding top 2.5%: OR=1.19358, 95% CI 1.01823-1.39912, P=0.0293862.
- Creatinine-normalized MCOP: OR=1.244, 95% CI 1.07474-1.43992, P=0.0037899.

## Interpretation

The continuous association is not explained by the top 1%-2.5% exposure tail or by the choice between creatinine covariate adjustment and creatinine-ratio normalization. However, the global interaction test indicates cycle heterogeneity, so the result is not yet a uniformly replicated effect across cycles. The RCS test does not support a stable nonlinear peak. Keep the candidate at yellow-green and prioritize an independent cohort rather than mechanistic expansion.

## 输出文件

- mcop_crc_phase2_per_cycle.csv
- mcop_crc_phase2_cycle_interaction.csv
- mcop_crc_phase2_spline.csv
- mcop_crc_phase2_spline_curve.csv
- mcop_crc_phase2_tail_exclusion.csv
- mcop_crc_phase2_creatinine_normalized.csv

RCS 使用 5th/35th/65th/95th percentile knots。非线性检验同时给出 design-based Wald F 近似和 chi-square 近似；主要参考 F 近似。
Creatinine-normalized 模型使用 log2(MCOP)-log2(creatinine)，因此不再把原始 creatinine_log2 作为同一模型协变量重复加入。

Run timestamp (UTC): 2026-08-22T11:36:17.261951+00:00