# MCOP-CRC Phase 2C：cycle heterogeneity audit

本轮只解释 MCOP cycle heterogeneity，不更换候选、不增加机制分析。所有人口学和 CRC 结构描述使用 MCOP 有效且 CRC outcome 为 CRC 或 cancer-free 的 primary population；模型型 quartile case counts 使用与前一轮相同的完整协变量 complete-case frame。

## 主要发现

- 2011-2012 的 primary MCOP median=16 ng/mL，Q1-Q3=6.6-45.35；CRC cases=11，CRC median age=73.
- 2011-2012 MCOP above LOD=100.00%，codebook LLOD=0.2 ng/mL；因此 2011-2012 的反向 OR 不是由最高 LOD/censoring 直接解释。
- 在与主模型相同的 complete-case frame 中，2011-2012 CRC cases=10，病例 MCOP median=9.8 ng/mL，对照 median=16.1 ng/mL，病例/对照 median ratio=0.609；这与该周期未调整 OR<1 的方向一致。
- Pooled primary MCOP median=7.8 ng/mL.
- 公开代码本的方法描述均属于 HPLC-ESI-MS/MS 类平台，但这不能证明跨周期完全没有校准/批次尺度差异。可见的 documented assay-scale change 是 MCOP LLOD：2005-2008 为 0.7 ng/mL，2009-2012 为 0.2，2013-2018 为 0.3。
- CDC explicitly states no lab method, equipment, or site changes for 2013-2014, 2015-2016, and 2017-2018; the 2011-2012 manual is Method 6306.04 and explicitly maps MCOP to the di-isononyl phthalate axis.
- 因此当前最合理的解释不是单一 LOD 故障，而是 2011-2012 暴露分布整体偏高、该周期仅 10 个 complete-case CRC cases，且病例与对照的暴露排序反向；年龄、种族、吸烟和肌酐构成差异可能进一步改变调整后效应。仅凭当前数据不能把其中任何一个因素定为唯一原因。

## Quartile case counts

| Quartile | N | CRC cases | Controls | Unweighted CRC % | Weighted CRC prevalence |
| --- | ---: | ---: | ---: | ---: | ---: |
| Q1 | 2523 | 16 | 2507 | 0.634% | 0.0038324 |
| Q2 | 2447 | 14 | 2433 | 0.572% | 0.00545015 |
| Q3 | 2479 | 21 | 2458 | 0.847% | 0.00857738 |
| Q4 | 2487 | 19 | 2468 | 0.764% | 0.00650035 |

Quartile complete-case frame: N=9936, CRC cases=70; log2 cutpoints=-2.8365;1.76553;2.96347;4.41278;11.6827.

## Files

- mcop_crc_phase2_cycle_heterogeneity_summary.csv
- mcop_crc_phase2_assay_lod_audit.csv
- mcop_crc_phase2_quartile_case_counts.csv

官方依据：[NHANES 2005-2006 PHTHTE_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PHTHTE_D.htm)、[2011-2012 PHTHTE_G](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PHTHTE_G.htm)、[2011-2012 laboratory manual](https://wwwn.cdc.gov/nchs/data/nhanes/public/2011/labmethods/phthte_g_met.pdf)、[2013-2014 PHTHTE_H](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/PHTHTE_H.htm)、[2015-2016 PHTHTE_I](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PHTHTE_I.htm)、[2017-2018 PHTHTE_J](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm)。

Run timestamp (UTC): 2026-08-25T15:32:09.169140+00:00