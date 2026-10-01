# GALAXY 开发、独立公共验证与动态风险更新

## 分析目标与时间定义

开发队列为 GALAXY 2024；Chen、COSMOS 和 ColonAiQ 在模型冻结后分别验证。结局为临床复发，不含 ctDNA molecular recurrence。临床变量限定为分期、年龄、性别；GALAXY 患者级补充表没有肿瘤部位字段，因此核心模型不含 site。MRD 使用阳性/阴性，不拼接 VAF 或不同 assay 的连续值。

**主分析：名义L6 landmark。** GALAXY 原文将术后160–200天定义为6-month timepoint。本研究把该研究定义的 ctDNA 状态按名义术后6个月 day 182.6 作为 L6 输入，纳入 day 182.6 时无临床复发者，预测 **(182.6天, 365.24天]** 内临床复发。此为临床主分析。由于补充表2没有 GALAXY 患者级采血日期，day 182.6 是名义 landmark，采血实际日在160–200天内；该时间精度作为局限报告。day 200 起点的分析另列为时间窗敏感性分析。

L6按各研究自己的临床节点定义，不统一转换成同一个采血日：GALAXY采用原文160–200天研究定义窗（补充表无个体采血日，主分析风险起点按名义day182.6）；Chen按原文预设的术后6个月采样节点，在实际数据中选取4.5–7.5个月内最接近6个月的一次；COSMOS使用原始补充表中的scheduled month-6列；ColonAiQ按每3个月监测方案，选取4.5–7.5个月内最接近术后6个月的一次。Chen与ColonAiQ容差取季度采样间隔的一半，实际采样日作为该患者的风险起点。复发结局从各自L6评估之后计算至术后12个月。L3仍按70–112天、MRD window按14–<70天整理。研究采样方案分别见[GALAXY论文](https://doi.org/10.1038/s41591-024-03254-6)、[Chen等](https://d-nb.info/1241319898/34)、[COSMOS论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11443202/)及[ColonAiQ论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10119774/)。实际纳入人数及所选L6采样月份见 `mrd_window_harmonization_audit.csv`。

## GALAXY 开发

匹配的三时间点风险集：N=1081，(182.6天,12月]复发=103。ridge logistic 的惩罚系数 λ=1；年龄以开发集均值66.32岁和标准差11.29标准化。Clinical-only 与 Clinical+MRD 在完全相同患者上比较。MRD模型变量为 L6 当前阳性及既往已观测 MRD 状态中任一次阳性；GALAXY 可用的既往列为 MRD-window/L3。外部队列按其实际 L6 前采样历史派生同一“ever observed positive”概念，检测频率不同会影响该特征的观察机会。

| 模型 | apparent AUC | optimism-corrected AUC | apparent AP | corrected AP | apparent Brier | corrected Brier |
|---|---:|---:|---:|---:|---:|---:|
| Clinical-only | 0.658 | 0.645 | 0.141 | 0.131 | 0.084 | 0.084 |
| Clinical+MRD | 0.840 | 0.832 | 0.566 | 0.550 | 0.055 | 0.056 |

患者级 bootstrap 重抽样 1000 次。校准截距、斜率及OOB校准分箱见相应表；完整回归系数与标准化参数已冻结至 `outputs/models/galaxy_frozen_ridge_models.json`。没有使用外部验证结果调参、重校准或选阈值。

### GALAXY校准

| 模型 | apparent截距 | optimism-corrected截距 (95% bootstrap CI) | apparent斜率 | optimism-corrected斜率 (95% bootstrap CI) |
|---|---:|---:|---:|---:|
| Clinical-only | 0.195 | 0.055 (-0.599, 1.000) | 1.094 | 1.024 (0.689, 1.457) |
| Clinical+MRD | 0.176 | 0.128 (-0.308, 0.555) | 1.091 | 1.063 (0.890, 1.235) |

### day200 次要时间窗敏感性分析

day200结果来自此前完整运行并保留在原分析目录的1000次bootstrap；它作为时间窗敏感性结果，不替代day182.6主分析。两个landmark的风险集不同，数值用于评估起点选择的稳健性。

| 分析 | N | events | 模型 | corrected AUC (95% bootstrap CI) | corrected AP (95% bootstrap CI) | corrected Brier (95% bootstrap CI) |
|---|---:|---:|---|---:|---:|---:|
| 主分析 day182.6 | 1081 | 103 | Clinical-only | 0.645 (0.598, 0.692) | 0.131 (0.091, 0.164) | 0.084 (0.071, 0.098) |
| 主分析 day182.6 | 1081 | 103 | Clinical+MRD | 0.832 (0.791, 0.878) | 0.550 (0.446, 0.657) | 0.056 (0.045, 0.067) |
| 敏感性 day200 | 1070 | 92 | Clinical-only | 0.634 (0.586, 0.680) | 0.116 (0.077, 0.148) | 0.077 (0.064, 0.090) |
| 敏感性 day200 | 1070 | 92 | Clinical+MRD | 0.812 (0.762, 0.860) | 0.489 (0.375, 0.598) | 0.055 (0.044, 0.066) |

## 冻结模型的公共验证

表内顺序为 AUC / average precision / Brier；校准截距和斜率见 `external_validation_metrics.csv`。

| 队列 | 验证角色 | L6采样月，中位数（范围） | N | events | Clinical-only | Clinical+MRD |
|---|---|---:|---:|---:|---:|---:|
| COSMOS | 优先独立验证 | 6.000 (6.000–6.000) | 300 | 15 | 0.752 / 0.200 / 0.045 | 0.851 / 0.382 / 0.038 |
| Chen | 优先独立验证 | 6.374 (4.731–7.392) | 128 | 7 | 0.281 / 0.041 / 0.053 | 0.728 / 0.504 / 0.040 |
| ColonAiQ | 探索性（事件少） | 6.400 (4.600–7.500) | 85 | 7 | 0.617 / 0.271 / 0.078 | 0.756 / 0.452 / 0.074 |
| Pooled priority external | 优先队列合并 | 按各自L6节点 | 428 | 22 | 0.601 / 0.152 / 0.048 | 0.807 / 0.365 / 0.038 |

外部验证按COSMOS和Chen作为优先队列报告；ColonAiQ因可用事件少仅作探索性结果。纳入者须在各自L6评估前未复发、有可用既往MRD历史，并能判定术后12个月结局；事件要求发生在L6评估后、且不迟于术后12个月。各患者的L6到12个月预测时长会随实际评估日期略有不同，这保留了真实研究采样节奏。未做验证集特征选择、调参或重校准。MRD增量的配对bootstrap CI见 `clinical_vs_mrd_increment.csv`。样本量/事件少的队列其calibration slope与AP会不稳定，应结合N、events和区间解释。

| 队列 | 模型 | 校准截距 | 校准斜率 |
|---|---|---:|---:|
| COSMOS | Clinical-only | 0.656 | 1.432 |
| COSMOS | Clinical+MRD | -0.730 | 0.875 |
| Chen | Clinical-only | NA | NA |
| Chen | Clinical+MRD | 0.613 | 1.068 |
| ColonAiQ | Clinical-only | 6.041 | 2.512 |
| ColonAiQ | Clinical+MRD | -0.843 | 0.576 |

## GALAXY L6→L12 主轨迹结果

以下独立主表按**当前L6状态 + L6之前已观察到的MRD历史**划分持续阴性、转阳、清除和持续阳性。事件为名义day 182.6之后至术后12个月内临床复发；“平均预测风险/中位预测风险”来自GALAXY主模型的患者级bootstrap OOB预测。既往“无阳性”指纳入模型的既往可观测结果均无阳性，不代表未检测时间段的真实阴性。

| L6→L12轨迹 | 当前L6 MRD | 既往MRD历史 | N | 6–12月复发 | 观察复发率 (95% Wilson CI) | 平均OOB预测风险 | 中位OOB预测风险 |
|---|---|---|---:|---:|---:|---:|---:|
| 持续阴性 | 阴性 | 既往未观察到阳性 | 933 | 38 | 4.1% (3.0%–5.5%) | 4.5% | 5.2% |
| 转阳 | 阳性 | 既往未观察到阳性 | 22 | 15 | 68.2% (47.3%–83.6%) | 54.9% | 57.3% |
| 清除 | 阴性 | 既往至少一次阳性 | 80 | 14 | 17.5% (10.7%–27.3%) | 16.4% | 18.1% |
| 持续阳性 | 阳性 | 既往至少一次阳性 | 46 | 36 | 78.3% (64.4%–87.7%) | 77.3% | 81.8% |

可下载的独立汇总表为 `outputs/tables/galaxy_l6_l12_trajectory_main_table.csv`；患者级预测仍只保存在本地 `outputs/local_only/`。

## 分子记忆

在 L6 ctDNA 阴性者中，既往组定义为 L6 前至少一次已观测阳性；从未阳性组定义为所有已观测既往结果均阴性。

| 队列 | 既往状态 | N | 复发 | 观察复发率 | 平均预测风险 |
|---|---|---:|---:|---:|---:|
| GALAXY | never_positive | 933 | 38 | 4.1% | 4.5% |
| GALAXY | previously_positive_now_negative | 80 | 14 | 17.5% | 16.4% |
| Chen | never_positive | 115 | 2 | 1.7% | 1.9% |
| Chen | previously_positive_now_negative | 6 | 2 | 33.3% | 5.2% |
| ColonAiQ | never_positive | 59 | 2 | 3.4% | 1.9% |
| ColonAiQ | previously_positive_now_negative | 16 | 2 | 12.5% | 5.8% |
| COSMOS | never_positive | 271 | 6 | 2.2% | 3.2% |
| COSMOS | previously_positive_now_negative | 15 | 1 | 6.7% | 12.2% |

GALAXY 的分子记忆在 OOB 估计中有清楚分层；验证队列的“既往阳性、当前阴性”人数/事件很少，不能据此宣称已跨队列复制。

## 局限

1. GALAXY 的主分析将研究定义的160–200天状态映射到名义day 182.6 L6 landmark；实际采血日不可见，因此不能确认每个状态都在该名义日期前获得。day 200起点另作敏感性分析。
2. GALAXY 缺少患者级3/6月采血日，不能完成严格≤6月敏感性分析，也不能用真实日期拟合患者特异的动态生存模型。
3. GALAXY 没有个体 site，因此临床底模使用 stage、age、sex；无法检验加入部位后的效能。
4. 化疗时间、CEA、影像随访等没有在所有公共队列同定义，未进入公共 core model；模型估计关联与预测，不解释 MRD 清除的因果或治疗反应效应。
5. COSMOS 非复发者随访采用补充表 RFS/censor 月数；Chen/ColonAiQ 非复发者以可确认的随访采样时间判断12月结局。早期截尾或竞争事件不作为对照。

## 文件

系数、指标、增量、转移、分子记忆、缺失与时间窗审计表位于 `outputs/tables/`；指定图形位于 `outputs/figures/`。患者级轨迹及真实病例示例只写入 `outputs/local_only/`，不得提交到 GitHub。
