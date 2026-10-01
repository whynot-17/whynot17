# GALAXY 开发、独立公共验证与动态风险更新

## 分析目标与时间定义

开发队列为 GALAXY 2024；Chen、COSMOS 和 ColonAiQ 在模型冻结后分别验证。结局为临床复发，不含 ctDNA molecular recurrence。临床变量限定为分期、年龄、性别；GALAXY 患者级补充表没有肿瘤部位字段，因此核心模型不含 site。MRD 使用阳性/阴性，不拼接 VAF 或不同 assay 的连续值。

**时间泄漏保护：** GALAXY 的 study-defined 6个月采样窗口为术后160–200天，但表2没有个体采血日。主分析将共同 landmark 放在第200天，以确保整段 GALAXY 采样窗口已结束；纳入者需在第200天仍未复发，预测结局为 **(200天, 365.24天]** 临床复发。这是对临床 L6→12m 问题的保守近似，并非严格的(182.6天,12月]。对 GALAXY 严格 ≤182.6天敏感性分析无法识别，不能臆造采血日。

其他队列的主分析状态取各自160–200天窗口内最接近182.6天的实测结果；L3 取70–112天，MRD window 取14–<70天。所有结局使用同一第200天起点，且事件必须晚于起点。窗口和可用人数见 `mrd_window_harmonization_audit.csv`。

## GALAXY 开发

匹配的三时间点风险集：N=1070，(200天,12月]复发=92。ridge logistic 的惩罚系数 λ=1；年龄以开发集均值66.32岁和标准差11.34标准化。Clinical-only 与 Clinical+MRD 在完全相同患者上比较。MRD模型变量为 L6 当前阳性及既往已观测 MRD 状态中任一次阳性；GALAXY 可用的既往列为 MRD-window/L3。外部队列按其实际 L6 前采样历史派生同一“ever observed positive”概念，检测频率不同会影响该特征的观察机会。

| 模型 | apparent AUC | optimism-corrected AUC | apparent AP | corrected AP | apparent Brier | corrected Brier |
|---|---:|---:|---:|---:|---:|---:|
| Clinical-only | 0.650 | 0.634 | 0.127 | 0.116 | 0.077 | 0.077 |
| Clinical+MRD | 0.821 | 0.812 | 0.507 | 0.489 | 0.054 | 0.055 |

患者级 bootstrap 重抽样 1000 次。校准截距、斜率及OOB校准分箱见相应表；完整回归系数与标准化参数已冻结至 `outputs/models/galaxy_frozen_ridge_models.json`。没有使用外部验证结果调参、重校准或选阈值。

## 冻结模型的公共验证

表内顺序为 AUC / average precision / Brier；校准截距和斜率见 `external_validation_metrics.csv`。

| 队列 | N | events | Clinical-only | Clinical+MRD |
|---|---:|---:|---:|---:|
| Chen | 79 | 5 | 0.170 / 0.045 / 0.061 | 0.630 / 0.459 / 0.048 |
| COSMOS | 295 | 10 | 0.721 / 0.070 / 0.033 | 0.849 / 0.216 / 0.034 |
| ColonAiQ | 28 | 1 | 0.278 / 0.048 / 0.035 | 0.204 / 0.043 / 0.049 |
| Pooled external | 402 | 16 | 0.521 / 0.057 / 0.038 | 0.731 / 0.198 / 0.038 |

外部队列按各自窗口、相同结局起点与同一套冻结系数评估；未做验证集特征选择、调参或重校准。MRD 增量的配对 bootstrap CI 见 `clinical_vs_mrd_increment.csv`。评估数据来自既有公开补充材料，样本量/事件少的队列其 calibration slope 与 AP 会不稳定，应结合 N、events 和区间解释。

## 分子记忆

在 L6 ctDNA 阴性者中，既往组定义为 L6 前至少一次已观测阳性；从未阳性组定义为所有已观测既往结果均阴性。

| 队列 | 既往状态 | N | 复发 | 观察复发率 | 平均预测风险 |
|---|---|---:|---:|---:|---:|
| GALAXY | never_positive | 933 | 38 | 4.1% | 4.5% |
| GALAXY | previously_positive_now_negative | 79 | 13 | 16.5% | 15.6% |
| Chen | never_positive | 72 | 2 | 2.8% | 1.8% |
| Chen | previously_positive_now_negative | 4 | 1 | 25.0% | 4.9% |
| ColonAiQ | never_positive | 19 | 1 | 5.3% | 1.9% |
| ColonAiQ | previously_positive_now_negative | 7 | 0 | 0.0% | 5.5% |
| COSMOS | never_positive | 270 | 5 | 1.9% | 3.2% |
| COSMOS | previously_positive_now_negative | 15 | 1 | 6.7% | 11.7% |

GALAXY 的分子记忆在 OOB 估计中有清楚分层；验证队列的“既往阳性、当前阴性”人数/事件很少，不能据此宣称已跨队列复制。

## 局限

1. 第200天起点使时间泄漏风险最低，但排除了第182.6–200天复发；结论不能直接称为精确6月至12月复发风险。
2. GALAXY 缺少患者级3/6月采血日，不能完成严格≤6月敏感性分析，也不能用真实日期拟合患者特异的动态生存模型。
3. GALAXY 没有个体 site，因此临床底模使用 stage、age、sex；无法检验加入部位后的效能。
4. 化疗时间、CEA、影像随访等没有在所有公共队列同定义，未进入公共 core model；模型估计关联与预测，不解释 MRD 清除的因果或治疗反应效应。
5. COSMOS 非复发者随访采用补充表 RFS/censor 月数；Chen/ColonAiQ 非复发者以可确认的随访采样时间判断12月结局。早期截尾或竞争事件不作为对照。

## 文件

系数、指标、增量、转移、分子记忆、缺失与时间窗审计表位于 `outputs/tables/`；指定图形位于 `outputs/figures/`。患者级轨迹及真实病例示例只写入 `outputs/local_only/`，不得提交到 GitHub。
