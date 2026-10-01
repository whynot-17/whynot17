# Dynamic survival sensitivity 与时间泄漏审计

本轮完成了固定结局窗口的 same-horizon 更新，不将其称作连续时间 joint model 或 time-dependent Cox。

GALAXY 原文将术后160–200天称为6-month timepoint。本轮主分析按用户指定，以名义day 182.6作为 L6 landmark，预测至术后365.24天。补充表没有患者级3/6月采血日期，不能把个体时间变 covariate 放在确切采血时点，也不能检验采血与复发的逐患者先后顺序。因此本轮**不拟合 time-dependent Cox 或 landmark supermodel**，不插值采血日。day 200 作为次要敏感性分析，用于检验等待完整L6窗口结束的影响。

## 主要泄漏检查

- GALAXY主分析采用研究定义的160–200天L6状态和名义day182.6风险起点；每个外部队列按自己的study-defined/clinically designated L6节点使用对应结果，没有跨队列统一的日数截止。
- 各队列分别排除L6评估前或当日已复发者；预测结局是L6评估后至术后365.24天的临床复发。GALAXY个体采血日不可见，故只能使用名义day182.6；该研究窗内的实际事件与采血先后关系不能逐患者核验。
- 复发晚于术后365.24天只作为12月前无复发对照；非复发对照须有至少365.24天可确认随访。
- 未观测 ctDNA 不编码为阴性；临床结局早期截尾或随访不足者排除于拟合/验证集。
- 外部队列未用于特征选择、惩罚参数选择、阈值选择或重校准。
- 患者级同一结局窗风险更新用同一人群与配对OOB风险，避免将不同预测时距的风险直接相减。

## 留待未来院内数据

要估计真正患者特异的 `P(T_recurrence > u | MRD history through t)`，需院内队列同时保留确切 surgery date、每次采血日期、复发日期及末次临床无复发评估日期；届时优先用 landmark supermodel 或 time-dependent Cox，并报告 time-dependent AUC、C-index、Brier/IBS 与校准。当前公开GALAXY表的时间精度不支持这些指标，因此本轮不报告。
