# Dynamic survival sensitivity 与时间泄漏审计

本轮完成了固定结局窗口的 same-horizon 更新，不将其称作连续时间 joint model 或 time-dependent Cox。

GALAXY 补充表提供DFS天数和研究标注的3月、6月 ctDNA 状态，但没有患者级3/6月采血日期。3月窗口为70–112天、6月窗口为160–200天；无法诚实地将患者级时间变 covariate 放在确切采血时点，也无法评估采血前后事件顺序。因此本轮**不拟合 time-dependent Cox 或 landmark supermodel**，不插值采血日。严格≤182.6天的 Galaxy 分析同样不可识别。主模型改用200天共同风险起点的保守 landmark，并明确为 L6-window-to-12m 近似。

## 主要泄漏检查

- L6 MRD 输入只来自 GALAXY 表中的研究定义6月列，外部队列只取160–200天内实测状态。
- 主风险集排除第200天及以前复发者；外部队列的采样日有记录，所选状态必须早于复发事件。
- 目标窗口为(200,365.24]天；复发早于或等于200天的患者不进入风险集；复发晚于365.24天只作为12月前无复发对照；非复发对照须有至少365.24天可确认随访。
- 未观测 ctDNA 不编码为阴性；临床结局早期截尾或随访不足者排除于拟合/验证集。
- 外部队列未用于特征选择、惩罚参数选择、阈值选择或重校准。
- 患者级同一结局窗风险更新用同一人群与配对OOB风险，避免将不同预测时距的风险直接相减。

## 留待未来院内数据

要估计真正患者特异的 `P(T_recurrence > u | MRD history through t)`，需院内队列同时保留确切 surgery date、每次采血日期、复发日期及末次临床无复发评估日期；届时优先用 landmark supermodel 或 time-dependent Cox，并报告 time-dependent AUC、C-index、Brier/IBS 与校准。当前公开GALAXY表的时间精度不支持这些指标，因此本轮不报告。
