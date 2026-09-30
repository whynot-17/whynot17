# CRC MRD 纵向预测 digital twin：三队列纵向模型

**日期：**2026-09-30  
**项目：**E:/crc_mrd_digital_twin

## Digital twin 定义

预测型 digital twin 用患者基线临床特征初始化纵向疾病状态，并在获得新的术后 MRD 检测后更新。在每个 landmark，只使用该时点及此前可获得的信息预测后续复发风险。

形式化表示为 Z_t = f(X_baseline, MRD_{≤t})，Risk_t = P(T_recurrence ∈ (t,H_t] | Z_t)。本轮仅纳入可跨队列统一的 stage、site 和 assay-agnostic 二分类 MRD 状态；ACT 给药时间、CEA 及肿瘤分子基线无法在公开队列间可靠统一，未纳入模型。

## 风险集与时间窗

- L3：只使用采样时间≤3月的信息，预测术后 (3,12] 月复发。
- L6：只使用采样时间≤6月的信息，预测术后 (6,12] 月复发。
- L12：只使用采样时间≤12月的信息，预测术后 (12,24] 月复发。
- 三个风险对应不同的未来时间窗。MRD 特征包括当前状态及既往检测是否曾阳性。缺失检测不编码为阴性。
- 每位患者特征所用 MRD 时间均通过断言检查，必须≤landmark；L12 的任何>12月 MRD 检测均会过滤，若有泄漏脚本将中止。

## L12 结局可判定性与队列纳入

L12 结局审计见 reports/12m_landmark_24m_outcome_audit_zh.md。ColonAiQ 有19例在 (12,24] 月复发，但12月风险患者中没有可确认随访至24月的未复发对照；123人的24月结局未知，因此整个 ColonAiQ 队列从 L12 模型和 LOCO 排除。L12 主模型仅使用 Chen 与 COSMOS：N=382、事件=35（Chen 101人/16事件；COSMOS 281人/19事件）。L3/L6 仍使用 Chen、ColonAiQ、COSMOS 三队列。L12 按队列整体排除，不会从结局不完整的队列单独拼入事件病例。

## 模型与验证

低维 ridge logistic，λ=1.0。每个时间点分别拟合 clinical-only 和 clinical+MRD 模型；逐一留出整个队列。L3/L6 为三折 LOCO；L12 因 ColonAiQ 结局无法判定而只有两个队列，每轮用一个队列训练、另一个队列留出。这是有限的 internal-external transportability assessment，不是独立外部验证。

## 完整 LOCO 患者分层 bootstrap 95% CI

在每个队列内以患者为单位有放回抽样；每次对 clinical-only、clinical+MRD 模型以及每个 LOCO fold 重新拟合。同一 bootstrap replicate 用于配对差异。重复次数 1000，随机种子 20260930。差值均为 MRD 模型减 clinical-only；ΔBrier 为负表示 MRD 模型 Brier 更低。区间为 percentile 95% CI，成功/失败次数见 outputs/tables/landmark_bootstrap_ci.csv。AUC、AP、Brier 差值未进行额外假设检验。

| Landmark | N | 事件 | AUC clinical | AUC + MRD | ΔAUC | Brier clinical | Brier + MRD | ΔBrier |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| L3 | 658 | 60 | 0.638 [0.558, 0.699] | 0.828 [0.489, 0.872] | +0.190 [-0.140, +0.270] | 0.081 [0.064, 0.100] | 0.064 [0.052, 0.083] | -0.017 [-0.028, -0.004] |
| L6 | 642 | 40 | 0.698 [0.611, 0.753] | 0.830 [0.457, 0.875] | +0.133 [-0.230, +0.180] | 0.056 [0.041, 0.072] | 0.046 [0.035, 0.065] | -0.010 [-0.015, +0.000] |
| L12 | 382 | 35 | 0.487 [0.329, 0.728] | 0.603 [0.365, 0.770] | +0.116 [-0.083, +0.193] | 0.085 [0.062, 0.117] | 0.082 [0.059, 0.113] | -0.003 [-0.011, +0.006] |

## LOCO 点估计和校准

| Landmark | 模型 | 留出队列 | N | 事件 | AUC | Brier | 平均预测风险 | 校准截距 | 校准斜率 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| L3 | clinical_only | COSMOS | 330 | 19 | 0.729 | 0.054 | 0.106 | 0.606 | 1.696 |
| L3 | clinical_only | Chen | 181 | 15 | 0.629 | 0.077 | 0.119 | -1.044 | 0.655 |
| L3 | clinical_only | ColonAiQ | 147 | 26 | 0.668 | 0.148 | 0.081 | 1.120 | 1.084 |
| L3 | clinical_plus_MRD | COSMOS | 330 | 19 | 0.794 | 0.049 | 0.065 | nan | nan |
| L3 | clinical_plus_MRD | Chen | 181 | 15 | 0.799 | 0.061 | 0.091 | 0.072 | 1.105 |
| L3 | clinical_plus_MRD | ColonAiQ | 147 | 26 | 0.872 | 0.103 | 0.137 | 0.979 | 1.438 |
| L6 | clinical_only | COSMOS | 326 | 15 | 0.734 | 0.042 | 0.064 | 0.370 | 1.316 |
| L6 | clinical_only | Chen | 179 | 13 | 0.647 | 0.066 | 0.070 | -0.787 | 0.650 |
| L6 | clinical_only | ColonAiQ | 137 | 12 | 0.768 | 0.075 | 0.064 | nan | nan |
| L6 | clinical_plus_MRD | COSMOS | 326 | 15 | 0.849 | 0.035 | 0.062 | -0.268 | 1.071 |
| L6 | clinical_plus_MRD | Chen | 179 | 13 | 0.771 | 0.057 | 0.055 | 0.388 | 1.017 |
| L6 | clinical_plus_MRD | ColonAiQ | 137 | 12 | 0.904 | 0.057 | 0.087 | 0.599 | 1.346 |
| L12 | clinical_only | COSMOS | 281 | 19 | 0.693 | 0.066 | 0.143 | 1.087 | 2.164 |
| L12 | clinical_only | Chen | 101 | 16 | 0.620 | 0.138 | 0.067 | 0.463 | 0.794 |
| L12 | clinical_plus_MRD | COSMOS | 281 | 19 | 0.730 | 0.064 | 0.145 | -0.430 | 1.307 |
| L12 | clinical_plus_MRD | Chen | 101 | 16 | 0.676 | 0.132 | 0.066 | 1.434 | 1.161 |

校准截距和斜率在每个留出队列单独计算；若某折近似分离或不可识别则记为 NaN。汇总性能应与各留出队列表现一起解读，尤其是 L12 两队列的小样本结果。

## Digital twin 风险轨迹

对同时进入 L3、L6、L12 LOCO 的患者，保存 cross-fitted Risk_3m → Risk_6m → Risk_12m。三个风险各自对应不同的未来时间窗，图和数据表均按窗口标注。

| 截至12月 MRD 状态模式 | 三时点人数 | 12–24月事件 | 中位 Risk3 | 中位 Risk6 | 中位 Risk12 |
|---|---:|---:|---:|---:|---:|
| intermittent_mixed | 6 | 3 | 0.073 | 0.170 | 0.478 |
| negative_to_positive_conversion | 4 | 1 | 0.060 | 0.039 | 0.163 |
| persistent_negative | 338 | 23 | 0.031 | 0.019 | 0.100 |
| persistent_positive | 3 | 1 | 0.521 | 0.681 | 0.561 |
| positive_to_negative_clearance | 18 | 5 | 0.073 | 0.121 | 0.271 |
| single_negative | 12 | 2 | 0.027 | 0.013 | 0.029 |
| single_positive | 1 | 0 | 0.228 | 0.117 | 0.058 |

典型患者图展示实际 MRD 采样、landmark 预测和临床复发时间。患者级 ID 与轨迹仅保存在本地 outputs/local_only/，患者轨迹 SVG 也仅保存在 outputs/figures/digital_twin_patient_trajectories.svg，不会提交到 GitHub。若某种模式在可判定病例中不存在，图中不会虚构病例。

## 解释边界

本工作属于 multicohort development + internal-external validation。LOCO 评估跨公开队列的 transportability，不等于独立院外前瞻性验证；clinical utility 尚未通过临床决策实验或前瞻性研究验证。独立 external validation 仍预留给院内队列。L12 仅覆盖 Chen/COSMOS 两个来源，估计不确定性较大，应视为探索性分析。

## 输出文件

- outputs/tables/l12_cohort_audit.csv
- outputs/tables/l12_loco_metrics.csv
- outputs/tables/l12_loco_by_heldout_cohort.csv
- outputs/tables/landmark_bootstrap_ci.csv
- outputs/tables/digital_twin_trajectory_summary.csv
- outputs/figures/landmark_auc_comparison.svg
- outputs/figures/landmark_brier_comparison.svg
- 患者级轨迹数据和典型病例仅保存在 outputs/local_only/，患者轨迹 SVG 也不会提交到 GitHub。
