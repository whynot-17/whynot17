# 三队列 CRC MRD：18个月复发动态预测

**日期：**2026-09-30  
**队列：**Chen、ColonAiQ、COSMOS；Tie 不纳入。

## 18个月结局和 ColonAiQ 可用性

终点定义为每个 landmark 后、术后18个月内临床复发。阳性为 recurrence time ∈ (landmark, 18]；阴性为 recurrence time >18，或无复发记录且逐患者 follow-up proxy ≥18个月。仅用 landmark 当时及此前 MRD，缺失检测不作阴性填补。

| Landmark | 队列 | 有 MRD 历史且结局可判定 N | 事件 | 对照 | 风险集中18月结局未知 |
|---:|---|---:|---:|---:|---:|
| 3 | Chen | 155 | 25 | 130 | 68 |
| 6 | Chen | 153 | 23 | 130 | 65 |
| 12 | Chen | 140 | 10 | 130 | 26 |
| 3 | ColonAiQ | 70 | 40 | 30 | 174 |
| 6 | ColonAiQ | 56 | 26 | 30 | 158 |
| 12 | ColonAiQ | 55 | 16 | 39 | 87 |
| 3 | COSMOS | 326 | 23 | 303 | 6 |
| 6 | COSMOS | 322 | 19 | 303 | 6 |
| 12 | COSMOS | 307 | 4 | 303 | 4 |

ColonAiQ 在三个 landmark 均有阳性、阴性病例可用：L3 为 70 人（40/30），L6 为 56 人（26/30），L12 为 55 人（16/39）。因此可进入三队列18个月模型。ColonAiQ 的阴性观察证据使用补充表中末次术后采血时间作为随访代理，并非逐患者影像学随访时间，解释时需考虑该限制。

三队列 L12 风险集的 MRD 可判定样本合计 502 人、30 例事件。Chen 为 140/10，ColonAiQ 为 55/16，COSMOS 为 307/4。每个 LOCO 训练折均至少有10例事件；但 COSMOS 留出折只有4例事件，L12 的 COSMOS 特异性能估计会很不稳定，故该部分为探索性。未因事件少而事后移除 COSMOS。

## 预测设置与无信息泄漏

- L3 使用 ≤3月 MRD，预测 (3,18] 月复发；L6 使用 ≤6月 MRD，预测 (6,18] 月复发；L12 使用 ≤12月 MRD，预测 (12,18] 月复发。
- 两个模型分别为 clinical-only（stage、site）与 clinical + MRD（再加入最近一次 MRD 状态、当前样本前曾阳性）。
- ridge logistic，λ=1.0；不加入 cohort 指示变量，以便留出队列预测。
- 所有患者逐一断言 MRD 特征时间不超过对应 landmark。

## LOCO 表现与完整 bootstrap 95% CI

在每个队列×结局类别层内以患者为单位有放回抽样，每次重新拟合所有 LOCO 折和两类模型。AUC/AP 为三个留出队列折内指标的宏平均，Brier 为全部 out-of-fold 预测的汇总均值；每个 landmark 1000次，随机种子 20260930。差值为 clinical + MRD 减 clinical-only。

| Landmark | N | 事件 | AUC clinical | AUC + MRD [95% CI] | ΔAUC [95% CI] | AP clinical | AP + MRD | ΔAP [95% CI] | Brier clinical | Brier + MRD | ΔBrier [95% CI] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| L3 | 551 | 88 | 0.674 | 0.780 [0.713, 0.832] | 0.105 [0.053, 0.157] | 0.337 | 0.552 | 0.214 [0.091, 0.258] | 0.148 | 0.102 | -0.045 [-0.063, -0.027] |
| L6 | 531 | 68 | 0.687 | 0.784 [0.710, 0.844] | 0.097 [0.043, 0.154] | 0.285 | 0.531 | 0.246 [0.086, 0.294] | 0.117 | 0.096 | -0.021 [-0.034, -0.006] |
| L12 | 502 | 30 | 0.642 | 0.773 [0.611, 0.848] | 0.131 [0.012, 0.251] | 0.182 | 0.397 | 0.215 [0.024, 0.337] | 0.063 | 0.058 | -0.005 [-0.012, 0.003] |

## 各留出队列表现

| Landmark | 模型 | 留出队列 | N | 事件 | AUC | AP | Brier | 校准截距 | 校准斜率 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| L3 | clinical_only | COSMOS | 326 | 23 | 0.728 | 0.162 | 0.096 | -0.871 | 1.822 |
| L3 | clinical_only | Chen | 155 | 25 | 0.649 | 0.221 | 0.134 | -0.635 | 0.789 |
| L3 | clinical_only | ColonAiQ | 70 | 40 | 0.646 | 0.629 | 0.419 | 1.918 | 0.787 |
| L3 | clinical_plus_MRD | COSMOS | 326 | 23 | 0.789 | 0.408 | 0.066 | nan | nan |
| L3 | clinical_plus_MRD | Chen | 155 | 25 | 0.767 | 0.490 | 0.103 | 0.053 | 1.048 |
| L3 | clinical_plus_MRD | ColonAiQ | 70 | 40 | 0.783 | 0.757 | 0.268 | 1.491 | 0.902 |
| L6 | clinical_only | COSMOS | 322 | 19 | 0.731 | 0.148 | 0.074 | -0.731 | 1.644 |
| L6 | clinical_only | Chen | 153 | 23 | 0.661 | 0.187 | 0.122 | -0.358 | 0.825 |
| L6 | clinical_only | ColonAiQ | 56 | 26 | 0.668 | 0.521 | 0.351 | 1.857 | 0.910 |
| L6 | clinical_plus_MRD | COSMOS | 322 | 19 | 0.860 | 0.477 | 0.060 | nan | nan |
| L6 | clinical_plus_MRD | Chen | 153 | 23 | 0.744 | 0.440 | 0.105 | 0.563 | 1.154 |
| L6 | clinical_plus_MRD | ColonAiQ | 56 | 26 | 0.747 | 0.677 | 0.279 | 1.193 | 0.772 |
| L12 | clinical_only | COSMOS | 307 | 4 | 0.712 | 0.145 | 0.025 | 1.649 | 3.212 |
| L12 | clinical_only | Chen | 140 | 10 | 0.669 | 0.093 | 0.065 | 2.291 | 1.979 |
| L12 | clinical_only | ColonAiQ | 55 | 16 | 0.545 | 0.308 | 0.270 | 0.040 | 0.270 |
| L12 | clinical_plus_MRD | COSMOS | 307 | 4 | 0.932 | 0.395 | 0.025 | nan | nan |
| L12 | clinical_plus_MRD | Chen | 140 | 10 | 0.708 | 0.300 | 0.059 | 0.962 | 1.256 |
| L12 | clinical_plus_MRD | ColonAiQ | 55 | 16 | 0.679 | 0.496 | 0.244 | 0.976 | 0.619 |

## 动态风险轨迹

同时具有 L3、L6、L12 cross-fitted 风险预测的 L12 患者数为 491。三个时点均预测到同一个术后18个月终点，风险随新 MRD 信息更新。

| MRD 状态模式 | 人数 | (12,18]月事件 | 中位 Risk3 | 中位 Risk6 | 中位 Risk12 |
|---|---:|---:|---:|---:|---:|
| intermittent_mixed | 9 | 5 | 0.216 | 0.407 | 0.515 |
| negative_to_positive_conversion | 6 | 1 | 0.152 | 0.130 | 0.155 |
| persistent_negative | 418 | 12 | 0.106 | 0.095 | 0.073 |
| persistent_positive | 5 | 1 | 0.769 | 0.900 | 0.630 |
| positive_to_negative_clearance | 32 | 7 | 0.256 | 0.407 | 0.228 |
| single_negative | 17 | 1 | 0.058 | 0.043 | 0.032 |
| single_positive | 4 | 1 | 0.356 | 0.219 | 0.053 |

患者级明细、预测、典型病例和患者轨迹图仅保存在 outputs/local_only/，未提交至 GitHub。

## 结果解释

该模型可在三个公开队列中计算，ColonAiQ 在18个月终点下可纳入。但 ColonAiQ 的随访判定来自末次采血时间代理，且 L12 总事件只有30例、COSMOS 留出折仅4例事件。结果属于多队列开发和 internal-external LOCO 评估，不是独立外部验证或临床效用证明。

## 输出

- outputs/tables/l18_cohort_outcome_audit.csv
- outputs/tables/l18_loco_metrics.csv
- outputs/tables/l18_loco_by_heldout_cohort.csv
- outputs/tables/l18_landmark_bootstrap_ci.csv
- outputs/tables/l18_trajectory_summary.csv
- outputs/figures/l18_auc_comparison.svg
- outputs/figures/l18_brier_comparison.svg
- 患者级输出保存在本地 outputs/local_only/
