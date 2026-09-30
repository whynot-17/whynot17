# 多队列术后 CRC 12 个月复发动态预测：第一版模型结果

**日期：**2026-09-30  
**分析根目录：**`E:\crc_mrd_digital_twin`  
**结局：**术后12个月内临床复发（Y12）；landmark 模型在该时点仍无复发者中预测至12个月。

## 分析队列核对

| 队列 | Y12可判定 N | Y12事件 | 3月风险集 N | 3月后至12月事件 | 6月风险集 N | 6月后至12月事件 |
|---|---:|---:|---:|---:|---:|---:|
| Chen | 181 | 15 | 181 | 15 | 179 | 13 |
| ColonAiQ | 178 | 36 | 147 | 26 | 137 | 12 |
| COSMOS | 332 | 21 | 330 | 19 | 326 | 15 |
| Tie | 96 | 15 | 96 | 15 | 93 | 12 |
| **合计** | **787** | **87** | **754** | **75** | **735** | **52** |

脚本在拟合前逐队列核对了既有审计计数；若人数或事件数不符会中止。Henriksen 因没有逐人复发日期未进入 Y12 模型。

## 模型和验证

采用带 L2 收缩的低维 logistic 模型（惩罚参数 λ=1.0），避免在87个 Y12 事件上拟合过多参数。共同临床特征为 III/IV期与直肠原发；不同 assay 统一为二分类 MRD。主模型暂不纳入 ACT：补充表未给精确治疗时间，无法证明它在术前或 landmark 时已可用。

- **M0 baseline clinical：**所有 Y12 可判定者，预测 ≤12月复发。
- **L3 clinical only vs clinical + MRD history：**在术后3月仍无复发且截至该时点有可判读 MRD 历史者中，预测 (3,12] 月复发。
- **L6 clinical only vs clinical + MRD history：**同理预测 (6,12] 月复发。
- MRD 摘要为截至 landmark 最新状态，以及其前一次及更早样本中是否曾阳性。缺失采样绝不赋值为阴性；不使用 landmark 后/复发后的状态。
- 采用 leave-one-cohort-out (LOCO)：每轮整队列留出。留出队列的来源效应不可从训练折估计，因此 LOCO 传输模型不将 cohort ID 当预测器。另提供含 cohort 固定截距的全队列系数，仅适用于这四个已知来源；未来新医院需要独立校准。

不同 landmark 的风险集和预测结局不同，故 M0 全队列 AUC 与 L3/L6 AUC 不作直接的显著性或增益比较。同一 landmark 内，clinical-only 与 clinical+MRD 使用完全相同的患者和标签，描述性差异如下；尚未计算正式差异检验或置信区间。

| 配对比较 | ΔAUC | ΔAP | ΔBrier（负值较好）|
|---|---:|---:|---:|
| 3月风险集：MRD vs 临床 | AUC +0.182 | AP +0.209 | Brier -0.014 |
| 6月风险集：MRD vs 临床 | AUC +0.118 | AP +0.240 | Brier -0.008 |

## LOCO 汇总指标

AUC 为 ROC AUC；AP 为 average precision；校准截距/斜率来自交叉拟合预测的 logistic 校准回归。

| 模型 | N | 事件 | AUC | AP | Brier | 平均预测风险 | 校准截距 | 校准斜率 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| M0_baseline_clinical | 787 | 87 | 0.608 | 0.142 | 0.097 | 0.120 | -0.546 | 0.759 |
| L3_clinical_only | 754 | 75 | 0.626 | 0.137 | 0.087 | 0.110 | -0.188 | 0.965 |
| L3_clinical_plus_MRD_history | 754 | 75 | 0.808 | 0.346 | 0.073 | 0.099 | 0.223 | 1.122 |
| L6_clinical_only | 735 | 52 | 0.679 | 0.109 | 0.063 | 0.074 | 0.416 | 1.208 |
| L6_clinical_plus_MRD_history | 735 | 52 | 0.797 | 0.349 | 0.055 | 0.075 | 0.057 | 1.063 |

## 各留出队列表现

| 模型 | 留出队列 | N | 事件 | AUC | Brier |
|---|---|---:|---:|---:|---:|
| M0_baseline_clinical | COSMOS | 332 | 21 | 0.715 | 0.059 |
| M0_baseline_clinical | Chen | 181 | 15 | 0.629 | 0.077 |
| M0_baseline_clinical | ColonAiQ | 178 | 36 | 0.667 | 0.167 |
| M0_baseline_clinical | Tie | 96 | 15 | 0.500 | 0.133 |
| L3_clinical_only | COSMOS | 330 | 19 | 0.729 | 0.054 |
| L3_clinical_only | Chen | 181 | 15 | 0.629 | 0.076 |
| L3_clinical_only | ColonAiQ | 147 | 26 | 0.668 | 0.146 |
| L3_clinical_only | Tie | 96 | 15 | 0.500 | 0.132 |
| L3_clinical_plus_MRD_history | COSMOS | 330 | 19 | 0.794 | 0.048 |
| L3_clinical_plus_MRD_history | Chen | 181 | 15 | 0.799 | 0.061 |
| L3_clinical_plus_MRD_history | ColonAiQ | 147 | 26 | 0.872 | 0.105 |
| L3_clinical_plus_MRD_history | Tie | 96 | 15 | 0.653 | 0.128 |
| L6_clinical_only | COSMOS | 326 | 15 | 0.734 | 0.042 |
| L6_clinical_only | Chen | 179 | 13 | 0.647 | 0.066 |
| L6_clinical_only | ColonAiQ | 137 | 12 | 0.768 | 0.075 |
| L6_clinical_only | Tie | 93 | 12 | 0.500 | 0.112 |
| L6_clinical_plus_MRD_history | COSMOS | 326 | 15 | 0.849 | 0.035 |
| L6_clinical_plus_MRD_history | Chen | 179 | 13 | 0.771 | 0.057 |
| L6_clinical_plus_MRD_history | ColonAiQ | 137 | 12 | 0.904 | 0.057 |
| L6_clinical_plus_MRD_history | Tie | 93 | 12 | 0.586 | 0.117 |

小队列的事件数有限，单个留出队列指标会有较大抽样不确定性；本版为方法学开发结果，不能视为临床验证。

## 决策曲线

`outputs/tables/multicohort_decision_curve.csv` 提供风险阈值 5%、10%、15%、20%、30% 下的模型净获益，并与 treat-all / treat-none 基准比较。决策曲线依赖风险校准；对跨队列偏移明显的模型，应先校准再讨论临床阈值。

## 关键限制与解读

1. 可判定病例来自公开补充表的结局/末次随访字段，并沿用已审计的保守判定；结局未知者不作为阴性。
2. ColonAiQ 末次术后 ctDNA 样本仅作为无复发患者12个月观察证据的保守代理；它不等同于正式临床随访时间。
3. COSMOS 334人分析集按补充表重建：排除 R2，并按已发表纵向分析样本覆盖规则重建。该规则涉及纵向样本可用性，结果可能受观察/入组机制影响。
4. Tie 的术后 ctDNA 仅有 POM1 结果，公开表没有精确抽血日或 6月样本；6月模型在该队列仍使用 POM1 历史，没有假造新一次检测。
5. 仅4个队列，assay、来源、分期范围、随访和抽样规则同时变化。ColonAiQ 是唯一的甲基化 assay 队列，因而 assay 与 cohort 完全混杂，不能从这4组资料中独立估计 assay 效应。LOCO 是压力测试，不等同于独立前瞻性外部验证；AUC/校准估计仍需不确定区间和院内独立队列复核。
6. 单个留出队列事件数较少，且部分预测在该队列内几乎为常数；遇到校准回归不可识别或近似分离时，逐队列校准截距/斜率留空，不将不稳定极值作为结果。
7. 本版使用小型 ridge logistic 而不是复杂机器学习；当前主要价值是确认统一定义、可复现提取与跨队列验证流程是否成立。

## 生成文件

- `scripts/03_multicohort_landmark_model.py`：完整解析、风险集构建、LOCO、校准、Brier、AUC/AP、决策曲线脚本。
- `outputs/tables/multicohort_model_sample_counts.csv`：逐队列聚合人数核对。
- `outputs/tables/multicohort_loco_metrics.csv`：LOCO 汇总指标。
- `outputs/tables/multicohort_loco_by_heldout_cohort.csv`：逐留出队列指标。
- `outputs/tables/multicohort_decision_curve.csv`：净获益数据。
- `outputs/tables/multicohort_development_coefficients.csv`：全队列开发系数。
- 患者级交叉拟合预测仅保存在本机 `outputs/local_only/`，不会纳入版本库。
