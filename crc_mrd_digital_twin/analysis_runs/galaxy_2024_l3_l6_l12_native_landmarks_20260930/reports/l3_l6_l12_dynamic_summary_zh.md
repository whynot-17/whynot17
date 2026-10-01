# 队列特异节点的 L3→L6→L12 动态分析

## 时间点定义

- **GALAXY**：L3采用论文定义的术后70–112天窗口；L6采用160–200天窗口，因补充表无患者级采血日期，L6→L12结局风险起点仍按预先确定的day182.6。
- **COSMOS**：L3和L6分别采用补充表中scheduled month-3和month-6节点；结局起点为患者的L6评估。
- **ColonAiQ**：按季度监测记录选择最接近术后3个月、且早于L6的实际样本作为L3；L6按最近6个月节点选择。每个患者的实际L6样本日作为结局起点。
- **Chen**：原文术后ctDNA采样在day3–7后从术后6个月开始，此后每3个月；没有研究定义的3个月MRD节点，因此不纳入L3→L6分析。少数非计划早期采样不重新标作L3。

各队列结局为：患者在其L6评估前无复发，并统计L6评估之后至术后12个月的临床复发。只有能够确认12个月结局者进入轨迹率与模型评价；未知ctDNA不编码为阴性。采样方案参见[GALAXY论文](https://doi.org/10.1038/s41591-024-03254-6)、[Chen队列原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8130394/)、[COSMOS论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11443202/)和[ColonAiQ论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10119774/)。

## L3/L6可用性与12个月结局风险集

| 队列 | L3节点 | 有L3 | 有L6 | 两节点均有 | L6后结局可判定 | 其中复发 |
|---|---|---:|---:|---:|---:|---:|
| GALAXY | published study-defined 70-112-day L3 window; patient-specific date absent | 1708 | 1396 | 1307 | 1115 | 106 |
| Chen | not available: published postoperative ctDNA schedule starts at month 6 | 0 | 154 | 0 | 0 | 0 |
| ColonAiQ | quarterly surveillance; closest observed test within 1.5 months of month 3 and before selected L6 | 61 | 137 | 61 | 40 | 3 |
| COSMOS | scheduled month-3 status column; exact 3-month node | 281 | 303 | 281 | 278 | 14 |

## L3→L6转移与L6后复发

| 队列 | 轨迹 | N | L6后至12月复发 | 复发率 (95% Wilson CI) | L3采样月中位数 | L6起点月中位数 |
|---|---|---:|---:|---:|---:|---:|
| GALAXY | L3阴→L6阴 | 1027 | 45 | 4.4% (3.3%–5.8%) | 研究定义窗 | 6.000 |
| GALAXY | L3阴→L6阳 | 30 | 21 | 70.0% (52.1%–83.3%) | 研究定义窗 | 6.000 |
| GALAXY | L3阳→L6阴 | 18 | 8 | 44.4% (24.6%–66.3%) | 研究定义窗 | 6.000 |
| GALAXY | L3阳→L6阳 | 40 | 32 | 80.0% (65.2%–89.5%) | 研究定义窗 | 6.000 |
| ColonAiQ | L3阴→L6阴 | 30 | 1 | 3.3% (0.6%–16.7%) | 3.700 | 6.700 |
| ColonAiQ | L3阴→L6阳 | 5 | 1 | 20.0% (3.6%–62.4%) | 3.900 | 6.900 |
| ColonAiQ | L3阳→L6阴 | 3 | 0 | 0.0% (0.0%–56.1%) | 3.800 | 6.100 |
| ColonAiQ | L3阳→L6阳 | 2 | 1 | 50.0% (9.5%–90.5%) | 3.250 | 5.500 |
| COSMOS | L3阴→L6阴 | 263 | 7 | 2.7% (1.3%–5.4%) | 3.000 | 6.000 |
| COSMOS | L3阴→L6阳 | 4 | 3 | 75.0% (30.1%–95.4%) | 3.000 | 6.000 |
| COSMOS | L3阳→L6阴 | 3 | 0 | 0.0% (0.0%–56.1%) | 3.000 | 6.000 |
| COSMOS | L3阳→L6阳 | 8 | 4 | 50.0% (21.5%–78.5%) | 3.000 | 6.000 |
| Pooled public (descriptive) | L3阴→L6阴 | 1320 | 53 | 4.0% (3.1%–5.2%) | 研究定义窗 | 队列特异 |
| Pooled public (descriptive) | L3阴→L6阳 | 39 | 25 | 64.1% (48.4%–77.3%) | 研究定义窗 | 队列特异 |
| Pooled public (descriptive) | L3阳→L6阴 | 24 | 8 | 33.3% (18.0%–53.3%) | 研究定义窗 | 队列特异 |
| Pooled public (descriptive) | L3阳→L6阳 | 50 | 37 | 74.0% (60.4%–84.1%) | 研究定义窗 | 队列特异 |

合并行只作描述性汇总；具体分析仍以各队列分层结果为主。不同队列assay、病例构成和采样节奏不同，不能把该合并比例解释成单一目标人群风险。

## same-horizon 预测更新

Clinical-only、Clinical+L3 history和Clinical+L6 history均在相同的L6后至12个月风险集上比较。GALAXY报告患者级bootstrap的OOB预测；外部队列使用冻结的GALAXY系数，不重拟合、不校准。Chen因缺少L3节点不参与这组比较。

| 队列 | 模型 | N | 事件 | AUC | AP | Brier |
|---|---|---:|---:|---:|---:|---:|
| GALAXY | Clinical-only | 1081 | 103 | 0.611 | 0.115 | 0.084 |
| GALAXY | Clinical+L3 history | 1081 | 103 | 0.753 | 0.383 | 0.064 |
| GALAXY | Clinical+L6 history | 1081 | 103 | 0.789 | 0.482 | 0.056 |
| COSMOS | Clinical-only | 278 | 14 | 0.732 | 0.199 | 0.046 |
| COSMOS | Clinical+L3 history | 278 | 14 | 0.764 | 0.250 | 0.047 |
| COSMOS | Clinical+L6 history | 278 | 14 | 0.838 | 0.368 | 0.039 |
| ColonAiQ | Clinical-only | 40 | 3 | 0.505 | 0.388 | 0.071 |
| ColonAiQ | Clinical+L3 history | 40 | 3 | 0.649 | 0.424 | 0.066 |
| ColonAiQ | Clinical+L6 history | 40 | 3 | 0.712 | 0.530 | 0.066 |
| Pooled external (descriptive) | Clinical-only | 318 | 17 | 0.684 | 0.177 | 0.049 |
| Pooled external (descriptive) | Clinical+L3 history | 318 | 17 | 0.743 | 0.239 | 0.050 |
| Pooled external (descriptive) | Clinical+L6 history | 318 | 17 | 0.813 | 0.363 | 0.042 |

| 队列 | 更新比较 | 指标 | 更新后−更新前 (95% paired bootstrap CI) |
|---|---|---|---:|
| GALAXY | Clinical-only -> Clinical+L3 history | auc | +0.142 (0.088, 0.198) |
| GALAXY | Clinical-only -> Clinical+L3 history | average_precision | +0.268 (0.180, 0.370) |
| GALAXY | Clinical-only -> Clinical+L3 history | brier | -0.020 (-0.028, -0.011) |
| GALAXY | Clinical+L3 history -> Clinical+L6 history | auc | +0.036 (0.002, 0.073) |
| GALAXY | Clinical+L3 history -> Clinical+L6 history | average_precision | +0.100 (0.015, 0.191) |
| GALAXY | Clinical+L3 history -> Clinical+L6 history | brier | -0.009 (-0.016, -0.002) |
| GALAXY | Clinical-only -> Clinical+L6 history | auc | +0.178 (0.121, 0.236) |
| GALAXY | Clinical-only -> Clinical+L6 history | average_precision | +0.367 (0.279, 0.482) |
| GALAXY | Clinical-only -> Clinical+L6 history | brier | -0.028 (-0.039, -0.020) |
| COSMOS | Clinical-only -> Clinical+L3 history | auc | +0.033 (-0.011, 0.095) |
| COSMOS | Clinical-only -> Clinical+L3 history | average_precision | +0.051 (-0.136, 0.304) |
| COSMOS | Clinical-only -> Clinical+L3 history | brier | +0.001 (-0.014, 0.013) |
| COSMOS | Clinical+L3 history -> Clinical+L6 history | auc | +0.073 (0.001, 0.187) |
| COSMOS | Clinical+L3 history -> Clinical+L6 history | average_precision | +0.118 (0.013, 0.269) |
| COSMOS | Clinical+L3 history -> Clinical+L6 history | brier | -0.008 (-0.019, 0.002) |
| COSMOS | Clinical-only -> Clinical+L6 history | auc | +0.106 (0.013, 0.226) |
| COSMOS | Clinical-only -> Clinical+L6 history | average_precision | +0.169 (-0.032, 0.432) |
| COSMOS | Clinical-only -> Clinical+L6 history | brier | -0.007 (-0.023, 0.009) |
| ColonAiQ | Clinical-only -> Clinical+L3 history | auc | +0.144 (-0.051, 0.513) |
| ColonAiQ | Clinical-only -> Clinical+L3 history | average_precision | +0.036 (-0.003, 0.144) |
| ColonAiQ | Clinical-only -> Clinical+L3 history | brier | -0.005 (-0.028, 0.008) |
| ColonAiQ | Clinical+L3 history -> Clinical+L6 history | auc | +0.063 (0.000, 0.256) |
| ColonAiQ | Clinical+L3 history -> Clinical+L6 history | average_precision | +0.106 (0.000, 0.460) |
| ColonAiQ | Clinical+L3 history -> Clinical+L6 history | brier | -0.001 (-0.056, 0.043) |
| ColonAiQ | Clinical-only -> Clinical+L6 history | auc | +0.207 (-0.051, 0.692) |
| ColonAiQ | Clinical-only -> Clinical+L6 history | average_precision | +0.142 (-0.003, 0.563) |
| ColonAiQ | Clinical-only -> Clinical+L6 history | brier | -0.005 (-0.074, 0.044) |
| Pooled external (descriptive) | Clinical-only -> Clinical+L3 history | auc | +0.059 (-0.011, 0.141) |
| Pooled external (descriptive) | Clinical-only -> Clinical+L3 history | average_precision | +0.062 (-0.100, 0.271) |
| Pooled external (descriptive) | Clinical-only -> Clinical+L3 history | brier | +0.000 (-0.012, 0.012) |
| Pooled external (descriptive) | Clinical+L3 history -> Clinical+L6 history | auc | +0.070 (0.003, 0.175) |
| Pooled external (descriptive) | Clinical+L3 history -> Clinical+L6 history | average_precision | +0.124 (0.030, 0.265) |
| Pooled external (descriptive) | Clinical+L3 history -> Clinical+L6 history | brier | -0.007 (-0.019, 0.003) |
| Pooled external (descriptive) | Clinical-only -> Clinical+L6 history | auc | +0.130 (0.026, 0.259) |
| Pooled external (descriptive) | Clinical-only -> Clinical+L6 history | average_precision | +0.186 (0.016, 0.451) |
| Pooled external (descriptive) | Clinical-only -> Clinical+L6 history | brier | -0.007 (-0.023, 0.007) |

上述同一结局窗比较检验在L3获得的信息基础上，加入L6状态是否进一步更新风险排序与误差。外部队列事件数有限，尤其ColonAiQ，AUC/AP及其差值区间应谨慎解读。L6 landmark前复发者被排除，因此本分析回答的是“L6时仍未复发者未来至12个月的风险”，不代表从L3起至12个月的总复发风险。

详细数据见 `outputs/tables/l3_l6_cohort_eligibility_audit.csv`、`l3_l6_l12_transition_summary.csv`、`l3_l6_l12_model_performance.csv` 与 `l3_l6_l12_model_update_increment.csv`。患者级预测文件仅保存于本地 `outputs/local_only/`。
