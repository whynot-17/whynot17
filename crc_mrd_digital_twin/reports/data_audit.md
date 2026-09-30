# 数据审计报告

## 范围与数据来源

本报告记录 Chen 等人发表的 II-III 期结直肠癌前瞻性队列的第 0-1 阶段数据获取与审计结果（[论文 DOI](https://doi.org/10.1186/s13045-021-01089-z)）。从 Springer Nature Figshare 获取了与论文关联的全部 18 项补充材料，18/18 个文件均通过发布方提供的 MD5 校验。文件编号、DOI、网址、字节数及校验和保存在 data_raw/figshare_manifest.csv。

患者级数据来源包括：[补充文件 1](https://springernature.figshare.com/articles/journal_contribution/Additional_file_1_of_Postoperative_circulating_tumor_DNA_as_markers_of_recurrence_risk_in_stages_II_to_III_colorectal_cancer/14608859)（Table S1，DOCX）、[补充文件 6](https://springernature.figshare.com/articles/dataset/Additional_file_6_of_Postoperative_circulating_tumor_DNA_as_markers_of_recurrence_risk_in_stages_II_to_III_colorectal_cancer/14608874)（Dataset S1，原发肿瘤变异）、[补充文件 7](https://springernature.figshare.com/articles/dataset/Additional_file_7_of_Postoperative_circulating_tumor_DNA_as_markers_of_recurrence_risk_in_stages_II_to_III_colorectal_cancer/14608877)（Dataset S2，连续 ctDNA）及[补充文件 9](https://springernature.figshare.com/articles/dataset/Additional_file_9_of_Postoperative_circulating_tumor_DNA_as_markers_of_recurrence_risk_in_stages_II_to_III_colorectal_cancer/14608883)（Dataset S3，连续 CEA）。其余补充材料包括方法、基因面板、图和汇总结果。未对缺失的患者级数据进行填补。

## 队列与样本核对

- Table S1 中共有 **240** 名患者。
- 有复发时间记录者 **32** 人（II 期 10 人；III 期 22 人）。
- 连续 ctDNA 观测 **1,290** 条，覆盖 240 名患者。
- 每位患者采血次数（含 P1）：中位数 **5**，四分位距 **3-8**，范围 **2-9**。
- CEA 观测 **970** 条，覆盖 239 名患者。术前 CEA 有 231 条（231 名患者）；9 名患者缺少术前 CEA，其中 1 名患者任何时间都没有 CEA 结果。
- 原发肿瘤变异记录 **3,990** 条，覆盖 240 名患者。每位患者变异数中位数 **6**（四分位距 5-9；范围 1-327）。这些是肿瘤组织记录，不是血浆检测结果。

明细见 outputs/tables/patient_draw_counts.csv、samplepoint_summary.csv、tumor_variant_summary.csv 和 replication_checks.csv。

## 数据字典与字段可用性

字段级机器可读字典见 outputs/tables/data_dictionary.csv。Table S1 含 Patient_No、Age、Sex、PTL、Stage、Hist_type、Hist_grade、LV_invasion、N_invasion、MSI、ACT、R_time 和 R_treatment。R_time 是从手术到首次记录复发的月数；空值不代表末次随访时间。

Dataset S2 仅含 Patient_No、Time_postop（月）、Sample_point（P1-P9）及二元 ctDNA_status。没有唯一样本编号、血浆 ctDNA 定量值、血浆 VAF、血浆突变数或逐样本检测 QC。Dataset S1 是原发肿瘤突变级数据，含 DP 和 AD；这些测序读数不能解释为血浆 ctDNA 数值。Dataset S3 含 Patient_No、Time_postop（月）和 CEA_value。补充方法定义 CEA 正常范围为 0.00-5.00 ng/mL；本次审计将 >5.00 ng/mL 派生为升高。

没有单独的 T 分期和 N 分期、手术绝对日期、辅助化疗方案与时间、死亡状态/日期/死因，以及逐患者删失或末次随访时间。ACT 仅记录是/否。采样时间以十进制月表示，并非日期或天数。P1 标注为术前采血，但数值时间记录为 0；解释时须同时参考 Sample_point 与数值时间。

## 纵向结构与 MRD

各采样点的数量、阳性率及实际采样时间分布见 outputs/tables/samplepoint_summary.csv。所有记录中的阳性率为：P1 154/240（64.2%）；P2 20/240（8.3%）；P3 8/164（4.9%）；P4 9/140（6.4%）；P5 10/131（7.6%）；P6 8/113（7.1%）；P7 7/104（6.7%）；P8 2/79（2.5%）；P9 1/79（1.3%）。该表也报告已知复发前观测数及复发当日或之后的记录数。复发当日或之后的 ctDNA 属于结局发生后的信息，不可纳入预测未来复发的模型。

P1 为术前采血，P2 为早期术后采血。应保留实际记录的十进制月时间，不要将方案采样点换算为假定日期。转换频数见 outputs/tables/MRD_transition_frequencies.csv；分母仅包含相邻的已观测采血。

| 已观测历史范围 | 状态转换 | 配对数 |
|---|---|---:|
| 所有相邻的已观测采血 | 阴性→阴性 | 825 / 1050 |
| 所有相邻的已观测采血 | 阴性→阳性 | 26 / 1050 |
| 所有相邻的已观测采血 | 阳性→阴性 | 160 / 1050 |
| 所有相邻的已观测采血 | 阳性→阳性 | 39 / 1050 |
| 复发前相邻的术后采血 | 阴性→阴性 | 734 / 786 |
| 复发前相邻的术后采血 | 阴性→阳性 | 22 / 786 |
| 复发前相邻的术后采血 | 阳性→阴性 | 19 / 786 |
| 复发前相邻的术后采血 | 阳性→阳性 | 11 / 786 |

未观测到的采血绝不按 MRD 阴性处理。治疗与轨迹交叉计数使用 ACT 是/否字段；复发患者仅纳入复发前观测：

| ACT | 已观测术后轨迹 | 患者数 |
|---|---|---:|
| 否 | 末次观测转阳 | 4 |
| 否 | 持续阴性 | 42 |
| 否 | 末次观测转阴 | 1 |
| 否 | 仅一次术后观测：阴性 | 17 |
| 否 | 仅一次术后观测：阳性 | 2 |
| 是 | 间歇阴性；首末次均阳性 | 2 |
| 是 | 间歇阳性；首末次均阴性 | 7 |
| 是 | 末次观测转阳 | 8 |
| 是 | 持续阴性 | 118 |
| 是 | 持续阳性 | 3 |
| 是 | 末次观测转阴 | 8 |
| 是 | 仅一次术后观测：阴性 | 24 |
| 是 | 仅一次术后观测：阳性 | 4 |

术后采样间隔不规则。在至少有 3 次术后采血的患者中，患者内相邻采血间隔变异系数的中位数为 0.36（四分位距 0.30-0.45）；168 人中有 156 人超过描述性标记 CV 0.20。该标记不是临床阈值。按复发组别统计的采样覆盖率与采血次数见 outputs/tables/sample_coverage_by_recurrence.csv 和 draw_counts_by_recurrence.csv。复发患者采血次数中位数为 4（四分位距 3-6；n=32）；未记录复发患者中位数为 5（四分位距 3-8；n=208）。

## 复发、随访与信息性缺失

在复发患者中，R_time 中位数为 12.59 个月（四分位距 8.23-16.23）。这是复发患者的事件时间，不是全队列随访时间。未提供未记录复发患者的逐患者删失时间，也没有死亡字段。论文报告队列随访中位数为 27.4 个月（95% CI 26.2-28.5）；补充材料不足以独立重建复发组与未复发组的随访时长。

复发患者的复发前末次 MRD 观测，以及从最早一次术后 ctDNA 阳性到复发的探索性间隔，见 outputs/tables/recurrence_sample_timing_audit.csv。复发前末次采血距事件中位数为 2.89 个月（四分位距 1.36-6.00）。22 名复发患者可计算阳性至复发的代理间隔（均值 7.82、中位数 7.10 个月）。这不等同于论文定义的提前检出时间，因为缺少 ACT 日期，且论文使用根治性治疗后检出 ctDNA 的时间。此代理指标排除了复发当日或之后的样本；原始数据有 24 条观测发生在复发当日或之后，不能进入预测历史。

观测强度可能具有信息性。复发组与未复发组的后期访视覆盖率不同，计划采样并未完整实现，且未复发患者的删失/随访时长未知。现有文件无法区分漏访、提前失访、治疗、死亡、撤回或行政性删失。缺失访视未被编码为阴性；本次未拟合观测过程模型。

## 第 0-1 阶段结论

公开数据支持患者级临床表、带相对时间的二元连续 ctDNA、原发肿瘤变异记录及数值型连续 CEA 数据，但无法提供所有计划字段。血浆 MRD 定量值、逐样本 QC、详细 ACT 时间、个体删失时间及死亡数据均不可用。
