# GALAXY 2024 患者级 MRD 数据可用性审计

**项目：** `E:/crc_mrd_digital_twin`  
**审计日期：** 2026-09-30  
**数据：** Nature Medicine 2024 Supplementary Tables 1–3 PDF；仅使用其中 Supplementary Table 2，原始 PDF 留在本地 `data_raw/`，不提交 GitHub。

## 六项核心数字

| 指标 | 人数 | 口径 |
|---|---:|---|
| 总患者数 | 2240 | 表2中连续编号且患者行解析成功 |
| 有 Recurred 标志和 DFS 天数 | 2240 | 两字段均非缺失；其中任意随访复发 514 例 |
| 有 study-defined 3-month ctDNA 状态 | 1708 | 阳性 183；阴性 1525 |
| 有 study-defined 6-month ctDNA 状态 | 1396 | 阳性 119；阴性 1277 |
| 同时有 3m 与 6m 状态 | 1307 | 两个名义时间点均有阳性/阴性结果 |
| 其中 12m 复发终点可判定 | 1164 | 复发事件 155；明确非复发 1009 |

另外，**MRD window + 3m + 6m 三个状态均齐全 1261 人**，可构造基线 MRD 历史状态。

## 按当前项目终点定义估计可用风险集

采用术后天数作为近似边界：3个月 = 91.3 天，6个月 = 182.6 天，12个月 = 365.2 天。12个月复发终点按复发日期判定：复发≤365.24天为事件；复发>365.24天或无复发且 DFS≥365.24天为明确非复发；无复发但 DFS<365.24天为未知。死亡在复发终点中作为竞争事件/未知，不编码为无复发。

| 分析集 | n | 12个月内复发 | 明确无复发 | 注 |
|---|---:|---:|---:|---|
| L3：有3m状态、Y12可判定且DFS≥91.3天 | 1516 | 251 | 1265 | 名义3m landmark |
| L6：有6m状态、Y12可判定且DFS≥182.6天 | 1185 | 113 | 1072 | 名义6m landmark |
| 3m→6m 配对轨迹，L3 风险集 | 1162 | 153 | 1009 | 两状态都齐全 |
| 3m→6m 配对轨迹，L6 风险集 | 1115 | 106 | 1009 | 两状态都齐全 |
| MRD window+3m+6m 完整，L6 风险集 | 1081 | 103 | 978 | 具备完整既往状态 |

现有三队列模型 L3 为 658 人/60 事件，L6 为 642 人/40 事件，L3→L6 配对为 638 人。若采用 GALAXY 名义时间点，L3 端可与既有集拼接至约 **2174 人**；L6 端约 **1827 人**；配对 L6 轨迹约 **1753 人**。这些只是可用性规模估算，尚未训练或合并模型。

## 关键时间对齐限制

论文方法定义 3-month ctDNA 为术后 **70–112 天**，6-month ctDNA 为 **160–200 天**，但 Supplementary Table 2 只提供两个状态列，没有每位患者的 3m/6m 实际采血日。MRD-window 日数有单独字段。严格的术后≤3.0月（91.3天）和≤6.0月（182.6天）可用人数因此**无法从公开表精确识别**：两个采样窗都跨越我们固定的landmark边界。上表按 DFS 日数筛选的是结局与风险时间的可用性近似，不能替代采血日审计。

所以 GALAXY 非常值得加入，但应先二选一：

1. 取得逐患者 3m/6m 采血日，再按现有严格≤3/≤6月规则重审；或
2. 将 landmark 预先改为研究定义的 70–112天 / 160–200天采样窗，并把 Chen、ColonAiQ、COSMOS 也重做同一时间规则。

在完成这个决定前，不把这些行拼进当前冻结 L3/L6 预测，不声称严格跨队列验证。论文提到的 1,664 人是同时有 MRD-window 与后续 surveillance-window 结果的子集，并非 3m+6m trajectory 的样本数。

## 数据与抽取校验

表2解析出 2,240 个唯一连续患者行；已校验 MRD-window 阳性/阴性状态共 2110 人、3m状态 1708 人、6m状态 1396 人。没有患者行或 ID 写入汇总输出。

原始补充 PDF：`E:/crc_mrd_digital_twin/data_raw/galaxy_2024_supplementary_tables_1_3.pdf`  
审计脚本：`scripts/08_galaxy_patient_level_audit.py`  
汇总表：`outputs/tables/galaxy_patient_level_availability_audit.csv`、`galaxy_timepoint_availability.csv`、`galaxy_nominal_mrd_state_transitions.csv`

## 来源

- [Nature Medicine 2024 论文](https://www.nature.com/articles/s41591-024-03254-6)
- [补充表 PDF](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41591-024-03254-6/MediaObjects/41591_2024_3254_MOESM1_ESM.pdf)
