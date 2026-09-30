# 公开/可申请的术后 CRC serial ctDNA/MRD 队列审计

**检索截止：2026-09-30**  
**目的：**在正式训练模型前，盘点能否取得患者级临床结局及纵向 ctDNA/MRD 结果；重点判断其是否支持术后 12 个月动态复发预测。  
**范围说明：**这是基于主要原始论文、补充材料和数据可用性声明的定向审计，不是 PRISMA 系统综述。这里的“公开”指可下载的患者级/样本级表格；“申请”指论文明确允许向作者、研究委员会或数据平台提交申请。能下载补充表不自动代表允许任何二次使用，须逐篇核对许可。

## 结论

1. **当前最值得先整理的可下载队列有 4 批：**
   - Chen：240 人、1,290 份纵向 ctDNA 样本；补充 Dataset S2/S3 分别有患者级 ctDNA 和 CEA 结果。
   - 丹麦/西班牙：最终分析 160 人、1,204 份血浆 ctDNA 样本；公开补充表提供逐患者纵向结果。
   - 中国 ColonAiQ 甲基化队列：299 人、1,228 份血样；补充 eTable 3 标明涵盖全部样本的 ctDNA 检测结果。
   - COSMOS-CRC-01：334 人、1,902 份纵向样本；Supplementary Table 2 及数据声明提供处理后的 ctDNA 检出状态和患者级临床资料。
2. **Tie 2019 的 96 人队列适合作为短程动态补充**：公开逐患者术后/化疗后检测结果和结局，但没有同等密度的季度监测。
3. **不要把 GALAXY 直接放进开发集。**尽管 Supplementary Table 2 含去标识逐患者资料与 ctDNA 窗口结果，CIRCULATE-Japan 的声明把访问用途限定为科学核验，并明确禁止二次使用。除非拿到允许模型开发的书面许可，否则应排除出训练及 LOCO 验证。
4. **不能把文章报告的队列人数简单相加当成可用模型样本量。**Reinert 2019 与 Henriksen 2022 至少有 77 人重叠；GALAXY 的多个子研究也不是独立队列。各队列 12 个月分析人数还要在下载后逐患者按观察时间、复发事件和 landmark 资格重新计算。
5. **12 个月结局尤其要防止把随访不足的未复发者错标为对照。**例如 Chen 补充表缺少未复发者逐患者末次随访/删失时间；只可把有证据显示无复发随访达到 12 个月的患者作为 12-month control。3 月和 6 月预测还要排除 landmark 前已经复发者，并确认 12 月前结局状态。

## 数据访问分层

### A. 可下载患者级/样本级纵向结果：优先核验和整理

| 队列 | 样本规模与时间结构 | 已公开的患者级信息 | 适用性与限制 |
|---|---|---|---|
| Chen 等，2021，中国 II–III 期 | 240 人；1,290 份血浆；术后早期及约每 3 个月监测 | Supplementary Dataset S2：逐次 ctDNA 状态；S3：CEA；另有临床特征与肿瘤变异表。[论文及补充材料](https://pmc.ncbi.nlm.nih.gov/articles/PMC8130394/) | 核心候选。ctDNA 是二元结果；缺少未复发患者逐人末次随访时间，不能把 R_time 空值当 12 月阴性。项目已有的 Chen 数据审计报告记录了这些字段限制。 |
| Henriksen 等，2022，丹麦/西班牙 III 期 | 入组 168 人，最终分析 160 人；1,204 份血浆；随访期约每 3 个月 | 多份公开补充工作簿含逐患者 ctDNA 动态和临床结局。[论文与补充表](https://pmc.ncbi.nlm.nih.gov/articles/PMC9401484/) | 最强的季度纵向候选之一。论文指出其中 77 人曾在早期研究中报告；与 Reinert 2019 队列做患者级去重，不能按 160+125 相加。 |
| Mo 等，2023，ColonAiQ 多中心甲基化 | 299 人；1,228 份血样；术后、治疗中及约每 3 个月，最长至 2 年 | Supplement 1 的 eTable 3 标为全部 1,228 份样本的 ctDNA 结果，含相对采样时间与患者结局信息。[论文与补充表](https://pmc.ncbi.nlm.nih.gov/articles/PMC10119774/) | 值得纳入多 assay 开发。为甲基化 qPCR 二元状态，不宜直接合并 assay 强度。中位随访约 21 个月、范围 8–27 个月，12 月对照仍须按逐患者实际观察确认。 |
| Nakamura 等，2024，COSMOS-CRC-01 | 334 人；1,902 份纵向 surveillance 样本；tissue-free epigenomic assay | Supplementary Table 2 为临床与检测数据；作者声明处理后的 ctDNA detected/not-detected 及患者级临床数据见补充材料。[论文与补充表](https://pmc.ncbi.nlm.nih.gov/articles/PMC11443202/) | 高价值外部 assay 队列。只用可分享的二元处理结果，原始测序数据不公开；将其作为独立研究来源并保留 cohort 标签。 |
| Tie 等，2019，澳大利亚 III 期 | 96 人；174 份 serial plasma；主要为术后及化疗后两个窗口 | eTable 3 提供逐患者术后/化疗后 ctDNA、CEA、复发结局。[论文及补充表](https://jamanetwork.com/journals/jamaoncology/fullarticle/2752788) | 可补 clearance/持久阳性/清除后复发等短程轨迹；季度 surveillance 信息少，不应视作高密度纵向队列。 |
| Reinert 等，2019，丹麦 I–III 期 | 125 人；795 份样本；术后、治疗及 surveillance 阶段 | eTable 5 列出 125 人详细 ctDNA 动态。[论文及补充表](https://pmc.ncbi.nlm.nih.gov/articles/PMC6512280/) | 仅在 Henriksen 2022 数据有缺项时补充；与 Henriksen 2022 有至少 77 人已发表重叠，先用中心/入组期/采样时间和可用 ID 去重。 |

### B. 有患者级材料，但需许可或不能直接当训练数据

| 队列 | 可用规模 | 访问状态 | 建议 |
|---|---:|---|---|
| CIRCULATE-Japan / GALAXY，2024 | 主队列 2,240 人、13,429 份样本；补充表含 MRD、3/6 月及 surveillance 结果 | 文件可见，但声明限定科学核验用途并禁止二次使用。[论文与数据声明](https://pmc.ncbi.nlm.nih.gov/articles/PMC11564113/) | 先书面询问是否允许独立模型开发；未获准前排除。2026 年 GALAXY tissue-free assay 验证是该注册人群的子集，不能另算独立外部队列。 |
| UK TRACC Part B | 214 人、约 997 份连续血样 | 补充材料以分组表和图为主；患者级资料需向作者合理申请。[论文与补充材料](https://pmc.ncbi.nlm.nih.gov/articles/PMC11325146/) | 若希望加入英国 tissue-free assay，按 12 月动态预测方案向作者申请患者×样本长表及事件/末次随访。 |
| α-CORRECT，2025 | 124 人进入分析；1,029 个 ctDNA 结果；中位随访 4.8 年 | Supplementary Information 有患者 swimmer plot；Data Availability Statement 说明患者级数据可向通讯作者申请。[论文与数据声明](https://onlinelibrary.wiley.com/doi/10.1002/jso.27989) | 强申请候选：纵向每 3 个月至 3 年，之后每 6 个月至 5 年。公开图形不能替代机器可读的逐次检测表。 |
| FIND III 期试验，2026 | 随机试验；mITT 分析 584 人 | 季度甲基化 ctDNA surveillance 有公开论文，患者级数据需按研究数据共享流程申请。[论文](https://ascopubs.org/doi/10.1200/JCO-25-03009) | 可作为请求获得后的多中心敏感性/外部验证；试验治疗策略会影响结局与采样，需预先定义 estimand。 |
| ALTAIR，CIRCULATE-Japan | 243 人左右，ctDNA 阳性富集人群 | 个体数据需申请/审查；人群来自 GALAXY/CIRCULATE-Japan 体系。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13375530/) | ctDNA 阳性富集且与 GALAXY 人群关联；不作为独立普通术后队列，也不能与 GALAXY 重复计数。 |
| JGAS000590 日本小队列 | 38 人；31 人有纵向样本 | 原始数据在日本 JGA 仓库，需按其流程申请访问。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9909342/) | 规模小、使用条件需先确认，作为外部敏感性候选。 |

### C. 公开论文/图表可用，但当前不是结构化核心数据

- **Sasaki 等，2025 dPCR 密集监测队列：**52 人，867 个检测时间点、1,526 份 plasma aliquots，平均每人 16.7 个时间点；论文补充材料主要提供逐患者轨迹图，尚未确认有完整机器可读患者×样本表。它非常值得向作者询问原始长表，但不应从图上数字化后冒充原始数据。[论文与补充材料](https://pmc.ncbi.nlm.nih.gov/articles/PMC11932597/)
- **Jin 等，2021 mqMSP：**82 人、182 份样本，但只有 19 人获得额外纵向随访样本；论文展示患者轨迹图并公开部分测序数据 accession，未发现完整 patient-by-visit 检测长表。[论文与补充材料](https://pmc.ncbi.nlm.nih.gov/articles/PMC7865146/)
- **Tie 等，2016 II 期队列：**230 人、1,046 份血浆；公开材料的纵向部分主要覆盖复发病例，不能把复发病例的密集随访当成全队列数据。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC5346159/)
- **Parikh 等，2021 plasma-only 队列：**约 103 人、252 份 prospectively collected plasma，但公开补充表并未确认包含完整逐患者、逐时间点检测长表；先向作者申请或仅用于文献层面比较。[论文与补充材料](https://aacrjournals.org/clincancerres/article/27/20/5586/671719/Minimal-Residual-Disease-Detection-using-a-Plasma)
- **中国 309 人 NGS 队列：**包含 I–IV 期，术后样本约 573 份、平均每人 1.9 份；发表材料以汇总分析为主，未确认有完整患者级纵向数据文件。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10822076/)
- **DYNAMIC-III：**约 968 人，主要 ctDNA 测量是术后 5–6 周单一 landmark；更适合 landmark 结果的外部对照，不足以支持多轮动态更新；数据另有申请要求。[论文](https://www.nature.com/articles/s41591-025-04030-w)
- **肝转移切除（CRLM）队列：**可用于 metastatic-resection 专项敏感性分析，但既往 IV 期病史及围手术期治疗和局限期 II–III 期目标人群不同，应与主队列分层，不能直接混为同一风险集。[代表性前瞻研究](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1003620)

## 对 12 个月动态模型的直接影响

建议把预测目标写成 landmark 条件风险，而不是将所有患者静态标注为“12 个月复发/未复发”：

- 在术后 landmark \(L\in\{3,6\}\) 月，限定患者截至 \(L\) 尚未复发；用截至该时点的 MRD/CEA 历史估计 \(P(T_{rec}\leq 12\text{月}\mid T_{rec}>L,\;History(L))\)。
- 复发病例应有可靠的手术至复发时间；对照必须能证明至少随访至术后 12 月且期间未复发。若队列只有总体中位随访、没有患者级随访终点，不能将该患者标成对照。
- 保留真实的采血日/术后相对时间、检测方法、治疗前后状态、队列来源及缺测标记。缺测不是 MRD 阴性；MRD assay 的“阳性/阴性”含义可 harmonize，定量单位和阈值不能假设跨 assay 等价。
- 预先统一 recurrence endpoint、手术起始日、ACT 状态/结束日与末次无事件随访；对“复发事件”与 DFS/RFS 包含死亡的复合终点作明确区分。
- 用多个数据源 pooled development，按 cohort 做 leave-one-cohort-out；避免随机患者切分造成同一中心/同一 assay 同时出现在训练与验证中。每轮独立按完整队列留出，报告 cohort-specific calibration 与 discrimination。
- 在尚未申请到受限数据之前，可先从 Chen、Henriksen、ColonAiQ、COSMOS 建长表；Tie 2019 用于 clearance/两阶段敏感性；最后用院内数据做真正独立 external validation。

**实施顺序：**先对 5 个可下载队列逐一读取表头和患者ID结构，再计算各自 12 月事件数、达到 12 月观察的对照数，以及 3/6 月 landmark 的在险数；这一步之后再定 pooled/LOCO 的实际样本配置。