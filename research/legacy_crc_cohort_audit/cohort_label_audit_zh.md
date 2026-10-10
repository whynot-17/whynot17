# 三套旧 CRC 队列临床标签审计

**审计日期：** 2026-10-10
**范围：** 只审计 GSE39582、GSE28702、GSE122246 的公开样本规模、治疗、结局标签、随访与既往研究。没有下载或打开表达矩阵；没有进行差异表达、通路富集、机器学习或细胞比例分析。统计单位按患者/样本对应患者，而不是细胞。

## 结论先行

- **GSE39582：有条件保留为主要临床发现队列。** 它有术后辅助化疗状态和复发随访，适合检验“治疗前原发瘤表达的候选细胞状态是否与接受辅助治疗后的复发相关”。但文献中的 stage II–III 化疗子集从 **142 到约 210 例**不等，具体可用复发事件数也没有在公开汇总表中按“stage II–III 且接受化疗”交叉列出。正式分析前必须固定一个临床注释版本并在人级别重新核实分母、事件数和时间字段。
- **GSE28702：保留为 FOLFOX 应答的方向性验证。** 83 名不可切除 CRC 患者有 42 应答、41 不应答；肿瘤样本均在 mFOLFOX6 前采集，治疗 4 个周期后用 CT 评估。它不是术后 MRD/复发队列，且样本混合原发瘤和转移灶，不能把应答标签解释成“残留细胞”。
- **GSE122246：只作低权重的外部方向核对。** 24 名 stage III 患者全部接受辅助化疗，按极端结局抽取 12 名 5 年内未复发者和 12 名 3 年内复发者。它不是连续入组队列，也不适合多变量调整；其原研究已经报告 PNN、KCNQ1OT1 和一组差异转录信号。
- 三队列都没有纵向 ctDNA/MRD 测量，也没有化疗后的配对肿瘤样本。它们只能为 MRD/残留细胞假说提供**基线状态—后续临床结局的间接证据**，不能证明某种细胞真的在化疗后存活下来。bulk 芯片/RNA-seq 也不能直接证明肿瘤、巨噬细胞和 CAF 的空间共存或邻近关系。

## 队列逐项审计

### GSE39582：可作发现队列，但先锁定病例和事件口径

- GEO 的 GSE39582 Series 共 **585 个样本**，其中原论文描述为 443 个发现集结肠癌、123 个验证集结肠癌和 19 个非肿瘤黏膜；常见后续文献将 **566 个肿瘤表达谱**作为患者肿瘤队列。2024 年一篇辅助化疗签名研究汇总 stage II **264** 例、stage III **205** 例（合计 469）；较早的复发标志物重分析则报告 **460** 个 stage II–III 且有复发信息的肿瘤谱，若再要求复发信息和突变资料则为 **417**。这些数字反映不同过滤/资料完备性口径，不能直接互换。[GEO GSE39582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582)；[Marisa et al. 2013 原研究](https://pubmed.ncbi.nlm.nih.gov/23700391/)；[stage II–III 复发及突变信息的筛选示例](https://pmc.ncbi.nlm.nih.gov/articles/PMC5862619/)
- 结局与随访：2024 年汇总表对 stage II–III 子集报告中位随访 **48 个月**（SD 40.4，范围 0–201），复发 **139**、未复发 **322**、结局缺失 **8**；这组复发数对应全 stage II–III 汇总，不是“接受化疗者”内部的事件数。其他已发表清理版本也使用不同分母和结局口径，例如全 RFS 集的中位随访 43 个月、另一临床注释版本报告 52 个月。[RFS 分析示例](https://pmc.ncbi.nlm.nih.gov/articles/PMC4741632/)；[另一临床注释版本](https://www.dovepress.com/prognostic-value-of-the-combination-of-microsatellite-instability-and--peer-reviewed-fulltext-article-CMAR)。目标分析需要的是“stage II–III + 治疗前肿瘤 + 明确辅助治疗 + 可用 RFS 时间/事件”，该交叉集合的事件数尚不能从目前可查的汇总表可靠确定。[Xu et al. 2024，Table 1 与方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC11384724/)
- 辅助治疗病例数有口径差异：一个 2024 分析在 GSE39582 中使用 **142** 个接受 5-FU-based 辅助治疗的原发肿瘤；其他 stage II–III 治疗结局分析分别使用约 **202** 或 **210** 个接受治疗且符合其结局/协变量筛选的病例。治疗定义、stage、缺失值处理和随访可用性需要逐项核对，当前不应把任何一个数当成统一的最终样本数。[Xu et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11384724/)；[scRNA 整合签名研究中的治疗子集](https://pmc.ncbi.nlm.nih.gov/articles/PMC9604003/)；[按治疗状态分析 RFS 的研究](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.620811/full)
- 临床变量：公开研究可使用年龄、性别、AJCC stage、肿瘤部位/左右侧、MSI/MMR、KRAS/BRAF 等信息，并有辅助治疗状态。CMS 并非 GEO Series 页面上已核实的统一临床标签；后续若纳入 CMS，需先确定一致的表达分类器和病例映射，不能把 Marisa 的六类 C1–C6 直接叫作 CMS。[GSE39582 原始六亚型研究](https://pubmed.ncbi.nlm.nih.gov/23700391/)；[GSE39582 的临床变量和复发分析示例](https://pmc.ncbi.nlm.nih.gov/articles/PMC7951964/)
- 旧知识撞车：该队列早已用于分子亚型/预后、CMS 与辅助治疗获益、复发风险签名及微环境特征分析；2024 年又纳入 GSE39582 建立并验证 stage II–III 辅助化疗 **10 基因获益签名**。所以创新点不能只是“复发相关表达分数”或给旧预后基因换一个 MRD 名称；需要预先定义细胞来源明确的状态假说，再检验其跨结局一致性。[CMS 与化疗结局分析](https://ascopubs.org/doi/10.1200/PO.17.00241)；[2024 年 10 基因签名](https://pmc.ncbi.nlm.nih.gov/articles/PMC11384724/)
- 生物学解释边界：这是切除时/辅助化疗前原发肿瘤状态与后续复发的关系，不是化疗后残留肿瘤的直接采样。即使某状态与复发相关，也只能称为“与后续失败相关的基线状态”。

### GSE28702：直接 FOLFOX 应答标签，但临床场景不同

- GEO 记录 **83 名不可切除 CRC 患者**，42 responders、41 non-responders；原设计按入组时间分为训练集 54 和测试集 29。[GEO GSE28702](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28702)
- 原研究说明样本均在 mFOLFOX6 前采集，患者此前未接受化疗或放疗；病例由 **56 个原发 CRC 和 27 个转移灶**组成（肝 23、肺 1、腹膜 3）。均接受 mFOLFOX6，4 个周期后以 CT 按该研究使用的 RECIST 规则评估 CR/PR/SD/PD，再归纳应答与不应答。原作者已经用随机森林筛选应答基因并建立 15-probe-set 预测器。[Tsuji et al. 原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC3251854/)
- 结局质量：有明确、接近治疗过程的影像学反应标签，是三队列中唯一适合直接核对 FOLFOX 不应答的队列；但没有术后 RFS/ctDNA/MRD 标签。原论文也报告了依据 RF 预测分组的生存曲线，不等于公开了标准化的患者级复发时间结局。
- 混杂/调整：必须至少记录并敏感性检查取材为原发还是转移、转移部位、肿瘤部位和同期治疗。全队列均不可切除，常规 stage 可变性有限；当前公开 Series 摘要没有证实患者级 MSI/CMS 标签可用，不能承诺对 MSI/CMS 调整。样本量只有 83，建议预先指定一个主要效应检验，避免再次训练新分类器。
- 撞车：数据最初就是随机森林 FOLFOX 应答研究，之后又被其他 FOLFOX 基因预测和机器学习工作重复使用；也有研究用它评估 Treg/免疫浸润与反应。因此它是有用但并非未触碰的验证集，候选分数必须在看本队列结果前冻结。[Tsuji et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3251854/)；[后续 FOLFOX 基因分析](https://pmc.ncbi.nlm.nih.gov/articles/PMC6149834/)；[Treg 与治疗反应分析](https://pmc.ncbi.nlm.nih.gov/articles/PMC8899991/)

### GSE122246：极小型、极端结局 RNA-seq 复核

- GEO 有 **24 个 stage III 肿瘤样本**；原研究设计为接受辅助化疗后的极端结局队列。公开原始研究摘要显示不良组为化疗后 3 年内复发（n=12），良好组为化疗后至少 5 年无复发（n=12）；常规治疗包括 5-FU/capecitabine，部分方案含 oxaliplatin。它不是所有连续入组病例，也没有中间结局患者。[GEO GSE122246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122246)；[Mini et al. 原研究摘要](https://pubmed.ncbi.nlm.nih.gov/30973654/)；[两组的具体极端结局定义和样本数](https://www.annalsofoncology.org/article/S0923-7534%2819%2956692-3/fulltext)
- 原研究已报告 108 个差异基因，并提出/验证 PNN 与 KCNQ1OT1 和辅助化疗后 DFS 的关系；其结论属于既有分子机制/标志物结果，不能把这些已知信号包装成新的残留细胞发现。[Mini et al. 2019](https://pubmed.ncbi.nlm.nih.gov/30973654/)
- 可做：只对事先锁定的候选细胞状态分数做方向性核对，并报告每组 12 人及不确定性。不可做：在 24 人中训练模型、同时调整 stage/site/MSI/CMS、或把极端抽样标签当作完整生存队列来估计总体风险。该队列 stage 固定为 III，但公开记录/原摘要没有证明可用于统一调整的 MSI/CMS/左右侧临床表已公开。

## 对两个预注册问题的审计判断

| 问题 | 目前判断 | 可支持的结论边界 |
|---|---|---|
| 同一候选状态是否同时与辅助化疗后复发及 FOLFOX 不应答相关？ | **可以做有限的跨队列方向性检验。** GSE39582 提供辅助治疗后复发结局；GSE28702 提供 4 周期后的影像学应答；GSE122246 可作小样本极端结局复核。 | 三种终点、疾病阶段、样本位置和平台不同；结果一致是跨队列间接支持，不是 ctDNA 持续阳性或 MRD 机制证明。GSE122246 由于已发表/极端选样，不是完全独立的新机制发现。 |
| 控制 stage、部位、MSI/CMS 后是否仍有关？ | **GSE39582 最有可能，但事件数和协变量缺失必须在人级别核对；GSE28702 及 GSE122246 不足以完整回答。** | GSE28702 是不可切除病例且原发/转移混合，MSI/CMS 未在当前公开摘要中证实；GSE122246 n=24 不适合多变量模型。CMS 需预先确定分类方法并由表达重算/可靠映射。 |

## 建议的 go/no-go 条件

1. **先做标签整合，不先碰表达矩阵。** 为 GSE39582 固定一个表型表版本，建立 GEO GSM—患者 ID—stage—取材部位—辅助治疗方案—RFS 时间—复发事件—年龄/性别—MSI/MMR—CMS 映射；输出完整纳入/排除流程及化疗子集事件数。对照已发表 n=142、202、210 的差异，解释所用口径。
2. 如果 GSE39582 的化疗组能提供足够 RFS 事件，作为主要临床发现；早期复发阈值应预先锁定，优先用完整的时间-事件资料，并把长期无复发者定义为有足够随访的删失病例。不要把“没有记录到复发”直接等于长期无复发。
3. GSE28702 仅用于一个冻结候选分数与 response/non-response 的人级别比较，预先按原发/转移取材做敏感性检查；不再在此拟合/优化机器学习预测器。
4. GSE122246 仅作为小样本方向性补充。任何显著性分析都需呈现宽置信区间和 12/12 的组别选择限制。
5. 若不能取得可追溯的患者级 clinical annotation，或 GSE39582 的“治疗组复发事件数”低到无法支持预设模型，则**暂停表达分析**。只有临床标签审计通过后再启动第一轮表达评分；现阶段不进入空间组学/湿实验。

## 资料来源

1. [GSE39582 GEO 记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582)；[Marisa et al., 2013 原研究](https://pubmed.ncbi.nlm.nih.gov/23700391/)。
2. [GSE39582 stage II–III 复发/突变样本筛选](https://pmc.ncbi.nlm.nih.gov/articles/PMC5862619/)；[2024 年 GSE39582 10-gene ACT/RFS 分析](https://pmc.ncbi.nlm.nih.gov/articles/PMC11384724/)；[化疗人群的独立已发表筛选示例](https://pmc.ncbi.nlm.nih.gov/articles/PMC9604003/)。
3. [GSE28702 GEO 记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28702)；[Tsuji et al., 2011/2012 原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC3251854/)。
4. [GSE122246 GEO 记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122246)；[Mini et al., 2019 原研究](https://pubmed.ncbi.nlm.nih.gov/30973654/)；[原始会议摘要中的两组定义](https://www.annalsofoncology.org/article/S0923-7534%2819%2956692-3/fulltext)。
5. [GSE39582 的 CMS/辅助化疗既往分析](https://ascopubs.org/doi/10.1200/PO.17.00241)；[GSE28702 后续预测基因分析](https://pmc.ncbi.nlm.nih.gov/articles/PMC6149834/)；[GSE28702 Treg/应答分析](https://pmc.ncbi.nlm.nih.gov/articles/PMC8899991/)。

---

**本审计的决定：** 保留 GSE39582 + GSE28702 做第一轮候选状态的临床标签可行性验证；GSE122246 降为可选补充。进入任何表达层分析前，先解决 GSE39582 stage II–III 治疗组的病例映射、具体复发事件数和随访截尾口径。结果均是临床间接证据，不命名为 MRD signature，也不声称直接识别了化疗后存活的细胞。
