# 环境毒理学与 CRC：Phase 1 分析协议

## 当前状态

本文件记录已经审定的 **CTD × GeneCards 环境化学物富集筛选方案**。

当前只完成方法学设计，尚未执行 CTD 化学物筛选、GeneCards 全量导出或 Top 20 结果计算。因此本文件不把候选化学物、P 值或 FDR 写成已观察到的结果。

## 研究问题

在不预先指定镉、砷、PFAS 等热门化学物的前提下，哪些环境相关化学物的 CTD 直接化学物–基因靶标集合，显著富集于 GeneCards 的 CRC 相关基因集合？

本阶段的 estimand 是：

> **环境化学物的 CTD 人类直接 chemical–gene target set 是否相对于 CTD 可观测人类基因背景，富集于 GeneCards CRC-associated gene set。**

这是一项候选发现和机制假设生成分析，不等价于证明某种暴露会导致 CRC。

## 1. GeneCards CRC 基因集

### 查询

- 查询词：`"colorectal cancer"`
- 主分析：按 Relevance Score 排名的 top 1000
- 敏感性分析：top 500、top 2000
- 记录：GeneCards snapshot date、查询词、rank、symbol、gene type、Relevance Score、Knowledge score（若导出可用）

GeneCards 当前 exact-phrase 结果约为 17,924 个基因；Relevance Score 表示查询匹配强度，不被解释为因果效应量或疾病风险效应量。

### ID 规范化

优先使用 NCBI Gene ID 或 HGNC ID 合并 CTD 与 GeneCards，gene symbol 仅用于展示。

- 主分析：保留可映射到 CTD 的人类基因
- 敏感性分析：仅保留 protein-coding genes
- 报告未映射基因数量

定义：

\[
G_{CRC,K}=\text{GeneCards top-}K\text{ genes after ID harmonization}
\]

## 2. CTD 化学物–基因数据

### 主数据表

仅使用 CTD 的直接、人工整理的：

> Chemical–Gene–Interaction–Organism–Evidence

不纳入：

- chemical–disease associations
- gene–disease associations
- CTD inference networks
- tetramers 或其他计算推断关系
- 仅由 CRC 疾病查询得到的关系

### 物种规则

主分析限定：`Homo sapiens`。

CTD 的跨物种关系不与人类关系混合。将其他物种通过 ortholog 映射到人类基因作为独立敏感性分析，并单独标记映射来源。

对每个化学物按 unique chemical ID × human gene ID 去重，而不是按 interaction 行数计数。

对每个 chemical 保留：

- CTD Chemical ID
- chemical name、CAS 或其他标识
- gene ID、gene symbol
- interaction type
- action
- organism
- evidence statement count
- distinct PMID count

## 3. 环境化学物纳入规则

化学物分类在查看富集结果前预先完成。

### Core environmental toxicants

- heavy metals and metalloids
- pesticides
- PFAS
- phthalates
- bisphenols
- PAHs
- VOCs
- flame retardants
- persistent organic pollutants
- industrial chemicals
- environmental endocrine disruptors

### Gray-zone chemicals

- pharmaceuticals detected in the environment
- endogenous metabolites
- dietary compounds
- nutrients
- common laboratory reagents

主分析只纳入 core environmental toxicants；gray-zone chemicals 作为敏感性分析。治疗药物、内源性分子和普通营养素不因名称判断，而依据预先建立的 CTD Chemical ID 分类表执行。

化学物家族、母体化合物、盐型和具体 congeners 需要显式记录，避免同一暴露实体被重复计数。

## 4. 背景集与富集检验

不能使用所有人类基因作为默认背景。背景定义为：

\[
U=\text{all unique human genes observed in eligible CTD direct chemical–gene data}
\]

令：

\[
G=G_{CRC,K}\cap U
\]

对于化学物 \(i\)：

\[
T_i=\text{unique human CTD target genes for chemical }i
\]

四格表为：

|  | CRC gene | Non-CRC gene |
|---|---:|---:|
| Chemical target | \(a=|T_i\cap G|\) | \(b=|T_i|-a\) |
| Non-target | \(c=|G|-a\) | \(d=|U|-a-b-c\) |

报告：

\[
OR_i=\frac{ad}{bc}
\]

\[
ER_i=\frac{a/|T_i|}{|G|/|U|}
\]

并执行 one-sided Fisher exact / hypergeometric enrichment test：

\[
P_i=P(X\geq a\mid N=|U|,K=|G|,n=|T_i|)
\]

对全部符合纳入规则的化学物统一进行 BH-FDR 校正。若四格表存在零格，保留 Fisher exact 结果，并使用连续性修正计算 OR 置信区间。

建议将 \(|T_i|\geq20\) 且 \(|T_i\cap G|\geq5\) 作为主排序的稳定性门槛；所有未过门槛的化学物仍保留在完整结果表中。

## 5. GeneCards 加权 overlap

Raw Relevance Score 的加和作为描述性结果，不假设该分数是线性测量尺度。

主加权指标采用 rank-based weight：

\[
w_g=\frac{1}{\log_2(rank_g+1)}
\]

\[
W_i=\sum_{g\in T_i\cap G}w_g
\]

同时报告：

- unweighted overlap
- rank-weighted overlap sum
- mean rank weight among overlapping genes
- overlap 中的 top-ranked genes

加权 overlap 用于生物学优先级排序，不替代 Fisher/FDR 的统计证据。

## 6. 研究热度与 CTD degree 偏倚

GeneCards 和 CTD 都受到文献可见度影响。TP53、APC、KRAS 等热门基因更容易在两个数据库中同时出现；研究较多的化学物也天然拥有更多 CTD 靶基因。

因此，除 Fisher 主分析外，执行 degree-matched permutation sensitivity analysis：

1. 在 CTD 人类基因宇宙中随机抽取与 \(|G|\) 相同的基因集。
2. 按基因的 CTD chemical degree 匹配随机集与 GeneCards CRC 集。
3. 重新计算每个化学物的 overlap 和富集统计量。
4. 报告经验 P 值及其与 Fisher/FDR 结论的一致性。

如果 GeneCards Knowledge score 可用，再进行 Knowledge score 匹配的敏感性分析。

## 7. 核心输出表

主表字段：

| Field | Definition |
|---|---|
| chemical | CTD chemical name |
| ctd_chemical_id | Stable CTD chemical identifier |
| chemical_class | Predefined environmental class |
| n_ctd_human_genes | \(|T_i|\) |
| crc_overlap | \(|T_i\cap G|\) |
| enrichment_ratio | \(ER_i\) |
| odds_ratio | \(OR_i\) |
| fisher_p | One-sided enrichment P |
| bh_fdr | BH-adjusted P across all tested chemicals |
| rank_weighted_overlap | \(W_i\) |
| n_evidence_statements | CTD evidence count |
| n_pmids | Distinct PMID count |
| top_overlap_genes | Highest-ranked overlapping genes |

结果排序采用双层结构：

1. 统计层：FDR、OR、overlap
2. 生物学优先级层：rank-weighted overlap、核心基因构成、CTD evidence/PMID 数量

## 8. Phase 1 不做的事情

第一轮不进行独立文献检索、不根据热门程度手工提升镉、砷、PFAS 等化学物，也不对 Top 20 做机制叙述。

但保留 CTD PMID 和 evidence 字段，用于可追溯性和后续 Phase 2 的文献撞车审查。

Phase 2 再回答：

> Top 20 中哪些已经被 CRC 研究充分覆盖，哪些具有较强富集但 CRC 专项研究较少？

## 9. 解释边界

本分析可以支持：

> 某环境化学物的 CTD 直接靶基因集合与 CRC-associated gene set 存在统计富集。

本分析不能单独支持：

- 该化学物会导致 CRC
- 该化学物在真实人群中的暴露足以产生该效应
- 富集基因是该化学物的特异性作用靶点
- 化学物对 CRC 具有治疗或预防效果

## 数据与引用

- GeneCards CRC exact-phrase search: <https://www.genecards.org/search/results?q=colorectal+cancer>
- GeneCards term-search definition: <https://docs.genecards.org/genecards/guide/search/term-search>
- CTD data resource and update description: <https://academic.oup.com/nar/article/51/D1/D1257/6725767>
- CTD downloads: <http://ctdbase.org/downloads/>

正式运行时必须额外记录 GeneCards 版本/导出日期、CTD 数据版本/下载日期、化学物分类表版本、ID 映射版本、随机种子和软件环境。
