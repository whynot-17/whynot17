# 环境毒理学与 CRC：Phase 1 分析协议

## 当前状态

本文件记录已经审定的 **CTD × GeneCards 环境化学物富集筛选方案**。

当前只完成方法学设计，尚未执行 CTD 化学物筛选、GeneCards 全量导出或 Top 20 结果计算。因此本文件不把候选化学物、P 值或 FDR 写成已观察到的结果。

## 研究问题

在不预先指定镉、砷、PFAS 等热门化学物的前提下，哪些环境相关化学物的 CTD 直接化学物–基因相互作用集合，显著富集于 GeneCards 的 CRC 相关基因集合？

本阶段的 estimand 是：

> **环境化学物的 CTD 人类直接 chemical-interacting gene set 是否相对于同一候选化学物集合形成的 CTD 人类基因背景，富集于 GeneCards CRC-associated gene set。**

这是一项候选发现和机制假设生成分析，不等价于证明某种暴露会导致 CRC。

## 1. GeneCards CRC 基因集

### 查询

- 主分析：`[disorders] "colorectal cancer"`，按 Relevance Score 排名的 top 1000
- 敏感性分析：全域 `"colorectal cancer"` 的 top 500、top 1000、top 2000
- 记录：GeneCards snapshot date、查询词、rank、symbol、gene type、Relevance Score、Knowledge score（若导出可用）

GeneCards 支持用 `[scope] term` 将查询限制到 Disorders section，并通过当前过滤结果页 export。主分析优先使用该 scoped query，降低 GeneCards 与 CTD chemical–gene 文献整合之间的循环污染；全域查询保留为敏感性分析。Relevance Score 表示查询匹配强度，不被解释为因果效应量或疾病风险效应量。

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

对每个化学物按 unique chemical ID × human gene ID 去重，而不是按 interaction 行数计数。CTD interaction 可以是 increases expression、decreases activity 等关系，因此下文统一称为 **chemical-interacting gene set**，不把所有关系都称为直接结合靶点。

对每个 chemical 保留：

- CTD Chemical ID
- chemical name、CAS 或其他标识
- gene ID、gene symbol
- interaction type
- action
- organism
- `n_raw_interaction_rows`
- `n_unique_pmids`
- `n_unique_chemical_gene_pairs`

## 3. 环境化学物纳入规则

化学物分类在查看富集结果前预先完成。分类优先使用 CTD chemical vocabulary 的 `ParentIDs`、`TreeNumbers` 和 `ParentTreeNumbers` 的 MeSH hierarchy；不按化学物名称人工挑选。当前 CTD bulk chemical file 的字段不包含 DrugBankID，因此 DrugBank 只能作为可用时的辅助排除字段，不能假定本次下载中存在该列。当前实现同时使用预先冻结的 CTD Definition/MESHSynonyms/CTDCuratedSynonyms 药物语义正则和 DrugCentral 的 CAS/INN 参照作为自动化 fallback；包含药物后代的过宽 MeSH 树根不单独作为环境类别。

### Core environmental toxicants

- heavy metals and metalloids
- pesticides
- PFAS
- phthalates
- bisphenols
- PAHs（MeSH 层级后再要求 PubChem 分子式为纯 C/H，作为结构级 guard）
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

主分析只纳入 core environmental toxicants；gray-zone chemicals 作为敏感性分析。主分析分类表冻结为 CTD Chemical ID、MeSH tree path 和分类规则，药物/内源性/营养素排除也只依据冻结的结构化字段执行，不因 Top 20 结果调整。

当前实现使用的 core MeSH path 规则包括：

- heavy metals：`D01.268.556`、`D01.552.544`
- pesticides：`D27.720.031.700`、`D27.888.723`
- PFAS/perfluoro chemicals：`D02.455.526.510.435`
- phthalates：`D02.241.223.805`
- bisphenols：`D02.455.426.559.389.657.110`
- PAHs：`D02.455.426.559.847`、`D04.615`
- VOCs：`D02.974`
- flame retardants：`D27.720.361`
- persistent organic pollutants：`D27.888.284.429`
- endocrine disruptors：`D27.505.696.353`、`D27.888.141`
- environmental pollutants：`D27.888.284`

化学物家族、母体化合物、盐型和具体 congeners 需要显式记录，避免同一暴露实体被重复计数。

## 4. 背景集与富集检验

不能使用所有人类基因作为默认背景。背景定义为：

\[
U_{core}=\bigcup_{i\in core} I_i
\]

其中 `core` 是最终纳入的 core environmental toxicant 集合，且：

\[
U_{allCTD}=\text{all unique human genes observed in all eligible CTD direct chemical–gene data}
\]

主分析使用 \(U_{core}\)，\(U_{allCTD}\) 作为背景集敏感性分析。这样被排除的药物、内源性分子或营养素不会进入主分析的 background frequency。

令：

\[
G=G_{CRC,K}\cap U_{core}
\]

对于化学物 \(i\)：

\[
I_i=\text{unique human CTD interacting genes for chemical }i
\]

四格表为：

|  | CRC gene | Non-CRC gene |
|---|---:|---:|
|  | CRC gene | Non-CRC gene |
|---|---:|---:|
| Present in \(I_i\) | \(a=|I_i\cap G|\) | \(b=|I_i|-a\) |
| Absent from \(I_i\) | \(c=|G|-a\) | \(d=|U_{core}|-a-b-c\) |

报告：

\[
OR_i=\frac{ad}{bc}
\]

\[
ER_i=\frac{a/|I_i|}{|G|/|U_{core}|}
\]

并执行 one-sided Fisher exact / hypergeometric enrichment test：

\[
P_i=P(X\geq a\mid N=|U_{core}|,K=|G|,n=|I_i|)
\]

对全部符合纳入规则的化学物统一进行 BH-FDR 校正。若四格表存在零格，保留 Fisher exact 结果，并使用连续性修正计算 OR 置信区间。

建议将 \(|I_i|\geq20\) 且 \(|I_i\cap G|\geq5\) 作为主排序的稳定性门槛；所有未过门槛的化学物仍保留在完整结果表中。

## 5. GeneCards 加权 overlap

Raw Relevance Score 的加和作为描述性结果，不假设该分数是线性测量尺度。

主加权指标采用 rank-based weight：

\[
w_g=\frac{1}{\log_2(rank_g+1)}
\]

\[
W_i=\sum_{g\in I_i\cap G}w_g
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

1. 在 \(U_{core}\) CTD 人类基因宇宙中随机抽取与 \(|G|\) 相同的基因集。
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
| n_ctd_human_genes | \(|I_i|\) |
| crc_overlap | \(|I_i\cap G|\) |
| enrichment_ratio | \(ER_i\) |
| odds_ratio | \(OR_i\) |
| fisher_p | One-sided enrichment P |
| bh_fdr | BH-adjusted P across all tested chemicals |
| rank_weighted_overlap | \(W_i\) |
| n_raw_interaction_rows | Raw CTD interaction rows |
| n_unique_pmids | Distinct PMID count |
| n_unique_chemical_gene_pairs | Unique chemical × gene pairs |
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
