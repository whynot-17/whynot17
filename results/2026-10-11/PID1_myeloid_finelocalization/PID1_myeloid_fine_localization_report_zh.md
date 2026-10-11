# PID1 在 CRC 肿瘤髓系细胞中的精细定位与功能状态分析

日期：2026-10-11  
结果目录：`results/2026-10-11/PID1_myeloid_finelocalization/`

## 结论摘要

1. **PID1在单核细胞和巨噬细胞中都可检出，但没有在两队列中复现一个独立的 FCGR3A⁺ 单核细胞簇。** FCGR3A/MS4A7 信号主要和 C1QC⁺ 巨噬细胞群共现；因此不能把它解释为已经解析出的 FCGR3A⁺ 单核细胞特异表达。
2. **最稳定的亚群差别是 C1QC⁺ 巨噬细胞中 PID1 较低、炎症状态髓系细胞中较高。** GSE132465患者级配对比较支持 C1QC⁺ 与炎症巨噬的差异；GSE200997也呈相同的细胞检出率方向，但样本数和炎症巨噬细胞数偏少，患者级比较未达FDR显著。
3. **单核细胞整体高于巨噬细胞的证据有限。** GSE132465在每个患者两类细胞均≥10时方向和校正后检验支持单核细胞较高；提高到≥20后未通过BH-FDR，≥50时不显著。GSE200997在≥10时同向但未通过FDR，较高门槛可配对的样本很少。
4. **募集、炎症、干扰素、糖酵解、脂质/胆固醇和抗原呈递程序没有显示可跨队列重复的PID1关联。** GSE132465中PID1与糖酵解分数的调整后关联通过BH-FDR，但GSE200997未复现；患者级相关检验也未通过FDR。
5. **DPT拟时序只作探索。** 两队列都可计算连续的图排序，但PID1与排序的负向关联只在GSE132465显著，GSE200997未显著；这不是谱系追踪，也不能证明分化方向或因果机制。

本报告只使用两个GEO队列提供的 processed raw UMI 矩阵与原始样本/细胞注释重新分析；没有读取旧项目的派生结果或缓存。没有运行机器学习。

## 数据与质量控制

GSE132465使用原始注释中的肿瘤髓系细胞作为候选。GSE200997原始注释没有髓系细胞类型标签，因此先在全部肿瘤细胞的原始UMI上用核心髓系和竞争谱系marker筛选；PTPRC单独不作为髓系证据。各队列分别处理，没有合并不同队列的注释。

过滤条件：UMI≥500、检测基因≥200、线粒体比例≤20%，并按样本对UMI与基因数使用4 MAD离群值审查。每样本在至少80个QC通过的髓系候选细胞时运行Scrublet。跨谱系marker共检出保存在逐细胞审计表，不单独当作双细胞删除；清楚的T/NK簇单列为非髓系，不进入单核/巨噬对照。

| 队列 | 注释细胞 | 肿瘤细胞 | 髓系候选 | QC通过 | Scrublet调用 | 保留髓系/相关细胞 | Leiden簇数 | 患者/样本单位 | 运行Scrublet单位 | 跳过Scrublet单位（候选<80） |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE132465 | 63689 | 47285 | 6400 | 6394 | 56 | 6338 | 13 | 23 | 16 | 7 |
| GSE200997 | 49859 | 31586 | 804 | 789 | 72 | 717 | 9 | 15 | 4 | 11 |

GSE132465有6个候选细胞未通过QC，56个QC通过细胞被Scrublet标记，最终保留6338/6400。GSE200997有15个未通过QC，72个Scrublet调用，最终保留717/804。GSE200997只有4/15个样本达到Scrublet运行门槛；这点限制了双细胞判断覆盖率。每个细胞的QC指标、Scrublet分数与运行状态、宽/窄marker证据均见对应 `*_cell_QC_doublet_marker_audit.tsv.gz`。宽marker冲突标记不是删除数。

两队列独立使用Scanpy 1.12.1重新归一化、log1p、选择2000个Seurat风格高变基因（按样本批次可用时进行样本内选择）、PCA、邻居图、Leiden 0.6聚类和UMAP。细胞状态由不含PID1的marker程序分数结合每簇差异marker核对。聚类marker p值仅用于描述性注释，不作为患者级差异检验。

## 亚群注释、细胞数和患者/样本覆盖

以下 `n_units` 在GSE132465代表患者，在GSE200997代表肿瘤样本。少量细胞或较少单位覆盖的群体应视为探索性。T/NK簇与未解析髓系细胞仍保存在注释表中以供审查；`cellchat_eligible=否`表示先不要作为细胞通讯的明确髓系亚群输入。

| 队列 | 亚群 | 细胞数 | 患者/样本数 |
| --- | --- | --- | --- |
| GSE132465 | DC_like | 348 | 22 |
| GSE132465 | DC_pDC_like | 275 | 22 |
| GSE132465 | Macrophage_C1QC | 2061 | 23 |
| GSE132465 | Macrophage_IFN | 163 | 10 |
| GSE132465 | Macrophage_SPP1_TREM2 | 243 | 16 |
| GSE132465 | Macrophage_inflammatory | 1800 | 23 |
| GSE132465 | Monocyte_CD14_classical | 519 | 17 |
| GSE132465 | Monocyte_CD14_classical_inflammatory | 755 | 23 |
| GSE132465 | Non_myeloid_T_NK_contamination | 174 | 22 |
| GSE200997 | DC_like | 127 | 13 |
| GSE200997 | DC_mature_migratory | 10 | 6 |
| GSE200997 | DC_pDC_like | 14 | 7 |
| GSE200997 | Macrophage_C1QC | 253 | 14 |
| GSE200997 | Macrophage_inflammatory | 28 | 6 |
| GSE200997 | Monocyte_CD14_classical | 122 | 12 |
| GSE200997 | Monocyte_CD14_classical_inflammatory | 126 | 12 |
| GSE200997 | Myeloid_unresolved | 37 | 5 |

关键marker核验：

- CD14经典/炎症单核细胞由FCN1、VCAN、S100A8/A9/A12及CD14支撑。
- C1QC巨噬细胞以C1QA/B/C、APOC1、GPNMB、CD163/VSIG4等为特征；GSE132465 C1QC簇中FCGR3A检出率约81%，GSE200997 C1QC簇约68%。这支持FCGR3A信号主要位于巨噬细胞样群，不能据此重命名为FCGR3A⁺单核细胞。
- SPP1/TREM2、炎症和IFN巨噬亚群在GSE132465均可辨认；GSE200997没有可稳定复现的SPP1/TREM2簇，炎症巨噬群也只有28个细胞、6个样本。
- pDC由LILRA4/GZMB/IL3RA等支撑，cDC样群有CD1C/CLEC10A，迁移成熟DC有LAMP3/CCR7/FSCN1。GSE132465一个簇以CD3D/E、TRAC/TRBC为主，明确标为T/NK污染。

逐簇分数、患者数和前12个差异marker见 [`PID1_cluster_annotation_with_top_markers_all_cohorts.tsv`](PID1_cluster_annotation_with_top_markers_all_cohorts.tsv) 与队列级 `*_cluster_rank_markers.tsv.gz`。

## PID1检出与亚群定位

PID1检出率按该亚群所有保留细胞中原始UMI是否大于0计算；表达均值为每细胞log1p归一化值。图中的患者/样本点仅显示该单位该亚群至少有10个细胞，完整分母与逐单位结果保存在表格中。

| 队列 | 亚群 | 细胞数 | 患者/样本数 | PID1阳性细胞数 | PID1细胞检出率（%） | 至少检出1个PID1阳性细胞的单位数 | 可用于CellChat |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE132465 | DC_like | 348 | 22 | 96 | 27.586 | 18 | 是 |
| GSE132465 | DC_pDC_like | 275 | 22 | 1 | 0.364 | 1 | 是 |
| GSE132465 | Macrophage_C1QC | 2061 | 23 | 63 | 3.057 | 19 | 是 |
| GSE132465 | Macrophage_IFN | 163 | 10 | 10 | 6.135 | 3 | 是 |
| GSE132465 | Macrophage_SPP1_TREM2 | 243 | 16 | 18 | 7.407 | 3 | 是 |
| GSE132465 | Macrophage_inflammatory | 1800 | 23 | 346 | 19.222 | 22 | 是 |
| GSE132465 | Monocyte_CD14_classical | 519 | 17 | 69 | 13.295 | 12 | 是 |
| GSE132465 | Monocyte_CD14_classical_inflammatory | 755 | 23 | 102 | 13.510 | 17 | 是 |
| GSE132465 | Non_myeloid_T_NK_contamination | 174 | 22 | 21 | 12.069 | 8 | 否 |
| GSE200997 | DC_like | 127 | 13 | 17 | 13.386 | 6 | 是 |
| GSE200997 | DC_mature_migratory | 10 | 6 | 0 | 0.000 | 0 | 是 |
| GSE200997 | DC_pDC_like | 14 | 7 | 0 | 0.000 | 0 | 是 |
| GSE200997 | Macrophage_C1QC | 253 | 14 | 19 | 7.510 | 5 | 是 |
| GSE200997 | Macrophage_inflammatory | 28 | 6 | 8 | 28.571 | 3 | 是 |
| GSE200997 | Monocyte_CD14_classical | 122 | 12 | 22 | 18.033 | 4 | 是 |
| GSE200997 | Monocyte_CD14_classical_inflammatory | 126 | 12 | 27 | 21.429 | 6 | 是 |
| GSE200997 | Myeloid_unresolved | 37 | 5 | 5 | 13.514 | 2 | 否 |

**患者级单核细胞与巨噬细胞比较：**每个单位先按该亚群细胞数加权汇总log1p PID1，再在同时有单核和巨噬细胞的患者/样本内配对。表中差值为“单核均值−巨噬均值”；p值用Wilcoxon配对检验，并在相同细胞数门槛内对各亚群比较进行BH校正。

| 队列 | 每单位每类最低细胞数 | 配对单位数 | 配对差值中位数（单核−巨噬） | Wilcoxon p | BH-FDR |
| --- | --- | --- | --- | --- | --- |
| GSE132465 | 10 | 18 | 0.034 | 0.012 | 0.036 |
| GSE132465 | 20 | 14 | 0.034 | 0.025 | 0.074 |
| GSE132465 | 50 | 9 | 0.031 | 0.203 | 0.500 |
| GSE200997 | 10 | 5 | 0.191 | 0.062 | 0.188 |
| GSE200997 | 20 | 2 | 0.237 | 0.500 | 0.500 |
| GSE200997 | 50 | 1 | 0.283 | — | — |

GSE132465的≥10细胞门槛结果为18个配对患者，14/18名患者差值为正，BH-FDR=0.036；≥20时14个患者，FDR=0.074；≥50时9个患者，FDR=0.50。GSE200997在≥10时只有5个配对样本，5/5同向，但FDR=0.188；≥20时2个，≥50时1个，不能据此建立稳定差异。

GSE132465中更清楚的定位是**炎症巨噬细胞高于C1QC巨噬细胞**：

| 队列 | 每亚群每单位最低细胞数 | 配对单位数 | 配对差值中位数（C1QC−炎症巨噬） | Wilcoxon p | BH-FDR |
| --- | --- | --- | --- | --- | --- |
| GSE132465 | 10 | 20 | -0.182 | 3.81e-06 | 5.72e-05 |
| GSE132465 | 20 | 17 | -0.177 | 3.05e-05 | 4.58e-04 |
| GSE132465 | 50 | 10 | -0.191 | 0.002 | 0.020 |
| GSE200997 | 10 | 1 | -0.387 | — | — |
| GSE200997 | 20 | 0 | — | — | — |
| GSE200997 | 50 | 0 | — | — | — |

在GSE132465，这一差异在≥10、≥20、≥50细胞门槛下均通过BH-FDR；GSE200997在炎症巨噬亚群中仅28个细胞、6个样本，没有FDR显著的细分亚群配对比较。两队列的总体细胞检出率均显示炎症髓系亚群高于C1QC巨噬细胞，但跨队列患者级证据仍有限。

## 单核细胞内功能程序关联

六类程序分数均不含PID1：趋化募集、炎症、干扰素、糖酵解、脂质/胆固醇和抗原呈递。细胞级模型以log1p归一化PID1为因变量、程序分数为自变量，纳入样本固定效应、log10总UMI和线粒体比例，并按患者/样本聚类稳健标准误；每队列六个程序内BH校正。另以每单位至少10个单核细胞计算患者级Spearman相关。

| 队列 | PID1-free程序 | 患者/样本数 | 每SD程序分数的回归系数 | p | BH-FDR |
| --- | --- | --- | --- | --- | --- |
| GSE132465 | recruitment_chemotaxis | 23 | -0.065 | 0.040 | 0.063 |
| GSE132465 | inflammation | 23 | -0.052 | 0.026 | 0.063 |
| GSE132465 | interferon_response | 23 | -0.067 | 0.042 | 0.063 |
| GSE132465 | glycolysis | 23 | 0.168 | 1.02e-10 | 6.15e-10 |
| GSE132465 | lipid_cholesterol | 23 | -0.014 | 0.832 | 0.832 |
| GSE132465 | antigen_presentation | 23 | -0.063 | 0.172 | 0.206 |
| GSE200997 | recruitment_chemotaxis | 14 | -0.014 | 0.914 | 0.921 |
| GSE200997 | inflammation | 14 | 0.038 | 0.385 | 0.770 |
| GSE200997 | interferon_response | 14 | -0.117 | 0.275 | 0.770 |
| GSE200997 | glycolysis | 14 | 0.018 | 0.921 | 0.921 |
| GSE200997 | lipid_cholesterol | 14 | 0.135 | 0.269 | 0.770 |
| GSE200997 | antigen_presentation | 14 | 0.016 | 0.866 | 0.921 |

GSE132465只有糖酵解程序的细胞级调整关联通过BH-FDR（β=0.168，FDR=6.15×10⁻¹⁰）；同一程序患者级相关不显著（18个单位，ρ=0.135，FDR=0.848），且GSE200997细胞级关联未复现（FDR=0.921）。GSE200997患者级炎症和糖酵解相关各有名义p=0.042，但仅6个单位，BH-FDR均为0.125。两队列没有可重复的募集、炎症、IFN、脂质/胆固醇或抗原呈递关联。

## 单核到巨噬状态的DPT探索

DPT仅在单核和巨噬细胞均≥200个、至少8个单位、至少3个状态且邻接图最大连通分量≥90%时运行。根取CD14经典单核细胞中最接近该群PCA中位数的细胞。DPT是一种表达相似性排序，不是谱系追踪。

| 队列 | 拟时序状态 | 单核细胞数 | 巨噬细胞数 | 患者/样本数 | 运行条件/状态 |
| --- | --- | --- | --- | --- | --- |
| GSE132465 | exploratory_DPT_computed | 1274 | 4267 | 23 | n_monocytes=1274; n_macrophages=4267; units=23; graph_components=1; largest_component_fraction=1.000 |
| GSE200997 | exploratory_DPT_computed | 248 | 281 | 15 | n_monocytes=248; n_macrophages=281; units=15; graph_components=1; largest_component_fraction=1.000 |

| 队列 | 状态 | 细胞数 | 患者/样本数 | PID1与DPT的回归系数 | p | 解释限制 |
| --- | --- | --- | --- | --- | --- | --- |
| GSE132465 | exploratory_DPT_computed | 5541 | 23 | -0.313 | 1.84e-09 | DPT is an expression ordering, not lineage tracing |
| GSE200997 | exploratory_DPT_computed | 529 | 15 | -0.330 | 0.071 | DPT is an expression ordering, not lineage tracing |

GSE132465中按上述排序的PID1系数为−0.313（p=1.84×10⁻⁹）；GSE200997为−0.330（p=0.071）。只有第一队列有统计支持。不能把这个排序或系数解释成PID1随真实单核→巨噬分化必然下降。

## 与已发表Nature Cancer PID1机制的范围核对

2026年的Nature Cancer研究已报道PID1在肿瘤相关巨噬细胞中的作用，并提出PID1相关脂质代谢/氧固醇机制：PID1缺失改变LDLR介导的LDL摄取和胆固醇/ROS，生成5α,6α-EC与7β-OHC，抑制mTOR–STAT6信号并推动巨噬细胞转向抗肿瘤状态；论文也用人CRC数据集（包括GSE146771）支撑TAM中的PID1表达现象。[Nature Cancer原文](https://www.nature.com/articles/s43018-026-01189-0)；[PubMed记录](https://pubmed.ncbi.nlm.nih.gov/42303776/)。

本分析的增量是用GSE132465与GSE200997进行独立的单细胞亚群定位和患者/样本级复核，重点看到C1QC与炎症髓系状态的PID1差别。它没有测定LDL摄取、胆固醇、氧固醇、ROS、mTOR/STAT6或巨噬细胞功能，也没有重复论文的干预实验；因此不构成新的机制证据。单细胞表达与程序分数的相关性均按观察性关联解释，不作因果推断。

## 图件与可复用文件

- UMAP：[`GSE132465_PID1_myeloid_UMAP.png`](GSE132465_PID1_myeloid_UMAP.png)、[`GSE200997_PID1_myeloid_UMAP.png`](GSE200997_PID1_myeloid_UMAP.png)
- DotPlot：[`GSE132465_PID1_myeloid_DotPlot.png`](GSE132465_PID1_myeloid_DotPlot.png)、[`GSE200997_PID1_myeloid_DotPlot.png`](GSE200997_PID1_myeloid_DotPlot.png)
- ViolinPlot：[`GSE132465_PID1_myeloid_ViolinPlot.png`](GSE132465_PID1_myeloid_ViolinPlot.png)、[`GSE200997_PID1_myeloid_ViolinPlot.png`](GSE200997_PID1_myeloid_ViolinPlot.png)
- 患者/样本检出率：[`GSE132465_PID1_detection_rate.png`](GSE132465_PID1_detection_rate.png)、[`GSE200997_PID1_detection_rate.png`](GSE200997_PID1_detection_rate.png)
- DPT图：[`GSE132465_PID1_DPT_exploration.png`](GSE132465_PID1_DPT_exploration.png)、[`GSE200997_PID1_DPT_exploration.png`](GSE200997_PID1_DPT_exploration.png)
- 逐细胞CellChat注释：[`GSE132465_PID1_CellChat_myeloid_annotations.tsv.gz`](GSE132465_PID1_CellChat_myeloid_annotations.tsv.gz)、[`GSE200997_PID1_CellChat_myeloid_annotations.tsv.gz`](GSE200997_PID1_CellChat_myeloid_annotations.tsv.gz)。字段含细胞、患者/样本、状态、PID1值、QC、Scrublet、marker冲突和 `cellchat_eligible`。未解析髓系与明确非髓系污染标为不可直接用于CellChat。
- 患者级PID1值：[`PID1_patient_sample_subtype_summary_all_cohorts.tsv`](PID1_patient_sample_subtype_summary_all_cohorts.tsv)
- 10/20/50细胞阈值敏感性：[`PID1_unit_level_threshold_sensitivity_all_cohorts.tsv`](PID1_unit_level_threshold_sensitivity_all_cohorts.tsv)
- 程序模型与单位相关：[`PID1_monocyte_program_adjusted_models_all_cohorts.tsv`](PID1_monocyte_program_adjusted_models_all_cohorts.tsv)、[`PID1_program_coverage_and_unit_correlations_all_cohorts.tsv`](PID1_program_coverage_and_unit_correlations_all_cohorts.tsv)
- QC/双细胞审计：`*_cell_QC_doublet_marker_audit.tsv.gz`、`*_QC_flow_by_sample.tsv`、`*_broad_myeloid_marker_screen.tsv.gz`
- 重跑参数和来源：[`analysis_manifest.json`](analysis_manifest.json)

## 解释边界

- GSE200997只剩248个单核细胞、281个巨噬细胞；更高细胞数阈值下可配对样本不足，亚群差异的阴性结果不能视为等效性证明。
- GSE200997的Scrublet覆盖4/15个样本；不足80个候选髓系细胞的样本没有运行Scrublet。宽marker共表达留在审计表中，可能含环境RNA/吞噬RNA，也可能含残余异型双细胞。
- PID1在C1QC巨噬细胞中的检出率低，按患者分层后稀疏；表达阳性率变化可能受细胞数量与样本组成影响。
- marker程序和聚类注释是本次数据中的可复核标签；队列差异较大时没有强行合并。CellChat前仍应按患者与每样本最小细胞数再次过滤。
