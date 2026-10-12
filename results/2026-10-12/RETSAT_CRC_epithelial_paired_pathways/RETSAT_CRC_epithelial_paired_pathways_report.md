# RETSAT：CRC恶性上皮定位、患者级配对与独立通路线索

## 结论摘要

- **细胞定位：** GSE132465中RETSAT主要落在肿瘤上皮细胞：26.4%肿瘤上皮细胞检出RETSAT，非上皮细胞为3.0%（检出率比8.8倍）；83.7%的肿瘤样本RETSAT阳性细胞为上皮细胞。GSE200997的标记基因上皮样细胞检出率7.3%，非上皮细胞2.4%（3.1倍）。
- **配对表达：** 四个队列中RETSAT的肿瘤−正常中位差均为负。GSE132465 10对时未达显著（Δ=-0.67, P=0.16）；GSE200997 7对（Δ=-1.32, P=0.0156），GSE44076 98对和GSE32323 17对均显著下降。
- **视黄醇代谢：** 排除RETSAT自身后，两个bulk配对队列的视黄醇代谢分数均下降（GSE44076 Δ=-0.45, P=3.4e-16; GSE32323 Δ=-0.43, P=0.00038），且患者配对的RETSAT差值与通路差值正相关（ρ分别为0.56、0.55）。但是该通路总分在两套单细胞肿瘤上皮伪bulk中都未显著，故当前证据支持bulk层面的独立视黄醇代谢线索，尚未证明它是肿瘤上皮特异改变。
- **铁死亡：** KEGG铁死亡转录分数在GSE132465肿瘤上皮伪bulk与GSE44076 bulk中上升；GSE200997及GSE32323不显著。通路成员方向有升有降，整体缺少跨队列一致性，不能据此推断铁死亡增强或减弱。

## 研究设计与统计

只读取本项目work目录下的原始GEO表达矩阵、原始细胞注释及平台探针注释；未读取历史分析结果。GSE132465采用作者给出的Cell_type标签；将肿瘤样本的上皮细胞称为肿瘤上皮候选群，未在本分析中进行inferCNV/拷贝数恶性判定。该队列含23位CRC患者肿瘤样本及10位配对正常黏膜样本。[GSE132465 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)。

GSE200997的本地细胞注释没有通用细胞类型字段，因此以EPCAM、KRT8、KRT18、KRT19、KRT20中至少2个基因检出定义上皮样细胞；配对肿瘤样本细胞数均≥49，正常样本均≥161。该操作性注释用于敏感性定位，不等同于CNV恶性鉴定。[GSE200997 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997)。

单细胞按患者/样本与细胞类型汇总原始UMI为伪bulk，RETSAT按总UMI计算log2(CPM+0.5)。配对患者为统计单位，不把单个细胞当作独立重复。所有单细胞样本上皮细胞数要求≥10；另报告≥20细胞敏感性分析。GSE44076与GSE32323使用GEO提供的log2 RMA处理表达值；GSE44076保留98对CRC与相邻正常，排除50名健康供者；GSE32323保留17对患者样本，排除细胞系。[GSE44076 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44076)，[GSE32323 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32323)。

视黄醇代谢和铁死亡基因集分别来自KEGG人类通路[hsa00830](https://www.kegg.jp/pathway/hsa00830)与[hsa04216](https://www.kegg.jp/pathway/hsa04216)。视黄醇集合原有69个成员，其中RETSAT被剔除，余68个；通路分数是各队列每个基因跨配对样本标准化后Z分数的均值。所有配对差值为肿瘤−正常。配对检验使用双侧Wilcoxon signed-rank；三项主要指标（RETSAT、视黄醇分数、铁死亡分数）在队列内BH校正；通路成员在各队列/通路内BH校正。RETSAT与通路的关联用患者配对差值做Spearman相关。

## 细胞定位

| 队列 | 上皮定义 | 肿瘤上皮细胞 | RETSAT阳性上皮细胞 | 上皮检出率 | 非上皮检出率 | 上皮/非上皮率比 |
|---|---|---:|---:|---:|---:|---:|
| GSE132465 | 作者注释Epithelial cells | 17469 | 4613 | 26.41% | 3.01% | 8.79 |
| GSE200997 | ≥2/5上皮标志物 | 6814 | 499 | 7.32% | 2.39% | 3.06 |

GSE132465各细胞类型的RETSAT检出比例见图；上皮细胞不是唯一RETSAT阳性类型，但在肿瘤样本中显著富集。GSE200997提示方向一致，但其标记分类是本次重建，证据等级较低。

## 患者配对结果

| 队列 | 指标 | 配对数 | 中位Δ | Δ为正比例 | Wilcoxon P | 队列内BH q |
|---|---|---:|---:|---:|---:|---:|
| GSE132465 | RETSAT_log2expr | 10 | -0.6706 | 20.0% | 0.16 | 0.24 |
| GSE132465 | Retinol_metabolism_score | 10 | -0.002812 | 50.0% | 0.492 | 0.492 |
| GSE132465 | Ferroptosis_score | 10 | 0.5796 | 100.0% | 0.00195 | 0.00586 |
| GSE200997 | RETSAT_log2expr | 7 | -1.319 | 0.0% | 0.0156 | 0.0469 |
| GSE200997 | Retinol_metabolism_score | 7 | -0.1948 | 14.3% | 0.156 | 0.234 |
| GSE200997 | Ferroptosis_score | 7 | 0.05543 | 57.1% | 0.812 | 0.812 |
| GSE44076 | RETSAT_log2expr | 98 | -2.076 | 2.0% | 1.28e-17 | 3.84e-17 |
| GSE44076 | Retinol_metabolism_score | 98 | -0.4469 | 7.1% | 3.38e-16 | 5.06e-16 |
| GSE44076 | Ferroptosis_score | 98 | 0.0501 | 57.1% | 0.00568 | 0.00568 |
| GSE32323 | RETSAT_log2expr | 17 | -1.183 | 5.9% | 0.000107 | 0.00032 |
| GSE32323 | Retinol_metabolism_score | 17 | -0.4312 | 5.9% | 0.000381 | 0.000572 |
| GSE32323 | Ferroptosis_score | 17 | 0.1824 | 64.7% | 0.225 | 0.225 |

RETSAT在四个队列的配对中位变化方向均为下降。GSE132465主单细胞配对分析10对，P=0.160；移除每侧上皮细胞<20的配对后为8对，RETSAT中位差仍为负且P=0.0234（见敏感性表）。

### 单细胞每对上皮细胞数敏感性

| 队列 | 最低细胞数/样本 | 配对数 | RETSAT Δ | P | 视黄醇分数 Δ | P | 铁死亡分数 Δ | P |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GSE132465 | 10 | 10 | -0.671 | 0.16 | -0.003 | 0.492 | 0.580 | 0.00195 |
| GSE132465 | 20 | 8 | -0.853 | 0.0234 | -0.021 | 0.641 | 0.608 | 0.00781 |
| GSE200997 | 10 | 7 | -1.319 | 0.0156 | -0.195 | 0.156 | 0.055 | 0.812 |
| GSE200997 | 20 | 7 | -1.319 | 0.0156 | -0.195 | 0.156 | 0.055 | 0.812 |

### RETSAT变化与独立通路分数变化的患者级相关

| 队列 | 通路 | 配对数 | Spearman ρ | P | 队列内BH q |
|---|---|---:|---:|---:|---:|
| GSE132465 | 视黄醇代谢 | 10 | -0.079 | 0.829 | 0.829 |
| GSE132465 | 铁死亡 | 10 | -0.345 | 0.328 | 0.656 |
| GSE200997 | 视黄醇代谢 | 7 | 0.286 | 0.535 | 0.702 |
| GSE200997 | 铁死亡 | 7 | 0.179 | 0.702 | 0.702 |
| GSE44076 | 视黄醇代谢 | 98 | 0.557 | 2.61e-09 | 5.21e-09 |
| GSE44076 | 铁死亡 | 98 | 0.286 | 0.00434 | 0.00434 |
| GSE32323 | 视黄醇代谢 | 17 | 0.554 | 0.021 | 0.0421 |
| GSE32323 | 铁死亡 | 17 | 0.179 | 0.492 | 0.492 |

## 独立视黄醇代谢与铁死亡线索

两个bulk队列共同有17个视黄醇代谢通路成员在肿瘤中下降且各自队列内FDR<0.05：ADH1A, ADH1B, ADH1C, ADH6, ALDH1A1, AOX1, CYP26B1, CYP2C18, CYP2C8, CYP3A5, DGAT1, DHRS3, DHRS9, RDH13, RDH5, UGT2A3, UGT2B28。通路评分中没有RETSAT本身，因此这一结果是RETSAT之外成员提供的转录层面支持；不过两套单细胞上皮伪bulk均未显示视黄醇总分显著下降，提示需要将其解释为bulk组织层面的线索，不能直接归因于恶性上皮的内在代谢改变。

铁死亡通路评分在GSE132465上皮伪bulk与GSE44076中显著上升，但GSE200997和GSE32323未显著。成员既有上升也有下降，例如铁摄取/转运、脂质重塑、抗氧化与铁储存相关成员方向不一。KEGG转录成员集不直接测量脂质过氧化或细胞死亡活性，因此目前只能说有部分、非稳定的转录关联线索，不能声称铁死亡增强或受RETSAT调控。

## 文件

- RETSAT_GSE132465_celltype_localization.png：GSE132465肿瘤细胞类型定位及RETSAT阳性率。
- RETSAT_paired_cohort_signals.png：四队列患者配对表达与通路分数。
- cohort_paired_summary.csv：四队列配对检验及BH校正。
- paired_delta_correlations.csv：患者级RETSAT差值与通路分数差值相关。
- pathway_member_paired_tests.csv：单个通路成员配对效应及FDR。
- pathway_gene_coverage.csv：每个队列进入通路评分的有效基因数。
- singlecell_epithelial_paired_cell_counts.csv、singlecell_cell_count_sensitivity.csv：单细胞配对细胞数及敏感性分析。
- singlecell_pseudobulk_all_groups_*.csv、*_paired_expression_and_scores.csv：样本级伪bulk表达数据。
- 分析脚本：analysis/retsat_crc_epithelial_paired_pathways.py；本次输出整理脚本：analysis/finalize_retsat_crc_outputs.py。

## 复现

在项目根目录依次运行：python analysis/retsat_crc_epithelial_paired_pathways.py，然后运行 python analysis/finalize_retsat_crc_outputs.py。前者从work中的本地原始矩阵重算；后者重建细胞数敏感性、有效通路覆盖、图和报告。

## 边界

GSE44076/GSE32323为bulk组织，细胞组成和分化状态差异可能影响代谢通路分数；单细胞患者数较少，GSE200997上皮亚群标签依赖本次的标志基因定义。RETSAT与通路变化相关不代表因果。GSE21510未纳入主bulk分析，因为原始元数据显示其肿瘤样本为LCM、配对正常为homogenized，取材处理区室不匹配。
