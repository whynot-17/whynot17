# PID1 与上皮–髓系亚群通讯：逐患者/样本 CellChat 分析

日期：2026-10-11  
结果目录：`results/2026-10-11/PID1_CellChat_subtype_consistency/`

## 结论摘要

- 上皮→C1QC巨噬细胞的 APP–CD74 在两个队列均较常复现：GSE132465 20/21（95.2%），GSE200997 7/9（77.8%）；MIF–CD74/CD44 分别为 20/21（95.2%） 和 6/9（66.7%）。GSE200997细分髓系样本数偏少，应作为敏感性支持。
- 在 GSE132465 同患者配对比较中，APP–CD74概率在 C1QC 巨噬细胞高于炎症巨噬细胞：配对n=19，中位差=0.02036，Wilcoxon p=0.000196，BH-FDR=0.0204。相反，本次标签中PID1亚群中位均值为 C1QC 0.0225、炎症巨噬细胞 0.198，因此该通讯差异不支持“PID1高表达状态对应更强APP–CD74”的简单解释。
- 同一比较中 MIF–CD74/CD44 的亚群概率差异不显著：配对n=19，中位差=0.003512，Wilcoxon p=0.446，BH-FDR=0.62。控制亚群身份、患者/样本和细胞数后，APP–CD74与MIF–CD74/CD44的PID1关联 BH-FDR 分别为 0.83 和 0.915；全部身份调整检验中 BH-FDR<0.05 为 0/132。
- 因此，当前数据支持部分上皮–髓系预测通讯在细分亚群后仍可见，但**没有证据证明这些通讯与PID1独立相关**。尤其GSE200997炎症髓系样本覆盖有限，不足以验证其亚群差异。

## 分析问题

沿用原逐患者 CellChat 的输入表达矩阵和肿瘤上皮标签，将旧版宽泛的 `Monocyte/Macrophage/Other_myeloid` 替换为本次 PID1 分析中独立确定的单核/巨噬细胞亚群。只分析肿瘤上皮↔单核/巨噬细胞，不纳入CAF、DC、未解析髓系或T/NK污染。两队列分开分析，单位分别是 GSE132465 患者和 GSE200997 肿瘤样本。

本次重用的表达矩阵是由原始 GEO UMI 数据生成的 CellChat 配体/受体信号基因面板矩阵（各队列1471和1349个基因，均含 APP、MIF、CD74、CD44、PID1），不是全转录组。上皮标签沿用原分析；髓系标签完全换成本次新注释。对照旧分析时仅使用你指定的 CellChat 文件，不读取其他历史分析结果。

最终纳入的表达矩阵/标签对齐细胞为 GSE132465 23,010 个（17,469个肿瘤上皮+5,541个新版单核/巨噬细胞）和 GSE200997 7,400 个（6,871个上皮+529个新版单核/巨噬细胞）；GSE200997有一个上皮样本没有达到任何单核/巨噬亚群的分析门槛。

旧版宽分组的关键结果已从你指定的旧版复现表核对：上皮→单核细胞 APP–CD74 为 GSE132465 19/19、GSE200997 11/12；上皮→巨噬细胞为22/23、11/12。MIF–CD74/CD44 在上皮→单核细胞为19/19、8/12；上皮→巨噬细胞为22/23、7/12。下面的新结果将这些宽标签拆成当前可辨认的单核/巨噬状态，并重新计算各状态的分母和网络。

## CellChat与统计方法

- CellChat 2.2.0.9001，`CellChatDB.human`；每名患者/样本独立建网，`triMean`、`population.size=FALSE`、100次置换。
- 与旧结果相同，CellChat内部过滤要求群体细胞数 **>10**。另按目标亚群每单位≥20和≥50个细胞重新计算覆盖率敏感性；边概率保持同一次患者/样本模型估计，仅改变纳入分母。
- 每个单位内先将置换p值下限设为1/101，再在每个发送者–接收者组合内做BH。`q<0.05`的跨患者比例是复现度描述，不是患者级因果或临床推断检验。
- 亚群间比较以同一患者/样本内相同上皮、不同髓系状态的 CellChat 概率做配对 Wilcoxon；未返回的边按0记为“本次网络未预测到”，不是证明生物学信号绝对不存在。
- PID1关联检验分两层：固定髓系亚群后，跨患者/样本回归通信概率与该亚群PID1均值，并调整亚群和上皮细胞数；以及加入患者/样本固定效应和亚群固定效应的合并模型。后者检验在控制患者背景、亚群身份及细胞数后，患者内不同亚群的PID1差异是否仍与边概率相关。模型仅对多单位重复的候选边运行，并对检验做BH校正。
- 本次所有具备上皮和至少一个目标亚群的单位均成功建网：CellChat运行异常 0 个；异常单位不纳入统计分母。

## 每个亚群的患者/样本覆盖

| 队列 | 亚群 | >10细胞单位数 | ≥20细胞单位数 | ≥50细胞单位数 | 中位细胞数 | PID1中位均值 | PID1中位检出率 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE132465 | Macrophage_C1QC | 21 | 21 | 13 | 88 | 0.0225 | 0.0286 |
| GSE132465 | Macrophage_IFN | 3 | 2 | 2 | 66 | 0.147 | 0.0909 |
| GSE132465 | Macrophage_SPP1_TREM2 | 4 | 3 | 2 | 49.5 | 0.00714 | 0.0082 |
| GSE132465 | Macrophage_inflammatory | 19 | 17 | 10 | 63 | 0.198 | 0.152 |
| GSE132465 | Monocyte_CD14_classical | 9 | 7 | 5 | 56 | 0.208 | 0.146 |
| GSE132465 | Monocyte_CD14_classical_inflammatory | 17 | 12 | 4 | 26 | 0.14 | 0.087 |
| GSE200997 | Macrophage_C1QC | 9 | 4 | 0 | 19 | 0.0348 | 0.0526 |
| GSE200997 | Macrophage_inflammatory | 1 | 0 | 0 | 17 | 0.497 | 0.353 |
| GSE200997 | Monocyte_CD14_classical | 2 | 1 | 1 | 43.5 | 0.251 | 0.159 |
| GSE200997 | Monocyte_CD14_classical_inflammatory | 3 | 2 | 1 | 33 | 0.249 | 0.211 |

> GSE200997没有患者ID，以样本作单位；不少细分状态覆盖很少，未达到门槛的亚群不作“未通讯”的结论。

## 预先关注的 APP/MIF 轴

下表报告 APP–CD74、MIF–CD74或其受体复合体、MIF–CD44 在不同亚群和两个方向中的复现。分母是该上皮–目标亚群具备足够细胞数的患者/样本。

| 队列 | 方向 | 髓系亚群 | 配体_受体 | 合格单位数 | BH-q<0.05单位数 | 复现比例 | 全部合格单位的中位概率 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | APP_CD74 | 21 | 20 | 0.952 | 0.0953 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | MIF_CD74_CD44 | 21 | 20 | 0.952 | 0.152 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | MIF_CD74_CXCR4 | 21 | 17 | 0.81 | 0.0675 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_IFN | APP_CD74 | 3 | 3 | 1 | 0.0558 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_IFN | MIF_CD74_CD44 | 3 | 3 | 1 | 0.133 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_IFN | MIF_CD74_CXCR4 | 3 | 1 | 0.333 | 0.0653 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_SPP1_TREM2 | APP_CD74 | 4 | 4 | 1 | 0.116 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_SPP1_TREM2 | MIF_CD74_CD44 | 4 | 4 | 1 | 0.108 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_SPP1_TREM2 | MIF_CD74_CXCR4 | 4 | 4 | 1 | 0.088 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | APP_CD74 | 19 | 18 | 0.947 | 0.0736 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | MIF_CD74_CD44 | 19 | 18 | 0.947 | 0.127 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | MIF_CD74_CXCR4 | 19 | 15 | 0.789 | 0.0476 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical | APP_CD74 | 9 | 9 | 1 | 0.0678 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical | MIF_CD74_CD44 | 9 | 9 | 1 | 0.129 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical | MIF_CD74_CXCR4 | 9 | 5 | 0.556 | 0.0486 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | APP_CD74 | 17 | 16 | 0.941 | 0.0705 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | MIF_CD74_CD44 | 17 | 16 | 0.941 | 0.143 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | MIF_CD74_CXCR4 | 17 | 2 | 0.118 | 0 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_C1QC | APP_CD74 | 21 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_C1QC | MIF_CD74_CD44 | 21 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_IFN | APP_CD74 | 3 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_IFN | MIF_CD74_CD44 | 3 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_SPP1_TREM2 | APP_CD74 | 4 | 0 | 0 | 0.00492 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_SPP1_TREM2 | MIF_CD74_CD44 | 4 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_inflammatory | APP_CD74 | 19 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_inflammatory | MIF_CD74_CD44 | 19 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Monocyte_CD14_classical | APP_CD74 | 9 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Monocyte_CD14_classical | MIF_CD74_CD44 | 9 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Monocyte_CD14_classical_inflammatory | APP_CD74 | 17 | 0 | 0 | 0 |
| GSE132465 | Myeloid_to_epithelial | Monocyte_CD14_classical_inflammatory | MIF_CD74_CD44 | 17 | 0 | 0 | 0 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | APP_CD74 | 9 | 7 | 0.778 | 0.116 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | MIF_CD74_CD44 | 9 | 6 | 0.667 | 0.13 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | MIF_CD74_CXCR4 | 9 | 2 | 0.222 | 0 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_inflammatory | APP_CD74 | 1 | 1 | 1 | 0.0523 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_inflammatory | MIF_CD74_CD44 | 1 | 1 | 1 | 0.282 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_inflammatory | MIF_CD74_CXCR4 | 1 | 1 | 1 | 0.0976 |
| GSE200997 | Epithelial_to_myeloid | Monocyte_CD14_classical | APP_CD74 | 2 | 2 | 1 | 0.123 |
| GSE200997 | Epithelial_to_myeloid | Monocyte_CD14_classical | MIF_CD74_CD44 | 2 | 1 | 0.5 | 0.0612 |
| GSE200997 | Epithelial_to_myeloid | Monocyte_CD14_classical | MIF_CD74_CXCR4 | 2 | 1 | 0.5 | 0.0948 |
| GSE200997 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | APP_CD74 | 3 | 3 | 1 | 0.0582 |
| GSE200997 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | MIF_CD74_CD44 | 3 | 3 | 1 | 0.244 |
| GSE200997 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | MIF_CD74_CXCR4 | 3 | 2 | 0.667 | 0.102 |
| GSE200997 | Myeloid_to_epithelial | Macrophage_C1QC | APP_CD74 | 9 | 0 | 0 | 0 |
| GSE200997 | Myeloid_to_epithelial | Macrophage_C1QC | MIF_CD74_CD44 | 9 | 0 | 0 | 0 |
| GSE200997 | Myeloid_to_epithelial | Macrophage_inflammatory | APP_CD74 | 1 | 0 | 0 | 0 |
| GSE200997 | Myeloid_to_epithelial | Macrophage_inflammatory | MIF_CD74_CD44 | 1 | 0 | 0 | 0 |
| GSE200997 | Myeloid_to_epithelial | Monocyte_CD14_classical | APP_CD74 | 2 | 0 | 0 | 0 |
| GSE200997 | Myeloid_to_epithelial | Monocyte_CD14_classical | MIF_CD74_CD44 | 2 | 0 | 0 | 0 |
| GSE200997 | Myeloid_to_epithelial | Monocyte_CD14_classical_inflammatory | APP_CD74 | 3 | 0 | 0 | 0 |
| GSE200997 | Myeloid_to_epithelial | Monocyte_CD14_classical_inflammatory | MIF_CD74_CD44 | 3 | 0 | 0 | 0 |

## 其他重复出现的配体–受体边

以下按每亚群至少5个合格单位、至少2个单位内BH-q<0.05筛出；展示排序靠前的边，不作为独立机制验证。

| 队列 | 方向 | 髓系亚群 | 配体_受体 | 合格单位数 | BH-q<0.05单位数 | 复现比例 | 中位概率 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | LAMB3_CD44 | 21 | 21 | 1 | 0.0232 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | LAMB3_CD44 | 19 | 19 | 1 | 0.0286 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | LAMB3_CD44 | 17 | 17 | 1 | 0.0321 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical | MIF_CD74_CD44 | 9 | 9 | 1 | 0.129 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical | APP_CD74 | 9 | 9 | 1 | 0.0678 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical | LAMB3_CD44 | 9 | 9 | 1 | 0.0336 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | MIF_CD74_CD44 | 21 | 20 | 0.952 | 0.152 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | APP_CD74 | 21 | 20 | 0.952 | 0.0953 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | APP_TREM2_TYROBP | 21 | 20 | 0.952 | 0.0373 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | MIF_CD74_CD44 | 19 | 18 | 0.947 | 0.127 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | APP_CD74 | 19 | 18 | 0.947 | 0.0736 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | MIF_CD74_CD44 | 17 | 16 | 0.941 | 0.143 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | APP_CD74 | 9 | 7 | 0.778 | 0.116 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | MIF_CD74_CD44 | 9 | 6 | 0.667 | 0.13 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | APP_TREM2_TYROBP | 9 | 6 | 0.667 | 0.0118 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | MDK_ITGA4_ITGB1 | 9 | 5 | 0.556 | 0.000675 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | CD99_CD99 | 9 | 4 | 0.444 | 0 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | CD99_PILRA | 9 | 4 | 0.444 | 0 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | MDK_LRP1 | 9 | 4 | 0.444 | 0 |
| GSE200997 | Myeloid_to_epithelial | Macrophage_C1QC | CD99_CD99 | 9 | 4 | 0.444 | 0 |
| GSE200997 | Myeloid_to_epithelial | Macrophage_C1QC | LGALS9_CD44 | 9 | 2 | 0.222 | 0.00197 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | LAMB3_CD44 | 9 | 2 | 0.222 | 0 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | LTE4-DPEP1_CYSLTR1 | 9 | 2 | 0.222 | 0 |
| GSE200997 | Epithelial_to_myeloid | Macrophage_C1QC | MIF_CD74_CXCR4 | 9 | 2 | 0.222 | 0 |

## 同患者亚群对照

重点查看经典与炎症单核细胞，以及 C1QC 与炎症巨噬细胞。差值为亚群A减亚群B；概率缺失按0。表中只显示配对单位数≥5且边在至少两个单位中达到组内BH阈值、或属于预先关注的APP/MIF轴的比较。

| 队列 | 方向 | 配体_受体 | 亚群A | 亚群B | 配对单位数 | 中位概率差（A−B） | Wilcoxon p | BH-FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE132465 | Epithelial_to_myeloid | APP_CD74 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.0204 | 0.000196 | 0.0204 |
| GSE132465 | Epithelial_to_myeloid | APP_TREM2_TYROBP | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.0316 | 0.000196 | 0.0204 |
| GSE132465 | Epithelial_to_myeloid | HLA-DRA_CD4 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00579 | 0.000438 | 0.0204 |
| GSE132465 | Epithelial_to_myeloid | LAMA3_CD44 | Macrophage_C1QC | Macrophage_inflammatory | 19 | -0.00146 | 0.000438 | 0.0204 |
| GSE132465 | Epithelial_to_myeloid | LAMB3_CD44 | Macrophage_C1QC | Macrophage_inflammatory | 19 | -0.00539 | 0.000132 | 0.0204 |
| GSE132465 | Epithelial_to_myeloid | LAMC2_CD44 | Macrophage_C1QC | Macrophage_inflammatory | 19 | -0.00255 | 0.000438 | 0.0204 |
| GSE132465 | Epithelial_to_myeloid | LGALS9_CD44 | Macrophage_C1QC | Macrophage_inflammatory | 19 | -0.00151 | 0.000196 | 0.0204 |
| GSE132465 | Epithelial_to_myeloid | LGALS9_HAVCR2 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.000899 | 0.000196 | 0.0204 |
| GSE132465 | Myeloid_to_epithelial | PLAU_PLAUR | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.0209 | 0.000155 | 0.0204 |
| GSE132465 | Myeloid_to_epithelial | TNFSF12_TNFRSF12A | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00336 | 0.000293 | 0.0204 |
| GSE132465 | Epithelial_to_myeloid | MDK_LRP1 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00268 | 0.000982 | 0.025 |
| GSE132465 | Epithelial_to_myeloid | CD99_CD99 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00479 | 0.00123 | 0.0286 |
| GSE132465 | Myeloid_to_epithelial | CD99_CD99 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00479 | 0.00123 | 0.0286 |
| GSE132465 | Myeloid_to_epithelial | PPIA_BSG | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.0144 | 0.00123 | 0.0286 |
| GSE132465 | Epithelial_to_myeloid | ANXA1_FPR3 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00352 | 0.00222 | 0.0466 |
| GSE132465 | Epithelial_to_myeloid | ANXA1_FPR1 | Macrophage_C1QC | Macrophage_inflammatory | 19 | -2.79e-05 | 0.00585 | 0.0769 |
| GSE132465 | Epithelial_to_myeloid | CCL3_CCR1 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.000226 | 0.00444 | 0.0769 |
| GSE132465 | Epithelial_to_myeloid | GAS6_AXL | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00036 | 0.00506 | 0.0769 |
| GSE132465 | Epithelial_to_myeloid | GAS6_MERTK | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0 | 0.00769 | 0.0769 |
| GSE132465 | Epithelial_to_myeloid | HLA-DPB1_CD4 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00018 | 0.00506 | 0.0769 |
| GSE132465 | Epithelial_to_myeloid | HLA-DQB1_CD4 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.000921 | 0.00506 | 0.0769 |
| GSE132465 | Epithelial_to_myeloid | HLA-DRB1_CD4 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.000196 | 0.00506 | 0.0769 |
| GSE132465 | Myeloid_to_epithelial | CD55_ADGRE5 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0 | 0.00769 | 0.0769 |
| GSE132465 | Myeloid_to_epithelial | LGALS9_P4HB | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.0127 | 0.00523 | 0.0769 |
| GSE132465 | Myeloid_to_epithelial | RETN_CAP1 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0 | 0.00769 | 0.0769 |
| GSE132465 | Myeloid_to_epithelial | SPP1_CD44 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00314 | 0.00836 | 0.0826 |
| GSE132465 | Myeloid_to_epithelial | LGALS9_CD44 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0.00259 | 0.00861 | 0.0841 |
| GSE132465 | Epithelial_to_myeloid | CCL3L3_CCR1 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0 | 0.0117 | 0.0903 |
| GSE132465 | Epithelial_to_myeloid | COL1A1_CD44 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0 | 0.0117 | 0.0903 |
| GSE132465 | Epithelial_to_myeloid | Desmosterol-DHCR7_NR1H3 | Macrophage_C1QC | Macrophage_inflammatory | 19 | 0 | 0.0117 | 0.0903 |

## PID1是否独立关联通讯

“亚群内模型”固定一个髓系状态，在该状态各患者/样本间回归边概率与PID1均值，并调整髓系和上皮细胞数。身份调整模型再加入患者/样本与髓系亚群固定效应。共运行 452 个亚群内模型和 132 个身份调整模型；BH后显著数分别为 0 和 0。

### 固定亚群后的跨患者/样本关联

| 队列 | 方向 | 髓系亚群 | 配体_受体 | 单位数 | 正概率单位数 | PID1系数（每SD） | p | BH-FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE200997 | Myeloid_to_epithelial | Macrophage_C1QC | PPIA_BSG | 9 | 2 | 0.0437 | 0.00463 | 0.167 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | PPIA_BSG | 17 | 5 | 0.0094 | 0.49 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_inflammatory | APP_SORL1 | 19 | 3 | -9.61e-05 | 0.138 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Monocyte_CD14_classical_inflammatory | THBS1_CD47 | 17 | 11 | 0.00516 | 0.349 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_inflammatory | APP_CD74 | 19 | 1 | 0.000138 | 0.396 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_C1QC | APP_CD74 | 21 | 1 | -0.000668 | 0.441 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | VEGFA_FLT1 | 17 | 1 | -0.000139 | 0.362 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical | VEGFA_FLT1 | 9 | 1 | 7.78e-05 | 0.655 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | VEGFA_FLT1 | 19 | 2 | -0.000158 | 0.144 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | THBS1_CD47 | 17 | 3 | 0.000751 | 0.307 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_C1QC | THBS1_ITGA3_ITGB1 | 21 | 2 | 2.71e-05 | 0.445 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | THBS1_CD47 | 19 | 2 | -0.000151 | 0.438 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | THBS1_CD47 | 21 | 3 | 0.000392 | 0.126 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical_inflammatory | THBS1_CD36 | 17 | 2 | 0.000584 | 0.432 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Monocyte_CD14_classical | THBS1_CD36 | 9 | 1 | -0.000104 | 0.716 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | THBS1_CD36 | 19 | 1 | -0.000207 | 0.691 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | THBS1_CD36 | 21 | 2 | 0.000678 | 0.342 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_inflammatory | SPP1_ITGA5_ITGB1 | 19 | 5 | -0.000256 | 0.273 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | SPP1_ITGA5_ITGB1 | 21 | 3 | 0.000111 | 0.457 | 0.964 |
| GSE132465 | Epithelial_to_myeloid | Macrophage_C1QC | SEMA3C_PLXND1 | 21 | 5 | 3.98e-05 | 0.661 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_C1QC | APP_TNFRSF21 | 21 | 7 | 0.000109 | 0.51 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_C1QC | AREG_EGFR | 21 | 2 | -3.7e-05 | 0.24 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Macrophage_inflammatory | AREG_EGFR | 19 | 3 | 0.000532 | 0.201 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Monocyte_CD14_classical | AREG_EGFR | 9 | 1 | -0.000115 | 0.469 | 0.964 |
| GSE132465 | Myeloid_to_epithelial | Monocyte_CD14_classical_inflammatory | CD99_PILRB | 17 | 3 | -3.76e-05 | 0.636 | 0.964 |

### 同时调整亚群身份和患者/样本后的关联

| 队列 | 方向 | 配体_受体 | 单位数 | 亚群数 | 正概率单位数 | PID1调整系数（每SD） | p | BH-FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE132465 | Epithelial_to_myeloid | ANXA1_FPR1 | 21 | 6 | 38 | 0.00111 | 0.439 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | CD99_CD99 | 21 | 6 | 38 | -0.00223 | 0.468 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | CD55_ADGRE5 | 21 | 6 | 36 | 0.0011 | 0.349 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | CD46_JAG1 | 21 | 6 | 10 | -2.81e-05 | 0.472 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | CD209_CEACAM1 | 21 | 6 | 9 | 0.000284 | 0.109 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | AREG_EGFR_ERBB2 | 21 | 6 | 10 | 0.000257 | 0.313 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | AREG_EGFR | 21 | 6 | 10 | 0.000171 | 0.293 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | APP_TNFRSF21 | 21 | 6 | 13 | -0.000265 | 0.395 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | APP_SORL1 | 21 | 6 | 12 | -0.000201 | 0.196 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | APP_CD74 | 21 | 6 | 5 | -0.0022 | 0.324 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | VEGFA_FLT1 | 21 | 6 | 4 | -3.57e-05 | 0.356 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | THBS1_CD47 | 21 | 6 | 11 | -0.00037 | 0.388 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | SAA1_FPR2 | 21 | 6 | 3 | 0.000337 | 0.364 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | PROS1_MERTK | 21 | 6 | 3 | -3.12e-05 | 0.439 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | PGE2-PTGES3_PTGER4 | 21 | 6 | 36 | -0.00144 | 0.469 | 0.83 |
| GSE132465 | Myeloid_to_epithelial | TNFSF14_LTBR | 21 | 6 | 9 | 0.000154 | 0.448 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | PGE2-PTGES2_PTGER4 | 21 | 6 | 34 | -0.000429 | 0.424 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | PGE2-PTGES2_PTGER2 | 21 | 6 | 14 | 0.000213 | 0.289 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | NAMPT_ITGA5_ITGB1 | 21 | 6 | 21 | -0.000737 | 0.356 | 0.83 |
| GSE132465 | Epithelial_to_myeloid | MDK_NCL | 21 | 6 | 43 | 0.000832 | 0.468 | 0.83 |

表格中的系数代表观察性关联。即使某条通讯边和PID1在统计上相关，也不能证明PID1调控配体/受体或实际细胞通讯。若身份调整后关联消失，而某亚群边更稳定，则更符合“亚群身份/细胞组成可解释该差异”的结果；若关联仍在，也只能视为后续实验的候选线索。

## 结果文件

- [`PID1_CellChat_APP_MIF_candidate_recurrence.tsv`](PID1_CellChat_APP_MIF_candidate_recurrence.tsv)：重点通路按亚群、方向、阈值的复现率。
- [`PID1_CellChat_pair_eligibility_model_units.tsv`](PID1_CellChat_pair_eligibility_model_units.tsv)：实际进入统计分母的逐患者/样本亚群覆盖。
- [`PID1_CellChat_edge_recurrence_by_subtype_threshold.tsv.gz`](PID1_CellChat_edge_recurrence_by_subtype_threshold.tsv.gz)：所有配体–受体边与10/20/50细胞覆盖敏感性。
- [`PID1_CellChat_paired_subtype_edge_contrasts.tsv.gz`](PID1_CellChat_paired_subtype_edge_contrasts.tsv.gz)：同患者/样本亚群配对比较。
- [`PID1_CellChat_PID1_association_within_subtype.tsv`](PID1_CellChat_PID1_association_within_subtype.tsv)：亚群内部的PID1关联模型。
- [`PID1_CellChat_PID1_association_adjusted_for_subtype_identity.tsv`](PID1_CellChat_PID1_association_adjusted_for_subtype_identity.tsv)：同时调整亚群身份与患者/样本的模型。
- [`GSE132465_PID1_CellChat_patient_sample_interactions.tsv.gz`](GSE132465_PID1_CellChat_patient_sample_interactions.tsv.gz)、[`GSE200997_PID1_CellChat_patient_sample_interactions.tsv.gz`](GSE200997_PID1_CellChat_patient_sample_interactions.tsv.gz)：逐单位全边与概率、p/q值。
- [`GSE132465_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz`](GSE132465_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz)、[`GSE200997_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz`](GSE200997_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz)：本次CellChat实际纳入的逐细胞标签。
- `APP_MIF_subtype_recurrence.png`：预先关注轴的逐亚群复现图。
- 若存在可估计的身份调整模型，`PID1_adjusted_communication_associations.png` 汇总其系数与区间。

## 解释限制

细胞通讯由配体/受体表达模式预测，不能证明蛋白结合、分泌或方向性信号已发生。GSE200997的上皮标签来自先前的多marker推断，且细分髓系状态单位覆盖较少；它是敏感性分析。CellChat概率也会受表达检测率、每亚群细胞数及样本组成影响；因此报告纳入细胞数、阈值敏感性和身份调整模型。此次分析不支持因果机制结论，也不扩展到机器学习。
