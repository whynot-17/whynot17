# Phase 5：扰动签名与 OXA-R reversal score

## 结论先行

本轮结果**不支持**把 Leflunomide 或 Teriflunomide 直接升级为“salvage-low/DHODH-high OXA-R 亚型的选择性药物”。两者的反转分数处于中等水平，且相对于 matched-negative OXA-R 状态的选择性很弱。Bortezomib 的反转分数和 DHODH-network score 反而更高，但它只是机制阳性对照，且 proteasome/oxaliplatin 组合已有临床探索，因此不能作为新药物发现结论。

当前最稳妥的判断是：

> **DHODH compensation 的生物学假设仍可保留，但 Leflunomide/Teriflunomide 的扰动签名没有显示出足够的 DHODH-subtype selectivity。Phase 5 尚未筛出可以直接进入湿实验的 Drug X。**

## 数据与评分

使用六个 acquired OXA-R CRC cell-line contrasts：

`GSE77932|HCT116`、`GSE77932|DLD1`、`GSE42387|HCT116`、`GSE42387|HT29`、`GSE42387|LoVo`、`GSE119603|HCT116`。

扰动数据来自 SigCom/LINCS L1000 characteristic-direction signatures，优先选取 HT29、24 h 签名：Leflunomide 12 个、Teriflunomide 12 个、Bortezomib 24 个。Bortezomib 被限制为代表性子集，避免因为沉积签名数量较多而获得机械性优势。Meldonium 没有进入主排序，也没有使用同靶点药物替代。

定义：

`RRS = (1 − Spearman(OXA-R Δ, drug perturbation coefficient)) / 2`

越高表示药物扰动越倾向于反向 OXA-R 状态。全局分数使用疾病共识的 top 250 up + top 250 down 基因；网络分数使用 Phase 4 的补偿性网络；去 DHODH 分析只从该网络中删除 DHODH。

## 主要结果

### 全局 OXA-R reversal

基于各药物的 HT29 扰动签名中位数：

| Drug | Global RRS | Global target-network RRS | Per-signature RRS median |
|---|---:|---:|---:|
| Bortezomib | 0.619 | 0.645 | 0.527 |
| Leflunomide | 0.609 | 0.492 | 0.522 |
| Teriflunomide | 0.602 | 0.484 | 0.524 |

药物之间差异不大；逐签名中位数几乎重叠。因此不能把 0.619 versus 0.609 解读为可靠的药物优先级差异。

还要注意，L1000 signature 与全转录组 OXA-R 共识的实际重叠在全局评分中约为 65 个基因；因此这些数字是方向性证据，不是高维全转录组的完整验证。

### 亚型结果

| Drug | salvage-low/DHODH-high | salvage-low/RRM2-high | UPR-low/ERAD-high |
|---|---:|---:|---:|
| Leflunomide | 0.591 | 0.448 | 0.617 |
| Teriflunomide | 0.548 | 0.508 | 0.603 |
| Bortezomib | 0.634 | 0.597 | 0.611 |

Leflunomide/Teriflunomide 并没有在 DHODH 亚型上形成清晰的第一名优势；Leflunomide 的最高分反而出现在 UPR-low/ERAD-high 亚型。

### Matched-negative selectivity

| Drug | Target subtype | Target RRS | Matched-negative RRS | Difference |
|---|---|---:|---:|---:|
| Leflunomide | salvage-low/DHODH-high | 0.591 | 0.549 | +0.042 |
| Teriflunomide | salvage-low/DHODH-high | 0.548 | 0.504 | +0.044 |
| Bortezomib | salvage-low/DHODH-high | 0.634 | 0.438 | +0.196 |

预设的“DHODH inhibitor 在 DHODH-high 亚型选择性增强”目前没有得到强支持。Leflunomide/Teriflunomide 的 +0.04 级别差异更像弱趋势，不能作为进入湿实验的核心依据。

### Remove-DHODH test

| Drug | Network RRS | Remove-DHODH RRS | Change |
|---|---:|---:|---:|
| Leflunomide | 0.455 | 0.457 | +0.003 |
| Teriflunomide | 0.466 | 0.467 | +0.002 |
| Bortezomib | 0.599 | 0.596 | −0.003 |

删除 DHODH 后排序基本不变。这个结果说明当前 network score 不是由一个 DHODH 单点驱动；但它也意味着“Leflunomide 命中 DHODH”并没有在扰动层面形成很强的网络特异性。

### Leave-one-dataset-out

在 salvage-low/DHODH-high 亚型中，Leflunomide 的 leave-one-dataset-out RRS 为 0.447–0.578，Teriflunomide 为 0.493–0.562，Bortezomib 为 0.538–0.630。方向没有完全崩溃，但波动足以说明目前只能称为探索性结果。

UPR-low/ERAD-high 亚型主要由 GSE42387 与 GSE119603 构成；去掉 GSE119603 后只剩一个模型，因此不能把该亚型的稳定性写成跨数据集验证。

## Novelty audit

“Leflunomide/Teriflunomide 在 CRC 中是全新作用”这个表述现在不能使用。

- Leflunomide 已经被用于 CRC 肝转移/缺氧生长与 DHODH 抑制研究，并且对 5-FU-resistant CRC 的转移定植有过机制和动物证据：[Yamaguchi et al., eLife](https://pmc.ncbi.nlm.nih.gov/articles/PMC7299340/)。
- Teriflunomide 已有 CRC 免疫检查点调节研究；DHODH 抑制剂也已经在 5-FU-resistant CRC 中被研究：[5-FU-resistant CRC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11152977/)、[teriflunomide/PD-1–PD-L1 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13002968/)。
- DHODH 与 OXA resistance 的关系也已经出现直接机制证据：NSUN2–DHODH 轴通过 ferroptosis resistance 促进 CRC OXA resistance：[PubMed 41808866](https://pubmed.ncbi.nlm.nih.gov/41808866/)。
- Bortezomib 与 oxaliplatin 的联合已经进入晚期实体瘤临床一期探索：[PubMed 18376220](https://pubmed.ncbi.nlm.nih.gov/18376220/)。

因此，真正可能成立的新颖点只能收窄为：

> **跨多个 acquired OXA-R 模型的数据整合，识别 salvage/UPR/EMT 组合状态，并检验现有药物扰动是否对某一 OXA-R 亚型具有选择性反转。**

这不是“发现了一个从未用于 CRC 的药”，而是“提出了一个 OXA-R 分层和药物匹配框架”。正式投稿前仍需做 PubMed、Embase、ClinicalTrials.gov 和专利的组合检索。

## 当前决策

1. **Leflunomide/Teriflunomide：降级为 computational comparator，不进入主湿实验候选。** 除非后续加入更多 perturbation 数据或 pharmacogenomic sensitivity 后，DHODH 亚型选择性明显增强。
2. **Bortezomib：保留为机制阳性对照，不作为 novelty lead。** 它可用于确认 proteostasis/NF-κB/proteasome 方向是否具有可重复的反转信号。
3. **Meldonium：正式放入 appendix/failed candidate。** 本轮不再投入主分析资源。
4. 下一步不是继续解释 Leflunomide，而是扩大“已上市、非肿瘤、具有真实 HT29/CRC perturbation signature”的候选池，并要求候选同时满足：subtype selectivity、leave-one-dataset-out 稳定、matched-negative 失败、以及独立药敏数据支持。

本轮 perturbation signature 来自 HT29，不是 OXA-resistant HT29；因此它回答的是“药物扰动是否反向 OXA-R 状态”，还没有回答“药物在 OXA-R 细胞中是否真的更敏感”。后者必须由 DepMap/PRISM/GDSC 药敏或独立实验数据补上。

## 输出文件

- [drug_rrs_global_subtype_and_remove_dhodh.csv](drug_rrs_global_subtype_and_remove_dhodh.csv)
- [drug_subtype_selectivity.csv](drug_subtype_selectivity.csv)
- [per_model_drug_rrs.csv](per_model_drug_rrs.csv)
- [per_signature_global_rrs.csv](per_signature_global_rrs.csv)
- [leave_one_dataset_out.csv](leave_one_dataset_out.csv)
- [perturbation_signature_metadata.csv](perturbation_signature_metadata.csv)
