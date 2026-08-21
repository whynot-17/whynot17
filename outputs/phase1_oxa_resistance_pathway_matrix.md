# OXA-resistant CRC 的跨模型 pathway × dataset 矩阵

## 当前结论

在 6 个“细胞获得性 OXA-resistance”对照中，FAO 不是共识性状态；更稳定的是一个组合状态：

1. **pyrimidine metabolism score 下降：6/6 个细胞模型方向一致**；
2. **carnitine-entry score 下降：5/6**，但它仍然只是转录层面的 transport/entry proxy，不等于已经证明 carnitine dependency；
3. **UPR/ER-stress score 下降：5/6**；
4. **EMT score 上升：5/6**；
5. **apoptosis gene-set score 上升：5/6**，这更可能表示 apoptotic priming/stress-related transcription，而不是“耐药细胞更容易发生功能性凋亡”。

NRF2、ferroptosis-resistance、OXPHOS、TNF/NFκB、ABC transport、glycolysis、lipid/FA metabolism 和 drug metabolism 大多为 4/6 同向，属于**中等稳定、需要继续验证的模块**。GSH/redox、ROS、DNA repair、autophagy、TGFβ、IL6/JAK/STAT3、purine metabolism 只有 3/6 同向，当前不能称作跨模型共识。

因此，目前最合理的 working model 不是“FAO-high OXA-resistant CRC”，而是：

> **OXA-resistant CRC 的跨模型核心更接近：EMT/状态转移 + pyrimidine-metabolism remodeling + UPR attenuation；carnitine-entry reduction 是一个有潜力连接到 Meldonium 的代谢脆弱性候选，而不是已经确立的 universal state。**

## 主要 pathway × dataset matrix

矩阵中的数值是：

\[
\Delta pathway = mean(score_{OXA-R}) - mean(score_{parental})
\]

其中每个数据集内部先对每个基因做 gene-wise z-score，再对该 pathway 中可测基因取均值。正值表示 OXA-R 侧更高，负值表示更低。它是跨平台可比较的标准化 effect-size proxy，不是原始表达量，也不是 pathway activation 的直接实验测量。

完整 CSV：

- [cell-line pathway × dataset matrix](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase1_pathway_by_dataset_matrix_cell_lines.csv)
- [primary matrix including xenograft](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase1_pathway_by_dataset_matrix_primary.csv)
- [all-context long table](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase1_pathway_effects_long_all_contexts.csv)
- [cell-line stability summary](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase1_pathway_stability_summary_cell_lines.csv)
- [all-primary stability summary](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase1_pathway_stability_summary_primary.csv)

主分析包含：GSE77932 的 HCT116/DLD1、GSE42387 的 HCT116/HT29/LoVo、GSE119603 的 HCT116；另把 GSE124808 的 HCT116 resistant xenograft 单独放入 primary matrix，但它的 parental 样本数为 1，不能与细胞三重复等权解释。数据来自对应的 [GEO series records](https://www.ncbi.nlm.nih.gov/geo/)。

## 跨细胞模型稳定性摘要

| pathway | 同向比例 | 方向 | median Δ | 解释 |
|---|---:|---|---:|---|
| pyrimidine metabolism | 6/6 | ↓ | -0.128 | 当前最稳定的转录模块，但需排除增殖/细胞周期混杂 |
| carnitine entry | 5/6 | ↓ | -0.105 | 对 Meldonium 最直接的候选连接；仍需依赖性与 rescue 验证 |
| UPR/ER stress | 5/6 | ↓ | -0.057 | 状态性变化，需区分真正 UPR attenuation 与平台/生长差异 |
| EMT | 5/6 | ↑ | +0.066 | 跨模型最稳定的耐药伴随状态之一 |
| apoptosis | 5/6 | ↑ | +0.065 | 只能解读为 apoptosis-related transcriptional state |
| NRF2 response | 4/6 | ↑ | +0.266 | 效应较大但一致性不足，适合做 subgroup/interaction |
| ferroptosis resistance | 4/6 | ↑ | +0.190 | 与 NRF2 同步，但不是 universal |
| OXPHOS | 4/6 | ↑ | +0.124 | 方向存在异质性，不能直接作为统一代谢依赖 |
| ABC transport | 4/6 | ↑ | +0.081 | 中等稳定；需避免将表达上升直接等同于药物外排功能 |
| GSH/redox | 3/6 | — | +0.039 | 当前不支持跨模型共识 |
| DNA repair | 3/6 | — | +0.033 | 当前不支持跨模型共识 |
| FAO mitochondrial | 3/6 | — | -0.068 | 明确不是稳定的 FAO-high state |
| TGFβ / IL6-JAK-STAT3 | 3/6 | — | +0.017 / -0.005 | 需要模型分层，而非总平均 |

## 外部验证层

- **GSE30011：OXA-resistant vs OXA-sensitive 的横断面细胞系**。它不是获得性耐药对照，因此不能放入主稳定性统计；它作为外部方向验证。当前结果：GSH/redox -0.310、carnitine entry -0.229、FAO mitochondrial -0.077、EMT +0.231、DNA repair +0.237、ROS +0.324。它支持“carnitine-entry 下降”和“EMT/ROS/DNA-repair 方向上升”的可能性，但不能单独证明因果。
- **GSE83129：患者 OXA non-responder vs responder**。这是临床 response layer，不是细胞获得性耐药层。当前结果：ferroptosis-resistance -0.270、pyrimidine metabolism +0.246、EMT +0.163、UPR -0.140、OXPHOS -0.108；方向与细胞模型并不完全一致，提示患者 response 与体外 acquired resistance 不是同一个状态。

外部数据也来自 [GEO](https://www.ncbi.nlm.nih.gov/geo/)。Hallmark 和 Reactome gene sets 来自 [MSigDB 2026.1.Hs release](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/)；MSigDB 的 gene-set 文件格式说明见 [官方文档](https://docs.gsea-msigdb.org/GSEA/Data_Formats/)。

## 对 Meldonium 课题的直接影响

这一步以后，Meldonium 不应再围绕“抑制 CPT1A”来写。更稳妥的生信假说是：

> **Meldonium 可能优先作用于一个在多数 OXA-R 细胞模型中降低的 carnitine-entry axis，并进一步影响 FAO/lipid redox buffering；但它是否能逆转真正的 OXA-R phenotype，必须依赖后续药物扰动数据。**

下一步应做三件事：

1. 用独立的 OXA-R 转录组或单细胞/公共药物扰动数据确认 `carnitine_entry` 的方向和可重复性；
2. 把每个 OXA-R 模型按 `EMT-high / NRF2-ferroptosis-high / carnitine-low` 做 subgroup，而不是把所有模型硬平均；
3. 在 Meldonium 的 drug-perturbation signature 中检验：它是否使 `carnitine_entry`、EMT、NRF2/ferroptosis 等模块发生反向改变，并以反向一致性构建 reversal score。

## 重要限制

这张矩阵回答的是“哪些转录模块跨模型方向更稳定”，不回答“哪些模块是耐药的必要因果节点”。尤其是 apoptosis、EMT、UPR、pyrimidine metabolism 都可能受增殖速度、培养条件、克隆选择和平台差异影响。真正进入论文主结论前，至少要增加：独立数据集复现、基因集替换敏感性分析、留一模型交叉验证，以及临床 response layer 的独立解释。
