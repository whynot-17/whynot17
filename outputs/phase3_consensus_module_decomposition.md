# Phase 3：Consensus-module decomposition + vulnerability inference

## 结论先行

这一步没有把原来的结论推翻，但把它精确化了：

> **“pyrimidine metabolism↓”并不是 de novo pyrimidine synthesis 全面关闭，而更像是 catabolism/salvage 下降、interconversion 部分重排；UPR↓更偏向 PERK/ATF4 与 ATF6/ER-proteostasis 支路；EMT↑主要由 mesenchymal gain 支撑，epithelial loss 不够稳定。**

因此，目前最值得筛的不是“抑制整个 pyrimidine/UPR/EMT pathway”的药，而是：

> **OXA-R 在关闭或削弱某些输入/周转支路后，是否被迫依赖仍然升高的补偿节点。**

## 1. Pyrimidine decomposition

| 子模块 | 6 个细胞模型方向 | median Δ | 当前解释 |
|---|---:|---:|---|
| 原始 pyrimidine composite | 0/6 上升，6/6 下降 | -0.128 | 这是最稳定的表型，但原始集合主要由 salvage/catabolism 构成 |
| pyrimidine catabolism | 1/6 上升，5/6 下降 | -0.136 | 是原始下降信号的主要来源之一 |
| pyrimidine salvage | 2/6 上升，4/6 下降 | -0.062 | 不是完全关闭，而是发生重排 |
| de novo core：CAD–DHODH–UMPS–CTPS1/2 | 3/6 上升，3/6 下降 | +0.014 | 不支持“de novo synthesis 全面下降” |
| interconversion | 4/6 上升，2/6 下降 | +0.275 | 可能是耐药状态下的补偿/再分配层 |

基因层面最值得注意的是：

- **DHODH：5/6 上升**；
- **RRM2：5/6 上升**；
- **DCTD：5/6 上升**；
- **PUDP：3/3 下降，但可测模型数较少**；
- **TYMP、UCK1、TK1：多数模型下降**；
- **DCK：5/6 上升，提示 salvage 并非整体关闭，而是选择性重编程。**

这意味着目前不能把 DHODH 直接降级。相反，逻辑已经从“DHODH 低表达，所以抑制它没有意义”变成：

> **DHODH/CTPS/RRM2 等仍然维持或升高，可能是 OXA-R 对低 salvage/低 catabolism 状态的补偿依赖候选。**

但这仍然只是 vulnerability hypothesis，不是因果证据。下一步应检验 OXA-R 是否比 parental 对 DHODH、CAD、CTPS 或 RRM2 抑制更敏感，以及 uridine/cytidine rescue 是否能够救回表型。

## 2. UPR decomposition

| 支路 | 6 个细胞模型方向 | median Δ | 当前解释 |
|---|---:|---:|---|
| Hallmark UPR total | 1/6 上升，5/6 下降 | -0.057 | 总体下降，但不是所有模型一致 |
| PERK–eIF2A–ATF4 | 2/6 上升，4/6 下降 | -0.312 | 有下降趋势，效应量较大但模型间异质性明显 |
| IRE1–XBP1 | 3/6 上升，3/6 下降 | +0.080 | 不是稳定下降支路；仅靠总 XBP1 RNA 也不能判断 XBP1 splicing |
| ATF6/ER-proteostasis | 1/6 上升，5/6 下降 | -0.096 | 方向上最接近稳定下降，但效应受 DLD1 影响 |

同时，部分 ERAD/proteostasis genes 反而偏上：**DERL1、MANF、EDEM1、HSP90B1** 在多数模型中上升。这更像：

> **UPR transcriptional response 被压低，但残余 ER quality-control/ERAD 负荷可能被保留或增强。**

所以不应直接选择 UPR inhibitor。更合理的脆弱性方向是测试：

- proteostasis overload；
- ERAD/HSP90/proteasome stress；
- PERK/ATF4 下降后对额外蛋白折叠压力的耐受性。

## 3. EMT decomposition

| 子模块 | 6 个细胞模型方向 | median Δ | 当前解释 |
|---|---:|---:|---|
| EMT Hallmark | 5/6 上升 | +0.066 | 跨模型较稳定 |
| mesenchymal gain | 5/6 上升 | +0.076 | 主要支撑 EMT 结论 |
| epithelial state | 2/6 上升，4/6 下降 | -0.057 | 有 epithelial loss 趋势，但远不如 mesenchymal gain 稳定 |

基因层面：

- **FN1：6/6 上升**；
- **TJP1：6/6 下降**；
- **MMP2、ITGB1：5/6 上升**；
- VIM、CDH1、CDH2、SNAI1 等方向并不统一。

因此 EMT 不能简单写成“所有 epithelial markers 下降 + 所有 mesenchymal markers 上升”。目前更准确的表述是：

> **OXA-R 更稳定地获得 mesenchymal/adhesion-remodeling program，而 epithelial loss 是模型依赖的。**

## 4. 三个模块是否共存于同一批模型？

用原始矩阵中的三个模块定义：

- EMT Δ > 0；
- pyrimidine Δ < 0；
- UPR Δ < 0。

结果是 **4/6 个细胞模型满足三者的方向性共存**：GSE119603-HCT116、GSE42387-HCT116、GSE42387-HT29、GSE77932-HCT116。

两个例外很有价值：

- **GSE42387-LoVo**：pyrimidine↓、UPR↓，但 EMT Hallmark↓；
- **GSE77932-DLD1**：EMT↑、pyrimidine↓，但 UPR↑。

因此可以提出一个 composite state，但暂时不能说是全部 OXA-R 的单一 universal state。模型间 Spearman 相关也提示：EMT 与 pyrimidine 为中等反向关系，而 EMT 与 UPR 几乎不相关；三者更像“部分共现的状态组合”，不是一个完全耦合的单轴程序。

## 5. 是否只是增殖变慢造成的 pyrimidine↓？

目前不能直接排除。pyrimidine composite 与 proliferation/E2F/G2M/MYC score 在 6 个模型中有中等至较强相关趋势，但 proliferation score 本身并没有 6/6 下降；它在 4/6 模型中反而上升。因此当前最稳妥的判断是：

> **pyrimidine↓ 可能包含 proliferation/dNTP-demand 成分，但不能被简单解释为“所有耐药细胞都变慢”。**

另一个有利于该信号的结果是：原始 pyrimidine gene set 与 proliferation gene set 只有 3 个直接重叠基因；去掉这 3 个重叠基因后，pyrimidine score 仍然是 6/6 下降，而且 median Δ 从 -0.128 加强到 -0.209。这说明它不是 gene-set 直接重叠造成的数学假象，但仍不能排除耐药状态整体重编程或生长动力学造成的间接混杂。

下一步应使用：

1. 以 E2F/G2M score 为协变量的 residual pyrimidine score；
2. 只用 metabolic genes、排除 cell-cycle-overlapping genes 的敏感性分析；
3. 对 parental/OXA-R 做 matched growth-rate 或 doubling-time 校正；
4. 独立 OXA-R 数据集验证 residual signal。

## 6. 当前的 vulnerability inference

| 观察到的状态 | 不应直接做的事 | 更合理的候选脆弱性 |
|---|---|---|
| salvage/catabolism↓ | 直接抑制整个 pyrimidine pathway | 依赖仍维持/升高的 de novo/interconversion compensation；优先验证 DHODH、CTPS、RRM2，并做 nucleoside rescue |
| DHODH、RRM2、DCTD 多数模型↑ | 直接宣布 DHODH 是靶点 | 作为 OXA-R-selective dependency 候选，先做 DepMap/药敏/CRISPR 和 rescue |
| PERK/ATF6 下降 + ERAD genes 上升 | 直接使用 UPR inhibitor | proteostasis overload、ERAD/HSP90/proteasome stress sensitivity |
| mesenchymal gain 稳定 | 直接找“EMT inhibitor” | 用 EMT-reversal signature 筛非肿瘤上市药，并要求同时逆转 OXA-R phenotype |

## 7. 对课题主线的更新

当前最值得保留的主线不是：

> Meldonium → FAO↓ → OXA resistance reversal

而是：

> **OXA-R → mesenchymal gain + pyrimidine salvage/catabolism remodeling + partial UPR attenuation → compensatory metabolic/proteostasis vulnerability → drug reversal screening**

Meldonium 目前应降为候选库中的一个待验证药物，而不是先验主角。

完整结果文件：

- [module effects](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase3_module_effects_primary.csv)
- [module decomposition summary](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase3_module_decomposition_summary_cell_lines.csv)
- [gene-level decomposition](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase3_gene_decomposition_summary_cell_lines.csv)
- [model-level EMT/pyrimidine/UPR state](C:/Users/21634/Documents/Codex/2026-08-16/meldonium-crc/outputs/phase3_emt_pyrimidine_upr_state_by_model.csv)

本分析是探索性、转录层面的证据。它用于决定下一轮筛选空间，不等于证明任何药物或节点的必要性。
