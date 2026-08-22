# MCOP–CRC Phase 2E：DINP-axis molecular validation

## 结论先行

本轮不是重新筛化学物，也不是把 MCOP 事后包装成 CTD 的直接发现物。冻结的逻辑是：

> **CTD 发现 DINP/MiNP 轴的分子桥接 → NHANES 以 MCOP 作为 DINP 轴尿液 biomarker 出现 CRC association → 机制验证先审计桥接证据的可重复性。**

- **MiNP**：24 个 human interacting genes，5 个 GeneCards Disorders CRC overlap；Phase 1 primary OR=10.1，BH-FDR=0.00346。
- **MCOP**：19 个 human interacting genes，2 个 overlap；Phase 1 primary BH-FDR=0.284，因此不能写成 CTD 已直接发现 MCOP–CRC 分子桥。
- MiNP 的 5 个 CRC overlap gene 中，4 个有 co-treatment evidence，1 个有 single-chemical statement；共处理依赖性是本桥接目前最大的限制。
- 本轮本地 Hallmark/Reactome ORA 使用 **U_core=22,786 个 CTD 可检测 human genes** 作背景。全局 BH-FDR<0.05 且 overlap≥2 的 exploratory rows：**293**（142 个不同 term；MiNP query 为 74 rows/51 terms）；Reactome 层级 term 高度冗余，因此不把这个计数当作独立机制数量。

## 1. 分析边界与输入冻结

| 元素 | 冻结定义 |
|---|---|
| CTD interaction | Homo sapiens；统计按 unique ChemicalID × GeneID；PMID 只用于证据审计 |
| CRC gene set | GeneCards Disorders-scoped export；与 U_core 相交后 585 genes；Phase 1 primary row 为 k=1000 |
| U_core | 最终 core environmental toxicants 的所有 CTD human interacting genes；本轮复核为 22,786 genes |
| 分子比较 | MiNP、MCOP、DINP parent、MBzP；不重新筛选候选 |
| pathway | local MSigDB Hallmark 2026.1 与 Reactome 2026.1 GMT；ORA exploratory；不使用全基因组背景 |

## 2. MiNP 与 MCOP 的桥接对照

| 化学物 | 角色 | CTD human genes | CRC overlap | unique PMIDs | co-treatment genes | single-chemical genes | Phase 1 BH-FDR |
|---|---|---:|---:|---:|---:|---:|---:|
| MiNP | DINP-axis molecular discovery | 24 | 5 | 4 | 14 | 9 | 0.00346 |
| MCOP | NHANES urinary biomarker | 19 | 2 | 2 | 0 | 18 | 0.284 |
| DINP_parent | parent DINP comparator | 86 | 4 | 16 | 56 | 30 | 0.449 |
| MBzP | historical human-validation comparator | 87 | 18 | 20 | 31 | 32 | 8.99e-10 |

### MiNP CRC overlap genes

`BAX;CASP8;CDKN2A;MIR141;PPARG`

逐基因证据见 `mcop_phase2e_molecular_bridge_qc.csv`。关键限制是：BAX、CASP8、MIR141、CDKN2A 的 MiNP CTD 记录来自同一个 phthalate co-treatment study/PMID，而 PPARG 才有单化学物的 binding/activity 记录。因此，这个桥接支持“DINP/MiNP 轴存在 CRC-relevant molecular evidence”，但目前不支持“MiNP 单独通过 5 个 CRC genes 直接驱动该桥接”。

## 3. Co-treatment sensitivity

去掉所有带 co-treatment 标记的 CTD rows 后，重新使用同一个 U_core 与 GeneCards CRC background 计算的是**未进行跨化学物重排的敏感性 P 值**，不是新的 confirmatory FDR。

| 化学物 | no-co-treatment genes | CRC overlap | OR | Fisher P（未调整） | overlap genes |
|---|---:|---:|---:|---:|---|
| MiNP | 10 | 1 | 4.22 | 0.229 | PPARG |
| MCOP | 19 | 2 | 4.48 | 0.0844 | CHASERR;NEAT1 |
| DINP_parent | 30 | 3 | 4.23 | 0.0409 | PPARG;RELA;STAT3 |
| MBzP | 59 | 14 | 12.1 | 2.13e-10 | CCND1;CDKN2A;MALAT1;MIR126;MIR127;MIR143;MIR221;MIR25;MIR320A;MIR328;MIR451A;PPARG;TP53;ZFAS1 |

MiNP 的结果在去除共处理证据后明显变弱，说明 Phase 1 的强富集不能被解释为一个完全独立、纯 MiNP 单化学物的 CRC gene signature。它仍可作为轴级别机制假设，但需要单化学物实验或独立暴露组学来确认。

## 4. Pathway ORA

本轮 ORA 只用于定位可检验的方向，不用于把 24 个 CTD genes 变成机制定论。对于 MiNP，4/5 CRC overlap gene 受共处理证据影响，且 gene list 很短；因此 pathway 结果必须按 exploratory 处理。

结果文件：`mcop_phase2e_pathway_ora.csv`。所有可测试通路（包括 0 overlap）均纳入 BH；全局 BH-FDR 作为跨 query 与 Hallmark+Reactome 的保守参考。

MiNP 的领先方向（只作候选机制）：

最稳定的方向是核受体/PPAR 相关转录调控；这个方向在去掉 co-treatment 后仍由 PPARA、PPARG、PPARD、NR1H2/NR1H3、NR1I2/NR1I3 等基因驱动，但它是 **DINP-axis level** 的候选机制，不等于 CRC-specific pathway 已被验证。Reactome 的 SUMOylation、脂质代谢和凋亡 terms 之间存在明显基因集重叠。
- `REACTOME_NUCLEAR_RECEPTOR_TRANSCRIPTION_PATHWAY`：query=MiNP_no_cotreatment_genes，overlap=7，genes=NR1H2;NR1H3;NR1I2;NR1I3;PPARA;PPARD;PPARG，global BH-FDR=4.38e-13。
- `REACTOME_NUCLEAR_RECEPTOR_TRANSCRIPTION_PATHWAY`：query=MiNP_all_human_interacting_genes，overlap=8，genes=ESR1;NR1H2;NR1H3;NR1I2;NR1I3;PPARA;PPARD;PPARG，global BH-FDR=2.65e-12。
- `REACTOME_SUMOYLATION`：query=MiNP_all_human_interacting_genes，overlap=9，genes=CDKN2A;DNMT1;ESR1;NCOA1;NR1H2;NR1H3;NR1I2;PPARA;PPARG，global BH-FDR=5.17e-10。
- `REACTOME_SUMOYLATION_OF_INTRACELLULAR_RECEPTORS`：query=MiNP_all_human_interacting_genes，overlap=6，genes=ESR1;NR1H2;NR1H3;NR1I2;PPARA;PPARG，global BH-FDR=8.71e-10。
- `REACTOME_SUMOYLATION_OF_INTRACELLULAR_RECEPTORS`：query=MiNP_no_cotreatment_genes，overlap=5，genes=NR1H2;NR1H3;NR1I2;PPARA;PPARG，global BH-FDR=1.32e-09。
- `REACTOME_SUMOYLATION`：query=MiNP_no_cotreatment_genes，overlap=6，genes=NCOA1;NR1H2;NR1H3;NR1I2;PPARA;PPARG，global BH-FDR=6.4e-08。
- `REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA`：query=MiNP_no_cotreatment_genes，overlap=5，genes=NCOA1;NR1H2;NR1H3;PPARA;PPARG，global BH-FDR=7.06e-07。
- `REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS`：query=MiNP_all_human_interacting_genes，overlap=5，genes=BAX;BCL2;CASP3;CASP8;CDKN2A，global BH-FDR=2e-06。

## 5. 当前机制判定

### 支持

- MiNP 的 Phase 1 富集在冻结的 CTD×GeneCards 流程中可复核；
- MiNP overlap 中包含 PPARG 的 single-chemical binding/activity evidence；
- MCOP 人群信号可以合理描述为 DINP-axis biomarker signal，而不是 MCOP 已被 CTD 分子桥接直接发现；
- 机制验证的下一步应围绕 PPAR/核受体、氧化应激/凋亡和表观遗传调控提出可检验假设，而不是继续扩大无先验通路清单。

### 限制

- MiNP CRC overlap 的主要证据集中在共处理研究；
- MCOP 本身 CTD overlap 很弱（2 genes，FDR≈0.284），不能作为 MCOP-specific molecular proof；
- ORA 的输入列表很短且来自数据库交集，显著性不等于暴露因果；
- NHANES 仍是横断面人群发现，WHI 尚未获得或分析真实 biospecimen。

## 文件

- `mcop_phase2e_molecular_candidate_summary.csv`：四个化学物的桥接与证据汇总
- `mcop_phase2e_molecular_bridge_qc.csv`：CRC overlap gene-level evidence QC
- `mcop_phase2e_molecular_evidence_long.csv`：CTD interaction 原始证据长表
- `mcop_phase2e_pathway_ora.csv`：custom-background Hallmark/Reactome exploratory ORA
- `mcop_phase2e_figure_bridge_evidence.png`：桥接大小与证据构成图

## 可复现性

- 运行时间（UTC）：2026-08-22T13:09:22.611663+00:00
- 脚本：`work\scripts\mcop_phase2e_molecular_validation.py`
- CTD chemicals SHA256：`9e4b642c8716140d30a9376d6b2229acb81ee48a8b106a4f41ed06b29894ef6c`
- CTD interactions SHA256：`05e1b0d2d93bb33f72659e6ed2d590304cdb58da9538c787ad51fc6624d0e055`
- GeneCards Disorders SHA256：`79bff2923c36b02f8f7307a36604890cf3de73e2e17488ccb6f60aa392a87362`
- Phase 1 ranked output SHA256：`609239674c641a307db49bfbf7392a392dcfea6ea50b6bda0d477ed2266a2ce3`

**最终判断：Phase 2E 通过的是“继续做 DINP-axis molecular validation”的门，不是“MCOP 已有一对一 CTD 机制证明”的门。**
