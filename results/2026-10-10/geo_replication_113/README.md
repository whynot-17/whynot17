# 113 个 DINP–CRC 候选基因在四个 GEO 队列中的复现

- 分析日期：2026-10-10
- 目标基因：来自本项目已保存的 113 基因交集明细。
- 复现阈值：每队列以全平台唯一映射基因做 BH 校正；FDR < 0.05 且 |log2FC| > 1 记为该队列复现。
- 稳定复现定义：至少 3/4 队列达到上述阈值，且达到阈值的队列方向一致。
- 独立来源敏感性分析：GSE21510 与 GSE32323 的样本标题含有相同患者编号；因此保守地将二者合并为一个来源组，再加 GSE44076、GSE41258，共 3 个来源组。独立来源复现要求三个来源组均至少有一个队列达到显著性与效应量阈值，且同向。
- 差异方向：log2FC = 肿瘤 − 正常；正值为肿瘤上调，负值为肿瘤下调。
- 探针处理：各平台按 NCBI GEO 平台注释将探针映射到 gene symbol；多探针基因按样本取中位数；未做别名映射。
- GEO 系列矩阵已归一化；GSE41258 的 PLIER 信号按 log2(intensity) 转换，其余三个矩阵保留提交者提供的 log 表达尺度。
- GSE21510 的 LCM 肿瘤与全组织样本分开归一化，且没有 LCM 正常对照；主比较仅用 homogenized 肿瘤/正常。
- GSE44076 的 50 个健康黏膜样本不纳入主肿瘤对癌旁分析；仅使用 98 个肿瘤与 98 个配对癌旁样本。
- GSE41258 限定 included in analysis=Yes 的 Primary Tumor 与 Normal Colon；重复/配对患者按 patient ID 聚类计算稳健标准误。

## 队列样本数

| GEO 队列 | 矩阵总样本 | 纳入肿瘤 | 纳入正常 | 可配对患者 |
|---|---:|---:|---:|---:|
| GSE21510 | 148 | 19 | 25 | 1 |
| GSE32323 | 44 | 17 | 17 | 17 |
| GSE44076 | 246 | 98 | 98 | 98 |
| GSE41258 | 390 | 182 | 53 | 44 |

## 稳定复现结果

- 稳定复现基因数：**23 / 113**。
- 肿瘤上调方向：6 个；肿瘤下调方向：17 个。
- 基因：BMP2, INHBA, MUC2, NR3C1, PCK1, PLAC8, SCD, TIMP1, UBE2C, AOC3, CPT2, DDIT4, DES, ELOVL6, FABP1, FHL1, IL18, RARRES2, RETSAT, RGS2, THRB, TPSAB1, ZP3。
- 按重叠患者合并为 3 个来源组后，仍在全部 3 组同向复现：**15 / 113**。
- 三来源组复现基因：BMP2, INHBA, MUC2, NR3C1, PCK1, PLAC8, SCD, TIMP1, UBE2C, AOC3, DES, FABP1, FHL1, RGS2, TPSAB1。
- GSE21510 与 GSE32323 共享患者编号：肿瘤 17 个、正常 1 个；逐例清单见 GEO_GSE21510_GSE32323_patient_id_overlap.tsv。

## 输出文件

- GEO_cohort_sample_counts.tsv：各队列样本纳入数与分析方法。
- GEO_platform_gene_coverage.tsv：113 个基因的平台探针覆盖率及测试基因数。
- 113_gene_per_cohort_DE_results.tsv：113 基因在每个队列的 log2FC、P、FDR 与方向。
- 113_gene_stability_summary.tsv：每个目标基因的跨队列复现次数与方向。
- stable_genes_replicated_in_3_or_more_cohorts.tsv：达到预设稳定标准的基因。
- stable_genes_replicated_in_3_source_groups.tsv：按重叠患者保守合并后，在三个来源组均复现的基因。
- GEO_GSE21510_GSE32323_patient_id_overlap.tsv：两个 GEO 系列之间共享患者编号的样本清单。
- source_checksums.json：下载矩阵与注释来源 SHA-256。

## 复现边界

平台注释的基因 symbol 来自 GEO/GPL；本分析按 symbol 精确映射，不做当前 HGNC 别名回溯。队列间采用与样本设计相符的配对检验或按患者聚类的稳健检验；候选基因是否进入稳定列表还同时要求效应量阈值和同向复现。GSE21510 与 GSE32323 患者编号重叠提示其并非完全独立验证；因此同时给出按四个 GEO 系列统计的结果和将这两个系列合并后的保守来源组结果。
