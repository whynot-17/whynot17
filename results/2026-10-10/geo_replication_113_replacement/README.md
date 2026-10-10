# 113 个 DINP×CRC 交集基因在替换后的四个 GEO 队列中的复现

- 分析日期：2026-10-10
- 查询集：本项目 113 个 DINP×TCGA-GDC FDR DEG 交集基因；表达矩阵及平台注释从 NCBI GEO 重新下载。
- 队列：GSE32323、GSE44076、GSE41258、GSE106582。
- 每队列复现阈值沿用此前口径：BH-FDR < 0.05 且 |log2FC| > 1；log2FC 为 Tumor−Normal。
- 主复现集合要求至少 3/4 队列达到阈值，且显著队列方向一致；严格交集另按 4/4 列出。
- 多探针以 gene symbol 聚合取表达中位数；BH-FDR 在该平台全部可映射基因上计算。
- GSE32323、GSE44076、GSE106582 使用完整配对的 paired t-test。GSE41258 有患者重复切片，使用患者聚类稳健标准误；44 个完整配对患者和其他符合纳入条件的样本均进入模型。

## 队列样本与配对数

| 队列 | 平台 | 肿瘤样本 | 正常/黏膜样本 | 完整配对 | 实际使用配对 |
|---|---|---:|---:|---:|---:|
| GSE32323 | GPL570 | 17 | 17 | 17 | 17 |
| GSE44076 | GPL13667 | 98 | 98 | 98 | 98 |
| GSE41258 | GPL96 | 182 | 53 | 44 | 44 |
| GSE106582 | GPL10558 | 77 | 117 | 68 | 68 |

- 至少 3/4 队列同向复现：**19 个基因**。
- 严格 4/4 同向复现：**9 个基因**。

### 至少 3/4 队列复现

| 基因 | 队列数 | 达标队列 | 方向 |
|---|---:|---|---|
| AOC3 | 3 | GSE32323;GSE41258;GSE44076 | Tumor_down |
| BMP2 | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_down |
| DDIT4 | 3 | GSE106582;GSE32323;GSE44076 | Tumor_up |
| DES | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_down |
| FABP4 | 3 | GSE106582;GSE41258;GSE44076 | Tumor_down |
| FHL1 | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_down |
| INHBA | 3 | GSE32323;GSE41258;GSE44076 | Tumor_up |
| MUC2 | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_down |
| NR3C1 | 3 | GSE32323;GSE41258;GSE44076 | Tumor_down |
| PCK1 | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_down |
| PDK4 | 3 | GSE106582;GSE41258;GSE44076 | Tumor_down |
| PID1 | 3 | GSE106582;GSE41258;GSE44076 | Tumor_down |
| PLAC8 | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_down |
| RETSAT | 3 | GSE106582;GSE32323;GSE44076 | Tumor_down |
| RGS2 | 3 | GSE32323;GSE41258;GSE44076 | Tumor_down |
| SCD | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_up |
| TIMP1 | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_up |
| TPSAB1 | 3 | GSE32323;GSE41258;GSE44076 | Tumor_down |
| UBE2C | 4 | GSE106582;GSE32323;GSE41258;GSE44076 | Tumor_up |

### 严格 4/4 交集

| 基因 | 方向 |
|---|---|
| BMP2 | Tumor_down |
| DES | Tumor_down |
| FHL1 | Tumor_down |
| MUC2 | Tumor_down |
| PCK1 | Tumor_down |
| PLAC8 | Tumor_down |
| SCD | Tumor_up |
| TIMP1 | Tumor_up |
| UBE2C | Tumor_up |

## 数据来源与结果文件

- [GSE32323](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32323)
- [GSE44076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44076)
- [GSE41258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41258)
- [GSE106582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106582)
- 113_gene_per_cohort_DE_results.tsv：各基因在队列中的效应、P 值、BH-FDR 和方向。
- 113_gene_stability_summary.tsv：逐基因复现情况汇总。
- stable_genes_replicated_in_3_or_more_of_4_cohorts.tsv 和 stable_genes_replicated_in_all_4_cohorts.tsv：3/4 及 4/4 交集。
- GEO_cohort_sample_counts.tsv、GEO_sample_selection_audit.tsv：队列样本数、配对数和逐样本筛选记录。
- GEO_platform_gene_coverage.tsv：候选基因在各平台的探针覆盖。
- analysis_metadata.json、source_checksums.json：分析参数、来源 URL 和输入文件 SHA-256。
- 可复现脚本：analysis/validate_113_across_geo_replacement.py。

## 解释边界

配对检验仅使用完整肿瘤/正常患者配对；GSE41258 的聚类稳健模型将患者作为聚类单位。不同平台的探针注释和信号尺度可能影响效应量；跨队列复现不代表因果关系或临床诊断性能。
