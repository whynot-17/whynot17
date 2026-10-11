# 11 个主分析候选基因的单细胞定位

## 输入与标签来源

- 分析基因：ABHD6, ACAT1, BMP2, FABP4, FHL1, INHBA, PID1, SCD, SCP2, TIMP1, UBE2C。这些是主分析中通过 BH-FDR 的 11 个候选基因。
- 只读取本项目 `work/single_cell_GSE132465_GSE200997/` 中的两套 GEO 原始 UMI 矩阵和配套注释；未读取历史结果文件。
- 表达量按每细胞 UMI 总量归一到 10,000 后 log1p。
- GSE132465 使用 GEO 原始 `Cell_type` / `Cell_subtype` 注释。上皮以 CMS1–4 分层；成纤维细胞限于 stromal 中的 Myofibroblasts、Stromal 1、2、3。原注释的内皮、周细胞和平滑肌不当作成纤维细胞。
- GSE200997 没有细胞类型字段。按上皮/成纤维等经典 marker 均值推断宽泛谱系；上皮 CMS 使用 GEO 的逐细胞 `prediction` 字段。成纤维细胞按 myCAF、iCAF、ECM_CAF、apCAF 预设 marker 模块中最高的基因标准分数赋予描述性状态。该状态是模块推断，不等同于作者验证的细胞标签或独立聚类。
- 脂质代谢分数为预先列出的脂质合成、去饱和、延长、转运、氧化与储存基因的跨细胞标准分数均值，并排除 SCD 本身。

## 统计单位

单细胞表达仅作定位与描述；配对亚群比较以样本为统计单位，要求每个样本/亚群至少 10 个细胞、至少 3 个配对样本。两两比较的 P 值在每个队列、区室、指标内做 BH 校正。

## 结果文件

- `tumor_lineage_expression_summary.tsv`：所有肿瘤细胞宽谱系的 11 基因表达及检测率。
- `subpopulation_expression_summary.tsv`：上皮 CMS / CAF 亚群的表达及检测率。
- `top_population_by_gene.tsv`：每基因在宽谱系、上皮亚群、成纤维亚群中表达均值最高的群体。
- `sample_subpopulation_means.tsv`、`tumor_lineage_sample_means.tsv`：样本层面均值。
- `paired_subtype_tests.tsv`：11 个候选基因及模块的样本配对比较。
- `tumor_cell_gene_expression.tsv.gz`：所有肿瘤细胞的细胞类型标签和 11 基因表达。
- `cell_subpopulation_expression.tsv.gz`：目标上皮/成纤维细胞的细胞级表达、亚群与模块分数。
- `primary11_tumor_lineage_dotplot.svg`：跨队列宽谱系点图；点大小代表检测率，颜色代表该基因在细胞群间的标准化均值。
- `analysis_metadata.json`：原始输入、矩阵行数、面板基因覆盖情况。

## 解读限制

GSE132465 的原始亚群标签与 GSE200997 的推断标签来源不同，因此跨队列一致性仅作为方向性支持。GSE200997 的 CAF 状态标签是预设 marker 模块评分，不应写成无监督聚类结果。低细胞数亚群的表达与检验需谨慎解释。


## 本轮结果摘要

- 两套矩阵均检测到全部 11 个候选基因；肿瘤细胞数：GSE132465 47,285，GSE200997 31,586。
- 两队列宽谱系最高表达均指向基质/成纤维细胞：BMP2, FABP4, FHL1, INHBA, TIMP1。
- 两队列宽谱系最高表达均指向上皮细胞：SCD, UBE2C。
- 两队列宽谱系最高表达均指向髓系细胞：PID1。
- 宽谱系最高均值不等同于疾病特异性或统计学显著；群体表达均值、检测率和细胞/样本数请结合 TSV 解读。

### 配对亚群比较通过 BH-FDR 的结果

以下为样本配对的亚群均值比较；A−B 为正表示 A 亚群表达更高。FDR 在每队列、区室、指标内校正。

| 队列 | 区室 | 基因 | 亚群 A | 亚群 B | 配对样本数 | A−B 中位数 | BH-FDR |
|---|---|---|---|---|---:|---:|---:|
| GSE132465 | Tumor_epithelial | ACAT1 | CMS2 | CMS3 | 14 | 0.0910 | 0.001099 |
| GSE132465 | Tumor_epithelial | UBE2C | CMS2 | CMS3 | 14 | 0.1410 | 0.001099 |
| GSE132465 | Tumor_epithelial | SCP2 | CMS1 | CMS3 | 7 | -0.1972 | 0.02344 |
| GSE132465 | Tumor_epithelial | SCP2 | CMS2 | CMS3 | 14 | -0.0861 | 0.02344 |
| GSE132465 | Tumor_epithelial | UBE2C | CMS1 | CMS3 | 7 | 0.1942 | 0.02344 |
| GSE132465 | Tumor_epithelial | BMP2 | CMS2 | CMS3 | 14 | -0.0282 | 0.02563 |
