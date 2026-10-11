# SCD / PID1 单细胞亚群定位分析

## 输入与标签来源

- 只读取本项目 `work/single_cell_GSE132465_GSE200997/` 中的两套 GEO 原始 UMI 矩阵和配套注释；未读取历史结果文件。
- 表达量按每细胞 UMI 总量归一到 10,000 后 log1p。
- GSE132465 使用 GEO 原始 `Cell_type` / `Cell_subtype` 注释。上皮以 CMS1–4 分层；成纤维细胞限于 stromal 中的 Myofibroblasts、Stromal 1、2、3。原注释的内皮、周细胞和平滑肌不当作成纤维细胞。
- GSE200997 没有细胞类型字段。按上皮/成纤维等经典 marker 均值推断宽泛谱系，仅保留 marker 判定可信的肿瘤细胞；上皮 CMS 使用 GEO 的逐细胞 `prediction` 字段。成纤维细胞按 myCAF、iCAF、ECM_CAF、apCAF 预设 marker 模块中最高的基因标准分数赋予描述性状态。该状态是模块推断，不等同于作者验证的细胞标签或独立聚类。
- 脂质代谢分数为预先列出的脂质合成、去饱和、延长、转运、氧化与储存基因的跨细胞标准分数均值，并排除 SCD 本身。

## 统计单位

细胞级结果用于描述；亚群间检验先按样本汇总每个亚群的均值，再对同一肿瘤样本内的亚群做配对 Wilcoxon 检验。每个亚群/样本至少 10 个细胞，且至少 3 个配对样本才检验；P 值在同一队列、区室、指标的亚群两两比较内做 BH 校正。

## 结果文件

- `subpopulation_expression_summary.tsv`：每亚群细胞数、覆盖样本数、SCD/PID1 表达和模块分数。
- `sample_subpopulation_means.tsv`：用于复核独立样本层面的均值。
- `paired_subtype_tests.tsv`：样本配对检验。
- `cell_subpopulation_expression.tsv.gz`：纳入分析的细胞及其亚群、目标表达和模块分数。
- `analysis_metadata.json`：原始输入、矩阵行数、面板基因覆盖情况。

## 解读限制

GSE132465 的原始亚群标签与 GSE200997 的推断标签来源不同，因此跨队列一致性仅作为方向性支持。GSE200997 的 CAF 状态标签是预设 marker 模块评分，不应写成无监督聚类结果。低细胞数亚群的表达与检验需谨慎解释。
