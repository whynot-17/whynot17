# 16个基因单细胞队列的样本级配对检验

## 统计单位与配对

以样本/患者为独立统计单位；每个样本、细胞类型、基因先汇总为该细胞类型所有细胞的平均 `log1p(CP10K)`，随后进行双侧精确配对 Wilcoxon 符号秩检验。未把细胞当作独立重复。仅使用该患者该细胞类型在肿瘤和正常样本中均有细胞的完整配对；缺失细胞类型不当作零表达。

- GSE132465：从 GEO 注释中的 `Patient`、`Class`、`Sample` 字段识别10对患者样本。
- GSE200997：从注释中的 `B_cacN`/`T_cacN`样本编号识别7对正常/肿瘤患者样本；此对应也与GEO样本记录中的患者编号一致。
- 共检验 256 个可配对的基因×细胞类型×队列组合；全局 BH-FDR 在全部组合上校正。另附每队列内BH校正作为敏感性结果。

## 结果

名义 `p < 0.05` 的比较有 47 个；全局 BH 校正后 `q < 0.05` 有 13 个；分队列 BH 校正后 `q < 0.05` 有 21 个。

全局 BH 校正后显著的结果：

| 队列 | 基因 | 细胞类型 | 完整配对数 | 肿瘤−正常均值差（log1p CP10K） | 精确 p | 全局 BH q |
|---|---|---|---:|---:|---:|---:|
| GSE132465 | FABP4 | B cells | 10 | -0.0064 | 0.00195312 | 0.0384615 |
| GSE132465 | FHL1 | B cells | 10 | -0.0260 | 0.00195312 | 0.0384615 |
| GSE132465 | UBE2C | B cells | 10 | 0.0448 | 0.00195312 | 0.0384615 |
| GSE132465 | SCD | Epithelial cells | 10 | 0.4368 | 0.00195312 | 0.0384615 |
| GSE132465 | TIMP1 | Epithelial cells | 10 | 0.9304 | 0.00195312 | 0.0384615 |
| GSE132465 | ACAT1 | Stromal cells | 10 | -0.0775 | 0.00195312 | 0.0384615 |
| GSE132465 | FHL1 | Stromal cells | 10 | -0.5348 | 0.00195312 | 0.0384615 |
| GSE132465 | INHBA | Stromal cells | 10 | 0.7355 | 0.00195312 | 0.0384615 |
| GSE132465 | TIMP1 | Stromal cells | 10 | 1.0076 | 0.00195312 | 0.0384615 |
| GSE132465 | ABHD6 | T cells | 10 | 0.0080 | 0.00195312 | 0.0384615 |
| GSE132465 | FHL1 | T cells | 10 | -0.0163 | 0.00195312 | 0.0384615 |
| GSE132465 | INHBA | T cells | 10 | 0.0085 | 0.00195312 | 0.0384615 |
| GSE132465 | UBE2C | T cells | 10 | 0.0505 | 0.00195312 | 0.0384615 |

## 主分析：GSE132465

若将GSE132465作为主分析、GSE200997作为敏感性探索，GSE132465队列内BH校正后有 21 个显著的基因×细胞类型比较；去重后是 11 个不同基因（共16个算法基因）。

| 基因 | 显著细胞类型 | 显著比较数 | 最小 p | 最小队列内 BH q |
|---|---|---:|---:|---:|
| ABHD6 | T cells | 1 | 0.00195312 | 0.0144231 |
| ACAT1 | Stromal cells | 1 | 0.00195312 | 0.0144231 |
| BMP2 | Epithelial cells | 1 | 0.00390625 | 0.0220588 |
| FABP4 | B cells, Stromal cells | 2 | 0.00195312 | 0.0144231 |
| FHL1 | B cells, Stromal cells, T cells | 3 | 0.00195312 | 0.0144231 |
| INHBA | Stromal cells, T cells | 2 | 0.00195312 | 0.0144231 |
| PID1 | Stromal cells | 1 | 0.00976562 | 0.0446429 |
| SCD | Epithelial cells, T cells | 2 | 0.00195312 | 0.0144231 |
| SCP2 | T cells | 1 | 0.00976562 | 0.0446429 |
| TIMP1 | B cells, Epithelial cells, Stromal cells | 3 | 0.00195312 | 0.0144231 |
| UBE2C | B cells, Epithelial cells, Stromal cells, T cells | 4 | 0.00195312 | 0.0144231 |

队列内检验计数（数字是基因×细胞类型比较，不是唯一基因数）：

| 队列 | 名义 p<0.05 | 全局BH q<0.05 | 队列内BH q<0.05 |
|---|---:|---:|---:|
| GSE132465 | 34 | 13 | 21 |
| GSE200997 | 13 | 0 | 0 |

全局BH的13个显著比较去重后为 8 个基因，全部来自GSE132465。GSE200997作为敏感性探索，有13个名义 `p<0.05` 比较，但没有通过队列内或全局BH校正的比较；部分细胞类型完整配对数仅2–4，检验效能有限。

## 解读

统计显著表示样本层面表达差异与零差异不一致，不等同于DINP机制证据。GSE200997细胞类型为标志基因推断标签；两队列的细胞类型标签并非完全等价。效应值为平均 `log1p(CP10K)` 差值，适合描述方向和尺度，不应直接解释为倍数变化。

## 文件

- `sample_level_paired_tests.tsv`：所有组合的配对数、方向、精确p值与BH校正q值。
- `analysis/single_cell_sample_level_pvalues.py`：可复现分析脚本。

队列来源：[GSE132465 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)；[GSE200997 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997)。
## 主分析基因与跨队列单细胞一致信号的交集

按前面报告的跨队列同细胞区室、肿瘤/正常方向一致信号定义，第二组有 8 个唯一基因：`FHL1, INHBA, PCK1, PID1, SCD, TIMP1, TPM1, UBE2C`。它与主分析GSE132465中BH-FDR通过的11个唯一基因相交，得到 **6 个唯一基因**：`FHL1, INHBA, PID1, SCD, TIMP1, UBE2C`。TIMP1在上皮和基质/成纤维两个区室均有一致信号，但交集按基因去重只计1个。

| 基因 | 跨队列对应细胞类型 | 肿瘤−正常方向 | GSE132465主分析BH q | GSE200997方向差值（均值 log1p CP10K） |
|---|---|---|---:|---:|
| TIMP1 | Epithelial cells ↔ Epithelial | up | 0.01442 | 0.2959 |
| TIMP1 | Stromal cells ↔ Fibroblast | up | 0.01442 | 1.1773 |
| SCD | Epithelial cells ↔ Epithelial | up | 0.01442 | 0.0678 |
| UBE2C | Epithelial cells ↔ Epithelial | up | 0.02206 | 0.0715 |
| INHBA | Stromal cells ↔ Fibroblast | up | 0.01442 | 0.5624 |
| FHL1 | Stromal cells ↔ Fibroblast | down | 0.01442 | -0.2314 |
| PID1 | Stromal cells ↔ Fibroblast | down | 0.04464 | -0.2883 |

GSE200997作为敏感性探索时，上述方向一致仅表示方向/定位支持，不代表该队列通过自身BH-FDR。逐区室配对效应与两队列p/q值见 `primary11_scRNA_concordance_intersection.tsv`；复现交集的脚本为 `analysis/intersect_primary11_with_scRNA_concordance.py`。
