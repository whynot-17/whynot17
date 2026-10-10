# GSE41258 44 对患者配对分析与复现名单复核

- 分析日期：2026-10-10
- 本次只重算 GSE41258；其他三个队列沿用本项目上一轮的结果。
- GSE41258/GPL96 Series Matrix 和平台注释来自 NCBI GEO；按 included in analysis=Yes 选 Primary Tumor 与 Normal Colon。
- 对同一患者、同一组织的重复芯片先取表达均值，再对 44 个完整配对做 paired t-test；PLIER 信号按 log2 转换。
- 多 probe 按 gene symbol 取表达中位数；BH-FDR 在全部可映射基因上重新计算。阈值保持 FDR<0.05 且 |log2FC|>1。

## 原 GSE41258 模型的 29 个显著候选

- 新配对分析仍达到阈值：29/29。
- 同时达到阈值且方向与原模型一致：29/29。
- 原模型与配对模型的 29 个 log2FC 高度相关：Pearson r=0.995，Spearman ρ=0.977；log2FC 绝对差异中位数=0.081。
- 未达到阈值：无。
- 新配对模型达到阈值、但不在原 29 个中的候选：INHBB。

## 原 19 个 3/4 基因

- 原名单仍保留：19/19。
- 移出：无。
- 新增：无。
- 更新名单（19 个）：AOC3, BMP2, DDIT4, DES, FABP4, FHL1, INHBA, MUC2, NR3C1, PCK1, PDK4, PID1, PLAC8, RETSAT, RGS2, SCD, TIMP1, TPSAB1, UBE2C。

## 原 9 个 4/4 基因

- 仍满足 4/4：9/9。
- 不再满足：无。
- 更新 4/4 名单（9 个）：BMP2, DES, FHL1, MUC2, PCK1, PLAC8, SCD, TIMP1, UBE2C。

## 文件

- GSE41258_113_candidates_paired_DE.tsv：候选基因配对检验结果。
- original_model_29_candidates_rechecked.tsv：原 29 个候选逐项比较。
- original_19_status_after_paired_reanalysis.tsv、original_9_four_of_four_status.tsv：原名单复核。
- updated_stable_genes_3_of_4.tsv、updated_stable_genes_4_of_4.tsv：更新后的名单。
- updated_113_gene_per_cohort_results.tsv、updated_113_gene_stability_summary.tsv：替换单队列后完整汇总。
- GSE41258_all_mapped_genes_paired_DE.tsv：全平台基因的 P 值与 BH-FDR。
- GSE41258_44_pair_sample_counts.tsv：44 对的患者内重复芯片数。
- analysis_summary.json：关键数字和输入数据校验和。
- 可复现脚本：analysis/recheck_gse41258_44_paired.py。

数据记录：[NCBI GEO GSE41258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41258)。
