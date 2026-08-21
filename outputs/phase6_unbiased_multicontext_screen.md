# Phase 6：无偏 approved non-oncology perturbation screen

## 结论

第一轮无偏反筛没有发现可以直接定义为 Drug X 的候选。可解释药物面板中共有 28 个药物、501 个候选签名记录，最终成功评分 442 个签名，覆盖 HT29、HCC515、A549、MCF7 和 PC3 五类背景。

最重要的结果不是某个药排名第一，而是：

> **在多背景稳健性和 subtype selectivity 同时约束下，ERAD/proteostasis subtype 没有出现 RRS ≥0.60 且 selectivity ≥0.10 的非肿瘤药物。**

因此目前不能把任意一个“最高分药物”直接推进湿实验。

## 评分结果

| Drug | Multi-context RRS | ERAD-subtype RRS | ERAD selectivity | Contexts | Confidence |
|---|---:|---:|---:|---:|---|
| Atorvastatin | 0.535 | 0.530 | −0.005 | 5 | supported |
| Sulfasalazine | 0.527 | 0.522 | −0.005 | 4 | supported |
| Topiramate | 0.520 | 0.521 | +0.001 | 4 | supported |
| Amlodipine | 0.520 | 0.520 | −0.001 | 5 | supported |
| Teriflunomide | 0.510 | 0.518 | +0.007 | 4 | supported |
| Bortezomib | 0.515 | 0.508 | −0.008 | 5 | supported |

Atorvastatin 只是本轮 ERAD-subtype 分数最高的充分覆盖候选，但其 selectivity 是负值，不能解释为 ERAD-specific reversal。Minocycline 的 DHODH-subtype RRS 为 0.593，但只有 2 个签名/2 个背景，属于 insufficient-contexts，不能升级。

Leflunomide 与 Teriflunomide 在多背景数据中仍然没有形成优势：Teriflunomide 的 ERAD RRS 为 0.518，DHODH RRS 为 0.462；Leflunomide 的 ERAD RRS 为 0.487。

## Cross-context consistency

同一药物在不同细胞背景的 RRS 离散度已经被单独计算。大多数充分覆盖药物的 context-level SD 约为 0.02–0.05；Ivermectin、Allopurinol、Disulfiram、Celecoxib、Tacrolimus 等背景差异更大，只能作为 exploratory candidates。

这一步排除了“一个细胞背景高分、其他背景接近零”的伪阳性，但也提示：当前扰动 signature 与 OXA-R signature 的直接反向匹配整体偏弱，不能把 0.52–0.54 当作强药效预测。

## 首轮药物排序的解释

- Atorvastatin：脂质代谢/心血管药物背景；多背景稳定，但没有 ERAD selectivity。
- Sulfasalazine：炎症性肠病/免疫炎症背景；分数中等，没有 subtype enrichment。
- Topiramate：神经系统药物背景；ERAD selectivity 接近零。
- Amlodipine：心血管药物背景；多背景稳定，但没有明显反转优势。
- Teriflunomide：保留为 DHODH comparator，不是主候选。
- Bortezomib：保留为 proteostasis positive control，不是 novelty lead。

这些适应症分类是首轮候选注释，不应替代正式 FDA/ChEMBL indication audit。

## Novelty quick audit

本轮最高的充分覆盖候选也不能直接包装成“CRC 新药”：

- Atorvastatin 已有 CRC chemoresistance/高糖环境研究，并且已有 statin–chemotherapy 相关工作；[代表性研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC11887183/)、[CRC cholesterol vulnerability](https://pmc.ncbi.nlm.nih.gov/articles/PMC10311014/)。
- Sulfasalazine 已有 CRC stemness/metastasis 研究，且存在与 5-FU/oxaliplatin 联合的临床试验登记：[NCT06134388](https://clinicaltrials.gov/study/NCT06134388)。
- 因此本轮的“高分”候选仍然缺乏“CRC 基本空白 + OXA-R subtype selective”的组合条件。

## 为什么仍然不能结束 Phase 6

本轮面板仍然是“可解释药名优先”的第一轮，而不是整个 LINCS 匿名 BRD 化合物宇宙。HT29 unrestricted metadata page 中大量记录是 BRD 化合物 ID，不能直接标记为 approved drug；本轮选择先保证药物身份可追溯，再扩大到 ChEMBL/FDA 确认后的完整上市药物集合。

此外，L1000 signature 不是 OXA-resistant cell 的 signature。它只能回答药物扰动是否像 OXA-R 状态的反向改变，不能直接证明 OXA-R 细胞药敏增加。下一步必须加入 PRISM、DepMap、GDSC 或 CTRP 的药敏/依赖数据。

## 当前决策

1. 不把 Atorvastatin、Sulfasalazine、Topiramate 或 Amlodipine 直接推进湿实验；它们只是首轮可解释候选。
2. Bortezomib 继续作为 proteostasis positive control。
3. DHODH inhibitors 继续作为 pharmacological comparator。
4. 下一轮优先做两个增强：
   - 用 FDA/ChEMBL 核验后的完整 approved non-oncology universe，避免 curated list 偏倚；
   - 将 RRS 与独立药敏数据合并，要求“转录反转 + OXA-R selective sensitivity”同时成立。

## 输出文件

- [drug_subtype_rrs_ranked.csv](drug_subtype_rrs_ranked.csv)
- [per_signature_subtype_rrs.csv](per_signature_subtype_rrs.csv)
- [candidate_signature_metadata_ht29.csv](candidate_signature_metadata_ht29.csv)
- [Phase 5 结果报告](phase5_perturbation_reversal_score.md)
