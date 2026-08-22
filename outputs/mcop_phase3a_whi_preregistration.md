# MCOP–CRC Phase 3A：WHI 外部复制预注册冻结版

## Material Passport

- Origin Skill: academic-research-suite + experiment-agent
- Origin Mode: plan
- Origin Date: 2026-08-22
- Verification Status: UNVERIFIED — WHI biospecimen/data access pending
- Version Label: `mcop_whi_phase3a_prereg_v1.0`

## 研究问题与假设

在 WHI BMD subcohort 中，诊断前尿 MCOP 浓度是否与随后发生的 incident invasive colorectal cancer 相关？

方向性假设：MCOP 每翻倍对应 CRC odds 增加，主效应 OR > 1。

## 设计与研究对象

- Design：WHI prospective nested case-control replication study。
- Cases：随访期间首次发生的 adjudicated invasive CRC；排除 baseline 已有 CRC 或其他影响入组资格的癌症记录。
- Controls：由 WHI 提供的 cancer-free risk-set matched controls；每个 matched set 恰好 1 个 case，固定为 1:1 或 1:2 两种预先定义设计之一。
- Exposure timing：尿样采集必须早于 CRC 诊断；记录 sample-to-diagnosis interval。
- Primary population：所有符合上述条件且 MCOP、肌酐和协变量可用者。

## 暴露定义

- Primary exposure：`log2(MCOP_ng_mL)`；效应解释为 MCOP 每翻倍 OR。
- MCOP 命名兼容：MCOP、MCiOP、mCIOP、mono(carboxyisooctyl) phthalate。
- Non-detect：使用 assay-specific LLOQ/LOD 在实验室盲态下冻结；主分析采用 `LLOQ / sqrt(2)`，同时报告检出率。
- Urinary dilution：主模型加入 `log2(creatinine_mg_dL)`；不把 creatinine-normalized MCOP 作为主暴露。
- Repeated urine：若 baseline 与 year-3 均可用，次要分析使用两个 `log2(MCOP)` 的算术平均；不得在看到 CRC 结果后选择时间点。

## 结局与协变量

- Outcome：incident invasive CRC，整体作为主结局。
- Secondary outcomes：colon-only、rectum-only（若病例数允许，仅作描述性/方向性分析）。
- Prespecified covariates：age、BMI、smoking、alcohol、physical activity、SES、sex、race、urinary creatinine、assay batch。
- Continuous covariates：在分析数据中 z-score；MCOP 保持原始 log2 单位，因此主 OR 不改变。
- Assay batch：作为分类固定效应；若 batch 与 matched set 完全共线，则记录为 non-estimable，不进行事后替代性调节。
- 条件 logistic 只能估计 matched set 内有变异的协变量；任何在所有 matched sets 内均无变异的预设协变量都会被脚本自动记录并标记为 non-estimable，不进行事后替代性调节。

## 主模型与敏感性分析

主模型固定为 conditional logistic regression，以 matched set 为 stratum：

`CRC ~ log2(MCOP) + age + BMI + smoking + alcohol + physical_activity + SES + sex + race + log2(creatinine) + assay_batch`

预先定义：

1. Primary：所有合格的 prediagnostic samples；
2. Lag ≥2 years：sample-to-diagnosis interval ≥2 年；
3. Lag ≥5 years：sample-to-diagnosis interval ≥5 年；
4. MCOP quartiles：仅作为非线性/剂量结构描述，不替代连续主模型；
5. Repeated-urine mean：baseline/year-3 log2 MCOP 平均值；
6. 结果按 case 数、matched-set 数、检测率、缺失率和 assay batch 报告。

## 统计判定

- Primary alpha：0.05，two-sided Wald test；
- 报告 OR、95% CI、P、matched-set 数、病例数和 complete-case N；
- 不以单一 P<0.05 作为外部复制唯一标准；预先关注方向、效应量、CI 是否与中等正关联相容，以及 lag/repeated urine 是否同向；
- 不进行新的化学物筛选、TCGA/PPI/GO、全基因组机制扩展或结果驱动的多模型挑选。

## Power 情景

在 WHI 样本数确认前，使用 NHANES Phase 2 complete-case MCOP/covariate 分布作 planning proxy，模拟 CRC cases = 181、150、120、100、80；controls/case = 1、2；目标 OR = 1.15、1.20、1.25、1.30。每个情景采用 1000 次 conditional-logistic simulation，alpha=0.05，two-sided。

这些数字仅用于申请前规划，不能替代 WHI 实际 MCOP 方差、缺失率、尿样剩余量和 assay performance。

## Access gate

真实 WHI 分析开始前必须确认：

1. CRC cases 中有可用 prediagnostic urine 的人数；
2. 每个时间点的剩余体积与冻融/消耗记录；
3. 可匹配 controls 数量及 matched-set 结构；
4. AS458/AS498 等既往 ancillary study 的样本消耗和 assay overlap；
5. MCOP assay 的 LLOQ、批次、盲法 QC 和检测成功率。

本文件在 WHI access 确认前冻结；任何改变主暴露、主结局、匹配设计或主要协变量的修改都必须产生新版本并记录理由。
