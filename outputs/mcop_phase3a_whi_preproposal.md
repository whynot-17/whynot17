# WHI MCOP–CRC 外部复制：Phase 3A 一页 pre-proposal

## 项目摘要

本项目拟在 Women’s Health Initiative（WHI）Bone Mineral Density subcohort 中开展 MCOP 与 incident colorectal cancer 的外部复制。研究目标是把 NHANES 中的横断面关联转化为“诊断前尿 MCOP → 随后 incident CRC”的前瞻性证据。

WHI AS498 已有 181 名 adjudicated CRC cases，样本平均在诊断前约 7.2 年采集；WHI AS458 已在相关 WHI 尿样中使用 CDC 方法测量包含 MCOP 在内的 13 种 urinary phthalate metabolites。因此，本研究可以直接复用既有 WHI 尿样处理和 MCOP assay 经验，但仍需确认 CRC 病例的残余尿量和可匹配对照数。

## 申请目的

请求 WHI/CCC 协助完成样本可行性查询，并申请对可用 CRC nested case-control 样本进行 MCOP assay。首轮只申请 MCOP；若样本量和剩余体积允许，可同步测 MCNP、MiNP 作为 DINP/DIDP exposure-axis QC，但不改变 MCOP 主假设。

## 需要 WHI 确认的字段

1. BMD subcohort 中 incident CRC cases 的尿样可用数；
2. 每个病例的 urine collection time point、sample-to-diagnosis interval 和剩余体积；
3. risk-set matched cancer-free controls 的可用数与 set structure；
4. 样本是否已被 AS458/AS498 或其他 ancillary studies 消耗；
5. 允许的剩余尿样体积、冻融记录、批次信息和可用协变量；
6. MCOP assay 的 LLOQ、QC 规则和预计检测成功率。

## 预先冻结的分析

主暴露为 `log2(MCOP)`，主结局为 incident invasive CRC，主模型为 conditional logistic regression，按 WHI matched set 分层。协变量包括 age、BMI、smoking、alcohol、physical activity、SES、sex、race、urinary creatinine 和 assay batch。预先定义 lag ≥2 年、lag ≥5 年及重复尿样平均值分析。

## 成功标准

不把 P<0.05 作为唯一标准。主要看：

- MCOP–CRC 方向是否为正；
- OR 是否与 1.15–1.30 的中等正关联相容；
- 95% CI 是否排除强烈反向效应；
- lag 和 repeated-urine 分析是否同向；
- 结果是否不依赖单一 assay batch 或少数病例。

## 可行性与资源

WHI 官方流程要求 sponsoring investigator、Ancillary Study Committee 审核和独立经费；应先使用 WHI Query Builder / specimen-results 页面完成样本核对，再提交正式 AS proposal。项目在 biospecimen access 确认前不进行真实 WHI 回归，不使用 TCGA/PPI/GO 作为替代终点。

## 相关入口

- [WHI AS498](https://www.whi.org/study/498)
- [WHI AS458](https://www.whi.org/study/458)
- [WHI ancillary study application path](https://www.whi.org/plan-a-study)
- [WHI specimen results](https://www.whi.org/datasets/specimen-results)
