# MCOP/DINP 独立 CRC 复制队列扫描

日期：2026-08-22

## 结论先行

首选 **Women’s Health Initiative（WHI）Bone Mineral Density subcohort**，不是把 CHMS 当作第二个 NHANES，也不是继续在横断面病例样本里重复 MCOP。

WHI 已经具备：

- 前诊断尿样和长期随访；
- adjudicated incident CRC；
- 181 名 CRC 病例及配对对照的代谢组样本框架；
- 同一 WHI 资源已经完成包含 MCOP 在内的 13 种尿 phthalate 检测。

因此，最现实的 Phase 3A 是：**向 WHI 申请 CRC nested case-control 的残余尿样，按既有 CDC phthalate 方法测 MCOP，做预先注册的外部复制。** 这不是现成公开分析，样本量和残余尿量必须先向 WHI Query Builder/CCC 确认。

## 候选排序

| 候选 | 设计 | MCOP/尿样证据 | CRC 结局证据 | 判定 |
| --- | --- | --- | --- | --- |
| WHI BMD subcohort | 前瞻性 nested case-control | WHI 已在相关子样本测过 13 种 phthalates，包括 MCOP；重复尿样可用 | 181 incident CRC cases，样本平均在诊断前 7.2 年 | **首选** |
| CHMS | 人群生物监测调查 | 加拿大官方 dashboard 已测 MCIOP/MCOP，2018–2019 检出率 87.4% | 公开资源是暴露分布，不是可直接链接的 incident CRC 队列 | 暴露校准，不是主复制 |
| PLCO | 前瞻性癌症筛查队列 | 公开 CDAS 资料强调 blood、buccal、tumor specimens；未确认有可申请尿 MCOP 样本 | CRC incidence 很强 | 备用，先问 NCI 是否有尿样 |
| 台湾 KMU | 临床横断面病例-对照 | 测的是 MEHP，不是 MCOP | 122 CRC、19 adenoma、80 healthy | 机制/方向参考，不算 MCOP 复制 |
| WHI AS458 / LIBCSP | 前瞻性或病例-对照癌症样本 | MCOP 已测 | 结局主要是 breast cancer，不是 CRC | 方法学先例，不是 CRC 复制 |

## 为什么 WHI 是第一优先级

WHI 的 AS498 已公开描述了同一 BMD subcohort 中的代谢物-癌症研究：758 名 adjudicated cancer women，其中 181 名 CRC，另有 758 名对照；样本平均在诊断前 7.2 年收集。该研究实际使用了血清和尿液样本，但没有把 MCOP 作为 CRC 暴露分析对象。

WHI 的 AS458 则明确使用前瞻性 baseline/year-3 尿样，测量 13 种 urinary phthalate metabolites；官方项目页和已发表论文均列出 MCOP。也就是说，MCOP 的实验方法、质量控制和 WHI 样本处理已经有现成先例，不需要从头开发方法。

MCOP 的命名申请时要同时写：`MCOP`、`MCiOP`、`mCIOP`、`mono(carboxyisooctyl) phthalate`。IARC Exposome-Explorer 将其定义为 DINP 的 secondary metabolite，并且当前 MCOP 页面没有列出 cancer-association study，说明 CRC 复制仍具有新颖性，但不能把数据库空表写成“绝对没有文献”。

## 建议的 Phase 3A 方案

### 1. 先做 WHI 可行性查询

必须先得到以下四个数字：

1. BMD subcohort 中 incident CRC cases 的尿样可用数；
2. 每名病例诊断前可用尿样的时间点和剩余体积；
3. 可匹配的 cancer-free controls 数量；
4. 这些样本是否已被 AS458/AS498 消耗或存在 assay overlap。

WHI 官方流程要求寻找 sponsoring investigator，并由 Ancillary Study Committee 审核新增 biospecimen assay；申请前应使用 Query Builder 和 specimen-results 页面核对样本可得性。

### 2. 主分析预注册

- Outcome：incident invasive CRC；排除 baseline 已患癌者；
- Exposure：MCOP，连续 `log2(MCOP)` 为主；
- 主要效应：每翻倍 OR；
- 次要效应：预先定义四分位、colon/rectum、lag ≥2 年和 lag ≥5 年；
- 样本设计：年龄、入组时间、种族等匹配因素下的 conditional logistic regression；
- 协变量：BMI、吸烟、饮酒、体力活动、教育/SES、饮食相关变量、尿肌酐和 assay batch；
- 若有 baseline 与 year-3 尿样：优先报告单次 baseline，重复尿样平均值作为预先定义的稳定性分析；
- 不把 MiNP/MCNP/总 phthalate burden 作为新的主要假设，最多作为暴露轴 QC。

### 3. 判定标准

外部复制不应再使用“P<0.05 才算成功”的单一门槛。建议在分析前冻结：

- 方向：OR/HR > 1；
- 效应：点估计是否落在 1.15–1.25 附近或更高；
- 95% CI 是否与中等正关联相容；
- single-sample 与 repeated-sample 结果是否同向；
- 结果是否在 lag 分析中保持，而不是只在诊断临近的尿样中出现。

181 个 CRC 病例足以做有价值的方向性复制，但对 OR 1.15–1.25 这种小效应不能预设有很高 power；正式 power 计算要等实际 MCOP 分布、匹配比例、缺失率和可用病例数确认后再做。

## 不建议现在做的事

- 不再把 NHANES 的更多分层/变换当成独立验证；
- 不把 CHMS 的 MCOP 分布直接当 CRC 结局；
- 不用台湾 MEHP 横断面结果替代 MCOP；
- 不因为 WHI 的 181 例规模有限就改回另一个没有尿样的 cohort；
- 不先做 TCGA、PPI 或湿实验。

## 证据与入口

- [WHI AS498：Metabolite Predictors of Breast and Colorectal Cancer Risk](https://www.whi.org/study/498)
- [WHI AS458：Phthalate metabolites and breast cancer risk](https://www.whi.org/study/458)
- [WHI Plan an Ancillary Study](https://www.whi.org/plan-a-study)
- [WHI specimen results](https://www.whi.org/datasets/specimen-results)
- [WHI phthalate biomarker paper](https://pubmed.ncbi.nlm.nih.gov/30629220/)
- [Canadian Biomonitoring Dashboard: MCIOP](https://health-infobase.canada.ca/biomonitoring/?mx=urine&v=LAB_MCIO)
- [PLCO Cancer Data Access System](https://cdas.cancer.gov/plco/)
- [Taiwan urinary MEHP–CRC study](https://pubs.acs.org/doi/10.1021/acs.jafc.1c00953)
- [IARC Exposome-Explorer MCOP entry](https://exposome-explorer.iarc.fr/compounds/1455)

## 当前状态

**推荐：WHI BMD subcohort → MCOP assay → incident CRC external replication。**

当前仍是 feasibility-ready，不是 access-confirmed。下一步不是继续分析 NHANES，而是准备一页 WHI sponsor/ASC pre-proposal，并先拿到上述四个样本可得性数字。
