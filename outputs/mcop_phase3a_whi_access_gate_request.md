# MCOP–CRC Phase 3A：WHI access gate 请求包

## 用途

这是 biospecimen/data access 前的可行性查询，不是正式 WHI 回归结果，也不是在请求先行释放个体数据。目标是确认 MCOP assay 是否值得提交正式 ancillary study proposal。

WHI 官方流程要求：需要 sponsoring investigator；涉及使用/产生 biospecimen 数据的 ancillary study 由 Ancillary Study Committee（ASC）审查，并需要独立经费。Query Builder 和 specimen resources 用于前置可行性核对。

## 必须返回的四个数字

请按 WHI 的真实 matched-set definition 返回：

1. `n_crc_urine_eligible`：incident invasive CRC 病例中，诊断前尿样仍可用、剩余体积足够且有可追溯采集日期者；
2. `control_ratio_and_n_sets`：可匹配 controls 数量、实际 1:1 或 1:2 结构，以及 matched-set 数；
3. `prediagnostic_interval_summary`：sample-to-diagnosis interval 的 median、IQR、范围，以及 lag ≥2 年和 lag ≥5 年可用病例数；
4. `residual_volume_and_assay_feasibility`：剩余尿量、冻融/消耗记录、MCOP assay 所需体积、预计检测成功率、LLOQ/LOD、批次和 QC 条件。

## Query Builder / CCC 字段清单

### 病例和结局

- WHI component/cohort：BMD subcohort / AS498 cancer case-control framework；
- incident invasive CRC：首次诊断日期、结肠/直肠部位、adjudication status；
- baseline cancer exclusion status；
- case identifier 与 matched-set identifier（去标识化即可）。

### 尿样和 assay 可行性

- urine collection time point/date；
- sample-to-diagnosis interval；
- urine specimen availability and residual volume；
- freeze–thaw count / prior ancillary consumption；
- whether MCOP/MCiOP/mCIOP has already been measured；
- assay method, LLOQ/LOD, expected detection rate, batch/QC information；
- whether one assay can be run blinded to case/control status.

### 对照和 matching definition

- cancer-free risk-set controls available per case；
- exact matching factors and matching ratio；
- whether age, race, clinic/region, enrollment period or other factors were used to construct the sets；
- whether any proposed adjustment variable is constant within sets.

### 预设协变量可用性

- age, BMI, smoking, alcohol, physical activity, SES；
- race；
- urinary creatinine；
- assay batch；
- missingness and collection-time availability for each field.

## 预先决策门槛

以 OR per MCOP doubling = 1.25 的女性 proxy simulation 为规划依据：

- `≥150` analyzable CRC cases：GREEN；
- `120–149`：值得正式申请，优先争取 1:2 controls；
- `100–119`：仅在 1:2 controls、样本质量和 lag 分布良好时继续；
- `<80`：统计风险明显上升，需重新评估是否仍作为主复制队列。

这些门槛不替代 WHI 实测 power；最终 power 需用真实 MCOP 方差、缺失率和 assay success rate 更新。

## 可直接发送的询问信

**Subject:** Feasibility query: prediagnostic urinary MCOP assay in WHI incident colorectal cancer nested case-control samples

We are preparing a WHI ancillary study to evaluate prediagnostic urinary mono(carboxyisooctyl) phthalate (MCOP/MCiOP) in relation to incident invasive colorectal cancer. The proposed analysis is a prospective nested case-control replication with one case per WHI matched set and conditional logistic regression. MCOP will be analyzed as log2 concentration, with urinary creatinine and pre-specified covariates handled according to the final WHI matched-set definition.

Before submitting a full ASC proposal, could WHI/CCC advise on the feasibility of the following for the AS498/BMD-related CRC case-control resource?

1. Number of incident invasive CRC cases with sufficient residual prediagnostic urine;
2. Number and structure of available risk-set matched controls, including the exact matching factors;
3. Sample-to-diagnosis interval summary and numbers meeting lag ≥2 years and lag ≥5 years;
4. Residual volume, freeze–thaw/prior-use history, MCOP assay feasibility, LLOQ/LOD, expected detection rate, batch/QC requirements and estimated assay cost;
5. Availability and missingness of BMI, smoking, alcohol, physical activity, SES, race and urinary creatinine.

We understand that this feasibility query does not constitute data access and that a formal ancillary study proposal, sponsoring investigator, ASC review and separate funding are required before any new assay or analysis is initiated.

## 用户需要完成的动作

1. 用本人/机构 WHI 账号进入 Query Builder 或联系 WHI Help Desk/CCC；
2. 确认 sponsoring investigator 和机构信息；
3. 将上面询问信发送给 CCC/AS498 相关联系人，并按 WHI 要求提交 Query Builder 结果；
4. 把返回的四个数字、matching factors 和 assay 条件贴回本项目；
5. 我将据此更新正式 sample-size memo、最终模型变量表和 ASC proposal，不会在 access 前运行真实 WHI 回归。

## 官方入口

- [WHI Plan an Ancillary Study](https://www.whi.org/plan-a-study)
- [WHI specimen results and data dictionaries](https://www.whi.org/datasets/specimen-results)
- [WHI AS498](https://www.whi.org/study/498)
