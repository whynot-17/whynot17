# Phase 7A：GDSC phenotype-direct OXA pharmacogenomic validation

## 结论

本轮在可提供 Oxaliplatin 的 GDSC CRC 细胞系中，按细胞系内平均 LN_IC50 进行分析；共输出 272 个 dataset-drug 组合。LN_IC50 越高代表越不敏感。

这一步验证的是 OXA-R-like phenotype，而不是表达谱定义的 ERAD/DHODH 亚型。因此它是 Phase 7 的独立药敏锚点，不能替代 subtype projection。

## 判定规则

- OXA-R-like：CRC 细胞系 Oxaliplatin LN_IC50 高于该数据集 CRC 中位数；OXA-S-like：不高于中位数。
- collateral sensitivity：OXA-R-like 组的候选药 LN_IC50 更低，即 `delta_R_minus_S < 0`；Spearman OXA/药物敏感性相关也应为负。
- 预设复制要求：若两个数据集都含 Oxaliplatin，则要求 GDSC1 与 GDSC2 同方向；若数据集没有 Oxaliplatin，则标为不可复制，不把缺失当阴性。

## 当前结果

数据可用性：

- GDSC1: Oxaliplatin not found

按‘两个 GDSC 数据集方向一致’的严格盲筛门槛，当前没有通过者；这不是证明没有药物效应，而是说明不能把单一数据集的负相关当作可靠 Drug X。

GDSC2/可用数据中，方向最偏向 OXA-R collateral sensitivity 的前 5 个药物为：
- WEHI-539：Spearman r=-0.273，BH-p=0.080，R−S median LN_IC50=-1.186，n=45
- UMI-77：Spearman r=-0.105，BH-p=0.531，R−S median LN_IC50=-0.490，n=41
- Navitoclax：Spearman r=-0.088，BH-p=0.575，R−S median LN_IC50=-0.618，n=46
- Dasatinib：Spearman r=-0.079，BH-p=0.619，R−S median LN_IC50=0.133，n=45
- Acetalax：Spearman r=-0.005，BH-p=0.973，R−S median LN_IC50=-0.585，n=45

已测到的预先指定候选均未呈现 collateral-sensitivity 方向：
- Bortezomib：Spearman r=0.320，BH-p=0.044，R−S median LN_IC50=0.401
- Leflunomide：Spearman r=0.436，BH-p=0.006，R−S median LN_IC50=0.392

### 预先指定候选

候选药物的 GDSC 覆盖和方向见 `phase7_gdsc_named_candidate_validation.csv`。若候选不在 GDSC1/2，则记为 pharmacogenomic unavailable，而不是阴性。

## 重要限制

1. GDSC1/GDSC2 是药敏关联数据，不是 acquired OXA-R 与 parental 的配对实验；只能支持 phenotype-level prioritization。
2. 细胞系之间的 lineage、增殖速度与普遍药敏会造成混杂；下一步若有独立表达谱，应对 subtype score 做连续建模并控制 growth-related covariates。
3. 本轮不能证明 ERAD-high 或 DHODH-high 选择性，除非表达谱投影成功并与药敏矩阵在同一细胞系 ID 上可靠合并。

## 文件

- `phase7_gdsc_crc_oxaliplatin_phenotype.csv`
- `phase7_gdsc_drug_oxa_association_by_dataset.csv`
- `phase7_gdsc_cross_dataset_replication.csv`
- `phase7_gdsc_named_candidate_validation.csv`
- `phase7_gdsc_top_negative_oxa_associations.csv`
- `phase7_gdsc_blind_collateral_sensitivity_gate.csv`