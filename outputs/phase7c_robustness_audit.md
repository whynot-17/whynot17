# Phase 7C robustness audit

本轮只审计 Phase 7C 已有 NES，不重新筛药。

## 规则

- trajectory leave-one-out：去掉任意一条 trajectory 后，模块方向仍需完全一致。
- GEO dataset leave-one-out：分别去掉 GSE77932、GSE42387 或 GSE119603 后，剩余 trajectory 方向仍需完全一致。
- 重点模块同时要求完整分析中 |median NES| ≥ 0.75 且至少 2/6 trajectory 的 FDR q≤0.25。
- 高重叠 gene sets 用 overlap coefficient ≥0.60 聚类，防止同一生物学被 Reactome/CORUM/co-essentiality 重复计数。

重点模块数：290；冗余簇数：140；同时通过两类 leave-out 的模块数：101。

## 同时通过两类 leave-out 的模块

| Module | Collection | Direction | full consistency | median NES | FDR≤0.25 hits |
|---|---|---|---:|---:|---:|
| COESSENTIALITY::1400::biological process:endocytic recycling | COESSENTIALITY | positive | 1.00 | 1.13 | 2 |
| COESSENTIALITY::1255::biological process:DNA replication initiation | COESSENTIALITY | negative | 1.00 | -0.95 | 2 |
| REACTOME::REACTOME_FORMATION_OF_DEFINITIVE_ENDODERM | REACTOME | negative | 1.00 | -1.08 | 2 |
| COESSENTIALITY::132::biological process:DNA recombination | COESSENTIALITY | negative | 1.00 | -1.15 | 2 |
| REACTOME::REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS | REACTOME | negative | 1.00 | -1.24 | 2 |
| REACTOME::REACTOME_PURINE_SALVAGE | REACTOME | negative | 1.00 | -1.37 | 2 |
| COESSENTIALITY::1787::molecular function:transcription factor activity, RNA polymerase II proximal promoter sequence-specific DNA binding | COESSENTIALITY | positive | 1.00 | 1.32 | 1 |
| COESSENTIALITY::1625::biological process:COPII vesicle coating | COESSENTIALITY | positive | 1.00 | 1.15 | 1 |
| COESSENTIALITY::1224::biological process:inorganic anion transport | COESSENTIALITY | positive | 1.00 | 1.14 | 1 |
| COESSENTIALITY::238::biological process:response to nutrient levels | COESSENTIALITY | positive | 1.00 | 0.99 | 1 |
| COESSENTIALITY::1773::biological process:inorganic anion transport | COESSENTIALITY | positive | 1.00 | 0.98 | 1 |
| COESSENTIALITY::1240:: | COESSENTIALITY | positive | 1.00 | 0.94 | 1 |
| COESSENTIALITY::914::biological process:ATP metabolic process | COESSENTIALITY | positive | 1.00 | 0.80 | 1 |
| CORUM::6664::STAGA complex, SPT3-linked | CORUM | negative | 1.00 | -0.78 | 1 |
| CORUM::1495::PID complex | CORUM | negative | 1.00 | -0.81 | 1 |
| REACTOME::REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION_TERMINATION | REACTOME | negative | 1.00 | -0.83 | 1 |
| REACTOME::REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING | REACTOME | negative | 1.00 | -0.97 | 1 |
| REACTOME::REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT | REACTOME | negative | 1.00 | -0.97 | 1 |
| REACTOME::REACTOME_PROCESSING_OF_CAPPED_INTRONLESS_PRE_MRNA | REACTOME | negative | 1.00 | -1.03 | 1 |
| REACTOME::REACTOME_SARS_COV_1_TARGETS_HOST_INTRACELLULAR_SIGNALLING_AND_REGULATORY_PATHWAYS | REACTOME | negative | 1.00 | -1.09 | 1 |
| CORUM::742::eIF3 complex (EIF3S6, EIF3S5, EIF3S4, EIF3S3, EIF3S6IP, EIF3S2, EIF3S9, EIF3S12,  EIF3S10, EIF3S8,  EIF3S1, EIF3S7) | CORUM | negative | 1.00 | -1.14 | 1 |
| CORUM::1097::eIF3 complex (EIF3S6, EIF3S5, EIF3S4, EIF3S3, EIF3S6IP, EIF3S2, EIF3S9, EIF3S12,  EIF3S10, EIF3S8,  EIF3S1, EIF3S7, PCID1) | CORUM | negative | 1.00 | -1.34 | 1 |
| COESSENTIALITY::1403::cellular component:keratin filament | COESSENTIALITY | negative | 1.00 | -1.47 | 1 |
| COESSENTIALITY::979:: | COESSENTIALITY | positive | 1.00 | 1.37 | 0 |
| COESSENTIALITY::1218::biological process:intra-Golgi vesicle-mediated transport | COESSENTIALITY | positive | 1.00 | 1.27 | 0 |
| REACTOME::REACTOME_KINESINS | REACTOME | positive | 1.00 | 1.17 | 0 |
| COESSENTIALITY::906:: | COESSENTIALITY | positive | 1.00 | 1.11 | 0 |
| COESSENTIALITY::150::biological process:type I interferon signaling pathway | COESSENTIALITY | positive | 1.00 | 1.09 | 0 |
| CORUM::6304::DDX11-Ctf18-RFC complex | CORUM | positive | 1.00 | 1.09 | 0 |
| COESSENTIALITY::1541:: | COESSENTIALITY | positive | 1.00 | 1.09 | 0 |
| COESSENTIALITY::1493::biological process:thyroid gland development | COESSENTIALITY | positive | 1.00 | 1.09 | 0 |
| COESSENTIALITY::1644::biological process:ATP metabolic process | COESSENTIALITY | positive | 1.00 | 1.06 | 0 |
| REACTOME::REACTOME_FATTY_ACIDS | REACTOME | positive | 1.00 | 1.04 | 0 |
| COESSENTIALITY::190::biological process:ATP metabolic process | COESSENTIALITY | positive | 1.00 | 1.04 | 0 |
| CORUM::443::BP-SMAD complex | CORUM | positive | 1.00 | 1.03 | 0 |
| COESSENTIALITY::1607:: | COESSENTIALITY | positive | 1.00 | 1.03 | 0 |
| COESSENTIALITY::195::biological process:ATP metabolic process | COESSENTIALITY | positive | 1.00 | 1.02 | 0 |
| REACTOME::REACTOME_COPI_DEPENDENT_GOLGI_TO_ER_RETROGRADE_TRAFFIC | REACTOME | positive | 1.00 | 1.02 | 0 |
| COESSENTIALITY::83::biological process:ATP metabolic process | COESSENTIALITY | positive | 1.00 | 1.02 | 0 |
| REACTOME::REACTOME_BIOSYNTHESIS_OF_EPA_DERIVED_SPMS | REACTOME | positive | 1.00 | 1.01 | 0 |

## 解释

- 通过 leave-out 只说明信号不是由单一 trajectory 或单一 GEO 数据集完全驱动；仍不能证明 acquired OXA-R causality。
- 通过 redundancy collapse 后才把模块簇作为下一阶段药物映射单位。
- 若某模块只在完整数据中显著、但 leave-out 失败，应降为 exploratory，不进入主候选。
