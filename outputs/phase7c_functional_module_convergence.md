# Phase 7C：OXA-R trajectory-conditioned functional-module convergence

## 分析定义

- 输入：6 条 parental→OXA-R trajectory 在 80 个 DepMap CRC 模型上的完整 gene-level vulnerability ranking。
- 方法：对 Reactome、Hallmark、curated modules、CORUM protein complexes 和 public co-essentiality modules 做 preranked GSEA。
- 正方向：NES > 0 表示更像该 OXA-R trajectory 的模型对模块内基因总体更依赖。
- 这一阶段不做药物筛选，也不把单个基因命中写成机制结论。

## Gene-set universe

{"REACTOME": 1818, "HALLMARK": 50, "CURATED": 22, "CORUM": 1592, "COESSENTIALITY": 1878}

## Universal convergence candidates

| Module | Collection | + trajectories | - trajectories | consistency | median NES | FDR≤0.25 hits |
|---|---:|---:|---:|---:|---:|---:|
| COESSENTIALITY::979:: | COESSENTIALITY | 6 | 0 | 1.00 | 1.37 | 0 |
| COESSENTIALITY::1787::molecular function:transcription factor activity, RNA polymerase II proximal promoter sequence-specific DNA binding | COESSENTIALITY | 6 | 0 | 1.00 | 1.32 | 1 |
| COESSENTIALITY::1218::biological process:intra-Golgi vesicle-mediated transport | COESSENTIALITY | 6 | 0 | 1.00 | 1.27 | 0 |
| REACTOME::REACTOME_KINESINS | REACTOME | 6 | 0 | 1.00 | 1.17 | 0 |
| COESSENTIALITY::1625::biological process:COPII vesicle coating | COESSENTIALITY | 6 | 0 | 1.00 | 1.15 | 1 |
| COESSENTIALITY::1224::biological process:inorganic anion transport | COESSENTIALITY | 6 | 0 | 1.00 | 1.14 | 1 |
| COESSENTIALITY::1400::biological process:endocytic recycling | COESSENTIALITY | 6 | 0 | 1.00 | 1.13 | 2 |
| COESSENTIALITY::906:: | COESSENTIALITY | 6 | 0 | 1.00 | 1.11 | 0 |
| COESSENTIALITY::150::biological process:type I interferon signaling pathway | COESSENTIALITY | 6 | 0 | 1.00 | 1.09 | 0 |
| CORUM::6304::DDX11-Ctf18-RFC complex | CORUM | 6 | 0 | 1.00 | 1.09 | 0 |
| COESSENTIALITY::1541:: | COESSENTIALITY | 6 | 0 | 1.00 | 1.09 | 0 |
| COESSENTIALITY::1493::biological process:thyroid gland development | COESSENTIALITY | 6 | 0 | 1.00 | 1.09 | 0 |
| COESSENTIALITY::1644::biological process:ATP metabolic process | COESSENTIALITY | 6 | 0 | 1.00 | 1.06 | 0 |
| REACTOME::REACTOME_FATTY_ACIDS | REACTOME | 6 | 0 | 1.00 | 1.04 | 0 |
| COESSENTIALITY::190::biological process:ATP metabolic process | COESSENTIALITY | 6 | 0 | 1.00 | 1.04 | 0 |
| CORUM::443::BP-SMAD complex | CORUM | 6 | 0 | 1.00 | 1.03 | 0 |
| COESSENTIALITY::1607:: | COESSENTIALITY | 6 | 0 | 1.00 | 1.03 | 0 |
| COESSENTIALITY::195::biological process:ATP metabolic process | COESSENTIALITY | 6 | 0 | 1.00 | 1.02 | 0 |
| REACTOME::REACTOME_COPI_DEPENDENT_GOLGI_TO_ER_RETROGRADE_TRAFFIC | REACTOME | 6 | 0 | 1.00 | 1.02 | 0 |
| COESSENTIALITY::83::biological process:ATP metabolic process | COESSENTIALITY | 6 | 0 | 1.00 | 1.02 | 0 |

## Interpretation guardrails

- `universal_consistent_5of6` 和 `universal_consistent_4of6` 是方向一致性描述，不等于 acquired OXA-R causality。
- CORUM/co-essentiality 模块用于把分散的 gene-level signals 聚合到功能层；它们不能替代 paired parental/OXA-R CRISPR。
- subtype patterns are exploratory cluster contrasts and require independent models or paired screens for confirmation.
