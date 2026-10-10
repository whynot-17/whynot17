# 23 个稳定复现基因的 GO/KEGG 富集与 PPI 网络分析

- 日期：2026-10-10
- 物种：Homo sapiens（taxon 9606）；基因 ID：HGNC gene symbol。
- 查询集：23 个在至少 3/4 GEO accession 达到此前复现标准的基因；背景集：本项目 113 个 DINP×CRC 交集基因。查询集全部包含于背景集。
- GO/KEGG：使用 g:Profiler g:GOSt，来源 GO:BP、GO:MF、GO:CC、KEGG；custom background（113 symbols）；g:SCS 多重检验校正；显著阈值为校正后 p < 0.05。ORA 不区分肿瘤上调/下调方向。
- STRING PPI：STRING v12.5，人类；功能关联网络与物理相互作用网络分开导出；置信分数阈值 700/1000；STRING PPI enrichment 使用 113 基因背景中成功映射到 STRING 的蛋白 ID。功能关联包含直接物理证据以外的已知/预测功能联系，不能将每条边都解释为直接结合。
- PPI 网络节点仅为 23 个查询基因；add_nodes=0，不添加外部互作蛋白。网络指标在指定 score 阈值下计算；degree/betweenness 是拓扑描述，不代表因果或调控方向。

## 显著富集项数量

| 数据源 | g:SCS < 0.05 项数 |
|---|---:|
| GO:BP | 0 |
| GO:MF | 0 |
| GO:CC | 0 |
| KEGG | 0 |

BH-FDR 结果作为敏感性分析单独报告；与 g:SCS 分开，不混作主分析。

| 数据源 | BH-FDR < 0.05 项数 | 最小 BH-FDR |
|---|---:|---:|
| GO:BP | 0 | 0.451 |
| GO:MF | 0 | 0.562 |
| GO:CC | 0 | 0.783 |
| KEGG | 0 | 0.695 |

## GO 代表性条目

无 GO 条目达到 g:SCS < 0.05。

## KEGG 显著通路

无 KEGG 通路达到 g:SCS < 0.05。

### 最靠前但未达显著的条目（仅作探索）

以下按 BH-FDR 排序；全部仍未达到 0.05，不作为显著富集结论。

| 数据源 | 条目 | BH-FDR p | 重叠基因 |
|---|---|---:|---|
| GO:BP | G1/S transition of mitotic cell cycle (GO:0000082) | 0.451 | FHL1;INHBA |
| GO:BP | G2/M transition of mitotic cell cycle (GO:0000086) | 0.451 | FHL1 |
| GO:BP | protein polyubiquitination (GO:0000209) | 0.451 | UBE2C |
| GO:BP | mitotic cell cycle (GO:0000278) | 0.451 | FHL1;INHBA;UBE2C |
| GO:BP | retinoid metabolic process (GO:0001523) | 0.451 | RARRES2;RETSAT |
| GO:CC | ubiquitin ligase complex (GO:0000151) | 0.783 | UBE2C |
| GO:CC | nuclear ubiquitin ligase complex (GO:0000152) | 0.783 | UBE2C |
| GO:CC | extracellular region (GO:0005576) | 0.783 | BMP2;DES;FABP1;IL18;INHBA;MUC2;PCK1;PLAC8;RARRES2;TIMP1;TPSAB1;ZP3 |
| GO:CC | basement membrane (GO:0005604) | 0.783 | TIMP1 |
| GO:CC | extracellular space (GO:0005615) | 0.783 | BMP2;DES;FABP1;IL18;INHBA;PCK1;RARRES2;TIMP1;TPSAB1;ZP3 |
| GO:MF | magnesium ion binding (GO:0000287) | 0.562 | PCK1 |
| GO:MF | G-protein alpha-subunit binding (GO:0001965) | 0.562 | RGS2 |
| GO:MF | protease binding (GO:0002020) | 0.562 | TIMP1 |
| GO:MF | chromatin binding (GO:0003682) | 0.562 | FABP1;NR3C1;PLAC8;THRB |
| GO:MF | catalytic activity (GO:0003824) | 0.562 | AOC3;CPT2;ELOVL6;PCK1;RETSAT;RGS2;SCD;TPSAB1;UBE2C |
| KEGG | Glycolysis / Gluconeogenesis (KEGG:00010) | 0.695 | PCK1 |
| KEGG | Citrate cycle (TCA cycle) (KEGG:00020) | 0.695 | PCK1 |
| KEGG | Fatty acid elongation (KEGG:00062) | 0.695 | ELOVL6 |
| KEGG | Glycine, serine and threonine metabolism (KEGG:00260) | 0.695 | AOC3 |
| KEGG | Tyrosine metabolism (KEGG:00350) | 0.695 | AOC3 |

## PPI 网络摘要

- functional（score≥0.700）：23 个查询节点、3 条边，平均度 0.26，连通分量 20 个；分量规模 [3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]。
- physical（score≥0.700）：23 个查询节点、1 条边，平均度 0.09，连通分量 22 个；分量规模 [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]。
- functional（score≥0.400）：23 个查询节点、16 条边，平均度 1.39，连通分量 13 个；分量规模 [10, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]。
- physical（score≥0.400）：23 个查询节点、1 条边，平均度 0.09，连通分量 22 个；分量规模 [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]。
- score≥0.400 功能网络的度最高节点：FABP1 (degree=5), PCK1 (degree=5), CPT2 (degree=4), ELOVL6 (degree=4), SCD (degree=4).
- score≥0.400 功能网络最大连通分量：CPT2, DDIT4, ELOVL6, FABP1, IL18, NR3C1, PCK1, RGS2, SCD, TIMP1.
- score≥0.700 功能网络边：SCD–ELOVL6 (0.982); FABP1–CPT2 (0.723); PCK1–SCD (0.714).
- score≥0.400 物理网络边：SCD–ELOVL6 (0.900).
- STRING PPI enrichment（自定义背景，score≥0.700）：{"number_of_nodes": 23, "number_of_edges": 3, "average_node_degree": 0.26, "local_clustering_coefficient": 0.174, "expected_number_of_edges": 2, "p_value": 0.3732}
- STRING PPI enrichment（自定义背景，score≥0.400）：{"number_of_nodes": 23, "number_of_edges": 16, "average_node_degree": 1.39, "local_clustering_coefficient": 0.357, "expected_number_of_edges": 13, "p_value": 0.2821}

## 文件清单

- GO_KEGG_all_results.tsv / GO_KEGG_significant_gSCS.tsv：完整结果与显著项。
- GO_KEGG_all_results_FDR_sensitivity.tsv / GO_KEGG_significant_FDR_sensitivity.tsv：BH-FDR 敏感性分析。
- GO_KEGG_top_ranked_non_significant.tsv：每个数据源 BH-FDR 最靠前的未显著条目。
- GO_significant_representative_terms.tsv：GO 重叠条目聚类后的代表项。
- KEGG_significant_gSCS.tsv：KEGG 显著通路单独表。
- PPI_edges_*_score_ge_0.700.tsv 与 PPI_edges_*_score_ge_0.400.tsv：STRING 边表，可导入 Cytoscape。
- PPI_functional_node_metrics.tsv、PPI_physical_node_metrics.tsv：节点度、中心性、分量与模块。
- PPI_*_network*.graphml：图网络文件，可导入 Cytoscape/Gephi。
- STRING_PPI_enrichment_custom_background*.tsv：STRING PPI enrichment 数值，使用 113 基因自定义背景。
- STRING_mapping_background_113.tsv、STRING_mapping_query_23.tsv：STRING ID 映射。
- STRING_request_parameters.json：查询及背景的 STRING ID 与参数。
- query_23_genes.txt、background_113_genes.txt：本次实际使用的输入集。
- gprofiler_request_g_SCS.json、gprofiler_response_g_SCS.json、gprofiler_request_fdr.json、gprofiler_response_fdr.json、analysis_metadata.json：API 参数、原始响应与版本/统计元数据。
- 复现脚本：analysis/enrichment_ppi_23_vs_113.py。

## 解释边界

这是基于阈值化基因集的过度代表分析（ORA），结果表示注释条目的基因在 23 个基因中相对 113 个背景基因富集，不表示通路被激活或抑制。GO 各项存在层级和重叠，代表表按基因交集进行简化。STRING 功能关联网络汇总多种证据来源；应区分功能关联与物理相互作用。
