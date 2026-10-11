# SCENIC：PID1相关巨噬细胞调控网络（患者级分析）

## 设计与样本量

以GSE132465推断GRN并进行motif pruning/AUCell；GSE200997仅用同一组primary regulons打分，未单独推断网络。只纳入Inflammatory与C1QC巨噬细胞，并要求当前注释的QC通过、非Scrublet doublet、无marker lineage conflict。
GSE132465：1256个细胞（C1QC 451；炎症型 805），12位患者两亚群均≥10个细胞。≥5/≥20阈值下分别有17/8位配对患者。
GSE200997：198个细胞（C1QC 173；炎症型 25），仅1个样本两亚群均≥10个细胞（≥5：1；≥20：0），外部队列只能作探索性方向核对。
GSE132465识别并打分140个regulons；采用患者/样本为统计单位，主比较为每位患者内炎症型减C1QC的AUCell差，Wilcoxon并对regulons做BH校正。

## 主要结果

在≥10细胞/亚群的患者级配对分析中，BH-FDR<0.05的差异regulons有57/140个；炎症型高/低方向分别按患者差值统计，见`GSE132465_paired_regulon_stats_min10_cells.tsv`及逐患者差值表。
显著且炎症型活性较高的前五个为ETS2(+), MXD1(+), REL(+), IRF7(+), NFKB1(+)；显著且C1QC较高的前五个为FOSB(+), FOS(+), STAT1(+), JUNB(+), JUN(+)。
细胞数阈值敏感性：≥5、≥10、≥20时分别有63、57、48个regulons达到BH-FDR<0.05；其中41个在三种阈值下均显著，见`GSE132465_regulon_threshold_sensitivity.tsv`。
患者内ΔPID1与Δregulon AUCell的相关分析中，最强的是ETS1(+)（Spearman ρ=0.846，原始P=0.0005211，140个regulon校正后FDR=0.073）；目前没有相关性通过BH-FDR 0.05。
motif-pruned regulon中将PID1列为候选靶基因的有3个：NR3C1, SP3, ZNF580。其中可结合细胞状态差异看作优先候选，但其患者级PID1相关性仍需以表格中的效应量和FDR判断，不能视为调控已被证明。
GSE200997在1个可配对样本中可核对140个regulon方向，其中99个与主队列方向一致；由于仅1个样本达到每亚群≥10个细胞，不将方向一致率称为跨患者复现。

## 解释边界

SCENIC输出的是表达共变驱动的候选TF regulon及motif支持的调控关系，不等同于直接结合或因果调控。患者数较少，尤其外部队列亚群细胞稀疏；细胞数阈值敏感性、患者方向一致率与FDR应共同解读。不同患者/样本作为统计单位，未把细胞当作独立重复。
本次只依赖GSE132465/GSE200997原始UMI矩阵、当前PID1髓系精细注释以及官方人类TF/motif-ranking资源；不读取旧CellChat预测或历史SCENIC结果。

## 资源与软件

pySCENIC 0.12.1：GRNBoost2 → cisTarget motif pruning → AUCell；人类v10 clustered hg38 RefSeq R80 10kb upstream/downstream ranking database与匹配motif注释。资源说明见Aerts lab官方资源库及pySCENIC文档。
- https://github.com/aertslab/pySCENIC
- https://resources.aertslab.org/cistarget/databases/
- https://resources.aertslab.org/cistarget/motif2tf/

关键文件：`GSE132465_PID1_related_candidate_regulons.tsv`、`candidate_regulatory_network_edges_top30.tsv`、`primary_external_regulon_direction_comparison.tsv`及患者级活动/统计表。
