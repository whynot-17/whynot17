# CNP0004138 队列复用分析报告

## 结论摘要
- scRNA 患者内配对显示，治疗后上皮样细胞比例没有可检出的系统性改变；但上皮身份、WNT/干性、增殖、缺氧应激和 EMT 程序分数总体下降，其中 IFN/ISG 分数未见改变。由于公开条码没有 XFZ 患者编号到 TRG 应答组的映射，以上不能回答这些残留上皮状态是否在 NR 特异富集，也不能把上皮样群直接称为 CNV 确认的残留恶性细胞。
- CAF 样群在配对术后样本中增加，CAF_PI16/SLIT2/ADAM28 程序上升；CAF_BMP4/FAP/MMP1 程序总体下降。髓系总比例变化不显著；NLRP3-like 均分向下但不显著，C1QA-like 向上但经多重校正未达显著，SPP1-like 分数上升。所有 scRNA 配对统计以患者为单位（27 对）。
- 27 位配对患者的术后样本级共丰度中，上皮 EMT 与 CAF_FAP 的 Spearman 相关为 rho=0.79（原定义模块）；去掉两边重叠的 MIR4435-2HG 后仍为 rho=0.66，p=0.00017。该结果支持患者间共变化，不证明细胞直接互作。
- 四张术后空间切片呈现 CR 与 NR 的方向性差异：两张 NR 切片的上皮身份高分 bin 周围，CAF_BMP4/FAP/MMP1、PI16/SLIT2、NLRP3-like 与 C1QA-like 等程序均高于该切片其他有效 bin；CR 切片中相应方向多为零或负。空间结果是每张切片内的局部邻近支持，每组只有两位患者，不构成患者层面的应答组推断。

## 数据与统计单位
- 原文队列为 29 位患者，27 份治疗前、29 份治疗后，27 位可配对患者；有 TRG 的患者为 CR 7、PR 15、NR 5。空间组学只覆盖 4 个术后样本（CR 两张、NR 两张；无 PR）。原文将恶性上皮细胞通过 CNV 与聚类识别。来源：[Qin et al., Cell Reports Medicine (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10591051/)。
- CNGBdb 公布矩阵含 378,344 个条码，按患者×阶段识别出 61 个样本组：PRET 27、POSTT 29、POSTN 5；POSTN 是邻近正常样本，不参加肿瘤配对比较。MatrixMarket 中 159 个技术文库已合并到 61 个患者阶段样本。
- 应用论文所述检出基因数（200–6000）和线粒体比例（≤20%）阈值后，矩阵中 378,344 个条码全部通过（每细胞检出 201–5,999 个基因）；这比原文报告的 333,285 个高质量细胞多 45,059 个。公开矩阵与原文最终细胞集合未能完全对齐，可能涉及额外筛选/版本差异；本报告保留此差异，不宣称复现作者最终 QC 集合。
- 细胞仅用于估计每个样本的细胞比例与程序分数；统计检验的 n 是患者数，配对前后比较均使用一位患者一个 PRET 与一个 POSTT 值。治疗后共丰度相关也用患者级汇总（n=27）。

## 1. 化疗后保留的上皮样状态
- 过滤后 POSTT 中识别出 36,114 个上皮样标记细胞（29 个治疗后样本）。这是表达标志定义的 epithelial-like 计数，不是 CNV 确认的残留癌细胞数。
- 上皮样细胞占样本全部 QC 细胞比例：配对差中位数 -0.002，n=27，Wilcoxon p=0.897，BH q=0.966。
- 上皮身份分数：配对差中位数 -0.181，n=27，Wilcoxon p=0.00325，BH q=0.00649。
- EMT：配对差中位数 -0.032，n=27，Wilcoxon p=0.0245，BH q=0.0404；WNT/干性：配对差中位数 -0.116，n=27，Wilcoxon p=7.99e-06，BH q=5.59e-05；增殖：配对差中位数 -0.069，n=27，Wilcoxon p=0.00271，BH q=0.00583；缺氧/应激：配对差中位数 -0.049，n=27，Wilcoxon p=0.00461，BH q=0.00861；IFN/ISG：配对差中位数 0.019，n=27，Wilcoxon p=1，BH q=1。BH 校正范围为对应配对指标集合。
- 因 XFZ—TRG 对照缺失，无法判断上述任何细胞状态是否在 NR 比 CR/PR 更突出。不能用 ST-CR/NR 标签给 scRNA 病人反向补组。

## 2. 髓系/巨噬细胞与 CAF 状态
- 髓系细胞总体比例：配对差中位数 -0.000，n=27，Wilcoxon p=0.934，BH q=0.969。巨噬细胞候选程序的每样本平均分数：NLRP3-like（配对差中位数 -0.205，n=27，Wilcoxon p=0.123，BH q=0.181）；C1QA-like（配对差中位数 0.242，n=27，Wilcoxon p=0.0463，BH q=0.072）；SPP1-like（配对差中位数 0.250，n=27，Wilcoxon p=0.00645，BH q=0.0113）。这是连续 marker-module 分数，不等同于已经验证的巨噬细胞亚群比例。
- CAF 样群比例：配对差中位数 0.069，n=27，Wilcoxon p=4.09e-05，BH q=0.000191。候选状态程序的配对变化：PI16（配对差中位数 0.474，n=27，Wilcoxon p=7.99e-06，BH q=5.59e-05）；SLIT2（配对差中位数 0.408，n=27，Wilcoxon p=0.000301，BH q=0.0012）；ADAM28（配对差中位数 0.188，n=27，Wilcoxon p=0.00153，BH q=0.00389）；BMP4（配对差中位数 -0.203，n=27，Wilcoxon p=0.00169，BH q=0.00394）；FAP（配对差中位数 -0.265，n=27，Wilcoxon p=0.000427，BH q=0.0015）；MMP1（配对差中位数 -0.341，n=27，Wilcoxon p=1.35e-05，BH q=7.54e-05）。
- 术后患者级共丰度中，EMT—CAF_FAP 关系最强（原签名 rho=0.79，q=5.4e-05）；剔除两边共用的 MIR4435-2HG 后仍为 rho=0.66，p=0.00017。髓系信号较弱：上皮缺氧/应激与 C1QA-like 为 rho=0.49、q=0.044；上皮身份与 NLRP3-like 为 rho=0.46、q=0.060。其余相关见 `post_patient_state_cooccurrence.csv`；样本级共丰度不代表同一空间位点或因果作用。

## 3. 空间邻近
- 对每张切片，在 panel-positive 的 bins 中将上皮身份分数最高四分位定义为 epithelial-high 区；对每个这样的 bin 取最近 8 个 bin，计算候选髓系/CAF/EMT 程序均分，与该切片其余 panel-positive bins 的邻居分数比较。计数矩阵只解析固定的标志基因 panel；为统一 bin 深度，采用该 panel 的总计数作切片内归一化，不是全转录组 UMI 归一化。
- “epithelial-high”仅是上皮标志富集区，不代表存活肿瘤或癌细胞定位；空间测量单位是约 25 μm 的 bin（非单细胞分割）。邻近统计使用切片内标签置换；置换 p 值是 bin-level 局部关联，不是患者层面的 CR 对 NR 显著性检验。
- 上皮高分区邻居差值（NR1、NR2；随后 CR1、CR2）分别为：CAF_BMP4 0.218, 0.100; -0.104, -0.176；CAF_PI16 0.042, 0.024; -0.047, -0.059；CAF_SLIT2 0.059, 0.030; -0.079, -0.079；NLRP3-like 0.177, 0.080; -0.018, 0.032；C1QA-like 0.035, 0.033; -0.006, 0.005。FAP、MMP1 在两张 NR 切片均为正（FAP 0.135/0.092；MMP1 0.417/0.095），而 CR 接近零或为负。
- 这种方向一致性可作为后续靶向复核的空间候选证据，但病例数仅 CR=2、NR=2，且没有术前切片，不可外推到整个应答组或全体患者。

## 分析定义与可复现文件
- 上皮状态：预设 EMT、IFN/ISG、增殖、缺氧应激、WNT/干性、上皮身份基因集；分数为每细胞 CP10K 后 log1p 表达的基因均分，再汇总至样本/患者。
- 巨噬/髓系与 CAF 状态：依原文描述的 marker 及常见标记建立预设程序；属于候选程序复核，不重跑作者 Seurat 聚类、inferCNV、NMF 或原文图像解卷积。
- 配对检验：双侧 Wilcoxon signed-rank，针对一组样本级指标作 BH 校正；术后程序共丰度以 Spearman 相关并按相关矩阵作 BH 校正。
- 复跑：先执行 `src/analyze_reuse.py`，再执行 `src/summarize_outputs.py`。
- 主要文件：`results/paired_patient_tests.csv`、`results/sample_level_state_summaries.csv`、`results/post_patient_state_cooccurrence.csv`、`results/overlap_free_cooccurrence_sensitivity.csv`、`results/spatial_neighbor_permutation.csv`、`results/spatial_group_direction_summary.csv`、`results/cell_counts_audit.csv`、`results/scRNA_qc_analysis_audit.json`；图件在 `figures/`。

## 补做 NR 分层所缺的数据
需要 Table S1 中 XFZ 编号到 CR/PR/NR 的对照表。作者补充 PDF 的公开下载在本次会话遇到 reCAPTCHA/Cloudflare 验证，我没有尝试绕过验证；用户确认手头没有映射表。因此本次不输出 scRNA 的 NR 特异统计。如果之后把 Table S1 PDF 或一张 XFZ—TRG 对照表放进该 D 盘项目目录，可以在 27 位配对患者上补做 NR/PR/CR 分层模型及患者级效应量。
