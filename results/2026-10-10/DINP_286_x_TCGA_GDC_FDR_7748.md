# DINP 286 × TCGA GDC 7,748 FDR DEG 交集

- 计算日期：2026-10-10
- DINP 输入：D:\DINP_CRC_Frozen_Archive_20260918\outputs\DINP_target_master.csv（286 行，286 个唯一 gene_symbol）
- DEG 输入：D:\Codex\projectless-workspaces\2026-09-08\dinp-crc-replication-b24a399\outputs\dinp_crc_dehp_replication\TCGA_GDC_current_DEG_FDR_threshold.tsv（7748 行，7748 个唯一 gene；来自 FDR 阈值结果）
- 匹配规则：去除首尾空格后按 gene symbol 精确匹配，大小写不敏感；未做别名映射或模糊匹配。
- 交集：**113 个基因**（肿瘤上调 31 个、肿瘤下调 82 个）
- DEG 输入表 fdr_hit 非 TRUE 的行数：0
- DINP 输入 SHA-256：3a1a8d6367330ddafa4bd909434136b9a322ebd8e1770c01f2b1e173a9dba016
- DEG 输入 SHA-256：7975288bc71e9ed80f71f75f8c64f5b32c3cc2e727297dbc820dc808f9b47ebc

交集明细见同目录 DINP_286_x_TCGA_GDC_FDR_7748.tsv。方向定义为 logFC 正值 = Tumor_up，负值 = Tumor_down。
