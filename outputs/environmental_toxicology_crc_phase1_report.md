# 环境毒理学与 CRC：Phase 1 CTD × GeneCards 结果

## 运行状态

本报告由冻结的 CTD chemical hierarchy、GeneCards Disorders-scoped CRC 查询、Fisher/BH-FDR 和 degree-matched permutation 流程生成。第一轮未进行独立 CRC 文献撞车审查。

- CTD 人类 interaction rows: 1,327,615
- CTD unique chemical–gene pairs: 828,394
- Core environmental chemicals with human interactions: 267
- Core chemical classes observed: 14
- GeneCards mapped inputs: GeneCards_Anywhere, GeneCards_Disorders

## Primary top 20

排序预先固定为 BH-FDR 升序、log2(OR) 降序、rank-weighted overlap 降序、overlap 降序；没有按化学物名称或研究热门程度人工挑选。

| unfiltered_rank | ChemicalName | chemical_class | n_ctd_human_genes | crc_overlap | odds_ratio | enrichment_ratio | bh_fdr | rank_weighted_overlap | n_unique_pmids | top_overlap_genes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | butylbenzyl phthalate | phthalates | 214 | 28 | 5.94989 | 5.09632 | 4.42339e-10 | 4.29823 | 54 | PTEN;ATM;TP53;CDH1;EGFR;PDGFRA;AKT1;MIR19A;CTNNB1;PPARG |
| 2 | mono-benzyl phthalate | phthalates | 87 | 18 | 10.1827 | 8.05871 | 8.9924e-10 | 2.50881 | 20 | TP53;CASP8;MALAT1;PPARG;BAX;CCND1;MIR126;MIR143;MIR221;MIR141 |
| 3 | nickel monoxide | heavy_metals | 45 | 11 | 12.4942 | 9.52122 | 1.19078e-06 | 1.68486 | 13 | ATM;CDH1;AKT1;PIK3R1;SMO;BAX;JUN;MEG3;HOTAIRM1;IL1B |
| 4 | Copper | heavy_metals | 3542 | 138 | 1.70479 | 1.51755 | 1.15754e-05 | 21.86 | 284 | MSH6;MLH1;MSH2;EPCAM;CHEK2;POLD1;APC;BRAF;SMAD4;NTHL1 |
| 5 | monobutyl phthalate | phthalates | 84 | 13 | 7.08387 | 6.02804 | 1.16178e-05 | 1.89952 | 26 | TP53;BRCA1;CASP8;BUB1B;NFE2L2;PPARG;BAX;MIR141;NBR2;SPRY4-AS1 |
| 6 | Silver | heavy_metals | 873 | 48 | 2.316 | 2.1416 | 2.99908e-05 | 7.40401 | 46 | MSH6;TP53;PALB2;RNF43;SMAD7;MIR19A;FGFR3;NFE2L2;SOX9;SOD2 |
| 7 | mono-(2-ethylhexyl)phthalate | phthalates | 1118 | 56 | 2.10713 | 1.95101 | 5.67665e-05 | 8.25762 | 94 | CHEK2;POLD1;NTHL1;PTEN;TP53;PIK3CA;PLA2G2A;SCG5;AKT1;CASP8 |
| 8 | ammonium 2,3,3,3-tetrafluoro-2-(heptafluoropropoxy)-propanoate | pfas_perfluoro | 175 | 17 | 4.17555 | 3.78376 | 9.70542e-05 | 2.81123 | 14 | ATM;KRAS;TP53;CDH1;EGFR;NRAS;PDGFRA;CTNNB1;HRAS;PPARG |
| 9 | chromium hexavalent ion | heavy_metals | 1632 | 72 | 1.85704 | 1.7184 | 0.000139045 | 12.8751 | 85 | MSH6;MLH1;MSH2;EPCAM;POLE;CHEK2;POLD1;ATM;TP53;BRCA2 |
| 10 | 9,10-Dimethyl-1,2-benzanthracene | pahs | 79 | 11 | 6.23752 | 5.42348 | 0.000144482 | 1.79038 | 38 | TP53;CDH1;BRCA1;CTNNB1;HRAS;PPARG;MYC;SOD2;BAX;CCND1 |
| 11 | monoethyl phthalate | phthalates | 56 | 9 | 7.36503 | 6.25989 | 0.00028584 | 1.21198 | 14 | CASP8;PPARG;BAX;MIR34A;MIR141;MIR192;MIR25;SPRY4-AS1;CDKN2A |
| 12 | Zinc | heavy_metals | 1755 | 74 | 1.76775 | 1.64235 | 0.000365368 | 10.7922 | 255 | APC;BRAF;PTEN;ATM;PMS1;TP53;BRCA2;PALB2;CDH1;BRCA1 |
| 13 | Volatile Organic Compounds | vocs | 1909 | 77 | 1.68528 | 1.57108 | 0.00101193 | 11.7603 | 6 | MLH1;PMS2;POLE;MUTYH;NTHL1;ATM;PMS1;CTNNA1;CDH1;BRCA1 |
| 14 | lead nitrate | heavy_metals | 68 | 9 | 5.86388 | 5.1552 | 0.00111182 | 1.37406 | 17 | TP53;RB1;EGFR;CTNNB1;BAX;NEAT1;MIR139;DNMT3A;RELA |
| 15 | cobalt ferrite | heavy_metals | 10 | 4 | 25.4676 | 15.5802 | 0.00133234 | 0.643693 | 3 | TP53;CASP8;BAX;RELA |
| 16 | monomethyl phthalate | phthalates | 10 | 4 | 25.4676 | 15.5802 | 0.00133234 | 0.505182 | 4 | PPARG;NEAT1;TINCR;IL1B |
| 17 | mono(2-ethyl-5-hydroxyhexyl) phthalate | phthalates | 42 | 7 | 7.6699 | 6.49174 | 0.00137702 | 0.926512 | 21 | IGF2;PTPRJ;SOD2;MEG3;SOCS2-AS1;SPRY4-AS1;TUG1 |
| 18 | gallium arsenide | heavy_metals | 57 | 8 | 6.26803 | 5.46673 | 0.00145683 | 1.17041 | 1 | TP53;AKT1;CTNNB1;MYC;JUN;NFKB1;SMAD3;MLC1 |
| 19 | CPS 49 | phthalates | 11 | 4 | 21.8284 | 14.1638 | 0.0017273 | 0.560225 | 1 | AKT1;CASP8;STAT3;XIAP |
| 20 | tetrathiomolybdate | heavy_metals | 112 | 11 | 4.19326 | 3.82549 | 0.00195899 | 1.67765 | 12 | TP53;EGFR;AKT1;PIK3R1;NFE2L2;CP;SOD2;BAX;IL2RA;CD40 |

## 主分析摘要

- Primary tested chemicals: 267
- Primary FDR < 0.05: 52
- Primary stable candidates (n_interacting ≥ 20 and overlap ≥ 5): 69

## GeneCards scope stability

| configuration_a | configuration_b | top20_overlap | top20_jaccard |
| --- | --- | --- | --- |
| GeneCards_Disorders_top500 | GeneCards_Disorders_top1000 | 15 | 0.6 |
| GeneCards_Disorders_top500 | GeneCards_Disorders_top2000 | 15 | 0.6 |
| GeneCards_Disorders_top500 | GeneCards_Anywhere_top500 | 11 | 0.37931 |
| GeneCards_Disorders_top500 | GeneCards_Anywhere_top1000 | 11 | 0.37931 |
| GeneCards_Disorders_top500 | GeneCards_Anywhere_top2000 | 7 | 0.212121 |
| GeneCards_Disorders_top1000 | GeneCards_Disorders_top2000 | 20 | 1 |
| GeneCards_Disorders_top1000 | GeneCards_Anywhere_top500 | 10 | 0.333333 |
| GeneCards_Disorders_top1000 | GeneCards_Anywhere_top1000 | 11 | 0.37931 |
| GeneCards_Disorders_top1000 | GeneCards_Anywhere_top2000 | 8 | 0.25 |
| GeneCards_Disorders_top2000 | GeneCards_Anywhere_top500 | 10 | 0.333333 |
| GeneCards_Disorders_top2000 | GeneCards_Anywhere_top1000 | 11 | 0.37931 |
| GeneCards_Disorders_top2000 | GeneCards_Anywhere_top2000 | 8 | 0.25 |
| GeneCards_Anywhere_top500 | GeneCards_Anywhere_top1000 | 14 | 0.538462 |
| GeneCards_Anywhere_top500 | GeneCards_Anywhere_top2000 | 9 | 0.290323 |
| GeneCards_Anywhere_top1000 | GeneCards_Anywhere_top2000 | 14 | 0.538462 |

## 解释边界

结果表示 CTD chemical-interacting gene set 与 GeneCards CRC-associated gene set 的富集关系，不证明真实人群暴露因果、CRC 发病风险或治疗效果。CTD PMID 字段仅用于审计与证据成熟度，第一轮未读文献。

## Reproducibility

- Run timestamp (UTC): 2026-08-22T06:38:36.119958+00:00
- Random seed: 20260822
- Degree-matched permutations: 1000
- Raw CTD archives are excluded from Git by repository .gitignore.
