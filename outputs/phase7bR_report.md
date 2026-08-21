# Phase 7B-R：OXA-resistant CRC trajectory-conditioned dependency mapping

## 结论先行

本轮使用 DepMap 23Q4 的 56 个 CRC 模型、6 条 acquired OXA-R trajectory，但按 HCT116、DLD1、HT29、LoVo 四个独立生物学背景聚合。主模型为 top 250+250 weighted directional score；HCT116/DLD1/HT29/LoVo 自身 DepMap cell line 在对应 trajectory 中排除。
最终判定：**B: strong but non-universal dependency identified**。Tier1=0，Tier2=1105，Tier3=0。

## 1. 这六条 trajectory 是否其实高度相似？

状态分数相关性见 `phase7bR_state_score_correlation.csv`。主 weighted top250 中，HCT116 三条 trajectory 的两两相关为：GSE77932|HCT116 vs GSE42387|HCT116: rho=0.629; GSE77932|HCT116 vs GSE119603|HCT116: rho=0.092; GSE42387|HCT116 vs GSE119603|HCT116: rho=-0.032。相关性若不高，不把它解释为失败，而是保留异质轨迹并依赖背景级收敛。
Tier2=1105 是按预注册的宽松发现阈值得到的候选层，不等于 1105 个已验证靶点；Tier2 本身不要求 meta FDR 或每个背景都有 resampling support。因此本轮没有产生可以直接进入药物筛选的 universal target，下一步若继续应先在独立功能/药敏数据中收缩候选。

## 2. HCT116 三条轨迹是否只是重复实验？

不是独立生物学背景；本分析将其作为同一 HCT116 background，并用 median 聚合。它们仍保留在 trajectory-level 表中，用于审计轨迹间一致性和 leave-one-dataset/trajectory 稳健性。

## 3–5. 功能依赖收敛层级

Tier1 universal: 0；Tier2 strong: 1105；Tier3 subtype: 0。主排名前 20：

```text
 rank     gene primary_tier  median_background_vulnerability_rho  n_positive_backgrounds  n_resampling_supported_backgrounds  meta_q_value
    1 C16ORF86 Tier2_strong                             0.319697                       3                                   1      0.231863
    2 HS3ST3B1 Tier2_strong                             0.283983                       3                                   1      0.269461
    3     TJP1 Tier2_strong                             0.275866                       3                                   1      0.064374
    4    FOXF1 Tier2_strong                             0.269769                       3                                   1      0.370157
    5   CFAP90 Tier2_strong                             0.264250                       3                                   0      0.328570
    6  ZFP36L1 Tier2_strong                             0.263384                       3                                   2      0.233174
    7  ADAMTS9 Tier2_strong                             0.263095                       3                                   2      0.159394
    8   DIAPH3 Tier2_strong                             0.257612                       3                                   1      0.269461
    9    SPRY1 Tier2_strong                             0.253319                       4                                   0      0.509125
   10     SOD2 Tier2_strong                             0.252381                       3                                   0      0.231863
   11 ATP6V0E2 Tier2_strong                             0.250974                       3                                   1      0.367239
   12     MOB4 Tier2_strong                             0.250866                       3                                   2      0.064374
   13 SLC25A14 Tier2_strong                             0.248846                       3                                   1      0.334178
   14   PPP6R2 Tier2_strong                             0.248629                       3                                   1      0.231863
   15     MXD4 Tier2_strong                             0.246501                       4                                   0      0.482538
   16  SLC52A3 Tier2_strong                             0.243831                       4                                   0      0.529312
   17     AFDN Tier2_strong                             0.242280                       4                                   2      0.360875
   18   IFNA10 Tier2_strong                             0.241414                       3                                   2      0.328570
   19     JDP2 Tier2_strong                             0.240620                       3                                   1      0.233174
   20    KDM5B Tier2_strong                             0.239899                       3                                   1      0.237826
```

## 6. 预设机制节点审计（不参与排名）

```text
   gene  rank primary_tier  median_background_vulnerability_rho  n_positive_backgrounds
    FN1  1057 Tier2_strong                             0.121681                       3
HSP90B1  1129         none                             0.120310                       2
  ITGB1  1882         none                             0.098160                       3
   GPX4  2237         none                             0.090368                       2
   CPT2  2656         none                             0.081818                       2
  PSMB5  3195         none                             0.072403                       2
  EDEM1  4608         none                             0.050108                       3
SLC22A5  6714         none                             0.024459                       3
  DERL1  7886         none                             0.010606                       2
SLC7A11  7931         none                             0.010101                       2
    VCP 10007         none                            -0.011941                       2
  CPT1A 12279         none                            -0.040224                       1
  DHODH 13392         none                            -0.055339                       1
   RRM2 14476         none                            -0.072872                       2
```

## 7. 稳健性审计

LOBO 结果见 `phase7bR_leave_one_background_out.csv`；HCT116 留出后的正向候选数（median rho > 0）为 9156。LODO 全局 rank Spearman 稳定性：{'GSE119603': 0.9800271835969915, 'GSE42387': 0.7938646719447395, 'GSE77932': 0.7264383024616113}。
signature size 100/250/500 与 rank score 结果见 `phase7bR_signature_size_sensitivity.csv`。协变量敏感性见 `phase7bR_covariate_sensitivity.csv`；adjusted-vs-raw vulnerability rho 的绝对差中位数为 0.041。
置换和 Bootstrap 仅对主 top200/trajectory 候选执行，各 1000 次，结果见 `phase7bR_bootstrap_permutation_results.csv`。

## 研究边界

本轮没有执行 Phase 7C/8，没有把预设机制节点、药物机制或 LINCS 反向签名用于主排名。Meldonium 不作强行挽救；只有后续出现 CPT2/SLC22A5/BBOX1 等功能依赖才可重新进入候选层。原始 DepMap/GEO 文件保持本地，不纳入 GitHub。

## 复现

```powershell
python work/scripts/phase7bR_robust_dependency_mapping.py
```