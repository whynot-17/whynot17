# Step 6 T2D robustness and exposure-cluster audit

Generated (UTC): 2026-08-26T11:41:29.965166+00:00

## Scope and firewall

- Frozen primary family: **29 tests**; no FDR was recomputed or narrowed.
- FDR-positive tests audited uniformly: **14**.
- Primary reproductions within absolute log-OR difference <=1e-8: **14/14**.
- Exposure clusters under |cycle-adjusted Spearman rho| >= 0.7 and pairwise N >= 500: **11**.
- Deterministic robust-FDR candidates under the locked rubric: **13**.
- Robust-FDR candidate list: **URXUUR, URXUMO, URXUTU, URXUPB, URXMIB, URXUBA, URXP02, URXUSN, URXMHH, URXMOH, URXUSR, URXECP, LBXPFHS**.
- FDR-positive test(s) with a stability downgrade: **URXCOP**.
- Heterogeneity concerns (H0/H1): **URXUTU Pinteraction=0.044226; URXECP Pinteraction=0.0613155**.
- High-correlation edges: **URXECP-URXMHH |rho|=0.933; URXECP-URXMOH |rho|=0.942; URXMHH-URXMOH |rho|=0.979; URXUBA-URXUSR |rho|=0.774**.
- No GeneCards, disease-specific CTD, transcriptomics, literature, or mechanistic analysis was performed.

## Candidate audit

| Biomarker | OR | 95% CI | P | q (29) | LOCO same | Cycle same | Pinteraction | Tail max | Cluster | Priority |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|
| URXUUR | 1.165 | 1.099-1.234 | 7.127e-07 | 2.067e-05 | 8/8 | 7/8 | Pinteraction=0.209331 | 0.0175 | cluster_11 | robust_fdr_candidate |
| URXUMO | 1.19 | 1.107-1.279 | 4.725e-06 | 6.851e-05 | 10/10 | 10/10 | Pinteraction=0.867852 | 0.02271 | cluster_7 | robust_fdr_candidate |
| URXUTU | 1.135 | 1.069-1.206 | 5.146e-05 | 0.0004974 | 10/10 | 9/10 | Pinteraction=0.044226 | 0.03321 | cluster_10 | robust_fdr_candidate |
| URXUPB | 0.8496 | 0.777-0.929 | 0.0004238 | 0.003073 | 10/10 | 8/10 | Pinteraction=0.806086 | 0.007843 | cluster_8 | robust_fdr_candidate |
| URXMIB | 1.096 | 1.036-1.159 | 0.001602 | 0.009289 | 9/9 | 9/9 | Pinteraction=0.623442 | 0.01472 | cluster_4 | robust_fdr_candidate |
| URXUBA | 0.9229 | 0.8756-0.9727 | 0.003016 | 0.01458 | 10/10 | 8/10 | Pinteraction=0.214931 | 0.006585 | cluster_6 | robust_fdr_candidate |
| URXP02 | 1.085 | 1.025-1.148 | 0.005067 | 0.02099 | 8/8 | 6/8 | Pinteraction=0.162626 | 0.02294 | cluster_5 | robust_fdr_candidate |
| URXUSN | 1.125 | 1.035-1.223 | 0.006557 | 0.02377 | 4/4 | 4/4 | Pinteraction=0.372659 | 0.007578 | cluster_9 | robust_fdr_candidate |
| URXMHH | 1.054 | 1.013-1.098 | 0.01032 | 0.03324 | 9/9 | 8/9 | Pinteraction=0.274796 | 0.03948 | cluster_3 | robust_fdr_candidate |
| URXCOP | 1.051 | 1.009-1.093 | 0.01604 | 0.04587 | 7/7 | 4/7 | Pinteraction=0.613251 | 0.03381 | cluster_2 | fdr_candidate_with_instability |
| URXMOH | 1.052 | 1.009-1.097 | 0.0174 | 0.04587 | 9/9 | 8/9 | Pinteraction=0.245085 | 0.029 | cluster_3 | robust_fdr_candidate |
| URXUSR | 0.8865 | 0.8021-0.9798 | 0.01935 | 0.04653 | 3/3 | 3/3 | Pinteraction=0.17195 | 0.0348 | cluster_6 | robust_fdr_candidate |
| URXECP | 1.054 | 1.008-1.103 | 0.02205 | 0.04653 | 8/8 | 6/8 | Pinteraction=0.0613155 | 0.02346 | cluster_3 | robust_fdr_candidate |
| LBXPFHS | 0.9338 | 0.8806-0.9902 | 0.02246 | 0.04653 | 8/8 | 6/8 | Pinteraction=0.634027 | 0.007916 | cluster_1 | robust_fdr_candidate |

## Interpretation

The audit distinguishes primary FDR support from stability and exposure dependence. A high-correlation cluster does not invalidate a primary test and was not used to alter multiplicity. Conversely, a robust-FDR label is a prioritization aid and does not establish temporality or causality.
