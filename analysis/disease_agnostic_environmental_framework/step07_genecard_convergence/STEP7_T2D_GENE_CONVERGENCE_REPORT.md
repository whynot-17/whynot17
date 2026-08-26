# Step 7 — T2D-specific CTD × GeneCards biological convergence

- Status: **complete_primary_genecards_convergence**
- Frozen exposure clusters: **11**
- CTD human chemical–gene pairs represented: **7,551**
- Cluster CTD gene memberships summed over clusters: **7,379**

## GeneCards input

- Primary ordinary query: `type 2 diabetes mellitus`; complete public result: **20,554 rows**.
- Multiple testing is corrected across the 11 frozen exposure clusters in this single primary family.
- An earlier overly restrictive Disorders-scoped exact-phrase query was deprecated during method auditing and is retained solely for provenance; it does not contribute to the formal Step 7 results.

## Primary full-list result

- Primary q < 0.05 clusters: **5** — cluster_11, cluster_2, cluster_5, cluster_6, cluster_8.
- Minimum primary q: **0.0001214**.
- Step 6 robustness + Step 7 primary convergence: **4 Tier A clusters**.

## Cluster-level primary summary

| Cluster | Biomarker(s) | CTD genes | GeneCards overlap | OR | q |
|---|---|---:|---:|---:|---:|
| cluster_6 | URXUBA;URXUSR | 919 | 754 | 1.47 | 0.0001214 |
| cluster_2 | URXCOP | 104 | 93 | 2.61 | 0.003172 |
| cluster_5 | URXP02 | 65 | 60 | 3.69 | 0.003172 |
| cluster_8 | URXUPB | 3016 | 2357 | 1.18 | 0.009762 |
| cluster_11 | URXUUR | 245 | 204 | 1.54 | 0.01409 |
| cluster_4 | URXMIB | 49 | 42 | 1.83 | 0.1435 |
| cluster_7 | URXUMO | 9 | 9 | inf | 0.1435 |
| cluster_1 | LBXPFHS | 497 | 388 | 1.09 | 0.3267 |
| cluster_9 | URXUSN | 4 | 4 | inf | 0.422 |
| cluster_10 | URXUTU | 3 | 3 | inf | 0.4955 |
| cluster_3 | URXECP;URXMHH;URXMOH | 2468 | 1878 | 0.948 | 0.8152 |

## Interpretation boundary

The GeneCards analysis is post-firewall biological prioritization. It does not modify the 29-test T2D screen, the 14 FDR-positive tests, the 13 robustness-supported tests, or the 11 exposure clusters. CTD chemical–gene associations and GeneCards disease associations are not causal evidence.
