# Step 7 — T2D-specific CTD × GeneCards biological convergence

- Status: **complete_two_gene_sets**
- Frozen exposure clusters: **11**
- CTD human chemical–gene pairs represented: **7,551**
- Cluster CTD gene memberships summed over clusters: **7,379**

## GeneCards input sets

- Primary ordinary query: `type 2 diabetes mellitus`; complete public result: **20,554 rows**.
- High-specificity sensitivity: `[Disorders] "type 2 diabetes mellitus"`; complete public result: **111 rows**.
- Primary analysis uses the complete ordinary-query list; the 111-row exact Disorders result is not used as the primary set.
- Primary full-list and strict-set enrichment are each corrected across the same 11-cluster family.

## Primary full-list result

- Primary q < 0.05 clusters: **5** — cluster_11, cluster_2, cluster_5, cluster_6, cluster_8.
- Minimum primary q: **0.0001214**.

## Strict exact-phrase sensitivity

- Strict q < 0.05 clusters: **0** — none.
- Minimum strict q: **0.1514**.

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
