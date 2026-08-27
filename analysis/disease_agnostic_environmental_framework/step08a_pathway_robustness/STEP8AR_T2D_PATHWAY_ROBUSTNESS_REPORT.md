# Step 8A-R pathway enrichment robustness audit

Generated (UTC): 2026-08-27T16:57:35.951636+00:00

## Scope and boundary

This is an additive robustness audit of the frozen Step 8A pathway analysis. It does not replace or alter the original g:Profiler results, the 1,647 globally significant terms, the 321 reduced modules, or the 32 compact representatives.

- Frozen Tier-A axes: **4**.
- Frozen all-axis background: **6,076 genes**.
- Null replicates per axis and null type: **1,000**.
- Null types: **gene-size matched** and **annotation-burden matched**.
- Independent annotation snapshots: **GOA human, Reactome GMT, KEGG REST exports**.
- Term filter for the audit: **5–5000 genes after restriction to the frozen background**.
- ORA is one-sided hypergeometric with BH correction within each axis/source audit family; empirical null P values are descriptive robustness metrics, not new discovery claims.

## Observed axis-level audit

| Axis | Query genes | Reactome significant terms | Concordant themes (>=2 sources) | Reactome xenobiotic/CYP null P |
|---|---:|---:|---:|---:|
| cluster_11 | 245 | 17 | 0 | 1 |
| cluster_5 | 65 | 53 | 4 | 0.000999 |
| cluster_6 | 919 | 114 | 3 | 0.000999 |
| cluster_8 | 3016 | 43 | 1 | 1 |

## Interpretation

A pathway theme is treated as more credible when it recurs across independent annotation resources, exceeds matched-null expectations, and remains interpretable after removal of the most recurrent driver genes. This audit is not a causal exposure-to-pathway test and does not establish pathway activation or mediation of T2D.

## Key robustness findings

- **cluster_5 xenobiotic/CYP:** the theme recurred in GO:BP, Reactome, and KEGG. In Reactome, the observed xenobiotic/CYP term count was 3; the empirical P was **0.000999** under both gene-size matching and annotation-burden matching. The observed Reactome count of 53 significant terms was also beyond the 1,000 gene-size-matched null replicates (empirical P **0.000999**). Annotation-burden matching produced a larger null median for total significant terms (28), so the total term count is not interpreted in isolation; the theme-specific result remained extreme.
- **cluster_5 driver audit:** the top recurrence-ranked genes were TP53, CCND1, BAX, FOS, and JUN. After removing all five, the xenobiotic/CYP theme remained represented by 9 GO:BP, 3 Reactome, and 3 KEGG significant terms, supporting persistence beyond a single recurring driver set.
- **cluster_6:** the large raw term count was partly sensitive to annotation burden (for example, Reactome annotation-matched null median 17 versus 114 observed significant terms). It is therefore treated as a broad, less specific convergence pattern rather than as evidence proportional to the number of significant terms.
- **cluster_8 and cluster_11:** these axes showed resource-specific pathway structure rather than a cross-resource xenobiotic/CYP pattern; they remain supporting axes and are not promoted by this audit.

The empirical P values above are permutation-based descriptive robustness metrics with a lower resolution of 1/(1,000+1); they are not additional BH-corrected disease-discovery claims.

## Driver analysis

Driver recurrence counts how often each frozen query gene appears in significant term intersections across the three annotation resources. Leave-driver-out results remove the top 1, 3, or 5 recurrence-ranked genes and re-run the same ORA engine.
