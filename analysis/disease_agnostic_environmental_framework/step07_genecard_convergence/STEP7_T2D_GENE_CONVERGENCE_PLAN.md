# Step 7 — T2D-specific CTD × GeneCards convergence

## Objective

Test all 11 exposure clusters carried forward from the frozen T2D Step 6 audit against a disease-specific T2D gene set. Step 7 is a post-firewall biological prioritization step; it must not change the 29-test screen, the 14 FDR-positive tests, the 13 robust tests, or the 11-cluster definition.

## Frozen inputs

- `step06_t2d_robustness/t2d_exposure_clusters.csv`
- `step06_t2d_robustness/t2d_robustness_results.csv`
- `hypothesis_unit_audit/step4_test_chemical_membership.csv`
- CTD human chemical–gene interaction export, with unique `ChemicalID × GeneID` pairs
- A complete ranked GeneCards result from the ordinary T2D query `type 2 diabetes mellitus`, retaining rank, symbol, and score fields, as the sole formal gene set. The source may be the official public results table or official export; acquisition mode and checksum must be recorded.

The CRC GeneCards export is explicitly prohibited as a substitute. If the T2D input is absent, the script runs a CTD-side preflight only and records the GeneCards analysis as blocked; it never fabricates a T2D gene set from web snippets or another disease.

## Exposure-unit rule

The unit of analysis is the frozen Step 6 exposure cluster. Chemical IDs and parent/proxy annotations are inherited from Step 4 exactly as recorded. Parent relationships are not inferred from chemical names, and no candidate is removed after seeing T2D gene overlap.

## CTD rule

Restrict to `Organism == Homo sapiens`. Deduplicate interactions by `ChemicalID × GeneID`; retain raw interaction-row count and unique PMID count only as evidence-audit fields. Cluster genes are the union of deduplicated human CTD genes over all mapped chemical IDs in that cluster. The primary background is the union of all human CTD genes represented by the 11 frozen clusters.

## GeneCards rule

The formal gene set is the complete ranked result of the ordinary query `type 2 diabetes mellitus` (currently 20,554 captured rows). K=100, 500, 1,000, and 2,000 are descriptive rank cutoffs; the full list is the primary analysis. The previously audited exact `[Disorders] "type 2 diabetes mellitus"` result (111 rows) and synonym-OR query are deprecated preflight artifacts only; neither is a sensitivity analysis or a formal input. The primary query, acquisition mode, retrieval timestamp, GeneCards build when available, and checksum must be recorded. No CRC-specific GeneCards set, CRC outcome, or T2D outcome P value is used to construct the gene set.

## Enrichment and prioritization

For each cluster and each available K, calculate the one-sided hypergeometric enrichment of the cluster CTD gene set in the T2D GeneCards set using the 11-cluster CTD union as background. Report overlap genes, overlap size, odds ratio, raw P, and BH-FDR across the 11 clusters at each K. Also report a rank-weighted overlap using the supplied GeneCards rank; this is descriptive and does not replace the enrichment test.

Final flagship axes are selected only after all 11 clusters have been analyzed, using a transparent joint table of epidemiologic strength/robustness and T2D-specific biological convergence. No post hoc literature or pathway search is used to rescue a cluster.

## Required outputs

- cluster-to-chemical mapping audit
- cluster CTD evidence and gene table
- T2D GeneCards input audit and checksum
- cluster × GeneCards overlap/enrichment table for the primary set
- overlap-gene table
- analysis manifest and report, including missing-input or blocked status when applicable
