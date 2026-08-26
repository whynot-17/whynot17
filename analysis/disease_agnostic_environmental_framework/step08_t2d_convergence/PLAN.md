# Step 8 — T2D pathway convergence

## Scope

Step 8 begins after the locked Step 7 result. It analyzes only the four Tier A
exposure axes carried forward from Step 7:

```text
cluster_6, cluster_5, cluster_8, cluster_11
```

Step 5 epidemiology, Step 6 robustness labels, the 11-cluster definition, and
the Step 7 primary GeneCards input are frozen. No CRC/MCOP mechanism output is
an input to this T2D analysis.

## Primary pathway input and background

For each Tier A cluster, the query is its frozen Homo sapiens CTD gene union
from `step07_genecard_convergence/t2d_cluster_ctd_gene_membership.csv`.
The statistical background is the union of all 11 frozen cluster CTD gene
sets, not the whole human genome. CTD associations are descriptive exposure-
gene evidence; they are not treated as causal or directional expression data.

## Primary ORA

Run one-sided over-representation analysis separately for each Tier A axis
against these prespecified sources:

- Gene Ontology Biological Process (`GO:BP`)
- Reactome (`REAC`)
- KEGG (`KEGG`)

Queries use the g:Profiler official API with `domain_scope=custom`, the frozen
11-cluster background, human symbols, and `all_results=true`. The complete raw
responses are retained with the g:Profiler version, database date, query,
background, and API options.

The raw hypergeometric P value is reconstructed from the returned effective
domain size, term size, query size, and intersection size. A single global
Benjamini–Hochberg correction is then applied across all returned
axis × source × term tests. The g:Profiler source-level adjusted value is
retained as a secondary audit field, not used to replace the global correction.

## Interpretation boundary

Pathway ORA is direction-agnostic and does not establish pathway activation,
exposure causality, or mediation of T2D. Terms are interpreted only with their
overlap size, effective background, leading/intersection genes, and redundancy
context. Very broad or tiny terms are not promoted solely because of a small P
value.

## Exit criteria for the pathway stage

- All 4 Tier A axes queried successfully;
- query and background ID mapping audited;
- GO:BP, Reactome, and KEGG database metadata captured;
- full results retained, with raw P and global BH-FDR;
- no pathway is promoted without adequate overlap and an explicit axis-level
  evidence table.

Transcriptomic directionality and network context are separate follow-up
stages. They require independent expression or interaction inputs and must not
be inferred from this ORA alone.

## Stage 2 — redundancy reduction and module summary

The 1,647 terms passing the frozen global BH-FDR threshold are retained in
full.  A separate descriptive layer assigns terms to modules using the
captured g:Profiler parent/ancestor structure.  Terms covering more than 25%
of the effective background are not allowed to bridge modules, and a guarded
ancestor-set similarity threshold is used only when a direct significant
parent is absent.  The module representative is selected deterministically
from overlap size, global q value, and term specificity; it is not a new
statistical test.

For the compact axis summary, at most eight eligible representatives are
retained per Tier A axis.  Representatives require at least three overlapping
genes and a term size no greater than 25% of the effective background.  All
other significant terms and all module assignments remain available for
audit.  Cross-source similarity in the compact summary is lexical only and is
not treated as semantic or mechanistic evidence.
