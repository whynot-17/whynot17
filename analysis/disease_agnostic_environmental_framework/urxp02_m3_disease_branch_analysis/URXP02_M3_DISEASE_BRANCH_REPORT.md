# URXP02 M3 disease-branch analysis

Generated 2026-08-30T07:51:15.205537+00:00. This package analyzes the frozen M1b 828-gene universe against the M2 thyroid and hypertension disease branches.

## Scope and safeguards

- Custom-background ORA was run separately for thyroid-specific (30), hypertension-specific (251), and shared (189) genes.
- Enrichment sources were GO Biological Process, Reactome, and KEGG using g:Profiler g:GOSt with the complete 828-gene expanded universe as background and g:SCS correction.
- STRING functional networks were retrieved for the full 219-gene thyroid branch, 440-gene hypertension branch, and 189-gene shared core at combined score >=700; no added interactors were requested.
- Modules use weighted greedy modularity communities; centrality is descriptive and calculated on the retrieved network.
- Shared hub candidates require exact 2-NAP evidence (human or experimental), at least two evidence sources, and top-10% network centrality by any of degree, betweenness, or PageRank. This is a transparent gate, not a composite score.
- No NHANES model, pathway-to-disease causal claim, sex-specific molecular claim, tissue/cell mapping, or figures were produced.

## Branch evidence composition

- **thyroid-specific**: 30 genes; exact 2-NAP human 0; exact 2-NAP experimental 0; parent naphthalene 30; multi-source (>=2) 0.
- **hypertension-specific**: 251 genes; exact 2-NAP human 26; exact 2-NAP experimental 4; parent naphthalene 244; multi-source (>=2) 11.
- **shared**: 189 genes; exact 2-NAP human 34; exact 2-NAP experimental 2; parent naphthalene 185; multi-source (>=2) 17.

## Enrichment audit

- **thyroid-specific**: status OK; retrieved 1871 terms; g:SCS-significant terms 0; query genes 30; effective background is API-reported per term.
- **hypertension-specific**: status OK; retrieved 6508 terms; g:SCS-significant terms 0; query genes 251; effective background is API-reported per term.
- **shared**: status OK; retrieved 6952 terms; g:SCS-significant terms 148; query genes 189; effective background is API-reported per term.

## PPI/network audit

- **thyroid-branch**: 219 network nodes, 1058 edges, 33 components, largest component 186, 39 modules; STRING PPI-enrichment p=0.0.
- **hypertension-branch**: 440 network nodes, 2360 edges, 52 components, largest component 387, 61 modules; STRING PPI-enrichment p=0.0.
- **shared-core**: 189 network nodes, 978 edges, 23 components, largest component 166, 29 modules; STRING PPI-enrichment p=0.0.

## Shared-core hub candidates

The shared-core evidence/network gate identifies **4** priority candidates: TP53, JUN, CASP3, ESR1.

The complete ranked shared centrality table is in `07_shared_hub_candidates.csv`; it includes genes that fail one or more evidence/centrality gates so that exclusions remain auditable.
