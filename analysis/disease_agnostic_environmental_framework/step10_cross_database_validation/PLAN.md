# Step 10 — Cross-database robustness and annotation-bias stress test

## Purpose

The main Step 10 question is whether candidate evidence remains structured when
the knowledge source is changed. It is **not** a request to repeat the same
candidate across more diseases. The legacy randomized multi-disease analysis is
retained as a supplementary transportability demonstration in
`step10_randomized_replication_expanded/`.

The primary estimand is source sensitivity:

> Does a candidate or exposure axis retain comparable evidence, rank, and
> biological direction under prespecified replacement or perturbation of the
> underlying knowledge source?

This step cannot establish causality, temporal ordering, or independence of
databases.

## Frozen boundaries

1. The environmental actionability universe and the post-firewall candidate
   labels are read-only inputs. No Step 5–9 result, FDR family, robustness tier,
   or candidate rank may be changed by this step.
2. The CTD annotation-burden audit is descriptive. CTD chemical–gene row count,
   unique human chemical–gene pairs, PMIDs, interaction actions, and vocabulary
   fields are not treated as effect sizes or as selection gates.
3. T2D/CRC outcome labels may be carried as provenance metadata for an audit,
   but may not be used to choose a replacement database, tune a threshold, or
   re-rank candidates.
4. Each source must have an eligibility and entity-resolution rule frozen before
   its results are inspected. Missing evidence is not equivalent to negative
   evidence.
5. Sources with different evidence semantics are not pooled into one edge list.
   For example, a CTD expression interaction, an experimental binding record,
   and a LINCS perturbational response are separate evidence types and require
   source-specific summaries.

## Analysis layers

### A. Chemical-to-gene source replacement

CTD remains the reference source for the current audit. Potential independent
sources are evaluated separately, subject to public availability and compatible
human gene identifiers:

- ChEMBL and BindingDB: experimentally observed compound–protein binding or
  activity evidence, not interchangeable with CTD expression interactions.
- LINCS/L1000: perturbational expression response; summarize gene-direction
  signatures rather than pretending they are direct target edges.
- Other toxicogenomics or curated toxicology resources may be added only after
  their evidence semantics, release/version, species scope, and identifier
  mapping are recorded.

### B. Disease-gene source replacement

GeneCards remains the source used in the post-firewall T2D convergence analysis.
Independent disease-gene sources are tested as separate sensitivity layers:
DisGeNET, Open Targets, GWAS Catalog, OMIM, and ClinGen where the disease
definition and gene-level evidence are sufficiently explicit. The ordinary
GeneCards T2D result is the canonical input for the existing analysis; the
historical 111-gene scoped preflight is deprecated and excluded.

### C. Mechanistic corroboration

Mechanistic evidence is summarized independently through transcriptomic
perturbation, human tissue expression, CRISPR dependency, proteomics, or other
well-defined resources. These are corroboration layers, not interchangeable
replacements for chemical–gene edges.

## Annotation-burden audit

The first executed component compares the 134 post-firewall T2D candidate
chemicals with the 2,042-chemical actionability universe and its 409 actionable
chemicals. CTD interactions are restricted to human records and deduplicated at
chemical ID × GeneID (GeneSymbol fallback only when GeneID is absent). Raw rows,
unique PMIDs, interaction actions, and vocabulary annotations are reported
separately.

The audit reports global and actionable-background percentiles, within-class
percentiles when class labels are available, and descriptive candidate-versus-
noncandidate comparisons. A top-decile burden flag is an audit descriptor only;
it is not a candidate filter.

## Cross-database comparison metrics

For every source, before examining results, define the comparable unit and
report:

- entity-resolution yield and source-specific eligible record counts;
- top-k retention and overlap/Jaccard at fixed k;
- Spearman/Kendall rank concordance where a source provides a rankable score;
- rank-biased overlap or weighted overlap for long-tailed rankings;
- gene/pathway direction concordance where direction is defined by the source;
- missingness and unresolvable-identifier rates;
- source release, download date, species, evidence scope, and full input hash.

No single composite score is primary. If a composite summary is later useful,
its weights and missing-data handling must be frozen before calculation and
reported alongside all component metrics.

## Null and stress tests

The following are planned only after the relevant source inputs are complete:

- chemical-class- and gene-set-size-matched nulls;
- CTD annotation-degree/PMID-burden-matched nulls;
- degree-preserving network perturbation where a network is available;
- identifier-remapping and source-dropout sensitivity;
- leave-one-source-out consensus summaries.

The null tests answer whether observed cross-source retention exceeds what is
expected from annotation density, set size, or network degree. They do not turn
database concordance into causal evidence.

## Gating and outputs

The main Step 10 gate is passed only when at least two independent replacement
sources have complete provenance, predeclared eligibility rules, and a valid
comparison on the same candidate unit. Until then, the status is
`audit_complete_external_validation_pending`.

Current executed outputs are the CTD burden tables, summary, manifest, and audit
report in this directory. The old multi-disease Step 10R outputs remain
unchanged and are supplementary rather than the primary external-validation
claim.
