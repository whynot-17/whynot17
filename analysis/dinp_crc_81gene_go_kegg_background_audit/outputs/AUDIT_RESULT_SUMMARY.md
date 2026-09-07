# DINP–CRC 81-gene background audit: observed result

The run completed against four statistical universes using the same 81-gene
query. The primary DINP exposure universe contains 86 genes; the CRC disease
union contains 15,885 genes; the legacy exposure∪CRC union contains 15,890
genes. The query is therefore 81/86 (94.2%) of the DINP universe.

## Global result

The primary DINP-background analysis returned 11 terms with global
BH_FDR_all_GO_KEGG < 0.05, mostly broad ontology terms such as molecular
function, binding, biological regulation, and response to stimulus. The
corresponding counts were 534 terms for the CRC union, 534 for the legacy
union, and 634 for the annotated-genome descriptive background.

None of the requested pathway families passed the primary global FDR filter.
The selected terms below are representative values from
`key_pathway_background_comparison.csv`:

| Term | DINP 86-gene raw p / global FDR | CRC union raw p / global FDR | Legacy union raw p / global FDR | Genome raw p / global FDR |
|---|---:|---:|---:|---:|
| Arachidonic acid metabolism | 0.316 / 1.000 | 4.08e-10 / 4.34e-7 | 4.06e-10 / 4.31e-7 | 1.03e-8 / 3.93e-6 |
| Prostaglandin metabolic process | 0.232 / 1.000 | 5.10e-15 / 2.50e-11 | 5.06e-15 / 2.49e-11 | 1.48e-15 / 4.83e-12 |
| Prostaglandin biosynthetic process | 0.367 / 1.000 | 1.85e-11 / 2.80e-8 | 1.84e-11 / 2.79e-8 | 9.11e-12 / 1.02e-8 |
| Prostaglandin E receptor activity | 0.493 / 1.000 | 2.67e-12 / 5.03e-9 | 2.66e-12 / 4.94e-9 | 7.52e-13 / 1.14e-9 |
| Inflammatory response | 0.0175 / 1.000 | 9.61e-14 / 2.70e-10 | 9.46e-14 / 2.66e-10 | 6.41e-16 / 2.52e-12 |
| PPAR signaling pathway | 0.493 / 1.000 | 3.08e-5 / 4.28e-3 | 3.07e-5 / 4.26e-3 | 1.94e-4 / 1.63e-2 |

## Interpretation

Under the prespecified rule, these pathway signals should not be described as
DINP–CRC convergence-specific enrichment. They are features of the 81-gene
set and/or the DINP-related biology that become highly significant when the
null universe is expanded to the CRC or genome scale.

This audit does not invalidate the PGE2–PTGER4 expression analysis, but it
removes the enrichment result as independent support for a CRC-specific
convergence mechanism. The PGE2–PTGER4 axis should therefore be framed as a
separate expression-level hypothesis requiring donor-level and external
cohort validation.

The custom-background analysis follows g:Profiler's documented custom-domain
workflow; the API reports an effective annotated domain after identifier
mapping, while the supplied primary list contains 86 DINP symbols.
