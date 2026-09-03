# DINP--CRC bulk WGCNA report

## Scope

This analysis places the frozen 81-gene DINP--CRC intersection and the seven macrophage driver candidates in the public GSE39582 GPL570 CRC bulk expression data. WGCNA is used as exploratory co-expression convergence evidence; it does not establish that DINP causes any module or trait.

## Input audit

- Samples: 585 total; 566 tumor and 19 non-tumoral.
- Probe sets: 54675; collapsed gene-level matrix: 21755 genes.
- WGCNA input: 5000 genes selected by MAD, with all available frozen target/driver genes forced into the input.
- Frozen 81-gene list present on GPL570: 75/81; seven-driver list present: 6/7.
- Duplicate probes were collapsed per gene by retaining the probe with the largest across-sample MAD.
- The processed GEO matrix was used as supplied; no outcome-driven gene filtering was applied.

## Network analyses

- All samples: n=585; 9 non-grey modules; soft power=9.
- Tumor-only: n=566; 9 non-grey modules; soft power=9.
- The all-sample analysis enables a descriptive tumor-status module-trait comparison; the tumor-only analysis is the less status-dominated sensitivity network.

## Highest target-set enrichment modules

### All samples

- turquoise: 19 target genes; Fisher P=0.05463; BH-FDR=0.4916
- brown: 7 target genes; Fisher P=0.7437; BH-FDR=    1
- yellow: 4 target genes; Fisher P=0.9195; BH-FDR=    1
- red: 2 target genes; Fisher P=0.9905; BH-FDR=    1
- green: 2 target genes; Fisher P=0.9933; BH-FDR=    1

### Tumor-only

- turquoise: 19 target genes; Fisher P=0.05304; BH-FDR=0.4773
- green: 5 target genes; Fisher P=0.8572; BH-FDR=    1
- black: 3 target genes; Fisher P=0.9361; BH-FDR=    1
- yellow: 4 target genes; Fisher P=0.9445; BH-FDR=    1
- red: 2 target genes; Fisher P=0.9915; BH-FDR=    1

## Interpretation boundary

A target gene set being enriched in a CRC co-expression module supports transcriptomic/network convergence in this dataset. It is not a causal exposure-to-gene test, and module membership/hubness should be treated as prioritization rather than mechanistic proof. Independent bulk or single-cell replication remains necessary.
