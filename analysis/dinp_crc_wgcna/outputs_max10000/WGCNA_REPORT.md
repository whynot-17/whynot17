# DINP--CRC bulk WGCNA report

## Scope

This analysis builds an unbiased CRC co-expression network and overlays the frozen 81-gene DINP--CRC intersection and seven macrophage driver candidates. WGCNA is used as exploratory co-expression convergence evidence; it does not establish that DINP causes any module or trait.

## Input audit

- Samples: 585 total; 566 tumor and 19 non-tumoral.
- Probe sets: 54675; collapsed gene-level matrix: 21755 genes.
- WGCNA input: 10000 genes selected by MAD only; frozen target/driver genes were not forced into network construction.
- Frozen 81-gene list present on GPL570: 75/81; selected by MAD: 41; seven-driver list present: 6/7; selected by MAD: 5.
- Macrophage/myeloid proxy: 14/15 fixed markers available; the resulting score is an abundance/state proxy, not cell deconvolution.
- Duplicate probes were collapsed per gene by retaining the probe with the largest across-sample MAD.
- The processed GEO matrix was used as supplied; no outcome-driven gene filtering was applied.

## Network analyses

- All samples: n=585; 22 non-grey modules; soft power=10.
- Tumor-only: n=566; 19 non-grey modules; soft power=10.
- The all-sample analysis enables a descriptive tumor-status module-trait comparison; the tumor-only analysis is the less status-dominated sensitivity network.
- Macrophage-associated module correlations are reported in macrophage_module_association_all_samples.csv and macrophage_module_association_tumor_only.csv.

## Highest target-set enrichment modules

### All samples

- turquoise: 14 target genes; Fisher P=0.0003427; BH-FDR=0.00754
- greenyellow: 3 target genes; Fisher P=0.1389; BH-FDR=    1
- magenta: 3 target genes; Fisher P=0.3257; BH-FDR=    1
- salmon: 1 target genes; Fisher P=0.671; BH-FDR=    1
- black: 2 target genes; Fisher P=0.7019; BH-FDR=    1

### Tumor-only

- turquoise: 13 target genes; Fisher P=0.001368; BH-FDR=0.026
- red: 6 target genes; Fisher P=0.04546; BH-FDR=0.4319
- green: 4 target genes; Fisher P=0.3186; BH-FDR=    1
- cyan: 1 target genes; Fisher P=0.4207; BH-FDR=    1
- salmon: 1 target genes; Fisher P=0.5199; BH-FDR=    1

## Interpretation boundary

Target-set enrichment is calculated against the natural WGCNA input and uses only target genes selected by the prespecified MAD filter; all frozen targets remain visible in the overlay audit. Target/driver module membership or hubness should be treated as prioritization rather than mechanistic proof. Independent bulk or single-cell replication remains necessary.
