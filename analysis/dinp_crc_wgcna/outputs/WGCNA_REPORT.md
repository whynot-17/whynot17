# DINP--CRC bulk WGCNA report

## Scope

This analysis builds an unbiased CRC co-expression network and overlays the frozen 81-gene DINP--CRC intersection and seven macrophage driver candidates. WGCNA is used as exploratory co-expression convergence evidence; it does not establish that DINP causes any module or trait.

## Input audit

- Samples: 585 total; 566 tumor and 19 non-tumoral.
- Probe sets: 54675; collapsed gene-level matrix: 21755 genes.
- WGCNA input: 5000 genes selected by MAD only; frozen target/driver genes were not forced into network construction.
- Frozen 81-gene list present on GPL570: 75/81; selected by MAD: 33; seven-driver list present: 6/7; selected by MAD: 4.
- Macrophage/myeloid proxy: 14/15 fixed markers available; the resulting score is an abundance/state proxy, not cell deconvolution.
- Duplicate probes were collapsed per gene by retaining the probe with the largest across-sample MAD.
- The processed GEO matrix was used as supplied; no outcome-driven gene filtering was applied.

## Network analyses

- All samples: n=585; 8 non-grey modules; soft power=8.
- Tumor-only: n=566; 10 non-grey modules; soft power=8.
- The all-sample analysis enables a descriptive tumor-status module-trait comparison; the tumor-only analysis is the less status-dominated sensitivity network.
- Macrophage-associated module correlations are reported in macrophage_module_association_all_samples.csv and macrophage_module_association_tumor_only.csv.

## Highest target-set enrichment modules

### All samples

- turquoise: 14 target genes; Fisher P=0.0006273; BH-FDR=0.005018
- green: 8 target genes; Fisher P=0.02146; BH-FDR=0.08584
- red: 3 target genes; Fisher P=0.6296; BH-FDR=    1
- yellow: 1 target genes; Fisher P=0.9792; BH-FDR=    1
- black: 0 target genes; Fisher P=    1; BH-FDR=    1

### Tumor-only

- turquoise: 13 target genes; Fisher P=0.002023; BH-FDR=0.02023
- green: 5 target genes; Fisher P=0.1512; BH-FDR=0.756
- blue: 5 target genes; Fisher P=0.2786; BH-FDR=0.9286
- magenta: 1 target genes; Fisher P=0.8281; BH-FDR=    1
- pink: 1 target genes; Fisher P=0.879; BH-FDR=    1

## Interpretation boundary

Target-set enrichment is calculated against the natural WGCNA input and uses only target genes selected by the prespecified MAD filter; all frozen targets remain visible in the overlay audit. Target/driver module membership or hubness should be treated as prioritization rather than mechanistic proof. Independent bulk or single-cell replication remains necessary.
