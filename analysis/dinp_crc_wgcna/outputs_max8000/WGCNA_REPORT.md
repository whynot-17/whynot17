# DINP--CRC bulk WGCNA report

## Scope

This analysis builds an unbiased CRC co-expression network and overlays the frozen 81-gene DINP--CRC intersection and seven macrophage driver candidates. WGCNA is used as exploratory co-expression convergence evidence; it does not establish that DINP causes any module or trait.

## Input audit

- Samples: 585 total; 566 tumor and 19 non-tumoral.
- Probe sets: 54675; collapsed gene-level matrix: 21755 genes.
- WGCNA input: 8000 genes selected by MAD only; frozen target/driver genes were not forced into network construction.
- Frozen 81-gene list present on GPL570: 75/81; selected by MAD: 38; seven-driver list present: 6/7; selected by MAD: 5.
- Macrophage/myeloid proxy: 14/15 fixed markers available; the resulting score is an abundance/state proxy, not cell deconvolution.
- Duplicate probes were collapsed per gene by retaining the probe with the largest across-sample MAD.
- The processed GEO matrix was used as supplied; no outcome-driven gene filtering was applied.

## Network analyses

- All samples: n=585; 15 non-grey modules; soft power=10.
- Tumor-only: n=566; 16 non-grey modules; soft power=9.
- The all-sample analysis enables a descriptive tumor-status module-trait comparison; the tumor-only analysis is the less status-dominated sensitivity network.
- Macrophage-associated module correlations are reported in macrophage_module_association_all_samples.csv and macrophage_module_association_tumor_only.csv.

## Highest target-set enrichment modules

### All samples

- blue: 14 target genes; Fisher P=0.0003816; BH-FDR=0.005724
- greenyellow: 3 target genes; Fisher P=0.218; BH-FDR=    1
- tan: 1 target genes; Fisher P=0.3375; BH-FDR=    1
- black: 3 target genes; Fisher P=0.3549; BH-FDR=    1
- purple: 2 target genes; Fisher P=0.4955; BH-FDR=    1

### Tumor-only

- turquoise: 14 target genes; Fisher P=0.0004635; BH-FDR=0.007416
- green: 5 target genes; Fisher P=0.1383; BH-FDR=    1
- red: 3 target genes; Fisher P=0.4979; BH-FDR=    1
- pink: 1 target genes; Fisher P=0.9191; BH-FDR=    1
- yellow: 1 target genes; Fisher P=0.9476; BH-FDR=    1

## Interpretation boundary

Target-set enrichment is calculated against the natural WGCNA input and uses only target genes selected by the prespecified MAD filter; all frozen targets remain visible in the overlay audit. Target/driver module membership or hubness should be treated as prioritization rather than mechanistic proof. Independent bulk or single-cell replication remains necessary.
