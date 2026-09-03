# DINP--CRC bulk WGCNA

This analysis builds an unbiased CRC bulk co-expression network and then
overlays the frozen 81-gene DINP--CRC intersection and the seven macrophage
driver candidates.

## Input

- Public processed GEO series matrix: GSE39582, platform GPL570.
- Platform annotation: GPL570 annotation table.
- Frozen target list: `analysis/dinp_crc_81gene_singlecell_localization/outputs/input_81_genes.csv`.
- Macrophage driver list: `analysis/dinp_crc_81gene_macrophage_driver_decomposition/outputs/macrophage_driver_candidates.csv`.
- The raw GEO matrix and annotation are stored on D: under
  `D:/whynot17/work/wgcna_crc/raw/` and are intentionally not committed to the
  repository.

## Analysis

`run_dinp_crc_wgcna.R` performs:

1. GEO metadata and expression-table parsing;
2. GPL570 probe-to-gene mapping, with duplicate probes collapsed by largest
   across-sample MAD;
3. selection of the 5,000 highest-MAD genes by a prespecified variance filter;
   frozen target/driver genes are not forced into the network;
4. signed, bicor WGCNA with a transparent soft-threshold rule, minimum module
   size 30, and eigengene merge cut height 0.25;
5. an all-sample network (tumor-status trait included) and a tumor-only network
   sensitivity analysis;
6. module-trait correlations, BH-FDR by trait, a fixed macrophage/myeloid
   marker-score association, and target/driver overlay tables;
7. target-set enrichment against the natural WGCNA input, with module-level
   BH-FDR calculated jointly across non-grey modules.

The frozen target and driver sets are overlaid after network construction. The
target enrichment test therefore uses only target/driver genes that happen to
be selected by the natural MAD filter; genes absent from that input remain in
the overlay audit with missing module membership rather than being silently
discarded.

The macrophage score is the within-subset mean z-score of the fixed marker set
`CD68, LST1, TYROBP, AIF1, FCER1G, CTSS, C1QA, C1QB, C1QC, MS4A7, LILRB1,
LGALS3, CD14, CTSB, SPI1`. It is a bulk myeloid abundance/state proxy and is
not a substitute for cell deconvolution.

The primary 5,000-gene run is accompanied by natural-MAD sensitivity runs at
8,000 and 10,000 genes in `outputs_max8000/` and `outputs_max10000/`.
`WGCNA_VARIANCE_FILTER_AUDIT.md` summarizes the prespecified cross-threshold
comparison. Module colors are local labels and are not compared as if they
were stable identities across independently built networks.

## Interpretation boundary

WGCNA is exploratory co-expression convergence evidence. It does not establish
that DINP causes a module, does not identify a causal target, and does not
replace independent transcriptomic or experimental validation. The grey bin is
unassigned and is excluded from the module-level enrichment multiplicity
family.

The expression source is the public NCBI GEO GSE39582 record:
<https://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582+>
