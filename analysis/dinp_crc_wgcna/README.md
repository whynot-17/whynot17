# DINP--CRC bulk WGCNA

This analysis places the frozen 81-gene DINP--CRC intersection and the seven
macrophage driver candidates in CRC bulk expression data.

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
3. selection of the 5,000 highest-MAD genes while forcing available frozen
   target/driver genes into the input;
4. signed, bicor WGCNA with a transparent soft-threshold rule, minimum module
   size 30, and eigengene merge cut height 0.25;
5. an all-sample network (tumor-status trait included) and a tumor-only network
   sensitivity analysis;
6. module-trait correlations, BH-FDR by trait, target-set enrichment, and
   target/module membership tables.

## Interpretation boundary

WGCNA is exploratory co-expression convergence evidence. It does not establish
that DINP causes a module, does not identify a causal target, and does not
replace independent transcriptomic or experimental validation. The grey bin is
unassigned and is excluded from the module-level enrichment multiplicity
family.

The expression source is the public NCBI GEO GSE39582 record:
<https://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582+>
