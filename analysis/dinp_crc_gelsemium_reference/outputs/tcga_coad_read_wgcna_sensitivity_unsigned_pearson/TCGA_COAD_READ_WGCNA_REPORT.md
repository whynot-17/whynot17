# TCGA-COAD + TCGA-READ full-transcriptome WGCNA

## Scope

A unsigned, pearson WGCNA was built from the top 8000 genes selected by full-transcriptome MAD after gene-symbol deduplication. The fresh 41-gene DINP–CRC set was overlaid after network construction; target genes were not used to select or force network genes.

- Samples: **380** primary tumors (288 COAD, 92 READ).
- Network genes: **8000**; fresh 41-gene targets present as network genes: **16/41**.
- Targets projected to module eigengenes only: **25/41**.
- Non-grey modules: **4**; target-enriched modules by module-level BH-FDR<0.05: **0**.
- Soft-threshold power: **2** (first power with unsigned R2 >= 0.80 and mean connectivity > 1).

## Interpretation boundary

WGCNA describes co-expression structure in CRC bulk tissue. It does not establish DINP exposure, direct chemical binding, causality, or direction of regulation. Module projections for targets absent from the variance-filtered network are descriptive and are not counted in Fisher enrichment. Bulk modules may also reflect tissue composition; macrophage claims require independent cell-type support.

## Main outputs

- `tcga_coad_read_wgcna_gene_module_membership.csv`: module labels and eigengene correlations for all network genes.
- `tcga_coad_read_wgcna_41_gene_overlay.csv`: exact target overlay and projection-only results.
- `tcga_coad_read_wgcna_target_module_enrichment.csv`: correct Fisher ORA using the observed network universe.
- `tcga_coad_read_wgcna_module_trait_correlations.csv`: module eigengene associations with the 41-gene score and cohort.
- `tcga_coad_read_wgcna_module_summary.csv`: module sizes, target overlap, and trait associations.
- `tcga_coad_read_wgcna_sample_dendrogram.pdf`, `tcga_coad_read_wgcna_module_dendrogram.pdf`, and `tcga_coad_read_wgcna_module_trait_heatmap.pdf`: QC/overview figures.

## Reproducibility

- Expression scale: Xena-delivered log2(TPM+0.001).
- Network: unsigned pearson; power=2; minModuleSize=30; mergeCutHeight=0.25; deepSplit=2.
- R: R version 4.6.1 (2026-06-24 ucrt); WGCNA: 1.74; dynamicTreeCut: 1.63.1; fastcluster: 1.3.0.
- No automatic sample removal was performed; sample clustering is saved for inspection.
- The 41-gene set was not used for network-gene selection.
