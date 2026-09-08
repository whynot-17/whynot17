# DINP–CRC TCGA/GEO transcriptomic validation

## Scope

The frozen input was the 97-gene `DINP_CRC_overlap.csv` list. The analysis tests tumor-versus-normal expression only; it does not infer DINP causality and does not perform MR, survival, PPI, enrichment, TCGA mutation, or single-cell analysis.

## Statistical design

- GEO paired cohorts and the TCGA matched-patient sensitivity use two-sided paired Wilcoxon signed-rank tests.
- All-sample TCGA and TCGA-versus-GTEx contrasts use two-sided Mann–Whitney U tests.
- BH FDR is computed separately within each contrast across the 97 target genes.
- Effect direction is tumor median minus normal median on the delivered within-dataset expression scale; positive means tumor-high.
- Cross-dataset summary reports direction/sign consistency and does not pool platform-specific expression scales into a causal meta-analysis.

## Cohorts

- `TCGA-COAD+GTEx-colon`: tumor=288, normal=41, GTEx-colon-normal=308, paired=NA; Toil/Xena TcgaTargetGtex_rsem_gene_tpm.
- `GSE74602`: tumor=30, normal=30, GTEx-colon-normal=NA, paired=30; 30 paired colorectal tumor and normal samples; Illumina humanRef-8 v2.0.
- `GSE10950`: tumor=24, normal=24, GTEx-colon-normal=NA, paired=24; 24 paired colon tumor and normal mucosa samples; Illumina humanRef-8 v2.0.
- `GSE156355`: tumor=6, normal=6, GTEx-colon-normal=NA, paired=6; 6 paired colorectal tumor and adjacent normal tissues; Agilent human microarray.

## Contrast-level results

- `GSE10950`: 92/97 genes measured; 66 FDR<0.05; tumor-high=36; normal-high=30.
- `GSE156355`: 94/97 genes measured; 0 FDR<0.05; tumor-high=0; normal-high=0.
- `GSE74602`: 92/97 genes measured; 64 FDR<0.05; tumor-high=21; normal-high=43.
- `TCGA-COAD`: 95/97 genes measured; 77 FDR<0.05; tumor-high=36; normal-high=41.
- `TCGA-COAD_paired`: 95/97 genes measured; 58 FDR<0.05; tumor-high=28; normal-high=30.
- `TCGA-COAD_vs_GTEx-colon`: 95/97 genes measured; 85 FDR<0.05; tumor-high=45; normal-high=40.

## Genes with at least two significant independent primary validations

- `ACACA`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `ACLY`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `ACSL4`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `AGTR1`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `ANGPT2`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `BMP2`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `CCR7`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `CD36`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `CEBPB`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `CNR1`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `COPS5`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `CYP3A4`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `DUSP1`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `EPAS1`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `ESR1`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `FABP1`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `FAS`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `FASN`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `FDPS`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `G6PD`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `GATA2`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `GRN`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `HSPA9`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `IL18`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `IL1A`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `IL1B`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `LDLR`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `LEPR`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.
- `MKI67`: 3 independent primary datasets; consensus=tumor_high; concordance=1.00.
- `NLRP3`: 3 independent primary datasets; consensus=normal_high; concordance=1.00.

## Interpretation boundary

A tumor-versus-normal difference validates CRC disease-state expression association for the target, not DINP exposure response. Directional discordance across cohorts should be retained as heterogeneity rather than collapsed into a single claim. TCGA internal normal and GTEx colon normal are reported separately because normal-reference choice can change the contrast.

## Source links

- UCSC Toil/Xena hub: https://toil.xenahubs.net; expression dataset `TcgaTargetGtex_rsem_gene_tpm`; phenotype dataset `TcgaTargetGTEX_phenotype.txt`.
- GSE74602: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74602
- GSE10950: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10950
- GSE156355: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156355

## Files

- `DINP_CRC_TCGA_GEO_expression_validation.csv`: gene-level contrast results.
- `DINP_CRC_TCGA_GEO_cross_dataset_summary.csv`: cross-contrast directional/significance summary.
- `DINP_CRC_TCGA_GEO_sample_manifest.csv`: sample groups and pairing metadata.
- `DINP_CRC_TCGA_GEO_target_expression_long.csv`: target-gene expression values retained for audit.
- `DINP_CRC_TCGA_GEO_expression_validation_manifest.json`: data and run manifest.
