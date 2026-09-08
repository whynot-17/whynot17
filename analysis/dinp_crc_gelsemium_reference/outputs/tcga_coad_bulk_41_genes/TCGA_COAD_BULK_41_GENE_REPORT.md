# TCGA-COAD bulk validation of the fresh 41-gene DINP–CRC intersection

## Status

This is a new bulk expression check using only the fresh 41 genes reconstructed from the published Gelsemium elegans–CRC GeneCards supplement and the frozen DINP multi-source matrix. It does not reuse the legacy 81-gene or 9-gene analyses.

- COAD primary tumour samples: **288**
- COAD solid-tissue normal samples: **41**
- Genes queried: **41**
- Independent tumor–normal tests with BH-FDR < 0.05: **31 / 41**
- Positive median shifts (tumor > normal): **16 / 41**
- Patient-matched pairs available: **26**

## Interpretation boundary

The analysis tests whether the 41-gene DINP–CRC intersection is expressed differently in TCGA COAD tumour tissue than in COAD solid-tissue normal tissue. It is a disease-state and tissue-context validation only; it does not establish DINP exposure, direct target binding, causality, or temporal direction.

The primary inferential family is the 41 queried genes, with BH-FDR applied across all 41 independent Mann–Whitney tests. Patient-paired results, when available, are reported separately and are not used to rescue or re-rank the independent analysis.

## Data source

- UCSC Toil Xena hub: `https://toil.xenahubs.net`
- Expression dataset: `TcgaTargetGtex_rsem_gene_tpm`
- Phenotype dataset: `TcgaTargetGTEX_phenotype.txt`
- Disease field: `primary disease or tissue == Colon Adenocarcinoma`
- Groups: `_sample_type == Primary Tumor` versus `_sample_type == Solid Tissue Normal`
- Expression values are analyzed on the Xena-delivered gene-TPM scale; no cross-platform normalization is applied.

## Outputs

- `tcga_coad_sample_manifest.csv`: complete Xena phenotype mapping and inclusion flag.
- `tcga_coad_expression_41_genes.csv`: downloaded expression matrix for the scoped COAD samples.
- `tcga_coad_independent_gene_stats.csv`: independent tumor–normal comparison with 41-gene BH-FDR.
- `tcga_coad_paired_gene_stats.csv`: available patient-matched sensitivity comparison.
- `tcga_coad_bulk_manifest.json`: source, input hash and run metadata.

## Caveat

COAD solid-tissue normal availability is limited in TCGA. If the normal reference is small, the estimate is treated as a precision-limited tissue comparison rather than a definitive universal CRC direction. The other TCGA CRC subtype, GTEx colon, and single-cell analyses are separate contexts and are not silently pooled into this COAD-only run.
