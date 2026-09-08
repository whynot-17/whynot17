# DINP–CRC: DEHP-paper method replication

Template paper: Zhu et al., *Discover Oncology* 2025, "Computational analysis of DEHP's oncogenic role in colorectal cancer".

## Stage 1 — transcriptomic preprocessing

Run from repository root:

```bash
Rscript work/dinp_crc_dehp_replication/01_preprocess_crc_transcriptomes.R
Rscript work/dinp_crc_dehp_replication/02_tcga_deg_and_geo_ready.R
```

### Datasets replicated from the paper

- TCGA-COADREAD: discovery / CRC-vs-normal differential expression.
  - Paper reports 41 normal + 476 tumor.
  - Pipeline uses UCSC Xena `TCGA.COADREAD.sampleMap/HiSeqV2` and labels TCGA sample-type codes 11=normal, 01=primary tumor.
  - Current sample counts are audited rather than silently forced to match the paper.
- GSE32323 / GPL570: 17 paired CRC + adjacent non-cancer tissues; cell-line/5-aza samples excluded.
- GSE21510 / GPL570: 25 normal + 123 CRC samples.
- GEO merge expected from paper: 42 normal + 140 tumor.

### Preprocessing choices

- Both GEO cohorts are GPL570 and deposited as RMA/log2-normalized expression.
- Probe IDs are mapped to HGNC gene symbols with `hgu133plus2.db`.
- Multiple probes per gene are collapsed by median.
- Two merged GEO matrices are retained:
  1. `GEO_merged_uncorrected.tsv.gz` — closest to the paper's described merge.
  2. `GEO_merged_combat.tsv.gz` — batch-corrected version for robust ML.
- TCGA DEGs reproduce the paper threshold: nominal P < 0.05 and |log2FC| > 1.
- An FDR < 0.05 version is also saved for sensitivity analysis.

## Stage 2 — DINP target intersection

Compound manifest: `DINP_compound_manifest.tsv`

- DINP CAS: 28553-12-0
- PubChem CID: 590836
- Target sources to mirror the paper: SwissTargetPrediction + ChEMBL

After DINP targets are exported and standardized to gene symbols, intersect:

`DINP targets ∩ CRC disease targets (GeneCards/OMIM) ∩ TCGA DEGs`

The resulting common genes then enter STRING PPI, followed by LASSO / SVM / Random Forest on the merged GEO matrix, matching the DEHP paper's downstream structure.

## Important replication note

The DEHP paper does not clearly document batch correction for the merged GEO cohorts and uses a nominal P-value DEG threshold. Both the literal-style and more defensible alternatives are therefore retained rather than silently changing the method.
