# DINP–CRC sex-stratified molecular convergence

This follow-up tests whether the frozen 81-gene accessible DINP–CRC
intersection has different expression states in male versus female TCGA CRC
primary tumors.

The analysis uses the public UCSC Toil Xena harmonized
`TcgaTargetGtex_rsem_gene_tpm` expression dataset and its phenotype table.
The primary comparison is male versus female among TCGA colon and rectal
adenocarcinoma primary tumors. A covariate-adjusted OLS sensitivity includes
COAD/READ disease type; expression-level results are not interpreted as a
sex-specific epidemiologic effect of MCOP.

Run from the repository root:

```powershell
C:\Users\21634\anaconda3\python.exe analysis/dinp_crc_sex_convergence/run_dinp_crc_sex_convergence.py
```

The 81-gene input is read from the upstream three-source convergence output;
no upstream intersection, source evidence, or exposure result is modified.
