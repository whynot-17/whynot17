# DINP–CRC male-versus-female molecular convergence

Generated: `2026-09-02T17:04:11.977698+00:00`

## Frozen interpretation

This analysis compares the frozen 81-gene accessible DINP–CRC intersection between male and female TCGA CRC primary tumors. It is a tumor molecular-state comparison, not evidence that MCOP has a male-specific epidemiologic effect.

## Cohort

- TCGA CRC primary tumors analyzed: **377** (205 male, 172 female).
- Tissue composition: **286 COAD** and **91 READ**.

| Sex | Tumor type | N |
|---|---|---:|
| Female | COAD | 130 |
| Female | READ | 42 |
| Male | COAD | 156 |
| Male | READ | 49 |

## Main result

- Frozen gene family: **81 genes**.
- Adjusted OLS nominal P<0.05: **4 genes**; adjusted-OLS BH-FDR<0.05: **0 genes**.
- Mann–Whitney BH-FDR<0.05: **0 genes**.
- The adjusted OLS effect is male minus female, with COAD/READ included as a tissue covariate and HC3-robust standard errors.

### Top adjusted-OLS results

| Rank | Gene | Male−female beta | 95% CI | P | BH-FDR |
|---:|---|---:|---|---:|---:|
| 1 | `ABCC4` | 0.3305 | 0.0999 to 0.561 | 0.004966 | 0.2093 |
| 2 | `PTGES` | 0.3948 | 0.1181 to 0.6715 | 0.005168 | 0.2093 |
| 3 | `SIRT3` | 0.1182 | 0.01409 to 0.2224 | 0.02608 | 0.7041 |
| 4 | `SPP2` | 0.3616 | 0.004087 to 0.7192 | 0.04744 | 0.806 |
| 5 | `DKK1` | 0.6142 | -0.03454 to 1.263 | 0.06351 | 0.806 |
| 6 | `DDIT4` | 0.1902 | -0.01125 to 0.3916 | 0.06423 | 0.806 |
| 7 | `ATG5` | 0.08615 | -0.006924 to 0.1792 | 0.06965 | 0.806 |
| 8 | `HSPA1A` | 0.2645 | -0.04998 to 0.5789 | 0.09926 | 0.8091 |
| 9 | `PTGS2` | 0.3338 | -0.07602 to 0.7436 | 0.1104 | 0.8091 |
| 10 | `PTGES3` | 0.07833 | -0.01834 to 0.175 | 0.1123 | 0.8091 |

## Boundaries

- Gene-level tests are descriptive molecular-state analyses; they do not prove sex-specific exposure susceptibility.
- The 81-gene family was inherited unchanged from the upstream three-source convergence output.
- No NHANES MCOP×sex interaction statistic is reinterpreted here.
- Expression values use the Xena-delivered scale; no new exposure or disease outcome was used to redefine the 81-gene family.

## Files

- `tcga_crc_sex_gene_results.csv`: all 81 gene-level results.
- `tcga_crc_sex_top_genes.csv`: top 20 adjusted-OLS rows.
- `tcga_crc_81_gene_expression.csv`: expression and phenotype data used for the cohort.
- `tcga_crc_sex_sample_audit.csv` and `tcga_crc_sex_group_counts.csv`: sample QC.
- `sex_convergence_manifest.json`: frozen inputs and provenance.

## Reproducibility

- Script: `analysis/dinp_crc_sex_convergence/run_dinp_crc_sex_convergence.py`
- Xena hub: `https://toil.xenahubs.net`
- Expression dataset: `TcgaTargetGtex_rsem_gene_tpm`
- Phenotype dataset: `TcgaTargetGTEX_phenotype.txt`
- Upstream 81-gene intersection SHA-256: `305265981b0239b2d83265c08716403ca0598c3c2a340d5449c29ff66a092dec`

**Status: completed as a sex-stratified CRC molecular-state audit; no epidemiologic sex-specific MCOP claim is made.**
