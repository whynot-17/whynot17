# Phase 2G local-H5AD completion audit

Generated: 2026-08-24T11:18:36.631277+00:00

**Completion QC: PASS**

## Provenance and integrity

- Official source H5AD bytes: **3,759,324,000**; exact expected size: **True**.
- C106 9-gene raw-expression equivalence: cells **3,392/3,392**, nnz **3,628/3,628**, count sum **5230/5230**.
- Successful online-cache donor manifests with exact source-H5AD cell counts: **25/25**.
- Eligible epithelial cells: **43,801**; pseudobulk cell sum: **43,801**; conservation: **True**.
- Frozen target universe present: **3,817/3,826**; missing: `CCL4L1, CHLSN, DLEU1, GPX1, H4C6, HSP90AA2P, MTRNR2L8, PRSS2, PTTG3P`.

## Frozen PPAR/NR reproduction

- Phase 2G median paired tumor-normal delta: **-0.418601**; P=**4.29e-07**; n=**36**.
- Phase 2F frozen delta: **-0.418601**; absolute difference: **5.55e-17**; P exact: **True**.
- Cell-level PPAR quartiles define state labels only; donor-level PPAR inference reuses the Phase 2F score standardized in the full epithelial dataset context.

## New state and regulator results

Top paired tumor-normal state contrasts:

| feature | median_delta_tumor_minus_normal | p_value | BH_FDR | direction |
| --- | --- | --- | --- | --- |
| enterocyte_differentiation | -0.64073 | 5.9226e-08 | 1.5991e-06 | down |
| intestinal_epithelial_differentiation | -0.54457 | 3.7692e-07 | 5.0885e-06 | down |
| UPR | 0.43952 | 1.2978e-06 | 1.168e-05 | up |
| stress_like_epithelial | 0.40594 | 1.936e-05 | 0.00012565 | up |
| MYC_targets_V1 | 0.49468 | 2.3269e-05 | 0.00012565 | up |
| OXPHOS | -0.53136 | 3.9719e-05 | 0.00017874 | down |
| Fatty_acid_metabolism | -0.2325 | 0.00018572 | 0.00062681 | down |
| stemness | 0.41547 | 0.00018572 | 0.00062681 | up |

Top donor-level PPAR/state associations:

| group | state | n_donors | spearman_rho | P | BH_FDR |
| --- | --- | --- | --- | --- | --- |
| pooled | NR1I3 | 72 | 0.8231 | 7.2155e-19 | 3.1027e-17 |
| pooled | enterocyte_differentiation | 72 | 0.78819 | 2.0916e-16 | 4.497e-15 |
| normal | NR1I3 | 36 | 0.9021 | 5.8259e-14 | 2.5051e-12 |
| pooled | intestinal_epithelial_differentiation | 72 | 0.71281 | 2.1622e-12 | 3.0992e-11 |
| pooled | NR1H3 | 72 | 0.62623 | 4.0234e-09 | 4.3252e-08 |
| pooled | MYC | 72 | -0.59649 | 3.2183e-08 | 2.3865e-07 |
| pooled | PPARA | 72 | 0.59598 | 3.33e-08 | 2.3865e-07 |
| pooled | Fatty_acid_metabolism | 72 | 0.55865 | 3.4025e-07 | 2.0901e-06 |

Top regulator-activity contrasts:

| regulator | comparison | n_pairs | activity_delta | P | BH_FDR |
| --- | --- | --- | --- | --- | --- |
| NR1H2 | PPAR_low_vs_high | 39 | -0.60385 | 3.638e-12 | 5.0932e-11 |
| TCF7L2 | tumor_vs_normal | 36 | 1.0415 | 2.9104e-11 | 4.0745e-10 |
| PPARG | tumor_vs_normal | 36 | -1.6205 | 5.8208e-11 | 4.0745e-10 |
| JUN | tumor_vs_normal | 36 | 1.342 | 2.9104e-10 | 1.3582e-09 |
| RELA | tumor_vs_normal | 36 | 1.5476 | 7.276e-10 | 2.5466e-09 |
| PPARG | PPAR_low_vs_high | 39 | -0.55968 | 2.3283e-09 | 1.6298e-08 |
| STAT3 | tumor_vs_normal | 36 | 1.1138 | 8.9349e-09 | 2.5018e-08 |
| NR1I2 | tumor_vs_normal | 36 | -0.74404 | 1.6697e-07 | 3.8959e-07 |

- Tumor×PPAR interaction minimum BH-FDR: **0.499**; no interaction program passes multiplicity control.
- Directly supported regulatory anchors under frozen evidence tags: **RELA, STAT3**.

## Boundary

The CRC epithelial disease-state bridge is now executable from the official source H5AD and passes integrity QC. The analysis remains associative: it does not establish DINP/MCOP exposure as the cause of the PPAR/NR-low state or prove mediation of the epidemiologic association.
