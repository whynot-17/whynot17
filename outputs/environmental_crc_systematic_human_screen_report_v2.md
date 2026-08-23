# Systematic NHANES human screen v2

Generated: 2026-08-23T16:53:30.056378+00:00

## Frozen scope

The complete actionability audit covers 267 original chemicals. 87 chemical rows satisfy the permissive outcome-blinded rule and collapse to 15 unique exposure-axis/primary-analyte tests. Every one of these 15 axes was passed through the same complex-survey logistic model; BH-FDR was calculated over all finite P values from this unified screen.

Model: `CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR`; urinary analytes additionally include `log2(creatinine)`. Serum/blood analytes do not receive creatinine adjustment.

Finite fitted P values: 15/15. Survey weights were taken from the selected analyte's own NHANES file and divided by the number of included cycles for that analyte. The existing validated Taylor-style PSU/strata implementation was used.

## Results

|   screen_rank | exposure_axis                                      | primary_biomarker   |   eligible_chemical_count |   analytic_n |   crc_cases |       OR |   CI_low |   CI_high |          P |    BH_FDR | status                 |
|--------------:|:---------------------------------------------------|:--------------------|--------------------------:|-------------:|------------:|---------:|---------:|----------:|-----------:|----------:|:-----------------------|
|             1 | PFAS exposure axis                                 | LBXPFHS             |                         1 |         5322 |          32 | 0.624414 | 0.470702 |  0.828322 | 0.00146008 | 0.0219012 | converged_with_warning |
|             2 | DINP-related exposure axis                         | URXCOP              |                         3 |         9936 |          70 | 1.24551  | 1.07753  |  1.43968  | 0.00331136 | 0.0248352 | ok                     |
|             3 | phthalate exposure axis                            | URXMOH              |                         1 |        12603 |          97 | 1.20368  | 1.0226   |  1.41683  | 0.0261383  | 0.130691  | ok                     |
|             4 | DEHP-related exposure axis;phthalate exposure axis | URXECP              |                         2 |        11224 |          85 | 1.21093  | 0.987355 |  1.48513  | 0.0658442  | 0.246916  | ok                     |
|             5 | phthalate exposure axis                            | URXMEP              |                         1 |        13775 |         104 | 0.932921 | 0.834064 |  1.0435   | 0.222586   | 0.542938  | ok                     |
|             6 | PFAS exposure axis                                 | LBXPFNA             |                         1 |         5322 |          32 | 0.805854 | 0.555426 |  1.1692   | 0.250749   | 0.542938  | converged_with_warning |
|             7 | PFAS exposure axis                                 | LBXPFDE             |                         1 |         5322 |          32 | 0.716091 | 0.385047 |  1.33175  | 0.286125   | 0.542938  | ok                     |
|             8 | DEHP-related exposure axis                         | URXMHP              |                         1 |        13779 |         104 | 1.09293  | 0.926513 |  1.28924  | 0.289567   | 0.542938  | ok                     |
|             9 | phthalate exposure axis                            | URXMBP              |                         1 |        13779 |         104 | 0.895411 | 0.690063 |  1.16186  | 0.403443   | 0.622786  | ok                     |
|            10 | bisphenol exposure axis                            | URXBPH              |                        17 |         5751 |          45 | 1.06132  | 0.918099 |  1.22688  | 0.41519    | 0.622786  | converged_with_warning |
|            11 | phthalate exposure axis                            | URXMZP              |                         2 |        13779 |         104 | 0.954508 | 0.800546 |  1.13808  | 0.601792   | 0.752778  | ok                     |
|            12 | PAH exposure axis                                  | URXP10              |                        39 |         8265 |          70 | 1.07212  | 0.817917 |  1.40532  | 0.610624   | 0.752778  | ok                     |
|            13 | phthalate exposure axis                            | URXMIB              |                         1 |        12603 |          97 | 0.965853 | 0.821817 |  1.13513  | 0.671228   | 0.752778  | ok                     |
|            14 | PAH exposure axis                                  | URXP04              |                         5 |         8290 |          70 | 0.954216 | 0.748394 |  1.21664  | 0.702592   | 0.752778  | ok                     |
|            15 | PAH exposure axis                                  | URXP02              |                        11 |         8301 |          71 | 1.00297  | 0.778007 |  1.29299  | 0.981543   | 0.981543  | ok                     |

## MCOP and MiNP interpretation

MCOP (`URXCOP`) is included as one of the prespecified eligible axes. MiNP (`URXMNP`) remains a separate DINP molecular nomination but is not a primary human-screen axis because direct MiNP detectability fails the frozen D>=1 gate; it was not removed or merged into MCOP.

This screen is an epidemiologic association scan, not causal evidence. Candidate selection/eligibility was frozen before reading these OR/P/FDR results.
