# DINP cross-tissue perturbation consensus

This report is generated from public GEO expression matrices and explicit sample metadata. The legacy 81-gene list is post-hoc annotation only.

## Dataset status

| dataset_id | tissue | species | contrast_id | dose | timepoint | n_treated | n_control | expression_scale | analysis_method | quality | status | source_url | exclusion_reason | n_tumor | n_normal | method | n_pairs | platform | design | n_genes_tested |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE158473 | ovary | mouse | GSE158473_0month_DiNP_20ugkgd | 20 ug/kg/day | 0 month | 3 | 3 | log2 CPM | Welch t-test | secondary_processed_expression | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158473 | NA | NA | NA | NA | NA | NA | NA | NA |
| GSE158473 | ovary | mouse | GSE158473_0month_DiNP_100ugkgd | 100 ug/kg/day | 0 month | 3 | 3 | log2 CPM | Welch t-test | secondary_processed_expression | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158473 | NA | NA | NA | NA | NA | NA | NA | NA |
| GSE158473 | ovary | mouse | GSE158473_3months_DiNP_20ugkgd | 20 ug/kg/day | 3 months | 3 | 3 | log2 CPM | Welch t-test | secondary_processed_expression | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158473 | NA | NA | NA | NA | NA | NA | NA | NA |
| GSE158473 | ovary | mouse | GSE158473_3months_DiNP_100ugkgd | 100 ug/kg/day | 3 months | 3 | 3 | log2 CPM | Welch t-test | secondary_processed_expression | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158473 | NA | NA | NA | NA | NA | NA | NA | NA |
| GSE313258 | liver | mouse | GSE313258_20weeks_DiNP_1.5mgkgd | 1.5 mg/kg/day | 20 weeks | 13 | 14 | TPM (log2 transformed for test) | Welch t-test | secondary_processed_expression | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313258 | NA | NA | NA | NA | NA | NA | NA | NA |
| GSE313258 | liver | mouse | GSE313258_20weeks_DiNP_15mgkgd | 15 mg/kg/day | 20 weeks | 13 | 14 | TPM (log2 transformed for test) | Welch t-test | secondary_processed_expression | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313258 | NA | NA | NA | NA | NA | NA | NA | NA |
| GSE313258 | liver | mouse | GSE313258_20weeks_DiNP_150mgkgd | 150 mg/kg/day | 20 weeks | 14 | 14 | TPM (log2 transformed for test) | Welch t-test | secondary_processed_expression | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313258 | NA | NA | NA | NA | NA | NA | NA | NA |
| THYROID_ORGANOID | thyroid organoid | mouse |  |  |  |  |  |  |  |  | UNAVAILABLE |  | No unique DINP GEO accession confirmed from the source record; excluded rather than guessed. | NA | NA | NA | NA | NA | NA | NA |
| GSE44076 | NA | NA | NA | NA | NA | NA | NA | log2 RMA-like | NA | NA | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44076 | NA | 98 | 98 | paired | 98 | GPL13667 | paired adjacent-normal vs primary tumor | 20217 |
| GSE37364 | NA | NA | NA | NA | NA | NA | NA | log2(x+1) probe intensity | NA | NA | PASS | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37364 | NA | 27 | 38 | welch | 0 | GPL570 | unpaired normal mucosa vs CRC biopsy | 22834 |

## Analysis decisions

- GSE158473 is processed logCPM and GSE313258 is processed TPM; no raw-count DESeq2 call was fabricated. Welch t-tests on the native processed scale are marked secondary-quality.
- Dose and time contrasts are kept separate. Independent-dataset consensus counts GSE158473 and GSE313258 once each, never each dose as an independent tissue.
- Mouse-to-human conversion uses g:Orth, which retrieves Ensembl-backed canonical orthologues. Unmapped genes remain auditable.
- Pathway output uses a transparent preranked rank-sum implementation over locally available Hallmark and Reactome GMTs. KEGG/GO:BP are recorded as unavailable when no local gene-set file is present.

## DINP consensus

| gene | n_independent_datasets | n_significant_datasets | direction | all_directions_same | median_log2FC | min_FDR | meta_log2FC | meta_SE | meta_p | strict_supported | direction_supported | label |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EHHADH | 2 | 1 | up | True | 0.97342 | 3.3893e-07 | 1.0023 | 0.56747 | 0.077362 | False | True | DINP |
| ACOX1 | 2 | 1 | up | True | 0.39868 | 1.6395e-06 | 0.41152 | 0.24271 | 0.089975 | False | True | DINP |
| ANGPTL4 | 2 | 1 | up | True | 0.87303 | 1.0761e-05 | 0.88627 | 0.59786 | 0.13823 | False | True | DINP |
| KIAA1671 | 2 | 1 | up | True | 0.4381 | 1.5792e-05 | 0.55479 | 0.17308 | 0.0013492 | False | True | DINP |
| ACOT1 | 2 | 1 | up | True | 1.2157 | 2.1217e-05 | 1.4402 | 0.23455 | 8.2426e-10 | False | True | DINP |
| CYP4A11 | 2 | 1 | up | True | 1.3002 | 2.1217e-05 | 1.3059 | 0.43299 | 0.0025608 | False | True | DINP |
| HSD17B11 | 2 | 1 | up | True | 0.46208 | 2.1217e-05 | 0.47416 | 0.2783 | 0.088419 | False | True | DINP |
| HSDL2 | 2 | 1 | up | True | 0.26038 | 2.1217e-05 | 0.27452 | 0.15086 | 0.06881 | False | True | DINP |
| ECI1 | 2 | 1 | up | True | 0.43969 | 0.00010811 | 0.47788 | 0.06121 | 5.8481e-15 | False | True | DINP |
| BACE1 | 2 | 1 | up | True | 0.37664 | 0.00012396 | 0.47084 | 0.06236 | 4.3423e-14 | False | True | DINP |
| PDP2 | 2 | 1 | up | True | 0.38235 | 0.00040754 | 0.40156 | 0.22102 | 0.069241 | False | True | DINP |
| HMGCL | 2 | 1 | up | True | 0.20341 | 0.00082267 | 0.21946 | 0.068859 | 0.0014372 | False | True | DINP |
| RAB4B | 2 | 1 | down | True | -0.4411 | 0.00090744 | -0.45932 | 0.20167 | 0.02275 | False | True | DINP |
| SLC25A22 | 2 | 1 | down | True | -0.4083 | 0.00095057 | -0.41176 | 0.20523 | 0.044821 | False | True | DINP |
| CERS2 | 2 | 1 | up | True | 0.23904 | 0.001576 | 0.26317 | 0.0372 | 1.5016e-12 | False | True | DINP |
| FBF1 | 2 | 1 | up | True | 0.68241 | 0.0017429 | 0.66603 | 0.62035 | 0.28298 | False | True | DINP |
| POLR2M | 2 | 1 | up | True | 0.4003 | 0.0018709 | 0.39447 | 0.2995 | 0.1878 | False | True | DINP |
| RPL35 | 2 | 1 | up | True | 0.45656 | 0.0019094 | 0.5053 | 0.078 | 9.2848e-11 | False | True | DINP |
| RPL4 | 2 | 1 | up | True | 0.49068 | 0.0019852 | 0.49736 | 0.3369 | 0.13987 | False | True | DINP |
| ETNK1 | 2 | 1 | down | True | -0.32323 | 0.0022573 | -0.36236 | 0.052901 | 7.3981e-12 | False | True | DINP |
| ALDH9A1 | 2 | 1 | up | True | 0.17642 | 0.0024197 | 0.18602 | 0.13595 | 0.17122 | False | True | DINP |
| TMEM134 | 2 | 1 | up | True | 0.20694 | 0.0026869 | 0.20478 | 0.12922 | 0.11302 | False | True | DINP |
| CALU | 2 | 1 | down | True | -0.36461 | 0.0031139 | -0.36651 | 0.20557 | 0.07461 | False | True | DINP |
| PIN4 | 2 | 1 | up | True | 0.79134 | 0.0031139 | 0.81285 | 0.463 | 0.079157 | False | True | DINP |
| RCC1 | 2 | 1 | down | True | -0.5338 | 0.0031198 | -0.53384 | 0.40215 | 0.18435 | False | True | DINP |
| LLGL2 | 2 | 1 | up | True | 0.35633 | 0.0031984 | 0.35863 | 0.13314 | 0.0070697 | False | True | DINP |
| FAU | 2 | 1 | up | True | 0.3429 | 0.0036544 | 0.34205 | 0.052609 | 7.9437e-11 | False | True | DINP |
| CBX1 | 2 | 1 | down | True | -0.23888 | 0.0036669 | -0.24019 | 0.11856 | 0.042779 | False | True | DINP |
| SCAMP3 | 2 | 1 | up | True | 0.30148 | 0.0036669 | 0.2989 | 0.15203 | 0.049291 | False | True | DINP |
| UBA1 | 2 | 1 | up | True | 0.2691 | 0.0036669 | 0.26701 | 0.11832 | 0.024025 | False | True | DINP |

## CRC expression consensus

| gene | n_independent_datasets | n_significant_datasets | direction | all_directions_same | median_log2FC | min_FDR | meta_log2FC | meta_SE | meta_p | strict_supported | direction_supported | label |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CA7 | 2 | 2 | down | True | -4.4634 | 1.3716e-54 | -4.6129 | 0.11949 | 0 | True | True | CRC |
| ETV4 | 2 | 2 | up | True | 3.4911 | 3.0535e-54 | 3.4518 | 0.27797 | 2.0863e-35 | True | True | CRC |
| GUCA2B | 2 | 2 | down | True | -6.5522 | 2.3331e-53 | -6.548 | 0.17944 | 1.5465e-291 | True | True | CRC |
| CLDN1 | 2 | 2 | up | True | 3.5525 | 3.3188e-53 | 3.5255 | 0.26285 | 5.0959e-41 | True | True | CRC |
| ABCA8 | 2 | 2 | down | True | -1.6635 | 6.9159e-53 | -1.6344 | 0.35101 | 3.2177e-06 | True | True | CRC |
| CDH3 | 2 | 2 | up | True | 4.8994 | 1.0528e-52 | 4.8847 | 1.8667 | 0.0088773 | True | True | CRC |
| GTF2IRD1 | 2 | 2 | up | True | 1.9361 | 1.4563e-52 | 1.9427 | 0.43672 | 8.6491e-06 | True | True | CRC |
| UGP2 | 2 | 2 | down | True | -1.4366 | 3.4008e-52 | -1.4552 | 0.17378 | 5.5702e-17 | True | True | CRC |
| BEST4 | 2 | 2 | down | True | -4.3679 | 8.8021e-52 | -4.2048 | 0.11866 | 4.5828e-275 | True | True | CRC |
| AQP8 | 2 | 2 | down | True | -7.0341 | 1.1471e-51 | -6.7473 | 0.30334 | 1.3026e-109 | True | True | CRC |
| PPM1H | 2 | 2 | up | True | 1.7295 | 4.8848e-51 | 1.7324 | 0.84337 | 0.039965 | True | True | CRC |
| FOXQ1 | 2 | 2 | up | True | 4.5118 | 1.2719e-50 | 4.5223 | 0.65669 | 5.7212e-12 | True | True | CRC |
| TMEM206 | 2 | 2 | up | True | 1.1396 | 1.8266e-50 | 1.142 | 0.73374 | 0.11961 | True | True | CRC |
| GLTP | 2 | 2 | down | True | -1.3447 | 6.2899e-50 | -1.357 | 0.23342 | 6.1126e-09 | True | True | CRC |
| ASCL2 | 2 | 2 | up | True | 3.1264 | 2.5908e-49 | 3.2095 | 0.498 | 1.1573e-10 | True | True | CRC |
| GREM2 | 2 | 2 | down | True | -2.7889 | 2.6219e-49 | -2.8318 | 0.55876 | 4.0217e-07 | True | True | CRC |
| TMEM100 | 2 | 2 | down | True | -2.8962 | 2.9299e-49 | -2.9446 | 0.52871 | 2.5553e-08 | True | True | CRC |
| EPB41L3 | 2 | 2 | down | True | -2.0593 | 1.8643e-48 | -2.1238 | 0.066868 | 2.2074e-221 | True | True | CRC |
| SLC6A6 | 2 | 2 | up | True | 1.3207 | 2.0267e-48 | 1.2982 | 0.34913 | 0.00020052 | True | True | CRC |
| IL6R | 2 | 2 | down | True | -1.5497 | 3.1284e-48 | -1.4702 | 0.093612 | 1.388e-55 | True | True | CRC |
| ABCG2 | 2 | 2 | down | True | -4.0972 | 3.3654e-48 | -3.9909 | 0.1654 | 1.2544e-128 | True | True | CRC |
| XPOT | 2 | 2 | up | True | 1.3955 | 4.7982e-48 | 1.4358 | 0.05715 | 2.6952e-139 | True | True | CRC |
| RDH5 | 2 | 2 | down | True | -1.8184 | 7.1483e-48 | -1.7825 | 0.32744 | 5.2157e-08 | True | True | CRC |
| NONO | 2 | 2 | up | True | 0.87235 | 1.5921e-47 | 0.87693 | 0.27668 | 0.0015272 | True | True | CRC |
| MTHFD1L | 2 | 2 | up | True | 1.6667 | 2.2931e-47 | 1.6686 | 0.51182 | 0.001114 | True | True | CRC |
| GFRA2 | 2 | 2 | down | True | -1.7414 | 3.1293e-47 | -1.7307 | 0.77248 | 0.025064 | True | True | CRC |
| ENC1 | 2 | 2 | up | True | 1.5621 | 4.1498e-47 | 1.5762 | 0.28016 | 1.8433e-08 | True | True | CRC |
| SLC4A4 | 2 | 2 | down | True | -3.5002 | 5.323e-47 | -3.4073 | 0.36252 | 5.5064e-21 | True | True | CRC |
| ADH1B | 2 | 2 | down | True | -3.637 | 5.843e-47 | -3.6727 | 0.82498 | 8.5134e-06 | True | True | CRC |
| SLCO4A1 | 2 | 2 | up | True | 1.5625 | 9.4493e-47 | 1.5915 | 0.052593 | 3.784e-201 | True | True | CRC |

## Direction-aware DINP × CRC convergence

| gene | n_independent_datasets_DINP | n_significant_datasets_DINP | direction_DINP | all_directions_same_DINP | median_log2FC_DINP | min_FDR_DINP | meta_log2FC_DINP | meta_SE_DINP | meta_p_DINP | strict_supported_DINP | direction_supported_DINP | label_DINP | n_independent_datasets_CRC | n_significant_datasets_CRC | direction_CRC | all_directions_same_CRC | median_log2FC_CRC | min_FDR_CRC | meta_log2FC_CRC | meta_SE_CRC | meta_p_CRC | strict_supported_CRC | direction_supported_CRC | label_CRC | same_direction | opposite_direction | legacy_81_annotation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EHHADH | 2 | 1 | up | True | 0.97342 | 3.3893e-07 | 1.0023 | 0.56747 | 0.077362 | False | True | DINP | 2 | 2 | down | True | -0.94953 | 1.7189e-19 | -0.89978 | 0.1644 | 4.4177e-08 | True | True | CRC | False | True | False |
| ACOX1 | 2 | 1 | up | True | 0.39868 | 1.6395e-06 | 0.41152 | 0.24271 | 0.089975 | False | True | DINP | 2 | 2 | down | True | -1.0873 | 2.2206e-29 | -1.0899 | 0.40515 | 0.0071457 | True | True | CRC | False | True | False |
| ANGPTL4 | 2 | 1 | up | True | 0.87303 | 1.0761e-05 | 0.88627 | 0.59786 | 0.13823 | False | True | DINP | 2 | 1 | up | True | 1.0244 | 2.0431e-08 | 1.0089 | 1.0191 | 0.32221 | False | True | CRC | True | False | False |
| KIAA1671 | 2 | 1 | up | True | 0.4381 | 1.5792e-05 | 0.55479 | 0.17308 | 0.0013492 | False | True | DINP | 2 | 2 | down | True | -0.46725 | 8.8973e-15 | -0.4333 | 0.073208 | 3.2438e-09 | True | True | CRC | False | True | False |
| ACOT1 | 2 | 1 | up | True | 1.2157 | 2.1217e-05 | 1.4402 | 0.23455 | 8.2426e-10 | False | True | DINP | 2 | 2 | down | True | -0.5508 | 8.4662e-12 | -0.56028 | 0.059656 | 5.8981e-21 | True | True | CRC | False | True | False |
| CYP4A11 | 2 | 1 | up | True | 1.3002 | 2.1217e-05 | 1.3059 | 0.43299 | 0.0025608 | False | True | DINP | 2 | 0 | down | True | -0.037026 | 0.001669 | -0.053722 | 0.015809 | 0.0006783 | False | True | CRC | False | True | False |
| HSD17B11 | 2 | 1 | up | True | 0.46208 | 2.1217e-05 | 0.47416 | 0.2783 | 0.088419 | False | True | DINP | 2 | 2 | down | True | -1.09 | 4.3762e-29 | -1.0932 | 0.057757 | 6.8406e-80 | True | True | CRC | False | True | False |
| HSDL2 | 2 | 1 | up | True | 0.26038 | 2.1217e-05 | 0.27452 | 0.15086 | 0.06881 | False | True | DINP | 2 | 0 | down | False | -0.041935 | 0.12842 | -0.032352 | 0.1189 | 0.78556 | False | False | CRC | False | True | False |
| ECI1 | 2 | 1 | up | True | 0.43969 | 0.00010811 | 0.47788 | 0.06121 | 5.8481e-15 | False | True | DINP | 1 | 1 | down | True | -0.3984 | 0.00041428 | NA | NA | NA | False | False | CRC | False | True | False |
| BACE1 | 2 | 1 | up | True | 0.37664 | 0.00012396 | 0.47084 | 0.06236 | 4.3423e-14 | False | True | DINP | 2 | 0 | up | True | 0.19932 | 5.3069e-06 | 0.22229 | 0.043047 | 2.4204e-07 | False | True | CRC | True | False | False |
| PDP2 | 2 | 1 | up | True | 0.38235 | 0.00040754 | 0.40156 | 0.22102 | 0.069241 | False | True | DINP | 2 | 0 | up | True | 0.14289 | 0.0018456 | 0.13109 | 0.03476 | 0.00016241 | False | True | CRC | True | False | False |
| HMGCL | 2 | 1 | up | True | 0.20341 | 0.00082267 | 0.21946 | 0.068859 | 0.0014372 | False | True | DINP | 2 | 2 | down | True | -0.93993 | 3.8555e-29 | -0.92774 | 0.047948 | 2.0739e-83 | True | True | CRC | False | True | False |
| RAB4B | 2 | 1 | down | True | -0.4411 | 0.00090744 | -0.45932 | 0.20167 | 0.02275 | False | True | DINP | 2 | 2 | down | True | -0.85514 | 2.214e-30 | -0.8664 | 0.15218 | 1.2457e-08 | True | True | CRC | True | False | False |
| SLC25A22 | 2 | 1 | down | True | -0.4083 | 0.00095057 | -0.41176 | 0.20523 | 0.044821 | False | True | DINP | 2 | 2 | up | True | 0.48634 | 4.2529e-25 | 0.5086 | 0.077931 | 6.742e-11 | True | True | CRC | False | True | False |
| CERS2 | 2 | 1 | up | True | 0.23904 | 0.001576 | 0.26317 | 0.0372 | 1.5016e-12 | False | True | DINP | 1 | 0 | up | True | 0.18048 | 0.046449 | NA | NA | NA | False | False | CRC | True | False | False |
| FBF1 | 2 | 1 | up | True | 0.68241 | 0.0017429 | 0.66603 | 0.62035 | 0.28298 | False | True | DINP | 2 | 0 | up | True | 0.088088 | 0.28778 | 0.040089 | 0.064849 | 0.53645 | False | True | CRC | True | False | False |
| POLR2M | 2 | 1 | up | True | 0.4003 | 0.0018709 | 0.39447 | 0.2995 | 0.1878 | False | True | DINP | 1 | 0 | down | True | -0.12698 | 0.19112 | NA | NA | NA | False | False | CRC | False | True | False |
| RPL35 | 2 | 1 | up | True | 0.45656 | 0.0019094 | 0.5053 | 0.078 | 9.2848e-11 | False | True | DINP | 2 | 1 | up | True | 0.27776 | 4.9165e-27 | 0.28022 | 0.20805 | 0.17802 | False | True | CRC | True | False | False |
| RPL4 | 2 | 1 | up | True | 0.49068 | 0.0019852 | 0.49736 | 0.3369 | 0.13987 | False | True | DINP | 2 | 1 | up | True | 0.22792 | 4.972e-18 | 0.23895 | 0.11364 | 0.035494 | False | True | CRC | True | False | False |
| ETNK1 | 2 | 1 | down | True | -0.32323 | 0.0022573 | -0.36236 | 0.052901 | 7.3981e-12 | False | True | DINP | 2 | 2 | down | True | -0.91994 | 6.4744e-14 | -0.95671 | 0.1598 | 2.1404e-09 | True | True | CRC | True | False | False |
| ALDH9A1 | 2 | 1 | up | True | 0.17642 | 0.0024197 | 0.18602 | 0.13595 | 0.17122 | False | True | DINP | 2 | 1 | down | True | -0.2863 | 6.5155e-06 | -0.27931 | 0.11894 | 0.018861 | False | True | CRC | False | True | False |
| TMEM134 | 2 | 1 | up | True | 0.20694 | 0.0026869 | 0.20478 | 0.12922 | 0.11302 | False | True | DINP | 2 | 0 | up | True | 0.086167 | 0.099859 | 0.060134 | 0.031205 | 0.053972 | False | True | CRC | True | False | False |
| CALU | 2 | 1 | down | True | -0.36461 | 0.0031139 | -0.36651 | 0.20557 | 0.07461 | False | True | DINP | 2 | 2 | up | True | 1.2697 | 9.1177e-20 | 1.2606 | 0.12909 | 1.5863e-22 | True | True | CRC | False | True | False |
| PIN4 | 2 | 1 | up | True | 0.79134 | 0.0031139 | 0.81285 | 0.463 | 0.079157 | False | True | DINP | 2 | 0 | up | True | 0.066286 | 0.029747 | 0.078626 | 0.043627 | 0.07151 | False | True | CRC | True | False | False |
| RCC1 | 2 | 1 | down | True | -0.5338 | 0.0031198 | -0.53384 | 0.40215 | 0.18435 | False | True | DINP | 2 | 2 | up | True | 0.84021 | 1.0214e-24 | 0.84698 | 0.14426 | 4.3282e-09 | True | True | CRC | False | True | False |
| LLGL2 | 2 | 1 | up | True | 0.35633 | 0.0031984 | 0.35863 | 0.13314 | 0.0070697 | False | True | DINP | 2 | 1 | down | True | -0.34773 | 0.0002312 | -0.343 | 0.10618 | 0.0012366 | False | True | CRC | False | True | False |
| FAU | 2 | 1 | up | True | 0.3429 | 0.0036544 | 0.34205 | 0.052609 | 7.9437e-11 | False | True | DINP | 2 | 0 | up | True | 0.08328 | 4.4368e-07 | 0.09355 | 0.049071 | 0.056596 | False | True | CRC | True | False | False |
| CBX1 | 2 | 1 | down | True | -0.23888 | 0.0036669 | -0.24019 | 0.11856 | 0.042779 | False | True | DINP | 2 | 2 | up | True | 0.47736 | 2.6448e-16 | 0.48416 | 0.10599 | 4.9294e-06 | True | True | CRC | False | True | False |
| SCAMP3 | 2 | 1 | up | True | 0.30148 | 0.0036669 | 0.2989 | 0.15203 | 0.049291 | False | True | DINP | 2 | 1 | up | True | 0.35892 | 5.2026e-15 | 0.35894 | 0.14048 | 0.010613 | False | True | CRC | True | False | False |
| UBA1 | 2 | 1 | up | True | 0.2691 | 0.0036669 | 0.26701 | 0.11832 | 0.024025 | False | True | DINP | 2 | 0 | up | False | 0.049176 | 0.0020976 | 0.088895 | 0.030962 | 0.0040904 | False | False | CRC | True | False | False |
| ADIPOR1 | 2 | 1 | up | True | 0.18833 | 0.0037174 | 0.19489 | 0.10978 | 0.075846 | False | True | DINP | 2 | 1 | up | True | 0.19897 | 8.9282e-08 | 0.20502 | 0.10719 | 0.055785 | False | True | CRC | True | False | False |
| GC | 2 | 1 | up | True | 0.80411 | 0.0037174 | 0.53198 | 0.42581 | 0.21153 | False | True | DINP | 2 | 0 | down | True | -0.055056 | 0.66581 | -0.011404 | 0.021889 | 0.60236 | False | True | CRC | False | True | False |
| SEC62 | 2 | 1 | up | True | 0.40562 | 0.0037174 | 0.4581 | 0.22872 | 0.045193 | False | True | DINP | 2 | 0 | down | False | -0.038692 | 0.081299 | -0.056317 | 0.03871 | 0.14571 | False | False | CRC | False | True | False |
| LPIN2 | 2 | 1 | up | True | 0.51516 | 0.0037493 | 0.50108 | 0.43507 | 0.24944 | False | True | DINP | 2 | 2 | down | True | -0.38037 | 1.1357e-09 | -0.40469 | 0.052996 | 2.234e-14 | True | True | CRC | False | True | False |
| OPLAH | 2 | 1 | up | True | 0.46446 | 0.0037493 | 0.43686 | 0.077425 | 1.6775e-08 | False | True | DINP | 2 | 0 | down | True | -0.19796 | 0.010263 | -0.2006 | 0.060713 | 0.0009532 | False | True | CRC | False | True | False |
| CPT2 | 2 | 1 | up | True | 0.33942 | 0.0037992 | 0.34103 | 0.05038 | 1.2959e-11 | False | True | DINP | 2 | 2 | down | True | -1.3048 | 1.1132e-32 | -1.3126 | 0.058086 | 4.6643e-113 | True | True | CRC | False | True | False |
| GNAS | 2 | 1 | up | True | 0.24727 | 0.003936 | 0.24246 | 0.17161 | 0.1577 | False | True | DINP | 2 | 1 | up | True | 0.1766 | 7.5469e-08 | 0.1767 | 0.17177 | 0.30362 | False | True | CRC | True | False | False |
| TAX1BP1 | 2 | 1 | up | True | 0.22148 | 0.003936 | 0.22082 | 0.116 | 0.056954 | False | True | DINP | 2 | 0 | up | True | 0.053522 | 0.10306 | 0.050296 | 0.044711 | 0.26063 | False | True | CRC | True | False | False |
| ALDOB | 2 | 1 | up | True | 0.53216 | 0.004562 | 0.52142 | 0.091706 | 1.302e-08 | False | True | DINP | 2 | 1 | up | True | 0.38418 | 3.4892e-05 | 0.46646 | 0.23411 | 0.046317 | False | True | CRC | True | False | False |
| RXRA | 2 | 1 | up | True | 0.23582 | 0.004562 | 0.22973 | 0.19564 | 0.2403 | False | True | DINP | 2 | 2 | down | True | -0.53055 | 4.9283e-19 | -0.55895 | 0.045389 | 7.5442e-35 | True | True | CRC | False | True | False |
| THRA | 2 | 1 | down | True | -0.3125 | 0.004562 | -0.30959 | 0.2609 | 0.23539 | False | True | DINP | 2 | 0 | up | True | 0.06837 | 0.41318 | 0.025333 | 0.036855 | 0.49185 | False | True | CRC | False | True | False |
| YTHDF2 | 2 | 1 | down | True | -0.32248 | 0.004562 | -0.31729 | 0.26869 | 0.23765 | False | True | DINP | 2 | 0 | up | True | 0.14039 | 6.2794e-06 | 0.14772 | 0.038875 | 0.00014482 | False | True | CRC | False | True | False |
| SLC25A20 | 2 | 1 | up | True | 0.20472 | 0.0047227 | 0.1985 | 0.079827 | 0.012895 | False | True | DINP | 2 | 2 | down | True | -1.3714 | 2.1749e-35 | -1.4124 | 0.062282 | 7.4139e-114 | True | True | CRC | False | True | False |
| IFNAR2 | 2 | 1 | up | True | 0.33226 | 0.0047693 | 0.33422 | 0.055027 | 1.2497e-09 | False | True | DINP | 2 | 0 | up | True | 0.18787 | 6.8262e-06 | 0.1787 | 0.033324 | 8.205e-08 | False | True | CRC | True | False | False |
| RAPGEF1 | 2 | 1 | up | True | 0.28829 | 0.0047693 | 0.28228 | 0.21165 | 0.1823 | False | True | DINP | 2 | 0 | down | False | -0.033512 | 6.7214e-08 | -0.062128 | 0.087658 | 0.47848 | False | False | CRC | False | True | False |
| BCL2L13 | 2 | 1 | up | True | 0.25107 | 0.0048064 | 0.25044 | 0.10665 | 0.018858 | False | True | DINP | 2 | 1 | down | True | -0.22541 | 3.5552e-06 | -0.21967 | 0.088395 | 0.01295 | False | True | CRC | False | True | False |
| ECH1 | 2 | 1 | up | True | 0.31855 | 0.0048064 | 0.31459 | 0.20993 | 0.134 | False | True | DINP | 2 | 2 | down | True | -0.82018 | 4.2073e-17 | -0.81283 | 0.13754 | 3.4269e-09 | True | True | CRC | False | True | False |
| NTRK2 | 2 | 1 | down | True | -0.70841 | 0.0048064 | -0.78053 | 0.4127 | 0.058587 | False | True | DINP | 2 | 2 | down | True | -0.98888 | 5.6666e-09 | -0.97254 | 0.57914 | 0.093096 | True | True | CRC | True | False | False |
| ARL1 | 2 | 1 | up | True | 0.35012 | 0.0049241 | 0.34101 | 0.25353 | 0.17861 | False | True | DINP | 2 | 1 | up | True | 0.21937 | 2.2261e-09 | 0.24917 | 0.03542 | 1.9981e-12 | False | True | CRC | True | False | False |
| COX20 | 2 | 1 | up | True | 0.34977 | 0.0049241 | 0.36858 | 0.14657 | 0.011915 | False | True | DINP | 1 | 0 | up | True | 0.15277 | 0.06222 | NA | NA | NA | False | False | CRC | True | False | False |

## Legacy post-hoc audit

| gene | n_independent_datasets | n_significant_datasets | direction | all_directions_same | median_log2FC | min_FDR | meta_log2FC | meta_SE | meta_p | strict_supported | direction_supported | label | legacy_81_annotation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RXRA | 2 | 1 | up | True | 0.23582 | 0.004562 | 0.22973 | 0.19564 | 0.2403 | False | True | DINP | True |
| SIRT2 | 2 | 0 | up | True | 0.15681 | 0.019083 | 0.15196 | 0.068347 | 0.026191 | False | True | DINP | True |
| CYP2A6 | 2 | 0 | up | True | 1.0315 | 0.080131 | 0.76778 | 0.22314 | 0.00058015 | False | True | DINP | True |
| DDIT4 | 2 | 0 | down | True | -0.55838 | 0.082725 | -0.53808 | 0.14256 | 0.00016043 | False | True | DINP | True |
| ESR1 | 2 | 0 | down | True | -0.31499 | 0.10509 | -0.32743 | 0.092997 | 0.00043008 | False | True | DINP | True |
| SUSD2 | 2 | 0 | up | True | 0.58463 | 0.12163 | 0.58606 | 0.10652 | 3.7531e-08 | False | True | DINP | True |
| NR3C1 | 2 | 0 | up | True | 0.14936 | 0.22457 | 0.15139 | 0.050806 | 0.0028848 | False | True | DINP | True |
| ADAMTS1 | 2 | 0 | down | True | -0.32132 | 0.24351 | -0.28086 | 0.086824 | 0.0012171 | False | True | DINP | True |
| NR1I2 | 2 | 0 | up | True | 0.37562 | 0.27586 | 0.28658 | 0.099869 | 0.00411 | False | True | DINP | True |
| HSPA1A | 2 | 0 | down | True | -0.81466 | 0.29987 | -0.7436 | 0.47447 | 0.11707 | False | True | DINP | True |
| PPARD | 2 | 0 | down | True | -0.14029 | 0.36474 | -0.11821 | 0.077254 | 0.12599 | False | True | DINP | True |
| CPS1 | 2 | 0 | up | True | 0.42325 | 0.44687 | 0.35742 | 0.1545 | 0.020702 | False | True | DINP | True |
| ATF4 | 2 | 0 | down | True | -0.21957 | 0.46051 | -0.20991 | 0.072452 | 0.0037641 | False | True | DINP | True |
| DUSP1 | 2 | 0 | down | True | -0.31298 | 0.64146 | -0.26561 | 0.19693 | 0.1774 | False | True | DINP | True |
| SNORA65 | 2 | 0 | down | True | -0.16402 | 0.64146 | -0.16655 | 0.13858 | 0.22945 | False | True | DINP | True |
| PTGES | 2 | 0 | down | True | -0.25199 | 0.69186 | -0.17506 | 0.16314 | 0.28325 | False | True | DINP | True |
| SLCO2A1 | 2 | 0 | up | True | 0.21516 | 0.69186 | 0.20319 | 0.12982 | 0.11754 | False | True | DINP | True |
| PLA2G4A | 2 | 0 | down | True | -0.33714 | 0.71689 | -0.20958 | 0.11425 | 0.066585 | False | True | DINP | True |
| PTGES2 | 2 | 0 | up | True | 0.16487 | 0.71689 | 0.11092 | 0.088433 | 0.20972 | False | True | DINP | True |
| PTGS1 | 2 | 0 | down | True | -0.42965 | 0.71689 | -0.36606 | 0.29132 | 0.20891 | False | True | DINP | True |
| PTGER2 | 2 | 0 | up | True | 0.22812 | 0.83612 | 0.15759 | 0.15177 | 0.29913 | False | True | DINP | True |
| CXCL13 | 2 | 0 | up | True | 0.37891 | 0.87233 | 0.22193 | 0.22043 | 0.31404 | False | True | DINP | True |
| PPARA | 2 | 1 | down | False | -0.017974 | 0.010974 | 0.027154 | 0.60096 | 0.96396 | False | False | DINP | True |
| SIRT3 | 2 | 1 | up | False | 0.1454 | 0.022666 | 0.17139 | 0.27838 | 0.53811 | False | False | DINP | True |
| SPP2 | 1 | 1 | up | True | 0.3128 | 0.038056 | NA | NA | NA | False | False | DINP | True |
| RXRB | 2 | 1 | down | False | -0.16219 | 0.044602 | -0.16526 | 0.25919 | 0.52372 | False | False | DINP | True |
| AHSG | 2 | 0 | up | False | 0.15859 | 0.060229 | 0.093186 | 0.50512 | 0.85363 | False | False | DINP | True |
| STAT3 | 2 | 0 | down | False | -0.04897 | 0.062096 | -0.066792 | 0.26716 | 0.80258 | False | False | DINP | True |
| C4BPA | 1 | 0 | down | True | -0.20206 | 0.065473 | NA | NA | NA | False | False | DINP | True |
| SIRT5 | 2 | 0 | up | False | 0.073857 | 0.11726 | 0.075643 | 0.29821 | 0.79976 | False | False | DINP | True |
| TIMP2 | 2 | 0 | down | False | -0.19964 | 0.17586 | -0.18005 | 0.36266 | 0.61956 | False | False | DINP | True |
| PARP10 | 2 | 0 | up | False | 0.024785 | 0.23764 | 0.015315 | 0.31531 | 0.96126 | False | False | DINP | True |
| HSD3B2 | 2 | 0 | up | False | 0.082583 | 0.28341 | 0.081362 | 0.30213 | 0.7877 | False | False | DINP | True |
| HPGD | 2 | 0 | up | False | 0.0657 | 0.31019 | 0.12291 | 0.21986 | 0.57614 | False | False | DINP | True |
| SQSTM1 | 2 | 0 | up | False | 0.0013478 | 0.31103 | -0.026648 | 0.1751 | 0.87904 | False | False | DINP | True |
| MMP9 | 2 | 0 | up | False | 0.31129 | 0.31622 | 0.2814 | 0.76525 | 0.71308 | False | False | DINP | True |
| LUM | 2 | 0 | down | False | -0.077363 | 0.32986 | -0.097665 | 0.4339 | 0.82191 | False | False | DINP | True |
| RELA | 2 | 0 | down | False | -0.069157 | 0.37848 | -0.052799 | 0.13303 | 0.69144 | False | False | DINP | True |
| SIRT1 | 2 | 0 | up | False | 0.038175 | 0.39744 | 0.057121 | 0.25567 | 0.82321 | False | False | DINP | True |
| IGFBP1 | 1 | 0 | down | True | -0.60356 | 0.50265 | NA | NA | NA | False | False | DINP | True |
| ATG5 | 2 | 0 | up | False | 0.047717 | 0.52965 | 0.053024 | 0.090018 | 0.55583 | False | False | DINP | True |
| TIMP1 | 2 | 0 | up | False | 0.14004 | 0.53807 | 0.038213 | 0.76141 | 0.95997 | False | False | DINP | True |
| NR1I3 | 1 | 0 | up | True | 0.37361 | 0.56074 | NA | NA | NA | False | False | DINP | True |
| RGS2 | 2 | 0 | down | False | -0.0076796 | 0.57552 | -0.051797 | 0.29031 | 0.8584 | False | False | DINP | True |
| HSPA1L | 2 | 0 | down | False | -0.1811 | 0.64146 | 0.0075063 | 0.17384 | 0.96556 | False | False | DINP | True |
| BECN1 | 2 | 0 | up | False | 0.031414 | 0.65159 | 0.0008499 | 0.11288 | 0.99399 | False | False | DINP | True |
| DPT | 2 | 0 | down | False | -0.041554 | 0.67286 | -0.0078883 | 0.3779 | 0.98335 | False | False | DINP | True |
| MMP2 | 2 | 0 | up | False | 0.019076 | 0.67286 | 0.034169 | 0.42631 | 0.93612 | False | False | DINP | True |
| PPARG | 2 | 0 | down | False | -0.019519 | 0.67286 | -0.009962 | 0.1749 | 0.95458 | False | False | DINP | True |
| PTGER4 | 2 | 0 | up | False | 0.019611 | 0.67286 | -0.030449 | 0.13433 | 0.82068 | False | False | DINP | True |

## Overlap statistics

```json
{
  "dinp_mapped_gene_universe": 15027,
  "crc_gene_universe": 26219,
  "jointly_testable_universe": 14194,
  "dinp_strict_n": 0,
  "crc_strict_n": 4980,
  "strict_same_direction_overlap_n": 0,
  "hypergeom_p_strict_overlap": null
}
```

## Interpretation

Strict support requires at least two independent DINP datasets, significance in each contributing dataset, and no direction conflict. With two confirmed DINP tissues, this is a deliberately conservative gate; an empty strict list is a valid result and calls for external DINP replication rather than threshold relaxation.
