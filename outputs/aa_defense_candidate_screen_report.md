# AA-routing proxy versus candidate ferroptosis-defense axes

## Scope

This is a targeted transcriptomic screen of AIFM2/FSP1, GCH1, SCD, MBOAT1, MBOAT2, GPX4, SLC7A11 and DHODH in GSE39582 and TCGA-COAD.
It tests sidedness, proxy-high/low differences and continuous coupling. It does not measure tissue AA, lipid flux, ferroptotic death or genetic dependency.

## AA-high definition

The primary AA-routing proxy is the within-cohort mean z-score of PLA2G4A, PTGS2 and PTGES. The sensitivity proxy adds ALOX5 and ALOX5AP.
Because patient-level AA lipidomics are unavailable, AA-high/low means high/low transcriptomic routing proxy, not measured AA concentration.

## Primary right-sided results

| Cohort | Candidate | Left n | Right n | Right−left mean z | Welch P | FDR |
|---|---|---:|---:|---:|---:|---:|
| GSE39582 | AIFM2 | 342 | 224 | 0.070 | 0.416 | 0.416 |
| GSE39582 | GCH1 | 342 | 224 | 0.511 | 2.19e-09 | 1.75e-08 |
| GSE39582 | SCD | 342 | 224 | -0.123 | 0.152 | 0.174 |
| GSE39582 | MBOAT1 | 342 | 224 | 0.297 | 0.000483 | 0.00129 |
| GSE39582 | MBOAT2 | 342 | 224 | 0.162 | 0.0691 | 0.0922 |
| GSE39582 | GPX4 | 342 | 224 | -0.171 | 0.0487 | 0.0779 |
| GSE39582 | SLC7A11 | 342 | 224 | 0.427 | 1.28e-06 | 5.12e-06 |
| GSE39582 | DHODH | 342 | 224 | -0.213 | 0.014 | 0.0281 |
| TCGA-COAD | AIFM2 | 135 | 200 | -0.055 | 0.627 | 0.627 |
| TCGA-COAD | GCH1 | 135 | 200 | 0.528 | 1.29e-06 | 1.03e-05 |
| TCGA-COAD | SCD | 135 | 200 | -0.135 | 0.22 | 0.289 |
| TCGA-COAD | MBOAT1 | 135 | 200 | 0.315 | 0.00394 | 0.0105 |
| TCGA-COAD | MBOAT2 | 135 | 200 | 0.127 | 0.253 | 0.289 |
| TCGA-COAD | GPX4 | 135 | 200 | -0.165 | 0.143 | 0.229 |
| TCGA-COAD | SLC7A11 | 135 | 200 | 0.342 | 0.00176 | 0.00703 |
| TCGA-COAD | DHODH | 135 | 200 | -0.169 | 0.128 | 0.229 |

## Primary coupling within right-sided CRC

| Cohort | Candidate | n | Spearman ρ | P | FDR |
|---|---|---:|---:|---:|---:|
| GSE39582 | AIFM2 | 224 | 0.125 | 0.0621 | 0.0828 |
| GSE39582 | GCH1 | 224 | 0.272 | 3.64e-05 | 0.000291 |
| GSE39582 | SCD | 224 | -0.026 | 0.701 | 0.801 |
| GSE39582 | MBOAT1 | 224 | 0.169 | 0.0111 | 0.0222 |
| GSE39582 | MBOAT2 | 224 | 0.128 | 0.0557 | 0.0828 |
| GSE39582 | GPX4 | 224 | -0.006 | 0.933 | 0.933 |
| GSE39582 | SLC7A11 | 224 | 0.174 | 0.00916 | 0.0222 |
| GSE39582 | DHODH | 224 | -0.203 | 0.00225 | 0.00899 |
| TCGA-COAD | AIFM2 | 200 | -0.231 | 0.00102 | 0.00272 |
| TCGA-COAD | GCH1 | 200 | 0.208 | 0.00311 | 0.00622 |
| TCGA-COAD | SCD | 200 | -0.172 | 0.015 | 0.02 |
| TCGA-COAD | MBOAT1 | 200 | 0.137 | 0.0522 | 0.0522 |
| TCGA-COAD | MBOAT2 | 200 | 0.375 | 4.47e-08 | 1.79e-07 |
| TCGA-COAD | GPX4 | 200 | -0.152 | 0.0316 | 0.0361 |
| TCGA-COAD | SLC7A11 | 200 | 0.196 | 0.00532 | 0.00851 |
| TCGA-COAD | DHODH | 200 | -0.379 | 3.18e-08 | 1.79e-07 |

## Primary AA-routing proxy high versus low

The split analysis is descriptive and complements the continuous correlation; it should not replace the continuous test.

| Cohort | Scope | Candidate | High n | Low n | High−low mean z | Welch P | FDR |
|---|---|---|---:|---:|---:|---:|---:|
| GSE39582 | all_side_known | AIFM2 | 283 | 283 | 0.046 | 0.588 | 0.659 |
| GSE39582 | all_side_known | GCH1 | 283 | 283 | 0.432 | 2.13e-07 | 1.32e-06 |
| GSE39582 | all_side_known | SCD | 283 | 283 | -0.037 | 0.659 | 0.659 |
| GSE39582 | all_side_known | MBOAT1 | 283 | 283 | 0.413 | 6.92e-07 | 1.85e-06 |
| GSE39582 | all_side_known | MBOAT2 | 283 | 283 | 0.173 | 0.0393 | 0.0629 |
| GSE39582 | all_side_known | GPX4 | 283 | 283 | -0.062 | 0.46 | 0.613 |
| GSE39582 | all_side_known | SLC7A11 | 283 | 283 | 0.425 | 3.3e-07 | 1.32e-06 |
| GSE39582 | all_side_known | DHODH | 283 | 283 | -0.216 | 0.0102 | 0.0204 |
| GSE39582 | right | AIFM2 | 142 | 82 | 0.312 | 0.0226 | 0.0453 |
| GSE39582 | right | GCH1 | 142 | 82 | 0.555 | 2.96e-05 | 0.000236 |
| GSE39582 | right | SCD | 142 | 82 | -0.108 | 0.43 | 0.491 |
| GSE39582 | right | MBOAT1 | 142 | 82 | 0.425 | 0.00133 | 0.00533 |
| GSE39582 | right | MBOAT2 | 142 | 82 | 0.209 | 0.165 | 0.22 |
| GSE39582 | right | GPX4 | 142 | 82 | -0.079 | 0.593 | 0.593 |
| GSE39582 | right | SLC7A11 | 142 | 82 | 0.415 | 0.00375 | 0.01 |
| GSE39582 | right | DHODH | 142 | 82 | -0.296 | 0.0309 | 0.0494 |
| TCGA-COAD | all_side_known | AIFM2 | 168 | 167 | -0.352 | 0.00119 | 0.00318 |
| TCGA-COAD | all_side_known | GCH1 | 168 | 167 | 0.276 | 0.0113 | 0.015 |
| TCGA-COAD | all_side_known | SCD | 168 | 167 | -0.185 | 0.09 | 0.09 |
| TCGA-COAD | all_side_known | MBOAT1 | 168 | 167 | 0.260 | 0.0172 | 0.0196 |
| TCGA-COAD | all_side_known | MBOAT2 | 168 | 167 | 0.527 | 9.53e-07 | 7.63e-06 |
| TCGA-COAD | all_side_known | GPX4 | 168 | 167 | -0.295 | 0.00681 | 0.0136 |
| TCGA-COAD | all_side_known | SLC7A11 | 168 | 167 | 0.278 | 0.0107 | 0.015 |
| TCGA-COAD | all_side_known | DHODH | 168 | 167 | -0.384 | 0.000394 | 0.00158 |
| TCGA-COAD | right | AIFM2 | 112 | 88 | -0.402 | 0.00407 | 0.0109 |
| TCGA-COAD | right | GCH1 | 112 | 88 | 0.371 | 0.0087 | 0.0174 |
| TCGA-COAD | right | SCD | 112 | 88 | -0.215 | 0.136 | 0.136 |
| TCGA-COAD | right | MBOAT1 | 112 | 88 | 0.227 | 0.119 | 0.136 |
| TCGA-COAD | right | MBOAT2 | 112 | 88 | 0.567 | 7e-05 | 0.00056 |
| TCGA-COAD | right | GPX4 | 112 | 88 | -0.302 | 0.0301 | 0.0401 |
| TCGA-COAD | right | SLC7A11 | 112 | 88 | 0.344 | 0.017 | 0.0271 |
| TCGA-COAD | right | DHODH | 112 | 88 | -0.514 | 0.000259 | 0.00104 |

## Interpretation rule

A candidate becomes interesting only if its direction is reasonably reproducible across cohorts and its continuous coupling is positive within right-sided tumors. A positive RNA association is still not proof of AA-induced defense or dependency.
The next step, if a candidate survives this screen, is a context-specific functional test rather than another expression-only expansion.

## Sensitivity analysis

The expanded AA-routing proxy results are saved in the machine-readable tables. Adding ALOX5/ALOX5AP is interpreted cautiously because these genes can be strongly contributed by myeloid cells in bulk CRC.

## Provenance

- GSE39582: {'n_samples_in_matrix': 585, 'matrix_rows_seen': 54675, 'matched_probe_rows': 27, 'missing_genes': [], 'probe_rows_used': {'AIFM2': 2, 'GCH1': 1, 'SCD': 4, 'MBOAT1': 1, 'MBOAT2': 2, 'GPX4': 1, 'SLC7A11': 3, 'DHODH': 3, 'PLA2G4A': 1, 'PTGS2': 2, 'PTGES': 2, 'ALOX5': 4, 'ALOX5AP': 1}, 'metadata_keys': ['Sample_channel_count', 'Sample_characteristics_ch1', 'Sample_characteristics_ch1.10', 'Sample_characteristics_ch1.11', 'Sample_characteristics_ch1.12', 'Sample_characteristics_ch1.13', 'Sample_characteristics_ch1.14', 'Sample_characteristics_ch1.15', 'Sample_characteristics_ch1.16', 'Sample_characteristics_ch1.17', 'Sample_characteristics_ch1.18', 'Sample_characteristics_ch1.19', 'Sample_characteristics_ch1.2', 'Sample_characteristics_ch1.20', 'Sample_characteristics_ch1.21', 'Sample_characteristics_ch1.22', 'Sample_characteristics_ch1.23', 'Sample_characteristics_ch1.24', 'Sample_characteristics_ch1.25', 'Sample_characteristics_ch1.26', 'Sample_characteristics_ch1.27', 'Sample_characteristics_ch1.28', 'Sample_characteristics_ch1.29', 'Sample_characteristics_ch1.3', 'Sample_characteristics_ch1.30', 'Sample_characteristics_ch1.31', 'Sample_characteristics_ch1.32', 'Sample_characteristics_ch1.33', 'Sample_characteristics_ch1.4', 'Sample_characteristics_ch1.5', 'Sample_characteristics_ch1.6', 'Sample_characteristics_ch1.7', 'Sample_characteristics_ch1.8', 'Sample_characteristics_ch1.9', 'Sample_contact_address', 'Sample_contact_city', 'Sample_contact_country', 'Sample_contact_department', 'Sample_contact_email', 'Sample_contact_institute', 'Sample_contact_laboratory', 'Sample_contact_name', 'Sample_contact_zip/postal_code', 'Sample_data_processing', 'Sample_data_row_count', 'Sample_description', 'Sample_description.10', 'Sample_description.11', 'Sample_description.12', 'Sample_description.13', 'Sample_description.14', 'Sample_description.15', 'Sample_description.2', 'Sample_description.3', 'Sample_description.4', 'Sample_description.5', 'Sample_description.6', 'Sample_description.7', 'Sample_description.8', 'Sample_description.9', 'Sample_extract_protocol_ch1', 'Sample_geo_accession', 'Sample_growth_protocol_ch1', 'Sample_hyb_protocol', 'Sample_label_ch1', 'Sample_label_protocol_ch1', 'Sample_last_update_date', 'Sample_molecule_ch1', 'Sample_organism_ch1', 'Sample_platform_id', 'Sample_scan_protocol', 'Sample_source_name_ch1', 'Sample_status', 'Sample_submission_date', 'Sample_supplementary_file', 'Sample_taxid_ch1', 'Sample_title', 'Sample_treatment_protocol_ch1', 'Sample_type'], 'probe_counts': {'PTGES': 2, 'ALOX5AP': 1, 'MBOAT2': 2, 'PTGS2': 2, 'SCD': 4, 'GPX4': 1, 'AIFM2': 2, 'MBOAT1': 1, 'DHODH': 3, 'PLA2G4A': 1, 'ALOX5': 4, 'SLC7A11': 3, 'GCH1': 1}, 'n_side_known_primary': 566, 'n_analyzed': 566, 'side_counts': {'left': 342, 'right': 224}}
- TCGA-COAD: {'path': 'C:\\Users\\21634\\Documents\\Codex\\2026-08-31\\1-janssen-kp-et-al-extrinsic\\work\\data\\tcga_aa_defense_candidate_expression_uqfpkm.tsv', 'source': 'cached GDC gene_expression/values response', 'n_cases': 458, 'n_genes': 13, 'missing_genes': [], 'n_analyzed': 335, 'side_counts': {'right': 200, 'left': 135}}
- Expression values were converted to within-cohort z-scores before score construction.
- Statistical unit: bulk tumor sample/case.
- FDR is controlled separately within each cohort/comparison family across the eight pre-specified candidates.

## Files

- `aa_defense_candidate_sample_scores.csv`
- `aa_defense_candidate_sidedness.csv`
- `aa_defense_candidate_proxy_splits.csv`
- `aa_defense_candidate_coupling.csv`
- `aa_defense_candidate_screen.png`
