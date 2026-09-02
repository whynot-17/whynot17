# Minimal CRC spatial PUFA/AA–neuronal–ferroptosis validation

## Frozen scope

- Dataset: GSE280315 P2CRC Visium HD, one stage III-B CRC section.
- Scores are within-section z-scored means of log1p(CPM+1) expression; they are transcriptional proxies, not lipid concentrations or measured flux.
- Neoplasm bins touching any non-neoplasm pathologist-labeled bin in the 8-neighbor grid are labeled `invasive_front_proxy`; other Neoplasm bins are `tumor_core`.
- `neural_containing` means at least two neuronal genes detected in the bin.
- The fixed ferroptosis-stress score overlaps with PUFA pressure through ACSL4/LPCAT3/ALOX15; `ferroptosis_stress_nonoverlap_score` (TFRC/POR/CYB5R1) is reported as a sensitivity analysis.
- Requested models: `neuronal_score ~ pufa_aa_pressure_score + tumor_region` and `ferroptosis_stress_score ~ pufa_aa_pressure_score * neural_containing + tumor_region`, with HC3 standard errors.

## Data audit

- H5 bins: **545,913**; analyzed non-zero-UMI bins: **545,850**; annotation overlap: **545,913**; coordinate overlap: **545,913**.
- Stage: **III-B**; early-versus-late stage comparison: **not estimable in this one-section run**.
- Neural-containing bins (≥2 neuronal genes): **66** (0.01%).

## Genes used

- pufa_incorporation: ACSL4, LPCAT3, AGPAT3
- aa_routing: PLA2G4A, PTGS2, PTGES, ALOX5, ALOX15
- neuronal: ELAVL3, ELAVL4, SNAP25, UCHL1, SYP, PHOX2B
- ferroptosis_stress: TFRC, ACSL4, LPCAT3, ALOX15, POR, CYB5R1
- ferroptosis_defense: SLC7A11, GPX4, AIFM2, GCH1
- ferroptosis_stress_nonoverlap: TFRC, POR, CYB5R1

## Region summary

| region | bins | PUFA/AA pressure | neuronal | ferroptosis stress | non-overlap stress | defense |
|---|---:|---:|---:|---:|---:|---:|
| adjacent_non_tumor | 65064 | -0.062 | -0.022 | -0.099 | -0.129 | -0.101 |
| connective_tissue | 107327 | -0.095 | -0.019 | -0.143 | -0.195 | -0.107 |
| invasive_front_proxy | 4437 | -0.036 | -0.019 | -0.049 | -0.065 | -0.035 |
| muscularis_smooth_muscle | 51234 | -0.110 | -0.021 | -0.167 | -0.224 | -0.140 |
| outside | 1 | -0.143 | -0.025 | -0.197 | -0.257 | -0.176 |
| tumor_core | 310402 | 0.019 | -0.020 | 0.066 | 0.100 | 0.026 |
| vascular | 7385 | -0.128 | -0.024 | -0.182 | -0.242 | -0.161 |

## Descriptive correlations

| analysis | subset | n | Spearman rho | p |
|---|---|---:|---:|---:|
| PUFA/AA pressure vs neuronal score | all tissue bins | 545850 | 0.0101 | 8.426e-14 |
| PUFA/AA pressure vs neuronal score | neural-containing bins | 66 | 0.2728 | 0.02669 |
| PUFA/AA pressure vs ferroptosis stress | neural-containing bins | 66 | 0.593 | 1.556e-07 |
| PUFA/AA pressure vs ferroptosis stress (non-overlap sensitivity) | neural-containing bins | 66 | -0.0462 | 0.7126 |

## Requested exploratory models

| model | term | estimate | HC3 p | n | R² |
|---|---|---:|---:|---:|---:|
| neuronal_score_model | Intercept | -0.02181 | 0 | 545850 | 0.0003008 |
| neuronal_score_model | C(tumor_region)[T.connective_tissue] | 0.002788 | 3.524e-20 | 545850 | 0.0003008 |
| neuronal_score_model | C(tumor_region)[T.invasive_front_proxy] | 0.003084 | 0.004492 | 545850 | 0.0003008 |
| neuronal_score_model | C(tumor_region)[T.muscularis_smooth_muscle] | 0.0007223 | 0.04748 | 545850 | 0.0003008 |
| neuronal_score_model | C(tumor_region)[T.outside] | -0.002739 | 8.172e-06 | 545850 | 0.0003008 |
| neuronal_score_model | C(tumor_region)[T.tumor_core] | 0.00164 | 4.612e-13 | 545850 | 0.0003008 |
| neuronal_score_model | C(tumor_region)[T.vascular] | -0.001977 | 3.703e-08 | 545850 | 0.0003008 |
| neuronal_score_model | pufa_aa_pressure_score | 0.002186 | 2.849e-09 | 545850 | 0.0003008 |
| ferroptosis_stress_interaction_model | Intercept | -0.06152 | 0 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | C(tumor_region)[T.connective_tissue] | -0.02496 | 6.514e-122 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | C(tumor_region)[T.invasive_front_proxy] | 0.03385 | 8.446e-16 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | C(tumor_region)[T.muscularis_smooth_muscle] | -0.03999 | 2.522e-291 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | C(tumor_region)[T.outside] | -0.05012 | 0.003326 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | C(tumor_region)[T.tumor_core] | 0.1164 | 0 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | C(tumor_region)[T.vascular] | -0.04369 | 5.362e-201 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | pufa_aa_pressure_score | 0.5983 | 0 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | neural_containing | 0.06669 | 0.03757 | 545850 | 0.3289 |
| ferroptosis_stress_interaction_model | pufa_aa_pressure_score:neural_containing | 0.0323 | 0.8181 | 545850 | 0.3289 |
| ferroptosis_stress_nonoverlap_interaction_model | Intercept | -0.1206 | 0 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | C(tumor_region)[T.connective_tissue] | -0.06097 | 3.277e-249 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | C(tumor_region)[T.invasive_front_proxy] | 0.06068 | 8.176e-17 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | C(tumor_region)[T.muscularis_smooth_muscle] | -0.08873 | 0 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | C(tumor_region)[T.outside] | -0.1165 | 0.1613 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | C(tumor_region)[T.tumor_core] | 0.2178 | 0 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | C(tumor_region)[T.vascular] | -0.1036 | 0 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | pufa_aa_pressure_score | 0.138 | 0 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | neural_containing | 0.07036 | 0.1636 | 545850 | 0.08294 |
| ferroptosis_stress_nonoverlap_interaction_model | pufa_aa_pressure_score:neural_containing | -0.2492 | 0.0003204 | 545850 | 0.08294 |

## Interpretation guardrails

This is a one-section spatial screen. Bin-level p-values are descriptive and can be strongly anti-conservative because neighboring bins are not independent biological replicates. A positive association would motivate multi-patient spatial replication and/or direct neuronal and lipid-peroxidation measurements; it does not establish ENS depletion, AA flux, ferroptosis, causality or stage progression.
