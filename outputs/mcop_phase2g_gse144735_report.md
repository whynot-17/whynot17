# Phase 2G alternative primary — GSE144735 complete-matrix epithelial state analysis

Generated: 2026-08-24T00:33:28.485398+00:00

## Data path

- Primary alternative dataset: GSE144735 GEO processed natural-log TPM matrix; live CELLxGENE Census/TileDB raw query was not used.
- Cells in epithelial annotation: **6,168**; primary tumor/normal epithelial cells: **3,356**; paired donors: **6**.
- The matrix was streamed row-wise and only the frozen state/regulator gene universe was retained; the complete raw/processed matrix is not committed to the repository.
- Primary contrast excludes Border epithelial cells; Border cells are retained only in the cell audit and subtype localization.

## Frozen boundaries

- PPAR/NR core: **PPARA, PPARD, PPARG, NR1I2, NR1I3, NR1H2, NR1H3**; no result-driven gene editing.
- Unbiased state universe: **27** programs; present gene counts: `{"EMT": 198, "E2F_targets": 195, "G2M_checkpoint": 189, "MYC_targets_V1": 193, "MYC_targets_V2": 57, "Hypoxia": 194, "Inflammatory_response": 199, "TNF_NFkB": 198, "IL6_JAK_STAT3": 87, "ROS": 47, "OXPHOS": 184, "Fatty_acid_metabolism": 156, "Cholesterol_homeostasis": 73, "Apoptosis": 160, "p53": 193, "UPR": 108, "Glycolysis": 198, "IFN_alpha": 95, "IFN_gamma": 196, "TGF_beta": 54, "WNT_beta_catenin": 42, "intestinal_epithelial_differentiation": 17, "stemness": 10, "secretory_differentiation": 9, "enterocyte_differentiation": 10, "goblet_program": 7, "stress_like_epithelial": 10}`.
- Cell-level PPAR-low/high is defined by the bottom/top quartile of the fixed 7-gene PPAR/NR score; donor-level `PPAR_NR_score` validation reuses the frozen Phase 2F nine-gene pseudobulk score (7 core genes plus RELA/STAT3 denominator).
- State and regulator expression comes from the complete processed log-TPM matrix; donor-level inference uses donor-by-group mean expression because GEO supplies processed expression rather than raw full-matrix counts.
- Targeted state/regulator expression is not a genome-wide DE claim.

## Primary alternative results

- PPAR/NR paired result: `[{'feature_type': 'state', 'feature': 'PPAR_NR_score', 'analysis_level': 'GSE144735 donor-level paired epithelial analysis', 'n_paired_donors': 6, 'mean_delta_tumor_minus_normal': -0.1057725637355198, 'median_delta_tumor_minus_normal': -0.31168562963879776, 'p_value': 0.6875, 'direction': 'down', 'donor_consistency': 0.6666666666666666, 'BH_FDR': 0.6875}]`.
- State programs with estimable tumor-normal paired summaries: **28**.
- Largest absolute state contrast: **enterocyte_differentiation (median Δ=-0.461; BH-FDR=0.384)**; no state is called discovery-positive without multiplicity control.
- DoRothEA regulator summaries: **14** regulators; A-C weighted ULM with minimum five observed targets.
- Minimum regulator-activity BH-FDR across tested summaries: **0.068**; minimum tumor×PPAR interaction BH-FDR: **0.943**. These small-sample values are descriptive, not mechanistic confirmation.
- This alternative primary analysis completes the prespecified state, interaction, subtype, and regulator modules for one independent dataset, but the six-patient sample is underpowered for stable inference and cannot supply Census-style leave-one-dataset-out validation.

## Interpretation

The result tests whether DINP-axis molecular candidates converge on a CRC epithelial state; it does not establish DINP/MCOP causality or mediation. The Phase 2F Census epithelial/myeloid compartment result remains the prior validated support, and the myeloid opposite-direction result is not erased.

## Verdict

**PARTIALLY — alternative primary analysis completed, but independent-dataset stability and causal exposure-to-state direction remain unresolved.**

## Outputs

- `mcop_phase2g_gse144735_cell_state_scores.csv`
- `mcop_phase2g_gse144735_ppar_low_high_de.csv`
- `mcop_phase2g_gse144735_pathway_state_scores.csv`
- `mcop_phase2g_gse144735_state_correlations.csv`
- `mcop_phase2g_gse144735_donor_level_validation.csv`
- `mcop_phase2g_gse144735_subtype_localization.csv`
- `mcop_phase2g_gse144735_tumor_ppar_interaction.csv`
- `mcop_phase2g_gse144735_regulator_activity.csv`
- `mcop_phase2g_gse144735_regulatory_anchor_ranking.csv`
- `mcop_phase2g_gse144735_external_validation.csv`
- `mcop_phase2g_gse144735_bridge_evidence_table.csv`
- `mcop_phase2g_gse144735_donor_state_pseudobulk.csv`
