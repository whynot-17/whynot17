# Phase 2G — DINP/MCOP–CRC epithelial state convergence and regulatory anchoring

Generated: 2026-08-24T11:13:05.045813+00:00

## Frozen analysis boundaries

- Census release: **2025-11-08**; primary data filter enforced.
- Expression access path: **official source H5AD, streamed locally**.
- Primary dataset selected from Phase 2F matched epithelial audit: **16023185-de21-4c0d-a9c8-73abdd52d142**; paired donors queried=36; epithelial cell types=8.
- PPAR/NR core fixed as: **PPARA, PPARD, PPARG, NR1I2, NR1I3, NR1H2, NR1H3**; no result-driven gene editing.
- State universe: **27** programs; gene-set sizes present in queried expression: `{"EMT": 199, "E2F_targets": 200, "G2M_checkpoint": 199, "MYC_targets_V1": 200, "MYC_targets_V2": 58, "Hypoxia": 200, "Inflammatory_response": 200, "TNF_NFkB": 200, "IL6_JAK_STAT3": 87, "ROS": 49, "OXPHOS": 200, "Fatty_acid_metabolism": 158, "Cholesterol_homeostasis": 74, "Apoptosis": 160, "p53": 200, "UPR": 113, "Glycolysis": 200, "IFN_alpha": 97, "IFN_gamma": 200, "TGF_beta": 54, "WNT_beta_catenin": 42, "intestinal_epithelial_differentiation": 17, "stemness": 10, "secretory_differentiation": 10, "enterocyte_differentiation": 10, "goblet_program": 8, "stress_like_epithelial": 10}`.
- PPAR-low/high is defined at cell level as bottom/top quartile; tercile and median labels are retained as sensitivity labels. Inference uses donor-level pseudobulk.
- Expression DE is targeted to the frozen state/regulator gene universe, not genome-wide. This limitation is explicit and is not called a genome-wide DE result.

## Primary answer

- Paired epithelial PPAR/NR tumor-normal result: [{'median_delta_tumor_minus_normal': -0.4186007613938946, 'p_value': 4.2922329157590866e-07}].
- Directly supported candidate anchors under the frozen evidence tags: **RELA, STAT3**.
- The most defensible interpretation remains a CRC epithelial disease-state convergence. It does not establish that DINP/MCOP causes the state or that any regulator mediates the epidemiologic association.

## State discovery

All prespecified state programs were scored before ranking. The state-correlation table reports donor-level Spearman associations with PPAR/NR; no EMT, stemness, inflammatory or metabolic state was selected in advance as the expected winner.

## Regulatory activity and anchor boundary

DoRothEA confidence levels A–C with decoupler ULM were used for the listed candidate regulators. Activity shifts are descriptive donor-level contrasts. Anchor tiers are evidence tags, not a subjective numeric total score.

## Opposite-direction compartment result

The Phase 2F myeloid PPAR/NR increase is retained as a visible localization result; this Phase 2G epithelial analysis does not reinterpret bulk tissue or erase the compartment contrast.

## Final verdict

**PARTIALLY** — the analysis tests whether epithelial PPAR/NR suppression defines a reproducible CRC state and whether a regulatory bridge is plausible. A positive disease-state convergence can be supported, but the DINP/MCOP → PPAR/NR arrow remains untested; GSE144735 is small and directional, and no causal perturbation is asserted.

## Output files

- `mcop_phase2g_epithelial_state_scores.csv`
- `mcop_phase2g_ppar_low_high_de.csv`
- `mcop_phase2g_pathway_state_scores.csv`
- `mcop_phase2g_state_correlations.csv`
- `mcop_phase2g_donor_level_validation.csv`
- `mcop_phase2g_subtype_localization.csv`
- `mcop_phase2g_tumor_ppar_interaction.csv`
- `mcop_phase2g_regulator_activity.csv`
- `mcop_phase2g_regulatory_anchor_ranking.csv`
- `mcop_phase2g_external_validation.csv`
- `mcop_phase2g_bridge_evidence_table.csv`
- `mcop_phase2g_donor_state_pseudobulk.csv`
