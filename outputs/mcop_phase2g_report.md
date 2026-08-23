# Phase 2G — DINP/MCOP–CRC epithelial state convergence and regulatory anchoring

Generated: 2026-08-23T18:49:04.552035+00:00

## Runtime status

**PARTIAL / RUNTIME-BLOCKED:** the pinned WSL CELLxGENE Census raw-expression query repeatedly terminated at the TileDB `query.X('raw').tables()` stage, including a validated 9-gene minimal query and a single-donor probe. No Census primary state/DE/regulator result is claimed from this fallback.

- Existing Phase 2F donor-level epithelial output retained for dataset `16023185-de21-4c0d-a9c8-73abdd52d142`.
- GSE144735 full target-universe state scores and existing TCGA paired output were recomputed locally as external support only; the PPAR/NR external row reuses the frozen Phase 2F score definition for exact comparability.
- The 7-gene PPAR/NR core was not redefined; no cell-level PPAR-low/high state was fabricated from donor-level rows.

## What is and is not concluded

- Retained Phase 2F paired epithelial result: PPAR/NR median tumor-minus-normal Δ = −0.419 (36 donors; P = 4.29×10⁻⁷); RELA/STAT3 Δ = +1.167 (P = 1.08×10⁻⁷). Myeloid PPAR/NR remains opposite-direction (Δ = +0.610; 35 donors; P = 7.97×10⁻⁹).
- GSE144735 directionally concordant but underpowered: PPAR/NR median Δ = −0.312 (6 paired patients; P = 0.6875), using the frozen Phase 2F core score.

The available Phase 2F evidence continues to support an epithelial PPAR/NR disease-state signal, but Phase 2G cannot answer the prespecified unbiased state-discovery, tumor×PPAR interaction, or DoRothEA regulator-activity questions until primary Census target-universe expression is readable. Therefore the final verdict is **PARTIALLY**, not YES.

## Outputs

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
