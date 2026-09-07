# PGE2–PTGER4 robustness re-analysis

## Scope

This re-analysis keeps `analysis/dinp_crc_pge2_ptger4_axis/` unchanged and tests whether its source → receiver conclusion depends on score definition, PTGES3, pooled cell-level z-scoring, sparse donor groups, winner-takes-all ranking, or the degenerate PTGER4-high threshold. No CellChat, spatial, machine learning, docking, MD, clustering, or new macrophage annotation was performed.

## Inputs and primary design

- Local GSE144735 CRC cache: 27,414 cells and 6 donors (KUL01, KUL19, KUL21, KUL28, KUL30, KUL31).
- Primary unit: donor × `analysis_group` pseudobulk, using mean normalized log1p expression, detection fraction, expressing-cell mean, and `log1p(mean(expm1(log1p expression)))` pseudobulk expression.
- Primary source filter: tumor donor × group with **≥20 cells**. Threshold sensitivity: ≥10, ≥20, and ≥30 cells.
- Scores A–C are donor-group gene-wise z-score means. Score D is the bottleneck-aware geometric mean: `geometric_mean(PLA2G4A/q95, max(PTGS1,PTGS2)/q95, max(PTGES,PTGES2)/q95)`, with each component clipped to [0,1]. Score D therefore requires upstream, cyclooxygenase, and terminal synthase support together.
- Score C is the original six-gene definition and includes PTGES3 only as a sensitivity analysis. Scores A, B, and D exclude PTGES3.
- The prior pooled cell-level global z-score is retained in the per-donor-group tables as `cell_global_sensitivity_mean`; it is not the primary evidence.

## Required questions

1. **Is mast-cell first robust?** Pooled mast cells remain a high signal, but donor tumor mast-cell counts are: **KUL01:8, KUL19:31, KUL21:27, KUL28:10, KUL30:1, KUL31:1**. At n≥20, mast-cell source ranks are therefore coverage- and sparsity-sensitive; see `mast_cell_donor_counts.csv`, `cell_count_threshold_sensitivity.csv`, and the count-versus-score figure.
2. **Who is highest without PTGES3?** At n≥20, the highest mean donor-group source is **fibroblast_like** under Score B and **C1QC_like_TAM** under Score D. Score C is the PTGES3-inclusive sensitivity; threshold-specific donor calls are in `source_rank_consistency.csv`.
3. **Most stable score:** **Score_B** has the highest mean pairwise Spearman concordance of donor-group ranks at n≥20 (0.767); full concordance values are in `final_robustness_summary.csv`.
4. **Does SPP1-like TAM remain high?** Score A: **median 2.00; top-2 1/1 eligible (1/6 total)**; Score B: **median 2.00; top-2 1/1 eligible (1/6 total)**; Score D: **median 3.00; top-2 0/1 eligible (0/6 total)**. The eligible-donor denominator is essential: at n≥20, SPP1-like TAM has adequate coverage in only one tumor donor, so its apparent high rank cannot be generalized to six donors. This directly tests the PTGES3 removal and bottleneck score.
5. **Does C1QC remain PTGER4-enriched?** The pooled receiver ranking is **C1QC_like_TAM**; C1QC's pooled PTGER4 mean is **0.191** with detection **0.232**. Donor pseudobulk median PTGER4 rank is **2.00**; C1QC is higher than SPP1 in **2/4** donors with both subtypes (**2/6** overall), and higher than the other-subtype mean in **3/5** valid donors. Full Δ1/Δ2 contrasts are in `c1qc_ptger4_donor_contrasts.csv` and pooled values are in `ptger4_receiver_pooled_summary.csv`.
6. **Why are donors discordant?** The tables separate rank-rule effects from sparse coverage. Mast-cell donor counts are low in several donors; C1QC/SPP1 state counts also vary. Dropout is reflected by detection and expressing-cell means. Remaining differences after count filtering are consistent with biological heterogeneity, but six donors cannot distinguish it decisively.
7. **Axis evidence level:** source = **WEAK SUPPORT**; receiver = **WEAK SUPPORT**; combined descriptive axis = **WEAK SUPPORT**. Axis 1/2/3 are ranking sensitivities only and are not communication probabilities.

## Interpretation

`PTGER4-high fraction` is not reported as an independent endpoint here because the previous macrophage q75 was zero and it collapsed to detection fraction. The receiver analysis instead uses PTGER4 mean, detection fraction, expressing-cell mean, pseudobulk PTGER4, and donor rank.

### Recommended interpretation

The data **suggests only weakly** a potential PGE2–PTGER4 macrophage response axis at the expression level. The pooled C1QC PTGER4 enrichment is reproducible as a pooled description, while donor-level source and receiver ranks must be read with sparse-group coverage and heterogeneity visible. This does not demonstrate extracellular PGE2 production, receptor binding, spatial proximity, or functional communication.

### Recommended next step

larger single-cell cohort replication plus spatial validation before extending the axis.

## Files

- `pge2_score_definitions.csv`: exact Score A–D definitions and formula.
- `pseudobulk_donor_group_expression.csv`: donor-group means, detection, expressing-cell means, and pseudobulk expression.
- `pge2_score_A_by_donor_group.csv` through `pge2_score_D_by_donor_group.csv`: score tables with old global-z sensitivity.
- `source_rank_consistency.csv`, `source_score_contrasts.csv`, `cell_count_threshold_sensitivity.csv`: source robustness and SPP1 contrasts.
- `ptger4_receiver_by_donor.csv`, `ptger4_receiver_pooled_summary.csv`, `c1qc_ptger4_donor_contrasts.csv`: receiver metrics, pooled description, and donor contrasts.
- `pge2_ptger4_axis_sensitivity.csv`: three descriptive source → C1QC axis scores.
- `final_robustness_summary.csv`: numeric evidence ledger and final support grades.
- `figures/`: PNG, PDF, and SVG figures requested in the protocol.
