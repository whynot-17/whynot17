# PGE2–PTGER4 robustness re-analysis

## Scope

This re-analysis keeps `analysis/dinp_crc_pge2_ptger4_axis/` unchanged and tests whether its source → receiver conclusion depends on score definition, PTGES3, pooled cell-level z-scoring, sparse donor groups, winner-takes-all ranking, or the degenerate PTGER4-high threshold. No CellChat, spatial, machine learning, docking, MD, clustering, or new macrophage annotation was performed.

## Inputs and primary design

- Local GSE178341 CRC cache: 370,115 cells and 62 donors (C103, C104, C105, C106, C107, C109, C110, C111, C112, C113, C114, C115, C116, C118, C119, C122, C123, C124, C125, C126, C129, C130, C132, C133, C134, C135, C136, C137, C138, C139, C140, C142, C143, C144, C145, C146, C147, C149, C150, C151, C152, C153, C154, C155, C156, C157, C158, C159, C160, C161, C162, C163, C164, C165, C166, C167, C168, C169, C170, C171, C172, C173).
- Primary unit: donor × `analysis_group` pseudobulk, using mean normalized log1p expression, detection fraction, expressing-cell mean, and `log1p(mean(expm1(log1p expression)))` pseudobulk expression.
- Primary source filter: tumor donor × group with **≥20 cells**. Threshold sensitivity: ≥10, ≥20, and ≥30 cells.
- Scores A–C are donor-group gene-wise z-score means. Score D is the bottleneck-aware geometric mean: `geometric_mean(PLA2G4A/q95, max(PTGS1,PTGS2)/q95, max(PTGES,PTGES2)/q95)`, with each component clipped to [0,1]. Score D therefore requires upstream, cyclooxygenase, and terminal synthase support together.
- Score C is the original six-gene definition and includes PTGES3 only as a sensitivity analysis. Scores A, B, and D exclude PTGES3.
- The prior pooled cell-level global z-score is retained in the per-donor-group tables as `cell_global_sensitivity_mean`; it is not the primary evidence.

## Required questions

1. **Is mast-cell first robust?** Pooled mast cells remain a high signal, but donor tumor mast-cell counts are: **C103:28, C104:6, C105:6, C106:2, C107:63, C109:26, C110:7, C111:26, C112:63, C113:24, C114:47, C115:2, C116:5, C118:11, C122:15, C123:60, C124:375, C125:120, C126:19, C129:105, C130:21, C132:81, C133:53, C134:5, C135:39, C136:57, C137:15, C138:51, C139:151, C140:33, C142:92, C143:24, C144:22, C145:22, C146:25, C147:47, C149:90, C150:7, C152:3, C153:100, C154:14, C155:58, C156:29, C157:49, C158:21, C159:39, C160:33, C161:11, C162:225, C163:38, C165:21, C166:5, C168:8, C169:7, C170:60, C171:68, C172:69, C173:12**. At n≥20, mast-cell source ranks are therefore coverage- and sparsity-sensitive; see `mast_cell_donor_counts.csv`, `cell_count_threshold_sensitivity.csv`, and the count-versus-score figure.
2. **Who is highest without PTGES3?** At n≥20, the highest mean donor-group source is **mast_cells** under Score B and **C1QC_like_TAM** under Score D. Score C is the PTGES3-inclusive sensitivity; threshold-specific donor calls are in `source_rank_consistency.csv`.
3. **Most stable score:** **Score_B** has the highest mean pairwise Spearman concordance of donor-group ranks at n≥20 (0.741); full concordance values are in `final_robustness_summary.csv`.
4. **Does SPP1-like TAM remain high?** Score A: **median 5.00; top-2 2/56 eligible (2/62 total)**; Score B: **median 5.00; top-2 9/56 eligible (9/62 total)**; Score D: **median 3.00; top-2 20/56 eligible (20/62 total)**. The eligible-donor denominator is essential: at n≥20, SPP1-like TAM is covered unevenly across donors, so its apparent high rank must be interpreted with the eligible-donor denominator and cannot be generalized to all 62 donors. This directly tests the PTGES3 removal and bottleneck score.
5. **Does C1QC remain PTGER4-enriched?** The pooled receiver ranking is **resident_like**; C1QC's pooled PTGER4 mean is **0.225** with detection **0.256**. Donor pseudobulk median PTGER4 rank is **2.00**; C1QC is higher than SPP1 in **32/62** donors with both subtypes (**32/62** overall), and higher than the other-subtype mean in **34/62** valid donors. Full Δ1/Δ2 contrasts are in `c1qc_ptger4_donor_contrasts.csv` and pooled values are in `ptger4_receiver_pooled_summary.csv`.
6. **Why are donors discordant?** The tables separate rank-rule effects from sparse coverage. Mast-cell donor counts are low in several donors; C1QC/SPP1 state counts also vary. Dropout is reflected by detection and expressing-cell means. Remaining differences after count filtering are consistent with biological heterogeneity, but 62 donors cannot distinguish it decisively.
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
