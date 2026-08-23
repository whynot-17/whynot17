# Source-data map

## Figure 1

- `figure1_primary_screen.csv`: all 267 chemicals from the frozen primary GeneCards-disorders screen, including screen rank, enrichment odds ratio, BH-FDR, CRC overlap count and degree-matched permutation values
- `figure1_panelA_workflow_nodes.csv`: discovery-universe workflow nodes
- `figure1_panelC_translation_links.csv`: MiNP/DINP-to-MCOP translation relations
- `figure1_panelD_study_roadmap.csv`: discovery, human biomonitoring and biological-interpretation stages

## Figure 2

- `figure2_primary_python_vs_r.csv`: pooled R and independent Python survey estimates
- `figure2_loco.csv`: pooled and seven leave-one-cycle-out estimates
- `figure2_per_cycle.csv`: seven cycle-specific estimates, case counts and convergence status
- `figure2_cycle_interaction.csv`: global MCOP-by-cycle interaction test

## Figure 3

- `figure3_sensitivity.csv`: prespecified population, diagnosis-timing, tail and creatinine analyses
- `figure3_coexposure.csv`: pairwise co-exposure-adjusted MCOP estimates
- `figure3_rcs_curve_with_ci.csv`: survey-weighted restricted-cubic-spline curve with pointwise 95% CI
- `figure3_rcs_curve_with_ci_metadata.json`: knots, reference, design degrees of freedom and Wald-test metadata
- `figure3_weighted_quartiles.csv`: survey-weighted quartile estimates and trend test

The spline display is restricted to the survey-weighted 5th–95th percentiles. Knots are the weighted 5th, 35th, 65th and 95th percentiles; the reference is the weighted median (7.9 ng/mL).

## Figure 4

- `figure4_bulk_scores.csv` and `figure4_bulk_sample_manifest.csv`: TCGA score observations and matched-patient identifiers
- `figure4_tcga_paired_summary.csv`: TCGA paired summary statistics
- `figure4_census_donor_scores.csv`: Census donor-compartment paired pseudobulk scores
- `figure4_census_paired_summary.csv`: Census compartment-level paired summaries
- `figure4_gse144735_scores.csv`: six matched epithelial donor scores
- `figure4_gse144735_paired_summary.csv`: GSE144735 paired summary statistics

All displayed points come from these tables. The only positional jitter is deterministic and affects vertical display coordinates, not data values.

## Figure 5

- `figure5_evidence_nodes.csv`: discovery, human-biomonitoring and CRC-state evidence nodes with evidence status
- `figure5_evidence_links.csv`: observed solid relations and the hypothetical dashed mechanistic bridge
