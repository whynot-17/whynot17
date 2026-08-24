# Frozen figure-build lock

Date: 2026-08-24
Branch: `phase2f-compartment-external-replication`
Status: FIGURE CONSTRUCTION ONLY — NO NEW EXPLORATORY ANALYSIS

This file implements the frozen five-figure architecture from `work/environmental_crc_final_analysis_freeze_20260824.md`.

## Figure 1 — Data-first molecular discovery
Purpose: establish unbiased discovery and preserve the distinction between database association and causality.

Main visual logic:
1. CTD human chemical–gene associations + GeneCards CRC genes.
2. Evidence filtering / enrichment / degree-matched permutation.
3. Highlight DINP/MiNP-related molecular nomination without implying that MCOP was a CTD molecular hit.

Do not use causal arrows from CTD association to CRC.

## Figure 2 — Outcome-blinded actionability prioritization
Primary message:
> 267 environmental chemicals were filtered by prespecified human-actionability gates to 87 eligible chemical–biomarker mappings, corresponding to 15 unique NHANES biomarker tests.

Frozen sequential counts:
- 267 total
- 259 E-valid
- 135 E+X interpretable
- 134 +B biomarker available
- 127 +D detectable
- 124 +C adequate coverage
- 87 +T human-testable
- 27 strict-eligible

Visual architecture:
A. left-to-right or radial attrition architecture using the true sequential counts above.
B. a 15-test biomarker ring/grid showing biomarker, biological matrix, eligible chemical count, and mapping confidence.
C. identity callout showing MiNP / DINP parent / MCOP as distinct entities.

Mandatory wording:
`87 eligible chemical–biomarker mappings corresponded to 15 unique NHANES biomarker tests.`

Forbidden shorthand:
`87 chemicals collapsed into 15 equivalent exposures.`

Shared biomarker mappings must not imply parent-chemical equivalence.

## Figure 3 — Systematic human screening
Primary message:
> Across the 15 unique biomarker tests, two passed BH-FDR<0.05, and the DINP-related MCOP test subsequently showed the strongest overall robustness profile.

Use all 15 tests in one integrated effect landscape rather than a generic forest plot if possible.

Frozen FDR-supported tests:
- LBXPFHS: OR 0.6244, BH-FDR 0.02190, Tier B
- URXCOP / MCOP: OR 1.2455, BH-FDR 0.02484, Robust Tier A

MCOP must not be visually selected before the full 15-test screen is displayed.

## Figure 4 — Robustness / heterogeneity / exposure shape
Primary message:
> MCOP remained positive under the prespecified robustness audit and LOCO analyses despite significant between-cycle heterogeneity.

Show robustness fingerprint, LOCO, cycle heterogeneity, and exposure-shape/sensitivity results without hiding the discordant 2011–2012 cycle.

Do not call seven cycles independent replications.

## Figure 5 — CRC epithelial disease-state convergence
Status: YELLOW — disease-state convergence, not causal mechanism.

Primary message:
> CRC shows epithelial-state-specific PPAR/nuclear-receptor and metabolic-differentiation remodeling, while the DINP-related exposure-to-epithelial-state bridge remains untested.

Required elements:
- frozen PPAR/NR core down in paired CRC epithelium
- independent KEGG/Reactome/Hallmark/regulon definitions predominantly reproduce the down direction
- enterocyte-like state: PPAR/NR down
- secretory-like state: modestly up
- inflammatory/stress remodeling displayed as a parallel CRC disease-state program, not downstream of PPAR loss
- DINP/MiNP environmental bridge drawn dashed/dotted and labeled candidate / untested

Forbidden arrow:
`PPAR/NR down -> RELA/STAT3 up`

## Figure-wide aesthetic lock
- High-density journal-style information architecture; avoid default bioinformatics box-and-arrow layouts.
- Main conclusion first; secondary QC details visually subordinate.
- Restrained palette and consistent visual grammar across all five figures.
- Real data points/effect sizes should dominate where available.
- Dashed/dotted styling reserved for hypotheses or untested links.
- Database associations, computed results, observational associations, and hypotheses must be visually distinguishable.

## Source locks
Figure 2 source files:
- `outputs/environmental_crc_267_actionability_flow.csv`
- `outputs/environmental_crc_267_human_testable_candidates.csv`
- `outputs/environmental_crc_267_actionability_matrix_v2.csv`

Figure 3 source file:
- `outputs/environmental_crc_15axis_robustness_scorecard.csv`

Figure 5 source files:
- `outputs/figure5_final_report.md`
- `outputs/figure5_ppar_definition_comparison.csv`
- `outputs/figure5_within_state_ppar_analysis.csv`
- `outputs/figure5_state_correlation_matrix.csv`
- `outputs/figure5_evidence_tier_lock.csv`

No values may be manually substituted from memory when an authoritative frozen source exists.
