# Frozen-primary directional fork scan

This audit reads only the frozen 29 x 14 primary result grid (`406` fixed interaction family). It does not read robustness results, refit models, recalculate FDR, or perform mechanism, literature, figure, spline, quartile, or phenome-wide analyses.

## Pre-specified directional rule

For each exposure, endpoint A was required to satisfy `male beta > 0` and `female beta <= 0`. Endpoint B was required to satisfy `female beta > 0` and `male beta <= 0`. Every A x B combination was retained. Because the frozen interaction estimand is `female beta - male beta`, A additionally requires a negative interaction coefficient (male-enhanced) and B a positive interaction coefficient (female-enhanced); these direction checks were recorded, not used to cherry-pick rows.

## Results

- The frozen grid contained 406/406 estimable rows, 29 exposures, and 14 outcomes.
- The sex-specific sign rules generated 88 A -> B combinations across 16 exposures; 13 exposures generated no such combination.
- Requiring the corresponding pooled interaction directions leaves 80 interaction-supported directional pairs across 15 exposures; 14 exposures have no interaction-supported pair. The remaining 8 sign-only combinations are retained in the CSV with `strict_directional_pair = False` and are not treated as reciprocal candidates.
- No interaction-supported pair had both endpoints at fixed-406 BH q < 0.05 (0); 4 supported pairs had at least one endpoint at q < 0.05. These q-values are the existing fixed-family values, not a new pair-level FDR.
- URXP25 and URXUBA are not strict reciprocal candidates because neither has a male-only-positive endpoint A in the frozen primary grid.

The 4 pairs with at least one endpoint at the existing fixed-406 q < 0.05 threshold are listed descriptively below; this is not a new candidate-ranking rule:

- LBXPFNA: myocardial_infarction (male-enhanced endpoint, interaction beta -0.335385, q 0.0209695) -> any_cancer_history (female-enhanced endpoint, interaction beta 0.076549, q 0.6310196).
- LBXPFNA: myocardial_infarction (male-enhanced endpoint, interaction beta -0.335385, q 0.0209695) -> emphysema (female-enhanced endpoint, interaction beta 0.001917, q 0.9912109).
- URXCOP: emphysema (male-enhanced endpoint, interaction beta -0.136270, q 0.5719261) -> T2D (female-enhanced endpoint, interaction beta 0.134692, q 0.0223701).
- URXUCD: thyroid_disease (male-enhanced endpoint, interaction beta -0.007598, q 0.9680362) -> arthritis (female-enhanced endpoint, interaction beta 0.140688, q 0.0147906).

## Interpretation boundary

The scan identifies directionally compatible patterns, not confirmed reciprocal disease splits. In particular, the absence of any pair with both endpoints at q < 0.05 means that no exposure currently meets a two-endpoint fixed-FDR criterion. Any follow-up robustness or focused validation should be pre-specified separately and must not be inferred from this scan alone.
