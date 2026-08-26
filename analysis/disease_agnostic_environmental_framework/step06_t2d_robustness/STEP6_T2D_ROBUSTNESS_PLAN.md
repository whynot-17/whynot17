# Step 6 T2D robustness and exposure-cluster audit lock

This lock is applied to the frozen first-round T2D screen after the outcome
firewall. It does not alter the 29-test primary screen or its BH-FDR family.

## Scope

- Include exactly the tests marked `FDR_supported=True` in the frozen
  `step05_t2d_screen/t2d_primary_29_tests.csv` output.
- Apply the same modules and thresholds to every included test.
- Do not rerun, narrow, or replace the 29-test BH-FDR correction.
- Do not query GeneCards, disease-specific CTD, transcriptomic resources, or
  literature in this stage.

## Uniform audit modules

1. Reproduce the frozen primary model from the same assay-specific exposure
   files, cycles, weights, outcome frame, and covariates.
2. Leave one included NHANES cycle out at a time.
3. Fit each included cycle separately and test exposure-by-cycle interaction.
4. Exclude the upper 1% and 2.5% of the exposure distribution.
5. For urinary tests, use log2(exposure) - log2(urinary creatinine) as a
   normalization sensitivity; serum/blood tests are not applicable.
6. Run age >=40 restriction.
7. Run sex-stratified models and a formal exposure-by-sex interaction.
8. Run LOD/detectable-only sensitivity when the frozen registry indicates that
   fewer than 90% of measured values are above LOD.
9. Build a descriptive exposure correlation and cluster audit across all
   FDR-positive tests, using pairwise complete log2 values and a cycle-adjusted
   Spearman correlation sensitivity. This audit is not an outcome-based
   selection step.

## Frozen interpretation rules

- LOCO direction: L2 if all estimable LOCO estimates preserve the pooled
  direction and all 95% CIs exclude 1; L1 if all preserve direction but some
  CIs cross 1; L0 if any estimate reverses direction.
- Cycle concordance: C2 if >=80% of estimable cycle-specific estimates preserve
  the pooled direction; C1 if 60-79%; C0 if <60%.
- Heterogeneity tag: H2 for interaction P>=0.10; H1 for 0.05<=P<0.10; H0
  for P<0.05. H is reported and is not a deletion gate.
- Tail stability: T2 if both tail exclusions preserve direction and the
  maximum absolute log-OR change is <=0.25; T1 if direction is preserved but
  the change is larger; T0 for direction instability or no estimable result.
- Technical status: A2 if all applicable fits converge; A1 if there are only
  localized convergence warnings; A0 if any applicable fit fails or the
  warning is persistent.
- A candidate is labelled `robust_fdr_candidate` only when it is FDR-positive
  and has L>=1, C>=1, T>=1, and A>=1. This label prioritizes follow-up; it does
  not change the primary q-value or prove causality.

## Exposure clustering rule

- Pairwise correlations are calculated without using T2D status, P values, or
  effect estimates.
- A high-correlation edge is defined as absolute cycle-adjusted Spearman
  |rho| >=0.70 with pairwise N >=500.
- Exposure clusters are connected components of that prespecified graph. Tests
  without a high-correlation edge remain singleton clusters.
- Cluster membership is reported as a dependence/interpretation diagnostic,
  not as a reason to delete a primary test or reassign FDR.

## Output firewall

The audit may prioritize candidates for later work, but it must not perform
mechanistic analysis or reinterpret a T2D association as causal.
