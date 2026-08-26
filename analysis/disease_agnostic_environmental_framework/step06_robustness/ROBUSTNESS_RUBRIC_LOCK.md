# Step 6 robustness rubric lock

This file freezes the robustness rubric before inspecting the Step 6 results.
It is the two-signal adaptation of the historical environmental CRC
robustness audit. It does not alter the Step 5 CRC screen, its 29-test
multiple-testing family, or the list of FDR-supported signals.

## Frozen scope

- The audit includes exactly the two signals supported by the frozen Step 5
  outcome-aware screen: `URXCOP` (MCOP) and `LBXPFHS` (PFHS).
- Both signals receive the identical audit modules, thresholds, score fields,
  and tier logic. No module is added, removed, or redefined after looking at
  either signal's audit result.
- The Step 5 BH-FDR denominator remains **29**, including the one frozen test
  whose CRC model was not estimable. Step 6 never recomputes FDR on a
  robustness subset.
- The primary estimator remains the Step 5 Python survey-weighted
  Taylor/PSU-sandwich logistic estimator, with analyte-specific laboratory
  weight, cycle-specific PSU and strata, and the frozen covariate model:
  `CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR`; urine
  tests additionally include `log2(creatinine)`.

## Required audit modules

Each signal is evaluated with the same modules:

1. Primary screen reproduction and fixed Step 5 FDR status.
2. Event/sample-size context and design diagnostics.
3. Leave-one-cycle-out (LOCO) fits.
4. Cycle-specific fits and global exposure-by-cycle interaction test.
5. Diagnosis-timing exclusions for CRC diagnoses within 1, 2, and 5 years.
6. Exclusion of the upper 1% and 2.5% of the exposure distribution.
7. LOD/detectability sensitivity using the frozen registry flags.
8. Urinary creatinine-normalized sensitivity where applicable.
9. Age >=40 restriction.
10. Sex-stratified fits and a formal exposure-by-sex interaction test.
11. Pairwise co-exposure adjustment between the two supported signals, using
    the target signal's survey weight and complete cases for both exposures.

Cycle-specific estimates are descriptive precision-limited components; the
global interaction test is the prespecified heterogeneity statistic. A
co-exposure model does not replace the primary model and does not alter the
Step 5 screen.

## Frozen component scores

- **F (false-discovery support):** F2 = Step 5 BH-FDR < 0.05; F1 = nominal
  P < 0.05; F0 otherwise.
- **L (LOCO direction):** L2 = all estimable LOCO estimates preserve the
  pooled direction and all LOCO 95% CIs exclude 1; L1 = all preserve the
  pooled direction but at least one CI crosses 1; L0 = direction instability.
  If fewer than three cycles are available, L is not applicable.
- **C (cycle concordance):** C2 = >=80% of estimable cycle-specific estimates
  preserve the pooled direction; C1 = 60–79%; C0 = <60%.
- **H (heterogeneity tag):** H2 = exposure-by-cycle interaction P >= 0.10;
  H1 = 0.05 <= P < 0.10; H0 = P < 0.05. H is reported, not used as a hard
  deletion gate.
- **D (diagnosis timing):** D2 = all applicable timing exclusions preserve
  direction and the maximum absolute log-OR change is <= 0.25; D1 = direction
  preserved with a larger change; D0 = direction instability or no estimable
  timing result.
- **T (upper tail):** T2/T1/T0 use the identical rule as D for top 1% and
  top 2.5% exclusion.
- **A (algorithmic/technical status):** A2 = no applicable fit warning;
  A1 = localized `converged_with_warning` with finite estimable primary and
  sensitivity results; A0 = any fit failure/non-estimability in an applicable
  audit component or a persistent warning (warning in >=75% of applicable
  fits, or >=50% when the primary fit itself warns).
- **E (event support):** E2 = >=60 CRC cases; E1 = 30–59; E0 = <30.

The attenuation threshold for D/T is a maximum absolute change in log(OR) of
**0.25**. An extreme survey-weight ratio above **100** is recorded as a
diagnostic warning; it does not silently delete a signal.

## Frozen tier logic

- **Robust Tier A:** F2, L >= 1, C >= 1, D >= 1, T >= 1, and A >= 1.
- **Tier B:** not Tier A, but F >= 1, E >= 1, and at least one of L/C/D/T
  is >= 1.
- **Exploratory:** all remaining combinations.

The Tier A definition intentionally does not use H as a hard gate, matching
the historical rubric. H0 remains a prominently reported heterogeneity
concern. Event support E is always reported and contributes to Tier B; it is
not silently substituted for any robustness component.

## Interpretation firewall

This is a robustness/association audit of cross-sectional NHANES data. It does
not establish temporality or causality and does not eliminate reverse
causation or survivor bias. MCOP is the urine biomarker used for the DINP-axis
human analysis; Step 6 must not relabel it as direct molecular proof that MCOP
itself caused the CRC-associated state.
