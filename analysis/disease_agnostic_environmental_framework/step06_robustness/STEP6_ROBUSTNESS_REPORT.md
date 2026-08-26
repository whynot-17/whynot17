# Step 6 robustness audit report

Generated (UTC): 2026-08-26T10:32:50.673829+00:00

## Frozen scope

This audit was run after locking `ROBUSTNESS_RUBRIC_LOCK.md`. It applies the same modules, thresholds, estimator, and tier logic to exactly the two Step 5 FDR-supported signals: URXCOP (MCOP) and LBXPFHS (PFHS).

- Step 5 test family: **29 frozen tests**; denominator unchanged and not recomputed here.
- Step 6 audited signals: **2**.
- Primary rerun reproduction: **2/2** within absolute log-OR difference <=1e-8.

## Primary estimates and scorecard

| Biomarker | OR | 95% CI | P | BH-FDR (29) | Fingerprint | Tier |
|---|---:|---|---:|---:|---|---|
| LBXPFHS | 0.6244 | 0.4707–0.8283 | 0.00146 | 0.04234 | F2 | L2 | C2 | H1 | D2 | T2 | A0 | E1 | Tier B |
| URXCOP | 1.246 | 1.078–1.44 | 0.003311 | 0.04801 | F2 | L2 | C2 | H0 | D2 | T2 | A1 | E2 | Robust Tier A |

Robust Tier A signals under the frozen rubric: **URXCOP**.
Persistent technical-warning signals (A0): **LBXPFHS**.

## Audit interpretation

- `H` is a reported exposure-by-cycle heterogeneity tag, not a hard deletion gate.
- PFHS and MCOP are evaluated in the same audit. A `converged_with_warning` status is retained and interpreted as a technical warning, never silently upgraded to clean convergence.
- The pairwise co-exposure model is secondary and uses the target biomarker's own survey weight; it does not replace the primary model or change the Step 5 FDR result.
- Sex-specific estimates and the formal exposure-by-sex interaction are reported descriptively/secondarily; they are not used to redefine the primary signal.
- MCOP is interpreted as the urinary biomarker for the DINP-related exposure axis. This audit does not establish causality or prove that MCOP itself is a direct CRC mechanism.

## Frozen component definitions

The complete definitions and thresholds are in `ROBUSTNESS_RUBRIC_LOCK.md`; in particular, the attenuation threshold is 0.25 absolute log(OR), and Tier A is F2 + L>=1 + C>=1 + D>=1 + T>=1 + A>=1.
