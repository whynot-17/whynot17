# Environmental CRC 15-axis systematic robustness audit

Generated: 2026-08-25T15:31:46.850592+00:00

## Frozen scope

All 15 axes were read automatically from the frozen primary screen. The same primary model and the same robustness modules were applied to every axis. The primary BH-FDR remains the original 15-axis correction; no robustness subset was used to recompute significance.

Primary reproduction QC: 15/15 axes reproduced their stored OR to absolute log-OR difference <=1e-8.

## Robustness scorecard

| Axis | Biomarker | OR | BH-FDR (15-axis) | Fingerprint | Tier |
|---|---|---:|---:|---|---|
| DEHP-related exposure axis;phthalate exposure axis | URXECP | 1.211 | 0.2469 | F0 | L1 | C1 | H0 | D2 | T2 | A1 | E2 | Exploratory |
| PFAS exposure axis | LBXPFDE | 0.7161 | 0.5429 | F0 | L1 | C0 | H1 | D2 | T2 | A1 | E1 | Exploratory |
| PFAS exposure axis | LBXPFNA | 0.8059 | 0.5429 | F0 | L1 | C1 | H0 | D2 | T2 | A0 | E1 | Exploratory |
| phthalate exposure axis | URXMEP | 0.9329 | 0.5429 | F0 | L1 | C0 | H0 | D2 | T2 | A1 | E2 | Exploratory |
| DEHP-related exposure axis | URXMHP | 1.093 | 0.5429 | F0 | L1 | C1 | H0 | D2 | T2 | A1 | E2 | Exploratory |
| bisphenol exposure axis | URXBPH | 1.061 | 0.6228 | F0 | L1 | C1 | H2 | D2 | T2 | A0 | E1 | Exploratory |
| phthalate exposure axis | URXMBP | 0.8954 | 0.6228 | F0 | L0 | C1 | H0 | D2 | T2 | A1 | E2 | Exploratory |
| phthalate exposure axis | URXMIB | 0.9659 | 0.7528 | F0 | L0 | C1 | H0 | D2 | T2 | A1 | E2 | Exploratory |
| phthalate exposure axis | URXMZP | 0.9545 | 0.7528 | F0 | L1 | C1 | H0 | D2 | T2 | A1 | E2 | Exploratory |
| PAH exposure axis | URXP04 | 0.9542 | 0.7528 | F0 | L0 | C1 | H2 | D2 | T2 | A1 | E2 | Exploratory |
| PAH exposure axis | URXP10 | 1.072 | 0.7528 | F0 | L0 | C0 | H2 | D2 | T2 | A1 | E2 | Exploratory |
| PAH exposure axis | URXP02 | 1.003 | 0.9815 | F0 | L0 | C0 | H2 | D0 | T2 | A2 | E2 | Exploratory |
| DINP-related exposure axis | URXCOP | 1.246 | 0.02484 | F2 | L2 | C2 | H0 | D2 | T2 | A1 | E2 | Robust Tier A |
| PFAS exposure axis | LBXPFHS | 0.6244 | 0.0219 | F2 | L2 | C2 | H1 | D2 | T2 | A0 | E1 | Tier B |
| phthalate exposure axis | URXMOH | 1.204 | 0.1307 | F1 | L1 | C1 | H0 | D2 | T2 | A1 | E2 | Tier B |

Robust Tier A axes: **1**.

## Required interpretation questions

- FDR-supported primary axes: **2**; the original 15-axis correction is unchanged.
- Axes with all available LOCO estimates directionally concordant: **10**; axes with direction instability: **5**.
- Axes with significant exposure×cycle heterogeneity (P<0.05): **9**.
- Axes retaining direction after all diagnosis-timing exclusions: **14**.
- Axes with any fit warning or non-estimable audit component: **14**; axes with persistent warning/technical concern (A0): **3**.

### Warning interpretation

Convergence warnings were retained rather than hidden. For warning-only fits, the Newton-IRLS routine returned finite estimable coefficients and sandwich variance but did not reach the configured tolerance within the iteration limit; the warning is therefore algorithmic and does not by itself imply coefficient failure. The frozen primary estimates were independently reproduced (15/15; absolute log-OR difference <=1e-8).
- Persistent warning axes: **LBXPFHS, LBXPFNA, URXBPH**.
- Localized warning axes: **LBXPFDE, URXCOP, URXECP, URXMBP, URXMEP, URXMHP, URXMIB, URXMOH, URXMZP, URXP04, URXP10**.
- A persistent warning is treated as an A0 technical concern in the scorecard; it is not converted into a positive finding or silently removed. Warning repetition across secondary fits is explicitly recorded in `environmental_crc_15axis_model_warnings.csv`.

MCOP/URXCOP is evaluated in exactly the same scorecard as PFHS/LBXPFHS and the other 13 axes. MiNP remains a distinct molecular DINP nominee and is not silently converted into MCOP; it is outside the 15-axis human screen because its direct detectability failed the frozen actionability gate.

## Prespecified tag definitions

- F: F2=15-axis BH-FDR<0.05; F1=nominal P<0.05; F0 otherwise.
- L: L2=same direction and all LOCO CIs exclude 1; L1=same direction with some CIs crossing 1; L0=direction instability.
- C: C2>=80% same-direction cycle estimates; C1=60–79%; C0<60%.
- H: H2=Pinteraction>=0.10; H1=0.05–<0.10; H0<0.05.
- D/T: D2/T2 preserve direction with maximum absolute log-OR change <=0.25; D1/T1 preserve direction with greater attenuation; D0/T0 direction instability.
- A: A2 no warning; A1 localized warning with primary/sensitivity fits otherwise estimable; A0 persistent warning (>=75% of applicable fits, or >=50% when the primary itself warns) or any fit failure.
- E: E2>=60 CRC cases; E1=30–59; E0<30.
- Robust Tier A was frozen as F2, L>=1, C>=1, D>=1, T>=1, A>=1; H is reported as a penalty/evidence tag, not a hard deletion gate.

This is a robustness/association audit of cross-sectional NHANES data and does not establish causality or eliminate reverse causation/survivor bias.
