# Environmental chemicals → CRC — Final analysis freeze

Date: 2026-08-24
Status: **ANALYSIS FROZEN — MANUSCRIPT/FIGURE CONSTRUCTION ONLY**
Branch: `phase2f-compartment-external-replication`
Freeze reference HEAD at lock initiation: `b0aa8c5774401d2238b70e3ec232f05ab6be9ed8`

This document supersedes earlier exploratory figure/manuscript architecture where it conflicts with the current 267-chemical actionability framework or the Figure 5 PPAR contradiction audit.

## 1. Paper-level claim lock

This paper is a **data-first environmental epidemiology + toxicogenomic prioritization + human disease-state convergence study**.

It is **not** a causal-mechanism paper and **not** a network-pharmacology paper.

Frozen central claim:

> An outcome-blinded data-first screen of environmental chemicals identified a DINP-related axis that was human-testable through urinary MCOP. Across 15 unique NHANES biomarker tests, the DINP-related MCOP test was FDR-supported and the sole Robust Tier A signal after a uniform robustness audit. Independent CRC transcriptomic analyses showed epithelial-state-specific PPAR/nuclear-receptor remodeling, but the exposure-to-epithelial-state bridge remains untested and must not be presented as causal.

Forbidden upgrades without new prospective or experimental evidence:
- DINP causes CRC.
- MCOP causes CRC.
- DINP causes PPAR/NR suppression in CRC epithelium.
- PPAR/NR mediates the MCOP–CRC association.
- RELA/STAT3 is downstream of PPAR loss in the current human data.

## 2. Discovery/actionability freeze

### 2.1 Starting universe
- 267 environmental chemicals/candidates.
- These are chemicals/candidate entities, not 267 biomarkers.

### 2.2 Outcome firewall
Candidate prioritization is outcome-blinded.

Candidate-specific CRC OR, P, CI, FDR, LOCO, and cycle-specific CRC effects must not enter pre-outcome actionability eligibility.

### 2.3 Frozen actionability gates

- `E`: valid/interpretable entity identity.
- `X`: interpretable human exposure axis mapping.
- `B`: human biomarker mapping available.
- `D`: detectability of the selected primary biomarker.
  - D0: <50% above LOD
  - D1: 50% to <90%
  - D2: >=90%
- `C`: NHANES cycle coverage.
  - C0: <=2 cycles
  - C1: 3–4 cycles
  - C2: >=5 cycles
- `T`: predefined survey-model testability.
  - T2: complete N>=5000 and CRC cases>=60 with required survey/covariate infrastructure
  - T1: complete N>=500 and CRC cases>=20 with required infrastructure
  - T0: otherwise
- `M`: molecular evidence level; does not determine human-testability eligibility.
- `N`: novelty/manual-review annotation; does not determine main eligibility.

Frozen eligibility definitions:

`permissive = E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1`

`moderate = E=1 & X=1 & B=1 & D>=1 & C=2 & T>=1`

`strict = E=1 & X=1 & B=1 & D=2 & C=2 & T=2`

### 2.4 Frozen attrition

- 267 total core chemicals
- 259 E-valid
- 135 E+X interpretable
- 134 E+X+B biomarker available
- 127 +D detectable
- 124 +C adequate coverage
- 87 +T human-testable chemical–biomarker mappings
- 87 moderate-eligible mappings
- 27 strict-eligible mappings

Important wording lock:

> **87 eligible chemical–biomarker mappings corresponded to 15 unique NHANES biomarker tests.**

Do not write:
- “87 chemicals biologically collapsed into 15 equivalent exposures.”
- “15 exposure axes” without qualification when the mapping is family/proxy based.

The 15 tests are the unique biomarker/model units actually screened in NHANES. Shared biomarker mapping does not imply chemical equivalence or parent-specific attribution.

## 3. Chemical identity lock

### MiNP
- Molecular DINP-related nominee.
- Must remain chemically distinct from MCOP.
- Latest v2 direct biomarker detectability: ~40.7% above LOD, D0.
- Therefore fails the direct detectability gate for primary human testing.

### MCOP
- `URXCOP` human urinary biomarker/test used for the DINP-related human exposure axis.
- Must not be described as the same molecule as MiNP.
- Human-testable and robustly associated with prevalent CRC.
- MCOP is an exposure biomarker, not a demonstrated molecular perturbagen of CRC epithelium.

### DINP parent
- Parent/exposure-axis concept.
- Molecular evidence is not identical to MiNP or MCOP evidence.

Frozen manuscript relationship wording:

> The data-first screen nominated a DINP-related molecular/exposure axis; urinary MCOP provided the viable human biomarker test for this axis.

## 4. Human screen freeze

The human statistical screening unit is the **unique NHANES biomarker test**, not every mapped parent chemical.

Frozen primary model:

`CRC ~ log2(exposure) + age + sex + race/ethnicity + BMI + smoking + PIR`

For urinary biomarkers add `log2(creatinine)`.

BH-FDR is controlled across all 15 unique biomarker tests.

Frozen key screen results:

- DINP-related / `URXCOP` (MCOP): OR ~1.246, BH-FDR ~0.02484.
- PFAS / `LBXPFHS`: OR ~0.624, BH-FDR ~0.02190.
- Other 13 tests do not pass 15-test BH-FDR<0.05.

Do not select MCOP merely because it was positive; it advances by the prespecified uniform robustness audit.

## 5. Robustness freeze

Frozen robustness framework uses F/L/C/H/D/T/A/E evidence tags.

### MCOP / DINP-related test
- Fingerprint: `F2 | L2 | C2 | H0 | D2 | T2 | A1 | E2`
- Final status: **Robust Tier A**
- Sole Robust Tier A test.

### PFHS
- FDR-supported but persistent technical-warning/event-count limitations.
- Final status: Tier B.

### Important heterogeneity lock
MCOP has statistically significant between-cycle heterogeneity.

Allowed wording:

> MCOP showed the strongest overall robustness profile despite significant between-cycle heterogeneity.

Forbidden wording:
- homogeneous across cycles
- replicated in seven independent cohorts
- consistent in every cycle

LOCO analyses are robustness checks on the pooled NHANES analysis, not seven independent external replications.

## 6. Primary human estimate lock

Primary complete-case N: 9,936
CRC cases: 70
NHANES cycles: 2005–2006 through 2017–2018

Primary R `survey::svyglm` estimate:
- OR per doubling of MCOP = 1.2455068
- 95% CI = 1.0773085–1.4399655
- design-df P = 0.0033113

The association is cross-sectional/prevalent.

Mandatory limitations:
- temporality is not established
- reverse causation remains possible
- survivor bias remains possible
- spot urine introduces exposure measurement error
- only ~70 complete-case CRC cases
- cycle heterogeneity is significant
- residual confounding remains possible

Diagnosis-timing exclusions, tail exclusions, creatinine sensitivity and co-exposure models may reduce selected alternative explanations but do not establish causality.

## 7. Figure 5 mechanism freeze

Figure 5 status: **YELLOW — disease-state convergence, not causal mechanism**.

### 7.1 Definition audit result
The frozen PPAR/NR result survives independent definitions and is predominantly DOWN in paired CRC epithelium.

Key locked examples:
- Frozen 7-gene PPAR/NR core: median tumor-normal delta = -0.419, FDR ~9.3e-07.
- PPAR receptor-only module: down.
- NR partner module: down.
- KEGG PPAR signaling: down.
- Reactome PPAR-alpha lipid regulation: down.
- Peroxisomal lipid metabolism: down.
- Mitochondrial fatty-acid beta oxidation: down.
- Hallmark fatty-acid metabolism: down.
- Enterocyte metabolic differentiation: down.
- DoRothEA PPARA activity: down.
- DoRothEA PPARG activity: strongly down.

Therefore the current PPAR/NR disease-state result is not an artifact of one custom score.

### 7.2 State-specificity lock
The overall epithelial result reflects **mixed composition change + within-state remodeling**.

Within-state audit:
- enterocyte-like: PPAR/NR significantly down
- secretory-like: PPAR/NR modestly up
- other epithelial: underpowered/not interpretable

Allowed wording:

> CRC is associated with epithelial-state-specific PPAR/nuclear-receptor remodeling, including suppression in the enterocyte-like compartment alongside epithelial compositional redistribution.

Forbidden wording:
- universal PPAR suppression in all CRC epithelial states
- pan-epithelial PPAR shutdown

### 7.3 Inflammation/RELA/STAT3 lock
PPAR/NR loss and inflammatory programs are parallel disease-state changes in the current human data.

Do not draw or write:

`PPAR/NR down -> RELA/STAT3 up`

The donor-delta audit does not support a simple inverse mediation/coupling model. RELA/STAT3 regulon-delta associations with PPAR/NR delta are not significant.

Allowed representation:

`CRC -> PPAR/NR/differentiation remodeling`

and separately

`CRC -> inflammatory/stress remodeling`

with association/convergence indicated where justified, not directional mediation.

### 7.4 Exposure bridge lock
The DINP/MiNP toxicology literature supports general nuclear-receptor/PPAR plausibility but does not reproduce or establish the CRC epithelial state.

Non-colon liver/ovary/hepatocyte evidence must not be used as the primary mechanistic bridge.

The link

`DINP/MCOP exposure -> CRC epithelial PPAR/NR state`

remains **UNTESTED**.

If shown in Figure 5, use dashed/dotted/hypothetical representation and explicit wording such as “candidate bridge” or “causal link untested.”

## 8. Collision/novelty freeze

Targeted collision audit did not identify an exact prior complete chain combining:

DINP/MCOP human exposure + CRC human association + epithelial single-cell PPAR/NR state remodeling.

Allowed wording:

> Targeted search did not identify an exact prior study integrating this complete chain.

Forbidden wording:
- first ever
- world-first
- unprecedented

unless a future formal systematic novelty review justifies that language.

## 9. Final main-figure architecture lock

### Figure 1 — Data-first molecular discovery
Purpose: establish unbiased environmental-chemical discovery and distinguish database association from causality.

### Figure 2 — Outcome-blinded actionability prioritization
Frozen narrative:

267 chemicals -> 87 human-testable chemical–biomarker mappings -> 15 unique NHANES biomarker tests

Must visually distinguish direct/metabolite mappings from family/proxy mappings.

### Figure 3 — Systematic human screening and primary MCOP association
Show all 15 tests and the two FDR-supported tests, with MCOP highlighted only after the uniform screen.

### Figure 4 — Robustness / heterogeneity / exposure-shape
Show why MCOP becomes the sole Robust Tier A signal while retaining cycle heterogeneity and limitations.

### Figure 5 — CRC epithelial disease-state convergence
Primary message:

> Epithelial-state-specific PPAR/NR and metabolic-differentiation remodeling in CRC, with a clearly untested dashed DINP-related environmental bridge.

Figure 5 is not a causal DAG.

## 10. Analysis stop rule

From this freeze onward, **no new exploratory analysis is allowed** unless it satisfies at least one of the following:

1. verifies a number already intended for the manuscript;
2. resolves a reproducibility/QC discrepancy;
3. addresses a specific predictable reviewer concern that cannot be answered from frozen outputs;
4. corrects a factual or coding error discovered during figure/manuscript construction.

A new analysis must be explicitly labeled `UNFREEZE_REQUEST` before execution and document:
- exact question
- why frozen outputs are insufficient
- prespecified method
- expected manuscript decision affected

If none applies, do not run it.

## 11. What remains allowed

Allowed without unfreezing:
- figure rendering and layout refinement
- source-table formatting
- manuscript drafting
- source-map construction
- reference verification
- numeric cross-checking without re-estimation
- supplementary-table assembly from frozen outputs
- wording and causal-boundary audits

Not allowed without unfreezing:
- changing eligibility thresholds
- trying additional covariates because a result looks better
- choosing new biomarker mappings using CRC outcomes
- adding new exposure candidates after viewing CRC associations
- rerunning Figure 5 with alternative gene sets solely to optimize significance
- adding new public datasets solely to strengthen the desired story
- subgroup fishing

## 12. Immediate manuscript-build tasks

1. Update the old Phase 2I manuscript lock wherever it conflicts with this file.
2. Replace all “15 exposure axes” shorthand with “15 unique NHANES biomarker tests” or “biomarker-defined exposure tests” as appropriate.
3. Update figure numbering/architecture to the frozen five-figure structure above.
4. Update Figure 5 wording to state-specific PPAR/NR remodeling and remove any PPAR->RELA/STAT3 causal arrow.
5. Build `outputs/manuscript/manuscript_source_map.csv` mapping every main-text numerical claim to an exact frozen output.
6. Build the manuscript skeleton from the frozen results only.
7. Perform a final MiNP/MCOP/DINP identity audit across all main and supplementary materials.

## 13. Final status

**Core scientific analysis: FROZEN.**

The project now moves from discovery/analysis to:

`FIGURE BUILD -> MANUSCRIPT BUILD -> CONSISTENCY AUDIT -> SUBMISSION PACKAGE`

Any future analytical expansion requires an explicit unfreeze decision.
