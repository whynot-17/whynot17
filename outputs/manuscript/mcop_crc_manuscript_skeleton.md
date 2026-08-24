# Manuscript skeleton — Environmental chemical prioritization, urinary MCOP, and CRC

Status: MANUSCRIPT BUILD v1
Analysis status: FROZEN
Freeze authority: `work/environmental_crc_final_analysis_freeze_20260824.md`

## Working title

**Urinary MCOP, a biomarker of DINP exposure, is associated with prevalent colorectal cancer in NHANES: a data-first environmental screening and transcriptomic convergence study**

## Central claim

An outcome-blinded data-first screen of environmental chemicals identified a DINP-related axis that was human-testable through urinary MCOP. Across 15 unique NHANES biomarker tests, the DINP-related MCOP test was FDR-supported and the sole Robust Tier A signal after a uniform robustness audit. Independent CRC transcriptomic analyses showed epithelial-state-specific PPAR/nuclear-receptor remodeling, while the exposure-to-epithelial-state bridge remained untested.

---

# Abstract

## Background
Environmental chemical exposures may contribute to colorectal carcinogenesis, but candidate selection is often hypothesis-led and vulnerable to post hoc prioritization. We developed an outcome-blinded, data-first framework linking toxicogenomic nomination, human biomarker testability, population-level screening, and CRC disease-state transcriptomics.

## Methods
We first evaluated 267 environmental chemicals using CTD human chemical–gene associations and GeneCards CRC genes, followed by an outcome-blinded actionability framework incorporating entity validity, interpretable exposure mapping, human biomarker availability, biomarker detectability, NHANES cycle coverage, and survey-model testability. Eligible chemical–biomarker mappings were reduced to unique NHANES biomarker tests before CRC outcome screening. Survey-weighted logistic regression was used for the human screen. The leading signal was subjected to a uniform robustness audit, and independent CRC transcriptomic datasets were used to evaluate epithelial PPAR/nuclear-receptor disease-state remodeling.

## Results
Of 267 starting chemicals, 87 chemical–biomarker mappings met the predefined human-testability criteria and corresponded to 15 unique NHANES biomarker tests. Two tests passed BH-FDR <0.05: the DINP-related urinary MCOP test and serum PFHS. After a prespecified robustness audit, MCOP was the sole Robust Tier A signal. In the frozen primary NHANES model (N=9,936; 70 CRC cases), each doubling of urinary MCOP was associated with higher odds of prevalent CRC (OR=1.2455, 95% CI 1.0773–1.4400; design-df P=0.00331). The pooled association remained positive in leave-one-cycle-out analyses, although between-cycle heterogeneity was significant. In paired CRC epithelial transcriptomic analyses, the prespecified PPAR/nuclear-receptor core was lower in tumor epithelium (median tumor-normal delta=-0.419; FDR≈9.3×10^-7), with independent KEGG, Reactome, metabolic, and regulon definitions showing predominantly concordant suppression. Within-state analyses indicated suppression in enterocyte-like epithelium but modest elevation in secretory-like epithelium, supporting state-specific remodeling rather than universal epithelial suppression.

## Conclusions
A data-first, outcome-blinded framework identified urinary MCOP as the strongest human-testable DINP-related CRC signal in this analysis. The epidemiologic association is cross-sectional and does not establish temporality or causality. CRC transcriptomics independently support epithelial-state-specific PPAR/nuclear-receptor remodeling as disease-state convergence, but direct DINP/MCOP-to-epithelial-state causality remains untested.

---

# Introduction

## Paragraph 1 — Clinical and environmental context
- CRC remains a major global cancer burden.
- Established hereditary and lifestyle factors do not account for all inter-individual risk variation.
- Environmental chemicals are biologically plausible contributors, but human evidence is uneven and often fragmented across toxicology, biomonitoring, and epidemiology.

## Paragraph 2 — Methodological gap
- Conventional environmental-cancer studies frequently start with a preselected chemical or class.
- Candidate-first designs can create literature-driven selection and obscure how many alternatives were considered.
- A stronger design is to nominate candidates before outcome inspection and then require independent human biomarker feasibility.

## Paragraph 3 — DINP-related context
- DINP is a high-volume phthalate exposure of current toxicologic interest.
- Human carcinogenic evidence remains limited.
- MiNP and MCOP are distinct DINP-related metabolites/biomarkers and must not be treated as chemically interchangeable.

## Paragraph 4 — Study objective
We therefore constructed an outcome-blinded multistage framework to: (1) screen environmental chemicals using toxicogenomic CRC relevance; (2) evaluate human biomonitoring actionability without using candidate-specific CRC outcomes; (3) systematically screen the resulting unique NHANES biomarker tests for association with prevalent CRC; (4) apply the same robustness audit to all leading tests; and (5) assess whether independent CRC transcriptomic data show disease-state programs compatible with the nominated toxicologic axis.

---

# Methods

## Study design overview
Describe the five-stage design corresponding to Figures 1–5. Explicitly state that candidate-specific CRC outcome estimates were excluded from actionability prioritization.

## Phase 1: Toxicogenomic discovery
- CTD Homo sapiens chemical–gene interaction universe.
- GeneCards CRC gene set.
- Fisher/hypergeometric enrichment, BH-FDR, degree-matched permutation sensitivity.
- Database associations are treated as associations, not causal targets.

## Outcome-blinded actionability prioritization
Starting universe: 267 chemicals/candidate entities.

Frozen gates:
- E: valid/interpretable entity identity.
- X: interpretable human exposure mapping.
- B: human biomarker mapping available.
- D: biomarker detectability: D0 <50%, D1 50–<90%, D2 ≥90% above LOD.
- C: NHANES coverage: C0 ≤2 cycles, C1 3–4, C2 ≥5.
- T: survey-model testability: T2 if complete N≥5000 and CRC cases≥60 with required infrastructure; T1 if complete N≥500 and CRC cases≥20; T0 otherwise.

Eligibility:
`permissive = E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1`

`moderate = E=1 & X=1 & B=1 & D>=1 & C=2 & T>=1`

`strict = E=1 & X=1 & B=1 & D=2 & C=2 & T=2`

State explicitly that multiple chemicals can map to the same NHANES biomarker test and that shared mapping does not imply chemical equivalence or parent-specific attribution.

## NHANES systematic human screen
- Statistical unit: 15 unique NHANES biomarker tests.
- Primary model: `CRC ~ log2(exposure) + age + sex + race/ethnicity + BMI + smoking + PIR`.
- Urinary biomarkers additionally adjust for log2 urinary creatinine.
- BH-FDR controlled across all 15 tests.
- Complex survey design retained.

## Primary MCOP analysis
- NHANES 2005–2006 through 2017–2018.
- Adults age ≥20 years.
- CRC cases versus cancer-free controls.
- Complete-case N=9,936; CRC cases=70.
- Primary exposure: log2(URXCOP), interpreted as OR per doubling of urinary MCOP.
- Cycle-specific phthalate subsample weight divided by seven; cycle-unique strata and PSU identifiers.
- Primary estimator: R `survey::svyglm`.

## Robustness analyses
Prespecified analyses retained in the manuscript:
- leave-one-cycle-out;
- age ≥40 years;
- diagnosis-timing exclusions (<1, <2, <5 years);
- top 1% and 2.5% exposure exclusions;
- creatinine-normalized exposure;
- pairwise co-exposure adjustment;
- sex-specific estimates with formal interaction;
- exposure-shape assessment.

## CRC transcriptomic disease-state analyses
Describe TCGA, paired single-cell epithelial analysis, GSE144735 direction replication, PPAR/NR definition audit, within-state analyses, and donor-level program correlations. Make clear these datasets contain CRC disease-state evidence, not DINP perturbation evidence.

## Statistical principles
- Donor/patient-level paired inference where applicable.
- No single-cell pseudo-replication for main inference.
- Multiple-testing correction reported where prespecified.
- No post-freeze exploratory analysis unless formally unfreezed for QC/reviewer concern.

---

# Results

## Result 1 — Outcome-blinded prioritization reduces 267 chemicals to 15 unique human biomarker tests
The starting universe contained 267 environmental chemicals/candidate entities. Sequential actionability filtering retained 259 entity-valid candidates, 135 with interpretable exposure mappings, 134 with a human biomarker mapping, 127 with adequate detectability, 124 with sufficient cycle coverage, and 87 with predefined human survey-model testability. These 87 eligible chemical–biomarker mappings corresponded to 15 unique NHANES biomarker tests. Twenty-seven mappings met the strict D2/C2/T2 definition.

MiNP remained a molecular DINP-related nominee but failed the direct human detectability gate in the updated audit (~40.7% above LOD, D0). MCOP (URXCOP), a distinct urinary DINP-related biomarker, had high detectability and adequate cycle coverage and therefore entered the systematic human screen. This transition reflects biomarker actionability rather than chemical equivalence.

## Result 2 — Systematic screening identifies two FDR-supported biomarker tests
Across the 15 unique NHANES biomarker tests, two passed BH-FDR <0.05. The DINP-related urinary MCOP test was positively associated with prevalent CRC (screen OR≈1.246; BH-FDR≈0.02484), whereas serum PFHS showed an inverse association (OR≈0.624; BH-FDR≈0.02190). The remaining 13 tests did not pass the 15-test FDR threshold.

A uniform post-screen robustness framework classified MCOP as the sole Robust Tier A test (`F2 | L2 | C2 | H0 | D2 | T2 | A1 | E2`). PFHS remained Tier B because of persistent technical/event-count limitations.

## Result 3 — Urinary MCOP is positively associated with prevalent CRC in the frozen primary survey model
The frozen complete-case primary analysis included 9,936 participants, including 70 prevalent CRC cases. In R `survey::svyglm`, each doubling of urinary MCOP was associated with higher odds of prevalent CRC (OR=1.2455068, 95% CI 1.0773085–1.4399655; design-df P=0.0033113). Independent Python Taylor-sandwich implementation reproduced the point estimate to numerical precision.

This result is an association with prevalent CRC and does not establish whether higher MCOP preceded cancer development.

## Result 4 — The MCOP association is robust overall but heterogeneous across NHANES cycles
All seven leave-one-cycle-out pooled re-estimates remained above OR=1 with 95% CIs excluding 1, indicating that no single cycle alone accounted for the pooled association. Prespecified age, diagnosis-timing, upper-tail, creatinine-handling, and co-exposure analyses retained the overall positive direction. However, the MCOP×cycle interaction was statistically significant, with 2011–2012 the major discordant cycle. Therefore, the pooled signal should be described as robust to cycle omission despite significant between-cycle heterogeneity, not as homogeneous replication across seven independent cohorts.

## Result 5 — CRC epithelium shows state-specific PPAR/nuclear-receptor remodeling
In paired CRC epithelial analyses, the frozen seven-gene PPAR/nuclear-receptor core was lower in tumor-derived epithelium (median tumor-normal delta=-0.4186; FDR≈9.3×10^-7). The direction was independently reproduced by receptor-only and NR-partner modules, KEGG PPAR signaling, Reactome PPAR-alpha lipid regulation, peroxisomal lipid metabolism, mitochondrial fatty-acid beta oxidation, Hallmark fatty-acid metabolism, enterocyte metabolic differentiation, and DoRothEA PPARA/PPARG activity.

The signal was not uniform across epithelial states. Within-state analysis showed significant PPAR/NR suppression in enterocyte-like epithelium (n=24 paired donors; median delta≈-0.191; FDR≈5.26×10^-4) but modest elevation in secretory-like epithelium (n=27; median delta≈+0.068; FDR≈0.018). These findings support mixed epithelial compositional redistribution plus within-state regulatory remodeling rather than universal PPAR suppression.

Donor-level analyses did not support a simple PPAR-low → RELA/STAT3-high mediation model. PPAR/NR and inflammatory/stress programs should therefore be interpreted as parallel CRC disease-state changes. Toxicology evidence for DINP/MiNP provides general nuclear-receptor/PPAR plausibility, but the direct DINP/MCOP → CRC epithelial-state link remains untested.

---

# Discussion

## Principal finding
This study used a data-first, outcome-blinded selection framework to move from a broad environmental chemical universe to a limited set of human-testable biomarker models. The DINP-related urinary MCOP test emerged as the sole Robust Tier A signal after systematic human screening and uniform robustness assessment.

## Why the human result matters
The MCOP association was identified only after the same 15-test screen and therefore should not be presented as a candidate chosen because it produced a favorable outcome. The result was robust to multiple predefined sensitivities but was cross-sectional and heterogeneous across cycles.

## Interpretation of the DINP/MiNP/MCOP relationship
MiNP, MCOP, and DINP are related but non-equivalent chemical entities. The molecular screen nominated a DINP-related molecular/exposure axis; MCOP provided the viable human biomarker test. MCOP itself should not be described as the molecular perturbagen nominated by CTD.

## Disease-state convergence rather than causal mechanism
Independent CRC transcriptomics showed that PPAR/nuclear-receptor remodeling is strongly detectable in the epithelial disease state and survives multiple pathway definitions. However, this evidence does not demonstrate that DINP or MCOP caused the observed epithelial state. The state-specific results further argue against a simplistic pan-epithelial shutdown model.

## Strengths
- broad data-first discovery rather than single-chemical nomination;
- explicit outcome firewall before human screening;
- unique biomarker-test level FDR control;
- complex-survey NHANES implementation;
- uniform robustness scoring of competing signals;
- transparent cycle heterogeneity;
- paired donor/patient transcriptomic inference;
- explicit separation of exposure evidence, disease-state evidence, and untested mechanism.

## Limitations
- prevalent/cross-sectional CRC outcome;
- temporality and reverse causation cannot be excluded;
- survivor bias is possible;
- spot urine exposure measurement error;
- only approximately 70 complete-case CRC cases;
- statistically significant cycle heterogeneity;
- residual confounding;
- biomarker-to-parent specificity varies across environmental chemicals;
- MiNP/MCOP/DINP are related but non-equivalent;
- no prospective external MCOP–CRC replication;
- no direct colon/intestinal epithelial DINP perturbation transcriptome closing the exposure-to-state bridge.

## Future work
Priority future studies are prospective or prediagnostic validation of DINP-related biomarkers and controlled intestinal epithelial/organoid perturbation experiments that directly test whether DINP-related exposure can reproduce the state-specific PPAR/NR remodeling observed in CRC.

---

# Conclusion

An outcome-blinded data-first framework identified urinary MCOP as the strongest robust DINP-related human biomarker signal associated with prevalent CRC in NHANES. Independent CRC transcriptomics demonstrated epithelial-state-specific PPAR/nuclear-receptor remodeling, providing disease-state convergence but not exposure-to-mechanism causality. These findings motivate prospective epidemiologic replication and direct intestinal epithelial exposure experiments.

---

# Figure architecture lock

- **Figure 1:** Data-first molecular discovery.
- **Figure 2:** Outcome-blinded actionability prioritization: 267 → 87 eligible chemical–biomarker mappings → 15 unique NHANES biomarker tests.
- **Figure 3:** Systematic 15-test human screen and primary MCOP association.
- **Figure 4:** MCOP robustness, cycle heterogeneity, and exposure-shape analyses.
- **Figure 5:** CRC epithelial state-specific PPAR/NR disease-state convergence with dashed, explicitly untested DINP-related environmental bridge.

# Writing prohibitions

Do not use causal terms such as `causes`, `drives`, `mediates`, or `mechanism of` for the DINP/MCOP–CRC relationship. Do not describe LOCO cycles as independent replication cohorts. Do not describe the Figure 5 disease-state program as proof of DINP perturbation. Do not merge MiNP, MCOP, and DINP identities.
