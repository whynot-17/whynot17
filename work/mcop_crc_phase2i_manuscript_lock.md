# Phase 2I — MCOP/DINP–CRC manuscript lock

Date: 2026-08-23
Status: FROZEN FOR MANUSCRIPT BUILD
Branch: `phase2f-compartment-external-replication`

## 1. Paper-level scientific claim

This paper is **not** a network-toxicology paper and is **not** a causal-mechanism paper.

The central claim is:

> A data-first toxicogenomic screen nominated a DINP/MiNP exposure axis, and urinary MCOP — a widely used human biomarker of DINP exposure — was positively associated with prevalent colorectal cancer in seven NHANES cycles using complex-survey methods, with strong robustness to leave-one-cycle-out, alternative creatinine handling, tail exclusion, diagnosis-timing exclusions, and co-exposure adjustment. CRC transcriptomic analyses independently identify an epithelial-centered PPAR/nuclear-receptor program as a candidate mechanistic bridge, but do not establish DINP→PPAR→CRC causality.

Do not strengthen this claim later without new prospective or experimental evidence.

## 2. Frozen primary human analysis

### Population
- NHANES 2005–2006 through 2017–2018 only: 7 two-year cycles.
- Adults age >=20.
- Primary comparator: CRC cases versus cancer-free controls (`MCQ220=2`).
- Complete-case primary N = 9,936.
- CRC cases = 70.

### Exposure
- Primary: `log2(URXCOP)`.
- Interpretation: OR per doubling of urinary MCOP.
- MCOP is used as a human biomarker of DINP exposure; do not call MCOP identical to MiNP or a molecularly equivalent CTD exposure.

### Complex survey design
- Cycle-specific phthalate subsample weight divided by 7.
- Cycle-unique strata and PSU identifiers.
- Primary complete-case singleton strata = 0.
- Standard implementation: R `survey::svyglm`.
- Python Taylor-sandwich implementation is an independent numerical replication, not the sole manuscript estimator.

### Frozen covariates
- age
- sex
- race/ethnicity
- BMI
- smoking
- PIR
- log2 urinary creatinine

Do not add/remove covariates based on P values.

### Frozen primary result
R `survey::svyglm`:
- OR per MCOP doubling = 1.2455068
- 95% CI = 1.0773085–1.4399655
- standard P = 0.0033958
- design-df P = 0.0033113

Python check:
- OR = 1.2455068
- 95% CI = 1.0775254–1.4396756
- P = 0.0033114

Relative Python-vs-R logOR change = 4.877e-11%.

## 3. Frozen robustness hierarchy

Primary robustness analyses to retain:

1. Leave-one-cycle-out (all seven pooled re-estimates remain OR>1 and 95% CI excludes 1).
2. Age >=40.
3. Exclude CRC diagnosis <1 year.
4. Exclude CRC diagnosis <2 years.
5. Exclude CRC diagnosis <5 years.
6. Exclude top 1% MCOP.
7. Exclude top 2.5% MCOP.
8. Creatinine-normalized exposure sensitivity.
9. Pairwise adjustment for MEHHP, MEOHP, MECPP, MBzP and predefined non-MCOP phthalate burden.
10. Sex-specific estimates plus formal MCOP×sex interaction.

Do not promote male-only significance as effect modification because the formal interaction is not significant.

## 4. Known weaknesses that MUST remain visible

### Cross-sectional/prevalent CRC
NHANES exposure is measured at examination, after CRC diagnosis for prevalent cases. Timing exclusions reduce simple recent-diagnosis explanations but do not establish temporality or eliminate reverse causation.

### Cycle heterogeneity
The previously detected MCOP×cycle interaction is significant and must not be hidden. The paper should report that pooled association is robust to LOCO while cycle-specific slopes are heterogeneous, with 2011–2012 being the major discordant cycle.

### Small CRC case count
Only ~70 complete-case CRC cases. Avoid excessive subgrouping and avoid overinterpreting quartile categories.

### Quartiles
Survey-weighted Q4 vs Q1 is not statistically significant and P-trend is ~0.051. Quartiles are secondary descriptive analyses only.

### RCS
Survey-weighted RCS shows a significant overall exposure association but no evidence of nonlinearity. This supports the prespecified continuous log2 exposure as the primary parameterization; it should not be marketed as a dose-response proof of causality.

## 5. Molecular discovery language lock

### What Phase 1 actually nominated
- CTD/GeneCards screen nominated **MiNP/DINP exposure axis**.
- MiNP primary enrichment was strong but much of the CRC overlap was driven by a multi-phthalate co-treatment record.
- Pure/single-chemical evidence is substantially weaker.

### MCOP
- MCOP is the viable NHANES DINP biomarker.
- MCOP itself was not a strong CTD CRC molecular hit.

Allowed wording:
> The molecular screen nominated a DINP/MiNP exposure axis, motivating evaluation of MCOP as a human biomarker of DINP exposure.

Forbidden wording:
> CTD identified MCOP as a CRC-causing chemical.

## 6. Mechanism language lock

### Primary candidate module
PPAR/nuclear receptor score:
- PPARA
- PPARD
- PPARG
- NR1I2
- NR1I3
- NR1H2
- NR1H3

### Disease-state evidence
- TCGA paired bulk: PPAR/NR score lower in tumor.
- Census paired epithelial: strong PPAR/NR suppression in tumor-derived epithelium.
- GSE144735 paired epithelial: same direction but underpowered (6 pairs; not significant).
- Non-epithelial compartments are not uniformly suppressed; myeloid shows the opposite direction.

Interpretation:
> CRC-associated PPAR/NR remodeling appears epithelial-centered rather than a universal whole-tissue program.

### Secondary inflammatory module
RELA/STAT3 may be shown separately.

### Do not use
- The 9-gene DINP composite as a primary mechanism score.
- Statements that DINP causes PPAR suppression in human CRC epithelium.
- Statements that PPAR is a novel DINP mechanism.

Allowed wording:
> PPAR/nuclear-receptor remodeling is a candidate mechanistic bridge between DINP-related toxicology and CRC epithelial state.

## 7. Main figure architecture

### Figure 1 — Data-first discovery and study design
Goal: show the paper was not generated by starting from a fashionable chemical.

Panels:
A. CTD human chemical–gene universe and GeneCards CRC gene set workflow.
B. Full ranked screen highlighting MiNP at rank 24 while retaining DINP parent and MBzP as contextual anchors.
C. High-level multistage prioritization and translation: rank-24 MiNP signal -> DINP-related exposure axis -> MCOP as the human DINP biomarker. Display MiNP 27.4% versus MCOP 98.4% above LOD across seven cycles.

Detailed Top-30 gate-by-gate triage is reserved for a dedicated candidate-prioritization figure. Do not imply that MiNP was the top-ranked molecular hit, that molecular rank alone determined advancement, or that CTD directly nominated MCOP.

Do not show PPI or molecular docking.

### Figure 2 — Primary NHANES MCOP–CRC association
Panels:
A. Main R `survey::svyglm` result with N=9,936 / 70 CRC cases.
B. Python vs R numerical replication.
C. Leave-one-cycle-out forest plot (7 estimates + pooled estimate).
D. Per-cycle estimates may be shown in a separate panel to transparently display heterogeneity.

Primary visual message:
> The pooled association is not driven by a single NHANES cycle.

### Figure 3 — Robustness and exposure-shape analysis
Panels:
A. Forest plot: age>=40, diagnosis timing exclusions, top-tail exclusions, creatinine-normalized exposure.
B. Co-exposure-adjusted MCOP estimates.
C. Survey-weighted RCS curve.
D. Quartiles as a small secondary panel/inset only; explicitly label P-trend and wide CIs.

Do not headline quartile results.

### Figure 4 — CRC epithelial PPAR/NR disease-state validation
Panels:
A. TCGA matched tumor-normal paired score.
B. Census matched donor epithelial score.
C. GSE144735 paired epithelial direction replication, clearly showing n=6.
D. Compartment comparison (epithelial, endothelial, fibroblast, myeloid) to show epithelial specificity / opposing myeloid behavior.

Do not present this figure as DINP perturbation evidence.

### Figure 5 — Integrated evidence model
A restrained schematic, not a causal DAG claiming proof.

Show three evidence layers:
1. Molecular discovery: MiNP/DINP axis nominated.
2. Human association: urinary MCOP positively associated with prevalent CRC.
3. CRC disease state: epithelial PPAR/NR remodeling.

Use dashed or explicitly labeled hypothetical connector for:
DINP exposure -> epithelial PPAR/NR remodeling.

The figure must visually distinguish:
- observed association
- disease-state validation
- unproven mechanistic bridge

## 8. Main tables

### Table 1
NHANES analytic population characteristics, weighted where appropriate. Include CRC vs cancer-free controls.

### Table 2
Primary and prespecified sensitivity MCOP estimates.

### Supplementary tables
- chemical-screen top candidates and classification
- exact CTD/GeneCards audit
- all LOCO estimates
- per-cycle estimates and interaction test
- sex analyses
- diagnosis-timing exclusions
- tail exclusions
- creatinine sensitivity
- co-exposure models
- weighted/unweighted quartile cutpoints
- weighted/unweighted RCS knots
- TCGA/scRNA score definitions and paired results

## 9. Manuscript Results order

### Results 1
Data-first toxicogenomic screen nominates a DINP/MiNP exposure axis.

### Results 2
Urinary MCOP is positively associated with prevalent CRC in NHANES.

### Results 3
The MCOP association is robust to complex-survey reimplementation and leave-one-cycle-out analyses.

### Results 4
The association persists across predefined sensitivity analyses but exhibits cycle heterogeneity.

### Results 5
CRC transcriptomic data identify epithelial-centered PPAR/nuclear-receptor remodeling as a candidate mechanistic bridge.

Do not put mechanism before the human association.

## 10. Discussion structure

1. Principal finding: urinary MCOP–CRC association in a nationally representative complex-survey framework.
2. Why this matters in the context of limited human DINP cancer evidence.
3. Relation to prior phthalate–CRC epidemiology without claiming the broader field is empty.
4. Direct DINP colon toxicology as biological plausibility.
5. PPAR/NR as candidate bridge, with explicit recognition that DINP-to-PPAR direction is not directly validated in colon exposure data.
6. Strengths: data-first discovery, standard survey replication, LOCO robustness, exposure-shape audit, co-exposure sensitivity, multi-layer CRC transcriptomic localization.
7. Limitations: prevalent/cross-sectional outcome, reverse causation, 70 cases, cycle heterogeneity, biomarker specificity, MiNP-to-MCOP bridge, lack of prospective replication and direct colon DINP perturbation transcriptomics.
8. Future work: prediagnostic prospective cohort and/or intestinal epithelial/organoid DINP perturbation.

## 11. Working title options

Preferred neutral title:
> Urinary MCOP, a biomarker of DINP exposure, is associated with prevalent colorectal cancer in NHANES: a data-first toxicogenomic and transcriptomic study

Alternative shorter title:
> A data-first analysis links the DINP exposure axis to colorectal cancer through urinary MCOP and epithelial nuclear-receptor remodeling

Do not use "causes", "drives", "mediates", or "mechanism of" in the title.

## 12. Immediate execution tasks for Codex

1. Create `outputs/manuscript/`.
2. Generate publication-ready data tables for Figures 1–4 from frozen outputs; no re-estimation except formatting/plot preparation.
3. Generate Figure 2 and Figure 3 first because the human association is the paper's strongest evidence.
4. Then build Figure 4 with paired-patient/donor points visible, not only boxplots.
5. Build Figure 1 only from audited Phase 1 outputs and avoid implying MCOP was directly nominated by CTD.
6. Build Figure 5 last as a restrained evidence schematic.
7. Create `outputs/manuscript/mcop_crc_manuscript_skeleton.md` with Introduction/Methods/Results/Discussion headings and numeric results inserted from frozen files.
8. Create `outputs/manuscript/manuscript_source_map.csv` mapping every main-text number/claim to its exact source output file.
9. Do not run new exploratory analyses during Phase 2I unless a reproducibility or figure-integrity problem is discovered.

## 13. Stop rule

Phase 2I is manuscript construction, not hypothesis generation.

If a new analysis is proposed, it must answer one of only three questions:
- Does it verify a manuscript number?
- Does it address a predictable reviewer concern?
- Does it correct a reproducibility problem?

If none applies, do not run it.
