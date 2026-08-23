# Figure 1–5 Master Blueprint v1

Status: FROZEN FOR V1 FIGURE PRODUCTION
Target manuscript style: Environment International submission, with visual hierarchy inspired by high-impact Nature/Cell/Science figure design principles.

## Global figure principles

1. Figure 1 explains the entire study design and logic in one glance.
2. Figures 2–5 follow **one figure = one message**.
3. Visual hierarchy must be stronger than decorative styling: clear focal point, restraint, whitespace, limited colors, consistent typography, and real data points whenever possible.
4. Avoid generic bioinformatics aesthetics: rainbow heatmaps, large PPI webs, crowded GO bar plots, decorative molecular docking graphics, or excessive gradients.
5. Maintain fixed visual identities across all figures:
   - DINP / MCOP / exposure axis: one consistent accent family.
   - CRC / tumor: one consistent darker identity.
   - Normal / control: neutral gray family.
   - Hypothesized mechanistic bridge: lighter treatment and dashed connectors.
6. Main statistical figures must be generated directly from frozen result files; no manual recreation of numeric estimates.
7. Do not visually imply causal evidence where only association or disease-state convergence has been demonstrated.

---

# Figure 1 — Data-first prioritization of the DINP exposure axis and overall study design

## One-line purpose
Allow an editor or reviewer to understand the entire paper within ~5 seconds: unbiased environmental discovery -> DINP/MiNP prioritization -> MCOP human biomonitoring -> CRC association -> candidate biological convergence.

## Core message
**A DINP-related exposure axis was identified by a data-first environmental screen, translated to MCOP as a human biomarker, and then evaluated in NHANES and CRC transcriptomic datasets.**

## Layout
Horizontal four-part / funnel-like multi-panel composition: A -> B -> C -> D.

### Panel A — Environmental chemical discovery universe

Purpose: show that DINP was not preselected.

Content:
- 267 core environmental chemicals.
- CTD human chemical–gene interactions.
- GeneCards CRC-associated genes.
- Fisher/hypergeometric enrichment + BH-FDR.
- Degree-matched permutation.

Plot type:
- restrained discovery funnel / structured workflow diagram.

On-figure labels:
- `267 environmental chemicals`
- `CTD human chemical–gene interactions`
- `GeneCards CRC genes`
- `Enrichment + degree-matched permutation`

Do not include detailed ORs/P values here.

### Panel B — Ranked screen and MiNP/DINP prioritization

Purpose: show how the DINP/MiNP axis emerged from the wider chemical screen.

Plot type:
- horizontal ranked lollipop plot preferred over bubble plot.

Suggested x-axis:
- `-log10(BH-FDR)` or another frozen significance/enrichment metric chosen consistently from Phase 1 output.

Suggested encoding:
- point size = CRC overlap count, if visually clean.
- most chemicals in neutral treatment.
- MiNP/DINP strongly but tastefully highlighted.
- MBzP may be lightly marked as a later human-validation counterexample.

Do not imply MiNP was the absolute top chemical; communicate prioritization within the audited shortlist.

### Panel C — Translation bridge: MiNP/DINP -> MCOP

Purpose: explain why the molecular nomination and population biomarker are not the same molecule.

Core logic:
`MiNP/DINP molecular nomination -> MCOP as a human biomarker of DINP exposure`

Short supporting labels:
- `MiNP molecular evidence nominated a DINP-related exposure axis`
- `MCOP selected for NHANES because of broad availability and high detectability`

Plot type:
- clean exposure-biomarker translation diagram.

Critical boundary:
- do not state or visually imply `CTD discovered MCOP`.

### Panel D — Overall study roadmap

Three stages:
1. Discovery: CTD/GeneCards -> MiNP/DINP axis.
2. Human biomonitoring: NHANES 2005–2018, N=9,936, CRC=70, urinary MCOP.
3. Biological interpretation: TCGA paired bulk, CELLxGENE donor-level single-cell analyses, GSE144735 external epithelial dataset.

Use wording such as `transcriptomic validation` and `candidate mechanistic bridge`, not `mechanism proven`.

## Figure 1 legend core sentence
A data-first environmental screen prioritized a DINP-related exposure axis from 267 chemicals using CTD–GeneCards CRC enrichment. MCOP was then selected as a human biomarker of DINP exposure for NHANES evaluation, followed by transcriptomic and single-cell analyses to identify candidate CRC-associated biological convergence.

---

# Figure 2 — Urinary MCOP is positively associated with prevalent colorectal cancer across seven NHANES cycles

## One-line purpose
Show the strongest human result immediately: positive MCOP–CRC association, independently reproduced by standard survey implementation, and robust to leaving out any single NHANES cycle.

## Core message
**Urinary MCOP is positively associated with prevalent CRC, and the pooled result is neither an implementation artifact nor driven by one NHANES cycle.**

## Layout
Four panels with strong hierarchy: A largest, C second largest, B compact validation panel, D transparency panel.

### Panel A — Hero primary result

Primary model:
- N=9,936.
- CRC=70.
- exposure = log2 urinary MCOP.
- interpretation = OR per doubling.
- R `survey::svyglm`: OR=1.2455068, 95% CI=1.0773085–1.4399655, standard P=0.0033958, design-df P=0.0033113.

Plot type:
- hero single-estimate forest / compact estimate card.

Visual focus:
- `OR 1.25 per doubling of MCOP` should be immediately visible.
- show 95% CI and P with restraint.

### Panel B — Independent implementation check

Compare:
- R `survey::svyglm`.
- independent Python Taylor-sandwich implementation.

Python result:
- OR=1.2455068.
- 95% CI=1.0775254–1.4396756.
- P=0.0033114.

Plot type:
- mini forest / paired estimate comparison.

Message:
- virtually identical point estimates.

### Panel C — Leave-one-cycle-out robustness

Rows:
- pooled overall.
- drop 2005–06.
- drop 2007–08.
- drop 2009–10.
- drop 2011–12.
- drop 2013–14.
- drop 2015–16.
- drop 2017–18.

Frozen LOCO estimates:
- drop 2005–06: OR 1.196759, 95% CI 1.027085–1.394464, P 0.021804.
- drop 2007–08: OR 1.264449, 95% CI 1.086107–1.472077, P 0.002852.
- drop 2009–10: OR 1.217756, 95% CI 1.052295–1.409233, P 0.008734.
- drop 2011–12: OR 1.334879, 95% CI 1.121840–1.588374, P 0.001379.
- drop 2013–14: OR 1.260707, 95% CI 1.076562–1.476350, P 0.004470.
- drop 2015–16: OR 1.233493, 95% CI 1.052823–1.445168, P 0.009953.
- drop 2017–18: OR 1.216737, 95% CI 1.038627–1.425391, P 0.015681.

Plot type:
- horizontal forest plot.

Message:
- all pooled re-estimates remain >1 and 95% CIs exclude 1.

Do not label these as seven independent replications.

### Panel D — Per-cycle estimates and heterogeneity

Plot type:
- compact forest / strip plot.

Purpose:
- transparently display cycle-specific heterogeneity, including the discordant 2011–12 estimate.

Statistical note to show in legend or small annotation:
- MCOP×cycle heterogeneity test P approximately 0.006.

This panel should not visually dominate Panel A or C.

## Figure 2 legend core sentence
In complex-survey analyses of NHANES 2005–2018, each doubling of urinary MCOP was associated with 24.6% higher odds of prevalent CRC. The estimate was numerically replicated using an independent R survey implementation and remained positive in all leave-one-cycle-out analyses, although cycle-specific heterogeneity was present.

---

# Figure 3 — The MCOP–CRC association is robust across predefined sensitivity analyses and co-exposure adjustment

## One-line purpose
Address the predictable reviewer objections in one figure: recent diagnosis, extreme exposure values, urinary dilution handling, co-exposure, and functional form.

## Core message
**The positive MCOP–CRC association persists across major prespecified sensitivity analyses; categorical dose-response is less monotonic, whereas continuous/RCS modeling supports an overall positive association without clear nonlinearity.**

## Layout
A large grouped forest panel + a smaller co-exposure forest + a clean RCS panel. Quartiles are optional small inset or Supplementary.

### Panel A — Main robustness forest grouped by reviewer concern

Group 1: Population / diagnosis-timing sensitivity
- age >=40: OR 1.221445, 95% CI 1.048222–1.423293, P 0.01084.
- exclude diagnosis <1 year: OR 1.24123, 95% CI 1.0624–1.4501, P 0.00691.
- exclude diagnosis <2 years: OR 1.26294, 95% CI 1.0713–1.4889, P 0.00586.
- exclude diagnosis <5 years: OR 1.26451, 95% CI 1.0586–1.5105, P 0.01013.

Group 2: Exposure-distribution / dilution sensitivity
- exclude top 1%: OR 1.270787, 95% CI 1.08307–1.49104, P 0.003647.
- exclude top 2.5%: OR 1.193576, 95% CI 1.01823–1.39912, P 0.029386.
- creatinine-normalized exposure: OR 1.243998, 95% CI 1.074737–1.439917, P 0.003790.

Plot type:
- grouped horizontal forest plot.

### Panel B — Co-exposure adjustment

MCOP after pairwise adjustment for:
- MEHHP: OR 1.222318, 95% CI 1.05741–1.41294, P 0.007071.
- MEOHP: OR 1.218929, 95% CI 1.05457–1.40891, P 0.007839.
- MECPP: use frozen current output value.
- MBzP: use frozen current output value.
- predefined non-MCOP phthalate burden: use frozen current output value.

Codex must pull exact current values from the frozen co-exposure output before plotting; do not hard-code from memory if source file differs.

Plot type:
- compact forest plot.

Message:
- MCOP estimate remains positive after other phthalate adjustment; do not claim mixture effects are fully eliminated.

### Panel C — Survey-weighted restricted cubic spline

Use frozen weighted knots and final survey-weighted model output.

Display:
- adjusted OR curve.
- 95% CI ribbon.
- reference OR=1 line.
- rug or exposure-density marks only if visually clean.

Annotate:
- overall F-test P ≈ 0.00233.
- nonlinear F-test P ≈ 0.365.

Message:
- significant overall association; no evidence of nonlinearity.
- do not state that linearity is proven.

### Panel D — Quartile analysis (optional main inset; otherwise Supplementary)

Survey-weighted category estimates:
- Q2 vs Q1 OR 0.8575, 95% CI 0.3769–1.9509.
- Q3 vs Q1 OR 1.6076, 95% CI 0.6711–3.8511.
- Q4 vs Q1 OR 1.3061, 95% CI 0.6464–2.6389.
- P-trend 0.051321.

Message:
- categorical estimates are imprecise and not monotonic.
- primary inference remains the continuous per-doubling model.

## Figure 3 legend core sentence
The positive association between urinary MCOP and prevalent CRC remained directionally and statistically robust after excluding recent CRC diagnoses, trimming the upper exposure tail, alternative creatinine handling, and pairwise co-exposure adjustment. Weighted spline analysis supported a positive overall association without clear evidence of nonlinearity; categorical quartile analyses were less monotonic and remained secondary.

---

# Figure 4 — CRC-associated PPAR/nuclear-receptor remodeling is epithelial-centered and accompanied by inflammatory activation

## One-line purpose
Show where the candidate biology is localized rather than merely reporting another bulk pathway result.

## Core message
**CRC-associated PPAR/NR suppression localizes primarily to tumor-derived epithelium, while myeloid cells show the opposite PPAR/NR direction; epithelial RELA/STAT3 activation occurs in parallel.**

## Layout
A: TCGA matched tissue; B: Census paired epithelial; C: compartment localization; D: external GSE144735 directional replication.

### Panel A — TCGA matched tumor-normal PPAR/NR score

Data:
- 32 patient-matched tumor-normal pairs.
- PPAR/NR median tumor-minus-normal delta = -0.533.
- Wilcoxon P = 7.404e-5.

Plot type:
- paired dot/line plot with all patient pairs visible.

Do not replace with a boxplot alone.

### Panel B — Census paired epithelial donors

Primary paired dataset:
- 36 paired donors.

Scores:
- PPAR/NR median delta = -0.419, P = 4.29e-7.
- RELA/STAT3 median delta = +1.167, P = 1.08e-7.
- 9-gene composite median delta ≈ +0.011, P = 0.636 (not a primary panel unless needed as a small negative-control inset).

Plot type:
- two side-by-side paired donor plots, PPAR/NR and RELA/STAT3.

Visual emphasis:
- PPAR/NR down and RELA/STAT3 up as distinct modules.

### Panel C — Compartment-specific paired deltas

Compartments and frozen paired results:
- epithelial: PPAR/NR delta -0.419, P 4.29e-7.
- endothelial: PPAR/NR delta -0.019, P 0.609.
- fibroblast: PPAR/NR delta +0.172, P 0.504.
- myeloid: PPAR/NR delta +0.610, P 7.97e-9.

Optional secondary inflammatory values can be included only if clarity remains high.

Plot type:
- paired-delta summary plot / dot-box hybrid.

Message:
- PPAR/NR behavior is compartment-specific, not a uniform whole-tissue program.
- this helps explain TCGA-vs-GTEx reference dependence.

### Panel D — Independent GSE144735 epithelial direction check

Data:
- 6 matched patients.
- PPAR/NR median tumor-minus-normal delta = -0.312.
- Wilcoxon P = 0.688.

Plot type:
- paired patient plot with all six pairs visible.

Annotation:
- `n=6; directionally concordant, underpowered`.

Do not call this statistically significant independent validation.

## Figure 4 legend core sentence
Bulk and single-cell analyses identified lower PPAR/nuclear-receptor signaling in CRC tumor epithelium relative to matched normal epithelium, accompanied by increased RELA/STAT3 activity in the paired Census epithelial analysis. Myeloid cells showed the opposite PPAR/NR direction, supporting compartment-specific remodeling rather than a uniform whole-tissue program. GSE144735 was directionally concordant but underpowered.

---

# Figure 5 — Integrated evidence model linking a DINP-related exposure axis, urinary MCOP, and epithelial PPAR/nuclear-receptor remodeling in CRC

## One-line purpose
Close the paper by distinguishing observed evidence from the unproven mechanistic bridge.

## Core message
**Three evidence layers converge — data-driven DINP-axis nomination, human MCOP–CRC association, and CRC epithelial PPAR/NR remodeling — but DINP exposure -> epithelial PPAR/NR remodeling remains a hypothesis requiring direct validation.**

## Layout
Three clean evidence blocks from left to right.

### Layer 1 — Data-driven discovery
- CTD/GeneCards.
- MiNP/DINP exposure axis nominated.

### Layer 2 — Human biomonitoring
- urinary MCOP.
- NHANES 2005–2018.
- OR 1.246 per doubling.
- robustness across LOCO and major sensitivity analyses.

### Layer 3 — CRC biological state
- epithelial PPAR/NR down.
- epithelial RELA/STAT3 up.
- compartment-specific remodeling.

## Connector rules

Use solid connectors only for directly observed evidence relationships.

Use a dashed connector for:
`DINP exposure -> epithelial PPAR/NR remodeling`

Label the dashed connector:
`candidate mechanistic bridge`

Do not use language such as `causes`, `drives`, `mediates`, or `mechanism established`.

## Figure 5 legend core sentence
The study provides three converging layers of evidence: data-driven prioritization of a DINP-related exposure axis, a robust human biomonitoring association between urinary MCOP and prevalent CRC, and CRC-associated epithelial PPAR/nuclear-receptor remodeling. The causal bridge between DINP exposure and epithelial remodeling remains hypothetical and requires prospective and experimental validation.

---

# Production order

## First production batch
1. Figure 2 — strongest human result.
2. Figure 3 — reviewer-facing robustness.
3. Figure 4 — biological localization.

## Second production batch
4. Figure 1 — manuscript-level study overview.
5. Figure 5 — integrated synthesis.

---

# Output requirements for V1 figures

For each figure generate:
- editable vector PDF.
- SVG.
- high-resolution PNG (>=600 dpi where appropriate).
- individual panel files in addition to assembled multi-panel figure.
- underlying plotting table used for each panel.
- exact script used to generate the figure.

All figures must:
- use the same typography.
- use consistent line widths and panel labels.
- use one restrained manuscript-wide palette.
- retain vector text in PDF/SVG.
- display real paired patient/donor points when such data exist.
- avoid excessive legends and redundant numeric labels.

# Visual review questions for V1

Before V2, evaluate:
1. Can the main message of each figure be identified within 5 seconds?
2. Which panel looks most like a generic bioinformatics figure?
3. Which elements should be removed rather than added?
4. Is the visual hierarchy comparable to a high-impact journal figure?
5. Do Figures 1–5 clearly belong to one manuscript and one design system?

# Final scientific boundary

The figures may support:
- data-driven nomination of a DINP-related exposure axis;
- positive association of urinary MCOP with prevalent CRC;
- robustness of the human association;
- CRC-associated, epithelial-centered PPAR/NR remodeling.

The figures must not claim:
- DINP causes CRC;
- MCOP was directly nominated by the CTD screen;
- DINP causes PPAR/NR suppression in human CRC epithelium;
- GSE144735 provides statistically significant replication.
