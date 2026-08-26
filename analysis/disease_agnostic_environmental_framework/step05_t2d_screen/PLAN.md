# T2D Disease Plug-in Plan

## Purpose

Use type 2 diabetes (T2D) as the first new disease plug-in for the already frozen disease-agnostic environmental screening framework.

The purpose is not to alter the front-end environmental test set or to tune the analysis toward any particular chemical. The purpose is to test whether the outcome-firewalled framework can produce a statistically and biologically coherent downstream disease demonstration once a disease outcome is introduced after the environmental test set has been frozen.

## Non-negotiable firewall

The following components are inherited unchanged from the existing disease-agnostic framework:

1. CTD-derived environmental universe and deterministic exclusions.
2. Human/NHANES biomarker mapping.
3. Human actionability filtering.
4. Collapse at unique NHANES biomarker-test level.
5. The frozen set of exactly 29 NHANES biomarker tests.
6. Test-specific assay files, analytic cycle coverage, and survey/subsample weights.
7. Exposure transformations used in the rebuilt assay-specific Step 5 architecture.
8. Primary multiplicity family size: 29 tests.
9. Benjamini-Hochberg correction across all 29 primary T2D association tests.

No disease information, T2D association result, GeneCards disease gene set, disease-specific CTD interaction evidence, transcriptomic result, or literature finding may be used to add, remove, or re-rank the 29 tests before the primary T2D screen is complete.

## Disease plug-in

### Primary outcome

Construct a T2D outcome from NHANES diabetes data without using any exposure laboratory file to define the participant frame.

The implementation must distinguish clearly between:

- diagnosed diabetes,
- probable diabetes identified from objective glycemic measurements when available,
- non-diabetic controls,
- indeterminate participants,
- likely type 1 diabetes / early-onset insulin-dependent cases when identifiable.

The exact operational definition must be documented before inspecting any exposure-T2D association result.

The primary T2D definition should prioritize a clinically defensible adult T2D phenotype rather than a maximally broad diabetes label. Objective laboratory information may be used to improve disease classification, but it must not depend on exposure values.

### Primary comparison

T2D cases versus participants classified as non-diabetic controls.

Participants with ambiguous diabetes status should not be silently assigned to controls.

## Primary model

Retain the rebuilt Step 5 assay-specific survey architecture.

For each of the 29 frozen biomarker tests:

T2D ~ log2(exposure) + age + sex + race/ethnicity + BMI + smoking + PIR

For urinary biomarkers additionally adjust for log2(urinary creatinine), using the same creatinine handling principle as the rebuilt CRC screen.

Use the test-specific NHANES laboratory/subsample survey weight divided by the number of included cycles for that assay, consistent with the existing rebuilt screen.

The primary model must be identical across all 29 tests except for the prespecified urine creatinine term.

## Primary multiplicity analysis

- Number of primary hypotheses: 29.
- Procedure: Benjamini-Hochberg FDR.
- Primary conventional threshold: q < 0.05.
- Also report exact q-values for every biomarker.
- Do not change the denominator according to observed results.
- Do not use parent-chemical collapse as the primary analysis.

A secondary exploratory label of q < 0.10 may be reported descriptively, but it must not replace the q < 0.05 primary multiplicity assessment.

## Required primary outputs

Produce at minimum:

1. `t2d_primary_29_tests.csv`
   - test_id
   - biomarker
   - exposure axis
   - matrix
   - analytic N
   - T2D cases
   - controls
   - OR per exposure doubling
   - 95% CI
   - P value
   - BH-FDR q value
   - convergence/status
   - cycle coverage
   - survey weight used

2. `t2d_primary_ranked.csv`
   - all 29 tests sorted by P and q value.

3. `t2d_outcome_qc.csv`
   - cycle-specific participant counts
   - T2D cases
   - controls
   - excluded/indeterminate diabetes classifications
   - key diabetes variables available in each cycle.

4. `t2d_merge_audit.csv`
   - test-specific exposure rows
   - outcome-frame rows
   - merged rows
   - complete-case N
   - complete-case T2D cases.

5. `t2d_analysis_manifest.json`
   - frozen test-set hash
   - script hashes
   - input data paths/hashes where feasible
   - model specification
   - FDR denominator
   - timestamps
   - outcome-definition text.

6. `t2d_screen_summary.md`
   - number estimable
   - number nominal P<0.05
   - number q<0.10
   - number q<0.05
   - strongest signals
   - warnings
   - no biological interpretation at this stage.

## Mandatory QC before interpretation

The screen is not considered valid until all of the following are checked:

1. Each test uses its own assay-specific laboratory file rather than a shared convenience frame.
2. T2D outcome/covariate frames are built independently of exposure laboratory files.
3. Cycle coverage matches each frozen test registry entry.
4. Test-specific survey/subsample weights are used correctly.
5. 29 frozen tests enter the FDR family even if a test is non-estimable; any handling of non-estimable tests must be explicitly documented and must not silently shrink the prespecified family.
6. T2D case counts are plausible across cycles.
7. Urinary biomarkers receive the prespecified creatinine adjustment.
8. No downstream biological database is queried before the primary epidemiologic results and candidate set are frozen.

## Candidate progression rule

After the primary T2D screen is complete, freeze a candidate set before robustness analysis.

Evidence tiers should be reported rather than collapsed into a single binary label:

- Tier A screening evidence: BH-FDR q < 0.05.
- Tier B suggestive screening evidence: 0.05 <= q < 0.10.
- Nominal candidate: P < 0.05 but q >= 0.10.

All candidates admitted to the robustness stage must be evaluated using the same robustness rubric. No candidate-specific sensitivity analysis may be invented to rescue one biomarker.

## Robustness stage (only after primary screen freeze)

For every carried-forward candidate, apply the same stress tests where technically applicable:

- leave-one-cycle-out analysis,
- cycle-specific estimates,
- exposure-by-cycle heterogeneity,
- upper-tail exclusion / influential-value sensitivity,
- alternative creatinine handling for urinary biomarkers,
- age restriction / reverse-causation sensitivity,
- sex-stratified estimates and interaction,
- co-exposure adjustment within chemically coherent groups where prespecified and technically valid,
- convergence and numerical reproducibility audit.

Robustness does not replace FDR. It addresses stability and prioritization among screening candidates.

## Biological convergence stage (only after robustness candidate set is frozen)

Only after the epidemiologic candidate set is fixed:

1. map biomarker to the relevant chemical/parent exposure conservatively;
2. obtain candidate-chemical CTD-associated genes;
3. obtain T2D GeneCards disease genes using a documented query/version/threshold;
4. calculate overlap/enrichment;
5. test pathway/network convergence;
6. seek independent transcriptomic or other molecular support.

Database associations must not be described as causal evidence.

## Interpretation rules

Allowed:

- `The outcome-firewalled framework identified/prioritized ...`
- `X met the 5% FDR threshold in the T2D demonstration.`
- `X met a suggestive 10% FDR threshold but not the conventional 5% threshold.`
- `X was the strongest nominal screening signal and showed subsequent robustness.`

Not allowed:

- claiming the disease was prospectively preregistered if it was not;
- claiming the entire framework is outcome-blinded after the disease plug-in step;
- changing thresholds after viewing results and presenting them as prespecified;
- dropping unhelpful diseases/tests/results without recording the development history;
- using GeneCards/CTD/transcriptomics to retroactively define the primary epidemiologic test family.

## Stop/go criteria for this proof-of-concept

After primary T2D screening:

### Strong GO
At least one q < 0.05 signal with adequate analytic N/events and no major QC failure.

### Moderate GO
No q < 0.05, but a clearly separated q < 0.10 candidate with strong effect precision and adequate event counts. Continue only as an explicitly exploratory demonstration.

### Weak / STOP for T2D showcase
No q < 0.10 and no convincing candidate separation, or major outcome/assay/weight/QC problems. Do not alter the 29-test front end to improve the result.

## Immediate execution task

Implement a new assay-specific T2D Step 5 by reusing the validated architecture of:

`analysis/disease_agnostic_environmental_framework/step05_crc_screen_rebuilt/run_step05_rebuilt.py`

Do not overwrite CRC outputs.

The first implementation milestone is strictly:

1. define and audit the NHANES T2D outcome across all required cycles;
2. run the exact frozen 29 tests;
3. calculate BH-FDR with denominator 29;
4. generate the required QC and result tables;
5. stop before robustness, GeneCards, CTD disease-specific biology, transcriptomics, literature interpretation, or manuscript rewriting.
