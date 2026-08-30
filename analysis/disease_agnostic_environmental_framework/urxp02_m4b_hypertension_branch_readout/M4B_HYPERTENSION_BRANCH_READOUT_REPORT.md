# M4b — hypertension-branch tissue readout

## Scope

This is a readout of the frozen M4 hypertension-branch results only. The fixed family remains 440 genes × 9 prespecified GTEx tissues = 3,960 tests; no model was refit, no new multiple-testing correction was introduced, and no single-cell or mechanism analysis was performed.

## Main result

Across the seven hypertension-relevant tissues (three arteries, two heart tissues, kidney cortex, and adrenal gland), there are **65 fixed-FDR hits**: **24 female-biased** and **41 male-biased**. This is a mixed localized pattern, not a coherent female-biased hypertension-branch shift.

- Artery - Tibial: 27 FDR hits (10 female-biased, 17 male-biased).
- Artery - Aorta: 20 FDR hits (9 female-biased, 11 male-biased).
- Heart - Atrial Appendage: 6 FDR hits (2 female-biased, 4 male-biased).
- Heart - Left Ventricle: 6 FDR hits (1 female-biased, 5 male-biased).

Kidney cortex has **0 fixed-FDR hits**, despite 127 tests meeting the prespecified absolute-effect flag. The kidney GTEx sample size is much smaller than the large artery tissues, so this is a power/context limitation rather than evidence of a kidney null mechanism.

## Recurring genes

The repeated-gene pattern is descriptive only; recurrence was not assigned a new P value:

- PLA2G5: 4 tissues — Adrenal Gland; Artery - Aorta; Artery - Coronary; Artery - Tibial (male_only).
- CES1: 3 tissues — Adrenal Gland; Artery - Aorta; Artery - Tibial (female_only).
- BCL2: 3 tissues — Artery - Aorta; Artery - Tibial; Heart - Left Ventricle (mixed).
- SFN: 2 tissues — Artery - Aorta; Heart - Left Ventricle (female_only).
- SERPINE1: 2 tissues — Artery - Aorta; Artery - Tibial (male_only).
- UCHL1: 2 tissues — Heart - Atrial Appendage; Heart - Left Ventricle (male_only).
- CLCNKB: 2 tissues — Artery - Aorta; Artery - Tibial (female_only).
- DRD1: 2 tissues — Artery - Coronary; Artery - Tibial (male_only).
- PTH1R: 2 tissues — Artery - Aorta; Artery - Tibial (female_only).

The clearest repeated female-biased gene is CES1 (adrenal, aorta, tibial). PLA2G5 repeats across adrenal and arterial tissues with a male-biased direction. BCL2 is repeated but directionally mixed across tissues. These observations do not establish disease causality.

## Module context

FDR hits occur in **14 of 61 fixed hypertension-branch STRING modules**. The largest descriptive concentrations are:

- HYP_M001: 17 unique hit genes across 5 tissues; 7 female-biased / 14 male-biased.
- HYP_M002: 9 unique hit genes across 4 tissues; 1 female-biased / 9 male-biased.
- HYP_M003: 11 unique hit genes across 6 tissues; 6 female-biased / 10 male-biased.
- HYP_M004: 5 unique hit genes across 3 tissues; 4 female-biased / 3 male-biased.
- HYP_M005: 1 unique hit genes across 1 tissues; 1 female-biased / 0 male-biased.
- HYP_M006: 1 unique hit genes across 1 tissues; 0 female-biased / 1 male-biased.

This is a module-membership overlay on already-estimated gene×tissue results, not a new module enrichment test.

## Relation to the NHANES female-hypertension phenotype

The relevant-tissue readout provides some female-biased hits in aorta, tibial artery, heart, and adrenal tissue, but female-biased hits are fewer than male-biased hits and kidney has no fixed-FDR hits. Therefore M4b does **not** support a simple, globally female-biased molecular hypertension program. It supports localized and directionally mixed sex-by-tissue differences. The kidney result should be treated as unresolved because of limited GTEx female sample size, motivating targeted—but explicitly exploratory—kidney validation in a donor-balanced single-cell dataset.

## Files and guardrails

- `01_hypertension_branch_fdr_hits.csv`: all fixed-FDR hits in the seven relevant tissues.
- `02_hypertension_branch_tissue_summary.csv`: tissue-level and descriptive organ-group counts, including thyroid/liver context rows.
- `03_hypertension_branch_gene_recurrence.csv`: repeated-gene overlay across relevant tissues.
- `04_hypertension_branch_module_context.csv`: descriptive overlay onto fixed M3 STRING modules.

No figures, no new FDR family, no candidate ranking, and no causal or mechanistic interpretation were added.

Generated UTC: 2026-08-30T13:05:28.125134+00:00
