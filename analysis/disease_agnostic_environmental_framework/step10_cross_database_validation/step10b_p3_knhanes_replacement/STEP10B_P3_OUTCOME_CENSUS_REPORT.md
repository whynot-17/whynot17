# Step 10B-P3 — KNHANES disease-feasibility census

Generated: `2026-08-28T12:47:13Z`  Status: **pre-association outcome feasibility census complete; provisional freeze only**

## Decision

The predeclared feasibility rubric provisionally freezes three eligible directly defined disease/clinical-outcome candidates before any exposure–outcome association is inspected:

1. **Hypertension** — score 12/12;
2. **Diabetes / glycemic disease outcome** — score 11/12;
3. **Dyslipidemia / hypercholesterolemia** — score 11/12.

This is a **data-feasibility freeze**, not a claim that KNHANES is ready for immediate replication. The authorized raw files and codebooks must still confirm exact fields, same-person exposure/outcome overlap, missingness, harmonized cycle windows, and current weight/strata/PSU variables.

## What was and was not used

- Candidate selection used only outcome definition reproducibility, expected event/phenotype density, same-person component availability, core covariate availability, survey-design documentation, and official access feasibility.
- No exposure–outcome model, P value, FDR, effect estimate, or association direction was read.
- The frozen 29-test exposure crosswalk was carried forward only as a common constraint: **0/29 exact public matches and 2/29 related blood-matrix mismatches** in the prior audit. This exposure overlap was not used to rank outcomes.
- Current case counts are not calculable without registered raw data and are therefore not imputed from publications.

## Candidate census

| Outcome | Class | Score | Role |
|---|---|---:|---|
| Hypertension | clinical disease outcome | 12/12 | primary disease candidate |
| Diabetes / glycemic disease outcome | clinical disease outcome | 11/12 | primary disease candidate |
| Dyslipidemia / hypercholesterolemia | clinical phenotype outcome | 11/12 | primary clinical-outcome candidate |
| Metabolic syndrome | derived clinical phenotype | 11/12 | reserve derived phenotype |
| Chronic kidney disease | clinical disease outcome | 10/12 | reserve disease candidate |
| Obesity / adiposity phenotype | non-disease phenotype contrast | 12/12 | contrast only; not eligible for primary disease freeze |
| Liver disease / liver injury phenotype | conditional clinical phenotype | 9/12 | reserve conditional outcome |

Obesity scored highly but is retained only as a non-disease phenotype contrast. Metabolic syndrome remains a reserve derived phenotype even though its component variables are plausible, because the primary freeze excludes multicomponent derived outcomes. CKD and liver disease remain reserve/conditional outcomes because renal/liver definitions, synchronized components, or codebook-level availability require additional confirmation.

## Why the three provisional candidates are feasible

- **Hypertension:** direct blood-pressure examination plus medication/diagnosis information gives a reproducible outcome pathway and high expected event density. KNHANES blood-pressure measurement and hypertension definitions are documented in published analyses.
- **Diabetes:** glycemic laboratory measures plus diagnosis/medication fields are documented, but T2D-specific coding, type-1 exclusion, and fasting-subsample handling remain pending.
- **Dyslipidemia:** fasting lipid measurements and medication information are documented, but the historical threshold/definition changes require a harmonized cycle window before modeling.

The KDCA survey overview describes health examinations and interviews covering obesity, hypertension, diabetes, dyslipidemia, kidney/liver disease, smoking, drinking, physical activity, and socioeconomic/health information. The official raw-data route is registration/consent based; no personal information was entered and no microdata were downloaded in this census. [KDCA survey overview](https://kdca.go.kr/eng/4428/subview.do) · [official raw-data record](https://www.data.go.kr/tcs/dss/selectFileDataDetailView.do?publicDataPk=15076556)

Published KNHANES precedents support the feasibility of measured hypertension definitions, glycemic/diabetes components, and fasting lipid outcomes, but they are not imported as our external estimates: [hypertension precedent](https://pmc.ncbi.nlm.nih.gov/articles/PMC4661365/) · [diabetes/glycemia precedent](https://pmc.ncbi.nlm.nih.gov/articles/PMC3678002/) · [dyslipidemia precedent](https://pmc.ncbi.nlm.nih.gov/articles/PMC12488789/)

## Gate before association analysis

The three frozen outcomes may enter an external population-replacement analysis only after authorized KNHANES raw data/codebook review confirms:

1. exact outcome variables and coding;
2. exposure laboratory file and outcome components in the same persons;
3. analytic N, event counts, missingness, and cycle harmonization;
4. correct survey weight, strata, and PSU variables;
5. a prespecified outcome definition and exclusions.

Until then, the correct status is **provisional feasibility freeze, access/codebook confirmation pending**.
