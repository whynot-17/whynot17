# Step 10B-P3 — KNHANES disease-feasibility census

## Objective

Identify one to three KNHANES outcomes that can be frozen **before any exposure–outcome association is inspected**, so that KNHANES can serve as an external population-replacement demonstration without outcome shopping.

This is a feasibility census, not an association analysis. It does not import published effect estimates, inspect project P values, or rank outcomes by any exposure result.

## Frozen inputs

- The 29-test environmental exposure family from Step 4 remains unchanged.
- The prior KNHANES exact crosswalk is reused only to report the common exposure-side constraint: 0/29 exact public matches and 2/29 related blood-matrix mismatches.
- Exposure overlap is not used to rank diseases because it is the same frozen constraint for every candidate outcome.

## Candidate outcomes

The census covers common KNHANES health outcomes with clear measurement or derivation pathways:

1. hypertension;
2. diabetes / glycemic disease outcome;
3. dyslipidemia / hypercholesterolemia;
4. metabolic syndrome;
5. chronic kidney disease (CKD);
6. obesity (included as a phenotype contrast, not a primary disease outcome);
7. liver disease (conditional reserve candidate).

## Feasibility-only scoring

Each candidate receives 0–2 points on six pre-association dimensions:

- definition reproducibility;
- expected event/phenotype density;
- same-person availability of outcome components;
- core covariate availability;
- survey-design documentation;
- access feasibility.

The maximum is 12. Exposure overlap is recorded separately and is never included in the score. The primary provisional freeze selects up to three highest-scoring **directly defined disease/clinical-outcome** candidates. Derived multicomponent phenotypes (for example, metabolic syndrome) and non-disease phenotypes (for example, obesity) remain explicitly out of the primary disease freeze, even when their feasibility score is high. The final primary analysis remains conditional on authorized raw-file/codebook confirmation.

## Required access confirmation before association models

For each frozen outcome, authorized KNHANES files must confirm:

- exact variable names and coding;
- outcome component co-observation with the exposure laboratory file;
- missingness and analytic sample counts;
- weight, strata, and PSU fields for the relevant examination/laboratory subsample;
- harmonized cycle window;
- prespecified case/control exclusions.

Until those checks pass, `case_count_now` is recorded as `not_calculable_without_registered_raw_data` and no association model may be run.

## Sources

- KDCA survey overview: https://kdca.go.kr/eng/4428/subview.do
- Official raw-data access record: https://www.data.go.kr/tcs/dss/selectFileDataDetailView.do?publicDataPk=15076556
- Hypertension measurement/definition precedent: https://pmc.ncbi.nlm.nih.gov/articles/PMC4661365/
- Diabetes and glycemic-variable precedent: https://pmc.ncbi.nlm.nih.gov/articles/PMC3678002/
- Dyslipidemia definition and KNHANES laboratory precedent: https://pmc.ncbi.nlm.nih.gov/articles/PMC12488789/
