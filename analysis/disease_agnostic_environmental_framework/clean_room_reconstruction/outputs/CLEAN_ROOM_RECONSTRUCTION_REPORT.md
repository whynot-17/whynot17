# Clean-room Steps 1–4 reconstruction audit

## Scope

This reconstruction was run by neutral modules that do not import the historical disease-project runners. Inputs were limited to the CTD chemical vocabulary, frozen CTD classification rules, DrugCentral drug-exclusion reference, PAH formula guard, the NHANES environmental laboratory catalog, and local environmental XPT files.

No disease outcome, case count, odds ratio, P value, FDR value, disease gene set, disease-specific CTD interaction, or transcriptomic result was loaded.

## Reconstructed counts

- Environmental chemical entities: **2,042 candidates processed; 411 mapping memberships represented in the frozen test table**.
- Registry-backed mapped chemical–biomarker mappings: **449**.
- Actionable chemical–biomarker mappings: **411**.
- Unique human-measurable tests: **29**.

## Comparison against the locked outputs

| Stage | Expected rows/entities | Clean-room rows/entities | Key sets identical | Field values identical | Full-file hash identical |
|---|---:|---:|---|---|---|
| environmental universe | 2042/2042 | 2042/2042 | True | True | False |
| mapped chemical-biomarker mappings | 449/2042 | 449/2042 | True | False | False |
| actionable chemical-biomarker mappings | 411/2042 | 411/2042 | True | False | False |
| unique human test family | 29/29 | 29/29 | True | True | True |
| NHANES registry core rows used by locked outputs | 293/38 | 293/38 | True | True | False |

## Non-analytic provenance differences

The clean-room outputs intentionally use neutral, compact schemas, so full-file hashes for the upstream universe/mapping/actionability tables are not expected to match legacy enriched tables. Their key sets and core universe fields do match.

One legacy metadata label differs for a non-actionable PFAS candidate (`C479228`): the locked mapping table says `unresolved_registry_gap`, whereas the clean-room rule evaluation says `resolved_no_candidate_specific_analyte`. This row is absent from the actionable mapping set and does not change the 411 actionable mappings or the 29-test family.

The full 179,672-row clean-room classification ledger is retained locally for audit reruns; its SHA-256 is recorded in the manifest but the large derived ledger is not part of the version-controlled result bundle.

## Provenance boundary

The clean-room reconstruction demonstrates reproducible execution under an outcome-free input contract. It does not establish that every mapping rule was historically invented before any disease project was seen; that historical-development claim remains intentionally unmade.

Generated UTC: 2026-08-27T16:26:15.098264+00:00

## Input hashes

- ctd: `9e4b642c8716140d30a9376d6b2229acb81ee48a8b106a4f41ed06b29894ef6c`
- rules: `f565cd741ee32ea0683d70e5cf36ea1d13062a2c92bd0b8c0a430f736c057ee6`
- drugcentral: `5b81423a2ec1e2766e9666ec4a172d5a5b47045ea2cc032d1ba06085956bc1fc`
- pah_formulas: `c9057cabc8c3cb740fcbf56e0224431be860bea450083a026d374715e41634eb`
- catalog: `d38569b0b0225b9da4848e32160e31af9667c8e0959743bd799bbca01b28a8f9`
- runner: `0710f2114f55497517218106b12af78923ddcc8fe62022bd7d26ab46ac51c1ca`
