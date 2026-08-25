# MCOP–CRC Supplement Data Revalidation Report

Overall audit status: **PASS**

This report records an independent rerun of the frozen analysis components before Supplement/manuscript submission. It validates computation and data handling; it does not upgrade the cross-sectional NHANES association to causal evidence.

## Audit outcome

- Checks: 11 total; 10 PASS; 1 INFO; 0 FAIL.
- NHANES primary: complete-case N=9,936, CRC cases=70; MCOP per-doubling OR reproduced at approximately 1.2455.
- Independent R `survey::svyglm`: direction and CI conclusion agree with Python; relative logOR change is approximately 4.9×10^-11%.
- 15-axis screen: 15 unique axes, two BH-FDR-supported tests, one Robust Tier A axis.
- 267-chemical actionability matrix: 267 rows and gate/disposition fields reproduced; outcome-firewall audit retained.
- Phase2G: official H5AD local source QC passed for the frozen paired epithelial analysis; staged live Census cache is explicitly not treated as complete evidence.

## Reproducibility boundaries

- The NHANES analysis is cross-sectional and uses current urinary MCOP with prevalent/previous CRC information. Results are associations, not proof of DINP or MCOP causation.
- MCOP is the human biomarker used for the DINP-related exposure axis. The CTD/GeneCards molecular screen nominated the DINP/MiNP axis; it is not a one-to-one CTD mechanism for MCOP.
- The transcriptomic evidence supports CRC epithelial PPAR/NR suppression and RELA/STAT3 activation as disease-state observations. The exposure-to-state arrow remains untested and is not written as mediation or causality.
- The live staged Census run stopped before all 36 donors were successfully cached (35/36 before C136 failure); the final local Phase2G result uses the verified official H5AD source and its local equivalence QC.

## Files

- `logs/reproducibility_audit_results.csv` — machine-readable checks.
- `logs/data_inventory.csv` — input/output/script hashes and sizes.
- `logs/analysis_environment.txt` — software and design metadata.
- `supplementary/MCOP_CRC_Supplement_Tables.xlsx` — Tables S1–S9.

## Submission recommendation

The validated computational components are suitable for Supplement assembly. Before submission, preserve the exact software versions, keep the outcome-firewall language, and retain the Phase2G H5AD source/hash and the R survey cross-check in the audit trail.

## Check details

- **Phase1-2 primary MCOP rerun vs frozen** — PASS: N, events, OR, CI and P reproduced within 1e-8
- **Independent R survey::svyglm replication** — PASS: R OR=1.2455068; Python OR=1.2455068; relative logOR change=4.88e-11%
- **15-axis human screen rerun vs frozen** — PASS: rows=15; all OR/CI/P/FDR values reproduced
- **15-axis screen cardinality** — PASS: n_unique_axes=15
- **15-axis FDR-supported tests** — PASS: n_BH_FDR_lt_0.05=2
- **15-axis robustness rerun vs frozen** — PASS: rows=15; max primary OR difference=3.77e-15
- **Robust Tier A cardinality** — PASS: Robust Tier A=1
- **267-chemical actionability rerun vs frozen** — PASS: rows=267; all gate/disposition fields reproduced
- **Outcome firewall audit** — PASS: Eligibility/actionability selection audited without using CRC OR/P/CI
- **Phase2G official H5AD source QC** — PASS: eligible cells=43801; paired donors=36; source equivalence=True
- **Phase2G staged Census cache completeness** — INFO: staged success=35/36; failures=['C136']; not used as final H5AD source
