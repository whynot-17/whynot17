# Blood Cd versus urinary Cd × sex → OA

## Scope

This focused biomarker-discrepancy audit uses the same eight NHANES cycles (2003–2004 through 2017–2018), the same adult OA definition (`MCQ160A=1` plus cycle-mapped OA subtype versus `MCQ160A=2` controls), the same age/race/ethnicity/PIR/smoking/cycle fixed effects, and the same survey-weighted interaction framework. It compares `LBXBCD` blood cadmium with frozen `URXUCD` urinary cadmium. It does not add outcomes, create a new FDR family, or perform mechanism analysis.

Urinary Cd retains the established sex-specific urinary-creatinine adjustment (`log2(URXUCR)` plus its female interaction). Blood Cd uses the corresponding NHANES MEC exam weight (`WTMEC2YR`) and does not receive a urinary-dilution term. Urine Cd uses the analyte-specific urine weight (`WTSA2YR`). Each 2-year weight is divided by eight cycles for the pooled fit.

## Model results

| biomarker | interaction β (female−male) | SE | P | male β | female β | analytic N | OA cases | cycles |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| urinary_cadmium | 0.17474 | 0.07473 | 0.020975 | -0.24978 | -0.091685 | 10606 | 1309 | 2003-2004;2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 |
| blood_cadmium | 0.061358 | 0.054253 | 0.26026 | -0.12491 | -0.17865 | 26857 | 3288 | 2003-2004;2005-2006;2007-2008;2009-2010;2011-2012;2013-2014;2015-2016;2017-2018 |

## Direct readout

Urinary Cd has a positive pooled interaction (β=0.17474, P=0.020975), whereas blood Cd has a smaller positive interaction whose confidence interval includes zero (β=0.061358, P=0.26026). Thus this same-cycle analysis supports a biomarker-dependent difference in the Cd × sex signal, not a simple replication of one common OA association.

The descriptive sex-stratified slopes are negative for both sexes under both biomarker-specific models (urinary Cd: male β=-0.24978, female β=-0.091685; blood Cd: male β=-0.12491, female β=-0.17865). These slopes are descriptive and, because nuisance coefficients can differ between stratified and pooled fits, are not used to override the formal pooled interaction.

## Cycle and overlap audit

All 8 planned cycles were available for both exposure matrices after supplementing the local package with the corresponding official public blood-metal XPT files. The OA-eligible blood/urine SEQN overlap is shown in `02_biomarker_cycle_availability.csv`; biomarker-specific analytic N can differ because laboratory availability and missingness differ.

## Interpretation boundary

The formal comparison is the pooled Cd × female interaction, not a comparison of one significant and one non-significant sex-specific slope. A positive interaction means the female slope is higher than the male slope under that biomarker-specific model; it does not establish a positive OA risk, causality, or protection. The blood and urine models address different biological measurement windows and have different laboratory weights, so any discrepancy is interpreted as biomarker-dependent sex heterogeneity, not as proof that one matrix is superior.

Run timestamp (UTC): 2026-08-30T14:46:21.226704+00:00
