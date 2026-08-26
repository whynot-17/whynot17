# Disease-Agnostic Environmental Biomarker Discovery Framework

## Pre-disease protocol lock: Steps 1–4

This workspace separates the environment-exposure screening engine from any disease plug-in. Steps 1–4 are outcome-blinded and must be completed before any disease outcome, disease gene set, or disease-specific molecular result is loaded.

### Frozen rules

1. The environmental universe is defined from the CTD chemical vocabulary using the frozen CTD MeSH TreeNumber/ParentTreeNumber branches in `work/environmental_toxicology_crc_phase1/chemical_class_rules.json`.
2. Drug-like entries are excluded only by the prespecified CTD/DrugCentral semantic rules in that file; no disease information is used.
3. Chemical-to-biomarker mapping uses the CDC NHANES laboratory catalog and downloaded laboratory XPT files. Parent chemicals, metabolites, family proxies, and analytes remain distinct until Step 4.
4. Detectability and cycle coverage are infrastructure properties. They are assessed without reading cancer or other disease outcomes.
5. A chemical–biomarker mapping is actionable when chemical/analyte identity, human-exposure interpretability, NHANES availability, detectability, cycle coverage, and survey-design infrastructure all pass the fixed gates in `actionability_rules.json`.
6. Step 4 collapses only duplicate human tests. The mapping ledger remains chemical–biomarker level.
7. The resulting unique test set and its size are frozen in `PRE_DISEASE_TESTSET_LOCK.json`. The planned downstream multiplicity denominator is the frozen unique-test count.

### Explicit firewall

The Step 1–4 runner must not load CRC status, case counts, odds ratios, P values, FDR values, GeneCards disease genes, CTD × disease-gene overlap, PPAR results, or any historical candidate effect estimate. Historical results may be compared descriptively only after this lock and never used to change eligibility.

`DISEASE INFORMATION USED IN STEPS 1–4: NO`
