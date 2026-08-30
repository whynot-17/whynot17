# URXP02 / 2-NAP expanded molecular universe (M1b)

Generated 2026-08-30T07:21:33.162465+00:00. This package maximizes evidence coverage while keeping exact 2-naphthol separate from parent naphthalene. It does not re-prove the NHANES association, infer causality, or make sex-specific molecular claims.

## Evidence tiers

- **Tier A:** exact 2-naphthol human evidence.
- **Tier B:** exact 2-naphthol experimental/non-human evidence.
- **Tier C:** parent naphthalene support; never labelled exact 2-NAP.

The exact CTD identity is `C028405` (2-naphthol; synonyms 2-hydroxynaphthalene/2-NAP). Parent naphthalene is `C031721`. PubChem and ChEMBL identifiers are recorded in every long-table row.

## Size audit

- Exact 2-NAP genes: **110**
- Parent naphthalene genes: **784**
- Unified evidence rows: **2165**
- CTD exact rows: **21**; PubChem exact assay rows: **419**; ChEMBL exact activity rows: **84**
- Parent CTD rows: **662**; parent PubChem assay rows: **862**; parent ChEMBL activity rows: **117**

See `06_molecular_universe_audit.csv` for the full metric table.

## Source and coverage notes

CTD records retain interaction actions, interaction sentences, species, and PubMed identifiers. PubChem BioAssay records retain activity outcome/value, assay identifiers, target GeneID, and NCBI-mapped species/symbol when available. ChEMBL records retain target organism, assay type, standard measurement, target identifier, and component accession/symbol. GEO/BioStudies exact-term search results and unintegrated resource status are recorded in the audit/manifest; no exact expression dataset was silently converted into a gene list.

The current expanded universe is therefore suitable for downstream pathway/network/tissue/cell work, but the small exact-2-NAP human evidence base and the assay-heavy Tier C records should be kept visibly separate in later analyses.

## Restrictions honored

No disease intersection, pathway enrichment, tissue/cell analysis, figures, or new NHANES association model was run.
