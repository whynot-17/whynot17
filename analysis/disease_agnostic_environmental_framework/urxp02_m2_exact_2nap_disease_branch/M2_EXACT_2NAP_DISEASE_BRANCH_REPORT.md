# M2-RESET — exact 2-NAP disease branching

Generated UTC: 2026-08-30T13:35:38.192085+00:00

## Primary rule

> **Parent naphthalene-only evidence is excluded from the primary molecular universe.**

The primary universe is constructed only from M1b rows with chemical `2-naphthol` and the exact identity `CTD:C028405; PubChem CID:8663; ChEMBL:CHEMBL14126` (2-naphthol / 2-hydroxynaphthalene / 2-NAP). The filter yields 524 exact evidence rows and exactly **110 unique gene symbols**. No gene enters because of parent naphthalene evidence alone.

The 110 symbols are retained as supplied by M1b, including source labels that may not be canonical HGNC symbols. They are evidence-linked entities, not 110 proven causal targets of 2-NAP.

## Disease definitions

The original M2 definitions are frozen and reused without broadening:

- Thyroid disease: `MONDO_0003240` (`enableIndirect=false` Open Targets direct associations).
- Hypertension: `HP_0000822` (`enableIndirect=false` Open Targets direct associations).

No thyroid cancer, thyroid hormone traits, blood-pressure GWAS traits, cardiovascular surrogates, or kidney surrogates were introduced.

## Exact-2-NAP branch counts

| Branch | Genes |
|---|---:|
| Thyroid-specific (A − B) | **0** |
| Hypertension-specific (B − A) | **29** |
| Shared (A ∩ B) | **34** |
| Neither | **47** |
| **Total** | **110** |

The mandatory partition identity passes: **0 + 29 + 34 + 47 = 110**.

There is no exact-2-NAP thyroid-specific gene under these frozen disease lists. The exact universe therefore has a hypertension-specific component (29), a shared component (34), and a sizeable neither component (47), but no exact thyroid-only branch.

## Evidence characterization

Evidence flags in the branch summary are recomputed from exact rows only. `exact_experimental_supported` follows the M1b convention of any exact row whose species is not `Homo sapiens` (including source rows with unknown species labels). Evidence roles distinguish `gene_to_chemical_metabolism` from `direct_binding`, `bioassay_target`, and other exact interactions. A metabolism enzyme that changes 2-NAP abundance is not treated as a gene whose expression/activity is changed by 2-NAP.

These records preserve provenance and direction where supplied; they do not establish causality. In particular, this analysis does not call all 110 genes direct targets.

## Comparison with the old parent-expanded 828-gene M2

- Old shared core: **189 → 34 exact shared genes**; **155/189 (82.01%)** are removed under the exact identity restriction.
- Old thyroid intersection: 219 → 34 exact genes (15.53% retained).
- Old hypertension intersection: 440 → 63 exact genes (14.32% retained).
- Old thyroid-specific branch: 30 → **0** exact genes.
- Old hypertension-specific branch: 251 → 29 exact genes (11.55% retained).

Therefore the prior broad 828-gene shared-core interpretation **does not remain valid as the primary exact-2-NAP interpretation**. A 34-gene exact shared subset remains, but the broad 189-gene parent-expanded shared core has collapsed and the exact branch structure is not reciprocal thyroid-vs-hypertension.

## Stop rule

This reset stops after exact-universe construction, disease intersections, branch classification, evidence characterization, and old-vs-exact audit. No enrichment, GO/KEGG/Reactome, PPI, STRING, modules, GTEx, sex-DE, single-cell analysis, immune infiltration, figures, or new NHANES analysis was run.
