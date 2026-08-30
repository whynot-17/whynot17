# URXP02 / 2-NAP mechanism gene discovery (M1)

Generated 2026-08-30T07:01:45.269904+00:00. This is a gene-level discovery package only; it does not establish causality, sex-specific molecular mechanisms, or tissue/cell localization.

## Chemical identity and exposure universe

The exact CTD entry is **2-naphthol** (`C028405`; CAS 135-19-3). CTD curated synonyms include 2-hydroxynaphthalene and 2-NAP. Near matches (1-naphthol, 2-naphthyl sulfate, parent naphthalene, and unrelated derivatives) were excluded. The human CTD interaction universe contains **9 genes**.

## Disease gene sources

Disease branches were kept independent. Open Targets direct target–disease associations (`enableIndirect=false`) returned 2079 genes for `MONDO_0003240` (thyroid gland disorder) and 4979 genes for `HP_0000822` (Hypertension). These are ranked association resources, not causal gene lists.

## Intersections

- 2-NAP ∩ thyroid: **3** genes
- 2-NAP ∩ hypertension: **7** genes
- thyroid-specific (A−B): **0**
- hypertension-specific (B−A): **4**
- shared (A∩B): **3**

Full rows and component evidence are in `04`–`07`.

## Enrichment and candidate use

GO Biological Process, KEGG, and Reactome were queried separately for each branch with the **9-gene human CTD 2-NAP universe as the custom background**. The thyroid-specific branch contains no genes, and g:Profiler returned no terms for the four-gene hypertension-specific branch or the three-gene shared branch under this custom background. Thus no pathway enrichment is claimed in M1. The outputs preserve the null/empty results and g:SCS settings; no expected pathway was forced. Because the enrichment response does not provide a defensible per-gene attribution table, `11_cell_mapping_candidate_genes.csv` reports set-level pathway support only and does not claim sex specificity.

Next phase may map these transparent candidates to tissues/cells. It must test, rather than assume, any male/female molecular divergence.

## Provenance

- CTD chemical–gene interactions: `C:\Users\21634\Documents\Codex\2026-08-22\non\work\whynot17\work\environmental_toxicology_crc_phase1\data\CTD_chem_gene_ixns.tsv.gz` (SHA-256 `05e1b0d2d93bb33f72659e6ed2d590304cdb58da9538c787ad51fc6624d0e055`).
- Disease associations: Open Targets Platform GraphQL `https://api.platform.opentargets.org/api/v4/graphql`; direct associations only.
- Enrichment: g:Profiler g:GOSt `https://biit.cs.ut.ee/gprofiler/api/gost/profile/`; custom background, sources GO:BP/KEGG/REAC, g:SCS.
- No NHANES association or exposure–outcome result was read or refit by this M1 script.
