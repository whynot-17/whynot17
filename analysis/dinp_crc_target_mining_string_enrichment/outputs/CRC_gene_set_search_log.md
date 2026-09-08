# CRC-related gene set construction audit log

- Generated locally: `2026-09-08T20:26:21.556514+08:00`
- Scope: CRC disease gene set construction plus DINP-master intersection only.
- Downstream intentionally not run: PPI, GO, KEGG, CRC-only overlap extensions, MR, docking, and single-cell analyses.

## Method replicated

The DEP–CRC paper searched GeneCards with `colorectal cancer`, restricted to Homo sapiens, retained GeneCards Relevance score >=10, then unioned the GeneCards, OMIM, and TTD results after deduplication. This implementation uses the paper's public S6 supplementary table as the reproducible snapshot of those post-filter source lists, because the live GeneCards score table was not retrievable in this run and OMIM's API requires an access key.

- Article: [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0343038](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0343038)
- Article supplementary S6 file used: [10.1371/journal.pone.0343038.s011](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0343038.s011&type=supplementary)
- GeneCards live access attempt: HTTP 403 from the public search endpoint; no live score values were fabricated.
- OMIM live access attempt: API key not available in the workspace; no live OMIM rows were fabricated.
- TTD landing page resolves to the current TTD site, but the article snapshot was used for exact method replication rather than silently mixing time points.

## Input and normalization

- DINP input: `C:\Users\21634\Documents\Codex\2026-09-08\codex-dinp-target-mining-crc-mr\outputs\DINP_target_master.csv`; unique symbols in the current master: **286** (the user message referred to 285; the file actually contains 286).
- DINP master SHA-256: `3a1a8d6367330ddafa4bd909434136b9a322ebd8e1770c01f2b1e173a9dba016`
- CRC supplementary SHA-256: `185403c695e13865ca72342c6c5431cd668e35d85f347717723534033249d8c7`
- Source columns used: GeneCards, OMIM, TTD.
- Final CRC output contains only unique human gene symbols resolved against the NCBI Homo sapiens Gene_info snapshot in `work/Homo_sapiens.gene_info.gz`.
- Unique canonical symbols were matched first, then unique NCBI synonyms; a short, explicit alias map was applied only for unambiguous TTD protein aliases (for example HER2→ERBB2, PD-L1→CD274, TrkA→NTRK1). Ambiguous/unmapped names were excluded from the final HGNC-symbol set and listed below.
- GeneCards numeric scores are not present in the S6 table, so the CSV records the paper's filter rule but does not invent per-gene scores.

## Count reconciliation

| Quantity | Count |
|---|---:|
| Article-reported GeneCards rows | 2,651 |
| Article-reported OMIM rows | 31 |
| Article-reported TTD rows | 42 |
| Article-reported CRC union | 2,785 |
| S6 parsed unique GeneCards raw names | 2743 |
| S6 parsed unique OMIM raw names | 35 |
| S6 parsed unique TTD raw names | 101 |
| S6 raw source-column union | 2786 |
| Final normalized human CRC symbols | 2754 |
| DINP master unique symbols | 286 |
| DINP∩CRC overlap | **97** |

The S6 source columns do not exactly reconcile to the prose counts: the parsed source-column counts are larger, and a few TTD entries are protein/disease labels rather than unique gene symbols. This discrepancy is retained rather than hidden. The final set is based on explicit source-column values, human-gene normalization, and transparent exclusion of names without a unique human gene mapping.

## Normalization status counts

- `canonical`: 2845
- `manual_alias`: 17
- `synonym`: 7
- `unmapped`: 10

## Excluded raw names

These names were present in the S6 source columns but had no unique human gene-symbol mapping and were excluded from the final CRC CSV:

- `GeneCards`: `LOC654780` — unmapped; no unique human gene-symbol mapping
- `GeneCards`: `MT-CO1` — unmapped; no unique human gene-symbol mapping
- `GeneCards`: `MT-CYB` — unmapped; no unique human gene-symbol mapping
- `GeneCards`: `MT-ND1` — unmapped; no unique human gene-symbol mapping
- `GeneCards`: `MT-TT` — unmapped; no unique human gene-symbol mapping
- `TTD`: `COVID-19` — unmapped; no unique human gene-symbol mapping
- `TTD`: `Candi TMP1` — unmapped; no unique human gene-symbol mapping
- `TTD`: `Ftase` — unmapped; no unique human gene-symbol mapping
- `TTD`: `GD2` — unmapped; no unique human gene-symbol mapping
- `TTD`: `hDNA` — unmapped; no unique human gene-symbol mapping

## Deliverables

- `CRC_gene_set_GeneCards10_OMIM_TTD.csv`: normalized CRC gene set with per-source presence and raw source names.
- `DINP_CRC_overlap.csv`: exact gene-symbol intersection with the current DINP target master, including DINP evidence metadata.
- This log: method, access limitations, hashes, count reconciliation, and exclusions.

## Interpretation boundary

The overlap is a screening intersection defined by a literature-derived CRC gene set and the existing DINP target master. It is not evidence that DINP directly binds every overlapping protein, and it is not yet an MR instrument/target list. Review the source partition and evidence grade in the overlap CSV before deciding what can enter MR.
