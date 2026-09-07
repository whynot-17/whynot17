# Expanded DINP--CRC ORA and GeneCards background audit

Audit status: **PASS**

## Executive findings

- The sensitivity input is 93 gene symbols, but g:Profiler maps it to 100 unique Ensembl IDs.
- The complete one-to-many expansion is `SNORA72` to 8 Ensembl IDs; the other sensitivity symbols map one-to-one and none fail.
- Fixed GMT term membership, N/K/n/k metrics, historical local raw P values, and g:Profiler adjusted P values were checked source-by-source.
- g:Profiler `p_value` is a multiple-testing-adjusted value; the local raw hypergeometric P is compared to a source-wise BH reconstruction rather than to that adjusted field directly.

## Domain accounting

| Background | Input symbols | Canonical Ensembl IDs | GO:BP annotated N | Reactome annotated N |
|---|---:|---:|---:|---:|
| primary_genecards881 | 881 | 868 | 569 | 387 |
| sensitivity_dinp93 | 93 | 100 | 99 | 80 |

## Source-wise recheck

| Background | Source | API rows | API N | Fixed GMT metrics | Historical raw P | API adjusted P reproduced |
|---|---|---:|---:|---|---|---|
| primary_genecards881 | GO:BP | 2005 | 868 | PASS | PASS | PASS |
| primary_genecards881 | REAC | 287 | 868 | PASS | PASS | PASS |
| sensitivity_dinp93 | GO:BP | 2005 | 100 | PASS | PASS | PASS |
| sensitivity_dinp93 | REAC | 287 | 100 | PASS | PASS | PASS |

## GeneCards RelevanceScore >= 10 audit

- Source: `archived ordinary GeneCards CRC top-2000 reference`
- Rows/ranks: 2000 / ranks 1--2000; unique symbols: 2000
- Score >=10: **881**
- Rank 881: `NDRG4` score 10.0
- Rank 882: `FAS` score 9.9
- Rank 2000: `EIF2B5` score 3.1
- Within-file score monotonicity: **PASS**

Within the archived ranked top-2000 source, the threshold crosses from 10.0 at rank 881 to 9.9 at rank 882, and the scores are nonincreasing through rank 2000 (score 3.1). This supports 881 as the complete >=10 set under the source's globally sorted-ranking assumption; it is not evidence that a full login-restricted export was obtained.

## Files

- `local_vs_gprofiler_ora.csv`: one row per API-returned GO:BP/Reactome term, with custom-domain and annotated-domain local calculations.
- `canonical_id_mapping.csv`: symbol-to-Ensembl mappings used in the local calculation.
- `step_gene_background_ora_audit.json`: complete machine-readable audit.
- `genecards_881_threshold_audit.json`: GeneCards threshold proof.

## Boundary

The audit validates the historical `domain_scope=custom` implementation and separately reports an annotated-domain alternative. The official static g:Profiler GMT snapshot used here does not include KEGG; KEGG is therefore outside this fixed-mapping audit and its historical output is untouched. The audit does not silently rewrite the canonical enrichment results or claim that the archived top-2000 file is a full GeneCards export.
