# DINP--CRC restart using the published 1,893-gene CRC set

## Frozen inputs

- Source: Que et al., BMC Complementary Medicine and Therapies (2021), DOI `10.1186/s12906-021-03273-7`.
- CRC source: official Supplementary Table S2, normalized to 1,893 unique gene symbols.
- DINP source: frozen CTD + ToxCast/Tox21 + T3DB matrix, 93 unique gene symbols.
- Symbol-level intersection: **41 genes**.
- No previous 881-gene GeneCards background or 18-gene query was reused.

## Enrichment design

GO:BP, KEGG, and Reactome were queried separately with g:Profiler using `domain_scope=custom` and `all_results=true`.
The 1,893-gene CRC set is the primary source-defined background; the 93-gene DINP set is a sensitivity background.
g:Profiler may report an effective domain smaller or larger than the input symbol count after canonical Ensembl mapping; those values are retained in the output and manifest.

## Results by background

| Background | Input genes | Effective domain(s) | GO:BP significant | KEGG significant | Reactome significant |
|---|---:|---|---:|---:|---:|
| gelsemium_crc1893 | 1893 | 1875 | 202 | 13 | 50 |
| dinp93 | 93 | 100 | 128 | 0 | 1 |

## Files

- `crc_genecards_1893_normalized.csv`: converted article Table S2.
- `dinp_crc_intersection_41.csv`: 41-gene DINP--CRC intersection with DINP source flags.
- `enrichment_all_backgrounds.csv`: GO:BP/KEGG/Reactome results for both backgrounds.
- `manifest.json`: hashes, API provenance, and effective-domain accounting.
