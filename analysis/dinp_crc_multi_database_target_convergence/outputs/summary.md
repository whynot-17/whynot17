# DINP–CRC multi-database target convergence summary

Generated: `2026-09-02T16:00:36.209510+00:00`

## Design

Exposure sources were kept separate: CTD, EPA CompTox/ToxCast, and ChEMBL. Disease sources were kept separate: GeneCards, DisGeNET, and Open Targets. Source support counts are descriptive and do not represent a merged biological truth score.

Parent DINP was fixed to CTD `C012125`, CASRN `28553-12-0`, and DTXSID `DTXSID4022521`. No related phthalate was substituted.

## Counts

| Quantity | Count |
|---|---:|
| `CTD_human_genes` | 86 |
| `ToxCast_human_genes` | 0 |
| `ChEMBL_human_genes` | 0 |
| `DINP_exposure_union_genes` | 86 |
| `GeneCards_crc_genes` | 2000 |
| `DisGeNET_crc_genes` | 0 |
| `OpenTargets_crc_genes` | 15611 |
| `CRC_disease_union_genes` | 15885 |
| `DINP_CRC_intersection` | 81 |
| `high_confidence_intersection` | 0 |

## Pairwise and three-way overlaps

Overlap counts are reported separately within the exposure side and disease side; an unavailable source contributes an empty set only as a computational placeholder and is not interpreted as biological absence.

| Exposure-side overlap | Genes |
|---|---:|
| `CTD__ToxCast` | 0 |
| `CTD__ChEMBL` | 0 |
| `ToxCast__ChEMBL` | 0 |
| `CTD__ToxCast__ChEMBL` | 0 |

| Disease-side overlap | Genes |
|---|---:|
| `GeneCards__DisGeNET` | 0 |
| `GeneCards__OpenTargets` | 1726 |
| `DisGeNET__OpenTargets` | 0 |
| `GeneCards__DisGeNET__OpenTargets` | 0 |

## Bilateral-support overview

The following is a descriptive ordering of genes present in the all-source intersection, ordered by exposure-side support, then disease-side support, then total source coverage. It is not a biological truth score and is not used to infer causality.

| Gene | Exposure support | Disease support | Total source coverage |
|---|---:|---:|---:|
| `ABCC4` | 1 | 2 | 3 |
| `ATF4` | 1 | 2 | 3 |
| `ATG5` | 1 | 2 | 3 |
| `BECN1` | 1 | 2 | 3 |
| `CA9` | 1 | 2 | 3 |
| `CCL20` | 1 | 2 | 3 |
| `CRP` | 1 | 2 | 3 |
| `CXCR4` | 1 | 2 | 3 |
| `CYP2A6` | 1 | 2 | 3 |
| `DKK1` | 1 | 2 | 3 |
| `ESR1` | 1 | 2 | 3 |
| `HPGD` | 1 | 2 | 3 |
| `HSPA1A` | 1 | 2 | 3 |
| `IGFBP1` | 1 | 2 | 3 |
| `MAP1LC3B` | 1 | 2 | 3 |
| `MMP2` | 1 | 2 | 3 |
| `MMP9` | 1 | 2 | 3 |
| `NEAT1` | 1 | 2 | 3 |
| `NR1I2` | 1 | 2 | 3 |
| `PGR` | 1 | 2 | 3 |

## Interpretation boundaries

- CTD includes 166 human DINP rows in the archived interaction file; co-treatment rows are retained and flagged in the source record table.
- The EPA live CTX bioactivity request is recorded as unavailable when no API key is provided. Public bulk release availability is recorded separately; this is not treated as biological absence.
- ChEMBL returned no reliable exact parent DINP molecule match in the name search; related phthalates were not queried as substitutes.
- DisGeNET API access was not authorized in this run; its zero is not a biological negative.
- GeneCards is the archived ordinary CRC top-2000 export. The historical strict scoped CRC file is not used in the primary matrix.
- Open Targets uses the source-native disease concept resolved from the exact search hit `MONDO_0005575` (colorectal cancer).

## Files

- `dinp_exposure_gene_matrix.csv`: exposure-side source-preserving matrix.
- `crc_gene_matrix.csv`: disease-side source-preserving matrix.
- `dinp_crc_intersection.csv`: all genes supported on both sides by at least one accessible source.
- `high_confidence_intersection.csv`: prespecified >=2 exposure sources and >=2 disease sources.
- `intersection_subsets.csv`: requested descriptive subsets.
- `source_manifest.json` and `source_records/`: source versions, hashes, requests, and raw/flattened records.
- `source_records/hgnc_symbol_audit.csv`: HGNC approval/status audit for the exposure/intersection symbol universe.
