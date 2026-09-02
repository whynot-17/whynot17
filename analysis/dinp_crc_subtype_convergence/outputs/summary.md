# DINP–CRC anatomical subtype convergence summary

Generated: `2026-09-02T16:37:00.052666+00:00`

## Frozen interpretation

The 81-gene accessible DINP–CRC intersection was localized using source-native Open Targets anatomical CRC concepts. This is not a patient-level right-versus-left CRC analysis.

## Input counts

- DINP exposure union: **86 genes**
- General CRC union: **15885 genes**
- Frozen DINP–CRC intersection: **81 genes**

## Open Targets concepts

| Label | Concept | ID | Unique targets | Exact expected concept |
|---|---|---|---:|---|
| `right_ascending` | ascending colon cancer | `MONDO_0002238` | 23 | True |
| `left_sigmoid_strict` | sigmoid colon cancer | `MONDO_0001464` | 64 | True |
| `left_rectosigmoid_sensitivity` | rectosigmoid carcinoma | `MONDO_0002424` | 27 | True |

## 81-gene localization

| Pattern | Genes |
|---|---:|
| `both` | 0 |
| `right_only` | 0 |
| `left_only` | 2 |
| `neither` | 79 |

## Exploratory enrichment

| Comparison | Target genes in background | 81-gene overlap | Fisher OR | P | BH-FDR |
|---|---:|---:|---:|---:|---:|
| `right_ascending` | 23 | 0 | 0 | 1 | 1 |
| `left_sigmoid_strict` | 64 | 2 | 6.428 | 0.04224 | 0.06336 |
| `left_expanded_sigmoid_or_rectosigmoid` | 64 | 2 | 6.428 | 0.04224 | 0.06336 |

## Boundaries

- Anatomical anchor concepts are not equivalent to a clinical right/left-sided CRC phenotype label.
- Source support is descriptive; no heterogeneous source evidence was collapsed into a biological truth score.
- The upstream 81-gene intersection was not modified; these outputs are a subtype-localization follow-up only.

## Files

- `dinp_crc_subtype_gene_support.csv`: 81-gene source-preserving anatomical support table.
- `dinp_crc_subtype_enrichment.csv`: exploratory enrichment calculations.
- `subtype_manifest.json` and `source_records/`: concept resolution, target rows, requests, and hashes.
