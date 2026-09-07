# GO:BP semantic cleaning of the 99 primary significant terms

Generated (UTC): 2026-09-07T14:14:32.873245+00:00

## Scope

This is an additive GO-specific redundancy audit. It does not change the frozen 18-gene query, 881-gene GeneCards background, original P values/FDR values, or the original 106-term output. Reactome terms were not included in this GO cleaning.

## Method

The official `go-basic.obo` is used to recover the GO `is_a` DAG. Terms were ordered from more specific to broader (ontology depth, then term size, then original global FDR). A term was folded into a retained representative only when it was an enriched ancestor/descendant pair with query-gene Jaccard >= **0.50**. Non-ancestor terms were never folded based only on shared query genes, because a small query set can make unrelated terms look similar. All mappings are recorded term by term.

- GO:BP significant terms: **99**
- Retained GO representatives: **43**
- Folded GO terms: **56**
- GO theme families among representatives: **9**

## Clean theme families

| Family | Theme | Retained representatives | Representative sources | Example representative | Union genes |
|---|---|---:|---:|---|---:|
| family_01 | Lipid / fatty-acid / eicosanoid metabolism | 11 | 1 | lipid biosynthetic process | 15 |
| family_02 | Inflammation / defense / stimulus response | 4 | 1 | inflammatory response | 15 |
| family_03 | Chemical / oxygen / abiotic stress response | 6 | 1 | cellular response to chemical stimulus | 16 |
| family_04 | Apoptosis / programmed cell death | 2 | 1 | apoptotic process | 12 |
| family_05 | miRNA transcription / metabolism | 2 | 1 | negative regulation of miRNA transcription | 5 |
| family_06 | Signaling / receptor / hormone response | 4 | 1 | intracellular receptor signaling pathway | 11 |
| family_07 | Matrix remodeling / cell migration | 1 | 1 | regulation of cell migration | 9 |
| family_08 | Development / reproduction / circulation | 10 | 1 | regulation of cell population proliferation | 17 |
| family_09 | General metabolism / biosynthesis / storage | 3 | 1 | regulation of biological quality | 18 |

## Cleaned driver contribution

Driver counts below use retained GO representatives only. They are descriptive repeated-membership counts, not gene-level tests.

| Rank | Gene | Retained GO terms |
|---:|---|---:|
| 1 | PPARD | 35 |
| 2 | PPARG | 32 |
| 3 | PTGS2 | 32 |
| 4 | SIRT1 | 31 |
| 5 | RELA | 25 |
| 6 | ESR1 | 24 |
| 7 | MMP9 | 24 |
| 8 | MMP2 | 23 |
| 9 | STAT3 | 20 |
| 10 | BECN1 | 20 |

## Boundary

A cleaned representative set is not a new enrichment family and does not make the underlying 99 GO terms independent. It is intended to support readable figures/tables and prevent broad GO ancestors from being counted as separate biological stories. The original 99-term results remain the statistical record.
