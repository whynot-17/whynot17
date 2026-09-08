# DINP–CRC GO/KEGG enrichment de-redundancy and ranking log

- Generated locally: `2026-09-08T20:55:15.486754+08:00`
- Input: `DINP_CRC_GO_KEGG_enrichment.csv`; SHA-256 `7704caadf2e86091c478110a38596a0b8a5fbf8f88b25aa8512dd55173eef8f4`
- STRING-mapped input proteins used for enrichment: `96`
- Original STRING enrichment table was not modified.

## Reduction method

1. Exact duplicates were removed within each analysis category and term ID.
2. Terms were sorted deterministically by STRING FDR ascending, then input hit-gene count descending, then background term size ascending, then term ID.
3. Within each category, a greedy representative was retained. A later term was collapsed when its input-gene set had Jaccard similarity >= `0.5` OR overlap coefficient >= `0.8` with an already retained representative.
4. The retained representative is the most significant term under this deterministic order. Each collapsed term keeps its representative term, overlap metrics, cluster size, and cluster membership.
5. Ranking fields: `reduced_rank` ranks retained representatives within GO-BP, GO-MF, GO-CC, or KEGG; `neg_log10_fdr`, `gene_coverage_percent`, and `rank_score` provide additional sorting context.

This is gene-set-overlap redundancy reduction, not a GO-DAG semantic-similarity algorithm. It is explicit and reproducible from the STRING result table alone; use the full ranked file if a semantic-similarity/REVIGO-style reduction is preferred later.

## Count reduction

| Category | Input terms | Representatives | Collapsed | Significant representatives (FDR<0.05) |
|---|---:|---:|---:|---:|
| GO_BP | 1228 | 49 | 1179 | 49 |
| GO_CC | 41 | 9 | 32 | 9 |
| GO_MF | 65 | 4 | 61 | 4 |
| KEGG | 129 | 48 | 81 | 48 |

## Representative outputs

### GO_BP

| Rank | Term | Description | Hit genes | FDR |
|---:|---|---|---:|---:|
| 1 | `GO:0070887` | Cellular response to chemical stimulus | 70 | 6.12e-36 |
| 2 | `GO:0051128` | Regulation of cellular component organization | 42 | 1.25e-12 |
| 3 | `GO:1901362` | Organic cyclic compound biosynthetic process | 30 | 1.25e-12 |
| 4 | `GO:0008202` | Steroid metabolic process | 17 | 2.24e-12 |
| 5 | `GO:0050790` | Regulation of catalytic activity | 40 | 2.98e-11 |
| 6 | `GO:0002683` | Negative regulation of immune system process | 19 | 4.54e-11 |
| 7 | `GO:0150077` | Regulation of neuroinflammatory response | 9 | 4.54e-11 |
| 8 | `GO:0009894` | Regulation of catabolic process | 25 | 6.69e-10 |
| 9 | `GO:0050708` | Regulation of protein secretion | 13 | 2.39e-08 |
| 10 | `GO:0001525` | Angiogenesis | 14 | 4.8e-08 |

### GO_CC

| Rank | Term | Description | Hit genes | FDR |
|---:|---|---|---:|---:|
| 1 | `GO:0005615` | Extracellular space | 46 | 2.15e-09 |
| 2 | `GO:0005667` | Transcription regulator complex | 15 | 2.73e-05 |
| 3 | `GO:0012505` | Endomembrane system | 45 | 0.0002 |
| 4 | `GO:0009897` | External side of plasma membrane | 12 | 0.0002 |
| 5 | `GO:0045121` | Membrane raft | 8 | 0.0119 |
| 6 | `GO:0005829` | Cytosol | 43 | 0.0127 |
| 7 | `GO:0005739` | Mitochondrion | 19 | 0.0219 |
| 8 | `GO:0097413` | Lewy body | 2 | 0.0365 |
| 9 | `GO:0034399` | Nuclear periphery | 5 | 0.0389 |

### GO_MF

| Rank | Term | Description | Hit genes | FDR |
|---:|---|---|---:|---:|
| 1 | `GO:0005515` | Protein binding | 78 | 3.14e-15 |
| 2 | `GO:1901363` | Heterocyclic compound binding | 51 | 0.00042 |
| 3 | `GO:0005324` | Long-chain fatty acid transporter activity | 3 | 0.0115 |
| 4 | `GO:0001618` | Virus receptor activity | 4 | 0.043 |

### KEGG

| Rank | Term | Description | Hit genes | FDR |
|---:|---|---|---:|---:|
| 1 | `hsa05200` | Pathways in cancer | 26 | 1.45e-16 |
| 2 | `hsa05321` | Inflammatory bowel disease | 12 | 1.41e-13 |
| 3 | `hsa05418` | Fluid shear stress and atherosclerosis | 14 | 8.91e-13 |
| 4 | `hsa04932` | Non-alcoholic fatty liver disease | 14 | 3.24e-12 |
| 5 | `hsa04920` | Adipocytokine signaling pathway | 11 | 8.96e-12 |
| 6 | `hsa04060` | Cytokine-cytokine receptor interaction | 16 | 5.85e-11 |
| 7 | `hsa03320` | PPAR signaling pathway | 10 | 4.96e-10 |
| 8 | `hsa05332` | Graft-versus-host disease | 8 | 1.47e-09 |
| 9 | `hsa05164` | Influenza A | 12 | 1.91e-09 |
| 10 | `hsa04933` | AGE-RAGE signaling pathway in diabetic complications | 10 | 2.95e-09 |

## Files

- `DINP_CRC_GO_KEGG_reduced_ranked.csv`: all original GO/KEGG terms with representative assignment, rank, overlap metrics, and cluster membership.
- `DINP_CRC_GO_KEGG_representatives.csv`: compact nonredundant representative table.
- `DINP_CRC_GO_representatives.csv`: GO-only representatives.
- `DINP_CRC_KEGG_representatives.csv`: KEGG-only representatives.

## Interpretation boundary

Ranking is an evidence/prioritization order based on STRING FDR and input coverage; it is not a pathway activation score. GO terms are hierarchical and overlapping, so representatives should be treated as topic-level summaries rather than independent discoveries.
