# Primary 106-term redundancy reduction and 18-gene driver contribution

Generated (UTC): 2026-09-07T14:04:50.852390+00:00

## Frozen scope

- Input: the existing primary GO:BP/KEGG/Reactome enrichment output only.
- Query: the frozen 18-gene DINP–CRC intersection; no genes were added or removed.
- Background: the archived GeneCards CRC high-relevance set with RelevanceScore >= 10 (881 genes).
- Primary term universe: exactly 106 terms with global source-family BH-FDR < 0.05.
- The sensitivity background (93 genes; 0 globally significant terms) was not mixed into this reduction.

## Term-to-gene evidence

Term membership was reconstructed from a separate g:Profiler request with `no_evidences=false`. The response provides one evidence list per query gene; non-empty evidence denotes membership. All 106 terms passed the validation that the evidence-derived overlap count equals the stored `intersection_size`.

## Reduction method

Average-linkage agglomeration used a fixed similarity threshold of **0.58**. Pairwise similarity was 0.65 × query-gene Jaccard + 0.20 × normalized term-token Jaccard + 0.10 × parent-set Jaccard + 0.05 × direct parent/child linkage. The result is a descriptive redundancy reduction, not a new enrichment test and not a causal model.

- 106 significant terms -> **42 fine-grained redundancy clusters** -> **10 broad theme families**.
- Sources represented: GO:BP, REAC.

## Broad theme families

The family layer is the main interpretation view. It is a transparent, disjoint grouping of terms using the fixed term-name rules in the script; the fine-grained similarity clusters remain available in `primary_theme_clusters.csv` and `primary_cluster_term_membership.csv`.

| Family | Theme | Terms | Sources | Representative term | Union genes |
|---|---|---:|---:|---|---:|
| family_01 | Lipid / fatty-acid / eicosanoid metabolism | 31 | 2 | fatty acid metabolic process | 15 |
| family_02 | Inflammation / defense / stimulus response | 17 | 1 | inflammatory response | 18 |
| family_03 | Chemical / oxygen / abiotic stress response | 10 | 1 | cellular response to chemical stimulus | 16 |
| family_04 | Apoptosis / programmed cell death | 7 | 1 | regulation of apoptotic process | 12 |
| family_05 | miRNA transcription / metabolism | 6 | 1 | regulation of miRNA metabolic process | 6 |
| family_06 | Signaling / receptor / hormone response | 9 | 2 | intracellular receptor signaling pathway | 13 |
| family_07 | Matrix remodeling / cell migration | 2 | 2 | Activation of Matrix Metalloproteinases | 9 |
| family_08 | Development / reproduction / circulation | 18 | 1 | regulation of cell population proliferation | 18 |
| family_09 | General metabolism / biosynthesis / storage | 5 | 1 | regulation of biological quality | 18 |
| family_10 | Ontology umbrella / root term | 1 | 1 | REACTOME root term | 17 |

## Fine-grained redundancy clusters

| Cluster | Theme | Terms | Sources | Representative term | Union genes |
|---|---|---:|---:|---|---:|
| theme_01 | Inflammation / defense / stimulus response | 2 | 1 | inflammatory response | 11 |
| theme_02 | Inflammation / defense / stimulus response | 2 | 1 | defense response | 12 |
| theme_03 | Chemical / oxygen / abiotic stress response | 6 | 2 | cellular response to chemical stimulus | 18 |
| theme_04 | Development / reproduction / circulation | 2 | 1 | regulation of cell population proliferation | 14 |
| theme_05 | Lipid / fatty-acid / eicosanoid metabolism | 6 | 1 | fatty acid metabolic process | 8 |
| theme_06 | Lipid / fatty-acid / eicosanoid metabolism | 2 | 1 | lipid metabolic process | 9 |
| theme_07 | Inflammation / defense / stimulus response | 2 | 1 | regulation of response to external stimulus | 12 |
| theme_08 | Development / reproduction / circulation | 7 | 1 | regulation of multicellular organismal process | 16 |
| theme_09 | Development / reproduction / circulation | 1 | 1 | embryo implantation | 5 |
| theme_10 | Lipid / fatty-acid / eicosanoid metabolism | 4 | 1 | fatty acid biosynthetic process | 5 |
| theme_11 | Development / reproduction / circulation | 3 | 1 | female pregnancy | 6 |
| theme_12 | Signaling / receptor / hormone response | 1 | 1 | intracellular receptor signaling pathway | 7 |
| theme_13 | General metabolism / biosynthesis / storage | 1 | 1 | regulation of biological quality | 12 |
| theme_14 | Development / reproduction / circulation | 4 | 1 | positive regulation of response to stimulus | 15 |
| theme_15 | Inflammation / defense / stimulus response | 2 | 1 | response to stress | 15 |
| theme_16 | General metabolism / biosynthesis / storage | 1 | 1 | catabolic process | 11 |
| theme_17 | miRNA transcription / metabolism | 6 | 1 | regulation of miRNA metabolic process | 6 |
| theme_18 | Lipid / fatty-acid / eicosanoid metabolism | 1 | 1 | Metabolism of lipids | 6 |
| theme_19 | Lipid / fatty-acid / eicosanoid metabolism | 1 | 1 | Fatty acid metabolism | 4 |
| theme_20 | Development / reproduction / circulation | 2 | 1 | blood circulation | 6 |
| theme_21 | Signaling / receptor / hormone response | 4 | 1 | regulation of signal transduction | 13 |
| theme_22 | Lipid / fatty-acid / eicosanoid metabolism | 6 | 1 | prostaglandin biosynthetic process | 3 |
| theme_23 | Chemical / oxygen / abiotic stress response | 2 | 2 | cellular response to UV-A | 3 |
| theme_24 | Lipid / fatty-acid / eicosanoid metabolism | 2 | 1 | long-chain fatty acid metabolic process | 4 |
| theme_25 | Lipid / fatty-acid / eicosanoid metabolism | 3 | 1 | regulation of lipid storage | 4 |
| theme_26 | Chemical / oxygen / abiotic stress response | 5 | 1 | response to hypoxia | 7 |
| theme_27 | Signaling / receptor / hormone response | 1 | 1 | hormone-mediated signaling pathway | 5 |
| theme_28 | Chemical / oxygen / abiotic stress response | 2 | 1 | response to lipid | 11 |
| theme_29 | Lipid / fatty-acid / eicosanoid metabolism | 2 | 1 | unsaturated fatty acid metabolic process | 4 |
| theme_30 | Signaling / receptor / hormone response | 1 | 1 | negative regulation of intracellular signal transduction | 8 |
| theme_31 | Inflammation / defense / stimulus response | 1 | 1 | response to abiotic stimulus | 10 |
| theme_32 | Apoptosis / programmed cell death | 7 | 1 | regulation of apoptotic process | 12 |
| theme_33 | Lipid / fatty-acid / eicosanoid metabolism | 1 | 1 | steroid metabolic process | 5 |
| theme_34 | Inflammation / defense / stimulus response | 4 | 1 | response to external biotic stimulus | 9 |
| theme_35 | Chemical / oxygen / abiotic stress response | 1 | 1 | cellular response to chemical stress | 6 |
| theme_36 | General metabolism / biosynthesis / storage | 1 | 1 | primary metabolic process | 16 |
| theme_37 | Signaling / receptor / hormone response | 1 | 1 | regulation of apoptotic signaling pathway | 7 |
| theme_38 | Lipid / fatty-acid / eicosanoid metabolism | 1 | 1 | cholesterol metabolic process | 3 |
| theme_39 | Matrix remodeling / cell migration | 1 | 1 | regulation of cell migration | 9 |
| theme_40 | Lipid / fatty-acid / eicosanoid metabolism | 1 | 1 | lipid localization | 5 |
| theme_41 | Signaling / receptor / hormone response | 1 | 1 | Nuclear Receptor transcription pathway | 3 |
| theme_42 | Lipid / fatty-acid / eicosanoid metabolism | 2 | 1 | Biosynthesis of DHA-derived SPMs | 3 |

## Theme-family driver contribution

For each family, genes are ranked by the number of significant terms in that family containing the gene. A tie is resolved by the overall contribution rank. `family_top_driver` is descriptive and does not denote a gene-level test.

| Family | Theme | Top contributing genes (gene:term count) |
|---|---|---|
| family_01 | Lipid / fatty-acid / eicosanoid metabolism | PTGS2:26; SIRT1:24; PTGS1:23; CYP1A2:21; CYP3A4:18; PPARD:15 |
| family_02 | Inflammation / defense / stimulus response | PTGS2:17; PPARD:17; PPARG:17; RELA:17; MMP9:17; ESR1:16 |
| family_03 | Chemical / oxygen / abiotic stress response | PPARG:9; SIRT1:9; BECN1:9; PTGS2:8; PPARD:8; MMP2:8 |
| family_04 | Apoptosis / programmed cell death | PTGS2:7; PPARD:7; SIRT1:7; ESR1:7; RELA:7; MMP9:7 |
| family_05 | miRNA transcription / metabolism | PPARD:6; PPARG:6; ESR1:6; RELA:6; STAT3:4; NEAT1:2 |
| family_06 | Signaling / receptor / hormone response | PPARG:9; SIRT1:8; ESR1:8; RELA:8; PPARD:7; PTGS2:6 |
| family_07 | Matrix remodeling / cell migration | MMP9:2; MMP2:2; TIMP1:2; PTGS2:1; PPARD:1; PPARG:1 |
| family_08 | Development / reproduction / circulation | PTGS2:18; PPARD:18; MMP2:18; MMP9:16; TIMP1:16; ESR1:15 |
| family_09 | General metabolism / biosynthesis / storage | PPARD:5; SIRT1:5; PPARG:4; CYP1A2:4; CYP3A4:4; PTGS2:3 |
| family_10 | Ontology umbrella / root term | PTGS2:1; PPARD:1; PPARG:1; SIRT1:1; ESR1:1; RELA:1 |

## Broadest recurrent drivers

These are ranked by the number of the 106 significant terms containing each gene; this is a contribution count, not an independent gene-level significance claim.

| Rank | Gene | Significant terms | Fine clusters | Theme families | GO:BP | Reactome |
|---:|---|---:|---:|---:|---:|---:|
| 1 | PTGS2 | 87 | 32 | 9 | 82 | 5 |
| 2 | PPARD | 85 | 33 | 10 | 81 | 4 |
| 3 | PPARG | 77 | 30 | 10 | 74 | 3 |
| 4 | SIRT1 | 74 | 27 | 9 | 73 | 1 |
| 5 | ESR1 | 61 | 22 | 9 | 59 | 2 |
| 6 | RELA | 58 | 21 | 9 | 57 | 1 |
| 7 | MMP9 | 56 | 21 | 9 | 54 | 2 |
| 8 | BECN1 | 53 | 19 | 8 | 52 | 1 |

## Interpretation boundary

The 106 terms are not 106 independent biological findings. Most are related GO descendants/ancestors or share the same 18 query genes. The reduced clusters, theme families, and driver counts summarize redundancy and repeated support; they do not add evidence beyond the frozen enrichment analysis. Cluster/family labels are transparent descriptive labels and should not be read as pathway activation or exposure causality.
