# 41-gene DINP evidence and target-prioritization audit

Generated (UTC): 2026-09-07T16:34:27.724991+00:00

## Scope

This audit decomposes the fresh 41-gene DINP–CRC intersection into CRC reference relevance, DINP source support, and GO/Reactome/KEGG pathway recurrence. It does not change the intersection, recompute enrichment, or introduce a hidden composite score.

- Intersection genes: **41**
- Genes supported by at least two DINP source families (CTD, CompTox, T3DB): **1**
- Genes classified as strong cross-resource pathway anchors: **18**

## DINP source-family support

ToxCast and Tox21 are shown separately and also collapsed into one CompTox family flag for source-family counting. The raw source count is retained and is not interpreted as independent biological truth.

| DINP source-family count | Genes |
|---:|---:|
| 1 | 40 |
| 3 | 1 |

## Why RXRA ranks first

RXRA has CRC article rank **1244**, CRC relevance score **14.84**, DINP support from **CTD** (1 source families), and pathway recurrence of **47** terms across **3** resources. Its pathway recurrence is a breadth proxy, not an annotation-density correction; this audit does not claim that RXRA is the strongest direct DINP target.

## Prioritization decomposition

The original ranking is pathway-context prioritization. It is best read together with the independent DINP evidence profile and CRC relevance columns below.

| Rank | Gene | CRC rank | CRC score | DINP families | DINP profile | GO | Reactome | KEGG | Pathway anchor |
|---:|---|---:|---:|---:|---|---:|---:|---:|---|
| 1 | **RXRA** | 1244 | 14.84 | 1 | CTD | 15 | 28 | 4 | strong_cross-resource_pathway_anchor |
| 2 | **CYP2C9** | 1195 | 15.28 | 1 | CompTox | 17 | 15 | 8 | strong_cross-resource_pathway_anchor |
| 3 | **PPARG** | 67 | 59.83 | 1 | CTD | 50 | 16 | 3 | strong_cross-resource_pathway_anchor |
| 4 | **CYP1A2** | 935 | 17.98 | 1 | CompTox | 21 | 16 | 6 | strong_cross-resource_pathway_anchor |
| 5 | **CYP3A4** | 664 | 22.25 | 1 | CompTox | 20 | 12 | 6 | strong_cross-resource_pathway_anchor |
| 6 | **PPARA** | 1086 | 16.42 | 1 | CTD | 44 | 19 | 2 | strong_cross-resource_pathway_anchor |
| 7 | **PTGS2** | 143 | 46.06 | 1 | CTD | 34 | 11 | 4 | strong_cross-resource_pathway_anchor |
| 8 | **AKR1C3** | 1852 | 10.28 | 1 | CTD | 30 | 10 | 1 | strong_cross-resource_pathway_anchor |
| 9 | **PPARD** | 696 | 21.45 | 1 | CTD | 41 | 9 | 2 | strong_cross-resource_pathway_anchor |
| 10 | **HPGD** | 290 | 34.93 | 1 | CTD | 19 | 9 | 1 | strong_cross-resource_pathway_anchor |
| 11 | **PTGS1** | 517 | 26.08 | 1 | CTD | 11 | 8 | 2 | strong_cross-resource_pathway_anchor |
| 12 | **PTGES2** | 1392 | 13.61 | 1 | CTD | 7 | 6 | 1 | strong_cross-resource_pathway_anchor |
| 13 | **RELA** | 165 | 44.02 | 1 | CTD | 27 | 12 | 3 | strong_cross-resource_pathway_anchor |
| 14 | **STAT3** | 52 | 64.86 | 1 | CTD | 29 | 9 | 2 | strong_cross-resource_pathway_anchor |
| 15 | **ESR1** | 49 | 65.77 | 1 | CTD | 25 | 9 | 1 | strong_cross-resource_pathway_anchor |
| 16 | **MMP9** | 198 | 41.44 | 1 | CTD | 19 | 5 | 2 | strong_cross-resource_pathway_anchor |
| 17 | **SQSTM1** | 720 | 21.05 | 1 | CTD | 15 | 9 | 1 | cross-resource_pathway_supported |
| 18 | **MMP2** | 176 | 43.13 | 1 | CTD | 21 | 6 | 1 | cross-resource_pathway_supported |
| 19 | **PGR** | 563 | 24.73 | 1 | CTD | 14 | 12 | 1 | strong_cross-resource_pathway_anchor |
| 20 | **HSPA1A** | 1408 | 13.46 | 1 | CTD | 24 | 4 | 1 | strong_cross-resource_pathway_anchor |

## Interpretation boundary

RXRA, PPARA, PPARG, PPARD, PTGS2, and related genes should be described as cross-resource pathway anchors or prioritized candidates only after this decomposition. A high GO/Reactome/KEGG recurrence can reflect annotation breadth. Direct DINP binding, causal mediation, and in-vivo relevance require independent evidence and are not established here.

## Files

- `dinp_crc_41_gene_evidence_matrix.csv`: full source/relevance/pathway decomposition for all 41 genes.
- `prioritization_decomposition.csv`: compact table for reviewer-facing rank decomposition.
