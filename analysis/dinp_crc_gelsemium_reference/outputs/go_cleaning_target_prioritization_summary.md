# Fresh DINP--CRC GO cleaning and target prioritization

Generated (UTC): 2026-09-07T16:25:40.815145+00:00

## Frozen scope

This analysis uses only the published Supplementary Table S2 CRC reference set, the frozen 93-gene DINP multi-source set, and their fresh 41-gene symbol-level intersection. It does not reuse the old 881-gene background or 18-gene query.

- Fresh query genes: **41**
- Significant GO:BP terms before cleaning: **202**
- Significant Reactome terms: **50**
- Significant KEGG terms: **13**
- Retained GO representatives: **97**
- Folded GO terms: **105**

## GO cleaning

GO:BP terms were ordered by GO-DAG depth and reduced only for ancestor/descendant pairs with query-gene Jaccard >= **0.50**. Non-ancestor terms were not collapsed solely because they shared genes.

The original enrichment P values and FDR values remain unchanged. The cleaned representatives are an interpretation/readability layer, not a new statistical family.

## Cross-resource themes

GO representatives, Reactome, and KEGG were assigned to the same frozen term-name themes. Theme overlap is descriptive and does not treat GO, Reactome, and KEGG as independent tests. The unassigned/other bucket is reported for completeness but is excluded from positive target prioritization.

| Theme | GO reps | Reactome | KEGG | GO–Reactome shared genes | GO–KEGG shared genes | 3-source shared genes | Union genes |
|---|---:|---:|---:|---:|---:|---:|---:|
| Lipid / fatty-acid / eicosanoid metabolism | 19 | 8 | 2 | 14 | 13 | 8 | 30 |
| Chemical / xenobiotic / oxygen-stress response | 6 | 3 | 4 | 4 | 9 | 4 | 27 |
| Inflammation / immune / stimulus response | 12 | 4 | 0 | 9 | 0 | 0 | 38 |
| Signaling / receptor / hormone response | 10 | 7 | 1 | 22 | 4 | 4 | 30 |
| Matrix remodeling / adhesion / migration | 0 | 1 | 0 | 0 | 0 | 0 | 4 |
| Cell death / cell cycle / genome maintenance | 2 | 0 | 0 | 0 | 0 | 0 | 21 |
| General metabolism / biosynthesis | 14 | 7 | 2 | 21 | 3 | 3 | 35 |
| Other / mixed biological processes | 34 | 20 | 4 | 36 | 18 | 18 | 40 |

## Target prioritization

The ranking is lexicographic, not a fitted composite score: three-source same-theme support, two-source same-theme support, number of supporting resources, Reactome recurrence, GO recurrence, KEGG recurrence, then theme breadth and gene symbol. This prioritizes reproducible cross-resource anchors without claiming causality.

| Rank | Gene | GO | Reactome | KEGG | Resources | Same-theme 2-source | Same-theme labels |
|---:|---|---:|---:|---:|---:|---:|---|
| 1 | **RXRA** | 15 | 28 | 4 | 3 | 4 | Lipid / fatty-acid / eicosanoid metabolism; Chemical / xenobiotic / oxygen-stress response; Signaling / receptor / hormone response; General metabolism / biosynthesis |
| 2 | **CYP2C9** | 17 | 15 | 8 | 3 | 3 | Lipid / fatty-acid / eicosanoid metabolism; Chemical / xenobiotic / oxygen-stress response; General metabolism / biosynthesis |
| 3 | **PPARG** | 50 | 16 | 3 | 3 | 3 | Lipid / fatty-acid / eicosanoid metabolism; Signaling / receptor / hormone response; General metabolism / biosynthesis |
| 4 | **CYP1A2** | 21 | 16 | 6 | 3 | 3 | Lipid / fatty-acid / eicosanoid metabolism; Chemical / xenobiotic / oxygen-stress response; General metabolism / biosynthesis |
| 5 | **CYP3A4** | 20 | 12 | 6 | 3 | 3 | Lipid / fatty-acid / eicosanoid metabolism; Chemical / xenobiotic / oxygen-stress response; General metabolism / biosynthesis |
| 6 | **PPARA** | 44 | 19 | 2 | 3 | 4 | Lipid / fatty-acid / eicosanoid metabolism; Chemical / xenobiotic / oxygen-stress response; Signaling / receptor / hormone response; General metabolism / biosynthesis |
| 7 | **PTGS2** | 34 | 11 | 4 | 3 | 4 | Lipid / fatty-acid / eicosanoid metabolism; Chemical / xenobiotic / oxygen-stress response; Inflammation / immune / stimulus response; General metabolism / biosynthesis |
| 8 | **AKR1C3** | 30 | 10 | 1 | 3 | 3 | Lipid / fatty-acid / eicosanoid metabolism; Signaling / receptor / hormone response; General metabolism / biosynthesis |
| 9 | **PPARD** | 41 | 9 | 2 | 3 | 3 | Lipid / fatty-acid / eicosanoid metabolism; Signaling / receptor / hormone response; General metabolism / biosynthesis |
| 10 | **HPGD** | 19 | 9 | 1 | 3 | 2 | Lipid / fatty-acid / eicosanoid metabolism; General metabolism / biosynthesis |
| 11 | **PTGS1** | 11 | 8 | 2 | 3 | 2 | Lipid / fatty-acid / eicosanoid metabolism; General metabolism / biosynthesis |
| 12 | **PTGES2** | 7 | 6 | 1 | 3 | 2 | Lipid / fatty-acid / eicosanoid metabolism; General metabolism / biosynthesis |
| 13 | **RELA** | 27 | 12 | 3 | 3 | 5 | Lipid / fatty-acid / eicosanoid metabolism; Chemical / xenobiotic / oxygen-stress response; Inflammation / immune / stimulus response; Signaling / receptor / hormone response; General metabolism / biosynthesis |
| 14 | **STAT3** | 29 | 9 | 2 | 3 | 5 | Lipid / fatty-acid / eicosanoid metabolism; Chemical / xenobiotic / oxygen-stress response; Inflammation / immune / stimulus response; Signaling / receptor / hormone response; General metabolism / biosynthesis |
| 15 | **ESR1** | 25 | 9 | 1 | 3 | 3 | Chemical / xenobiotic / oxygen-stress response; Signaling / receptor / hormone response; General metabolism / biosynthesis |

## Boundary

A recurrent gene is a pathway-context candidate, not a proven DINP target. Database overlap, pathway recurrence, and target prioritization do not establish direct binding, directionality, causality, or in-vivo mediation. Any docking or structural follow-up must remain a separate hypothesis-generating layer.

## Files

- `go_bp_cleaning_1893/`: GO:BP term-level cleaning and representatives.
- `go_reactome_kegg_aligned_1893/`: descriptive cross-resource theme alignment and drivers.
- `target_prioritization_1893/`: current 41-gene target ranking and cross-resource shortlist.
