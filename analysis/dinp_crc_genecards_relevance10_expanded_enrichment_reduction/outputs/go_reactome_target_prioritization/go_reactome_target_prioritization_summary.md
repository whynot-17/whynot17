# GO–Reactome target prioritization

Generated (UTC): 2026-09-07T14:57:56.184799+00:00

## Boundary

This is an evidence-ordering layer over the existing 18-gene DINP–CRC intersection. No enrichment P values or FDR values were recalculated. Database recurrence is not causality, and structural/docking availability is not evidence of in-vivo binding.

Biological priority is ordered lexicographically by: (1) number of same-theme GO–Reactome concordances, (2) presence in both sources, (3) Reactome term recurrence, (4) GO representative-term recurrence, and (5) total theme breadth. This avoids a hidden weighted composite score.

## Cross-resource biological anchors

- 18 total genes; **10** have GO and Reactome support in the same theme.
- Shared GO–Reactome themes: **3** — lipid/eicosanoid metabolism, nuclear-receptor signaling, and matrix remodeling.

| Rank | Gene | GO | Reactome | Same-theme concordance | Role |
|---:|---|---:|---:|---:|---|
| 1 | **PPARD** | 35 | 3 | 2 | direct lipid/eicosanoid pathway context; nuclear-receptor pathway context |
| 2 | **PPARG** | 32 | 2 | 2 | direct lipid/eicosanoid pathway context; nuclear-receptor pathway context |
| 3 | **PTGS2** | 32 | 4 | 1 | direct lipid/eicosanoid pathway context |
| 4 | **CYP1A2** | 13 | 4 | 1 | direct lipid/eicosanoid pathway context |
| 5 | **CYP3A4** | 10 | 3 | 1 | direct lipid/eicosanoid pathway context |
| 6 | **PTGS1** | 12 | 2 | 1 | direct lipid/eicosanoid pathway context |
| 7 | **ESR1** | 24 | 1 | 1 | nuclear-receptor pathway context |
| 8 | **MMP9** | 24 | 1 | 1 | matrix-remodeling pathway context |

### Interpretation

- **PPARD / PPARG** are the strongest multi-theme cross-resource anchors because they recur in lipid/eicosanoid and nuclear-receptor themes.
- **PTGS2** is the strongest direct lipid/eicosanoid recurrence anchor, with support across all four informative lipid/eicosanoid Reactome terms and extensive GO recurrence.
- **MMP9 / MMP2 / TIMP1** form the matrix-remodeling cross-resource branch; MMP9 additionally has macrophage and STRING-network support.
- **CYP1A2 / CYP3A4 / PTGS1** are direct lipid/eicosanoid pathway-context candidates, but their prioritization is pathway-contextual rather than proof of a DINP-specific target interaction.

## Structural follow-up branch

The structural branch is reported separately because the existing macrophage/STRING workflow covers seven candidates and is not the same universe as the 18-gene GO–Reactome alignment.

| Follow-up rank | Gene | Macrophage driver rank | STRING/network rank | Biological source status | Role |
|---:|---|---:|---:|---|---|
| 1 | **MMP9** | 2.0 | 1 | GO+Reactome_same_theme | network bridge candidate |
| 2 | **STAT3** | 4.0 | 2 | GO_only | network bridge candidate |
| 3 | **TIMP1** | 3.0 | 3 | GO+Reactome_same_theme | supporting/context |
| 4 | **CXCR4** | 7.0 | 4 | GO_only | supporting/context |
| 5 | **NEAT1** | 1.0 | 5 | GO_only | non-protein state/regulatory node |
| 6 | **PTGER4** | 5.0 | 5 | outside_aligned_18_gene_universe | direct prostaglandin node |
| 7 | **PTGES3** | 6.0 | 5 | outside_aligned_18_gene_universe | direct prostaglandin node |

## Recommended use

For the biological mechanism narrative, carry forward the cross-resource anchors first: **PPARD/PPARG–PTGS2**, with the **MMP9/MMP2/TIMP1 matrix branch** as a complementary tumor-microenvironment axis.

For protein-level follow-up, retain MMP9 and STAT3 as network-bridge candidates and PTGER4/PTGES3 as direct prostaglandin-context candidates. NEAT1 remains a non-protein state node and is not docking-eligible.

These are prioritization recommendations only. They do not establish that DINP binds any target or that any target mediates the epidemiologic association.

## Files

- `go_reactome_target_prioritization.csv`: all 18 genes with aligned evidence and follow-up annotations.
- `go_reactome_cross_resource_target_shortlist.csv`: genes with same-theme GO–Reactome concordance.
- `go_reactome_structural_followup_prioritization.csv`: macrophage/network/structure follow-up branch.
