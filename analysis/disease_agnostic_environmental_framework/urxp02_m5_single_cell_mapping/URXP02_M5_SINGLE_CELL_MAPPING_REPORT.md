# URXP02 M5 single-cell mapping

Generated 2026-08-30T11:52:55.225913+00:00. This is a focused donor-sex/cell-type context audit, not a causal mechanism analysis.

## Scope and design

- Six primary frozen M4 handoff genes were analysed with five prespecified secondary audit genes: JUN, TP53, AHR, CYP19A1, NFE2L2, NFKB1, ESR1, ESR2, AR, THRB, CASP3.
- Existing M3 module membership was also mapped using the union of 470 module genes (total expression-test gene universe: 470); no PPI was rerun.
- Public human CELLxGENE Discover h5ad datasets were selected before inspecting gene-level results: thyroid cells from Human Cell Landscape, adult kidney from the mature kidney atlas, and mesenteric arterial cells from an endothelium-enriched dataset.
- Curated `cell_type` labels were retained verbatim. Explicit keyword classes were used only to define the prespecified thyroid, endothelial, vascular smooth-muscle/pericyte, fibroblast, immune, and renal-epithelial contexts.
- Donor-level means were the inferential unit. Cells were not treated as independent biological replicates. Donor-level cell fractions are a separate composition diagnostic and are not population abundance estimates.
- Expression values use the deterministic transform recorded per row; no clustering, reannotation, disease screen, figure, or causal claim was added.
- Expression FDR family: all prespecified candidate-gene × target curated-cell-type × dataset tests (15510 planned rows; non-estimable rows have FDR=1). Composition diagnostics use a separate family (33 rows).

## Dataset and donor audit

- **thyroid / Construction of a human cell landscape at single-cell level**: 0 male and 12647 female selected cells; 0 male and 2 female donors; 456/470 analysis genes found.
- **kidney / Mature kidney dataset: full**: 13950 male and 26318 female selected cells; 6 male and 7 female donors; 469/470 analysis genes found.
- **vascular / scRNA-seq data analysis of endothelium-enriched mesenteric arterial tissues from human donors**: 7403 male and 3840 female selected cells; 3 male and 1 female donors; 468/470 analysis genes found.

## Results

- Expression contexts with FDR <0.05: **0**.
- Cell-composition contexts with FDR <0.05: **0**.
- M3 module activity contexts mapped: **4212** (2709 donor-supported tests); module-activity FDR hits: **0**.
- The smallest unadjusted expression P was **0.0007523** for **ANGPTL4** in **kidney / kidney interstitial fibroblast**; its fixed-family FDR was **1**.
- The largest absolute donor-level expression contrast was **2.437** for **FOS** in **kidney / podocyte**; it is retained as a descriptive estimate, not a validated sex-specific effect.
- The thyroid slice contained **0 male and 2 female donors**, so every thyroid sex contrast is non-estimable and carries FDR=1. The arterial slice contained only **3 male and 1 female donor**, so its null result is low-powered.
- Because the thyroid dataset has no male donors, the prespecified cross-tissue working pattern cannot be classified as supportive in this M5 resource set; the cross-tissue CSV records it as non-estimable rather than imputing a direction.
- These are localized donor-sex contrasts. They do not demonstrate that a cell-type contrast mediates the NHANES URXP02 phenotype.

## Interpretation guardrail

Any FDR hit is a context-specific candidate for follow-up. A lack of a hit is not evidence of no biology when donor numbers are small, especially in the four-donor arterial dataset. The analysis therefore distinguishes within-cell expression differences from donor-level captured-cell composition and does not combine them into a score.

## Files

The CSVs preserve the original curated labels, donor counts, effect estimates, raw P/FDR, expression transform, dataset accession, and explicit mapping rules. No figures were generated.
