# Figure 1–5 V1 QA report

## Automated preflight

- Final plotting-source validator: 18 PASS, 2 reviewed WARN, 0 FAIL.
- Reviewed WARN — TIFF absent: the requested delivery formats were PNG, PDF and SVG. The assembled figures are supplied as editable vector PDFs; the 600 dpi PNGs are review/PPT files, not substitutes for vector submission artwork.
- Reviewed WARN — width not inferred by static parser: the source explicitly defines `FIGURE_WIDTH_IN = 190 / 25.4`. Independent PDF page-size inspection confirmed 190.00 mm width for all five assembled figures.
- PDF text audit: 24 PDFs inspected (five assembled figures and 19 standalone panels), minimum glyph 5.1 pt, zero text runs below the 5 pt floor.
- SVG audit: all five assembled SVGs contain live `<text>` nodes.

## Dimensions

| Figure | PDF dimensions | PNG dimensions | PNG resolution |
|---|---:|---:|---:|
| Figure 1 | 190.00 × 95.25 mm | 4,488 × 2,250 px | 600 dpi |
| Figure 2 | 190.00 × 135.89 mm | 4,488 × 3,210 px | 600 dpi |
| Figure 3 | 190.00 × 146.05 mm | 4,488 × 3,450 px | 600 dpi |
| Figure 4 | 190.00 × 149.86 mm | 4,488 × 3,540 px | 600 dpi |
| Figure 5 | 190.00 × 91.95 mm | 4,488 × 2,172 px | 600 dpi |

## Data-integrity checks

- Figure 1B uses all 267 chemicals from the frozen primary Phase 1 screen (`GeneCards_Disorders`, `gene_cards_k = 1000`, `U_core`), rather than a visually selected top subset. The source table contains 52 BH-FDR-significant chemicals and 69 stability-retained candidates. MiNP is rank 24 and is not depicted as the top hit.
- Figure 2 reproduces N = 9,936, CRC = 70 and the frozen R/Python estimates. Seven LOCO and seven cycle-specific records are asserted in code; all LOCO lower confidence limits exceed 1.
- Figure 3 reads all requested frozen sensitivity, co-exposure, spline and weighted-quartile records. Log-scale axes are guarded against nonpositive values.
- Figure 4 asserts 32 TCGA pairs, 36 Census epithelial pairs, compartment pair counts of 36/33/31/35 and six GSE144735 pairs. Every paired observation is displayed; positional jitter is deterministic and changes no data value.
- Figure 5 distinguishes observed evidence from the hypothetical bridge in both source tables and connector style.

## Statistical display definitions

| Figure/panels | Biological/statistical unit | Center | Interval/spread | Test |
|---|---|---|---|---|
| 1B | Chemical | Enrichment odds ratio/rank | BH-FDR; degree-matched FDR for MiNP | Fisher/hypergeometric enrichment with BH correction; degree-matched permutation |
| 2A–D | NHANES participant within complex survey | Odds ratio per MCOP doubling | 95% CI | Survey-weighted logistic regression; global cycle interaction Wald test in 2D |
| 3A–B | NHANES participant within complex survey | Odds ratio per MCOP doubling | 95% CI | Survey-weighted logistic regression |
| 3C | NHANES participant within complex survey | RCS odds ratio vs weighted median | Pointwise 95% CI | Overall and nonlinear survey Wald F tests |
| 3D | NHANES participant within complex survey | Quartile odds ratio vs Q1 | 95% CI | Survey-weighted trend test |
| 4A/B/D | Matched patient or donor | Raw pairs, group medians and median paired delta | Raw pair distribution | Two-sided Wilcoxon signed-rank test |
| 4C | Matched Census donor | Median paired delta | IQR, 1.5-IQR whiskers and raw deltas | Two-sided Wilcoxon signed-rank test |

No technical-replicate aggregates are presented. No multiplicity correction is applied to prespecified robustness or paired disease-state panels unless explicitly shown in the source data.

## Rendered visual audit

All 19 panels and five assembled figures were inspected after final rendering. Panel titles, labels, data marks, confidence intervals and annotations are fully visible at assembled scale. Standalone panel exports temporarily hide neighboring axes, preventing fragments of adjacent panels from entering cropped files. The teal exposure family, burgundy tumor family, gray normal/control family and muted orange inflammatory accent remain consistent across Figures 1–5. Direct labels and marker/fill differences make the interpretation non-dependent on color alone.

No microscopy, photographic, blot or gel content is present; therefore crop, contrast, gamma, scale-bar and image-stitching checks are not applicable.

