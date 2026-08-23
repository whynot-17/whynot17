# Academic Figure Skill QA Report — Figure 1 v2

## Contract

- Core conclusion: a data-first environmental screen prioritized a DINP-related exposure axis, translated it to urinary MCOP for NHANES evaluation, and connected the human association to CRC transcriptomic localization without claiming a causal DINP-to-CRC mechanism.
- Archetype: schematic-led horizontal composite with Panel B as the quantitative anchor.
- Target: Nature/Cell/Science-style double-column figure, 183 mm wide.
- Export: RGB PDF vector master, 300 dpi PNG preview, editable-text SVG.

## Pass 0 — Anti-pattern scan

- **PASS** AP-0: typography, color, and export baselines are present in the production script.
- **PASS** AP-1/AP-2: no default qualitative palette, jet, rainbow, or HSV colormap.
- **PASS** AP-3: plot spines are restricted to the quantitative Panel B axes; schematic panels are axis-free.
- **PASS** AP-4: no legend or annotation occludes data.
- **PASS** AP-5/AP-6: no 3D rendering, pie chart, or decorative chart type.
- **PASS** AP-7: causal language is explicitly avoided; the MCOP translation boundary and future WHI stage are labelled.

## Pass 1 — Code compliance

- **PASS** CL-1: smallest visible text is approximately 5.15 pt; primary labels are larger.
- **PASS** CL-2: figure is created at 183 mm width and 88 mm height.
- **PASS** CL-3: PDF TrueType embedding and editable SVG text are configured.
- **PASS** CL-4: all source paths are explicit and traceable.
- **PASS** CL-5: Panel B uses all 267 source rows; no downsampling or silent row removal.
- **PASS** CL-6: source-data validation checks 267 chemicals, 52 BH-FDR-significant chemicals, and 69 stable candidates.
- **PASS** CL-7: statistics/reproducibility report is generated beside the figure.

## Pass 2 — Visual logic and data integrity

- **PASS** VI-1: Panel B carries the strongest visual weight and the MiNP/DINP/MBzP anchors are directly labelled.
- **PASS** VI-2: Panel A establishes the unbiased screen; Panel C prevents MiNP/DINP–MCOP molecule conflation; Panel D separates completed evidence from future replication.
- **PASS** VI-3: blue exposure, red CRC state, grey background, and dashed future bridge have consistent semantic roles.
- **PASS** VI-4: no panel duplicates the same quantitative question.
- **PASS** VI-5: the 267-point ranked screen is finite and non-empty; the BH-FDR threshold is visible.
- **PASS** VI-6: text was shortened and reflowed after the first render to prevent title, footer, and roadmap overflow.

## Pass 3 — Rendered visual verification

- **PASS** VV-1: no visible data, label, or card occlusion in the 300 dpi PNG preview.
- **PASS** VV-2: panel edges and spacing are aligned; the horizontal narrative reads A → B → C → D.
- **PASS** VV-3: labels remain legible at the rendered double-column size.
- **PASS** VV-4: colors are distinguishable and remain interpretable without relying on a red–green pair.
- **PASS** VV-5: PNG dimensions are 2208 × 1045 pixels (300 dpi equivalent for the requested canvas); non-background pixel fraction is 18.3%, indicating visible content in all panels without saturation.

## Verdict

**READY — Figure 1 v2 passes the publication-oriented QA checks.**

The figure is a study-design/logic figure, not a causal model. The future stage is shown only as prospective replication, with the cohort identity reserved for the manuscript legend; the DINP-to-epithelial PPAR/NR relationship remains visually subordinate and unproven.
