# Academic Figure Skill QA — Figures 3–5 redesign v4

Target: Nature-family double-column, 183 mm, RGB. Backend: Python/matplotlib.

## Figure contract and anti-redundancy

- Figure 3: pooled association, LOCO stability, and cycle heterogeneity answer three distinct questions.
- Figure 4: spline shape, prespecified stress-test fingerprint, and pairwise co-exposure specificity are non-redundant.
- Figure 5: cross-platform paired replication, within-donor module coupling, and compartment localization are non-redundant.
- No conventional forest plot or repeated spaghetti panel is used.

## Code/data checks

- PASS: mandatory typography, palette, and export baselines included.
- PASS: 183-mm dimensions; PDF/SVG vector masters and 300-dpi PNG previews.
- PASS: exact estimates and 95% CIs imported from frozen source data.
- PASS: no downsampling; every paired donor/patient is plotted in Figure 5A/B.
- PASS: nulls, confidence encoding, sample sizes, and exact P values are defined on-figure or in the statistics report.
- PASS: color is redundant with shape, arrows, labels, and direction.

## Visual review checklist

- Check panel-label alignment, title clearance, tick-label legibility, data occlusion, and cropped text in the rendered PNGs.
- Verify Figure 3 radial labels remain readable at manuscript scale.
- Verify Figure 4 ridgelines remain visually separated and the right-side OR labels fit inside the panel.
- Verify Figure 5 small-n GSE144735 points remain individually visible.
