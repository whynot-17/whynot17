# Academic Figure Skill QA Report — Figure 1 v2

## Contract

- Core conclusion: MiNP ranked 24th in the hypothesis-agnostic molecular screen but advanced through multistage actionability prioritization to a human-testable DINP exposure axis represented by urinary MCOP.
- Archetype: schematic-led horizontal composite with Panel B as the quantitative anchor.
- Target: Nature-family double-column figure, 183 mm wide.
- Export: RGB PDF vector master, 300 dpi PNG preview, editable-text SVG.

## Four-pass QA

- **PASS — Anti-pattern scan:** restrained semantic colors, no decorative chart types, no causal connector, no default rainbow palette.
- **PASS — Code compliance:** mandatory typography/color/export baselines retained; all 267 chemicals used; minimum text size is 5.0 pt.
- **PASS — Visual logic:** MiNP is labelled as rank 24; DINP parent is labelled rank 107 with BH-FDR 0.449; MCOP is explicitly a biomarker translation rather than a direct CTD hit.
- **PASS — Translation evidence:** Panel C reads the seven-cycle MiNP and MCOP above-LOD percentages directly from the frozen Phase 2A audit.
- **PASS — Anti-redundancy:** the roadmap/future-replication panel was removed; detailed Top-30 triage is reserved for a dedicated candidate-prioritization figure.
- **PASS — Statistical terminology:** MiNP is labelled `degree-matched FDR 0.036`; the raw degree-matched empirical P is not conflated with the BH-adjusted value.
- **PASS — Card containment:** rendered text bounding boxes were checked programmatically and all registered card text remains inside its frame.
- **PASS — Evidence boundary:** multistage prioritization is presented as selection logic, not as proof that MiNP was the top-ranked molecular hit.

## Verdict

**READY — Figure 1 passes the publication-oriented scientific and overflow QA checks.**
