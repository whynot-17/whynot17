# Figure 1 academic v2 — statistics and reproducibility

## Scientific question
How did a rank-24 MiNP molecular signal advance to a human-testable DINP exposure axis?

## Quantitative panel
- Panel B uses all 267 frozen Phase 1 chemicals; no rows were downsampled or dropped.
- Primary ranking metric: `-log10(BH-FDR)` from the GeneCards Disorders-scoped, `gene_cards_k=1000`, `U_core` screen.
- Frozen screen checks: 52 chemicals with BH-FDR < 0.05; 69 stable candidates.
- Panel B is descriptive; it does not introduce a new statistical test or re-rank candidates. Molecular rank is explicitly separated from final actionability.

## Schematic panels
- Panel A: workflow nodes from `figure1_panelA_workflow_nodes.csv`; no quantitative effect estimate is encoded.
- Panel C summarizes the prespecified prioritization dimensions at a high level; the detailed Top-30 tournament is reserved for a dedicated candidate-prioritization figure.
- Panel C: translation boundary is explicitly labelled as biomarker translation; MCOP is not represented as a direct CTD nomination.
- Panel C detectability values are read from the Phase 2A NHANES biomarker audit: MiNP 27.4% and MCOP 98.4% above LOD across seven cycles.
- The previous roadmap/future-replication panel was removed as redundant; prospective replication belongs in the Discussion rather than this figure.

## Source traceability
- `outputs\manuscript\figures\source_data\figure1_primary_screen.csv` — all 267 rows — Panel B ranked screen.
- `outputs\manuscript\figures\source_data\figure1_panelA_workflow_nodes.csv` — 4 workflow nodes — Panel A workflow labels.
- `outputs\nhanes_dinp_phase2a_audit_summary.csv` — MiNP and MCOP rows — Panel C seven-cycle detectability.

## Export
- Vector master: PDF.
- Preview: RGB PNG at 300 dpi.
- Editable text: SVG export.
