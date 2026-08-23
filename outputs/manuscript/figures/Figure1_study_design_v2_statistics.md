# Figure 1 academic v2 — statistics and reproducibility

## Scientific question
How did the data-first environmental screen lead to DINP-axis biomonitoring and the subsequent CRC biological-state analyses?

## Quantitative panel
- Panel B uses all 267 frozen Phase 1 chemicals; no rows were downsampled or dropped.
- Primary ranking metric: `-log10(BH-FDR)` from the GeneCards Disorders-scoped, `gene_cards_k=1000`, `U_core` screen.
- Frozen screen checks: 52 chemicals with BH-FDR < 0.05; 69 stable candidates.
- Panel B is descriptive; it does not introduce a new statistical test or re-rank candidates.

## Schematic panels
- Panel A: workflow nodes from `figure1_panelA_workflow_nodes.csv`; no quantitative effect estimate is encoded.
- Panel C: translation boundary is explicitly labelled as biomarker translation; MCOP is not represented as a direct CTD nomination.
- Panel D: completed evidence layers are solid; the future prospective replication stage is dashed and contains no result. The specific cohort identity is reserved for the manuscript legend.

## Source traceability
- `outputs\manuscript\figures\source_data\figure1_primary_screen.csv` — all 267 rows — Panel B ranked screen.
- `outputs\manuscript\figures\source_data\figure1_panelA_workflow_nodes.csv` — 4 workflow nodes — Panel A workflow labels.
- `outputs\manuscript\figures\source_data\figure1_panelD_study_roadmap.csv` — 3 frozen stages — Panel D completed roadmap.
- `outputs/mcop_crc_phase2h_primary_reanalysis.csv` — frozen primary survey result — Panel D NHANES evidence label.

## Export
- Vector master: PDF.
- Preview: RGB PNG at 300 dpi.
- Editable text: SVG export.
