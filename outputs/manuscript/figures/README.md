# Figure 1–5 V1 production package

This directory contains the five assembled manuscript figures defined by the frozen Figure 1–5 Master Blueprint v1, together with individual panels, plotting source tables, exact Python scripts and QA documentation.

## Assembled figures

- `Figure1_study_design_v1.{png,pdf,svg}`
- `Figure2_nhanes_primary_v1.{png,pdf,svg}`
- `Figure3_robustness_v1.{png,pdf,svg}`
- `Figure4_ppar_singlecell_v1.{png,pdf,svg}`
- `Figure5_integrated_model_v1.{png,pdf,svg}`

The `panels/` directory contains 19 separately exported panels in the same three formats. The `source_data/` directory contains the exact plotting tables used by each quantitative panel.

## Rebuild

From the repository root:

```powershell
python work/scripts/make_manuscript_figures_v1.py `
  --repo-root . `
  --source-dir outputs/manuscript/figures/source_data `
  --output-dir outputs/manuscript/figures
```

`make_manuscript_figures_v1.py` is copied into this delivery directory for auditability. `derive_manuscript_rcs_ci_v1.py` reproduces the survey-weighted spline curve and confidence interval from the frozen analysis root.

## Production specification

- Final assembled width: 190 mm (Elsevier double-column target).
- PNG: 600 dpi RGB review files, 4,488 pixels wide.
- PDF: single-page editable vector artwork with embedded TrueType text.
- SVG: editable text nodes retained.
- Font stack: Arial, Helvetica, DejaVu Sans, sans-serif.
- Minimum audited PDF glyph: 5.1 pt.
- Captions and figure-level titles are kept outside the artwork in `figure_legends_v1.md`.

For submission, use the vector PDF as the primary artwork file. PNG is retained for rapid review and slide/PPT import; SVG is the editable master.

