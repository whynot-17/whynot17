# Focused PGE2 -> PTGER4 axis

This module tests one narrow molecular hypothesis in the local GSE144735 CRC
single-cell cache:

`PGE2 synthesis-high source group -> PTGER4-high macrophage receiver`

It calculates a six-gene PGE2 synthesis score (`PLA2G4A`, `PTGS1`, `PTGS2`,
`PTGES`, `PTGES2`, `PTGES3`) for every annotated cell, then compares
`PTGER4` mean expression, detection fraction, and PTGER4-high-cell fraction
across the four frozen macrophage states from the preceding module. Tumor
donor-level aggregates retain the six-donor structure and expose sparse
donor/state combinations.

## Run

```powershell
E:\chatgpt\sc_env\Scripts\python.exe `
  analysis/dinp_crc_pge2_ptger4_axis/run_pge2_ptger4_axis.py
```

The default inputs are the cached matrix and annotation under `E:\mcop` and
the prior macrophage subtype output in this repository. Override them with
`--matrix`, `--annotation`, `--subtype`, or `--out` when needed.

## Outputs

- `outputs/PGE2_PTGER4_REPORT.md`: concise interpretation with pooled and donor-level checks.
- `outputs/pge2_synthesis_by_group.csv`: all-cell group ranking and six-gene components.
- `outputs/ptger4_macrophage_subtype_statistics.csv`: PTGER4 mean, detection and high-cell fractions.
- `outputs/pge2_ptger4_axis_donor_summary.csv`: source -> C1QC-like donor consistency.
- `outputs/pge2_ptger4_top_calls_by_donor.csv`: per-donor top source and receiver calls.
- `outputs/figures/`: PNG, PDF and SVG figures.

The score is a descriptive expression summary. It does not establish
extracellular PGE2 production, receptor binding, spatial proximity, or true
cell-cell communication, so conclusions should use “supports” or “suggests”.
No ordinary CellChat run is required for this module.
