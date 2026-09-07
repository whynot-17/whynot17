# PGE2–PTGER4 robustness re-analysis

This module stress-tests the focused `PGE2 synthesis source -> PTGER4 receiver`
hypothesis without changing the prior module, clustering, or macrophage
annotation. The primary unit is donor × cell-group pseudobulk. It compares:

- Score A: `PLA2G4A + PTGS2 + PTGES`;
- Score B: `PLA2G4A + PTGS1 + PTGS2 + PTGES + PTGES2` (no `PTGES3`);
- Score C: the original six-gene score, retained as sensitivity only;
- Score D: a bottleneck-aware geometric mean of upstream, cyclooxygenase, and terminal-synthase components.

Source ranks are evaluated at donor-group minimum cell thresholds of 10, 20,
and 30, with 20 as the primary threshold. Receiver metrics are PTGER4 mean,
detection, expressing-cell mean, pseudobulk PTGER4, and donor rank. The old
pooled cell-level z-score is retained only as a sensitivity column.

## Run

```powershell
E:\chatgpt\sc_env\Scripts\python.exe `
  analysis/dinp_crc_pge2_ptger4_axis_robustness/run_pge2_ptger4_robustness.py
```

The default input is the target-gene cell-score table from the preceding
module, with raw matrix and annotation paths recorded and hashed in
`outputs/input_manifest.json`. No CellChat, spatial, ML, docking, MD,
clustering, or new macrophage annotation is performed.

The final interpretation is in `outputs/PGE2_PTGER4_ROBUSTNESS_REPORT.md`.
All axis values are descriptive ranking sensitivities and do not establish
extracellular PGE2 communication.
