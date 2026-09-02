# DINP–CRC multi-database target convergence

This analysis implements the frozen 3+3 source design:

- Exposure side: CTD, EPA CompTox/ToxCast, ChEMBL.
- CRC side: the ordinary GeneCards CRC top-2000 export already archived in the
  Phase 1 data directory, DisGeNET, and Open Targets.

The analysis keeps source-specific evidence separate.  A source that is not
accessible is recorded as `unavailable` or `conditional`; it is never treated
as biological negative evidence, and a related phthalate is never substituted
for parent DINP.

Run from the repository root with:

```powershell
C:\Users\21634\anaconda3\python.exe analysis/dinp_crc_multi_database_target_convergence/run_dinp_crc_multi_database_target_convergence.py
```

Outputs are written to `outputs/`.  The raw CTD interaction snapshot is not
duplicated; its original path and SHA-256 are recorded in
`outputs/source_manifest.json`.
