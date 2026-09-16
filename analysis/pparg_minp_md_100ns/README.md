# PPARG–MiNP 100 ns MD analysis

This module stores the complete analysis of the 100 ns PPARG–MiNP production trajectory (`seed20260917`). The trajectory was analyzed after protein-backbone fitting with whole-ligand, aromatic-core, and two branch RMSD/COM metrics; pocket-distance and residency metrics; 5 ns contact summaries; residue contact occupancy; ligand and pocket RMSF; final-30-ns core clustering; and representative structures.

The raw production DCD remains on the local E: drive at `E:\chatgpt\pparg_minp_md\run_20260915_seed20260917\production.dcd` and is intentionally not copied into Git. The checked-in CSV/JSON/PDB files are the reproducible result tables and representative structures.

## Re-run

```powershell
& 'E:\chatgpt\pparg_md_env\python.exe' -u .\analyze_minp_complete.py
& 'E:\chatgpt\pparg_md_env\python.exe' -u .\plot_minp_complete.py
```

The report is in `MINP_PPARG_MD_COMPLETE_REPORT.md`. Figures are under `figures/`; numeric tables and structures are under `outputs/`.

## Main result

The ligand undergoes immediate relaxation, then a brief excursion around 75.50 ns, and returns to a later shifted ensemble. The aromatic core remains internally rigid, total contacts remain persistent, the original six-atom contact retention is low because the contact geometry changes, and the PPARG pocket does not show structural collapse.
