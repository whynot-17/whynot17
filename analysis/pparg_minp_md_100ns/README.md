# PPARG–MiNP 100 ns MD analysis

This module stores the complete analysis of the 100 ns PPARG–MiNP production trajectory (`seed20260917`). The trajectory was analyzed after protein-backbone fitting with whole-ligand, aromatic-core, and two branch RMSD/COM metrics; pocket-distance and residency metrics; 5 ns contact summaries; residue contact occupancy; ligand and pocket RMSF; final-30-ns core clustering; and representative structures.

The raw production DCD remains on the local E: drive at `E:\chatgpt\pparg_minp_md\run_20260915_seed20260917\production.dcd` and is intentionally not copied into Git. The checked-in CSV/JSON/PDB files are the reproducible result tables and representative structures.

## Re-run

```powershell
& 'E:\chatgpt\pparg_md_env\python.exe' -u .\analyze_minp_complete.py
& 'E:\chatgpt\pparg_md_env\python.exe' -u .\plot_minp_complete.py
```

The report is in `MINP_PPARG_MD_COMPLETE_REPORT.md`. Figures are under `figures/`; numeric tables and structures are under `outputs/`.

## RG/SASA audit

`calculate_rg_sasa.py` samples the complete production DCD every 0.1 ns and writes `outputs/rg_sasa_0_100ns.csv`, `outputs/rg_sasa_window_summary.csv`, and `outputs/rg_sasa_0_100ns_summary.json`. Protein and ligand radius of gyration are mass-weighted over heavy atoms. SASA uses FreeSASA's Lee–Richards algorithm with a 1.4 Å probe; protein, ligand, and complex values are reported, together with isolated-vs-complex buried surface areas and SASA of the initial pocket atom set. The raw DCD remains on the E: drive.

```powershell
& 'E:\chatgpt\pparg_md_env\python.exe' -u .\calculate_rg_sasa.py
```

## XVG RMSD exports

`export_minp_rmsd_xvg.py` exports GROMACS-style XVG traces from the completed 100 ns trajectory. Both traces use the initial PDB as reference and fit each frame on the protein backbone; the complex trace contains protein plus ligand heavy atoms, and the protein trace contains protein heavy atoms. Values are in nm and time is in ns.

```powershell
& 'E:\\chatgpt\\pparg_md_env\\python.exe' -u .\\export_minp_rmsd_xvg.py
```

`export_minp_metrics_xvg.py` additionally exports `MiNP_rmsd_0_100ns.xvg`, `MiNP_ligand_rmsf_0_100ns.xvg`, `MiNP_pocket_ca_rmsf_0_100ns.xvg`, `MiNP_sasa_0_100ns.xvg`, `MiNP_rg_0_100ns.xvg`, and `MiNP_hbond_0_100ns.xvg`. RMSD/Rg/SASA/H-bond traces keep their native sampling; RMSF files are profiles over ligand atoms or pocket residues.

```powershell
& 'E:\\chatgpt\\pparg_md_env\\python.exe' -u .\\export_minp_metrics_xvg.py
```

`export_minp_gyrate_xvg.py` is the authoritative gyrate export. It writes total Rg plus X/Y/Z components for both protein and ligand; `MiNP_gyrate_0_100ns.xvg` therefore has nine numeric columns including time.

## Main result

The ligand undergoes immediate relaxation, then a brief excursion around 75.50 ns, and returns to a later shifted ensemble. The aromatic core remains internally rigid, total contacts remain persistent, the original six-atom contact retention is low because the contact geometry changes, and the PPARG pocket does not show structural collapse.
