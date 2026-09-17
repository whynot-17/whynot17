# DINP–PPARG vs MiNP–PPARG same-method 100 ns audit

This module reanalyzes the completed DINP–PPARG and MiNP–PPARG production DCDs with one code path and one set of definitions. Both trajectories are protein-backbone fitted to their corresponding initial system PDB. Whole-ligand, core, branch RMSD/COM, core-to-pocket minimum distance, initial core-contact retention, protein/pocket RMSD, contact chemistry, pocket residency, component RMSF, and last-30-ns core clustering are reported.

The contact cutoff is 4.5 Å with periodic minimum-image distances. Contact chemistry is sampled every 0.1 ns and summarized in 5 ns bins. The audit uses the initial system PDB as the common docking-to-trajectory reference. Earlier DINP QC used the first production frame as its reference; those values should not be compared directly with this table.

MiNP components are disjoint: aromatic core heavy-atom indices 12–17, long ester/alkyl branch A indices 0–11, and short carboxyl branch B indices 18–20. DINP uses core indices 0–5 and two 12-atom ester/alkyl branches (6–17 and 18–29).

Run from the repository module with the E: drive environment:

```powershell
& 'E:\chatgpt\pparg_md_env\python.exe' -u .\run_same_method_audit.py
```

The combined tables are `outputs/comparison_audit_table.csv`, `outputs/contact_comparison_5ns.csv`, and `outputs/contact_window_summary.csv`. Per-system trajectories, residue occupancies, RMSF, clustering, and JSON summaries are under `outputs/DINP/` and `outputs/MiNP/`. The interpretation is structural plausibility only; it does not establish experimental affinity or causality.
