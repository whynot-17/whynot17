# DINP–PPARG vs MiNP–PPARG same-method 100 ns audit

Both systems were reanalyzed from their production DCD with the same code and definitions. Protein-backbone fitting uses the corresponding initial system PDB; ligand RMSD/COM values are relative to that initial pose. Contact cutoffs are 4.5 Å with PBC-aware distances and contact summaries are sampled every 0.1 ns and binned every 5 ns.

MiNP component definitions are disjoint and chemically explicit: core indices 12–17, long ester/alkyl branch A indices 0–11, and short carboxyl branch B indices 18–20. This replaces the earlier overlapping branch-B bookkeeping used in the preliminary MiNP QC.

The earlier DINP report used the first production frame as its RMSD reference. This audit deliberately uses the corresponding initial system PDB for both ligands, so the values below are a docking-to-trajectory comparison on one baseline. The per-frame tables preserve the exact definitions and make the baseline explicit.

## 80–100 ns comparison

| metric | DINP | MiNP |
|---|---:|---:|
| whole RMSD (Å) | 5.735 ± 0.295 | 7.614 ± 0.328 |
| core RMSD (Å) | 5.350 ± 0.727 | 5.000 ± 0.753 |
| branch A RMSD (Å) | 6.029 ± 0.455 | 8.593 ± 0.316 |
| branch B RMSD (Å) | 5.578 ± 0.420 | 7.697 ± 0.614 |
| whole COM (Å) | 1.964 ± 0.616 | 5.687 ± 0.426 |
| core COM (Å) | 5.285 ± 0.721 | 4.609 ± 0.799 |
| pocket min distance (Å) | 3.504 ± 0.245 | 3.800 ± 0.238 |
| initial core-contact retention | 0.288 ± 0.126 | 0.261 ± 0.125 |
| protein backbone RMSD (Å) | 2.898 ± 0.121 | 2.629 ± 0.121 |
| pocket Cα RMSD (Å) | 1.358 ± 0.249 | 0.907 ± 0.172 |

## Contact summary

| window | DINP total residues | MiNP total residues | DINP hydrophobic | MiNP hydrophobic | DINP H-bond proxy | MiNP H-bond proxy |
|---|---:|---:|---:|---:|---:|---:|
| 0–10 ns | 20.48 | 18.74 | 12.43 | 10.99 | 0.71 | 0.29 |
| 10–50 ns | 20.20 | 17.80 | 12.42 | 10.84 | 0.70 | 0.09 |
| 50–80 ns | 21.34 | 17.89 | 12.76 | 12.56 | 0.76 | 0.06 |
| 80–100 ns | 19.45 | 18.01 | 13.20 | 12.15 | 0.59 | 0.05 |

## Interpretation

- DINP initial core-pocket retention: 0.998 residency; farthest sampled core distance 4.99 Å at 14.33 ns.
- MiNP initial core-pocket retention: 0.957 residency; farthest sampled core distance 7.14 Å at 72.77 ns.
- DINP last-30-ns core clustering: 4 clusters, dominant occupancy 0.917.
- MiNP last-30-ns core clustering: 7 clusters, dominant occupancy 0.497.

The detailed audit tables are in `outputs/`. Geometric hydrogen bonds are a distance/angle proxy and do not replace an energy decomposition. The comparison evaluates structural persistence and rearrangement, not experimental affinity.
