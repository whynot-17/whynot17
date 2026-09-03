# DINP–CRC Stage-2 docking summary

DINP PubChem CID: 590836
RDKit conformers generated: 50
Ligand minimization: MMFF94

Docking is performed only when a co-crystal-like pocket is detected and both ligand and receptor PDBQT preparation succeed.
Targets without an auditable pocket are intentionally skipped rather than subjected to blind docking.

| Target | PDB | Pocket ligand | Docking | Best Vina affinity (kcal/mol) |
|---|---|---|---|---:|
| PTGER4 | 9JQZ | A1ECR | skipped |  |
| CXCR4 | 3ODU | ITD | ok | -5.499 |
| MMP9 | 6ESM | B9Z | ok | -3.995 |
| STAT3 | 6NJS | KQV | ok | -4.022 |

The MD shortlist is a computational prioritization output only. Manual inspection of chain identity, receptor completeness, pocket occupancy, protonation, cofactors and membrane context is required before production MD.
