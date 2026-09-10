# 8CH8 native-ligand self-redocking QC

## Scope

This QC uses the co-crystallized 8CH8 ligand **ULC** (the PXR/NR1I2 ligand in chain A, residue 501), not DINP. Native ULC heavy-atom coordinates were extracted from `inputs/8CH8.pdb` and combined with the RCSB Chemical Component Dictionary topology before Meeko preparation.

## Protocol

- Receptor: `prepared/NR1I2_8CH8_rigid.pdbqt`
- Native ligand: ULC, chain A, residue 501 from `inputs/8CH8.pdb`
- Box: `prepared/NR1I2_8CH8_box.txt` (24 × 24 × 24 Å)
- Engine: AutoDock Vina 1.2.7
- Seed: `20260909`
- Exhaustiveness: `32`
- CPU: `8`
- Number of modes: `20`
- Energy range: `5` kcal/mol
- Pre-specified receptor-frame heavy-atom RMSD pass threshold: **≤ 2.0 Å**

## Result

| metric | value |
|---|---:|
| top-ranked affinity | -7.904 kcal/mol |
| top-ranked receptor-frame heavy-atom RMSD | 8.736 Å |
| top-ranked Kabsch-fit heavy-atom RMSD | 1.847 Å |
| top-ranked pose QC | FAIL |
| lowest-RMSD mode | 5 |
| lowest receptor-frame heavy-atom RMSD | 7.010 Å |
| lowest-RMSD pose affinity | -7.578 kcal/mol |
| number of poses with RMSD ≤2 Å | 0 |

The top-ranked pose is considered protocol-valid only if its receptor-frame heavy-atom RMSD is ≤2.0 Å. A low RMSD in a lower-ranked mode is reported as supportive but does not replace the top-ranked criterion.

## Files

- `8CH8_native_ULC.pdb`: extracted native ligand coordinates
- `8CH8_native_ULC.sdf`: native-coordinate ligand with CCD topology
- `8CH8_native_ULC.pdbqt`: prepared docking ligand
- `8CH8_native_redocked_ULC.pdbqt`: Vina poses
- `8CH8_native_redocking_rmsd.csv`: affinity and RMSD for each returned mode
- `8CH8_native_redocking.log`: complete Vina output
- `../run_8CH8_native_redocking.py`: reproducible script

This QC does not perform molecular dynamics and does not validate CEBPB as a docking target. It only tests the PXR/8CH8 docking protocol used for PXR-DINP structural interpretation.
