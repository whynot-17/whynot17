# 4XLD native-ligand self-redocking QC

## Scope

This is a protocol QC run using the co-crystallized 4XLD ligand **BRL** (rosiglitazone), not DINP. The native BRL heavy-atom coordinates were extracted from 4XLD and combined with the RCSB Chemical Component Dictionary topology before ligand preparation.

## Protocol

- Receptor: `prepared/PPARG_4XLD_rigid.pdbqt`
- Native ligand: BRL, chain A, residue 502 from `inputs/4XLD.pdb`
- Box: `prepared/PPARG_4XLD_box.txt`
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
| top-ranked affinity | -8.816 kcal/mol |
| top-ranked receptor-frame heavy-atom RMSD | 7.158 Å |
| top-ranked Kabsch-fit heavy-atom RMSD | 1.674 Å |
| top-ranked pose QC | FAIL |
| lowest-RMSD mode | 5 |
| lowest receptor-frame heavy-atom RMSD | 1.270 Å |
| lowest-RMSD pose affinity | -8.260 kcal/mol |

The top-ranked pose is considered protocol-valid only if its receptor-frame heavy-atom RMSD is ≤2.0 Å. A low RMSD in a lower-ranked mode is reported as supportive but does not replace the top-ranked criterion.

## Files

- `4XLD_native_BRL.pdb`: extracted native ligand coordinates
- `4XLD_native_BRL.sdf`: native-coordinate ligand with CCD topology
- `4XLD_native_BRL.pdbqt`: prepared docking ligand
- `4XLD_native_redocked_BRL.pdbqt`: Vina poses
- `4XLD_native_redocking_rmsd.csv`: affinity and RMSD for each returned mode
- `4XLD_native_redocking.log`: complete Vina output
- `run_4XLD_native_redocking.py`: reproducible script

This QC does not perform molecular dynamics and does not validate CEBPB as a docking target. It only tests the PPARG/4XLD docking protocol used for the subsequent PPARG-DINP MD starting pose.
