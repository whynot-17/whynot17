# 4XLD focused-box native redocking QC (exhaustiveness 128)

## Purpose

This is a focused follow-up to the original 4XLD BRL self-redocking QC. It keeps the same PPARG receptor, BRL native-coordinate ligand, Vina seed, and scoring protocol, but uses a ligand-centered pocket box and 4× higher exhaustiveness. The aim is to test whether the native-like pose can move into the top-ranked/top-three poses.

## Box definition

The box is centered on the native BRL heavy-atom bounding-box midpoint. Each axis uses the native ligand span plus an 8 Å total margin, rounded upward to practical dimensions. This gives a focused box of `14` × `20` × `16` Å, compared with the original 24 × 24 × 24 Å box.

- center: (57.2365, 29.3710, 3.7130) Å
- native ligand span: (4.783, 11.176, 6.656) Å

## Protocol

- Receptor: `../prepared/PPARG_4XLD_rigid.pdbqt`
- Native ligand: BRL/rosiglitazone, chain A residue 502 from `../inputs/4XLD.pdb`
- Engine: AutoDock Vina 1.2.7
- Seed: `20260909`
- Exhaustiveness: **128** (original: 32)
- CPU: `8`
- Modes: `20`
- Energy range: `5` kcal/mol
- Pre-specified native-redocking pass threshold: receptor-frame heavy-atom RMSD ≤ `2.0` Å

## Result

| metric | value |
|---|---:|
| top-ranked affinity | -8.517 kcal/mol |
| top-ranked receptor-frame RMSD | 2.346 Å |
| top-ranked Kabsch-fit RMSD | 1.621 Å |
| native-like pose in top1 | no |
| native-like pose in top3 | no |
| lowest-RMSD mode | 15 |
| lowest receptor-frame RMSD | 1.358 Å |
| lowest-RMSD pose affinity | -7.720 kcal/mol |
| number of poses with RMSD ≤2 Å | 2 |

The focused/high-exhaustiveness run is considered successful for ranking if a pose with RMSD ≤2 Å reaches top1 or top3. The result is reported separately from the original 24 Å-box run and does not alter the original DINP docking results.

## Files

- `PPARG_4XLD_focused_native_box.txt`: focused box parameters
- `4XLD_native_redocked_BRL_focused_exh128.pdbqt`: focused/high-exhaustiveness poses
- `4XLD_native_redocking_focused_exh128_rmsd.csv`: affinity and RMSD for all modes
- `4XLD_native_redocking_focused_exh128.log`: complete Vina output
- `../run_4XLD_native_redocking_focused.py`: reproducible script
