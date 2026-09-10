# PPARG/4XLD BRL structural overlap and core RMSD audit

## Scope

This audit compares the 4XLD crystal BRL ligand with the focused-box/exhaustiveness-128 PPARG native-redocking poses **Top1, Top2, and Mode15**. It addresses whether the Top1 receptor-frame RMSD of 2.346 Å reflects a local flexible-tail displacement or a whole-ligand flip/rearrangement.

## Coordinate treatments

Two coordinate views are retained:

1. **Receptor-frame view:** crystal and docked coordinates are compared without moving the docked pose. This preserves the ligand's position/orientation relative to the fixed PPARG receptor and is the primary native-redocking QC view.
2. **Full-heavy-atom aligned view:** each docked pose is superposed onto the crystal using a least-squares Kabsch fit over all 25 BRL heavy atoms. Core RMSD is then calculated on the transformed coordinates using the fixed core below. This separates ligand-shape/core preservation from receptor-frame placement.

## Pre-declared BRL core definition

The core is fixed before inspecting the pose-specific values and contains **19 heavy atoms**:

- thiazolidinedione pharmacophore: `S1, C2, N3, C4, C5, O2, O4`;
- phenyl ring: `C7, C8, C9, C10, C11, C12`;
- 2-pyridyl ring: `N18, C17, C19, C20, C21, C22`.

The excluded **6 heavy atoms** are the linker/side-chain atoms: `C6, O13, C14, C15, N16, C16`. They are not removed selectively per pose; the same fixed atom list is used for Crystal, Top1, Top2, and Mode15.

## Results

See `4XLD_BRL_full_and_core_RMSD.csv` for the machine-readable table. `full_aligned_rmsd_A` is the all-heavy-atom fit RMSD. `core_aligned_rmsd_A` is the core RMSD after that same all-heavy-atom fit. `excluded_linker_sidechain_aligned_rmsd_A` is provided only as a diagnostic for localization of deviations.

The PNG uses common axis limits across all panels. Gray is crystal BRL; the colored trace is the docked pose after full-heavy-atom alignment. Filled colored atoms are the pre-defined core, while open colored atoms are the excluded linker/side-chain atoms.

## Structural interpretation rule

- Low full-aligned and core RMSD with higher excluded-atom RMSD supports a preserved core with a localized flexible-tail/linker deviation.
- Core RMSD above 2 Å after full-heavy-atom fitting is treated as evidence against a preserved rigid core and is compatible with core rearrangement or an alternative orientation.
- Receptor-frame RMSD remains separate because full-atom alignment can make a displaced pose look geometrically similar while hiding its position relative to PPARG.

## Outputs

- `4XLD_BRL_crystal_top1_top2_mode15_receptor_frame_overlay.pdb`: raw receptor-frame overlap for molecular viewers
- `4XLD_BRL_crystal_top1_top2_mode15_full_heavy_aligned_overlay.pdb`: all-heavy-atom aligned overlap for shape/core inspection
- `4XLD_BRL_crystal_top1_top2_mode15_full_alignment_overlay.png`: common-scale XY/XZ visual overlap
- `4XLD_BRL_full_and_core_RMSD.csv`: full/core RMSD and interpretation table
- `../analyze_4XLD_BRL_overlap_core.py`: reproducible analysis script

This analysis does not change the docking ranking, does not perform MD, and does not redefine the previously specified ≤2 Å receptor-frame QC threshold.
