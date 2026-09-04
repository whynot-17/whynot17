# PTGER4–DINP PPM Pose-Transfer QC

- Overall pose-transfer QC: `PASS`
- Selected DINP affinity: `-7.429 kcal/mol`
- Timestamp (UTC): `2026-09-04T03:18:11.573881+00:00`

## Provenance and protocol identity

- Source receptor: `/mnt/d/Codex/projectless-workspaces/2026-08-22/non/work/whynot17/analysis/dinp_crc_structural_pipeline/outputs/stage2_refine_and_rescue/ptger4_rescue/PTGER4_9JQZ_chainA_clean.pdb`
- PPM-oriented receptor: `/mnt/d/Codex/projectless-workspaces/2026-08-22/non/work/whynot17/analysis/dinp_crc_structural_pipeline/outputs/ppm_local/9JQZ_PTGER4_PPM_oriented.pdb`
- Source DINP pose: `/mnt/d/Codex/projectless-workspaces/2026-08-22/non/work/whynot17/analysis/dinp_crc_structural_pipeline/outputs/stage2_refine_and_rescue/ptger4_rescue/DINP_PTGER4_9JQZ_rescue_pose.pdbqt`
- Stage-2 audit: `/mnt/d/Codex/projectless-workspaces/2026-08-22/non/work/whynot17/analysis/dinp_crc_structural_pipeline/outputs/stage2_refine_and_rescue/audit.json`
- PTGER4 control audit: `/mnt/d/Codex/projectless-workspaces/2026-08-22/non/work/whynot17/analysis/dinp_crc_structural_pipeline/outputs/ptger4_control_docking/audit.json`
- Same prepared docking receptor: `PASS`
- Same pocket/protocol: `PASS`
- Selected Vina affinity from pose remark: `-7.429 kcal/mol`

## A. Receptor mapping and rigid fit

- Mapping key: `chain + residue number + insertion code + residue name + atom name`
- Common receptor atoms: `2209`
- Unmatched source atoms: `0`
- Unmatched PPM atoms: `0`
- Common Cα atoms / RMSD after fit: `281` / `0.000491874065 Å` (threshold `<0.01 Å`)
- Common backbone atoms / RMSD after fit: `1124` / `0.000488286545 Å` (threshold `<0.01 Å`)
- Common heavy atoms / RMSD after fit: `2209` / `0.000492461955 Å` (threshold `<0.01 Å`)
- Receptor pre-fit RMSD: `224.457 Å`
- Receptor post-fit all-atom RMSD: `0.000492461955 Å`
- Receptor transform QC: `PASS`

## B. Ligand integrity

- DINP atoms before/after: `30` / `30`
- DINP heavy atoms before/after: `30` / `30`
- Ligand internal pairwise-distance RMSD (pre-serialization rigid transform): `1.15700456e-14 Å` (threshold `<1e-5 Å`)
- Ligand internal geometry changed: `NO`
- Ligand processing: no redocking, minimization, embedding, conformer regeneration, coordinate reordering, or bond-order inference.
- PDB serialization note: source PDBQT uses repeated element-only atom names (`C`/`O`); output PDB uses unique names (`C01`, `O08`, ...) so BioPython/MDAnalysis retain every atom. The original names are preserved in the REMARK mapping and transform JSON.

## C. Receptor–ligand relative pose preservation

- Pocket definition: `24 receptor residues with any DINP heavy atom within 5.0 Å in the original pose`
- Minimum ligand–pocket heavy-atom distance before/after: `3.25154` / `3.25154 Å`
- Ligand COM → pocket COM distance before/after: `1.09714` / `1.09714 Å`
- Ligand–pocket pairwise distance-matrix RMSD: `1.43070801e-14 Å` (threshold `<1e-4 Å`)
- Ten nearest contact distances before: `[3.2515, 3.2725, 3.2739, 3.3068, 3.3758, 3.377, 3.3981, 3.4166, 3.4247, 3.4267]`
- Ten nearest contact distances after: `[3.2515, 3.2725, 3.2739, 3.3068, 3.3758, 3.377, 3.3981, 3.4166, 3.4247, 3.4267]`
- Ligand–pocket relative pose preserved: `PASS`

## D. PPM membrane sanity check

- Parsed PPM DUM boundary records: `495`
- Membrane center Z: `0.000 Å`
- Hydrophobic slab / boundary Z: `-16.200 to 16.200 Å`
- Half-thickness inferred from DUM records: `16.200 Å`
- Transformed DINP ligand COM Z: `12.769 Å`
- DINP near receptor pocket: `YES`
- Interpretation: geometric sanity check only; no membrane-binding energy or biological conclusion is inferred.

## E. Output re-read checks

- Usable protein-only complex BioPython reread: `PASS (2239)`
- Usable protein-only complex MDAnalysis reread: `PASS (2239 atoms)`
- Native-DUM complex BioPython reread: `PASS`
- Native-DUM complex MDAnalysis reread: `PASS`
- PPM native DUM information preserved: `YES` in `9JQZ_PTGER4_PPM_oriented.pdb` and `PTGER4_DINP_PPM_oriented_complex_native_dum.pdb`
- The primary complex omits only PPM DUM annotation records so it is readable by standard protein/ligand tools; receptor coordinates are unchanged.

## Boundary

This task performed coordinate transfer and QC only. No docking, OpenMM, membrane construction, POPC insertion, or molecular dynamics was run.
