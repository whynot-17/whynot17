# PPARG–DINP 100 ns MD (seed 20260909)

## Run

- Complex: PPARG (PDB 4XLD) with DINP
- Engine: OpenMM, OpenCL platform, mixed precision
- Temperature: 310 K
- Production length: 100 ns (50,000,000 steps; 2 fs timestep)
- Ensemble settings: Amber14-all/TIP3P, 0.15 M ionic strength, 1 nm padding, MonteCarloBarostat at 1 bar
- Ligand parameterization: OpenFF 2.2.1 with RDKit Gasteiger charges
- The production run was resumed from the preserved checkpoint at approximately 49.16 ns and completed normally.

## Versioned artifacts

- `final.pdb`: final 100 ns snapshot
- `md_completion.json`: completion status and final step
- `md_setup_manifest.json`, `resume_manifest.json`, `system_build_manifest.json`: run and system provenance
- `production_state.csv`: compact 10 ps progress/energy trace
- `qc_0_100ns/`: whole-trajectory structural QC
- `qc_50_80ns/`: stable middle/late segment QC requested during the run
- `qc_80_100ns/`: terminal-segment QC
- `qc_components_0_100ns/`, `qc_components_50_80ns/`, and `qc_components_80_100ns/`: whole-ligand, aromatic-core, and branch-resolved QC

The raw `production.dcd` trajectory is approximately 6 GB and remains at `E:\chatgpt\pparg_dinp_md_run_20260909\seed20260909\production.dcd`; it is intentionally excluded from GitHub.

## QC interpretation

The protein backbone remains in a stable range over the run (0–100 ns mean RMSD 2.58 Å; 80–100 ns mean 2.86 Å). The 50–80 ns segment is comparatively well behaved: DINP heavy-atom RMSD 2.19 Å, COM displacement 1.09 Å, pocket contact retention 71.8%, and nearest pocket distance 3.06 Å.

The terminal 80–100 ns segment shows increased ligand reorientation/partial loosening (DINP RMSD 3.36 Å, COM displacement 1.58 Å, contact retention 66.3%; 95–100 ns contact retention falls to about 60%). The nearest pocket distance remains close to 3.10 Å, so this QC supports continued pocket association but does not establish experimental affinity or causality. Flexible-chain motion may contribute to the ligand RMSD increase.

## Core and branch-resolved QC

DINP heavy atoms were partitioned by the verified ligand atom order: the six phthalate aromatic-ring atoms form the core; each ester group plus its isononyl chain forms branch A or branch B. The A/B labels identify the two atom-order branches and do not imply different chemical identities.

The 50–80 ns reference segment has core RMSD 1.50 Å, branch A RMSD 2.29 Å, branch B RMSD 2.29 Å, core COM displacement 1.40 Å, and core contact occupancy 62.0%. In the 80–100 ns terminal segment, core RMSD remains lower than the whole-ligand value at 1.83 Å, while branch A rises to 4.49 Å and branch B remains 2.44 Å. Core COM displacement is 1.74 Å and core contact occupancy is 56.6%; during 95–100 ns, these are 2.08 Å, 4.81 Å, 2.67 Å, 2.00 Å, and 45.3%, respectively.

This identifies branch A rearrangement as the primary source of the late whole-ligand RMSD increase. The concurrent reduction in core-contact occupancy warrants a conservative interpretation: DINP remains pocket-associated by the nearest-distance and whole-pocket-contact metrics, but the terminal segment is not as stable as 50–80 ns.

All RMSD/COM/contact measurements are structural stability checks relative to the initial production frame and use PBC-aware ligand imaging.
