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

The raw `production.dcd` trajectory is approximately 6 GB and remains at `E:\chatgpt\pparg_dinp_md_run_20260909\seed20260909\production.dcd`; it is intentionally excluded from GitHub.

## QC interpretation

The protein backbone remains in a stable range over the run (0–100 ns mean RMSD 2.58 Å; 80–100 ns mean 2.86 Å). The 50–80 ns segment is comparatively well behaved: DINP heavy-atom RMSD 2.19 Å, COM displacement 1.09 Å, pocket contact retention 71.8%, and nearest pocket distance 3.06 Å.

The terminal 80–100 ns segment shows increased ligand reorientation/partial loosening (DINP RMSD 3.36 Å, COM displacement 1.58 Å, contact retention 66.3%; 95–100 ns contact retention falls to about 60%). The nearest pocket distance remains close to 3.10 Å, so this QC supports continued pocket association but does not establish experimental affinity or causality. Flexible-chain motion may contribute to the ligand RMSD increase.

All RMSD/COM/contact measurements are structural stability checks relative to the initial production frame and use PBC-aware ligand imaging.
