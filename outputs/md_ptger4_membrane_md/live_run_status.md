# PTGER4–DINP MD run status

Started 2026-09-04 on the E: drive.

- Environment: `E:\chatgpt\md_env` (Python 3.11, OpenMM 8.6, MDAnalysis 2.10, PDBFixer 1.12, Open Babel wheel)
- Platform: OpenCL on NVIDIA GeForce RTX 4060 Laptop GPU, mixed precision
- System: 73,768 atoms; POPC membrane, 310 K, 0.15 M NaCl
- Target: 100 ns production, 2 fs timestep
- Source checkpoint: resumed from the latest completed checkpoint at 24.1 ns after the desktop window was closed
- Background process: PID 215624
- Output directory: `E:\chatgpt\ptger4_membrane_md_20260904`
- Live progress at this update: approximately 50.3 ns at 62.7 ns/day; about 19 hours remain to 100 ns

The runner writes a DCD frame and CSV row every 10 ps and a checkpoint every 100 ps. The 50 ns snapshot shows stable protein Cα RMSD, retained ligand contacts and no nonfinite coordinates. Stable trajectories indicate computational structural plausibility only and do not establish experimental binding or causal biology.
