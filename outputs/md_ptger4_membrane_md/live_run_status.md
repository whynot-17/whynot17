# PTGER4–DINP MD run status

Started 2026-09-04 on the E: drive.

- Environment: `E:\chatgpt\md_env` (Python 3.11, OpenMM 8.6, MDAnalysis 2.10, PDBFixer 1.12, Open Babel wheel)
- Platform: OpenCL on NVIDIA GeForce RTX 4060 Laptop GPU, mixed precision
- System: 73,768 atoms; POPC membrane, 310 K, 0.15 M NaCl
- Target: 100 ns production, 2 fs timestep
- Source checkpoint: resumed from the latest completed checkpoint at 24.1 ns after the desktop window was closed
- Background process: PID 215624 (resume with `run_ptger4_existing_system_md.py --production-ns 100 --resume` if interrupted)
- Output directory: `E:\chatgpt\ptger4_membrane_md_20260904`

The runner writes a DCD frame and CSV row every 10 ps and a checkpoint every 100 ps. After restart the speed readout is transient and is settling back toward approximately 63 ns/day; the estimated wall time for the remaining production is about 36 hours. Stable trajectories indicate computational structural plausibility only and do not establish experimental binding or causal biology.
