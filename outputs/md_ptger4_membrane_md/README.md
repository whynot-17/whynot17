# PTGER4–DINP membrane MD results

This directory records milestone QC for the PTGER4–DINP membrane production run used in the Meldonium/DINP–CRC structural pipeline.

## Run

- System: 73,768 atoms, POPC membrane, 310 K, 0.15 M NaCl
- Integrator: Langevin middle, 2 fs timestep
- Platform: OpenMM OpenCL, mixed precision, NVIDIA GeForce RTX 4060 Laptop GPU
- Target: 100 ns production
- Local output: `E:\chatgpt\ptger4_membrane_md_20260904`
- At the time of this commit the live run had reached approximately 50.3 ns at about 62.7 ns/day; production continues toward 100 ns.

## Milestone QC

| Milestone | Protein Cα RMSD vs minimized | DINP RMSD (PBC corrected) | DINP COM displacement | 4.5 Å ligand–protein contacts | Interpretation |
|---|---:|---:|---:|---:|---|
| 1.6 ns checkpoint | 1.32 Å | 3.91 Å | — | 112 (initial 115) | Early numerical stability; ligand rearrangement |
| 10.0 ns frame | 1.68 Å | 4.32 Å | 0.47 Å | 132 | Protein stable; ligand remains near the starting site |
| 25.0 ns frame | 1.86 Å | 4.00 Å | 0.56 Å | — | Stable protein and retained ligand location |
| 35.0 ns frame | 1.71 Å | 4.05 Å | 0.73 Å | 113 | No sustained contact loss or clear dissociation |
| 50.0 ns frame | 1.54 Å | 4.66 Å | 0.33 Å | 120 | Protein stable; contacts retained; ligand remains near the starting region |

For the trajectory through 50.23 ns, protein Cα RMSD averaged 1.63 Å (maximum 2.04 Å). PBC-corrected DINP RMSD averaged 4.18 Å (maximum 5.01 Å), while the ligand COM displacement averaged 0.73 Å (maximum 2.83 Å). The 4.5 Å contact count averaged 122 (range 94–146); the exact 50 ns frame has 120 contacts. Temperature after the initial 0.5 ns was 310.29 ± 1.10 K and density was 1.0263 ± 0.0018 g/mL.

These results support computational structural plausibility and show ligand conformational flexibility with retention near the starting binding region. They do not establish experimental affinity, residence time, or causal DINP biology. Longer production and, preferably, independent replicas are needed for those questions.

## Files

- `status_25ns.json`, `status_35ns.json`, `status_50ns.json`: milestone status records
- `qc_10ns_vs_minimized.json`, `qc_35ns_vs_minimized.json`, `qc_50ns_vs_minimized.json`: PBC-corrected trajectory summaries
- `qc_10ns_vs_minimized.csv`, `qc_35ns_vs_minimized.csv`, `qc_50ns_vs_minimized.csv`: per-frame RMSD and COM series
- `production_state_log_to_10ns.csv`, `production_state_log_to_35ns.csv`, `production_state_log_to_50ns.csv`: OpenMM thermodynamic logs
- `qc_snapshot_1p6ns.json`, `qc_10ns_pbc_corrected_qc.json`, `frame_50ns.json`, `thermo_50ns.json`, `contacts_50ns.json`: frame and thermodynamic QC records

Large binary trajectories (`.dcd`), checkpoints (`.chk`) and coordinate snapshots are retained on the E: drive and intentionally excluded from this repository commit.
