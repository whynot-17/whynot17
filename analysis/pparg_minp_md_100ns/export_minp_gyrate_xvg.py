"""Export GROMACS-style total and X/Y/Z radius-of-gyration components."""

from __future__ import annotations

import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np


ROOT = Path(r"E:\chatgpt\whynot17\analysis\pparg_minp_md_100ns")
OUT = ROOT / "outputs"
TOPOLOGY = Path(r"E:\chatgpt\pparg_minp_md\system\seed20260917\initial.pdb")
TRAJECTORY = Path(r"E:\chatgpt\pparg_minp_md\run_20260915_seed20260917\production.dcd")
STRIDE = 10  # 0.1 ns from the 10 ps reporter interval


def components(xyz, masses):
    center = np.average(xyz, axis=0, weights=masses)
    centered = xyz - center
    axis = np.sqrt(np.average(centered * centered, axis=0, weights=masses))
    total = float(np.sqrt(np.sum(axis * axis)))
    return total, float(axis[0]), float(axis[1]), float(axis[2])


def main():
    u = mda.Universe(str(TOPOLOGY), str(TRAJECTORY))
    protein = u.select_atoms("protein and not name H*")
    ligand = u.select_atoms("resname UNK and not name H*")
    protein_masses = np.asarray(protein.masses, dtype=float)
    ligand_masses = np.asarray(ligand.masses, dtype=float)
    rows = []
    for frame_index, ts in enumerate(u.trajectory):
        if frame_index % STRIDE:
            continue
        protein_xyz = protein.positions.copy()
        ligand_xyz = ligand.positions.copy()
        box = np.asarray(ts.dimensions[:3], dtype=float)
        if np.all(box > 0):
            ligand_xyz += np.rint((protein_xyz.mean(axis=0) - ligand_xyz.mean(axis=0)) / box) * box
        rows.append((float(ts.time / 1000.0), components(protein_xyz, protein_masses), components(ligand_xyz, ligand_masses)))

    path = OUT / "MiNP_gyrate_0_100ns.xvg"
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# MiNP–PPARG completed 100 ns MD; seed20260917; heavy atoms; mass-weighted Rg\n")
        handle.write("# Columns: time, protein Rg, protein RgX, protein RgY, protein RgZ, ligand Rg, ligand RgX, ligand RgY, ligand RgZ\n")
        handle.write('@ title "MiNP–PPARG radius of gyration (total and XYZ components)"\n')
        handle.write('@ xaxis label "Time (ns)"\n')
        handle.write('@ yaxis label "Rg (nm)"\n')
        legends = ["Protein Rg", "Protein RgX", "Protein RgY", "Protein RgZ", "Ligand Rg", "Ligand RgX", "Ligand RgY", "Ligand RgZ"]
        for i, legend in enumerate(legends):
            handle.write(f'@ s{i} legend "{legend}"\n')
        for time_ns, protein_vals, ligand_vals in rows:
            vals = protein_vals + ligand_vals
            handle.write(f"{time_ns:.5f} " + " ".join(f"{x / 10.0:.8f}" for x in vals) + "\n")

    manifest_path = OUT / "gyrate_xvg_manifest.json"
    manifest = {
        "source": str(TRAJECTORY),
        "sampling": "every 0.1 ns",
        "groups": {"protein": "protein heavy atoms", "ligand": "MiNP heavy atoms"},
        "columns": legends,
        "units": {"time": "ns", "rg": "nm"},
        "n_samples": len(rows),
        "output": path.name,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
