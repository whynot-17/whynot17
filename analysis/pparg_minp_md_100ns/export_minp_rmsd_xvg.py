"""Export GROMACS-style XVG RMSD traces for the completed MiNP 100 ns run."""

from __future__ import annotations

import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align


ROOT = Path(r"E:\chatgpt\whynot17\analysis\pparg_minp_md_100ns")
OUT = ROOT / "outputs"
TOPOLOGY = Path(r"E:\chatgpt\pparg_minp_md\system\seed20260917\initial.pdb")
TRAJECTORY = Path(r"E:\chatgpt\pparg_minp_md\run_20260915_seed20260917\production.dcd")


def apply_left_rotation(coords, rotation, mobile_center, reference_center):
    v = coords - mobile_center
    out = np.empty_like(v)
    out[:, 0] = v[:, 0] * rotation[0, 0] + v[:, 1] * rotation[0, 1] + v[:, 2] * rotation[0, 2] + reference_center[0]
    out[:, 1] = v[:, 0] * rotation[1, 0] + v[:, 1] * rotation[1, 1] + v[:, 2] * rotation[1, 2] + reference_center[1]
    out[:, 2] = v[:, 0] * rotation[2, 0] + v[:, 1] * rotation[2, 1] + v[:, 2] * rotation[2, 2] + reference_center[2]
    return out


def rmsd(coords, reference):
    delta = coords - reference
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def write_xvg(path, title, legend, rows):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"# MiNP–PPARG completed 100 ns MD; reference={TOPOLOGY}; fit=protein backbone\n")
        handle.write("# Heavy atoms only; RMSD values are in nm; time is in ns.\n")
        handle.write(f'@ title "{title}"\n')
        handle.write('@ xaxis label "Time (ns)"\n')
        handle.write('@ yaxis label "RMSD (nm)"\n')
        handle.write(f'@ s0 legend "{legend}"\n')
        for time_ns, value_A in rows:
            handle.write(f"{time_ns:.5f} {value_A / 10.0:.8f}\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    u = mda.Universe(str(TOPOLOGY), str(TRAJECTORY))
    ref = mda.Universe(str(TOPOLOGY))
    bb = u.select_atoms("protein and backbone")
    protein_heavy = u.select_atoms("protein and not name H*")
    ligand_heavy = u.select_atoms("resname UNK and not name H*")
    ref_bb = ref.select_atoms("protein and backbone")
    ref_protein_heavy = ref.select_atoms("protein and not name H*")
    ref_ligand_heavy = ref.select_atoms("resname UNK and not name H*")
    if len(bb) != 1052 or len(protein_heavy) != 2112 or len(ligand_heavy) != 21:
        raise RuntimeError(f"Unexpected selections: backbone={len(bb)}, protein={len(protein_heavy)}, ligand={len(ligand_heavy)}")

    ref_bb_xyz = ref_bb.positions.copy()
    ref_protein_xyz = ref_protein_heavy.positions.copy()
    ref_ligand_xyz = ref_ligand_heavy.positions.copy()
    ref_bb_center = ref_bb_xyz.mean(axis=0)
    ref_complex_xyz = np.vstack([ref_protein_xyz, ref_ligand_xyz])
    complex_rows = []
    protein_rows = []
    for frame_index, ts in enumerate(u.trajectory):
        mobile_bb = bb.positions.copy()
        mobile_center = mobile_bb.mean(axis=0)
        rotation, _ = align.rotation_matrix(mobile_bb - mobile_center, ref_bb_xyz - ref_bb_center)
        box = np.asarray(ts.dimensions[:3], dtype=float)
        ligand_xyz = ligand_heavy.positions.copy()
        if np.all(box > 0):
            ligand_xyz += np.rint((mobile_center - ligand_xyz.mean(axis=0)) / box) * box
        aligned_protein = apply_left_rotation(protein_heavy.positions.copy(), rotation, mobile_center, ref_bb_center)
        aligned_ligand = apply_left_rotation(ligand_xyz, rotation, mobile_center, ref_bb_center)
        aligned_complex = np.vstack([aligned_protein, aligned_ligand])
        time_ns = float(ts.time / 1000.0)
        complex_rows.append((time_ns, rmsd(aligned_complex, ref_complex_xyz)))
        protein_rows.append((time_ns, rmsd(aligned_protein, ref_protein_xyz)))

    complex_path = OUT / "MiNP_complex_rmsd_0_100ns.xvg"
    protein_path = OUT / "MiNP_protein_rmsd_0_100ns.xvg"
    write_xvg(complex_path, "MiNP–PPARG complex RMSD", "Complex heavy atoms", complex_rows)
    write_xvg(protein_path, "MiNP–PPARG protein RMSD", "Protein heavy atoms", protein_rows)
    manifest = {
        "trajectory": str(TRAJECTORY),
        "reference": str(TOPOLOGY),
        "fit_group": "protein backbone",
        "complex_group": "protein and ligand heavy atoms",
        "protein_group": "protein heavy atoms",
        "n_frames": len(complex_rows),
        "time_ns": [complex_rows[0][0], complex_rows[-1][0]],
        "outputs": [complex_path.name, protein_path.name],
    }
    (OUT / "rmsd_xvg_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
