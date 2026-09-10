from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
}
WATER_ION_RESNAMES = {"HOH", "WAT", "TIP3", "SOL", "NA", "CL", "K", "CA", "MG"}


def fit_kabsch(mobile: np.ndarray, reference: np.ndarray):
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    p = mobile - mobile_center
    q = reference - reference_center
    covariance = p.T @ q
    u, _, vt = np.linalg.svd(covariance)
    d = np.eye(3)
    d[2, 2] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ d @ vt
    aligned = p @ rotation + reference_center
    return aligned, mobile_center, rotation


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def main():
    parser = argparse.ArgumentParser(description="Analyze a PPARG-DINP trajectory segment without stopping the MD run.")
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--start-ns", type=float, default=0.0)
    parser.add_argument("--end-ns", type=float, default=10.0)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    universe = mda.Universe(str(args.topology), str(args.trajectory))
    protein = universe.select_atoms("protein")
    backbone = universe.select_atoms("protein and backbone")
    if protein.n_atoms == 0 or backbone.n_atoms == 0:
        raise RuntimeError("Protein or backbone selection is empty")
    candidates = [
        residue for residue in universe.residues
        if residue.resname not in PROTEIN_RESNAMES | WATER_ION_RESNAMES
        and len(residue.atoms) >= 20
    ]
    if not candidates:
        raise RuntimeError("Could not identify the DINP residue")
    ligand_residue = max(candidates, key=lambda residue: len(residue.atoms))
    ligand = ligand_residue.atoms
    ligand_heavy = ligand.select_atoms("not name H*")
    protein_heavy = protein.select_atoms("not name H*")
    if ligand_heavy.n_atoms < 20:
        raise RuntimeError(f"DINP heavy-atom selection unexpectedly small: {ligand_heavy.n_atoms}")

    universe.trajectory[0]
    reference_backbone = backbone.positions.copy()
    reference_ligand = ligand_heavy.positions.copy()
    initial_distances = distance_array(ligand_heavy.positions, protein_heavy.positions)
    initial_min_by_atom = initial_distances.min(axis=0)
    pocket_resids = sorted(set(int(resid) for resid in protein_heavy.resids[initial_min_by_atom <= 4.5]))
    if not pocket_resids:
        raise RuntimeError("No initial protein pocket atoms within 4.5 A of DINP")
    pocket_mask = np.isin(protein_heavy.resids, pocket_resids)
    pocket_atoms = protein_heavy[pocket_mask]
    pocket_resids_array = np.asarray(pocket_resids, dtype=int)

    rows = []
    end_ps = args.end_ns * 1000.0
    start_ps = args.start_ns * 1000.0
    if start_ps < 0 or end_ps <= start_ps:
        raise ValueError("Require 0 <= start-ns < end-ns")
    for frame_index, ts in enumerate(universe.trajectory):
        if ts.time > end_ps:
            break
        if ts.time < start_ps:
            continue
        # DCDReporter may wrap the protein and ligand into different periodic
        # images.  Translate DINP to the nearest image of the protein before
        # fitting or measuring contacts, otherwise a harmless box crossing
        # appears as an artificial ~80 A jump.
        box = np.asarray(ts.dimensions[:3], dtype=float)
        ligand_positions = ligand_heavy.positions.copy()
        if np.all(box > 0):
            shift = np.rint((ligand_positions.mean(axis=0) - backbone.positions.mean(axis=0)) / box)
            ligand_positions -= shift * box
        aligned_backbone, mobile_center, rotation = fit_kabsch(backbone.positions, reference_backbone)
        aligned_ligand = (ligand_positions - mobile_center) @ rotation + reference_backbone.mean(axis=0)
        protein_rmsd = rmsd(aligned_backbone, reference_backbone)
        ligand_rmsd = rmsd(aligned_ligand, reference_ligand)
        ligand_com = aligned_ligand.mean(axis=0)
        reference_com = reference_ligand.mean(axis=0)
        com_displacement = float(np.linalg.norm(ligand_com - reference_com))
        distances = distance_array(ligand_positions, pocket_atoms.positions, box=ts.dimensions)
        nearest = float(distances.min())
        retention = []
        for resid in pocket_resids_array:
            atom_mask = pocket_atoms.resids == resid
            retention.append(float(distances[:, atom_mask].min() <= 4.5))
        rows.append({
            "frame": frame_index,
            "time_ps": float(ts.time),
            "protein_backbone_rmsd_A": protein_rmsd,
            "ligand_heavy_rmsd_A": ligand_rmsd,
            "ligand_com_displacement_A": com_displacement,
            "pocket_contact_retention": float(np.mean(retention)),
            "nearest_pocket_distance_A": nearest,
        })
    if not rows:
        raise RuntimeError("No frames found in requested segment")

    tag = f"{args.start_ns:g}_{args.end_ns:g}"
    csv_path = args.out_dir / f"pparg_dinp_{tag}ns_metrics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    times = np.asarray([r["time_ps"] for r in rows])
    arrays = {key: np.asarray([r[key] for r in rows], dtype=float) for key in rows[0] if key not in {"frame", "time_ps"}}
    first_mask = times <= min(times[0] + 1000.0, times[-1])
    last_mask = times >= max(times[-1] - 1000.0, times[0])

    def stats(values):
        return {
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "first_1ns_mean": float(np.mean(values[first_mask])),
            "last_1ns_mean": float(np.mean(values[last_mask])),
        }

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "topology": str(args.topology),
        "trajectory": str(args.trajectory),
        "segment": f"{args.start_ns:g}-{args.end_ns:g} ns of available production frames",
        "n_atoms": universe.atoms.n_atoms,
        "n_frames": len(rows),
        "last_time_ps": float(times[-1]),
        "ligand_resname": ligand_residue.resname,
        "ligand_resid": int(ligand_residue.resid),
        "ligand_atoms": ligand.n_atoms,
        "ligand_heavy_atoms": ligand_heavy.n_atoms,
        "initial_pocket_resids": pocket_resids,
        "initial_pocket_residue_count": len(pocket_resids),
        "metrics": {key: stats(values) for key, values in arrays.items()},
        "interpretation_boundary": "This is a structural stability check; it does not establish experimental affinity or causality.",
    }
    summary_path = args.out_dir / f"pparg_dinp_{tag}ns_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(4, 1, figsize=(9, 11), sharex=True)
        x = times / 1000.0
        axes[0].plot(x, arrays["protein_backbone_rmsd_A"], lw=0.8)
        axes[0].set_ylabel("Protein BB RMSD (A)")
        axes[1].plot(x, arrays["ligand_heavy_rmsd_A"], lw=0.8, color="#c0392b")
        axes[1].set_ylabel("DINP RMSD (A)")
        axes[2].plot(x, arrays["ligand_com_displacement_A"], lw=0.8, color="#8e44ad")
        axes[2].set_ylabel("DINP COM shift (A)")
        axes[3].plot(x, arrays["pocket_contact_retention"], lw=0.8, color="#16803c")
        axes[3].set_ylabel("Pocket contact retention")
        axes[3].set_xlabel("Production time (ns)")
        axes[3].set_ylim(-0.02, 1.02)
        fig.suptitle(f"PPARG-DINP MD: {args.start_ns:g}–{args.end_ns:g} ns")
        fig.tight_layout()
        fig.savefig(args.out_dir / f"pparg_dinp_{tag}ns_metrics.png", dpi=180)
        plt.close(fig)
    except Exception as exc:
        summary["plot_warning"] = repr(exc)
        summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
