"""Calculate MiNP–PPARG radius-of-gyration and SASA metrics from the 100 ns DCD.

The raw trajectory remains on the local E: drive.  This script writes compact
sampled tables to the repository module and is intentionally independent of
the RMSD/contact analysis.
"""

from __future__ import annotations

import csv
import json
import warnings
from pathlib import Path

import freesasa
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array


ROOT = Path(r"E:\chatgpt\whynot17\analysis\pparg_minp_md_100ns")
OUT = ROOT / "outputs"
TOPOLOGY = Path(r"E:\chatgpt\pparg_minp_md\system\seed20260917\initial.pdb")
TRAJECTORY = Path(r"E:\chatgpt\pparg_minp_md\run_20260915_seed20260917\production.dcd")
SAMPLE_STRIDE = 10  # production DCD is 10 ps/frame; report every 0.1 ns
CUTOFF_A = 4.5
CORE_INDICES = np.arange(12, 18, dtype=int)  # disjoint MiNP aromatic core


def _metadata(atoms):
    """Return FreeSASA metadata for one static atom group."""
    names = [str(x)[:4] for x in atoms.names]
    resnames = [str(x)[:3] for x in atoms.resnames]
    # The source PDB lacks residue numbers; resindices provide stable IDs.
    resnums = [int(x) + 1 for x in atoms.resindices]
    chains = [str(x)[:1] or "A" for x in atoms.segids]
    return names, resnames, resnums, chains


def _structure(meta, xyz):
    names, resnames, resnums, chains = meta
    structure = freesasa.Structure()
    structure.addAtoms(
        names,
        resnames,
        resnums,
        chains,
        xyz[:, 0].astype(float).tolist(),
        xyz[:, 1].astype(float).tolist(),
        xyz[:, 2].astype(float).tolist(),
    )
    return structure


def _sasa(meta, xyz):
    result = freesasa.calc(_structure(meta, xyz))
    return np.fromiter((result.atomArea(i) for i in range(result.nAtoms())), dtype=float)


def _mass_weighted_rg(xyz, masses):
    center = np.average(xyz, axis=0, weights=masses)
    return float(np.sqrt(np.average(np.sum((xyz - center) ** 2, axis=1), weights=masses)))


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _summary(rows):
    fields = [
        "protein_rg_A",
        "ligand_rg_A",
        "protein_sasa_complex_A2",
        "protein_sasa_isolated_A2",
        "protein_buried_by_ligand_A2",
        "ligand_sasa_complex_A2",
        "complex_sasa_A2",
        "ligand_sasa_isolated_A2",
        "ligand_buried_A2",
        "initial_pocket_sasa_complex_A2",
    ]
    windows = {
        "0_10ns": (0.0, 10.0),
        "10_50ns": (10.0, 50.0),
        "50_80ns": (50.0, 80.0),
        "80_100ns": (80.0, 100.0),
        "0_100ns": (0.0, 100.0),
    }
    output = {"n_samples": len(rows), "windows": {}}
    for label, (lo, hi) in windows.items():
        subset = [r for r in rows if lo <= r["time_ns"] < hi or (label == "0_100ns" and lo <= r["time_ns"] <= hi)]
        output["windows"][label] = {
            "n_samples": len(subset),
            **{
                field: {
                    "mean": float(np.mean([r[field] for r in subset])),
                    "sd": float(np.std([r[field] for r in subset], ddof=1)) if len(subset) > 1 else 0.0,
                    "min": float(np.min([r[field] for r in subset])),
                    "max": float(np.max([r[field] for r in subset])),
                }
                for field in fields
            },
        }
    return output


def main():
    warnings.filterwarnings("ignore", message="PDB file is missing resid information")
    freesasa.setVerbosity(freesasa.nowarnings)
    OUT.mkdir(parents=True, exist_ok=True)

    u = mda.Universe(str(TOPOLOGY), str(TRAJECTORY))
    protein = u.select_atoms("protein and not name H*")
    ligand = u.select_atoms("resname UNK and not name H*")
    if len(ligand) != 21:
        raise RuntimeError(f"Expected 21 MiNP heavy atoms, found {len(ligand)}")

    combined = protein + ligand
    protein_meta = _metadata(protein)
    ligand_meta = _metadata(ligand)
    combined_meta = _metadata(combined)
    masses_protein = np.asarray(protein.masses, dtype=float)
    masses_ligand = np.asarray(ligand.masses, dtype=float)
    protein_n = len(protein)

    # Initial pocket is defined from the six-atom MiNP core at the initial pose.
    initial_protein = protein.positions.copy()
    initial_ligand = ligand.positions.copy()
    initial_core = initial_ligand[CORE_INDICES]
    initial_dist = distance_array(initial_protein, initial_core)
    pocket_indices = np.where(initial_dist.min(axis=1) <= CUTOFF_A)[0]
    if len(pocket_indices) == 0:
        raise RuntimeError("No initial pocket atoms found")

    rows = []
    for frame_index, ts in enumerate(u.trajectory):
        if frame_index % SAMPLE_STRIDE:
            continue
        protein_xyz = protein.positions.copy()
        ligand_xyz = ligand.positions.copy()
        box = np.asarray(ts.dimensions[:3], dtype=float)
        if np.all(box > 0):
            protein_center = protein_xyz.mean(axis=0)
            ligand_center = ligand_xyz.mean(axis=0)
            ligand_xyz += np.rint((protein_center - ligand_center) / box) * box
        combined_xyz = np.vstack([protein_xyz, ligand_xyz])

        areas = _sasa(combined_meta, combined_xyz)
        protein_isolated = _sasa(protein_meta, protein_xyz)
        ligand_isolated = _sasa(ligand_meta, ligand_xyz)
        protein_complex = float(areas[:protein_n].sum())
        protein_iso = float(protein_isolated.sum())
        ligand_complex = float(areas[protein_n:].sum())
        ligand_iso = float(ligand_isolated.sum())
        rows.append(
            {
                "time_ns": float(ts.time / 1000.0),
                "protein_rg_A": _mass_weighted_rg(protein_xyz, masses_protein),
                "ligand_rg_A": _mass_weighted_rg(ligand_xyz, masses_ligand),
                "protein_sasa_complex_A2": protein_complex,
                "protein_sasa_isolated_A2": protein_iso,
                "protein_buried_by_ligand_A2": float(protein_iso - protein_complex),
                "ligand_sasa_complex_A2": ligand_complex,
                "complex_sasa_A2": float(areas.sum()),
                "ligand_sasa_isolated_A2": ligand_iso,
                "ligand_buried_A2": float(ligand_iso - ligand_complex),
                "initial_pocket_sasa_complex_A2": float(areas[pocket_indices].sum()),
            }
        )
        if len(rows) % 100 == 0:
            print(f"sampled {len(rows)} frames ({rows[-1]['time_ns']:.2f} ns)", flush=True)

    csv_path = OUT / "rg_sasa_0_100ns.csv"
    window_csv_path = OUT / "rg_sasa_window_summary.csv"
    summary_path = OUT / "rg_sasa_0_100ns_summary.json"
    manifest_path = OUT / "rg_sasa_manifest.json"
    _write_csv(csv_path, rows)
    summary = _summary(rows)
    window_fields = [
        "protein_rg_A", "ligand_rg_A", "protein_sasa_complex_A2",
        "protein_sasa_isolated_A2", "protein_buried_by_ligand_A2",
        "ligand_sasa_complex_A2", "complex_sasa_A2", "ligand_sasa_isolated_A2",
        "ligand_buried_A2", "initial_pocket_sasa_complex_A2",
    ]
    window_rows = []
    for window, metrics in summary["windows"].items():
        row = {"window": window, "n_samples": metrics["n_samples"]}
        for field in window_fields:
            row[f"{field}_mean"] = metrics[field]["mean"]
            row[f"{field}_sd"] = metrics[field]["sd"]
        window_rows.append(row)
    _write_csv(window_csv_path, window_rows)
    summary.update(
        {
            "method": "Mass-weighted protein/ligand radius of gyration plus FreeSASA Lee-Richards SASA",
            "sasa_atom_selection": "protein and resname UNK heavy atoms; hydrogens excluded",
            "freesasa_probe_radius_A": 1.4,
            "pocket_definition": "protein heavy atoms within 4.5 A of initial MiNP core indices 12-17",
            "sample_stride_frames": SAMPLE_STRIDE,
            "topology": str(TOPOLOGY),
            "trajectory": str(TRAJECTORY),
            "n_protein_heavy_atoms": int(len(protein)),
            "n_ligand_heavy_atoms": int(len(ligand)),
            "n_initial_pocket_atoms": int(len(pocket_indices)),
        }
    )
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    manifest = {
        "script": "analysis/pparg_minp_md_100ns/calculate_rg_sasa.py",
        "outputs": [csv_path.name, window_csv_path.name, summary_path.name],
        "raw_trajectory_not_committed": True,
        "summary": summary,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(summary["windows"]["80_100ns"], indent=2))


if __name__ == "__main__":
    main()
