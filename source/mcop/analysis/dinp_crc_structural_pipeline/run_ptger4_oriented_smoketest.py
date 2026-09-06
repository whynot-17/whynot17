#!/usr/bin/env python3
"""Prepare an OPM/PPM-oriented PTGER4–DINP complex and run a short membrane-MD smoke test.

This script fixes an important coordinate issue: when PTGER4 is re-oriented by OPM/PPM,
the DINP docking pose must receive the *same rigid-body transform* before membrane setup.
It therefore:
  1) aligns the original 9JQZ PTGER4 coordinates to a user-supplied oriented PTGER4 PDB;
  2) applies that transform to the validated DINP docking pose;
  3) audits the alignment and ligand placement;
  4) optionally builds POPC + 0.15 M NaCl and runs a short CPU smoke test.

Default smoke test: 1 ns equilibration + 5 ns production at 310 K, 2 fs.
The full 100 ns run should only start after this smoke test passes manual QC.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
from Bio.PDB import PDBParser, Superimposer

import run_ptger4_membrane_md as md

ROOT = Path(__file__).resolve().parent
REFINE = ROOT / "outputs" / "stage2_refine_and_rescue"
CONTROL = ROOT / "outputs" / "ptger4_control_docking"
OUT = ROOT / "outputs" / "ptger4_membrane_smoketest"


def save_json(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def find_existing(paths: Iterable[Path]) -> Path:
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError("None of the expected files exists:\n" + "\n".join(map(str, paths)))


def ca_map(pdb: Path) -> Dict[Tuple[str, int, str], object]:
    structure = PDBParser(QUIET=True).get_structure(pdb.stem, str(pdb))
    out = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " " or "CA" not in res:
                    continue
                key = (chain.id, int(res.id[1]), res.resname.strip())
                out[key] = res["CA"]
        break
    return out


def derive_transform(reference_pdb: Path, oriented_pdb: Path):
    ref = ca_map(reference_pdb)
    ori = ca_map(oriented_pdb)
    common = sorted(set(ref) & set(ori))
    if len(common) < 100:
        raise RuntimeError(
            f"Only {len(common)} common CA atoms found. Expected >=100; check chain/residue numbering."
        )
    moving = [ref[k] for k in common]       # original 9JQZ
    fixed = [ori[k] for k in common]        # OPM/PPM-oriented target
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    rot, tran = sup.rotran
    return np.asarray(rot, float), np.asarray(tran, float), float(sup.rms), len(common)


def transform_pdbqt(src: Path, dst: Path, rot: np.ndarray, tran: np.ndarray):
    lines: List[str] = []
    transformed = 0
    for line in src.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            try:
                xyz = np.array([
                    float(line[30:38]), float(line[38:46]), float(line[46:54])
                ])
                new = np.dot(xyz, rot) + tran
                line = (
                    line[:30]
                    + f"{new[0]:8.3f}{new[1]:8.3f}{new[2]:8.3f}"
                    + line[54:]
                )
                transformed += 1
            except ValueError:
                pass
        lines.append(line)
    if transformed == 0:
        raise RuntimeError("No ligand coordinates were transformed from PDBQT")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return transformed


def ligand_centroid_from_pdbqt(path: Path) -> np.ndarray:
    coords = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                pass
    if not coords:
        raise RuntimeError("No ligand atoms parsed")
    return np.asarray(coords, float).mean(axis=0)


def protein_centroid(pdb: Path) -> np.ndarray:
    coords = []
    s = PDBParser(QUIET=True).get_structure("p", str(pdb))
    for atom in s.get_atoms():
        if atom.element != "H":
            coords.append(atom.coord)
    return np.asarray(coords, float).mean(axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oriented-pdb", required=True, type=Path,
                    help="OPM/PPM-oriented PTGER4 PDB; membrane normal must be Z")
    ap.add_argument("--run", action="store_true",
                    help="Actually build membrane and run smoke test; otherwise preparation/audit only")
    ap.add_argument("--equilibration-ns", type=float, default=1.0)
    ap.add_argument("--production-ns", type=float, default=5.0)
    ap.add_argument("--platform", default="CPU", choices=["CPU", "CUDA", "OpenCL", "HIP", "auto"])
    ap.add_argument("--temperature-k", type=float, default=310.0)
    args = ap.parse_args()

    if not args.oriented_pdb.exists():
        raise FileNotFoundError(args.oriented_pdb)
    OUT.mkdir(parents=True, exist_ok=True)
    md.prepend_runtime_path()

    reference = find_existing([
        REFINE / "ptger4_rescue" / "PTGER4_9JQZ_chainA_clean.pdb",
        CONTROL / "PTGER4_9JQZ_chainA_clean.pdb",
    ])
    docked = find_existing([
        CONTROL / "DINP_PTGER4_pose.pdbqt",
        CONTROL / "DINP_PTGER4_9JQZ_pose.pdbqt",
        REFINE / "ptger4_rescue" / "DINP_PTGER4_9JQZ_rescue_pose.pdbqt",
    ])

    rot, tran, ca_rmsd, n_ca = derive_transform(reference, args.oriented_pdb)
    transformed_pose = OUT / "DINP_PTGER4_oriented_pose.pdbqt"
    n_lig = transform_pdbqt(docked, transformed_pose, rot, tran)

    lig_cent = ligand_centroid_from_pdbqt(transformed_pose)
    prot_cent = protein_centroid(args.oriented_pdb)
    centroid_distance = float(np.linalg.norm(lig_cent - prot_cent))

    orientation_audit = {
        "reference_protein": str(reference),
        "oriented_protein": str(args.oriented_pdb),
        "source_docking_pose": str(docked),
        "transformed_docking_pose": str(transformed_pose),
        "matched_CA_atoms": n_ca,
        "protein_alignment_CA_RMSD_A": ca_rmsd,
        "ligand_atoms_transformed": n_lig,
        "rotation_matrix": rot.tolist(),
        "translation_A": tran.tolist(),
        "oriented_ligand_centroid_A": lig_cent.tolist(),
        "oriented_protein_centroid_A": prot_cent.tolist(),
        "ligand_to_protein_centroid_distance_A": centroid_distance,
        "membrane_normal_expected": "Z",
        "manual_qc_required": [
            "confirm OPM/PPM membrane orientation is valid",
            "confirm DINP remains in the A1ECR-defined orthosteric pocket after transform",
            "inspect protein/lipid clashes after membrane construction",
        ],
    }
    save_json(orientation_audit, OUT / "orientation_audit.json")

    deps = md.dependency_report()
    save_json(deps, OUT / "dependency_report.json")

    if not args.run:
        print(json.dumps({"status": "orientation_prepared", **orientation_audit}, indent=2, ensure_ascii=False))
        return

    required = ["openmm", "openff_toolkit", "openmmforcefields", "pdbfixer"]
    missing = [x for x in required if not deps.get(x)]
    if missing or not deps.get("obabel"):
        raise RuntimeError(f"MD dependency check failed: missing={missing}, obabel={deps.get('obabel')}")

    # Repair the already-oriented protein; rigid-body orientation is preserved.
    prepared = OUT / "PTGER4_oriented_prepared.pdb"
    repair = md.prepare_protein(args.oriented_pdb, prepared)
    save_json(repair, OUT / "protein_repair_audit.json")

    ligand_dir = OUT / "ligand"
    ligand_sdf, ligand_pdb = md.convert_ligand(str(deps["obabel"]), transformed_pose, ligand_dir)

    built = md.build_system(
        prepared,
        ligand_pdb,
        ligand_sdf,
        OUT,
        membrane_padding_nm=1.0,
        ionic_strength_m=0.15,
        temperature_k=args.temperature_k,
    )
    save_json({
        "n_atoms": built["n_atoms"],
        "barostat": built["barostat"],
        "residue_counts": built["residue_counts"],
        "built_pdb": built["built_pdb"],
        "system_xml": built["system_xml"],
    }, OUT / "build_audit.json")

    runinfo = md.run_md(
        built=built,
        outdir=OUT,
        temperature_k=args.temperature_k,
        timestep_fs=2.0,
        equilibration_ns=args.equilibration_ns,
        production_ns=args.production_ns,
        report_ps=10.0,
        checkpoint_ps=100.0,
        platform_name=args.platform,
        resume=False,
    )
    save_json(runinfo, OUT / "smoketest_run_audit.json")

    qc = md.basic_trajectory_analysis(Path(str(built["built_pdb"])), Path(runinfo["dcd"]), OUT)
    save_json(qc, OUT / "smoketest_trajectory_qc.json")

    summary = [
        "# PTGER4-DINP oriented membrane MD smoke test",
        "",
        f"- OPM/PPM alignment CA RMSD: {ca_rmsd:.4f} A ({n_ca} CA atoms)",
        f"- Equilibration: {args.equilibration_ns:.3f} ns",
        f"- Production smoke test: {args.production_ns:.3f} ns",
        f"- Platform: {runinfo['platform']}",
        f"- System atoms: {built['n_atoms']}",
        f"- Trajectory QC status: {qc.get('status')}",
        "",
        "Do not start the 100 ns production run until membrane placement, lipid clashes,",
        "protein stability, and DINP pocket retention are manually inspected.",
    ]
    (OUT / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
