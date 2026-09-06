#!/usr/bin/env python3
"""Run a short, staged PTGER4--DINP membrane-MD smoke test.

This runner is intentionally separate from the long-production launcher.  It
uses the PPM-oriented PTGER4 coordinates and the DINP pose after the validated
rigid-body transfer, then performs:

    minimization -> restrained equilibration -> unrestrained equilibration
    -> 1 ns production

The output is a technical stability/QC record only.  It is not evidence of
experimental affinity, residence time, or causal biology.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent
PPM = ROOT / "outputs" / "ppm_local"
TRANSFER = ROOT / "outputs" / "ppm_pose_transfer"
REFINE = ROOT / "outputs" / "stage2_refine_and_rescue"
OUT = ROOT / "outputs" / "ptger4_membrane_smoke_1ns"
ORIENTED = PPM / "9JQZ_PTGER4_PPM_oriented.pdb"
COMPLEX = TRANSFER / "PTGER4_DINP_PPM_oriented_complex.pdb"
SOURCE_RECEPTOR = REFINE / "ptger4_rescue" / "PTGER4_9JQZ_chainA_clean.pdb"
SOURCE_POSE = REFINE / "ptger4_rescue" / "DINP_PTGER4_9JQZ_rescue_pose.pdbqt"


def save_json(obj: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def require_files(paths: Iterable[Path]) -> None:
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError("Required input(s) missing:\n" + "\n".join(missing))


def sanitize_openbabel_pdb(raw_path: Path, clean_path: Path) -> int:
    """Keep coordinates and CONECT records, remove nonstandard MODEL metadata."""

    kept = []
    atom_count = 0
    for line in raw_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        record = line[:6].strip()
        if record in {"ATOM", "HETATM"}:
            atom_count += 1
            try:
                serial = int(line[6:11])
            except ValueError as exc:
                raise RuntimeError(f"Invalid ligand atom serial in {raw_path}: {line!r}") from exc
            element = line[76:78].strip() or line[12:16].strip()[:1] or "C"
            # OpenBabel writes element-only atom names (C/O), which OpenMM
            # collapses as duplicate atoms within the UNL residue.  Give each
            # atom a unique PDB name while leaving serials, elements, bonds,
            # and all coordinate columns untouched.
            unique_name = f"{element[0].upper()}{serial:03d}"[-4:]
            kept.append(line[:12] + f"{unique_name:>4}" + line[16:])
        elif record in {"CONECT", "TER", "END"}:
            kept.append(line)
    if atom_count == 0 or not any(x.startswith("CONECT") for x in kept):
        raise RuntimeError("Sanitized DINP PDB lacks coordinates or CONECT topology")
    if kept[-1].strip() != "END":
        kept.append("END")
    clean_path.parent.mkdir(parents=True, exist_ok=True)
    clean_path.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return atom_count


def import_pipeline_modules():
    # The existing modules contain the validated PPM pose-transfer and force-
    # field setup code.  This runner adds only staged equilibration and QC.
    import run_ptger4_membrane_md as md
    import run_ptger4_oriented_smoketest as transfer

    return md, transfer


def finite_energy(simulation, unit) -> float:
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    positions_nm = np.asarray(
        state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float
    )
    if not np.isfinite(energy) or not np.isfinite(positions_nm).all():
        raise RuntimeError("Non-finite energy or coordinates detected")
    if np.max(np.abs(positions_nm)) > 1000.0:
        raise RuntimeError("Coordinate explosion detected during MD")
    return energy


def write_state(simulation, path: Path, PDBFile) -> float:
    from openmm import unit

    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with path.open("w", encoding="utf-8") as fh:
        PDBFile.writeFile(simulation.topology, state.getPositions(), fh, keepIds=True)
    return float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))


def add_protein_heavy_atom_restraint(system, modeller, initial_positions, openmm, unit):
    """Add a global-strength restraint that can be switched off in-context."""

    force = openmm.CustomExternalForce(
        "0.5*restraint_k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
    )
    force.addGlobalParameter("restraint_k", 1000.0)
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)

    indices: List[int] = []
    for atom in modeller.topology.atoms():
        # Chain A is the single PTGER4 receptor chain in the PPM input.  The
        # ligand is chain L, and addMembrane adds lipids/solvent separately.
        if atom.residue.chain.id != "A":
            continue
        if atom.element is not None and atom.element.symbol == "H":
            continue
        xyz = initial_positions[atom.index].value_in_unit(unit.nanometer)
        force.addParticle(atom.index, [float(xyz[0]), float(xyz[1]), float(xyz[2])])
        indices.append(atom.index)
    if len(indices) < 1000:
        raise RuntimeError(
            f"Only {len(indices)} receptor heavy atoms were restrained; expected a full GPCR"
        )
    force_index = system.addForce(force)
    return force_index, len(indices)


def stage_steps(ns: float, timestep_fs: float) -> int:
    return int(round(ns * 1_000_000.0 / timestep_fs))


def run_stage(simulation, steps: int, log_path: Path, report_steps: int, total_steps: int) -> Dict[str, object]:
    from openmm.app import StateDataReporter

    if steps <= 0:
        return {"steps": 0, "log": None}
    simulation.reporters.append(
        StateDataReporter(
            str(log_path),
            report_steps,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            remainingTime=True,
            totalSteps=total_steps,
            separator=",",
        )
    )
    t0 = time.time()
    simulation.step(steps)
    elapsed = time.time() - t0
    simulation.reporters.clear()
    return {"steps": steps, "log": str(log_path), "elapsed_seconds": elapsed}


def kabsch_row(moving: np.ndarray, reference: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return row-vector rotation/translation mapping moving onto reference."""

    cm = moving.mean(axis=0)
    cr = reference.mean(axis=0)
    x = moving - cm
    y = reference - cr
    u, _, vt = np.linalg.svd(x.T @ y)
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1.0
        r = u @ vt
    t = cr - cm @ r
    fitted = moving @ r + t
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - reference) ** 2, axis=1))))
    return r, t, rmsd


def select_ligand(u):
    ligand = u.select_atoms("resname DIN UNL LIG DINP")
    if ligand.n_atoms == 0:
        excluded = "protein or resname POPC HOH WAT TIP3 NA CL Na+ Cl-"
        ligand = u.select_atoms(f"not ({excluded})")
    if ligand.n_atoms == 0:
        raise RuntimeError("No DINP ligand atoms found in the built trajectory")
    return ligand


def pocket_residues_from_reference():
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    ref = mda.Universe(str(COMPLEX))
    ligand = ref.select_atoms("resname DIN UNL LIG DINP")
    protein = ref.select_atoms("protein")
    if ligand.n_atoms == 0 or protein.n_atoms == 0:
        raise RuntimeError("Could not identify protein/ligand in pose-transfer complex")
    keys = []
    for residue in protein.residues:
        distances = distance_array(residue.atoms.positions, ligand.positions)
        if float(distances.min()) <= 5.0:
            keys.append((str(residue.segid), int(residue.resid), str(residue.resname)))
    if len(keys) < 5:
        raise RuntimeError(f"Only {len(keys)} reference pocket residues found")
    return keys


def trajectory_qc(topology_pdb: Path, dcd: Path, production_csv: Path, outdir: Path) -> Dict[str, object]:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    u = mda.Universe(str(topology_pdb), str(dcd))
    ligand = select_ligand(u)
    protein_ca = u.select_atoms("protein and name CA")
    if protein_ca.n_atoms < 200:
        raise RuntimeError(f"Only {protein_ca.n_atoms} protein C-alpha atoms in trajectory")
    pocket_keys = pocket_residues_from_reference()
    pocket_resids = [x[1] for x in pocket_keys]
    pocket = u.select_atoms("protein and resid " + " ".join(str(x) for x in pocket_resids))
    if pocket.n_atoms == 0:
        raise RuntimeError("Reference pocket residues were not found in built topology")

    # The first production frame is the reference for all dynamic RMSDs.
    u.trajectory[0]
    ca_ref = protein_ca.positions.copy()
    lig_ref = ligand.positions.copy()
    n_frames = len(u.trajectory)
    rows = []
    contact_counts = {int(resid): 0 for resid in pocket_resids}
    thickness = []
    areas = []
    max_abs_coord = 0.0
    nonfinite = False

    headgroups = u.select_atoms("resname POPC and name P")
    if headgroups.n_atoms == 0:
        headgroups = u.select_atoms("resname POPC and name PO4")

    for ts in u.trajectory:
        coords = np.asarray(u.atoms.positions, dtype=float)
        if not np.isfinite(coords).all():
            nonfinite = True
        max_abs_coord = max(max_abs_coord, float(np.nanmax(np.abs(coords))))
        if ts.dimensions is None or not np.isfinite(ts.dimensions).all():
            nonfinite = True
        else:
            areas.append(float(ts.dimensions[0] * ts.dimensions[1]))

        r, t, ca_rmsd = kabsch_row(protein_ca.positions, ca_ref)
        lig_fit = ligand.positions @ r + t
        lig_rmsd = float(np.sqrt(np.mean(np.sum((lig_fit - lig_ref) ** 2, axis=1))))
        lig_com = ligand.positions.mean(axis=0)
        pocket_com = pocket.positions.mean(axis=0)
        pocket_distance = float(np.linalg.norm(lig_com - pocket_com))
        per_residue = {}
        for residue in pocket.residues:
            d = distance_array(ligand.positions, residue.atoms.positions)
            min_d = float(d.min())
            per_residue[int(residue.resid)] = min_d
            if min_d <= 4.5:
                contact_counts[int(residue.resid)] += 1
        rows.append((float(ts.time), ca_rmsd, lig_rmsd, pocket_distance))

        if headgroups.n_atoms >= 10:
            z = headgroups.positions[:, 2]
            center = float(np.median(z))
            upper = z[z > center]
            lower = z[z < center]
            if len(upper) >= 3 and len(lower) >= 3:
                thickness.append(float(np.mean(upper) - np.mean(lower)))

    arr = np.asarray(rows, dtype=float)
    np.savetxt(
        outdir / "trajectory_metrics.csv",
        arr,
        delimiter=",",
        header="time_ps,protein_ca_rmsd_A,dinp_rmsd_after_protein_fit_A,ligand_pocket_com_distance_A",
        comments="",
    )
    with (outdir / "pocket_contact_occupancy.csv").open("w", encoding="utf-8") as fh:
        fh.write("resid,resname,contact_cutoff_A,occupancy\n")
        for segid, resid, resname in pocket_keys:
            fh.write(f"{resid},{resname},4.5,{contact_counts.get(resid, 0) / max(1, n_frames):.6f}\n")

    log_values = []
    if production_csv.exists():
        with production_csv.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                parsed = {}
                for key, val in row.items():
                    if key is None or val is None:
                        continue
                    try:
                        parsed[key] = float(val)
                    except ValueError:
                        pass
                log_values.append(parsed)
    log_finite = True
    for row in log_values:
        if not all(np.isfinite(x) for x in row.values()):
            log_finite = False

    qc = {
        "status": "ok" if (not nonfinite and log_finite and max_abs_coord < 1000.0) else "failed",
        "n_frames": n_frames,
        "protein_ca_atoms": int(protein_ca.n_atoms),
        "ligand_atoms": int(ligand.n_atoms),
        "pocket_residues": [
            {"resid": resid, "resname": resname, "segid": segid} for segid, resid, resname in pocket_keys
        ],
        "protein_ca_rmsd_mean_A": float(arr[:, 1].mean()),
        "protein_ca_rmsd_last_A": float(arr[-1, 1]),
        "protein_ca_rmsd_max_A": float(arr[:, 1].max()),
        "dinp_rmsd_after_protein_fit_mean_A": float(arr[:, 2].mean()),
        "dinp_rmsd_after_protein_fit_last_A": float(arr[-1, 2]),
        "dinp_rmsd_after_protein_fit_max_A": float(arr[:, 2].max()),
        "ligand_pocket_com_distance_mean_A": float(arr[:, 3].mean()),
        "ligand_pocket_com_distance_last_A": float(arr[-1, 3]),
        "ligand_pocket_com_distance_min_A": float(arr[:, 3].min()),
        "ligand_pocket_com_distance_max_A": float(arr[:, 3].max()),
        "membrane_headgroup_selection": "POPC:P or POPC:PO4",
        "membrane_phosphate_leaflet_separation_mean_A": float(np.mean(thickness)) if thickness else None,
        "membrane_phosphate_leaflet_separation_min_A": float(np.min(thickness)) if thickness else None,
        "membrane_phosphate_leaflet_separation_max_A": float(np.max(thickness)) if thickness else None,
        "membrane_area_mean_A2": float(np.mean(areas)) if areas else None,
        "membrane_area_min_A2": float(np.min(areas)) if areas else None,
        "membrane_area_max_A2": float(np.max(areas)) if areas else None,
        "coordinate_nonfinite": bool(nonfinite),
        "production_log_rows": len(log_values),
        "production_log_finite": bool(log_finite),
        "max_absolute_coordinate_A": max_abs_coord,
        "explosion_check": bool(nonfinite or not log_finite or max_abs_coord >= 1000.0),
    }
    save_json(qc, outdir / "trajectory_qc.json")
    return qc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--production-ns", type=float, default=1.0)
    ap.add_argument("--restrained-equilibration-ns", type=float, default=0.25)
    ap.add_argument("--unrestrained-equilibration-ns", type=float, default=0.25)
    ap.add_argument("--timestep-fs", type=float, default=2.0)
    ap.add_argument("--temperature-k", type=float, default=310.0)
    ap.add_argument("--ionic-strength-m", type=float, default=0.15)
    ap.add_argument("--barostat-frequency", type=int, default=500,
                    help="Membrane barostat attempt interval in MD steps; 500 is frozen for this CPU smoke test")
    ap.add_argument("--platform", default="CPU", choices=["CPU", "Reference"])
    ap.add_argument("--report-ps", type=float, default=10.0)
    ap.add_argument("--checkpoint-ps", type=float, default=100.0)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if args.production_ns != 1.0:
        raise ValueError("This smoke-test runner is frozen to exactly 1.0 ns production")
    require_files([ORIENTED, COMPLEX, SOURCE_RECEPTOR, SOURCE_POSE])
    OUT.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("OPENMM_CPU_THREADS", "4")

    md, transfer = import_pipeline_modules()
    md.prepend_runtime_path()
    # In WSL, OpenBabel is installed beside the selected Python interpreter,
    # while the legacy launcher only knows the Windows-side runtime path.
    # Prefer the active environment so the dependency audit reflects the
    # executable that will actually be used.
    active_bin = Path(sys.executable).resolve().parent
    os.environ["PATH"] = os.pathsep.join([str(active_bin), os.environ.get("PATH", "")])
    deps = md.dependency_report()
    save_json(deps, OUT / "dependency_report.json")
    required = ["openmm", "openff_toolkit", "openmmforcefields", "pdbfixer"]
    missing = [x for x in required if not deps.get(x)]
    if missing or not deps.get("obabel"):
        raise RuntimeError(f"MD dependency check failed: missing={missing}, obabel={deps.get('obabel')}")
    if args.dry_run:
        print(json.dumps({"status": "dry_run_ok", "production_ns": 1.0, "platform": args.platform, "dependencies": deps}, indent=2))
        return

    from openmm import unit
    import openmm
    from openmm.app import PDBFile, Simulation

    # Recompute the pose transfer from the validated source files, rather than
    # using the old pre-PPM pose.  The QC JSON is saved for provenance.
    rot, tran, ca_rmsd, n_ca = transfer.derive_transform(SOURCE_RECEPTOR, ORIENTED)
    oriented_pose = OUT / "DINP_PTGER4_PPM_oriented_pose.pdbqt"
    n_lig = transfer.transform_pdbqt(SOURCE_POSE, oriented_pose, rot, tran)
    save_json({
        "source_receptor": str(SOURCE_RECEPTOR),
        "oriented_receptor": str(ORIENTED),
        "source_pose": str(SOURCE_POSE),
        "oriented_pose": str(oriented_pose),
        "matched_ca_atoms": n_ca,
        "receptor_alignment_rmsd_A": ca_rmsd,
        "ligand_atoms_transformed": n_lig,
        "interpretation": "Rigid-body coordinate transfer only; no redocking or ligand minimization.",
    }, OUT / "pose_transfer_audit.json")

    ligand_dir = OUT / "ligand"
    ligand_dir.mkdir(parents=True, exist_ok=True)
    # The validated pose-transfer PDB is the coordinate source.  OpenBabel's
    # PDB writer preserves the coordinates but may emit an empty MODEL field
    # from PDBQT metadata, which OpenMM rejects.  We therefore use OpenBabel
    # only for the chemical SDF and use the already QC-passed PPM-oriented PDB
    # for coordinates.
    ligand_sdf = ligand_dir / "DINP_docked_pose.sdf"
    cp = md.run([str(deps["obabel"]), str(oriented_pose), "-O", str(ligand_sdf), "-f", "1", "-l", "1"])
    (ligand_dir / "obabel_sdf.log").write_text(cp.stdout + "\n" + cp.stderr, encoding="utf-8")
    if cp.returncode != 0 or not ligand_sdf.exists() or ligand_sdf.stat().st_size == 0:
        raise RuntimeError("OpenBabel SDF conversion failed")
    raw_ligand_pdb = ligand_dir / "DINP_docked_pose_openbabel_raw.pdb"
    cp = md.run([str(deps["obabel"]), str(oriented_pose), "-O", str(raw_ligand_pdb), "-f", "1", "-l", "1"])
    (ligand_dir / "obabel_pdb.log").write_text(cp.stdout + "\n" + cp.stderr, encoding="utf-8")
    if cp.returncode != 0 or not raw_ligand_pdb.exists() or raw_ligand_pdb.stat().st_size == 0:
        raise RuntimeError("OpenBabel PDB conversion failed")
    ligand_pdb = ligand_dir / "DINP_PPM_oriented_coordinates_openmm.pdb"
    sanitize_openbabel_pdb(raw_ligand_pdb, ligand_pdb)
    # Preserve the independently QC-passed coordinate file alongside the
    # OpenMM-ready bonded topology input for direct comparison.
    shutil.copy2(TRANSFER / "DINP_PPM_oriented.pdb", ligand_dir / "DINP_PPM_oriented_coordinate_reference.pdb")
    prepared = OUT / "PTGER4_PPM_oriented_prepared.pdb"
    repair = md.prepare_protein(ORIENTED, prepared)
    save_json(repair, OUT / "protein_repair_audit.json")

    built = md.build_system(
        prepared,
        ligand_pdb,
        ligand_sdf,
        OUT,
        membrane_padding_nm=1.0,
        ionic_strength_m=args.ionic_strength_m,
        temperature_k=args.temperature_k,
        barostat_frequency=args.barostat_frequency,
    )
    save_json({
        "n_atoms": built["n_atoms"],
        "barostat": built["barostat"],
        "residue_counts": built["residue_counts"],
        "built_pdb": built["built_pdb"],
        "system_xml": built["system_xml"],
    }, OUT / "build_audit.json")

    restraint_index, restrained_atoms = add_protein_heavy_atom_restraint(
        built["system"], built["modeller"], built["modeller"].positions, openmm, unit
    )
    (OUT / "system_with_restraint.xml").write_text(
        openmm.XmlSerializer.serialize(built["system"]), encoding="utf-8"
    )

    integrator = openmm.LangevinMiddleIntegrator(
        args.temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        args.timestep_fs * unit.femtoseconds,
    )
    platform = openmm.Platform.getPlatformByName(args.platform)
    simulation = Simulation(built["modeller"].topology, built["system"], integrator, platform)
    simulation.context.setPositions(built["modeller"].positions)
    # A bounded 1 ns smoke-test minimization: enough to remove initial clashes
    # while keeping the CPU-only preflight tractable for this large membrane
    # system.  The finite-energy/coordinate gate below must still pass before
    # either equilibration stage is allowed to start.
    min_energy = simulation.minimizeEnergy(maxIterations=1000)
    min_potential = finite_energy(simulation, unit)
    write_state(simulation, OUT / "minimized.pdb", PDBFile)

    simulation.context.setVelocitiesToTemperature(args.temperature_k * unit.kelvin, 20260904)
    restrained_steps = stage_steps(args.restrained_equilibration_ns, args.timestep_fs)
    unrestrained_steps = stage_steps(args.unrestrained_equilibration_ns, args.timestep_fs)
    production_steps = stage_steps(args.production_ns, args.timestep_fs)
    report_steps = max(1, stage_steps(args.report_ps / 1000.0, args.timestep_fs))
    checkpoint_steps = max(1, stage_steps(args.checkpoint_ps / 1000.0, args.timestep_fs))

    restrained = run_stage(
        simulation,
        restrained_steps,
        OUT / "restrained_equilibration.csv",
        report_steps,
        restrained_steps,
    )
    restrained["potential_energy_kJ_mol"] = finite_energy(simulation, unit)
    write_state(simulation, OUT / "restrained_equilibrated.pdb", PDBFile)

    simulation.context.setParameter("restraint_k", 0.0)
    unrestrained = run_stage(
        simulation,
        unrestrained_steps,
        OUT / "unrestrained_equilibration.csv",
        report_steps,
        unrestrained_steps,
    )
    unrestrained["potential_energy_kJ_mol"] = finite_energy(simulation, unit)
    write_state(simulation, OUT / "unrestrained_equilibrated.pdb", PDBFile)

    # Production reporters are installed only after both equilibration stages.
    from openmm.app import CheckpointReporter, DCDReporter, StateDataReporter

    dcd = OUT / "production_1ns.dcd"
    production_csv = OUT / "production_1ns.csv"
    checkpoint = OUT / "production_1ns.chk"
    simulation.reporters.append(DCDReporter(str(dcd), report_steps))
    simulation.reporters.append(
        StateDataReporter(
            str(production_csv),
            report_steps,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            remainingTime=True,
            totalSteps=production_steps,
            separator=",",
        )
    )
    simulation.reporters.append(CheckpointReporter(str(checkpoint), checkpoint_steps))
    t0 = time.time()
    simulation.step(production_steps)
    elapsed = time.time() - t0
    simulation.reporters.clear()
    production_energy = finite_energy(simulation, unit)
    write_state(simulation, OUT / "final.pdb", PDBFile)
    simulation.saveCheckpoint(str(checkpoint))

    qc = trajectory_qc(Path(str(built["built_pdb"])), dcd, production_csv, OUT)
    speed = args.production_ns / elapsed * 86400.0 if elapsed > 0 else None
    run_audit = {
        "status": "completed",
        "platform": args.platform,
        "openmm_version": deps.get("openmm_version"),
        "production_ns": args.production_ns,
        "restrained_equilibration_ns": args.restrained_equilibration_ns,
        "unrestrained_equilibration_ns": args.unrestrained_equilibration_ns,
        "timestep_fs": args.timestep_fs,
        "temperature_k": args.temperature_k,
        "ionic_strength_m": args.ionic_strength_m,
        "barostat_frequency_steps": args.barostat_frequency,
        "restrained_protein_heavy_atoms": restrained_atoms,
        "minimization_potential_energy_kJ_mol": min_potential,
        "restrained_equilibration": restrained,
        "unrestrained_equilibration": unrestrained,
        "production_steps": production_steps,
        "production_elapsed_seconds": elapsed,
        "production_speed_ns_per_day": speed,
        "production_final_potential_energy_kJ_mol": production_energy,
        "files": {"dcd": str(dcd), "log": str(production_csv), "checkpoint": str(checkpoint)},
        "trajectory_qc_status": qc["status"],
        "interpretation_boundary": "Smoke-test stability is computational structural plausibility only; it does not establish binding affinity, residence time, or causal DINP biology.",
    }
    save_json(run_audit, OUT / "smoketest_1ns_run_audit.json")

    summary = [
        "# PTGER4--DINP 1 ns membrane-MD smoke test",
        "",
        "## Frozen protocol",
        "",
        f"- Receptor: PPM-oriented 9JQZ PTGER4 (`{ORIENTED.name}`)",
        f"- Ligand: DINP pose after validated rigid-body transfer (`{oriented_pose.name}`)",
        "- Membrane: POPC; 0.15 M NaCl; 310 K; 1 bar; PME",
        f"- Membrane barostat attempt interval: every {args.barostat_frequency} steps ({args.barostat_frequency * args.timestep_fs / 1000.0:.3f} ps)",
        f"- Workflow: minimization -> {args.restrained_equilibration_ns:.3f} ns restrained -> {args.unrestrained_equilibration_ns:.3f} ns unrestrained -> {args.production_ns:.3f} ns production",
        f"- Platform: {args.platform}; timestep {args.timestep_fs:.1f} fs",
        "",
        "## Technical QC",
        "",
        f"- Production wall time: {elapsed:.1f} s",
        f"- Observed production speed: {speed:.3f} ns/day" if speed is not None else "- Observed production speed: unavailable",
        f"- Protein C-alpha RMSD (mean / max): {qc['protein_ca_rmsd_mean_A']:.3f} / {qc['protein_ca_rmsd_max_A']:.3f} A",
        f"- DINP RMSD after protein alignment (mean / max): {qc['dinp_rmsd_after_protein_fit_mean_A']:.3f} / {qc['dinp_rmsd_after_protein_fit_max_A']:.3f} A",
        f"- DINP--pocket COM distance (mean / range): {qc['ligand_pocket_com_distance_mean_A']:.3f} / {qc['ligand_pocket_com_distance_min_A']:.3f}--{qc['ligand_pocket_com_distance_max_A']:.3f} A",
        f"- Phosphate leaflet separation (mean / range): {qc['membrane_phosphate_leaflet_separation_mean_A']:.3f} / {qc['membrane_phosphate_leaflet_separation_min_A']:.3f}--{qc['membrane_phosphate_leaflet_separation_max_A']:.3f} A" if qc["membrane_phosphate_leaflet_separation_mean_A"] is not None else "- Phosphate leaflet separation: unavailable",
        f"- Membrane XY area (mean / range): {qc['membrane_area_mean_A2']:.1f} / {qc['membrane_area_min_A2']:.1f}--{qc['membrane_area_max_A2']:.1f} A2" if qc["membrane_area_mean_A2"] is not None else "- Membrane XY area: unavailable",
        f"- NaN/explosion check: {'FAIL' if qc['explosion_check'] else 'PASS'}",
        f"- Overall trajectory QC: **{qc['status'].upper()}**",
        "",
        "## Boundary",
        "",
        "This smoke test asks whether the membrane system is numerically and geometrically tractable for a longer controlled simulation. It does not prove experimental binding, pharmacology, residence time, or a DINP->PTGER4->CRC causal mechanism.",
    ]
    (OUT / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
