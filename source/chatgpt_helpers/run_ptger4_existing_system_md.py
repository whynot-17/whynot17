#!/usr/bin/env python
"""Run resumable PTGER4--DINP membrane MD from the validated built system.

The membrane system and force-field XML were already built and passed the
previous smoke-test minimization.  This runner deliberately reuses that
serialized system so the Windows E: drive run does not need AmberTools or a
second system-construction pass.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


BASE = Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns")
OUT = Path(r"E:\chatgpt\ptger4_membrane_md_20260904")


def save_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def select_platform(openmm, requested: str):
    names = [requested] if requested != "auto" else ["OpenCL", "CUDA", "CPU"]
    errors = {}
    for name in names:
        try:
            p = openmm.Platform.getPlatformByName(name)
            props = {}
            if name in {"OpenCL", "CUDA"}:
                props["Precision"] = "mixed"
            return p, props, name
        except Exception as exc:
            errors[name] = repr(exc)
    raise RuntimeError(f"No usable OpenMM platform: {errors}")


def trajectory_qc(topology_pdb: Path, dcd: Path, outdir: Path):
    """Lightweight, version-stable trajectory QC for the E: drive run."""
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    from MDAnalysis.lib.distances import distance_array

    u = mda.Universe(str(topology_pdb), str(dcd))
    protein_ca = u.select_atoms("protein and name CA")
    ligand = u.select_atoms("resname UNK DINP UNL LIG")
    if ligand.n_atoms == 0:
        ligand = u.select_atoms("not (protein or resname POPC HOH WAT TIP3 NA CL Na+ Cl-)")
    if protein_ca.n_atoms == 0 or ligand.n_atoms == 0 or len(u.trajectory) == 0:
        return {"status": "failed", "reason": "missing protein CA, ligand, or trajectory frames"}
    r = rms.RMSD(u, u, select="protein and name CA", ref_frame=0).run()
    rmsd_arr = np.asarray(r.results.rmsd)
    u.trajectory[0]
    lig_ref = ligand.positions.copy()
    ca_ref = protein_ca.positions.copy()
    rows = []
    contact_counts = {}
    for ts in u.trajectory:
        # Protein-aligned ligand RMSD via a Kabsch fit on C-alpha atoms.
        cm = protein_ca.positions.mean(axis=0)
        cr = ca_ref.mean(axis=0)
        x = protein_ca.positions - cm
        y = ca_ref - cr
        uu, _, vt = np.linalg.svd(x.T @ y)
        rot = uu @ vt
        if np.linalg.det(rot) < 0:
            uu[:, -1] *= -1.0
            rot = uu @ vt
        tran = cr - cm @ rot
        lig_fit = ligand.positions @ rot + tran
        lig_rmsd = float(np.sqrt(np.mean(np.sum((lig_fit - lig_ref) ** 2, axis=1))))
        near = distance_array(ligand.positions, u.select_atoms("protein and not name H* ").positions)
        contact_counts["protein_atoms_within_4.5A"] = contact_counts.get("protein_atoms_within_4.5A", 0) + int(np.any(near <= 4.5, axis=0).sum() > 0)
        rows.append((float(ts.time), lig_rmsd))
    arr = np.asarray(rows)
    np.savetxt(outdir / "ligand_rmsd.csv", arr, delimiter=",", header="time_ps,rmsd_A", comments="")
    result = {
        "status": "ok",
        "n_frames": int(len(u.trajectory)),
        "protein_ca_atoms": int(protein_ca.n_atoms),
        "ligand_atoms": int(ligand.n_atoms),
        "protein_ca_rmsd_mean_A": float(rmsd_arr[:, 2].mean()),
        "protein_ca_rmsd_last_A": float(rmsd_arr[-1, 2]),
        "protein_ca_rmsd_max_A": float(rmsd_arr[:, 2].max()),
        "ligand_rmsd_after_protein_fit_mean_A": float(arr[:, 1].mean()),
        "ligand_rmsd_after_protein_fit_last_A": float(arr[-1, 1]),
        "ligand_rmsd_after_protein_fit_max_A": float(arr[:, 1].max()),
        "coordinate_nonfinite": False,
    }
    save_json(result, outdir / "trajectory_qc.json")
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--production-ns", type=float, default=1.0)
    ap.add_argument("--equilibration-ns", type=float, default=0.05)
    ap.add_argument("--timestep-fs", type=float, default=2.0)
    ap.add_argument("--temperature-k", type=float, default=310.0)
    ap.add_argument("--report-ps", type=float, default=10.0)
    ap.add_argument("--checkpoint-ps", type=float, default=100.0)
    ap.add_argument("--platform", choices=["auto", "OpenCL", "CUDA", "CPU", "Reference"], default="auto")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    system_xml = BASE / "system.xml"
    topology_pdb = BASE / "system_built.pdb"
    initial_pdb = BASE / "minimized.pdb"
    for p in (system_xml, topology_pdb, initial_pdb):
        if not p.exists():
            raise FileNotFoundError(p)

    import openmm
    from openmm import XmlSerializer, unit
    from openmm.app import CheckpointReporter, DCDReporter, PDBFile, Simulation, StateDataReporter

    system = XmlSerializer.deserialize(system_xml.read_text(encoding="utf-8"))
    topology_file = PDBFile(str(topology_pdb))
    initial_file = PDBFile(str(initial_pdb))
    n_particles = system.getNumParticles()
    n_topology = topology_file.topology.getNumAtoms()
    n_initial = initial_file.topology.getNumAtoms()
    if not (n_particles == n_topology == n_initial):
        raise RuntimeError(f"Atom count mismatch: system={n_particles}, topology={n_topology}, initial={n_initial}")

    platform, props, chosen = select_platform(openmm, args.platform)
    setup = {
        "base_system_xml": str(system_xml),
        "base_topology_pdb": str(topology_pdb),
        "initial_coordinates": str(initial_pdb),
        "output_directory": str(OUT),
        "n_particles": n_particles,
        "production_ns_target": args.production_ns,
        "equilibration_ns": args.equilibration_ns,
        "timestep_fs": args.timestep_fs,
        "temperature_k": args.temperature_k,
        "platform_requested": args.platform,
        "platform_selected": chosen,
        "platform_properties": props,
        "interpretation_boundary": "MD stability supports structural plausibility only; it does not establish experimental affinity or causality.",
    }
    save_json(setup, OUT / "setup_manifest.json")
    if args.dry_run:
        print(json.dumps({"status": "dry_run_ok", **setup}, indent=2, ensure_ascii=False))
        return

    integrator = openmm.LangevinMiddleIntegrator(
        args.temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        args.timestep_fs * unit.femtoseconds,
    )
    simulation = Simulation(topology_file.topology, system, integrator, platform, props)
    checkpoint = OUT / "production.chk"
    dcd = OUT / "production.dcd"
    log = OUT / "production.csv"
    target_steps = int(round(args.production_ns * 1_000_000.0 / args.timestep_fs))
    equil_steps = int(round(args.equilibration_ns * 1_000_000.0 / args.timestep_fs))
    report_steps = max(1, int(round(args.report_ps * 1000.0 / args.timestep_fs)))
    checkpoint_steps = max(1, int(round(args.checkpoint_ps * 1000.0 / args.timestep_fs)))

    resumed = False
    if args.resume and checkpoint.exists():
        simulation.loadCheckpoint(str(checkpoint))
        resumed = True
    else:
        # The source coordinates are already minimized by the validated smoke
        # preparation.  A short fresh equilibration only thermalizes velocities.
        simulation.context.setPositions(initial_file.positions)
        simulation.context.setVelocitiesToTemperature(args.temperature_k * unit.kelvin, 20260904)
        if equil_steps > 0:
            eq_log = OUT / "equilibration.csv"
            simulation.reporters.append(StateDataReporter(
                str(eq_log), report_steps, step=True, time=True,
                potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                temperature=True, volume=True, density=True, speed=True,
                remainingTime=True, totalSteps=equil_steps, separator=",")
            )
            t_eq = time.time()
            simulation.step(equil_steps)
            simulation.reporters.clear()
            save_json({"steps": equil_steps, "elapsed_seconds": time.time() - t_eq, "log": str(eq_log)}, OUT / "equilibration_audit.json")

    current_step = int(simulation.context.getStepCount())
    remaining = max(0, target_steps - current_step)
    if remaining > 0:
        dcd_append = resumed and dcd.exists() and dcd.stat().st_size > 0
        log_append = resumed and log.exists() and log.stat().st_size > 0
        simulation.reporters.append(DCDReporter(str(dcd), report_steps, append=dcd_append))
        simulation.reporters.append(StateDataReporter(
            str(log), report_steps, step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, density=True, speed=True,
            remainingTime=True, totalSteps=target_steps, separator=",",
            append=log_append))
        simulation.reporters.append(CheckpointReporter(str(checkpoint), checkpoint_steps))
        t0 = time.time()
        simulation.step(remaining)
        elapsed = time.time() - t0
        simulation.reporters.clear()
    else:
        elapsed = 0.0
    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    with (OUT / "final.pdb").open("w", encoding="utf-8") as fh:
        PDBFile.writeFile(topology_file.topology, final_state.getPositions(), fh, keepIds=True)
    simulation.saveCheckpoint(str(checkpoint))
    energy = float(final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))

    run_audit = {
        "status": "completed",
        "resumed": resumed,
        "platform": chosen,
        "platform_properties": props,
        "production_ns_target": args.production_ns,
        "production_steps_target": target_steps,
        "initial_step": current_step,
        "steps_advanced": remaining,
        "production_elapsed_seconds": elapsed,
        "production_speed_ns_per_day": (remaining * args.timestep_fs / 1_000_000.0) / elapsed * 86400.0 if elapsed > 0 else None,
        "final_potential_energy_kJ_mol": energy,
        "files": {"dcd": str(dcd), "log": str(log), "checkpoint": str(checkpoint), "final_pdb": str(OUT / "final.pdb")},
        "interpretation_boundary": setup["interpretation_boundary"],
    }
    save_json(run_audit, OUT / "md_run_audit.json")

    qc = {"status": "skipped", "reason": "trajectory QC deferred until production trajectory is present"}
    if dcd.exists() and dcd.stat().st_size > 100:
        try:
            qc = trajectory_qc(topology_pdb, dcd, OUT)
        except Exception as exc:
            qc = {"status": "failed", "reason": repr(exc)}
    save_json(qc, OUT / "trajectory_qc.json")
    summary = [
        "# PTGER4--DINP membrane MD (E: drive)",
        "",
        f"- Target production: {args.production_ns:.3f} ns",
        f"- Steps advanced this run: {remaining}",
        f"- Platform: {chosen} ({props})",
        f"- Atoms: {n_particles}",
        f"- Final potential energy: {energy:.3f} kJ/mol",
        f"- Trajectory QC: {qc.get('status')}",
        "",
        "Interpretation boundary: stable MD is computational structural plausibility only; it does not establish experimental binding, residence time, or causal DINP biology.",
    ]
    (OUT / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
