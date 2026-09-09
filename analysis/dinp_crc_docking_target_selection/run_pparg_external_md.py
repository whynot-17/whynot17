from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def save_json(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def choose_platform(openmm, requested):
    names = [requested] if requested != "auto" else ["CUDA", "OpenCL", "CPU"]
    errors = {}
    for name in names:
        try:
            platform = openmm.Platform.getPlatformByName(name)
            props = {}
            if name in {"CUDA", "OpenCL"}:
                props["Precision"] = "mixed"
            return platform, props, name
        except Exception as exc:
            errors[name] = repr(exc)
    raise RuntimeError(f"No usable OpenMM platform: {errors}")


def main():
    parser = argparse.ArgumentParser(description="Run external-machine PPARG-DINP OpenMM production MD from a validated system.xml.")
    parser.add_argument("--system-xml", required=True, type=Path)
    parser.add_argument("--topology-pdb", required=True, type=Path)
    parser.add_argument("--initial-pdb", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--production-ns", type=float, default=100.0)
    parser.add_argument("--equilibration-ns", type=float, default=0.1)
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--temperature-k", type=float, default=310.0)
    parser.add_argument("--friction-per-ps", type=float, default=1.0)
    parser.add_argument("--report-ps", type=float, default=10.0)
    parser.add_argument("--checkpoint-ps", type=float, default=100.0)
    parser.add_argument("--platform", choices=["auto", "CUDA", "OpenCL", "CPU", "Reference"], default="auto")
    parser.add_argument("--seed", type=int, default=20260909)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    for path in (args.system_xml, args.topology_pdb, args.initial_pdb):
        if not path.exists():
            raise FileNotFoundError(path)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    import openmm
    from openmm import XmlSerializer, unit
    from openmm.app import CheckpointReporter, DCDReporter, PDBFile, Simulation, StateDataReporter

    system = XmlSerializer.deserialize(args.system_xml.read_text(encoding="utf-8"))
    topology_file = PDBFile(str(args.topology_pdb))
    initial_file = PDBFile(str(args.initial_pdb))
    n_system = system.getNumParticles()
    n_topology = topology_file.topology.getNumAtoms()
    n_initial = initial_file.topology.getNumAtoms()
    if not (n_system == n_topology == n_initial):
        raise RuntimeError(f"Atom count mismatch: system={n_system}, topology={n_topology}, initial={n_initial}")

    platform, platform_properties, selected_platform = choose_platform(openmm, args.platform)
    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "system_xml": str(args.system_xml),
        "topology_pdb": str(args.topology_pdb),
        "initial_pdb": str(args.initial_pdb),
        "output_dir": str(args.output_dir),
        "n_particles": n_system,
        "production_ns": args.production_ns,
        "equilibration_ns": args.equilibration_ns,
        "timestep_fs": args.timestep_fs,
        "temperature_k": args.temperature_k,
        "random_seed": args.seed,
        "platform_requested": args.platform,
        "platform_selected": selected_platform,
        "platform_properties": platform_properties,
        "interpretation_boundary": "MD stability supports structural plausibility only; it does not establish experimental affinity or causality.",
    }
    save_json(metadata, args.output_dir / "md_setup_manifest.json")
    if args.dry_run:
        print(json.dumps({"status": "dry_run_ok", **metadata}, indent=2, ensure_ascii=False))
        return

    integrator = openmm.LangevinMiddleIntegrator(
        args.temperature_k * unit.kelvin,
        args.friction_per_ps / unit.picosecond,
        args.timestep_fs * unit.femtoseconds,
    )
    simulation = Simulation(topology_file.topology, system, integrator, platform, platform_properties)
    checkpoint = args.output_dir / "production.chk"
    dcd = args.output_dir / "production.dcd"
    state_csv = args.output_dir / "production_state.csv"
    target_steps = int(round(args.production_ns * 1_000_000.0 / args.timestep_fs))
    equil_steps = int(round(args.equilibration_ns * 1_000_000.0 / args.timestep_fs))
    report_steps = max(1, int(round(args.report_ps * 1000.0 / args.timestep_fs)))
    checkpoint_steps = max(1, int(round(args.checkpoint_ps * 1000.0 / args.timestep_fs)))

    if args.resume and checkpoint.exists():
        simulation.loadCheckpoint(str(checkpoint))
    else:
        simulation.context.setPositions(initial_file.positions)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(args.temperature_k * unit.kelvin, args.seed)
        if equil_steps:
            simulation.step(equil_steps)

    simulation.reporters.clear()
    simulation.reporters.append(DCDReporter(str(dcd), report_steps, enforcePeriodicBox=True))
    simulation.reporters.append(StateDataReporter(str(state_csv), report_steps, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, progress=True, remainingTime=True, speed=True, totalSteps=target_steps, separator=","))
    simulation.reporters.append(CheckpointReporter(str(checkpoint), checkpoint_steps))
    remaining = max(0, target_steps - simulation.currentStep)
    simulation.step(remaining)
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    with (args.output_dir / "final.pdb").open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(topology_file.topology, state.getPositions(), handle)
    save_json({"status": "completed", "current_step": simulation.currentStep, "target_steps": target_steps, "platform": selected_platform}, args.output_dir / "md_completion.json")
    print(json.dumps({"status": "completed", "current_step": simulation.currentStep, "target_steps": target_steps, "output_dir": str(args.output_dir)}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
