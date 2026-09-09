from __future__ import annotations

import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import openmm


HERE = Path(__file__).resolve().parent


def available(name):
    return bool(importlib.util.find_spec(name))


def main():
    platforms = []
    for index in range(openmm.Platform.getNumPlatforms()):
        platform = openmm.Platform.getPlatform(index)
        item = {"name": platform.getName(), "properties": platform.getPropertyNames()}
        if platform.getName() in {"OpenCL", "CUDA"}:
            try:
                system = openmm.System()
                system.addParticle(1.0 * openmm.unit.dalton)
                integrator = openmm.VerletIntegrator(0.001 * openmm.unit.picoseconds)
                context = openmm.Context(system, integrator, platform)
                item["device_name"] = platform.getPropertyValue(context, "DeviceName")
                item["platform_name"] = platform.getPropertyValue(context, "OpenCLPlatformName") if platform.getName() == "OpenCL" else ""
                del context
            except Exception as exc:
                item["probe_error"] = repr(exc)
        platforms.append(item)
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "openmm_version": openmm.version.version,
        "platforms": platforms,
        "python_packages": {
            "openmm": available("openmm"),
            "pdbfixer": available("pdbfixer"),
            "rdkit": available("rdkit"),
            "openff_toolkit": available("openff"),
            "openmmforcefields": available("openmmforcefields"),
            "mdtraj": available("mdtraj"),
        },
        "md_readiness": "not_ready_for_production",
        "reason": "Protein force field is available through OpenMM, but a validated ligand force-field/charge pipeline for DINP is not installed; PPARG docking pose consensus also requires inspection before MD.",
    }
    (HERE / "DINP_CRC_PPARG_MD_runtime_probe.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
