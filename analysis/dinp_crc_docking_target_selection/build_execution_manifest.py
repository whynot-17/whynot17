from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
EXCLUDED_PACKAGE_FILES = {
    "inputs/6NX1.pdb",  # superseded PXR structure retained only as local prep history
    "prepared/NR1I2_6NX1_protein_h.pdb",  # superseded local prep artifact
}


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    summary_path = HERE / "docking_runs" / "DINP_CRC_docking_summary.csv"
    with summary_path.open("r", encoding="utf-8-sig", newline="") as handle:
        summary = list(csv.DictReader(handle))
    package_files = []
    for path in sorted(HERE.rglob("*")):
        relative_path = str(path.relative_to(HERE)).replace("\\", "/")
        if not path.is_file() or path.name.startswith("test_") or path.suffix == ".pyc" or "__pycache__" in path.parts or ".venv" in str(path) or relative_path in EXCLUDED_PACKAGE_FILES:
            continue
        package_files.append({"path": relative_path, "sha256": sha256(path), "bytes": path.stat().st_size})
    manifest = {
        "analysis": "DINP_CRC_docking_target_selection_and_vina_screen",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "primary_target": "PPARG",
        "support_targets": ["PPARA", "RXRA", "NR1I2", "ESR1"],
        "downstream_mediator_excluded": ["CEBPB"],
        "docking_engine": "AutoDock Vina 1.2.7",
        "docking_settings": {"exhaustiveness": 32, "num_modes": 20, "energy_range_kcal_mol": 5, "cpu": 8, "pparg_seeds": [20260909, 20260910, 20260911]},
        "docking_summary": summary,
        "md_status": "deferred_to_external_machine_by_design",
        "md_note": "MD is intentionally outside the current delivery scope; all docking inputs and poses are included for transfer to the user's external machine.",
        "package_files": package_files,
    }
    (HERE / "DINP_CRC_docking_execution_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
