from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from run_4XLD_native_redocking import (
    BOX as ORIGINAL_BOX,
    CPU,
    ENERGY_RANGE,
    EXHAUSTIVENESS,
    IDEAL_SDF,
    LOG as ORIGINAL_LOG,
    NATIVE_PDBQT,
    NATIVE_SDF,
    NUM_MODES,
    PDB,
    PREP,
    PYTHON,
    RECEPTOR,
    RMSD_PASS_A,
    SEED,
    VINA,
    native_brl_atoms,
    parse_index_map,
    parse_pdbqt_models,
    parse_vina_modes,
    rmsd,
)


HERE = Path(__file__).resolve().parent
OUT = HERE / "native_redocking_4XLD_focused_exh128"
FOCUSED_BOX = OUT / "PPARG_4XLD_focused_native_box.txt"
REDOCKED_PDBQT = OUT / "4XLD_native_redocked_BRL_focused_exh128.pdbqt"
LOG = OUT / "4XLD_native_redocking_focused_exh128.log"
RMSD_CSV = OUT / "4XLD_native_redocking_focused_exh128_rmsd.csv"
REPORT = OUT / "4XLD_native_redocking_focused_exh128_QC_report.md"
MANIFEST = OUT / "4XLD_native_redocking_focused_exh128_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_focused_box(native_atoms):
    xyz = np.asarray([[a["x"], a["y"], a["z"]] for a in native_atoms], dtype=float)
    lower = xyz.min(axis=0)
    upper = xyz.max(axis=0)
    center = (lower + upper) / 2.0
    span = upper - lower
    # Native-ligand span plus an 8-A total margin, rounded upward to practical values.
    size = np.ceil(span + 8.0).astype(int)
    size = np.maximum(size, np.array([14, 20, 16]))
    FOCUSED_BOX.write_text(
        f"center_x = {center[0]:.4f}\n"
        f"center_y = {center[1]:.4f}\n"
        f"center_z = {center[2]:.4f}\n"
        f"size_x = {size[0]:.4f}\n"
        f"size_y = {size[1]:.4f}\n"
        f"size_z = {size[2]:.4f}\n",
        encoding="utf-8",
    )
    return center, span, size


def run_vina():
    command = [
        str(VINA),
        "--config", str(FOCUSED_BOX),
        "--receptor", str(RECEPTOR),
        "--ligand", str(NATIVE_PDBQT),
        "--out", str(REDOCKED_PDBQT),
        "--exhaustiveness", "128",
        "--num_modes", str(NUM_MODES),
        "--energy_range", str(ENERGY_RANGE),
        "--cpu", str(CPU),
        "--seed", str(SEED),
        "--verbosity", "1",
    ]
    completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace")
    text = (completed.stdout or "") + "\n" + (completed.stderr or "")
    LOG.write_text(text, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"Vina failed with exit code {completed.returncode}\n{text}")
    return command, text


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for required in (PDB, IDEAL_SDF, NATIVE_SDF, NATIVE_PDBQT, RECEPTOR, ORIGINAL_BOX, PYTHON, VINA):
        if not required.exists():
            raise FileNotFoundError(required)

    native_atoms = native_brl_atoms(PDB)
    center, span, size = write_focused_box(native_atoms)
    command, vina_text = run_vina()
    modes = parse_vina_modes(vina_text)
    models = parse_pdbqt_models(REDOCKED_PDBQT)
    index_map = parse_index_map(REDOCKED_PDBQT)
    native_xyz = [[a["x"], a["y"], a["z"]] for a in native_atoms]
    native_elements = [a["element"] for a in native_atoms]

    rows = []
    for mode in modes:
        mapped = []
        for atom in models[mode["mode"]]:
            input_idx = index_map.get(atom["serial"])
            if input_idx is None:
                raise AssertionError(f"Missing index map for PDBQT serial {atom['serial']}")
            if input_idx <= len(native_atoms):
                mapped.append((input_idx, atom))
        mapped.sort(key=lambda pair: pair[0])
        if [native_elements[idx - 1] for idx, _ in mapped] != native_elements:
            raise AssertionError(f"Mode {mode['mode']} does not map to the 25 native BRL heavy atoms")
        direct, fitted = rmsd(native_xyz, [atom["xyz"] for _, atom in mapped])
        rows.append({
            **mode,
            "receptor_frame_heavy_atom_rmsd_A": direct,
            "kabsch_fit_heavy_atom_rmsd_A": fitted,
            "top3": "yes" if mode["mode"] <= 3 else "no",
            "pass_2A": "yes" if direct <= RMSD_PASS_A else "no",
        })

    with RMSD_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    top = rows[0]
    top3 = rows[:3]
    native_like = [row for row in rows if row["receptor_frame_heavy_atom_rmsd_A"] <= RMSD_PASS_A]
    best = min(rows, key=lambda row: row["receptor_frame_heavy_atom_rmsd_A"])
    top3_native_like = [row for row in top3 if row["receptor_frame_heavy_atom_rmsd_A"] <= RMSD_PASS_A]
    report = f"""# 4XLD focused-box native redocking QC (exhaustiveness 128)

## Purpose

This is a focused follow-up to the original 4XLD BRL self-redocking QC. It keeps the same PPARG receptor, BRL native-coordinate ligand, Vina seed, and scoring protocol, but uses a ligand-centered pocket box and 4× higher exhaustiveness. The aim is to test whether the native-like pose can move into the top-ranked/top-three poses.

## Box definition

The box is centered on the native BRL heavy-atom bounding-box midpoint. Each axis uses the native ligand span plus an 8 Å total margin, rounded upward to practical dimensions. This gives a focused box of `{size[0]}` × `{size[1]}` × `{size[2]}` Å, compared with the original 24 × 24 × 24 Å box.

- center: ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f}) Å
- native ligand span: ({span[0]:.3f}, {span[1]:.3f}, {span[2]:.3f}) Å

## Protocol

- Receptor: `../prepared/PPARG_4XLD_rigid.pdbqt`
- Native ligand: BRL/rosiglitazone, chain A residue 502 from `../inputs/4XLD.pdb`
- Engine: AutoDock Vina 1.2.7
- Seed: `{SEED}`
- Exhaustiveness: **128** (original: 32)
- CPU: `{CPU}`
- Modes: `{NUM_MODES}`
- Energy range: `{ENERGY_RANGE}` kcal/mol
- Pre-specified native-redocking pass threshold: receptor-frame heavy-atom RMSD ≤ `{RMSD_PASS_A:.1f}` Å

## Result

| metric | value |
|---|---:|
| top-ranked affinity | {top['affinity_kcal_mol']:.3f} kcal/mol |
| top-ranked receptor-frame RMSD | {top['receptor_frame_heavy_atom_rmsd_A']:.3f} Å |
| top-ranked Kabsch-fit RMSD | {top['kabsch_fit_heavy_atom_rmsd_A']:.3f} Å |
| native-like pose in top1 | {'yes' if top['receptor_frame_heavy_atom_rmsd_A'] <= 2.0 else 'no'} |
| native-like pose in top3 | {'yes' if top3_native_like else 'no'} |
| lowest-RMSD mode | {best['mode']} |
| lowest receptor-frame RMSD | {best['receptor_frame_heavy_atom_rmsd_A']:.3f} Å |
| lowest-RMSD pose affinity | {best['affinity_kcal_mol']:.3f} kcal/mol |
| number of poses with RMSD ≤2 Å | {len(native_like)} |

The focused/high-exhaustiveness run is considered successful for ranking if a pose with RMSD ≤2 Å reaches top1 or top3. The result is reported separately from the original 24 Å-box run and does not alter the original DINP docking results.

## Files

- `PPARG_4XLD_focused_native_box.txt`: focused box parameters
- `4XLD_native_redocked_BRL_focused_exh128.pdbqt`: focused/high-exhaustiveness poses
- `4XLD_native_redocking_focused_exh128_rmsd.csv`: affinity and RMSD for all modes
- `4XLD_native_redocking_focused_exh128.log`: complete Vina output
- `../run_4XLD_native_redocking_focused.py`: reproducible script
"""
    REPORT.write_text(report, encoding="utf-8")
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "assay": "4XLD focused-box native-ligand self-redocking",
        "ligand": {"resname": "BRL", "description": "co-crystallized rosiglitazone", "chain": "A", "resseq": 502},
        "box": {
            "derivation": "native BRL heavy-atom bounding-box midpoint; native span plus 8 A total margin; rounded upward",
            "center_A": [float(x) for x in center],
            "native_span_A": [float(x) for x in span],
            "size_A": [int(x) for x in size],
        },
        "parameters": {
            "engine": "AutoDock Vina 1.2.7",
            "seed": SEED,
            "exhaustiveness": 128,
            "original_exhaustiveness": EXHAUSTIVENESS,
            "num_modes": NUM_MODES,
            "energy_range_kcal_mol": ENERGY_RANGE,
            "cpu": CPU,
            "rmsd_pass_threshold_A": RMSD_PASS_A,
        },
        "result": {
            "top_affinity_kcal_mol": top["affinity_kcal_mol"],
            "top_rmsd_A": top["receptor_frame_heavy_atom_rmsd_A"],
            "native_like_in_top1": bool(top["receptor_frame_heavy_atom_rmsd_A"] <= RMSD_PASS_A),
            "native_like_in_top3": bool(top3_native_like),
            "best_mode": best["mode"],
            "best_rmsd_A": best["receptor_frame_heavy_atom_rmsd_A"],
            "native_like_pose_count": len(native_like),
        },
        "provenance": {
            "receptor_sha256": sha256(RECEPTOR),
            "native_ligand_pdbqt_sha256": sha256(NATIVE_PDBQT),
            "original_box_sha256": sha256(ORIGINAL_BOX),
            "command": [str(x) for x in command],
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest["result"], indent=2))


if __name__ == "__main__":
    main()
