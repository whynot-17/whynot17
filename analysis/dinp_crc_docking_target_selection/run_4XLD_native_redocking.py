from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
PREP = HERE / "prepared"
OUT = HERE / "native_redocking_4XLD"
PYTHON = Path(__file__).resolve().parents[2] / "work" / ".venv_docking" / "Scripts" / "python.exe"
VINA = Path(__file__).resolve().parents[2] / "work" / "bin" / "vina_1.2.7_win.exe"
PDB = INPUTS / "4XLD.pdb"
IDEAL_SDF = OUT / "BRL_ideal.sdf"
NATIVE_SDF = OUT / "4XLD_native_BRL.sdf"
NATIVE_PDBQT = OUT / "4XLD_native_BRL.pdbqt"
RECEPTOR = PREP / "PPARG_4XLD_rigid.pdbqt"
BOX = PREP / "PPARG_4XLD_box.txt"
REDOCKED_PDBQT = OUT / "4XLD_native_redocked_BRL.pdbqt"
LOG = OUT / "4XLD_native_redocking.log"

SEED = 20260909
EXHAUSTIVENESS = 32
NUM_MODES = 20
ENERGY_RANGE = 5
CPU = 8
RMSD_PASS_A = 2.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def native_brl_atoms(path: Path):
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip() != "BRL":
            continue
        if line[21].strip() != "A" or line[22:26].strip() != "502":
            continue
        element = line[76:78].strip() or re.sub(r"[^A-Za-z]", "", line[12:16]).upper()[:1]
        atoms.append(
            {
                "name": line[12:16].strip(),
                "element": element.upper(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            }
        )
    if len(atoms) != 25:
        raise AssertionError(f"Expected 25 BRL heavy atoms in 4XLD, found {len(atoms)}")
    if any(atom["element"] == "H" for atom in atoms):
        raise AssertionError("Native BRL extraction unexpectedly contains hydrogen")
    return atoms


def write_native_sdf(native_atoms):
    supplier = Chem.SDMolSupplier(str(IDEAL_SDF), removeHs=False, sanitize=True)
    mol = supplier[0] if supplier and supplier[0] is not None else None
    if mol is None:
        raise ValueError(f"Could not parse {IDEAL_SDF}")
    heavy = [atom for atom in mol.GetAtoms() if atom.GetSymbol().upper() != "H"]
    if len(heavy) != len(native_atoms):
        raise AssertionError(f"CCD heavy atom count {len(heavy)} != native count {len(native_atoms)}")
    native_elements = [atom["element"] for atom in native_atoms]
    sdf_elements = [atom.GetSymbol().upper() for atom in heavy]
    if native_elements != sdf_elements:
        raise AssertionError(f"Native/CCD heavy atom order mismatch: {native_elements} != {sdf_elements}")
    conformer = mol.GetConformer()
    for atom, native in zip(heavy, native_atoms):
        conformer.SetAtomPosition(atom.GetIdx(), (native["x"], native["y"], native["z"]))
    writer = Chem.SDWriter(str(NATIVE_SDF))
    writer.write(mol)
    writer.close()
    return mol


def run_command(command, log_path: Path):
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    text = (completed.stdout or "") + "\n" + (completed.stderr or "")
    log_path.write_text(text, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(map(str, command))}\n{text}")
    return text


def parse_vina_modes(text: str):
    rows = []
    pattern = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$")
    in_table = False
    for line in text.splitlines():
        if "-----+------------+----------+----------" in line:
            in_table = True
            continue
        if not in_table:
            continue
        match = pattern.match(line)
        if match:
            rows.append(
                {
                    "mode": int(match.group(1)),
                    "affinity_kcal_mol": float(match.group(2)),
                    "vina_rmsd_lb_A": float(match.group(3)),
                    "vina_rmsd_ub_A": float(match.group(4)),
                }
            )
        elif rows and line.strip() and not re.match(r"^\s*\d+", line):
            break
    if not rows:
        raise ValueError("Could not parse Vina mode table from log")
    return rows


def parse_pdbqt_models(path: Path):
    models = {}
    current = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            current = int(line[5:].strip())
            models[current] = []
        elif line.startswith(("ATOM", "HETATM")):
            if current is None:
                current = 1
                models.setdefault(current, [])
            atom_name = line[12:16].strip()
            atom_type = line[77:79].strip().upper()
            element = re.sub(r"[^A-Z]", "", atom_type)[:1]
            if not element:
                element = re.sub(r"[^A-Z]", "", atom_name.upper())[:1]
            models[current].append(
                {
                    "serial": int(line[6:11]),
                    "name": atom_name,
                    "element": element,
                    "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                }
            )
    if not models:
        raise ValueError(f"No atom models parsed from {path}")
    return models


def parse_index_map(path: Path):
    """Return PDBQT atom serial -> input SDF atom index (both 1-based)."""
    mapping = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("REMARK INDEX MAP"):
            continue
        values = line.split()[3:]
        if len(values) % 2:
            raise AssertionError(f"Malformed INDEX MAP line: {line}")
        for i in range(0, len(values), 2):
            input_idx = int(values[i])
            pdbqt_idx = int(values[i + 1])
            mapping[pdbqt_idx] = input_idx
    if not mapping:
        raise AssertionError(f"No REMARK INDEX MAP found in {path}")
    return mapping


def rmsd(native_xyz, dock_xyz):
    native = np.asarray(native_xyz, dtype=float)
    docked = np.asarray(dock_xyz, dtype=float)
    if native.shape != docked.shape:
        raise AssertionError(f"RMSD coordinate shape mismatch: {native.shape} != {docked.shape}")
    direct = float(np.sqrt(np.mean(np.sum((docked - native) ** 2, axis=1))))
    native_centered = native - native.mean(axis=0)
    docked_centered = docked - docked.mean(axis=0)
    covariance = docked_centered.T @ native_centered
    u, _, vh = np.linalg.svd(covariance)
    rotation = u @ vh
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vh
    fitted = docked_centered @ rotation
    fitted_rmsd = float(np.sqrt(np.mean(np.sum((fitted - native_centered) ** 2, axis=1))))
    return direct, fitted_rmsd


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for required in (PDB, IDEAL_SDF, RECEPTOR, BOX, PYTHON, VINA):
        if not required.exists():
            raise FileNotFoundError(required)

    native_atoms = native_brl_atoms(PDB)
    write_native_sdf(native_atoms)

    prepare_command = [
        str(PYTHON),
        "-m",
        "meeko.cli.mk_prepare_ligand",
        "-i",
        str(NATIVE_SDF),
        "-o",
        str(NATIVE_PDBQT),
        "--add_index_map",
    ]
    preparation_log = OUT / "meeko_prepare_native_BRL.log"
    run_command(prepare_command, preparation_log)

    vina_command = [
        str(VINA),
        "--config",
        str(BOX),
        "--receptor",
        str(RECEPTOR),
        "--ligand",
        str(NATIVE_PDBQT),
        "--out",
        str(REDOCKED_PDBQT),
        "--exhaustiveness",
        str(EXHAUSTIVENESS),
        "--num_modes",
        str(NUM_MODES),
        "--energy_range",
        str(ENERGY_RANGE),
        "--cpu",
        str(CPU),
        "--seed",
        str(SEED),
        "--verbosity",
        "1",
    ]
    vina_text = run_command(vina_command, LOG)
    mode_rows = parse_vina_modes(vina_text)

    native_xyz = [[atom["x"], atom["y"], atom["z"]] for atom in native_atoms]
    native_elements = [atom["element"] for atom in native_atoms]
    index_map = parse_index_map(NATIVE_PDBQT)
    models = parse_pdbqt_models(REDOCKED_PDBQT)
    rmsd_rows = []
    for row in mode_rows:
        model_atoms = models.get(row["mode"], [])
        mapped_atoms = []
        for atom in model_atoms:
            input_idx = index_map.get(atom["serial"])
            if input_idx is None:
                raise AssertionError(f"Mode {row['mode']} atom serial {atom['serial']} missing from INDEX MAP")
            if input_idx <= len(native_atoms):
                mapped_atoms.append((input_idx, atom))
        mapped_atoms.sort(key=lambda pair: pair[0])
        model_elements = [native_elements[input_idx - 1] for input_idx, _ in mapped_atoms]
        if model_elements != native_elements:
            raise AssertionError(
                f"Mode {row['mode']} atom mapping/elements do not match native BRL: {model_elements} != {native_elements}"
            )
        docked_xyz = [atom["xyz"] for _, atom in mapped_atoms]
        direct, fitted = rmsd(native_xyz, docked_xyz)
        rmsd_rows.append(
            {
                **row,
                "receptor_frame_heavy_atom_rmsd_A": direct,
                "kabsch_fit_heavy_atom_rmsd_A": fitted,
                "top_ranked_pass_2A": "yes" if row["mode"] == 1 and direct <= RMSD_PASS_A else "no",
                "pass_2A": "yes" if direct <= RMSD_PASS_A else "no",
            }
        )

    rmsd_csv = OUT / "4XLD_native_redocking_rmsd.csv"
    fields = list(rmsd_rows[0].keys())
    with rmsd_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rmsd_rows)

    top = rmsd_rows[0]
    best_rmsd = min(rmsd_rows, key=lambda row: row["receptor_frame_heavy_atom_rmsd_A"])
    top_pass = top["receptor_frame_heavy_atom_rmsd_A"] <= RMSD_PASS_A
    report = f"""# 4XLD native-ligand self-redocking QC

## Scope

This is a protocol QC run using the co-crystallized 4XLD ligand **BRL** (rosiglitazone), not DINP. The native BRL heavy-atom coordinates were extracted from 4XLD and combined with the RCSB Chemical Component Dictionary topology before ligand preparation.

## Protocol

- Receptor: `prepared/PPARG_4XLD_rigid.pdbqt`
- Native ligand: BRL, chain A, residue 502 from `inputs/4XLD.pdb`
- Box: `prepared/PPARG_4XLD_box.txt`
- Engine: AutoDock Vina 1.2.7
- Seed: `{SEED}`
- Exhaustiveness: `{EXHAUSTIVENESS}`
- CPU: `{CPU}`
- Number of modes: `{NUM_MODES}`
- Energy range: `{ENERGY_RANGE}` kcal/mol
- Pre-specified receptor-frame heavy-atom RMSD pass threshold: **≤ {RMSD_PASS_A:.1f} Å**

## Result

| metric | value |
|---|---:|
| top-ranked affinity | {top['affinity_kcal_mol']:.3f} kcal/mol |
| top-ranked receptor-frame heavy-atom RMSD | {top['receptor_frame_heavy_atom_rmsd_A']:.3f} Å |
| top-ranked Kabsch-fit heavy-atom RMSD | {top['kabsch_fit_heavy_atom_rmsd_A']:.3f} Å |
| top-ranked pose QC | {'PASS' if top_pass else 'FAIL'} |
| lowest-RMSD mode | {best_rmsd['mode']} |
| lowest receptor-frame heavy-atom RMSD | {best_rmsd['receptor_frame_heavy_atom_rmsd_A']:.3f} Å |
| lowest-RMSD pose affinity | {best_rmsd['affinity_kcal_mol']:.3f} kcal/mol |

The top-ranked pose is considered protocol-valid only if its receptor-frame heavy-atom RMSD is ≤2.0 Å. A low RMSD in a lower-ranked mode is reported as supportive but does not replace the top-ranked criterion.

## Files

- `4XLD_native_BRL.pdb`: extracted native ligand coordinates
- `4XLD_native_BRL.sdf`: native-coordinate ligand with CCD topology
- `4XLD_native_BRL.pdbqt`: prepared docking ligand
- `4XLD_native_redocked_BRL.pdbqt`: Vina poses
- `4XLD_native_redocking_rmsd.csv`: affinity and RMSD for each returned mode
- `4XLD_native_redocking.log`: complete Vina output
- `run_4XLD_native_redocking.py`: reproducible script

This QC does not perform molecular dynamics and does not validate CEBPB as a docking target. It only tests the PPARG/4XLD docking protocol used for the subsequent PPARG-DINP MD starting pose.
"""
    (OUT / "4XLD_native_redocking_QC_report.md").write_text(report, encoding="utf-8")

    native_pdb = OUT / "4XLD_native_BRL.pdb"
    native_lines = [
        line
        for line in PDB.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("HETATM") and line[17:20].strip() == "BRL" and line[21].strip() == "A" and line[22:26].strip() == "502"
    ]
    native_pdb.write_text("\n".join(native_lines) + "\nEND\n", encoding="utf-8")

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "assay": "4XLD native-ligand self-redocking",
        "ligand": {"resname": "BRL", "description": "co-crystallized rosiglitazone", "chain": "A", "resseq": 502},
        "inputs": {
            "pdb": {"path": str(PDB), "sha256": sha256(PDB)},
            "ideal_sdf": {"path": str(IDEAL_SDF), "sha256": sha256(IDEAL_SDF)},
            "receptor_pdbqt": {"path": str(RECEPTOR), "sha256": sha256(RECEPTOR)},
            "box": {"path": str(BOX), "sha256": sha256(BOX)},
        },
        "parameters": {
            "engine": "AutoDock Vina 1.2.7",
            "seed": SEED,
            "exhaustiveness": EXHAUSTIVENESS,
            "num_modes": NUM_MODES,
            "energy_range_kcal_mol": ENERGY_RANGE,
            "cpu": CPU,
            "rmsd_pass_threshold_A": RMSD_PASS_A,
        },
        "result": {
            "top_affinity_kcal_mol": top["affinity_kcal_mol"],
            "top_receptor_frame_rmsd_A": top["receptor_frame_heavy_atom_rmsd_A"],
            "top_ranked_pass": bool(top_pass),
            "lowest_rmsd_mode": best_rmsd["mode"],
            "lowest_rmsd_A": best_rmsd["receptor_frame_heavy_atom_rmsd_A"],
        },
        "outputs": {
            "native_pdb": str(native_pdb),
            "native_sdf": str(NATIVE_SDF),
            "native_pdbqt": str(NATIVE_PDBQT),
            "redocked_pdbqt": str(REDOCKED_PDBQT),
            "rmsd_csv": str(rmsd_csv),
            "log": str(LOG),
            "report": str(OUT / "4XLD_native_redocking_QC_report.md"),
        },
    }
    (OUT / "4XLD_native_redocking_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest["result"], indent=2))


if __name__ == "__main__":
    main()
