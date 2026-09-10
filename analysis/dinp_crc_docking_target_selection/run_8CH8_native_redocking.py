from __future__ import annotations

import csv
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem

from run_4XLD_native_redocking import (
    CPU,
    ENERGY_RANGE,
    EXHAUSTIVENESS,
    PYTHON,
    RMSD_PASS_A,
    SEED,
    VINA,
    parse_index_map,
    parse_pdbqt_models,
    parse_vina_modes,
    rmsd,
    sha256,
)


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
PREP = HERE / "prepared"
OUT = HERE / "native_redocking_8CH8_ULC"
PDB = INPUTS / "8CH8.pdb"
IDEAL_SDF = OUT / "ULC_ideal.sdf"
NATIVE_SDF = OUT / "8CH8_native_ULC.sdf"
NATIVE_PDBQT = OUT / "8CH8_native_ULC.pdbqt"
RECEPTOR = PREP / "NR1I2_8CH8_rigid.pdbqt"
BOX = PREP / "NR1I2_8CH8_box.txt"
REDOCKED_PDBQT = OUT / "8CH8_native_redocked_ULC.pdbqt"
LOG = OUT / "8CH8_native_redocking.log"

NUM_MODES = 20


def native_ulc_atoms(path: Path):
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip() != "ULC":
            continue
        if line[21].strip() != "A" or line[22:26].strip() != "501":
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
    if len(atoms) != 23:
        raise AssertionError(f"Expected 23 ULC heavy atoms in 8CH8, found {len(atoms)}")
    if any(atom["element"] == "H" for atom in atoms):
        raise AssertionError("Native ULC extraction unexpectedly contains hydrogen")
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
        raise AssertionError(f"Native/CCD ULC heavy atom order mismatch: {native_elements} != {sdf_elements}")
    conformer = mol.GetConformer()
    for atom, native in zip(heavy, native_atoms):
        conformer.SetAtomPosition(atom.GetIdx(), (native["x"], native["y"], native["z"]))
    writer = Chem.SDWriter(str(NATIVE_SDF))
    writer.write(mol)
    writer.close()


def run_command(command, log_path: Path):
    completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace")
    text = (completed.stdout or "") + "\n" + (completed.stderr or "")
    log_path.write_text(text, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(map(str, command))}\n{text}")
    return text


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for required in (PDB, IDEAL_SDF, RECEPTOR, BOX, PYTHON, VINA):
        if not required.exists():
            raise FileNotFoundError(required)

    native_atoms = native_ulc_atoms(PDB)
    write_native_sdf(native_atoms)

    preparation_command = [
        str(PYTHON), "-m", "meeko.cli.mk_prepare_ligand",
        "-i", str(NATIVE_SDF), "-o", str(NATIVE_PDBQT), "--add_index_map",
    ]
    run_command(preparation_command, OUT / "meeko_prepare_native_ULC.log")

    vina_command = [
        str(VINA), "--config", str(BOX), "--receptor", str(RECEPTOR),
        "--ligand", str(NATIVE_PDBQT), "--out", str(REDOCKED_PDBQT),
        "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
        "--energy_range", str(ENERGY_RANGE), "--cpu", str(CPU),
        "--seed", str(SEED), "--verbosity", "1",
    ]
    vina_text = run_command(vina_command, LOG)
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
                raise AssertionError(f"Mode {mode['mode']} atom serial {atom['serial']} missing from INDEX MAP")
            if input_idx <= len(native_atoms):
                mapped.append((input_idx, atom))
        mapped.sort(key=lambda pair: pair[0])
        if [native_elements[idx - 1] for idx, _ in mapped] != native_elements:
            raise AssertionError(f"Mode {mode['mode']} does not map to the 23 native ULC heavy atoms")
        direct, fitted = rmsd(native_xyz, [atom["xyz"] for _, atom in mapped])
        rows.append({
            **mode,
            "receptor_frame_heavy_atom_rmsd_A": direct,
            "kabsch_fit_heavy_atom_rmsd_A": fitted,
            "top_ranked_pass_2A": "yes" if mode["mode"] == 1 and direct <= RMSD_PASS_A else "no",
            "pass_2A": "yes" if direct <= RMSD_PASS_A else "no",
        })

    rmsd_csv = OUT / "8CH8_native_redocking_rmsd.csv"
    with rmsd_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    native_pdb = OUT / "8CH8_native_ULC.pdb"
    native_lines = [
        line for line in PDB.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("HETATM") and line[17:20].strip() == "ULC" and line[21].strip() == "A" and line[22:26].strip() == "501"
    ]
    native_pdb.write_text("\n".join(native_lines) + "\nEND\n", encoding="utf-8")

    top = rows[0]
    best = min(rows, key=lambda row: row["receptor_frame_heavy_atom_rmsd_A"])
    native_like = [row for row in rows if row["receptor_frame_heavy_atom_rmsd_A"] <= RMSD_PASS_A]
    report = f"""# 8CH8 native-ligand self-redocking QC

## Scope

This QC uses the co-crystallized 8CH8 ligand **ULC** (the PXR/NR1I2 ligand in chain A, residue 501), not DINP. Native ULC heavy-atom coordinates were extracted from `inputs/8CH8.pdb` and combined with the RCSB Chemical Component Dictionary topology before Meeko preparation.

## Protocol

- Receptor: `prepared/NR1I2_8CH8_rigid.pdbqt`
- Native ligand: ULC, chain A, residue 501 from `inputs/8CH8.pdb`
- Box: `prepared/NR1I2_8CH8_box.txt` (24 × 24 × 24 Å)
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
| top-ranked pose QC | {'PASS' if top['receptor_frame_heavy_atom_rmsd_A'] <= RMSD_PASS_A else 'FAIL'} |
| lowest-RMSD mode | {best['mode']} |
| lowest receptor-frame heavy-atom RMSD | {best['receptor_frame_heavy_atom_rmsd_A']:.3f} Å |
| lowest-RMSD pose affinity | {best['affinity_kcal_mol']:.3f} kcal/mol |
| number of poses with RMSD ≤2 Å | {len(native_like)} |

The top-ranked pose is considered protocol-valid only if its receptor-frame heavy-atom RMSD is ≤2.0 Å. A low RMSD in a lower-ranked mode is reported as supportive but does not replace the top-ranked criterion.

## Files

- `8CH8_native_ULC.pdb`: extracted native ligand coordinates
- `8CH8_native_ULC.sdf`: native-coordinate ligand with CCD topology
- `8CH8_native_ULC.pdbqt`: prepared docking ligand
- `8CH8_native_redocked_ULC.pdbqt`: Vina poses
- `8CH8_native_redocking_rmsd.csv`: affinity and RMSD for each returned mode
- `8CH8_native_redocking.log`: complete Vina output
- `../run_8CH8_native_redocking.py`: reproducible script

This QC does not perform molecular dynamics and does not validate CEBPB as a docking target. It only tests the PXR/8CH8 docking protocol used for PXR-DINP structural interpretation.
"""
    report_path = OUT / "8CH8_native_redocking_QC_report.md"
    report_path.write_text(report, encoding="utf-8")

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "assay": "8CH8 native-ligand self-redocking",
        "target": {"gene_symbol": "NR1I2", "protein_name": "pregnane X receptor", "pdb_id": "8CH8"},
        "ligand": {"resname": "ULC", "chain": "A", "resseq": 501, "heavy_atom_count": 23},
        "parameters": {
            "engine": "AutoDock Vina 1.2.7", "seed": SEED,
            "exhaustiveness": EXHAUSTIVENESS, "num_modes": NUM_MODES,
            "energy_range_kcal_mol": ENERGY_RANGE, "cpu": CPU,
            "rmsd_pass_threshold_A": RMSD_PASS_A,
        },
        "result": {
            "top_affinity_kcal_mol": top["affinity_kcal_mol"],
            "top_receptor_frame_rmsd_A": top["receptor_frame_heavy_atom_rmsd_A"],
            "top_ranked_pass": bool(top["receptor_frame_heavy_atom_rmsd_A"] <= RMSD_PASS_A),
            "lowest_rmsd_mode": best["mode"],
            "lowest_rmsd_A": best["receptor_frame_heavy_atom_rmsd_A"],
            "native_like_pose_count": len(native_like),
        },
        "provenance": {
            "pdb_sha256": sha256(PDB), "ideal_sdf_sha256": sha256(IDEAL_SDF),
            "receptor_sha256": sha256(RECEPTOR), "box_sha256": sha256(BOX),
            "command": [str(x) for x in vina_command],
        },
    }
    (OUT / "8CH8_native_redocking_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest["result"], indent=2))


if __name__ == "__main__":
    main()
