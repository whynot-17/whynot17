#!/usr/bin/env python3
"""CXCR4 control-docking validation for the DINP structural pipeline.

Purpose
-------
Validate the CXCR4 docking protocol on PDB 3ODU by using the co-crystallized
small-molecule ligand IT1t (PDB residue ID ITD) as a positive-control ligand.
The script performs two complementary analyses under the SAME receptor/pocket
settings used for DINP:

1) IT1t redocking into the 3ODU pocket.
2) Side-by-side Vina affinity comparison: DINP vs IT1t.

Optional QC:
- Heavy-atom RMSD between the best redocked IT1t pose and the crystallographic
  IT1t pose after graph-based atom mapping with RDKit.

Important interpretation boundary
---------------------------------
A good redocking result validates the docking protocol for this structure and
pocket; it does NOT prove that DINP binds CXCR4 in vivo. Vina scores are
comparative computational estimates, not experimental binding free energies.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
except Exception as exc:  # pragma: no cover
    raise SystemExit("RDKit is required") from exc

ROOT = Path(__file__).resolve().parent
STAGE2 = ROOT / "outputs" / "stage2_docking"
OUT = ROOT / "outputs" / "cxcr4_control_docking"
RUNTIME = Path(r"D:\whynot17\dinp_stage2_runtime")

TARGET = "CXCR4"
PDB_ID = "3ODU"
POCKET_RESNAME = "ITD"


def run(cmd: List[str], cwd: Optional[Path] = None, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        capture_output=True,
        check=check,
    )


def which_any(names: List[str]) -> Optional[str]:
    for name in names:
        p = shutil.which(name)
        if p:
            return p
    return None


def prepend_runtime_to_path(runtime: Path) -> None:
    candidates = [
        runtime,
        runtime / "Scripts",
        runtime / "Library" / "bin",
        runtime / "bin",
    ]
    existing = os.environ.get("PATH", "")
    parts = [str(p) for p in candidates if p.exists()]
    if parts:
        os.environ["PATH"] = os.pathsep.join(parts + [existing])


def load_json(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def save_json(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, ensure_ascii=False)


def parse_vina_best(stdout: str) -> Optional[float]:
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0].isdigit():
            try:
                return float(fields[1])
            except ValueError:
                continue
    return None


def extract_itd_from_pdb(pdb_path: Path, out_pdb: Path, chain: Optional[str] = None) -> int:
    """Extract crystallographic ITD HETATM records from 3ODU."""
    lines = []
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        ch = line[21].strip()
        if resname != POCKET_RESNAME:
            continue
        if chain and ch != chain:
            continue
        lines.append(line)
    if not lines:
        raise RuntimeError(f"No {POCKET_RESNAME} HETATM records found in {pdb_path}")
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return len(lines)


def ligand_center_and_box_from_pdb(pdb_path: Path, padding: float = 12.0, min_size: float = 22.0, max_size: float = 34.0) -> Dict:
    coords = []
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords.append((x, y, z))
    if not coords:
        raise RuntimeError("No ligand coordinates available for pocket box")
    center = [sum(p[i] for p in coords) / len(coords) for i in range(3)]
    radius = max(math.dist(center, p) for p in coords)
    side = max(min_size, 2.0 * radius + padding)
    side = min(side, max_size)
    return {"center": center, "size": [side, side, side], "radius": radius}


def prepare_ligand_pdbqt(input_ligand: Path, out_pdbqt: Path) -> Tuple[bool, str, str]:
    """Prefer Meeko; fall back to OpenBabel."""
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    meeko = which_any(["mk_prepare_ligand.py", "mk_prepare_ligand"])
    if meeko:
        cp = run([meeko, "-i", str(input_ligand), "-o", str(out_pdbqt)], check=False)
        if cp.returncode == 0 and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            return True, "Meeko", cp.stdout + "\n" + cp.stderr

    obabel = which_any(["obabel", "babel"])
    if obabel:
        cp = run([
            obabel,
            str(input_ligand),
            "-O", str(out_pdbqt),
            "-xh",
            "--partialcharge", "gasteiger",
        ], check=False)
        if cp.returncode == 0 and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            return True, "OpenBabel", cp.stdout + "\n" + cp.stderr

    return False, "none", "No working ligand-preparation backend"


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    box: Dict,
    out_pose: Path,
    log_path: Path,
    exhaustiveness: int,
    num_modes: int,
    seed: int,
) -> Dict:
    vina = which_any(["vina", "vina.exe"])
    if not vina:
        return {"status": "failed", "reason": "vina executable not found"}

    c = box["center"]
    s = box["size"]
    cmd = [
        vina,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(c[0]),
        "--center_y", str(c[1]),
        "--center_z", str(c[2]),
        "--size_x", str(s[0]),
        "--size_y", str(s[1]),
        "--size_z", str(s[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--out", str(out_pose),
    ]
    cp = run(cmd, check=False)
    log_path.write_text(cp.stdout + "\n" + cp.stderr, encoding="utf-8")
    best = parse_vina_best(cp.stdout)
    return {
        "status": "ok" if cp.returncode == 0 and best is not None else "failed",
        "returncode": cp.returncode,
        "best_affinity_kcal_mol": best,
        "pose": str(out_pose),
        "log": str(log_path),
        "command": cmd,
    }


def convert_pose_pdbqt_to_sdf(pdbqt: Path, sdf: Path) -> bool:
    obabel = which_any(["obabel", "babel"])
    if not obabel:
        return False
    cp = run([obabel, str(pdbqt), "-O", str(sdf)], check=False)
    return cp.returncode == 0 and sdf.exists() and sdf.stat().st_size > 0


def convert_pdb_to_sdf(pdb: Path, sdf: Path) -> bool:
    obabel = which_any(["obabel", "babel"])
    if not obabel:
        return False
    cp = run([obabel, str(pdb), "-O", str(sdf)], check=False)
    return cp.returncode == 0 and sdf.exists() and sdf.stat().st_size > 0


def compute_best_pose_rmsd(crystal_sdf: Path, docked_sdf: Path) -> Dict:
    """Graph-map docked IT1t to crystal IT1t and report heavy-atom RMSD.

    This is a QC metric only. If atom typing/connectivity conversion prevents a
    valid mapping, RMSD is reported as unavailable rather than fabricated.
    """
    crystal = Chem.SDMolSupplier(str(crystal_sdf), removeHs=False)[0]
    docked = Chem.SDMolSupplier(str(docked_sdf), removeHs=False)[0]
    if crystal is None or docked is None:
        return {"rmsd_available": False, "reason": "RDKit could not parse SDF"}

    crystal_noh = Chem.RemoveHs(crystal)
    docked_noh = Chem.RemoveHs(docked)

    match = docked_noh.GetSubstructMatch(crystal_noh)
    reverse = False
    if not match:
        match2 = crystal_noh.GetSubstructMatch(docked_noh)
        if not match2:
            return {"rmsd_available": False, "reason": "No graph-based heavy-atom mapping"}
        reverse = True
        atom_map = [(i, j) for i, j in enumerate(match2)]
        prb = crystal_noh
        ref = docked_noh
    else:
        atom_map = [(i, j) for i, j in enumerate(match)]
        prb = docked_noh
        ref = crystal_noh

    try:
        rmsd = float(rdMolAlign.AlignMol(prb, ref, atomMap=atom_map))
    except Exception as exc:
        return {"rmsd_available": False, "reason": f"Alignment failed: {exc}"}

    return {
        "rmsd_available": True,
        "heavy_atom_rmsd_A": rmsd,
        "atom_count": len(atom_map),
        "mapping_reversed": reverse,
        "redocking_qc_pass_rmsd_lt_2A": bool(rmsd < 2.0),
    }


def find_existing_cxcr4_receptor() -> Path:
    candidates = [
        STAGE2 / "docking" / "CXCR4" / "CXCR4_3ODU_receptor.pdbqt",
        STAGE2 / "docking" / "CXCR4" / "receptor.pdbqt",
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            return p
    raise RuntimeError("Existing production CXCR4 receptor PDBQT not found")


def find_existing_3odu_pdb() -> Path:
    candidates = [
        STAGE2 / "receptors" / "3ODU.pdb",
        ROOT / "outputs" / "stage2_docking" / "receptors" / "3ODU.pdb",
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            return p
    raise RuntimeError("3ODU.pdb from Stage-2 not found")


def find_existing_dinp_pdbqt() -> Path:
    candidates = [
        STAGE2 / "ligand" / "DINP.pdbqt",
        STAGE2 / "ligand" / "DINP_best_3d.pdbqt",
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            return p
    raise RuntimeError("Prepared DINP PDBQT from production docking not found")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runtime", type=Path, default=RUNTIME)
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--num-modes", type=int, default=20)
    ap.add_argument("--seed", type=int, default=590836)
    ap.add_argument("--chain", default="A")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    prepend_runtime_to_path(args.runtime)

    audit: Dict = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "target": TARGET,
        "pdb_id": PDB_ID,
        "positive_control": "IT1t",
        "pdb_residue_id": POCKET_RESNAME,
        "runtime": str(args.runtime),
        "exhaustiveness": args.exhaustiveness,
        "num_modes": args.num_modes,
        "seed": args.seed,
    }

    receptor = find_existing_cxcr4_receptor()
    pdb = find_existing_3odu_pdb()
    dinp_pdbqt = find_existing_dinp_pdbqt()

    crystal_itd_pdb = OUT / "IT1t_crystal_ITD.pdb"
    n_atoms = extract_itd_from_pdb(pdb, crystal_itd_pdb, chain=args.chain)
    box = ligand_center_and_box_from_pdb(crystal_itd_pdb)
    audit["crystal_it1t_atom_records"] = n_atoms
    audit["docking_box"] = box

    # Prepare IT1t as the positive-control docking ligand.
    it1t_pdbqt = OUT / "IT1t_control.pdbqt"
    prep_ok, prep_method, prep_log = prepare_ligand_pdbqt(crystal_itd_pdb, it1t_pdbqt)
    (OUT / "IT1t_prepare.log").write_text(prep_log, encoding="utf-8")
    audit["it1t_prep"] = {"ok": prep_ok, "method": prep_method}
    if not prep_ok:
        save_json(audit, OUT / "audit.json")
        raise SystemExit("IT1t ligand preparation failed")

    # Re-dock IT1t.
    it1t_result = run_vina(
        receptor,
        it1t_pdbqt,
        box,
        OUT / "IT1t_redocked.pdbqt",
        OUT / "IT1t_redocking.log",
        args.exhaustiveness,
        args.num_modes,
        args.seed,
    )
    audit["it1t_redocking"] = it1t_result

    # Re-run DINP with the exact same receptor and IT1t-defined pocket box for
    # a strict within-protocol affinity comparison.
    dinp_result = run_vina(
        receptor,
        dinp_pdbqt,
        box,
        OUT / "DINP_same_box.pdbqt",
        OUT / "DINP_same_box.log",
        args.exhaustiveness,
        args.num_modes,
        args.seed,
    )
    audit["dinp_same_box"] = dinp_result

    # Optional redocking RMSD QC.
    rmsd_result = {"rmsd_available": False, "reason": "not attempted"}
    crystal_sdf = OUT / "IT1t_crystal.sdf"
    docked_sdf = OUT / "IT1t_redocked.sdf"
    if it1t_result.get("status") == "ok":
        if convert_pdb_to_sdf(crystal_itd_pdb, crystal_sdf) and convert_pose_pdbqt_to_sdf(OUT / "IT1t_redocked.pdbqt", docked_sdf):
            rmsd_result = compute_best_pose_rmsd(crystal_sdf, docked_sdf)
        else:
            rmsd_result = {"rmsd_available": False, "reason": "OpenBabel conversion unavailable/failed"}
    audit["redocking_rmsd_qc"] = rmsd_result

    it1t_aff = it1t_result.get("best_affinity_kcal_mol")
    dinp_aff = dinp_result.get("best_affinity_kcal_mol")
    delta = None
    if it1t_aff is not None and dinp_aff is not None:
        delta = float(dinp_aff - it1t_aff)

    comparison = {
        "target": TARGET,
        "pdb_id": PDB_ID,
        "pocket_reference": "3ODU co-crystallized IT1t (ITD)",
        "same_receptor": True,
        "same_box": True,
        "same_vina_parameters": True,
        "IT1t_best_affinity_kcal_mol": it1t_aff,
        "DINP_best_affinity_kcal_mol": dinp_aff,
        "DINP_minus_IT1t_kcal_mol": delta,
        "redocking_heavy_atom_rmsd_A": rmsd_result.get("heavy_atom_rmsd_A"),
        "redocking_qc_pass_rmsd_lt_2A": rmsd_result.get("redocking_qc_pass_rmsd_lt_2A"),
    }

    with (OUT / "affinity_comparison.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(comparison.keys()))
        writer.writeheader()
        writer.writerow(comparison)

    save_json(audit, OUT / "audit.json")

    lines = [
        "# CXCR4 control-docking validation",
        "",
        f"- Target: {TARGET} ({PDB_ID})",
        "- Positive control: IT1t (PDB ligand ID ITD)",
        "- DINP and IT1t were docked to the same prepared receptor and the same IT1t-defined pocket box.",
        f"- IT1t best Vina affinity: {it1t_aff if it1t_aff is not None else 'NA'} kcal/mol",
        f"- DINP best Vina affinity: {dinp_aff if dinp_aff is not None else 'NA'} kcal/mol",
        f"- DINP − IT1t score difference: {delta if delta is not None else 'NA'} kcal/mol",
        f"- IT1t redocking heavy-atom RMSD: {rmsd_result.get('heavy_atom_rmsd_A', 'NA')} Å",
        f"- RMSD <2 Å QC pass: {rmsd_result.get('redocking_qc_pass_rmsd_lt_2A', 'NA')}",
        "",
        "Interpretation: the IT1t redocking serves as a protocol QC/positive control. A more negative IT1t score than DINP would indicate that DINP is computationally weaker than the canonical co-crystallized antagonist under identical docking conditions. This comparison does not establish in-vivo binding.",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    if it1t_result.get("status") != "ok" or dinp_result.get("status") != "ok":
        raise SystemExit("Control docking did not fully complete; inspect audit.json/logs")

    print(json.dumps(comparison, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
