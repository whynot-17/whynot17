#!/usr/bin/env python3
"""Refine CXCR4 control docking and rescue PTGER4 docking.

Goals
-----
1. CXCR4 protocol refinement:
   * Use the same 3ODU receptor and IT1t/ITD pocket.
   * Sweep a prespecified set of box sizes, Vina exhaustiveness values, and seeds.
   * Redock IT1t for every setting.
   * Compute heavy-atom RMSD to the crystal pose.
   * Freeze the best protocol only if RMSD < 2.0 A; otherwise retain the best
     achieved setting as exploratory and report QC failure.
   * Re-dock DINP using the frozen/best setting, without tuning on DINP score.

2. PTGER4 rescue:
   * Convert 9JQZ mmCIF to a clean PDB retaining protein atoms for the selected chain.
   * Use the Stage-2 pocket audit A1ECR pocket coordinates.
   * Prepare the rescued receptor with Meeko.
   * Dock DINP with the same Vina defaults used in the main production workflow.

This script never selects CXCR4 parameters by DINP affinity. Parameter selection is
based only on positive-control IT1t redocking RMSD to avoid outcome-driven tuning.
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
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from Bio.PDB import MMCIFParser, PDBIO, Select

ROOT = Path(__file__).resolve().parent
STAGE2 = ROOT / "outputs" / "stage2_docking"
CTRL = ROOT / "outputs" / "cxcr4_control_docking"
OUT = ROOT / "outputs" / "stage2_refine_and_rescue"
OUT.mkdir(parents=True, exist_ok=True)

RUNTIME = Path(r"D:\whynot17\dinp_stage2_runtime")


def save_json(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def prepend_runtime_path():
    candidates = [
        RUNTIME,
        RUNTIME / "Scripts",
        RUNTIME / "bin",
        RUNTIME / "Library" / "bin",
    ]
    existing = [str(p) for p in candidates if p.exists()]
    os.environ["PATH"] = os.pathsep.join(existing + [os.environ.get("PATH", "")])


def which(names: List[str]) -> Optional[str]:
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def run(cmd: List[str], cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True,
                          capture_output=True, check=False)


def parse_vina_best(text: str) -> Optional[float]:
    for line in text.splitlines():
        f = line.split()
        if len(f) >= 2 and f[0].isdigit():
            try:
                return float(f[1])
            except ValueError:
                pass
    return None


def vina_dock(vina: str, receptor: Path, ligand: Path, center: List[float],
              size: List[float], exhaustiveness: int, num_modes: int,
              seed: int, out_pose: Path, log: Path) -> Dict:
    cmd = [
        vina, "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes),
        "--seed", str(seed), "--out", str(out_pose),
    ]
    cp = run(cmd)
    log.write_text(cp.stdout + "\n" + cp.stderr, encoding="utf-8")
    return {
        "returncode": cp.returncode,
        "affinity": parse_vina_best(cp.stdout),
        "pose": str(out_pose),
        "command": cmd,
    }


def pdbqt_to_sdf(obabel: str, pdbqt: Path, sdf: Path) -> bool:
    cp = run([obabel, str(pdbqt), "-O", str(sdf)])
    return cp.returncode == 0 and sdf.exists() and sdf.stat().st_size > 0


def mol_from_sdf(path: Path):
    supp = Chem.SDMolSupplier(str(path), removeHs=True)
    for m in supp:
        if m is not None:
            return m
    return None


def heavy_atom_rmsd(reference_sdf: Path, docked_sdf: Path) -> Optional[float]:
    ref = mol_from_sdf(reference_sdf)
    dock = mol_from_sdf(docked_sdf)
    if ref is None or dock is None or ref.GetNumAtoms() != dock.GetNumAtoms():
        return None
    try:
        return float(rdMolAlign.GetBestRMS(ref, dock))
    except Exception:
        return None


def heavy_atom_rmsd_pdbqt(reference_sdf: Path, docked_pdbqt: Path) -> Optional[float]:
    """Map the first Vina pose to the reference graph without OpenBabel.

    OpenBabel can convert PDBQT to SDF but may lose the bond orders needed for
    RDKit graph matching.  Meeko/Vina writes a SMILES plus ``SMILES IDX`` atom
    map, which lets us place the docked heavy-atom coordinates on the original
    reference graph and calculate a reproducible RMSD directly.
    """
    ref = mol_from_sdf(reference_sdf)
    if ref is None:
        return None

    smiles = None
    smiles_to_pdb: Dict[int, int] = {}
    coordinates: Dict[int, Tuple[float, float, float]] = {}
    saw_model = False
    model_finished = False
    for line in docked_pdbqt.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("REMARK SMILES IDX"):
            fields = line.split()[3:]
            if len(fields) % 2 == 0:
                for i in range(0, len(fields), 2):
                    smiles_to_pdb[int(fields[i])] = int(fields[i + 1])
            continue
        if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX"):
            smiles = line.split("REMARK SMILES ", 1)[1].strip()
            continue
        if line.startswith("MODEL"):
            if saw_model:
                break
            saw_model = True
            continue
        if saw_model and line.startswith("ENDMDL"):
            model_finished = True
            break
        if saw_model and line.startswith(("ATOM", "HETATM")):
            try:
                serial = int(line[6:11])
                atom_type = line[76:78].strip() if len(line) >= 78 else line[12:16].strip()
                if atom_type.upper().startswith("H"):
                    continue
                coordinates[serial] = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            except (ValueError, IndexError):
                continue
    if not saw_model or not model_finished or not smiles or not smiles_to_pdb:
        return None

    query = Chem.MolFromSmiles(smiles)
    if query is None:
        return None
    match = ref.GetSubstructMatch(query)
    if not match or len(match) != query.GetNumAtoms():
        return None

    docked = Chem.Mol(ref)
    docked.RemoveAllConformers()
    conformer = Chem.Conformer(docked.GetNumAtoms())
    atom_map = []
    for smiles_idx, ref_idx in enumerate(match, start=1):
        pdb_serial = smiles_to_pdb.get(smiles_idx)
        if pdb_serial is None or pdb_serial not in coordinates:
            return None
        conformer.SetAtomPosition(ref_idx, coordinates[pdb_serial])
        atom_map.append((ref_idx, ref_idx))
    docked.AddConformer(conformer, assignId=True)
    try:
        return float(rdMolAlign.AlignMol(docked, ref, atomMap=atom_map))
    except Exception:
        return None


def find_existing(paths: List[Path]) -> Path:
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError("None of the expected files exist:\n" + "\n".join(map(str, paths)))


def load_cxcr4_inputs() -> Dict[str, Path]:
    receptor = find_existing([
        STAGE2 / "docking" / "CXCR4" / "CXCR4_3ODU_receptor.pdbqt",
        CTRL / "CXCR4_3ODU_receptor.pdbqt",
    ])
    dinp = find_existing([
        STAGE2 / "ligand" / "DINP.pdbqt",
        CTRL / "DINP.pdbqt",
    ])
    it1t = find_existing([
        CTRL / "IT1t_control.pdbqt",
        CTRL / "IT1t.pdbqt",
        CTRL / "ITD.pdbqt",
        CTRL / "IT1t_crystal.pdbqt",
    ])
    crystal_sdf = find_existing([
        CTRL / "IT1t_crystal_ITD.converted.sdf",
        CTRL / "IT1t_crystal.sdf",
        CTRL / "ITD_crystal.sdf",
        CTRL / "IT1t.sdf",
    ])
    return {"receptor": receptor, "dinp": dinp, "it1t": it1t, "crystal_sdf": crystal_sdf}


def cxcr4_refinement(vina: str, obabel: str, pocket_audit: Dict, args) -> Dict:
    inputs = load_cxcr4_inputs()
    cx = pocket_audit["CXCR4"] if "CXCR4" in pocket_audit else pocket_audit.get("targets", {}).get("CXCR4")
    if not cx:
        raise RuntimeError("CXCR4 pocket not found in pocket audit")
    box = cx.get("proposed_box") or cx.get("box")
    if not box:
        raise RuntimeError("CXCR4 proposed box missing")
    center = [float(x) for x in box["center"]]
    base_size = [float(x) for x in box["size"]]

    # Prespecified refinement grid. Selection uses IT1t RMSD only.
    size_scales = [0.70, 0.80, 0.90, 1.00]
    exhaustiveness_values = [32, 64, 128]
    seeds = [590836, 3001, 3002]

    rows = []
    refdir = OUT / "cxcr4_refinement"
    refdir.mkdir(parents=True, exist_ok=True)
    run_id = 0
    for scale in size_scales:
        size = [max(16.0, x * scale) for x in base_size]
        for ex in exhaustiveness_values:
            for seed in seeds:
                run_id += 1
                stem = f"it1t_s{scale:.2f}_e{ex}_seed{seed}"
                pose = refdir / f"{stem}.pdbqt"
                log = refdir / f"{stem}.log"
                sdf = refdir / f"{stem}.sdf"
                resumed_from_cache = False
                cached_affinity = None
                if pose.exists() and log.exists():
                    cached_affinity = parse_vina_best(log.read_text(encoding="utf-8", errors="ignore"))
                if cached_affinity is not None:
                    # A valid pose/log pair is enough to resume the Vina
                    # calculation.  SDF is only a fallback RMSD route; the
                    # primary RMSD implementation reads the PDBQT directly.
                    res = {
                        "returncode": 0,
                        "affinity": cached_affinity,
                        "pose": str(pose),
                        "command": ["RESUMED_FROM_CACHE"],
                    }
                    resumed_from_cache = True
                else:
                    res = vina_dock(vina, inputs["receptor"], inputs["it1t"], center, size,
                                    ex, args.num_modes, seed, pose, log)
                rmsd = None
                if res["returncode"] == 0 and pose.exists():
                    rmsd = heavy_atom_rmsd_pdbqt(inputs["crystal_sdf"], pose)
                    if rmsd is None and pdbqt_to_sdf(obabel, pose, sdf):
                        rmsd = heavy_atom_rmsd(inputs["crystal_sdf"], sdf)
                rows.append({
                    "run_id": run_id,
                    "size_scale": scale,
                    "size_x": size[0], "size_y": size[1], "size_z": size[2],
                    "exhaustiveness": ex,
                    "seed": seed,
                    "it1t_affinity_kcal_mol": res["affinity"],
                    "it1t_rmsd_A": rmsd,
                    "qc_pass_lt2A": bool(rmsd is not None and rmsd < 2.0),
                    "resumed_from_cache": resumed_from_cache,
                    "pose": str(pose),
                })
                # Checkpoint the grid after every completed run so a later
                # interruption never discards already finished combinations.
                with (refdir / "refinement_grid.csv").open("w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                    w.writeheader(); w.writerows(rows)

    with (refdir / "refinement_grid.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    valid = [r for r in rows if r["it1t_rmsd_A"] is not None]
    if not valid:
        raise RuntimeError("No CXCR4 redocking RMSD could be computed")
    valid.sort(key=lambda r: (r["it1t_rmsd_A"], -(r["it1t_affinity_kcal_mol"] or 999)))
    best = valid[0]
    frozen = dict(best)
    frozen["selection_rule"] = "minimum IT1t redocking heavy-atom RMSD; DINP score not used"
    frozen["protocol_qc_pass"] = bool(best["it1t_rmsd_A"] < 2.0)
    save_json(frozen, refdir / "frozen_protocol.json")

    # DINP is run only once under the frozen/best control-derived protocol.
    size = [best["size_x"], best["size_y"], best["size_z"]]
    dinp_pose = refdir / "DINP_frozen_protocol.pdbqt"
    dinp_log = refdir / "DINP_frozen_protocol.log"
    dinp_res = vina_dock(vina, inputs["receptor"], inputs["dinp"], center, size,
                         int(best["exhaustiveness"]), args.num_modes, int(best["seed"]),
                         dinp_pose, dinp_log)
    return {
        "best_control_protocol": frozen,
        "dinp_affinity_kcal_mol": dinp_res["affinity"],
        "dinp_pose": str(dinp_pose),
        "n_grid_runs": len(rows),
    }


class ProteinChainSelect(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return 1 if chain.id == self.chain_id else 0
    def accept_residue(self, residue):
        # Keep only standard amino-acid-like polymer residues; drop waters/ligands.
        return 1 if residue.id[0] == " " else 0


def rescue_ptger4(vina: str, pocket_audit: Dict, args) -> Dict:
    pt = pocket_audit["PTGER4"] if "PTGER4" in pocket_audit else pocket_audit.get("targets", {}).get("PTGER4")
    if not pt:
        raise RuntimeError("PTGER4 pocket not found in pocket audit")
    box = pt.get("proposed_box") or pt.get("box")
    if not box:
        raise RuntimeError("PTGER4 proposed box missing")
    chain = "A"
    src_cif = find_existing([
        STAGE2 / "receptors" / "9JQZ.cif",
        ROOT / "outputs" / "receptors" / "9JQZ.cif",
        ROOT / "outputs" / "structures" / "9JQZ.cif",
    ])
    rescue_dir = OUT / "ptger4_rescue"
    rescue_dir.mkdir(parents=True, exist_ok=True)
    clean_pdb = rescue_dir / "PTGER4_9JQZ_chainA_clean.pdb"
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("9JQZ", str(src_cif))
    io = PDBIO(); io.set_structure(structure); io.save(str(clean_pdb), ProteinChainSelect(chain))

    mkrec = which(["mk_prepare_receptor.py", "mk_prepare_receptor.exe", "mk_prepare_receptor"])
    if not mkrec:
        return {"status": "skipped", "reason": "Meeko receptor-prep executable not found", "clean_pdb": str(clean_pdb)}
    receptor_pdbqt = rescue_dir / "PTGER4_9JQZ_rescued_receptor.pdbqt"
    # The clean rescue input is a PDB.  Use Meeko's --read_pdb path rather
    # than -i/--read_with_prody, which is reserved for the optional ProDy
    # parser and unnecessarily makes this rescue depend on ProDy.
    attempts = [
        [
            mkrec, "--read_pdb", str(clean_pdb), "-p", str(receptor_pdbqt),
            "--delete_bad_res", "--default_altloc", "A",
        ],
        [mkrec, "--read_pdb", str(clean_pdb), "-p", str(receptor_pdbqt)],
    ]
    prep_logs = []
    prep_ok = False
    for cmd in attempts:
        cp = run(cmd)
        prep_logs.append({"cmd": cmd, "returncode": cp.returncode, "stdout": cp.stdout, "stderr": cp.stderr})
        if cp.returncode == 0 and receptor_pdbqt.exists() and receptor_pdbqt.stat().st_size > 0:
            prep_ok = True
            break
    save_json(prep_logs, rescue_dir / "receptor_prep_attempts.json")
    if not prep_ok:
        return {"status": "skipped", "reason": "rescued PDB still failed receptor preparation", "clean_pdb": str(clean_pdb)}

    dinp = find_existing([STAGE2 / "ligand" / "DINP.pdbqt", CTRL / "DINP.pdbqt"])
    center = [float(x) for x in box["center"]]
    size = [float(x) for x in box["size"]]
    pose = rescue_dir / "DINP_PTGER4_9JQZ_rescue_pose.pdbqt"
    log = rescue_dir / "DINP_PTGER4_9JQZ_rescue_vina.log"
    res = vina_dock(vina, receptor_pdbqt, dinp, center, size,
                    args.ptger4_exhaustiveness, args.num_modes, args.ptger4_seed,
                    pose, log)
    return {
        "status": "ok" if res["returncode"] == 0 else "failed",
        "affinity_kcal_mol": res["affinity"],
        "clean_pdb": str(clean_pdb),
        "receptor_pdbqt": str(receptor_pdbqt),
        "pose": str(pose),
        "pocket_center": center,
        "pocket_size": size,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--num-modes", type=int, default=20)
    ap.add_argument("--ptger4-exhaustiveness", type=int, default=32)
    ap.add_argument("--ptger4-seed", type=int, default=590836)
    args = ap.parse_args()

    prepend_runtime_path()
    vina = which(["vina.exe", "vina"])
    obabel = which(["obabel.exe", "obabel", "babel.exe", "babel"])
    if not vina:
        raise SystemExit("AutoDock Vina not found in PATH/runtime")
    if not obabel:
        raise SystemExit("OpenBabel not found in PATH/runtime")

    pocket_audit_path = STAGE2 / "pocket_audit.json"
    if not pocket_audit_path.exists():
        raise SystemExit(f"Missing pocket audit: {pocket_audit_path}")
    pocket_audit = load_json(pocket_audit_path)

    audit = {
        "vina": vina,
        "obabel": obabel,
        "cxcr4": None,
        "ptger4": None,
    }
    audit["cxcr4"] = cxcr4_refinement(vina, obabel, pocket_audit, args)
    audit["ptger4"] = rescue_ptger4(vina, pocket_audit, args)
    save_json(audit, OUT / "audit.json")

    lines = [
        "# Stage-2 docking refinement + rescue summary",
        "",
        "## CXCR4 control-derived refinement",
        f"- Best IT1t redocking RMSD: {audit['cxcr4']['best_control_protocol']['it1t_rmsd_A']:.3f} Å",
        f"- Protocol QC pass (<2 Å): {audit['cxcr4']['best_control_protocol']['protocol_qc_pass']}",
        f"- IT1t affinity under selected protocol: {audit['cxcr4']['best_control_protocol']['it1t_affinity_kcal_mol']} kcal/mol",
        f"- DINP affinity under the same frozen protocol: {audit['cxcr4']['dinp_affinity_kcal_mol']} kcal/mol",
        "- Selection of the protocol used only IT1t redocking RMSD, never DINP affinity.",
        "",
        "## PTGER4 rescue",
        f"- Status: {audit['ptger4'].get('status')}",
        f"- DINP affinity: {audit['ptger4'].get('affinity_kcal_mol')} kcal/mol",
        "",
        "Interpretation boundary: docking remains computational structural plausibility evidence. "
        "MD should be run only after manual inspection of the selected pose/receptor and, for CXCR4, preferably after protocol QC passes.",
    ]
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("Completed. Results:", OUT)


if __name__ == "__main__":
    main()
