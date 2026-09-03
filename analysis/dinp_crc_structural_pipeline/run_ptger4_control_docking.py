#!/usr/bin/env python3
"""PTGER4 positive-control redocking and DINP comparison for 9JQZ.

Purpose
-------
Validate the PTGER4 docking protocol before interpreting the strong DINP score.
The script:
1) reuses the rescued PTGER4 9JQZ receptor;
2) extracts the co-crystallized ligand nearest the previously defined PTGER4 pocket;
3) prepares that ligand with Meeko;
4) redocks the co-crystal ligand into the same PTGER4 pocket;
5) computes heavy-atom RMSD versus the crystal pose;
6) docks DINP with the exact same frozen receptor / box / Vina settings;
7) writes an affinity comparison and protocol-QC summary.

Important: the protocol is NOT tuned on DINP affinity.
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
from typing import Dict, List, Optional, Tuple

from Bio.PDB import MMCIFParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import rdMolAlign

ROOT = Path(__file__).resolve().parent
STAGE2 = ROOT / "outputs" / "stage2_docking"
REFINE = ROOT / "outputs" / "stage2_refine_and_rescue"
OUT = ROOT / "outputs" / "ptger4_control_docking"
OUT.mkdir(parents=True, exist_ok=True)
RUNTIME = Path(r"D:\whynot17\dinp_stage2_runtime")


def run(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, text=True, capture_output=True, check=False)


def save_json(obj, path: Path):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def which(names: List[str]) -> Optional[str]:
    for name in names:
        p = shutil.which(name)
        if p:
            return p
    return None


def prepend_runtime_path():
    candidates = [RUNTIME, RUNTIME / "Scripts", RUNTIME / "bin", RUNTIME / "Library" / "bin"]
    existing = [str(p) for p in candidates if p.exists()]
    os.environ["PATH"] = os.pathsep.join(existing + [os.environ.get("PATH", "")])


def find_existing(paths: List[Path]) -> Path:
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError("None of these expected files exists:\n" + "\n".join(map(str, paths)))


def parse_vina_best(text: str) -> Optional[float]:
    for line in text.splitlines():
        f = line.split()
        if len(f) >= 2 and f[0].isdigit():
            try:
                return float(f[1])
            except ValueError:
                pass
    return None


def vina_dock(vina: str, receptor: Path, ligand: Path, center: List[float], size: List[float],
              exhaustiveness: int, num_modes: int, seed: int, out_pose: Path, log: Path) -> Dict:
    cmd = [
        vina, "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes),
        "--seed", str(seed), "--out", str(out_pose),
    ]
    cp = run(cmd)
    log.write_text(cp.stdout + "\n" + cp.stderr, encoding="utf-8")
    return {"returncode": cp.returncode, "affinity": parse_vina_best(cp.stdout), "command": cmd}


def load_ptger4_pocket() -> Tuple[List[float], List[float]]:
    audits = [
        STAGE2 / "pocket_audit.json",
        STAGE2 / "stage2_pocket_audit.json",
        REFINE / "audit.json",
    ]
    for path in audits:
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        candidates = []
        if isinstance(data, dict):
            candidates.append(data.get("PTGER4"))
            candidates.append(data.get("targets", {}).get("PTGER4") if isinstance(data.get("targets"), dict) else None)
            candidates.append(data.get("ptger4"))
        for pt in candidates:
            if not isinstance(pt, dict):
                continue
            box = pt.get("proposed_box") or pt.get("box")
            if box and "center" in box and "size" in box:
                return [float(x) for x in box["center"]], [float(x) for x in box["size"]]
    # Fallback: recover center/size from rescue summary audit if present.
    audit = REFINE / "audit.json"
    if audit.exists():
        data = json.loads(audit.read_text(encoding="utf-8"))
        pt = data.get("ptger4_rescue", {})
        if "pocket_center" in pt and "pocket_size" in pt:
            return [float(x) for x in pt["pocket_center"]], [float(x) for x in pt["pocket_size"]]
    raise RuntimeError("PTGER4 pocket center/size not found in existing Stage-2 outputs")


class ResidueSelect(Select):
    def __init__(self, chain_id: str, residue_id):
        self.chain_id = chain_id
        self.residue_id = residue_id
    def accept_chain(self, chain):
        return int(chain.id == self.chain_id)
    def accept_residue(self, residue):
        return int(residue.id == self.residue_id)


def residue_centroid(residue) -> Tuple[float, float, float]:
    pts = []
    for atom in residue.get_atoms():
        if atom.element and atom.element.upper() == "H":
            continue
        pts.append(atom.coord)
    if not pts:
        raise ValueError("No heavy atoms")
    n = float(len(pts))
    return (sum(float(p[0]) for p in pts)/n, sum(float(p[1]) for p in pts)/n, sum(float(p[2]) for p in pts)/n)


def dist(a, b) -> float:
    return math.sqrt(sum((float(a[i]) - float(b[i]))**2 for i in range(3)))


def extract_nearest_cocrystal_ligand(cif: Path, center: List[float], out_pdb: Path) -> Dict:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("9JQZ", str(cif))
    hits = []
    excluded = {"HOH", "WAT", "NA", "CL", "K", "CA", "MG", "ZN", "SO4", "PO4", "GOL", "EDO"}
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0]
                resname = residue.get_resname().strip()
                if hetflag == " " or resname in excluded:
                    continue
                heavy = [a for a in residue.get_atoms() if not (a.element and a.element.upper() == "H")]
                if len(heavy) < 6:
                    continue
                c = residue_centroid(residue)
                hits.append({
                    "chain": chain.id,
                    "resname": resname,
                    "residue_id": residue.id,
                    "heavy_atoms": len(heavy),
                    "centroid": list(c),
                    "distance_to_pocket_center_A": dist(c, center),
                    "residue_obj": residue,
                })
    if not hits:
        raise RuntimeError("No plausible co-crystallized ligand found in 9JQZ")
    hits.sort(key=lambda x: (x["distance_to_pocket_center_A"], -x["heavy_atoms"]))
    best = hits[0]
    io = PDBIO(); io.set_structure(structure)
    io.save(str(out_pdb), ResidueSelect(best["chain"], best["residue_id"]))
    report = {k: v for k, v in best.items() if k not in {"residue_obj", "residue_id"}}
    report["residue_id"] = [str(x) for x in best["residue_id"]]
    report["all_candidates"] = [
        {k: v for k, v in h.items() if k not in {"residue_obj", "residue_id"}} for h in hits[:10]
    ]
    return report


def convert_pdb_to_sdf(obabel: str, pdb: Path, sdf: Path) -> bool:
    cp = run([obabel, str(pdb), "-O", str(sdf)])
    (OUT / "crystal_ligand_conversion.log").write_text(cp.stdout + "\n" + cp.stderr, encoding="utf-8")
    return cp.returncode == 0 and sdf.exists() and sdf.stat().st_size > 0


def prepare_ligand(mklig: str, sdf: Path, pdbqt: Path) -> bool:
    attempts = [
        [mklig, "-i", str(sdf), "-o", str(pdbqt)],
        [mklig, "--input", str(sdf), "--output_filename", str(pdbqt)],
    ]
    logs = []
    for cmd in attempts:
        cp = run(cmd)
        logs.append({"cmd": cmd, "returncode": cp.returncode, "stdout": cp.stdout, "stderr": cp.stderr})
        if cp.returncode == 0 and pdbqt.exists() and pdbqt.stat().st_size > 0:
            save_json(logs, OUT / "ligand_prep_attempts.json")
            return True
    save_json(logs, OUT / "ligand_prep_attempts.json")
    return False


def mol_from_sdf(path: Path):
    supp = Chem.SDMolSupplier(str(path), removeHs=True)
    for m in supp:
        if m is not None:
            return m
    return None


def heavy_atom_rmsd_pdbqt(reference_sdf: Path, docked_pdbqt: Path) -> Optional[float]:
    ref = mol_from_sdf(reference_sdf)
    if ref is None:
        return None
    smiles = None
    smiles_to_pdb: Dict[int, int] = {}
    coordinates: Dict[int, Tuple[float, float, float]] = {}
    saw_model = False
    for line in docked_pdbqt.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("REMARK SMILES IDX"):
            fields = line.split()[3:]
            if len(fields) % 2 == 0:
                for i in range(0, len(fields), 2):
                    smiles_to_pdb[int(fields[i])] = int(fields[i+1])
        elif line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX"):
            smiles = line.split("REMARK SMILES ", 1)[1].strip()
        elif line.startswith("MODEL"):
            if saw_model:
                break
            saw_model = True
        elif saw_model and line.startswith("ENDMDL"):
            break
        elif saw_model and line.startswith(("ATOM", "HETATM")):
            try:
                serial = int(line[6:11])
                coordinates[serial] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except Exception:
                pass
    if not smiles or not smiles_to_pdb or not coordinates:
        return None
    query = Chem.MolFromSmiles(smiles)
    if query is None:
        return None
    match = ref.GetSubstructMatch(query)
    if not match or len(match) != query.GetNumAtoms():
        return None
    docked = Chem.Mol(ref); docked.RemoveAllConformers()
    conf = Chem.Conformer(docked.GetNumAtoms())
    atom_map = []
    for smiles_idx, ref_idx in enumerate(match, start=1):
        serial = smiles_to_pdb.get(smiles_idx)
        if serial not in coordinates:
            return None
        conf.SetAtomPosition(ref_idx, coordinates[serial])
        atom_map.append((ref_idx, ref_idx))
    docked.AddConformer(conf, assignId=True)
    try:
        return float(rdMolAlign.AlignMol(docked, ref, atomMap=atom_map))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--seed", type=int, default=590836)
    ap.add_argument("--num-modes", type=int, default=20)
    args = ap.parse_args()

    prepend_runtime_path()
    vina = which(["vina.exe", "vina"])
    obabel = which(["obabel.exe", "obabel"])
    mklig = which(["mk_prepare_ligand.py", "mk_prepare_ligand.exe", "mk_prepare_ligand"])
    if not vina or not obabel or not mklig:
        raise RuntimeError(f"Missing runtime tool(s): vina={vina}, obabel={obabel}, mk_prepare_ligand={mklig}")

    receptor = find_existing([
        REFINE / "ptger4_rescue" / "PTGER4_9JQZ_rescued_receptor.pdbqt",
    ])
    dinp = find_existing([
        STAGE2 / "ligand" / "DINP.pdbqt",
        ROOT / "outputs" / "cxcr4_control_docking" / "DINP.pdbqt",
    ])
    cif = find_existing([
        STAGE2 / "receptors" / "9JQZ.cif",
        ROOT / "outputs" / "receptors" / "9JQZ.cif",
        ROOT / "outputs" / "structures" / "9JQZ.cif",
    ])
    center, size = load_ptger4_pocket()

    crystal_pdb = OUT / "PTGER4_9JQZ_cocrystal_ligand.pdb"
    crystal_sdf = OUT / "PTGER4_9JQZ_cocrystal_ligand.sdf"
    control_pdbqt = OUT / "PTGER4_control_ligand.pdbqt"
    ligand_info = extract_nearest_cocrystal_ligand(cif, center, crystal_pdb)
    if not convert_pdb_to_sdf(obabel, crystal_pdb, crystal_sdf):
        raise RuntimeError("Could not convert PTGER4 crystal ligand PDB to SDF")
    if not prepare_ligand(mklig, crystal_sdf, control_pdbqt):
        raise RuntimeError("Could not prepare PTGER4 co-crystal ligand with Meeko")

    ctrl_pose = OUT / "PTGER4_control_redocked.pdbqt"
    ctrl_log = OUT / "PTGER4_control_redocking.log"
    ctrl = vina_dock(vina, receptor, control_pdbqt, center, size, args.exhaustiveness,
                     args.num_modes, args.seed, ctrl_pose, ctrl_log)
    rmsd = heavy_atom_rmsd_pdbqt(crystal_sdf, ctrl_pose) if ctrl["returncode"] == 0 else None

    dinp_pose = OUT / "DINP_PTGER4_same_protocol.pdbqt"
    dinp_log = OUT / "DINP_PTGER4_same_protocol.log"
    dinp_res = vina_dock(vina, receptor, dinp, center, size, args.exhaustiveness,
                         args.num_modes, args.seed, dinp_pose, dinp_log)

    rows = [
        {"ligand": f"co-crystal control ({ligand_info['resname']})", "affinity_kcal_mol": ctrl["affinity"], "rmsd_A": rmsd},
        {"ligand": "DINP", "affinity_kcal_mol": dinp_res["affinity"], "rmsd_A": ""},
    ]
    with (OUT / "affinity_comparison.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    qc = bool(rmsd is not None and rmsd < 2.0)
    delta = None
    if ctrl["affinity"] is not None and dinp_res["affinity"] is not None:
        delta = float(dinp_res["affinity"] - ctrl["affinity"])
    audit = {
        "target": "PTGER4",
        "structure": "9JQZ",
        "control_ligand": ligand_info,
        "receptor": str(receptor),
        "pocket_center": center,
        "pocket_size": size,
        "exhaustiveness": args.exhaustiveness,
        "seed": args.seed,
        "num_modes": args.num_modes,
        "control_affinity_kcal_mol": ctrl["affinity"],
        "dinp_affinity_kcal_mol": dinp_res["affinity"],
        "dinp_minus_control_kcal_mol": delta,
        "control_redocking_rmsd_A": rmsd,
        "qc_pass_lt2A": qc,
        "protocol_selection": "no DINP-based tuning; same receptor, pocket and Vina settings for control and DINP",
    }
    save_json(audit, OUT / "audit.json")

    summary = [
        "# PTGER4 control-docking validation",
        "",
        "- Target: PTGER4 (9JQZ)",
        f"- Co-crystal control selected automatically: {ligand_info['resname']} (chain {ligand_info['chain']})",
        f"- Control best Vina affinity: {ctrl['affinity']} kcal/mol",
        f"- DINP best Vina affinity under identical protocol: {dinp_res['affinity']} kcal/mol",
        f"- DINP − control score difference: {delta} kcal/mol",
        f"- Control redocking heavy-atom RMSD: {rmsd} Å",
        f"- RMSD <2 Å QC pass: {qc}",
        "",
        "Interpretation: only compare DINP with the PTGER4 co-crystal control within this same receptor/pocket/protocol. "
        "A more negative docking score does not prove stronger experimental binding; redocking validates pose reproduction, not biological affinity.",
    ]
    (OUT / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
