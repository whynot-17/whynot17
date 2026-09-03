#!/usr/bin/env python3
"""Stage 2: DINP 3D generation, pocket-guided receptor preparation and AutoDock Vina docking.

This script is designed to follow the Stage-1 structural inventory already generated in
analysis/dinp_crc_structural_pipeline/outputs/structures/selected_structures.json.

Primary targets are PTGER4, CXCR4, MMP9 and STAT3. The code:
  1) builds an ensemble of DINP 3D conformers with RDKit and MMFF/UFF minimization;
  2) downloads the selected RCSB structure in PDB or mmCIF format;
  3) identifies a source-structure ligand pocket from non-polymer HETATM residues when present;
  4) writes pocket proposals for manual audit;
  5) prepares ligand/receptor PDBQT using Meeko/MGLTools/OpenBabel if available;
  6) runs AutoDock Vina only when the pocket and PDBQT preparation pass QC;
  7) emits a ranked docking summary and provenance log.

Important: an automatically detected co-crystal pocket is a prioritization aid, not proof that
DINP binds the target. Results should be manually inspected before MD.
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
from typing import Dict, Iterable, List, Optional, Tuple

import requests

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except Exception as exc:  # pragma: no cover
    raise SystemExit("RDKit is required for Stage 2 ligand generation") from exc

try:
    from Bio.PDB import PDBParser, MMCIFParser, PDBIO
except Exception as exc:  # pragma: no cover
    raise SystemExit("Biopython is required for receptor/pocket parsing") from exc


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "stage2_docking"
LIGDIR = OUT / "ligand"
RECDIR = OUT / "receptors"
DOCKDIR = OUT / "docking"

PRIMARY_TARGETS = ["PTGER4", "CXCR4", "MMP9", "STAT3"]
EXCLUDE_HET = {
    "HOH", "WAT", "DOD", "NA", "CL", "K", "CA", "MG", "ZN", "MN", "FE", "CU",
    "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "ACT", "FMT", "ACE", "TRS", "MES",
}


def run(cmd: List[str], cwd: Optional[Path] = None, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True, capture_output=True, check=check)


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def save_json(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, ensure_ascii=False)


def pubchem_canonical_smiles(cid: int) -> str:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    props = r.json()["PropertyTable"]["Properties"][0]
    return props["ConnectivitySMILES"] if "ConnectivitySMILES" in props else props["CanonicalSMILES"]


def build_ligand_ensemble(smiles: str, n_conf: int = 50, seed: int = 590836) -> Dict:
    LIGDIR.mkdir(parents=True, exist_ok=True)
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        raise RuntimeError("Could not parse DINP SMILES")
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.4
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not ids:
        raise RuntimeError("RDKit failed to generate DINP 3D conformers")

    mmff_ok = AllChem.MMFFHasAllMoleculeParams(mol)
    scores: List[Tuple[int, float, str]] = []
    for cid in ids:
        if mmff_ok:
            try:
                AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=2000)
                props = AllChem.MMFFGetMoleculeProperties(mol)
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
                e = float(ff.CalcEnergy())
                method = "MMFF94"
            except Exception:
                AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=2000)
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
                e = float(ff.CalcEnergy())
                method = "UFF"
        else:
            AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=2000)
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            e = float(ff.CalcEnergy())
            method = "UFF"
        scores.append((cid, e, method))

    scores.sort(key=lambda x: x[1])
    best_id, best_e, method = scores[0]
    sdf = LIGDIR / "DINP_best_3d.sdf"
    writer = Chem.SDWriter(str(sdf))
    mol.SetProp("PubChemCID", "590836")
    mol.SetProp("3D_method", f"ETKDGv3+{method}")
    mol.SetProp("energy", str(best_e))
    writer.write(mol, confId=best_id)
    writer.close()

    pdb = LIGDIR / "DINP_best_3d.pdb"
    Chem.MolToPDBFile(mol, str(pdb), confId=best_id)
    return {
        "smiles": smiles,
        "conformers_generated": len(ids),
        "force_field": method,
        "best_energy": best_e,
        "sdf": str(sdf),
        "pdb": str(pdb),
    }


def fetch_rcsb_structure(pdb_id: str) -> Tuple[Path, str]:
    RECDIR.mkdir(parents=True, exist_ok=True)
    for ext, fmt in [("pdb", "pdb"), ("cif", "mmcif")]:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.{ext}"
        r = requests.get(url, timeout=90)
        if r.ok and len(r.content) > 1000:
            path = RECDIR / f"{pdb_id.upper()}.{ext}"
            path.write_bytes(r.content)
            return path, fmt
    raise RuntimeError(f"Could not download {pdb_id} from RCSB")


def parse_structure(path: Path, fmt: str):
    if fmt == "pdb":
        return PDBParser(QUIET=True).get_structure(path.stem, str(path))
    return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))


def candidate_ligands(structure) -> List[Dict]:
    hits = []
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0]
                resname = residue.resname.strip().upper()
                if not hetflag.startswith("H_") or resname in EXCLUDE_HET:
                    continue
                atoms = [a for a in residue.get_atoms() if a.element != "H"]
                if len(atoms) < 6:
                    continue
                coords = [a.coord for a in atoms]
                center = [float(sum(c[i] for c in coords) / len(coords)) for i in range(3)]
                maxdist = max(math.dist(center, [float(x) for x in c]) for c in coords)
                hits.append({
                    "resname": resname,
                    "chain": chain.id,
                    "resseq": int(residue.id[1]),
                    "heavy_atoms": len(atoms),
                    "center": center,
                    "radius": maxdist,
                })
        break
    hits.sort(key=lambda x: (x["heavy_atoms"], x["radius"]), reverse=True)
    return hits


def box_from_ligand(hit: Dict, min_size: float = 22.0, padding: float = 12.0) -> Dict:
    side = max(min_size, 2.0 * float(hit["radius"]) + padding)
    side = min(side, 34.0)
    return {"center": hit["center"], "size": [side, side, side], "source_ligand": hit}


def choose_binary(names: Iterable[str]) -> Optional[str]:
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def prepare_ligand_pdbqt(sdf: Path, out_pdbqt: Path) -> Tuple[bool, str]:
    meeko = choose_binary(["mk_prepare_ligand.exe", "mk_prepare_ligand.py", "mk_prepare_ligand"])
    if meeko:
        cp = run([meeko, "-i", str(sdf), "-o", str(out_pdbqt)], check=False)
        if cp.returncode == 0 and out_pdbqt.exists():
            return True, "Meeko mk_prepare_ligand.py"
    obabel = choose_binary(["obabel", "babel"])
    if obabel:
        cp = run([obabel, str(sdf), "-O", str(out_pdbqt), "-xh", "--partialcharge", "gasteiger"], check=False)
        if cp.returncode == 0 and out_pdbqt.exists():
            return True, "OpenBabel"
    return False, "No working Meeko/OpenBabel ligand-preparation executable"


def prepare_receptor_pdbqt(receptor: Path, out_pdbqt: Path) -> Tuple[bool, str]:
    meeko = choose_binary(["mk_prepare_receptor.exe", "mk_prepare_receptor.py", "mk_prepare_receptor"])
    if meeko:
        prep_input = receptor
        if receptor.suffix.lower() == ".cif":
            # Meeko's --read_pdb route is more reproducible here than relying
            # on an optional ProDy mmCIF parser. Preserve the original mmCIF
            # and create a local conversion only for receptor preparation.
            prep_input = receptor.with_suffix(".converted.pdb")
            try:
                structure = MMCIFParser(QUIET=True).get_structure(receptor.stem, str(receptor))
                io = PDBIO()
                io.set_structure(structure)
                io.save(str(prep_input))
            except Exception:
                prep_input = receptor
        if prep_input.suffix.lower() == ".pdb":
            cmd = [
                meeko, "--read_pdb", str(prep_input), "-p", str(out_pdbqt),
                "--delete_bad_res", "--default_altloc", "A",
            ]
        else:
            cmd = [
                meeko, "-i", str(prep_input), "-p", str(out_pdbqt),
                "--delete_bad_res", "--default_altloc", "A",
            ]
        cp = run(cmd, check=False)
        if cp.returncode == 0 and out_pdbqt.exists():
            return True, "Meeko mk_prepare_receptor (delete_bad_res; default_altloc=A)"
    prep = choose_binary(["prepare_receptor4.py"])
    if prep and receptor.suffix.lower() == ".pdb":
        cp = run([prep, "-r", str(receptor), "-o", str(out_pdbqt), "-A", "hydrogens"], check=False)
        if cp.returncode == 0 and out_pdbqt.exists():
            return True, "MGLTools prepare_receptor4.py"
    return False, "No working receptor-preparation executable; mmCIF may require conversion/cleanup"


def run_vina(receptor: Path, ligand: Path, box: Dict, out_pose: Path, log_path: Path,
             exhaustiveness: int, num_modes: int) -> Dict:
    vina = choose_binary(["vina.exe", "vina"])
    if not vina:
        return {"status": "skipped", "reason": "vina executable not found"}
    c = box["center"]
    s = box["size"]
    cmd = [
        vina, "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(c[0]), "--center_y", str(c[1]), "--center_z", str(c[2]),
        "--size_x", str(s[0]), "--size_y", str(s[1]), "--size_z", str(s[2]),
        "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes),
        "--out", str(out_pose),
    ]
    cp = run(cmd, check=False)
    log_path.write_text(cp.stdout + "\n" + cp.stderr, encoding="utf-8")
    best = None
    for line in cp.stdout.splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0].isdigit():
            try:
                best = float(fields[1])
                break
            except ValueError:
                pass
    return {"status": "ok" if cp.returncode == 0 else "failed", "returncode": cp.returncode, "best_affinity_kcal_mol": best}


def load_stage1_selections(path: Path) -> Dict[str, str]:
    data = load_json(path)
    # Be permissive about the exact Stage-1 schema.
    result = {}
    if isinstance(data, dict):
        for gene, val in data.items():
            if isinstance(val, str):
                result[gene] = val
            elif isinstance(val, dict):
                for key in ("pdb_id", "pdb", "structure_id", "selected_pdb"):
                    if val.get(key):
                        result[gene] = str(val[key])
                        break
    elif isinstance(data, list):
        for row in data:
            if isinstance(row, dict):
                gene = row.get("gene") or row.get("gene_symbol") or row.get("target")
                pdb = row.get("pdb_id") or row.get("pdb") or row.get("structure_id")
                if gene and pdb:
                    result[str(gene)] = str(pdb)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected", type=Path, default=ROOT / "outputs" / "structures" / "selected_structures.json")
    ap.add_argument("--config", type=Path, default=ROOT / "config.json")
    ap.add_argument("--targets", nargs="*", default=PRIMARY_TARGETS)
    ap.add_argument("--cid", type=int, default=590836)
    ap.add_argument("--conformers", type=int, default=50)
    ap.add_argument("--dry-run", action="store_true", help="Prepare/audit but do not execute Vina")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    cfg = load_json(args.config) if args.config.exists() else {}
    selected = load_stage1_selections(args.selected)
    missing = [g for g in args.targets if g not in selected]
    if missing:
        raise SystemExit(f"Missing Stage-1 selected structures for: {', '.join(missing)}")

    audit = {"timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "targets": args.targets, "events": []}
    smiles = pubchem_canonical_smiles(args.cid)
    ligand = build_ligand_ensemble(smiles, n_conf=args.conformers, seed=args.cid)
    audit["ligand"] = ligand
    ligand_pdbqt = LIGDIR / "DINP.pdbqt"
    ligand_ok, ligand_method = prepare_ligand_pdbqt(Path(ligand["sdf"]), ligand_pdbqt)
    audit["ligand_pdbqt"] = {"ok": ligand_ok, "method": ligand_method, "path": str(ligand_pdbqt)}

    rows = []
    pockets = {}
    for gene in args.targets:
        pdb_id = selected[gene]
        rec_path, fmt = fetch_rcsb_structure(pdb_id)
        structure = parse_structure(rec_path, fmt)
        lig_hits = candidate_ligands(structure)
        pocket = box_from_ligand(lig_hits[0]) if lig_hits else None
        pockets[gene] = {"pdb_id": pdb_id, "format": fmt, "ligand_candidates": lig_hits[:10], "proposed_box": pocket}

        target_dir = DOCKDIR / gene
        target_dir.mkdir(parents=True, exist_ok=True)
        receptor_pdbqt = target_dir / f"{gene}_{pdb_id}_receptor.pdbqt"
        rec_ok, rec_method = prepare_receptor_pdbqt(rec_path, receptor_pdbqt)
        result = {"status": "not_run"}
        if pocket is None:
            result = {"status": "skipped", "reason": "no co-crystal-like non-polymer ligand detected; manual pocket required"}
        elif not ligand_ok:
            result = {"status": "skipped", "reason": ligand_method}
        elif not rec_ok:
            result = {"status": "skipped", "reason": rec_method}
        elif args.dry_run:
            result = {"status": "skipped", "reason": "--dry-run"}
        else:
            result = run_vina(
                receptor_pdbqt,
                ligand_pdbqt,
                pocket,
                target_dir / f"DINP_{gene}_vina.pdbqt",
                target_dir / "vina.log",
                int(cfg.get("vina_exhaustiveness", 32)),
                int(cfg.get("vina_num_modes", 20)),
            )

        rows.append({
            "gene": gene,
            "pdb_id": pdb_id,
            "format": fmt,
            "pocket_source": pocket["source_ligand"]["resname"] if pocket else "",
            "pocket_chain": pocket["source_ligand"]["chain"] if pocket else "",
            "receptor_prep_ok": rec_ok,
            "receptor_prep_method": rec_method,
            "ligand_prep_ok": ligand_ok,
            "docking_status": result.get("status"),
            "best_affinity_kcal_mol": result.get("best_affinity_kcal_mol"),
            "note": result.get("reason", ""),
        })

    save_json(pockets, OUT / "pocket_audit.json")
    save_json(audit, OUT / "audit.json")
    with (OUT / "docking_summary.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # Select structurally successful targets only; do not auto-promote by affinity alone.
    ok = [r for r in rows if r["docking_status"] == "ok" and r["best_affinity_kcal_mol"] is not None]
    ok.sort(key=lambda r: float(r["best_affinity_kcal_mol"]))
    md_candidates = ok[:2]
    save_json(md_candidates, OUT / "md_candidate_shortlist.json")

    lines = [
        "# DINP–CRC Stage-2 docking summary",
        "",
        f"DINP PubChem CID: {args.cid}",
        f"RDKit conformers generated: {ligand['conformers_generated']}",
        f"Ligand minimization: {ligand['force_field']}",
        "",
        "Docking is performed only when a co-crystal-like pocket is detected and both ligand and receptor PDBQT preparation succeed.",
        "Targets without an auditable pocket are intentionally skipped rather than subjected to blind docking.",
        "",
        "| Target | PDB | Pocket ligand | Docking | Best Vina affinity (kcal/mol) |",
        "|---|---|---|---|---:|",
    ]
    for r in rows:
        aff = "" if r["best_affinity_kcal_mol"] is None else str(r["best_affinity_kcal_mol"])
        lines.append(f"| {r['gene']} | {r['pdb_id']} | {r['pocket_source']} | {r['docking_status']} | {aff} |")
    lines += [
        "",
        "The MD shortlist is a computational prioritization output only. Manual inspection of chain identity, receptor completeness, pocket occupancy, protonation, cofactors and membrane context is required before production MD.",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Stage 2 complete. Outputs: {OUT}")


if __name__ == "__main__":
    main()
