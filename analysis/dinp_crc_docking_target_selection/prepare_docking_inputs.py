from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem
from openmm.app import PDBFile, PDBxFile


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUT = HERE / "inputs"
PREP = HERE / "prepared"
VINA = ROOT / "work" / "bin" / "vina_1.2.7_win.exe"
MEEKO_RECEPTOR = ROOT / "work" / ".venv_docking" / "Scripts" / "mk_prepare_receptor.exe"
MEEKO_LIGAND = ROOT / "work" / ".venv_docking" / "Scripts" / "mk_prepare_ligand.exe"


TARGETS = {
    "PPARG": {"pdb": "4XLD", "co_ligand": "BRL", "box_size": [24, 24, 24]},
    "PPARA": {"pdb": "6KB1", "co_ligand": "T4T", "box_size": [24, 24, 24]},
    "RXRA": {"pdb": "1FBY", "co_ligand": "9CR", "box_size": [22, 22, 22]},
    "NR1I2": {"pdb": "8CH8", "co_ligand": "ULC", "box_size": [24, 24, 24]},
    "ESR1": {"pdb": "1XPC", "co_ligand": "AIT", "box_size": [22, 22, 22]},
}


def ligand_center(pdb_path: Path, residue_name: str):
    xyz = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == residue_name:
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not xyz:
        raise RuntimeError(f"No co-crystallized ligand {residue_name} in {pdb_path}")
    return np.asarray(xyz, dtype=float).mean(axis=0)


def prepare_ligand():
    source = INPUT / "DINP_CID590836_2d.sdf"
    supplier = Chem.SDMolSupplier(str(source), removeHs=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise RuntimeError("Could not read DINP PubChem SDF")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260909
    params.pruneRmsThresh = 0.5
    params.numThreads = 0
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params))
    if not conf_ids:
        raise RuntimeError("RDKit could not generate a 3-D DINP conformer")
    energies = []
    for conf_id in conf_ids:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        ff.Initialize()
        ff.Minimize(maxIts=2000)
        energies.append((float(ff.CalcEnergy()), conf_id))
    energy, chosen = min(energies)
    best_conf = Chem.Conformer(mol.GetConformer(chosen))
    mol.RemoveAllConformers()
    mol.AddConformer(best_conf, assignId=True)
    mol.SetProp("PUBCHEM_CID", "590836")
    mol.SetProp("CAS", "28553-12-0; 68515-48-0")
    mol.SetProp("CONFORMER_METHOD", "RDKit ETKDGv3 + UFF; representative single conformer")
    mol.SetProp("CONFORMER_SCREEN_LOWEST_UFF_ENERGY", f"{energy:.6f}")
    mol.SetProp("ISOMER_SCOPE", "representative PubChem bis(7-methyloctyl) phthalate; DINP commercial mixture not fully represented")
    sdf_out = PREP / "DINP_CID590836_representative_3d.sdf"
    pdb_out = PREP / "DINP_CID590836_representative_3d.pdb"
    writer = Chem.SDWriter(str(sdf_out))
    writer.write(mol)
    writer.close()
    Chem.MolToPDBFile(mol, str(pdb_out))
    return {"sdf": str(sdf_out), "pdb": str(pdb_out), "n_atoms_with_h": mol.GetNumAtoms(), "conformer_energy_uff": float(energy)}


def prepare_receptor(gene, spec):
    raw = INPUT / f"{spec['pdb']}.pdb"
    fixer = PDBFixer(filename=str(raw))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    # Do not add unresolved terminal/missing loop residues automatically; the
    # deposited LBD is used as the structural template for docking.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # Leave protonation to Meeko. PDBFixer's explicit hydrogen naming can be
    # misread for a few partially resolved LBD residues by receptor perception.
    fixed = PREP / f"{gene}_{spec['pdb']}_protein_h.pdb"
    with fixed.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    center = ligand_center(raw, spec["co_ligand"])
    base = PREP / f"{gene}_{spec['pdb']}"
    command = [
            str(MEEKO_RECEPTOR),
            "--read_pdb",
            str(fixed),
            "-o",
            str(base),
            "--write_pdbqt",
            str(base) + "_rigid.pdbqt",
            "--box_center",
            *(f"{x:.3f}" for x in center),
            "--box_size",
            *(str(int(x)) for x in spec["box_size"]),
            "--write_vina_box",
            str(base) + "_box.txt",
            "--delete_bad_res",
            "--forgive_extra_bonds",
            "--default_altloc",
            "A",
        ]
    receptor_input = fixed
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError:
        # Some deposited LBDs contain locally unresolved/chemically ambiguous
        # residues that PDBFixer reconstructs in a way Meeko cannot perceive.
        # For support docking only, fall back to the deposited PDB and retain
        # the warning in the manifest; the primary PPARG preparation must not
        # use this fallback silently.
        if gene == "PPARG":
            raise
        command[command.index(str(fixed))] = str(raw)
        receptor_input = raw
        subprocess.run(command, check=True)
    return {
        "gene_symbol": gene,
        "pdb_id": spec["pdb"],
        "co_crystal_ligand": spec["co_ligand"],
        "protein_pdb": str(receptor_input),
        "receptor_pdbqt": str(base) + "_rigid.pdbqt",
        "box_file": str(base) + "_box.txt",
        "box_center_A": [round(float(x), 3) for x in center],
        "box_size_A": spec["box_size"],
    }


def main():
    PREP.mkdir(parents=True, exist_ok=True)
    if not MEEKO_RECEPTOR.exists() or not MEEKO_LIGAND.exists():
        raise FileNotFoundError("Meeko executables are not installed")
    ligand = prepare_ligand()
    subprocess.run([str(MEEKO_LIGAND), "-i", ligand["sdf"], "-o", str(PREP / "DINP_CID590836.pdbqt"), "--add_index_map"], check=True)
    receptors = [prepare_receptor(gene, spec) for gene, spec in TARGETS.items()]
    manifest = {"ligand": ligand, "receptors": receptors, "status": "prepared_for_docking", "cebp_b_excluded": True}
    (PREP / "docking_input_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
