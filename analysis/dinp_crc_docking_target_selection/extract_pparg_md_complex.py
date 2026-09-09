from __future__ import annotations

from pathlib import Path


HERE = Path(__file__).resolve().parent
PREP = HERE / "prepared"
RUNS = HERE / "docking_runs"


def element_from_pdbqt(atom_type: str) -> str:
    atom_type = atom_type.upper()
    if atom_type in {"A", "C", "CA"}:
        return "C"
    if atom_type.startswith("N"):
        return "N"
    if atom_type.startswith("O"):
        return "O"
    if atom_type.startswith("S"):
        return "S"
    if atom_type.startswith("P"):
        return "P"
    if atom_type.startswith("CL"):
        return "Cl"
    if atom_type.startswith("F"):
        return "F"
    if atom_type.startswith("BR"):
        return "Br"
    if atom_type.startswith("I"):
        return "I"
    if atom_type.startswith("H"):
        return "H"
    return atom_type[:2].title()


def write_complex(seed: int):
    pose = RUNS / f"PPARG_4XLD_DINP_seed{seed}.pdbqt"
    receptor = PREP / "PPARG_4XLD_protein_h.pdb"
    output = PREP / f"PPARG_4XLD_DINP_seed{seed}_mode1_complex.pdb"
    lines = []
    for line in receptor.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "TER")):
            lines.append(line)
    model = False
    atom_index = sum(line.startswith("ATOM") for line in lines) + 1
    ligand_atom_index = 1
    for line in pose.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            model = line[10:14].strip() in {"", "1"}
            continue
        if line.startswith("ENDMDL") and model:
            break
        if not model or not line.startswith("ATOM"):
            continue
        element = element_from_pdbqt(line.split()[-1])
        out = list(line.ljust(80))
        out[0:6] = list("HETATM")
        out[6:11] = list(f"{atom_index:5d}")
        ligand_atom_name = f"{element:>2s}{ligand_atom_index:02d}"
        out[12:16] = list(ligand_atom_name[:4])
        out[17:20] = list("DNP")
        out[21] = "Z"
        out[22:26] = list("   1")
        out[76:78] = list(f"{element:>2s}")
        lines.append("".join(out).rstrip())
        atom_index += 1
        ligand_atom_index += 1
    lines.extend(["TER", "END"])
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def main():
    for seed in (20260909, 20260911):
        print(write_complex(seed))


if __name__ == "__main__":
    main()
