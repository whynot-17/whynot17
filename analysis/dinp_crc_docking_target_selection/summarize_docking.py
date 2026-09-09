from __future__ import annotations

import csv
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
PREP = HERE / "prepared"
RUNS = HERE / "docking_runs"


def read_csv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def mean_sd(values):
    if not values:
        return "", ""
    mean = sum(values) / len(values)
    sd = math.sqrt(sum((x - mean) ** 2 for x in values) / len(values)) if len(values) > 1 else 0.0
    return round(mean, 4), round(sd, 4)


def parse_ligand_pose(path: Path):
    atoms = []
    in_model = False
    saw_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            saw_model = True
            in_model = line[10:14].strip() in {"", "1"}
            continue
        if saw_model and line.startswith("ENDMDL"):
            break
        if (not saw_model or in_model) and line.startswith(("ATOM", "HETATM")):
            atom_type = line.split()[-1]
            if atom_type.upper().startswith("H"):
                continue
            atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms


def parse_receptor_atoms(path: Path):
    residues = defaultdict(list)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        element = line[76:78].strip() or atom_name[:1]
        if element.upper().startswith("H"):
            continue
        key = f"{line[17:20].strip()}{line[21].strip()}:{line[22:26].strip()}"
        residues[key].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return residues


def contact_audit(ligand_path: Path, receptor_path: Path):
    ligand = parse_ligand_pose(ligand_path)
    receptor = parse_receptor_atoms(receptor_path)
    contacts = []
    for residue, atoms in receptor.items():
        minimum = min(
            math.sqrt((lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2)
            for lx, ly, lz in ligand
            for rx, ry, rz in atoms
        )
        if minimum <= 4.5:
            contacts.append((residue, round(minimum, 3)))
    return sorted(contacts, key=lambda item: item[1])


def rmsd(a, b):
    if len(a) != len(b) or not a:
        return ""
    return round(math.sqrt(sum((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2 for (x1, y1, z1), (x2, y2, z2) in zip(a, b)) / len(a)), 4)


def main():
    runs = read_csv(RUNS / "DINP_CRC_docking_run_summary.csv")
    modes = read_csv(RUNS / "DINP_CRC_docking_modes.csv")
    by_target = defaultdict(list)
    for row in runs:
        if row["status"] == "ok":
            by_target[row["gene_symbol"]].append(float(row["best_affinity_kcal_mol"]))
    summary = []
    for gene, values in sorted(by_target.items(), key=lambda kv: min(kv[1])):
        m, sd = mean_sd(values)
        summary.append(
            {
                "gene_symbol": gene,
                "n_successful_runs": len(values),
                "best_affinity_kcal_mol": min(values),
                "mean_top1_affinity_kcal_mol": m,
                "sd_top1_affinity_kcal_mol": sd,
                "top1_range_kcal_mol": round(max(values) - min(values), 4),
                "interpretation": "descriptive only; Vina scores are not directly comparable across different proteins/pockets",
            }
        )
    with (RUNS / "DINP_CRC_docking_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = list(summary[0].keys()) if summary else []
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)

    # Contact audit uses the top pose from each PPARG seed. The receptor is
    # fixed, so coordinate-frame comparison is valid for this same-structure
    # seed sensitivity check.
    contact_rows = []
    pparg_pose_paths = []
    for row in runs:
        if row["gene_symbol"] != "PPARG" or row["status"] != "ok":
            continue
        pose = RUNS / f"PPARG_4XLD_DINP_seed{row['seed']}.pdbqt"
        pparg_pose_paths.append(pose)
        receptor = PREP / "PPARG_4XLD_protein_h.pdb"
        for residue, distance in contact_audit(pose, receptor):
            contact_rows.append({"gene_symbol": "PPARG", "pdb_id": "4XLD", "seed": row["seed"], "residue": residue, "min_heavy_atom_distance_A": distance})
    with (RUNS / "PPARG_DINP_top_pose_contact_audit.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["gene_symbol", "pdb_id", "seed", "residue", "min_heavy_atom_distance_A"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(contact_rows)

    poses = [parse_ligand_pose(path) for path in pparg_pose_paths]
    pairwise = []
    for i in range(len(poses)):
        for j in range(i + 1, len(poses)):
            pairwise.append({"seed_a": pparg_pose_paths[i].stem.split("seed")[-1], "seed_b": pparg_pose_paths[j].stem.split("seed")[-1], "ligand_pose_rmsd_A": rmsd(poses[i], poses[j])})
    with (RUNS / "PPARG_DINP_seed_pose_rmsd.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["seed_a", "seed_b", "ligand_pose_rmsd_A"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(pairwise)

    report = "# DINP docking execution summary\n\n"
    report += f"Generated: {datetime.now(timezone.utc).isoformat()}\n\n"
    report += "## Target decision\n\n"
    report += "PPARG remains the only primary target for downstream 100 ns MD. CEBPB was not docked. PPARA, RXRA, NR1I2 and ESR1 were run as support docking candidates.\n\n"
    report += "## Vina results\n\n"
    report += "| Target | Successful runs | Best top-1 score (kcal/mol) | Mean top-1 | SD |\n|---|---:|---:|---:|---:|\n"
    for row in summary:
        report += f"| {row['gene_symbol']} | {row['n_successful_runs']} | {row['best_affinity_kcal_mol']} | {row['mean_top1_affinity_kcal_mol']} | {row['sd_top1_affinity_kcal_mol']} |\n"
    report += "\n"
    report += "PPARG top-1 scores across three seeds were tightly clustered, which supports sampling repeatability for this fixed receptor/pocket setup. Cross-target ranking is descriptive only: different proteins, pocket sizes, protonation states and structural preparations make raw Vina energies non-comparable. The more negative ESR1 score therefore does not override the PPARG mechanistic selection.\n\n"
    report += "## Pose/contact audit\n\n"
    if pairwise:
        report += "PPARG same-structure seed pose RMSDs (unrestrained coordinate RMSD of the top pose) were: " + ", ".join(f"{r['seed_a']} vs {r['seed_b']}: {r['ligand_pose_rmsd_A']} Å" for r in pairwise) + ".\n\n"
    top_contacts = defaultdict(list)
    for row in contact_rows:
        top_contacts[row["residue"]].append(float(row["min_heavy_atom_distance_A"]))
    ranked_contacts = sorted(((res, min(ds), len(ds)) for res, ds in top_contacts.items()), key=lambda x: x[1])[:20]
    report += "Closest PPARG top-pose protein contacts (≤4.5 Å in at least one seed):\n\n"
    report += "| Residue | Minimum distance (Å) | Seeds with contact |\n|---|---:|---:|\n"
    for residue, distance, n in ranked_contacts:
        report += f"| {residue} | {distance:.3f} | {n} |\n"
    report += "\n"
    report += "## Interpretation boundary\n\n"
    report += "Docking supplies a pose-plausibility hypothesis only. It does not prove that DINP binds PPARG, does not resolve the commercial DINP isomer mixture, and does not prove PPARG causes CEBPB induction. MD is intentionally deferred to an external machine; before starting there, inspect the PPARG poses and validate the DINP ligand parameters.\n"
    (HERE / "DINP_CRC_docking_results_report.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
