"""Export the completed MiNP 100 ns analysis tables as GROMACS-style XVG files."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(r"E:\chatgpt\whynot17\analysis\pparg_minp_md_100ns")
OUT = ROOT / "outputs"


def read_csv(name):
    with (OUT / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_time_xvg(path, title, y_label, legends, rows, columns, scale=1.0):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# MiNP–PPARG completed 100 ns MD; seed20260917; reference=initial PDB\n")
        handle.write("# Time is in ns.\n")
        handle.write(f'@ title "{title}"\n')
        handle.write('@ xaxis label "Time (ns)"\n')
        handle.write(f'@ yaxis label "{y_label}"\n')
        for index, legend in enumerate(legends):
            handle.write(f'@ s{index} legend "{legend}"\n')
        for row in rows:
            values = [float(row[column]) * scale for column in columns]
            handle.write(f"{float(row['time_ns']):.5f} " + " ".join(f"{x:.8f}" for x in values) + "\n")


def main():
    metrics = read_csv("trajectory_metrics_0_100ns_complete.csv")
    sasa = read_csv("rg_sasa_0_100ns.csv")
    contacts = read_csv("contact_samples_0_100ns.csv")
    if len(metrics) != 9992 or len(sasa) != 1000 or len(contacts) != 1000:
        raise RuntimeError(f"Unexpected row counts: RMSD={len(metrics)}, SASA/Rg={len(sasa)}, H-bond={len(contacts)}")

    # Add the derived complex RMSD to the same CSV rows without recomputing the trajectory.
    # This is the atom-count-weighted RMSD after the same protein-backbone fit.
    complex_rows = []
    for row in metrics:
        n_protein, n_ligand = 2112.0, 21.0
        complex_A = ((n_protein * float(row["protein_backbone_rmsd_A"]) ** 2 + n_ligand * float(row["whole_ligand_rmsd_A"]) ** 2) / (n_protein + n_ligand)) ** 0.5
        copy = dict(row)
        copy["_complex_rmsd_A"] = complex_A
        complex_rows.append(copy)
    write_time_xvg(
        OUT / "MiNP_rmsd_0_100ns.xvg",
        "MiNP–PPARG RMSD components",
        "RMSD (nm)",
        ["Complex heavy atoms", "Protein heavy atoms", "Whole ligand", "Core", "Branch A", "Branch B"],
        complex_rows,
        ["_complex_rmsd_A", "protein_backbone_rmsd_A", "whole_ligand_rmsd_A", "core_rmsd_A", "branch_A_rmsd_A", "branch_B_rmsd_A"],
        scale=0.1,
    )

    write_time_xvg(
        OUT / "MiNP_sasa_0_100ns.xvg",
        "MiNP–PPARG SASA",
        "SASA (nm^2)",
        ["Protein SASA in complex", "Protein isolated SASA", "Ligand SASA in complex", "Ligand isolated SASA", "Complex SASA", "Ligand buried SASA"],
        sasa,
        ["protein_sasa_complex_A2", "protein_sasa_isolated_A2", "ligand_sasa_complex_A2", "ligand_sasa_isolated_A2", "complex_sasa_A2", "ligand_buried_A2"],
        scale=0.01,
    )
    write_time_xvg(
        OUT / "MiNP_rg_0_100ns.xvg",
        "MiNP–PPARG radius of gyration",
        "Rg (nm)",
        ["Protein", "Ligand"],
        sasa,
        ["protein_rg_A", "ligand_rg_A"],
        scale=0.1,
    )
    write_time_xvg(
        OUT / "MiNP_hbond_0_100ns.xvg",
        "MiNP–PPARG hydrogen-bond proxy",
        "Count",
        ["Geometric H-bond proxy"],
        contacts,
        ["hbond_count"],
        scale=1.0,
    )

    ligand_rmsf = read_csv("ligand_heavy_rmsf.csv")
    ligand_internal = read_csv("ligand_internal_rmsf.csv")
    with (OUT / "MiNP_ligand_rmsf_0_100ns.xvg").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# MiNP–PPARG completed 100 ns MD; per ligand heavy atom; RMSF in nm\n")
        handle.write('@ title "MiNP ligand RMSF"\n@ xaxis label "Ligand heavy-atom index"\n@ yaxis label "RMSF (nm)"\n')
        handle.write('@ s0 legend "Protein-aligned RMSF"\n@ s1 legend "Core-internal RMSF"\n')
        for row, internal in zip(ligand_rmsf, ligand_internal):
            handle.write(f"{int(row['atom_index_1based'])} {float(row['rmsf_A'])/10.0:.8f} {float(internal['internal_rmsf_A'])/10.0:.8f} # {row['atom_name']} {row['part']}\n")

    pocket_rmsf = read_csv("pocket_ca_rmsf.csv")
    with (OUT / "MiNP_pocket_ca_rmsf_0_100ns.xvg").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# MiNP–PPARG completed 100 ns MD; initial-pocket C-alpha RMSF; RMSF in nm\n")
        handle.write('@ title "MiNP initial-pocket C-alpha RMSF"\n@ xaxis label "Pocket residue index"\n@ yaxis label "RMSF (nm)"\n@ s0 legend "Pocket C-alpha RMSF"\n')
        for index, row in enumerate(pocket_rmsf, start=1):
            handle.write(f"{index} {float(row['rmsf_A'])/10.0:.8f} # {row['residue']}\n")

    manifest = {
        "source": "completed MiNP seed20260917 100 ns analysis tables",
        "outputs": [
            "MiNP_rmsd_0_100ns.xvg", "MiNP_ligand_rmsf_0_100ns.xvg", "MiNP_pocket_ca_rmsf_0_100ns.xvg",
            "MiNP_sasa_0_100ns.xvg", "MiNP_rg_0_100ns.xvg", "MiNP_hbond_0_100ns.xvg",
        ],
        "units": {"time": "ns", "rmsd": "nm", "rmsf": "nm", "sasa": "nm^2", "rg": "nm", "hbond": "count"},
    }
    (OUT / "minp_metrics_xvg_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
