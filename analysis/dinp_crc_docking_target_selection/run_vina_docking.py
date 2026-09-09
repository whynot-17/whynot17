from __future__ import annotations

import csv
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
PREP = HERE / "prepared"
OUT = HERE / "docking_runs"
VINA = Path(__file__).resolve().parents[2] / "work" / "bin" / "vina_1.2.7_win.exe"
LIGAND = PREP / "DINP_CID590836.pdbqt"

RUNS = [
    ("PPARG", "4XLD", 20260909),
    ("PPARG", "4XLD", 20260910),
    ("PPARG", "4XLD", 20260911),
    ("PPARA", "6KB1", 20260909),
    ("RXRA", "1FBY", 20260909),
    ("NR1I2", "8CH8", 20260909),
    ("ESR1", "1XPC", 20260909),
]


def parse_modes(text):
    rows = []
    # Vina table: mode | affinity | rmsd l.b. | rmsd u.b.
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
                    "rmsd_lb_A": float(match.group(3)),
                    "rmsd_ub_A": float(match.group(4)),
                }
            )
        elif rows and line.strip() and not re.match(r"^\s*\d+", line):
            break
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not VINA.exists():
        raise FileNotFoundError(VINA)
    if not LIGAND.exists():
        raise FileNotFoundError(LIGAND)
    result_rows = []
    run_rows = []
    for gene, pdb_id, seed in RUNS:
        receptor = PREP / f"{gene}_{pdb_id}_rigid.pdbqt"
        config = PREP / f"{gene}_{pdb_id}_box.txt"
        stem = f"{gene}_{pdb_id}_DINP_seed{seed}"
        out_pdbqt = OUT / f"{stem}.pdbqt"
        log_path = OUT / f"{stem}.log"
        if not receptor.exists() or not config.exists():
            run_rows.append({"gene_symbol": gene, "pdb_id": pdb_id, "seed": seed, "status": "missing_input", "best_affinity_kcal_mol": ""})
            continue
        if log_path.exists() and out_pdbqt.exists():
            log_text = log_path.read_text(encoding="utf-8", errors="replace")
            completed_returncode = 0
        else:
            command = [
            str(VINA),
            "--config", str(config),
            "--receptor", str(receptor),
            "--ligand", str(LIGAND),
            "--out", str(out_pdbqt),
            "--exhaustiveness", "32",
            "--num_modes", "20",
            "--energy_range", "5",
            "--cpu", "8",
            "--seed", str(seed),
            "--verbosity", "1",
            ]
            completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace")
            log_text = (completed.stdout or "") + "\n" + (completed.stderr or "")
            log_path.write_text(log_text, encoding="utf-8")
            completed_returncode = completed.returncode
        modes = parse_modes(log_text)
        status = "ok" if completed_returncode == 0 and modes else f"failed_exit_{completed_returncode}"
        best = modes[0]["affinity_kcal_mol"] if modes else ""
        run_rows.append({"gene_symbol": gene, "pdb_id": pdb_id, "seed": seed, "status": status, "best_affinity_kcal_mol": best})
        for mode in modes:
            result_rows.append({"gene_symbol": gene, "pdb_id": pdb_id, "seed": seed, "output_pdbqt": str(out_pdbqt), "log_file": str(log_path), **mode})
    run_fields = ["gene_symbol", "pdb_id", "seed", "status", "best_affinity_kcal_mol"]
    with (OUT / "DINP_CRC_docking_run_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=run_fields)
        writer.writeheader()
        writer.writerows(run_rows)
    result_fields = ["gene_symbol", "pdb_id", "seed", "output_pdbqt", "log_file", "mode", "affinity_kcal_mol", "rmsd_lb_A", "rmsd_ub_A"]
    with (OUT / "DINP_CRC_docking_modes.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=result_fields)
        writer.writeheader()
        writer.writerows(result_rows)
    (OUT / "docking_execution_metadata.txt").write_text(
        f"generated_utc={datetime.now(timezone.utc).isoformat()}\n"
        "engine=AutoDock Vina 1.2.7 Windows binary\n"
        "ligand=PubChem CID 590836 representative DINP conformer\n"
        "exhaustiveness=32\nnum_modes=20\nenergy_range_kcal_mol=5\ncpu=8\n",
        encoding="utf-8",
    )
    print(f"Completed {len(run_rows)} docking runs; successful={sum(r['status'] == 'ok' for r in run_rows)}")


if __name__ == "__main__":
    main()
