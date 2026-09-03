#!/usr/bin/env python3
"""Production launcher for DINP docking on the local D: runtime.

This wrapper does NOT change the scientific logic in run_stage2_docking.py. It only:
  1) discovers the local runtime installed under D:/whynot17/dinp_stage2_runtime;
  2) prepends likely executable directories to PATH;
  3) verifies that Vina and at least one ligand/receptor preparation route are visible;
  4) runs Stage-2 without --dry-run;
  5) audits the produced docking summary and fails loudly if no real affinity was produced.

No affinity or MD shortlist is fabricated when external tools fail.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
STAGE2 = ROOT / "run_stage2_docking.py"
OUT = ROOT / "outputs" / "stage2_docking"
DEFAULT_RUNTIME = Path(r"D:\whynot17\dinp_stage2_runtime")


def add_runtime_to_path(runtime: Path) -> list[str]:
    candidates = [
        runtime,
        runtime / "Scripts",
        runtime / "bin",
        runtime / "Library" / "bin",
        runtime / "Library" / "usr" / "bin",
    ]
    existing = [str(p) for p in candidates if p.exists()]
    os.environ["PATH"] = os.pathsep.join(existing + [os.environ.get("PATH", "")])
    return existing


def which_any(names: list[str]) -> str | None:
    for name in names:
        hit = shutil.which(name)
        if hit:
            return hit
    return None


def tool_audit() -> dict:
    return {
        "python": sys.executable,
        "vina": which_any(["vina.exe", "vina"]),
        "obabel": which_any(["obabel.exe", "obabel", "babel.exe", "babel"]),
        "mk_prepare_ligand": which_any([
            "mk_prepare_ligand.exe", "mk_prepare_ligand.py", "mk_prepare_ligand"
        ]),
        "mk_prepare_receptor": which_any([
            "mk_prepare_receptor.exe", "mk_prepare_receptor.py", "mk_prepare_receptor"
        ]),
        "prepare_receptor4": which_any([
            "prepare_receptor4.py", "prepare_receptor4"
        ]),
    }


def read_real_docking_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    rows = []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            status = str(row.get("status", "")).strip().lower()
            affinity_raw = str(
                row.get("best_affinity_kcal_mol", row.get("best_vina_affinity_kcal_mol", ""))
            ).strip()
            try:
                affinity = float(affinity_raw)
            except Exception:
                affinity = None
            if status == "ok" and affinity is not None:
                row["_affinity"] = affinity
                rows.append(row)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runtime", type=Path, default=DEFAULT_RUNTIME)
    ap.add_argument("--targets", nargs="*", default=["PTGER4", "CXCR4", "MMP9", "STAT3"])
    ap.add_argument("--conformers", type=int, default=50)
    ap.add_argument("--selected", type=Path, default=ROOT / "outputs" / "structures" / "selected_structures.json")
    ap.add_argument("--config", type=Path, default=ROOT / "config.json")
    args = ap.parse_args()

    if not STAGE2.exists():
        raise SystemExit(f"Missing Stage-2 script: {STAGE2}")
    if not args.runtime.exists():
        raise SystemExit(f"Runtime directory not found: {args.runtime}")

    path_entries = add_runtime_to_path(args.runtime)
    tools = tool_audit()

    # Production docking requires Vina and at least one ligand preparation route.
    if not tools["vina"]:
        raise SystemExit("AutoDock Vina is not visible on PATH after adding the runtime directory.")
    if not (tools["mk_prepare_ligand"] or tools["obabel"]):
        raise SystemExit("Neither Meeko ligand preparation nor OpenBabel is visible on PATH.")
    if not (tools["mk_prepare_receptor"] or tools["prepare_receptor4"]):
        raise SystemExit("No receptor preparation executable is visible on PATH.")

    OUT.mkdir(parents=True, exist_ok=True)
    audit_path = OUT / "production_runtime_audit.json"
    audit = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime": str(args.runtime),
        "path_entries_added": path_entries,
        "tools": tools,
        "targets": args.targets,
        "conformers": args.conformers,
        "stage2_script": str(STAGE2),
        "command": None,
        "returncode": None,
        "real_docking_rows": 0,
    }

    cmd = [
        sys.executable,
        str(STAGE2),
        "--selected", str(args.selected),
        "--config", str(args.config),
        "--conformers", str(args.conformers),
        "--targets", *args.targets,
    ]
    audit["command"] = cmd

    cp = subprocess.run(cmd, cwd=str(ROOT), text=True)
    audit["returncode"] = cp.returncode

    summary_csv = OUT / "docking_summary.csv"
    real_rows = read_real_docking_rows(summary_csv)
    audit["real_docking_rows"] = len(real_rows)
    audit["real_docking"] = [
        {
            "target": r.get("target") or r.get("gene"),
            "pdb": r.get("pdb_id") or r.get("pdb"),
            "best_affinity_kcal_mol": r.get("_affinity"),
        }
        for r in sorted(real_rows, key=lambda x: x["_affinity"])
    ]
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")

    if cp.returncode != 0:
        raise SystemExit(cp.returncode)
    if not real_rows:
        raise SystemExit(
            "Stage-2 completed but produced no real Vina affinities. "
            f"Inspect {summary_csv} and {audit_path}; do not create an MD shortlist yet."
        )

    print("\nReal docking completed:")
    for row in sorted(real_rows, key=lambda x: x["_affinity"]):
        target = row.get("target") or row.get("gene") or "unknown"
        print(f"  {target}: {row['_affinity']:.3f} kcal/mol")
    print(f"\nRuntime audit: {audit_path}")
    print(f"Docking summary: {summary_csv}")


if __name__ == "__main__":
    main()
