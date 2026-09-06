"""Staged recovery probes for CELLxGENE Census raw expression.

The script intentionally separates runtime recovery from Phase 2G inference:

1. one known paired donor and the frozen nine-gene minimal set;
2. the same donor and the frozen Phase 2G target universe;
3. only after 1 and 2 succeed, one donor at a time for the complete paired
   epithelial cohort.

Every raw-expression table is written immediately to a cache directory. A
failed later table or donor therefore cannot erase earlier successful work.
Run in the pinned WSL environment.
"""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from contextlib import nullcontext
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs"
CACHE = ROOT / "work" / "mcop_phase2g" / "census_staged_cache"
CENSUS_VERSION = "2025-11-08"
CENSUS_URI = "s3://cellxgene-census-public-us-west-2/cell-census/2025-11-08/soma/"
ORGANISM = "homo_sapiens"
DATASET_ID = "16023185-de21-4c0d-a9c8-73abdd52d142"
DEFAULT_DONOR = "C106"
MINIMAL_GENES = ["PPARA", "PPARD", "PPARG", "NR1I2", "NR1I3", "NR1H2", "NR1H3", "RELA", "STAT3"]


def quote(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def in_filter(field: str, values: list[str]) -> str:
    return f"{field} in [{', '.join(quote(x) for x in values)}]"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def manifest_path(stage: int, donor: str) -> Path:
    return CACHE / f"stage{stage}_{donor}_manifest.json"


def load_target_universe() -> list[str]:
    sys.path.insert(0, str(ROOT / "work" / "scripts"))
    import mcop_phase2g_epithelial_state_analysis as phase2g

    state_sets, _ = phase2g.build_state_sets()
    genes = set().union(*state_sets.values()) | set(phase2g.PPAR_NR_GENES) | set(phase2g.ANCHOR_CANDIDATES)
    network_path = ROOT / "work" / "mcop_phase2g" / "dorothea_raw.tsv"
    if network_path.exists():
        network = pd.read_csv(network_path, sep="\t")
        if "target_genesymbol" in network.columns:
            genes |= set(network["target_genesymbol"].dropna().astype(str))
    return sorted(str(g) for g in genes if str(g) and str(g) != "nan")


def paired_donors() -> list[str]:
    path = OUTPUT / "mcop_phase2f_singlecell_donor_scores.csv"
    scores = pd.read_csv(path, usecols=["dataset_id", "donor_id", "group", "compartment"])
    scores = scores.loc[
        scores["dataset_id"].astype(str).eq(DATASET_ID)
        & scores["compartment"].eq("epithelial")
    ]
    tumor = set(scores.loc[scores["group"].eq("tumor"), "donor_id"].astype(str))
    normal = set(scores.loc[scores["group"].eq("normal"), "donor_id"].astype(str))
    return sorted(tumor & normal)


def query_filter(donor: str) -> str:
    return f"is_primary_data == True and dataset_id == {quote(DATASET_ID)} and donor_id == {quote(donor)}"


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def run_one_donor(
    stage: int,
    donor: str,
    genes: list[str],
    force: bool = False,
    census_handle=None,
) -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    manifest_file = manifest_path(stage, donor)
    if manifest_file.exists() and not force:
        cached = json.loads(manifest_file.read_text(encoding="utf-8"))
        if cached.get("status") == "success":
            print(f"[stage {stage}] cache hit donor={donor}", flush=True)
            return cached

    stage_dir = CACHE / f"stage{stage}_{donor}"
    stage_dir.mkdir(parents=True, exist_ok=True)
    filter_text = query_filter(donor)
    payload = {
        "stage": stage,
        "donor": donor,
        "dataset_id": DATASET_ID,
        "census_version": CENSUS_VERSION,
        "census_uri": CENSUS_URI,
        "value_filter": filter_text,
        "n_requested_genes": len(genes),
        "requested_genes": genes,
        "started_utc": utc_now(),
        "status": "running",
        "raw_tables_completed": 0,
        "raw_rows_cached": 0,
    }
    write_json(manifest_file, payload)
    try:
        import cellxgene_census
        import tiledbsoma as soma

        print(f"[stage {stage}] opening Census {CENSUS_VERSION}; donor={donor}; genes={len(genes)}", flush=True)
        census_context = (
            cellxgene_census.open_soma(uri=CENSUS_URI)
            if census_handle is None
            else nullcontext(census_handle)
        )
        with census_context as census:
            experiment = census["census_data"][ORGANISM]
            with experiment.axis_query(
                measurement_name="RNA",
                obs_query=soma.AxisQuery(value_filter=filter_text),
                var_query=soma.AxisQuery(value_filter=in_filter("feature_name", genes)),
            ) as query:
                obs_frames = [table.to_pandas() for table in query.obs(column_names=["dataset_id", "donor_id", "disease", "cell_type", "tissue", "is_primary_data"])]
                obs = pd.concat(obs_frames, ignore_index=True) if obs_frames else pd.DataFrame()
                obs.to_csv(stage_dir / "obs.csv", index=False)
                var_frames = [table.to_pandas() for table in query.var(column_names=["feature_name", "feature_id"])]
                var = pd.concat(var_frames, ignore_index=True) if var_frames else pd.DataFrame()
                var.to_csv(stage_dir / "var.csv", index=False)
                payload["n_obs"] = int(len(obs))
                payload["n_var_returned"] = int(len(var))
                payload["obs_joinids_n"] = int(len(query.obs_joinids()))
                payload["var_joinids_n"] = int(len(query.var_joinids()))
                write_json(manifest_file, payload)
                print(f"[stage {stage}] metadata returned obs={len(obs):,}; var={len(var):,}; entering raw tables", flush=True)
                for table_number, table in enumerate(query.X("raw").tables(), start=1):
                    frame = table.to_pandas()
                    chunk_path = stage_dir / f"raw_chunk_{table_number:05d}.csv.gz"
                    frame.to_csv(chunk_path, index=False, compression="gzip")
                    payload["raw_tables_completed"] = table_number
                    payload["raw_rows_cached"] = int(payload["raw_rows_cached"] + len(frame))
                    write_json(manifest_file, payload)
                    print(f"[stage {stage}] cached table={table_number}; rows={len(frame):,}; total={payload['raw_rows_cached']:,}", flush=True)
        payload["status"] = "success"
        payload["finished_utc"] = utc_now()
        write_json(manifest_file, payload)
        print(f"[stage {stage}] SUCCESS donor={donor}; tables={payload['raw_tables_completed']}", flush=True)
    except Exception as exc:
        payload["status"] = "failed"
        payload["finished_utc"] = utc_now()
        payload["error_type"] = type(exc).__name__
        payload["error"] = str(exc)
        payload["traceback"] = traceback.format_exc()
        write_json(manifest_file, payload)
        print(f"[stage {stage}] FAILED donor={donor}: {type(exc).__name__}: {exc}", flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=int, choices=[1, 2, 3], required=True)
    parser.add_argument("--donor", default=DEFAULT_DONOR)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    if args.stage == 1:
        result = run_one_donor(1, args.donor, MINIMAL_GENES, force=args.force)
        raise SystemExit(0 if result.get("status") == "success" else 2)
    if args.stage == 2:
        stage1 = manifest_path(1, args.donor)
        if not stage1.exists() or json.loads(stage1.read_text(encoding="utf-8")).get("status") != "success":
            raise SystemExit("Stage 2 blocked: Stage 1 has not succeeded for this donor.")
        result = run_one_donor(2, args.donor, load_target_universe(), force=args.force)
        raise SystemExit(0 if result.get("status") == "success" else 2)
    stage2 = manifest_path(2, args.donor)
    if not stage2.exists() or json.loads(stage2.read_text(encoding="utf-8")).get("status") != "success":
        raise SystemExit("Stage 3 blocked: Stage 2 has not succeeded for the probe donor.")
    donors = paired_donors()
    print(f"[stage 3] eligible paired epithelial donors={len(donors)}", flush=True)
    genes = load_target_universe()
    import cellxgene_census

    failures = []
    success_count = 0
    print("[stage 3] opening one shared Census session for all remaining donors", flush=True)
    with cellxgene_census.open_soma(uri=CENSUS_URI) as census:
        for index, donor in enumerate(donors, start=1):
            print(f"[stage 3] donor {index}/{len(donors)}: {donor}", flush=True)
            result = run_one_donor(
                3,
                donor,
                genes,
                force=args.force,
                census_handle=census,
            )
            if result.get("status") != "success":
                failures.append(donor)
                break
            success_count += 1
    summary = {"stage": 3, "n_eligible_donors": len(donors), "n_success_before_failure": success_count, "failures": failures, "finished_utc": utc_now()}
    write_json(CACHE / "stage3_summary.json", summary)
    raise SystemExit(0 if not failures else 2)


if __name__ == "__main__":
    main()
