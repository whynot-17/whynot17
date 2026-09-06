"""Freeze the Step 10-R disease panel before loading exposure results."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


OUT_DIR = Path(__file__).resolve().parent
RANDOM_SEED = 20260828
PANEL_SIZE = 5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pool", type=Path, default=OUT_DIR / "step10r_eligible_disease_pool.csv")
    parser.add_argument("--outdir", type=Path, default=OUT_DIR)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--panel-size", type=int, default=PANEL_SIZE)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    pool = pd.read_csv(args.pool, dtype=str, keep_default_na=False)
    eligible = pool.loc[pool["eligible_for_randomization"].eq("True")].sort_values("disease_id").reset_index(drop=True)
    if len(eligible) < args.panel_size:
        raise ValueError(f"Expanded pool has {len(eligible)} eligible outcomes, fewer than requested {args.panel_size}")
    rng = random.Random(args.seed)
    selected_ids = rng.sample(eligible["disease_id"].tolist(), k=args.panel_size)
    selected = eligible.set_index("disease_id").loc[selected_ids].reset_index()
    selected.insert(0, "randomization_order", range(1, len(selected) + 1))
    selected["random_seed"] = args.seed
    selection_status = "frozen_step10r_full_pool_coverage" if args.panel_size == len(eligible) else "frozen_step10r_primary_panel"
    selected["selection_status"] = selection_status
    selected.to_csv(args.outdir / "step10r_randomized_disease_panel.csv", index=False)

    generated = datetime.now(timezone.utc).isoformat()
    lock = {
        "lock_type": "STEP10R_RANDOMIZATION_LOCK",
        "generated_utc": generated,
        "random_seed": args.seed,
        "requested_panel_size": args.panel_size,
        "frozen_panel_size": len(selected),
        "eligible_pool_n": len(eligible),
        "eligible_pool_sha256": sha256(args.pool),
        "selection_algorithm": "random.Random(seed).sample(sorted(eligible disease_id), k)",
        "analysis_mode": "full_eligible_pool_coverage" if args.panel_size == len(eligible) else "randomized_primary_panel",
        "association_results_loaded_before_randomization": False,
        "exposure_values_loaded_before_randomization": False,
        "replacement_rule": "Only R1 outcome non-reconstructible, R2 <50% technically estimable tests, R3 no compatible survey design, or R4 invalid outcome definition; no result-driven replacement.",
        "selected_disease_ids": selected["disease_id"].tolist(),
        "output_sha256": sha256(args.outdir / "step10r_randomized_disease_panel.csv"),
    }
    (args.outdir / "STEP10R_RANDOMIZATION_LOCK.json").write_text(json.dumps(lock, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report = [
        "# Step 10-R randomization lock",
        "",
        f"Generated (UTC): {generated}",
        "",
        f"- Expanded eligible pool: **{len(eligible)}** outcomes.",
        f"- Requested panel: **{args.panel_size}**; frozen panel: **{len(selected)}**.",
        f"- Random seed: **{args.seed}**.",
        "- Randomization occurred before exposure values and association results were loaded.",
        "",
        "## Frozen panel",
        "",
        "| Order | Disease ID | Disease | Component | Variable | Cases | Cycles |",
        "|---:|---|---|---|---|---:|---:|",
    ]
    for _, row in selected.iterrows():
        report.append(f"| {row['randomization_order']} | {row['disease_id']} | {row['disease_name']} | {row['source_component']} | {row['source_variable(s)']} | {row['case_count_pooled']} | {row['n_cycles']} |")
    (args.outdir / "STEP10R_RANDOMIZATION_LOCK.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
