"""Run the frozen 29-test screen for the frozen Step 10-R disease panel.

The statistical estimator and exposure reader are reused from the completed
Step 10 implementation. Only the outcome-frame reader is extended so that
the frozen panel can use MCQ and BPQ questionnaire files. Results are written
to the separate Step 10-R directory and never overwrite the original Step 10
outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
STEP10_DIR = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step10_cross_disease_replication"
OUT_DIR = Path(__file__).resolve().parent
DEFAULT_TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
DEFAULT_REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_PANEL = OUT_DIR / "step10r_randomized_disease_panel.csv"
FDR_DENOMINATOR = 29


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def build_outcome_frame(model, data_dir: Path, variable: str, cycles: list[str]) -> tuple[pd.DataFrame, list[dict]]:
    """Read a fixed questionnaire outcome without loading any exposure data."""
    component = "BPQ" if variable.upper().startswith("BPQ") else "MCQ"
    specs = {str(spec["cycle"]): (idx, spec) for idx, spec in enumerate(model.CYCLES)}
    frames: list[pd.DataFrame] = []
    audit: list[dict] = []
    for cycle in cycles:
        if cycle not in specs:
            raise ValueError(f"Unknown cycle: {cycle}")
        idx, _ = specs[cycle]
        source_path = data_dir / f"{cycle}_{component}.XPT"
        if not source_path.exists():
            raise FileNotFoundError(source_path)
        demo_path = data_dir / f"{cycle}_DEMO.XPT"
        bmx_path = data_dir / f"{cycle}_BMX.XPT"
        smq_path = data_dir / f"{cycle}_SMQ.XPT"
        source = model.read_xpt(source_path)
        if variable not in source.columns:
            raise ValueError(f"Outcome variable {variable} missing in {source_path.name}")
        core = model.derive_demo(model.read_xpt(demo_path), idx)
        core = core.merge(model.derive_bmx(model.read_xpt(bmx_path)), on="SEQN", how="left", validate="one_to_one")
        core = core.merge(model.derive_smoking(model.read_xpt(smq_path)), on="SEQN", how="left", validate="one_to_one")
        frame = source[["SEQN", variable]].merge(core, on="SEQN", how="inner", validate="one_to_one")
        value = numeric(frame[variable])
        frame["outcome"] = value.where(value.isin([1, 2])).map({1: 1.0, 2: 0.0})
        frame["cycle"] = cycle
        frame["cycle_index"] = idx
        frame["adult"] = numeric(frame["age"]).ge(20)
        frames.append(frame)
        adult = frame["adult"]
        valid = adult & frame["outcome"].notna()
        audit.append({
            "cycle": cycle,
            "outcome_component": component,
            "outcome_variable": variable,
            "source_file": source_path.name,
            "outcome_source_rows": int(len(source)),
            "merged_core_rows": int(len(frame)),
            "adult_rows": int(adult.sum()),
            "valid_binary_rows": int(valid.sum()),
            "case_rows": int((valid & frame["outcome"].eq(1)).sum()),
            "control_rows": int((valid & frame["outcome"].eq(0)).sum()),
        })
    out = pd.concat(frames, ignore_index=True)
    return out.loc[out["adult"]].copy(), audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", type=Path, default=DEFAULT_TESTS)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--outdir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tests = pd.read_csv(args.tests, dtype=str, keep_default_na=False)
    registry = pd.read_csv(args.registry, low_memory=False)
    panel = pd.read_csv(args.panel, dtype=str, keep_default_na=False).sort_values("randomization_order")
    if len(tests) != FDR_DENOMINATOR or tests["test_id"].nunique() != FDR_DENOMINATOR:
        raise ValueError("Step 10-R requires exactly the frozen 29-test family")
    if panel.empty:
        raise ValueError("Step 10-R frozen panel is empty")

    # The bundled workspace runtime supplies numpy/pandas, while the shared
    # survey estimator also requires scipy. Keep this analysis dependency
    # local and untracked.
    runtime = ROOT / "work" / "step08a_pathway_robustness" / "runtime"
    if runtime.exists():
        sys.path.insert(0, str(runtime))
    screen = load_module(STEP10_DIR / "run_step10d_cross_disease_screen.py", "step10r_screen")
    model = screen.load_statistical_model(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "step10r_survey_model")
    reader = screen.load_module(FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py", "step10r_exposure_reader")
    model.DATA_DIR = args.data_dir
    screen.build_outcome_frame = build_outcome_frame

    summaries = []
    for _, panel_row in panel.iterrows():
        order = int(panel_row["randomization_order"])
        _, _, summary = screen.run_one_disease(panel_row, tests, registry, model, reader, args.data_dir, args.outdir, order)
        summaries.append(summary)

    summary_df = pd.DataFrame(summaries).sort_values("disease_id").reset_index(drop=True)
    summary_df.to_csv(args.outdir / "step10r_cross_disease_replication_summary.csv", index=False)
    generated = datetime.now(timezone.utc).isoformat()
    eligible_pool_n = len(pd.read_csv(args.panel.parent / "step10r_eligible_disease_pool.csv", dtype=str, keep_default_na=False))
    if len(panel) == eligible_pool_n:
        panel_description = "the full eligible outcome pool (all outcomes in the pre-randomization registry)"
        analysis_mode = "additive full-pool coverage audit"
    else:
        panel_description = "a randomized five-outcome panel"
        analysis_mode = "randomized primary-panel audit"
    lines = [
        "# Step 10-R expanded cross-disease replication report",
        "",
        f"Generated (UTC): {generated}",
        "",
        f"The outcome pool was expanded using a fixed questionnaire/module inventory before exposure values and association results were loaded. The frozen 29-test environmental family was then screened separately within {panel_description}.",
        "",
        f"- Analysis mode: **{analysis_mode}**.",
        f"- Expanded eligible pool: **{eligible_pool_n}** outcomes (from the pre-randomization registry).",
        f"- Selected panel: **{len(panel)}** outcomes.",
        f"- Frozen exposure tests per outcome: **{FDR_DENOMINATOR}**.",
        f"- BH-FDR denominator per outcome: **{FDR_DENOMINATOR}**.",
        "- No outcome was replaced using its association results.",
        "",
        "## Replication summary",
        "",
        "| Disease | Component | Variable | Pooled cases | Estimable tests | Nominal P<0.05 | FDR<0.05 | Branch | Technical warnings |",
        "|---|---|---|---:|---:|---:|---:|---|---:|",
    ]
    for _, row in summary_df.iterrows():
        p = panel.loc[panel["disease_id"].eq(row["disease_id"])].iloc[0]
        lines.append(f"| {row['disease_name']} | {p['source_component']} | {p['source_variable(s)']} | {row['pooled_cases']} | {row['estimable_tests']}/{FDR_DENOMINATOR} | {row['nominal_positive_tests']} | {row['fdr_positive_tests']} | {row['branch']} | {row['technical_warning_count']} |")
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "This expanded panel is a replication/transportability audit of the outcome-firewalled test family. It does not establish causal effects, temporal ordering, or independence among disease outcomes.",
    ]
    report_path = args.outdir / "STEP10R_CROSS_DISEASE_REPLICATION_REPORT.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    manifest = {
        "lock_type": "STEP10R_CROSS_DISEASE_SCREEN",
        "generated_utc": generated,
        "tests_path": str(args.tests),
        "tests_sha256": sha256(args.tests),
        "registry_path": str(args.registry),
        "registry_sha256": sha256(args.registry),
        "panel_path": str(args.panel),
        "panel_sha256": sha256(args.panel),
        "randomization_lock_sha256": sha256(args.outdir / "STEP10R_RANDOMIZATION_LOCK.json"),
        "expanded_outcome_pool_n": int(len(pd.read_csv(args.panel.parent / "step10r_eligible_disease_pool.csv", dtype=str, keep_default_na=False))),
        "frozen_panel_n": int(len(panel)),
        "frozen_test_count": FDR_DENOMINATOR,
        "fdr_denominator_per_disease": FDR_DENOMINATOR,
        "model_specification": "disease ~ log2(exposure) + age + sex + race + BMI + smoking + PIR; urine tests additionally + log2(creatinine)",
        "association_results_generated_after_randomization_lock": True,
        "outcome_source_components": sorted(panel["source_component"].unique().tolist()),
        "output_hashes": {name: sha256(args.outdir / name) for name in ["step10r_cross_disease_replication_summary.csv", "STEP10R_CROSS_DISEASE_REPLICATION_REPORT.md"]},
    }
    (args.outdir / "STEP10R_SCREEN_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
