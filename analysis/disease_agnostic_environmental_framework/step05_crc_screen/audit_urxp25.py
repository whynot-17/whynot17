"""Diagnostic and isolated reanalysis for the URXP25 Step 5 failure.

This script is deliberately scoped to URXP25.  It does not edit the frozen
29-test screen or recompute BH-FDR.  It contrasts the existing Step 5
phthalate-shaped participant frame with a cycle-matched NHANES outcome and
covariate frame, then fits URXP25 only on the corrected frame.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
DEFAULT_TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
DEFAULT_REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_OUT = FRAMEWORK / "step05_crc_screen"


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


def read_universal_cycle(model, spec: dict, cycle_index: int, data_dir: Path) -> pd.DataFrame:
    """Build an outcome/covariate frame without requiring phthalate lab data."""
    model.DATA_DIR = data_dir
    cycle = spec["cycle"]
    paths = {
        "mcq": data_dir / f"{cycle}_MCQ.XPT",
        "demo": data_dir / f"{cycle}_DEMO.XPT",
        "bmx": data_dir / f"{cycle}_BMX.XPT",
        "smq": data_dir / f"{cycle}_SMQ.XPT",
        "alq": data_dir / f"{cycle}_ALQ.XPT",
        "diq": data_dir / f"{cycle}_DIQ.XPT",
        "paq": data_dir / f"{cycle}_PAQ.XPT",
        "creatinine": data_dir / f"{cycle}_ALB_CR.XPT",
    }
    frames = {key: model.read_xpt(path) for key, path in paths.items()}
    merged = model.cancer_flags(frames["mcq"], frames["demo"])
    merged = merged.merge(model.derive_demo(frames["demo"], cycle_index), on="SEQN", how="inner", validate="one_to_one")
    for supplement in (
        model.derive_bmx(frames["bmx"]),
        model.derive_smoking(frames["smq"]),
        model.derive_alcohol(frames["alq"]),
        model.derive_diabetes(frames["diq"]),
        model.derive_activity(frames["paq"]),
        frames["creatinine"][["SEQN", "URXUCR"]].drop_duplicates("SEQN"),
    ):
        merged = merged.merge(supplement, on="SEQN", how="left", validate="one_to_one")
    merged["cycle"] = cycle
    merged["cycle_index"] = cycle_index
    merged["creatinine_log2"] = np.log2(numeric(merged["URXUCR"]).where(numeric(merged["URXUCR"]) > 0))
    merged["years_since_crc"] = merged["age"] - merged["crc_diagnosis_age"]
    return merged


def population_complete(population: pd.DataFrame) -> pd.DataFrame:
    required = [
        "outcome", "axis_log2", "age", "bmi", "pir", "sex", "race",
        "smoking", "pooled_weight", "psu", "strata", "creatinine_log2",
    ]
    complete = population.dropna(subset=required).copy()
    return complete[complete["pooled_weight"].gt(0)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", type=Path, default=DEFAULT_TESTS)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    step5 = load_module(args.outdir / "run_step05_crc_screen.py", "urxp25_step5")
    model = load_module(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "urxp25_model")
    tests = pd.read_csv(args.tests, dtype=str, keep_default_na=False)
    registry = pd.read_csv(args.registry, low_memory=False)
    test_row = tests.loc[tests["variable"].eq("URXP25")].iloc[0]
    exposure, source = step5.read_test_exposure(test_row, registry)

    # The pre-repair Step 5 frame is retained only as an audit comparator.
    old_harmonized, _ = step5.load_harmonized(model, args.data_dir)
    old_merge = exposure.merge(old_harmonized, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
    old_population = model.population_frames(old_merge)["CRC_vs_cancer_free"]
    old_complete = old_population.iloc[0:0]

    specs = {spec["cycle"]: (index, spec) for index, spec in enumerate(model.CYCLES)}
    universal_frames = []
    for cycle in source["cycles"]:
        index, spec = specs[cycle]
        universal_frames.append(read_universal_cycle(model, spec, index, args.data_dir))
    universal = pd.concat(universal_frames, ignore_index=True)
    repaired_merge = exposure.merge(universal, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
    repaired_population = model.population_frames(repaired_merge)["CRC_vs_cancer_free"]
    repaired_population["pooled_weight"] = repaired_population["test_weight"] / len(source["cycles"])
    repaired_complete = population_complete(repaired_population)

    fit = model.fit_survey_logistic(
        repaired_population,
        ["axis_log2", "age", "bmi", "pir", "creatinine_log2"],
        ["sex", "race", "smoking"],
        exposure_name="axis_log2",
        levels=model.LEVELS,
    )
    old_step5 = pd.read_csv(args.outdir / "full_29_test_crc_screen.csv", low_memory=False).loc[lambda x: x["variable"].eq("URXP25")].iloc[0]

    cycle_rows = []
    for cycle in source["cycles"]:
        e = exposure.loc[exposure["cycle"].eq(cycle)]
        u = universal.loc[universal["cycle"].eq(cycle)]
        old_u = old_harmonized.loc[old_harmonized["cycle"].eq(cycle)]
        z = e.merge(u, on=["SEQN", "cycle"], how="inner", validate="one_to_one")
        z_old = e.merge(old_u, on=["SEQN", "cycle"], how="inner", validate="one_to_one")
        cycle_rows.append(
            {
                "cycle": cycle,
                "source_exposure_rows": int(len(e)),
                "source_nonmissing_exposure": int(e["exposure_raw"].notna().sum()),
                "universal_outcome_covariate_rows": int(len(u)),
                "old_step5_phthalate_frame_rows": int(len(old_u)),
                "repaired_merge_rows": int(len(z)),
                "old_step5_merge_rows": int(len(z_old)),
                "repaired_crc_vs_cancer_free_rows": int(len(model.population_frames(z)["CRC_vs_cancer_free"])),
                "old_crc_vs_cancer_free_rows": int(len(model.population_frames(z_old)["CRC_vs_cancer_free"])),
                "repaired_complete_case_rows": int(len(population_complete(model.population_frames(z)["CRC_vs_cancer_free"].assign(pooled_weight=z["test_weight"] / len(source["cycles"])))) if len(z) else 0),
            }
        )
    cycle_df = pd.DataFrame(cycle_rows)

    result = {
        "variable": "URXP25",
        "chemical_mapping": str(test_row["chemical_names"]),
        "source_cycles": ";".join(source["cycles"]),
        "source_registry_status": source.get("status"),
        "source_registry_rows": int(source.get("n_raw", 0) or 0),
        "source_exposure_rows": int(len(exposure)),
        "source_nonmissing_exposure": int(exposure["exposure_raw"].notna().sum()),
        "old_step5_harmonized_rows": int(len(old_harmonized)),
        "old_step5_merge_rows": int(len(old_merge)),
        "old_step5_analytic_n": int(old_step5["N"]),
        "old_step5_status": str(old_step5["status"]),
        "repaired_universal_frame_rows": int(len(universal)),
        "repaired_merge_rows": int(len(repaired_merge)),
        "repaired_complete_case_n": int(len(repaired_complete)),
        "repaired_crc_cases": int(repaired_complete["outcome"].sum()) if len(repaired_complete) else 0,
        "repaired_controls": int(len(repaired_complete) - repaired_complete["outcome"].sum()) if len(repaired_complete) else 0,
        "repaired_status": str(fit.get("status")),
        "repaired_beta": fit.get("beta", np.nan),
        "repaired_SE": fit.get("SE", np.nan),
        "repaired_OR": fit.get("OR", np.nan),
        "repaired_CI_low": fit.get("CI_low", np.nan),
        "repaired_CI_high": fit.get("CI_high", np.nan),
        "repaired_P": fit.get("P", np.nan),
        "repaired_design_df": fit.get("design_df", np.nan),
        "repaired_message": fit.get("message", ""),
        "repaired_model": "CRC ~ log2(URXP25) + age + sex + race + BMI + smoking + PIR + log2(creatinine)",
        "repair_scope": "isolated URXP25 diagnostic/reanalysis; frozen 29-test screen and BH-FDR unchanged",
    }
    result_df = pd.DataFrame([result])
    result_df.to_csv(args.outdir / "URXP25_REANALYSIS.csv", index=False)
    cycle_df.to_csv(args.outdir / "URXP25_DIAGNOSTIC_BY_CYCLE.csv", index=False)

    report = [
        "# URXP25 diagnostic audit",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Finding",
        "",
        "URXP25 did not fail because it had too few CRC cases. The Step 5 outcome/covariate frame was rebuilt from the MBzP/phthalate laboratory file (`PHTHTE_*`), while URXP25 is measured in PAH files (`PAH_H/I/J`). The participant IDs in the three URXP25 source files had zero intersection with that phthalate-shaped frame.",
        "",
        f"- Registry/exposure source rows: **{len(exposure):,}**; non-missing positive exposure values: **{int(exposure['exposure_raw'].notna().sum()):,}**.",
        f"- Pre-repair Step 5 merge rows: **{len(old_merge):,}**; analytic N: **{int(old_step5['N'])}**; status: **{old_step5['status']}**.",
        f"- Corrected cycle-matched outcome/covariate merge rows: **{len(repaired_merge):,}**; complete-case N: **{len(repaired_complete):,}**; CRC cases: **{int(repaired_complete['outcome'].sum()) if len(repaired_complete) else 0}**.",
        "",
        "## Isolated corrected reanalysis",
        "",
        f"Using the same Step 5 estimator and covariate model on the corrected PAH-compatible frame: **OR={fit.get('OR', np.nan):.6g}**, 95% CI **{fit.get('CI_low', np.nan):.6g}–{fit.get('CI_high', np.nan):.6g}**, P **{fit.get('P', np.nan):.6g}**, status **{fit.get('status', '')}**.",
        "",
        "This isolated result is not inserted into `full_29_test_crc_screen.csv`, is not assigned a BH-FDR, and does not alter any of the other 28 tests.",
        "",
        "## Important scope note",
        "",
        "The zero-intersection diagnosis exposes a broader Step 5 architecture issue: the existing harmonized frame is phthalate-subsample shaped, so non-phthalate tests may not have been evaluated in their own laboratory subsamples. This audit intentionally does not rerun those tests; before the 29-test screen is treated as final, each non-phthalate test should be checked with its own cycle-matched outcome/covariate frame.",
    ]
    (args.outdir / "URXP25_DIAGNOSTIC_AUDIT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = {
        "lock_type": "URXP25_DIAGNOSTIC_REANALYSIS",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "URXP25 only",
        "does_not_modify": ["full_29_test_crc_screen.csv", "crc_bh_fdr_29_tests.csv", "other 28 tests"],
        "tests_input": str(args.tests),
        "tests_input_sha256": sha256(args.tests),
        "registry_input": str(args.registry),
        "registry_input_sha256": sha256(args.registry),
        "model_script": str(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"),
        "model_script_sha256": sha256(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"),
        "reanalysis": result,
        "outputs": {
            "URXP25_REANALYSIS.csv": sha256(args.outdir / "URXP25_REANALYSIS.csv"),
            "URXP25_DIAGNOSTIC_BY_CYCLE.csv": sha256(args.outdir / "URXP25_DIAGNOSTIC_BY_CYCLE.csv"),
            "URXP25_DIAGNOSTIC_AUDIT.md": sha256(args.outdir / "URXP25_DIAGNOSTIC_AUDIT.md"),
        },
    }
    (args.outdir / "URXP25_REANALYSIS_LOCK.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
