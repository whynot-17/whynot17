"""Step 6: frozen, uniform robustness audit for the two Step 5 FDR signals.

The Step 5 screen is the source of truth for the primary estimates and the
29-test BH-FDR family.  This script only reruns and stress-tests the two
signals supported by that screen (URXCOP/MCOP and LBXPFHS/PFHS); it never
recomputes or narrows the Step 5 multiplicity family.
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
STEP5 = FRAMEWORK / "step05_crc_screen"
DEFAULT_TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
DEFAULT_REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_OUT = FRAMEWORK / "step06_robustness"
RUBRIC = DEFAULT_OUT / "ROBUSTNESS_RUBRIC_LOCK.md"
FDR_DENOMINATOR = 29
TARGET_VARIABLES = ["URXCOP", "LBXPFHS"]
ATTENUATION_LOG_THRESHOLD = 0.25
EXTREME_WEIGHT_RATIO_THRESHOLD = 100.0


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


def finite(value: object) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def status_of(item: object) -> str | None:
    if isinstance(item, dict):
        value = item.get("fit_status", item.get("status"))
        return str(value) if value is not None else None
    return None


def dataframe_statuses(frame: pd.DataFrame) -> list[str]:
    if frame is None or frame.empty:
        return []
    column = "fit_status" if "fit_status" in frame.columns else "status" if "status" in frame.columns else None
    if column is None:
        return []
    return [str(x) for x in frame[column].tolist() if pd.notna(x)]


def axis_from_test(test_row: pd.Series) -> pd.Series:
    return pd.Series(
        {
            "axis_key": f"{test_row['matrix']}|{test_row['variable']}",
            "exposure_axis": str(test_row.get("exposure_axes", test_row.get("biomarker", test_row["variable"]))),
            "primary_biomarker": str(test_row["variable"]),
            "biological_matrix": str(test_row["matrix"]),
            "eligible_chemical_count": int(float(test_row.get("mapping_count", 0) or 0)),
            "eligible_chemical_ids": str(test_row.get("chemical_ids", "")),
            "eligible_chemical_names": str(test_row.get("chemical_names", "")),
            "mapping_confidence": "frozen Step 4 mapping",
        }
    )


def build_population(test_row: pd.Series, registry: pd.DataFrame, harmonized: pd.DataFrame, model, step5):
    exposure, source = step5.read_test_exposure(test_row, registry)
    if exposure.empty:
        return pd.DataFrame(), source
    merged = exposure.merge(
        harmonized,
        on=["SEQN", "cycle"],
        how="inner",
        validate="one_to_one",
        suffixes=("", "_outcome"),
    )
    return model.population_frames(merged)["CRC_vs_cancer_free"], source


def extra_fit(frame: pd.DataFrame, model, robust, urine: bool, label: str, categorical: list[str] | None = None, continuous_extra: list[str] | None = None, exposure_name: str = "axis_log2", include_creatinine: bool = True) -> tuple[dict, dict]:
    categorical = list(categorical or ["sex", "race", "smoking"])
    continuous_extra = list(continuous_extra or [])
    continuous = [exposure_name, "age", "bmi", "pir", *continuous_extra]
    if urine and include_creatinine:
        continuous.append("creatinine_log2")
    levels = model.LEVELS
    fit = model.fit_survey_logistic(frame, continuous, categorical, exposure_name=exposure_name, levels=levels)
    diagnostics = robust.diagnostic_context(frame, model, exposure_name, urine, categorical, [*continuous_extra, *( ["creatinine_log2"] if urine and include_creatinine else [])], levels)
    row = {
        "analysis": label,
        "exposure_variable": exposure_name,
        **robust.fit_clean(fit),
        **diagnostics,
        "fit_status": fit.get("status", "not_estimable"),
        "warning_message": fit.get("message", "") if fit.get("status") != "ok" else "",
        "warning_type": "none" if fit.get("status") == "ok" else ("convergence_warning" if fit.get("status") == "converged_with_warning" else "fit_failure_or_non_estimable"),
        "analytic_n": fit.get("N", np.nan),
        "crc_cases": fit.get("CRC_N", np.nan),
        "controls": fit.get("Control_N", np.nan),
    }
    return row, fit


def run_age40(axis: pd.Series, population: pd.DataFrame, model, robust, urine: bool) -> dict:
    frame = population.loc[population["age"].ge(40)].copy()
    row, _ = robust.fit_axis_model(frame, model, "axis_log2", urine, "age_ge_40")
    return robust.add_meta(row, axis)


def run_sex_specific(axis: pd.Series, population: pd.DataFrame, model, robust, urine: bool) -> list[dict]:
    rows = []
    for sex in ["Female", "Male"]:
        frame = population.loc[population["sex"].eq(sex)].copy()
        row, _ = extra_fit(frame, model, robust, urine, f"sex_specific_{sex}", categorical=["race", "smoking"])
        row.update({"sex_group": sex})
        rows.append(robust.add_meta(row, axis))
    return rows


def run_sex_interaction(axis: pd.Series, population: pd.DataFrame, model, robust, stability, urine: bool) -> tuple[dict, dict]:
    work = population.copy()
    work["axis_x_sex_Male"] = work["axis_log2"] * work["sex"].eq("Male").astype(float)
    row, fit = robust.fit_axis_model(
        work,
        model,
        "axis_log2",
        urine,
        "sex_interaction",
        continuous_extra=["axis_x_sex_Male"],
    )
    test = stability.wald_test(fit, ["axis_x_sex_Male"])
    row.update(
        {
            "interaction_term": "axis_x_sex_Male",
            "interaction_P": test.get("P_F", np.nan),
            "interaction_P_chi2": test.get("P_chi2", np.nan),
            "interaction_test": "Wald F test",
            "interaction_df_num": test.get("df_num", np.nan),
            "interaction_df_denom": test.get("df_denom", np.nan),
        }
    )
    return robust.add_meta(row, axis), fit


def run_pairwise_coexposure(axis: pd.Series, population: pd.DataFrame, coexposure: pd.DataFrame, model, robust, urine: bool, coex_variable: str) -> dict:
    coex = coexposure[["SEQN", "cycle", "axis_log2"]].rename(columns={"axis_log2": "coexposure_log2"}).drop_duplicates(["SEQN", "cycle"])
    work = population.merge(coex, on=["SEQN", "cycle"], how="inner", validate="one_to_one")
    row, _ = robust.fit_axis_model(work, model, "axis_log2", urine, "pairwise_coexposure", continuous_extra=["coexposure_log2"])
    row.update({"coexposure_variable": coex_variable, "coexposure_complete_n": int(len(work))})
    return robust.add_meta(row, axis)


def collect_warning_statuses(primary: dict, loco: pd.DataFrame, cycle: pd.DataFrame, hetero: dict, timing: pd.DataFrame, tail: pd.DataFrame, lod: dict, creat: dict, age40: dict, sex_rows: list[dict], interaction: dict, coexposure: dict) -> list[str]:
    statuses = []
    for item in [primary, hetero, lod, creat, age40, interaction, coexposure, *sex_rows]:
        value = status_of(item)
        if value and value != "not_applicable":
            statuses.append(value)
    for frame in [loco, cycle, timing, tail]:
        statuses.extend(x for x in dataframe_statuses(frame) if x != "not_applicable")
    return statuses


def build_warning(axis: pd.Series, primary_row: dict, primary_reproduction_diff: float, statuses: list[str], robust) -> dict:
    warning_count = sum(status == "converged_with_warning" for status in statuses)
    failure_count = sum(status not in {"ok", "converged_with_warning"} for status in statuses)
    fraction = float(warning_count / len(statuses)) if statuses else 0.0
    persistent = bool(warning_count and (fraction >= 0.75 or (primary_row.get("fit_status") == "converged_with_warning" and fraction >= 0.50)))
    if failure_count:
        reason = "At least one applicable audit fit was non-estimable or failed; inspect module-level fit_status."
        impact = "Not fully assessed because at least one applicable audit component was non-estimable or failed."
    elif warning_count:
        reason = "Finite estimable fits were returned, but Newton-IRLS did not reach the configured tolerance within the iteration limit for some audit fits."
        impact = "The Step 6 primary rerun was compared with the frozen Step 5 estimate; no material discrepancy is permitted."
    else:
        reason = "No convergence or estimability warning detected across applicable fits."
        impact = "Not applicable; all applicable fits were reported as converged."
    distribution = ";".join(f"{status}={statuses.count(status)}" for status in sorted(set(statuses)))
    row = {
        "fit_warning_count": warning_count,
        "fit_failure_count": failure_count,
        "warning_count": warning_count + failure_count,
        "n_fit_evaluations": len(statuses),
        "warning_fraction": fraction,
        "warning_persistent": persistent,
        "primary_fit_status": primary_row.get("fit_status"),
        "primary_warning_type": primary_row.get("warning_type"),
        "primary_warning_message": primary_row.get("warning_message", ""),
        "warning_status_distribution": distribution,
        "warning_technical_reason": reason,
        "coefficient_impact_assessment": impact,
        "warning_repeated_in_sensitivity": bool(warning_count > 0),
        "primary_reproduction_abs_logOR_difference": primary_reproduction_diff,
        "warning_assessment": "no meaningful warning" if not warning_count and not failure_count else "persistent convergence warning" if persistent else "localized convergence warning" if not failure_count else "persistent technical/non-estimable issue",
        "audit_status": "pass" if not failure_count and not persistent else "warning" if not failure_count else "technical_concern",
    }
    return robust.add_meta(row, axis)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", type=Path, default=DEFAULT_TESTS)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if not RUBRIC.exists():
        raise FileNotFoundError(f"Frozen rubric is required before running Step 6: {RUBRIC}")

    step5 = load_module(STEP5 / "run_step05_crc_screen.py", "step06_step5")
    robust = load_module(ROOT / "work" / "scripts" / "environmental_crc_15axis_robustness_audit.py", "step06_robustness_core")
    model = load_module(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "step06_survey_model")
    stability = load_module(ROOT / "work" / "scripts" / "mcop_crc_phase2_stability.py", "step06_stability")

    tests = pd.read_csv(args.tests, dtype=str, keep_default_na=False)
    registry = pd.read_csv(args.registry, low_memory=False)
    full_screen = pd.read_csv(STEP5 / "full_29_test_crc_screen.csv", low_memory=False)
    fdr_screen = pd.read_csv(STEP5 / "crc_bh_fdr_29_tests.csv", low_memory=False)
    if len(tests) != FDR_DENOMINATOR or tests["test_id"].nunique() != FDR_DENOMINATOR:
        raise ValueError("Step 6 requires the immutable 29-test Step 4 test set")
    if len(full_screen) != FDR_DENOMINATOR or len(fdr_screen) != FDR_DENOMINATOR:
        raise ValueError("Step 6 requires the immutable 29-test Step 5 screen outputs")
    supported = full_screen.loc[full_screen["variable"].isin(TARGET_VARIABLES) & full_screen["FDR_supported"].astype(str).str.lower().eq("true")].copy()
    if supported["variable"].tolist() != TARGET_VARIABLES:
        supported = supported.sort_values("variable")
    if len(supported) != 2 or set(supported["variable"]) != set(TARGET_VARIABLES):
        raise AssertionError("Step 5 supported-signal scope is not exactly URXCOP + LBXPFHS")

    model.DATA_DIR = args.data_dir
    harmonized, source_manifest = step5.load_harmonized(model, args.data_dir)
    populations: dict[str, pd.DataFrame] = {}
    exposures: dict[str, pd.DataFrame] = {}
    sources: dict[str, dict] = {}
    for _, test_row in tests.loc[tests["variable"].isin(TARGET_VARIABLES)].iterrows():
        exposure, source = step5.read_test_exposure(test_row, registry)
        exposures[str(test_row["variable"])] = exposure
        sources[str(test_row["variable"])] = source
        populations[str(test_row["variable"])], _ = build_population(test_row, registry, harmonized, model, step5)

    all_summary, all_loco, all_cycle, all_hetero, all_sensitivity, all_warnings = [], [], [], [], [], []
    primary_qc, source_qc = [], []
    score_rows, fingerprint_rows = [], []

    for _, screen_row in supported.sort_values("variable").iterrows():
        variable = str(screen_row["variable"])
        test_row = tests.loc[tests["variable"].eq(variable)].iloc[0]
        axis = axis_from_test(test_row)
        population = populations[variable]
        urine = str(test_row["matrix"]).lower() == "urine"
        if population.empty:
            raise RuntimeError(f"Empty CRC population for {variable}")

        primary_row, primary_raw = robust.primary_fit(axis, population, model, urine)
        stored_or = float(screen_row["OR"])
        rerun_or = float(primary_row["OR"])
        reproduction_diff = abs(float(np.log(rerun_or / stored_or)))
        primary_qc.append({
            "variable": variable,
            "step5_OR": stored_or,
            "step6_rerun_OR": rerun_or,
            "absolute_logOR_difference": reproduction_diff,
            "step5_P": float(screen_row["P"]),
            "step5_BH_FDR_29": float(screen_row["BH_FDR"]),
            "step5_status": str(screen_row["status"]),
            "step6_status": str(primary_row.get("fit_status")),
        })
        if not finite(reproduction_diff) or reproduction_diff > 1e-8:
            raise AssertionError(f"Step 6 primary rerun does not reproduce Step 5 for {variable}: {reproduction_diff}")

        pooled_or = float(primary_row["OR"])
        loco, loco_summary = robust.run_loco(axis, population, model, urine, pooled_or)
        cycle, cycle_summary = robust.run_cycle_specific(axis, population, model, urine, pooled_or)
        hetero, hetero_summary = robust.run_heterogeneity(axis, population, model, urine)
        timing = robust.run_timing(axis, population, model, urine, pooled_or)
        tail = robust.run_tail(axis, population, model, urine, pooled_or)
        axis_detect = registry.loc[registry["variable"].eq(variable)].copy()
        lod = robust.run_lod(axis, population, model, urine, pooled_or, axis_detect)
        creat = robust.run_creatinine(axis, population, model, urine, pooled_or)
        age40 = run_age40(axis, population, model, robust, urine)
        sex_rows = run_sex_specific(axis, population, model, robust, urine)
        sex_interaction, sex_interaction_raw = run_sex_interaction(axis, population, model, robust, stability, urine)
        other = "LBXPFHS" if variable == "URXCOP" else "URXCOP"
        coexposure = run_pairwise_coexposure(axis, population, exposures[other], model, robust, urine, other)

        statuses = collect_warning_statuses(primary_row, loco, cycle, hetero, timing, tail, lod, creat, age40, sex_rows, sex_interaction, coexposure)
        warning = build_warning(axis, primary_row, reproduction_diff, statuses, robust)
        score = robust.scorecard(axis, primary_raw, screen_row, loco_summary, cycle_summary, hetero_summary, timing, tail, creat, warning)
        score["primary_BH_FDR_29tests"] = score.pop("primary_BH_FDR_15axis")
        score["step5_fdr_denominator"] = FDR_DENOMINATOR
        score["sex_interaction_P"] = sex_interaction.get("interaction_P", np.nan)
        score["coexposure_variable"] = other
        score_rows.append(score)
        fingerprint_rows.append({key: score.get(key) for key in ["axis_key", "primary_biomarker", "primary_OR", "primary_P", "primary_BH_FDR_29tests", "F", "L", "C", "H", "D", "T", "A", "E", "robustness_fingerprint", "robustness_tier"]})

        all_summary.append({**score, **{f"primary_{key}": value for key, value in primary_row.items() if key not in {"axis_key", "exposure_axis", "primary_biomarker", "biological_matrix"}}})
        all_loco.append(loco)
        all_cycle.append(cycle)
        all_hetero.append(hetero)
        sensitivity = pd.concat(
            [
                timing.assign(sensitivity_domain="diagnosis_timing"),
                tail.assign(sensitivity_domain="upper_tail"),
                pd.DataFrame([lod]).assign(sensitivity_domain="lod"),
                pd.DataFrame([creat]).assign(sensitivity_domain="creatinine"),
                pd.DataFrame([age40]).assign(sensitivity_domain="age_ge_40"),
                pd.DataFrame(sex_rows).assign(sensitivity_domain="sex_specific"),
                pd.DataFrame([sex_interaction]).assign(sensitivity_domain="sex_interaction"),
                pd.DataFrame([coexposure]).assign(sensitivity_domain="pairwise_coexposure"),
            ],
            ignore_index=True,
            sort=False,
        )
        all_sensitivity.append(sensitivity)
        all_warnings.append(pd.DataFrame([warning]))
        source_qc.append({
            "variable": variable,
            "source_status": sources[variable].get("status"),
            "source_cycles": ";".join(sources[variable].get("cycles", [])),
            "source_n": sources[variable].get("n_raw", 0),
            "primary_population_n": len(population),
            "primary_population_crc_cases": int(population["outcome"].sum()),
            "primary_population_controls": int(len(population) - population["outcome"].sum()),
            "source_rows": len(sources[variable].get("source_rows", [])),
        })

    summary_df = pd.DataFrame(all_summary)
    loco_df = pd.concat(all_loco, ignore_index=True)
    cycle_df = pd.concat(all_cycle, ignore_index=True)
    hetero_df = pd.DataFrame(all_hetero)
    sensitivity_df = pd.concat(all_sensitivity, ignore_index=True)
    warnings_df = pd.concat(all_warnings, ignore_index=True)
    score_df = pd.DataFrame(score_rows)
    fingerprint_df = pd.DataFrame(fingerprint_rows)
    primary_qc_df = pd.DataFrame(primary_qc)
    source_qc_df = pd.DataFrame(source_qc)

    outputs = {
        "robustness_results.csv": summary_df,
        "loco_results.csv": loco_df,
        "cycle_specific_results.csv": cycle_df,
        "cycle_heterogeneity.csv": hetero_df,
        "sensitivity_results.csv": sensitivity_df,
        "robustness_fingerprints.csv": fingerprint_df,
        "model_warnings.csv": warnings_df,
        "primary_reproduction_qc.csv": primary_qc_df,
        "source_population_qc.csv": source_qc_df,
    }
    for name, frame in outputs.items():
        frame.to_csv(args.outdir / name, index=False)

    robust_a = score_df.loc[score_df["robustness_tier"].eq("Robust Tier A"), "primary_biomarker"].astype(str).tolist()
    persistent = warnings_df.loc[warnings_df["warning_persistent"].astype(bool), "primary_biomarker"].astype(str).tolist()
    report_lines = [
        "# Step 6 robustness audit report",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Frozen scope",
        "",
        "This audit was run after locking `ROBUSTNESS_RUBRIC_LOCK.md`. It applies the same modules, thresholds, estimator, and tier logic to exactly the two Step 5 FDR-supported signals: URXCOP (MCOP) and LBXPFHS (PFHS).",
        "",
        f"- Step 5 test family: **{FDR_DENOMINATOR} frozen tests**; denominator unchanged and not recomputed here.",
        f"- Step 6 audited signals: **{len(score_df)}**.",
        f"- Primary rerun reproduction: **{len(primary_qc_df)}/{len(primary_qc_df)}** within absolute log-OR difference <=1e-8.",
        "",
        "## Primary estimates and scorecard",
        "",
        "| Biomarker | OR | 95% CI | P | BH-FDR (29) | Fingerprint | Tier |",
        "|---|---:|---|---:|---:|---|---|",
    ]
    for _, row in score_df.sort_values("primary_BH_FDR_29tests").iterrows():
        report_lines.append(
            f"| {row['primary_biomarker']} | {float(row['primary_OR']):.4g} | {float(row['primary_CI_low']):.4g}–{float(row['primary_CI_high']):.4g} | {float(row['primary_P']):.4g} | {float(row['primary_BH_FDR_29tests']):.4g} | {row['robustness_fingerprint']} | {row['robustness_tier']} |"
        )
    report_lines += [
        "",
        f"Robust Tier A signals under the frozen rubric: **{', '.join(robust_a) if robust_a else 'none'}**.",
        f"Persistent technical-warning signals (A0): **{', '.join(persistent) if persistent else 'none'}**.",
        "",
        "## Audit interpretation",
        "",
        "- `H` is a reported exposure-by-cycle heterogeneity tag, not a hard deletion gate.",
        "- PFHS and MCOP are evaluated in the same audit. A `converged_with_warning` status is retained and interpreted as a technical warning, never silently upgraded to clean convergence.",
        "- The pairwise co-exposure model is secondary and uses the target biomarker's own survey weight; it does not replace the primary model or change the Step 5 FDR result.",
        "- Sex-specific estimates and the formal exposure-by-sex interaction are reported descriptively/secondarily; they are not used to redefine the primary signal.",
        "- MCOP is interpreted as the urinary biomarker for the DINP-related exposure axis. This audit does not establish causality or prove that MCOP itself is a direct CRC mechanism.",
        "",
        "## Frozen component definitions",
        "",
        "The complete definitions and thresholds are in `ROBUSTNESS_RUBRIC_LOCK.md`; in particular, the attenuation threshold is 0.25 absolute log(OR), and Tier A is F2 + L>=1 + C>=1 + D>=1 + T>=1 + A>=1.",
    ]
    (args.outdir / "STEP6_ROBUSTNESS_REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest = {
        "lock_type": "ROBUST_HUMAN_SIGNAL_LOCK",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "scope": TARGET_VARIABLES,
        "fdr_denominator": FDR_DENOMINATOR,
        "step5_screen": str(STEP5 / "full_29_test_crc_screen.csv"),
        "step5_screen_sha256": sha256(STEP5 / "full_29_test_crc_screen.csv"),
        "step5_fdr": str(STEP5 / "crc_bh_fdr_29_tests.csv"),
        "step5_fdr_sha256": sha256(STEP5 / "crc_bh_fdr_29_tests.csv"),
        "step4_test_set": str(args.tests),
        "step4_test_set_sha256": sha256(args.tests),
        "registry": str(args.registry),
        "registry_sha256": sha256(args.registry),
        "model_script": str(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"),
        "model_script_sha256": sha256(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"),
        "rubric": str(RUBRIC),
        "rubric_sha256": sha256(RUBRIC),
        "estimator": "Python NHANES survey-weighted Taylor/PSU-sandwich logistic estimator reused from Step 5",
        "primary_reproduction_max_absolute_logOR_difference": float(primary_qc_df["absolute_logOR_difference"].max()),
        "primary_results": primary_qc_df.to_dict(orient="records"),
        "tier_counts": score_df["robustness_tier"].value_counts().to_dict(),
        "warning_assessment": warnings_df[["primary_biomarker", "audit_status", "warning_assessment", "warning_status_distribution"]].to_dict(orient="records"),
        "outputs": {name: {"path": str(args.outdir / name), "sha256": sha256(args.outdir / name)} for name in [*outputs, "STEP6_ROBUSTNESS_REPORT.md"]},
    }
    (args.outdir / "ROBUST_HUMAN_SIGNAL_LOCK.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
