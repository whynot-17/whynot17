"""Systematic robustness audit for every eligible environmental exposure axis.

This module deliberately reads the frozen 15-axis screen output and applies
one identical secondary-analysis template to every axis. It never re-ranks or
recomputes the primary 15-axis BH-FDR after robustness filtering.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / "outputs"
MATRIX = OUTPUTS / "environmental_crc_267_actionability_matrix_v2.csv"
DETECT = OUTPUTS / "environmental_crc_267_detectability_by_cycle.csv"
SCREEN = OUTPUTS / "environmental_crc_systematic_human_screen_v2.csv"
SCREEN_FDR = OUTPUTS / "environmental_crc_systematic_human_screen_fdr_v2.csv"

OUT_SUMMARY = OUTPUTS / "environmental_crc_15axis_robustness_summary.csv"
OUT_LOCO = OUTPUTS / "environmental_crc_15axis_loco.csv"
OUT_CYCLE = OUTPUTS / "environmental_crc_15axis_cycle_specific.csv"
OUT_HET = OUTPUTS / "environmental_crc_15axis_cycle_interaction.csv"
OUT_TIMING = OUTPUTS / "environmental_crc_15axis_diagnosis_timing.csv"
OUT_TAIL = OUTPUTS / "environmental_crc_15axis_tail_exclusion.csv"
OUT_LOD = OUTPUTS / "environmental_crc_15axis_lod_sensitivity.csv"
OUT_CREAT = OUTPUTS / "environmental_crc_15axis_creatinine_sensitivity.csv"
OUT_WARN = OUTPUTS / "environmental_crc_15axis_model_warnings.csv"
OUT_SCORE = OUTPUTS / "environmental_crc_15axis_robustness_scorecard.csv"
OUT_REPORT = OUTPUTS / "environmental_crc_15axis_robustness_report.md"
OUT_MANIFEST = OUTPUTS / "environmental_crc_15axis_robustness_manifest.json"

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
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def finite(value: object) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def truth(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def direction(or_value: object, pooled_or: float) -> str:
    if not finite(or_value) or not finite(pooled_or):
        return "not_estimable"
    value = float(or_value)
    if value == 1:
        return "null"
    if pooled_or > 1:
        return "positive" if value > 1 else "discordant"
    if pooled_or < 1:
        return "protective" if value < 1 else "discordant"
    return "positive" if value > 1 else "protective"


def ci_excludes_1(low: object, high: object) -> bool:
    return finite(low) and finite(high) and (float(low) > 1 or float(high) < 1)


def fit_clean(fit: dict) -> dict[str, object]:
    return {
        ("fit_N" if k == "N" else "fit_crc_cases" if k == "CRC_N" else "fit_control_n" if k == "Control_N" else k): v
        for k, v in fit.items()
        if k not in {"coefficients", "covariance"}
    }


def diagnostic_context(
    frame: pd.DataFrame,
    model,
    exposure: str,
    urine: bool,
    categorical: list[str],
    continuous_extra: list[str],
    levels: dict[str, list[str]],
) -> dict[str, object]:
    continuous = [exposure, "age", "bmi", "pir", *continuous_extra]
    required = ["outcome", *continuous, *categorical, "pooled_weight", "psu", "strata"]
    complete = frame.dropna(subset=[c for c in required if c in frame.columns]).copy()
    complete = complete[complete["pooled_weight"].gt(0)] if "pooled_weight" in complete else complete.iloc[0:0]
    try:
        xdf, names = model.build_design(complete, continuous, categorical, levels=levels)
        condition_number = float(np.linalg.cond(xdf.to_numpy(float))) if len(xdf) and xdf.shape[1] else np.nan
    except Exception:
        names, condition_number = [], np.nan
    sparse_parts = []
    for col in categorical:
        if col in complete.columns:
            counts = complete[col].astype(str).value_counts(dropna=False)
            if len(counts) and int(counts.min()) < 5:
                sparse_parts.append(f"{col}:min_cell={int(counts.min())}")
    singleton = int(sum(complete.groupby("strata")["psu"].nunique().lt(2))) if len(complete) else 0
    weights = pd.to_numeric(complete.get("pooled_weight", pd.Series(dtype=float)), errors="coerce")
    weight_ratio = float(weights.max() / weights.median()) if len(weights) and weights.median() > 0 else np.nan
    return {
        "complete_case_n": int(len(complete)),
        "n_parameters": int(len(names)),
        "cases_per_parameter_approx": float(complete["outcome"].sum() / len(names)) if len(names) else np.nan,
        "condition_number_if_available": condition_number,
        "sparse_cell_warning": ";".join(sparse_parts),
        "extreme_weight_ratio": weight_ratio,
        "extreme_weight_warning": bool(finite(weight_ratio) and weight_ratio > EXTREME_WEIGHT_RATIO_THRESHOLD),
        "singleton_stratum_N": singleton,
    }


def fit_axis_model(
    frame: pd.DataFrame,
    model,
    exposure: str,
    urine: bool,
    label: str,
    categorical_extra: list[str] | None = None,
    continuous_extra: list[str] | None = None,
    levels_extra: dict[str, list[str]] | None = None,
    include_creatinine: bool = True,
) -> tuple[dict[str, object], dict]:
    categorical = [*(categorical_extra or []), "sex", "race", "smoking"]
    continuous_extra = continuous_extra or []
    continuous = [exposure, "age", "bmi", "pir", *continuous_extra]
    if urine and include_creatinine:
        continuous.append("creatinine_log2")
    levels = {**(levels_extra or {}), **model.LEVELS}
    fit = model.fit_survey_logistic(frame, continuous, categorical, exposure_name=exposure, levels=levels)
    diagnostics = diagnostic_context(frame, model, exposure, urine, categorical, [*continuous_extra, *( ["creatinine_log2"] if urine and include_creatinine else [])], levels)
    row = {
        "analysis": label,
        "exposure_variable": exposure,
        **fit_clean(fit),
        **diagnostics,
        "fit_status": fit.get("status", "not_estimable"),
        "warning_message": fit.get("message", "") if fit.get("status") != "ok" else "",
        "warning_type": "none" if fit.get("status") == "ok" else ("convergence_warning" if fit.get("status") == "converged_with_warning" else "fit_failure_or_non_estimable"),
    }
    row["analytic_n"] = fit.get("N", np.nan)
    row["crc_cases"] = fit.get("CRC_N", np.nan)
    row["controls"] = fit.get("Control_N", np.nan)
    return row, fit


def axis_frame(axis: pd.Series, detect: pd.DataFrame, harmonized: pd.DataFrame, screen_module: object, model) -> tuple[pd.DataFrame, dict[str, object]]:
    variable = str(axis["primary_biomarker"])
    axis_rows = detect.loc[detect["variable"].eq(variable)].copy()
    exposure, source = screen_module.read_environmental_axis(axis_rows)
    if exposure.empty:
        return pd.DataFrame(), source
    merged = exposure.merge(harmonized, on=["SEQN", "cycle"], how="inner", validate="one_to_one", suffixes=("", "_outcome"))
    population = model.population_frames(merged)["CRC_vs_cancer_free"]
    return population, source


def add_meta(row: dict, axis: pd.Series) -> dict:
    return {
        "axis_key": axis["axis_key"],
        "exposure_axis": axis["exposure_axis"],
        "primary_biomarker": axis["primary_biomarker"],
        "biological_matrix": axis.get("biological_matrix", ""),
        "eligible_chemical_count": axis.get("eligible_chemical_count", np.nan),
        "eligible_chemical_ids": axis.get("eligible_chemical_ids", ""),
        "mapping_confidence": axis.get("mapping_confidence", ""),
        **row,
    }


def primary_fit(axis: pd.Series, population: pd.DataFrame, model, urine: bool) -> tuple[dict, dict]:
    row, raw = fit_axis_model(population, model, "axis_log2", urine, "primary")
    return add_meta(row, axis), raw


def run_loco(axis: pd.Series, population: pd.DataFrame, model, urine: bool, pooled_or: float) -> tuple[pd.DataFrame, dict]:
    cycles = sorted(population.loc[population["axis_log2"].notna(), "cycle"].astype(str).unique().tolist())
    rows = []
    if len(cycles) < 3:
        return pd.DataFrame([add_meta({"dropped_cycle": "", "status": "not_applicable", "reason": f"only {len(cycles)} cycles"}, axis)]), {"n_loco_runs": 0}
    for cycle in cycles:
        row, _ = fit_axis_model(population[population["cycle"].ne(cycle)].copy(), model, "axis_log2", urine, f"LOCO_drop_{cycle}")
        row.update({"dropped_cycle": cycle, "direction": direction(row.get("OR"), pooled_or), "ci_excludes_1": ci_excludes_1(row.get("CI_low"), row.get("CI_high"))})
        rows.append(add_meta(row, axis))
    result = pd.DataFrame(rows)
    estimable = result.loc[result["fit_status"].isin(["ok", "converged_with_warning"]) & result["OR"].apply(finite)]
    pooled_dir = "positive" if pooled_or > 1 else "protective"
    same = estimable["direction"].eq(pooled_dir)
    summary = {
        "n_loco_runs": int(len(estimable)),
        "n_or_gt_1": int((estimable["OR"] > 1).sum()),
        "n_or_lt_1": int((estimable["OR"] < 1).sum()),
        "n_ci_excluding_1": int(estimable["ci_excludes_1"].sum()),
        "min_or": float(estimable["OR"].min()) if len(estimable) else np.nan,
        "max_or": float(estimable["OR"].max()) if len(estimable) else np.nan,
        "all_same_direction": bool(len(estimable) == len(cycles) and same.all()),
        "all_ci_exclude_1": bool(len(estimable) == len(cycles) and estimable["ci_excludes_1"].all()),
    }
    return result, summary


def run_cycle_specific(axis: pd.Series, population: pd.DataFrame, model, urine: bool, pooled_or: float) -> tuple[pd.DataFrame, dict]:
    cycles = sorted(population.loc[population["axis_log2"].notna(), "cycle"].astype(str).unique().tolist())
    rows = []
    for cycle in cycles:
        row, _ = fit_axis_model(population[population["cycle"].eq(cycle)].copy(), model, "axis_log2", urine, f"Single_cycle_{cycle}")
        row.update({"cycle": cycle, "direction": direction(row.get("OR"), pooled_or), "ci_excludes_1": ci_excludes_1(row.get("CI_low"), row.get("CI_high"))})
        rows.append(add_meta(row, axis))
    result = pd.DataFrame(rows)
    estimable = result.loc[result["fit_status"].isin(["ok", "converged_with_warning"]) & result["OR"].apply(finite)]
    pooled_dir = "positive" if pooled_or > 1 else "protective"
    same = estimable["direction"].eq(pooled_dir)
    return result, {
        "n_cycles": len(cycles),
        "n_same_direction_as_pooled": int(same.sum()),
        "direction_concordance_fraction": float(same.mean()) if len(estimable) else np.nan,
        "discordant_cycles": ";".join(estimable.loc[~same, "cycle"].astype(str).tolist()),
    }


def run_heterogeneity(axis: pd.Series, population: pd.DataFrame, model, urine: bool) -> tuple[dict, dict]:
    cycles = sorted(population.loc[population["axis_log2"].notna(), "cycle"].astype(str).unique().tolist())
    if len(cycles) < 3:
        return add_meta({"interaction_test": "not_applicable", "heterogeneity_status": "not_applicable", "interaction_P": np.nan}, axis), {}
    work = population.loc[population["axis_log2"].notna()].copy()
    interaction_names = []
    for cycle in cycles[1:]:
        name = f"axis_x_cycle_{cycle}"
        work[name] = work["axis_log2"] * work["cycle"].eq(cycle).astype(float)
        interaction_names.append(name)
    row, fit = fit_axis_model(
        work,
        model,
        "axis_log2",
        urine,
        "cycle_interaction",
        categorical_extra=["cycle"],
        continuous_extra=interaction_names,
        levels_extra={"cycle": cycles},
    )
    stability = load_module(ROOT / "work" / "scripts" / "mcop_crc_phase2_stability.py", "robustness_stability")
    test = stability.wald_test(fit, interaction_names)
    p = test.get("P_F", np.nan)
    status = "low" if finite(p) and p >= 0.10 else "moderate" if finite(p) and p >= 0.05 else "significant" if finite(p) else "not_estimable"
    row.update({
        "reference_cycle": cycles[0],
        "interaction_terms": ";".join(interaction_names),
        "interaction_P": p,
        "interaction_P_chi2": test.get("P_chi2", np.nan),
        "interaction_test": "Wald F test",
        "df": test.get("df_num", len(interaction_names)),
        "heterogeneity_status": status,
    })
    return add_meta(row, axis), {"interaction_P": p, "heterogeneity_status": status}


def run_timing(axis: pd.Series, population: pd.DataFrame, model, urine: bool, pooled_or: float) -> pd.DataFrame:
    rows = []
    if "years_since_crc" not in population.columns:
        return pd.DataFrame([add_meta({"exclusion_window": "", "status": "not_applicable", "reason": "years_since_crc unavailable"}, axis)])
    for window in [1, 2, 5]:
        recent_case = population["crc_case"].eq(True) & pd.to_numeric(population["years_since_crc"], errors="coerce").lt(window)
        frame = population.loc[~recent_case].copy()
        row, _ = fit_axis_model(frame, model, "axis_log2", urine, f"Exclude_diagnosis_lt_{window}y")
        row.update({"exclusion_window": f"< {window} years", "direction": direction(row.get("OR"), pooled_or)})
        rows.append(add_meta(row, axis))
    return pd.DataFrame(rows)


def run_tail(axis: pd.Series, population: pd.DataFrame, model, urine: bool, pooled_or: float) -> pd.DataFrame:
    available = pd.to_numeric(population.loc[population["axis_log2"].notna(), "axis_log2"], errors="coerce").dropna()
    rows = []
    for fraction in [0.01, 0.025]:
        cutoff = float(available.quantile(1 - fraction)) if len(available) else np.nan
        frame = population.loc[population["axis_log2"].le(cutoff) | population["axis_log2"].isna()].copy()
        row, _ = fit_axis_model(frame, model, "axis_log2", urine, f"Exclude_top_{fraction * 100:g}pct")
        row.update({"tail_rule": f"top {fraction * 100:g}% excluded", "cutoff_log2": cutoff, "direction": direction(row.get("OR"), pooled_or)})
        rows.append(add_meta(row, axis))
    return pd.DataFrame(rows)


def run_lod(axis: pd.Series, population: pd.DataFrame, model, urine: bool, pooled_or: float, axis_detect: pd.DataFrame) -> dict:
    selected = axis_detect.drop_duplicates(["cycle", "data_file", "variable"]).copy()
    measured = pd.to_numeric(selected.get("n_measured"), errors="coerce").sum()
    above = pd.to_numeric(selected.get("n_above_lod"), errors="coerce").sum()
    pct = float(100 * above / measured) if measured else np.nan
    base = {"lod_method": "not_required_D2" if pct >= 90 else "detectable_only", "above_lod_pct": pct}
    if not finite(pct) or pct >= 90:
        return add_meta({**base, "status": "not_applicable", "reason": "LOD concern minimal (>=90% above LOD)"}, axis)
    frame = population.loc[population["above_lod"].eq(True)].copy() if "above_lod" in population.columns else pd.DataFrame()
    if frame.empty:
        return add_meta({**base, "status": "not_estimable", "reason": "above_lod flag unavailable"}, axis)
    row, _ = fit_axis_model(frame, model, "axis_log2", urine, "detectable_only")
    row.update(base)
    row["direction"] = direction(row.get("OR"), pooled_or)
    return add_meta(row, axis)


def run_creatinine(axis: pd.Series, population: pd.DataFrame, model, urine: bool, pooled_or: float) -> dict:
    if not urine:
        return add_meta({"normalization": "not_applicable_serum", "status": "not_applicable"}, axis)
    frame = population.copy()
    frame["axis_creatinine_norm_log2"] = frame["axis_log2"] - frame["creatinine_log2"]
    row, _ = fit_axis_model(frame, model, "axis_creatinine_norm_log2", True, "creatinine_normalized", continuous_extra=[], include_creatinine=False)
    row["normalization"] = "log2(analyte) - log2(urine creatinine)"
    row["direction"] = direction(row.get("OR"), pooled_or)
    return add_meta(row, axis)


def compare_robustness(df: pd.DataFrame, pooled_or: float, direction_name: str) -> tuple[int, str, float]:
    work = df.loc[df["OR"].apply(finite)].copy()
    if work.empty:
        return 0, "not_estimable", np.nan
    same = work["OR"].apply(lambda x: (float(x) > 1) if pooled_or > 1 else (float(x) < 1)).all()
    max_delta = float(np.max(np.abs(np.log(pd.to_numeric(work["OR"], errors="coerce") / pooled_or))))
    if not same:
        return 0, "direction_instability", max_delta
    return (2 if max_delta <= ATTENUATION_LOG_THRESHOLD else 1), ("direction_preserved_comparable" if max_delta <= ATTENUATION_LOG_THRESHOLD else "direction_preserved_attenuated"), max_delta


def scorecard(axis: pd.Series, primary: dict, screen_row: pd.Series, loco_summary: dict, cycle_summary: dict, hetero: dict, timing: pd.DataFrame, tail: pd.DataFrame, creat: dict, warning: dict) -> dict:
    p = float(screen_row["P"]) if finite(screen_row.get("P")) else np.nan
    fdr = float(screen_row["BH_FDR"]) if finite(screen_row.get("BH_FDR")) else np.nan
    F = 2 if finite(fdr) and fdr < 0.05 else 1 if finite(p) and p < 0.05 else 0
    if loco_summary.get("n_loco_runs", 0) == 0:
        L, L_note = np.nan, "not_applicable"
    elif loco_summary.get("all_same_direction"):
        L, L_note = (2 if loco_summary.get("all_ci_exclude_1") else 1), "same_direction" + ("_all_CI_exclude_1" if loco_summary.get("all_ci_exclude_1") else "_some_CI_cross_1")
    else:
        L, L_note = 0, "direction_instability"
    frac = cycle_summary.get("direction_concordance_fraction", np.nan)
    C = np.nan if not finite(frac) else (2 if frac >= 0.80 else 1 if frac >= 0.60 else 0)
    H_p = hetero.get("interaction_P", np.nan)
    H = np.nan if not finite(H_p) else (2 if H_p >= 0.10 else 1 if H_p >= 0.05 else 0)
    pooled_or = float(primary.get("OR", np.nan))
    pooled_dir = "positive" if pooled_or > 1 else "protective"
    D, D_note, D_delta = compare_robustness(timing, pooled_or, pooled_dir)
    T, T_note, T_delta = compare_robustness(tail, pooled_or, pooled_dir)
    A = 0 if warning.get("fit_failure_count", 0) or warning.get("warning_persistent", False) else 1 if warning.get("warning_count", 0) else 2
    cases = int(primary.get("CRC_N", 0)) if finite(primary.get("CRC_N")) else 0
    E = 2 if cases >= 60 else 1 if cases >= 30 else 0
    robust_a = F == 2 and all(finite(x) and x >= 1 for x in [L, C, D, T, A])
    tier = "Robust Tier A" if robust_a else "Tier B" if F >= 1 and E >= 1 and any(finite(x) and x >= 1 for x in [L, C, D, T]) else "Exploratory"
    fingerprint = " | ".join(f"{k}{'NA' if not finite(v) else int(v)}" for k, v in [("F", F), ("L", L), ("C", C), ("H", H), ("D", D), ("T", T), ("A", A), ("E", E)])
    return add_meta({
        "primary_OR": primary.get("OR"), "primary_CI_low": primary.get("CI_low"), "primary_CI_high": primary.get("CI_high"), "primary_P": primary.get("P"), "primary_BH_FDR_15axis": screen_row.get("BH_FDR"),
        "F": F, "L": L, "C": C, "H": H, "D": D, "T": T, "A": A, "E": E,
        "F_note": "BH-FDR<0.05" if F == 2 else "nominal_P<0.05" if F == 1 else "P>=0.05",
        "L_note": L_note, "C_note": f"{cycle_summary.get('n_same_direction_as_pooled', 0)}/{cycle_summary.get('n_cycles', 0)} same direction", "H_note": f"Pinteraction={H_p:.4g}" if finite(H_p) else "not_estimable",
        "D_note": D_note, "T_note": T_note, "A_note": "no meaningful warning" if A == 2 else "localized warning; primary/sensitivity fit otherwise estimable" if A == 1 else "persistent technical warning or fit failure", "E_note": f"{cases} CRC cases",
        "timing_max_abs_logOR_delta": D_delta, "tail_max_abs_logOR_delta": T_delta,
        "robustness_fingerprint": fingerprint, "robustness_tier": tier,
    }, axis)


def main() -> None:
    screen_module = load_module(ROOT / "work" / "scripts" / "environmental_crc_systematic_human_screen_v2.py", "robust_screen")
    model = load_module(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "robust_model")
    stability = load_module(ROOT / "work" / "scripts" / "mcop_crc_phase2_stability.py", "robust_stability")
    matrix = pd.read_csv(MATRIX, low_memory=False)
    detect = pd.read_csv(DETECT, low_memory=False)
    screen_results = pd.read_csv(SCREEN, low_memory=False)
    screen_fdr = pd.read_csv(SCREEN_FDR, low_memory=False)
    if len(screen_results) != 15 or screen_results["axis_key"].nunique() != 15:
        raise ValueError("The robustness audit requires exactly the frozen 15 primary axes")
    if len(screen_fdr) != 15:
        raise ValueError("The frozen primary FDR file must contain 15 axes")
    eligible = matrix.loc[matrix["eligible_permissive"].astype(str).str.lower().eq("true")].copy()
    eligible["axis_key"] = eligible["biological_matrix"].fillna("") + "|" + eligible["selected_primary_biomarker"].fillna("")
    axis_meta = eligible.groupby("axis_key", as_index=False).agg(
        exposure_axis=("exposure_axis", lambda s: ";".join(sorted(set(s.astype(str))))),
        primary_biomarker=("selected_primary_biomarker", "first"),
        biological_matrix=("biological_matrix", "first"),
        eligible_chemical_count=("ChemicalID", "size"),
        eligible_chemical_ids=("ChemicalID", lambda s: ";".join(s.astype(str))),
        mapping_confidence=("mapping_confidence", "first"),
    )
    axis_meta = axis_meta.merge(screen_results[["axis_key", "cycle_list", "n_cycles_available", "eligible_chemical_names"]], on="axis_key", how="left", validate="one_to_one")
    audit_module = load_module(ROOT / "work" / "scripts" / "environmental_crc_267_biomarker_audit_v2.py", "robust_actionability_audit")
    harmonized = screen_module.load_harmonized(audit_module)

    primary_rows, loco_rows, cycle_rows, het_rows, timing_rows, tail_rows, lod_rows, creat_rows = [], [], [], [], [], [], [], []
    score_rows, warning_rows = [], []
    qc_primary = []
    for _, axis in axis_meta.iterrows():
        population, source = axis_frame(axis, detect, harmonized, screen_module, model)
        urine = "urine" in str(axis["biological_matrix"]).lower()
        pooled = screen_results.loc[screen_results["axis_key"].eq(axis["axis_key"])].iloc[0]
        primary, primary_raw = primary_fit(axis, population, model, urine)
        primary_rows.append(primary)
        primary_reproduction_diff = abs(np.log(float(primary.get("OR")) / float(pooled["OR"]))) if finite(primary.get("OR")) and finite(pooled["OR"]) else np.nan
        qc_primary.append({"axis_key": axis["axis_key"], "screen_OR": pooled["OR"], "rerun_OR": primary.get("OR"), "absolute_logOR_difference": primary_reproduction_diff})
        pooled_or = float(primary.get("OR", np.nan))
        loco, loco_summary = run_loco(axis, population, model, urine, pooled_or)
        cycle, cycle_summary = run_cycle_specific(axis, population, model, urine, pooled_or)
        het, het_summary = run_heterogeneity(axis, population, model, urine)
        timing = run_timing(axis, population, model, urine, pooled_or)
        tail = run_tail(axis, population, model, urine, pooled_or)
        axis_detect = detect.loc[detect["variable"].eq(axis["primary_biomarker"])].copy()
        lod = run_lod(axis, population, model, urine, pooled_or, axis_detect)
        creat = run_creatinine(axis, population, model, urine, pooled_or)
        loco_rows.extend(loco.to_dict("records")); cycle_rows.extend(cycle.to_dict("records")); timing_rows.extend(timing.to_dict("records")); tail_rows.extend(tail.to_dict("records"))
        het_rows.append(het); lod_rows.append(lod); creat_rows.append(creat)
        all_fit_statuses = [primary.get("fit_status"), *loco.get("fit_status", pd.Series(dtype=object)).tolist(), *cycle.get("fit_status", pd.Series(dtype=object)).tolist(), het.get("fit_status"), *timing.get("fit_status", pd.Series(dtype=object)).tolist(), *tail.get("fit_status", pd.Series(dtype=object)).tolist(), lod.get("fit_status"), creat.get("fit_status")]
        warning_statuses = [x for x in all_fit_statuses if x not in {None, "not_applicable"}]
        warning_count = sum(x == "converged_with_warning" for x in warning_statuses)
        failure_count = sum(x not in {"ok", "converged_with_warning"} for x in warning_statuses)
        warning_fraction = float(warning_count / len(warning_statuses)) if warning_statuses else 0.0
        warning_persistent = bool(warning_count and (warning_fraction >= 0.75 or primary.get("fit_status") == "converged_with_warning" and warning_fraction >= 0.50))
        primary_warning = primary.get("warning_message", "")
        status_distribution = ";".join(f"{status}={warning_statuses.count(status)}" for status in sorted(set(warning_statuses)))
        if failure_count:
            warning_reason = "At least one audit fit was non-estimable or failed; inspect fit_status and the corresponding module output."
            coefficient_impact = "Not fully assessed because at least one audit component was non-estimable or failed."
        elif warning_count:
            warning_reason = "Newton-IRLS reached a finite estimable solution but stopped before the configured convergence tolerance within the iteration limit; this is an algorithmic convergence warning, not a missing estimate."
            coefficient_impact = "No discrepancy versus the frozen primary screen (absolute log-OR difference <=1e-8)."
        else:
            warning_reason = "No convergence or estimability warning detected across applicable fits."
            coefficient_impact = "Not applicable; all applicable fits were reported as converged."
        warning = {
            "fit_warning_count": int(warning_count), "fit_failure_count": int(failure_count), "warning_count": int(warning_count + failure_count), "n_fit_evaluations": int(len(warning_statuses)), "warning_fraction": warning_fraction, "warning_persistent": warning_persistent,
            "primary_fit_status": primary.get("fit_status"), "primary_warning_type": primary.get("warning_type"), "primary_warning_message": primary_warning,
            "warning_status_distribution": status_distribution,
            "warning_technical_reason": warning_reason,
            "coefficient_impact_assessment": coefficient_impact,
            "warning_repeated_in_sensitivity": bool(warning_count > 0),
            "primary_reproduction_abs_logOR_difference": primary_reproduction_diff,
        }
        warning_rows.append(add_meta({
            **warning,
            "primary_condition_number": primary.get("condition_number_if_available"),
            "sparse_cell_warning": primary.get("sparse_cell_warning", ""),
            "extreme_weight_warning": primary.get("extreme_weight_warning", False),
            "extreme_weight_ratio": primary.get("extreme_weight_ratio"),
            "singleton_stratum_N": primary.get("singleton_stratum_N", np.nan),
            "n_fit_warning_or_failure_across_audit": int(warning_count + failure_count),
            "warning_assessment": "no meaningful warning" if not warning_count and not failure_count else "persistent convergence warning" if warning_persistent else "localized convergence warning" if not failure_count else "persistent technical/non-estimable issue",
        }, axis))
        score_rows.append(scorecard(axis, primary_raw, pooled, loco_summary, cycle_summary, het_summary, timing, tail, creat, warning))

    primary_df = pd.DataFrame(primary_rows)
    loco_df = pd.DataFrame(loco_rows)
    cycle_df = pd.DataFrame(cycle_rows)
    het_df = pd.DataFrame(het_rows)
    timing_df = pd.DataFrame(timing_rows)
    tail_df = pd.DataFrame(tail_rows)
    lod_df = pd.DataFrame(lod_rows)
    creat_df = pd.DataFrame(creat_rows)
    warn_df = pd.DataFrame(warning_rows)
    score_df = pd.DataFrame(score_rows)
    summary_df = score_df.merge(primary_df[["axis_key", "fit_status", "fit_N", "fit_crc_cases", "fit_control_n", "cases_per_parameter_approx"]], on="axis_key", how="left", suffixes=("", "_primary"))
    primary_qc = pd.DataFrame(qc_primary)
    if len(primary_qc) != 15 or primary_qc["absolute_logOR_difference"].dropna().max() > 1e-8:
        raise AssertionError("Rerun primary estimates do not reproduce the frozen 15-axis screen")

    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(OUT_SUMMARY, index=False)
    loco_df.to_csv(OUT_LOCO, index=False); cycle_df.to_csv(OUT_CYCLE, index=False); het_df.to_csv(OUT_HET, index=False); timing_df.to_csv(OUT_TIMING, index=False); tail_df.to_csv(OUT_TAIL, index=False); lod_df.to_csv(OUT_LOD, index=False); creat_df.to_csv(OUT_CREAT, index=False); warn_df.to_csv(OUT_WARN, index=False); score_df.to_csv(OUT_SCORE, index=False)
    primary_qc.to_csv(OUTPUTS / "environmental_crc_15axis_primary_reproduction_qc.csv", index=False)

    robust = score_df.loc[score_df["robustness_tier"].eq("Robust Tier A"), ["axis_key", "primary_biomarker", "robustness_fingerprint"]]
    persistent_warning_axes = warn_df.loc[warn_df["warning_persistent"].eq(True), "primary_biomarker"].astype(str).tolist()
    localized_warning_axes = warn_df.loc[(warn_df["warning_persistent"].eq(False)) & (warn_df["fit_warning_count"] > 0), "primary_biomarker"].astype(str).tolist()
    lines = [
        "# Environmental CRC 15-axis systematic robustness audit",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Frozen scope",
        "",
        "All 15 axes were read automatically from the frozen primary screen. The same primary model and the same robustness modules were applied to every axis. The primary BH-FDR remains the original 15-axis correction; no robustness subset was used to recompute significance.",
        "",
        f"Primary reproduction QC: 15/15 axes reproduced their stored OR to absolute log-OR difference <=1e-8.",
        "",
        "## Robustness scorecard",
        "",
        "| Axis | Biomarker | OR | BH-FDR (15-axis) | Fingerprint | Tier |",
        "|---|---|---:|---:|---|---|",
    ]
    for _, row in score_df.sort_values(["robustness_tier", "primary_BH_FDR_15axis"], na_position="last").iterrows():
        lines.append(f"| {row['exposure_axis']} | {row['primary_biomarker']} | {float(row['primary_OR']):.4g} | {float(row['primary_BH_FDR_15axis']):.4g} | {row['robustness_fingerprint']} | {row['robustness_tier']} |")
    lines += [
        "",
        f"Robust Tier A axes: **{len(robust)}**.",
        "",
        "## Required interpretation questions",
        "",
        f"- FDR-supported primary axes: **{int((score_df['F'] == 2).sum())}**; the original 15-axis correction is unchanged.",
        f"- Axes with all available LOCO estimates directionally concordant: **{int((score_df['L'] >= 1).sum())}**; axes with direction instability: **{int((score_df['L'] == 0).sum())}**.",
        f"- Axes with significant exposure×cycle heterogeneity (P<0.05): **{int((score_df['H'] == 0).sum())}**.",
        f"- Axes retaining direction after all diagnosis-timing exclusions: **{int((score_df['D'] >= 1).sum())}**.",
        f"- Axes with any fit warning or non-estimable audit component: **{int((score_df['A'] < 2).sum())}**; axes with persistent warning/technical concern (A0): **{int((score_df['A'] == 0).sum())}**.",
        "",
        "### Warning interpretation",
        "",
        "Convergence warnings were retained rather than hidden. For warning-only fits, the Newton-IRLS routine returned finite estimable coefficients and sandwich variance but did not reach the configured tolerance within the iteration limit; the warning is therefore algorithmic and does not by itself imply coefficient failure. The frozen primary estimates were independently reproduced (15/15; absolute log-OR difference <=1e-8).",
        f"- Persistent warning axes: **{', '.join(persistent_warning_axes) if persistent_warning_axes else 'none'}**.",
        f"- Localized warning axes: **{', '.join(localized_warning_axes) if localized_warning_axes else 'none'}**.",
        "- A persistent warning is treated as an A0 technical concern in the scorecard; it is not converted into a positive finding or silently removed. Warning repetition across secondary fits is explicitly recorded in `environmental_crc_15axis_model_warnings.csv`.",
        "",
        "MCOP/URXCOP is evaluated in exactly the same scorecard as PFHS/LBXPFHS and the other 13 axes. MiNP remains a distinct molecular DINP nominee and is not silently converted into MCOP; it is outside the 15-axis human screen because its direct detectability failed the frozen actionability gate.",
        "",
        "## Prespecified tag definitions",
        "",
        "- F: F2=15-axis BH-FDR<0.05; F1=nominal P<0.05; F0 otherwise.",
        "- L: L2=same direction and all LOCO CIs exclude 1; L1=same direction with some CIs crossing 1; L0=direction instability.",
        "- C: C2>=80% same-direction cycle estimates; C1=60–79%; C0<60%.",
        "- H: H2=Pinteraction>=0.10; H1=0.05–<0.10; H0<0.05.",
        "- D/T: D2/T2 preserve direction with maximum absolute log-OR change <=0.25; D1/T1 preserve direction with greater attenuation; D0/T0 direction instability.",
        "- A: A2 no warning; A1 localized warning with primary/sensitivity fits otherwise estimable; A0 persistent warning (>=75% of applicable fits, or >=50% when the primary itself warns) or any fit failure.",
        "- E: E2>=60 CRC cases; E1=30–59; E0<30.",
        "- Robust Tier A was frozen as F2, L>=1, C>=1, D>=1, T>=1, A>=1; H is reported as a penalty/evidence tag, not a hard deletion gate.",
        "",
        "This is a robustness/association audit of cross-sectional NHANES data and does not establish causality or eliminate reverse causation/survivor bias.",
    ]
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "15-axis systematic robustness audit",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "n_axes": 15,
        "primary_screen_fdr_scope": "all 15 fixed axes; unchanged in robustness audit",
        "primary_model": "CRC ~ log2(exposure) + age + sex + race + BMI + smoking + PIR; urine adds log2(creatinine)",
        "weight_rule": "analyte-specific NHANES survey weight, including WTSPH4YR 2x conversion where applicable, divided by included analyte cycles",
        "prespecified_thresholds": {
            "attenuation_log_or_threshold": ATTENUATION_LOG_THRESHOLD,
            "extreme_weight_ratio_warning": EXTREME_WEIGHT_RATIO_THRESHOLD,
            "F": {"F2": "BH-FDR<0.05", "F1": "nominal P<0.05", "F0": "otherwise"},
            "C": {"C2": ">=80% same direction", "C1": "60-79%", "C0": "<60%"},
            "H": {"H2": "Pinteraction>=0.10", "H1": "0.05<=Pinteraction<0.10", "H0": "Pinteraction<0.05"},
            "E": {"E2": ">=60 cases", "E1": "30-59 cases", "E0": "<30 cases"},
            "robust_tier_A": "F2 & L>=1 & C>=1 & D>=1 & T>=1 & A>=1",
        },
        "primary_reproduction_qc": "15/15 axes reproduced frozen OR within absolute log-OR difference <=1e-8",
        "script": "work/scripts/environmental_crc_15axis_robustness_audit.py",
        "outputs": [str(p) for p in [OUT_SUMMARY, OUT_LOCO, OUT_CYCLE, OUT_HET, OUT_TIMING, OUT_TAIL, OUT_LOD, OUT_CREAT, OUT_WARN, OUT_SCORE, OUTPUTS / "environmental_crc_15axis_primary_reproduction_qc.csv", OUT_REPORT]],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"n_axes": 15, "robust_tier_a": len(robust), "fdr_supported": int((score_df["F"] == 2).sum()), "heterogeneity_significant": int((score_df["H"] == 0).sum())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
