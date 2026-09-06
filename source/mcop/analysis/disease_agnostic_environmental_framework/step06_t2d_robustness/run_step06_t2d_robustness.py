"""Uniform Step 6 audit for all FDR-positive tests in the frozen T2D screen.

The script is deliberately T2D-specific. It rebuilds the assay-specific
exposure populations against the frozen, outcome-independent T2D outcome
frame, then applies the same stability modules to every FDR-positive test.
Exposure correlations are descriptive and are computed without T2D status or
outcome association statistics.
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
STEP5 = FRAMEWORK / "step05_t2d_screen"
DEFAULT_TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
DEFAULT_REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_OUT = FRAMEWORK / "step06_t2d_robustness"
PLAN = DEFAULT_OUT / "STEP6_T2D_ROBUSTNESS_PLAN.md"
FDR_DENOMINATOR = 29
TAIL_FRACTIONS = (0.01, 0.025)
ATTENUATION_LOG_THRESHOLD = 0.25
CORRELATION_THRESHOLD = 0.70
CORRELATION_MIN_N = 500


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


def numeric(value: object) -> pd.Series:
    return pd.to_numeric(value, errors="coerce")


def fit_public(fit: dict) -> dict:
    return {key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}}


def ci_excludes_1(low: object, high: object) -> bool:
    return finite(low) and finite(high) and (float(low) > 1 or float(high) < 1)


def direction(or_value: object, pooled_or: float) -> str:
    if not finite(or_value) or not finite(pooled_or):
        return "not_estimable"
    value = float(or_value)
    if pooled_or > 1:
        return "same" if value > 1 else "discordant"
    if pooled_or < 1:
        return "same" if value < 1 else "discordant"
    return "same" if value == 1 else "discordant"


def add_meta(row: dict, test_row: pd.Series) -> dict:
    return {
        "test_id": str(test_row["test_id"]),
        "biomarker": str(test_row["biomarker"]),
        "variable": str(test_row["variable"]),
        "exposure_axis": str(test_row["exposure_axes"]),
        "matrix": str(test_row["matrix"]),
        **row,
    }


def fit_frame(frame: pd.DataFrame, model, test_row: pd.Series, label: str,
              urine: bool | None = None, extra_continuous: list[str] | None = None,
              extra_categorical: list[str] | None = None,
              levels_extra: dict[str, list[str]] | None = None,
              base_categorical: list[str] | None = None,
              exposure_name: str = "axis_log2",
              include_creatinine: bool = True) -> tuple[dict, dict]:
    urine = str(test_row["matrix"]).lower() == "urine" if urine is None else urine
    extra_continuous = list(extra_continuous or [])
    extra_categorical = list(extra_categorical or [])
    continuous = [exposure_name, "age", "bmi", "pir", *extra_continuous]
    if urine and include_creatinine:
        continuous.append("creatinine_log2")
    levels = {**(levels_extra or {}), **model.LEVELS}
    base_categorical = list(base_categorical or ["sex", "race", "smoking"])
    fit = model.fit_survey_logistic(
        frame,
        continuous,
        [*extra_categorical, *base_categorical],
        outcome="outcome",
        exposure_name=exposure_name,
        levels=levels,
    )
    row = {
        "analysis": label,
        "exposure_variable": exposure_name,
        **fit_public(fit),
        "fit_status": fit.get("status", "not_estimable"),
        "warning_message": fit.get("message", "") if fit.get("status") != "ok" else "",
        "analytic_n": fit.get("N", np.nan),
        "t2d_cases": fit.get("CRC_N", np.nan),
        "controls": fit.get("Control_N", np.nan),
    }
    return row, fit


def prepare_population(exposure: pd.DataFrame, outcome: pd.DataFrame, test_row: pd.Series) -> pd.DataFrame:
    if exposure.empty:
        return pd.DataFrame()
    merged = exposure.merge(
        outcome,
        on=["SEQN", "cycle"],
        how="inner",
        validate="one_to_one",
        suffixes=("", "_outcome"),
    )
    population = merged.loc[merged["t2d_eligible"]].copy()
    population["outcome"] = population["t2d_case"].astype(int)
    if str(test_row["matrix"]).lower() == "urine":
        population["creatinine_log2"] = np.log2(numeric(population["URXUCR"]).where(numeric(population["URXUCR"]).gt(0)))
    return population


def wald_test(fit: dict, coefficient_names: list[str], stability) -> dict:
    if fit.get("status") not in {"ok", "converged_with_warning"}:
        return {"status": "not_estimable", "P_F": np.nan, "P_chi2": np.nan}
    coefficients = fit.get("coefficients", {})
    covariance = fit.get("covariance")
    if covariance is None or any(name not in coefficients for name in coefficient_names):
        return {"status": "not_estimable", "P_F": np.nan, "P_chi2": np.nan}
    return stability.wald_test(fit, coefficient_names)


def audit_primary(population: pd.DataFrame, model, test_row: pd.Series, stored: pd.Series) -> tuple[dict, dict]:
    row, fit = fit_frame(population, model, test_row, "primary")
    stored_or = float(stored["OR_per_doubling"])
    rerun_or = float(row.get("OR", np.nan))
    diff = abs(float(np.log(rerun_or / stored_or))) if finite(rerun_or) and stored_or > 0 else np.nan
    row.update({"step5_or": stored_or, "step5_p": float(stored["P"]), "step5_bh_fdr": float(stored["BH_FDR"]), "absolute_log_or_difference": diff})
    return add_meta(row, test_row), fit


def run_loco(population: pd.DataFrame, model, test_row: pd.Series, pooled_or: float) -> pd.DataFrame:
    rows = []
    for cycle in sorted(population["cycle"].astype(str).unique()):
        row, _ = fit_frame(population.loc[population["cycle"].ne(cycle)].copy(), model, test_row, f"LOCO_drop_{cycle}")
        row.update({"dropped_cycle": cycle, "direction": direction(row.get("OR"), pooled_or), "ci_excludes_1": ci_excludes_1(row.get("CI_low"), row.get("CI_high"))})
        rows.append(add_meta(row, test_row))
    return pd.DataFrame(rows)


def run_cycle_specific(population: pd.DataFrame, model, test_row: pd.Series, pooled_or: float) -> pd.DataFrame:
    rows = []
    for cycle in sorted(population["cycle"].astype(str).unique()):
        row, _ = fit_frame(population.loc[population["cycle"].eq(cycle)].copy(), model, test_row, f"Single_cycle_{cycle}")
        row.update({"cycle": cycle, "direction": direction(row.get("OR"), pooled_or), "ci_excludes_1": ci_excludes_1(row.get("CI_low"), row.get("CI_high"))})
        rows.append(add_meta(row, test_row))
    return pd.DataFrame(rows)


def run_cycle_interaction(population: pd.DataFrame, model, test_row: pd.Series, stability) -> dict:
    cycles = sorted(population["cycle"].astype(str).unique())
    work = population.copy()
    interaction_names = []
    for cycle in cycles[1:]:
        name = f"axis_x_cycle_{cycle}"
        work[name] = work["axis_log2"] * work["cycle"].eq(cycle).astype(float)
        interaction_names.append(name)
    row, fit = fit_frame(
        work, model, test_row, "cycle_interaction",
        extra_continuous=interaction_names,
        extra_categorical=["cycle"],
        levels_extra={"cycle": cycles},
    )
    test = wald_test(fit, interaction_names, stability)
    row.update({
        "reference_cycle": cycles[0] if cycles else "",
        "interaction_terms": ";".join(interaction_names),
        "interaction_P_F": test.get("P_F", np.nan),
        "interaction_P_chi2": test.get("P_chi2", np.nan),
        "interaction_df_num": test.get("df_num", np.nan),
        "interaction_df_denom": test.get("df_denom", np.nan),
    })
    return add_meta(row, test_row)


def run_tail(population: pd.DataFrame, model, test_row: pd.Series, pooled_or: float) -> pd.DataFrame:
    available = numeric(population["axis_log2"]).dropna()
    rows = []
    for fraction in TAIL_FRACTIONS:
        cutoff = float(available.quantile(1 - fraction))
        frame = population.loc[population["axis_log2"].le(cutoff) | population["axis_log2"].isna()].copy()
        row, _ = fit_frame(frame, model, test_row, f"Exclude_top_{fraction * 100:g}pct")
        row.update({"tail_fraction": fraction, "cutoff_log2": cutoff, "direction": direction(row.get("OR"), pooled_or)})
        rows.append(add_meta(row, test_row))
    return pd.DataFrame(rows)


def run_creatinine(population: pd.DataFrame, model, test_row: pd.Series, pooled_or: float) -> dict:
    if str(test_row["matrix"]).lower() != "urine":
        return add_meta({"analysis": "creatinine_normalized", "fit_status": "not_applicable", "normalization": "not_applicable_serum"}, test_row)
    frame = population.copy()
    frame["axis_creatinine_norm_log2"] = frame["axis_log2"] - frame["creatinine_log2"]
    row, _ = fit_frame(frame, model, test_row, "creatinine_normalized", urine=True, exposure_name="axis_creatinine_norm_log2", include_creatinine=False)
    row.update({"normalization": "log2(analyte) - log2(urinary creatinine)", "direction": direction(row.get("OR"), pooled_or)})
    return add_meta(row, test_row)


def run_age40(population: pd.DataFrame, model, test_row: pd.Series, pooled_or: float) -> dict:
    row, _ = fit_frame(population.loc[population["age"].ge(40)].copy(), model, test_row, "age_ge_40")
    row["direction"] = direction(row.get("OR"), pooled_or)
    return add_meta(row, test_row)


def run_sex(population: pd.DataFrame, model, test_row: pd.Series, pooled_or: float) -> pd.DataFrame:
    rows = []
    for sex in ["Female", "Male"]:
        work = population.loc[population["sex"].eq(sex)].copy()
        row, _ = fit_frame(work, model, test_row, f"sex_{sex}", extra_categorical=[], base_categorical=["race", "smoking"])
        row.update({"sex_group": sex, "direction": direction(row.get("OR"), pooled_or)})
        rows.append(add_meta(row, test_row))
    return pd.DataFrame(rows)


def run_sex_interaction(population: pd.DataFrame, model, test_row: pd.Series, stability) -> dict:
    work = population.copy()
    work["axis_x_sex_Male"] = work["axis_log2"] * work["sex"].eq("Male").astype(float)
    row, fit = fit_frame(work, model, test_row, "sex_interaction", extra_continuous=["axis_x_sex_Male"])
    test = wald_test(fit, ["axis_x_sex_Male"], stability)
    row.update({"interaction_term": "axis_x_sex_Male", "interaction_P_F": test.get("P_F", np.nan), "interaction_P_chi2": test.get("P_chi2", np.nan)})
    return add_meta(row, test_row)


def run_lod(population: pd.DataFrame, model, test_row: pd.Series, registry: pd.DataFrame, pooled_or: float) -> dict:
    variable = str(test_row["variable"])
    rows = registry.loc[registry["variable"].eq(variable)].copy().drop_duplicates(["cycle", "data_file", "variable"])
    measured = numeric(rows.get("n_measured", pd.Series(dtype=float))).sum()
    above = numeric(rows.get("n_above_lod", pd.Series(dtype=float))).sum()
    pct = float(100 * above / measured) if measured else np.nan
    base = {"above_lod_pct": pct, "lod_threshold_rule": "sensitivity only if <90% above LOD"}
    if not finite(pct) or pct >= 90:
        return add_meta({**base, "analysis": "lod", "fit_status": "not_applicable", "reason": "LOD concern minimal"}, test_row)
    frame = population.loc[population["above_lod"].eq(True)].copy()
    row, _ = fit_frame(frame, model, test_row, "detectable_only")
    row.update({**base, "analysis": "lod", "direction": direction(row.get("OR"), pooled_or)})
    return add_meta(row, test_row)


def summarize_loco(loco: pd.DataFrame) -> dict:
    est = loco.loc[loco["fit_status"].isin(["ok", "converged_with_warning"]) & loco["OR"].apply(finite)].copy()
    same = est["direction"].eq("same")
    all_same = len(est) == len(loco) and bool(same.all())
    all_ci = all_same and bool(est["ci_excludes_1"].all())
    return {
        "loco_n": int(len(est)),
        "loco_same_n": int(same.sum()),
        "loco_direction_fraction": float(same.mean()) if len(est) else np.nan,
        "loco_all_same_direction": all_same,
        "loco_all_ci_exclude_1": all_ci,
        "loco_min_or": float(est["OR"].min()) if len(est) else np.nan,
        "loco_max_or": float(est["OR"].max()) if len(est) else np.nan,
    }


def summarize_cycle(cycle: pd.DataFrame) -> dict:
    est = cycle.loc[cycle["fit_status"].isin(["ok", "converged_with_warning"]) & cycle["OR"].apply(finite)].copy()
    same = est["direction"].eq("same")
    fraction = float(same.mean()) if len(est) else np.nan
    return {
        "cycle_n": int(len(est)),
        "cycle_same_n": int(same.sum()),
        "cycle_direction_fraction": fraction,
        "cycle_discordant": ";".join(est.loc[~same, "cycle"].astype(str).tolist()),
        "cycle_min_or": float(est["OR"].min()) if len(est) else np.nan,
        "cycle_max_or": float(est["OR"].max()) if len(est) else np.nan,
    }


def summarize_tail(tail: pd.DataFrame, pooled_or: float) -> dict:
    est = tail.loc[tail["fit_status"].isin(["ok", "converged_with_warning"]) & tail["OR"].apply(finite)].copy()
    direction_ok = est["direction"].eq("same")
    deltas = np.abs(np.log(numeric(est["OR"]) / pooled_or)) if len(est) else pd.Series(dtype=float)
    max_delta = float(deltas.max()) if len(deltas) else np.nan
    if len(est) < len(TAIL_FRACTIONS) or not direction_ok.all():
        score, note = 0, "direction_instability_or_not_estimable"
    elif max_delta <= ATTENUATION_LOG_THRESHOLD:
        score, note = 2, "direction_preserved_comparable"
    else:
        score, note = 1, "direction_preserved_attenuated"
    return {"tail_score": score, "tail_note": note, "tail_max_abs_log_or_delta": max_delta}


def score_candidate(primary: dict, loco: dict, cycle: dict, hetero: dict, tail: dict, statuses: list[str], screen: pd.Series, corr_cluster: dict) -> dict:
    fdr = float(screen["BH_FDR"])
    p = float(screen["P"])
    F = 2 if fdr < 0.05 else 1 if p < 0.05 else 0
    if loco["loco_n"] == 0 or not loco["loco_all_same_direction"]:
        L, L_note = 0, "direction_instability_or_not_estimable"
    elif loco["loco_all_ci_exclude_1"]:
        L, L_note = 2, "all_direction_same_and_all_CI_exclude_1"
    else:
        L, L_note = 1, "all_direction_same_some_CI_cross_1"
    frac = cycle["cycle_direction_fraction"]
    C = np.nan if not finite(frac) else 2 if frac >= 0.80 else 1 if frac >= 0.60 else 0
    H_p = hetero.get("interaction_P_F", np.nan)
    H = np.nan if not finite(H_p) else 2 if H_p >= 0.10 else 1 if H_p >= 0.05 else 0
    T = tail["tail_score"]
    applicable_statuses = [status for status in statuses if status != "not_applicable"]
    warning_n = sum(status == "converged_with_warning" for status in applicable_statuses)
    failures = sum(status not in {"ok", "converged_with_warning"} for status in applicable_statuses)
    A = 0 if failures else 1 if warning_n else 2
    robust = bool(F == 2 and L >= 1 and finite(C) and C >= 1 and T >= 1 and A >= 1)
    tier = "robust_fdr_candidate" if robust else "fdr_candidate_with_instability"
    return {
        **primary,
        **loco,
        **cycle,
        "P": p,
        "BH_FDR": fdr,
        "F": F,
        "L": L,
        "L_note": L_note,
        "C": C,
        "H": H,
        "H_note": f"Pinteraction={H_p:.6g}" if finite(H_p) else "not_estimable",
        "T": T,
        "T_note": tail["tail_note"],
        "tail_max_abs_log_or_delta": tail["tail_max_abs_log_or_delta"],
        "A": A,
        "technical_status_count": len(applicable_statuses),
        "technical_warning_count": warning_n,
        "technical_failure_count": failures,
        "cluster_id": corr_cluster.get("cluster_id", ""),
        "max_abs_cycle_adjusted_rho": corr_cluster.get("max_abs_cycle_adjusted_rho", np.nan),
        "cluster_n": corr_cluster.get("cluster_n", 1),
        "priority_tier": tier,
        "robustness_fingerprint": f"F{F}|L{L}|C{C if finite(C) else 'NA'}|H{H if finite(H) else 'NA'}|T{T}|A{A}",
    }


def pairwise_correlations(exposures: dict[str, pd.DataFrame], test_rows: dict[str, pd.Series]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict]]:
    variables = sorted(exposures)
    matrix = pd.DataFrame(np.eye(len(variables)), index=variables, columns=variables)
    edges = []
    for i, left in enumerate(variables):
        a = exposures[left][["SEQN", "cycle", "axis_log2"]].rename(columns={"axis_log2": "left_value"})
        for right in variables[i + 1:]:
            b = exposures[right][["SEQN", "cycle", "axis_log2"]].rename(columns={"axis_log2": "right_value"})
            joined = a.merge(b, on=["SEQN", "cycle"], how="inner")
            joined = joined.dropna(subset=["left_value", "right_value"])
            pooled_rho = joined["left_value"].corr(joined["right_value"], method="spearman") if len(joined) >= 3 else np.nan
            if len(joined):
                left_resid = joined["left_value"] - joined.groupby("cycle")["left_value"].transform("mean")
                right_resid = joined["right_value"] - joined.groupby("cycle")["right_value"].transform("mean")
                cycle_rho = left_resid.corr(right_resid, method="spearman") if len(joined) >= 3 else np.nan
            else:
                cycle_rho = np.nan
            matrix.loc[left, right] = matrix.loc[right, left] = cycle_rho if finite(cycle_rho) else np.nan
            edges.append({
                "left_variable": left,
                "right_variable": right,
                "left_biomarker": str(test_rows[left]["biomarker"]),
                "right_biomarker": str(test_rows[right]["biomarker"]),
                "pairwise_n": int(len(joined)),
                "pooled_spearman_rho": pooled_rho,
                "cycle_adjusted_spearman_rho": cycle_rho,
                "abs_cycle_adjusted_rho": abs(cycle_rho) if finite(cycle_rho) else np.nan,
                "high_correlation_edge": bool(finite(cycle_rho) and len(joined) >= CORRELATION_MIN_N and abs(cycle_rho) >= CORRELATION_THRESHOLD),
            })
    edge_df = pd.DataFrame(edges)
    # Connected components of the frozen high-correlation graph.
    parent = {v: v for v in variables}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        x, y = find(x), find(y)
        if x != y:
            parent[y] = x
    for row in edge_df.loc[edge_df["high_correlation_edge"]].itertuples():
        union(row.left_variable, row.right_variable)
    roots = {}
    clusters = {}
    for variable in variables:
        root = find(variable)
        roots.setdefault(root, len(roots) + 1)
        clusters[variable] = f"cluster_{roots[root]}"
    cluster_rows = [{"variable": v, "biomarker": str(test_rows[v]["biomarker"]), "exposure_axis": str(test_rows[v]["exposure_axes"]), "cluster_id": clusters[v]} for v in variables]
    cluster_df = pd.DataFrame(cluster_rows)
    metadata = {}
    for variable in variables:
        partner = edge_df.loc[(edge_df["left_variable"].eq(variable)) | (edge_df["right_variable"].eq(variable))].copy()
        metadata[variable] = {
            "cluster_id": clusters[variable],
            "cluster_n": int((cluster_df["cluster_id"] == clusters[variable]).sum()),
            "max_abs_cycle_adjusted_rho": float(partner["abs_cycle_adjusted_rho"].max()) if len(partner) else np.nan,
        }
    return edge_df, cluster_df, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", type=Path, default=DEFAULT_TESTS)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if not PLAN.exists():
        raise FileNotFoundError(PLAN)

    t2d = load_module(STEP5 / "run_t2d_screen.py", "step6_t2d_screen")
    model = load_module(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "step6_t2d_model")
    stability = load_module(ROOT / "work" / "scripts" / "mcop_crc_phase2_stability.py", "step6_t2d_stability")
    tests = pd.read_csv(args.tests, dtype=str, keep_default_na=False)
    registry = pd.read_csv(args.registry, low_memory=False)
    primary = pd.read_csv(STEP5 / "t2d_primary_29_tests.csv", low_memory=False)
    if len(tests) != FDR_DENOMINATOR or tests["test_id"].nunique() != FDR_DENOMINATOR:
        raise ValueError("Step 6 requires exactly the frozen 29-test set")
    if len(primary) != FDR_DENOMINATOR or primary["test_id"].nunique() != FDR_DENOMINATOR:
        raise ValueError("Step 6 requires exactly 29 Step 5 T2D results")
    supported = primary.loc[primary["FDR_supported"].astype(str).str.lower().eq("true")].copy()
    if supported.empty:
        raise ValueError("No FDR-positive T2D tests are available for Step 6")

    specs = {str(spec["cycle"]): (idx, spec) for idx, spec in enumerate(model.CYCLES)}
    needed_cycles = sorted({cycle for value in tests.loc[tests["test_id"].isin(supported["test_id"]), "cycles"] for cycle in str(value).split(";") if cycle})
    outcome_frames = {}
    outcome_qc = []
    for cycle in needed_cycles:
        idx, spec = specs[cycle]
        frame, qc = t2d.read_t2d_outcome_cycle(model, spec, idx, args.data_dir)
        outcome_frames[cycle] = frame
        outcome_qc.append(qc)

    test_rows = {str(row["test_id"]): row for _, row in tests.loc[tests["test_id"].isin(supported["test_id"])].iterrows()}
    populations = {}
    exposures = {}
    sources = {}
    for test_id, test_row in test_rows.items():
        exposure, source = t2d.read_test_exposure(t2d.load_module(FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py", "step6_exposure_reader"), test_row, registry)
        exposures[test_id] = exposure
        sources[test_id] = source
        outcome = pd.concat([outcome_frames[cycle] for cycle in source.get("cycles", [])], ignore_index=True)
        populations[test_id] = prepare_population(exposure, outcome, test_row)

    edge_df, cluster_df, cluster_meta = pairwise_correlations(exposures, test_rows)
    variables = sorted(exposures)
    corr_matrix = pd.DataFrame(np.eye(len(variables)), index=variables, columns=variables)
    for corr_row in edge_df.itertuples(index=False):
        corr_matrix.loc[corr_row.left_variable, corr_row.right_variable] = corr_row.cycle_adjusted_spearman_rho
        corr_matrix.loc[corr_row.right_variable, corr_row.left_variable] = corr_row.cycle_adjusted_spearman_rho
    cluster_df.to_csv(args.outdir / "t2d_exposure_clusters.csv", index=False)
    edge_df.to_csv(args.outdir / "t2d_exposure_pairwise_correlation.csv", index=False)
    corr_matrix.to_csv(args.outdir / "t2d_exposure_correlation_matrix.csv")

    all_loco, all_cycle, all_hetero, all_tail, all_creat, all_lod, all_age40, all_sex, all_sexint = [], [], [], [], [], [], [], [], []
    score_rows = []
    source_rows = []
    for _, screen_row in supported.sort_values("BH_FDR").iterrows():
        test_id = str(screen_row["test_id"])
        test_row = test_rows[test_id]
        population = populations[test_id]
        if population.empty:
            raise RuntimeError(f"Empty T2D population for {test_id}")
        primary_row, primary_fit = audit_primary(population, model, test_row, screen_row)
        pooled_or = float(primary_row["OR"])
        loco = run_loco(population, model, test_row, pooled_or)
        cycle = run_cycle_specific(population, model, test_row, pooled_or)
        hetero = run_cycle_interaction(population, model, test_row, stability)
        tail = run_tail(population, model, test_row, pooled_or)
        creat = run_creatinine(population, model, test_row, pooled_or)
        lod = run_lod(population, model, test_row, registry, pooled_or)
        age40 = run_age40(population, model, test_row, pooled_or)
        sex = run_sex(population, model, test_row, pooled_or)
        sexint = run_sex_interaction(population, model, test_row, stability)
        loco_summary = summarize_loco(loco)
        cycle_summary = summarize_cycle(cycle)
        tail_summary = summarize_tail(tail, pooled_or)
        statuses = [str(primary_row.get("fit_status")), *loco["fit_status"].astype(str).tolist(), *cycle["fit_status"].astype(str).tolist(), str(hetero.get("fit_status")), *tail["fit_status"].astype(str).tolist(), str(creat.get("fit_status")), str(lod.get("fit_status")), str(age40.get("fit_status")), *sex["fit_status"].astype(str).tolist(), str(sexint.get("fit_status"))]
        score = score_candidate(primary_row, loco_summary, cycle_summary, hetero, tail_summary, statuses, screen_row, cluster_meta[test_id])
        score_rows.append(score)
        all_loco.append(loco)
        all_cycle.append(cycle)
        all_hetero.append(hetero)
        all_tail.append(tail)
        all_creat.append(creat)
        all_lod.append(lod)
        all_age40.append(age40)
        all_sex.append(sex)
        all_sexint.append(sexint)
        source_rows.append({"test_id": test_id, "source_status": sources[test_id].get("status"), "source_cycles": ";".join(sources[test_id].get("cycles", [])), "source_n": sources[test_id].get("n_raw", 0), "population_n": len(population), "population_t2d_cases": int(population["outcome"].sum()), "population_controls": int(len(population) - population["outcome"].sum()), "source_rows": len(sources[test_id].get("source_rows", []))})
        if not finite(primary_row.get("absolute_log_or_difference")) or primary_row["absolute_log_or_difference"] > 1e-8:
            raise AssertionError(f"Primary T2D reproduction failed for {test_id}: {primary_row.get('absolute_log_or_difference')}")

    score_df = pd.DataFrame(score_rows).sort_values(["priority_tier", "BH_FDR", "P"]).reset_index(drop=True)
    score_df.to_csv(args.outdir / "t2d_robustness_results.csv", index=False)
    pd.concat(all_loco, ignore_index=True).to_csv(args.outdir / "t2d_loco_results.csv", index=False)
    pd.concat(all_cycle, ignore_index=True).to_csv(args.outdir / "t2d_cycle_specific_results.csv", index=False)
    pd.DataFrame(all_hetero).to_csv(args.outdir / "t2d_cycle_heterogeneity.csv", index=False)
    pd.concat(all_tail, ignore_index=True).to_csv(args.outdir / "t2d_tail_results.csv", index=False)
    pd.DataFrame(all_creat).to_csv(args.outdir / "t2d_creatinine_results.csv", index=False)
    pd.DataFrame(all_lod).to_csv(args.outdir / "t2d_lod_results.csv", index=False)
    pd.DataFrame(all_age40).to_csv(args.outdir / "t2d_age40_results.csv", index=False)
    pd.concat(all_sex, ignore_index=True).to_csv(args.outdir / "t2d_sex_results.csv", index=False)
    pd.DataFrame(all_sexint).to_csv(args.outdir / "t2d_sex_interaction.csv", index=False)
    pd.DataFrame(outcome_qc).to_csv(args.outdir / "t2d_step6_outcome_qc.csv", index=False)
    pd.DataFrame(source_rows).to_csv(args.outdir / "t2d_step6_source_population_qc.csv", index=False)

    robust_n = int(score_df["priority_tier"].eq("robust_fdr_candidate").sum())
    robust_names = score_df.loc[score_df["priority_tier"].eq("robust_fdr_candidate"), "biomarker"].astype(str).tolist()
    unstable_names = score_df.loc[score_df["priority_tier"].ne("robust_fdr_candidate"), "biomarker"].astype(str).tolist()
    heterogeneity_concerns = score_df.loc[score_df["H"].isin([0, 1]), ["biomarker", "H", "H_note"]].to_dict(orient="records")
    high_edges = edge_df.loc[edge_df["high_correlation_edge"]].copy()
    heterogeneity_text = "; ".join(f"{item['biomarker']} {item['H_note']}" for item in heterogeneity_concerns) if heterogeneity_concerns else "none"
    high_edges_text = "; ".join(f"{row.left_biomarker}-{row.right_biomarker} |rho|={abs(float(row.cycle_adjusted_spearman_rho)):.3f}" for row in high_edges.itertuples()) if len(high_edges) else "none"
    report = [
        "# Step 6 T2D robustness and exposure-cluster audit",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Scope and firewall",
        "",
        f"- Frozen primary family: **{FDR_DENOMINATOR} tests**; no FDR was recomputed or narrowed.",
        f"- FDR-positive tests audited uniformly: **{len(score_df)}**.",
        f"- Primary reproductions within absolute log-OR difference <=1e-8: **{len(score_df)}/{len(score_df)}**.",
        f"- Exposure clusters under |cycle-adjusted Spearman rho| >= {CORRELATION_THRESHOLD} and pairwise N >= {CORRELATION_MIN_N}: **{cluster_df['cluster_id'].nunique()}**.",
        f"- Deterministic robust-FDR candidates under the locked rubric: **{robust_n}**.",
        f"- Robust-FDR candidate list: **{', '.join(robust_names)}**.",
        f"- FDR-positive test(s) with a stability downgrade: **{', '.join(unstable_names) if unstable_names else 'none'}**.",
        f"- Heterogeneity concerns (H0/H1): **{heterogeneity_text}**.",
        f"- High-correlation edges: **{high_edges_text}**.",
        "- No GeneCards, disease-specific CTD, transcriptomics, literature, or mechanistic analysis was performed.",
        "",
        "## Candidate audit",
        "",
        "| Biomarker | OR | 95% CI | P | q (29) | LOCO same | Cycle same | Pinteraction | Tail max | Cluster | Priority |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for _, row in score_df.sort_values("BH_FDR").iterrows():
        fmt = lambda x: "NA" if pd.isna(x) else f"{float(x):.4g}"
        report.append(f"| {row['biomarker']} | {fmt(row.get('OR'))} | {fmt(row.get('CI_low'))}-{fmt(row.get('CI_high'))} | {fmt(row.get('P'))} | {fmt(row.get('BH_FDR'))} | {int(row.get('loco_same_n', 0))}/{int(row.get('loco_n', 0))} | {int(row.get('cycle_same_n', 0))}/{int(row.get('cycle_n', 0))} | {row.get('H_note')} | {fmt(row.get('tail_max_abs_log_or_delta'))} | {row.get('cluster_id')} | {row.get('priority_tier')} |")
    report += [
        "",
        "## Interpretation",
        "",
        "The audit distinguishes primary FDR support from stability and exposure dependence. A high-correlation cluster does not invalidate a primary test and was not used to alter multiplicity. Conversely, a robust-FDR label is a prioritization aid and does not establish temporality or causality.",
    ]
    (args.outdir / "t2d_robustness_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    outputs = [
        "t2d_robustness_results.csv", "t2d_loco_results.csv", "t2d_cycle_specific_results.csv", "t2d_cycle_heterogeneity.csv",
        "t2d_tail_results.csv", "t2d_creatinine_results.csv", "t2d_lod_results.csv", "t2d_age40_results.csv",
        "t2d_sex_results.csv", "t2d_sex_interaction.csv", "t2d_exposure_pairwise_correlation.csv", "t2d_exposure_correlation_matrix.csv",
        "t2d_exposure_clusters.csv", "t2d_step6_outcome_qc.csv", "t2d_step6_source_population_qc.csv", "t2d_robustness_report.md",
    ]
    manifest = {
        "lock_type": "T2D_STEP6_ROBUSTNESS_AND_EXPOSURE_CLUSTER_AUDIT",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "scope_test_count": int(len(score_df)),
        "scope_test_ids": sorted(score_df["test_id"].astype(str).tolist()),
        "fdr_denominator": FDR_DENOMINATOR,
        "step5_primary": str(STEP5 / "t2d_primary_29_tests.csv"),
        "step5_primary_sha256": sha256(STEP5 / "t2d_primary_29_tests.csv"),
        "step5_script": str(STEP5 / "run_t2d_screen.py"),
        "step5_script_sha256": sha256(STEP5 / "run_t2d_screen.py"),
        "plan": str(PLAN),
        "plan_sha256": sha256(PLAN),
        "estimator": "Python NHANES survey-weighted Taylor/PSU-sandwich logistic estimator reused from frozen T2D screen",
        "correlation_rule": {"method": "pairwise Spearman on log2 exposure; cycle-adjusted sensitivity demeans within cycle", "absolute_rho_threshold": CORRELATION_THRESHOLD, "minimum_pairwise_n": CORRELATION_MIN_N, "cluster_rule": "connected components of high-correlation graph"},
        "primary_reproduction_max_absolute_log_or_difference": float(score_df["absolute_log_or_difference"].max()),
        "robust_fdr_candidate_count": robust_n,
        "cluster_count": int(cluster_df["cluster_id"].nunique()),
        "outputs": {name: {"path": str(args.outdir / name), "sha256": sha256(args.outdir / name)} for name in outputs},
    }
    (args.outdir / "t2d_step6_analysis_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"audited_tests": len(score_df), "robust_fdr_candidates": robust_n, "clusters": int(cluster_df["cluster_id"].nunique()), "max_primary_reproduction_abs_log_or": float(score_df["absolute_log_or_difference"].max())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
