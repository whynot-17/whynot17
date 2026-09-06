"""Phase 2H: seven-cycle NHANES MCOP complex-survey audit.

This audit is deliberately separate from the historical MBzP analysis.  It
rebuilds the MCOP frame for the seven cycles with URXCOP measurements
(2005-06 through 2017-18), applies the CDC pooled-weight rule ``2-year / 7``,
re-runs the frozen Python Taylor-style fit, and independently reproduces the
primary model with R ``survey::svyglm``.

The script also reruns the frozen, targeted sensitivity set without changing
the outcome definition or covariate processing: LOCO, age >=40, sex and
interaction, diagnosis-timing exclusions, tail exclusion, creatinine-ratio
normalization, and pairwise co-exposure models.  Weighted quartile cutpoints
and weighted restricted-cubic-spline knots are reported as sensitivities; the
original unweighted cutpoints are retained for direct comparability.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_HARMONIZED = ROOT / "work" / "nhanes_phase2a" / "phase2b_harmonized.pkl"
DEFAULT_OUTPUT = ROOT / "outputs"
DEFAULT_RSCRIPT = Path("D:/CodexData/R-4.5.1/bin/Rscript.exe")
R_SCRIPT = ROOT / "work" / "scripts" / "mcop_crc_phase2h_survey_standard.R"

MCOP_CYCLES = [
    "2005-2006",
    "2007-2008",
    "2009-2010",
    "2011-2012",
    "2013-2014",
    "2015-2016",
    "2017-2018",
]
N_CYCLES = len(MCOP_CYCLES)
EXPECTED_WEIGHT_SOURCES = {
    cycle: ("WTSA2YR" if cycle == "2011-2012" else "WTSB2YR")
    for cycle in MCOP_CYCLES
}
PRIMARY_CONTINUOUS = ["mcop_log2", "age", "bmi", "pir", "creatinine_log2"]
PRIMARY_CATEGORICAL = ["sex", "race", "smoking"]
PRIMARY_REQUIRED = [
    "outcome",
    *PRIMARY_CONTINUOUS,
    *PRIMARY_CATEGORICAL,
    "pooled_weight",
    "psu",
    "strata",
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_components():
    phase2 = load_module(ROOT / "work" / "scripts" / "mcop_crc_phase2.py", "phase2h_mcop")
    model = phase2.load_validated_functions()
    paper = load_module(ROOT / "work" / "scripts" / "mcop_crc_phase2d_paper_audit.py", "phase2h_paper")
    stability = load_module(ROOT / "work" / "scripts" / "mcop_crc_phase2_stability.py", "phase2h_stability")
    return phase2, model, paper, stability


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def public_fit(fit: dict, analysis: str, cycles: list[str] | None = None) -> dict:
    return {
        "Analysis": analysis,
        "Cycles": ";".join(cycles or MCOP_CYCLES),
        **{key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}},
    }


def build_seven_cycle_frame(phase2, harmonized_path: Path, data_dir: Path) -> tuple[pd.DataFrame, list[dict]]:
    harmonized = pd.read_pickle(harmonized_path)
    mcop, input_manifest = phase2.read_mcop(data_dir)
    data = phase2.build_frame(harmonized, mcop)
    data = data[data["cycle"].isin(MCOP_CYCLES)].copy()
    if sorted(data["cycle"].dropna().unique().tolist()) != MCOP_CYCLES:
        raise RuntimeError(f"MCOP frame does not contain exactly the seven requested cycles: {sorted(data['cycle'].dropna().unique())}")
    data["pooled_weight"] = pd.to_numeric(data["phthalate_weight_base"], errors="coerce") / N_CYCLES
    data["pooled_weight_rule"] = f"cycle-specific phthalate subsample weight / {N_CYCLES}"
    data["weight_rule_audit"] = np.where(
        data["weight_source"].eq(data["cycle"].map(EXPECTED_WEIGHT_SOURCES)),
        "pass",
        "FAIL",
    )
    bad = data.loc[data["weight_rule_audit"].eq("FAIL"), ["cycle", "weight_source"]].drop_duplicates()
    if not bad.empty:
        raise RuntimeError(f"Unexpected phthalate weight source(s): {bad.to_dict('records')}")
    return data, input_manifest


def weight_source_audit(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cycle in MCOP_CYCLES:
        part = data[data["cycle"].eq(cycle)].copy()
        weights = pd.to_numeric(part["phthalate_weight_base"], errors="coerce")
        pooled = pd.to_numeric(part["pooled_weight"], errors="coerce")
        near_zero_base = weights.abs().lt(1e-12)
        near_zero_pooled = pooled.abs().lt(1e-12)
        rows.append(
            {
                "cycle": cycle,
                "weight_source": ";".join(sorted(part["weight_source"].dropna().astype(str).unique())),
                "expected_weight_source": EXPECTED_WEIGHT_SOURCES[cycle],
                "N": int(len(part)),
                "nonzero_base_weight_N": int(weights.gt(0).sum()),
                "nonzero_pooled_weight_N": int(pooled.gt(0).sum()),
                "near_zero_sas_missing_sentinel_N": int(near_zero_base.sum()),
                "usable_base_weight_N_gt_1e-12": int(weights.gt(1e-12).sum()),
                "usable_pooled_weight_N_gt_1e-12": int(pooled.gt(1e-12).sum()),
                "base_weight_sum": float(weights.sum()),
                "pooled_weight_sum": float(pooled.sum()),
                "base_weight_min": float(weights[weights.gt(0)].min()) if weights.gt(0).any() else np.nan,
                "base_weight_max": float(weights.max()) if weights.notna().any() else np.nan,
                "pooled_weight_min": float(pooled[pooled.gt(0)].min()) if pooled.gt(0).any() else np.nan,
                "pooled_weight_max": float(pooled.max()) if pooled.notna().any() else np.nan,
                "weight_rule": f"{EXPECTED_WEIGHT_SOURCES[cycle]} / 7",
            }
        )
    return pd.DataFrame(rows)


def complete_case_primary(primary: pd.DataFrame) -> pd.DataFrame:
    work = primary.dropna(subset=PRIMARY_REQUIRED).copy()
    return work[work["pooled_weight"].gt(0)].copy()


def design_unit_audit(primary: pd.DataFrame) -> pd.DataFrame:
    cc = complete_case_primary(primary)
    rows = []
    for (cycle, strata), group in cc.groupby(["cycle", "strata"], dropna=False, sort=True):
        n_psu = int(group["psu"].nunique())
        rows.append(
            {
                "scope": "cycle_stratum",
                "cycle": cycle,
                "strata": strata,
                "n_participants": int(len(group)),
                "crc_cases": int(group["outcome"].sum()),
                "psu_N": n_psu,
                "singleton_stratum": bool(n_psu < 2),
            }
        )
    for label, group in [("ALL_7_CYCLES", cc)]:
        rows.append(
            {
                "scope": "pooled",
                "cycle": label,
                "strata": np.nan,
                "n_participants": int(len(group)),
                "crc_cases": int(group["outcome"].sum()),
                "psu_N": int(group["psu"].nunique()),
                "strata_N": int(group["strata"].nunique()),
                "singleton_strata_N": int(sum(group.groupby("strata")["psu"].nunique().lt(2))),
                "singleton_stratum": np.nan,
            }
        )
    result = pd.DataFrame(rows)
    if "strata_N" not in result.columns:
        result["strata_N"] = np.nan
    if "singleton_strata_N" not in result.columns:
        result["singleton_strata_N"] = np.nan
    return result


def fit_primary(primary: pd.DataFrame, model, label: str, categorical: list[str] | None = None) -> dict:
    fit = model.fit_survey_logistic(
        primary,
        PRIMARY_CONTINUOUS,
        categorical or PRIMARY_CATEGORICAL,
        exposure_name="mcop_log2",
        levels=model.LEVELS,
    )
    return public_fit(fit, label)


def rerun_python(primary: pd.DataFrame, data: pd.DataFrame, model) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = [fit_primary(primary, model, "Primary_7_cycle_weight_div_7")]
    legacy = data.copy()
    legacy["pooled_weight"] = pd.to_numeric(legacy["phthalate_weight_base"], errors="coerce") / 10.0
    legacy_primary = model.population_frames(legacy)["CRC_vs_cancer_free"]
    rows.append(fit_primary(legacy_primary, model, "Legacy_weight_div_10_same_7_cycles"))
    rows.append(fit_primary(primary[primary["age"] >= 40].copy(), model, "Age_ge_40_7_cycle"))

    loco = []
    for cycle in MCOP_CYCLES:
        row = fit_primary(primary[primary["cycle"].ne(cycle)].copy(), model, f"LOCO_drop_{cycle}")
        row["Dropped_cycle"] = cycle
        loco.append(row)
    return pd.DataFrame(rows), pd.DataFrame(loco)


def weighted_quantile(values: pd.Series, weights: pd.Series, quantiles: list[float]) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    w = pd.to_numeric(weights, errors="coerce").to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[valid], w[valid]
    if len(x) == 0:
        return np.full(len(quantiles), np.nan)
    order = np.argsort(x)
    x, w = x[order], w[order]
    positions = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(quantiles, positions, x)


def fit_quartiles(primary: pd.DataFrame, model, cutpoints: np.ndarray, method: str) -> pd.DataFrame:
    cutpoints = np.unique(cutpoints[np.isfinite(cutpoints)])
    if len(cutpoints) < 5:
        return pd.DataFrame([{"method": method, "status": "not_estimable", "reason": "fewer than four distinct cutpoints"}])
    work = primary.copy()
    work["mcop_quartile"] = pd.cut(
        work["mcop_log2"],
        bins=[-np.inf, *cutpoints[1:-1], np.inf],
        labels=False,
        include_lowest=True,
    ) + 1
    cc = work.dropna(subset=PRIMARY_REQUIRED + ["mcop_quartile"]).copy()
    cc = cc[cc["pooled_weight"].gt(0)]
    trend_frame = cc.assign(qtrend=cc["mcop_quartile"].astype(float))
    trend = model.fit_survey_logistic(
        trend_frame,
        ["qtrend", "age", "bmi", "pir", "creatinine_log2"],
        PRIMARY_CATEGORICAL,
        exposure_name="qtrend",
        levels=model.LEVELS,
    )
    rows = []
    for q in [1, 2, 3, 4]:
        if q == 1:
            row = {
                "method": method,
                "Quartile": "Q1",
                "Reference": True,
                "N": len(cc),
                "CRC_N": int(cc["outcome"].sum()),
                "OR": 1.0,
                "CI_low": 1.0,
                "CI_high": 1.0,
                "P": np.nan,
                "status": "reference",
            }
        else:
            contrast = cc.copy()
            contrast["q_exposure"] = contrast["mcop_quartile"].eq(q).astype(float)
            fit = model.fit_survey_logistic(
                contrast,
                ["q_exposure", "age", "bmi", "pir", "creatinine_log2"],
                PRIMARY_CATEGORICAL,
                exposure_name="q_exposure",
                levels=model.LEVELS,
            )
            row = {"method": method, "Quartile": f"Q{q}", "Reference": False, **{k: v for k, v in fit.items() if k not in {"coefficients", "covariance"}}}
        row["P_trend"] = trend.get("P", np.nan)
        row["cutpoints_log2"] = ";".join(f"{value:.8g}" for value in cutpoints)
        rows.append(row)
    return pd.DataFrame(rows)


def rcs_sensitivity(primary: pd.DataFrame, model, stability, cutpoints: np.ndarray, method: str) -> dict:
    cc = complete_case_primary(primary)
    knots = np.unique(cutpoints[np.isfinite(cutpoints)])
    if len(knots) < 4:
        return {"method": method, "status": "not_estimable", "reason": "insufficient distinct knots"}
    basis = model.spline_basis(cc["mcop_log2"].to_numpy(float), knots)
    basis_cols = []
    for index in range(basis.shape[1]):
        name = f"mcop_spline_{index + 1}"
        cc[name] = basis[:, index]
        basis_cols.append(name)
    fit = model.fit_survey_logistic(
        cc,
        basis_cols + ["age", "bmi", "pir", "creatinine_log2"],
        PRIMARY_CATEGORICAL,
        exposure_name=basis_cols[0],
        levels=model.LEVELS,
    )
    overall = stability.wald_test(fit, basis_cols)
    nonlinear = stability.wald_test(fit, basis_cols[1:])
    return {
        "method": method,
        "status": fit.get("status"),
        "N": fit.get("N"),
        "CRC_N": fit.get("CRC_N"),
        "knots_log2": ";".join(f"{value:.8g}" for value in knots),
        "overall_P_chi2": overall.get("P_chi2"),
        "overall_P_F": overall.get("P_F"),
        "nonlinear_P_chi2": nonlinear.get("P_chi2"),
        "nonlinear_P_F": nonlinear.get("P_F"),
        "design_df": fit.get("design_df"),
    }


def weighted_quantile_sensitivity(primary: pd.DataFrame, model, stability) -> pd.DataFrame:
    exposure = primary[primary["mcop_log2"].notna() & primary["pooled_weight"].gt(0)].copy()
    unweighted_q = exposure["mcop_log2"].quantile([0, .25, .5, .75, 1]).to_numpy()
    weighted_q = weighted_quantile(exposure["mcop_log2"], exposure["pooled_weight"], [0, .25, .5, .75, 1])
    quartiles = pd.concat(
        [fit_quartiles(primary, model, unweighted_q, "unweighted_cutpoints"), fit_quartiles(primary, model, weighted_q, "survey_weighted_cutpoints")],
        ignore_index=True,
    )
    cc = complete_case_primary(primary)
    unweighted_knots = cc["mcop_log2"].quantile([.05, .35, .65, .95]).to_numpy()
    weighted_knots = weighted_quantile(cc["mcop_log2"], cc["pooled_weight"], [.05, .35, .65, .95])
    rcs = pd.DataFrame(
        [
            rcs_sensitivity(primary, model, stability, unweighted_knots, "RCS_unweighted_knots"),
            rcs_sensitivity(primary, model, stability, weighted_knots, "RCS_survey_weighted_knots"),
        ]
    )
    quartiles["analysis_family"] = "quartiles"
    rcs["analysis_family"] = "RCS"
    return pd.concat([quartiles, rcs], ignore_index=True, sort=False)


def run_standard_survey(primary: pd.DataFrame, model, outdir: Path, rscript: Path) -> tuple[pd.DataFrame, dict]:
    cc = complete_case_primary(primary)
    input_path = outdir.parent / "nhanes_phase2a" / "data" / ".mcop_phase2h_r_input.csv"
    input_path.parent.mkdir(parents=True, exist_ok=True)
    r_output = outdir / ".mcop_phase2h_r_output.csv"
    columns = ["outcome", "mcop_log2", "age", "bmi", "pir", "creatinine_log2", "sex", "race", "smoking", "pooled_weight", "psu", "strata"]
    cc[columns].to_csv(input_path, index=False)
    result = {"r_status": "not_run"}
    try:
        command = [str(rscript), str(R_SCRIPT), str(input_path), str(r_output)]
        completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
        result["r_exit_code"] = int(completed.returncode)
        result["r_stdout"] = completed.stdout[-4000:]
        result["r_stderr"] = completed.stderr[-4000:]
        if completed.returncode != 0 or not r_output.exists():
            result["r_status"] = "failed"
            return pd.DataFrame([result]), result
        standard = pd.read_csv(r_output).iloc[0].to_dict()
        py = fit_primary(primary, model, "Python_Taylor_primary")
        result = {
            "Analysis": "Python_vs_R_survey_primary",
            "python_status": py.get("status"),
            "python_N": py.get("N"),
            "python_CRC_N": py.get("CRC_N"),
            "python_beta": py.get("beta"),
            "python_SE": py.get("SE"),
            "python_OR": py.get("OR"),
            "python_CI_low": py.get("CI_low"),
            "python_CI_high": py.get("CI_high"),
            "python_P": py.get("P"),
            "python_design_df": py.get("design_df"),
            "python_PSU_N": py.get("PSU_N"),
            "python_strata_N": py.get("strata_N"),
            "r_status": standard.get("status"),
            "r_N": standard.get("N"),
            "r_CRC_N": standard.get("CRC_N"),
            "r_beta": standard.get("beta"),
            "r_SE": standard.get("SE"),
            "r_OR": standard.get("OR"),
            "r_CI_low": standard.get("CI_low"),
            "r_CI_high": standard.get("CI_high"),
            "r_P_standard": standard.get("P_standard"),
            "r_P_design_df": standard.get("P_design_df"),
            "r_design_df": standard.get("design_df"),
            "r_model_residual_df": standard.get("model_residual_df"),
            "r_PSU_N": standard.get("PSU_N"),
            "r_strata_N": standard.get("strata_N"),
            "r_singleton_strata_N": standard.get("singleton_strata_N"),
            "r_lonely_psu_option": standard.get("survey_lonely_psu_option"),
            "r_survey_package_version": standard.get("survey_package_version"),
        }
        result["absolute_logOR_difference"] = abs(float(result["r_beta"]) - float(result["python_beta"]))
        result["relative_logOR_change_pct"] = result["absolute_logOR_difference"] / abs(float(result["python_beta"])) * 100
        result["OR_direction_same"] = bool((float(result["r_beta"]) > 0) == (float(result["python_beta"]) > 0))
        result["CI_null_conclusion_same"] = bool((float(result["r_CI_low"]) < 1) == (float(result["python_CI_low"]) < 1) and (float(result["r_CI_high"]) > 1) == (float(result["python_CI_high"]) > 1))
        return pd.DataFrame([result]), result
    finally:
        for path in [input_path, r_output]:
            try:
                path.unlink()
            except FileNotFoundError:
                pass


def concat_sensitivities(primary: pd.DataFrame, model, paper, stability, outdir: Path) -> pd.DataFrame:
    wrapper = type("Module", (), {"model_module": model})()
    rows: list[pd.DataFrame] = []
    _, loco = rerun_python(primary, primary, model)
    loco["analysis_family"] = "LOCO"
    rows.append(loco)
    age40 = pd.DataFrame([fit_primary(primary[primary["age"] >= 40].copy(), model, "Age_ge_40_7_cycle")])
    age40["analysis_family"] = "age_restriction"
    rows.append(age40)
    sex, interaction = paper.run_sex_audit(primary, model)
    sex["analysis_family"] = "sex_specific"
    interaction["analysis_family"] = "sex_interaction"
    rows.extend([sex, interaction])
    timing, timing_cases = paper.run_timing_audit(primary, model)
    timing["analysis_family"] = "diagnosis_timing"
    rows.append(timing)
    timing_cases.to_csv(outdir / "mcop_crc_phase2h_diagnosis_timing_cases.csv", index=False)
    tails = stability.tail_exclusion(primary, wrapper)
    tails["analysis_family"] = "tail_exclusion"
    rows.append(tails)
    normalized = stability.creatinine_normalized(primary, wrapper)
    normalized["analysis_family"] = "creatinine_normalization"
    rows.append(normalized)
    coexposure = paper.run_coexposure_audit(primary, model)
    coexposure["analysis_family"] = "coexposure"
    rows.append(coexposure)
    return pd.concat(rows, ignore_index=True, sort=False)


def write_report(
    outdir: Path,
    data: pd.DataFrame,
    weights: pd.DataFrame,
    design_units: pd.DataFrame,
    primary_rows: pd.DataFrame,
    loco: pd.DataFrame,
    standard: pd.DataFrame,
    weighted_sensitivity: pd.DataFrame,
    input_manifest: list[dict],
    harmonized: Path,
    rscript: Path,
) -> dict:
    py = primary_rows.loc[primary_rows["Analysis"].eq("Primary_7_cycle_weight_div_7")].iloc[0].to_dict()
    legacy = primary_rows.loc[primary_rows["Analysis"].eq("Legacy_weight_div_10_same_7_cycles")].iloc[0].to_dict()
    standard_row = standard.iloc[0].to_dict() if not standard.empty else {}
    r_ok = standard_row.get("r_status") == "ok"
    rel_change = float(standard_row.get("relative_logOR_change_pct", np.nan)) if r_ok else np.nan
    same_direction = bool(standard_row.get("OR_direction_same", False)) if r_ok else False
    same_ci = bool(standard_row.get("CI_null_conclusion_same", False)) if r_ok else False
    if r_ok and float(standard_row.get("r_OR", np.nan)) > 1 and rel_change <= 10 and same_ci:
        decision = "GREEN"
    elif r_ok and float(standard_row.get("r_OR", np.nan)) > 1:
        decision = "YELLOW"
    else:
        decision = "RED"
    singleton = int(design_units.loc[design_units["scope"].eq("pooled"), "singleton_strata_N"].iloc[0])
    manifest = {
        "analysis": "MCOP-CRC Phase 2H seven-cycle complex-survey audit",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "cycles": MCOP_CYCLES,
        "n_cycles": N_CYCLES,
        "weight_rule": "cycle-specific phthalate subsample weight divided by 7",
        "expected_weight_sources": EXPECTED_WEIGHT_SOURCES,
        "outcome": "CRC type code 16 or 31 versus MCQ220=2 cancer-free controls",
        "primary_exposure": "log2(URXCOP)",
        "covariates": ["age", "sex", "race", "BMI", "smoking", "PIR", "log2(URXUCR)"],
        "python_method": "existing validated Newton-IRLS weighted logistic estimating equations with stratified PSU sandwich",
        "standard_method": "R survey::svyglm quasibinomial Taylor linearization",
        "rscript": str(rscript),
        "r_script_sha256": sha256(R_SCRIPT),
        "input_harmonized": {"path": str(harmonized), "sha256": sha256(harmonized)},
        "input_files": input_manifest,
        "decision": decision,
        "singleton_strata_primary_complete_case": singleton,
        "standard_survey_ran": r_ok,
    }
    (outdir / "mcop_crc_phase2h_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    lines = [
        "# MCOP–CRC Phase 2H：复杂抽样设计与标准 survey 软件复核",
        "",
        "## Scope",
        "",
        "本轮只重建 MCOP 覆盖的七个 NHANES 2-year cycle（2005–06 至 2017–18）。MBzP 的十-cycle 历史分析没有修改。病例定义、cancer-free controls、协变量处理和 MCOP 主暴露均沿用 Phase 2 冻结版本。",
        "",
        f"- Pooled weight rule: **cycle-specific phthalate subsample weight / {N_CYCLES}**.",
        f"- Complete-case primary N={int(py.get('N', 0))}, CRC cases={int(py.get('CRC_N', 0))}.",
        f"- Weight-source audit: **{'PASS' if (weights['weight_source'] == weights['expected_weight_source']).all() else 'FAIL'}**.",
        f"- Near-zero SAS missing-weight sentinels are recorded separately; primary complete-case modeling uses no such row.",
        f"- Primary complete-case singleton strata: **{singleton}**.",
        "",
        "## Python `/10` versus `/7` check",
        "",
        f"- `/7`: OR={float(py.get('OR', np.nan)):.8g}, 95% CI {float(py.get('CI_low', np.nan)):.8g}–{float(py.get('CI_high', np.nan)):.8g}, P={float(py.get('P', np.nan)):.8g}.",
        f"- Legacy `/10` on the same seven-cycle participants: OR={float(legacy.get('OR', np.nan)):.8g}, 95% CI {float(legacy.get('CI_low', np.nan)):.8g}–{float(legacy.get('CI_high', np.nan)):.8g}, P={float(legacy.get('P', np.nan)):.8g}.",
        "- Because the Python fitter normalizes all weights by their sample mean before fitting, `/10` versus `/7` is a common multiplicative rescaling and should not change beta, SE, OR or P.",
        "",
        "## Independent R `survey::svyglm` gate",
        "",
    ]
    if r_ok:
        lines += [
            f"- R `svyglm`: OR={float(standard_row['r_OR']):.8g}, 95% CI {float(standard_row['r_CI_low']):.8g}–{float(standard_row['r_CI_high']):.8g}, standard P={float(standard_row['r_P_standard']):.8g}; design-df P={float(standard_row['r_P_design_df']):.8g}.",
            f"- Python: beta={float(standard_row['python_beta']):.8g}, SE={float(standard_row['python_SE']):.8g}; R: beta={float(standard_row['r_beta']):.8g}, SE={float(standard_row['r_SE']):.8g}.",
            f"- Relative logOR change: **{rel_change:.4g}%**; direction same: **{same_direction}**; CI null conclusion same: **{same_ci}**.",
            f"- R design df={int(standard_row['r_design_df'])}; R model residual df={int(standard_row['r_model_residual_df'])}; Python design df={int(standard_row['python_design_df'])}.",
        ]
    else:
        lines += ["- R standard survey run: **FAILED**. See `mcop_crc_phase2h_python_vs_standard_survey.csv` for the exit code and stderr."]
    lines += [
        "",
        "## Frozen sensitivity audit",
        "",
        "The sensitivity CSV contains LOCO, age ≥40, sex-specific effects and formal interaction, diagnosis-age exclusions (<1/<2/<5 years), top-tail exclusions, creatinine normalization, and pairwise co-exposure models. No mechanistic analysis was run in Phase 2H.",
        "",
        "Weighted quartiles use the phthalate subsample weights only for cutpoint construction; the unweighted cutpoint analysis is retained unchanged. The RCS sensitivity reports both unweighted and survey-weighted 5th/35th/65th/95th-percentile knot sets.",
        "",
        "## Decision",
        "",
        f"Phase 2H decision: **{decision}**. The gate is GREEN only when R `svyglm` remains positive, the relative logOR change is ≤10%, and the CI conclusion agrees with the Python implementation. This remains a cross-sectional association audit, not causal evidence.",
        "",
        "## CDC design references",
        "",
        "- [NHANES weighting tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx)",
        "- [NHANES sample design and analysis](https://wwwn.cdc.gov/nchs/nhanes/tutorials/sampledesign.aspx)",
        "",
        "## Outputs",
        "",
        "- `mcop_crc_phase2h_weight_sources.csv`",
        "- `mcop_crc_phase2h_design_units.csv`",
        "- `mcop_crc_phase2h_python_vs_standard_survey.csv`",
        "- `mcop_crc_phase2h_primary_reanalysis.csv`",
        "- `mcop_crc_phase2h_sensitivity_reanalysis.csv`",
        "- `mcop_crc_phase2h_weighted_quantile_sensitivity.csv`",
        "- `mcop_crc_phase2h_survey_design_audit.md`",
    ]
    (outdir / "mcop_crc_phase2h_survey_design_audit.md").write_text("\n".join(lines), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", type=Path, default=DEFAULT_HARMONIZED)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rscript", type=Path, default=DEFAULT_RSCRIPT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    phase2, model, paper, stability = load_components()
    data, input_manifest = build_seven_cycle_frame(phase2, args.harmonized, args.data_dir)
    weights = weight_source_audit(data)
    primary = model.population_frames(data)["CRC_vs_cancer_free"].copy()
    design_units = design_unit_audit(primary)
    primary_rows, loco = rerun_python(primary, data, model)
    standard, standard_meta = run_standard_survey(primary, model, args.outdir, args.rscript)
    weighted_sensitivity = weighted_quantile_sensitivity(primary, model, stability)
    sensitivity = concat_sensitivities(primary, model, paper, stability, args.outdir)

    weights.to_csv(args.outdir / "mcop_crc_phase2h_weight_sources.csv", index=False)
    design_units.to_csv(args.outdir / "mcop_crc_phase2h_design_units.csv", index=False)
    standard.to_csv(args.outdir / "mcop_crc_phase2h_python_vs_standard_survey.csv", index=False)
    primary_rows.to_csv(args.outdir / "mcop_crc_phase2h_primary_reanalysis.csv", index=False)
    sensitivity.to_csv(args.outdir / "mcop_crc_phase2h_sensitivity_reanalysis.csv", index=False)
    weighted_sensitivity.to_csv(args.outdir / "mcop_crc_phase2h_weighted_quantile_sensitivity.csv", index=False)
    write_report(
        args.outdir,
        data,
        weights,
        design_units,
        primary_rows,
        loco,
        standard,
        weighted_sensitivity,
        input_manifest,
        args.harmonized,
        args.rscript,
    )
    print(json.dumps({"decision": json.loads((args.outdir / "mcop_crc_phase2h_manifest.json").read_text(encoding="utf-8"))["decision"], "primary": primary_rows.iloc[0].to_dict(), "standard": standard_meta}, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
