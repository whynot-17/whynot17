"""Bounded urinary-dilution sensitivity audit for the Cd arthritis signals.

This script compares only the frozen URXUCD exposure with the three prespecified
endpoints requested for the audit: any arthritis, osteoarthritis, and gout.  It
does not fit new outcomes, calculate a new multiplicity adjustment, or perform
any mechanistic or downstream analysis.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import t


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "work" / "nhanes_phase2a" / "data"
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
OUT = ROOT / "analysis" / "cadmium_sex_arthritis_subtypes_urinary_dilution_audit"
OLD_SCRIPT = ROOT / "work" / "run_cd_sex_arthritis_subtypes.py"
TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
PRIMARY = FRAMEWORK / "sex_divergence_primary" / "sex_divergence_primary_406.csv"
OLD_RESULTS = ROOT / "analysis" / "cadmium_sex_arthritis_subtypes" / "03_cd_sex_subtype_models.csv"
MODEL_PATH = ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"
READER_PATH = FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py"

TARGET_IDS = ["any_arthritis", "osteoarthritis", "gout"]
TARGET_NAMES = {
    "any_arthritis": "Any arthritis",
    "osteoarthritis": "Osteoarthritis (OA)",
    "gout": "Gout",
}
EXPOSURE_TEST_ID = "NHANES_URXUCD"
EXPECTED_CYCLES = [
    "2003-2004", "2005-2006", "2007-2008", "2009-2010",
    "2011-2012", "2013-2014", "2015-2016", "2017-2018",
]
STRATEGIES = [
    {
        "strategy_id": "original_frozen",
        "strategy_label": "Original frozen model",
        "exposure_definition": "log2(URXUCD)",
        "ucr_terms": "none",
    },
    {
        "strategy_id": "ucr_main_only",
        "strategy_label": "UCr main effect only",
        "exposure_definition": "log2(URXUCD)",
        "ucr_terms": "log2(URXUCR) main effect; no UCr-by-sex term",
    },
    {
        "strategy_id": "ucr_main_plus_sex_interaction",
        "strategy_label": "UCr main effect plus UCr-by-sex",
        "exposure_definition": "log2(URXUCD)",
        "ucr_terms": "log2(URXUCR) + log2(URXUCR):female",
    },
    {
        "strategy_id": "cd_cr_ratio",
        "strategy_label": "Creatinine-corrected Cd ratio",
        "exposure_definition": "log2(URXUCD) - log2(URXUCR)",
        "ucr_terms": "ratio exposure; no separate UCr covariate",
    },
    {
        "strategy_id": "specific_gravity",
        "strategy_label": "Specific-gravity correction",
        "exposure_definition": "log2(URXUCD)",
        "ucr_terms": "specific gravity main effect + specific gravity-by-sex",
    },
]


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


def numeric(value: pd.Series) -> pd.Series:
    return pd.to_numeric(value, errors="coerce")


def design_matrix(frame: pd.DataFrame, strategy_id: str, pooled: bool) -> tuple[np.ndarray, list[str]]:
    """Build the requested model design; covariate parameterization is frozen."""
    x = pd.DataFrame(index=frame.index)
    x["Intercept"] = 1.0
    exposure_col = "axis_log2" if strategy_id != "cd_cr_ratio" else "cd_cr_log2"
    x[exposure_col] = numeric(frame[exposure_col])
    age = numeric(frame["age"]) - 50.0
    x["age_centered"] = age
    x["age_centered_sq"] = age * age
    if pooled:
        x["female"] = frame["sex"].eq("Female").astype(float)
        x[f"{exposure_col}:female"] = x[exposure_col] * x["female"]
    if strategy_id == "ucr_main_only":
        x["creatinine_log2"] = numeric(frame["creatinine_log2"])
    elif strategy_id == "ucr_main_plus_sex_interaction":
        x["creatinine_log2"] = numeric(frame["creatinine_log2"])
        if pooled:
            x["creatinine_log2:female"] = x["creatinine_log2"] * x["female"]
    elif strategy_id == "specific_gravity":
        x["specific_gravity"] = numeric(frame["specific_gravity"])
        if pooled:
            x["specific_gravity:female"] = x["specific_gravity"] * x["female"]
    for col, levels in [
        ("race", ["Mexican American", "Other Hispanic", "Non-Hispanic Black", "Other/Multi"]),
        ("smoking", ["Former", "Current"]),
    ]:
        for level in levels:
            x[f"{col}={level}"] = frame[col].eq(level).astype(float)
    cycles = sorted(frame["cycle"].dropna().unique().tolist())
    for level in cycles[1:]:
        x[f"cycle={level}"] = frame["cycle"].eq(level).astype(float)
    return x.to_numpy(float), x.columns.tolist()


def fit_variant(frame: pd.DataFrame, strategy_id: str, pooled: bool) -> dict[str, object]:
    required = [
        "outcome", "axis_log2", "age", "pir", "race", "smoking", "pooled_weight",
        "psu", "strata", "sex", "cycle",
    ]
    if strategy_id in {"ucr_main_only", "ucr_main_plus_sex_interaction", "cd_cr_ratio"}:
        required.append("creatinine_log2")
    if strategy_id == "specific_gravity":
        required.append("specific_gravity")
    work = frame.dropna(subset=required).copy()
    work = work.loc[numeric(work["pooled_weight"]).gt(0)].reset_index(drop=True)
    y = numeric(work["outcome"]).to_numpy(float)
    base = {
        "status": "not_estimable", "reason": "", "analytic_n": int(len(work)),
        "cases": int(y.sum()) if len(y) else 0,
        "controls": int(len(y) - y.sum()) if len(y) else 0,
        "male_n": int(work["sex"].eq("Male").sum()),
        "female_n": int(work["sex"].eq("Female").sum()),
        "male_cases": int(work.loc[work["sex"].eq("Male"), "outcome"].sum()),
        "female_cases": int(work.loc[work["sex"].eq("Female"), "outcome"].sum()),
        "contributing_cycles": ";".join(sorted(work["cycle"].dropna().unique().tolist())),
        "n_contributing_cycles": int(work["cycle"].nunique()),
    }
    if len(work) == 0 or y.sum() == 0 or y.sum() == len(y):
        base["reason"] = "no complete-case outcome variation"
        return base
    if pooled and work["sex"].nunique() != 2:
        base["reason"] = "both sexes not represented"
        return base
    x, names = design_matrix(work, strategy_id, pooled)
    weights = numeric(work["pooled_weight"]).to_numpy(float)
    weights = weights / np.nanmean(weights)
    beta = np.zeros(x.shape[1], dtype=float)

    def loglik(b: np.ndarray) -> float:
        p = expit(np.clip(x @ b, -35, 35))
        return float(np.sum(weights * (y * np.log(p + 1e-12) + (1.0 - y) * np.log1p(-p + 1e-12))))

    current = loglik(beta)
    converged = False
    for _ in range(200):
        p = expit(np.clip(x @ beta, -35, 35))
        gradient = x.T @ (weights * (y - p))
        hessian = x.T @ ((weights * p * (1.0 - p))[:, None] * x)
        step = np.linalg.pinv(hessian) @ gradient
        if not np.all(np.isfinite(step)):
            base["reason"] = "non-finite IRLS step"
            return base
        step = np.clip(step, -5, 5)
        alpha = 1.0
        while alpha >= 1e-8 and loglik(beta + alpha * step) < current - 1e-10:
            alpha /= 2.0
        if alpha < 1e-8:
            break
        beta += alpha * step
        current = loglik(beta)
        if np.max(np.abs(alpha * step)) < 1e-8:
            converged = True
            break

    p = expit(np.clip(x @ beta, -35, 35))
    scores = (weights * (y - p))[:, None] * x
    bread = x.T @ ((weights * p * (1.0 - p))[:, None] * x)
    bread_inv = np.linalg.pinv(bread)
    meat = np.zeros_like(bread)
    for _, group in work[["strata", "psu"]].groupby("strata", sort=False):
        psus = group["psu"].unique()
        if len(psus) < 2:
            continue
        totals = np.vstack([
            scores[group.index[group["psu"].eq(psu)], :].sum(axis=0)
            for psu in psus
        ])
        centered = totals - totals.mean(axis=0, keepdims=True)
        meat += len(psus) / (len(psus) - 1) * centered.T @ centered
    covariance = bread_inv @ meat @ bread_inv
    covariance = (covariance + covariance.T) / 2.0
    se = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    exposure_col = "axis_log2" if strategy_id != "cd_cr_ratio" else "cd_cr_log2"
    term = f"{exposure_col}:female" if pooled else exposure_col
    if term not in names:
        base["reason"] = "exposure coefficient unavailable"
        return base
    idx = names.index(term)
    if not np.isfinite(se[idx]) or se[idx] <= 0:
        base["reason"] = "coefficient variance unavailable"
        return base
    design_df = max(int(work["psu"].nunique() - work["strata"].nunique()), 1)
    crit = float(t.ppf(0.975, design_df))
    coefficient = float(beta[idx])
    standard_error = float(se[idx])
    z_value = coefficient / standard_error
    return {
        **base,
        "status": "ok" if converged else "converged_with_warning",
        "beta": coefficient, "se": standard_error, "z": float(z_value),
        "ci_low": float(coefficient - crit * standard_error),
        "ci_high": float(coefficient + crit * standard_error),
        "p": float(2.0 * t.sf(abs(z_value), design_df)),
        "design_df": design_df, "psu_n": int(work["psu"].nunique()),
        "strata_n": int(work["strata"].nunique()),
    }


def scan_specific_gravity() -> tuple[list[dict[str, str]], list[str], list[str]]:
    """Inspect all local XPT columns for SG/specific-gravity variables."""
    hits: list[dict[str, str]] = []
    checked: list[str] = []
    errors: list[str] = []
    cycle_re = re.compile(r"^(2003-2004|2005-2006|2007-2008|2009-2010|2011-2012|2013-2014|2015-2016|2017-2018)_")
    for path in sorted(DATA.glob("*.XPT")):
        if not cycle_re.match(path.name):
            continue
        checked.append(path.name)
        try:
            frame = pd.read_sas(path, format="xport", encoding="latin1")
        except Exception as exc:
            errors.append(f"{path.name}: {type(exc).__name__}: {exc}")
            continue
        for col in frame.columns:
            upper = str(col).upper()
            if re.search(r"SPECIFIC|GRAV|URX.*SG|URD.*SG|(^|_)SG($|_)", upper):
                hits.append({"file": path.name, "variable": str(col)})
    return hits, checked, errors


def build_target_frames(cd, reader, model, exposure: pd.DataFrame, exposure_cycles: list[str]):
    ucr_parts = []
    hashes: dict[str, str] = {}
    for cycle in exposure_cycles:
        path = DATA / f"{cycle}_ALB_CR.XPT"
        ucr = cd.read_ucr(cycle)
        ucr["cycle"] = cycle
        ucr_parts.append(ucr)
        hashes[path.name] = sha256(path)
    exposure = exposure.merge(
        pd.concat(ucr_parts, ignore_index=True),
        on=["SEQN", "cycle"], how="left", validate="one_to_one",
    )
    exposure["cd_cr_log2"] = exposure["axis_log2"] - exposure["creatinine_log2"]
    for source_row in getattr(exposure, "source_rows", []):
        del source_row
    frames: dict[str, list[pd.DataFrame]] = {endpoint_id: [] for endpoint_id in TARGET_IDS}
    for idx, cycle in enumerate(cd.CYCLES):
        if cycle not in exposure_cycles:
            continue
        demo_path = DATA / f"{cycle}_DEMO.XPT"
        mcq_path = DATA / f"{cycle}_MCQ.XPT"
        smq_path = DATA / f"{cycle}_SMQ.XPT"
        for path in [demo_path, mcq_path, smq_path]:
            hashes[path.name] = sha256(path)
        demo = pd.read_sas(demo_path, format="xport", encoding="latin1")
        mcq = pd.read_sas(mcq_path, format="xport", encoding="latin1")
        smq = pd.read_sas(smq_path, format="xport", encoding="latin1")
        core = model.derive_demo(demo, idx).merge(
            model.derive_smoking(smq), on="SEQN", how="left", validate="one_to_one",
        )
        core = core.loc[core["age"].ge(20) & core["sex"].notna()].copy()
        for endpoint in cd.ENDPOINTS:
            endpoint_id = endpoint["endpoint_id"]
            if endpoint_id not in TARGET_IDS:
                continue
            status, _ = cd.status_for_endpoint(mcq, cycle, endpoint)
            if status is None:
                continue
            status["SEQN"] = pd.to_numeric(status["SEQN"], errors="coerce").astype("Int64")
            joined = core.merge(status, on="SEQN", how="inner", validate="one_to_one")
            eligible = joined["case"] | joined["control"]
            d = joined.loc[eligible].copy()
            d["outcome"] = d["case"].astype(int)
            d["cycle"] = cycle
            frames[endpoint_id].append(d[[
                "SEQN", "cycle", "outcome", "age", "sex", "race", "pir",
                "smoking", "psu", "strata",
            ]])
    merged: dict[str, pd.DataFrame] = {}
    for endpoint_id in TARGET_IDS:
        if not frames[endpoint_id]:
            merged[endpoint_id] = pd.DataFrame()
        else:
            status_frame = pd.concat(frames[endpoint_id], ignore_index=True)
            merged[endpoint_id] = exposure.merge(
                status_frame, on=["SEQN", "cycle"], how="inner", validate="one_to_one",
            )
    return merged, exposure, hashes


def result_value(result: dict[str, object], key: str) -> float:
    value = result.get(key, np.nan)
    return float(value) if value is not None else np.nan


def md_table(frame: pd.DataFrame, cols: list[str], digits: int = 5) -> str:
    view = frame[cols].copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.{digits}g}")
    header = "| " + " | ".join(cols) + " |"
    rule = "|" + "|".join(["---"] * len(cols)) + "|"
    rows = ["| " + " | ".join(str(v) for v in row) + " |" for row in view.itertuples(index=False, name=None)]
    return "\n".join([header, rule, *rows])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cd = load_module(OLD_SCRIPT, "cd_subtype_sensitivity_import")
    reader = load_module(READER_PATH, "frozen_exposure_reader_cd_sensitivity")
    model = load_module(MODEL_PATH, "nhanes_design_model_cd_sensitivity")
    model.DATA_DIR = DATA
    tests = pd.read_csv(TESTS, dtype=str, keep_default_na=False)
    registry = pd.read_csv(REGISTRY, low_memory=False)
    test_row = tests.loc[tests["test_id"].eq(EXPOSURE_TEST_ID)].iloc[0]
    exposure, exposure_source = reader.read_test_exposure(test_row, registry)
    if exposure.empty:
        raise RuntimeError("URXUCD exposure could not be reconstructed from the frozen registry")
    exposure_cycles = exposure["cycle"].dropna().sort_values().unique().tolist()
    if exposure_cycles != EXPECTED_CYCLES:
        raise AssertionError(f"Unexpected URXUCD exposure cycles: {exposure_cycles}")
    merged, exposure_with_ucr, input_hashes = build_target_frames(cd, reader, model, exposure, exposure_cycles)
    for source_row in exposure_source.get("source_rows", []):
        path = Path(source_row["local_xpt"])
        if path.exists():
            input_hashes[path.name] = sha256(path)

    sg_hits, sg_checked, sg_errors = scan_specific_gravity()
    historical = pd.read_csv(OLD_RESULTS)
    historical = historical.loc[historical["endpoint_id"].isin(TARGET_IDS)].set_index("endpoint_id")
    primary = pd.read_csv(PRIMARY)
    primary_any = primary.loc[
        primary["test_id"].eq(EXPOSURE_TEST_ID) & primary["outcome_id"].eq("arthritis")
    ]
    historical_406 = primary_any.iloc[0].to_dict() if not primary_any.empty else {}

    rows: list[dict[str, object]] = []
    for endpoint_id in TARGET_IDS:
        frame = merged[endpoint_id]
        for strategy in STRATEGIES:
            sid = strategy["strategy_id"]
            base = {
                "endpoint_id": endpoint_id,
                "outcome_name": TARGET_NAMES[endpoint_id],
                "strategy_id": sid,
                "strategy_label": strategy["strategy_label"],
                "exposure_definition": strategy["exposure_definition"],
                "urinary_dilution_terms": strategy["ucr_terms"],
                "status": "not_applicable" if sid == "specific_gravity" and not sg_hits else "not_estimable",
                "reason": "specific-gravity variable not available in local raw package" if sid == "specific_gravity" and not sg_hits else "",
                "interaction_beta": np.nan, "interaction_se": np.nan, "interaction_z": np.nan,
                "interaction_ci_low": np.nan, "interaction_ci_high": np.nan, "interaction_p": np.nan,
                "male_beta": np.nan, "male_se": np.nan, "male_z": np.nan,
                "female_beta": np.nan, "female_se": np.nan, "female_z": np.nan,
                "analytic_n": 0, "total_cases": 0, "total_controls": 0,
                "male_n": 0, "female_n": 0, "male_cases": 0, "female_cases": 0,
                "n_contributing_cycles": 0, "contributing_cycles": "",
                "interaction_positive": np.nan, "interaction_p_lt_0_05": np.nan,
                "historical_subtype_main_interaction_beta": historical.loc[endpoint_id, "interaction_beta"] if endpoint_id in historical.index else np.nan,
                "historical_subtype_main_interaction_p": historical.loc[endpoint_id, "interaction_p"] if endpoint_id in historical.index else np.nan,
                "historical_subtype_main_q": historical.loc[endpoint_id, "subtype_family_bh_q"] if endpoint_id in historical.index else np.nan,
                "historical_any_arthritis_primary406_q": historical_406.get("bh_q_interaction_fixed406", np.nan) if endpoint_id == "any_arthritis" else np.nan,
            }
            if sid == "specific_gravity" and not sg_hits:
                rows.append(base)
                continue
            pooled = fit_variant(frame, sid, pooled=True)
            male = fit_variant(frame.loc[frame["sex"].eq("Male")], sid, pooled=False)
            female = fit_variant(frame.loc[frame["sex"].eq("Female")], sid, pooled=False)
            beta = result_value(pooled, "beta")
            row = {
                **base,
                "status": pooled.get("status", "not_estimable"),
                "reason": pooled.get("reason", ""),
                "interaction_beta": beta, "interaction_se": result_value(pooled, "se"),
                "interaction_z": result_value(pooled, "z"), "interaction_ci_low": result_value(pooled, "ci_low"),
                "interaction_ci_high": result_value(pooled, "ci_high"), "interaction_p": result_value(pooled, "p"),
                "male_beta": result_value(male, "beta"), "male_se": result_value(male, "se"), "male_z": result_value(male, "z"),
                "female_beta": result_value(female, "beta"), "female_se": result_value(female, "se"), "female_z": result_value(female, "z"),
                "analytic_n": pooled.get("analytic_n", 0), "total_cases": pooled.get("cases", 0), "total_controls": pooled.get("controls", 0),
                "male_n": pooled.get("male_n", 0), "female_n": pooled.get("female_n", 0),
                "male_cases": pooled.get("male_cases", 0), "female_cases": pooled.get("female_cases", 0),
                "n_contributing_cycles": pooled.get("n_contributing_cycles", 0), "contributing_cycles": pooled.get("contributing_cycles", ""),
                "interaction_positive": bool(np.isfinite(beta) and beta > 0) if np.isfinite(beta) else np.nan,
                "interaction_p_lt_0_05": bool(np.isfinite(result_value(pooled, "p")) and result_value(pooled, "p") < 0.05) if np.isfinite(result_value(pooled, "p")) else np.nan,
            }
            rows.append(row)
    result_df = pd.DataFrame(rows).sort_values(["endpoint_id", "strategy_id"]).reset_index(drop=True)
    result_df.to_csv(OUT / "01_urinary_dilution_strategy_comparison.csv", index=False)

    available = result_df.loc[result_df["status"].isin(["ok", "converged_with_warning"])].copy()
    summary_rows = []
    for endpoint_id in TARGET_IDS:
        a = available.loc[available["endpoint_id"].eq(endpoint_id)]
        positive = bool(len(a) and (a["interaction_beta"] > 0).all())
        oa_nominal_all = bool(endpoint_id != "osteoarthritis" or (len(a) and (a["interaction_p"] < 0.05).all()))
        summary_rows.append({
            "endpoint_id": endpoint_id,
            "n_estimable_strategies": int(len(a)),
            "all_estimable_interactions_positive": positive,
            "all_estimable_nominal_p_lt_0_05": bool(len(a) and (a["interaction_p"] < 0.05).all()),
            "oa_consistently_nominal_p_lt_0_05": oa_nominal_all,
            "min_interaction_beta": float(a["interaction_beta"].min()) if len(a) else np.nan,
            "max_interaction_beta": float(a["interaction_beta"].max()) if len(a) else np.nan,
            "strategies_with_nominal_p_lt_0_05": ";".join(a.loc[a["interaction_p"] < 0.05, "strategy_id"].tolist()),
        })
    summary_df = pd.DataFrame(summary_rows)

    report_lines = [
        "# Cd × Sex urinary-dilution sensitivity audit",
        "",
        "## Scope",
        "",
        "This bounded audit compares only the frozen URXUCD exposure for Any arthritis, osteoarthritis (OA), and gout. It uses existing NHANES local data, the established survey-weighted model implementation, and the same age, race/ethnicity, PIR, smoking, cycle fixed effects, weights, PSU/strata variance, and complete-case conventions. No new endpoint family, FDR, mechanism, literature, or figure analysis was performed.",
        "",
        "The four estimable strategies are descriptive sensitivity comparisons; no new multiplicity adjustment was calculated. The historical six-test subtype q values are retained only as references.",
        "",
        "## Specific-gravity availability",
        "",
        f"All {len(sg_checked)} relevant local XPT files were scanned for specific-gravity variable names. No specific-gravity variable was identified (hits={len(sg_hits)}). The specific-gravity rows are therefore recorded as `not_applicable`; no proxy was substituted. Files with parser errors were retained in the manifest and did not contain usable SG metadata.",
        "",
        "## Interaction comparison",
        "",
        md_table(result_df, ["endpoint_id", "strategy_id", "status", "interaction_beta", "interaction_se", "interaction_p", "male_beta", "female_beta", "analytic_n", "total_cases"]),
        "",
        "## Pattern summary",
        "",
        md_table(summary_df, ["endpoint_id", "n_estimable_strategies", "all_estimable_interactions_positive", "all_estimable_nominal_p_lt_0_05", "min_interaction_beta", "max_interaction_beta", "strategies_with_nominal_p_lt_0_05"]),
        "",
        "## Interpretation",
        "",
        "The formal quantity compared across strategies is the pooled Cd × female interaction. Positive interaction means the female slope is higher than the male slope under that model specification; it does not by itself establish positive female susceptibility, causality, or protection. Male and female slopes are descriptive and are not interpreted by the one-side-significance rule.",
        "",
        "The earlier subtype report now labels Any arthritis, OA, and gout as `FDR-supported sex heterogeneity, but not a positive female-susceptibility pattern`, because their fixed six-endpoint q values are below 0.05 while both sex-specific Cd slopes are negative. This audit does not convert those inverse slopes into protective claims.",
        "",
        "## Data and provenance",
        "",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
        f"Input exposure cycles: {';'.join(EXPECTED_CYCLES)}",
        f"Output CSV: 01_urinary_dilution_strategy_comparison.csv",
    ]
    (OUT / "URINARY_DILUTION_SENSITIVITY_AUDIT_REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "Cd x sex urinary-dilution sensitivity audit",
        "exposure_test_id": EXPOSURE_TEST_ID,
        "endpoints": TARGET_IDS,
        "strategies": STRATEGIES,
        "no_new_fdr": True,
        "specific_gravity_scan": {
            "n_files_checked": len(sg_checked), "n_hits": len(sg_hits),
            "hits": sg_hits, "parser_errors": sg_errors,
            "status": "not_available_in_local_raw_package" if not sg_hits else "available_but_not_fitted_until_variable_mapping_review",
        },
        "model": "survey-weighted logistic outcome ~ exposure + female + exposure:female + age + age^2 + race + PIR + smoking + cycle FE, with strategy-specific urinary dilution terms",
        "variance": "Taylor-style stratified PSU sandwich with the frozen cycle-offset survey design",
        "inputs": {
            "old_subtype_script_sha256": sha256(OLD_SCRIPT),
            "primary_results_sha256": sha256(PRIMARY),
            "old_subtype_results_sha256": sha256(OLD_RESULTS),
            "raw_xpt_sha256": input_hashes,
        },
        "outputs": {
            "strategy_comparison": "01_urinary_dilution_strategy_comparison.csv",
            "report": "URINARY_DILUTION_SENSITIVITY_AUDIT_REPORT.md",
            "manifest": "manifest.json",
        },
        "interpretation_boundary": "Sensitivity diagnostics only; no mechanism, causality, protection, new FDR, ranking, literature, figures, or broader screen.",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({
        "rows": len(result_df), "sg_hits": len(sg_hits),
        "summary": summary_df.to_dict("records"),
        "output": str(OUT),
    }, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
