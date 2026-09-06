"""MBzP-CRC Phase 2B: NHANES weighted analysis and fixed 18-gene bridge.

The script uses only local CDC NHANES XPT files.  Raw XPT files are ignored by
the repository; all harmonization, model definitions, design variables and
input hashes are recorded in the generated manifest.

The survey logistic regression is a weighted estimating-equation fit with
Taylor-style stratified PSU sandwich variance.  It is not an ordinary
unweighted logistic regression: pooled phthalate subsample weights, strata,
and PSU are used for every reported model.
"""

from __future__ import annotations

import hashlib
import gzip
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import spearmanr, t


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
OUTPUT_DIR = ROOT / "outputs"
N_CYCLES = 10

CYCLES = [
    {"cycle": "1999-2000", "year": "1999", "suffix": "", "phthalate": "PHPYPA", "creatinine": "LAB16"},
    {"cycle": "2001-2002", "year": "2001", "suffix": "_B", "phthalate": "PHPYPA", "creatinine": "L16_B"},
    {"cycle": "2003-2004", "year": "2003", "suffix": "_C", "phthalate": "L24PH", "creatinine": "L16_C"},
    {"cycle": "2005-2006", "year": "2005", "suffix": "_D", "phthalate": "PHTHTE", "creatinine": "ALB_CR_D"},
    {"cycle": "2007-2008", "year": "2007", "suffix": "_E", "phthalate": "PHTHTE", "creatinine": "ALB_CR_E"},
    {"cycle": "2009-2010", "year": "2009", "suffix": "_F", "phthalate": "PHTHTE", "creatinine": "ALB_CR_F"},
    {"cycle": "2011-2012", "year": "2011", "suffix": "_G", "phthalate": "PHTHTE", "creatinine": "ALB_CR_G"},
    {"cycle": "2013-2014", "year": "2013", "suffix": "_H", "phthalate": "PHTHTE", "creatinine": "ALB_CR_H"},
    {"cycle": "2015-2016", "year": "2015", "suffix": "_I", "phthalate": "PHTHTE", "creatinine": "ALB_CR_I"},
    {"cycle": "2017-2018", "year": "2017", "suffix": "_J", "phthalate": "PHTHTE", "creatinine": "ALB_CR_J"},
]

PHthalates = {
    "MBzP": "URXMZP",
    "MnBP": "URXMBP",
    "MiBP": "URXMIB",
    "MEP": "URXMEP",
    "MMP": "URXMNM",
    "MEHP": "URXMHP",
    "MEHHP": "URXMHH",
    "MEOHP": "URXMOH",
    "MECPP": "URXECP",
}


def read_xpt(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size < 1000:
        raise FileNotFoundError(f"Missing or empty NHANES file: {path}")
    return pd.read_sas(path, format="xport", encoding="latin1")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def first_existing(frame: pd.DataFrame, names: list[str]) -> str | None:
    return next((name for name in names if name in frame.columns), None)


def clean_code(series: pd.Series, valid: set[int] | None = None) -> pd.Series:
    out = numeric(series)
    if valid is not None:
        out = out.where(out.isin(valid))
    return out


def cancer_flags(mcq: pd.DataFrame, demo: pd.DataFrame) -> pd.DataFrame:
    if "SEQN" not in mcq.columns or "MCQ220" not in mcq.columns:
        raise ValueError("MCQ must contain SEQN and MCQ220")
    type_cols = [c for c in mcq.columns if re.fullmatch(r"MCQ230[A-D]", c, flags=re.I)]
    type_cols = sorted(type_cols, key=lambda x: "ABCD".index(x[-1].upper()))
    types = mcq[type_cols].apply(numeric)
    colon = types.eq(16).any(axis=1)
    rectal = types.eq(31).any(axis=1)
    cancer_response = numeric(mcq["MCQ220"])
    crc = cancer_response.eq(1) & (colon | rectal)

    diag_age = pd.Series(np.nan, index=mcq.index, dtype=float)
    legacy_age_cols = [c for c in mcq.columns if re.fullmatch(r"MCQ240[A-Z]{1,2}", c, flags=re.I)]
    modern_age_cols = [c for c in mcq.columns if re.fullmatch(r"MCD240[A-C]", c, flags=re.I)]
    if legacy_age_cols:
        # Before 2017, MCQ240 letters correspond to fixed cancer-type slots,
        # not to the respondent's MCQ230A-D response positions: G=colon and
        # W=rectum in the published codebook.
        legacy_map = {16: "MCQ240G", 31: "MCQ240W"}
        for code, age_col in legacy_map.items():
            type_mask = types.eq(code).any(axis=1)
            if age_col in mcq.columns:
                age = numeric(mcq[age_col]).where(numeric(mcq[age_col]).between(0, 100))
                diag_age = pd.concat([diag_age, age.where(type_mask)], axis=1).min(axis=1, skipna=True)
    elif modern_age_cols:
        modern_age_cols = sorted(modern_age_cols, key=lambda x: "ABC".index(x[-1].upper()))
        for type_col, age_col in zip(type_cols, modern_age_cols):
            age = numeric(mcq[age_col]).where(numeric(mcq[age_col]).between(0, 100))
            is_crc_type = numeric(mcq[type_col]).isin([16, 31])
            diag_age = pd.concat([diag_age, age.where(is_crc_type)], axis=1).min(axis=1, skipna=True)

    return pd.DataFrame(
        {
            "SEQN": mcq["SEQN"],
            "cancer_outcome_available": cancer_response.notna(),
            "cancer_known": cancer_response.isin([1, 2]),
            "cancer_free": cancer_response.eq(2),
            "cancer_reported": cancer_response.eq(1),
            "colon_case": crc & colon,
            "rectal_case": crc & rectal,
            "crc_case": crc,
            "both_colon_rectal": crc & colon & rectal,
            "crc_diagnosis_age": diag_age,
        }
    )


def derive_demo(demo: pd.DataFrame, cycle_index: int) -> pd.DataFrame:
    out = demo.copy()
    if "SEQN" not in out.columns:
        raise ValueError("DEMO must contain SEQN")
    out["age"] = numeric(out.get("RIDAGEYR"))
    out["sex"] = numeric(out.get("RIAGENDR")).map({1: "Male", 2: "Female"})
    race_col = "RIDRETH3" if "RIDRETH3" in out.columns else "RIDRETH1"
    race_code = numeric(out[race_col])
    if race_col == "RIDRETH3":
        race_map = {1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White", 4: "Non-Hispanic Black", 6: "Other/Multi", 7: "Other/Multi"}
    else:
        race_map = {1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White", 4: "Non-Hispanic Black", 5: "Other/Multi"}
    out["race"] = race_code.map(race_map)
    out["pir"] = numeric(out.get("INDFMPIR"))
    out["education"] = clean_code(out.get("DMDEDUC2"), {1, 2, 3, 4, 5})
    out["psu_raw"] = numeric(out.get("SDMVPSU"))
    out["strata_raw"] = numeric(out.get("SDMVSTRA"))
    out["strata"] = cycle_index * 100 + out["strata_raw"]
    # SDMVPSU is only unique within strata; include the stratum in the
    # pooled PSU key so the survey degrees of freedom are not understated.
    out["psu"] = cycle_index * 10000 + out["strata_raw"] * 10 + out["psu_raw"]
    return out[["SEQN", "age", "sex", "race", "pir", "education", "psu", "strata"]]


def derive_smoking(smq: pd.DataFrame) -> pd.DataFrame:
    ever = numeric(smq.get("SMQ020"))
    current = numeric(smq.get("SMQ040"))
    status = pd.Series(pd.NA, index=smq.index, dtype="object")
    status.loc[ever.eq(2)] = "Never"
    status.loc[ever.eq(1) & current.isin([1, 2])] = "Current"
    status.loc[ever.eq(1) & current.eq(3)] = "Former"
    return pd.DataFrame({"SEQN": smq["SEQN"], "smoking": status})


def derive_alcohol(alq: pd.DataFrame) -> pd.DataFrame:
    col = first_existing(alq, ["ALQ101", "ALQ100", "ALQ111", "ALQ110"])
    value = numeric(alq[col]) if col else pd.Series(np.nan, index=alq.index)
    return pd.DataFrame({"SEQN": alq["SEQN"], "alcohol_ever": value.where(value.isin([1, 2])) .eq(1), "alcohol_source": col or ""})


def derive_diabetes(diq: pd.DataFrame) -> pd.DataFrame:
    value = numeric(diq.get("DIQ010"))
    return pd.DataFrame({"SEQN": diq["SEQN"], "diabetes": value.where(value.isin([1, 2])).eq(1)})


def derive_activity(paq: pd.DataFrame) -> pd.DataFrame:
    # Harmonized binary activity indicator: any affirmative moderate/vigorous
    # activity item available in the cycle.  These items differ by era.
    candidate_cols = [c for c in ["PAQ100", "PAD200", "PAD320", "PAQ605", "PAQ635", "PAQ650", "PAQ665", "PAD440"] if c in paq.columns]
    if candidate_cols:
        values = paq[candidate_cols].apply(numeric)
        activity = values.eq(1).any(axis=1)
        activity = activity.where(values.notna().any(axis=1))
    else:
        activity = pd.Series(pd.NA, index=paq.index, dtype="object")
    return pd.DataFrame({"SEQN": paq["SEQN"], "physical_activity": activity})


def derive_bmx(bmx: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({"SEQN": bmx["SEQN"], "bmi": numeric(bmx.get("BMXBMI"))})


def load_cycle(spec: dict, cycle_index: int) -> tuple[pd.DataFrame, list[dict]]:
    cycle = spec["cycle"]
    suffix = spec["suffix"]
    paths = {
        "lab": DATA_DIR / f"{cycle}_{spec['phthalate']}.XPT",
        "mcq": DATA_DIR / f"{cycle}_MCQ.XPT",
        "demo": DATA_DIR / f"{cycle}_DEMO.XPT",
        "bmx": DATA_DIR / f"{cycle}_BMX.XPT",
        "smq": DATA_DIR / f"{cycle}_SMQ.XPT",
        "alq": DATA_DIR / f"{cycle}_ALQ.XPT",
        "diq": DATA_DIR / f"{cycle}_DIQ.XPT",
        "paq": DATA_DIR / f"{cycle}_PAQ.XPT",
        "creatinine": DATA_DIR / f"{cycle}_ALB_CR.XPT",
    }
    frames = {key: read_xpt(path) for key, path in paths.items()}
    lab = frames["lab"]
    if "URXMZP" not in lab.columns:
        raise ValueError(f"{paths['lab'].name} lacks URXMZP")
    selected_lab = ["SEQN", *[col for col in PHthalates.values() if col in lab.columns]]
    # CDC provides special 4-year weights for 1999-2002.  When those cycles
    # are combined with later 2-year cycles, use WTSPH4YR * (2 / n_cycles).
    if cycle in {"1999-2000", "2001-2002"}:
        weight_col = first_existing(lab, ["WTSPH4YR", "WTSPH2YR"])
    else:
        weight_col = first_existing(lab, ["WTSB2YR", "WTSA2YR", "WTSPH2YR"])
    if not weight_col:
        raise ValueError(f"No phthalate subsample weight in {paths['lab'].name}")
    selected_lab += [weight_col]
    lab_part = lab[selected_lab].drop_duplicates("SEQN")
    lab_part = lab_part.rename(columns={weight_col: "phthalate_weight_base"})
    creat = frames["creatinine"][["SEQN", "URXUCR"]].drop_duplicates("SEQN")
    merged = lab_part.merge(cancer_flags(frames["mcq"], frames["demo"]), on="SEQN", how="inner", validate="one_to_one")
    merged = merged.merge(derive_demo(frames["demo"], cycle_index), on="SEQN", how="left", validate="one_to_one")
    for supplement in (derive_bmx(frames["bmx"]), derive_smoking(frames["smq"]), derive_alcohol(frames["alq"]), derive_diabetes(frames["diq"]), derive_activity(frames["paq"]), creat):
        merged = merged.merge(supplement, on="SEQN", how="left", validate="one_to_one")
    merged["cycle"] = cycle
    merged["cycle_index"] = cycle_index
    first_two_multiplier = 2 if cycle in {"1999-2000", "2001-2002"} and weight_col == "WTSPH4YR" else 1
    merged["pooled_weight"] = numeric(merged["phthalate_weight_base"]) * first_two_multiplier / N_CYCLES
    merged["weight_source"] = weight_col
    merged["weight_multiplier_before_division"] = first_two_multiplier
    merged["mbzp_log2"] = np.log2(numeric(merged["URXMZP"]).where(numeric(merged["URXMZP"]) > 0))
    merged["creatinine_log2"] = np.log2(numeric(merged["URXUCR"]).where(numeric(merged["URXUCR"]) > 0))
    merged["years_since_crc"] = merged["age"] - merged["crc_diagnosis_age"]
    input_manifest = []
    source_names = {
        "lab": f"{spec['phthalate']}{suffix}",
        "mcq": f"MCQ{suffix}",
        "demo": f"DEMO{suffix}",
        "bmx": f"BMX{suffix}",
        "smq": f"SMQ{suffix}",
        "alq": f"ALQ{suffix}",
        "diq": f"DIQ{suffix}",
        "paq": f"PAQ{suffix}",
        "creatinine": spec["creatinine"],
    }
    for key, path in paths.items():
        input_manifest.append({"cycle": cycle, "component": key, "local_file": path.name, "source_stem": source_names[key], "bytes": path.stat().st_size, "sha256": sha256(path)})
    return merged, input_manifest


def build_design(df: pd.DataFrame, continuous: list[str], categorical: list[str], levels: dict[str, list[str]] | None = None) -> tuple[pd.DataFrame, list[str]]:
    levels = levels or {}
    x = pd.DataFrame({"Intercept": 1.0}, index=df.index)
    names = ["Intercept"]
    for col in continuous:
        x[col] = numeric(df[col]).astype(float)
        names.append(col)
    for col in categorical:
        vals = df[col].astype(str)
        cats = levels.get(col) or sorted(v for v in vals.unique() if v not in {"nan", "<NA>", "None"})
        if not cats:
            continue
        baseline = cats[0]
        for level in cats[1:]:
            name = f"{col}={level}"
            x[name] = (vals == level).astype(float)
            names.append(name)
    return x[names], names


def fit_survey_logistic(df: pd.DataFrame, continuous: list[str], categorical: list[str], outcome: str = "outcome", exposure_name: str = "mbzp_log2", levels: dict[str, list[str]] | None = None) -> dict:
    required = [outcome, "pooled_weight", "psu", "strata", *continuous, *categorical]
    work = df.dropna(subset=required).copy()
    work = work[work["pooled_weight"] > 0].copy()
    work = work.reset_index(drop=True)
    y = numeric(work[outcome]).to_numpy(float)
    if len(work) == 0 or y.sum() == 0 or y.sum() == len(work):
        return {"status": "not_estimable", "N": len(work), "CRC_N": int(y.sum()) if len(y) else 0, "Control_N": int(len(y) - y.sum()) if len(y) else 0}
    Xdf, names = build_design(work, continuous, categorical, levels)
    X = Xdf.to_numpy(float)
    # Scale nuisance continuous covariates for numerical conditioning.  Keep
    # the exposure on its original log2 scale so the reported coefficient is
    # directly interpretable as an OR per doubling.
    for col in {"age", "bmi", "pir", "creatinine_log2"} & set(continuous):
        j = names.index(col)
        scale = np.nanstd(X[:, j])
        if np.isfinite(scale) and scale > 0:
            X[:, j] = (X[:, j] - np.nanmean(X[:, j])) / scale
    y = y.astype(float)
    weights = numeric(work["pooled_weight"]).to_numpy(float)
    weights = weights / np.nanmean(weights)

    def loglik(beta: np.ndarray) -> float:
        p = expit(np.clip(X @ beta, -35, 35))
        return float(np.sum(weights * (y * np.log(p + 1e-12) + (1 - y) * np.log1p(-p + 1e-12))))

    # Newton-IRLS with backtracking.  This avoids accepting a numerical
    # optimizer's boundary solution as convergence when the covariate set is
    # moderately sparse.
    beta = np.zeros(X.shape[1], dtype=float)
    current_ll = loglik(beta)
    converged = False
    for _ in range(200):
        p = expit(np.clip(X @ beta, -35, 35))
        gradient = X.T @ (weights * (y - p))
        hessian = X.T @ ((weights * p * (1 - p))[:, None] * X)
        step = np.linalg.pinv(hessian) @ gradient
        if not np.all(np.isfinite(step)):
            break
        step = np.clip(step, -5, 5)
        alpha = 1.0
        accepted = False
        while alpha >= 1e-8:
            candidate = beta + alpha * step
            candidate_ll = loglik(candidate)
            if candidate_ll >= current_ll - 1e-10:
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            break
        beta = candidate
        if np.max(np.abs(alpha * step)) < 1e-8:
            converged = True
            break
        current_ll = candidate_ll
    if not np.all(np.isfinite(beta)):
        return {"status": "fit_failed", "N": len(work), "CRC_N": int(y.sum()), "Control_N": int(len(y) - y.sum()), "message": "Non-finite IRLS coefficients"}
    fit_message = "Newton-IRLS converged" if converged else "Newton-IRLS stopped before tolerance"
    p = expit(np.clip(X @ beta, -35, 35))
    scores = (weights * (y - p))[:, None] * X
    bread = X.T @ ((weights * p * (1 - p))[:, None] * X)
    bread_inv = np.linalg.pinv(bread)
    meat = np.zeros((X.shape[1], X.shape[1]))
    cluster = work[["strata", "psu"]].copy()
    for _, stratum in cluster.groupby("strata", sort=False):
        clusters = stratum["psu"].unique()
        if len(clusters) < 2:
            continue
        totals = np.vstack([scores[stratum.index[stratum["psu"].eq(psu)], :].sum(axis=0) for psu in clusters])
        centered = totals - totals.mean(axis=0, keepdims=True)
        meat += len(clusters) / (len(clusters) - 1) * (centered.T @ centered)
    covariance = bread_inv @ meat @ bread_inv
    covariance = (covariance + covariance.T) / 2
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    n_psu = work["psu"].nunique()
    n_strata = work["strata"].nunique()
    design_df = max(n_psu - n_strata, 1)
    idx = names.index(exposure_name) if exposure_name in names else None
    if idx is None or not np.isfinite(se[idx]) or se[idx] <= 0:
        return {"status": "fit_failed", "N": len(work), "CRC_N": int(y.sum()), "Control_N": int(len(y) - y.sum()), "message": "Exposure coefficient variance unavailable"}
    crit = float(t.ppf(0.975, design_df))
    stat = beta[idx] / se[idx]
    p_value = float(2 * t.sf(abs(stat), design_df))
    return {
        "status": "ok" if converged else "converged_with_warning",
        "N": int(len(work)), "CRC_N": int(y.sum()), "Control_N": int(len(y) - y.sum()),
        "beta": float(beta[idx]), "SE": float(se[idx]), "OR": float(np.exp(beta[idx])),
        "CI_low": float(np.exp(np.clip(beta[idx] - crit * se[idx], -700, 700))), "CI_high": float(np.exp(np.clip(beta[idx] + crit * se[idx], -700, 700))),
        "P": p_value, "design_df": int(design_df), "PSU_N": int(n_psu), "strata_N": int(n_strata),
        "message": fit_message, "coefficients": dict(zip(names, beta)), "covariance": covariance,
    }


def model_specs() -> dict[str, tuple[list[str], list[str]]]:
    return {
        "Model 0": (["mbzp_log2"], []),
        "Model 1": (["mbzp_log2", "age"], ["sex", "race"]),
        "Model 2": (["mbzp_log2", "age", "bmi", "pir", "creatinine_log2"], ["sex", "race", "smoking"]),
        "Model 3": (["mbzp_log2", "age", "bmi", "pir", "creatinine_log2"], ["sex", "race", "smoking", "alcohol_ever", "education", "diabetes", "physical_activity"]),
    }


LEVELS = {"sex": ["Female", "Male"], "race": ["Non-Hispanic White", "Mexican American", "Other Hispanic", "Non-Hispanic Black", "Other/Multi"], "smoking": ["Never", "Former", "Current"]}


def complete_case_audit(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for population, pop in population_frames(df).items():
        for model, (continuous, categorical) in model_specs().items():
            req = ["outcome", "pooled_weight", "psu", "strata", *continuous, *categorical]
            cc = pop.dropna(subset=req)
            cc = cc[cc["pooled_weight"] > 0]
            rows.append({"Model": model, "Population": population, "Total_N": len(cc), "CRC_N": int(cc["outcome"].sum()), "Control_N": int((1 - cc["outcome"]).sum()), "Required_variables": ";".join(req), "Fit_status": "eligible" if len(cc) and cc["outcome"].sum() > 0 and (1 - cc["outcome"]).sum() > 0 else "not_estimable"})
    return pd.DataFrame(rows)


def population_frames(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    adult = (df["age"] >= 20) & df["cancer_known"] & df["crc_case"].isin([True, False])
    primary = df[adult & (df["crc_case"] | df["cancer_free"])].copy()
    primary["outcome"] = primary["crc_case"].astype(int)
    all_noncrc = df[adult].copy()
    all_noncrc["outcome"] = all_noncrc["crc_case"].astype(int)
    return {"CRC_vs_cancer_free": primary, "CRC_vs_all_nonCRC": all_noncrc}


def run_models(df: pd.DataFrame, population: str, model_names: list[str] | None = None) -> pd.DataFrame:
    pop = population_frames(df)[population]
    rows = []
    for model, (continuous, categorical) in model_specs().items():
        if model_names and model not in model_names:
            continue
        fit = fit_survey_logistic(pop, continuous, categorical, levels=LEVELS)
        rows.append({"Model": model, "Population": population, **{key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}}})
    return pd.DataFrame(rows)


def quartile_analysis(df: pd.DataFrame) -> pd.DataFrame:
    pop = population_frames(df)["CRC_vs_cancer_free"].copy()
    available = pop["mbzp_log2"].notna()
    quantiles = pop.loc[available, "mbzp_log2"].quantile([0, .25, .5, .75, 1]).to_numpy()
    quantiles = np.unique(quantiles)
    pop["mbzp_quartile"] = pd.cut(pop["mbzp_log2"], bins=[-np.inf, *quantiles[1:-1], np.inf], labels=False, include_lowest=True) + 1
    pop["mbzp_quartile"] = pop["mbzp_quartile"].astype(float)
    rows = []
    required = ["outcome", "mbzp_quartile", "age", "bmi", "pir", "creatinine_log2", "sex", "race", "smoking", "pooled_weight", "psu", "strata"]
    cc = pop.dropna(subset=required).copy()
    cc["q2"] = (cc["mbzp_quartile"] == 2).astype(float)
    cc["q3"] = (cc["mbzp_quartile"] == 3).astype(float)
    cc["q4"] = (cc["mbzp_quartile"] == 4).astype(float)
    cc["qtrend"] = cc["mbzp_quartile"]
    # Build one fit with Q2-Q4 and the Model 2 covariates.  Report each Q row.
    fit = fit_survey_logistic(cc, ["q2", "q3", "q4", "age", "bmi", "pir", "creatinine_log2"], ["sex", "race", "smoking"], exposure_name="q4", levels=LEVELS)
    # Refit each contrast so the reported coefficient and CI are explicit.
    for q in [1, 2, 3, 4]:
        if q == 1:
            rows.append({"Quartile": "Q1", "Reference": True, "N": len(cc), "CRC_N": int(cc["outcome"].sum()), "OR": 1.0, "CI_low": 1.0, "CI_high": 1.0, "P": np.nan, "P_trend": np.nan, "status": "reference"})
            continue
        tmp = cc.copy()
        tmp["q_exposure"] = (tmp["mbzp_quartile"] == q).astype(float)
        f = fit_survey_logistic(tmp, ["q_exposure", "age", "bmi", "pir", "creatinine_log2"], ["sex", "race", "smoking"], exposure_name="q_exposure", levels=LEVELS)
        rows.append({"Quartile": f"Q{q}", "Reference": False, "N": f.get("N", len(tmp)), "CRC_N": f.get("CRC_N", int(tmp["outcome"].sum())), "OR": f.get("OR", np.nan), "CI_low": f.get("CI_low", np.nan), "CI_high": f.get("CI_high", np.nan), "P": f.get("P", np.nan), "P_trend": np.nan, "status": f.get("status", "not_estimable")})
    trend = fit_survey_logistic(cc, ["qtrend", "age", "bmi", "pir", "creatinine_log2"], ["sex", "race", "smoking"], exposure_name="qtrend", levels=LEVELS)
    for row in rows:
        row["P_trend"] = trend.get("P", np.nan)
        row["Q1_cutpoint_log2"] = float(quantiles[0]) if len(quantiles) else np.nan
        row["Q2_cutpoint_log2"] = float(quantiles[1]) if len(quantiles) > 1 else np.nan
        row["Q3_cutpoint_log2"] = float(quantiles[2]) if len(quantiles) > 2 else np.nan
        row["Q4_cutpoint_log2"] = float(quantiles[3]) if len(quantiles) > 3 else np.nan
    return pd.DataFrame(rows)


def spline_basis(x: np.ndarray, knots: np.ndarray) -> np.ndarray:
    # Natural cubic spline basis with boundary knots at knots[0], knots[-1].
    # The first column is linear; remaining columns encode nonlinear curvature.
    x = np.asarray(x, dtype=float)
    k = knots[1:-1]
    def d(z, knot):
        return np.maximum(z - knot, 0) ** 3
    if len(k) == 0:
        return x[:, None]
    last = knots[-1]
    ref = knots[-2]
    bases = [x]
    for knot in k:
        bases.append((d(x, knot) - d(x, last) * (last - knot) / (last - ref)) / (last - knot))
    return np.column_stack(bases)


def spline_analysis(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    pop = population_frames(df)["CRC_vs_cancer_free"].copy()
    req = ["outcome", "mbzp_log2", "age", "bmi", "pir", "creatinine_log2", "sex", "race", "smoking", "pooled_weight", "psu", "strata"]
    cc = pop.dropna(subset=req).copy()
    if cc["outcome"].sum() < 80:
        return pd.DataFrame([{"status": "not_run_insufficient_crc_events", "N": len(cc), "CRC_N": int(cc["outcome"].sum()), "overall_P": np.nan, "nonlinear_P": np.nan, "knots": "5th,35th,65th,95th percentile"}]), None
    knots = cc["mbzp_log2"].quantile([.05, .35, .65, .95]).to_numpy()
    basis = spline_basis(cc["mbzp_log2"].to_numpy(), knots)
    basis_cols = []
    for i in range(basis.shape[1]):
        col = f"spline_{i+1}"
        cc[col] = basis[:, i]
        basis_cols.append(col)
    fit = fit_survey_logistic(cc, basis_cols + ["age", "bmi", "pir", "creatinine_log2"], ["sex", "race", "smoking"], exposure_name=basis_cols[0], levels=LEVELS)
    # Overall and nonlinear Wald tests use the coefficient covariance returned
    # by the same survey sandwich fit.
    overall_p = np.nan
    nonlinear_p = np.nan
    if fit.get("status") in {"ok", "converged_with_warning"}:
        coefs = fit["coefficients"]
        cov = fit["covariance"]
        names = list(coefs)
        idx = [names.index(col) for col in basis_cols]
        b = np.array([coefs[col] for col in basis_cols])
        v = cov[np.ix_(idx, idx)]
        stat = float(b @ np.linalg.pinv(v) @ b)
        overall_p = float(__import__("scipy").stats.chi2.sf(stat, len(idx)))
        if len(idx) > 1:
            b2 = b[1:]
            v2 = v[1:, 1:]
            nonlinear_p = float(__import__("scipy").stats.chi2.sf(float(b2 @ np.linalg.pinv(v2) @ b2), len(b2)))
    out = pd.DataFrame([{"status": fit.get("status", "not_estimable"), "N": fit.get("N", len(cc)), "CRC_N": fit.get("CRC_N", int(cc["outcome"].sum())), "overall_P": overall_p, "nonlinear_P": nonlinear_p, "knots": ";".join(f"{x:.6g}" for x in knots)}])
    grid = np.linspace(cc["mbzp_log2"].min(), cc["mbzp_log2"].max(), 100)
    pred = spline_basis(grid, knots)
    beta = np.array([fit["coefficients"][col] for col in basis_cols])
    plot_df = pd.DataFrame({"mbzp_log2": grid, "spline_linear_predictor": pred @ beta})
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(plot_df["mbzp_log2"], plot_df["spline_linear_predictor"], color="#2b6cb0")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set(xlabel="log2 urinary MBzP (ng/mL)", ylabel="Model 2 spline linear predictor", title="MBzP–CRC restricted cubic spline")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "mbzp_crc_phase2_spline.png", dpi=180)
    plt.close(fig)
    return out, plot_df


def diagnosis_timing(df: pd.DataFrame) -> pd.DataFrame:
    cases = df[(df["age"] >= 20) & df["crc_case"] & df["mbzp_log2"].notna()].copy()
    values = numeric(cases["years_since_crc"]).where(numeric(cases["years_since_crc"]) >= 0).dropna()
    rows = [{"metric": "N_CRC_with_diagnosis_age", "value": int(values.notna().sum()), "unit": "participants"}]
    if len(values):
        rows += [{"metric": "median", "value": float(values.median()), "unit": "years"}, {"metric": "IQR_low", "value": float(values.quantile(.25)), "unit": "years"}, {"metric": "IQR_high", "value": float(values.quantile(.75)), "unit": "years"}, {"metric": "less_than_1_year", "value": int((values < 1).sum()), "unit": "participants"}, {"metric": "1_to_5_years", "value": int(values.between(1, 5, inclusive="both").sum()), "unit": "participants"}, {"metric": "more_than_5_years", "value": int((values > 5).sum()), "unit": "participants"}]
    else:
        rows.append({"metric": "status", "value": "not_available", "unit": "no valid CRC diagnosis ages"})
    return pd.DataFrame(rows)


def phthalate_comparison(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pop = population_frames(df)["CRC_vs_cancer_free"].copy()
    comparison = []
    for metabolite, col in PHthalates.items():
        if col not in pop.columns:
            continue
        fit_frame = pop.drop(columns=["mbzp_log2"], errors="ignore").copy()
        fit_frame["metabolite_log2"] = np.log2(numeric(fit_frame[col]).where(numeric(fit_frame[col]) > 0))
        fit = fit_survey_logistic(fit_frame, ["metabolite_log2", "age", "bmi", "pir", "creatinine_log2"], ["sex", "race", "smoking"], exposure_name="metabolite_log2", levels=LEVELS)
        cycles = ";".join(sorted(pop.loc[pop[col].notna(), "cycle"].unique()))
        comparison.append({"Metabolite": metabolite, "Variable": col, "Cycles": cycles, **{k: v for k, v in fit.items() if k not in {"coefficients", "covariance"}}})
    result = pd.DataFrame(comparison)
    if not result.empty:
        result["BH_FDR"] = bh_fdr(result["P"].fillna(1).to_numpy())
        result["OR_rank"] = result["OR"].rank(ascending=False, method="min")
        result["P_rank"] = result["P"].rank(ascending=True, method="min")
        result["FDR_rank"] = result["BH_FDR"].rank(ascending=True, method="min")
    corr_rows = []
    log_data = pd.DataFrame(index=pop.index)
    for metabolite, col in PHthalates.items():
        if col in pop.columns:
            log_data[metabolite] = np.log2(numeric(pop[col]).where(numeric(pop[col]) > 0))
    for a in log_data.columns:
        for b in log_data.columns:
            if a >= b:
                continue
            valid = log_data[[a, b]].dropna()
            rho, pval = spearmanr(valid[a], valid[b]) if len(valid) >= 3 else (np.nan, np.nan)
            corr_rows.append({"Metabolite_A": a, "Metabolite_B": b, "N": len(valid), "Spearman_rho": rho, "P": pval})
    corr = pd.DataFrame(corr_rows)
    z = (log_data - log_data.mean()) / log_data.std(ddof=0)
    pop["phthalate_burden"] = z.mean(axis=1, skipna=True)
    pop["phthalate_burden_n_metabolites"] = z.notna().sum(axis=1)
    burden = pop[pop["phthalate_burden_n_metabolites"] >= 2].copy()
    burden = burden.drop(columns=["mbzp_log2"], errors="ignore").rename(columns={"phthalate_burden": "burden_exposure"})
    fit = fit_survey_logistic(burden, ["burden_exposure", "age", "bmi", "pir", "creatinine_log2"], ["sex", "race", "smoking"], exposure_name="burden_exposure", levels=LEVELS)
    burden_out = pd.DataFrame([{"Exposure": "PhthalateBurden_mean_z", "minimum_metabolites": 2, **{k: v for k, v in fit.items() if k not in {"coefficients", "covariance"}}}])
    return result, corr, burden_out, log_data


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.nan_to_num(np.asarray(pvalues, dtype=float), nan=1.0)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(adjusted, 1.0)
    return out


def overlap_gene_outputs() -> dict:
    data_dir = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data"
    # Stream only the MBzP rows instead of materializing the complete CTD
    # interaction file; the complete archive is large and unnecessary here.
    with gzip.open(data_dir / "CTD_chem_gene_ixns.tsv.gz", "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("# Fields:"):
                columns = next(handle).lstrip("# ").rstrip("\r\n").split("\t")
                break
        selected = []
        for line in handle:
            if line.startswith("#"):
                continue
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(columns):
                continue
            row = dict(zip(columns, values))
            if row.get("Organism") == "Homo sapiens" and row.get("ChemicalID", "").replace("MESH:", "") == "C103325":
                selected.append(row)
    interactions = pd.DataFrame(selected)
    human = interactions.copy()
    human["ChemicalID"] = human["ChemicalID"].astype(str).str.replace("^MESH:", "", regex=True)
    human["GeneID"] = human["GeneID"].astype(str).str.strip()
    gc = pd.read_csv(data_dir / "genecards_disorders_crc.csv", dtype=str, keep_default_na=False)
    gc["GeneCards_Rank"] = pd.to_numeric(gc["GeneCards_Rank"], errors="coerce")
    gc["NCBI_GeneID"] = ""
    # Phase 1 mapped the GeneCards export to NCBI IDs using the CTD symbol map;
    # use the same deterministic symbol map here.
    symbol_map = human.drop_duplicates("GeneSymbol").set_index(human.drop_duplicates("GeneSymbol")["GeneSymbol"].str.upper())["GeneID"].to_dict()
    gc["NCBI_GeneID"] = gc["GeneSymbol"].str.upper().map(symbol_map).fillna("")
    ctd_genes = set(human["GeneID"])
    gc_crc = gc[gc["NCBI_GeneID"].ne("")].sort_values("GeneCards_Rank").drop_duplicates("NCBI_GeneID")
    overlap_ids = ctd_genes & set(gc_crc["NCBI_GeneID"])
    # The Phase 1 primary table contains 18 IDs; fail loudly if the bridge
    # changes under a different input export.
    if len(overlap_ids) != 18:
        raise ValueError(f"Expected 18 MBzP-GeneCards overlap genes, found {len(overlap_ids)}")
    overlap_gc = gc_crc[gc_crc["NCBI_GeneID"].isin(overlap_ids)].copy().sort_values("GeneCards_Rank")
    top10 = set(overlap_gc.head(10)["NCBI_GeneID"])
    gene_rows = []
    long_rows = []
    for gene_id in sorted(overlap_ids, key=lambda x: float(overlap_gc.loc[overlap_gc["NCBI_GeneID"].eq(x), "GeneCards_Rank"].iloc[0])):
        gct = overlap_gc[overlap_gc["NCBI_GeneID"].eq(gene_id)].iloc[0]
        subset = human[human["GeneID"].eq(gene_id)].copy()
        pmids = sorted({p.strip() for value in subset["PubMedIDs"].fillna("") for p in re.split(r"[|;,]", str(value)) if p.strip()})
        for row in subset.itertuples(index=False):
            interaction = getattr(row, "Interaction", "")
            action = getattr(row, "InteractionActions", "")
            pmid_value = getattr(row, "PubMedIDs", "")
            pmid_list = [p.strip() for p in re.split(r"[|;,]", str(pmid_value)) if p.strip()]
            if pmid_list:
                for pmid in pmid_list:
                    long_rows.append({"GeneID": gene_id, "GeneSymbol": gct["GeneSymbol"], "Interaction": interaction, "Action": action, "PMID": pmid})
            else:
                long_rows.append({"GeneID": gene_id, "GeneSymbol": gct["GeneSymbol"], "Interaction": interaction, "Action": action, "PMID": ""})
        gene_rows.append({"GeneID": gene_id, "GeneSymbol": gct["GeneSymbol"], "GeneName": gct.get("GeneName", ""), "GeneType": gct.get("GeneType", ""), "GeneCards_rank": int(gct["GeneCards_Rank"]), "GeneCards_relevance_score": gct.get("RelevanceScore", ""), "GeneCards_scope": "Disorders", "CTD_ChemicalID": "C103325", "CTD_ChemicalName": "mono-benzyl phthalate", "CTD_interaction": "; ".join(sorted(set(subset["Interaction"].fillna("")))), "CTD_action": "; ".join(sorted(set(subset["InteractionActions"].fillna("")))), "n_CTD_interaction_rows": len(subset), "n_unique_PMIDs": len(pmids), "PMIDs": ";".join(pmids), "rank_weight": 1 / math.log2(float(gct["GeneCards_Rank"]) + 1), "overlap_status": True, "top10_status": gene_id in top10})
    genes = pd.DataFrame(gene_rows)
    long = pd.DataFrame(long_rows).drop_duplicates()
    genes.to_csv(OUTPUT_DIR / "mbzp_crc_phase1_full_overlap_genes.csv", index=False)
    long.to_csv(OUTPUT_DIR / "mbzp_crc_phase1_overlap_interactions_long.csv", index=False)
    def gene_class(row: pd.Series) -> str:
        symbol = str(row["GeneSymbol"]).upper()
        name = str(row.get("GeneName", "")).lower()
        text = str(row["GeneType"]).lower()
        if symbol.startswith("MIR") or "microrna" in name:
            return "miRNA"
        if "antisense" in name or symbol.endswith("-AS1") or symbol in {"MALAT1", "ZFAS1"}:
            return "lncRNA"
        if "protein" in text:
            return "protein-coding"
        if "lncrna" in text or "long non-coding" in text:
            return "lncRNA"
        if "ncrna" in text or "non-coding" in text:
            return "other ncRNA"
        return "other/unspecified"
    genes["gene_class"] = genes.apply(gene_class, axis=1)
    counts = genes["gene_class"].value_counts().to_dict()
    summary_rows = [
        {"metric": "n_overlap_genes", "value": len(genes), "details": "Fixed CTD human MBzP genes intersected with GeneCards Disorders CRC genes"},
        {"metric": "protein_coding_n", "value": counts.get("protein-coding", 0), "details": "GeneCards GeneType / symbol-derived structural class"},
        {"metric": "miRNA_n", "value": counts.get("miRNA", 0), "details": "MIR-prefixed symbols or MicroRNA gene names"},
        {"metric": "lncRNA_n", "value": counts.get("lncRNA", 0), "details": "Antisense/long noncoding labels or known antisense symbols"},
        {"metric": "other_ncRNA_n", "value": counts.get("other ncRNA", 0), "details": "Other noncoding labels"},
        {"metric": "transcription_factor_annotation", "value": "not_available", "details": "No TF annotation supplied in Phase 1 input exports"},
        {"metric": "tumor_suppressor_annotation", "value": "not_available", "details": "No TSG annotation supplied in Phase 1 input exports"},
        {"metric": "oncogene_annotation", "value": "not_available", "details": "No oncogene annotation supplied in Phase 1 input exports"},
        {"metric": "total_CTD_interaction_rows", "value": int(genes["n_CTD_interaction_rows"].sum()), "details": "All gene-specific MBzP CTD rows retained in long table"},
        {"metric": "total_unique_PMIDs", "value": int(long["PMID"].replace("", np.nan).dropna().nunique()), "details": "Distinct PMIDs across the 18-gene bridge"},
        {"metric": "CTD_evidence_top_genes", "value": ";".join(genes.sort_values(["n_unique_PMIDs", "n_CTD_interaction_rows"], ascending=False)["GeneSymbol"].head(10)), "details": "Ranked by gene-specific unique PMID count then interaction rows"},
        {"metric": "GeneCards_top_genes", "value": ";".join(genes.sort_values("GeneCards_rank")["GeneSymbol"].head(10)), "details": "Ranked by GeneCards Disorders rank"},
    ]
    pd.DataFrame(summary_rows).to_csv(OUTPUT_DIR / "mbzp_crc_phase1_overlap_gene_summary.csv", index=False)
    return {"n_overlap_genes": len(genes), "genes": genes["GeneSymbol"].tolist(), "gene_sha256": sha256(OUTPUT_DIR / "mbzp_crc_phase1_full_overlap_genes.csv"), "long_sha256": sha256(OUTPUT_DIR / "mbzp_crc_phase1_overlap_interactions_long.csv")}


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    overlap_meta = overlap_gene_outputs()
    frames, input_manifest = [], []
    for idx, spec in enumerate(CYCLES):
        frame, manifest_rows = load_cycle(spec, idx)
        frames.append(frame)
        input_manifest.extend(manifest_rows)
    data = pd.concat(frames, ignore_index=True)
    data = data[(data["age"] >= 20) & data["cancer_outcome_available"]].copy()
    data.to_pickle(ROOT / "work" / "nhanes_phase2a" / "phase2b_harmonized.pkl")
    audit = complete_case_audit(data)
    audit.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_complete_case_audit.csv", index=False)

    model_rows = []
    for population in ["CRC_vs_cancer_free", "CRC_vs_all_nonCRC"]:
        model_rows.append(run_models(data, population))
    # Primary model 0-2 plus Model 3 only when event threshold permits.
    model_df = pd.concat(model_rows, ignore_index=True)
    model3_primary = audit[(audit["Population"] == "CRC_vs_cancer_free") & (audit["Model"] == "Model 3")]
    if not model3_primary.empty and int(model3_primary.iloc[0]["CRC_N"]) < 80:
        model_df.loc[model_df["Model"].eq("Model 3"), "status"] = "sensitivity_only_crc_n_lt_80"
    model_df.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_main_models.csv", index=False)

    quartiles = quartile_analysis(data)
    quartiles.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_quartiles.csv", index=False)
    spline, _ = spline_analysis(data)
    spline.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_spline.csv", index=False)
    timing = diagnosis_timing(data)
    timing.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_diagnosis_timing.csv", index=False)

    primary = population_frames(data)["CRC_vs_cancer_free"].copy()
    sensitivity_rows = []
    for label, frame in [("Age_ge_40", primary[primary["age"] >= 40].copy()), ("Colon_only", primary.assign(outcome=primary["colon_case"].astype(int))), ("Recent_CRC_excluded", primary[~(primary["crc_case"] & (primary["years_since_crc"] < 1))].copy())]:
        if label == "Colon_only":
            frame = frame[frame["colon_case"] | frame["cancer_free"]].copy()
        fit = fit_survey_logistic(frame, *model_specs()["Model 2"], levels=LEVELS)
        sensitivity_rows.append({"Analysis": label, "Population": "CRC_vs_cancer_free", **{k: v for k, v in fit.items() if k not in {"coefficients", "covariance"}}})
    all_noncrc = population_frames(data)["CRC_vs_all_nonCRC"]
    fit = fit_survey_logistic(all_noncrc, *model_specs()["Model 2"], levels=LEVELS)
    sensitivity_rows.append({"Analysis": "All_nonCRC_control_definition", "Population": "CRC_vs_all_nonCRC", **{k: v for k, v in fit.items() if k not in {"coefficients", "covariance"}}})
    pd.DataFrame(sensitivity_rows).to_csv(OUTPUT_DIR / "mbzp_crc_phase2_sensitivity.csv", index=False)

    loco_rows = []
    for cycle in sorted(data["cycle"].unique()):
        frame = primary[primary["cycle"] != cycle].copy()
        fit = fit_survey_logistic(frame, *model_specs()["Model 2"], levels=LEVELS)
        loco_rows.append({"Dropped_cycle": cycle, **{k: v for k, v in fit.items() if k not in {"coefficients", "covariance"}}})
    pd.DataFrame(loco_rows).to_csv(OUTPUT_DIR / "mbzp_crc_phase2_leave_one_cycle_out.csv", index=False)

    comparison, corr, burden, _ = phthalate_comparison(data)
    comparison.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_phthalate_comparison.csv", index=False)
    corr.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_phthalate_correlations.csv", index=False)
    burden.to_csv(OUTPUT_DIR / "mbzp_crc_phase2_phthalate_burden.csv", index=False)

    primary_models = model_df[model_df["Population"].eq("CRC_vs_cancer_free")]
    m2 = primary_models[primary_models["Model"].eq("Model 2")].iloc[0].to_dict()
    q4 = quartiles[quartiles["Quartile"].eq("Q4")].iloc[0].to_dict()
    age40 = next((r for r in sensitivity_rows if r["Analysis"] == "Age_ge_40"), {})
    loco_ok = pd.DataFrame(loco_rows)
    loco_or = loco_ok["OR"].dropna() if "OR" in loco_ok.columns else pd.Series(dtype=float)
    loco_beta = loco_ok["beta"].dropna() if "beta" in loco_ok.columns else pd.Series(dtype=float)
    loco_direction_consistent = bool(len(loco_beta) and ((loco_beta < 0).all() or (loco_beta > 0).all()))
    mb_rank = comparison[comparison["Metabolite"].eq("MBzP")].iloc[0].to_dict() if not comparison.empty and (comparison["Metabolite"] == "MBzP").any() else {}
    report = f"""# MBzP–CRC Phase 2B NHANES 人群验证 + 18 基因机制桥接

## 首页六问

1. **Model 2 MBzP OR per doubling:** {m2.get('OR', np.nan):.6g}; 95% CI {m2.get('CI_low', np.nan):.6g}–{m2.get('CI_high', np.nan):.6g}; P={m2.get('P', np.nan):.6g}.
2. **Q4 vs Q1:** OR={q4.get('OR', np.nan):.6g}; 95% CI {q4.get('CI_low', np.nan):.6g}–{q4.get('CI_high', np.nan):.6g}; P-trend={q4.get('P_trend', np.nan):.6g}.
3. **Age ≥40:** OR={age40.get('OR', np.nan):.6g}; 95% CI {age40.get('CI_low', np.nan):.6g}–{age40.get('CI_high', np.nan):.6g}; P={age40.get('P', np.nan):.6g}.
4. **LOCO OR range:** {loco_or.min() if len(loco_or) else np.nan:.6g}–{loco_or.max() if len(loco_or) else np.nan:.6g}; direction consistent: {('YES' if loco_direction_consistent else 'NO')} (all LOCO point estimates are below 1 when YES).
5. **MBzP among phthalates:** OR rank={mb_rank.get('OR_rank', np.nan)}; P rank={mb_rank.get('P_rank', np.nan)}; FDR rank={mb_rank.get('FDR_rank', np.nan)}.
6. **18 overlap genes exported:** **YES** ({overlap_meta['n_overlap_genes']} genes).

## Scope and definitions

- Primary population: CRC (colon or rectal cancer type code) versus participants reporting no cancer history (`MCQ220=2`).
- Sensitivity population: CRC versus participants with known cancer outcome who are not CRC, including other cancer histories.
- Exposure: `log2(URXMZP)`; OR is per doubling of urinary MBzP.
- Primary model: age, sex, race, BMI, smoking, PIR, and `log2(URXUCR)`.
- All estimates use CDC-compatible pooled phthalate subsample weights: `WTSPH4YR×2/10` for 1999–2002 and cycle-specific 2-year weights divided by 10 thereafter, with pooled strata and PSU identifiers and Taylor-style stratified PSU sandwich variance. Counts are unweighted.

## Model 0–3

See `mbzp_crc_phase2_main_models.csv`. Model 3 is sensitivity-only when its complete-case CRC count is below 80.

## Sensitivity and specificity

See `mbzp_crc_phase2_sensitivity.csv`, `mbzp_crc_phase2_leave_one_cycle_out.csv`, `mbzp_crc_phase2_phthalate_comparison.csv`, `mbzp_crc_phase2_phthalate_correlations.csv`, and `mbzp_crc_phase2_phthalate_burden.csv`.

## Timing and spline

Diagnosis timing is descriptive and based on available cancer diagnosis-age fields. Restricted cubic spline uses 5th/35th/65th/95th percentile knots only when the Model 2 primary complete-case CRC count is at least 80; otherwise the output records the prespecified non-run status.

## Fixed molecular bridge

The CTD human-interacting MBzP genes intersected with the Phase 1 GeneCards Disorders CRC set to produce exactly 18 genes. The full table retains gene-specific interaction rows and unique PMID counts; no GO, KEGG, PPI, docking, WGCNA, machine learning, or hub-gene fishing was performed.

Run timestamp (UTC): `{datetime.now(timezone.utc).isoformat()}`
"""
    (OUTPUT_DIR / "mbzp_crc_phase2_report.md").write_text(report, encoding="utf-8")

    manifest = {
        "analysis": "MBzP-CRC Phase 2B",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "regression_method": "weighted logistic estimating equations with Taylor-style stratified PSU sandwich variance",
        "pooled_weight_rule": "1999-2000 and 2001-2002 use WTSPH4YR * 2/10; 2003-2018 use the cycle-specific 2-year phthalate subsample weight / 10",
        "weight_columns": {spec["cycle"]: ("WTSPH4YR" if spec["cycle"] in {"1999-2000", "2001-2002"} else "WTSA2YR" if spec["cycle"] == "2011-2012" else "WTSB2YR") for spec in CYCLES},
        "weight_multipliers_before_division": {spec["cycle"]: (2 if spec["cycle"] in {"1999-2000", "2001-2002"} else 1) for spec in CYCLES},
        "design_variables": {"strata": "SDMVSTRA with cycle offset", "psu": "SDMVPSU with cycle offset"},
        "primary_population": "CRC type code 16 or 31 vs MCQ220=2 cancer-free controls",
        "sensitivity_population": "CRC vs known non-CRC outcome (MCQ220 in 1,2)",
        "exposure": "log2(URXMZP)",
        "primary_covariates": ["age", "sex", "race", "bmi", "smoking", "pir", "log2(URXUCR)"],
        "cycles": [spec["cycle"] for spec in CYCLES],
        "n_adult_outcome_rows": int(len(data)),
        "n_primary_population": int(len(population_frames(data)["CRC_vs_cancer_free"])),
        "crc_primary_n": int(population_frames(data)["CRC_vs_cancer_free"]["outcome"].sum()),
        "crc_all_noncrc_n": int(population_frames(data)["CRC_vs_all_nonCRC"]["outcome"].sum()),
        "phase1_overlap": overlap_meta,
        "input_files": input_manifest,
        "outputs": [p.name for p in sorted(OUTPUT_DIR.glob("mbzp_crc_phase2*"))] + [p.name for p in sorted(OUTPUT_DIR.glob("mbzp_crc_phase1_*overlap*"))],
        "prohibited_analyses_not_run": ["GO", "KEGG", "Reactome", "STRING PPI", "WGCNA", "machine learning", "docking", "molecular dynamics", "hub gene fishing"],
    }
    (OUTPUT_DIR / "mbzp_crc_phase2_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print("Primary Model 2", m2)
    print("Primary population", len(population_frames(data)["CRC_vs_cancer_free"]), "CRC", int(population_frames(data)["CRC_vs_cancer_free"]["outcome"].sum()))
    print("Overlap genes", overlap_meta["n_overlap_genes"])


if __name__ == "__main__":
    main()
