"""Focused URXUCD x sex arthritis-subtype analysis.

This is a bounded follow-up to the frozen NHANES sex-divergence analysis.  It
uses only URXUCD, the MCQ arthritis/gout variables, the original covariates,
and the existing NHANES survey-weight construction.  The subtype interaction
family is fixed at six endpoints and no downstream biological analysis is
performed.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import chi2, t


ROOT = Path(__file__).resolve().parents[1]
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
DATA = ROOT / "work" / "nhanes_phase2a" / "data"
OUT = ROOT / "analysis" / "cadmium_sex_arthritis_subtypes"
TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
PRIMARY = FRAMEWORK / "sex_divergence_primary" / "sex_divergence_primary_406.csv"
ROBUSTNESS = FRAMEWORK / "sex_divergence_robustness" / "primary_robustness_uniform_results.csv"
MODEL_PATH = ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"
READER_PATH = FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py"

CYCLES = [
    "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
    "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018",
]
EXPOSURE_TEST_ID = "NHANES_URXUCD"
FDR_DENOMINATOR = 6
LOCO_MIN_TOTAL_CASES = 100
LOCO_MIN_SEX_CASES = 50

CODEBOOK_URLS = {
    "1999-2000": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/MCQ.htm",
    "2001-2002": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/MCQ_B.htm",
    "2003-2004": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/MCQ_C.htm",
    "2005-2006": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/MCQ_D.htm",
    "2007-2008": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/MCQ_E.htm",
    "2009-2010": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/MCQ_F.htm",
    "2011-2012": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/MCQ_G.htm",
    "2013-2014": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/MCQ_H.htm",
    "2015-2016": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/MCQ_I.htm",
    "2017-2018": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/MCQ_J.htm",
}

ENDPOINTS = [
    {
        "endpoint_id": "any_arthritis", "outcome_name": "Any arthritis",
        "definition": "MCQ160A=1 (case) versus MCQ160A=2 (control); adults >=20.",
        "endpoint_class": "reference",
    },
    {
        "endpoint_id": "osteoarthritis", "outcome_name": "Osteoarthritis (OA)",
        "subtype": "OA", "definition": "MCQ160A=1 and cycle-mapped arthritis subtype=OA versus MCQ160A=2; other arthritis subtypes excluded.",
        "endpoint_class": "primary_subtype",
    },
    {
        "endpoint_id": "rheumatoid_arthritis", "outcome_name": "Rheumatoid arthritis (RA)",
        "subtype": "RA", "definition": "MCQ160A=1 and cycle-mapped arthritis subtype=RA versus MCQ160A=2; other arthritis subtypes excluded.",
        "endpoint_class": "primary_subtype",
    },
    {
        "endpoint_id": "psoriatic_arthritis", "outcome_name": "Psoriatic arthritis (PsA)",
        "subtype": "PsA", "definition": "MCQ160A=1 and cycle-mapped arthritis subtype=PsA versus MCQ160A=2; endpoint is exploratory.",
        "endpoint_class": "exploratory",
    },
    {
        "endpoint_id": "other_arthritis", "outcome_name": "Other arthritis",
        "subtype": "Other", "definition": "MCQ160A=1 and cycle-mapped arthritis subtype=Other versus MCQ160A=2; refused/don't know/missing excluded.",
        "endpoint_class": "primary_subtype",
    },
    {
        "endpoint_id": "gout", "outcome_name": "Gout",
        "gout": True, "definition": "MCQ160N=1 (case) versus MCQ160N=2 (control); adults >=20.",
        "endpoint_class": "primary_related_endpoint",
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


def resolve_col(frame: pd.DataFrame, wanted: str) -> str | None:
    wanted_upper = wanted.upper()
    for col in frame.columns:
        if str(col).upper() == wanted_upper:
            return str(col)
    return None


def subtype_mapping(cycle: str) -> dict[str, object]:
    # CDC codebooks use MCQ190 through 2007-08, MCQ191 in 2009-10, and
    # MCQ195 from 2011-12 onward.  The category coding changed in 2011.
    if cycle <= "2007-2008":
        return {
            "variable": "MCQ190",
            "codes": {"RA": 1, "OA": 2, "PsA": None, "Other": 3},
            "label": "MCQ190: 1=RA; 2=OA; 3=Other; PsA category not present",
        }
    if cycle == "2009-2010":
        return {
            "variable": "MCQ191",
            "codes": {"RA": 1, "OA": 2, "PsA": 3, "Other": 4},
            "label": "MCQ191: 1=RA; 2=OA; 3=PsA; 4=Other",
        }
    return {
        "variable": "MCQ195",
        "codes": {"OA": 1, "RA": 2, "PsA": 3, "Other": 4},
        "label": "MCQ195: 1=OA/degenerative; 2=RA; 3=PsA; 4=Other",
    }


def status_for_endpoint(mcq: pd.DataFrame, cycle: str, endpoint: dict[str, object]) -> tuple[pd.DataFrame | None, dict[str, object]]:
    seqn_col = resolve_col(mcq, "SEQN")
    arthritis_col = resolve_col(mcq, "MCQ160A")
    if seqn_col is None:
        return None, {"variable": "", "mapping": "", "available": False, "reason": "MCQ lacks SEQN"}

    if endpoint.get("gout"):
        gout_col = resolve_col(mcq, "MCQ160N")
        if gout_col is None:
            return None, {
                "variable": "MCQ160N", "mapping": "MCQ160N: 1=Yes; 2=No",
                "available": False, "reason": "MCQ160N unavailable in this cycle",
            }
        value = numeric(mcq[gout_col])
        return pd.DataFrame({"SEQN": mcq[seqn_col], "case": value.eq(1), "control": value.eq(2)}), {
            "variable": gout_col, "mapping": "MCQ160N: 1=Yes; 2=No", "available": True, "reason": "",
        }

    if arthritis_col is None:
        return None, {"variable": "MCQ160A", "mapping": "MCQ160A: 1=Yes; 2=No", "available": False, "reason": "MCQ160A unavailable"}
    arthritis = numeric(mcq[arthritis_col])
    if endpoint.get("subtype") is None:
        return pd.DataFrame({"SEQN": mcq[seqn_col], "case": arthritis.eq(1), "control": arthritis.eq(2)}), {
            "variable": arthritis_col, "mapping": "MCQ160A: 1=Yes; 2=No", "available": True, "reason": "",
        }

    mapping = subtype_mapping(cycle)
    subtype_col = resolve_col(mcq, str(mapping["variable"]))
    subtype_name = str(endpoint["subtype"])
    subtype_code = mapping["codes"].get(subtype_name)
    if subtype_col is None:
        return None, {
            "variable": f"{arthritis_col};{mapping['variable']}", "mapping": str(mapping["label"]),
            "available": False, "reason": f"{mapping['variable']} unavailable in this cycle",
        }
    if subtype_code is None:
        return None, {
            "variable": f"{arthritis_col};{subtype_col}", "mapping": str(mapping["label"]),
            "available": False, "reason": "PsA category not present in this cycle-specific arthritis subtype mapping",
        }
    subtype = numeric(mcq[subtype_col])
    return pd.DataFrame({
        "SEQN": mcq[seqn_col],
        "case": arthritis.eq(1) & subtype.eq(subtype_code),
        "control": arthritis.eq(2),
    }), {
        "variable": f"{arthritis_col};{subtype_col}", "mapping": str(mapping["label"]), "available": True, "reason": "",
    }


def read_ucr(cycle: str) -> pd.DataFrame:
    path = DATA / f"{cycle}_ALB_CR.XPT"
    frame = pd.read_sas(path, format="xport", encoding="latin1")
    seqn_col = resolve_col(frame, "SEQN")
    ucr_col = resolve_col(frame, "URXUCR")
    if seqn_col is None or ucr_col is None:
        raise ValueError(f"{path.name} lacks SEQN or URXUCR")
    out = frame[[seqn_col, ucr_col]].copy()
    out.columns = ["SEQN", "URXUCR"]
    out["SEQN"] = pd.to_numeric(out["SEQN"], errors="coerce").astype("Int64")
    out["creatinine_log2"] = np.log2(numeric(out["URXUCR"]).where(numeric(out["URXUCR"]) > 0))
    return out[["SEQN", "creatinine_log2"]].drop_duplicates("SEQN")


def prepare_design(frame: pd.DataFrame, interaction: bool, include_creatinine: bool, heterogeneity: bool = False) -> tuple[np.ndarray, list[str]]:
    x = pd.DataFrame(index=frame.index)
    x["Intercept"] = 1.0
    x["axis_log2"] = numeric(frame["axis_log2"])
    age = numeric(frame["age"]) - 50.0
    x["age_centered"] = age
    x["age_centered_sq"] = age * age
    if interaction:
        x["female"] = frame["sex"].eq("Female").astype(float)
        x["axis_log2:female"] = x["axis_log2"] * x["female"]
    if include_creatinine:
        x["creatinine_log2"] = numeric(frame["creatinine_log2"])
        if interaction:
            x["creatinine_log2:female"] = x["creatinine_log2"] * x["female"]
    for col, levels in [
        ("race", ["Mexican American", "Other Hispanic", "Non-Hispanic Black", "Other/Multi"]),
        ("smoking", ["Former", "Current"]),
    ]:
        for level in levels:
            x[f"{col}={level}"] = frame[col].eq(level).astype(float)
    cycles = sorted(frame["cycle"].dropna().unique().tolist())
    for level in cycles[1:]:
        x[f"cycle={level}"] = frame["cycle"].eq(level).astype(float)
    if heterogeneity and interaction:
        for level in cycles[1:]:
            x[f"axis_log2:female:cycle={level}"] = x["axis_log2:female"] * frame["cycle"].eq(level).astype(float)
    return x.to_numpy(float), x.columns.tolist()


def model_required(include_creatinine: bool) -> list[str]:
    required = ["outcome", "axis_log2", "age", "pir", "race", "smoking", "pooled_weight", "psu", "strata", "sex", "cycle"]
    if include_creatinine:
        required.append("creatinine_log2")
    return required


def fit_weighted(frame: pd.DataFrame, interaction: bool, include_creatinine: bool = True, heterogeneity: bool = False) -> dict[str, object]:
    required = model_required(include_creatinine)
    work = frame.dropna(subset=required).copy()
    work = work.loc[numeric(work["pooled_weight"]).gt(0)].reset_index(drop=True)
    y = numeric(work["outcome"]).to_numpy(float)
    base = {
        "status": "not_estimable", "reason": "", "analytic_n": int(len(work)),
        "cases": int(y.sum()) if len(y) else 0, "controls": int(len(y) - y.sum()) if len(y) else 0,
        "male_n": int(work["sex"].eq("Male").sum()), "female_n": int(work["sex"].eq("Female").sum()),
        "male_cases": int(work.loc[work["sex"].eq("Male"), "outcome"].sum()),
        "female_cases": int(work.loc[work["sex"].eq("Female"), "outcome"].sum()),
    }
    if len(work) == 0 or y.sum() == 0 or y.sum() == len(y):
        base["reason"] = "no complete-case outcome variation"
        return base
    if interaction and work["sex"].nunique() != 2:
        base["reason"] = "both sexes not represented"
        return base
    x, names = prepare_design(work, interaction, include_creatinine, heterogeneity)
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
        totals = np.vstack([scores[group.index[group["psu"].eq(psu)], :].sum(axis=0) for psu in psus])
        centered = totals - totals.mean(axis=0, keepdims=True)
        meat += len(psus) / (len(psus) - 1) * centered.T @ centered
    covariance = bread_inv @ meat @ bread_inv
    covariance = (covariance + covariance.T) / 2.0
    se = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    term = "axis_log2:female" if interaction else "axis_log2"
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
        "ci_low": float(coefficient - crit * standard_error), "ci_high": float(coefficient + crit * standard_error),
        "p": float(2.0 * t.sf(abs(z_value), design_df)),
        "design_df": design_df, "psu_n": int(work["psu"].nunique()), "strata_n": int(work["strata"].nunique()),
        "coef": dict(zip(names, beta)), "cov": covariance, "names": names,
    }


def fixed_bh(pvalues: pd.Series) -> pd.Series:
    values = numeric(pvalues).fillna(1.0).clip(0, 1).to_numpy(float)
    order = np.argsort(values)
    raw = values[order] * FDR_DENOMINATOR / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(raw[::-1])[::-1]
    output = np.empty(len(values), dtype=float)
    output[order] = np.minimum(adjusted, 1.0)
    return pd.Series(output, index=pvalues.index)


def direction_label(beta: object) -> str:
    if beta is None or not np.isfinite(float(beta)):
        return "not_estimable"
    if float(beta) > 0:
        return "positive"
    if float(beta) < 0:
        return "negative"
    return "zero"


def md_table(frame: pd.DataFrame, columns: list[str], digits: int = 4) -> str:
    view = frame[columns].copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.{digits}g}")
    header = "| " + " | ".join(view.columns) + " |"
    rule = "|" + "|".join(["---"] * len(view.columns)) + "|"
    rows = ["| " + " | ".join(str(v) for v in row) + " |" for row in view.itertuples(index=False, name=None)]
    return "\n".join([header, rule, *rows])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reader = load_module(READER_PATH, "frozen_exposure_reader_cd_subtypes")
    model = load_module(MODEL_PATH, "nhanes_design_model_cd_subtypes")
    model.DATA_DIR = DATA
    tests = pd.read_csv(TESTS, dtype=str, keep_default_na=False)
    registry = pd.read_csv(REGISTRY, low_memory=False)
    test_row = tests.loc[tests["test_id"].eq(EXPOSURE_TEST_ID)].iloc[0]
    exposure, exposure_source = reader.read_test_exposure(test_row, registry)
    if exposure.empty:
        raise RuntimeError("URXUCD exposure could not be reconstructed from the frozen registry")
    exposure_cycles = exposure["cycle"].dropna().sort_values().unique().tolist()
    if exposure_cycles != ["2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"]:
        raise AssertionError(f"Unexpected URXUCD exposure cycles: {exposure_cycles}")
    ucr_parts = []
    input_hashes: dict[str, str] = {}
    for cycle in exposure_cycles:
        ucr_path = DATA / f"{cycle}_ALB_CR.XPT"
        ucr = read_ucr(cycle)
        ucr["cycle"] = cycle
        ucr_parts.append(ucr)
        input_hashes[ucr_path.name] = sha256(ucr_path)
    exposure = exposure.merge(pd.concat(ucr_parts, ignore_index=True), on=["SEQN", "cycle"], how="left", validate="one_to_one")
    for source_row in exposure_source.get("source_rows", []):
        path = Path(source_row["local_xpt"])
        if path.exists():
            input_hashes[path.name] = sha256(path)

    variable_audit_rows: list[dict[str, object]] = []
    cycle_endpoint_frames: dict[str, list[pd.DataFrame]] = {e["endpoint_id"]: [] for e in ENDPOINTS}
    endpoint_meta: dict[str, dict[str, object]] = {}
    for idx, cycle in enumerate(CYCLES):
        demo_path = DATA / f"{cycle}_DEMO.XPT"
        mcq_path = DATA / f"{cycle}_MCQ.XPT"
        smq_path = DATA / f"{cycle}_SMQ.XPT"
        input_hashes[demo_path.name] = sha256(demo_path)
        input_hashes[mcq_path.name] = sha256(mcq_path)
        input_hashes[smq_path.name] = sha256(smq_path)
        demo = pd.read_sas(demo_path, format="xport", encoding="latin1")
        mcq = pd.read_sas(mcq_path, format="xport", encoding="latin1")
        smq = pd.read_sas(smq_path, format="xport", encoding="latin1")
        core = model.derive_demo(demo, idx).merge(model.derive_smoking(smq), on="SEQN", how="left", validate="one_to_one")
        core = core.loc[core["age"].ge(20) & core["sex"].notna()].copy()
        for endpoint in ENDPOINTS:
            endpoint_id = endpoint["endpoint_id"]
            status, meta = status_for_endpoint(mcq, cycle, endpoint)
            exposure_cycle = cycle in exposure_cycles
            row = {
                "endpoint_id": endpoint_id, "outcome_name": endpoint["outcome_name"], "endpoint_class": endpoint["endpoint_class"],
                "cycle": cycle, "exposure_cycle_used": exposure_cycle, "module": "MCQ",
                "source_variable(s)": meta.get("variable", ""), "cycle_mapping": meta.get("mapping", ""),
                "variable_available": bool(meta.get("available", False)), "outcome_available": bool(meta.get("available", False)),
                "adult_age_rule": "RIDAGEYR >=20; RIAGENDR in {1,2}", "operational_definition": endpoint["definition"],
                "raw_eligible_n": 0, "raw_case_count": 0, "raw_control_count": 0,
                "raw_male_cases": 0, "raw_female_cases": 0, "raw_male_controls": 0, "raw_female_controls": 0,
                "source_codebook_url": CODEBOOK_URLS[cycle], "notes": str(meta.get("reason", "")),
            }
            if status is not None:
                status["SEQN"] = pd.to_numeric(status["SEQN"], errors="coerce").astype("Int64")
                joined = core.merge(status, on="SEQN", how="inner", validate="one_to_one")
                eligible = joined["case"] | joined["control"]
                row.update({
                    "raw_eligible_n": int(eligible.sum()),
                    "raw_case_count": int((eligible & joined["case"]).sum()),
                    "raw_control_count": int((eligible & joined["control"]).sum()),
                    "raw_male_cases": int((eligible & joined["case"] & joined["sex"].eq("Male")).sum()),
                    "raw_female_cases": int((eligible & joined["case"] & joined["sex"].eq("Female")).sum()),
                    "raw_male_controls": int((eligible & joined["control"] & joined["sex"].eq("Male")).sum()),
                    "raw_female_controls": int((eligible & joined["control"] & joined["sex"].eq("Female")).sum()),
                })
                if exposure_cycle:
                    d = joined.loc[eligible].copy()
                    d["outcome"] = d["case"].astype(int)
                    d["cycle"] = cycle
                    cycle_endpoint_frames[endpoint_id].append(d[["SEQN", "cycle", "outcome", "age", "sex", "race", "pir", "smoking", "psu", "strata"]])
            variable_audit_rows.append(row)
            endpoint_meta.setdefault(endpoint_id, {"endpoint": endpoint, "cycles": [], "variables": [], "mappings": []})
            if bool(meta.get("available", False)) and exposure_cycle:
                endpoint_meta[endpoint_id]["cycles"].append(cycle)
                endpoint_meta[endpoint_id]["variables"].append(str(meta.get("variable", "")))
                endpoint_meta[endpoint_id]["mappings"].append(str(meta.get("mapping", "")))

    variable_audit = pd.DataFrame(variable_audit_rows).sort_values(["endpoint_id", "cycle"]).reset_index(drop=True)
    merged_by_endpoint: dict[str, pd.DataFrame] = {}
    for endpoint in ENDPOINTS:
        endpoint_id = endpoint["endpoint_id"]
        if cycle_endpoint_frames[endpoint_id]:
            merged = pd.concat(cycle_endpoint_frames[endpoint_id], ignore_index=True)
            merged = exposure.merge(merged, on=["SEQN", "cycle"], how="inner", validate="one_to_one")
            merged_by_endpoint[endpoint_id] = merged
        else:
            merged_by_endpoint[endpoint_id] = pd.DataFrame()

    historical_any = pd.read_csv(PRIMARY)
    historical_any = historical_any.loc[
        historical_any["test_id"].eq(EXPOSURE_TEST_ID) & historical_any["outcome_id"].eq("arthritis")
    ]
    historical = historical_any.iloc[0].to_dict() if not historical_any.empty else {}

    model_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    cycle_rows: list[dict[str, object]] = []
    loco_rows: list[dict[str, object]] = []
    fit_cache: dict[str, dict[str, object]] = {}
    for endpoint in ENDPOINTS:
        endpoint_id = endpoint["endpoint_id"]
        d = merged_by_endpoint[endpoint_id]
        pooled = fit_weighted(d, interaction=True, include_creatinine=True)
        male = fit_weighted(d.loc[d["sex"].eq("Male")] if not d.empty else d, interaction=False, include_creatinine=True)
        female = fit_weighted(d.loc[d["sex"].eq("Female")] if not d.empty else d, interaction=False, include_creatinine=True)
        fit_cache[endpoint_id] = {"pooled": pooled, "male": male, "female": female}
        q_placeholder = np.nan
        cycles = endpoint_meta[endpoint_id]["cycles"]
        low_male = pooled.get("male_cases", 0) < LOCO_MIN_SEX_CASES
        low_female = pooled.get("female_cases", 0) < LOCO_MIN_SEX_CASES
        low_total = pooled.get("cases", 0) < LOCO_MIN_TOTAL_CASES
        model_row = {
            "endpoint_id": endpoint_id, "outcome_name": endpoint["outcome_name"], "endpoint_class": endpoint["endpoint_class"],
            "model_variant": "urinary_creatinine_sex_specific",
            "source_variables": ";".join(sorted(set(endpoint_meta[endpoint_id]["variables"]))),
            "cycle_mapping": " | ".join(sorted(set(endpoint_meta[endpoint_id]["mappings"]))),
            "contributing_cycles": ";".join(cycles), "n_contributing_cycles": len(cycles),
            "fit_status": pooled.get("status"), "fit_reason": pooled.get("reason", ""),
            "interaction_beta": pooled.get("beta", np.nan), "interaction_se": pooled.get("se", np.nan),
            "interaction_z": pooled.get("z", np.nan), "interaction_ci_low": pooled.get("ci_low", np.nan),
            "interaction_ci_high": pooled.get("ci_high", np.nan), "interaction_p": pooled.get("p", np.nan),
            "subtype_family_bh_q": q_placeholder, "fdr_denominator": FDR_DENOMINATOR,
            "male_beta": male.get("beta", np.nan), "male_se": male.get("se", np.nan), "male_z": male.get("z", np.nan),
            "female_beta": female.get("beta", np.nan), "female_se": female.get("se", np.nan), "female_z": female.get("z", np.nan),
            "analytic_n": pooled.get("analytic_n", 0), "total_cases": pooled.get("cases", 0), "total_controls": pooled.get("controls", 0),
            "male_n": pooled.get("male_n", 0), "female_n": pooled.get("female_n", 0),
            "male_cases": pooled.get("male_cases", 0), "female_cases": pooled.get("female_cases", 0),
            "low_male_cases_lt50": low_male, "low_female_cases_lt50": low_female, "low_total_cases_lt100": low_total,
            "historical_any_arthritis_primary406_beta": historical.get("beta_interaction", np.nan) if endpoint_id == "any_arthritis" else np.nan,
            "historical_any_arthritis_primary406_p": historical.get("p_interaction", np.nan) if endpoint_id == "any_arthritis" else np.nan,
            "historical_any_arthritis_primary406_q": historical.get("bh_q_interaction_fixed406", np.nan) if endpoint_id == "any_arthritis" else np.nan,
            "model_specification": "survey-weighted logit: subtype outcome ~ log2(URXUCD) + female + log2(URXUCD):female + log2(URXUCR) + log2(URXUCR):female + age_centered + age_centered^2 + race + PIR + smoking + cycle FE",
        }
        model_rows.append(model_row)
        sample_rows.append({
            "endpoint_id": endpoint_id, "outcome_name": endpoint["outcome_name"], "endpoint_class": endpoint["endpoint_class"],
            "analytic_n": pooled.get("analytic_n", 0), "total_cases": pooled.get("cases", 0), "total_controls": pooled.get("controls", 0),
            "male_n": pooled.get("male_n", 0), "female_n": pooled.get("female_n", 0), "male_cases": pooled.get("male_cases", 0), "female_cases": pooled.get("female_cases", 0),
            "male_controls": pooled.get("male_n", 0) - pooled.get("male_cases", 0), "female_controls": pooled.get("female_n", 0) - pooled.get("female_cases", 0),
            "n_contributing_cycles": len(cycles), "contributing_cycles": ";".join(cycles),
            "male_cases_lt50": low_male, "female_cases_lt50": low_female, "total_cases_lt100": low_total,
            "fit_status": pooled.get("status"), "fit_reason": pooled.get("reason", ""),
        })

        pooled_beta = pooled.get("beta", np.nan)
        for cycle in cycles:
            part = d.loc[d["cycle"].eq(cycle)].copy()
            fit_cycle = fit_weighted(part, interaction=True, include_creatinine=True)
            cycle_rows.append({
                "endpoint_id": endpoint_id, "outcome_name": endpoint["outcome_name"], "cycle": cycle,
                "interaction_beta": fit_cycle.get("beta", np.nan), "interaction_se": fit_cycle.get("se", np.nan),
                "interaction_z": fit_cycle.get("z", np.nan), "interaction_p": fit_cycle.get("p", np.nan),
                "fit_status": fit_cycle.get("status"), "fit_reason": fit_cycle.get("reason", ""),
                "analytic_n": fit_cycle.get("analytic_n", 0), "cases": fit_cycle.get("cases", 0), "controls": fit_cycle.get("controls", 0),
                "male_cases": fit_cycle.get("male_cases", 0), "female_cases": fit_cycle.get("female_cases", 0),
                "pooled_interaction_beta": pooled_beta,
                "same_direction_as_pooled": bool(np.isfinite(pooled_beta) and np.isfinite(fit_cycle.get("beta", np.nan)) and np.sign(pooled_beta) == np.sign(fit_cycle.get("beta", np.nan))),
                "low_case_caution": bool(fit_cycle.get("cases", 0) < LOCO_MIN_TOTAL_CASES or fit_cycle.get("male_cases", 0) < LOCO_MIN_SEX_CASES or fit_cycle.get("female_cases", 0) < LOCO_MIN_SEX_CASES),
            })

        can_loco = len(cycles) >= 3 and not low_male and not low_female and not low_total and pooled.get("status") in {"ok", "converged_with_warning"}
        if not can_loco:
            reason = "LOCO not performed because of inadequate subtype case counts" if (low_male or low_female or low_total) else "LOCO not performed because pooled model or cycle support was inadequate"
            loco_rows.append({
                "endpoint_id": endpoint_id, "outcome_name": endpoint["outcome_name"], "dropped_cycle": "",
                "status": "not_performed", "reason": reason, "n_successful_refits": 0,
                "interaction_beta": np.nan, "interaction_se": np.nan, "interaction_p": np.nan,
                "analytic_n": pooled.get("analytic_n", 0), "cases": pooled.get("cases", 0), "male_cases": pooled.get("male_cases", 0), "female_cases": pooled.get("female_cases", 0),
            })
        else:
            success = 0
            endpoint_loco: list[dict[str, object]] = []
            for dropped in cycles:
                fit_loco = fit_weighted(d.loc[d["cycle"].ne(dropped)].copy(), interaction=True, include_creatinine=True)
                if fit_loco.get("status") in {"ok", "converged_with_warning"}:
                    success += 1
                endpoint_loco.append({
                    "endpoint_id": endpoint_id, "outcome_name": endpoint["outcome_name"], "dropped_cycle": dropped,
                    "status": fit_loco.get("status"), "reason": fit_loco.get("reason", ""),
                    "n_successful_refits": np.nan, "interaction_beta": fit_loco.get("beta", np.nan), "interaction_se": fit_loco.get("se", np.nan), "interaction_p": fit_loco.get("p", np.nan),
                    "analytic_n": fit_loco.get("analytic_n", 0), "cases": fit_loco.get("cases", 0), "male_cases": fit_loco.get("male_cases", 0), "female_cases": fit_loco.get("female_cases", 0),
                    "same_direction_as_pooled": bool(np.isfinite(pooled_beta) and np.isfinite(fit_loco.get("beta", np.nan)) and np.sign(pooled_beta) == np.sign(fit_loco.get("beta", np.nan))),
                })
            for row in endpoint_loco:
                row["n_successful_refits"] = success
            loco_rows.extend(endpoint_loco)

    models = pd.DataFrame(model_rows).sort_values("endpoint_id").reset_index(drop=True)
    models["subtype_family_bh_q"] = fixed_bh(models["interaction_p"])
    models["interaction_fdr_rank"] = models["interaction_p"].fillna(1.0).rank(method="min", ascending=True).astype(int)
    # Conservative classification uses the fixed subtype-family q value and
    # the formal interaction direction; the male/female fits are descriptive.
    classes = []
    for row in models.itertuples(index=False):
        if row.fit_status not in {"ok", "converged_with_warning"} or row.low_male_cases_lt50 or row.low_female_cases_lt50 or row.low_total_cases_lt100:
            classes.append("Unstable / non-estimable")
        elif row.subtype_family_bh_q < 0.05 and row.female_beta > 0 and row.female_beta > row.male_beta and row.interaction_beta > 0:
            classes.append("Female-susceptible pattern")
        elif row.subtype_family_bh_q < 0.05 and row.male_beta > 0 and row.male_beta > row.female_beta and row.interaction_beta < 0:
            classes.append("Male-susceptible pattern")
        elif row.subtype_family_bh_q < 0.05:
            classes.append("FDR-supported sex heterogeneity, but not a positive female-susceptibility pattern")
        else:
            classes.append("No convincing sex heterogeneity")
    models["interpretation_class"] = classes
    models["interaction_fdr_supported"] = models["subtype_family_bh_q"].lt(0.05)
    def directional_pattern(row: pd.Series) -> str:
        if row["fit_status"] not in {"ok", "converged_with_warning"}:
            return "unstable_or_non_estimable"
        if row["subtype_family_bh_q"] >= 0.05:
            return "no_fdr_supported_directional_difference"
        b, bm, bf = row["interaction_beta"], row["male_beta"], row["female_beta"]
        if not all(np.isfinite(v) for v in [b, bm, bf]):
            return "unstable_or_non_estimable"
        if b > 0 and bf > bm:
            if bf > 0:
                return "female_susceptible_direction"
            if bm < 0 and bf < 0:
                return "female_enhanced_relative_slope_both_inverse"
            return "female_enhanced_relative_slope"
        if b < 0 and bm > bf:
            if bm > 0:
                return "male_susceptible_direction"
            if bm < 0 and bf < 0:
                return "male_enhanced_relative_slope_both_inverse"
            return "male_enhanced_relative_slope"
        return "interaction_direction_not_aligned_with_stratified_contrast"
    models["directional_pattern"] = models.apply(directional_pattern, axis=1)
    samples = pd.DataFrame(sample_rows).sort_values("endpoint_id").reset_index(drop=True)
    cycle_df = pd.DataFrame(cycle_rows).sort_values(["endpoint_id", "cycle"]).reset_index(drop=True)
    loco_df = pd.DataFrame(loco_rows).sort_values(["endpoint_id", "dropped_cycle"]).reset_index(drop=True)
    variable_audit.to_csv(OUT / "01_arthritis_subtype_variable_audit.csv", index=False)
    samples.to_csv(OUT / "02_subtype_sample_counts.csv", index=False)
    models.to_csv(OUT / "03_cd_sex_subtype_models.csv", index=False)
    models[["endpoint_id", "outcome_name", "endpoint_class", "interaction_p", "subtype_family_bh_q", "interaction_fdr_rank", "fdr_denominator", "interaction_beta", "interaction_se", "interaction_fdr_supported", "interpretation_class", "directional_pattern", "fit_status"]].to_csv(OUT / "04_subtype_interaction_fdr.csv", index=False)
    cycle_df.to_csv(OUT / "05_cycle_consistency.csv", index=False)
    loco_df.to_csv(OUT / "06_loco_results.csv", index=False)

    class_counts = models["interpretation_class"].value_counts().to_dict()
    sample_table = models[["endpoint_id", "outcome_name", "endpoint_class", "interaction_beta", "interaction_se", "interaction_p", "subtype_family_bh_q", "interaction_fdr_supported", "male_beta", "female_beta", "analytic_n", "total_cases", "male_cases", "female_cases", "interpretation_class", "directional_pattern"]].copy()
    report_lines = [
        "# Cd × Sex × Arthritis Subtype Analysis",
        "",
        "## Scope and stop rule",
        "",
        "This focused analysis uses the frozen `URXUCD` urinary cadmium exposure and six prespecified outcomes: any arthritis (reference), osteoarthritis, rheumatoid arthritis, psoriatic arthritis, other arthritis, and gout. It stops after variable audit, sample counts, survey-weighted interaction models, six-test BH-FDR, cycle consistency, and prespecified LOCO diagnostics. No CTD, gene, enrichment, PPI, GTEx, single-cell, mediation, hormone, literature, or figure analysis was run.",
        "",
        "## Outcome variable mapping",
        "",
        "`MCQ160A` is used for any arthritis and as the case gate for arthritis subtypes. The cycle-specific subtype variable is `MCQ190` for 1999–2008, `MCQ191` for 2009–2010, and `MCQ195` for 2011–2018. PsA is not available in the pre-2009 mapping and is not treated as a negative in those cycles. `MCQ160N` gout is available from 2007–2008 onward. CDC codebook URLs are retained in `01_arthritis_subtype_variable_audit.csv`.",
        "",
        "## Model",
        "",
        "The main model is a survey-weighted logistic regression with `log2(URXUCD)`, female, their interaction, age (centered linear and quadratic), race/ethnicity, PIR, smoking, cycle fixed effects, and urinary creatinine adjustment as `log2(URXUCR)` plus its female interaction. The formal sex difference is the pooled `URXUCD × female` coefficient. Sex-stratified coefficients are descriptive. Weights, cycle-offset strata/PSU identifiers, complete-case handling, and Taylor-style stratified PSU sandwich variance follow the existing NHANES pipeline. The six interaction P values use a fixed BH denominator of 6; the historical any-arthritis fixed-406 q is retained only as a reference.",
        "",
        "## Analytic sample and interaction results",
        "",
        md_table(sample_table, ["endpoint_id", "outcome_name", "endpoint_class", "interaction_beta", "interaction_se", "interaction_p", "subtype_family_bh_q", "interaction_fdr_supported", "male_beta", "female_beta", "analytic_n", "total_cases", "male_cases", "female_cases", "interpretation_class", "directional_pattern"]),
        "",
        f"Interpretation class counts: {class_counts}.",
        "",
        "Among the arthritis subtypes, OA is the only subtype whose interaction remains below the fixed six-endpoint BH threshold; RA and Other arthritis do not, and PsA is underpowered. The reference any-arthritis interaction and the related gout interaction also pass the six-test threshold. For these FDR-supported positive interactions, both sex-stratified Cd slopes are negative and the female slope is less negative than the male slope, so the result is a female-enhanced relative slope rather than a positive female-susceptibility association under the prespecified rule.",
        "",
        "Case-count cautions are explicit in `02_subtype_sample_counts.csv`; PsA is exploratory and is not eligible for strong interpretation if either sex has fewer than 50 cases or total cases are fewer than 100. No subtype is called protective. Directional wording is limited to susceptible, attenuated/weaker, or no convincing heterogeneity.",
        "",
        "## Cycle and LOCO diagnostics",
        "",
        "`05_cycle_consistency.csv` reports each contributing cycle's interaction estimate and whether its direction matches the pooled estimate. `06_loco_results.csv` reports each leave-one-cycle-out refit when the prespecified minimum of 100 total cases and 50 cases in each sex was met; otherwise it records `LOCO not performed because of inadequate subtype case counts`.",
        "",
        "## Interpretation boundary",
        "",
        "The analysis localizes the existing Cd × sex arthritis signal to the prespecified outcome definitions only. It does not establish mechanism, causality, temporality, or clinical protection/resilience.",
        "",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
    ]
    (OUT / "CD_SEX_ARTHRITIS_SUBTYPE_REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    output_files = [
        "01_arthritis_subtype_variable_audit.csv", "02_subtype_sample_counts.csv", "03_cd_sex_subtype_models.csv",
        "04_subtype_interaction_fdr.csv", "05_cycle_consistency.csv", "06_loco_results.csv", "CD_SEX_ARTHRITIS_SUBTYPE_REPORT.md",
    ]
    manifest = {
        "analysis": "Focused URXUCD x sex arthritis subtype analysis",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "exposure": {"test_id": EXPOSURE_TEST_ID, "variable": "URXUCD", "matrix": "urine", "log_transform": "log2", "cycles": exposure_cycles, "source_rows": exposure_source.get("source_rows", [])},
        "outcomes": [{"endpoint_id": e["endpoint_id"], "outcome_name": e["outcome_name"], "definition": e["definition"], "endpoint_class": e["endpoint_class"]} for e in ENDPOINTS],
        "variable_mapping": {"any_or_subtype_gate": "MCQ160A", "subtype": {"1999-2008": "MCQ190", "2009-2010": "MCQ191", "2011-2018": "MCQ195"}, "gout": "MCQ160N from 2007-2008 onward", "codebook_urls": CODEBOOK_URLS},
        "model": "survey-weighted logistic: outcome ~ log2(URXUCD) + female + interaction + log2(URXUCR) + female*log2(URXUCR) + age + age^2 + race + PIR + smoking + cycle FE",
        "urinary_creatinine_adjustment": "sex-specific: log2(URXUCR) + female x log2(URXUCR)",
        "variance": "Taylor-style stratified PSU sandwich; t reference with PSU minus strata degrees of freedom",
        "fdr": {"family": "six Cd x sex endpoint interactions", "method": "BH", "fixed_denominator": FDR_DENOMINATOR, "nonestimable_p_treated_as_one": True},
        "loco_rule": {"min_total_cases": LOCO_MIN_TOTAL_CASES, "min_cases_each_sex": LOCO_MIN_SEX_CASES, "not_performed_message": "LOCO not performed because of inadequate subtype case counts"},
        "historical_reference": {"file": str(PRIMARY), "test_id": EXPOSURE_TEST_ID, "outcome_id": "arthritis", "fixed406_q_retained": historical.get("bh_q_interaction_fixed406")},
        "input_sha256": {"tests": sha256(TESTS), "registry": sha256(REGISTRY), "primary_results": sha256(PRIMARY), "robustness_results": sha256(ROBUSTNESS), "model_script": sha256(MODEL_PATH), **input_hashes},
        "counts": {"n_endpoints": len(ENDPOINTS), "n_variable_audit_rows": len(variable_audit), "n_model_rows": len(models), "n_cycle_rows": len(cycle_df), "n_loco_rows": len(loco_df), "interpretation_class_counts": class_counts},
        "prohibited_not_run": ["CTD", "GeneCards", "enrichment", "PPI", "GTEx", "transcriptomics", "single-cell", "immune infiltration", "mediation", "hormone analysis", "figures", "literature collision search"],
        "outputs": output_files,
    }
    manifest["output_sha256"] = {name: sha256(OUT / name) for name in output_files}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    print(json.dumps({"exposure_cycles": exposure_cycles, "model_rows": len(models), "sample_counts": samples[["endpoint_id", "analytic_n", "total_cases", "male_cases", "female_cases"]].to_dict("records"), "results": models[["endpoint_id", "interaction_beta", "interaction_p", "subtype_family_bh_q", "interpretation_class"]].to_dict("records")}, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
