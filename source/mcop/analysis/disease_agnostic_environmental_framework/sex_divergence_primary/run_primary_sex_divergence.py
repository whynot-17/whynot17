"""Run the locked 29 x 14 primary NHANES sex-divergence analysis.

This program is intentionally limited to the primary package.  It produces no
candidate selection, literature/CTD query, robustness analysis, or narrative
interpretation.  The formal estimand is the pooled exposure-by-female term.
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
from scipy.special import expit
from scipy.stats import t


ROOT = Path(__file__).resolve().parents[3]
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
PLAN_DIR = FRAMEWORK / "sex_divergence_plan"
OUTCOME_DIR = FRAMEWORK / "outcome_inventory"
TESTS = FRAMEWORK / "step04_testset_freeze" / "unique_biomarker_test_set.csv"
REGISTRY = FRAMEWORK / "data_processed" / "detectability_registry_outcome_blinded.csv"
DATA = ROOT / "work" / "nhanes_phase2a" / "data"
OUT = Path(__file__).resolve().parent
N_INTERACTIONS = 406


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def build_design(work: pd.DataFrame, interaction: bool) -> tuple[np.ndarray, list[str]]:
    """Fixed primary parameterization; age is centered only for conditioning."""
    x = pd.DataFrame(index=work.index)
    x["Intercept"] = 1.0
    x["axis_log2"] = num(work["axis_log2"])
    age_c = num(work["age"]) - 50.0
    x["age_centered"] = age_c
    x["age_centered_sq"] = age_c * age_c
    if interaction:
        x["female"] = work["sex"].eq("Female").astype(float)
        x["axis_log2:female"] = x["axis_log2"] * x["female"]
    for col, levels in [
        ("race", ["Mexican American", "Other Hispanic", "Non-Hispanic Black", "Other/Multi"]),
        ("smoking", ["Former", "Current"]),
    ]:
        for level in levels:
            x[f"{col}={level}"] = work[col].eq(level).astype(float)
    # Test-specific available cycles are a fixed-effect nuisance term.
    cycle_levels = sorted(work["cycle"].dropna().unique().tolist())
    for level in cycle_levels[1:]:
        x[f"cycle={level}"] = work["cycle"].eq(level).astype(float)
    x["pir"] = num(work["pir"])
    names = x.columns.tolist()
    return x.to_numpy(float), names


def fit_logistic(work: pd.DataFrame, interaction: bool) -> dict:
    required = ["outcome", "axis_log2", "age", "pir", "race", "smoking", "pooled_weight", "psu", "strata", "sex", "cycle"]
    d = work.dropna(subset=required).copy()
    d = d.loc[d["pooled_weight"].gt(0)].reset_index(drop=True)
    y = num(d["outcome"]).to_numpy(float)
    base = {"analytic_n": int(len(d)), "cases": int(y.sum()) if len(y) else 0, "controls": int(len(y) - y.sum()) if len(y) else 0}
    if len(d) == 0 or y.sum() == 0 or y.sum() == len(y):
        return {**base, "status": "not_estimable", "reason": "no complete-case outcome variation"}
    if interaction and d["sex"].nunique() != 2:
        return {**base, "status": "not_estimable", "reason": "both sexes not represented"}
    x, names = build_design(d, interaction)
    weights = num(d["pooled_weight"]).to_numpy(float)
    weights = weights / np.nanmean(weights)
    beta = np.zeros(x.shape[1])

    def ll(b):
        p = expit(np.clip(x @ b, -35, 35))
        return float(np.sum(weights * (y * np.log(p + 1e-12) + (1 - y) * np.log1p(-p + 1e-12))))

    current = ll(beta)
    converged = False
    for _ in range(200):
        p = expit(np.clip(x @ beta, -35, 35))
        gradient = x.T @ (weights * (y - p))
        hessian = x.T @ ((weights * p * (1 - p))[:, None] * x)
        step = np.linalg.pinv(hessian) @ gradient
        if not np.all(np.isfinite(step)):
            return {**base, "status": "fit_failed", "reason": "non-finite IRLS step"}
        step = np.clip(step, -5, 5)
        alpha = 1.0
        while alpha >= 1e-8 and ll(beta + alpha * step) < current - 1e-10:
            alpha /= 2
        if alpha < 1e-8:
            break
        beta += alpha * step
        current = ll(beta)
        if np.max(np.abs(alpha * step)) < 1e-8:
            converged = True
            break
    p = expit(np.clip(x @ beta, -35, 35))
    scores = (weights * (y - p))[:, None] * x
    bread = x.T @ ((weights * p * (1 - p))[:, None] * x)
    bread_inv = np.linalg.pinv(bread)
    meat = np.zeros_like(bread)
    for _, group in d[["strata", "psu"]].groupby("strata", sort=False):
        psus = group["psu"].unique()
        if len(psus) < 2:
            continue
        totals = np.vstack([scores[group.index[group["psu"].eq(psu)], :].sum(axis=0) for psu in psus])
        centered = totals - totals.mean(axis=0, keepdims=True)
        meat += len(psus) / (len(psus) - 1) * centered.T @ centered
    cov = bread_inv @ meat @ bread_inv
    cov = (cov + cov.T) / 2
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    term = "axis_log2:female" if interaction else "axis_log2"
    idx = names.index(term)
    if not np.isfinite(se[idx]) or se[idx] <= 0:
        return {**base, "status": "fit_failed", "reason": "coefficient variance unavailable"}
    design_df = max(int(d["psu"].nunique() - d["strata"].nunique()), 1)
    crit = float(t.ppf(0.975, design_df))
    stat = float(beta[idx] / se[idx])
    return {**base, "status": "ok" if converged else "converged_with_warning", "reason": "",
            "beta": float(beta[idx]), "se": float(se[idx]), "ci_low": float(beta[idx] - crit * se[idx]),
            "ci_high": float(beta[idx] + crit * se[idx]), "p_value": float(2 * t.sf(abs(stat), design_df)),
            "design_df": design_df, "psu_n": int(d["psu"].nunique()), "strata_n": int(d["strata"].nunique())}


def fixed_bh(p: pd.Series) -> pd.Series:
    values = num(p).fillna(1.0).clip(0, 1).to_numpy(float)
    order = np.argsort(values)
    ranked = values[order] * N_INTERACTIONS / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(values))
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=p.index)


def outcome_frame(outcome_id: str, cycles: list[str], inventory, model) -> tuple[pd.DataFrame, dict[str, str]]:
    spec = next(row for row in inventory.REGISTRY if row[0] == outcome_id)
    frames, files = [], {}
    for cycle in cycles:
        index = inventory.CYCLES.index(cycle)
        demo_path, smq_path = DATA / f"{cycle}_DEMO.XPT", DATA / f"{cycle}_SMQ.XPT"
        demo, smq = model.read_xpt(demo_path), model.read_xpt(smq_path)
        core = model.derive_demo(demo, index).merge(model.derive_smoking(smq), on="SEQN", how="left", validate="one_to_one")
        status, source, reason = inventory.derive_status(spec[5], DATA, cycle)
        if status is None:
            # Availability was frozen outcome-by-cycle before this analysis.
            # For example, thyroid disease begins after 1999-2002.  A missing
            # outcome cycle is excluded from that outcome's fixed source set,
            # never replaced by a proxy or inferred from a different module.
            continue
        d = core.merge(status, on="SEQN", how="inner", validate="one_to_one")
        d = d.loc[d["age"].ge(20) & (d["case"] | d["control"])].copy()
        d["outcome"] = d["case"].astype(int)
        d["cycle"] = cycle
        frames.append(d[["SEQN", "cycle", "outcome", "age", "sex", "race", "pir", "smoking", "psu", "strata"]])
        files[f"{cycle}_DEMO.XPT"] = digest(demo_path)
        files[f"{cycle}_SMQ.XPT"] = digest(smq_path)
        for filename in source.split(";"):
            p = DATA / filename
            if p.exists():
                files[filename] = digest(p)
    return pd.concat(frames, ignore_index=True), files


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    model = load(ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py", "nhanes_design_model")
    reader = load(FRAMEWORK / "step05_crc_screen" / "run_step05_crc_screen.py", "frozen_exposure_reader")
    inventory = load(OUTCOME_DIR / "run_outcome_sex_audit.py", "frozen_outcome_inventory")
    model.DATA_DIR = DATA
    tests = pd.read_csv(TESTS, dtype=str, keep_default_na=False)
    registry = pd.read_csv(REGISTRY, low_memory=False)
    frozen = pd.read_csv(OUTCOME_DIR / "frozen_outcome_set_v1.csv")
    ontology = pd.read_csv(PLAN_DIR / "frozen_outcome_organ_system_ontology_v1.csv")
    outcomes = frozen.loc[frozen["selection_for_followup"].eq(True), "outcome_id"].tolist()
    if len(tests) != 29 or tests["test_id"].nunique() != 29 or len(outcomes) != 14 or len(ontology) != 14:
        raise ValueError("Frozen 29-test / 14-outcome / ontology lock verification failed")
    if len(tests) * len(outcomes) != N_INTERACTIONS:
        raise ValueError("Primary grid is not 406")
    exposure_cache = {}
    for _, row in tests.sort_values("test_id").iterrows():
        exposure_cache[str(row.test_id)] = reader.read_test_exposure(row, registry)
    exposure_input_files = {}
    for _, source in exposure_cache.values():
        for source_row in source.get("source_rows", []):
            path = Path(source_row["local_xpt"])
            if path.exists():
                exposure_input_files[path.name] = digest(path)
    outcome_cache, input_files = {}, {}
    needed_cycles = sorted({cycle for exposure, source in exposure_cache.values() for cycle in source.get("cycles", [])})
    for outcome_id in outcomes:
        outcome_cache[outcome_id], hashes = outcome_frame(outcome_id, needed_cycles, inventory, model)
        input_files.update(hashes)
    results, qc = [], []
    for _, test in tests.sort_values("test_id").iterrows():
        test_id = str(test.test_id)
        exposure, source = exposure_cache[test_id]
        for outcome_id in outcomes:
            ontology_row = ontology.loc[ontology["outcome_id"].eq(outcome_id)].iloc[0]
            base = {"test_id": test_id, "biomarker": str(test.biomarker), "variable": str(test.variable), "matrix": str(test.matrix),
                    "exposure_axis": str(test.exposure_axes), "outcome_id": outcome_id, "outcome_name": str(ontology_row.outcome_name),
                    "organ_system_id": str(ontology_row.organ_system_id), "organ_system_name": str(ontology_row.organ_system_name),
                    "within_system_weight": float(ontology_row.within_system_weight), "cycles": ";".join(source.get("cycles", []))}
            if exposure.empty:
                pooled = {"status": "not_estimable", "reason": source.get("reason", "empty exposure"), "analytic_n": 0, "cases": 0, "controls": 0}
                male = female = pooled.copy()
            else:
                outcome = outcome_cache[outcome_id].loc[lambda x: x["cycle"].isin(source["cycles"])]
                merged = exposure.merge(outcome, on=["SEQN", "cycle"], how="inner", validate="one_to_one")
                pooled = fit_logistic(merged, interaction=True)
                male = fit_logistic(merged.loc[merged["sex"].eq("Male")], interaction=False)
                female = fit_logistic(merged.loc[merged["sex"].eq("Female")], interaction=False)
            estimable = pooled.get("status") in {"ok", "converged_with_warning"}
            result = {**base, "estimable": estimable, "fit_status": pooled.get("status"), "fit_reason": pooled.get("reason", ""),
                      "beta_interaction": pooled.get("beta", np.nan), "se_interaction": pooled.get("se", np.nan),
                      "ci_low_interaction": pooled.get("ci_low", np.nan), "ci_high_interaction": pooled.get("ci_high", np.nan), "p_interaction": pooled.get("p_value", np.nan),
                      "pooled_n": pooled.get("analytic_n", 0), "pooled_cases": pooled.get("cases", 0), "pooled_controls": pooled.get("controls", 0),
                      "male_beta": male.get("beta", np.nan), "male_se": male.get("se", np.nan), "male_n": male.get("analytic_n", 0), "male_cases": male.get("cases", 0),
                      "female_beta": female.get("beta", np.nan), "female_se": female.get("se", np.nan), "female_n": female.get("analytic_n", 0), "female_cases": female.get("cases", 0),
                      "model_specification": "survey-weighted logit: exposure + female + exposure:female + age_centered + age_centered_sq + race + PIR + smoking + cycle fixed effects"}
            results.append(result)
            qc.append({**base, "estimable": estimable, "fit_status": pooled.get("status"), "reason": pooled.get("reason", ""),
                       "pooled_n": pooled.get("analytic_n", 0), "pooled_cases": pooled.get("cases", 0), "pooled_controls": pooled.get("controls", 0),
                       "male_n": male.get("analytic_n", 0), "male_cases": male.get("cases", 0), "female_n": female.get("analytic_n", 0), "female_cases": female.get("cases", 0),
                       "design_df": pooled.get("design_df", np.nan), "psu_n": pooled.get("psu_n", np.nan), "strata_n": pooled.get("strata_n", np.nan)})
    primary = pd.DataFrame(results).sort_values(["test_id", "outcome_id"]).reset_index(drop=True)
    primary["bh_q_interaction_fixed406"] = fixed_bh(primary["p_interaction"])
    qc_df = pd.DataFrame(qc).sort_values(["test_id", "outcome_id"]).reset_index(drop=True)
    summary = primary.assign(abs_d=lambda x: x["beta_interaction"].abs()).groupby(["test_id", "biomarker", "variable", "matrix", "exposure_axis"], as_index=False).agg(DLD=("abs_d", "mean"))
    system = primary.assign(abs_d=lambda x: x["beta_interaction"].abs()).groupby(["test_id", "biomarker", "variable", "organ_system_id", "organ_system_name"], as_index=False).apply(lambda x: pd.Series({"SD": np.average(x["abs_d"], weights=x["within_system_weight"]), "SL": np.average(x["beta_interaction"], weights=x["within_system_weight"])}), include_groups=False).reset_index()
    system = system.drop(columns=["level_5"], errors="ignore")
    osd = system.groupby("test_id", as_index=False)["SD"].mean().rename(columns={"SD": "OSD"})
    summary = summary.merge(osd, on="test_id", how="left").sort_values("test_id")
    primary.to_csv(args.outdir / "sex_divergence_primary_406.csv", index=False)
    qc_df.to_csv(args.outdir / "sex_divergence_estimability_qc.csv", index=False)
    summary.to_csv(args.outdir / "sex_divergence_exposure_summary.csv", index=False)
    system.sort_values(["test_id", "organ_system_id"]).to_csv(args.outdir / "sex_divergence_system_summary.csv", index=False)
    manifest = {"analysis": "LOCKED_PRIMARY_SEX_DIVERGENCE", "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "grid": {"exposures": 29, "outcomes": 14, "interactions": 406}, "bh_fixed_denominator": 406,
                "formal_estimand": "pooled survey-weighted exposure:female coefficient", "stratified_estimates": "descriptive only",
                "model_specification": "logit outcome ~ exposure + female + exposure:female + age_centered + age_centered_sq + race + PIR + smoking + cycle fixed effects",
                "weight_rule": "test-specific laboratory/subsample weight divided by number of contributing cycles; first 1999-2002 four-year weight multiplier retained by frozen reader",
                "variance": "Taylor-style stratified PSU sandwich variance with t reference using PSU minus strata degrees of freedom",
                "input_hashes": {"test_set": digest(TESTS), "outcome_set": digest(OUTCOME_DIR / "frozen_outcome_set_v1.csv"), "ontology": digest(PLAN_DIR / "frozen_outcome_organ_system_ontology_v1.csv"), "sap_v1_1": digest(PLAN_DIR / "SEX_DIVERGENCE_STATISTICAL_ANALYSIS_PLAN_v1.1.md"), "analysis_lock_v1_1": digest(PLAN_DIR / "SEX_DIVERGENCE_ANALYSIS_LOCK_v1.1.json"), "script": digest(Path(__file__))},
                "outcome_input_xpt_hashes": input_files, "exposure_input_xpt_hashes": exposure_input_files,
                "estimable_pairs": int(primary["estimable"].sum()), "nonestimable_pairs": int((~primary["estimable"]).sum()),
                "prohibited_not_run": ["candidate selection", "CTD", "literature search", "robustness analysis", "result interpretation"]}
    (args.outdir / "primary_analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
