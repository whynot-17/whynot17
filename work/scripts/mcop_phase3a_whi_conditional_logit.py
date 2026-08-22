"""Pre-specified conditional-logistic pipeline for WHI MCOP-CRC data.

This script is intentionally data-access agnostic. It validates a de-identified
WHI nested case-control CSV after biospecimen access is granted, then fits the
primary model and pre-specified lag/repeated-exposure analyses. It does not
download or infer WHI data.

Expected input columns are documented in outputs/mcop_phase3a_whi_preregistration.md.
Optional columns `mcop_detected` and `mcop_lloq_ng_mL` allow the frozen
LLOQ/sqrt(2) substitution to be applied inside the pipeline; otherwise
`mcop_ng_mL` must already contain the pre-specified non-detect replacement.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "outputs"
REQUIRED_COLUMNS = [
    "set_id", "crc", "mcop_ng_mL", "creatinine_mg_dL", "age", "bmi",
    "smoking", "alcohol", "physical_activity", "ses", "sex", "race",
    "assay_batch",
]
NUMERIC_COVARIATES = [
    "age", "bmi", "alcohol", "physical_activity", "ses", "log2_creatinine",
]
CATEGORICAL_COVARIATES = ["smoking", "sex", "race", "assay_batch"]


def validate_input(data: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    if not set(pd.to_numeric(data["crc"], errors="coerce").dropna().unique()).issubset({0, 1}):
        raise ValueError("crc must be coded 0/1")
    counts = data.groupby("set_id")["crc"].agg(["sum", "count"])
    if not (counts["sum"].eq(1) & counts["count"].ge(2)).all():
        raise ValueError("Each set_id must contain exactly one case and at least one control")


def prepare(data: pd.DataFrame) -> pd.DataFrame:
    work = data.copy()
    work["mcop_ng_mL"] = pd.to_numeric(work["mcop_ng_mL"], errors="coerce")
    work["creatinine_mg_dL"] = pd.to_numeric(work["creatinine_mg_dL"], errors="coerce")
    if {"mcop_detected", "mcop_lloq_ng_mL"}.issubset(work.columns):
        detected = work["mcop_detected"].astype("boolean")
        lloq = pd.to_numeric(work["mcop_lloq_ng_mL"], errors="coerce")
        work.loc[detected.eq(False), "mcop_ng_mL"] = lloq.loc[detected.eq(False)] / np.sqrt(2.0)
    work["log2_mcop"] = np.log2(work["mcop_ng_mL"].where(work["mcop_ng_mL"] > 0))
    work["log2_creatinine"] = np.log2(work["creatinine_mg_dL"].where(work["creatinine_mg_dL"] > 0))
    return work


def encode_design(frame: pd.DataFrame, exposure_column: str = "log2_mcop") -> tuple[pd.DataFrame, list[str]]:
    numeric = frame[NUMERIC_COVARIATES].copy()
    for column in numeric.columns:
        numeric[column] = pd.to_numeric(numeric[column], errors="coerce")
        sd = numeric[column].std(ddof=0)
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"Covariate has no usable variance: {column}")
        numeric[column] = (numeric[column] - numeric[column].mean()) / sd
    categorical = pd.get_dummies(
        frame[CATEGORICAL_COVARIATES].astype("string"),
        drop_first=True,
        dtype=float,
    )
    design = pd.concat(
        [pd.to_numeric(frame[exposure_column], errors="coerce").rename(exposure_column), numeric, categorical],
        axis=1,
    )
    return design.astype(float), list(design.columns)


def fit_conditional_logit(frame: pd.DataFrame, exposure_column: str = "log2_mcop") -> dict:
    design_frame, columns = encode_design(frame, exposure_column)
    work = frame[["set_id", "crc"]].join(design_frame)
    work = work.dropna().copy()
    within_set_range = work.groupby("set_id", sort=False)[columns].agg(lambda values: values.max() - values.min())
    estimable = within_set_range.gt(1e-12).any(axis=0)
    if not bool(estimable.get(exposure_column, False)):
        raise ValueError(f"Exposure has no within-set variation: {exposure_column}")
    dropped = [column for column in columns if column != exposure_column and not bool(estimable[column])]
    columns = [exposure_column] + [column for column in columns if column != exposure_column and bool(estimable[column])]
    work = work[["set_id", "crc"] + columns]
    grouped = list(work.groupby("set_id", sort=False))
    x_blocks = []
    case_positions = []
    for _, group in grouped:
        if int(group["crc"].sum()) != 1:
            raise ValueError("Each analysis set must contain one case")
        x_blocks.append(group[columns].to_numpy(dtype=float))
        case_positions.append(int(np.flatnonzero(group["crc"].to_numpy(dtype=int))[0]))
    sizes = {block.shape[0] for block in x_blocks}
    if len(sizes) != 1:
        raise ValueError("All sets must have the same size for this implementation")
    design = np.stack(x_blocks, axis=0)
    case_positions = np.asarray(case_positions, dtype=int)
    n_sets, n_members, n_parameters = design.shape
    case_x = design[np.arange(n_sets), case_positions]
    beta = np.zeros(n_parameters, dtype=float)

    def state(current: np.ndarray):
        eta = np.einsum("skp,p->sk", design, current)
        log_denominator = logsumexp(eta, axis=1)
        probabilities = np.exp(eta - log_denominator[:, None])
        mean_x = np.einsum("sk,skp->sp", probabilities, design)
        gradient = (case_x - mean_x).sum(axis=0)
        second = np.einsum("sk,skp,skq->spq", probabilities, design, design)
        hessian = -(second - np.einsum("sp,sq->spq", mean_x, mean_x)).sum(axis=0)
        log_likelihood = float((case_x @ current - log_denominator).sum())
        return log_likelihood, gradient, hessian

    converged = False
    for _ in range(100):
        log_likelihood, gradient, hessian = state(beta)
        if np.max(np.abs(gradient)) < 1e-7:
            converged = True
            break
        try:
            direction = -np.linalg.solve(hessian, gradient)
        except np.linalg.LinAlgError:
            direction = -np.linalg.pinv(hessian) @ gradient
        if not np.all(np.isfinite(direction)):
            break
        step = 1.0
        improved = False
        while step >= 1e-8:
            candidate = beta + step * direction
            candidate_ll = state(candidate)[0]
            if np.isfinite(candidate_ll) and candidate_ll >= log_likelihood - 1e-10:
                beta = candidate
                improved = True
                break
            step *= 0.5
        if not improved:
            break

    _, gradient, hessian = state(beta)
    try:
        covariance = np.linalg.inv(-hessian)
    except np.linalg.LinAlgError:
        covariance = np.linalg.pinv(-hessian)
    variance = float(covariance[0, 0])
    se = float(np.sqrt(variance)) if np.isfinite(variance) and variance > 0 else np.nan
    beta_exposure = float(beta[0])
    z_score = beta_exposure / se if np.isfinite(se) and se > 0 else np.nan
    p_value = float(2.0 * norm.sf(abs(z_score))) if np.isfinite(z_score) else np.nan
    ci_low = float(np.exp(beta_exposure - 1.96 * se)) if np.isfinite(se) else np.nan
    ci_high = float(np.exp(beta_exposure + 1.96 * se)) if np.isfinite(se) else np.nan
    return {
        "status": "converged" if converged and np.isfinite(se) else "converged_with_warning",
        "N": int(len(work)),
        "matched_sets": int(n_sets),
        "members_per_set": int(n_members),
        "CRC_cases": int(work["crc"].sum()),
        "exposure": exposure_column,
        "beta": beta_exposure,
        "SE": se,
        "OR_per_doubling": float(np.exp(beta_exposure)),
        "CI_low": ci_low,
        "CI_high": ci_high,
        "P": p_value,
        "gradient_max_abs": float(np.max(np.abs(gradient))),
        "covariate_columns": ";".join(columns[1:]),
        "non_estimable_covariates_dropped": ";".join(dropped),
    }


def analysis_frame(
    work: pd.DataFrame,
    lag_min_years: float | None = None,
    exposure_column: str = "log2_mcop",
) -> pd.DataFrame:
    frame = work.copy()
    if lag_min_years is not None:
        if "lag_years" not in frame.columns:
            raise ValueError("lag_years is required for lagged analyses")
        frame = frame[pd.to_numeric(frame["lag_years"], errors="coerce") >= lag_min_years].copy()
    required = REQUIRED_COLUMNS + [exposure_column, "log2_creatinine"]
    return frame.dropna(subset=required).copy()


def run_pipeline(input_path: Path, outdir: Path) -> pd.DataFrame:
    data = pd.read_csv(input_path)
    validate_input(data)
    work = prepare(data)
    outputs = []
    for label, lag in [("primary", None), ("lag_ge_2y", 2.0), ("lag_ge_5y", 5.0)]:
        frame = analysis_frame(work, lag)
        fit = fit_conditional_logit(frame)
        fit["analysis"] = label
        fit["lag_min_years"] = lag
        outputs.append(fit)

    if {"mcop_baseline_ng_mL", "mcop_year3_ng_mL"}.issubset(work.columns):
        baseline = pd.to_numeric(work["mcop_baseline_ng_mL"], errors="coerce")
        year3 = pd.to_numeric(work["mcop_year3_ng_mL"], errors="coerce")
        work["log2_mcop_repeated_mean"] = (
            np.log2(baseline.where(baseline > 0))
            + np.log2(year3.where(year3 > 0))
        ) / 2.0
        repeated = analysis_frame(work, exposure_column="log2_mcop_repeated_mean")
        fit = fit_conditional_logit(repeated, exposure_column="log2_mcop_repeated_mean")
        fit["analysis"] = "repeated_urine_mean"
        fit["lag_min_years"] = np.nan
        outputs.append(fit)

    result = pd.DataFrame(outputs)
    result.to_csv(outdir / "mcop_phase3a_whi_model_results.csv", index=False)
    manifest = {
        "analysis": "WHI MCOP-CRC Phase 3A conditional logistic pipeline",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path),
        "status": "requires WHI-accessed de-identified input CSV",
        "primary_exposure": "log2(MCOP), OR per doubling",
        "outcome": "incident invasive CRC, one case per matched set",
        "covariates": NUMERIC_COVARIATES + CATEGORICAL_COVARIATES,
        "outputs": ["mcop_phase3a_whi_model_results.csv"],
    }
    (outdir / "mcop_phase3a_whi_model_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    print(run_pipeline(args.input, args.outdir).to_string(index=False))


if __name__ == "__main__":
    main()
