from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import f, t


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def weighted_quantile(
    values: pd.Series, weights: pd.Series, quantiles: np.ndarray
) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    w = pd.to_numeric(weights, errors="coerce").to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[valid], w[valid]
    order = np.argsort(x)
    x, w = x[order], w[order]
    positions = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(quantiles, positions, x)


def wald_f_test(
    beta: np.ndarray, covariance: np.ndarray, indices: list[int], design_df: int
) -> float:
    sub_beta = beta[indices]
    sub_cov = covariance[np.ix_(indices, indices)]
    statistic = float(sub_beta @ np.linalg.pinv(sub_cov) @ sub_beta)
    df_num = len(indices)
    return float(f.sf(statistic / df_num, df_num, design_df))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    scripts_dir = args.analysis_root / "work" / "scripts"
    mcop_module = load_module(scripts_dir / "mcop_crc_phase2.py", "mcop_phase2")
    model_module = mcop_module.load_validated_functions()

    harmonized_path = (
        args.analysis_root / "work" / "nhanes_phase2a" / "phase2b_harmonized.pkl"
    )
    data_dir = args.analysis_root / "work" / "nhanes_phase2a" / "data"
    harmonized = pd.read_pickle(harmonized_path)
    mcop, _ = mcop_module.read_mcop(data_dir)
    data = mcop_module.build_frame(harmonized, mcop)
    primary = model_module.population_frames(data)["CRC_vs_cancer_free"]

    required = [
        "outcome",
        "mcop_log2",
        "age",
        "bmi",
        "pir",
        "creatinine_log2",
        "sex",
        "race",
        "smoking",
        "pooled_weight",
        "psu",
        "strata",
    ]
    cc = primary.dropna(subset=required).copy()
    cc = cc.loc[cc["pooled_weight"] > 0].copy()

    knot_quantiles = np.array([0.05, 0.35, 0.65, 0.95])
    knots = weighted_quantile(
        cc["mcop_log2"], cc["pooled_weight"], knot_quantiles
    )
    basis = model_module.spline_basis(cc["mcop_log2"].to_numpy(float), knots)
    basis_cols = []
    for index in range(basis.shape[1]):
        name = f"mcop_spline_{index + 1}"
        cc[name] = basis[:, index]
        basis_cols.append(name)

    fit = model_module.fit_survey_logistic(
        cc,
        basis_cols + ["age", "bmi", "pir", "creatinine_log2"],
        ["sex", "race", "smoking"],
        exposure_name=basis_cols[0],
        levels=model_module.LEVELS,
    )
    if fit.get("status") not in {"ok", "converged_with_warning"}:
        raise RuntimeError(f"Spline model failed: {fit}")

    names = list(fit["coefficients"])
    beta = np.array([fit["coefficients"][name] for name in names], dtype=float)
    covariance = np.asarray(fit["covariance"], dtype=float)
    indices = [names.index(name) for name in basis_cols]

    grid_log2 = np.linspace(knots[0], knots[-1], 240)
    grid_basis = model_module.spline_basis(grid_log2, knots)
    reference_log2 = float(
        weighted_quantile(
            cc["mcop_log2"], cc["pooled_weight"], np.array([0.5])
        )[0]
    )
    reference_basis = model_module.spline_basis(
        np.array([reference_log2], dtype=float), knots
    )[0]

    contrasts = np.zeros((len(grid_log2), len(names)), dtype=float)
    contrasts[:, indices] = grid_basis - reference_basis
    log_or = contrasts @ beta
    variance = np.einsum("ij,jk,ik->i", contrasts, covariance, contrasts)
    se = np.sqrt(np.maximum(variance, 0.0))
    critical = float(t.ppf(0.975, int(fit["design_df"])))
    overall_p_f = wald_f_test(
        beta, covariance, indices, int(fit["design_df"])
    )
    nonlinear_p_f = wald_f_test(
        beta, covariance, indices[1:], int(fit["design_df"])
    )

    curve = pd.DataFrame(
        {
            "mcop_log2": grid_log2,
            "mcop_ng_ml": np.power(2.0, grid_log2),
            "log_or_vs_median": log_or,
            "se_log_or": se,
            "or_vs_median": np.exp(log_or),
            "ci_low": np.exp(log_or - critical * se),
            "ci_high": np.exp(log_or + critical * se),
        }
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    curve.to_csv(args.output_csv, index=False)

    metadata = {
        "display_range": "survey-weighted 5th to 95th percentile of complete-case log2(MCOP)",
        "reference": "survey-weighted complete-case median log2(MCOP)",
        "reference_log2": reference_log2,
        "reference_ng_ml": float(2.0**reference_log2),
        "knots_log2": knots.tolist(),
        "knots_ng_ml": np.power(2.0, knots).tolist(),
        "knot_quantiles": knot_quantiles.tolist(),
        "N": int(fit["N"]),
        "CRC_N": int(fit["CRC_N"]),
        "design_df": int(fit["design_df"]),
        "critical_t_0.975": critical,
        "overall_P_F": overall_p_f,
        "nonlinear_P_F": nonlinear_p_f,
        "model_status": fit["status"],
    }
    args.output_json.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(metadata, ensure_ascii=False))


if __name__ == "__main__":
    main()
