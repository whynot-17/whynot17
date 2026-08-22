"""Planning power simulation for the WHI MCOP-CRC nested case-control study.

This is a planning simulation, not a WHI analysis. It uses the observed pooled
NHANES complete-case log2(MCOP) distribution and covariate composition as a
proxy because WHI biospecimen access is not yet confirmed. Each simulated
matched set contains one incident CRC case and either one or two controls;
case status is generated from a conditional-logistic data-generating model.
The reported power is two-sided Wald power at alpha=0.05.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_HARMONIZED = ROOT / "work" / "nhanes_phase2a" / "phase2b_harmonized.pkl"
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_OUTPUT = ROOT / "outputs"
CASE_COUNTS = [181, 150, 120, 100, 80]
MATCH_RATIOS = [1, 2]
TARGET_ORS = [1.15, 1.20, 1.25, 1.30]
NUMERIC_PROXY = ["age", "bmi", "pir", "creatinine_log2"]
CATEGORICAL_PROXY = ["sex", "race", "smoking"]


def load_mcop_module():
    path = ROOT / "work" / "scripts" / "mcop_crc_phase2.py"
    spec = importlib.util.spec_from_file_location("mcop_phase2_for_power", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load MCOP implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_proxy_frame(harmonized: Path, data_dir: Path) -> pd.DataFrame:
    mcop_module = load_mcop_module()
    model_module = mcop_module.load_validated_functions()
    data = mcop_module.build_frame(
        pd.read_pickle(harmonized),
        mcop_module.read_mcop(data_dir)[0],
    )
    primary = model_module.population_frames(data)["CRC_vs_cancer_free"].copy()
    required = [
        "mcop_log2", "age", "bmi", "pir", "creatinine_log2",
        "sex", "race", "smoking", "pooled_weight", "psu", "strata",
    ]
    return primary.dropna(subset=required).reset_index(drop=True)


def build_proxy_design(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, dict]:
    exposure = frame["mcop_log2"].to_numpy(dtype=float)
    exposure_center = float(exposure.mean())
    exposure = exposure - exposure_center

    numeric_parts = []
    numeric_meta = {}
    for col in NUMERIC_PROXY:
        values = frame[col].to_numpy(dtype=float)
        mean = float(values.mean())
        sd = float(values.std(ddof=0))
        if sd <= 0:
            raise ValueError(f"Proxy covariate has zero variance: {col}")
        numeric_parts.append((values - mean) / sd)
        numeric_meta[col] = {"mean": mean, "sd": sd}

    categorical = pd.get_dummies(
        frame[CATEGORICAL_PROXY].astype("string"),
        drop_first=True,
        dtype=float,
    )
    z = np.column_stack(numeric_parts + [categorical.to_numpy(dtype=float)])
    metadata = {
        "exposure_center_log2": exposure_center,
        "exposure_sd_log2": float(frame["mcop_log2"].std(ddof=1)),
        "numeric_standardization": numeric_meta,
        "categorical_dummy_columns": list(categorical.columns),
        "n_proxy_rows": int(len(frame)),
    }
    return exposure, z, metadata


def conditional_fit(x: np.ndarray, z: np.ndarray, case_position: np.ndarray) -> dict:
    """Fit conditional logistic regression for equal-size matched sets."""
    n_sets, n_members = x.shape
    design = np.concatenate([x[:, :, None], z], axis=2)
    n_parameters = design.shape[2]
    case_x = design[np.arange(n_sets), case_position]
    beta = np.zeros(n_parameters, dtype=float)

    def state(current: np.ndarray):
        eta = np.einsum("skp,p->sk", design, current)
        log_denominator = logsumexp(eta, axis=1)
        probabilities = np.exp(eta - log_denominator[:, None])
        mean_x = np.einsum("sk,skp->sp", probabilities, design)
        gradient = (case_x - mean_x).sum(axis=0)
        second = np.einsum("sk,skp,skq->spq", probabilities, design, design)
        hessian = -(second - np.einsum("sp,sq->spq", mean_x, mean_x)).sum(axis=0)
        log_likelihood = float(
            (np.einsum("skp,p->sk", case_x[:, None, :], current)[:, 0] - log_denominator).sum()
        )
        return log_likelihood, gradient, hessian

    converged = False
    final_hessian = None
    for _ in range(80):
        log_likelihood, gradient, hessian = state(beta)
        final_hessian = hessian
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

    if final_hessian is None:
        return {"converged": False, "beta": np.nan, "se": np.nan, "p": np.nan, "or": np.nan}
    _, gradient, final_hessian = state(beta)
    try:
        covariance = np.linalg.inv(-final_hessian)
    except np.linalg.LinAlgError:
        covariance = np.linalg.pinv(-final_hessian)
    variance = float(covariance[0, 0]) if covariance.size else np.nan
    se = float(np.sqrt(variance)) if np.isfinite(variance) and variance > 0 else np.nan
    beta_exposure = float(beta[0])
    z_score = beta_exposure / se if np.isfinite(se) and se > 0 else np.nan
    p_value = float(2.0 * norm.sf(abs(z_score))) if np.isfinite(z_score) else np.nan
    if np.isfinite(beta_exposure) and beta_exposure < 700:
        odds_ratio = float(np.exp(beta_exposure))
    else:
        odds_ratio = np.nan
    return {
        "converged": bool(converged and np.isfinite(se)),
        "beta": beta_exposure,
        "se": se,
        "p": p_value,
        "or": odds_ratio,
        "gradient_max_abs": float(np.max(np.abs(gradient))),
    }


def simulate_cell(
    x_pool: np.ndarray,
    z_pool: np.ndarray,
    n_cases: int,
    controls_per_case: int,
    target_or: float,
    n_sim: int,
    rng: np.random.Generator,
) -> dict:
    members = controls_per_case + 1
    beta_true = float(np.log(target_or))
    records = []
    for _ in range(n_sim):
        indices = rng.integers(0, len(x_pool), size=(n_cases, members))
        x = x_pool[indices]
        z = z_pool[indices]
        eta = beta_true * x
        probabilities = np.exp(eta - logsumexp(eta, axis=1)[:, None])
        random_values = rng.random(n_cases)
        case_position = (np.cumsum(probabilities, axis=1) < random_values[:, None]).sum(axis=1)
        fit = conditional_fit(x, z, case_position)
        records.append(fit)

    results = pd.DataFrame(records)
    converged = results[results["converged"]].copy()
    power = float((converged["p"] < 0.05).mean()) if not converged.empty else np.nan
    return {
        "n_cases": n_cases,
        "controls_per_case": controls_per_case,
        "total_n": int(n_cases * members),
        "target_or_per_doubling": target_or,
        "n_sim": n_sim,
        "n_converged": int(len(converged)),
        "convergence_rate": float(len(converged) / n_sim),
        "power_two_sided_alpha_0.05": power,
        "median_estimated_or": float(converged["or"].median()) if not converged.empty else np.nan,
        "q025_estimated_or": float(converged["or"].quantile(0.025)) if not converged.empty else np.nan,
        "q975_estimated_or": float(converged["or"].quantile(0.975)) if not converged.empty else np.nan,
        "median_se_log_or": float(converged["se"].median()) if not converged.empty else np.nan,
    }


def write_report(outdir: Path, rows: pd.DataFrame, proxy_meta: dict, n_sim: int, seed: int) -> None:
    lines = [
        "# MCOP–CRC Phase 3A：WHI power planning simulation",
        "",
        "本文件是 WHI biospecimen access 确认前的规划模拟，不是 WHI 实际回归结果。模拟使用 NHANES Phase 2 complete-case 的 MCOP/covariate 分布作为代理，模拟一个病例对应 1 或 2 个 matched controls，并用 conditional logistic likelihood 生成和拟合病例状态。",
        "",
        f"- Simulation seed: `{seed}`",
        f"- Replicates per cell: `{n_sim}`",
        "- Alpha: 0.05, two-sided Wald test",
        f"- NHANES proxy rows: `{proxy_meta['n_proxy_rows']}`",
        f"- Proxy log2(MCOP) SD: `{proxy_meta['exposure_sd_log2']:.6g}`",
        "- Target effect is OR per MCOP doubling",
        "",
        "## Interpretation",
        "",
        "Power is scenario planning only. It will change with the actual WHI MCOP distribution, residual urine availability, matching factors, missingness, assay batch structure and the number of analyzable CRC cases. The table should be used to decide whether the WHI sample can support directionally informative replication, not to claim that WHI is already analyzable.",
        "",
        "## Results",
        "",
        "| CRC cases | Controls/case | Total N | Target OR | Convergence | Power | Median estimated OR | 2.5–97.5% estimated OR |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, row in rows.iterrows():
        lines.append(
            f"| {int(row['n_cases'])} | {int(row['controls_per_case'])} | {int(row['total_n'])} | {row['target_or_per_doubling']:.2f} | {row['convergence_rate']:.3f} | {row['power_two_sided_alpha_0.05']:.3f} | {row['median_estimated_or']:.3f} | {row['q025_estimated_or']:.3f}–{row['q975_estimated_or']:.3f} |"
        )
    lines += [
        "",
        "## Limitations frozen before WHI access",
        "",
        "1. NHANES is an exposure-distribution proxy, not a substitute for WHI biomarker data.",
        "2. The simulation sets nuisance outcome coefficients to zero; covariates are retained in the fitted model and their empirical correlation with MCOP is preserved through the NHANES proxy rows.",
        "3. Matching is represented by conditional-logistic matched sets, but exact WHI matching-factor distributions are not yet known.",
        "4. This is not a sample-size justification until WHI confirms analyzable urine samples and actual MCOP assay performance.",
        "",
        "## Next gate",
        "",
        "Confirm the number of CRC cases with sufficient prediagnostic urine volume, matched controls, sample timing and assay availability in WHI Query Builder/CCC. Then rerun this simulation with the observed WHI exposure variance and missingness assumptions before finalizing the sample request.",
    ]
    (outdir / "mcop_phase3a_power_simulation_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", type=Path, default=DEFAULT_HARMONIZED)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--n-sim", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260822)
    args = parser.parse_args()
    if args.n_sim < 100:
        raise ValueError("Use at least 100 simulations per cell for planning output")
    args.outdir.mkdir(parents=True, exist_ok=True)

    proxy_frame = load_proxy_frame(args.harmonized, args.data_dir)
    x_pool, z_pool, proxy_meta = build_proxy_design(proxy_frame)
    rng = np.random.default_rng(args.seed)
    rows = []
    total_cells = len(CASE_COUNTS) * len(MATCH_RATIOS) * len(TARGET_ORS)
    cell = 0
    for n_cases in CASE_COUNTS:
        for controls_per_case in MATCH_RATIOS:
            for target_or in TARGET_ORS:
                cell += 1
                print(f"[{cell}/{total_cells}] cases={n_cases} controls_per_case={controls_per_case} target_or={target_or}")
                rows.append(
                    simulate_cell(
                        x_pool,
                        z_pool,
                        n_cases,
                        controls_per_case,
                        target_or,
                        args.n_sim,
                        rng,
                    )
                )
    result = pd.DataFrame(rows)
    result.to_csv(args.outdir / "mcop_phase3a_power_simulation.csv", index=False)
    write_report(args.outdir, result, proxy_meta, args.n_sim, args.seed)
    manifest = {
        "analysis": "MCOP-CRC Phase 3A planning power simulation",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "seed": args.seed,
        "n_sim_per_cell": args.n_sim,
        "case_counts": CASE_COUNTS,
        "controls_per_case": MATCH_RATIOS,
        "target_or_per_doubling": TARGET_ORS,
        "alpha": 0.05,
        "test": "two-sided Wald test from custom conditional logistic likelihood",
        "proxy": proxy_meta,
        "outputs": [
            "mcop_phase3a_power_simulation.csv",
            "mcop_phase3a_power_simulation_report.md",
        ],
    }
    (args.outdir / "mcop_phase3a_power_simulation_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
