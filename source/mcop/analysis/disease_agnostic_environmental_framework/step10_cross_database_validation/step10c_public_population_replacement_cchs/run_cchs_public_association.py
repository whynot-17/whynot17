"""Run the frozen CCHS 2022 public-PUMF demonstration.

The exposure/outcome pair was frozen in the preceding feasibility audit. This
script performs no candidate search and does not compare alternative outcomes.
It fits a weighted logistic model and obtains coefficient precision from the
1,000 supplied CCHS bootstrap replicate weights.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
import zipfile
from pathlib import Path

import numpy as np
from scipy.special import expit
from scipy.stats import norm


OUT = Path(__file__).resolve().parent
ARCHIVE = Path(r"D:\whynot17\public_sources\cchs_2022\2022_TXT.zip")
EXPECTED_ARCHIVE_SHA256 = "a116a01fb35cc3204a8b14e36dda70cf6722e75380f88b6d822e492765ac8c41"
N_BOOTSTRAPS = 1000
DATA_WIDTH = 385
BSW_WIDTH = 8029

FIELDS = {
    "EDDVH3": (39, 39),
    "DHHGAGE": (40, 40),
    "DHH_SEX": (41, 41),
    "HWTDGBCC": (57, 57),
    "CCC_80": (91, 91),
    "SMKDVSTY": (305, 306),
    "INCDGHH": (374, 374),
    "WTS_M": (378, 385),
}
VALID = {
    "DHHGAGE": {"2", "3", "4", "5"},
    "DHH_SEX": {"1", "2"},
    "HWTDGBCC": {"1", "2"},
    "CCC_80": {"1", "2"},
    "SMKDVSTY": {"01", "02", "03", "04", "05", "06"},
    "EDDVH3": {"1", "2", "3"},
    "INCDGHH": {"1", "2", "3", "4", "5"},
}


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def field(line: bytes, name: str) -> str:
    start, end = FIELDS[name]
    return line[start - 1 : end].decode("latin1").strip()


def analytic_record(line: bytes) -> tuple[np.ndarray, int, float] | None:
    if len(line) != DATA_WIDTH:
        raise RuntimeError(f"Unexpected master record width: {len(line)}")
    values = {name: field(line, name) for name in FIELDS}
    if values["DHHGAGE"] not in VALID["DHHGAGE"]:
        return None
    if values["CCC_80"] not in VALID["CCC_80"]:
        return None
    if values["SMKDVSTY"] not in VALID["SMKDVSTY"]:
        return None
    if any(values[name] not in VALID[name] for name in ["DHH_SEX", "HWTDGBCC", "EDDVH3", "INCDGHH"]):
        return None
    try:
        weight = float(values["WTS_M"])
    except ValueError:
        return None
    if not math.isfinite(weight) or weight <= 0:
        return None

    # Reference categories: age 18-34, male, normal/underweight, least
    # education, and lowest income. Exposure is current daily/occasional
    # smoker (01/02) versus not currently smoking (03-06).
    x = np.array([
        1.0,
        1.0 if values["SMKDVSTY"] in {"01", "02"} else 0.0,
        1.0 if values["DHHGAGE"] == "3" else 0.0,
        1.0 if values["DHHGAGE"] == "4" else 0.0,
        1.0 if values["DHHGAGE"] == "5" else 0.0,
        1.0 if values["DHH_SEX"] == "2" else 0.0,
        1.0 if values["HWTDGBCC"] == "2" else 0.0,
        1.0 if values["EDDVH3"] == "2" else 0.0,
        1.0 if values["EDDVH3"] == "3" else 0.0,
        1.0 if values["INCDGHH"] == "2" else 0.0,
        1.0 if values["INCDGHH"] == "3" else 0.0,
        1.0 if values["INCDGHH"] == "4" else 0.0,
        1.0 if values["INCDGHH"] == "5" else 0.0,
    ], dtype=float)
    y = 1 if values["CCC_80"] == "1" else 0
    return x, y, weight


def bootstrap_vector(raw: bytes) -> np.ndarray:
    line = raw.rstrip(b"\r\n")
    if len(line) != BSW_WIDTH:
        raise RuntimeError(f"Unexpected bootstrap record width: {len(line)}")
    # The official bsw_i.sas layout places BSW1 at character 30 and uses
    # 1,000 fixed-width 8.2 fields.
    return np.fromiter((float(line[29 + i * 8 : 37 + i * 8]) for i in range(N_BOOTSTRAPS)), dtype=np.float32, count=N_BOOTSTRAPS)


def load_analytic_data(z: zipfile.ZipFile, data_member: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    x_rows: list[np.ndarray] = []
    y_rows: list[int] = []
    weights: list[float] = []
    total_rows = 0
    malformed = 0
    with z.open(data_member) as handle:
        for raw in handle:
            total_rows += 1
            line = raw.rstrip(b"\r\n")
            if len(line) != DATA_WIDTH:
                malformed += 1
                continue
            record = analytic_record(line)
            if record is None:
                continue
            x, y, weight = record
            x_rows.append(x)
            y_rows.append(y)
            weights.append(weight)
    if not x_rows:
        raise RuntimeError("No complete analytic records")
    return np.vstack(x_rows), np.asarray(y_rows, dtype=float), np.asarray(weights, dtype=float), total_rows, malformed


def load_bootstrap_matrix(z: zipfile.ZipFile, data_member: str, bsw_member: str, n_analytic: int, path: Path) -> tuple[np.memmap, int, int]:
    matrix = np.memmap(path, dtype="float32", mode="w+", shape=(n_analytic, N_BOOTSTRAPS))
    analytic_index = 0
    total_rows = 0
    malformed = 0
    with z.open(data_member) as data_handle, z.open(bsw_member) as bsw_handle:
        for data_raw, bsw_raw in zip(data_handle, bsw_handle):
            total_rows += 1
            data_line = data_raw.rstrip(b"\r\n")
            if len(data_line) != DATA_WIDTH:
                malformed += 1
                continue
            if analytic_record(data_line) is None:
                continue
            matrix[analytic_index, :] = bootstrap_vector(bsw_raw)
            analytic_index += 1
    if analytic_index != n_analytic:
        raise RuntimeError(f"Analytic master/BSW alignment mismatch: {analytic_index} vs {n_analytic}")
    matrix.flush()
    return matrix, total_rows, malformed


def weighted_loglik(X: np.ndarray, y: np.ndarray, w: np.ndarray, beta: np.ndarray) -> float:
    eta = np.clip(X @ beta, -35, 35)
    # Stable Bernoulli log likelihood.
    return float(np.sum(w * (y * (-np.logaddexp(0.0, -eta)) + (1.0 - y) * (-np.logaddexp(0.0, eta)))))


def fit_weighted_logistic(X: np.ndarray, y: np.ndarray, weights: np.ndarray, max_iter: int = 80, tol: float = 1e-8) -> tuple[np.ndarray, bool, int]:
    mask = np.isfinite(weights) & (weights > 0)
    if mask.sum() <= X.shape[1] or np.unique(y[mask]).size < 2:
        return np.full(X.shape[1], np.nan), False, 0
    x = X[mask]
    yy = y[mask]
    w = weights[mask].astype(float)
    w /= np.mean(w)
    beta = np.zeros(x.shape[1], dtype=float)
    old_ll = weighted_loglik(x, yy, w, beta)
    for iteration in range(1, max_iter + 1):
        eta = np.clip(x @ beta, -35, 35)
        mu = expit(eta)
        variance = np.maximum(mu * (1.0 - mu), 1e-9)
        score = x.T @ (w * (yy - mu))
        hessian = x.T @ (x * (w * variance)[:, None])
        try:
            step = np.linalg.solve(hessian, score)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(hessian, score, rcond=None)[0]
        if not np.all(np.isfinite(step)):
            return np.full(x.shape[1], np.nan), False, iteration

        accepted = False
        for fraction in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125):
            trial = beta + fraction * step
            trial_ll = weighted_loglik(x, yy, w, trial)
            if trial_ll >= old_ll - 1e-10:
                beta = trial
                old_ll = trial_ll
                accepted = True
                break
        if not accepted:
            return np.full(x.shape[1], np.nan), False, iteration
        if np.max(np.abs(step)) * fraction < tol:
            return beta, True, iteration
    return beta, False, max_iter


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-replicates", type=int, default=N_BOOTSTRAPS, help="Diagnostic override; default is all 1,000 supplied replicates")
    args = parser.parse_args()
    if not 1 <= args.max_replicates <= N_BOOTSTRAPS:
        raise ValueError("--max-replicates must be between 1 and 1000")
    if sha256_path(ARCHIVE) != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError("CCHS archive SHA-256 mismatch")

    with zipfile.ZipFile(ARCHIVE) as z:
        names = z.namelist()
        data_member = next(n for n in names if n.endswith("/PUMF_MASTER_CCHS.txt"))
        bsw_member = next(n for n in names if n.endswith("/bsw.txt"))
        X, y, main_weights, master_rows, master_bad_width = load_analytic_data(z, data_member)

        temp_dir = Path(r"D:\whynot17\public_sources\cchs_2022\runtime")
        temp_dir.mkdir(parents=True, exist_ok=True)
        memmap_path = temp_dir / "cchs_2022_analytic_bsw_float32.memmap"
        bsw_matrix, bsw_rows, bsw_bad_width = load_bootstrap_matrix(z, data_member, bsw_member, X.shape[0], memmap_path)

    point_beta, point_converged, point_iter = fit_weighted_logistic(X, y, main_weights)
    if not point_converged:
        raise RuntimeError("Full-sample weighted logistic model did not converge")

    max_reps = args.max_replicates
    rep_betas = np.full((max_reps, X.shape[1]), np.nan, dtype=float)
    rep_iters = np.zeros(max_reps, dtype=int)
    converged = np.zeros(max_reps, dtype=bool)
    for i in range(max_reps):
        rep_betas[i], converged[i], rep_iters[i] = fit_weighted_logistic(X, y, np.asarray(bsw_matrix[:, i], dtype=float))
        if (i + 1) % 100 == 0:
            print(f"bootstrap replicates completed: {i + 1}/{max_reps}", flush=True)

    exposure_index = 1
    valid_reps = rep_betas[converged, exposure_index]
    if valid_reps.size < max(100, int(0.95 * max_reps)):
        raise RuntimeError(f"Too many failed bootstrap fits: {valid_reps.size}/{max_reps}")
    se = float(np.std(valid_reps, ddof=1))
    beta = float(point_beta[exposure_index])
    or_value = math.exp(beta)
    ci_low = math.exp(beta - 1.96 * se)
    ci_high = math.exp(beta + 1.96 * se)
    percentile_ci_low = float(math.exp(np.quantile(valid_reps, 0.025)))
    percentile_ci_high = float(math.exp(np.quantile(valid_reps, 0.975)))
    p_value = float(2.0 * norm.sf(abs(beta / se)))
    mean_rep_beta = float(np.mean(valid_reps))

    def weighted_mean(v: np.ndarray) -> float:
        return float(np.sum(v * main_weights) / np.sum(main_weights))

    primary_row = {
        "analysis_id": "CCHS2022_CURRENT_SMOKING_HYPERTENSION_PRIMARY",
        "dataset": "CCHS 2022 PUMF",
        "exposure": "SMKDVSTY current daily/occasional smoker (01/02) vs not currently smoking (03-06)",
        "outcome": "CCC_80 high blood pressure, yes (1) vs no (2)",
        "population": "adults, DHHGAGE groups 2-5",
        "covariates": "DHHGAGE categorical + DHH_SEX + HWTDGBCC + EDDVH3 + INCDGHH",
        "survey_weight": "WTS_M",
        "variance": "1000 supplied bootstrap replicates; empirical variance 1/(B-1) around replicate mean",
        "n_analytic": int(X.shape[0]),
        "unweighted_hypertension_cases": int(y.sum()),
        "unweighted_current_smokers": int(X[:, exposure_index].sum()),
        "weighted_hypertension_prevalence": weighted_mean(y),
        "weighted_current_smoker_prevalence": weighted_mean(X[:, exposure_index]),
        "beta_log_odds": beta,
        "se_bootstrap": se,
        "odds_ratio_current_smoker_vs_not_current": or_value,
        "ci95_low": ci_low,
        "ci95_high": ci_high,
        "bootstrap_percentile_ci95_low": percentile_ci_low,
        "bootstrap_percentile_ci95_high": percentile_ci_high,
        "p_value_normal_reference": p_value,
        "point_fit_converged": point_converged,
        "point_fit_iterations": point_iter,
        "bootstrap_replicates_requested": max_reps,
        "bootstrap_replicates_converged": int(converged.sum()),
        "bootstrap_replicate_mean_beta": mean_rep_beta,
        "association_was_used_for_candidate_selection": False,
        "candidate_was_frozen_before_modeling": True,
    }

    estimate_rows = []
    for i in range(max_reps):
        row = {
            "replicate": i + 1,
            "converged": bool(converged[i]),
            "iterations": int(rep_iters[i]),
            "beta_log_odds": None if not converged[i] else float(rep_betas[i, exposure_index]),
            "odds_ratio": None if not converged[i] else float(math.exp(rep_betas[i, exposure_index])),
        }
        estimate_rows.append(row)

    diagnostics = {
        "analysis_id": primary_row["analysis_id"],
        "master_rows": master_rows,
        "master_bad_width_rows": master_bad_width,
        "bootstrap_rows": bsw_rows,
        "bootstrap_bad_width_rows": bsw_bad_width,
        "analytic_n": int(X.shape[0]),
        "unweighted_cases": int(y.sum()),
        "unweighted_exposed": int(X[:, exposure_index].sum()),
        "point_fit_converged": point_converged,
        "bootstrap_requested": max_reps,
        "bootstrap_converged": int(converged.sum()),
        "bootstrap_failed": int((~converged).sum()),
        "bootstrap_se_ddof": 1,
        "archive_sha256": EXPECTED_ARCHIVE_SHA256,
        "association_was_used_for_selection": False,
    }

    with (OUT / "cchs_2022_smoking_hypertension_primary_result.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(primary_row))
        writer.writeheader()
        writer.writerow(primary_row)
    with (OUT / "cchs_2022_smoking_hypertension_bootstrap_estimates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(estimate_rows[0]))
        writer.writeheader()
        writer.writerows(estimate_rows)
    (OUT / "CCHS_PUBLIC_ASSOCIATION_QC_SUMMARY.json").write_text(json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8")
    (OUT / "CCHS_PUBLIC_ASSOCIATION_MANIFEST.json").write_text(json.dumps({
        "analysis_id": primary_row["analysis_id"],
        "status": "complete_frozen_demonstration",
        "files": [
            "run_cchs_public_association.py",
            "cchs_2022_smoking_hypertension_primary_result.csv",
            "cchs_2022_smoking_hypertension_bootstrap_estimates.csv",
            "CCHS_PUBLIC_ASSOCIATION_QC_SUMMARY.json",
            "CCHS_PUBLIC_ASSOCIATION_MANIFEST.json",
            "STEP10C_CCHS_PUBLIC_ASSOCIATION_REPORT.md",
        ],
        "raw_archive_excluded_from_git": True,
        "archive_sha256": EXPECTED_ARCHIVE_SHA256,
        "candidate_freeze_source": "CCHS_PUBLIC_AUDIT_QC_SUMMARY.json",
        "association_was_used_for_selection": False,
    }, indent=2) + "\n", encoding="utf-8")

    report = f"""# Step 10C — CCHS public association demonstration

## Frozen analysis

The feasibility audit was completed before this model was run. The sole frozen demonstration is:

> **Current smoking status (`SMKDVSTY`) → self-reported high blood pressure (`CCC_80`) among adults**

The model uses categorical adult age group, sex, adult BMI classification, education, and household income. Point estimation uses `WTS_M`; precision uses the 1,000 supplied CCHS bootstrap replicate weights (`BSW1–BSW1000`). The bootstrap standard error is the empirical standard deviation with `ddof=1` around the replicate mean. No CCHS-specific Fay multiplier is applied.

## Result

- Analytic N: **{X.shape[0]:,}**
- Unweighted high-blood-pressure cases: **{int(y.sum()):,}**
- Unweighted current smokers: **{int(X[:, exposure_index].sum()):,}**
- Current smoker versus not-currently-smoking OR: **{or_value:.3f}**
- 95% CI: **{ci_low:.3f}–{ci_high:.3f}**
- Bootstrap percentile 95% CI sensitivity: **{percentile_ci_low:.3f}–{percentile_ci_high:.3f}**
- Normal-reference P value: **{p_value:.6g}**
- Bootstrap fits converged: **{int(converged.sum())}/{max_reps}**

This is a source-native population-replacement demonstration, not a replication of the 29 NHANES urinary biomarker tests. It is also not a causal estimate: CCHS is cross-sectional and the exposure/outcome are contemporaneous survey measures.

## Firewall

The exposure–outcome pair was frozen in `CCHS_PUBLIC_AUDIT_QC_SUMMARY.json` before this script was run. No association estimate, P value, confidence interval, or FDR was used to choose the pair. No alternative disease was searched after seeing this result.

The raw CCHS archive and bootstrap matrix remain outside Git under `D:\\whynot17\\public_sources\\cchs_2022`; only aggregate result/provenance files are versioned.
"""
    (OUT / "STEP10C_CCHS_PUBLIC_ASSOCIATION_REPORT.md").write_text(report, encoding="utf-8")

    del bsw_matrix
    try:
        os.remove(memmap_path)
    except OSError:
        pass


if __name__ == "__main__":
    main()
