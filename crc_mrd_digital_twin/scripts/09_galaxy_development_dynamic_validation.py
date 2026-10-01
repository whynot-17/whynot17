from __future__ import annotations

"""GALAXY-led CRC recurrence model, fixed-horizon update, and public validation.

Inputs remain in the project data_raw directory. The script writes aggregate
outputs to --output-root; patient-level predictions are written only under
outputs/local_only and must not be published.
"""

import argparse
import csv
import importlib.util
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np
import pdfplumber
from docx import Document
from openpyxl import load_workbook
from pypdf import PdfReader

DAY_PER_MONTH = 365.2425 / 12.0
ORIGIN_DAY = 200.0  # End of GALAXY's published 160-200 day 6-month window.
HORIZON_DAY = 365.2425
L3_LOW, L3_HIGH, L3_TARGET = 70.0, 112.0, 91.3
MRD0_LOW, MRD0_HIGH, MRD0_TARGET = 14.0, 70.0, 28.0
L6_LOW, L6_HIGH, L6_TARGET = 160.0, 200.0, 182.6
L6_STRICT_DAY = 6.0 * DAY_PER_MONTH
BOOT_REPS = 1000
SEED = 20260930


def clean(x):
    return str(x or "").strip()


def num(x):
    try:
        if x is None or clean(x).lower() in {"", "na", "n/a", "none"}:
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def status01(x):
    s = clean(x).lower()
    if s in {"positive", "detected", "1"}:
        return 1
    if s in {"negative", "not detected", "0"}:
        return 0
    if isinstance(x, (int, float)) and x in (0, 1):
        return int(x)
    return None


def sex01(x):
    s = clean(x).lower()
    if s in {"male", "m", "man"}:
        return 1
    if s in {"female", "f", "woman"}:
        return 0
    return None


def load_legacy(root: Path):
    path = root / "scripts" / "03_multicohort_landmark_model.py"
    spec = importlib.util.spec_from_file_location("legacy_crc_parser_for_galaxy_model", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def enrich_chen(root: Path, mod):
    patients = mod.read_chen()
    table = Document(root / "data_raw" / "13045_2021_1089_MOESM1_ESM.docx").tables[0]
    headers = [c.text.strip() for c in table.rows[1].cells]
    i_id, i_age, i_sex = headers.index("Patient No."), headers.index("Age"), headers.index("Sex")
    baseline = {}
    for row in table.rows[2:]:
        vals = [c.text.strip() for c in row.cells]
        if vals and vals[i_id]:
            baseline[vals[i_id]] = {"age": num(vals[i_age]), "sex": vals[i_sex]}
    for p in patients:
        p.update(baseline.get(str(p["patient_id"]), {}))
    return patients


def enrich_coloniaiq(root: Path, mod):
    patients = mod.read_coloniaiq()
    reader = PdfReader(root / "data_raw" / "mo_colonaiq_supplement1.pdf")
    baseline = {}
    for page in reader.pages[5:19]:
        for line in (page.extract_text() or "").splitlines():
            tok = line.strip().split()
            if len(tok) >= 6 and re.fullmatch(r"P\d+", tok[0]) and tok[1] in ("Male", "Female"):
                baseline[tok[0]] = {"sex": tok[1], "age": num(tok[2])}
    for p in patients:
        p.update(baseline.get(str(p["patient_id"]), {}))
    return patients


def enrich_cosmos(root: Path, mod):
    patients = mod.read_cosmos()
    ws = load_workbook(root / "data_raw" / "cosmos_crc01_supp_table2.xlsx",
                       read_only=True, data_only=True).active
    baseline = {}
    for r in ws.iter_rows(min_row=4, max_col=23, values_only=True):
        if not r[0]:
            continue
        baseline[str(r[0])] = {
            "age": num(r[1]), "sex": r[2], "site": r[3], "stage": r[4],
            "rfs_months": num(r[13]),
        }
    for p in patients:
        b = baseline.get(str(p["patient_id"]), {})
        p.update(b)
        # Column 14 is RFS from surgery; for non-recurrence patients it is their
        # censor/follow-up time. The legacy parser left this blank for controls.
        if p.get("event_time") is None:
            p["followup_time"] = b.get("rfs_months")
    return patients


def read_galaxy(root: Path):
    path = root / "data_raw" / "galaxy_2024_supplementary_tables_1_3.pdf"
    out = []
    with pdfplumber.open(path) as pdf:
        if len(pdf.pages) != 46:
            raise RuntimeError(f"GALAXY supplement page count changed: {len(pdf.pages)}")
        for page in pdf.pages[2:45]:
            tables = page.extract_tables()
            if not tables:
                continue
            for row in tables[0]:
                if len(row) < 19 or not re.fullmatch(r"\s*\d+\s*", clean(row[0])):
                    continue
                if clean(row[1]) not in {"Male", "Female"}:
                    continue
                r = [clean(x) for x in row[:19]]
                rec = clean(r[4]).lower() == "yes"
                dfs = num(r[6])
                out.append({
                    "cohort": "GALAXY", "patient_id": f"GALAXY-{int(r[0]):04d}",
                    "age": num(r[2]), "sex": r[1], "stage": r[3], "site": None,
                    "event_time_days": dfs if rec else None,
                    "followup_days": None if rec else dfs,
                    "recurred": rec, "dfs_days": dfs,
                    "mrd0": status01(r[10]), "mrd0_day": num(r[11]),
                    "mrd3": status01(r[12]), "mrd6": status01(r[13]),
                    "mrd0_day_exact": True, "mrd3_day_exact": False,
                    "mrd6_day_exact": False,
                })
    ids = [p["patient_id"] for p in out]
    if len(out) != 2240 or len(set(ids)) != 2240:
        raise RuntimeError(f"Expected 2,240 unique GALAXY rows, found {len(out)}")
    return out


def source_patients(root: Path):
    mod = load_legacy(root)
    galaxy = read_galaxy(root)
    cohorts = {
        "GALAXY": galaxy,
        "Chen": enrich_chen(root, mod),
        "ColonAiQ": enrich_coloniaiq(root, mod),
        "COSMOS": enrich_cosmos(root, mod),
    }
    for cohort, patients in cohorts.items():
        for p in patients:
            p["cohort"] = cohort
            p["sex_male"] = sex01(p.get("sex"))
            p["age"] = num(p.get("age"))
            p["stage_high"] = mod.stage_high(p.get("stage"))
            if cohort != "GALAXY":
                p["event_time_days"] = (num(p.get("event_time")) * DAY_PER_MONTH
                                         if num(p.get("event_time")) is not None else None)
                p["followup_days"] = (num(p.get("followup_time")) * DAY_PER_MONTH
                                       if num(p.get("followup_time")) is not None else None)
                hist = []
                for t, s in p.get("history", []):
                    if num(t) is not None and status01(s) is not None and float(t) > 0:
                        hist.append((float(t) * DAY_PER_MONTH, int(status01(s))))
                p["history_days"] = hist
                p["mrd0_status"], p["mrd0_day"] = pick_window(hist, MRD0_LOW, MRD0_HIGH - 1e-6, MRD0_TARGET)
                p["l3_status"], p["l3_day"] = pick_window(hist, L3_LOW, L3_HIGH, L3_TARGET)
                p["l6_status"], p["l6_day"] = pick_window(hist, L6_LOW, L6_HIGH, L6_TARGET)
                p["mrd0_day_exact"] = p["mrd3_day_exact"] = p["mrd6_day_exact"] = True
                p["mrd3"] = p["l3_status"]
                p["mrd6"] = p["l6_status"]
                p["mrd0"] = p["mrd0_status"]
            else:
                p["mrd0_status"] = p["mrd0"]
                p["l3_status"] = p["mrd3"]
                p["l6_status"] = p["mrd6"]
                p["l3_day"] = p["l6_day"] = None
                p["history_days"] = []
    return cohorts


def pick_window(hist, low, high, target):
    cand = [(t, s) for t, s in hist if low <= t <= high]
    if not cand:
        return None, None
    t, s = min(cand, key=lambda z: (abs(z[0] - target), z[0]))
    return int(s), float(t)


def outcome_at_origin(p, origin=ORIGIN_DAY):
    et, fu = p.get("event_time_days"), p.get("followup_days")
    if et is not None:
        if et <= origin:
            return False, None
        if et <= HORIZON_DAY:
            return True, 1
        return True, 0
    if fu is not None and fu >= HORIZON_DAY:
        return True, 0
    if fu is not None and fu >= origin:
        return True, None
    return False, None


def transition(a, b):
    names = {(0, 0): "negative_to_negative", (0, 1): "negative_to_positive",
             (1, 0): "positive_to_negative", (1, 1): "positive_to_positive"}
    return names[(int(a), int(b))]


def analysis_records(cohorts):
    result = {c: [] for c in cohorts}
    validation = {c: [] for c in cohorts}
    trajectories = {c: [] for c in cohorts}
    audit = []
    timing = []
    for cohort, patients in cohorts.items():
        counts = defaultdict(int)
        for p in patients:
            counts["total"] += 1
            counts["clinical_complete"] += all(p.get(k) is not None for k in ("age", "sex_male", "stage_high"))
            for key in ("mrd0_status", "l3_status", "l6_status"):
                counts[key] += p.get(key) in (0, 1)
            complete = all(p.get(k) in (0, 1) for k in ("mrd0_status", "l3_status", "l6_status"))
            counts["complete_3timepoints"] += complete
            if cohort == "GALAXY":
                prior_values = [s for s in (p.get("mrd0_status"), p.get("l3_status")) if s in (0, 1)]
                pre_l3_values = [p.get("mrd0_status")] if p.get("mrd0_status") in (0, 1) else []
            else:
                l6_day = p.get("l6_day")
                prior_values = [s for t, s in p.get("history_days", [])
                                if l6_day is not None and t < l6_day and s in (0, 1)]
                l3_day = p.get("l3_day")
                pre_l3_values = [s for t, s in p.get("history_days", [])
                                 if l3_day is not None and t < l3_day and s in (0, 1)]
            counts["prior_mrd_history_before_l6"] += bool(prior_values)
            at_risk, y = outcome_at_origin(p)
            counts["event_free_at_200d"] += at_risk
            counts["known_outcome_by_12m"] += at_risk and y is not None
            counts["events_200d_to_12m"] += at_risk and y == 1
            counts["controls_through_12m"] += at_risk and y == 0
            if cohort != "GALAXY" and p.get("l6_day") is not None:
                counts["l6_sample_160_200d"] += 1
                counts["l6_sample_at_or_before_6mo"] += p["l6_day"] <= L6_STRICT_DAY
            if cohort == "GALAXY":
                counts["l6_sample_160_200d"] += p.get("l6_status") in (0, 1)
                counts["strict_l6_date_known"] += 0
            clinical_ok = all(p.get(k) is not None for k in ("age", "sex_male", "stage_high"))
            if at_risk and clinical_ok and complete:
                traj = dict(p)
                traj["y"] = int(y) if y is not None else None
                traj["ever_prior_l3"] = int(any(s == 1 for s in pre_l3_values))
                traj["ever_prior_l6"] = int(any(s == 1 for s in prior_values))
                traj["transition"] = transition(p["l3_status"], p["l6_status"])
                trajectories[cohort].append(traj)
            if (at_risk and y is not None and clinical_ok and p.get("l6_status") in (0, 1)
                    and prior_values):
                val = dict(p)
                val["y"] = int(y)
                val["ever_prior_l6"] = int(any(s == 1 for s in prior_values))
                val["ever_prior_l3"] = int(any(s == 1 for s in pre_l3_values)) if pre_l3_values else 0
                val["transition"] = (transition(p["l3_status"], p["l6_status"])
                                      if p.get("l3_status") in (0, 1) else "l3_unavailable")
                val["has_pre_l6_history"] = True
                validation[cohort].append(val)
                counts["validation_risk_set"] += 1
                counts["validation_events"] += int(y == 1)
            # Strict matched-history risk set used for model fitting and comparison.
            if (at_risk and y is not None and complete and clinical_ok):
                rec = dict(p)
                rec["y"] = int(y)
                rec["transition"] = transition(p["l3_status"], p["l6_status"])
                rec["ever_prior_l3"] = int(any(s == 1 for s in pre_l3_values))
                rec["ever_prior_l6"] = int(any(s == 1 for s in prior_values))
                result[cohort].append(rec)
        audit.append({"cohort": cohort, **counts})
        timing.append({
            "cohort": cohort,
            "L3_rule": "study-defined 70-112 day" if cohort == "GALAXY" else "closest status within 70-112 days",
            "L6_rule": "study-defined 160-200 day; individual date absent" if cohort == "GALAXY" else "closest status within 160-200 days",
            "L6_records_160_200": counts["l6_sample_160_200d"],
            "L6_records_at_or_before_6mo": "not identifiable; individual dates absent" if cohort == "GALAXY" else counts["l6_sample_at_or_before_6mo"],
            "risk_origin_day": ORIGIN_DAY,
            "outcome_window": "(200 days, 365.24 days] recurrence",
            "patient_specific_L6_date_available": "no" if cohort == "GALAXY" else "yes",
        })
    return result, validation, trajectories, audit, timing


def age_parameters(rows):
    a = np.asarray([r["age"] for r in rows], dtype=float)
    mu, sd = float(a.mean()), float(a.std(ddof=0))
    return {"age_mean": mu, "age_sd": sd if sd > 1e-9 else 1.0}


def feature_vector(r, variant, scale):
    base = [float(r["stage_high"]), (float(r["age"]) - scale["age_mean"]) / scale["age_sd"], float(r["sex_male"])]
    names = ["stage_III_IV", "age_z", "sex_male"]
    if variant in ("L3_before", "L6_after"):
        key = "l3_status" if variant == "L3_before" else "l6_status"
        prior = "ever_prior_l3" if variant == "L3_before" else "ever_prior_l6"
        base += [float(r[key]), float(r[prior])]
        names += ["current_mrd_positive_at_3m" if variant == "L3_before" else "current_mrd_positive_at_6m",
                  "ever_prior_mrd_positive"]
    return np.asarray(base, dtype=float), names


def fit_model(rows, variant, mod):
    scale = age_parameters(rows)
    xrows = [feature_vector(r, variant, scale)[0] for r in rows]
    x = np.column_stack([np.ones(len(rows)), np.vstack(xrows)])
    y = np.asarray([r["y"] for r in rows], dtype=int)
    beta = fit_ridge_logit(x, y, penalty=1.0)
    names = ["intercept"] + feature_vector(rows[0], variant, scale)[1]
    return {"variant": variant, "scale": scale, "beta": beta, "terms": names, "ridge_lambda": 1.0}


def predict(model, rows, mod):
    x = np.vstack([feature_vector(r, model["variant"], model["scale"])[0] for r in rows])
    return mod.sigmoid(np.column_stack([np.ones(len(rows)), x]) @ model["beta"])


def metric_values(y, p, mod):
    y = np.asarray(y, dtype=int); p = np.asarray(p, dtype=float)
    lp = np.log(np.clip(p, 1e-8, 1-1e-8) / np.clip(1-p, 1e-8, 1))
    ci, cs = float("nan"), float("nan")
    if len(np.unique(y)) == 2 and float(np.std(lp)) >= 1e-8:
        try:
            cb = fit_ridge_logit(np.column_stack([np.ones(len(lp)), lp]), y, penalty=1e-6)
            if np.all(np.isfinite(cb)) and np.max(np.abs(cb)) < 10:
                ci, cs = float(cb[0]), float(cb[1])
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            pass
    return {"auc": mod.rank_auc(y, p), "average_precision": mod.average_precision(y, p),
            "brier": float(np.mean((y-p)**2)), "calibration_intercept": ci,
            "calibration_slope": cs, "mean_predicted_risk": float(np.mean(p))}


def fit_ridge_logit(x, y, penalty=1.0, max_iter=300):
    """Damped Newton solver with an objective line search; intercept is unpenalized."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    if len(np.unique(y)) < 2:
        raise ValueError("Both outcome classes are required for logistic fitting")
    npar = x.shape[1]
    beta = np.zeros(npar, dtype=float)
    prevalence = float(np.clip(y.mean(), 1e-5, 1-1e-5))
    beta[0] = math.log(prevalence/(1-prevalence))
    pmat = np.eye(npar, dtype=float) * penalty
    pmat[0, 0] = 0.0

    def objective(b):
        z = x @ b
        return float(np.sum(np.logaddexp(0.0, z) - y*z) + .5*b @ pmat @ b)

    loss = objective(beta)
    for _ in range(max_iter):
        z = np.clip(x @ beta, -35, 35)
        prob = 1.0/(1.0+np.exp(-z))
        w = np.maximum(prob*(1-prob), 1e-9)
        grad = x.T @ (prob-y) + pmat @ beta
        hess = (x.T*w) @ x + pmat + np.eye(npar)*1e-8
        step = np.linalg.solve(hess, grad)
        if float(np.max(np.abs(grad))) < 1e-7:
            break
        alpha = 1.0
        accepted = False
        descent = float(grad @ step)
        while alpha >= 1e-10:
            candidate = beta - alpha*step
            new_loss = objective(candidate)
            if new_loss <= loss - 1e-4*alpha*descent:
                beta, loss, accepted = candidate, new_loss, True
                break
            alpha *= .5
        if not accepted or float(np.max(np.abs(alpha*step))) < 1e-8:
            break
    if not np.all(np.isfinite(beta)):
        raise FloatingPointError("Non-finite ridge logistic coefficients")
    return beta


def percentile_ci(values):
    a = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if not len(a):
        return None, None
    return float(np.quantile(a, .025)), float(np.quantile(a, .975))


def bootstrap_internal(rows, mod, reps=BOOT_REPS, seed=SEED):
    rng = np.random.default_rng(seed)
    y = np.asarray([r["y"] for r in rows], dtype=int)
    full = {v: fit_model(rows, v, mod) for v in ("clinical", "L3_before", "L6_after")}
    apparent = {v: metric_values(y, predict(full[v], rows, mod), mod) for v in ("clinical", "L6_after")}
    optimism = {v: defaultdict(list) for v in ("clinical", "L6_after")}
    boot_test = {v: defaultdict(list) for v in ("clinical", "L6_after")}
    oob_sum = {v: np.zeros(len(rows), dtype=float) for v in ("clinical", "L3_before", "L6_after")}
    oob_n = np.zeros(len(rows), dtype=int)
    successes = 0
    for b in range(reps):
        idx = rng.integers(0, len(rows), len(rows))
        boot = [rows[int(i)] for i in idx]
        unique = np.zeros(len(rows), dtype=bool); unique[np.unique(idx)] = True
        try:
            models = {v: fit_model(boot, v, mod) for v in ("clinical", "L3_before", "L6_after")}
            yb = np.asarray([r["y"] for r in boot], dtype=int)
            for v in ("clinical", "L6_after"):
                mb = metric_values(yb, predict(models[v], boot, mod), mod)
                mt = metric_values(y, predict(models[v], rows, mod), mod)
                for metric in ("auc", "average_precision", "brier", "calibration_intercept", "calibration_slope"):
                    if np.isfinite(mb[metric]) and np.isfinite(mt[metric]):
                        optimism[v][metric].append(mb[metric] - mt[metric])
                        boot_test[v][metric].append(mt[metric])
            oob = np.flatnonzero(~unique)
            if len(oob):
                for v in ("clinical", "L3_before", "L6_after"):
                    oob_sum[v][oob] += predict(models[v], [rows[int(i)] for i in oob], mod)
                oob_n[oob] += 1
            successes += 1
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
    if successes < reps * .95:
        raise RuntimeError(f"Only {successes}/{reps} bootstrap fits succeeded")
    oob_pred = {v: np.divide(oob_sum[v], oob_n, out=np.full(len(rows), np.nan), where=oob_n > 0)
                for v in ("clinical", "L3_before", "L6_after")}
    rows_out = []
    for v in ("clinical", "L6_after"):
        for metric, app in apparent[v].items():
            if metric not in optimism[v] or not optimism[v][metric]:
                continue
            opt = float(np.mean(optimism[v][metric]))
            corrected_dist = [app - x for x in optimism[v][metric]]
            lo, hi = percentile_ci(corrected_dist)
            rows_out.append({"model": "Clinical-only" if v == "clinical" else "Clinical+MRD",
                             "metric": metric, "apparent": app, "mean_bootstrap_optimism": opt,
                             "optimism_corrected": app - opt, "corrected_ci_low": lo,
                             "corrected_ci_high": hi, "bootstrap_refits": successes})
    return full, rows_out, oob_pred, oob_n, apparent


def calibration_bins(y, p, cohort, model, bins=10):
    order = np.argsort(p)
    chunks = np.array_split(order, min(bins, len(order)))
    out = []
    for i, ids in enumerate(chunks, 1):
        if not len(ids):
            continue
        out.append({"cohort": cohort, "model": model, "bin": i, "n": len(ids),
                    "events": int(np.asarray(y)[ids].sum()),
                    "mean_predicted_risk": float(np.mean(np.asarray(p)[ids])),
                    "observed_event_rate": float(np.mean(np.asarray(y)[ids]))})
    return out


def wilson(x, n):
    if not n:
        return None, None
    z = 1.959963984540054
    p = x / n
    den = 1 + z*z/n
    mid = (p + z*z/(2*n)) / den
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return max(0.0, mid-half), min(1.0, mid+half)


def bootstrap_group_medians(group, reps=BOOT_REPS, seed=SEED):
    rng = np.random.default_rng(seed + len(group))
    before = np.asarray([r["risk_before_oob"] for r in group], dtype=float)
    after = np.asarray([r["risk_after_oob"] for r in group], dtype=float)
    vals = np.empty((reps, 3), dtype=float)
    for b in range(reps):
        ix = rng.integers(0, len(group), len(group))
        vals[b] = [np.median(before[ix]), np.median(after[ix]), np.median(after[ix]-before[ix])]
    cis = [percentile_ci(vals[:, j]) for j in range(3)]
    return cis


def paired_delta_ci(y, p0, p1, mod, reps=BOOT_REPS, seed=SEED):
    rng = np.random.default_rng(seed + len(y))
    y = np.asarray(y, dtype=int); p0 = np.asarray(p0); p1 = np.asarray(p1)
    store = defaultdict(list); valid = 0
    for _ in range(reps):
        ix = rng.integers(0, len(y), len(y))
        if len(np.unique(y[ix])) < 2:
            continue
        a = metric_values(y[ix], p0[ix], mod); b = metric_values(y[ix], p1[ix], mod)
        for m in ("auc", "average_precision", "brier"):
            if np.isfinite(a[m]) and np.isfinite(b[m]):
                store[m].append(b[m]-a[m])
        valid += 1
    return {m: (*percentile_ci(store[m]), valid) for m in ("auc", "average_precision", "brier")}


def csv_write(path, rows):
    rows = list(rows)
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def svg_write(path, body, width=1100, height=700):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/>{body}</svg>', encoding="utf-8")


def esc(x):
    return (str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def make_figures(figdir, metrics, updates, memory, transitions, cases):
    figdir.mkdir(parents=True, exist_ok=True)
    body = '<text x="550" y="45" text-anchor="middle" font-size="26" font-family="Arial">GALAXY 主开发与冻结模型公共验证</text>'
    labels = ["GALAXY development", "Freeze coefficients", "Chen validation", "COSMOS validation", "ColonAiQ if compatible", "Future institutional external validation"]
    xs = [90, 440, 790, 790, 790, 1010]
    ys = [115, 115, 70, 160, 250, 115]
    for i, (lab, x, y) in enumerate(zip(labels, xs, ys)):
        body += f'<rect x="{x}" y="{y}" width="220" height="64" rx="12" fill="{("#dbeafe" if i<2 else "#e7f5ec")}" stroke="#334155"/><text x="{x+110}" y="{y+27}" text-anchor="middle" font-size="15" font-family="Arial">{esc(lab)}</text>'
    for y in [102, 192, 282]:
        body += f'<path d="M 1010 147 C 1050 {y}, 1040 {y}, 1010 {y}" fill="none" stroke="#64748b" stroke-width="2"/>'
    body += '<path d="M 310 147 L 430 147" stroke="#334155" stroke-width="3" marker-end="url(#a)"/><defs><marker id="a" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#334155"/></marker></defs>'
    svg_write(figdir / "study_design.svg", body, 1260, 390)

    cohorts = [r["cohort"] for r in metrics if r["model"] == "Clinical-only"]
    body = '<text x="550" y="35" text-anchor="middle" font-size="23" font-family="Arial">冻结模型的外部区分度、校准与误差</text>'
    panels = [("auc", "ROC-AUC", .5, 1), ("brier", "Brier score", 0, max(.25, max([float(r['brier']) for r in metrics if r.get('brier') not in (None, '')] or [.25]))), ("calibration_slope", "Calibration slope", 0, 2)]
    for pi, (metric, title, lo, hi) in enumerate(panels):
        ox=70+pi*355; oy=95; pw=300; ph=420
        body += f'<text x="{ox+pw/2}" y="{oy-15}" text-anchor="middle" font-size="17" font-family="Arial">{title}</text><line x1="{ox}" y1="{oy+ph}" x2="{ox+pw}" y2="{oy+ph}" stroke="#64748b"/><line x1="{ox}" y1="{oy}" x2="{ox}" y2="{oy+ph}" stroke="#64748b"/>'
        for ci, cohort in enumerate(cohorts):
            barw=26; gap=10; groupw=68; gx=ox+25+ci*groupw
            for mi, model in enumerate(("Clinical-only", "Clinical+MRD")):
                row=next((z for z in metrics if z["cohort"]==cohort and z["model"]==model),None)
                val=float(row[metric]) if row and row.get(metric) not in (None, "") and np.isfinite(float(row[metric])) else None
                if val is not None:
                    vv=min(hi,max(lo,val)); bh=(vv-lo)/(hi-lo)*ph
                    body += f'<rect x="{gx+mi*(barw+4)}" y="{oy+ph-bh}" width="{barw}" height="{bh}" fill="{("#64748b" if mi==0 else "#e07a5f")}"/><text x="{gx+mi*(barw+4)+barw/2}" y="{oy+ph+18}" text-anchor="middle" font-size="10" font-family="Arial">{model.split("-")[0]}</text>'
            body += f'<text x="{gx+30}" y="{oy+ph+42}" text-anchor="middle" font-size="11" font-family="Arial">{esc(cohort)}</text>'
        body += f'<text x="{ox-15}" y="{oy+ph/2}" transform="rotate(-90 {ox-15} {oy+ph/2})" text-anchor="middle" font-size="12" font-family="Arial">Metric</text>'
    svg_write(figdir / "external_validation_auc_brier.svg", body)

    body = '<text x="550" y="38" text-anchor="middle" font-size="22" font-family="Arial">同一(200天,12月]结局：L3信息到L6信息的风险更新（OOB）</text>'
    order=["negative_to_negative","negative_to_positive","positive_to_negative","positive_to_positive"]
    zh={"negative_to_negative":"阴→阴","negative_to_positive":"阴→阳","positive_to_negative":"阳→阴","positive_to_positive":"阳→阳"}
    cols={"negative_to_negative":"#238b67","negative_to_positive":"#e67e22","positive_to_negative":"#4c91c6","positive_to_positive":"#b82e38"}
    x0=170; y0=100; rowh=115; scale=750
    for i,state in enumerate(order):
        row=next((r for r in updates if r["transition"]==state),None)
        if not row: continue
        y=y0+i*rowh
        before=float(row["median_risk_before"]); after=float(row["median_risk_after"])
        body += f'<text x="45" y="{y+28}" font-size="16" font-family="Arial">{zh[state]} (n={row["n"]})</text><line x1="{x0}" y1="{y+24}" x2="{x0+before*scale}" y2="{y+24}" stroke="#64748b" stroke-width="7"/><line x1="{x0}" y1="{y+52}" x2="{x0+after*scale}" y2="{y+52}" stroke="{cols[state]}" stroke-width="7"/><text x="{x0+before*scale+8}" y="{y+28}" font-size="12" font-family="Arial">before {before:.1%}</text><text x="{x0+after*scale+8}" y="{y+56}" font-size="12" font-family="Arial">after {after:.1%}; Δ {float(row["median_delta_risk"]):+.1%}</text>'
    body += '<text x="75" y="82" font-size="12" font-family="Arial">灰色：L3可用信息；彩色：L6更新后；均针对同一预测窗口</text>'
    svg_write(figdir / "same_horizon_risk_update.svg", body)

    body = '<text x="550" y="38" text-anchor="middle" font-size="22" font-family="Arial">L6 MRD 阴性患者的分子记忆：从未阳性 vs 既往阳性</text>'
    memcohorts=sorted(set(r["cohort"] for r in memory))
    ybase=105
    for ci,c in enumerate(memcohorts):
        y=ybase+ci*115
        body += f'<text x="40" y="{y+22}" font-size="15" font-family="Arial">{esc(c)}</text>'
        for gi,g in enumerate(("never_positive","previously_positive_now_negative")):
            r=next((x for x in memory if x["cohort"]==c and x["memory_group"]==g),None)
            if r and int(r["n"]):
                x=215+gi*390
                obs=float(r["observed_recurrence_rate"]); pred=float(r["mean_predicted_risk"])
                body += f'<text x="{x}" y="{y+18}" font-size="12" font-family="Arial">{("Never positive" if gi==0 else "Prior +, current -")} n={r["n"]}</text><rect x="{x}" y="{y+30}" width="{obs*230}" height="18" fill="#d95d39"/><text x="{x+240}" y="{y+44}" font-size="12" font-family="Arial">actual {obs:.1%}</text><rect x="{x}" y="{y+58}" width="{pred*230}" height="18" fill="#457b9d"/><text x="{x+240}" y="{y+72}" font-size="12" font-family="Arial">predicted {pred:.1%}</text>'
    svg_write(figdir / "molecular_memory_validation.svg", body, 1150, max(360, ybase+len(memcohorts)*115))

    body = '<text x="650" y="35" text-anchor="middle" font-size="23" font-family="Arial">各队列3月→6月MRD状态转移与(200天,12月]复发</text>'
    coh=sorted(set(r["cohort"] for r in transitions))
    trans_order=order
    for ci,c in enumerate(coh):
        y=85+ci*125
        body += f'<text x="40" y="{y+18}" font-size="15" font-family="Arial">{esc(c)}</text>'
        for ti,t in enumerate(trans_order):
            r=next((x for x in transitions if x["cohort"]==c and x["transition"]==t),None)
            if not r or not int(r["n"]): continue
            x=180+ti*245
            rate=float(r["event_rate"])
            body += f'<text x="{x}" y="{y+17}" font-size="12" font-family="Arial">{zh[t]} n={r["n"]}</text><rect x="{x}" y="{y+28}" width="200" height="18" fill="#e5e7eb"/><rect x="{x}" y="{y+28}" width="{200*rate}" height="18" fill="{cols[t]}"/><text x="{x}" y="{y+66}" font-size="12" font-family="Arial">复发 {rate:.1%} ({r["events"]})</text>'
    svg_write(figdir / "dynamic_state_transition.svg", body, 1200, 85+len(coh)*125)

    # Local-only examples are actual records selected from the observed GALAXY cohort.
    body = '<text x="500" y="38" text-anchor="middle" font-size="22" font-family="Arial">GALAXY真实患者的代表性MRD风险轨迹（仅本地）</text>'
    for i,c in enumerate(cases):
        y=105+i*115
        body += f'<text x="45" y="{y+15}" font-size="14" font-family="Arial">{esc(c["transition_zh"])} | {esc(c["patient_id"])} | {esc(c["outcome_label"])}</text><line x1="100" y1="{y+44}" x2="{100+float(c["risk_before_update"])*750}" y2="{y+44}" stroke="#64748b" stroke-width="7"/><line x1="100" y1="{y+72}" x2="{100+float(c["risk_after_update"])*750}" y2="{y+72}" stroke="{c["color"]}" stroke-width="7"/><text x="880" y="{y+49}" font-size="12" font-family="Arial">{float(c["risk_before_update"]):.1%} → {float(c["risk_after_update"]):.1%}</text>'
    svg_write(Path(cases[0]["figure_path"]), body, 1020, 105+len(cases)*115)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=r"E:\crc_mrd_digital_twin")
    ap.add_argument("--output-root", required=True)
    ap.add_argument("--bootstrap", type=int, default=BOOT_REPS)
    args = ap.parse_args()
    root, outroot = Path(args.project_root), Path(args.output_root)
    tables, reports, figs, local, modelsdir = [outroot / d for d in
        ("outputs/tables", "reports", "outputs/figures", "outputs/local_only", "outputs/models")]
    for d in (tables, reports, figs, local, modelsdir): d.mkdir(parents=True, exist_ok=True)

    cohorts = source_patients(root)
    rows_by_cohort, validation_by_cohort, trajectories_by_cohort, audit, timing = analysis_records(cohorts)
    galaxy = rows_by_cohort["GALAXY"]
    if len(galaxy) < 100 or sum(r["y"] for r in galaxy) < 20:
        raise RuntimeError(f"Insufficient GALAXY L6-window training sample: n={len(galaxy)} events={sum(r['y'] for r in galaxy)}")
    mod = load_legacy(root)
    full_models, boot_rows, oob, oob_n, apparent = bootstrap_internal(galaxy, mod, args.bootstrap)

    # The frozen final models are fitted once on all eligible GALAXY development rows.
    # External cohorts do not enter this step in any way.
    beta_rows = []
    for variant, model in full_models.items():
        model_name = {"clinical":"Clinical-only", "L3_before":"Same-horizon before update", "L6_after":"Clinical+MRD / after update"}[variant]
        for term, coefficient in zip(model["terms"], model["beta"]):
            beta_rows.append({"model": model_name, "term": term, "coefficient": float(coefficient),
                              "odds_ratio_per_unit": float(math.exp(np.clip(coefficient, -30, 30))),
                              **model["scale"], "ridge_lambda": model["ridge_lambda"], "development_n": len(galaxy),
                              "development_events": int(sum(r["y"] for r in galaxy))})
    csv_write(tables / "galaxy_model_coefficients.csv", beta_rows)
    (modelsdir / "galaxy_frozen_ridge_models.json").write_text(json.dumps({
        "development_cohort":"GALAXY", "origin_day":ORIGIN_DAY, "horizon_day":HORIZON_DAY,
        "outcome":"clinical recurrence in (200 days, 365.24 days]", "ridge_lambda":1.0,
        "models":{k:{"terms":v["terms"],"coefficients":[float(x) for x in v["beta"]],"age_scale":v["scale"]} for k,v in full_models.items()},
        "feature_note":"GALAXY patient-level site is unavailable; clinical core uses stage, age, sex. MRD is binary/assay-agnostic. Ever-prior-positive summarizes observed pre-L6 statuses; prior testing density differs by cohort."
    }, ensure_ascii=False, indent=2), encoding="utf-8")

    boot_table = []
    for r in boot_rows:
        boot_table.append(r)
    # OOB predictions yield the same-horizon update estimates and development calibration curve.
    for i, r in enumerate(galaxy):
        r["risk_before_oob"] = float(oob["L3_before"][i])
        r["risk_after_oob"] = float(oob["L6_after"][i])
    calibration_rows = []
    for key,label in (("clinical","Clinical-only OOB"),
                      ("L3_before","Same-horizon L3 history OOB"),
                      ("L6_after","Clinical+MRD OOB")):
        valid=np.isfinite(oob[key])
        calibration_rows.extend(calibration_bins(np.asarray([r["y"] for r,v in zip(galaxy,valid) if v]),
                                                  oob[key][valid], "GALAXY", label))

    # External validation: apply frozen GALAXY coefficients with no refitting or recalibration.
    metric_rows, increment_rows, calibration_ext = [], [], []
    evaluated = {"GALAXY": galaxy}
    for cohort in ("Chen", "COSMOS", "ColonAiQ"):
        evaluated[cohort] = validation_by_cohort.get(cohort, [])
    all_ext_y, all_ext_c, all_ext_m = [], [], []
    ext_for_memory = {}
    for cohort, rs in evaluated.items():
        if not rs:
            continue
        y=np.asarray([r["y"] for r in rs],dtype=int)
        pc=predict(full_models["clinical"],rs,mod)
        pm=predict(full_models["L6_after"],rs,mod)
        ext_for_memory[cohort]=(rs,pm)
        for name,pred in (("Clinical-only",pc),("Clinical+MRD",pm)):
            met=metric_values(y,pred,mod)
            metric_rows.append({"cohort":cohort,"model":name,"n":len(rs),"events":int(y.sum()),
                                "event_rate":float(y.mean()),**met,"model_frozen":cohort!="GALAXY",
                                "validation_role":"development bootstrap" if cohort=="GALAXY" else "independent public validation"})
            calibration_ext.extend(calibration_bins(y,pred,cohort,name))
        deltas=paired_delta_ci(y,pc,pm,mod,BOOT_REPS,SEED+len(cohort))
        for m,(lo,hi,nvalid) in deltas.items():
            inc={"cohort":cohort,"n":len(rs),"events":int(y.sum()),"metric":m,
                 "clinical_only":metric_values(y,pc,mod)[m],"clinical_plus_mrd":metric_values(y,pm,mod)[m],
                 "delta_mrd_minus_clinical":metric_values(y,pm,mod)[m]-metric_values(y,pc,mod)[m],
                 "bootstrap_ci_low":lo,"bootstrap_ci_high":hi,"valid_bootstrap_replicates":nvalid}
            increment_rows.append(inc)
        if cohort != "GALAXY":
            all_ext_y.extend(y.tolist()); all_ext_c.extend(pc.tolist()); all_ext_m.extend(pm.tolist())
    if all_ext_y:
        y=np.asarray(all_ext_y,dtype=int); pc=np.asarray(all_ext_c); pm=np.asarray(all_ext_m)
        for name,pred in (("Clinical-only",pc),("Clinical+MRD",pm)):
            metric_rows.append({"cohort":"Pooled external","model":name,"n":len(y),"events":int(y.sum()),
                                "event_rate":float(y.mean()),**metric_values(y,pred,mod),"model_frozen":True,
                                "validation_role":"pooled descriptive external validation"})
        d=paired_delta_ci(y,pc,pm,mod)
        for m,(lo,hi,nvalid) in d.items():
            increment_rows.append({"cohort":"Pooled external","n":len(y),"events":int(y.sum()),"metric":m,
                "clinical_only":metric_values(y,pc,mod)[m],"clinical_plus_mrd":metric_values(y,pm,mod)[m],
                "delta_mrd_minus_clinical":metric_values(y,pm,mod)[m]-metric_values(y,pc,mod)[m],
                "bootstrap_ci_low":lo,"bootstrap_ci_high":hi,"valid_bootstrap_replicates":nvalid})
    calibration_rows.extend(calibration_ext)

    # Overall same-horizon performance uses paired OOB predictions, not apparent fits.
    same_y=np.asarray([r["y"] for r in galaxy],dtype=int)
    before_metrics=metric_values(same_y,oob["L3_before"],mod)
    after_metrics=metric_values(same_y,oob["L6_after"],mod)
    same_delta_ci=paired_delta_ci(same_y,oob["L3_before"],oob["L6_after"],mod,BOOT_REPS,SEED+777)
    same_perf_rows=[]
    for metric in ("auc","average_precision","brier"):
        lo,hi,nvalid=same_delta_ci[metric]
        same_perf_rows.append({"metric":metric,"before_L3_history_oob":before_metrics[metric],
            "after_L6_history_oob":after_metrics[metric],
            "delta_after_minus_before":after_metrics[metric]-before_metrics[metric],
            "delta_ci_low":lo,"delta_ci_high":hi,"valid_bootstrap_replicates":nvalid,
            "n":len(galaxy),"events":int(same_y.sum()),"same_outcome_window":"(200 days, 365.24 days]"})

    # Same-horizon patient-level update summary, using paired OOB predictions in GALAXY.
    update_rows=[]
    for tr in ("negative_to_negative","negative_to_positive","positive_to_negative","positive_to_positive"):
        g=[r for r in galaxy if r["transition"]==tr and np.isfinite(r["risk_before_oob"]) and np.isfinite(r["risk_after_oob"])]
        if not g:
            continue
        ci=bootstrap_group_medians(g,BOOT_REPS,SEED+len(tr))
        events=sum(r["y"] for r in g); lo,hi=wilson(events,len(g))
        update_rows.append({"transition":tr,"transition_zh":{"negative_to_negative":"阴→阴","negative_to_positive":"阴→阳","positive_to_negative":"阳→阴","positive_to_positive":"阳→阳"}[tr],
            "n":len(g),"events":events,"observed_recurrence_rate":events/len(g),"event_rate_ci_low":lo,"event_rate_ci_high":hi,
            "median_risk_before":float(np.median([r["risk_before_oob"] for r in g])),
            "median_risk_after":float(np.median([r["risk_after_oob"] for r in g])),
            "median_delta_risk":float(np.median([r["risk_after_oob"]-r["risk_before_oob"] for r in g])),
            "before_ci_low":ci[0][0],"before_ci_high":ci[0][1],"after_ci_low":ci[1][0],"after_ci_high":ci[1][1],
            "delta_ci_low":ci[2][0],"delta_ci_high":ci[2][1],"prediction_method":"paired patient-level bootstrap out-of-bag",
            "bootstrap_replicates":BOOT_REPS})

    # Descriptive transitions are shown for all cohorts on the same L6-window-end risk set.
    transition_rows=[]
    for cohort,rs in rows_by_cohort.items():
        for tr in ("negative_to_negative","negative_to_positive","positive_to_negative","positive_to_positive"):
            g=[r for r in rs if r["transition"]==tr]
            if not g: continue
            ev=sum(r["y"] for r in g); lo,hi=wilson(ev,len(g))
            transition_rows.append({"cohort":cohort,"transition":tr,"n":len(g),"events":ev,
                                    "event_rate":ev/len(g),"event_rate_ci_low":lo,"event_rate_ci_high":hi,
                                    "risk_origin_day":ORIGIN_DAY,"outcome_window":"(200 days, 365.24 days]"})

    # Molecular memory at current L6-negative status; predictions use OOB in GALAXY.
    memory_rows=[]
    memory_source = dict(validation_by_cohort)
    memory_source["GALAXY"] = rows_by_cohort["GALAXY"]
    for cohort,rs in memory_source.items():
        if not rs:
            continue
        predmap={r["patient_id"]:float(pred) for r,pred in zip(rs, predict(full_models["L6_after"],rs,mod))}
        if cohort=="GALAXY":
            predmap.update({r["patient_id"]:r["risk_after_oob"] for r in galaxy})
        for group,cond in (("never_positive",lambda r:r["ever_prior_l6"]==0),
                           ("previously_positive_now_negative",lambda r:r["ever_prior_l6"]==1)):
            g=[r for r in rs if r["l6_status"]==0 and cond(r)]
            if not g: continue
            ev=sum(r["y"] for r in g); lo,hi=wilson(ev,len(g))
            pr=[predmap[r["patient_id"]] for r in g if r["patient_id"] in predmap]
            memory_rows.append({"cohort":cohort,"memory_group":group,"n":len(g),"events":ev,
                "observed_recurrence_rate":ev/len(g),"observed_rate_ci_low":lo,"observed_rate_ci_high":hi,
                "mean_predicted_risk":float(np.mean(pr)) if pr else None,
                "median_predicted_risk":float(np.median(pr)) if pr else None,
                "prediction_source":"OOB bootstrap" if cohort=="GALAXY" else "frozen GALAXY model",
                "risk_origin_day":ORIGIN_DAY,"outcome_window":"(200 days, 365.24 days]"})

    # Patient-specific GALAXY update rows stay in outputs/local_only and are never pushed.
    dynamic_rows=[]; case_candidates=[]
    full_before=predict(full_models["L3_before"],galaxy,mod)
    full_after=predict(full_models["L6_after"],galaxy,mod)
    galaxy_trajectories=trajectories_by_cohort["GALAXY"]
    trajectory_before=predict(full_models["L3_before"],galaxy_trajectories,mod)
    trajectory_after=predict(full_models["L6_after"],galaxy_trajectories,mod)
    for i,r in enumerate(galaxy_trajectories):
        state=r["transition"]
        et=r.get("event_time_days")
        if et is not None:
            if et<=HORIZON_DAY: label="recurrence_by_12m"
            else: label="recurrence_after_12m"
        elif r.get("followup_days") is not None and r["followup_days"]>=HORIZON_DAY:
            label="no_recurrence_through_12m"
        else:
            label="censored_before_12m"
        d={"patient_id":r["patient_id"],"timepoint":"L6_study_window",
           "current_mrd":r["l6_status"],"ever_prior_positive":r["ever_prior_l6"],"state":state,
           "risk_before_update":float(trajectory_before[i]),"risk_after_update":float(trajectory_after[i]),
           "risk_delta":float(trajectory_after[i]-trajectory_before[i]),"recurrence_time_days":et or "",
           "recurrence_status":label,"followup_days":r.get("followup_days") or "",
           "endpoint_200d_to_12m":r["y"] if r.get("y") is not None else "","risk_origin_day":ORIGIN_DAY}
        dynamic_rows.append(d)
        if r.get("y") is not None:
            case_candidates.append({**d,"transition_zh":{"negative_to_negative":"阴→阴","negative_to_positive":"阴→阳","positive_to_negative":"阳→阴","positive_to_positive":"阳→阳"}[state],
                                    "outcome_label":label,"color":{"negative_to_negative":"#238b67","negative_to_positive":"#e67e22","positive_to_negative":"#4c91c6","positive_to_positive":"#b82e38"}[state],
                                    "figure_path":str(local/"galaxy_l6_l12"/"digital_twin_case_examples.svg")})
    selected_cases=[]
    for tr in ("negative_to_negative","negative_to_positive","positive_to_negative","positive_to_positive"):
        g=[x for x in case_candidates if x["state"]==tr]
        if g:
            med=float(np.median([x["risk_after_update"] for x in g]))
            selected_cases.append(min(g,key=lambda x:abs(x["risk_after_update"]-med)))

    # Save aggregate files and local-only patient trajectories.
    csv_write(tables/"galaxy_model_coefficients.csv",beta_rows)
    csv_write(tables/"galaxy_internal_bootstrap_metrics.csv",boot_table)
    csv_write(tables/"external_validation_metrics.csv",metric_rows)
    csv_write(tables/"clinical_vs_mrd_increment.csv",increment_rows)
    csv_write(tables/"same_horizon_update_summary.csv",update_rows)
    csv_write(tables/"same_horizon_model_performance.csv",same_perf_rows)
    csv_write(tables/"molecular_memory_by_cohort.csv",memory_rows)
    csv_write(tables/"dynamic_state_transition_summary.csv",transition_rows)
    csv_write(tables/"galaxy_cohort_missingness_and_sample_audit.csv",audit)
    csv_write(tables/"mrd_window_harmonization_audit.csv",timing)
    csv_write(tables/"galaxy_calibration_curve_points.csv",calibration_rows)
    csv_write(local/"patient_dynamic_predictions.csv",dynamic_rows)
    csv_write(local/"digital_twin_case_examples.csv",selected_cases)
    make_figures(figs,metric_rows,update_rows,memory_rows,transition_rows,selected_cases)

    def fmt(x):
        try: return f"{float(x):.3f}"
        except (ValueError,TypeError): return "NA"
    g_clin=next(r for r in metric_rows if r["cohort"]=="GALAXY" and r["model"]=="Clinical-only")
    g_mrd=next(r for r in metric_rows if r["cohort"]=="GALAXY" and r["model"]=="Clinical+MRD")
    ext_lines=[]
    for c in ("Chen","COSMOS","ColonAiQ","Pooled external"):
        a=next((r for r in metric_rows if r["cohort"]==c and r["model"]=="Clinical-only"),None)
        b=next((r for r in metric_rows if r["cohort"]==c and r["model"]=="Clinical+MRD"),None)
        if a and b:
            ext_lines.append(f"| {c} | {b['n']} | {b['events']} | {fmt(a['auc'])} / {fmt(a['average_precision'])} / {fmt(a['brier'])} | {fmt(b['auc'])} / {fmt(b['average_precision'])} / {fmt(b['brier'])} |")
    group_lines=[f"| {r['transition_zh']} | {r['n']} | {r['events']} | {r['observed_recurrence_rate']:.1%} | {r['median_risk_before']:.1%} | {r['median_risk_after']:.1%} | {r['median_delta_risk']:+.1%} ({r['delta_ci_low']:+.1%}, {r['delta_ci_high']:+.1%}) |" for r in update_rows]
    perf_lines=[f"| {r['metric']} | {r['before_L3_history_oob']:.3f} | {r['after_L6_history_oob']:.3f} | {r['delta_after_minus_before']:+.3f} ({r['delta_ci_low']:+.3f}, {r['delta_ci_high']:+.3f}) |" for r in same_perf_rows]
    memory_lines=[f"| {r['cohort']} | {r['memory_group']} | {r['n']} | {r['events']} | {r['observed_recurrence_rate']:.1%} | {r['mean_predicted_risk']:.1%} |" for r in memory_rows]
    report=f"""# GALAXY 开发、独立公共验证与动态风险更新

## 分析目标与时间定义

开发队列为 GALAXY 2024；Chen、COSMOS 和 ColonAiQ 在模型冻结后分别验证。结局为临床复发，不含 ctDNA molecular recurrence。临床变量限定为分期、年龄、性别；GALAXY 患者级补充表没有肿瘤部位字段，因此核心模型不含 site。MRD 使用阳性/阴性，不拼接 VAF 或不同 assay 的连续值。

**时间泄漏保护：** GALAXY 的 study-defined 6个月采样窗口为术后160–200天，但表2没有个体采血日。主分析将共同 landmark 放在第200天，以确保整段 GALAXY 采样窗口已结束；纳入者需在第200天仍未复发，预测结局为 **(200天, 365.24天]** 临床复发。这是对临床 L6→12m 问题的保守近似，并非严格的(182.6天,12月]。对 GALAXY 严格 ≤182.6天敏感性分析无法识别，不能臆造采血日。

其他队列的主分析状态取各自160–200天窗口内最接近182.6天的实测结果；L3 取70–112天，MRD window 取14–<70天。所有结局使用同一第200天起点，且事件必须晚于起点。窗口和可用人数见 `mrd_window_harmonization_audit.csv`。

## GALAXY 开发

匹配的三时间点风险集：N={len(galaxy)}，(200天,12月]复发={int(sum(r['y'] for r in galaxy))}。ridge logistic 的惩罚系数 λ=1；年龄以开发集均值{full_models['L6_after']['scale']['age_mean']:.2f}岁和标准差{full_models['L6_after']['scale']['age_sd']:.2f}标准化。Clinical-only 与 Clinical+MRD 在完全相同患者上比较。MRD模型变量为 L6 当前阳性及既往已观测 MRD 状态中任一次阳性；GALAXY 可用的既往列为 MRD-window/L3。外部队列按其实际 L6 前采样历史派生同一“ever observed positive”概念，检测频率不同会影响该特征的观察机会。

| 模型 | apparent AUC | optimism-corrected AUC | apparent AP | corrected AP | apparent Brier | corrected Brier |
|---|---:|---:|---:|---:|---:|---:|
| Clinical-only | {fmt(g_clin['auc'])} | {fmt(next(r['optimism_corrected'] for r in boot_rows if r['model']=='Clinical-only' and r['metric']=='auc'))} | {fmt(g_clin['average_precision'])} | {fmt(next(r['optimism_corrected'] for r in boot_rows if r['model']=='Clinical-only' and r['metric']=='average_precision'))} | {fmt(g_clin['brier'])} | {fmt(next(r['optimism_corrected'] for r in boot_rows if r['model']=='Clinical-only' and r['metric']=='brier'))} |
| Clinical+MRD | {fmt(g_mrd['auc'])} | {fmt(next(r['optimism_corrected'] for r in boot_rows if r['model']=='Clinical+MRD' and r['metric']=='auc'))} | {fmt(g_mrd['average_precision'])} | {fmt(next(r['optimism_corrected'] for r in boot_rows if r['model']=='Clinical+MRD' and r['metric']=='average_precision'))} | {fmt(g_mrd['brier'])} | {fmt(next(r['optimism_corrected'] for r in boot_rows if r['model']=='Clinical+MRD' and r['metric']=='brier'))} |

患者级 bootstrap 重抽样 {args.bootstrap} 次。校准截距、斜率及OOB校准分箱见相应表；完整回归系数与标准化参数已冻结至 `outputs/models/galaxy_frozen_ridge_models.json`。没有使用外部验证结果调参、重校准或选阈值。

## 冻结模型的公共验证

表内顺序为 AUC / average precision / Brier；校准截距和斜率见 `external_validation_metrics.csv`。

| 队列 | N | events | Clinical-only | Clinical+MRD |
|---|---:|---:|---:|---:|
{chr(10).join(ext_lines)}

外部队列按各自窗口、相同结局起点与同一套冻结系数评估；未做验证集特征选择、调参或重校准。MRD 增量的配对 bootstrap CI 见 `clinical_vs_mrd_increment.csv`。评估数据来自既有公开补充材料，样本量/事件少的队列其 calibration slope 与 AP 会不稳定，应结合 N、events 和区间解释。

## 分子记忆

在 L6 ctDNA 阴性者中，既往组定义为 L6 前至少一次已观测阳性；从未阳性组定义为所有已观测既往结果均阴性。

| 队列 | 既往状态 | N | 复发 | 观察复发率 | 平均预测风险 |
|---|---|---:|---:|---:|---:|
{chr(10).join(memory_lines)}

GALAXY 的分子记忆在 OOB 估计中有清楚分层；验证队列的“既往阳性、当前阴性”人数/事件很少，不能据此宣称已跨队列复制。

## 局限

1. 第200天起点使时间泄漏风险最低，但排除了第182.6–200天复发；结论不能直接称为精确6月至12月复发风险。
2. GALAXY 缺少患者级3/6月采血日，不能完成严格≤6月敏感性分析，也不能用真实日期拟合患者特异的动态生存模型。
3. GALAXY 没有个体 site，因此临床底模使用 stage、age、sex；无法检验加入部位后的效能。
4. 化疗时间、CEA、影像随访等没有在所有公共队列同定义，未进入公共 core model；模型估计关联与预测，不解释 MRD 清除的因果或治疗反应效应。
5. COSMOS 非复发者随访采用补充表 RFS/censor 月数；Chen/ColonAiQ 非复发者以可确认的随访采样时间判断12月结局。早期截尾或竞争事件不作为对照。

## 文件

系数、指标、增量、转移、分子记忆、缺失与时间窗审计表位于 `outputs/tables/`；指定图形位于 `outputs/figures/`。患者级轨迹及真实病例示例只写入 `outputs/local_only/`，不得提交到 GitHub。
"""
    (reports/"galaxy_development_external_validation_zh.md").write_text(report,encoding="utf-8")

    upreport=f"""# Same-horizon L3→L6 风险更新

## 目标

在同一组 GALAXY 患者中，两个模型预测完全相同的结局 **(200天, 365.24天] 临床复发**：before 只用临床信息及 L3/MRD-window history；after 增加 L6 ctDNA 与截至L6的既往阳性记忆。模型形式均为 ridge logistic、λ=1；同一病人风险集；组间风险使用1000次患者级 bootstrap 的 out-of-bag (OOB) 预测，以降低开发集表观拟合偏倚。

## 3月→6月状态分组

| 转移 | N | 复发 | 实际复发率 (95% Wilson CI) | before median risk | after median risk | Δrisk median (bootstrap 95% CI) |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(group_lines)}

总体OOB性能也在相同患者和相同结局窗比较：

| 指标 | L3信息 OOB | L6信息 OOB | Δ (95% paired bootstrap CI) |
|---|---:|---:|---:|
{chr(10).join(perf_lines)}

`Δrisk = risk_after - risk_before`。负→阳应上调、阳→阴应下调的方向以估计结果为准，不预设必须成立。稀疏的转换/清除亚组需谨慎解读。完整患者级风险列只保存在 `outputs/local_only/patient_dynamic_predictions.csv`。真实代表性患者SVG也只在 `outputs/local_only/`。

## 解释

这是固定未来结局窗口上的纵向风险更新，不是比较 L3→12 与 L6→12 两个不同长度预测窗。它支持回答“同一患者在看到新L6 MRD后，面对同一结局窗口的风险估计如何变化”。由于两模型均在GALAXY开发，OOB预测用于组别汇总，独立可迁移性由冻结模型的Chen/COSMOS/ColonAiQ验证另行判断。
"""
    (reports/"same_horizon_dynamic_update_zh.md").write_text(upreport,encoding="utf-8")

    dynreport=f"""# Dynamic survival sensitivity 与时间泄漏审计

本轮完成了固定结局窗口的 same-horizon 更新，不将其称作连续时间 joint model 或 time-dependent Cox。

GALAXY 补充表提供DFS天数和研究标注的3月、6月 ctDNA 状态，但没有患者级3/6月采血日期。3月窗口为70–112天、6月窗口为160–200天；无法诚实地将患者级时间变 covariate 放在确切采血时点，也无法评估采血前后事件顺序。因此本轮**不拟合 time-dependent Cox 或 landmark supermodel**，不插值采血日。严格≤182.6天的 Galaxy 分析同样不可识别。主模型改用200天共同风险起点的保守 landmark，并明确为 L6-window-to-12m 近似。

## 主要泄漏检查

- L6 MRD 输入只来自 GALAXY 表中的研究定义6月列，外部队列只取160–200天内实测状态。
- 主风险集排除第200天及以前复发者；外部队列的采样日有记录，所选状态必须早于复发事件。
- 目标窗口为(200,365.24]天；复发早于或等于200天的患者不进入风险集；复发晚于365.24天只作为12月前无复发对照；非复发对照须有至少365.24天可确认随访。
- 未观测 ctDNA 不编码为阴性；临床结局早期截尾或随访不足者排除于拟合/验证集。
- 外部队列未用于特征选择、惩罚参数选择、阈值选择或重校准。
- 患者级同一结局窗风险更新用同一人群与配对OOB风险，避免将不同预测时距的风险直接相减。

## 留待未来院内数据

要估计真正患者特异的 `P(T_recurrence > u | MRD history through t)`，需院内队列同时保留确切 surgery date、每次采血日期、复发日期及末次临床无复发评估日期；届时优先用 landmark supermodel 或 time-dependent Cox，并报告 time-dependent AUC、C-index、Brier/IBS 与校准。当前公开GALAXY表的时间精度不支持这些指标，因此本轮不报告。
"""
    (reports/"dynamic_model_sensitivity_zh.md").write_text(dynreport,encoding="utf-8")

    # Compact stdout enables a quick, independent review of denominators after the run.
    print(json.dumps({"origin_day":ORIGIN_DAY,"development_n":len(galaxy),"development_events":int(sum(r["y"] for r in galaxy)),
        "transition_complete_cohorts":{c:{"n":len(rs),"events":int(sum(r["y"] for r in rs))} for c,rs in rows_by_cohort.items()},
        "validation_cohorts":{c:{"n":len(rs),"events":int(sum(r["y"] for r in rs))} for c,rs in validation_by_cohort.items()},
        "updates":update_rows,"memory":memory_rows,"bootstrap_successes":args.bootstrap,
        "site_available_in_galaxy":False},ensure_ascii=False))


if __name__ == "__main__":
    main()
