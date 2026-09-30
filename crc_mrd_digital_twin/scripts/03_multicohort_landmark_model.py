from __future__ import annotations

"""Multi-cohort, 12-month postoperative CRC recurrence landmark modeling.

Reads the de-identified public supplements already stored under data_raw/.
Writes aggregate model diagnostics to outputs/tables and a Chinese report.
No patient-level analytic file is written to the repository.
"""
from collections import defaultdict
from pathlib import Path
import csv
import math
import re
import json

import numpy as np
from docx import Document
from openpyxl import load_workbook
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw"
OUT = ROOT / "outputs" / "tables"
REPORTS = ROOT / "reports"
OUT.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)
HORIZON = 12.0
LANDMARKS = (3.0, 6.0)
RIDGE_LAMBDA = 1.0


def num(value):
    try:
        if value is None or str(value).strip() in ("", "NA", "n/a", "N/A"):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def status01(value):
    s = str(value).strip().lower()
    if s in ("positive", "1", "detected"):
        return 1
    if s in ("negative", "0", "not detected"):
        return 0
    if isinstance(value, (int, float)) and value in (0, 1):
        return int(value)
    return None


def stage_high(stage):
    s = str(stage or "").strip().upper()
    return int(s.startswith("III") or s.startswith("IV"))


def rectal(site):
    return int("rect" in str(site or "").strip().lower())


def outcome(event_time, followup_time):
    if event_time is not None:
        return int(event_time <= HORIZON)
    if followup_time is not None and followup_time >= HORIZON:
        return 0
    return None


def make_patient(cohort, pid, event_time, followup_time, stage, site, histories, assay):
    return {
        "cohort": cohort,
        "patient_id": str(pid),
        "event_time": event_time,
        "followup_time": followup_time,
        "y12": outcome(event_time, followup_time),
        "stage_high": stage_high(stage),
        "rectal_primary": rectal(site),
        "assay": assay,
        "history": sorted(
            [(float(t), int(s)) for t, s in histories if t is not None and s in (0, 1)],
            key=lambda x: x[0],
        ),
    }


def read_chen():
    cohort = "Chen"
    base_table = Document(RAW / "13045_2021_1089_MOESM1_ESM.docx").tables[0]
    baseline = {}
    for row in base_table.rows[2:]:
        x = [c.text.strip() for c in row.cells]
        if not x or not x[0].startswith("P"):
            continue
        baseline[x[0]] = {
            "stage": x[4], "site": x[3], "event_time": num(x[11]),
        }

    wb = load_workbook(RAW / "13045_2021_1089_MOESM7_ESM.xlsx", read_only=True, data_only=True)
    ws = wb.active
    ctdna = defaultdict(list)
    last_observation = defaultdict(list)
    for row in ws.iter_rows(min_row=7, values_only=True):
        pid, t, point, raw_status = row[:4]
        if not pid:
            continue
        t = num(t)
        s = status01(raw_status)
        rank_match = re.search(r"(\d+)", str(point or ""))
        rank = int(rank_match.group(1)) if rank_match else 0
        if t is not None and rank >= 2:
            last_observation[pid].append(t)
            if s is not None:
                ctdna[pid].append((t, s))

    # The report's conservative endpoint confirmation also accepts postoperative CEA follow-up.
    ctab = load_workbook(RAW / "13045_2021_1089_MOESM9_ESM.xlsx", read_only=True, data_only=True).active
    for row in ctab.iter_rows(min_row=6, values_only=True):
        pid, t = row[:2]
        if pid and (t := num(t)) is not None and t > 0:
            last_observation[pid].append(t)

    patients = []
    for pid, b in baseline.items():
        follow = max(last_observation.get(pid, [0.0]))
        patients.append(make_patient(cohort, pid, b["event_time"], follow, b["stage"], b["site"], ctdna.get(pid, []), "tumor-informed binary ctDNA"))
    return patients


def read_coloniaiq():
    cohort = "ColonAiQ"
    reader = PdfReader(RAW / "mo_colonaiq_supplement1.pdf")
    baseline = {}
    # eTable 1: patient rows P1-P299 on printed pages 6-19.
    for page in reader.pages[5:19]:
        for line in (page.extract_text() or "").splitlines():
            tok = line.strip().split()
            if len(tok) >= 13 and re.fullmatch(r"P\d+", tok[0]) and tok[1] in ("Male", "Female"):
                pid = tok[0]
                event_time = num(tok[-1])
                baseline[pid] = {"stage": tok[4], "site": tok[5], "event_time": event_time}
    if len(baseline) != 299:
        raise RuntimeError(f"ColonAiQ eTable 1 解析人数应为299，实际为{len(baseline)}")

    histories = defaultdict(list)
    follow = defaultdict(list)
    # eTable 3: sample rows P1-P299 on printed pages 20-51; month 0 is preoperative.
    for page in reader.pages[19:51]:
        for line in (page.extract_text() or "").splitlines():
            tok = line.strip().split()
            if len(tok) < 4 or not re.fullmatch(r"P\d+", tok[0]):
                continue
            t = num(tok[1])
            s = status01(tok[-1])
            if t is None:
                continue
            if t > 0:
                follow[tok[0]].append(t)
                if s is not None:
                    histories[tok[0]].append((t, s))
    patients = []
    for pid, b in baseline.items():
        patients.append(make_patient(cohort, pid, b["event_time"], max(follow.get(pid, [0.0])), b["stage"], b["site"], histories.get(pid, []), "ColonAiQ methylation panel"))
    return patients


def read_cosmos():
    cohort = "COSMOS"
    ws = load_workbook(RAW / "cosmos_crc01_supp_table2.xlsx", read_only=True, data_only=True).active
    raw = [r for r in ws.iter_rows(min_row=4, max_col=23, values_only=True) if r[0]]
    schedule = (28.0 / 30.4375, 3.0, 6.0, 9.0, 12.0, 18.0, 24.0, 30.0)
    patients = []
    for r in raw:
        pid = str(r[0])
        resection = str(r[6] or "").strip()
        site = str(r[3] or "")
        stage = str(r[4] or "")
        et = num(r[13])
        recurrence_site = str(r[12] or "").strip().lower()
        event = recurrence_site not in ("", "n/a", "na", "none")
        event_time = et if event else None
        # Reconstruct the published 334-person longitudinal cohort: omit R2 and require
        # >=2 recorded scheduled cells, except a single available cell for relapse by month 6.
        recorded = sum(v is not None for v in r[14:22])
        if resection == "R2" or not (recorded >= 2 or (event and et is not None and et <= 6 and recorded >= 1)):
            continue
        hist = []
        follow_times = []
        for t, v in zip(schedule, r[14:22]):
            s = status01(v)
            if v is not None:
                follow_times.append(t)
            if s is not None:
                hist.append((t, s))
        # RFS is event time for recurrence and censor/follow-up time otherwise.
        followup = et if et is not None else None
        patients.append(make_patient(cohort, pid, event_time, followup, stage, site, hist, "tumor-informed binary ctDNA"))
    if len(patients) != 334:
        raise RuntimeError(f"COSMOS longitudinal analytic cohort 应为334，实际为{len(patients)}")
    return patients


def read_tie():
    cohort = "Tie"
    reader = PdfReader(RAW / "tie_stage3_supplement.pdf")
    rows = []
    current = None
    # eTable 3 spans printed pages 6-9. Join wrapped lines within each patient row.
    id_re = re.compile(r"^(LCRA?)\s+(\d+[a-z]?)\b", re.I)
    for page in reader.pages[5:9]:
        for line in (page.extract_text() or "").splitlines():
            s = line.strip()
            if id_re.match(s):
                if current:
                    rows.append(current)
                current = s
            elif current and s and not s.startswith(("©", "eTable", "Patient", "Tumor mutation", "Chemo", "Post-op", "Post-", "MAF", "Site(s)", "NA =")):
                current += " " + s
    if current:
        rows.append(current)

    patients = []
    for text in rows:
        im = id_re.match(text)
        if not im:
            continue
        pid = f"{im.group(1).upper()} {im.group(2)}"
        rest = text[im.end():]
        # First treatment token occurs after the variable-length mutation annotation.
        tm = re.search(r"\b(CAPOX|capecitabine|5FU|FOLFOX|Declined)\b", rest, re.I)
        if not tm:
            continue
        tail = rest[tm.end():].split()
        postop = next((status01(x) for x in tail if status01(x) is not None), None)
        et_matches = list(re.finditer(r"\b(Yes|No)\s+(\d+(?:\.\d+)?)\b", rest, re.I))
        if postop is None or not et_matches:
            continue
        em = et_matches[-1]
        event = em.group(1).lower() == "yes"
        time = float(em.group(2))
        # The supplement gives POM1 ctDNA collection at 4-10 weeks, without an exact date.
        # It is certainly before either landmark, so 1.5 months is only a sorting proxy.
        patients.append(make_patient(cohort, pid, time if event else None, time, "III", "Colon", [(1.5, postop)], "tumor-informed binary ctDNA"))
    if len(patients) != 96:
        raise RuntimeError(f"Tie eTable 3 解析人数应为96，实际为{len(patients)}")
    return patients


def feature_row(patient, landmark=None, with_mrd=False):
    row = [float(patient["stage_high"]), float(patient["rectal_primary"])]
    names = ["stage_III_IV", "rectal_primary"]
    if with_mrd:
        hist = [(t, s) for t, s in patient["history"] if t <= landmark]
        if not hist:
            raise ValueError("无 landmark 前 MRD 历史的患者进入了预测样本")
        row.extend([float(hist[-1][1]), float(any(s == 1 for _, s in hist[:-1]))])
        names.extend(["current_mrd_positive", "prior_mrd_positive_before_current"])
    return names, row


def sigmoid(z):
    z = np.clip(z, -35, 35)
    return 1.0 / (1.0 + np.exp(-z))


def fit_logit(x, y, penalty=RIDGE_LAMBDA, max_iter=100):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    p = x.shape[1]
    beta = np.zeros(p, dtype=float)
    beta[0] = math.log((y.mean() + 1e-4) / (1 - y.mean() + 1e-4))
    pen = np.eye(p) * penalty
    pen[0, 0] = 0.0
    for _ in range(max_iter):
        prob = sigmoid(x @ beta)
        w = np.maximum(prob * (1 - prob), 1e-8)
        grad = x.T @ (y - prob) - pen @ beta
        hess = (x.T * w) @ x + pen + np.eye(p) * 1e-9
        step = np.linalg.solve(hess, grad)
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break
    return beta


def rank_auc(y, p):
    y = np.asarray(y, dtype=int); p = np.asarray(p, dtype=float)
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(p, kind="mergesort")
    ranks = np.empty(len(p), dtype=float)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and p[order[j]] == p[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + 1 + j) / 2.0
        i = j
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def average_precision(y, p):
    y = np.asarray(y, dtype=int); p = np.asarray(p, dtype=float)
    n1 = int(y.sum())
    if n1 == 0:
        return float("nan")
    order = np.argsort(-p, kind="mergesort")
    ys = y[order]
    precision = np.cumsum(ys) / np.arange(1, len(ys) + 1)
    return float((precision * ys).sum() / n1)


def calibration(y, p):
    lp = np.log(np.clip(p, 1e-6, 1 - 1e-6) / np.clip(1 - p, 1e-6, 1))
    # With constant predictions (for example, Tie's constant clinical features) the
    # calibration slope is not identifiable. Perfect/semi-separation is also unstable.
    if float(np.std(lp)) < 1e-8:
        return float("nan"), float("nan")
    x = np.column_stack([np.ones(len(lp)), lp])
    try:
        b = fit_logit(x, y, penalty=1e-8, max_iter=200)
        if not np.all(np.isfinite(b)) or np.max(np.abs(b)) > 10:
            return float("nan"), float("nan")
        return float(b[0]), float(b[1])
    except (np.linalg.LinAlgError, FloatingPointError):
        return float("nan"), float("nan")


def score_metrics(y, p):
    y = np.asarray(y, dtype=int); p = np.asarray(p, dtype=float)
    ci, cs = calibration(y, p)
    return {
        "n": len(y), "events": int(y.sum()), "event_rate": float(y.mean()),
        "auc": rank_auc(y, p), "average_precision": average_precision(y, p),
        "brier": float(np.mean((y - p) ** 2)), "mean_predicted_risk": float(p.mean()),
        "calibration_intercept": ci, "calibration_slope": cs,
    }


def dca_rows(model, y, p):
    n = len(y); prevalence = float(np.mean(y)); rows = []
    for threshold in (0.05, 0.10, 0.15, 0.20, 0.30):
        treat = p >= threshold
        tp = int(np.sum(treat & (y == 1)))
        fp = int(np.sum(treat & (y == 0)))
        weight = threshold / (1 - threshold)
        nb = tp / n - fp / n * weight
        rows.append({"model": model, "threshold": threshold, "n": n, "events": int(y.sum()), "net_benefit_model": nb,
                     "net_benefit_treat_all": prevalence - (1 - prevalence) * weight,
                     "net_benefit_treat_none": 0.0, "n_flagged": int(treat.sum())})
    return rows


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def as_design(patients, landmark, with_mrd):
    names, _ = feature_row(patients[0], landmark, with_mrd)
    x = np.vstack([feature_row(p, landmark, with_mrd)[1] for p in patients])
    y = np.array([1 if p["event_time"] is not None and landmark < p["event_time"] <= HORIZON else 0 for p in patients], dtype=int)
    return names, x, y


def run_loco(label, patients, landmark, with_mrd):
    names, x, y = as_design(patients, landmark, with_mrd)
    pred = np.full(len(patients), np.nan)
    fold_metrics = []
    for heldout in sorted({p["cohort"] for p in patients}):
        test = np.array([p["cohort"] == heldout for p in patients])
        train = ~test
        xx = np.column_stack([np.ones(train.sum()), x[train]])
        beta = fit_logit(xx, y[train])
        pp = sigmoid(np.column_stack([np.ones(test.sum()), x[test]]) @ beta)
        pred[test] = pp
        met = score_metrics(y[test], pp)
        fold_metrics.append({"model": label, "held_out_cohort": heldout, **met, "ridge_lambda": RIDGE_LAMBDA,
                             "features": ";".join(names)})
    if np.isnan(pred).any():
        raise RuntimeError(f"LOCO prediction missing: {label}")
    return {"model": label, "landmark_month": "baseline" if landmark == 0 else int(landmark),
            "with_mrd": with_mrd, "patients": patients, "y": y, "p": pred,
            "pooled": score_metrics(y, pred), "folds": fold_metrics, "feature_names": names}


def main():
    # Primary analysis excludes Tie at the user request; parser remains for future sensitivity work.
    cohorts = [read_chen(), read_coloniaiq(), read_cosmos()]
    all_patients = [p for c in cohorts for p in c]
    cohort_names = ["Chen", "ColonAiQ", "COSMOS"]
    if len(all_patients) == 0:
        raise RuntimeError("未读取到任何患者")
    # Ensure each cohort has unique IDs and the three primary cohorts are represented.
    for c, name in zip(cohorts, cohort_names):
        ids = [p["patient_id"] for p in c]
        if len(ids) != len(set(ids)):
            raise RuntimeError(f"{name} 存在重复 patient ID")

    count_rows = []
    for c, name in zip(cohorts, cohort_names):
        known = [p for p in c if p["y12"] is not None]
        row = {"cohort": name, "raw_or_reconstructed_n": len(c), "y12_known_n": len(known),
               "y12_events": sum(p["y12"] == 1 for p in known),
               "y12_confirmed_final_nonrecurrence_controls": sum(p["event_time"] is None and p["followup_time"] is not None and p["followup_time"] >= 12 for p in known),
               "y12_late_recurrence_controls": sum(p["event_time"] is not None and p["event_time"] > 12 for p in known),
               "y12_unknown": len(c) - len(known)}
        for L in LANDMARKS:
            eligible = [p for p in known if (p["event_time"] is None or p["event_time"] > L)
                        and any(t <= L for t, _ in p["history"])]
            row[f"L{int(L)}_n"] = len(eligible)
            row[f"L{int(L)}_events"] = sum(p["event_time"] is not None and L < p["event_time"] <= 12 for p in eligible)
            row[f"L{int(L)}_MRD_history_n"] = sum(any(t <= L for t, _ in p["history"]) for p in c)
        count_rows.append(row)

    # Strictly check the previously audited denominators before fitting anything.
    expected = {
        "Chen": (181, 15, 181, 15, 179, 13),
        "ColonAiQ": (178, 36, 147, 26, 137, 12),
        "COSMOS": (332, 21, 330, 19, 326, 15),
    }
    for row in count_rows:
        want = expected[row["cohort"]]
        got = (row["y12_known_n"], row["y12_events"], row["L3_n"], row["L3_events"], row["L6_n"], row["L6_events"])
        if got != want:
            raise RuntimeError(f"{row['cohort']} 样本/事件数与审计不一致：实际{got}，预期{want}")
    if sum(r["y12_known_n"] for r in count_rows) != 691 or sum(r["y12_events"] for r in count_rows) != 72:
        raise RuntimeError("三队列主分析 Y12 样本数应为 691/72")
    if sum(r["L3_n"] for r in count_rows) != 658 or sum(r["L3_events"] for r in count_rows) != 60:
        raise RuntimeError("三队列 3月风险集应为 658/60")
    if sum(r["L6_n"] for r in count_rows) != 642 or sum(r["L6_events"] for r in count_rows) != 40:
        raise RuntimeError("三队列 6月风险集应为 642/40")
    write_csv(OUT / "multicohort_model_sample_counts.csv", count_rows)

    known_all = [p for p in all_patients if p["y12"] is not None]
    results = []
    # Baseline all-known-outcome cohort (Y12 case vs control).
    results.append(run_loco("M0_baseline_clinical", known_all, 0.0, False))
    # Matched comparisons within each landmark risk set.
    for L in LANDMARKS:
        eligible = [p for p in known_all if (p["event_time"] is None or p["event_time"] > L)
                    and any(t <= L for t, _ in p["history"])]
        results.append(run_loco(f"L{int(L)}_clinical_only", eligible, L, False))
        results.append(run_loco(f"L{int(L)}_clinical_plus_MRD_history", eligible, L, True))

    summary_rows, by_cohort_rows, decision_rows, prediction_rows = [], [], [], []
    for r in results:
        summary_rows.append({"model": r["model"], "landmark_month": r["landmark_month"], **r["pooled"],
                             "ridge_lambda": RIDGE_LAMBDA, "features": ";".join(r["feature_names"]),
                             "validation": "leave-one-cohort-out; held-out cohort effect not used"})
        by_cohort_rows.extend(r["folds"])
        decision_rows.extend(dca_rows(r["model"], r["y"], r["p"]))
        for p, y, pred in zip(r["patients"], r["y"], r["p"]):
            prediction_rows.append({"model": r["model"], "cohort": p["cohort"], "patient_id": p["patient_id"], "landmark_month": r["landmark_month"], "outcome": int(y), "predicted_risk": float(pred)})
    write_csv(OUT / "multicohort_loco_metrics.csv", summary_rows)
    write_csv(OUT / "multicohort_loco_by_heldout_cohort.csv", by_cohort_rows)
    write_csv(OUT / "multicohort_decision_curve.csv", decision_rows)
    # Patient-level cross-fitted predictions stay on the E: drive and are not committed.
    pred_dir = ROOT / "outputs" / "local_only"
    pred_dir.mkdir(parents=True, exist_ok=True)
    write_csv(pred_dir / "multicohort_loco_predictions.csv", prediction_rows)

    # Final all-cohort fits: no-cohort-effect transport model plus cohort fixed intercepts
    # for use only within these four known sources. These are development coefficients,
    # not an independently validated clinical score.
    coefficient_rows = []
    for r in results:
        patients = r["patients"]; y = r["y"]
        names, x, _ = as_design(patients, 0.0 if r["landmark_month"] == "baseline" else float(r["landmark_month"]), r["with_mrd"])
        beta = fit_logit(np.column_stack([np.ones(len(x)), x]), y)
        for name, val in zip(["intercept", *names], beta):
            coefficient_rows.append({"model": r["model"], "fit_scope": "pooled_no_cohort_effect", "coefficient": name, "estimate": float(val)})
        cohorts_present = sorted({p["cohort"] for p in patients})
        reference = cohorts_present[0]
        cohort_cols = [c for c in cohorts_present if c != reference]
        dummies = np.array([[float(p["cohort"] == c) for c in cohort_cols] for p in patients]) if cohort_cols else np.zeros((len(patients), 0))
        xcf = np.column_stack([x, dummies])
        bcf = fit_logit(np.column_stack([np.ones(len(xcf)), xcf]), y)
        cf_names = ["intercept", *names, *[f"cohort[{c}]_vs_{reference}" for c in cohort_cols]]
        for name, val in zip(cf_names, bcf):
            coefficient_rows.append({"model": r["model"], "fit_scope": f"known_cohorts_fixed_intercepts_ref_{reference}", "coefficient": name, "estimate": float(val)})
    write_csv(OUT / "multicohort_development_coefficients.csv", coefficient_rows)

    # Chinese report with aggregate findings and statistical limits.
    count_md = "\n".join(
        f"| {r['cohort']} | {r['y12_known_n']} | {r['y12_events']} | {r['L3_n']} | {r['L3_events']} | {r['L6_n']} | {r['L6_events']} |"
        for r in count_rows)
    metric_md = "\n".join(
        f"| {r['model']} | {r['n']} | {r['events']} | {r['auc']:.3f} | {r['average_precision']:.3f} | {r['brier']:.3f} | {r['mean_predicted_risk']:.3f} | {r['calibration_intercept']:.3f} | {r['calibration_slope']:.3f} |"
        for r in summary_rows)
    loco_md = "\n".join(
        f"| {r['model']} | {r['held_out_cohort']} | {r['n']} | {r['events']} | {r['auc']:.3f} | {r['brier']:.3f} |"
        for r in by_cohort_rows)
    lookup = {r["model"]: r for r in summary_rows}
    paired_md = "\n".join([
        f"| 3月风险集：MRD vs 临床 | AUC {lookup['L3_clinical_plus_MRD_history']['auc']-lookup['L3_clinical_only']['auc']:+.3f} | AP {lookup['L3_clinical_plus_MRD_history']['average_precision']-lookup['L3_clinical_only']['average_precision']:+.3f} | Brier {lookup['L3_clinical_plus_MRD_history']['brier']-lookup['L3_clinical_only']['brier']:+.3f} |",
        f"| 6月风险集：MRD vs 临床 | AUC {lookup['L6_clinical_plus_MRD_history']['auc']-lookup['L6_clinical_only']['auc']:+.3f} | AP {lookup['L6_clinical_plus_MRD_history']['average_precision']-lookup['L6_clinical_only']['average_precision']:+.3f} | Brier {lookup['L6_clinical_plus_MRD_history']['brier']-lookup['L6_clinical_only']['brier']:+.3f} |",
    ])
    report = f'''# 多队列术后 CRC 12 个月复发动态预测：主分析更新（三队列）

**日期：**2026-09-30  
**分析根目录：**`{ROOT}`  
**结局：**术后12个月内临床复发（Y12）；landmark 模型在该时点仍无复发者中预测至12个月。

## 分析队列核对

| 队列 | Y12可判定 N | Y12事件 | 3月风险集 N | 3月后至12月事件 | 6月风险集 N | 6月后至12月事件 |
|---|---:|---:|---:|---:|---:|---:|
{count_md}
| **合计** | **691** | **72** | **658** | **60** | **642** | **40** |

脚本在拟合前逐队列核对了既有审计计数；若人数或事件数不符会中止。Henriksen 因没有逐人复发日期未进入 Y12 模型；Tie 按本次指定从主分析和 LOCO 中移除，保留在原始审计材料中。

## 模型和验证

采用带 L2 收缩的低维 logistic 模型（惩罚参数 λ={RIDGE_LAMBDA}），避免在72个 Y12 事件上拟合过多参数。共同临床特征为 III/IV期与直肠原发；不同 assay 统一为二分类 MRD。主模型暂不纳入 ACT：补充表未给精确治疗时间，无法证明它在术前或 landmark 时已可用。

- **M0 baseline clinical：**所有 Y12 可判定者，预测 ≤12月复发。
- **L3 clinical only vs clinical + MRD history：**在术后3月仍无复发且截至该时点有可判读 MRD 历史者中，预测 (3,12] 月复发。
- **L6 clinical only vs clinical + MRD history：**同理预测 (6,12] 月复发。
- MRD 摘要为截至 landmark 最新状态，以及其前一次及更早样本中是否曾阳性。缺失采样绝不赋值为阴性；不使用 landmark 后/复发后的状态。
- 采用 leave-one-cohort-out (LOCO)：每轮整队列留出。留出队列的来源效应不可从训练折估计，因此 LOCO 传输模型不将 cohort ID 当预测器。另提供含 cohort 固定截距的全队列系数，仅适用于这三个已知来源；未来新医院需要独立校准。

不同 landmark 的风险集和预测结局不同，故 M0 全队列 AUC 与 L3/L6 AUC 不作直接的显著性或增益比较。同一 landmark 内，clinical-only 与 clinical+MRD 使用完全相同的患者和标签，描述性差异如下；尚未计算正式差异检验或置信区间。

| 配对比较 | ΔAUC | ΔAP | ΔBrier（负值较好）|
|---|---:|---:|---:|
{paired_md}

## LOCO 汇总指标

AUC 为 ROC AUC；AP 为 average precision；校准截距/斜率来自交叉拟合预测的 logistic 校准回归。

| 模型 | N | 事件 | AUC | AP | Brier | 平均预测风险 | 校准截距 | 校准斜率 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{metric_md}

## 各留出队列表现

| 模型 | 留出队列 | N | 事件 | AUC | Brier |
|---|---|---:|---:|---:|---:|
{loco_md}

小队列的事件数有限，单个留出队列指标会有较大抽样不确定性；本版为方法学开发结果，不能视为临床验证。

## 决策曲线

`outputs/tables/multicohort_decision_curve.csv` 提供风险阈值 5%、10%、15%、20%、30% 下的模型净获益，并与 treat-all / treat-none 基准比较。决策曲线依赖风险校准；对跨队列偏移明显的模型，应先校准再讨论临床阈值。

## 关键限制与解读

1. 可判定病例来自公开补充表的结局/末次随访字段，并沿用已审计的保守判定；结局未知者不作为阴性。
2. ColonAiQ 末次术后 ctDNA 样本仅作为无复发患者12个月观察证据的保守代理；它不等同于正式临床随访时间。
3. COSMOS 334人分析集按补充表重建：排除 R2，并按已发表纵向分析样本覆盖规则重建。该规则涉及纵向样本可用性，结果可能受观察/入组机制影响。
4. Tie 已按本次要求移出主分析。其公开表仅有 POM1 ctDNA、没有精确抽血日或6月样本，后续可作为单独的早期 MRD 敏感性分析。
5. 主分析仅含3个队列。ColonAiQ 是唯一的甲基化 assay 队列，因而 assay 与 cohort 完全混杂，不能独立估计 assay 效应。LOCO 是压力测试，不等同于独立前瞻性外部验证；AUC/校准估计仍需不确定区间和院内独立队列复核。
6. 单个留出队列事件数较少，且部分预测在该队列内几乎为常数；遇到校准回归不可识别或近似分离时，逐队列校准截距/斜率留空，不将不稳定极值作为结果。
7. 本版使用小型 ridge logistic 而不是复杂机器学习；当前主要价值是确认统一定义、可复现提取与跨队列验证流程是否成立。

## 生成文件

- `scripts/03_multicohort_landmark_model.py`：完整解析、风险集构建、LOCO、校准、Brier、AUC/AP、决策曲线脚本。
- `outputs/tables/multicohort_model_sample_counts.csv`：逐队列聚合人数核对。
- `outputs/tables/multicohort_loco_metrics.csv`：LOCO 汇总指标。
- `outputs/tables/multicohort_loco_by_heldout_cohort.csv`：逐留出队列指标。
- `outputs/tables/multicohort_decision_curve.csv`：净获益数据。
- `outputs/tables/multicohort_development_coefficients.csv`：全队列开发系数。
- 患者级交叉拟合预测仅保存在本机 `outputs/local_only/`，不会纳入版本库。
'''
    (REPORTS / "multicohort_12m_landmark_model_zh.md").write_text(report, encoding="utf-8")
    print(json.dumps({"cohort_counts": count_rows, "metrics": summary_rows, "report": str(REPORTS / "multicohort_12m_landmark_model_zh.md")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
