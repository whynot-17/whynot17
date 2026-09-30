from __future__ import annotations

"""Summarize frozen L3/L6 cross-fitted predictions as longitudinal MRD states.

This script does not fit or alter any prediction model. It reconstructs state
labels from the original public supplement parsers, checks landmark leakage,
joins the saved LOCO predictions, and writes aggregate outputs plus local-only
patient-level trajectories and examples.
"""

import csv
import html
import importlib.util
import itertools
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(r"E:\crc_mrd_digital_twin")
PRED_PATH = ROOT / "outputs" / "local_only" / "landmark_loco_patient_predictions.csv"
SCRIPT_DIR = ROOT / "scripts"
TABLE_DIR = ROOT / "outputs" / "tables"
FIG_DIR = ROOT / "outputs" / "figures"
LOCAL_DIR = ROOT / "outputs" / "local_only"
REPORT_DIR = ROOT / "reports"
BOOTSTRAP_REPS = 2000
SEED = 20260930
LANDMARKS = (3, 6)
COHORTS = ("Chen", "ColonAiQ", "COSMOS")

for d in (TABLE_DIR, FIG_DIR, LOCAL_DIR, REPORT_DIR):
    d.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location("crc_mrd_original_parser", SCRIPT_DIR / "03_multicohort_landmark_model.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

STATE_ZH = {
    "single_negative": "单次阴性",
    "single_positive": "单次阳性",
    "persistent_negative": "持续阴性",
    "persistent_positive": "持续阳性",
    "positive_to_negative_clearance": "阳性后转阴（清除）",
    "negative_to_positive_conversion": "阴性后转阳（转换）",
    "intermittent_mixed": "间歇/混合",
}
STATUS_ZH = {"negative_to_negative": "阴→阴", "negative_to_positive": "阴→阳",
             "positive_to_negative": "阳→阴", "positive_to_positive": "阳→阳"}
STATE_ORDER = list(STATE_ZH)
STATE_COLORS = {
    "single_negative": "#72b7a1", "single_positive": "#e98989",
    "persistent_negative": "#238b67", "persistent_positive": "#b82e38",
    "positive_to_negative_clearance": "#4c91c6", "negative_to_positive_conversion": "#e67e22",
    "intermittent_mixed": "#9276b5",
}


def write_csv(path: Path, rows):
    rows = list(rows)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def parse_float(x):
    if x is None or str(x).strip() == "":
        return None
    return float(x)


def state_pattern(hist):
    statuses = [int(s) for _, s in hist]
    if len(statuses) == 1:
        return "single_positive" if statuses[0] else "single_negative"
    if all(s == 0 for s in statuses):
        return "persistent_negative"
    if all(s == 1 for s in statuses):
        return "persistent_positive"
    if statuses[-1] == 0 and any(s == 1 for s in statuses[:-1]):
        return "positive_to_negative_clearance"
    if statuses[-1] == 1 and all(s == 0 for s in statuses[:-1]):
        return "negative_to_positive_conversion"
    return "intermittent_mixed"


def status_transition(a, b):
    return ("positive" if a else "negative") + "_to_" + ("positive" if b else "negative")


def qsummary(values):
    a = np.asarray(list(values), dtype=float)
    if not len(a):
        return {"mean": "", "median": "", "q1": "", "q3": ""}
    return {"mean": float(np.mean(a)), "median": float(np.median(a)),
            "q1": float(np.quantile(a, .25)), "q3": float(np.quantile(a, .75))}


def bootstrap_median_ci(values, seed, reps=BOOTSTRAP_REPS):
    a = np.asarray(list(values), dtype=float)
    if len(a) < 5:
        return "", ""
    rng = np.random.default_rng(seed)
    vals = np.empty(reps, dtype=float)
    for i in range(reps):
        vals[i] = np.median(a[rng.integers(0, len(a), len(a))])
    return tuple(float(x) for x in np.quantile(vals, [.025, .975]))


def bootstrap_median_difference(a, b, seed, reps=BOOTSTRAP_REPS):
    """Median(a)-median(b), independent patient samples, stratified by group."""
    a, b = np.asarray(list(a), float), np.asarray(list(b), float)
    if min(len(a), len(b)) < 5:
        return "", ""
    rng = np.random.default_rng(seed)
    vals = np.empty(reps, float)
    for i in range(reps):
        aa = a[rng.integers(0, len(a), len(a))]
        bb = b[rng.integers(0, len(b), len(b))]
        vals[i] = np.median(aa) - np.median(bb)
    return tuple(float(x) for x in np.quantile(vals, [.025, .975]))


def bootstrap_paired_median_ci(deltas, seed, reps=BOOTSTRAP_REPS):
    return bootstrap_median_ci(deltas, seed, reps)


def wilcoxon_signed_rank(deltas):
    d = np.asarray(list(deltas), float)
    d = d[np.abs(d) > 1e-12]
    n = len(d)
    if n < 5:
        return n, "", "descriptive_only_n<5"
    absd = np.abs(d)
    order = np.argsort(absd, kind="mergesort")
    ranks = np.empty(n, float)
    i = 0
    while i < n:
        j = i + 1
        while j < n and absd[order[j]] == absd[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + 1 + j) / 2.0
        i = j
    wplus = float(ranks[d > 0].sum())
    total = float(ranks.sum())
    wmin = min(wplus, total - wplus)
    if n <= 20:
        extreme = 0
        for signs in itertools.product((0, 1), repeat=n):
            wp = float(np.dot(ranks, signs))
            if min(wp, total - wp) <= wmin + 1e-10:
                extreme += 1
        p = extreme / (2 ** n)
        method = "exact_two_sided_signed_rank"
    else:
        mean = total / 2.0
        sd = math.sqrt(float(np.dot(ranks, ranks)) / 4.0)
        z = max(0.0, abs(wplus - mean) - .5) / sd if sd else 0.0
        p = math.erfc(z / math.sqrt(2.0))
        if p == 0.0:
            p = "<1e-300"
        method = "normal_approximation_with_continuity_correction"
    return n, (p if isinstance(p, str) else float(min(1.0, p))), method


def svg_header(width, height, title, subtitle=""):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            '<rect width="100%" height="100%" fill="white"/>',
            '<style>text{font-family:Arial,"Microsoft YaHei",sans-serif;fill:#202124}.title{font-size:23px;font-weight:bold}.sub{font-size:14px;fill:#555}.lab{font-size:14px}.tiny{font-size:12px;fill:#555}</style>',
            f'<text x="{width/2}" y="34" class="title" text-anchor="middle">{html.escape(title)}</text>',
            (f'<text x="{width/2}" y="59" class="sub" text-anchor="middle">{html.escape(subtitle)}</text>' if subtitle else "")]


def svg_alluvial(transitions, path):
    W, H, top, bottom = 1500, 850, 125, 770
    rows = [r for r in transitions]
    left_tot = Counter(); right_tot = Counter(); flows = Counter()
    for r in rows:
        key = (r["state3"], r["state6"])
        flows[key] += 1; left_tot[key[0]] += 1; right_tot[key[1]] += 1
    left_order = [s for s in STATE_ORDER if left_tot[s]]
    right_order = [s for s in STATE_ORDER if right_tot[s]]
    scale = (bottom - top) / max(1, len(rows))
    left_pos = {}; right_pos = {}
    y = top
    for s in left_order:
        left_pos[s] = [y, y]; y += left_tot[s] * scale
    y = top
    for s in right_order:
        right_pos[s] = [y, y]; y += right_tot[s] * scale
    out = svg_header(W, H, "L3 至 L6 MRD 数字孪生状态转移", "节点宽度和连线宽度表示人数；患者级 ID 未展示。")
    x0, x1, nw = 300, 1200, 22
    for (a, b), n in sorted(flows.items(), key=lambda x: (STATE_ORDER.index(x[0][0]), STATE_ORDER.index(x[0][1]))):
        h = n * scale
        sy = left_pos[a][1]; ey = right_pos[b][1]
        left_pos[a][1] += h; right_pos[b][1] += h
        c1 = sy + h / 2; c2 = ey + h / 2
        out.append(f'<path d="M{x0+nw},{c1} C{x0+340},{c1} {x1-340},{c2} {x1},{c2}" fill="none" stroke="{STATE_COLORS[a]}" stroke-opacity="0.34" stroke-width="{max(h,0.55)}"/>')
    for s in left_order:
        y0 = left_pos[s][0]; h = left_tot[s] * scale
        out.append(f'<rect x="{x0}" y="{y0}" width="{nw}" height="{h}" fill="{STATE_COLORS[s]}"/>')
        out.append(f'<text x="{x0-12}" y="{y0+h/2+5}" text-anchor="end" class="lab">{STATE_ZH[s]} n={left_tot[s]}</text>')
    for s in right_order:
        y0 = right_pos[s][0]; h = right_tot[s] * scale
        out.append(f'<rect x="{x1}" y="{y0}" width="{nw}" height="{h}" fill="{STATE_COLORS[s]}"/>')
        out.append(f'<text x="{x1+34}" y="{y0+h/2+5}" class="lab">{STATE_ZH[s]} n={right_tot[s]}</text>')
    out.append(f'<text x="{x0+10}" y="95" class="lab" text-anchor="middle">L3，未来复发窗 (3,12]</text>')
    out.append(f'<text x="{x1+10}" y="95" class="lab" text-anchor="middle">L6，未来复发窗 (6,12]</text>')
    out.append("</svg>")
    path.write_text("\n".join(x for x in out if x), encoding="utf-8")


def svg_spaghetti(pairs, path):
    W, H, left, right, top, bottom = 1080, 690, 95, 1035, 100, 590
    out = svg_header(W, H, "L3 → L6 cross-fitted 预测风险轨迹", "同一患者的两个预测风险对应不同未来时间窗，Δ风险不应解释为纯 MRD 效应。")
    def yy(v): return bottom - max(0, min(1, v)) * (bottom - top)
    colors = {"negative_to_negative": "#27845c", "negative_to_positive": "#d66b24",
              "positive_to_negative": "#327eb7", "positive_to_positive": "#b43b45"}
    for k in range(6):
        v = k / 5; y = yy(v)
        out.append(f'<line x1="{left}" y1="{y}" x2="{right}" y2="{y}" stroke="#e2e5e9"/>')
        out.append(f'<text x="{left-12}" y="{y+5}" text-anchor="end" class="lab">{v:.1f}</text>')
    for r in pairs:
        x1, x2 = left + 180, right - 180
        color = colors[r["status_transition"]]
        out.append(f'<line x1="{x1}" y1="{yy(r["risk3"])}" x2="{x2}" y2="{yy(r["risk6"])}" stroke="{color}" stroke-opacity="0.13" stroke-width="1.2"/>')
    out.append(f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#333"/><line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#333"/>')
    out.append(f'<text x="{left+180}" y="{bottom+38}" class="lab" text-anchor="middle">L3：预测 (3,12] 月复发风险</text>')
    out.append(f'<text x="{right-180}" y="{bottom+38}" class="lab" text-anchor="middle">L6：预测 (6,12] 月复发风险</text>')
    x, y = left+190, 650
    for key, color in colors.items():
        out.append(f'<line x1="{x}" y1="{y-5}" x2="{x+24}" y2="{y-5}" stroke="{color}" stroke-width="4"/><text x="{x+30}" y="{y}" class="tiny">{STATUS_ZH[key]}</text>')
        x += 200
    out.append("</svg>")
    path.write_text("\n".join(z for z in out if z), encoding="utf-8")


def svg_delta(summary, path):
    rows = [r for r in summary if r["level"] == "current_status_transition" and int(r["n"]) >= 5]
    W, H, left, right, top, bottom = 1050, 520, 220, 980, 100, 405
    max_abs = max([abs(float(r["median_delta_risk"])) for r in rows] + [0.03])
    max_abs = max(max_abs * 1.25, 0.05)
    out = svg_header(W, H, "按 L3→L6 当前 MRD 转移分层的风险更新", "点为中位数；线为患者 bootstrap 95% CI；病例窗随 landmark 缩短。")
    def xx(v): return left + (float(v) + max_abs) / (2 * max_abs) * (right-left)
    out.append(f'<line x1="{xx(0)}" y1="{top}" x2="{xx(0)}" y2="{bottom}" stroke="#555" stroke-dasharray="5,4"/>')
    for i, r in enumerate(rows):
        y = top + i * ((bottom-top) / max(1, len(rows)-1))
        lo, hi, med = (float(r[k]) for k in ("delta_ci_lower", "delta_ci_upper", "median_delta_risk"))
        key = r["group"]
        color = {"negative_to_negative":"#27845c","negative_to_positive":"#d66b24","positive_to_negative":"#327eb7","positive_to_positive":"#b43b45"}[key]
        out.append(f'<text x="{left-15}" y="{y+5}" class="lab" text-anchor="end">{STATUS_ZH[key]} n={r["n"]}</text>')
        out.append(f'<line x1="{xx(lo)}" y1="{y}" x2="{xx(hi)}" y2="{y}" stroke="{color}" stroke-width="3"/><circle cx="{xx(med)}" cy="{y}" r="7" fill="{color}"/>')
    for j in range(5):
        v = -max_abs + j * max_abs / 2; x = xx(v)
        out.append(f'<line x1="{x}" y1="{bottom}" x2="{x}" y2="{bottom+6}" stroke="#333"/><text x="{x}" y="{bottom+25}" text-anchor="middle" class="tiny">{v:+.2f}</text>')
    out.append(f'<text x="{(left+right)/2}" y="{H-25}" text-anchor="middle" class="lab">Δ风险 = Risk6 − Risk3</text>')
    out.append("</svg>")
    path.write_text("\n".join(z for z in out if z), encoding="utf-8")


def svg_memory(memory_rows, path):
    rows = [r for r in memory_rows if r["analysis"] == "all_current_negative_group" and r["level"] == "group"]
    keep = [r for r in rows if r["group"] in ("never_prior_positive", "ever_prior_positive") and int(r["n"]) > 0]
    W, H, left, right, top, bottom = 1070, 550, 280, 1000, 105, 435
    maxv = max([float(r["median_risk"]) for r in keep] + [.05]) * 1.25
    out = svg_header(W, H, "当前 MRD 阴性患者中的分子记忆", "描述性比较：当前阴性但既往曾阳性 vs 截至当前从未阳性；不作因果解释。")
    def xx(v): return left + float(v) / maxv * (right-left)
    groups = ["never_prior_positive", "ever_prior_positive"]
    labels = {"never_prior_positive":"从未阳性（当前阴性）", "ever_prior_positive":"既往阳性、当前阴性"}
    for j, L in enumerate((3, 6)):
        yy = top + j * 155
        out.append(f'<text x="{left-22}" y="{yy+5}" class="lab" text-anchor="end">L{L} 预测窗 ({L},12]</text>')
        for i, g in enumerate(groups):
            row = next((r for r in keep if int(r["landmark"]) == L and r["group"] == g), None)
            if not row: continue
            y = yy + (i - .5) * 52
            lo, hi, med = row["median_ci_lower"], row["median_ci_upper"], row["median_risk"]
            if lo != "" and hi != "":
                out.append(f'<line x1="{xx(lo)}" y1="{y}" x2="{xx(hi)}" y2="{y}" stroke="#597da5" stroke-width="3"/>')
            out.append(f'<circle cx="{xx(med)}" cy="{y}" r="7" fill="#597da5"/><text x="{left-8}" y="{y+5}" text-anchor="end" class="tiny">{labels[g]} n={row["n"]}</text>')
    for j in range(5):
        v = maxv*j/4; x=xx(v)
        out.append(f'<line x1="{x}" y1="{bottom}" x2="{x}" y2="{bottom+6}" stroke="#333"/><text x="{x}" y="{bottom+25}" text-anchor="middle" class="tiny">{v:.2f}</text>')
    out.append(f'<text x="{(left+right)/2}" y="{H-30}" class="lab" text-anchor="middle">cross-fitted 预测风险，中位数及患者 bootstrap 95% CI</text>')
    out.append("</svg>")
    path.write_text("\n".join(z for z in out if z), encoding="utf-8")


def svg_cases(cases, path):
    W, H = 1280, max(270, 145 * len(cases) + 105)
    out = svg_header(W, H, "MRD 状态转移与风险更新：患者级示例（仅本地）", "真实公开队列记录；只显示截至 L6 的 MRD 与冻结模型预测，不含未来信息作预测特征。")
    x0, x1 = 370, 1175
    def xx(t): return x0 + max(0, min(6, float(t))) / 6 * (x1-x0)
    for i, r in enumerate(cases):
        y = 115 + 145 * i
        out.append(f'<text x="28" y="{y-22}" class="lab">示例{i+1}：{html.escape(r["label"])} | {html.escape(r["cohort"])} / {html.escape(r["patient_id"])}</text>')
        out.append(f'<text x="28" y="{y+1}" class="tiny">{STATE_ZH[r["state3"]]} → {STATE_ZH[r["state6"]]}；Risk3={r["risk3"]:.1%} → Risk6={r["risk6"]:.1%}；(6,12]复发={r["event6"]}</text>')
        out.append(f'<line x1="{x0}" y1="{y+18}" x2="{x1}" y2="{y+18}" stroke="#666" stroke-width="2"/>')
        for month in range(7):
            x=xx(month);out.append(f'<line x1="{x}" y1="{y+13}" x2="{x}" y2="{y+24}" stroke="#555"/><text x="{x}" y="{y+44}" class="tiny" text-anchor="middle">{month}月</text>')
        for t, s in r["history"]:
            col = "#bd363e" if s else "#25855e"
            out.append(f'<circle cx="{xx(t)}" cy="{y+18}" r="7" fill="{col}" stroke="white" stroke-width="1"/>')
    out.append("</svg>")
    path.write_text("\n".join(z for z in out if z), encoding="utf-8")


def main():
    if not PRED_PATH.exists():
        raise FileNotFoundError(f"missing frozen prediction file: {PRED_PATH}")
    pred = read_csv(PRED_PATH)
    pred = [r for r in pred if r["model"] == "clinical_plus_MRD" and r["landmark"] in ("L3", "L6") and r["cohort"] in COHORTS]
    # Reuse the original cohort parsers so status/date semantics exactly match model development.
    parsed = {"Chen": mod.read_chen(), "ColonAiQ": mod.read_coloniaiq(), "COSMOS": mod.read_cosmos()}
    patients = {(c, p["patient_id"]): p for c, ps in parsed.items() for p in ps}
    rows_by_lm = {3: [], 6: []}
    local_long = []
    for r in pred:
        L = int(r["landmark"][1:])
        key = (r["cohort"], r["patient_id"])
        if key not in patients:
            raise RuntimeError("saved model prediction could not be matched to parsed public cohort record")
        p = patients[key]
        hist = [(float(t), int(s)) for t, s in p["history"] if float(t) <= L and s in (0, 1)]
        if not hist:
            raise RuntimeError("a saved clinical_plus_MRD prediction has no MRD history at its landmark")
        max_t = max(t for t, _ in hist)
        if max_t > L + 1e-10:
            raise RuntimeError("landmark leakage guard failed")
        if abs(max_t - float(r["max_feature_mrd_time"])) > 1e-7:
            raise RuntimeError("reconstructed MRD timing differs from the frozen prediction input")
        event_time = parse_float(r["event_time"])
        expected_y = int(event_time is not None and L < event_time <= 12)
        if expected_y != int(r["outcome"]):
            raise RuntimeError(f"outcome-window guard failed at L{L}")
        state = state_pattern(hist)
        current = hist[-1][1]
        prior_positive = int(any(s == 1 for _, s in hist[:-1]))
        row = {"cohort": r["cohort"], "patient_id": r["patient_id"], "landmark": L,
               "state": state, "current_mrd": current, "ever_prior_positive": prior_positive,
               "history": hist, "n_tests": len(hist), "outcome": int(r["outcome"]),
               "event_time": event_time, "followup_time": parse_float(r["followup_time"]),
               "risk": float(r["predicted_risk"]), "stage_high": int(r["stage_III_IV"]),
               "rectal": int(r["rectal_primary"])}
        rows_by_lm[L].append(row)
        local_long.append({"cohort": row["cohort"], "patient_id": row["patient_id"], "landmark": L,
                           "outcome_window": f"({L},12] months", "digital_state": state,
                           "digital_state_zh": STATE_ZH[state], "current_mrd": current,
                           "ever_prior_positive_before_current": prior_positive, "n_mrd_tests": len(hist),
                           "mrd_history_month_status": ";".join(f"{t:g}:{s}" for t, s in hist),
                           "predicted_risk": row["risk"], "future_recurrence": row["outcome"],
                           "recurrence_time_months": event_time})
    for L in LANDMARKS:
        if not rows_by_lm[L]:
            raise RuntimeError(f"no L{L} cross-fitted MRD prediction rows")

    by_key_lm = {L: {(r["cohort"], r["patient_id"]): r for r in rows_by_lm[L]} for L in LANDMARKS}
    common = sorted(set(by_key_lm[3]) & set(by_key_lm[6]))
    pairs = []
    for key in common:
        a, b = by_key_lm[3][key], by_key_lm[6][key]
        st = status_transition(a["current_mrd"], b["current_mrd"])
        pairs.append({"cohort": key[0], "patient_id": key[1], "state3": a["state"], "state6": b["state"],
                      "current3": a["current_mrd"], "current6": b["current_mrd"], "status_transition": st,
                      "risk3": a["risk"], "risk6": b["risk"], "delta": b["risk"]-a["risk"],
                      "event6": b["outcome"], "event_time": b["event_time"],
                      "ever_prior_positive3": a["ever_prior_positive"], "ever_prior_positive6": b["ever_prior_positive"],
                      "hist3": a["history"], "hist6": b["history"], "n_tests3": a["n_tests"], "n_tests6": b["n_tests"]})
    # Some patients first have an MRD draw after month 3, so an L6 prediction may
    # exist without an L3 prediction. Paired changes use only the two-landmark
    # intersection; each landmark's state summary retains its full prediction set.
    if len(pairs) > len(rows_by_lm[6]):
        raise RuntimeError("L3/L6 paired set exceeds the L6 prediction risk set")

    local_transitions = []
    for r in pairs:
        local_transitions.append({"cohort": r["cohort"], "patient_id": r["patient_id"],
                                  "state_L3": r["state3"], "state_L6": r["state6"],
                                  "current_MRD_L3": r["current3"], "current_MRD_L6": r["current6"],
                                  "status_transition": r["status_transition"], "risk_L3": r["risk3"],
                                  "risk_L6": r["risk6"], "delta_risk_L6_minus_L3": r["delta"],
                                  "outcome_window_L6": "(6,12] months", "recurrence_6_to_12m": r["event6"],
                                  "recurrence_time_months": r["event_time"]})
    write_csv(LOCAL_DIR / "patient_digital_state_trajectory.csv", local_long)
    write_csv(LOCAL_DIR / "patient_state_transition_records.csv", local_transitions)

    state_counts = []
    state_recurrence = []
    for L in LANDMARKS:
        for cohort in COHORTS:
            groupc = [r for r in rows_by_lm[L] if r["cohort"] == cohort]
            for state in STATE_ORDER:
                g = [r for r in groupc if r["state"] == state]
                if not g: continue
                ev = sum(r["outcome"] for r in g)
                pos = sum(r["current_mrd"] == 1 for r in g)
                state_counts.append({"landmark": f"L{L}", "cohort": cohort, "digital_state": state,
                                     "digital_state_zh": STATE_ZH[state], "n": len(g), "current_positive_n": pos,
                                     "current_negative_n": len(g)-pos, "events": ev,
                                     "event_rate": ev/len(g), "outcome_window": f"({L},12] months"})
        for state in STATE_ORDER:
            g = [r for r in rows_by_lm[L] if r["state"] == state]
            if g:
                ev = sum(r["outcome"] for r in g)
                state_recurrence.append({"landmark": f"L{L}", "digital_state": state, "digital_state_zh": STATE_ZH[state],
                                         "n": len(g), "events": ev, "event_rate": ev/len(g),
                                         "outcome_window": f"({L},12] months"})
    write_csv(TABLE_DIR / "digital_state_counts.csv", state_counts)
    write_csv(TABLE_DIR / "digital_state_recurrence_summary.csv", state_recurrence)

    flow_counts = Counter((r["state3"], r["state6"]) for r in pairs)
    transition_matrix = []
    for (s3, s6), n in sorted(flow_counts.items(), key=lambda x:(STATE_ORDER.index(x[0][0]), STATE_ORDER.index(x[0][1]))):
        g = [r for r in pairs if r["state3"] == s3 and r["state6"] == s6]
        ev = sum(r["event6"] for r in g)
        risks3 = qsummary(r["risk3"] for r in g); risks6 = qsummary(r["risk6"] for r in g)
        delta = qsummary(r["delta"] for r in g)
        transition_matrix.append({"state_L3": s3, "state_L3_zh": STATE_ZH[s3], "state_L6": s6,
                                  "state_L6_zh": STATE_ZH[s6], "n": n, "events_6_to_12m": ev,
                                  "recurrence_rate_6_to_12m": ev/n, "median_risk_L3": risks3["median"],
                                  "median_risk_L6": risks6["median"], "median_delta_risk_L6_minus_L3": delta["median"],
                                  "delta_q1": delta["q1"], "delta_q3": delta["q3"], "risk_window_L3": "(3,12]",
                                  "risk_window_L6": "(6,12]"})
    write_csv(TABLE_DIR / "digital_state_transition_matrix.csv", transition_matrix)

    # Paired Risk3->Risk6 changes by detailed state transition and by current MRD switch.
    groups = [("all_paired", "all_patients", pairs)]
    for (s3, s6), _ in flow_counts.items():
        groups.append(("state_transition", s3 + "_to_" + s6,
                       [r for r in pairs if r["state3"] == s3 and r["state6"] == s6]))
    status_keys = ["negative_to_negative", "negative_to_positive", "positive_to_negative", "positive_to_positive"]
    for s in status_keys:
        groups.append(("current_status_transition", s, [r for r in pairs if r["status_transition"] == s]))
    update_rows = []
    paired_stats = []
    for ix, (level, label, g) in enumerate(groups):
        if not g: continue
        d = [r["delta"] for r in g]
        qs = qsummary(d)
        lo, hi = bootstrap_paired_median_ci(d, SEED + 100 + ix) if len(g) >= 5 else ("", "")
        nz, pvalue, method = wilcoxon_signed_rank(d)
        ev = sum(r["event6"] for r in g)
        update_rows.append({"level": level, "group": label, "group_zh": (STATUS_ZH.get(label, label.replace("_to_", " → "))),
                            "n": len(g), "events_6_to_12m": ev, "event_rate_6_to_12m": ev/len(g),
                            "median_risk_L3": float(np.median([r["risk3"] for r in g])),
                            "median_risk_L6": float(np.median([r["risk6"] for r in g])),
                            "median_delta_risk": qs["median"], "delta_q1": qs["q1"], "delta_q3": qs["q3"],
                            "delta_ci_lower": lo, "delta_ci_upper": hi,
                            "risk_window_note": "Risk3 predicts (3,12]; Risk6 predicts (6,12]. Delta includes time advancement/risk-set change."})
        paired_stats.append({"level": level, "group": label, "n": len(g), "n_nonzero_pairs": nz,
                             "median_delta_risk": qs["median"], "delta_ci_lower": lo, "delta_ci_upper": hi,
                             "wilcoxon_two_sided_p": pvalue, "wilcoxon_method": method,
                             "inference_note": "Exploratory, no multiplicity adjustment; different prediction windows."})
    write_csv(TABLE_DIR / "risk_update_by_transition.csv", update_rows)
    write_csv(TABLE_DIR / "paired_risk_update_statistics.csv", paired_stats)

    # Current-negative molecular-memory comparison at each landmark.
    memory_rows = []
    for L in LANDMARKS:
        allrows = rows_by_lm[L]
        neg = [r for r in allrows if r["current_mrd"] == 0]
        memory_groups = {
            "never_prior_positive": [r for r in neg if r["ever_prior_positive"] == 0],
            "ever_prior_positive": [r for r in neg if r["ever_prior_positive"] == 1],
        }
        for j, (group, g) in enumerate(memory_groups.items()):
            if not g: continue
            risk = [r["risk"] for r in g]; lo, hi = bootstrap_median_ci(risk, SEED+300+L*10+j) if len(g)>=5 else ("", "")
            ev=sum(r["outcome"] for r in g); qs=qsummary(risk)
            memory_rows.append({"landmark": L, "outcome_window": f"({L},12] months", "analysis":"all_current_negative_group",
                                "level":"group", "group":group, "n":len(g), "events":ev, "event_rate":ev/len(g),
                                "median_risk":qs["median"], "risk_q1":qs["q1"], "risk_q3":qs["q3"],
                                "median_ci_lower":lo, "median_ci_upper":hi})
        a=memory_groups["ever_prior_positive"]; b=memory_groups["never_prior_positive"]
        diff_ci=bootstrap_median_difference([r["risk"] for r in a], [r["risk"] for r in b], SEED+500+L) if a and b else ("", "")
        riskdiff=(float(np.median([r["risk"] for r in a]))-float(np.median([r["risk"] for r in b]))) if a and b else ""
        rate_a=sum(r["outcome"] for r in a)/len(a) if a else ""
        rate_b=sum(r["outcome"] for r in b)/len(b) if b else ""
        memory_rows.append({"landmark":L,"outcome_window":f"({L},12] months","analysis":"all_current_negative_group",
                            "level":"contrast_ever_minus_never","group":"ever_prior_positive_minus_never_prior_positive",
                            "n":len(a)+len(b),"events":"","event_rate_difference":(rate_a-rate_b if a and b else ""),
                            "median_risk_difference":riskdiff,"difference_ci_lower":diff_ci[0],"difference_ci_upper":diff_ci[1],
                            "interpretation":"Predictive association only; not causal evidence."})
        strict = {
            "persistent_negative": [r for r in neg if r["state"] == "persistent_negative"],
            "positive_to_negative_clearance": [r for r in neg if r["state"] == "positive_to_negative_clearance"],
        }
        for j,(group,g) in enumerate(strict.items()):
            if not g: continue
            risk=[r["risk"] for r in g];lo,hi=bootstrap_median_ci(risk,SEED+700+L*10+j) if len(g)>=5 else ("","")
            ev=sum(r["outcome"] for r in g);qs=qsummary(risk)
            memory_rows.append({"landmark":L,"outcome_window":f"({L},12] months","analysis":"strict_persistent_negative_vs_clearance",
                                "level":"group","group":group,"n":len(g),"events":ev,"event_rate":ev/len(g),
                                "median_risk":qs["median"],"risk_q1":qs["q1"],"risk_q3":qs["q3"],
                                "median_ci_lower":lo,"median_ci_upper":hi})
        a=strict["positive_to_negative_clearance"]; b=strict["persistent_negative"]
        dci=bootstrap_median_difference([r["risk"] for r in a],[r["risk"] for r in b],SEED+900+L) if a and b else ("","")
        rd=(sum(r["outcome"] for r in a)/len(a)-sum(r["outcome"] for r in b)/len(b)) if a and b else ""
        md=(float(np.median([r["risk"] for r in a]))-float(np.median([r["risk"] for r in b]))) if a and b else ""
        memory_rows.append({"landmark":L,"outcome_window":f"({L},12] months","analysis":"strict_persistent_negative_vs_clearance",
                            "level":"contrast_clearance_minus_persistent_negative","group":"clearance_minus_persistent_negative",
                            "n":len(a)+len(b),"events":"","event_rate_difference":rd,"median_risk_difference":md,
                            "difference_ci_lower":dci[0],"difference_ci_upper":dci[1],
                            "interpretation":"Small strata are descriptive; predictive association only."})
    write_csv(TABLE_DIR / "current_negative_memory_effect.csv", memory_rows)

    # Select at most four observed L3/L6 trajectories for a local-only illustration.
    case_specs = [
        ("持续阴性", lambda r: r["state3"]=="persistent_negative" and r["state6"]=="persistent_negative", lambda r:r["risk6"]),
        ("阳性后清除", lambda r: r["state6"]=="positive_to_negative_clearance", lambda r:abs(r["delta"])),
        ("L3阴性→L6阳性", lambda r:r["current3"]==0 and r["current6"]==1, lambda r:r["risk6"]-r["risk3"]),
        ("持续/新发阳性", lambda r:r["current3"]==1 and r["current6"]==1, lambda r:r["risk6"]),
    ]
    cases=[]; used=set()
    for label, test, score in case_specs:
        eligible=[r for r in pairs if test(r)]
        eligible=[r for r in eligible if (r["cohort"],r["patient_id"]) not in used]
        if not eligible: continue
        if label=="持续阴性": chosen=min(eligible,key=score)
        elif label=="阳性后清除": chosen=min(eligible,key=score)
        else: chosen=max(eligible,key=score)
        used.add((chosen["cohort"],chosen["patient_id"]))
        hist=[x for x in chosen["hist6"] if x[0]<=6]
        cases.append({"label":label,"cohort":chosen["cohort"],"patient_id":chosen["patient_id"],
                      "state3":chosen["state3"],"state6":chosen["state6"],"risk3":chosen["risk3"],
                      "risk6":chosen["risk6"],"event6":chosen["event6"],"history":hist})
    svg_cases(cases, LOCAL_DIR / "digital_twin_case_examples.svg")

    svg_alluvial(pairs, FIG_DIR / "digital_state_sankey.svg")
    svg_delta(update_rows, FIG_DIR / "risk_delta_by_transition.svg")
    svg_memory(memory_rows, FIG_DIR / "molecular_memory_effect.svg")
    svg_spaghetti(pairs, LOCAL_DIR / "risk_update_spaghetti.svg")

    # Aggregate report; no patient-level identifiers or records are written here.
    total3, total6 = len(rows_by_lm[3]), len(rows_by_lm[6])
    ev3=sum(r["outcome"] for r in rows_by_lm[3]); ev6=sum(r["outcome"] for r in rows_by_lm[6])
    state_table=[]
    for L in LANDMARKS:
        for r in state_recurrence:
            if r["landmark"]==f"L{L}":
                state_table.append(f"| L{L} | {r['digital_state_zh']} | {r['n']} | {r['events']} | {r['event_rate']:.1%} | ({L},12] |")
    flow_table=[f"| {STATE_ZH[r['state_L3']]} | {STATE_ZH[r['state_L6']]} | {r['n']} | {r['events_6_to_12m']} | {r['recurrence_rate_6_to_12m']:.1%} | {r['median_risk_L3']:.3f} | {r['median_risk_L6']:.3f} | {r['median_delta_risk_L6_minus_L3']:+.3f} |" for r in transition_matrix]
    status_table=[]
    test_lookup={(r["level"],r["group"]):r for r in paired_stats}
    for r in update_rows:
        if r["level"]=="current_status_transition":
            lo = "" if r["delta_ci_lower"] == "" else f"{r['delta_ci_lower']:+.3f}"
            hi = "" if r["delta_ci_upper"] == "" else f"{r['delta_ci_upper']:+.3f}"
            test = test_lookup[(r["level"], r["group"])]
            p = test["wilcoxon_two_sided_p"]
            if p == "": ptxt = "未检验"
            elif isinstance(p, str): ptxt = p
            elif p < .001: ptxt = "<0.001"
            else: ptxt = f"{p:.3g}"
            status_table.append(f"| {STATUS_ZH[r['group']]} | {r['n']} | {r['events_6_to_12m']} | {r['median_risk_L3']:.3f} | {r['median_risk_L6']:.3f} | {r['median_delta_risk']:+.3f} | {lo}–{hi} | {ptxt} |")
    memory_table=[]
    for r in memory_rows:
        if r["level"].startswith("contrast"):
            L=r["landmark"]
            if r.get("median_risk_difference","")!="":
                analysis_zh = "当前阴性：既往阳性 − 从未阳性" if r["analysis"] == "all_current_negative_group" else "阳性后清除 − 持续阴性"
                md_pp = 100 * float(r["median_risk_difference"])
                lo = r.get("difference_ci_lower", ""); hi = r.get("difference_ci_upper", "")
                ci_txt = "未估计" if lo == "" or hi == "" else f"[{100*float(lo):+.1f}, {100*float(hi):+.1f}] 个百分点"
                rate_txt = "未估计" if r.get("event_rate_difference", "") == "" else f"{100*float(r['event_rate_difference']):+.1f} 个百分点"
                memory_table.append(f"| L{L} | {analysis_zh} | {r['n']} | {md_pp:+.1f} 个百分点 | {ci_txt} | {rate_txt} |")
    report = f"""# CRC MRD 状态转移型预测数字孪生分析

**项目：** `E:/crc_mrd_digital_twin`  
**分析范围：** Chen、ColonAiQ、COSMOS；Tie 不纳入。  
**分析对象：** 已冻结的三队列 L3/L6 leave-one-cohort-out cross-fitted 预测；本分析没有重拟合模型、改变超参数或更新既有 AUC。

## 分析设计

- L3 预测术后 (3,12] 月复发；L6 预测 (6,12] 月复发。仅使用该 landmark 及以前的 MRD 记录，时间泄漏检查全部通过。
- 同一患者 L3→L6 配对样本为 **{len(pairs)}** 人（L3/L6 预测患者交集）；L3 风险集 **{total3}** 人、{ev3} 个 (3,12] 月复发结局；L6 风险集 **{total6}** 人、{ev6} 个 (6,12] 月复发结局。部分患者首次 MRD 采样晚于术后3个月，虽进入 L6 但没有 L3 预测，故不进入配对更新。
- MRD 状态由原始队列解析脚本重建：单次阴/阳、持续阴/阳、阳性后清除、阴转阳、间歇/混合。当前状态和既往阳性另行编码。
- 风险更新定义为 Risk6 − Risk3。两者预测时间窗不同，变化同时包含新的 MRD 信息、时间推进及进入 L6 风险集的条件变化，不能视为单纯 MRD 的因果效应。

## L3 与 L6 的状态人数和随访期复发

事件率按各 landmark 的后续窗口计算；分层人数较少时仅作描述。

| Landmark | MRD 状态 | n | 后续复发 | 事件率 | 结局窗 |
|---|---|---:|---:|---:|---|
{chr(10).join(state_table)}

## L3→L6 数字状态转移

| L3 状态 | L6 状态 | n | (6,12]复发 | 事件率 | 中位 Risk3 | 中位 Risk6 | 中位 Δ风险 |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(flow_table)}

## 按当前 MRD 阴阳转移的配对风险更新

患者 bootstrap CI 和 signed-rank 检验是探索性摘要；小样本组不做推断。未校正多重比较。
bootstrap 区间以患者为单位重抽样，但条件于已冻结的预测值，不包含模型重拟合不确定性；因低维模型产生重复预测值，部分区间可能很窄或退化。

| 当前状态转移 | n | (6,12]复发 | 中位 Risk3 | 中位 Risk6 | 中位 Δ风险 | bootstrap 95% CI | 配对 Wilcoxon p |
|---|---:|---:|---:|---:|---:|---|---:|
{chr(10).join(status_table)}

## 当前 MRD 阴性中的分子记忆

比较当前阴性且既往曾阳性与截至当前从未阳性者的冻结模型预测风险；同时单列持续阴性与阳性后清除。该分析量化预测模型在给定历史下的分层，不证明 MRD 历史对复发的因果作用。区间按患者分组 bootstrap；不足5人的组不估计区间。

| Landmark | 分析 | 总 n | 预测风险中位差 | bootstrap 95% CI | 后续事件率差 |
|---|---|---:|---:|---:|---|
{chr(10).join(memory_table)}

## 图形和患者级文件

- `outputs/figures/digital_state_sankey.svg`：L3→L6 状态人数流向。
- `outputs/figures/risk_delta_by_transition.svg`：按当前 MRD 转移分层的风险更新。
- `outputs/figures/molecular_memory_effect.svg`：当前阴性人群的既往阳性记忆比较。
- `outputs/local_only/risk_update_spaghetti.svg`：匿名患者级配对风险轨迹，仅保存在本地。
- `outputs/local_only/digital_twin_case_examples.svg`：公开队列真实个案示意，仅保存在本地。
- `outputs/local_only/patient_digital_state_trajectory.csv` 和 `patient_state_transition_records.csv`：患者级数据，仅保存在本地，不上传 GitHub。

## 解释边界

本分析使用既有 cross-fitted 风险预测，描述动态状态、预测值变化及预测分层。配对的 Risk3 与 Risk6 并非相同终点窗口，因此其差异不是纯粹由 MRD 更新产生。分子记忆是预测关联，不是因果效应。状态组内事件率未经竞争风险调整，也不能直接作临床治疗决策。该分析不模拟治疗干预或反事实治疗获益；临床效用仍需独立前瞻性验证。
"""
    (REPORT_DIR / "mrd_state_digital_twin_analysis_zh.md").write_text(report, encoding="utf-8")

    # Machine-readable aggregate total audit without any patient identifiers.
    audit = [{"landmark": f"L{L}", "n": len(rows_by_lm[L]), "events": sum(r["outcome"] for r in rows_by_lm[L]),
              "cohort_n": ";".join(f"{c}:{sum(r['cohort']==c for r in rows_by_lm[L])}" for c in COHORTS),
              "prediction_source": str(PRED_PATH.name), "model_refit": "no"} for L in LANDMARKS]
    audit.append({"landmark":"L3_L6_paired", "n":len(pairs), "events":sum(r["event6"] for r in pairs),
                  "cohort_n":";".join(f"{c}:{sum(r['cohort']==c for r in pairs)}" for c in COHORTS),
                  "prediction_source":str(PRED_PATH.name),"model_refit":"no"})
    write_csv(TABLE_DIR / "digital_state_analysis_sample_audit.csv", audit)
    print(f"DONE | L3={total3}/{ev3} | L6={total6}/{ev6} | paired={len(pairs)}")
    print("state counts L3:", {r["digital_state"]:r["n"] for r in state_recurrence if r["landmark"]=="L3"})
    print("state counts L6:", {r["digital_state"]:r["n"] for r in state_recurrence if r["landmark"]=="L6"})
    print("status transitions:", dict(Counter(r["status_transition"] for r in pairs)))
    print("outputs: aggregate tables/report/3 aggregate SVG; all patient-level outputs remain in outputs/local_only/")


if __name__ == "__main__":
    main()
