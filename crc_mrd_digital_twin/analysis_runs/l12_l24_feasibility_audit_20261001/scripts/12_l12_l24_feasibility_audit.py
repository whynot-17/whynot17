from __future__ import annotations

"""Availability-only audit for a postoperative 12-to-24-month MRD stage.

No predictive model is fitted. Aggregate rows are public-safe; the patient-level
eligibility ledger stays under outputs/local_only/ on the E: project.
"""

import csv
import importlib.util
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pdfplumber

PROJECT = Path(r"E:\crc_mrd_digital_twin")
OUT = PROJECT / "analysis_runs" / "l12_l24_feasibility_audit_20261001"
PUBLIC_OUT = Path(r"C:\Users\ASUS\Documents\Codex\2026-09-29\files-pasted-by-the-user-mrd\outputs\crc_mrd_l12_l24_feasibility_20261001")
SCRIPT_TARGET = PROJECT / "scripts" / "12_l12_l24_feasibility_audit.py"
MONTH_DAYS = 365.2425 / 12
DAY12, DAY24 = 12 * MONTH_DAYS, 24 * MONTH_DAYS


def load_analysis():
    path = PROJECT / "scripts" / "10_galaxy_l3_l6_l12_cohort_nodes.py"
    spec = importlib.util.spec_from_file_location("l12_l24_audit_analysis10", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_patients():
    mod = load_analysis()
    cohorts = mod.source_patients(PROJECT)
    return mod, cohorts


def choose_l12(cohort, patient):
    history = patient.get("history_days", [])
    if cohort == "COSMOS":
        candidates = [(t, s) for t, s in history if abs(t / MONTH_DAYS - 12.0) < 1e-6]
    elif cohort == "Chen":
        candidates = [(t, s) for t, s in history if 10.5 <= t / MONTH_DAYS <= 13.5]
    else:
        return None, None
    if not candidates:
        return None, None
    t, s = min(candidates, key=lambda v: (abs(v[0] / MONTH_DAYS - 12.0), v[0]))
    return float(t), int(s)


def classify_outcome(patient, origin):
    event_time = patient.get("event_time_days")
    followup = patient.get("followup_days")
    if event_time is not None:
        if event_time <= origin:
            return None, "recurrence_at_or_before_L12"
        if event_time <= DAY24:
            return 1, "recurrence_after_L12_by_month24"
        return 0, "recurrence_after_month24"
    if followup is not None and followup >= DAY24:
        return 0, "confirmed_recurrence_free_through_month24"
    return None, "followup_censored_before_month24"


def audit_galaxy():
    path = PROJECT / "data_raw" / "galaxy_2024_supplementary_tables_1_3.pdf"
    rows = []
    with pdfplumber.open(path) as pdf:
        for page in pdf.pages[2:45]:
            for table in page.extract_tables() or []:
                for raw in table:
                    if (len(raw) < 19 or not re.fullmatch(r"\s*\d+\s*", str(raw[0] or ""))
                            or str(raw[1] or "").strip() not in {"Male", "Female"}):
                        continue
                    r = [str(x or "").strip() for x in raw[:19]]
                    pid = int(r[0])
                    recurred = r[4].lower() == "yes"
                    dfs_days = float(r[6]) if r[6] else None
                    surveillance = r[14].upper() if r[14] else "MISSING"
                    pos_date = float(r[17]) if r[17] and re.fullmatch(r"\d+(?:\.\d+)?", r[17]) else None
                    rows.append({"pid": pid, "recurred": recurred, "dfs_days": dfs_days,
                                 "surveillance_window": surveillance, "post_mrd_positive": r[16].upper(),
                                 "post_mrd_positive_day": pos_date})
    if len(rows) != 2240:
        raise RuntimeError(f"Expected 2240 GALAXY rows; found {len(rows)}")
    status_counts = Counter(r["surveillance_window"] for r in rows)
    positive_dates = [r["post_mrd_positive_day"] for r in rows if r["post_mrd_positive_day"] is not None]
    date_bins = Counter()
    for day in positive_dates:
        if day <= DAY12:
            date_bins["<=12m"] += 1
        elif day <= DAY24:
            date_bins["12-24m"] += 1
        else:
            date_bins[">24m"] += 1
    outcome = Counter()
    for r in rows:
        t = r["dfs_days"]
        if t is None:
            outcome["missing_DFS"] += 1
        elif r["recurred"]:
            if t <= DAY12:
                outcome["recurrence_by_12m"] += 1
            elif t <= DAY24:
                outcome["recurrence_12_to_24m"] += 1
            else:
                outcome["recurrence_after_24m"] += 1
        elif t >= DAY24:
            outcome["followup_to_24m_or_longer_no_recurrence"] += 1
        else:
            outcome["nonrecurrence_followup_under_24m"] += 1
    return rows, status_counts, date_bins, outcome


def audit_cohort(cohort, patients):
    counts = Counter()
    ledger = []
    groups = defaultdict(Counter)
    counts["total_patients"] = len(patients)
    for p in patients:
        counts["patients_with_any_ctdna_history"] += bool(p.get("history_days"))
        l12_day, l12_status = choose_l12(cohort, p)
        if l12_status is None:
            counts["no_classifiable_L12_status"] += 1
            continue
        counts["classifiable_L12_status"] += 1
        counts["L12_negative"] += int(l12_status == 0)
        counts["L12_positive"] += int(l12_status == 1)
        et = p.get("event_time_days")
        if et is not None and et <= l12_day:
            counts["recurrence_at_or_before_L12"] += 1
            outcome, reason = None, "recurrence_at_or_before_L12"
        else:
            outcome, reason = classify_outcome(p, l12_day)
        if reason == "followup_censored_before_month24":
            counts["followup_censored_before_month24"] += 1
        elif reason == "recurrence_at_or_before_L12":
            pass
        elif outcome == 1:
            counts["eligible_recurrence_12_to_24m"] += 1
        elif outcome == 0:
            counts["eligible_control_known_through_24m"] += 1
        if outcome in (0, 1):
            counts["eligible_known_L12_to_24m"] += 1
            groups["negative" if l12_status == 0 else "positive"]["N"] += 1
            groups["negative" if l12_status == 0 else "positive"]["events"] += int(outcome == 1)
        ledger.append({
            "cohort": cohort, "patient_id": str(p["patient_id"]),
            "l12_day": round(l12_day, 3), "l12_month": round(l12_day / MONTH_DAYS, 3),
            "l12_mrd": "positive" if l12_status else "negative",
            "event_time_months": round(et / MONTH_DAYS, 3) if et is not None else "",
            "followup_months": round(p["followup_days"] / MONTH_DAYS, 3) if p.get("followup_days") is not None else "",
            "outcome_12_to_24": outcome if outcome is not None else "unknown",
            "outcome_reason": reason,
        })
    return counts, groups, ledger


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    _, cohorts = load_patients()
    summary_rows, transition_rows, private_rows = [], [], []
    all_counts, all_groups = Counter(), defaultdict(Counter)
    for cohort in ("GALAXY", "COSMOS", "Chen"):
        if cohort == "GALAXY":
            patients, surveillance, date_bins, outcomes = audit_galaxy()
            status_available = 0  # No dedicated patient-level 12-month ctDNA field in Table S2.
            row = {
                "cohort": cohort, "total_patients": len(patients),
                "patients_with_classifiable_L12_MRD": status_available,
                "L12_ctDNA_negative": "not separately reported",
                "L12_ctDNA_positive": "not separately reported",
                "recurrence_by_12m": outcomes["recurrence_by_12m"],
                "recurrence_12_to_24m": outcomes["recurrence_12_to_24m"],
                "confirmed_no_recurrence_to_24m": outcomes["followup_to_24m_or_longer_no_recurrence"],
                "nonrecurrence_followup_under_24m": outcomes["nonrecurrence_followup_under_24m"],
                "known_L12_to_24_risk_set_with_MRD": 0,
                "eligible_events_12_to_24_with_MRD": 0,
                "note": "Table S2 reports broad surveillance-window status and positive-post-MRD date, not a per-patient 12-month status/date; cannot define L12 current state.",
            }
            extra = [
                {"cohort": cohort, "item": f"ctDNA_surveillance_window_{k}", "n": v} for k, v in sorted(surveillance.items())
            ] + [
                {"cohort": cohort, "item": f"first_post_MRD_positive_date_{k}", "n": v} for k, v in sorted(date_bins.items())
            ] + [
                {"cohort": cohort, "item": f"clinical_outcome_{k}", "n": v} for k, v in sorted(outcomes.items())
            ]
            transition_rows.extend(extra)
        else:
            counts, groups, ledger = audit_cohort(cohort, cohorts[cohort])
            private_rows.extend(ledger)
            row = {"cohort": cohort, **counts,
                   "known_L12_to_24_risk_set_with_MRD": counts["eligible_known_L12_to_24m"],
                   "eligible_events_12_to_24_with_MRD": counts["eligible_recurrence_12_to_24m"],
                   "note": "COSMOS uses its scheduled month-12 node; Chen uses closest observed serial ctDNA within 10.5-13.5 months, after requiring recurrence-free status at sample time."}
            for state in ("negative", "positive"):
                g = groups[state]
                transition_rows.append({"cohort": cohort, "item": f"L12_{state}_known_12_to_24_outcome_N", "n": g["N"]})
                transition_rows.append({"cohort": cohort, "item": f"L12_{state}_events_12_to_24", "n": g["events"]})
        summary_rows.append(row)

    out_tables = OUT / "outputs" / "tables"
    public_tables = PUBLIC_OUT / "tables"
    fields = sorted({k for r in summary_rows for k in r})
    write_csv(out_tables / "cohort_l12_l24_feasibility.csv", summary_rows, fields)
    write_csv(out_tables / "l12_mrd_state_event_counts.csv", transition_rows, ["cohort", "item", "n"])
    write_csv(OUT / "outputs" / "local_only" / "patient_l12_l24_eligibility_ledger.csv", private_rows,
              ["cohort", "patient_id", "l12_day", "l12_month", "l12_mrd", "event_time_months", "followup_months", "outcome_12_to_24", "outcome_reason"])

    refs = """# L12→L24 动态阶段可行性审计（2026-10-01）

## 目的

评估现有公开患者级数据能否支撑独立的 12 个月 MRD 状态到 12–24 个月临床复发阶段。该审计不拟合新模型、不改变已冻结的 L3→L6→L12 主分析。

结局规则为：患者在所选 L12 MRD 评估时尚未临床复发；事件须发生在该评估之后且不晚于术后 24 个月。无事件者需有资料确认观察/随访至少到 24 个月；更早删失不作为对照。COSMOS 使用原表 scheduled month-12 状态；Chen 使用其季度采样记录中 10.5–13.5 个月内最接近 12 个月、且先于临床复发的 ctDNA。Chen 的实际采血时间会作为 landmark；因此该队列的结局窗是该次实际评估至术后 24 个月。

## 可用性结果

| 队列 | 总患者 | 有可判定 L12 状态记录 | L12→L24 已知结局风险集 | 12–24 月复发事件 | 24 月无复发对照 | 判断 |
|---|---:|---:|---:|---:|---:|---|
""" + "\n".join(
        f"| {r['cohort']} | {r.get('total_patients', 0)} | {r.get('classifiable_L12_status', r.get('patients_with_classifiable_L12_MRD', 0))} | {r.get('eligible_known_L12_to_24m', r.get('known_L12_to_24_risk_set_with_MRD', 0))} | {r.get('eligible_recurrence_12_to_24m', r.get('eligible_events_12_to_24_with_MRD', 0))} | {r.get('eligible_control_known_through_24m', r.get('confirmed_no_recurrence_to_24m', 0))} | {('不可用于 L12 MRD 更新：缺少患者级 12 月节点' if r['cohort']=='GALAXY' else '仅适合作为探索/支持性阶段：风险集或事件数有限')} |"
        for r in summary_rows
    ) + """

“有可判定 L12 状态记录”包括在发生临床复发后才采到的 L12 附近样本；这些患者不进入 landmark 风险集。COSMOS 有 6 例、Chen 有 5 例复发发生在其 L12 采样之前。

GALAXY 全队列可重建临床结局，但不能与 L12 MRD 状态配对：术后 12 个月前复发 373 例，12–24 个月复发 125 例，24 个月后复发 16 例；704 例可确认无复发随访达到至少 24 个月，另有 1,022 例无复发记录但随访不到 24 个月，不能作为 24 月对照。

## L12 MRD 状态与 12–24 月结局

| 队列 | L12 阴性组 | L12 阳性组 | 风险集总计 |
|---|---:|---:|---:|
| COSMOS | 16/247 复发（6.5%） | 3/9（33.3%） | 256 人，19 事件 |
| Chen | 6/73（8.2%） | 6/6（100%） | 79 人，12 事件 |
| COSMOS + Chen（描述性合并） | 22/320（6.9%） | 9/15（60.0%） | 335 人，31 事件 |

合并率只作描述性展示；两队列的检测平台、采样密度和入组人群不同。Chen 的阳性组仅 6 人，COSMOS 阳性组仅 9 人；总事件 31 例不足以支撑复杂多变量模型或稳定的留一队列外部验证。若后续进入模型阶段，应把该阶段标为预先限定的探索/支持性分析，并限制自由度。

## 队列判断

### GALAXY

补充表 2 有独立的 MRD-window、3 月、6 月列，但其后只报告宽泛的 “ctDNA Surveillance Window” 状态和 post-MRD 阳性日期。它没有逐患者的 12 月阴/阳结果及相应采样日，因此不能将 surveillance-window 阴/阳可靠地等同于 L12 当前状态。原文报告的采样计划包括术后 4、12、24、36、48、72、96 周并持续至复发；48 周附近存在计划采样节点，但公开补充表未提供所有患者该节点的逐项结果/日期。日期缺失时，阳性转化日期只能部分标识阳性发生时间，不能确定无转化患者在 L12 的当前状态。

因此 GALAXY 可用于描述临床 12–24 月复发结局数，但**不能进入 L12 MRD 动态更新模型**。不得把 surveillance-window 状态重命名为 L12 MRD。

### COSMOS 与 Chen

COSMOS 补充表含 scheduled 12-、18-、24-month ctDNA 状态及 RFS 时间；Chen 的患者级 serial ctDNA 与 recurrence/follow-up 时间可用于按上述风险集规则核对。它们可以提供一个独立于 GALAXY 的 L12→L24 探索/验证阶段；是否适合作为正式主分析，需要根据风险集事件数、状态分层事件数和删失比例决定。该结论不回写或替换已冻结的主分析队列角色。

## 数据治理与文件

- `outputs/tables/cohort_l12_l24_feasibility.csv`：队列级可用性及结局计数。
- `outputs/tables/l12_mrd_state_event_counts.csv`：L12 MRD 状态分层的 N/事件数与 GALAXY surveillance 审计计数。
- `outputs/local_only/patient_l12_l24_eligibility_ledger.csv`：患者级节点与结局审计，只留在 E 盘本地，不提交 GitHub。
- `scripts/12_l12_l24_feasibility_audit.py`：生成此项可用性审计，不拟合模型。

公开来源：GALAXY [Nature Medicine 2024](https://doi.org/10.1038/s41591-024-03254-6)；COSMOS [公开研究与补充表](https://pmc.ncbi.nlm.nih.gov/articles/PMC11443202/)；Chen [公开研究与补充表](https://pmc.ncbi.nlm.nih.gov/articles/PMC8130394/)。
"""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "reports").mkdir(parents=True, exist_ok=True)
    (OUT / "reports" / "l12_l24_feasibility_zh.md").write_text(refs, encoding="utf-8")
    OUT_SCRIPTS = OUT / "scripts"
    OUT_SCRIPTS.mkdir(parents=True, exist_ok=True)
    SCRIPT_TARGET.parent.mkdir(parents=True, exist_ok=True)
    SCRIPT_TARGET.write_bytes(Path(__file__).read_bytes())
    (OUT_SCRIPTS / SCRIPT_TARGET.name).write_bytes(Path(__file__).read_bytes())
    PUBLIC_OUT.mkdir(parents=True, exist_ok=True)
    (PUBLIC_OUT / "reports").mkdir(parents=True, exist_ok=True)
    (PUBLIC_OUT / "tables").mkdir(parents=True, exist_ok=True)
    (PUBLIC_OUT / "reports" / "l12_l24_feasibility_zh.md").write_text(refs, encoding="utf-8")
    for p in (out_tables / "cohort_l12_l24_feasibility.csv", out_tables / "l12_mrd_state_event_counts.csv"):
        (PUBLIC_OUT / "tables" / p.name).write_bytes(p.read_bytes())
    print(json.dumps({"public_aggregate": summary_rows, "state_counts": transition_rows,
                      "galaxy_surveillance_status_counts": dict(surveillance),
                      "galaxy_post_mrd_positive_date_bins": dict(date_bins),
                      "galaxy_clinical_outcomes": dict(outcomes),
                      "project_output": str(OUT), "public_output": str(PUBLIC_OUT)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

