from __future__ import annotations

"""Patient-level feasibility audit of GALAXY 2024 Supplementary Table 2.

Only aggregate output is written. Raw PDF remains in data_raw; no row-level
data or patient identifiers are exported or committed.
"""

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import pdfplumber

ROOT = Path(r"E:\crc_mrd_digital_twin")
PDF = ROOT / "data_raw" / "galaxy_2024_supplementary_tables_1_3.pdf"
TABLES = ROOT / "outputs" / "tables"
REPORTS = ROOT / "reports"
DAY3 = 3 * 365.2425 / 12
DAY6 = 6 * 365.2425 / 12
DAY12 = 365.2425

TABLES.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)


def clean(x):
    return str(x or "").strip()


def status(x):
    s = clean(x).upper()
    return s if s in ("POSITIVE", "NEGATIVE") else None


def yesno(x):
    s = clean(x).lower()
    return True if s == "yes" else False if s == "no" else None


def number(x):
    try:
        return float(clean(x))
    except (TypeError, ValueError):
        return None


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_records():
    records = []
    with pdfplumber.open(PDF) as pdf:
        if len(pdf.pages) != 46:
            raise RuntimeError(f"Unexpected supplement page count: {len(pdf.pages)}")
        # Table 2 starts on printed page 2 (zero-based page index 2) and ends
        # before Supplementary Table 3 on printed page 45.
        for page in pdf.pages[2:45]:
            tables = page.extract_tables()
            if not tables:
                continue
            for row in tables[0]:
                if (len(row) >= 19 and re.fullmatch(r"\s*\d+\s*", clean(row[0]))
                        and clean(row[1]) in ("Male", "Female")):
                    r = [clean(x) for x in row[:19]]
                    records.append({
                        "row_number": int(r[0]),
                        "sex": r[1],
                        "stage": r[3],
                        "recurred": yesno(r[4]),
                        "dead": yesno(r[5]),
                        "dfs_days": number(r[6]),
                        "mrd0": status(r[10]),
                        "mrd0_day": number(r[11]),
                        "mrd3": status(r[12]),
                        "mrd6": status(r[13]),
                        "surveillance": status(r[14]),
                        "clearance": clean(r[15]),
                        "molecular_recurrence": yesno(r[16]),
                        "molecular_recurrence_day": number(r[17]),
                    })
    ids = [r["row_number"] for r in records]
    if len(records) != 2240 or len(set(ids)) != 2240 or sorted(ids) != list(range(1, 2241)):
        raise RuntimeError("GALAXY Supplementary Table 2 did not parse as exactly 2,240 unique sequential rows")
    return records


def y12(r):
    """Conservative recurrence-only 12m status, distinct from DFS composite."""
    t, rec = r["dfs_days"], r["recurred"]
    if t is None or rec is None:
        return None
    if rec is True:
        return int(t <= DAY12)
    if rec is False and t >= DAY12:
        return 0
    # No recurrence flag with DFS follow-up <12m is censored/competing-event unknown.
    return None


def eligible_at(r, landmark_day):
    y = y12(r)
    return y is not None and r["dfs_days"] is not None and r["dfs_days"] >= landmark_day


def pattern3(r):
    seq = [r[k] for k in ("mrd0", "mrd3", "mrd6")]
    if any(x is None for x in seq):
        return None
    if all(x == "NEGATIVE" for x in seq):
        return "persistent_negative"
    if all(x == "POSITIVE" for x in seq):
        return "persistent_positive"
    if seq[-1] == "NEGATIVE" and "POSITIVE" in seq[:-1]:
        return "positive_to_negative_clearance"
    if seq[-1] == "POSITIVE" and all(x == "NEGATIVE" for x in seq[:-1]):
        return "negative_to_positive_conversion"
    return "intermittent_mixed"


def main():
    if not PDF.exists():
        raise FileNotFoundError(PDF)
    records = read_records()
    totals = {
        "total_n": len(records),
        "known_recurrence_flag_and_dfs_days": sum(r["recurred"] is not None and r["dfs_days"] is not None for r in records),
        "recurrence_yes_any_followup": sum(r["recurred"] is True for r in records),
        "y12_known": sum(y12(r) is not None for r in records),
        "y12_recurrence_events": sum(y12(r) == 1 for r in records),
        "y12_confirmed_no_recurrence": sum(y12(r) == 0 for r in records),
        "y12_unknown_early_censor_or_competing_event": sum(y12(r) is None for r in records),
        "mrd_window_status_available": sum(r["mrd0"] is not None for r in records),
        "ctdna_3m_status_available": sum(r["mrd3"] is not None for r in records),
        "ctdna_6m_status_available": sum(r["mrd6"] is not None for r in records),
        "ctdna_3m_and_6m_both_available": sum(r["mrd3"] is not None and r["mrd6"] is not None for r in records),
        "mrd_window_3m_6m_complete": sum(all(r[k] is not None for k in ("mrd0", "mrd3", "mrd6")) for r in records),
    }
    # Study-defined timepoint status summaries (not exact specimen dates).
    status_rows = []
    for label, field, window in [
        ("MRD window", "mrd0", "2–10 weeks post-surgery; date field supplied"),
        ("3-month ctDNA", "mrd3", "study-defined 70–112-day window; individual draw date absent"),
        ("6-month ctDNA", "mrd6", "study-defined 160–200-day window; individual draw date absent"),
    ]:
        vals = [r[field] for r in records]
        status_rows.append({"timepoint": label, "nominal_window": window, "n_total": len(records),
                            "n_available": sum(v is not None for v in vals),
                            "n_positive": sum(v == "POSITIVE" for v in vals),
                            "n_negative": sum(v == "NEGATIVE" for v in vals),
                            "n_missing_or_unclassifiable": sum(v is None for v in vals),
                            "patient_specific_collection_date_available": "yes" if field == "mrd0" else "no"})
    write_csv(TABLES / "galaxy_timepoint_availability.csv", status_rows)

    both = [r for r in records if r["mrd3"] is not None and r["mrd6"] is not None]
    complete = [r for r in both if r["mrd0"] is not None]
    availability_rows = []
    def add(metric, n, denom=2240, events="", controls="", note=""):
        availability_rows.append({"metric": metric, "n": n, "denominator": denom,
                                  "percent_of_denominator": (n / denom if denom else ""),
                                  "known_12m_recurrence_events": events,
                                  "known_12m_nonrecurrence": controls,
                                  "definition_note": note})
    add("total_patients", len(records), note="Parsed patient rows; sequential row numbers were validated locally.")
    add("complete_recurrence_flag_and_DFS_days", totals["known_recurrence_flag_and_dfs_days"], note="Recurred yes/no and DFS days both nonmissing.")
    add("recurrence_yes_any_followup", totals["recurrence_yes_any_followup"], note="Recurred=Yes at any reported follow-up.")
    add("confirmed_12m_recurrence_endpoint", totals["y12_known"], events=totals["y12_recurrence_events"], controls=totals["y12_confirmed_no_recurrence"],
        note="Recurrence by 365.24d=event; recurrence after 365.24d or no recurrence with DFS≥365.24d=control; earlier no-recurrence DFS is unknown.")
    add("3m_ctDNA_status_available", totals["ctdna_3m_status_available"], note="Study-defined 3-month status; exact day unavailable.")
    add("6m_ctDNA_status_available", totals["ctdna_6m_status_available"], note="Study-defined 6-month status; exact day unavailable.")
    add("both_3m_and_6m_ctDNA_status_available", len(both), note="Nominal 3m and 6m status columns both classifiable.")
    add("MRD_window_plus_3m_plus_6m_complete", len(complete), note="All three statuses classifiable; enables baseline-history states.")

    for L, f, cutoff in [(3, "mrd3", DAY3), (6, "mrd6", DAY6)]:
        g = [r for r in records if r[f] is not None]
        risk = [r for r in g if eligible_at(r, cutoff)]
        events = sum(y12(r) == 1 for r in risk)
        controls = len(risk) - events
        add(f"L{L}_nominal_status_Y12_known_and_event_free_at_{cutoff:.1f}d", len(risk), events=events, controls=controls,
            note=f"Study-labelled {L}m status; requires known recurrence-only Y12 and DFS≥{cutoff:.1f}d. Exact draw date unavailable.")
    for label, g, cutoff in [("paired_3m6m", both, DAY6), ("complete_MRDwindow3m6m", complete, DAY6)]:
        risk = [r for r in g if eligible_at(r, cutoff)]
        events = sum(y12(r) == 1 for r in risk)
        controls = len(risk) - events
        add(f"{label}_L6_nominal_risk_set", len(risk), events=events, controls=controls,
            note=f"Both nominal timepoint statuses; known Y12 and DFS≥{cutoff:.1f}d. Actual collection dates not supplied.")
    write_csv(TABLES / "galaxy_patient_level_availability_audit.csv", availability_rows)

    transition_rows = []
    for a in ("NEGATIVE", "POSITIVE"):
        for b in ("NEGATIVE", "POSITIVE"):
            g = [r for r in both if r["mrd3"] == a and r["mrd6"] == b]
            risk = [r for r in g if eligible_at(r, DAY6)]
            ev = sum(y12(r) == 1 for r in risk)
            transition_rows.append({"transition_3m_to_6m": f"{a.lower()}_to_{b.lower()}",
                                    "transition_zh": {("NEGATIVE","NEGATIVE"):"阴→阴",("NEGATIVE","POSITIVE"):"阴→阳",("POSITIVE","NEGATIVE"):"阳→阴",("POSITIVE","POSITIVE"):"阳→阳"}[(a,b)],
                                    "n_with_both_statuses": len(g), "n_L6_nominal_risk_set_Y12_known": len(risk),
                                    "recurrence_events_6_to_12m": ev,
                                    "confirmed_no_recurrence_through_12m": len(risk)-ev,
                                    "note":"L6 risk set uses surgery-relative 182.6-day threshold; individual 6m specimen dates are absent."})
    pattern_counts = Counter(pattern3(r) for r in complete)
    for pattern, n in sorted(pattern_counts.items()):
        g = [r for r in complete if pattern3(r) == pattern]
        risk = [r for r in g if eligible_at(r, DAY6)]
        ev = sum(y12(r) == 1 for r in risk)
        transition_rows.append({"transition_3m_to_6m": "three_timepoint_state:" + pattern,
                                "transition_zh": pattern, "n_with_both_statuses": n,
                                "n_L6_nominal_risk_set_Y12_known": len(risk),
                                "recurrence_events_6_to_12m": ev,
                                "confirmed_no_recurrence_through_12m": len(risk)-ev,
                                "note":"State uses MRD-window, 3m and 6m statuses; study-defined timepoints only."})
    write_csv(TABLES / "galaxy_nominal_mrd_state_transitions.csv", transition_rows)

    l3 = [r for r in records if r["mrd3"] is not None and eligible_at(r, DAY3)]
    l6 = [r for r in records if r["mrd6"] is not None and eligible_at(r, DAY6)]
    paired_l3 = [r for r in both if eligible_at(r, DAY3)]
    paired_l6 = [r for r in both if eligible_at(r, DAY6)]
    complete_l6 = [r for r in complete if eligible_at(r, DAY6)]
    report = f"""# GALAXY 2024 患者级 MRD 数据可用性审计

**项目：** `E:/crc_mrd_digital_twin`  
**审计日期：** 2026-09-30  
**数据：** Nature Medicine 2024 Supplementary Tables 1–3 PDF；仅使用其中 Supplementary Table 2，原始 PDF 留在本地 `data_raw/`，不提交 GitHub。

## 六项核心数字

| 指标 | 人数 | 口径 |
|---|---:|---|
| 总患者数 | {len(records)} | 表2中连续编号且患者行解析成功 |
| 有 Recurred 标志和 DFS 天数 | {totals['known_recurrence_flag_and_dfs_days']} | 两字段均非缺失；其中任意随访复发 {totals['recurrence_yes_any_followup']} 例 |
| 有 study-defined 3-month ctDNA 状态 | {totals['ctdna_3m_status_available']} | 阳性 {sum(r['mrd3']=='POSITIVE' for r in records)}；阴性 {sum(r['mrd3']=='NEGATIVE' for r in records)} |
| 有 study-defined 6-month ctDNA 状态 | {totals['ctdna_6m_status_available']} | 阳性 {sum(r['mrd6']=='POSITIVE' for r in records)}；阴性 {sum(r['mrd6']=='NEGATIVE' for r in records)} |
| 同时有 3m 与 6m 状态 | {len(both)} | 两个名义时间点均有阳性/阴性结果 |
| 其中 12m 复发终点可判定 | {sum(y12(r) is not None for r in both)} | 复发事件 {sum(y12(r)==1 for r in both)}；明确非复发 {sum(y12(r)==0 for r in both)} |

另外，**MRD window + 3m + 6m 三个状态均齐全 {len(complete)} 人**，可构造基线 MRD 历史状态。

## 按当前项目终点定义估计可用风险集

采用术后天数作为近似边界：3个月 = {DAY3:.1f} 天，6个月 = {DAY6:.1f} 天，12个月 = {DAY12:.1f} 天。12个月复发终点按复发日期判定：复发≤365.24天为事件；复发>365.24天或无复发且 DFS≥365.24天为明确非复发；无复发但 DFS<365.24天为未知。死亡在复发终点中作为竞争事件/未知，不编码为无复发。

| 分析集 | n | 12个月内复发 | 明确无复发 | 注 |
|---|---:|---:|---:|---|
| L3：有3m状态、Y12可判定且DFS≥{DAY3:.1f}天 | {len(l3)} | {sum(y12(r)==1 for r in l3)} | {len(l3)-sum(y12(r)==1 for r in l3)} | 名义3m landmark |
| L6：有6m状态、Y12可判定且DFS≥{DAY6:.1f}天 | {len(l6)} | {sum(y12(r)==1 for r in l6)} | {len(l6)-sum(y12(r)==1 for r in l6)} | 名义6m landmark |
| 3m→6m 配对轨迹，L3 风险集 | {len(paired_l3)} | {sum(y12(r)==1 for r in paired_l3)} | {len(paired_l3)-sum(y12(r)==1 for r in paired_l3)} | 两状态都齐全 |
| 3m→6m 配对轨迹，L6 风险集 | {len(paired_l6)} | {sum(y12(r)==1 for r in paired_l6)} | {len(paired_l6)-sum(y12(r)==1 for r in paired_l6)} | 两状态都齐全 |
| MRD window+3m+6m 完整，L6 风险集 | {len(complete_l6)} | {sum(y12(r)==1 for r in complete_l6)} | {len(complete_l6)-sum(y12(r)==1 for r in complete_l6)} | 具备完整既往状态 |

现有三队列模型 L3 为 658 人/60 事件，L6 为 642 人/40 事件，L3→L6 配对为 638 人。若采用 GALAXY 名义时间点，L3 端可与既有集拼接至约 **{658+len(l3)} 人**；L6 端约 **{642+len(l6)} 人**；配对 L6 轨迹约 **{638+len(paired_l6)} 人**。这些只是可用性规模估算，尚未训练或合并模型。

## 关键时间对齐限制

论文方法定义 3-month ctDNA 为术后 **70–112 天**，6-month ctDNA 为 **160–200 天**，但 Supplementary Table 2 只提供两个状态列，没有每位患者的 3m/6m 实际采血日。MRD-window 日数有单独字段。严格的术后≤3.0月（{DAY3:.1f}天）和≤6.0月（{DAY6:.1f}天）可用人数因此**无法从公开表精确识别**：两个采样窗都跨越我们固定的landmark边界。上表按 DFS 日数筛选的是结局与风险时间的可用性近似，不能替代采血日审计。

所以 GALAXY 非常值得加入，但应先二选一：

1. 取得逐患者 3m/6m 采血日，再按现有严格≤3/≤6月规则重审；或
2. 将 landmark 预先改为研究定义的 70–112天 / 160–200天采样窗，并把 Chen、ColonAiQ、COSMOS 也重做同一时间规则。

在完成这个决定前，不把这些行拼进当前冻结 L3/L6 预测，不声称严格跨队列验证。论文提到的 1,664 人是同时有 MRD-window 与后续 surveillance-window 结果的子集，并非 3m+6m trajectory 的样本数。

## 数据与抽取校验

表2解析出 2,240 个唯一连续患者行；已校验 MRD-window 阳性/阴性状态共 {totals['mrd_window_status_available']} 人、3m状态 {totals['ctdna_3m_status_available']} 人、6m状态 {totals['ctdna_6m_status_available']} 人。没有患者行或 ID 写入汇总输出。

原始补充 PDF：`E:/crc_mrd_digital_twin/data_raw/galaxy_2024_supplementary_tables_1_3.pdf`  
审计脚本：`scripts/08_galaxy_patient_level_audit.py`  
汇总表：`outputs/tables/galaxy_patient_level_availability_audit.csv`、`galaxy_timepoint_availability.csv`、`galaxy_nominal_mrd_state_transitions.csv`

## 来源

- [Nature Medicine 2024 论文](https://www.nature.com/articles/s41591-024-03254-6)
- [补充表 PDF](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41591-024-03254-6/MediaObjects/41591_2024_3254_MOESM1_ESM.pdf)
"""
    (REPORTS / "galaxy_patient_level_audit_zh.md").write_text(report, encoding="utf-8")

    audit_table = [
        {"cohort":"GALAXY", "landmark":"L3_nominal", "n":len(l3), "events_to_12m":sum(y12(r)==1 for r in l3), "controls_to_12m":len(l3)-sum(y12(r)==1 for r in l3), "time_rule":"3m study window; exact sample dates unavailable"},
        {"cohort":"GALAXY", "landmark":"L6_nominal", "n":len(l6), "events_to_12m":sum(y12(r)==1 for r in l6), "controls_to_12m":len(l6)-sum(y12(r)==1 for r in l6), "time_rule":"6m study window; exact sample dates unavailable"},
        {"cohort":"GALAXY", "landmark":"L3_to_L6_paired_at_L6", "n":len(paired_l6), "events_to_12m":sum(y12(r)==1 for r in paired_l6), "controls_to_12m":len(paired_l6)-sum(y12(r)==1 for r in paired_l6), "time_rule":"3m+6m status; exact sample dates unavailable"},
    ]
    write_csv(TABLES / "galaxy_nominal_landmark_audit.csv", audit_table)
    print(f"DONE | total={len(records)} | dfs_fields={totals['known_recurrence_flag_and_dfs_days']} | ctDNA3={totals['ctdna_3m_status_available']} | ctDNA6={totals['ctdna_6m_status_available']} | paired={len(both)} | Y12_known_paired={sum(y12(r) is not None for r in both)} events={sum(y12(r)==1 for r in both)} | L3={len(l3)}/{sum(y12(r)==1 for r in l3)} | L6={len(l6)}/{sum(y12(r)==1 for r in l6)} | paired_L6={len(paired_l6)}/{sum(y12(r)==1 for r in paired_l6)}")


if __name__ == "__main__":
    main()
