from __future__ import annotations

"""Create frozen public-analysis figures without fitting or selecting a new model.

The script reproduces the GALAXY patient-level bootstrap OOB predictions from
the already frozen L3/L6 analysis, selects four anonymous real-patient examples,
and renders aggregate/public figures. A separate ID key is written only to the
project's local_only directory and must not be published.
"""

import csv
import html
import importlib.util
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


PROJECT = Path(r"E:\crc_mrd_digital_twin")
PRIMARY = PROJECT / "analysis_runs" / "galaxy_2024_l6_primary_182d6"
OUTPUT = PROJECT / "analysis_runs" / "public_analysis_freeze_20261001"
PUBLIC_OUT = Path(r"C:\Users\ASUS\Documents\Codex\2026-09-29\files-pasted-by-the-user-mrd\outputs\crc_mrd_public_analysis_freeze_20261001")
SCRIPT = PROJECT / "scripts" / "11_frozen_public_analysis_figures.py"
FONT = Path(r"C:\Windows\Fonts\msyh.ttc")
FONT_BOLD = Path(r"C:\Windows\Fonts\msyhbd.ttc")
REPS = 1000
SEED = 20260930

COLORS = {
    "ink": "#18303B", "muted": "#5E7078", "line": "#D8E1E4",
    "paper": "#FFFFFF", "soft": "#F4F7F7", "teal": "#16836B",
    "teal_light": "#E5F3EE", "coral": "#CB554B", "coral_light": "#FBECEA",
    "blue": "#477DA8", "blue_light": "#EAF1F7", "gold": "#D18A28",
    "gold_light": "#FBF2E3", "dark": "#24343B", "grey": "#A2B0B5",
}


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class Canvas:
    """Emit matching high-resolution PNG and editable SVG drawings."""

    def __init__(self, width: int, height: int):
        self.width, self.height = width, height
        self.image = Image.new("RGB", (width, height), COLORS["paper"])
        self.draw = ImageDraw.Draw(self.image)
        self.svg = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            f'<rect width="100%" height="100%" fill="{COLORS["paper"]}"/>',
            '<style>text{font-family:"Microsoft YaHei","Noto Sans CJK SC",Arial,sans-serif}</style>',
        ]

    def font(self, size, bold=False):
        path = FONT_BOLD if bold and FONT_BOLD.exists() else FONT
        return ImageFont.truetype(str(path), int(size))

    def rect(self, x, y, w, h, fill, radius=0, outline=None, width=1):
        if radius:
            self.draw.rounded_rectangle((x, y, x + w, y + h), radius=radius, fill=fill, outline=outline, width=width)
            tag = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" fill="{fill}"'
        else:
            self.draw.rectangle((x, y, x + w, y + h), fill=fill, outline=outline, width=width)
            tag = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}"'
        if outline:
            tag += f' stroke="{outline}" stroke-width="{width}"'
        self.svg.append(tag + '/>')

    def line(self, x1, y1, x2, y2, fill, width=2):
        self.draw.line((x1, y1, x2, y2), fill=fill, width=width)
        self.svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{fill}" stroke-width="{width}"/>')

    def circle(self, cx, cy, r, fill, outline=None, width=1):
        self.draw.ellipse((cx-r, cy-r, cx+r, cy+r), fill=fill, outline=outline, width=width)
        tag = f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}"'
        if outline:
            tag += f' stroke="{outline}" stroke-width="{width}"'
        self.svg.append(tag + '/>')

    def text(self, x, y, value, size=28, color=None, bold=False, anchor="left"):
        color = color or COLORS["ink"]
        font = self.font(size, bold)
        box = self.draw.textbbox((0, 0), str(value), font=font)
        tw = box[2] - box[0]
        if anchor == "middle":
            px = x - tw / 2
            svg_anchor = "middle"
        elif anchor == "right":
            px = x - tw
            svg_anchor = "end"
        else:
            px = x
            svg_anchor = "start"
        self.draw.text((px, y), str(value), fill=color, font=font)
        weight = "700" if bold else "400"
        self.svg.append(f'<text x="{x}" y="{y + size}" text-anchor="{svg_anchor}" font-size="{size}" font-weight="{weight}" fill="{color}">{html.escape(str(value))}</text>')

    def save(self, stem: Path):
        stem.parent.mkdir(parents=True, exist_ok=True)
        self.image.save(stem.with_suffix(".png"), dpi=(300, 300), optimize=True)
        self.svg.append("</svg>")
        stem.with_suffix(".svg").write_text("\n".join(self.svg), encoding="utf-8")


def pct(x, digits=1):
    return f"{100*x:.{digits}f}%"


def safe_ci(row, low, high):
    return float(row[low]), float(row[high])


def select_cases(rows, oob_before, oob_after):
    groups = {
        "negative_to_negative": (0, 0, "低风险维持", 0),
        "negative_to_positive": (0, 1, "转阳后风险上升", 1),
        "positive_to_negative": (1, 0, "转阴后风险下降", 0),
        "positive_to_positive": (1, 1, "持续阳性，高风险", 1),
    }
    chosen = []
    for code, (l3, l6, desc, preferred_y) in groups.items():
        candidates = []
        for i, row in enumerate(rows):
            if (row.get("l3_status") != l3 or row.get("l6_status") != l6
                    or not np.isfinite(oob_before[i]) or not np.isfinite(oob_after[i])):
                continue
            candidates.append((i, row))
        preferred = [(i, r) for i, r in candidates if r["y"] == preferred_y]
        if preferred:
            candidates = preferred
        target = float(np.median([oob_after[i] for i, _ in candidates]))
        if code == "positive_to_negative":
            nn = next((c for c in chosen if c["transition"] == "negative_to_negative"), None)
            if nn:
                elevated = [(i, r) for i, r in candidates if oob_after[i] > nn["risk_after"]]
                if elevated:
                    candidates = elevated
        i, row = min(candidates, key=lambda pair: abs(oob_after[pair[0]] - target))
        event_day = row.get("event_time_days") if row["y"] == 1 else None
        chosen.append({
            "case": chr(ord("A") + len(chosen)), "transition": code, "description": desc,
            "l3_status": l3, "l6_status": l6, "risk_before": float(oob_before[i]),
            "risk_after": float(oob_after[i]), "event": int(row["y"]),
            "event_day": float(event_day) if event_day is not None else None,
            "patient_id": str(row["patient_id"]), "risk_origin_day": 182.62125,
        })
    return chosen


def trajectory_figure(cases, out_stem):
    c = Canvas(2400, 1570)
    c.text(105, 52, "四类真实患者：MRD 轨迹如何更新同一未来结局窗风险", 48, COLORS["ink"], True)
    c.text(108, 119, "GALAXY 开发队列｜L3 信息 → 加入 L6 MRD → 预测 L6 后至术后 12 个月临床复发", 27, COLORS["muted"])
    c.text(108, 166, "患者级风险为冻结模型的 1,000 次 bootstrap OOB 预测；病例匿名展示，仅作轨迹示例。", 23, COLORS["muted"])

    c.rect(78, 240, 2244, 104, COLORS["soft"], radius=20)
    c.text(125, 270, "病例 / 结局", 27, COLORS["muted"], True)
    c.text(770, 270, "Z₃：L3 MRD", 27, COLORS["muted"], True, "middle")
    c.text(1130, 270, "Z₆：L6 MRD", 27, COLORS["muted"], True, "middle")
    c.text(1535, 270, "L3 时风险", 27, COLORS["muted"], True, "middle")
    c.text(1980, 270, "纳入 L6 后风险", 27, COLORS["muted"], True, "middle")

    transition_names = {
        "negative_to_negative": "阴 → 阴",
        "negative_to_positive": "阴 → 阳",
        "positive_to_negative": "阳 → 阴",
        "positive_to_positive": "阳 → 阳",
    }
    row_y = [430, 700, 970, 1240]
    for idx, (case, y) in enumerate(zip(cases, row_y)):
        bg = "#FFFFFF" if idx % 2 == 0 else "#FBFCFC"
        c.rect(78, y - 82, 2244, 225, bg, radius=18, outline=COLORS["line"], width=2)
        c.text(125, y - 49, f"病例 {case['case']}  ·  {transition_names[case['transition']]}", 31, COLORS["ink"], True)
        if case["event"]:
            outcome = "12 个月内复发"
            if case["event_day"]:
                outcome += f"（术后第 {case['event_day']:.0f} 天）"
            outcome_color = COLORS["coral"]
        else:
            outcome = "观察至 12 个月：未见复发"
            outcome_color = COLORS["teal"]
        c.text(125, y + 1, outcome, 24, outcome_color, True)

        for x, status in ((770, case["l3_status"]), (1130, case["l6_status"])):
            fill = COLORS["coral"] if status else COLORS["teal"]
            c.circle(x, y + 5, 43, fill)
            c.text(x, y - 16, "+" if status else "−", 43, "#FFFFFF", True, "middle")
        c.line(830, y + 5, 1058, y + 5, COLORS["grey"], 4)
        c.text(945, y - 27, "MRD 更新", 18, COLORS["muted"], False, "middle")

        for x, risk, fill in ((1535, case["risk_before"], COLORS["blue"]),
                              (1980, case["risk_after"], COLORS["gold"] if case["risk_after"] < case["risk_before"] else COLORS["coral"])):
            c.text(x, y - 48, pct(risk), 38, fill, True, "middle")
            c.rect(x - 126, y + 22, 252, 21, COLORS["line"], radius=10)
            c.rect(x - 126, y + 22, max(4, 252 * risk), 21, fill, radius=10)
        c.line(1677, y + 32, 1830, y + 32, COLORS["grey"], 4)
        c.text(1755, y - 12, "更新", 19, COLORS["muted"], True, "middle")
        if idx < 3:
            c.line(105, y + 159, 2295, y + 159, COLORS["line"], 1)

    c.text(115, 1500, "风险数值均指向同一目标：(名义 L6 day 182.6, 术后 day 365.24] 内临床复发。", 22, COLORS["muted"])
    c.text(2280, 1500, "L3 / L6 状态：− 阴性　+ 阳性", 22, COLORS["muted"], False, "right")
    c.save(out_stem)


def memory_figure(mem, out_stem):
    c = Canvas(2400, 1420)
    c.text(105, 52, "分子记忆：当前 L6 阴性，不代表风险相同", 50, COLORS["ink"], True)
    c.text(108, 122, "GALAXY 开发队列｜两组当前均为 ctDNA 阴性；按 L6 前是否曾观察到阳性分层", 27, COLORS["muted"])
    c.rect(85, 205, 2230, 105, COLORS["teal_light"], radius=18)
    c.text(120, 230, "+ → −  ≠  − → −", 39, COLORS["teal"], True)
    c.text(560, 247, "既往阳性、现已转阴者，L6 后至 12 月的复发率仍较高", 27, COLORS["ink"], True)

    # Group labels and counts.
    mem_by_group = {r["memory_group"]: r for r in mem}
    groups = []
    for key, label, color in (
        ("never_positive", "此前未观察到阳性 → 当前阴性", COLORS["teal"]),
        ("previously_positive_now_negative", "既往曾阳性 → 当前阴性", COLORS["coral"]),
    ):
        row = mem_by_group[key]
        groups.append({
            "name": label, "n": int(row["n"]), "events": int(row["events"]),
            "rate": float(row["observed_recurrence_rate"]),
            "low": float(row["observed_rate_ci_low"]), "high": float(row["observed_rate_ci_high"]),
            "mean": float(row["mean_predicted_risk"]), "median": float(row["median_predicted_risk"]),
            "color": color,
        })
    yvals = [520, 820]
    for g, y in zip(groups, yvals):
        c.text(130, y-70, g["name"], 31, COLORS["ink"], True)
        c.text(130, y-22, f"N={g['n']}　复发 {g['events']} 例", 25, COLORS["muted"])
        c.text(130, y+28, "当前状态：L6 ctDNA −", 23, COLORS["muted"])

    # Observed recurrence rate panel, x-axis 0-30%.
    x0, x1 = 885, 1535
    c.text((x0+x1)//2, 366, "观察复发率（Wilson 95% CI）", 26, COLORS["muted"], True, "middle")
    for tick in range(0, 31, 5):
        x = x0 + (x1-x0)*tick/30
        c.line(x, 430, x, 1010, COLORS["line"], 1)
        c.text(x, 1016, f"{tick}%", 18, COLORS["muted"], False, "middle")
    for g, y in zip(groups, yvals):
        lo, hi, rate = g["low"], g["high"], g["rate"]
        scale = (x1-x0)/0.30
        c.line(x0 + lo*scale, y, x0 + hi*scale, y, g["color"], 8)
        c.line(x0 + lo*scale, y-14, x0 + lo*scale, y+14, g["color"], 4)
        c.line(x0 + hi*scale, y-14, x0 + hi*scale, y+14, g["color"], 4)
        c.circle(x0 + rate*scale, y, 14, g["color"], outline="#FFFFFF", width=3)
        c.text(x0 + rate*scale + 26, y-48, f"{pct(rate)}  ({pct(lo)}–{pct(hi)})", 24, g["color"], True)

    # Mean OOB predicted risk panel, x-axis 0-25%.
    p0, p1 = 1640, 2240
    c.text((p0+p1)//2, 366, "平均 OOB 预测风险", 26, COLORS["muted"], True, "middle")
    for tick in range(0, 26, 5):
        x = p0 + (p1-p0)*tick/25
        c.line(x, 430, x, 1010, COLORS["line"], 1)
        c.text(x, 1016, f"{tick}%", 18, COLORS["muted"], False, "middle")
    for g, y in zip(groups, yvals):
        x = p0 + g["mean"]*(p1-p0)/0.25
        c.line(p0, y, x, y, g["color"], 9)
        c.circle(x, y, 15, g["color"], outline="#FFFFFF", width=3)
        c.text(x + 26, y-48, f"{pct(g['mean'])}  ·  中位数 {pct(g['median'])}", 24, g["color"], True)

    c.rect(130, 1132, 2115, 155, COLORS["soft"], radius=18)
    observed_diff = groups[1]["rate"] - groups[0]["rate"]
    mean_diff = groups[1]["mean"] - groups[0]["mean"]
    c.text(175, 1157, f"观察复发率绝对差：+{100*observed_diff:.1f} 个百分点", 29, COLORS["ink"], True)
    c.text(1260, 1157, f"平均 OOB 风险绝对差：+{100*mean_diff:.1f} 个百分点", 29, COLORS["ink"], True)
    c.text(175, 1210, "解释为预测分层与既往 ctDNA 阳性史有关；不是治疗效应或因果效应。", 22, COLORS["muted"])
    c.text(175, 1250, "“未观察到阳性”仅描述可用采样史，不能证明未采样时段真实阴性；组间差异未作独立验证队列推断。", 20, COLORS["muted"])
    c.save(out_stem)


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    PUBLIC_OUT.mkdir(parents=True, exist_ok=True)
    analysis = load_module(PROJECT / "scripts" / "10_galaxy_l3_l6_l12_cohort_nodes.py", "frozen_analysis10")
    cohorts = analysis.source_patients(PROJECT)
    records, _, _, _, _ = analysis.analysis_records(cohorts)
    galaxy = records["GALAXY"]
    if len(galaxy) != 1081 or sum(r["y"] for r in galaxy) != 103:
        raise RuntimeError(f"Frozen GALAXY risk set mismatch: N={len(galaxy)}, events={sum(r['y'] for r in galaxy)}")
    legacy = analysis.load_legacy(PROJECT)
    _, _, oob_pred, oob_n, _ = analysis.bootstrap_internal(galaxy, legacy, REPS, SEED)
    if np.any(oob_n == 0):
        raise RuntimeError("At least one patient has no OOB prediction; do not render case figure")
    cases = select_cases(galaxy, oob_pred["L3_before"], oob_pred["L6_after"])
    if len(cases) != 4:
        raise RuntimeError("Could not select exactly four real representative trajectories")

    # Figure values and anonymous case table can be public; identifier key stays local.
    public_cases = [{k: v for k, v in row.items() if k not in {"patient_id", "event_day"}} |
                    {"outcome": "recurrence_by_12m" if row["event"] else "no_recurrence_through_12m",
                     "recurrence_day_postoperative": row["event_day"]}
                    for row in cases]
    public_fields = ["case", "transition", "description", "l3_status", "l6_status", "risk_before", "risk_after", "event", "outcome", "recurrence_day_postoperative", "risk_origin_day"]
    local_key = [{"case": r["case"], "patient_id": r["patient_id"], "transition": r["transition"]} for r in cases]
    write_csv(OUTPUT / "outputs" / "tables" / "anonymous_case_trajectory_values.csv", public_cases, public_fields)
    write_csv(OUTPUT / "outputs" / "local_only" / "case_id_key.csv", local_key, ["case", "patient_id", "transition"])

    table_path = PRIMARY / "outputs" / "tables" / "molecular_memory_by_cohort.csv"
    with table_path.open(encoding="utf-8-sig", newline="") as f:
        memory_all = list(csv.DictReader(f))
    memory = [r for r in memory_all if r["cohort"] == "GALAXY"]
    if len(memory) != 2:
        raise RuntimeError("Expected two frozen GALAXY molecular-memory groups")
    (OUTPUT / "outputs" / "tables").mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT / "outputs" / "tables" / "galaxy_molecular_memory_values.csv", memory, list(memory[0].keys()))

    for base in (OUTPUT / "outputs" / "figures", PUBLIC_OUT / "figures"):
        trajectory_figure(cases, base / "galaxy_four_case_l3_l6_risk_update")
        memory_figure(memory, base / "galaxy_molecular_memory_main")
    for sub in ("reports", "outputs/tables", "tables"):
        (PUBLIC_OUT / sub).mkdir(parents=True, exist_ok=True)
    # Copy public aggregate tables into the user-facing deliverable folder.
    for filename in ("anonymous_case_trajectory_values.csv", "galaxy_molecular_memory_values.csv"):
        source = OUTPUT / "outputs" / "tables" / filename
        (PUBLIC_OUT / "tables" / filename).write_bytes(source.read_bytes())
    case_rows = "\n".join(
        f"| {r['case']} | {r['description']} | {pct(r['risk_before'])} | {pct(r['risk_after'])} | "
        f"{'12个月内复发' if r['event'] else '观察至12个月未复发'} |" for r in cases
    )
    report = f"""# 公共数据主分析冻结记录（2026-10-01）

## 冻结范围

- **GALAXY 2024**：开发队列；按原研究定义的 L3（70–112 天）和 L6（160–200 天）节点。
- **COSMOS**：主外部验证；使用研究提供的术后第 3、6 个月节点。
- **ColonAiQ**：支持性外部验证；按其纵向采样记录选取约 L3、L6 节点，事件数很少。
- **Chen**：仅作 L6 单节点验证；该研究的术后纵向采样从第 6 个月开始，因此不进入 L3→L6 更新分析。

主要估计目标保持不变：在同一 L6 后至术后 12 个月临床复发结局窗内，比较临床+L3 历史与再加入 L6 MRD 信息后的预测。GALAXY 的 L6 名义风险起点为术后 day 182.6；其补充资料没有逐患者采血日，研究定义的 MRD 状态来自 160–200 天窗口。外部队列按各自研究规定的 6 个月节点或实际评估时间处理，不强加跨队列同日阈值。

**冻结主结论：L3 MRD 信息已对未来复发有预测价值；在同一未来结局窗中，纳入新获得的 L6 MRD 后，预测表现进一步改善。** 本轮仅形成病例展示图和分子记忆图，没有新增候选算法、重新选模或依据外部队列调参。

## 冻结性能锚点

| 队列与角色 | 同一结局窗比较 | AUC | AP | Brier |
|---|---|---:|---:|---:|
| GALAXY，OOB 开发评估，N=1081，事件=103 | 临床+L3 → 临床+L3/L6 | 0.753 → 0.789 | 0.383 → 0.482 | 0.064 → 0.056 |
| COSMOS，主外部验证，N=278，事件=14 | 冻结模型的 L3 → L6 更新 | 0.764 → 0.838 | 0.250 → 0.368 | 0.047 → 0.039 |
| ColonAiQ，支持性验证，N=40，事件=3 | 冻结模型的 L3 → L6 更新 | 0.649 → 0.712 | 0.424 → 0.530 | 0.066 → 0.066 |

GALAXY 的 L3→L6 配对差值：AUC +0.036（95% bootstrap CI 0.002–0.073），AP +0.100（0.011–0.191），Brier −0.009（−0.016 至 −0.002）。COSMOS 和 ColonAiQ 的事件数有限，按队列角色解释，不将支持性结果写作确定性结论。

## 患者级数字孪生展示

主图展示四个真实 GALAXY 个案的 **Z₃ → Z₆ → 风险更新**。风险点由原冻结模型的 1,000 次 bootstrap OOB 预测生成，预测结局均为名义 day 182.6 至术后 day 365.24 的临床复发。图中仅用病例 A–D，不公开患者编号；编号映射留在项目 `outputs/local_only/`。

| 病例 | 轨迹 | L3 信息下风险 | 加入 L6 后风险 | 观察结局 |
|---|---|---:|---:|---|
{case_rows}

个案用于说明纵向风险如何更新，不替代队列层面的性能估计或验证。个案结局为观察结果，图中风险不能解释为个体因果效应。

## 分子记忆主结果

GALAXY 中两组 L6 ctDNA 当前均阴性。既往未观察到阳性者 N=933，38 例复发，观察复发率 4.1%（Wilson 95% CI 3.0%–5.5%），平均 OOB 预测风险 4.5%；既往曾观察到阳性、目前转阴者 N=80，14 例复发，观察复发率 17.5%（10.7%–27.3%），平均 OOB 风险 16.4%。观察复发率差为 +13.4 个百分点。

“未观察到阳性”只描述已有检测记录；未采样时段的状态未知。该结果支持在描述和预测中保留既往 MRD 历史，但不证明阳性清除的因果作用。其他外部队列中的既往阳性、当前阴性子组较小，本报告不据此宣称独立复制。

## 文件与数据治理

- `outputs/figures/galaxy_four_case_l3_l6_risk_update.png` 和 `.svg`：四个匿名个案的 L3→L6 风险更新。
- `outputs/figures/galaxy_molecular_memory_main.png` 和 `.svg`：当前 L6 阴性患者中的分子记忆主图。
- `outputs/tables/anonymous_case_trajectory_values.csv`、`galaxy_molecular_memory_values.csv`：不含病例编号的绘图数值。
- `outputs/local_only/case_id_key.csv`：病例编号映射，仅保存于 E 盘本地，不提交到 GitHub。
- `scripts/11_frozen_public_analysis_figures.py`：从冻结模型重建 OOB 示例分数并产出图形；不拟合新模型。

## 公开原始队列来源

- GALAXY：[Nature Medicine 2024](https://doi.org/10.1038/s41591-024-03254-6)
- COSMOS：[公开研究与补充资料](https://pmc.ncbi.nlm.nih.gov/articles/PMC11443202/)
- ColonAiQ：[公开研究与补充资料](https://pmc.ncbi.nlm.nih.gov/articles/PMC10119774/)
- Chen：[公开研究与补充资料](https://pmc.ncbi.nlm.nih.gov/articles/PMC8130394/)
"""
    (OUTPUT / "reports").mkdir(parents=True, exist_ok=True)
    (OUTPUT / "reports" / "public_analysis_freeze_zh.md").write_text(report, encoding="utf-8")
    (PUBLIC_OUT / "reports" / "public_analysis_freeze_zh.md").write_text(report, encoding="utf-8")
    (OUTPUT / "scripts").mkdir(parents=True, exist_ok=True)
    (OUTPUT / "scripts" / SCRIPT.name).write_bytes(Path(__file__).read_bytes())
    SCRIPT.parent.mkdir(parents=True, exist_ok=True)
    SCRIPT.write_bytes(Path(__file__).read_bytes())
    print(json.dumps({"galaxy_n": len(galaxy), "events": int(sum(r["y"] for r in galaxy)),
                      "case_examples": public_cases, "private_case_ids": local_key,
                      "output": str(OUTPUT), "public_output": str(PUBLIC_OUT)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

