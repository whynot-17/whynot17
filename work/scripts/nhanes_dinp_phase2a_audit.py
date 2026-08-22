"""Audit NHANES availability of DINP-axis urinary biomarkers for CRC analysis.

The audit deliberately stops before regression.  It links the existing age >=20
CRC-outcome harmonized frame to the NHANES PHTHTE laboratory files and reports
cycle coverage, analyte availability, CRC-case counts, and detection-limit
status for MiNP, MONP, and MCOP.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ANALYTES = {
    "MiNP": {"value": "URXMNP", "comment": "URDMNPLC"},
    "MONP": {"value": "URXMONP", "comment": "URDMONLC"},
    "MCOP": {"value": "URXCOP", "comment": "URDCOPLC"},
}


def cycle_from_path(path: Path) -> str:
    return path.name.split("_", 1)[0]


def read_phthe(path: Path) -> pd.DataFrame:
    try:
        return pd.read_sas(path, format="xport")
    except ValueError as exc:
        # The 1999-2004 URLs in the existing local cache are HTML/error files,
        # not XPORT data.  Keep them in the audit manifest and continue.
        raise ValueError(f"not_an_xport_file: {path.name}: {exc}") from exc


def audit_cycle(raw: pd.DataFrame, harmonized: pd.DataFrame, cycle: str) -> pd.DataFrame:
    needed = ["SEQN"]
    for spec in ANALYTES.values():
        needed.extend([spec["value"], spec["comment"]])
    present = [column for column in needed if column in raw.columns]
    raw = raw[present].copy()
    raw["SEQN"] = raw["SEQN"].astype("int64")

    # The harmonized CRC frame is already restricted to age >=20 and carries
    # the validated cancer outcome definition used in Phase 2B.
    frame = harmonized[[
        "SEQN", "cycle", "age", "cancer_outcome_available", "cancer_free", "crc_case"
    ]].copy()
    frame["SEQN"] = frame["SEQN"].astype("int64")
    frame = frame[frame["cycle"].eq(cycle)].copy()
    merged = frame.merge(raw, on="SEQN", how="inner", validate="one_to_one")

    rows: list[dict] = []
    for analyte, spec in ANALYTES.items():
        value_col = spec["value"]
        comment_col = spec["comment"]
        if value_col not in merged.columns:
            continue
        value = pd.to_numeric(merged[value_col], errors="coerce")
        comment = pd.to_numeric(merged[comment_col], errors="coerce")
        analytic = value.notna()
        outcome = merged["cancer_outcome_available"].fillna(False).astype(bool)
        eligible = analytic & outcome
        # In the XPORT files, the SAS numeric representation of code 0 is
        # read by pandas as the tiny missing-value sentinel (~5.4e-79), while
        # code 1 remains literal 1.  Therefore the robust interpretation is:
        # nonmissing and not 1 = at/above LOD; 1 = below LOD.
        below_lod = comment.eq(1)
        above_lod = comment.notna() & ~below_lod
        rows.append(
            {
                "cycle": cycle,
                "analyte": analyte,
                "value_variable": value_col,
                "lod_comment_variable": comment_col,
                "nhanes_phthe_rows": int(len(raw)),
                "age20_crc_frame_rows": int(len(frame)),
                "matched_rows": int(len(merged)),
                "analytic_value_n": int(analytic.sum()),
                "exposure_and_crc_outcome_n": int(eligible.sum()),
                "crc_cases": int(merged.loc[eligible, "crc_case"].fillna(False).sum()),
                "cancer_free_controls": int(merged.loc[eligible, "cancer_free"].fillna(False).sum()),
                "above_lod_n": int((above_lod & eligible).sum()),
                "above_lod_pct": float(100.0 * (above_lod & eligible).sum() / eligible.sum()) if eligible.sum() else None,
                "below_lod_n": int((below_lod & eligible).sum()),
                "missing_comment_n": int((comment.isna() & eligible).sum()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    harmonized = pd.read_pickle(args.harmonized)
    all_rows: list[pd.DataFrame] = []
    skipped: list[dict] = []
    files = sorted(args.data_dir.glob("*_PHTHTE.XPT"))
    for path in files:
        cycle = cycle_from_path(path)
        try:
            raw = read_phthe(path)
        except ValueError as exc:
            skipped.append({"file": path.name, "cycle": cycle, "reason": str(exc)})
            continue
        if "SEQN" not in raw.columns:
            skipped.append({"file": path.name, "cycle": cycle, "reason": "missing SEQN"})
            continue
        result = audit_cycle(raw, harmonized, cycle)
        if not result.empty:
            all_rows.append(result)

    by_cycle = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    by_cycle.to_csv(args.outdir / "nhanes_dinp_phase2a_audit_by_cycle.csv", index=False)

    if by_cycle.empty:
        summary = pd.DataFrame()
    else:
        grouped = []
        for analyte, group in by_cycle.groupby("analyte", sort=False):
            grouped.append(
                {
                    "analyte": analyte,
                    "value_variable": ";".join(sorted(group["value_variable"].unique())),
                    "cycles_with_file": ";".join(group["cycle"].tolist()),
                    "n_cycles": int(group["cycle"].nunique()),
                    "nhanes_phthe_rows": int(group["nhanes_phthe_rows"].sum()),
                    "age20_crc_frame_rows": int(group["age20_crc_frame_rows"].sum()),
                    "matched_rows": int(group["matched_rows"].sum()),
                    "analytic_value_n": int(group["analytic_value_n"].sum()),
                    "exposure_and_crc_outcome_n": int(group["exposure_and_crc_outcome_n"].sum()),
                    "crc_cases": int(group["crc_cases"].sum()),
                    "cancer_free_controls": int(group["cancer_free_controls"].sum()),
                    "above_lod_n": int(group["above_lod_n"].sum()),
                    "above_lod_pct": float(100.0 * group["above_lod_n"].sum() / group["exposure_and_crc_outcome_n"].sum())
                    if group["exposure_and_crc_outcome_n"].sum()
                    else None,
                    "below_lod_n": int(group["below_lod_n"].sum()),
                    "missing_comment_n": int(group["missing_comment_n"].sum()),
                }
            )
        summary = pd.DataFrame(grouped)
    summary.to_csv(args.outdir / "nhanes_dinp_phase2a_audit_summary.csv", index=False)

    manifest = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).as_posix()),
        "harmonized_input": str(args.harmonized),
        "data_dir": str(args.data_dir),
        "analytes": ANALYTES,
        "outcome_frame": "existing Phase 2B harmonized NHANES CRC frame, age >=20",
        "lod_rule": "NHANES comment code 0 = at or above LOD; 1 = below LOD",
        "official_source": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm",
        "skipped_files": skipped,
    }
    (args.outdir / "nhanes_dinp_phase2a_audit_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    lines = [
        "# DINP axis：NHANES Phase 2A availability audit",
        "",
        "本审计只回答能否进入 CRC 人群分析，不进行回归。分析框架为已有 Phase 2B 的 20 岁以上 CRC outcome frame 与 PHTHTE urine laboratory files 按 SEQN 连接。",
        "",
        "NHANES comment code 按官方定义解释：0 = at or above detection limit，1 = below detection limit。",
        "",
        "官方变量定义与 LOD 规则：[CDC NHANES PHTHTE_J codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm)。",
        "",
        "## Summary",
        "",
    ]
    if summary.empty:
        lines.append("没有可读的 PHTHTE XPORT 文件。")
    else:
        lines += [
            "| Analyte | Variable | Cycles | Exposure + CRC outcome | CRC cases | Above LOD | Above LOD % |",
            "| --- | --- | --- | ---: | ---: | ---: | ---: |",
        ]
        for row in summary.itertuples(index=False):
            lines.append(
                f"| {row.analyte} | {row.value_variable} | {row.cycles_with_file} | "
                f"{row.exposure_and_crc_outcome_n} | {row.crc_cases} | {row.above_lod_n} | "
                f"{row.above_lod_pct:.2f}% |"
            )
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "Exposure + CRC outcome 是可用于当前 CRC outcome 分析的最低数据条件，不等同于足够的统计 power。MONP 若只覆盖单一 NHANES cycle，病例数必须先过稀疏性门槛，再决定是否运行回归。",
        "",
        "跳过的文件和逐 cycle 明细见同目录 CSV/JSON。",
    ]
    (args.outdir / "nhanes_dinp_phase2a_audit_report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
