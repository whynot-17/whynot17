"""Clean-room construction of the NHANES environmental analyte registry."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def paired_flag(value_col: str, columns: set[str]) -> str:
    if value_col.startswith("URX"):
        candidate = "URD" + value_col[3:] + "LC"
        return candidate if candidate in columns else ""
    if value_col.startswith("LBX"):
        for candidate in ("LBD" + value_col[3:] + "L", "LBD" + value_col[3:] + "LC"):
            if candidate in columns:
                return candidate
    return ""


def build_registry(catalog: pd.DataFrame, xpt_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for item in catalog.itertuples(index=False):
        local = xpt_dir / f"{item.cycle}_{item.data_file}"
        if not local.exists() or local.stat().st_size < 1000:
            continue
        frame = pd.read_sas(local, format="xport", encoding="latin1")
        columns = set(frame.columns)
        weights = [column for column in frame.columns if column.startswith("WTS")]
        preferred = [column for column in weights if column in {"WTSA2YR", "WTSB2YR"}]
        weight_col = preferred[0] if preferred else (weights[0] if weights else "")
        if not weight_col:
            continue

        for value_col in frame.columns:
            if not (value_col.startswith("URX") or value_col.startswith("LBX")) or value_col == "URXUCR":
                continue
            flag_col = paired_flag(value_col, columns)
            values = pd.to_numeric(frame[value_col], errors="coerce")
            weight = pd.to_numeric(frame[weight_col], errors="coerce")
            valid = values.notna() & weight.gt(0)
            flag = pd.to_numeric(frame[flag_col], errors="coerce") if flag_col else pd.Series(np.nan, index=frame.index)
            above = valid & (flag.ne(1) if flag_col else True)
            rows.append({
                "cycle": item.cycle,
                "cycle_begin_year": item.cycle_begin_year,
                "laboratory_title": item.laboratory_title,
                "data_file": item.data_file,
                "data_url": item.data_url,
                "doc_url": item.doc_url,
                "local_xpt": str(local),
                "variable": value_col,
                "flag_variable": flag_col,
                "weight_variable": weight_col,
                "matrix": "urine" if value_col.startswith("URX") else "serum_or_blood",
                "n_measured": int(valid.sum()),
                "n_above_lod": int(above.sum()),
                "above_lod_pct": float(100 * above.sum() / valid.sum()) if valid.sum() else np.nan,
            })
    return pd.DataFrame(rows)
