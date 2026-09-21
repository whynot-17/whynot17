#!/usr/bin/env python3
"""Export the MiNP--PPARG FEL and its source trajectory for Origin.

Two different three-column files are produced deliberately:

1. Origin FEL XYZ: X, Y, G
   One row per XPM grid cell, equivalent to xpm2all.bsh -xyz.
2. Origin trajectory: time, RMSD, Rg
   One row per MD frame, retained from the original Excel source.

The XPM grid does not contain frame times; it is a spatial free-energy
surface generated after binning the time series.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


NUMBER_PATTERN = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def read_gromacs_xpm(xpm_file: Path):
    lines = xpm_file.read_text(encoding="utf-8", errors="ignore").splitlines()

    x_axis, y_axis = [], []
    for line in lines:
        if "x-axis:" in line:
            text = line.split("x-axis:", 1)[1].replace("*/", "")
            x_axis.extend(float(v) for v in re.findall(NUMBER_PATTERN, text))
        elif "y-axis:" in line:
            text = line.split("y-axis:", 1)[1].replace("*/", "")
            y_axis.extend(float(v) for v in re.findall(NUMBER_PATTERN, text))

    header_idx = None
    width = height = ncolors = cpp = None
    for i, line in enumerate(lines):
        match = re.search(r'"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"', line)
        if match:
            width, height, ncolors, cpp = map(int, match.groups())
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"未找到 XPM 头信息: {xpm_file}")

    symbol_to_value = {}
    color_lines = lines[header_idx + 1: header_idx + 1 + ncolors]
    for line in color_lines:
        quoted = re.findall(r'"([^"]*)"', line)
        value_match = re.search(r'/\*\s*"([^"]+)"\s*\*/', line)
        if quoted and value_match:
            symbol_to_value[quoted[0][:cpp]] = float(value_match.group(1))

    rows = []
    for line in lines[header_idx + 1 + ncolors:]:
        quoted = re.findall(r'"([^"]*)"', line)
        if not quoted:
            continue
        row_string = quoted[0]
        if len(row_string) != width * cpp:
            continue
        rows.append([
            symbol_to_value[row_string[j:j + cpp]]
            for j in range(0, len(row_string), cpp)
        ])
        if len(rows) == height:
            break

    z = np.asarray(rows, dtype=float)
    if z.shape != (height, width):
        raise ValueError(f"矩阵尺寸错误: {z.shape}; 预期 {(height, width)}")

    x = np.asarray(x_axis, dtype=float)
    y = np.asarray(y_axis, dtype=float)
    # GROMACS stores N+1 bin boundaries for an N-cell matrix.
    if len(x) == width + 1:
        x = (x[:-1] + x[1:]) / 2.0
    if len(y) == height + 1:
        y = (y[:-1] + y[1:]) / 2.0
    if len(x) != width or len(y) != height:
        raise ValueError(
            f"轴长度错误: x={len(x)}, y={len(y)}; 预期 {width}, {height}"
        )

    # XPM image rows run from top to bottom; restore increasing Y for Origin.
    return x, y, np.flipud(z)


def export_fel_xyz(xpm_file: Path, output_file: Path, mask_max: bool) -> int:
    x, y, z = read_gromacs_xpm(xpm_file)
    z_out = z.copy()
    max_energy = float(np.nanmax(z_out))
    if mask_max:
        z_out[z_out >= max_energy] = np.nan

    rows = []
    for iy, y_value in enumerate(y):
        for ix, x_value in enumerate(x):
            rows.append((x_value, y_value, z_out[iy, ix]))

    df = pd.DataFrame(rows, columns=[
        "Complex_RMSD_nm", "Protein_Rg_nm", "Free_energy_kJ_mol"
    ])
    # NaN is intentionally retained for unsampled/high-background cells so
    # Origin can treat them as missing rather than as a real energy value.
    df.to_csv(output_file, sep="\t", index=False, float_format="%.9g",
              na_rep="NaN")
    return len(df)


def export_trajectory(excel_file: Path, output_file: Path) -> int:
    # The frozen source has no header and uses Sheet2: time, complex RMSD, Rg.
    raw = pd.read_excel(excel_file, sheet_name="Sheet2", header=None)
    if raw.shape[1] < 3:
        raise ValueError(f"原始 Excel 少于三列: {excel_file}")
    raw = raw.iloc[:, :3].copy()
    raw.columns = ["Time_ns", "Complex_RMSD_nm", "Protein_Rg_nm"]
    for col in raw.columns:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
    raw = raw.dropna(subset=list(raw.columns)).reset_index(drop=True)
    raw.to_csv(output_file, sep="\t", index=False, float_format="%.9g")
    return len(raw)


def main() -> None:
    if len(sys.argv) not in {3, 4}:
        print(
            "用法: python export_origin_data.py <gibbs.xpm> <energy.xlsx> "
            "[output_dir]"
        )
        raise SystemExit(2)

    xpm_file = Path(sys.argv[1]).resolve()
    excel_file = Path(sys.argv[2]).resolve()
    output_dir = Path(sys.argv[3]).resolve() if len(sys.argv) == 4 else xpm_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    tag = xpm_file.stem
    masked = output_dir / f"{tag}_origin_fel_xyz_masked.txt"
    complete = output_dir / f"{tag}_origin_fel_xyz_all.txt"
    trajectory = output_dir / "energyxlsx_origin_time_rmsd_rg.txt"

    n_masked = export_fel_xyz(xpm_file, masked, mask_max=True)
    n_complete = export_fel_xyz(xpm_file, complete, mask_max=False)
    n_frames = export_trajectory(excel_file, trajectory)

    x, y, z = read_gromacs_xpm(xpm_file)
    print(f"FEL masked: {masked} ({n_masked} grid rows)")
    print(f"FEL complete: {complete} ({n_complete} grid rows)")
    print(f"Trajectory: {trajectory} ({n_frames} frames)")
    print(f"FEL X range: {x.min():.9g} to {x.max():.9g} nm")
    print(f"FEL Y range: {y.min():.9g} to {y.max():.9g} nm")
    print(f"FEL G range: {z.min():.9g} to {z.max():.9g} kJ/mol")


if __name__ == "__main__":
    main()
