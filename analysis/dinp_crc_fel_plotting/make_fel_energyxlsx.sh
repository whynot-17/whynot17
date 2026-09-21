#!/usr/bin/env bash
set -euo pipefail

# Excel -> GROMACS gmx sham 2D free-energy landscape
# Usage:
#   bash make_fel_energyxlsx.sh /mnt/c/Users/21634/Desktop/energy.xlsx
#
# Sheet2 (or $SHEET) is expected to contain:
#   column 1 = time
#   column 2 = complex RMSD
#   column 3 = Protein Rg

INPUT="${1:-/mnt/c/Users/21634/Desktop/energy.xlsx}"
SHEET="${SHEET:-Sheet2}"
OUTDIR="${OUTDIR:-$(pwd)/gmx_sham_energyxlsx_output}"
CLEAN="${OUTDIR}/rmsd_rg_for_sham.xvg"
OUTPUT="${OUTDIR}/gibbs.xpm"
EPS="${OUTDIR}/gibbs.eps"
PYTHON="${PYTHON:-python3}"
RUN_GMX="${RUN_GMX:-1}"
TSHAM_K="${TSHAM_K:-310}"

if [[ ! -f "$INPUT" ]]; then
  echo "Error: input workbook not found: $INPUT" >&2
  exit 1
fi

if ! command -v "$PYTHON" >/dev/null 2>&1; then
  echo "Error: Python interpreter not found: $PYTHON" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

echo "Input: $INPUT"
echo "Sheet: $SHEET"

"$PYTHON" - "$INPUT" "$SHEET" "$CLEAN" <<'PY'
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

def native_path(value: str) -> Path:
    # When WSL launches a Windows Python executable, convert /mnt/X/... to X:\\...
    if os.name == "nt" and value.startswith("/mnt/") and len(value) > 6:
        drive = value[5].upper()
        tail = value[7:].replace("/", "\\\\")
        return Path(f"{drive}:\\\\{tail}")
    return Path(value)

input_path = native_path(sys.argv[1])
sheet = sys.argv[2]
clean_path = native_path(sys.argv[3])

raw = pd.read_excel(input_path, sheet_name=sheet, header=None, usecols=[0, 1, 2])
if raw.shape[0] < 10 or raw.shape[1] != 3:
    raise SystemExit(f"Expected at least 10 rows and 3 columns, found {raw.shape}")

raw.columns = ["time", "complex_rmsd", "protein_rg"]
data = raw.apply(pd.to_numeric, errors="coerce")
if data.isna().any().any():
    raise SystemExit("Input contains missing or non-numeric values")
if not np.isfinite(data.to_numpy()).all():
    raise SystemExit("Input contains non-finite values")
if not data["time"].is_monotonic_increasing:
    raise SystemExit("Time column is not monotonically increasing")

with clean_path.open("w", encoding="utf-8", newline="\n") as handle:
    handle.write('@ title "MiNP-PPARG RMSD-Rg free-energy input"\n')
    handle.write('@ xaxis label "Complex RMSD"\n')
    handle.write('@ yaxis label "Protein Rg"\n')
    handle.write('@ TYPE xy\n')
    handle.write('@ s0 legend "Sheet2 coordinates"\n')
    for rmsd, rg in zip(data["complex_rmsd"], data["protein_rg"]):
        handle.write(f"{rmsd:.10g} {rg:.10g}\n")

print(f"Prepared 2D coordinates: X = complex RMSD; Y = Protein Rg")
print(f"Frames: {len(data)}")
print(f"RMSD range: {data['complex_rmsd'].min():.8g} to {data['complex_rmsd'].max():.8g}")
print(f"Rg range: {data['protein_rg'].min():.8g} to {data['protein_rg'].max():.8g}")
PY

if [[ "$RUN_GMX" == "0" ]]; then
  echo "Prepared input only (RUN_GMX=0): $CLEAN"
  exit 0
fi

if ! command -v gmx >/dev/null 2>&1; then
  echo "Error: gmx was not found in PATH." >&2
  echo "The two-column XVG was prepared successfully: $CLEAN" >&2
  echo "Install/load GROMACS first, then rerun this script." >&2
  exit 1
fi

read -r XMIN XMAX YMIN YMAX < <(
  awk '
    !/^[@#]/ && NF >= 2 {
      if (n == 0 || $1 < xmin) xmin = $1
      if (n == 0 || $1 > xmax) xmax = $1
      if (n == 0 || $2 < ymin) ymin = $2
      if (n == 0 || $2 > ymax) ymax = $2
      n++
    }
    END {
      dx = (xmax - xmin) * 1e-6
      dy = (ymax - ymin) * 1e-6
      if (dx == 0) dx = 1e-6
      if (dy == 0) dy = 1e-6
      print xmin - dx, xmax + dx, ymin - dy, ymax + dy
    }
  ' "$CLEAN"
)

echo "Running gmx sham..."
echo "Grid: 220 x 220"
echo "Temperature: ${TSHAM_K} K"
echo "X range: ${XMIN} to ${XMAX}"
echo "Y range: ${YMIN} to ${YMAX}"
gmx sham \
  -f "$CLEAN" \
  -ls "$OUTPUT" \
  -ngrid 220 220 1 \
  -xmin "$XMIN" "$YMIN" 0 \
  -xmax "$XMAX" "$YMAX" 1 \
  -tsham "$TSHAM_K" \
  -nlevels 30 \
  -notime

if [[ -f "$OUTPUT" ]]; then
  # gmx sham writes generic PC1/PC2 metadata for -notime input; replace only
  # those labels so the exported XPM/EPS document the actual CVs.
  sed -i 's/title:   "Gibbs Energy Landscape"/title:   "MiNP-PPARG Gibbs Energy Landscape"/' "$OUTPUT"
  sed -i 's/x-label: "PC1"/x-label: "Complex RMSD (nm)"/' "$OUTPUT"
  sed -i 's/y-label: "PC2"/y-label: "Protein Rg (nm)"/' "$OUTPUT"
  gmx xpm2ps \
    -f "$OUTPUT" \
    -o "$EPS" \
    -rainbow blue
  echo ""
  echo "Done!"
  echo "Prepared input: $CLEAN"
  echo "Generated: $OUTPUT"
  echo "Generated vector export: $EPS"
else
  echo "Failed to generate: $OUTPUT" >&2
  exit 1
fi
