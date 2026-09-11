from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


SERIES = (
    ("whole_ligand_rmsd_A", "Whole DINP", "#c0392b"),
    ("core_rmsd_A", "Aromatic core", "#2c3e50"),
    ("branch_A_rmsd_A", "Branch A", "#d35400"),
    ("branch_B_rmsd_A", "Branch B", "#16a085"),
)


def rolling_mean(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return values
    return np.convolve(values, np.ones(window, dtype=float) / window, mode="valid")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot whole-DINP, aromatic-core, branch-A, and branch-B RMSD on one axis."
    )
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--smooth-ns",
        type=float,
        default=0.5,
        help="Centered visual smoothing window in ns; raw traces remain visible (default: 0.5).",
    )
    args = parser.parse_args()
    if args.smooth_ns < 0:
        raise ValueError("smooth-ns must be non-negative")

    with args.metrics_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError("Metrics CSV contains no rows")
    missing = [key for key, _, _ in SERIES if key not in rows[0]]
    if missing:
        raise RuntimeError(f"Metrics CSV lacks component columns: {', '.join(missing)}")

    time_ns = np.asarray([float(row["time_ps"]) / 1000.0 for row in rows])
    if len(time_ns) > 1:
        dt_ns = float(np.median(np.diff(time_ns)))
    else:
        dt_ns = 0.0
    window = max(1, int(round(args.smooth_ns / dt_ns))) if dt_ns > 0 else 1

    import matplotlib.pyplot as plt

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(10, 5.4))
    for key, label, color in SERIES:
        values = np.asarray([float(row[key]) for row in rows])
        axis.plot(time_ns, values, color=color, alpha=0.13, lw=0.5)
        smoothed = rolling_mean(values, window)
        offset = (len(values) - len(smoothed)) // 2
        axis.plot(
            time_ns[offset : offset + len(smoothed)],
            smoothed,
            color=color,
            lw=1.8,
            label=label,
        )
    axis.set_title("PPARG–DINP MD: component-resolved RMSD (0–100 ns)")
    axis.set_xlabel("Production time (ns)")
    axis.set_ylabel("RMSD relative to 0 ns (Å)")
    axis.set_xlim(float(time_ns[0]), float(time_ns[-1]))
    axis.grid(alpha=0.25, lw=0.6)
    axis.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(args.output, dpi=220)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
