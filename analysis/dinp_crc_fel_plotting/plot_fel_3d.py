# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) 3D free-energy surface -> current frozen XYZ FEL grid -> native Python render
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from matplotlib import cm, colors
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


def read_grid(path: Path):
    df = pd.read_csv(path, sep="\t")
    x = np.sort(df["Complex_RMSD_nm"].unique())
    y = np.sort(df["Protein_Rg_nm"].unique())
    z = (df.pivot(index="Protein_Rg_nm", columns="Complex_RMSD_nm",
                  values="Free_energy_kJ_mol")
           .reindex(index=y, columns=x)
           .to_numpy(dtype=float))
    return df, x, y, z


def normalized_smooth(z: np.ndarray, x: np.ndarray, y: np.ndarray,
                      sigma: float = 3.2, min_support: float = 0.035):
    """Interpolate locally, then smooth only the sampled region.

    This is display-only interpolation. The source FEL grid is never changed.
    """
    valid = np.isfinite(z)
    if not np.any(valid):
        raise ValueError("FEL grid contains no finite energy values")

    X, Y = np.meshgrid(x, y)
    # Fill local gaps between observed FEL bins so the 3D mesh is not a set of
    # needle-like spikes. Linear interpolation avoids the overshoot seen with
    # cubic interpolation for this sparse energy surface.
    interpolated = griddata(
        (X[valid], Y[valid]), z[valid], (X, Y), method="linear"
    )
    support = gaussian_filter(valid.astype(float), sigma=2.0,
                              mode="nearest")
    keep = np.isfinite(interpolated) & (support > min_support)

    numerator = gaussian_filter(np.where(keep, interpolated, 0.0), sigma=sigma,
                                mode="nearest")
    denominator = gaussian_filter(keep.astype(float), sigma=sigma,
                                  mode="nearest")
    out = np.full_like(z, np.nan, dtype=float)
    final_keep = denominator > min_support
    out[final_keep] = numerator[final_keep] / denominator[final_keep]
    # Linear interpolation can produce tiny numerical overshoot only at the
    # boundary; keep the displayed surface within the observed FEL range.
    out[final_keep] = np.clip(out[final_keep], 0.0, float(np.nanmax(z)))
    return out, keep


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python plot_fel_3d.py gibbs_origin_fel_xyz_masked.txt")

    xyz_file = Path(sys.argv[1]).resolve()
    if not xyz_file.exists():
        raise FileNotFoundError(xyz_file)

    _df, x, y, z = read_grid(xyz_file)
    z_smooth, keep = normalized_smooth(z, x, y)

    # The surface uses the observed FEL energy scale after visualization-only
    # local interpolation and smoothing. The source FEL grid is not refit.
    zmin = float(np.nanmin(z_smooth))
    zmax = float(np.nanmax(z_smooth))
    norm = colors.Normalize(vmin=zmin, vmax=max(zmax, zmin + 0.01))
    cmap = colors.LinearSegmentedColormap.from_list(
        "dinp_fel",
        ["#762A83", "#2166AC", "#4393C3", "#1B7837", "#F1A340", "#B2182B"],
        N=256,
    )

    X, Y = np.meshgrid(x, y)
    Z = np.ma.masked_invalid(z_smooth)
    facecolors = cmap(norm(np.ma.filled(Z, 0.0)))
    facecolors[..., 3] = np.where(np.ma.getmaskarray(Z), 0.0, 0.96)

    fig = plt.figure(figsize=(183 / 25.4, 145 / 25.4), dpi=300)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    surf = ax.plot_surface(
        X, Y, Z,
        facecolors=facecolors,
        rstride=2,
        cstride=2,
        linewidth=0,
        antialiased=True,
        shade=False,
        alpha=0.96,
    )

    # A faint projected contour helps identify the dominant basin without
    # adding a second data source.
    levels = np.linspace(zmin, max(zmax, zmin + 0.01), 8)
    ax.contour(
        X, Y, Z,
        zdir="z",
        offset=zmin,
        levels=levels,
        colors="#6A5A74",
        linewidths=0.35,
        alpha=0.42,
    )

    ax.set_xlabel("Complex RMSD (nm)", labelpad=8)
    ax.set_ylabel("Protein Rg (nm)", labelpad=8)
    ax.set_zlabel("Free Energy (kJ/mol)", labelpad=8)
    ax.set_xlim(float(x.min()), float(x.max()))
    ax.set_ylim(float(y.min()), float(y.max()))
    ax.set_zlim(zmin, max(zmax, zmin + 0.01))
    ax.view_init(elev=31, azim=-125)
    ax.set_title("MiNP–PPARG 3D Gibbs free-energy landscape", pad=16,
                 color=BLACK, fontweight="bold")

    # Reduce visual weight of the 3D panes and keep the background white.
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.pane.set_edgecolor("#D9D4DE")
    ax.grid(False)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.66, pad=0.10, aspect=24)
    cbar.set_label("Free Energy (kJ/mol)", labelpad=8)
    cbar.ax.tick_params(labelsize=7, width=0.5, length=2.5)

    out_base = xyz_file.parent / "MiNP_PPARG_FEL_3D_surface"
    save_cns_figure(fig, str(out_base))
    fig.savefig(f"{out_base}.svg", bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"Input grid: {xyz_file}")
    print(f"Finite raw bins: {np.isfinite(z).sum()} / {z.size}")
    print(f"Rendered energy range: {zmin:.4f} to {zmax:.4f} kJ/mol")
    print(f"PDF: {out_base}.pdf")
    print(f"SVG: {out_base}.svg")
    print(f"PNG: {out_base}.png")


if __name__ == "__main__":
    main()
