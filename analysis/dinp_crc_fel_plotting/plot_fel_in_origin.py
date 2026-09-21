#!/usr/bin/env python3
"""Create an Origin contour plot from the frozen MiNP--PPARG FEL XYZ file."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import originpro as op
from scipy.ndimage import gaussian_filter


def main() -> None:
    if len(sys.argv) not in {2, 3}:
        print(
            "用法: python plot_fel_in_origin.py <fel_xyz.tsv> "
            "[CONTOUR.otpu]"
        )
        raise SystemExit(2)

    xyz_file = Path(sys.argv[1]).resolve()
    if not xyz_file.exists():
        raise FileNotFoundError(xyz_file)
    template = Path(
        sys.argv[2]
        if len(sys.argv) == 3
        else r"F:\origin\CONTOUR.otpu"
    ).resolve()
    if not template.exists():
        raise FileNotFoundError(
            f"Origin contour template not found: {template}. "
            "Pass its path as the second argument."
        )

    df = pd.read_csv(xyz_file, sep="\t")
    x = np.sort(df["Complex_RMSD_nm"].unique())
    y = np.sort(df["Protein_Rg_nm"].unique())
    z = (df.pivot(index="Protein_Rg_nm", columns="Complex_RMSD_nm",
                  values="Free_energy_kJ_mol")
           .reindex(index=y, columns=x)
           .to_numpy(dtype=float))

    # Keep the original XYZ data visible in the Origin project as a worksheet.
    # The matrix sheet is used for the actual contour plot.
    op.set_show(True)
    op.new()
    xyz_book = op.new_book('w', 'MiNP_PPARG_FEL_XYZ')
    xyz_sheet = xyz_book[0]
    xyz_sheet.from_df(df)
    xyz_sheet.name = 'FEL_XYZ'

    matrix = op.new_sheet('m', 'MiNP_PPARG_FEL_Matrix')
    # Display-only smoothing on sampled cells only.  If NaNs mark the
    # unsampled/background bins, use normalized convolution so the high-energy
    # background does not wash out the low-energy basin.
    valid = np.isfinite(z)
    if np.any(~valid):
        background = float(np.nanmax(z))
        numerator = gaussian_filter(np.where(valid, z, 0.0), sigma=2.0,
                                    mode='nearest')
        denominator = gaussian_filter(valid.astype(float), sigma=2.0,
                                      mode='nearest')
        z_plot = np.divide(numerator, denominator,
                           out=np.full_like(numerator, background),
                           where=denominator > 0.05)
        z_plot[denominator <= 0.05] = background
        # Visualization-only contrast stretch: preserve the ordering of
        # energies while restoring the low basin to the 0 kJ/mol end of the
        # color map after smoothing.
        sampled = denominator > 0.05
        lo = float(np.nanmin(z_plot[sampled]))
        hi = float(np.nanmax(z_plot[sampled]))
        if hi > lo:
            z_plot[sampled] = (z_plot[sampled] - lo) / (hi - lo) * background
    else:
        z_plot = gaussian_filter(z, sigma=2.0, mode='nearest')
    matrix.from_np(z_plot)
    matrix.xymap = (float(x.min()), float(x.max()),
                    float(y.min()), float(y.max()))
    # Do not let the matrix long name become an oversized Origin graph title.
    matrix.set_label(0, '', 'L')
    matrix.set_label(0, '', 'C')

    graph = op.new_graph('MiNP_PPARG_FEL_Origin', template=str(template))
    layer = graph[0]
    plot = layer.add_mplot(matrix, z=0, type='contour')
    # Conventional FEL palette: low free energy is blue/purple and high
    # free energy is yellow/red.
    plot.colormap = 'Rainbow Isolum.PAL'
    finite_z = z_plot[np.isfinite(z_plot)]
    # Leave a small headroom above the actual maximum so the top bin is not
    # rendered by Origin as an Above-Max white swatch.
    zmax = float(np.nanmax(finite_z)) + 0.01
    plot.zlevels = {
        'minors': 0,
        'levels': np.linspace(0.0, zmax, 30).tolist(),
    }

    # Origin contour settings: smooth cell transitions and do not fill
    # missing values with the top color of the palette.
    op.lt_exec('set %C -cs 1;')
    op.lt_exec('set %C -cm 0;')

    layer.axis('x').title = 'Complex RMSD (nm)'
    layer.axis('y').title = 'Protein Rg (nm)'
    layer.axis('x').set_limits(float(x.min()), float(x.max()))
    layer.axis('y').set_limits(float(y.min()), float(y.max()))
    layer.rescale('m')

    # Replace the template's placeholder color-scale title with a concise
    # scientific label.  In Origin contour templates the color-scale title
    # is controlled by SPECTRUM1.title$, not by the ordinary Text property.
    layer.lt_exec('SPECTRUM1.title$="Free Energy (kJ/mol)";')

    # The white swatch at the top of the color scale is Origin's
    # above-maximum bin.  Lock that bin to the same red end of the palette so
    # the legend is continuous rather than showing a misleading white square.
    layer.lt_exec('layer.cmap.colorAbove=1;')

    out_dir = xyz_file.parent
    graph.save_fig(str(out_dir / 'MiNP_PPARG_FEL_Origin_smoothed.pdf'), width=2400)
    graph.save_fig(str(out_dir / 'MiNP_PPARG_FEL_Origin_smoothed.png'), width=2400)
    op.save(str(out_dir / 'MiNP_PPARG_FEL_Origin_smoothed.opju'))

    print(f"Origin project: {out_dir / 'MiNP_PPARG_FEL_Origin_smoothed.opju'}")
    print(f"Origin PDF: {out_dir / 'MiNP_PPARG_FEL_Origin_smoothed.pdf'}")
    print(f"Origin PNG: {out_dir / 'MiNP_PPARG_FEL_Origin_smoothed.png'}")
    print(f"Origin template: {template}")
    print(f"grid={z.shape}; x={x.min():.9g}..{x.max():.9g}; "
          f"y={y.min():.9g}..{y.max():.9g}")


if __name__ == '__main__':
    main()
