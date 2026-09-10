from __future__ import annotations

import csv
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem

from run_4XLD_native_redocking import parse_index_map, parse_pdbqt_models


HERE = Path(__file__).resolve().parent
NATIVE_DIR = HERE / "native_redocking_4XLD"
FOCUSED_DIR = HERE / "native_redocking_4XLD_focused_exh128"
OUT = HERE / "brl_structure_overlap_core_audit"
NATIVE_PDB = NATIVE_DIR / "4XLD_native_BRL.pdb"
NATIVE_SDF = NATIVE_DIR / "4XLD_native_BRL.sdf"
FOCUSED_PDBQT = FOCUSED_DIR / "4XLD_native_redocked_BRL_focused_exh128.pdbqt"

ATOM_NAMES = [
    "S1", "O2", "O4", "O13", "N3", "N16", "N18", "C2", "C4", "C5",
    "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C14", "C15", "C16",
    "C17", "C19", "C20", "C21", "C22",
]

# Pre-declared BRL rigid pharmacophore/core. This is intentionally not tuned per pose.
# 1) thiazolidinedione ring plus its two carbonyl oxygens;
# 2) phenyl ring;
# 3) 2-pyridyl ring.
CORE_ATOM_NAMES = [
    "S1", "O2", "O4", "N3", "C2", "C4", "C5",
    "C7", "C8", "C9", "C10", "C11", "C12",
    "N18", "C17", "C19", "C20", "C21", "C22",
]
CORE_INDICES = np.array([ATOM_NAMES.index(name) for name in CORE_ATOM_NAMES], dtype=int)
FLEXIBLE_INDICES = np.array([i for i in range(len(ATOM_NAMES)) if i not in set(CORE_INDICES)], dtype=int)

POSES = [(1, "Top1", "#d55e00"), (2, "Top2", "#0072b2"), (15, "Mode15", "#009e73")]


def load_native_atoms():
    atoms = []
    for line in NATIVE_PDB.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        atoms.append(
            {
                "name": line[12:16].strip(),
                "element": (line[76:78].strip() or re.sub(r"[^A-Za-z]", "", line[12:16])).upper()[:1],
                "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            }
        )
    if [a["name"] for a in atoms] != ATOM_NAMES:
        raise AssertionError(f"Native atom names do not match the fixed BRL map: {[a['name'] for a in atoms]}")
    return atoms


def load_poses():
    models = parse_pdbqt_models(FOCUSED_PDBQT)
    index_map = parse_index_map(FOCUSED_PDBQT)
    poses = {}
    for mode, label, color in POSES:
        mapped = {}
        for atom in models[mode]:
            input_idx = index_map.get(atom["serial"])
            if input_idx is None:
                raise AssertionError(f"Mode {mode} atom serial {atom['serial']} missing from INDEX MAP")
            if input_idx <= len(ATOM_NAMES):
                mapped[input_idx - 1] = np.asarray(atom["xyz"], dtype=float)
        if sorted(mapped) != list(range(len(ATOM_NAMES))):
            raise AssertionError(f"Mode {mode} does not contain exactly the 25 BRL heavy atoms")
        poses[mode] = {"label": label, "color": color, "xyz": np.vstack([mapped[i] for i in range(len(ATOM_NAMES))])}
    return poses


def kabsch_fit(mobile: np.ndarray, reference: np.ndarray):
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    mobile_centered = mobile - mobile_center
    reference_centered = reference - reference_center
    covariance = mobile_centered.T @ reference_centered
    u, _, vh = np.linalg.svd(covariance)
    rotation = u @ vh
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vh
    transformed = mobile_centered @ rotation + reference_center
    return transformed, rotation, mobile_center, reference_center


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def hex_rgb(value: str):
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def get_font(size: int, bold: bool = False):
    candidates = [
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/segoeuib.ttf" if bold else "C:/Windows/Fonts/segoeui.ttf",
    ]
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def draw_projection(draw, reference, docked, core_indices, x_index, y_index, bounds, origin, size, color, title, subtitle):
    ox, oy = origin
    width, height = size
    xmin, xmax, ymin, ymax = bounds
    margin = 65

    def project(xyz):
        x = ox + margin + (xyz[:, x_index] - xmin) / (xmax - xmin) * (width - 2 * margin)
        y = oy + height - margin - (xyz[:, y_index] - ymin) / (ymax - ymin) * (height - 2 * margin)
        return np.column_stack([x, y])

    ref2 = project(reference)
    dock2 = project(docked)
    for i in range(len(ATOM_NAMES)):
        # Bonds are drawn below atoms; each bond is drawn from the fixed SDF topology.
        pass
    # Bond list is supplied through the closure in render_figure.
    bonds = draw_projection.bonds
    for a, b in bonds:
        draw.line([tuple(ref2[a]), tuple(ref2[b])], fill=(155, 155, 155), width=4)
    for a, b in bonds:
        draw.line([tuple(dock2[a]), tuple(dock2[b])], fill=hex_rgb(color), width=3)
    for i, (x, y) in enumerate(ref2):
        draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=(165, 165, 165), outline=(80, 80, 80))
    core_set = set(int(i) for i in core_indices)
    rgb = hex_rgb(color)
    for i, (x, y) in enumerate(dock2):
        if i in core_set:
            draw.ellipse((x - 5, y - 5, x + 5, y + 5), fill=rgb, outline=(20, 20, 20))
        else:
            draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=(255, 255, 255), outline=rgb, width=2)

    draw.rectangle((ox + margin, oy + margin, ox + width - margin, oy + height - margin), outline=(40, 40, 40), width=2)
    title_font = get_font(24, bold=True)
    small_font = get_font(17)
    draw.text((ox + 12, oy + 8), title, fill=(20, 20, 20), font=title_font)
    draw.text((ox + 12, oy + 38), subtitle, fill=(70, 70, 70), font=small_font)
    draw.text((ox + width / 2 - 40, oy + height - 43), "Å", fill=(40, 40, 40), font=small_font)
    draw.text((ox + 9, oy + height / 2 - 8), "Å", fill=(40, 40, 40), font=small_font)


def render_figure(reference, transformed, metrics, bonds):
    width, height = 1800, 1200
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_projection.bonds = bonds
    all_coords = [reference] + [item["aligned_xyz"] for item in transformed.values()]
    all_coords = np.vstack(all_coords)
    bounds = {}
    for x_index, y_index, key in [(0, 1, "xy"), (0, 2, "xz")]:
        low = all_coords[:, [x_index, y_index]].min(axis=0)
        high = all_coords[:, [x_index, y_index]].max(axis=0)
        pad = 1.5
        bounds[key] = (low[0] - pad, high[0] + pad, low[1] - pad, high[1] + pad)
    panel_w, panel_h = 600, 500
    positions = [(0, 0), (600, 0), (1200, 0), (0, 600), (600, 600), (1200, 600)]
    for row, (x_index, y_index, key) in enumerate([(0, 1, "xy"), (0, 2, "xz")]):
        for col, (mode, label, color) in enumerate(POSES):
            item = transformed[mode]
            row_data = metrics[mode]
            subtitle = (
                f"full RMSD {row_data['full_aligned_rmsd_A']:.2f} Å | "
                f"core RMSD {row_data['core_aligned_rmsd_A']:.2f} Å"
            )
            draw_projection(
                draw,
                reference,
                item["aligned_xyz"],
                CORE_INDICES,
                x_index,
                y_index,
                bounds[key],
                positions[row * 3 + col],
                (panel_w, panel_h),
                color,
                f"{label} · {key.upper()}",
                subtitle,
            )
    # Legend and interpretation key.
    legend_font = get_font(20, bold=True)
    small_font = get_font(17)
    draw.text((30, 1130), "Gray = crystal BRL; colored = docked pose after full-heavy-atom Kabsch alignment; filled colored atoms = pre-defined core; open colored atoms = excluded linker/side-chain atoms.", fill=(30, 30, 30), font=small_font)
    draw.text((30, 1100), "Structure overlap projections (top row XY, bottom row XZ; common axis limits across panels)", fill=(20, 20, 20), font=legend_font)
    return image


def write_overlay_pdb(path: Path, copies, labels, molecule):
    lines = ["REMARK 4XLD BRL structural overlap; coordinates are in receptor frame unless noted", "MODEL        1"]
    serial = 1
    serial_maps = []
    for copy_i, (label, coords) in enumerate(zip(labels, copies)):
        chain = "ABCD"[copy_i]
        serial_map = {}
        lines.append(f"REMARK COPY {chain} {label}")
        for atom_i, (name, xyz) in enumerate(zip(ATOM_NAMES, coords)):
            serial_map[atom_i] = serial
            element = re.sub(r"[^A-Za-z]", "", name).upper()[:1]
            lines.append(
                f"HETATM{serial:5d} {name:>4s} BRL {chain}{502:4d}    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}"
            )
            serial += 1
        lines.append("TER")
        serial_maps.append(serial_map)
    for serial_map in serial_maps:
        for bond in molecule.GetBonds():
            a = serial_map[bond.GetBeginAtomIdx()]
            b = serial_map[bond.GetEndAtomIdx()]
            lines.append(f"CONECT{a:5d}{b:5d}")
    lines.extend(["ENDMDL", "END"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for required in (NATIVE_PDB, NATIVE_SDF, FOCUSED_PDBQT):
        if not required.exists():
            raise FileNotFoundError(required)
    native_atoms = load_native_atoms()
    reference = np.asarray([a["xyz"] for a in native_atoms], dtype=float)
    poses = load_poses()
    transformed = {}
    metrics = {}
    for mode, item in poses.items():
        raw = item["xyz"]
        aligned, rotation, mobile_center, reference_center = kabsch_fit(raw, reference)
        full_direct = rmsd(raw, reference)
        full_aligned = rmsd(aligned, reference)
        core_direct = rmsd(raw[CORE_INDICES], reference[CORE_INDICES])
        core_aligned = rmsd(aligned[CORE_INDICES], reference[CORE_INDICES])
        flex_aligned = rmsd(aligned[FLEXIBLE_INDICES], reference[FLEXIBLE_INDICES])
        metrics[mode] = {
            "mode": mode,
            "label": item["label"],
            "full_receptor_frame_rmsd_A": full_direct,
            "full_aligned_rmsd_A": full_aligned,
            "core_receptor_frame_rmsd_A": core_direct,
            "core_aligned_rmsd_A": core_aligned,
            "excluded_linker_sidechain_aligned_rmsd_A": flex_aligned,
            "core_atom_count": int(len(CORE_INDICES)),
            "excluded_atom_count": int(len(FLEXIBLE_INDICES)),
            "mobile_centroid_minus_crystal_centroid_A": (mobile_center - reference_center).tolist(),
            "core_atom_names": CORE_ATOM_NAMES,
            "excluded_atom_names": [ATOM_NAMES[i] for i in FLEXIBLE_INDICES],
        }
        transformed[mode] = {**item, "aligned_xyz": aligned, "rotation": rotation}

    # Retrieve the original affinity values from the focused QC table if available.
    affinity_csv = FOCUSED_DIR / "4XLD_native_redocking_focused_exh128_rmsd.csv"
    affinities = {}
    if affinity_csv.exists():
        with affinity_csv.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                affinities[int(row["mode"])] = float(row["affinity_kcal_mol"])
    for mode in metrics:
        metrics[mode]["affinity_kcal_mol"] = affinities.get(mode)
        core = metrics[mode]["core_aligned_rmsd_A"]
        full = metrics[mode]["full_aligned_rmsd_A"]
        excluded = metrics[mode]["excluded_linker_sidechain_aligned_rmsd_A"]
        if core <= 2.0 and full <= 2.0 and excluded > core + 0.5:
            interpretation = "core preserved; deviation is concentrated in excluded linker/side-chain atoms"
        elif core <= 2.0 and full <= 2.0:
            interpretation = "core and overall ligand geometry preserved after full-heavy-atom alignment"
        elif core > 2.0:
            interpretation = "core geometry is not preserved; consistent with a core rearrangement/alternative orientation"
        else:
            interpretation = "mixed geometry; inspect receptor-frame and aligned overlays together"
        metrics[mode]["geometric_interpretation"] = interpretation

    # Fixed SDF topology is used for visual bond connectivity.
    molecule = Chem.RemoveHs(Chem.SDMolSupplier(str(NATIVE_SDF), removeHs=False, sanitize=True)[0])
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in molecule.GetBonds()]
    raw_copies = [reference] + [poses[mode]["xyz"] for mode, _, _ in POSES]
    aligned_copies = [reference] + [transformed[mode]["aligned_xyz"] for mode, _, _ in POSES]
    labels = ["Crystal BRL", "Top1", "Top2", "Mode15"]
    write_overlay_pdb(OUT / "4XLD_BRL_crystal_top1_top2_mode15_receptor_frame_overlay.pdb", raw_copies, labels, molecule)
    write_overlay_pdb(OUT / "4XLD_BRL_crystal_top1_top2_mode15_full_heavy_aligned_overlay.pdb", aligned_copies, labels, molecule)

    figure = render_figure(reference, transformed, metrics, bonds)
    figure_path = OUT / "4XLD_BRL_crystal_top1_top2_mode15_full_alignment_overlay.png"
    figure.save(figure_path, format="PNG", optimize=True)

    metrics_csv = OUT / "4XLD_BRL_full_and_core_RMSD.csv"
    fields = [
        "mode", "label", "affinity_kcal_mol", "full_receptor_frame_rmsd_A", "full_aligned_rmsd_A",
        "core_receptor_frame_rmsd_A", "core_aligned_rmsd_A", "excluded_linker_sidechain_aligned_rmsd_A",
        "core_atom_count", "excluded_atom_count", "mobile_centroid_minus_crystal_centroid_A",
        "geometric_interpretation",
    ]
    with metrics_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for mode in [1, 2, 15]:
            row = {key: metrics[mode][key] for key in fields}
            row["mobile_centroid_minus_crystal_centroid_A"] = json.dumps(row["mobile_centroid_minus_crystal_centroid_A"])
            writer.writerow(row)

    report = """# PPARG/4XLD BRL structural overlap and core RMSD audit

## Scope

This audit compares the 4XLD crystal BRL ligand with the focused-box/exhaustiveness-128 PPARG native-redocking poses **Top1, Top2, and Mode15**. It addresses whether the Top1 receptor-frame RMSD of 2.346 Å reflects a local flexible-tail displacement or a whole-ligand flip/rearrangement.

## Coordinate treatments

Two coordinate views are retained:

1. **Receptor-frame view:** crystal and docked coordinates are compared without moving the docked pose. This preserves the ligand's position/orientation relative to the fixed PPARG receptor and is the primary native-redocking QC view.
2. **Full-heavy-atom aligned view:** each docked pose is superposed onto the crystal using a least-squares Kabsch fit over all 25 BRL heavy atoms. Core RMSD is then calculated on the transformed coordinates using the fixed core below. This separates ligand-shape/core preservation from receptor-frame placement.

## Pre-declared BRL core definition

The core is fixed before inspecting the pose-specific values and contains **19 heavy atoms**:

- thiazolidinedione pharmacophore: `S1, C2, N3, C4, C5, O2, O4`;
- phenyl ring: `C7, C8, C9, C10, C11, C12`;
- 2-pyridyl ring: `N18, C17, C19, C20, C21, C22`.

The excluded **6 heavy atoms** are the linker/side-chain atoms: `C6, O13, C14, C15, N16, C16`. They are not removed selectively per pose; the same fixed atom list is used for Crystal, Top1, Top2, and Mode15.

## Results

See `4XLD_BRL_full_and_core_RMSD.csv` for the machine-readable table. `full_aligned_rmsd_A` is the all-heavy-atom fit RMSD. `core_aligned_rmsd_A` is the core RMSD after that same all-heavy-atom fit. `excluded_linker_sidechain_aligned_rmsd_A` is provided only as a diagnostic for localization of deviations.

The PNG uses common axis limits across all panels. Gray is crystal BRL; the colored trace is the docked pose after full-heavy-atom alignment. Filled colored atoms are the pre-defined core, while open colored atoms are the excluded linker/side-chain atoms.

## Structural interpretation rule

- Low full-aligned and core RMSD with higher excluded-atom RMSD supports a preserved core with a localized flexible-tail/linker deviation.
- Core RMSD above 2 Å after full-heavy-atom fitting is treated as evidence against a preserved rigid core and is compatible with core rearrangement or an alternative orientation.
- Receptor-frame RMSD remains separate because full-atom alignment can make a displaced pose look geometrically similar while hiding its position relative to PPARG.

## Outputs

- `4XLD_BRL_crystal_top1_top2_mode15_receptor_frame_overlay.pdb`: raw receptor-frame overlap for molecular viewers
- `4XLD_BRL_crystal_top1_top2_mode15_full_heavy_aligned_overlay.pdb`: all-heavy-atom aligned overlap for shape/core inspection
- `4XLD_BRL_crystal_top1_top2_mode15_full_alignment_overlay.png`: common-scale XY/XZ visual overlap
- `4XLD_BRL_full_and_core_RMSD.csv`: full/core RMSD and interpretation table
- `../analyze_4XLD_BRL_overlap_core.py`: reproducible analysis script

This analysis does not change the docking ranking, does not perform MD, and does not redefine the previously specified ≤2 Å receptor-frame QC threshold.
"""
    report_path = OUT / "4XLD_BRL_overlap_core_RMSD_report.md"
    report_path.write_text(report, encoding="utf-8")

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "native_pdb": str(NATIVE_PDB),
            "native_sdf": str(NATIVE_SDF),
            "focused_redocked_pdbqt": str(FOCUSED_PDBQT),
        },
        "poses": ["Top1", "Top2", "Mode15"],
        "core_atom_names": CORE_ATOM_NAMES,
        "excluded_atom_names": [ATOM_NAMES[i] for i in FLEXIBLE_INDICES],
        "alignment": "Kabsch least-squares fit using all 25 BRL heavy atoms; core RMSD computed after the same transform",
        "outputs": {
            "figure_png": str(figure_path),
            "rmsd_csv": str(metrics_csv),
            "report": str(report_path),
            "receptor_frame_overlay_pdb": str(OUT / "4XLD_BRL_crystal_top1_top2_mode15_receptor_frame_overlay.pdb"),
            "full_heavy_aligned_overlay_pdb": str(OUT / "4XLD_BRL_crystal_top1_top2_mode15_full_heavy_aligned_overlay.pdb"),
        },
        "results": metrics,
    }
    (OUT / "4XLD_BRL_overlap_core_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
