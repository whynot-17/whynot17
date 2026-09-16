from pathlib import Path
import csv
import json
import math
import warnings
from collections import Counter, defaultdict

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

RUN = Path(r"E:\chatgpt\pparg_minp_md\run_20260915_seed20260917")
TOPO = Path(r"E:\chatgpt\pparg_minp_md\system\seed20260917\initial.pdb")
DCD = RUN / "production.dcd"
CONTACT_CUTOFF = 4.5
SAMPLE_STRIDE = 10  # contact/H-bond sampling every 0.1 ns; RMSD/RMSF use every frame

HYDRO_RESNAMES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO", "CYS"}
HBOND_ELEMENTS = {"N", "O", "S"}


def apply_left_rotation(coords, rotation, mobile_center, reference_center):
    v = coords - mobile_center
    out = np.empty_like(v)
    out[:, 0] = v[:, 0] * rotation[0, 0] + v[:, 1] * rotation[0, 1] + v[:, 2] * rotation[0, 2] + reference_center[0]
    out[:, 1] = v[:, 0] * rotation[1, 0] + v[:, 1] * rotation[1, 1] + v[:, 2] * rotation[1, 2] + reference_center[1]
    out[:, 2] = v[:, 0] * rotation[2, 0] + v[:, 1] * rotation[2, 1] + v[:, 2] * rotation[2, 2] + reference_center[2]
    return out


def rmsd(a, b):
    d = a - b
    return float(np.sqrt(np.mean(np.sum(d * d, axis=1))))


def pair_distances(a, b, box):
    d = a[:, None, :] - b[None, :, :]
    if box is not None and np.all(box > 0):
        d -= np.rint(d / box) * box
    return np.sqrt(np.sum(d * d, axis=2))


def atom_element(atom):
    return str(getattr(atom, "element", "") or "").strip().upper()


def residue_key(atom):
    return f"{atom.resname}{int(atom.resid)}"


def build_donor_pairs(group, ref_positions):
    elems = [atom_element(a) for a in group]
    heavy = [i for i, e in enumerate(elems) if e in HBOND_ELEMENTS]
    hydrogens = [i for i, e in enumerate(elems) if e == "H"]
    pairs = []
    for hi in hydrogens:
        h_atom = group[hi]
        candidates = [i for i in heavy if group[i].resid == h_atom.resid]
        if not candidates:
            continue
        best = min(candidates, key=lambda i: float(np.sum((ref_positions[i] - ref_positions[hi]) ** 2)))
        if float(np.sqrt(np.sum((ref_positions[best] - ref_positions[hi]) ** 2))) <= 1.35:
            pairs.append((best, hi))
    return pairs


def acceptor_indices(group):
    out = []
    for i, atom in enumerate(group):
        e = atom_element(atom)
        if e not in HBOND_ELEMENTS:
            continue
        name = atom.name.strip().upper()
        res = atom.resname.strip().upper()
        if e == "N":
            # Keep histidine/aromatic nitrogens; exclude amide and protonated nitrogens.
            if res not in {"HIS", "HID", "HIE", "HIP", "ASN", "GLN", "TRP"}:
                continue
            if name in {"N", "NE", "NZ", "NH1", "NH2"} and res not in {"HIS", "HID", "HIE", "HIP", "TRP"}:
                continue
        out.append(i)
    return out


def hydrogen_bond_count(prot_pos, lig_pos, box, prot_donors, prot_acceptors, lig_donors, lig_acceptors):
    count = 0
    cos120 = -0.5
    lig_acc = lig_pos[lig_acceptors] if lig_acceptors else np.empty((0, 3))
    prot_acc = prot_pos[prot_acceptors] if prot_acceptors else np.empty((0, 3))
    for donor_i, h_i in prot_donors:
        if len(lig_acc) == 0:
            break
        d = prot_pos[donor_i]
        h = prot_pos[h_i]
        ha = lig_acc - h
        ha -= np.rint(ha / box) * box
        da = lig_acc - d
        da -= np.rint(da / box) * box
        ha_len = np.sqrt(np.sum(ha * ha, axis=1))
        da_len = np.sqrt(np.sum(da * da, axis=1))
        valid = (ha_len <= 2.5) & (da_len <= 3.5) & (ha_len > 1e-8)
        if np.any(valid):
            dh = d - h
            dh_len = math.sqrt(float(np.sum(dh * dh)))
            cosang = np.sum(ha[valid] * dh[None, :], axis=1) / (ha_len[valid] * dh_len)
            count += int(np.sum(cosang <= cos120))
    for donor_i, h_i in lig_donors:
        if len(prot_acc) == 0:
            break
        d = lig_pos[donor_i]
        h = lig_pos[h_i]
        ha = prot_acc - h
        ha -= np.rint(ha / box) * box
        da = prot_acc - d
        da -= np.rint(da / box) * box
        ha_len = np.sqrt(np.sum(ha * ha, axis=1))
        da_len = np.sqrt(np.sum(da * da, axis=1))
        valid = (ha_len <= 2.5) & (da_len <= 3.5) & (ha_len > 1e-8)
        if np.any(valid):
            dh = d - h
            dh_len = math.sqrt(float(np.sum(dh * dh)))
            cosang = np.sum(ha[valid] * dh[None, :], axis=1) / (ha_len[valid] * dh_len)
            count += int(np.sum(cosang <= cos120))
    return count


def format_pdb_coord(line, xyz):
    return line[:30] + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + line[54:]


def write_overlay(ref_pdb, models, out_path):
    protein_lines = []
    ligand_lines = []
    with ref_pdb.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM"):
                protein_lines.append(line.rstrip("\n"))
            elif line.startswith("HETATM") and line[17:20].strip() == "UNK":
                ligand_lines.append(line.rstrip("\n"))
    with out_path.open("w", encoding="utf-8") as f:
        for model_no, (label, coords) in enumerate(models, start=1):
            f.write(f"MODEL     {model_no:4d}    {label}\n")
            for line in protein_lines:
                f.write(line + "\n")
            for i, line in enumerate(ligand_lines):
                f.write(format_pdb_coord(line, coords[i]) + "\n")
            f.write("ENDMDL\n")
        f.write("END\n")


def aggregate_bins(rows, width=5.0):
    bins = defaultdict(list)
    for row in rows:
        start = math.floor(float(row["time_ns"]) / width) * width
        bins[round(start, 2)].append(row)
    fields = [
        "bin_start_ns", "bin_end_ns", "n_samples", "total_contact_residues_mean",
        "total_contact_pairs_mean", "hydrophobic_contact_residues_mean", "hydrophobic_contact_pairs_mean",
        "hbond_count_mean", "pocket_min_distance_mean_A", "pocket_min_distance_min_A",
        "whole_rmsd_mean_A", "core_rmsd_mean_A", "whole_com_mean_A", "core_com_mean_A",
        "protein_backbone_rmsd_mean_A", "pocket_ca_rmsd_mean_A", "initial_contact_retention_mean",
    ]
    output = []
    numeric = fields[3:]
    for start in sorted(bins):
        group = bins[start]
        out = {"bin_start_ns": start, "bin_end_ns": round(start + width, 2), "n_samples": len(group)}
        for key in numeric:
            source = "pocket_min_distance_A" if key == "pocket_min_distance_min_A" else key.replace("_mean", "")
            vals = [float(r[source]) for r in group]
            if key.endswith("_min_A"):
                out[key] = min(vals)
            else:
                out[key] = sum(vals) / len(vals)
        output.append(out)
    return output, fields


def main():
    warnings.filterwarnings("ignore")
    u = mda.Universe(str(TOPO), str(DCD))
    ref = mda.Universe(str(TOPO))
    bb = u.select_atoms("protein and backbone")
    prot_heavy = u.select_atoms("protein and not name H*")
    prot_all = u.select_atoms("protein")
    lig_heavy = u.select_atoms("resname UNK and not name H*")
    lig_all = u.select_atoms("resname UNK")
    ref_bb = ref.select_atoms("protein and backbone")
    ref_prot_heavy = ref.select_atoms("protein and not name H*")
    ref_prot_all = ref.select_atoms("protein")
    ref_lig_heavy = ref.select_atoms("resname UNK and not name H*")
    ref_lig_all = ref.select_atoms("resname UNK")
    if (len(bb), len(prot_heavy), len(lig_heavy)) != (1052, 2112, 21):
        raise RuntimeError(f"unexpected selections: backbone={len(bb)}, protein_heavy={len(prot_heavy)}, ligand_heavy={len(lig_heavy)}")

    ref_bb_xyz = ref_bb.positions.copy()
    ref_prot_heavy_xyz = ref_prot_heavy.positions.copy()
    ref_lig_heavy_xyz = ref_lig_heavy.positions.copy()
    ref_lig_all_xyz = ref_lig_all.positions.copy()
    ref_bb_center = ref_bb_xyz.mean(axis=0)
    ref_core_idx = np.arange(12, 18)
    ref_core_xyz = ref_lig_heavy_xyz[ref_core_idx]
    initial_d = pair_distances(ref_prot_heavy_xyz, ref_core_xyz, None)
    initial_nearest = initial_d.min(axis=1)
    initial_contact_mask = initial_nearest <= CONTACT_CUTOFF
    initial_keys = sorted({residue_key(atom) for atom, flag in zip(ref_prot_heavy, initial_contact_mask) if flag})
    pocket_mask = np.array([residue_key(atom) in initial_keys for atom in ref_prot_heavy], dtype=bool)
    pocket_ref_heavy = ref_prot_heavy_xyz[pocket_mask]
    pocket_ca_ref = ref.select_atoms("protein and name CA and resid " + " ".join(str(x.resid) for x in ref.select_atoms("protein and name CA") if residue_key(x) in initial_keys)).positions.copy()
    pocket_ca_traj = u.select_atoms("protein and name CA and resid " + " ".join(str(x.resid) for x in u.select_atoms("protein and name CA") if residue_key(x) in initial_keys))
    if len(pocket_ca_ref) == 0:
        pocket_ca_ref = ref.select_atoms("protein and name CA").positions.copy()

    prot_donors = build_donor_pairs(ref_prot_all, ref_prot_all.positions.copy())
    lig_donors = build_donor_pairs(ref_lig_all, ref_lig_all.positions.copy())
    prot_acceptors = acceptor_indices(ref_prot_all)
    lig_acceptors = acceptor_indices(ref_lig_all)
    prot_hydro_mask = np.array([atom_element(a) in {"C", "S"} and a.resname in HYDRO_RESNAMES for a in prot_heavy], dtype=bool)
    lig_hydro_mask = np.array([atom_element(a) in {"C", "S"} for a in lig_heavy], dtype=bool)
    prot_keys = [residue_key(a) for a in prot_heavy]
    lig_parts = ["branch_A"] * 9 + ["branch_B"] * 12
    lig_parts[12:18] = ["core"] * 6

    rows = []
    sample_rows = []
    residue_full = Counter()
    residue_late = Counter()
    residue_50 = Counter()
    lig_sum = np.zeros((len(lig_heavy), 3), dtype=float)
    lig_sq = np.zeros((len(lig_heavy), 3), dtype=float)
    lig_internal_sum = np.zeros((len(lig_heavy), 3), dtype=float)
    lig_internal_sq = np.zeros((len(lig_heavy), 3), dtype=float)
    ca_sum = np.zeros((len(pocket_ca_ref), 3), dtype=float)
    ca_sq = np.zeros((len(pocket_ca_ref), 3), dtype=float)
    representative_pool = {}
    last_aligned_all = None
    last_frame_index = None

    for frame_index, ts in enumerate(u.trajectory):
        bb_xyz = bb.positions.copy()
        bb_center = bb_xyz.mean(axis=0)
        rotation, fit_rmsd = align.rotation_matrix(bb_xyz - bb_center, ref_bb_xyz - ref_bb_center)
        box = np.asarray(ts.dimensions[:3], dtype=float)
        if not np.all(box > 0):
            box = np.array([81.072, 81.072, 81.072], dtype=float)
        lig_all_xyz = lig_all.positions.copy()
        shift = np.rint((bb_center - lig_all_xyz.mean(axis=0)) / box) * box
        lig_all_xyz = lig_all_xyz + shift
        aligned_all = apply_left_rotation(lig_all_xyz, rotation, bb_center, ref_bb_center)
        aligned_heavy = aligned_all[:len(lig_heavy)]
        aligned_core = aligned_heavy[ref_core_idx]
        current_ca = pocket_ca_traj.positions.copy()
        aligned_ca = None
        if len(current_ca) == len(pocket_ca_ref):
            aligned_ca = apply_left_rotation(current_ca, rotation, bb_center, ref_bb_center)
        else:
            aligned_ca = np.zeros_like(pocket_ca_ref)

        lig_sum += aligned_heavy
        lig_sq += aligned_heavy * aligned_heavy
        if len(aligned_ca) == len(pocket_ca_ref):
            ca_sum += aligned_ca
            ca_sq += aligned_ca * aligned_ca
        last_aligned_all = aligned_all.copy()
        last_frame_index = frame_index

        core_internal_rotation, _ = align.rotation_matrix(aligned_core - aligned_core.mean(axis=0), ref_core_xyz - ref_core_xyz.mean(axis=0))
        core_internal = apply_left_rotation(aligned_core, core_internal_rotation, aligned_core.mean(axis=0), ref_core_xyz.mean(axis=0))
        internal_heavy = apply_left_rotation(aligned_heavy, core_internal_rotation, aligned_core.mean(axis=0), ref_core_xyz.mean(axis=0))
        lig_internal_sum += internal_heavy
        lig_internal_sq += internal_heavy * internal_heavy
        current_pocket_d = np.sqrt(np.sum((pocket_ref_heavy[:, None, :] - aligned_core[None, :, :]) ** 2, axis=2))
        pocket_min = float(current_pocket_d.min())
        # Initial-contact retention is defined on the six aromatic core atoms, not on flexible tails.
        current_initial_d = pair_distances(prot_heavy.positions.copy(), lig_heavy.positions.copy()[ref_core_idx], box)
        current_nearest_initial = current_initial_d.min(axis=1)
        retention = float(np.mean(current_nearest_initial[initial_contact_mask] <= CONTACT_CUTOFF))
        rows.append({
            "time_ns": float(ts.time / 1000.0),
            "protein_backbone_rmsd_A": float(fit_rmsd),
            "whole_ligand_rmsd_A": rmsd(aligned_heavy, ref_lig_heavy_xyz),
            "core_rmsd_A": rmsd(aligned_core, ref_core_xyz),
            "branch_A_rmsd_A": rmsd(aligned_heavy[:9], ref_lig_heavy_xyz[:9]),
            "branch_B_rmsd_A": rmsd(aligned_heavy[9:21], ref_lig_heavy_xyz[9:21]),
            "whole_com_displacement_A": float(np.sqrt(np.sum((aligned_heavy.mean(axis=0) - ref_lig_heavy_xyz.mean(axis=0)) ** 2))),
            "core_com_displacement_A": float(np.sqrt(np.sum((aligned_core.mean(axis=0) - ref_core_xyz.mean(axis=0)) ** 2))),
            "core_internal_rmsd_A": rmsd(core_internal, ref_core_xyz),
            "pocket_min_distance_A": pocket_min,
            "initial_contact_retention": retention,
            "pocket_ca_rmsd_A": rmsd(aligned_ca, pocket_ca_ref) if len(aligned_ca) == len(pocket_ca_ref) else None,
        })

        if frame_index % SAMPLE_STRIDE == 0:
            # Contact chemistry is evaluated against the complete ligand;
            # initial-contact retention above is intentionally core-only.
            contact_d = pair_distances(prot_heavy.positions.copy(), lig_heavy.positions.copy(), box)
            contact_mask = contact_d <= CONTACT_CUTOFF
            contact_prot_indices = np.where(contact_mask.any(axis=1))[0]
            contact_keys = sorted({prot_keys[i] for i in contact_prot_indices})
            residue_full.update(contact_keys)
            if ts.time / 1000.0 >= 50.0:
                residue_50.update(contact_keys)
            if ts.time / 1000.0 >= 70.0:
                residue_late.update(contact_keys)
            hydro_d = contact_d[np.ix_(prot_hydro_mask, lig_hydro_mask)]
            hydro_mask = hydro_d <= CONTACT_CUTOFF
            hydro_prot_idx = np.where(hydro_mask.any(axis=1))[0]
            hydro_keys = sorted({prot_keys[i] for i in np.where(prot_hydro_mask)[0][hydro_prot_idx]}) if len(hydro_prot_idx) else []
            p_all = prot_all.positions.copy()
            l_all = lig_all.positions.copy()
            hbond_count = hydrogen_bond_count(p_all, l_all, box, prot_donors, prot_acceptors, lig_donors, lig_acceptors)
            sample = {
                "time_ns": float(ts.time / 1000.0),
                "total_contact_residues": len(contact_keys),
                "total_contact_pairs": int(contact_mask.sum()),
                "hydrophobic_contact_residues": len(hydro_keys),
                "hydrophobic_contact_pairs": int(hydro_mask.sum()),
                "hbond_count": hbond_count,
                "pocket_min_distance_A": pocket_min,
                "whole_rmsd_A": rows[-1]["whole_ligand_rmsd_A"],
                "core_rmsd_A": rows[-1]["core_rmsd_A"],
                "whole_com_A": rows[-1]["whole_com_displacement_A"],
                "core_com_A": rows[-1]["core_com_displacement_A"],
                "protein_backbone_rmsd_A": rows[-1]["protein_backbone_rmsd_A"],
                "pocket_ca_rmsd_A": rows[-1]["pocket_ca_rmsd_A"] or 0.0,
                "initial_contact_retention": retention,
                "contact_residues": ";".join(contact_keys),
                "hydrophobic_residues": ";".join(hydro_keys),
            }
            sample_rows.append(sample)
            if 70.0 <= sample["time_ns"] <= 100.0:
                representative_pool[frame_index] = (sample["time_ns"], aligned_all.copy(), aligned_core.copy())

    # Save full continuous metrics.
    full_fields = list(rows[0])
    with (RUN / "trajectory_metrics_0_100ns_complete.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=full_fields); w.writeheader(); w.writerows(rows)
    sample_fields = list(sample_rows[0])
    with (RUN / "contact_samples_0_100ns.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sample_fields); w.writeheader(); w.writerows(sample_rows)
    bins, bin_fields = aggregate_bins(sample_rows, width=5.0)
    with (RUN / "contacts_5ns_0_100ns.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=bin_fields); w.writeheader(); w.writerows(bins)

    # Residue occupancy, with initial key residues retained for comparison.
    n_full = len(sample_rows)
    n_late = sum(1 for r in sample_rows if r["time_ns"] >= 70.0)
    n_50 = sum(1 for r in sample_rows if r["time_ns"] >= 50.0)
    all_keys = sorted(set(residue_full) | set(initial_keys))
    occ_rows = []
    for key in all_keys:
        occ_rows.append({
            "residue": key,
            "initial_docking_contact": key in initial_keys,
            "occupancy_0_100": residue_full.get(key, 0) / n_full,
            "occupancy_50_100": residue_50.get(key, 0) / n_50 if n_50 else 0.0,
            "occupancy_70_100": residue_late.get(key, 0) / n_late if n_late else 0.0,
            "samples_contacted_0_100": residue_full.get(key, 0),
            "samples_contacted_70_100": residue_late.get(key, 0),
        })
    with (RUN / "residue_contact_occupancy.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(occ_rows[0])); w.writeheader(); w.writerows(occ_rows)

    # RMSF after protein-backbone alignment.
    lig_mean = lig_sum / len(rows)
    lig_mean_sq = lig_sq / len(rows)
    lig_rmsf = np.sqrt(np.maximum(0.0, (lig_mean_sq - lig_mean * lig_mean).sum(axis=1)))
    lig_rmsf_rows = []
    for i, atom in enumerate(ref_lig_heavy):
        lig_rmsf_rows.append({"atom_index_1based": i + 1, "atom_name": atom.name, "element": atom_element(atom), "part": lig_parts[i], "rmsf_A": float(lig_rmsf[i])})
    with (RUN / "ligand_heavy_rmsf.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lig_rmsf_rows[0])); w.writeheader(); w.writerows(lig_rmsf_rows)
    lig_internal_mean = lig_internal_sum / len(rows)
    lig_internal_mean_sq = lig_internal_sq / len(rows)
    lig_internal_rmsf = np.sqrt(np.maximum(0.0, (lig_internal_mean_sq - lig_internal_mean * lig_internal_mean).sum(axis=1)))
    lig_internal_rows = []
    for i, atom in enumerate(ref_lig_heavy):
        lig_internal_rows.append({"atom_index_1based": i + 1, "atom_name": atom.name, "element": atom_element(atom), "part": lig_parts[i], "internal_rmsf_A": float(lig_internal_rmsf[i])})
    with (RUN / "ligand_internal_rmsf.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lig_internal_rows[0])); w.writeheader(); w.writerows(lig_internal_rows)
    if len(pocket_ca_ref):
        ca_mean = ca_sum / len(rows)
        ca_mean_sq = ca_sq / len(rows)
        ca_rmsf = np.sqrt(np.maximum(0.0, (ca_mean_sq - ca_mean * ca_mean).sum(axis=1)))
        pocket_ca_atoms = ref.select_atoms("protein and name CA")
        pocket_ca_rows = []
        for i, atom in enumerate(pocket_ca_atoms):
            if residue_key(atom) in initial_keys:
                pos = list(pocket_ca_atoms.resids).index(atom.resid)
                pocket_ca_rows.append({"residue": residue_key(atom), "resid": int(atom.resid), "rmsf_A": float(ca_rmsf[len(pocket_ca_rows)])})
        with (RUN / "pocket_ca_rmsf.csv").open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(pocket_ca_rows[0])); w.writeheader(); w.writerows(pocket_ca_rows)

    # Pocket residency episodes at 0.1 ns sampling.
    residency = []
    # Use every 10 ps frame for residency; contact chemistry remains sampled at 0.1 ns.
    residency_rows = rows
    if residency_rows:
        dt = float(np.median(np.diff([r["time_ns"] for r in residency_rows])))
        start = None; prev = None
        for r in residency_rows:
            is_out = r["pocket_min_distance_A"] > CONTACT_CUTOFF
            if is_out and start is None:
                start = r["time_ns"]
            if not is_out and start is not None:
                residency.append({"start_ns": start, "end_ns": prev + dt, "duration_ns": prev + dt - start})
                start = None
            prev = r["time_ns"]
        if start is not None:
            residency.append({"start_ns": start, "end_ns": prev + dt, "duration_ns": prev + dt - start})
    residency_summary = {
        "definition": "inside original docking pocket when min distance from MiNP core to the six initial-contact residues is <=4.5 A",
        "sample_interval_ns": float(np.median(np.diff([r["time_ns"] for r in residency_rows]))) if len(residency_rows) > 1 else None,
        "fraction_inside": float(np.mean([r["pocket_min_distance_A"] <= CONTACT_CUTOFF for r in residency_rows])) if residency_rows else None,
        "fraction_outside": float(np.mean([r["pocket_min_distance_A"] > CONTACT_CUTOFF for r in residency_rows])) if residency_rows else None,
        "outside_total_sampled_ns": float(sum(x["duration_ns"] for x in residency)),
        "longest_outside_episode_ns": max((x["duration_ns"] for x in residency), default=0.0),
        "farthest_distance_A": max((r["pocket_min_distance_A"] for r in residency_rows), default=None),
        "farthest_time_ns": max(residency_rows, key=lambda r: r["pocket_min_distance_A"])["time_ns"] if residency_rows else None,
        "outside_episodes": residency,
    }
    (RUN / "pocket_residency_summary.json").write_text(json.dumps(residency_summary, indent=2), encoding="utf-8")

    # First excursion and late-window statistics.
    core_vals = [r["core_rmsd_A"] for r in rows]
    core_com_vals = [r["core_com_displacement_A"] for r in rows]
    first_core6 = next((r for r in rows if r["core_rmsd_A"] >= 6.0), None)
    max_50 = max((r for r in rows if r["time_ns"] >= 50.0), key=lambda r: r["core_com_displacement_A"])
    late_rows = [r for r in rows if r["time_ns"] >= 70.0]
    final20 = [r for r in rows if r["time_ns"] >= 80.0]
    def window_stats(group, key):
        vals = [r[key] for r in group]
        return {"mean": float(np.mean(vals)), "sd": float(np.std(vals)), "min": float(np.min(vals)), "max": float(np.max(vals))} if vals else None
    event_summary = {
        "first_core_rmsd_ge_6_A_time_ns": first_core6["time_ns"] if first_core6 else None,
        "first_core_rmsd_ge_6_A_value": first_core6["core_rmsd_A"] if first_core6 else None,
        "largest_core_com_after_50ns_time_ns": max_50["time_ns"],
        "largest_core_com_after_50ns_A": max_50["core_com_displacement_A"],
        "late_70_100": {k: window_stats(late_rows, k) for k in ("whole_ligand_rmsd_A", "core_rmsd_A", "whole_com_displacement_A", "core_com_displacement_A", "pocket_min_distance_A", "initial_contact_retention", "core_internal_rmsd_A")},
        "final_20_80_100": {k: window_stats(final20, k) for k in ("whole_ligand_rmsd_A", "core_rmsd_A", "whole_com_displacement_A", "core_com_displacement_A", "pocket_min_distance_A", "initial_contact_retention", "core_internal_rmsd_A")},
    }
    (RUN / "event_and_late_stability_summary.json").write_text(json.dumps(event_summary, indent=2), encoding="utf-8")

    # Last-30 ns clustering on protein-aligned core coordinates; threshold 2 A.
    cluster_items = [(idx, t, all_xyz, core_xyz) for idx, (t, all_xyz, core_xyz) in representative_pool.items()]
    clusters = []
    for item in cluster_items:
        assigned = False
        for c in clusters:
            center = c[0][3]
            if rmsd(item[3], center) <= 2.0:
                c.append(item); assigned = True; break
        if not assigned:
            clusters.append([item])
    clusters.sort(key=len, reverse=True)
    cluster_rows = []
    centroid_items = []
    for cid, cluster in enumerate(clusters, start=1):
        medoid = min(cluster, key=lambda x: sum(rmsd(x[3], y[3]) for y in cluster))
        centroid_items.append(medoid)
        cluster_rows.append({"cluster": cid, "n_samples": len(cluster), "occupancy": len(cluster) / len(cluster_items), "time_min_ns": min(x[1] for x in cluster), "time_max_ns": max(x[1] for x in cluster), "medoid_time_ns": medoid[1]})
    with (RUN / "last30ns_core_clusters.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cluster_rows[0]) if cluster_rows else ["cluster"]); w.writeheader(); w.writerows(cluster_rows)
    dominant = clusters[0] if clusters else []
    dominant_medoid = min(dominant, key=lambda x: sum(rmsd(x[3], y[3]) for y in dominant)) if dominant else None

    # Representative poses: initial, largest 50-80 ns core-COM excursion, final, dominant last30 medoid.
    intermediate = max((r for r in rows if 50.0 <= r["time_ns"] <= 80.0), key=lambda r: r["core_com_displacement_A"], default=rows[len(rows)//2])
    intermediate_idx = min(range(len(rows)), key=lambda i: abs(rows[i]["time_ns"] - intermediate["time_ns"]))
    final_coords = last_aligned_all
    pool_by_idx = {idx: (t, all_xyz, core_xyz) for idx, t, all_xyz, core_xyz in cluster_items}
    intermediate_coords = None
    for idx, (t, all_xyz, core_xyz) in representative_pool.items():
        if idx == intermediate_idx:
            intermediate_coords = all_xyz
            break
    if intermediate_coords is None:
        nearest_idx = min(representative_pool, key=lambda idx: abs(representative_pool[idx][0] - intermediate["time_ns"])) if representative_pool else None
        intermediate_coords = representative_pool[nearest_idx][1] if nearest_idx is not None else final_coords
    cluster_coords = dominant_medoid[2] if dominant_medoid else final_coords
    models = [("INITIAL", ref_lig_all_xyz), (f"REARRANGED_{intermediate['time_ns']:.2f}ns", intermediate_coords), ("FINAL_99.92ns", final_coords), ("CLUSTER_MEDOID_LAST30", cluster_coords)]
    write_overlay(TOPO, models, RUN / "representative_initial_rearranged_final_cluster_overlay.pdb")
    manifest = {
        "trajectory_frames": len(rows), "trajectory_time_max_ns": rows[-1]["time_ns"],
        "initial_docking_residues": initial_keys,
        "representative_intermediate_time_ns": intermediate["time_ns"],
        "cluster_threshold_core_rmsd_A": 2.0,
        "last30_cluster_count": len(clusters),
        "dominant_last30_cluster_occupancy": len(dominant) / len(cluster_items) if cluster_items else None,
        "files": {
            "metrics": str(RUN / "trajectory_metrics_0_100ns_complete.csv"),
            "contact_samples": str(RUN / "contact_samples_0_100ns.csv"),
            "contacts_5ns": str(RUN / "contacts_5ns_0_100ns.csv"),
            "residue_occupancy": str(RUN / "residue_contact_occupancy.csv"),
            "pocket_residency": str(RUN / "pocket_residency_summary.json"),
            "late_summary": str(RUN / "event_and_late_stability_summary.json"),
            "overlay": str(RUN / "representative_initial_rearranged_final_cluster_overlay.pdb"),
        },
    }
    (RUN / "complete_analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
