from __future__ import annotations

import csv
import json
import math
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.lib.distances import distance_array


ROOT = Path(r"E:\chatgpt\whynot17\analysis\dinp_vs_minp_pparg_md_audit")
OUT = ROOT / "outputs"
CUTOFF = 4.5
CONTACT_STRIDE = 10  # 0.1 ns for the 10 ps reporter interval
CLUSTER_STRIDE = 10

HYDRO_RESNAMES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO", "CYS"}
HBOND_ELEMENTS = {"N", "O", "S"}
PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

SYSTEMS = {
    "DINP": {
        "ligand_label": "DINP",
        "topology": Path(r"E:\chatgpt\pparg_dinp_md_run_20260909\seed20260909_system\system_topology.pdb"),
        "trajectory": Path(r"E:\chatgpt\pparg_dinp_md_run_20260909\seed20260909\production.dcd"),
        "core": np.arange(0, 6, dtype=int),
        "branch_A": np.arange(6, 18, dtype=int),
        "branch_B": np.arange(18, 30, dtype=int),
    },
    "MiNP": {
        "ligand_label": "MiNP",
        "topology": Path(r"E:\chatgpt\pparg_minp_md\system\seed20260917\initial.pdb"),
        "trajectory": Path(r"E:\chatgpt\pparg_minp_md\run_20260915_seed20260917\production.dcd"),
        # MiNP is asymmetric: the aromatic ring is 12:18, the long ester/alkyl arm is 0:12,
        # and the short carboxyl arm is 18:21. These are disjoint chemical components.
        "core": np.arange(12, 18, dtype=int),
        "branch_A": np.arange(0, 12, dtype=int),
        "branch_B": np.arange(18, 21, dtype=int),
    },
}


def atom_element(atom):
    return str(getattr(atom, "element", "") or "").strip().upper()


def residue_key(atom):
    return f"{atom.resname}{int(atom.resid)}"


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
        if float(np.linalg.norm(ref_positions[best] - ref_positions[hi])) <= 1.35:
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
            if res not in {"HIS", "HID", "HIE", "HIP", "ASN", "GLN", "TRP"}:
                continue
            if name in {"N", "NE", "NZ", "NH1", "NH2"} and res not in {"HIS", "HID", "HIE", "HIP", "TRP"}:
                continue
        out.append(i)
    return out


def hbond_proxy(prot_pos, lig_pos, box, prot_donors, prot_acceptors, lig_donors, lig_acceptors):
    count = 0
    cos120 = -0.5
    lig_acc = lig_pos[lig_acceptors] if lig_acceptors else np.empty((0, 3))
    prot_acc = prot_pos[prot_acceptors] if prot_acceptors else np.empty((0, 3))
    for donor_i, h_i in prot_donors:
        if len(lig_acc) == 0:
            break
        ha = lig_acc - prot_pos[h_i]
        da = lig_acc - prot_pos[donor_i]
        ha -= np.rint(ha / box) * box
        da -= np.rint(da / box) * box
        ha_len = np.linalg.norm(ha, axis=1)
        da_len = np.linalg.norm(da, axis=1)
        valid = (ha_len <= 2.5) & (da_len <= 3.5) & (ha_len > 1e-8)
        if np.any(valid):
            dh = prot_pos[donor_i] - prot_pos[h_i]
            dh_len = float(np.linalg.norm(dh))
            cosang = (ha[valid] @ dh) / (ha_len[valid] * dh_len)
            count += int(np.sum(cosang <= cos120))
    for donor_i, h_i in lig_donors:
        if len(prot_acc) == 0:
            break
        ha = prot_acc - lig_pos[h_i]
        da = prot_acc - lig_pos[donor_i]
        ha -= np.rint(ha / box) * box
        da -= np.rint(da / box) * box
        ha_len = np.linalg.norm(ha, axis=1)
        da_len = np.linalg.norm(da, axis=1)
        valid = (ha_len <= 2.5) & (da_len <= 3.5) & (ha_len > 1e-8)
        if np.any(valid):
            dh = lig_pos[donor_i] - lig_pos[h_i]
            dh_len = float(np.linalg.norm(dh))
            cosang = (ha[valid] @ dh) / (ha_len[valid] * dh_len)
            count += int(np.sum(cosang <= cos120))
    return count


def stats(values):
    values = np.asarray(values, dtype=float)
    return {"mean": float(np.mean(values)), "sd": float(np.std(values)), "min": float(np.min(values)), "max": float(np.max(values))}


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def bin_contacts(samples, width=5.0):
    bins = defaultdict(list)
    for row in samples:
        start = math.floor(row["time_ns"] / width) * width
        bins[round(start, 2)].append(row)
    rows = []
    for start in sorted(bins):
        group = bins[start]
        row = {"bin_start_ns": start, "bin_end_ns": round(start + width, 2), "n_samples": len(group)}
        for key in (
            "total_contact_residues", "total_contact_pairs", "hydrophobic_contact_residues",
            "hydrophobic_contact_pairs", "hbond_proxy", "core_contact_retention", "pocket_min_distance_A",
        ):
            row[f"{key}_mean"] = float(np.mean([x[key] for x in group]))
        rows.append(row)
    return rows


def cluster_last30(items, threshold=2.0):
    clusters = []
    for item in items:
        assigned = False
        for cluster in clusters:
            if rmsd(item[1], cluster[0][1]) <= threshold:
                cluster.append(item)
                assigned = True
                break
        if not assigned:
            clusters.append([item])
    clusters.sort(key=len, reverse=True)
    rows = []
    for i, cluster in enumerate(clusters, start=1):
        medoid = min(cluster, key=lambda x: sum(rmsd(x[1], y[1]) for y in cluster))
        rows.append({
            "cluster": i,
            "n_samples": len(cluster),
            "occupancy": len(cluster) / len(items) if items else 0.0,
            "time_min_ns": min(x[0] for x in cluster),
            "time_max_ns": max(x[0] for x in cluster),
            "medoid_time_ns": medoid[0],
        })
    return rows


def analyze_system(name, cfg):
    out = OUT / name
    out.mkdir(parents=True, exist_ok=True)
    u = mda.Universe(str(cfg["topology"]), str(cfg["trajectory"]))
    ref = mda.Universe(str(cfg["topology"]))
    bb = u.select_atoms("protein and backbone")
    prot_heavy = u.select_atoms("protein and not name H*")
    prot_all = u.select_atoms("protein")
    lig_all = u.select_atoms("resname UNK")
    lig_heavy = u.select_atoms("resname UNK and not name H*")
    ref_bb = ref.select_atoms("protein and backbone")
    ref_prot_heavy = ref.select_atoms("protein and not name H*")
    ref_prot_all = ref.select_atoms("protein")
    ref_lig_all = ref.select_atoms("resname UNK")
    ref_lig_heavy = ref.select_atoms("resname UNK and not name H*")
    if len(lig_heavy) != int(max(cfg["core"].max(), cfg["branch_A"].max(), cfg["branch_B"].max()) + 1):
        raise RuntimeError(f"{name}: unexpected ligand heavy-atom count {len(lig_heavy)}")
    if len(bb) != len(ref_bb) or len(prot_heavy) != len(ref_prot_heavy):
        raise RuntimeError(f"{name}: reference and trajectory selections differ")

    ref_bb_xyz = ref_bb.positions.copy()
    ref_prot_heavy_xyz = ref_prot_heavy.positions.copy()
    ref_lig_heavy_xyz = ref_lig_heavy.positions.copy()
    ref_core_xyz = ref_lig_heavy_xyz[cfg["core"]]
    ref_bb_center = ref_bb_xyz.mean(axis=0)
    ref_core_center = ref_core_xyz.mean(axis=0)

    initial_core_dist = distance_array(ref_prot_heavy_xyz, ref_core_xyz)
    initial_core_nearest = initial_core_dist.min(axis=1)
    initial_keys = sorted({residue_key(atom) for atom, flag in zip(ref_prot_heavy, initial_core_nearest <= CUTOFF) if flag})
    if not initial_keys:
        raise RuntimeError(f"{name}: no initial core-pocket contacts")
    prot_keys = [residue_key(atom) for atom in prot_heavy]
    pocket_indices = np.asarray([i for i, key in enumerate(prot_keys) if key in initial_keys], dtype=int)
    pocket_key_to_indices = {key: np.asarray([i for i, x in enumerate(prot_keys) if x == key], dtype=int) for key in initial_keys}
    ref_ca = ref.select_atoms("protein and name CA")
    traj_ca = u.atoms[ref_ca.indices]
    pocket_ca_indices = np.asarray([i for i, atom in enumerate(ref_ca) if residue_key(atom) in initial_keys], dtype=int)
    ref_pocket_ca = ref_ca.positions.copy()[pocket_ca_indices]

    prot_donors = build_donor_pairs(ref_prot_all, ref_prot_all.positions.copy())
    lig_donors = build_donor_pairs(ref_lig_all, ref_lig_all.positions.copy())
    prot_acceptors = acceptor_indices(ref_prot_all)
    lig_acceptors = acceptor_indices(ref_lig_all)
    prot_hydro_mask = np.asarray([atom_element(a) in {"C", "S"} and a.resname in HYDRO_RESNAMES for a in prot_heavy], dtype=bool)
    lig_hydro_mask = np.asarray([atom_element(a) in {"C", "S"} for a in lig_heavy], dtype=bool)

    rows = []
    samples = []
    residue_counter = Counter()
    residue_counter_50 = Counter()
    residue_counter_70 = Counter()
    lig_sum = np.zeros((len(lig_heavy), 3), dtype=float)
    lig_sq = np.zeros_like(lig_sum)
    internal_sum = np.zeros_like(lig_sum)
    internal_sq = np.zeros_like(lig_sum)
    ca_sum = np.zeros_like(ref_pocket_ca)
    ca_sq = np.zeros_like(ref_pocket_ca)
    cluster_items = []
    previous_time = None
    previous_rmsd = None
    largest_step = {"delta_A": -1.0, "time_ns": None}

    for frame_index, ts in enumerate(u.trajectory):
        time_ns = float(ts.time / 1000.0)
        box6 = np.asarray(ts.dimensions, dtype=float)
        if not np.all(box6[:3] > 0):
            box6 = np.asarray([81.0, 81.0, 81.0, 90.0, 90.0, 90.0], dtype=float)
        box = box6[:3]
        bb_xyz = bb.positions.copy()
        bb_center = bb_xyz.mean(axis=0)
        rotation, fit_rmsd = align.rotation_matrix(bb_xyz - bb_center, ref_bb_xyz - ref_bb_center)
        lig_xyz = lig_all.positions.copy()
        lig_xyz += np.rint((bb_center - lig_xyz.mean(axis=0)) / box) * box
        aligned_all = apply_left_rotation(lig_xyz, rotation, bb_center, ref_bb_center)
        aligned_heavy = aligned_all[:len(lig_heavy)]
        aligned_core = aligned_heavy[cfg["core"]]
        aligned_ca_all = apply_left_rotation(traj_ca.positions.copy(), rotation, bb_center, ref_bb_center)
        aligned_ca = aligned_ca_all[pocket_ca_indices]
        whole_rmsd = rmsd(aligned_heavy, ref_lig_heavy_xyz)
        core_rmsd = rmsd(aligned_core, ref_core_xyz)
        whole_com = float(np.linalg.norm(aligned_heavy.mean(axis=0) - ref_lig_heavy_xyz.mean(axis=0)))
        core_com = float(np.linalg.norm(aligned_core.mean(axis=0) - ref_core_center))
        core_internal_rotation, _ = align.rotation_matrix(aligned_core - aligned_core.mean(axis=0), ref_core_xyz - ref_core_xyz.mean(axis=0))
        internal_heavy = apply_left_rotation(aligned_heavy, core_internal_rotation, aligned_core.mean(axis=0), ref_core_center)
        core_internal = internal_heavy[cfg["core"]]
        core_internal_rmsd = rmsd(core_internal, ref_core_xyz)
        current_core_dist = distance_array(prot_heavy.positions.copy(), lig_heavy.positions.copy()[cfg["core"]], box=box6)
        current_core_min = current_core_dist.min(axis=1)
        core_retention = float(np.mean([current_core_min[pocket_key_to_indices[key]].min() <= CUTOFF for key in initial_keys]))
        pocket_min = float(current_core_min[pocket_indices].min())
        pocket_ca_rmsd = rmsd(aligned_ca, ref_pocket_ca)
        row = {
            "time_ns": time_ns,
            "protein_backbone_rmsd_A": float(fit_rmsd),
            "whole_ligand_rmsd_A": whole_rmsd,
            "core_rmsd_A": core_rmsd,
            "branch_A_rmsd_A": rmsd(aligned_heavy[cfg["branch_A"]], ref_lig_heavy_xyz[cfg["branch_A"]]),
            "branch_B_rmsd_A": rmsd(aligned_heavy[cfg["branch_B"]], ref_lig_heavy_xyz[cfg["branch_B"]]),
            "whole_com_displacement_A": whole_com,
            "core_com_displacement_A": core_com,
            "core_internal_rmsd_A": core_internal_rmsd,
            "pocket_min_distance_A": pocket_min,
            "initial_core_contact_retention": core_retention,
            "pocket_ca_rmsd_A": pocket_ca_rmsd,
        }
        rows.append(row)
        lig_sum += aligned_heavy
        lig_sq += aligned_heavy * aligned_heavy
        internal_sum += internal_heavy
        internal_sq += internal_heavy * internal_heavy
        ca_sum += aligned_ca
        ca_sq += aligned_ca * aligned_ca
        if previous_rmsd is not None and abs(core_rmsd - previous_rmsd) > largest_step["delta_A"]:
            largest_step = {"delta_A": abs(core_rmsd - previous_rmsd), "time_ns": time_ns}
        previous_rmsd = core_rmsd
        previous_time = time_ns

        if frame_index % CONTACT_STRIDE == 0:
            contact_d = distance_array(prot_heavy.positions.copy(), lig_heavy.positions.copy(), box=box6)
            contact_mask = contact_d <= CUTOFF
            contact_indices = np.where(contact_mask.any(axis=1))[0]
            contact_keys = sorted({prot_keys[i] for i in contact_indices})
            for key in contact_keys:
                residue_counter[key] += 1
                if time_ns >= 50.0:
                    residue_counter_50[key] += 1
                if time_ns >= 70.0:
                    residue_counter_70[key] += 1
            hydro_d = contact_d[np.ix_(prot_hydro_mask, lig_hydro_mask)]
            hydro_mask = hydro_d <= CUTOFF
            hydro_prot_indices = np.where(prot_hydro_mask)[0][np.where(hydro_mask.any(axis=1))[0]]
            hydro_keys = sorted({prot_keys[i] for i in hydro_prot_indices})
            hb = hbond_proxy(prot_all.positions.copy(), lig_all.positions.copy(), box, prot_donors, prot_acceptors, lig_donors, lig_acceptors)
            samples.append({
                "time_ns": time_ns,
                "total_contact_residues": len(contact_keys),
                "total_contact_pairs": int(contact_mask.sum()),
                "hydrophobic_contact_residues": len(hydro_keys),
                "hydrophobic_contact_pairs": int(hydro_mask.sum()),
                "hbond_proxy": hb,
                "core_contact_retention": core_retention,
                "pocket_min_distance_A": pocket_min,
            })
            if time_ns >= 70.0:
                cluster_items.append((time_ns, aligned_core.copy()))

    write_csv(out / "metrics_0_100ns.csv", rows)
    write_csv(out / "contacts_5ns.csv", bin_contacts(samples))
    total_samples = len(samples)
    occ_rows = []
    all_keys = sorted(set(residue_counter) | set(initial_keys))
    for key in all_keys:
        occ_rows.append({
            "residue": key,
            "initial_core_contact": key in initial_keys,
            "occupancy_0_100": residue_counter.get(key, 0) / total_samples,
            "occupancy_50_100": residue_counter_50.get(key, 0) / max(1, sum(1 for x in samples if x["time_ns"] >= 50.0)),
            "occupancy_70_100": residue_counter_70.get(key, 0) / max(1, sum(1 for x in samples if x["time_ns"] >= 70.0)),
        })
    write_csv(out / "residue_contact_occupancy.csv", occ_rows)

    aligned_rmsf = np.sqrt(np.maximum(0.0, ((lig_sq / len(rows)) - (lig_sum / len(rows)) ** 2).sum(axis=1)))
    internal_rmsf = np.sqrt(np.maximum(0.0, ((internal_sq / len(rows)) - (internal_sum / len(rows)) ** 2).sum(axis=1)))
    write_csv(out / "ligand_rmsf_components.csv", [
        {"atom_index_1based": i + 1, "atom_name": atom.name, "protein_aligned_rmsf_A": float(aligned_rmsf[i]), "core_aligned_internal_rmsf_A": float(internal_rmsf[i]), "component": "core" if i in cfg["core"] else "branch_A" if i in cfg["branch_A"] else "branch_B"}
        for i, atom in enumerate(ref_lig_heavy)
    ])
    pocket_ca_rmsf = np.sqrt(np.maximum(0.0, ((ca_sq / len(rows)) - (ca_sum / len(rows)) ** 2).sum(axis=1)))
    write_csv(out / "pocket_ca_rmsf.csv", [{"residue": residue_key(atom), "rmsf_A": float(pocket_ca_rmsf[i])} for i, atom in enumerate(ref_ca[pocket_ca_indices])])
    write_csv(out / "last30ns_core_clusters.csv", cluster_last30(cluster_items))

    def window(lo, hi):
        subset = [r for r in rows if lo <= r["time_ns"] < hi]
        return {key: stats([r[key] for r in subset]) for key in (
            "whole_ligand_rmsd_A", "core_rmsd_A", "branch_A_rmsd_A", "branch_B_rmsd_A",
            "whole_com_displacement_A", "core_com_displacement_A", "core_internal_rmsd_A",
            "pocket_min_distance_A", "initial_core_contact_retention", "protein_backbone_rmsd_A", "pocket_ca_rmsd_A",
        )}

    outside = [r for r in rows if r["pocket_min_distance_A"] > CUTOFF]
    inside_fraction = 1.0 - len(outside) / len(rows)
    outside_episodes = []
    start = None
    prev = None
    dt = float(np.median(np.diff([r["time_ns"] for r in rows])))
    for r in rows:
        is_out = r["pocket_min_distance_A"] > CUTOFF
        if is_out and start is None:
            start = r["time_ns"]
        if not is_out and start is not None:
            outside_episodes.append({"start_ns": start, "end_ns": prev + dt, "duration_ns": prev + dt - start})
            start = None
        prev = r["time_ns"]
    if start is not None:
        outside_episodes.append({"start_ns": start, "end_ns": prev + dt, "duration_ns": prev + dt - start})
    clusters = cluster_last30(cluster_items)
    max50 = max((r for r in rows if r["time_ns"] >= 50.0), key=lambda r: r["core_com_displacement_A"])
    first_core6 = next((r for r in rows if r["core_rmsd_A"] >= 6.0), None)
    summary = {
        "system": name,
        "ligand": cfg["ligand_label"],
        "topology": str(cfg["topology"]),
        "trajectory": str(cfg["trajectory"]),
        "n_frames": len(rows),
        "first_time_ns": rows[0]["time_ns"],
        "last_time_ns": rows[-1]["time_ns"],
        "ligand_heavy_atoms": len(lig_heavy),
        "component_indices_0based": {k: cfg[k].tolist() for k in ("core", "branch_A", "branch_B")},
        "initial_core_contact_residues": initial_keys,
        "initial_core_contact_residue_count": len(initial_keys),
        "largest_core_rmsd_step": largest_step,
        "first_core_rmsd_ge_6_A": {"time_ns": first_core6["time_ns"], "value_A": first_core6["core_rmsd_A"]} if first_core6 else None,
        "largest_core_com_after_50ns": {"time_ns": max50["time_ns"], "value_A": max50["core_com_displacement_A"]},
        "window_0_10ns": window(0.0, 10.0),
        "window_10_50ns": window(10.0, 50.0),
        "window_50_80ns": window(50.0, 80.0),
        "window_80_100ns": window(80.0, 100.1),
        "pocket_residency": {
            "definition": "core minimum distance to initial core-contact residues <=4.5 A",
            "fraction_inside": inside_fraction,
            "fraction_outside": 1.0 - inside_fraction,
            "outside_total_ns": float(sum(x["duration_ns"] for x in outside_episodes)),
            "longest_outside_episode_ns": max((x["duration_ns"] for x in outside_episodes), default=0.0),
            "farthest_distance_A": float(max(r["pocket_min_distance_A"] for r in rows)),
            "farthest_time_ns": float(max(rows, key=lambda r: r["pocket_min_distance_A"])["time_ns"]),
            "outside_episodes": outside_episodes,
        },
        "last30ns_clustering": {
            "threshold_core_rmsd_A": 2.0,
            "cluster_count": len(clusters),
            "dominant_occupancy": clusters[0]["occupancy"] if clusters else None,
            "second_occupancy": clusters[1]["occupancy"] if len(clusters) > 1 else None,
        },
        "rmsf_summary": {
            "protein_aligned_by_component": {
                "core_mean_A": float(np.mean(aligned_rmsf[cfg["core"]])),
                "branch_A_mean_A": float(np.mean(aligned_rmsf[cfg["branch_A"]])),
                "branch_B_mean_A": float(np.mean(aligned_rmsf[cfg["branch_B"]])),
            },
            "core_aligned_internal_by_component": {
                "core_mean_A": float(np.mean(internal_rmsf[cfg["core"]])),
                "branch_A_mean_A": float(np.mean(internal_rmsf[cfg["branch_A"]])),
                "branch_B_mean_A": float(np.mean(internal_rmsf[cfg["branch_B"]])),
            },
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return {"summary": summary, "rows": rows, "contacts": bin_contacts(samples), "occupancy": occ_rows}


def write_report(results):
    lines = [
        "# DINP–PPARG vs MiNP–PPARG same-method 100 ns audit",
        "",
        "Both systems were reanalyzed from their production DCD with the same code and definitions. Protein-backbone fitting uses the corresponding initial system PDB; ligand RMSD/COM values are relative to that initial pose. Contact cutoffs are 4.5 Å with PBC-aware distances and contact summaries are sampled every 0.1 ns and binned every 5 ns.",
        "",
        "MiNP component definitions are disjoint and chemically explicit: core indices 12–17, long ester/alkyl branch A indices 0–11, and short carboxyl branch B indices 18–20. This replaces the earlier overlapping branch-B bookkeeping used in the preliminary MiNP QC.",
        "",
        "The earlier DINP report used the first production frame as its RMSD reference. This audit deliberately uses the corresponding initial system PDB for both ligands, so the values below are a docking-to-trajectory comparison on one baseline. The per-frame tables preserve the exact definitions and make the baseline explicit.",
        "",
        "## 80–100 ns comparison",
        "",
        "| metric | DINP | MiNP |",
        "|---|---:|---:|",
    ]
    keys = [
        ("whole_ligand_rmsd_A", "whole RMSD (Å)"), ("core_rmsd_A", "core RMSD (Å)"),
        ("branch_A_rmsd_A", "branch A RMSD (Å)"), ("branch_B_rmsd_A", "branch B RMSD (Å)"),
        ("whole_com_displacement_A", "whole COM (Å)"), ("core_com_displacement_A", "core COM (Å)"),
        ("pocket_min_distance_A", "pocket min distance (Å)"), ("initial_core_contact_retention", "initial core-contact retention"),
        ("protein_backbone_rmsd_A", "protein backbone RMSD (Å)"), ("pocket_ca_rmsd_A", "pocket Cα RMSD (Å)"),
    ]
    for key, label in keys:
        vals = []
        for name in ("DINP", "MiNP"):
            s = results[name]["summary"]["window_80_100ns"][key]
            vals.append(f"{s['mean']:.3f} ± {s['sd']:.3f}")
        lines.append(f"| {label} | {vals[0]} | {vals[1]} |")
    lines += ["", "## Contact summary", "", "| window | DINP total residues | MiNP total residues | DINP hydrophobic | MiNP hydrophobic | DINP H-bond proxy | MiNP H-bond proxy |", "|---|---:|---:|---:|---:|---:|---:|"]
    for window_name, lo, hi in (("0–10 ns", 0.0, 10.0), ("10–50 ns", 10.0, 50.0), ("50–80 ns", 50.0, 80.0), ("80–100 ns", 80.0, 100.1)):
        vals = []
        for name in ("DINP", "MiNP"):
            subset = [r for r in results[name]["contacts"] if lo <= r["bin_start_ns"] < hi]
            vals.append((np.mean([r["total_contact_residues_mean"] for r in subset]), np.mean([r["hydrophobic_contact_residues_mean"] for r in subset]), np.mean([r["hbond_proxy_mean"] for r in subset])))
        lines.append(f"| {window_name} | {vals[0][0]:.2f} | {vals[1][0]:.2f} | {vals[0][1]:.2f} | {vals[1][1]:.2f} | {vals[0][2]:.2f} | {vals[1][2]:.2f} |")
    lines += ["", "## Interpretation", ""]
    d = results["DINP"]["summary"]; m = results["MiNP"]["summary"]
    lines.append(f"- DINP initial core-pocket retention: {d['pocket_residency']['fraction_inside']:.3f} residency; farthest sampled core distance {d['pocket_residency']['farthest_distance_A']:.2f} Å at {d['pocket_residency']['farthest_time_ns']:.2f} ns.")
    lines.append(f"- MiNP initial core-pocket retention: {m['pocket_residency']['fraction_inside']:.3f} residency; farthest sampled core distance {m['pocket_residency']['farthest_distance_A']:.2f} Å at {m['pocket_residency']['farthest_time_ns']:.2f} ns.")
    lines.append(f"- DINP last-30-ns core clustering: {d['last30ns_clustering']['cluster_count']} clusters, dominant occupancy {d['last30ns_clustering']['dominant_occupancy']:.3f}.")
    lines.append(f"- MiNP last-30-ns core clustering: {m['last30ns_clustering']['cluster_count']} clusters, dominant occupancy {m['last30ns_clustering']['dominant_occupancy']:.3f}.")
    lines += ["", "The detailed audit tables are in `outputs/`; figures are in `figures/`. Geometric hydrogen bonds are a distance/angle proxy and do not replace an energy decomposition. The comparison evaluates structural persistence and rearrangement, not experimental affinity.", ""]
    (ROOT / "DINP_vs_MINP_PPARG_MD_AUDIT_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    warnings.filterwarnings("ignore")
    results = {name: analyze_system(name, cfg) for name, cfg in SYSTEMS.items()}
    combined = []
    metric_keys = [
        "whole_ligand_rmsd_A", "core_rmsd_A", "branch_A_rmsd_A", "branch_B_rmsd_A",
        "whole_com_displacement_A", "core_com_displacement_A", "pocket_min_distance_A",
        "initial_core_contact_retention", "protein_backbone_rmsd_A", "pocket_ca_rmsd_A",
    ]
    for name, result in results.items():
        for window_name, lo, hi in (("0_10ns", 0.0, 10.0), ("10_50ns", 10.0, 50.0), ("50_80ns", 50.0, 80.0), ("80_100ns", 80.0, 100.1)):
            subset = [r for r in result["rows"] if lo <= r["time_ns"] < hi]
            row = {"system": name, "window": window_name, "n_frames": len(subset)}
            for key in metric_keys:
                s = stats([r[key] for r in subset])
                row[f"{key}_mean"] = s["mean"]; row[f"{key}_sd"] = s["sd"]; row[f"{key}_min"] = s["min"]; row[f"{key}_max"] = s["max"]
            combined.append(row)
    write_csv(OUT / "comparison_audit_table.csv", combined)
    contact_combined = []
    for name, result in results.items():
        for row in result["contacts"]:
            contact_combined.append({"system": name, **row})
    write_csv(OUT / "contact_comparison_5ns.csv", contact_combined)
    contact_windows = []
    for name, result in results.items():
        for window_name, lo, hi in (("0_10ns", 0.0, 10.0), ("10_50ns", 10.0, 50.0), ("50_80ns", 50.0, 80.0), ("80_100ns", 80.0, 100.1)):
            subset = [r for r in result["contacts"] if lo <= r["bin_start_ns"] < hi]
            out = {"system": name, "window": window_name, "n_bins": len(subset)}
            for key in ("total_contact_residues_mean", "total_contact_pairs_mean", "hydrophobic_contact_residues_mean", "hydrophobic_contact_pairs_mean", "hbond_proxy_mean", "core_contact_retention_mean", "pocket_min_distance_A_mean"):
                out[key] = float(np.mean([r[key] for r in subset])) if subset else None
            contact_windows.append(out)
    write_csv(OUT / "contact_window_summary.csv", contact_windows)
    write_report(results)
    manifest = {
        "method": "same-method DINP vs MiNP PPARG 100 ns audit",
        "systems": {name: {"topology": str(cfg["topology"]), "trajectory": str(cfg["trajectory"])} for name, cfg in SYSTEMS.items()},
        "outputs": [str(p) for p in sorted(OUT.rglob("*")) if p.is_file()],
        "figures": [str(p) for p in sorted(FIG.glob("*")) if p.is_file()],
        "report": str(ROOT / "DINP_vs_MINP_PPARG_MD_AUDIT_REPORT.md"),
    }
    (ROOT / "audit_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({name: result["summary"] for name, result in results.items()}, indent=2))


if __name__ == "__main__":
    main()
