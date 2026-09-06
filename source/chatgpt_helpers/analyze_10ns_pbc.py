from pathlib import Path
import json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

TOP = Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb")
DCD = Path(r"E:\chatgpt\qc_10ns\production_10ns.dcd")
OUT = Path(r"E:\chatgpt\qc_10ns")


def unwrap_group(pos, box):
    anchor = pos[0].copy()
    delta = pos - anchor
    delta -= box * np.round(delta / box)
    return anchor + delta


def kabsch(mobile, target):
    cm = mobile.mean(axis=0)
    ct = target.mean(axis=0)
    x = mobile - cm
    y = target - ct
    u, _, vt = np.linalg.svd(x.T @ y)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1.0
        rot = u @ vt
    tran = ct - cm @ rot
    return rot, tran


u = mda.Universe(str(TOP), str(DCD))
ca = u.select_atoms("protein and name CA")
lig = u.select_atoms("resname UNK DINP UNL LIG")
prot_heavy = u.select_atoms("protein and not name H*")
if not ca.n_atoms or not lig.n_atoms:
    raise RuntimeError("missing protein CA or ligand")

u.trajectory[0]
box0 = np.asarray(u.trajectory.ts.dimensions[:3], dtype=float)
ca_ref = ca.positions.copy().astype(float)
lig_ref = unwrap_group(lig.positions.copy().astype(float), box0)
ca_ref_com = ca_ref.mean(axis=0)
lig_ref_com = lig_ref.mean(axis=0)
ref_rel = lig_ref_com - ca_ref_com

rows = []
contact_counts = []
for ts in u.trajectory:
    box = np.asarray(ts.dimensions[:3], dtype=float)
    ca_cur = ca.positions.copy().astype(float)
    lig_cur = unwrap_group(lig.positions.copy().astype(float), box)
    # Choose the ligand image closest to the reference ligand/protein relative vector.
    rel = lig_cur.mean(axis=0) - ca_cur.mean(axis=0)
    rel -= box * np.round((rel - ref_rel) / box)
    lig_cur += (ca_cur.mean(axis=0) + rel) - lig_cur.mean(axis=0)

    rot, tran = kabsch(ca_cur, ca_ref)
    ca_fit = ca_cur @ rot + tran
    lig_fit = lig_cur @ rot + tran
    ca_rmsd = float(np.sqrt(np.mean(np.sum((ca_fit - ca_ref) ** 2, axis=1))))
    lig_rmsd = float(np.sqrt(np.mean(np.sum((lig_fit - lig_ref) ** 2, axis=1))))
    lig_com_disp = float(np.linalg.norm(lig_fit.mean(axis=0) - lig_ref_com))
    near = distance_array(lig_cur, prot_heavy.positions, box=ts.dimensions)
    contacts = int(np.any(near <= 4.5, axis=0).sum())
    rows.append((float(ts.time), ca_rmsd, lig_rmsd, lig_com_disp))
    contact_counts.append(contacts)

arr = np.asarray(rows)
result = {
    "status": "ok",
    "n_frames": int(len(arr)),
    "time_start_ps": float(arr[0, 0]),
    "time_end_ps": float(arr[-1, 0]),
    "protein_ca_rmsd_vs_first_frame_mean_A": float(arr[:, 1].mean()),
    "protein_ca_rmsd_vs_first_frame_last_A": float(arr[-1, 1]),
    "protein_ca_rmsd_vs_first_frame_max_A": float(arr[:, 1].max()),
    "ligand_rmsd_pbc_corrected_vs_first_frame_mean_A": float(arr[:, 2].mean()),
    "ligand_rmsd_pbc_corrected_vs_first_frame_last_A": float(arr[-1, 2]),
    "ligand_rmsd_pbc_corrected_vs_first_frame_max_A": float(arr[:, 2].max()),
    "ligand_com_displacement_pbc_corrected_mean_A": float(arr[:, 3].mean()),
    "ligand_com_displacement_pbc_corrected_last_A": float(arr[-1, 3]),
    "ligand_com_displacement_pbc_corrected_max_A": float(arr[:, 3].max()),
    "ligand_protein_contacts_4p5A_mean": float(np.mean(contact_counts)),
    "ligand_protein_contacts_4p5A_last": int(contact_counts[-1]),
    "ligand_protein_contacts_4p5A_min": int(np.min(contact_counts)),
    "ligand_protein_contacts_4p5A_max": int(np.max(contact_counts)),
    "coordinate_nonfinite": bool(not np.isfinite(arr).all()),
}
np.savetxt(OUT / "pbc_corrected_rmsd.csv", arr, delimiter=",", header="time_ps,protein_ca_rmsd_A,ligand_rmsd_A,ligand_com_displacement_A", comments="")
(OUT / "pbc_corrected_qc.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(result, ensure_ascii=False))
