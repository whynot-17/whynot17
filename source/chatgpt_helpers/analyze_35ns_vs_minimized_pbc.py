from pathlib import Path
import itertools
import json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

TOP = Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb")
REFPDB = Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb")
DCD = Path(r"E:\chatgpt\qc_35ns\production_35ns.dcd")
OUT = Path(r"E:\chatgpt\qc_35ns")


def unwrap_group(pos, box):
    anchor = pos[0].copy()
    delta = pos - anchor
    delta -= box * np.round(delta / box)
    return anchor + delta


def kabsch(mobile, target):
    cm = mobile.mean(axis=0)
    ct = target.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - cm).T @ (target - ct))
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1.0
        rot = u @ vt
    return rot, ct - cm @ rot


u = mda.Universe(str(TOP), str(DCD))
ref = mda.Universe(str(REFPDB))
ca = u.select_atoms("protein and name CA")
lig = u.select_atoms("resname UNK DINP UNL LIG")
ref_ca = ref.select_atoms("protein and name CA").positions.copy().astype(float)
ref_lig = ref.select_atoms("resname UNK DINP UNL LIG").positions.copy().astype(float)
ref_lig_com = ref_lig.mean(axis=0)
prot_heavy = u.select_atoms("protein and not name H*")
rows = []
contacts = []
for ts in u.trajectory:
    box = np.asarray(ts.dimensions[:3], dtype=float)
    ca_cur = ca.positions.copy().astype(float)
    lig_cur = unwrap_group(lig.positions.copy().astype(float), box)
    rot, tran = kabsch(ca_cur, ref_ca)
    ca_fit = ca_cur @ rot + tran
    # Pick the periodic image closest to the minimized reference ligand.
    best = None
    for n in itertools.product([-1, 0, 1], repeat=3):
        shifted = lig_cur + np.asarray(n, dtype=float) * box
        fit = shifted @ rot + tran
        score = float(np.linalg.norm(fit.mean(axis=0) - ref_lig_com))
        if best is None or score < best[0]:
            best = (score, fit)
    lig_fit = best[1]
    rows.append((float(ts.time),
                 float(np.sqrt(np.mean(np.sum((ca_fit - ref_ca) ** 2, axis=1)))),
                 float(np.sqrt(np.mean(np.sum((lig_fit - ref_lig) ** 2, axis=1)))),
                 float(np.linalg.norm(lig_fit.mean(axis=0) - ref_lig_com))))
    near = distance_array(lig_cur, prot_heavy.positions, box=ts.dimensions)
    contacts.append(int(np.any(near <= 4.5, axis=0).sum()))

arr = np.asarray(rows)
result = {
    "status": "ok",
    "n_frames": int(len(arr)),
    "time_start_ps": float(arr[0, 0]),
    "time_end_ps": float(arr[-1, 0]),
    "protein_ca_rmsd_vs_minimized_mean_A": float(arr[:, 1].mean()),
    "protein_ca_rmsd_vs_minimized_last_A": float(arr[-1, 1]),
    "protein_ca_rmsd_vs_minimized_max_A": float(arr[:, 1].max()),
    "ligand_rmsd_pbc_corrected_vs_minimized_mean_A": float(arr[:, 2].mean()),
    "ligand_rmsd_pbc_corrected_vs_minimized_last_A": float(arr[-1, 2]),
    "ligand_rmsd_pbc_corrected_vs_minimized_max_A": float(arr[:, 2].max()),
    "ligand_com_displacement_pbc_corrected_vs_minimized_mean_A": float(arr[:, 3].mean()),
    "ligand_com_displacement_pbc_corrected_vs_minimized_last_A": float(arr[-1, 3]),
    "ligand_com_displacement_pbc_corrected_vs_minimized_max_A": float(arr[:, 3].max()),
    "ligand_protein_contacts_4p5A_mean": float(np.mean(contacts)),
    "ligand_protein_contacts_4p5A_last": int(contacts[-1]),
    "ligand_protein_contacts_4p5A_min": int(np.min(contacts)),
    "ligand_protein_contacts_4p5A_max": int(np.max(contacts)),
    "coordinate_nonfinite": bool(not np.isfinite(arr).all()),
}
np.savetxt(OUT / "pbc_corrected_vs_minimized.csv", arr, delimiter=",", header="time_ps,protein_ca_rmsd_A,ligand_rmsd_A,ligand_com_displacement_A", comments="")
(OUT / "pbc_corrected_vs_minimized.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(result, ensure_ascii=False))

