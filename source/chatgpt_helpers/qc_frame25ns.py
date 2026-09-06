from pathlib import Path
import itertools, json
import numpy as np
import MDAnalysis as mda

TOP = Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb")
REF = Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb")
DCD = Path(r"E:\chatgpt\qc_25ns\production_25ns.dcd")

def kabsch(a, b):
    ca, cb = a.mean(0), b.mean(0)
    u, _, vt = np.linalg.svd((a-ca).T @ (b-cb))
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    return r, cb - ca @ r

def unwrap(x, box):
    a = x[0].copy()
    return a + (x-a) - box*np.round((x-a)/box)

u = mda.Universe(str(TOP), str(DCD))
ref = mda.Universe(str(REF))
ca = u.select_atoms('protein and name CA')
lig = u.select_atoms('resname UNK DINP UNL LIG')
rca = ref.select_atoms('protein and name CA').positions.copy().astype(float)
rlig = ref.select_atoms('resname UNK DINP UNL LIG').positions.copy().astype(float)
u.trajectory[2499]
ts = u.trajectory.ts
box = ts.dimensions[:3].astype(float)
xca = ca.positions.copy().astype(float)
xlig = unwrap(lig.positions.copy().astype(float), box)
r, t = kabsch(xca, rca)
best = None
for n in itertools.product([-1,0,1], repeat=3):
    fit = (xlig + np.asarray(n)*box) @ r + t
    score = float(np.linalg.norm(fit.mean(0) - rlig.mean(0)))
    if best is None or score < best[0]: best = (score, fit)
lf = best[1]
out = {
    'frame_index': int(ts.frame),
    'time_ps': float(ts.time),
    'protein_ca_rmsd_vs_minimized_A': float(np.sqrt(np.mean(np.sum(((xca@r+t)-rca)**2,axis=1)))),
    'ligand_rmsd_pbc_corrected_vs_minimized_A': float(np.sqrt(np.mean(np.sum((lf-rlig)**2,axis=1)))),
    'ligand_com_displacement_pbc_corrected_vs_minimized_A': float(np.linalg.norm(lf.mean(0)-rlig.mean(0))),
    'ligand_com': lf.mean(0).tolist(),
    'coordinate_nonfinite': bool(not np.isfinite(xca).all() or not np.isfinite(xlig).all()),
}
print(json.dumps(out, ensure_ascii=False))
