from pathlib import Path
import json, numpy as np, MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb")
DCD=Path(r"E:\chatgpt\qc_50ns\production.dcd")
u=mda.Universe(str(TOP),str(DCD)); lig=u.select_atoms('resname UNK DINP UNL LIG'); prot=u.select_atoms('protein and not name H*')
# choose closest frame to 50 ns
best=None
for ts in u.trajectory:
 d=abs(float(ts.time)-50000.0)
 if best is None or d<best[0]: best=(d,int(ts.frame))
u.trajectory[best[1]]; ts=u.trajectory.ts
near=distance_array(lig.positions,prot.positions,box=ts.dimensions)
mask=np.any(near<=4.5,axis=0)
out={'frame_index':int(ts.frame),'time_ps':float(ts.time),'protein_heavy_atoms_with_ligand_within_4p5A':int(mask.sum()),'ligand_atoms':int(lig.n_atoms),'protein_heavy_atoms':int(prot.n_atoms),'coordinate_nonfinite':bool(not np.isfinite(lig.positions).all() or not np.isfinite(prot.positions).all())}
Path(r"E:\chatgpt\qc_50ns\contacts_50ns.json").write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(out,ensure_ascii=False))
