import MDAnalysis as mda, numpy as np
from pathlib import Path
TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb"); REF=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb")
u=mda.Universe(str(TOP)); r=mda.Universe(str(REF)); a=u.select_atoms('resname UNK'); b=r.select_atoms('resname UNK'); print(a.n_atoms,b.n_atoms); print('core raw rmsd',np.sqrt(np.mean((a.positions[:6]-b.positions[:6])**2))); print('all raw',np.sqrt(np.mean((a.positions-b.positions)**2))); print('core com',a.positions[:6].mean(0),b.positions[:6].mean(0)); print('prot ca rmsd raw',np.sqrt(np.mean((u.select_atoms('protein and name CA').positions-r.select_atoms('protein and name CA').positions)**2)))
