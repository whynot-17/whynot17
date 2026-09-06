import MDAnalysis as mda, numpy as np
from MDAnalysis.lib.distances import distance_array
from pathlib import Path
TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb")
u=mda.Universe(str(TOP)); core=u.select_atoms('resname UNK and name C1x C2x C3x C4x C5x C6x'); prot=u.select_atoms('protein and not name H*'); d=distance_array(core.positions,prot.positions); mask=np.any(d<=4.5,axis=0); pairs=np.argwhere(d<=4.5); print(core.n_atoms,prot.n_atoms,'pocket_atoms',mask.sum(),'pairs',len(pairs),'min',d[:,mask].min());
for a in prot[mask]: print(a.resid,a.resname,a.name,round(float(d[:,a.index-prot.indices[0]].min()),3))
