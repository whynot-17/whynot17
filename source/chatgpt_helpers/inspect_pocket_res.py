import MDAnalysis as mda,numpy as np
from pathlib import Path
from MDAnalysis.lib.distances import distance_array
p=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb");u=mda.Universe(str(p)); core=u.select_atoms('resname UNK and name C1x C2x C3x C4x C5x C6x'); prot=u.select_atoms('protein and not name H*'); d=distance_array(core.positions,prot.positions); mask=np.any(d<=4.5,axis=0); groups={}
for j,a in enumerate(prot):
 if mask[j]: groups.setdefault((int(a.resid),a.resname),[]).append(j)
print('residues',len(groups)); print([(k,len(v)) for k,v in groups.items()]); print('pairs',int((d<=4.5).sum()))
