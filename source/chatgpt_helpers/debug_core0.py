import numpy as np,MDAnalysis as mda,itertools
from pathlib import Path
TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb"); REF=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb"); DCD=Path(r"E:\chatgpt\qc_50ns\production.dcd")
def k(a,b):
 ca,cb=a.mean(0),b.mean(0);u,_,vt=np.linalg.svd((a-ca).T@(b-cb));r=u@vt
 if np.linalg.det(r)<0:u[:,-1]*=-1;r=u@vt
 return r,cb-ca@r
def unwrap(x,box):a=x[0].copy();d=x-a;d-=box*np.round(d/box);return a+d
def rms(a,b):return np.sqrt(np.mean(np.sum((a-b)**2,axis=1)))
u=mda.Universe(str(TOP),str(DCD));r=mda.Universe(str(REF));lig=u.select_atoms('resname UNK');rl=r.select_atoms('resname UNK');ca=u.select_atoms('protein and name CA');rca=r.select_atoms('protein and name CA').positions.copy();names=list(lig.names);ci=np.array([names.index(n) for n in ['C1x','C2x','C3x','C4x','C5x','C6x']]);rc=rl.positions[:6].copy();rr=rl.positions.copy()
for ts in u.trajectory:
 if ts.frame in (0,1,2,10,100,999):
  box=ts.dimensions[:3];lu=unwrap(lig.positions.copy(),box);rp,tp=k(ca.positions.copy(),rca);best=None
  for sh in itertools.product([-1,0,1],repeat=3):
   img=lu+np.array(sh)*box;fit=img@rp+tp;score=np.linalg.norm(fit.mean(0)-rr.mean(0))
   if best is None or score<best[0]:best=(score,img,fit)
  print('fr',ts.frame,'t',ts.time,'box',box,'core raw com',lu[ci].mean(0),'core fit com',best[2][ci].mean(0),'ref core',rc.mean(0),'core rms',rms(best[2][ci],rc),'all rms',rms(best[2],rr),'shift score all',best[0])
