from pathlib import Path
import itertools,json,numpy as np,MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb")
REFPDB=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb")
DCD=Path(r"E:\chatgpt\qc_50ns\production.dcd"); OUT=Path(r"E:\chatgpt\qc_50ns"); OUT.mkdir(parents=True,exist_ok=True)
def unwrap_group(pos,box):
 a=pos[0].copy(); d=pos-a; d-=box*np.round(d/box); return a+d
def kabsch(mobile,target):
 cm,ct=mobile.mean(0),target.mean(0); u,_,vt=np.linalg.svd((mobile-cm).T@(target-ct)); r=u@vt
 if np.linalg.det(r)<0: u[:,-1]*=-1; r=u@vt
 return r,ct-cm@r
u=mda.Universe(str(TOP),str(DCD)); ref=mda.Universe(str(REFPDB)); ca=u.select_atoms('protein and name CA'); lig=u.select_atoms('resname UNK DINP UNL LIG'); ref_ca=ref.select_atoms('protein and name CA').positions.copy().astype(float); ref_lig=ref.select_atoms('resname UNK DINP UNL LIG').positions.copy().astype(float); ref_com=ref_lig.mean(0); prot=u.select_atoms('protein and not name H*')
rows=[]; contacts=[]
for ts in u.trajectory:
 box=np.asarray(ts.dimensions[:3],float); ca_cur=ca.positions.copy().astype(float); lig_cur=unwrap_group(lig.positions.copy().astype(float),box); rot,tran=kabsch(ca_cur,ref_ca); ca_fit=ca_cur@rot+tran; best=None
 for n in itertools.product([-1,0,1],repeat=3):
  fit=(lig_cur+np.asarray(n,float)*box)@rot+tran; score=float(np.linalg.norm(fit.mean(0)-ref_com))
  if best is None or score<best[0]: best=(score,fit)
 lf=best[1]; rows.append((float(ts.time),float(np.sqrt(np.mean(np.sum((ca_fit-ref_ca)**2,axis=1)))),float(np.sqrt(np.mean(np.sum((lf-ref_lig)**2,axis=1)))),float(np.linalg.norm(lf.mean(0)-ref_com)))); contacts.append(int(np.any(distance_array(lig_cur,prot.positions,box=ts.dimensions)<=4.5,axis=0).sum()))
a=np.asarray(rows); c=np.asarray(contacts); out={'status':'ok','n_frames':int(len(a)),'time_start_ps':float(a[0,0]),'time_end_ps':float(a[-1,0]),'protein_ca_rmsd_vs_minimized_mean_A':float(a[:,1].mean()),'protein_ca_rmsd_vs_minimized_last_A':float(a[-1,1]),'protein_ca_rmsd_vs_minimized_max_A':float(a[:,1].max()),'ligand_rmsd_pbc_corrected_vs_minimized_mean_A':float(a[:,2].mean()),'ligand_rmsd_pbc_corrected_vs_minimized_last_A':float(a[-1,2]),'ligand_rmsd_pbc_corrected_vs_minimized_max_A':float(a[:,2].max()),'ligand_com_displacement_pbc_corrected_vs_minimized_mean_A':float(a[:,3].mean()),'ligand_com_displacement_pbc_corrected_vs_minimized_last_A':float(a[-1,3]),'ligand_com_displacement_pbc_corrected_vs_minimized_max_A':float(a[:,3].max()),'ligand_protein_contacts_4p5A_mean':float(c.mean()),'ligand_protein_contacts_4p5A_last':int(c[-1]),'ligand_protein_contacts_4p5A_min':int(c.min()),'ligand_protein_contacts_4p5A_max':int(c.max()),'coordinate_nonfinite':bool(not np.isfinite(a).all())}; np.savetxt(OUT/'pbc_corrected_vs_minimized.csv',a,delimiter=',',header='time_ps,protein_ca_rmsd_A,ligand_rmsd_A,ligand_com_displacement_A',comments=''); (OUT/'pbc_corrected_vs_minimized.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(out,ensure_ascii=False))
