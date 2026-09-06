from pathlib import Path
import itertools,json,numpy as np,MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb"); REF=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb"); DCD=Path(r"E:\chatgpt\qc_50ns\production.dcd"); OUT=Path(r"E:\chatgpt\qc_dinp_minimal_0_50ns"); OUT.mkdir(parents=True,exist_ok=True)
def k(a,b):
 ca,cb=a.mean(0),b.mean(0);u,_,vt=np.linalg.svd((a-ca).T@(b-cb));r=u@vt
 if np.linalg.det(r)<0:u[:,-1]*=-1;r=u@vt
 return r,cb-ca@r
def unwrap(x,box):a=x[0].copy();d=x-a;d-=box*np.round(d/box);return a+d
def rms(a,b):return float(np.sqrt(np.mean(np.sum((a-b)**2,axis=1))))
u=mda.Universe(str(TOP),str(DCD));ref=mda.Universe(str(REF)); lig=u.select_atoms('resname UNK'); rlig=ref.select_atoms('resname UNK'); names=list(lig.names); ci=np.array([names.index(n) for n in ['C1x','C2x','C3x','C4x','C5x','C6x']]); rnames=list(rlig.names); rci=np.array([rnames.index(n) for n in ['C1x','C2x','C3x','C4x','C5x','C6x']]); rcore=rlig.positions[rci].copy(); ca=u.select_atoms('protein and name CA'); rca=ref.select_atoms('protein and name CA').positions.copy(); prot=u.select_atoms('protein and not name H*'); rprot_ag=ref.select_atoms('protein and not name H*'); rprot=rprot_ag.positions.copy(); dref=distance_array(rcore,rprot); pocket=np.any(dref<=4.5,axis=0); pair=np.argwhere(dref<=4.5); pocket_res=[]
for j,a in enumerate(rprot_ag):
 if pocket[j]:
  key=(int(a.resid),a.resname)
  if key not in pocket_res:pocket_res.append(key)
res_groups=[]
for key in pocket_res: res_groups.append(np.array([j for j,a in enumerate(rprot_ag) if (int(a.resid),a.resname)==key],int))
rows=[]
for ts in u.trajectory:
 if ts.time>50000.1: break
 box=np.asarray(ts.dimensions[:3],float); raw=lig.positions.copy().astype(float); lu=unwrap(raw,box); core=lu[ci]; rc,tr=k(core,rcore); internal=rms(core@rc+tr,rcore)
 dcur=distance_array(raw[ci],prot.positions,box=ts.dimensions); pairret=float(np.mean(dcur[pair[:,0],pair[:,1]]<=4.5)); resret=float(np.mean([np.any(dcur[:,g]<=4.5) for g in res_groups])); nearest=float(dcur[:,pocket].min())
 rows.append((float(ts.time),internal,pairret,resret,nearest))
a=np.asarray(rows); np.savetxt(OUT/'core_internal_and_pocket_residue_timeseries.csv',a,delimiter=',',header='time_ps,core_internal_rmsd_A,pair_retention_4p5A,residue_contact_fraction_4p5A,core_pocket_nearest_distance_A',comments='')
def summ(x):return {'n_frames':int(len(x)),'core_internal_rmsd_mean_A':float(x[:,1].mean()),'core_internal_rmsd_sd_A':float(x[:,1].std(ddof=1)),'core_internal_rmsd_min_A':float(x[:,1].min()),'core_internal_rmsd_max_A':float(x[:,1].max()),'pair_retention_mean':float(x[:,2].mean()),'pair_retention_sd':float(x[:,2].std(ddof=1)),'pair_retention_min':float(x[:,2].min()),'pair_retention_max':float(x[:,2].max()),'residue_contact_fraction_mean':float(x[:,3].mean()),'residue_contact_fraction_sd':float(x[:,3].std(ddof=1)),'residue_contact_fraction_min':float(x[:,3].min()),'residue_contact_fraction_max':float(x[:,3].max()),'nearest_distance_mean_A':float(x[:,4].mean()),'nearest_distance_min_A':float(x[:,4].min()),'nearest_distance_max_A':float(x[:,4].max()),'nearest_le_4p5_fraction':float(np.mean(x[:,4]<=4.5))}
out={'reference_pocket_residues':pocket_res,'reference_pocket_atoms':int(pocket.sum()),'reference_core_pocket_pairs':int(len(pair)),'overall_0_50ns':summ(a),'first_half_0_25ns':summ(a[a[:,0]<=25000.1]),'second_half_25_50ns':summ(a[a[:,0]>=25000.0]),'segment_40_50ns':summ(a[(a[:,0]>=40000)&(a[:,0]<=50000.1)])}
for ns in (10,25,40,50):
 i=int(np.argmin(abs(a[:,0]-ns*1000))); out[f'nearest_{ns}ns']={'time_ps':float(a[i,0]),'core_internal_rmsd_A':float(a[i,1]),'pair_retention':float(a[i,2]),'residue_contact_fraction':float(a[i,3]),'nearest_distance_A':float(a[i,4])}
(OUT/'core_internal_and_pocket_residue_summary.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(out,ensure_ascii=False))


