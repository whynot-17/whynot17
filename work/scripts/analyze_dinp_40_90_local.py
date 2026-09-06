from pathlib import Path
import itertools,json,numpy as np,MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb"); REF=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb"); DCD=Path(r"E:\chatgpt\qc_90ns\production.dcd"); OUT=Path(r"E:\chatgpt\qc_dinp_40_90ns"); OUT.mkdir(parents=True,exist_ok=True)
CUTOFF=4.5
def k(a,b):
 ca,cb=a.mean(0),b.mean(0);u,_,vt=np.linalg.svd((a-ca).T@(b-cb));r=u@vt
 if np.linalg.det(r)<0:u[:,-1]*=-1;r=u@vt
 return r,cb-ca@r
def unwrap(x,box):a=x[0].copy();d=x-a;d-=box*np.round(d/box);return a+d
def rms(a,b):return float(np.sqrt(np.mean(np.sum((a-b)**2,axis=1))))
# Load trajectory and topology.
u=mda.Universe(str(TOP),str(DCD)); ref=mda.Universe(str(REF)); lig=u.select_atoms('resname UNK'); rlig=ref.select_atoms('resname UNK'); names=list(lig.names); rn=list(rlig.names)
core_i=np.array([names.index(n) for n in ['C1x','C2x','C3x','C4x','C5x','C6x']]); b1_i=np.array([names.index(n) for n in ['C7x','O1x','O2x','C8x','C9x','C10x','C11x','C12x','C13x','C14x','C15x','C16x']]); b2_i=np.array([names.index(n) for n in ['C17x','O3x','O4x','C18x','C19x','C20x','C21x','C22x','C23x','C24x','C25x','C26x']]); rcore_min=rlig.positions[[rn.index(n) for n in ['C1x','C2x','C3x','C4x','C5x','C6x']]].copy(); rb1_min=rlig.positions[[rn.index(n) for n in ['C7x','O1x','O2x','C8x','C9x','C10x','C11x','C12x','C13x','C14x','C15x','C16x']]].copy(); rb2_min=rlig.positions[[rn.index(n) for n in ['C17x','O3x','O4x','C18x','C19x','C20x','C21x','C22x','C23x','C24x','C25x','C26x']]].copy(); rall_min=rlig.positions.copy(); rcom_min=rall_min.mean(0)
ca=u.select_atoms('protein and name CA'); prot=u.select_atoms('protein and not name H*'); rca_min=ref.select_atoms('protein and name CA').positions.copy(); rprot_min=ref.select_atoms('protein and not name H*').positions.copy()
# Find 40 ns reference frame and cache its local coordinates.
ref40=None
for ts in u.trajectory:
 if ref40 is None or abs(float(ts.time)-40000.0)<ref40[0]:
  raw=lig.positions.copy().astype(float); box=np.asarray(ts.dimensions[:3],float); lu=unwrap(raw,box); ref40=(abs(float(ts.time)-40000.0),int(ts.frame),float(ts.time),np.asarray(ts.dimensions,float).copy(),lu.copy(),ca.positions.copy().astype(float),prot.positions.copy().astype(float))
if ref40 is None: raise RuntimeError('no 40 ns frame')
_,f40,t40,box40,lu40,ca40,prot40=ref40; core40=lu40[core_i]; b140=lu40[b1_i]; b240=lu40[b2_i]; all40=lu40; com40=all40.mean(0)
# Define pocket from the 40 ns local state and from minimized reference for independent contact views.
d40=distance_array(core40,prot40,box=box40); pair40=np.argwhere(d40<=CUTOFF); pocket40=np.any(d40<=CUTOFF,axis=0); rprot40_ag=ref.select_atoms('protein and not name H*'); dmin=distance_array(rcore_min,rprot_min); pairmin=np.argwhere(dmin<=CUTOFF); pocketmin=np.any(dmin<=CUTOFF,axis=0); pocket40_keys=[]; pocketmin_keys=[]
for j,a in enumerate(rprot40_ag):
 if pocket40[j] and (int(a.resid),a.resname) not in pocket40_keys:pocket40_keys.append((int(a.resid),a.resname))
for j,a in enumerate(rprot40_ag):
 if pocketmin[j] and (int(a.resid),a.resname) not in pocketmin_keys:pocketmin_keys.append((int(a.resid),a.resname))
def groups_for(keys): return [np.array([j for j,a in enumerate(rprot40_ag) if (int(a.resid),a.resname)==key],int) for key in keys]
g40=groups_for(pocket40_keys); gmin=groups_for(pocketmin_keys)
rows=[]
for ts in u.trajectory:
 t=float(ts.time)
 if t<39999.9: continue
 if t>90000.1: break
 box=np.asarray(ts.dimensions[:3],float); raw=lig.positions.copy().astype(float); lu=unwrap(raw,box); rotp,tranp=k(ca.positions.copy().astype(float),ca40); best=None
 for sh in itertools.product([-1,0,1],repeat=3):
  image=lu+np.asarray(sh,dtype=float)*box; fit=image@rotp+tranp; score=float(np.linalg.norm(fit[core_i].mean(0)-core40.mean(0)))
  if best is None or score<best[0]: best=(score,image,fit)
 image,fit=best[1],best[2]; core=image[core_i]; b1=image[b1_i]; b2=image[b2_i]; core_pos=rms(fit[core_i],core40); rc,tr=k(core,core40); core_int=rms(core@rc+tr,core40); b1r=rms(b1@rc+tr,b140); b2r=rms(b2@rc+tr,b240); com40d=float(np.linalg.norm(fit.mean(0)-com40)); rotmin,tranmin=k(ca.positions.copy().astype(float),rca_min); bestmin=None
 for sh in itertools.product([-1,0,1],repeat=3):
  image_m=lu+np.asarray(sh,dtype=float)*box; fit_m=image_m@rotmin+tranmin; score_m=float(np.linalg.norm(fit_m.mean(0)-rcom_min))
  if bestmin is None or score_m<bestmin[0]: bestmin=(score_m,fit_m)
 commind=float(bestmin[0])
 dcur=distance_array(raw[core_i],prot.positions,box=ts.dimensions); pairret40=float(np.mean(dcur[pair40[:,0],pair40[:,1]]<=CUTOFF)) if len(pair40) else float('nan'); resret40=float(np.mean([np.any(dcur[:,g]<=CUTOFF) for g in g40])) if g40 else float('nan'); near40=float(dcur[:,pocket40].min()) if pocket40.any() else float('nan'); pairretmin=float(np.mean(dcur[pairmin[:,0],pairmin[:,1]]<=CUTOFF)) if len(pairmin) else float('nan'); resretmin=float(np.mean([np.any(dcur[:,g]<=CUTOFF) for g in gmin])) if gmin else float('nan'); nearmin=float(dcur[:,pocketmin].min()) if pocketmin.any() else float('nan')
 rows.append((t,core_pos,core_int,b1r,b2r,com40d,commind,pairret40,resret40,near40,pairretmin,resretmin,nearmin))
a=np.asarray(rows); np.savetxt(OUT/'dinp_40_90_local_timeseries.csv',a,delimiter=',',header='time_ps,core_rmsd_vs_40ns_after_protein_fit_A,core_internal_rmsd_vs_40ns_A,branch1_rmsd_vs_40ns_after_core_fit_A,branch2_rmsd_vs_40ns_after_core_fit_A,com_drift_vs_40ns_A,com_displacement_vs_minimized_A,pair_retention_vs_40ns,pocket_residue_fraction_vs_40ns,nearest_distance_to_40ns_pocket_A,pair_retention_vs_minimized,pocket_residue_fraction_vs_minimized,nearest_distance_to_minimized_pocket_A',comments='')
def summ(x):
 return {'n_frames':int(len(x)),'time_start_ns':float(x[0,0]/1000),'time_end_ns':float(x[-1,0]/1000),'core_rmsd_vs_40ns_mean_A':float(x[:,1].mean()),'core_rmsd_vs_40ns_sd_A':float(x[:,1].std(ddof=1)),'core_rmsd_vs_40ns_min_A':float(x[:,1].min()),'core_rmsd_vs_40ns_max_A':float(x[:,1].max()),'core_internal_rmsd_mean_A':float(x[:,2].mean()),'core_internal_rmsd_max_A':float(x[:,2].max()),'branch1_rmsd_vs_40ns_mean_A':float(x[:,3].mean()),'branch1_rmsd_vs_40ns_sd_A':float(x[:,3].std(ddof=1)),'branch1_rmsd_vs_40ns_max_A':float(x[:,3].max()),'branch2_rmsd_vs_40ns_mean_A':float(x[:,4].mean()),'branch2_rmsd_vs_40ns_sd_A':float(x[:,4].std(ddof=1)),'branch2_rmsd_vs_40ns_max_A':float(x[:,4].max()),'com_drift_vs_40ns_mean_A':float(x[:,5].mean()),'com_drift_vs_40ns_sd_A':float(x[:,5].std(ddof=1)),'com_drift_vs_40ns_max_A':float(x[:,5].max()),'com_displacement_vs_minimized_mean_A':float(x[:,6].mean()),'pair_retention_vs_40ns_mean':float(x[:,7].mean()),'pair_retention_vs_40ns_min':float(x[:,7].min()),'pair_retention_vs_40ns_max':float(x[:,7].max()),'pocket_residue_fraction_vs_40ns_mean':float(x[:,8].mean()),'pocket_residue_fraction_vs_40ns_min':float(x[:,8].min()),'pocket_residue_fraction_vs_40ns_max':float(x[:,8].max()),'nearest_distance_to_40ns_pocket_mean_A':float(x[:,9].mean()),'nearest_distance_to_40ns_pocket_min_A':float(x[:,9].min()),'nearest_distance_to_40ns_pocket_max_A':float(x[:,9].max()),'nearest_to_40ns_pocket_le_4p5_fraction':float(np.mean(x[:,9]<=4.5)),'pair_retention_vs_minimized_mean':float(x[:,10].mean()),'pocket_residue_fraction_vs_minimized_mean':float(x[:,11].mean()),'nearest_distance_to_minimized_pocket_mean_A':float(x[:,12].mean()),'nearest_to_minimized_pocket_le_4p5_fraction':float(np.mean(x[:,12]<=4.5))}
res={'reference_40ns_frame_index':int(f40),'reference_40ns_time_ps':float(t40),'reference_40ns_pocket_atoms':int(pocket40.sum()),'reference_40ns_core_pocket_pairs':int(len(pair40)),'reference_40ns_pocket_residues':pocket40_keys,'minimized_reference_pocket_atoms':int(pocketmin.sum()),'minimized_reference_core_pocket_pairs':int(len(pairmin)),'overall_40_90ns':summ(a),'first_half_40_55ns':summ(a[a[:,0]<=55000.1]),'second_half_55_90ns':summ(a[a[:,0]>=55000.0])}
for name,col in [('core_rmsd_vs_40ns',1),('branch1_rmsd_vs_40ns',3),('branch2_rmsd_vs_40ns',4),('com_drift_vs_40ns',5),('pair_retention_vs_40ns',7),('nearest_distance_to_40ns_pocket',9)]:res['overall_40_90ns'][name+'_slope_per_ns']=float(np.polyfit(a[:,0]/1000,a[:,col],1)[0])
res['interpretation']='Using the 40 ns frame as a local reference removes the 0-10 ns relocation from the stability test. A small core RMSD/COM drift and stationary pocket nearest distance/contact fraction over 40-90 ns indicate no ongoing departure; branch RMSDs quantify local flexible motion.'
(OUT/'status_40_70_local_reference.json').write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(res,ensure_ascii=False))


