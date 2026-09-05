from pathlib import Path
import itertools, json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb")
REF=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\minimized.pdb")
DCD=Path(r"E:\chatgpt\qc_50ns\production.dcd")
OUT=Path(r"E:\chatgpt\qc_dinp_minimal_0_50ns"); OUT.mkdir(parents=True,exist_ok=True)
CUTOFF=4.5

def kabsch(mobile,target):
    cm,ct=mobile.mean(0),target.mean(0)
    u,_,vt=np.linalg.svd((mobile-cm).T@(target-ct)); r=u@vt
    if np.linalg.det(r)<0: u[:,-1]*=-1; r=u@vt
    return r,ct-cm@r

def unwrap_group(pos,box):
    # Keep all ligand atoms in the same image as the first aromatic atom.
    a=pos[0].copy(); d=pos-a; d-=box*np.round(d/box); return a+d

def rmsd(a,b): return float(np.sqrt(np.mean(np.sum((a-b)**2,axis=1))))

u=mda.Universe(str(TOP),str(DCD)); ref=mda.Universe(str(REF))
lig=u.select_atoms('resname UNK DINP UNL LIG'); rlig=ref.select_atoms('resname UNK DINP UNL LIG')
if lig.n_atoms!=72: raise RuntimeError(f'Unexpected ligand atom count: {lig.n_atoms}')
# DINP atom naming in this system: six aromatic carbons, two ester/alkyl branches.
core_names=['C1x','C2x','C3x','C4x','C5x','C6x']
branch1_names=['C7x','O1x','O2x','C8x','C9x','C10x','C11x','C12x','C13x','C14x','C15x','C16x']
branch2_names=['C17x','O3x','O4x','C18x','C19x','C20x','C21x','C22x','C23x','C24x','C25x','C26x']
def indices(names,ag):
    return np.array([ag.indices[list(ag.names).index(n)] for n in names],dtype=int)
# Use indices relative to full ligand selection.
lig_name_to_i={name:i for i,name in enumerate(lig.names)}
core_i=np.array([lig_name_to_i[n] for n in core_names],int); b1_i=np.array([lig_name_to_i[n] for n in branch1_names],int); b2_i=np.array([lig_name_to_i[n] for n in branch2_names],int)
ref_name_to_i={name:i for i,name in enumerate(rlig.names)}
rcore=rlig.positions[ [ref_name_to_i[n] for n in core_names] ].copy().astype(float)
rbranch1=rlig.positions[[ref_name_to_i[n] for n in branch1_names]].copy().astype(float)
rbranch2=rlig.positions[[ref_name_to_i[n] for n in branch2_names]].copy().astype(float)
rall=rlig.positions.copy().astype(float); rcom=rall.mean(0)
ca=u.select_atoms('protein and name CA'); rca=ref.select_atoms('protein and name CA').positions.copy().astype(float)
prot=u.select_atoms('protein and not name H*'); rprot=ref.select_atoms('protein and not name H*').positions.copy().astype(float)
# Define the pocket from minimized-reference core contacts.
dref=distance_array(rcore,rprot)
pocket_mask=np.any(dref<=CUTOFF,axis=0); pair_idx=np.argwhere(dref<=CUTOFF)
if pair_idx.size==0: raise RuntimeError('Reference core has no pocket contacts')
rows=[]; n=0
for ts in u.trajectory:
    t=float(ts.time)
    if t>50000.1: break
    n+=1; box=np.asarray(ts.dimensions[:3],float)
    raw=lig.positions.copy().astype(float); lu=unwrap_group(raw,box)
    # Protein-fit frame: choose the ligand periodic image closest to the reference core.
    # This prevents a whole-box jump from appearing as a false core RMSD spike.
    rotp,tranp=kabsch(ca.positions.copy().astype(float),rca); best=None
    for sh in itertools.product([-1,0,1],repeat=3):
        image=lu+np.asarray(sh,dtype=float)*box
        fit=image@rotp+tranp
        score=float(np.linalg.norm(fit[core_i].mean(0)-rcore.mean(0)))
        if best is None or score<best[0]: best=(score,image,fit)
    image=best[1]; fit_lig=best[2]; core=image[core_i]; b1=image[b1_i]; b2=image[b2_i]
    core_p=fit_lig[core_i]; core_r=rmsd(core_p,rcore)
    # Core-fit branches: isolate flexible-chain rearrangement from whole-ligand translation/rotation.
    rotc,tranc=kabsch(core,rcore); b1_c=b1@rotc+tranc; b2_c=b2@rotc+tranc
    b1_r=rmsd(b1_c,rbranch1); b2_r=rmsd(b2_c,rbranch2)
    # Whole DINP COM displacement after protein fit, using the same closest periodic image.
    com_disp=float(np.linalg.norm(fit_lig.mean(0)-rcom))
    # Core-pocket contacts measured with minimum-image distances; pocket is fixed from reference.
    dcur=distance_array(raw[core_i],prot.positions,box=ts.dimensions)
    pair_d=dcur[pair_idx[:,0],pair_idx[:,1]]
    retention=float(np.mean(pair_d<=CUTOFF)); pocket_contact=float(np.mean(np.any(dcur[:,pocket_mask]<=CUTOFF,axis=0))); nearest=float(dcur[:,pocket_mask].min())
    rows.append((t,core_r,b1_r,b2_r,com_disp,retention,pocket_contact,nearest))
arr=np.asarray(rows,float)
np.savetxt(OUT/'dinp_minimal_timeseries.csv',arr,delimiter=',',header='time_ps,core_rmsd_after_protein_fit_A,branch1_rmsd_after_core_fit_A,branch2_rmsd_after_core_fit_A,ligand_com_displacement_A,core_pocket_pair_retention_4p5A,pocket_atom_contact_fraction_4p5A,core_pocket_nearest_distance_A',comments='')
def summary(a):
    return {'n_frames':int(len(a)),'time_start_ps':float(a[0,0]),'time_end_ps':float(a[-1,0]),'core_rmsd_mean_A':float(a[:,1].mean()),'core_rmsd_sd_A':float(a[:,1].std(ddof=1)),'core_rmsd_min_A':float(a[:,1].min()),'core_rmsd_max_A':float(a[:,1].max()),'branch1_rmsd_mean_A':float(a[:,2].mean()),'branch1_rmsd_sd_A':float(a[:,2].std(ddof=1)),'branch1_rmsd_min_A':float(a[:,2].min()),'branch1_rmsd_max_A':float(a[:,2].max()),'branch2_rmsd_mean_A':float(a[:,3].mean()),'branch2_rmsd_sd_A':float(a[:,3].std(ddof=1)),'branch2_rmsd_min_A':float(a[:,3].min()),'branch2_rmsd_max_A':float(a[:,3].max()),'com_displacement_mean_A':float(a[:,4].mean()),'com_displacement_sd_A':float(a[:,4].std(ddof=1)),'com_displacement_min_A':float(a[:,4].min()),'com_displacement_max_A':float(a[:,4].max()),'core_pocket_pair_retention_mean':float(a[:,5].mean()),'core_pocket_pair_retention_sd':float(a[:,5].std(ddof=1)),'core_pocket_pair_retention_min':float(a[:,5].min()),'core_pocket_pair_retention_max':float(a[:,5].max()),'pocket_atom_contact_fraction_mean':float(a[:,6].mean()),'pocket_atom_contact_fraction_min':float(a[:,6].min()),'pocket_atom_contact_fraction_max':float(a[:,6].max()),'core_pocket_nearest_distance_mean_A':float(a[:,7].mean()),'core_pocket_nearest_distance_min_A':float(a[:,7].min()),'core_pocket_nearest_distance_max_A':float(a[:,7].max())}
result={'definitions':{'core':'UNK heavy atoms C1x-C6x (aromatic ring)','branch1':'UNK heavy atoms C7x,O1x,O2x,C8x-C16x','branch2':'UNK heavy atoms C17x,O3x,O4x,C18x-C26x','pocket':'protein heavy atoms within 4.5 A of the minimized-reference aromatic core','contact_retention':'fraction of minimized-reference core-pocket atom pairs (4.5 A) still within 4.5 A','rmsd_alignment':'core RMSD after protein C-alpha fit; each branch RMSD after aromatic-core fit','com':'whole-DINP COM after protein fit and periodic-image correction'},'reference_pocket_atoms':int(pocket_mask.sum()),'reference_core_pocket_pairs':int(len(pair_idx)),'overall_0_50ns':summary(arr),'first_half_0_25ns':summary(arr[arr[:,0]<=25000.1]),'second_half_25_50ns':summary(arr[arr[:,0]>=25000.0])}
for ns in (0,10,25,40,50):
 i=int(np.argmin(abs(arr[:,0]-ns*1000.0))); result[f'nearest_{ns}ns']={'time_ps':float(arr[i,0]),'core_rmsd_A':float(arr[i,1]),'branch1_rmsd_A':float(arr[i,2]),'branch2_rmsd_A':float(arr[i,3]),'com_displacement_A':float(arr[i,4]),'core_pocket_pair_retention':float(arr[i,5]),'pocket_atom_contact_fraction':float(arr[i,6]),'core_pocket_nearest_distance_A':float(arr[i,7])}
# Simple drift estimates over 0-50 ns.
for name,col in [('core_rmsd',1),('branch1_rmsd',2),('branch2_rmsd',3),('com_displacement',4),('contact_retention',5),('nearest_distance',7)]: result['overall_0_50ns'][name+'_slope_per_ns']=float(np.polyfit(arr[:,0]/1000.0,arr[:,col],1)[0])
result['interpretation']='The minimal test should be read with the selected definitions: the aromatic core is a six-carbon ring, branches are the two ester/alkyl arms, and the pocket is defined from minimized-reference core contacts. A low core RMSD and COM displacement with fluctuating arm RMSDs and retained core-pocket contacts indicate pocket-bound flexible rearrangement.'
(OUT/'status_minimal_0_50ns.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False))
