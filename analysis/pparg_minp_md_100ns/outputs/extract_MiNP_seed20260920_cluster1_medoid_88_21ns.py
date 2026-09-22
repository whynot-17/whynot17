from pathlib import Path
import json, re
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

run = Path(r"E:\chatgpt\pparg_minp_md\run_20260920_seed20260920")
top = Path(r"E:\chatgpt\pparg_minp_md\system\seed20260917\system_topology.pdb")
ref_path = Path(r"E:\chatgpt\pparg_minp_md\system\seed20260917\initial.pdb")
dcd = run / "production.dcd"
template = Path(r"E:\chatgpt\pparg_minp_md\prepared\MiNP_CID110394.pdbqt")
target_ns = 88.21

u = mda.Universe(str(top), str(dcd))
ref = mda.Universe(str(ref_path))
print("loaded_universes", flush=True)
idx = 8820  # 0.01 ns first frame; 88.21 ns medoid maps to frame 8820
ts = None
bb_xyz = None
protein_xyz = None
ligand_xyz = None
frame_box = None
for _i, _ts in enumerate(u.trajectory):
    if _i == idx:
        ts = _ts.copy()
        actual_ns = float(_ts.time/1000.0)
        bb = u.select_atoms("protein and backbone")
        protein = u.select_atoms("protein")
        ligand = u.select_atoms("resname UNK")
        bb_xyz = bb.positions.copy()
        protein_xyz = protein.positions.copy()
        ligand_xyz = ligand.positions.copy()
        frame_box = np.asarray(_ts.dimensions[:3], dtype=float).copy()
        break
if ts is None:
    raise RuntimeError(f"Could not reach medoid frame {idx}")
print("reached_frame", idx, actual_ns, flush=True)

refbb = ref.select_atoms("protein and backbone")
refbb_xyz = refbb.positions.copy()
bb_center = bb_xyz.mean(axis=0)
refbb_center = refbb_xyz.mean(axis=0)
rotation, fit_rmsd = align.rotation_matrix(bb_xyz-bb_center, refbb_xyz-refbb_center)
print("aligned_protein", float(fit_rmsd), flush=True)

def apply_rot(x, mobile_center, reference_center):
    v=x-mobile_center
    out=np.empty_like(v)
    out[:,0]=v[:,0]*rotation[0,0]+v[:,1]*rotation[0,1]+v[:,2]*rotation[0,2]+reference_center[0]
    out[:,1]=v[:,0]*rotation[1,0]+v[:,1]*rotation[1,1]+v[:,2]*rotation[1,2]+reference_center[1]
    out[:,2]=v[:,0]*rotation[2,0]+v[:,1]*rotation[2,1]+v[:,2]*rotation[2,2]+reference_center[2]
    return out

box=frame_box
if not np.all(box>0):
    box=np.array([80.0,80.0,80.0])
print("box_ready", box, flush=True)
lig_raw=ligand_xyz
print("ligand_raw", ligand_xyz.shape, flush=True)
lig_shift=lig_raw + np.rint((bb_center-lig_raw.mean(axis=0))/box)*box
print("ligand_shift", flush=True)
ligand_aligned=apply_rot(lig_shift, bb_center, refbb_center)
print("aligned_ligand", ligand_aligned.shape, flush=True)

# Write a protein+ligand complex PDB manually in the reference-aligned frame.
# This avoids a Windows/PDBWriter crash on the very large solvated topology.
coord_map={int(i):xyz for i,xyz in zip(ligand.indices, ligand_aligned)}
selected=set(int(i) for i in protein.indices) | set(int(i) for i in ligand.indices)
ligand_indices=set(int(i) for i in ligand.indices)
pdb_path=run/"MiNP_seed20260920_cluster1_medoid_88.21ns_complex.pdb"
out_pdb=[f"REMARK MiNP-PPARG Cluster 1 medoid; actual {actual_ns:.5f} ns", f"REMARK Protein-backbone fit RMSD {float(fit_rmsd):.4f} A; aligned to initial.pdb"]
atom_i=-1
for line in ref_path.read_text(encoding="utf-8", errors="ignore").splitlines():
    if line.startswith(("ATOM", "HETATM")):
        atom_i += 1
        if atom_i in selected:
            if atom_i in ligand_indices:
                x,y,z=coord_map[atom_i]
                line=line[:30]+f"{x:8.3f}{y:8.3f}{z:8.3f}"+line[54:]
            out_pdb.append(line)
    elif line.startswith("TER") and atom_i < min(ligand.indices):
        out_pdb.append(line)
out_pdb.append("END")
pdb_path.write_text("\n".join(out_pdb)+"\n",encoding="utf-8")
print("wrote_pdb", pdb_path, flush=True)

# PDBQT template maps source ligand atom index -> PDBQT atom number.
pairs=[(13,1),(1,2),(14,3),(2,4),(15,5),(16,6),(17,7),(18,8),(19,9),(20,10),(21,11),(4,12),
       (3,13),(45,14),(10,15),(9,16),(7,17),(5,18),(6,19),(8,20),(11,21),(12,22)]
src_to_pdbqt={src:pq for src,pq in pairs}
pdbqt_to_src={pq:src for src,pq in pairs}
lines=template.read_text(encoding="utf-8").splitlines()
out=[]
out.append(f"REMARK MiNP-PPARG Cluster 1 medoid; target {target_ns:.2f} ns; actual {actual_ns:.5f} ns")
out.append(f"REMARK Protein-backbone fit RMSD {float(fit_rmsd):.4f} A; coordinates aligned to initial.pdb")
for line in lines:
    if line.startswith(("ATOM","HETATM")):
        pqnum=int(line[6:11])
        srcidx=pdbqt_to_src[pqnum]
        x,y,z=ligand_aligned[srcidx-1]
        line=line[:30]+f"{x:8.3f}{y:8.3f}{z:8.3f}"+line[54:]
    out.append(line)
pdbqt_path=run/"MiNP_seed20260920_cluster1_medoid_88.21ns_ligand.pdbqt"
pdbqt_path.write_text("\n".join(out)+"\n",encoding="utf-8")
print("wrote_pdbqt", pdbqt_path, flush=True)

meta={
 "cluster":1,
 "cluster_occupancy":0.8633333333333333,
 "target_time_ns":target_ns,
 "actual_time_ns":actual_ns,
 "trajectory_frame_index":int(idx),
 "protein_backbone_fit_rmsd_A":float(fit_rmsd),
 "complex_pdb":str(pdb_path),
 "ligand_pdbqt":str(pdbqt_path),
 "pdbqt_template":str(template),
 "pdb_atoms":int(len(selected)),
 "ligand_atoms_in_pdbqt":22,
}
(run/"MiNP_seed20260920_cluster1_medoid_88.21ns_metadata.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
print(json.dumps(meta,indent=2))

