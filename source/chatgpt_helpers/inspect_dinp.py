import MDAnalysis as mda
from pathlib import Path
TOP=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb")
u=mda.Universe(str(TOP))
for sel in ['resname UNK DINP UNL LIG','not (protein or resname POPC HOH WAT TIP3 NA CL Na+ Cl-)']:
 a=u.select_atoms(sel)
 print('SEL',sel,'n=',a.n_atoms)
 for x in a:
  print(x.index+1,x.resid,x.resname,x.name,x.type,x.element,*(f'{v:.3f}' for v in x.position))
 if a.n_atoms: break
