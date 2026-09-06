from pathlib import Path
p=Path(r"E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb")
adj={}
for line in p.read_text().splitlines():
 if line.startswith('CONECT'):
  vals=[]
  for i in range(6,len(line),5):
   s=line[i:i+5].strip()
   if s:
    try: vals.append(int(s))
    except: pass
  if vals: adj.setdefault(vals[0],set()).update(vals[1:])
for i in range(73697,73769): print(i, sorted(adj.get(i,set())))
