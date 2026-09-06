from pathlib import Path
import time, shutil, json
src=Path(r"E:\chatgpt\ptger4_membrane_md_20260904"); dst=Path(r"E:\chatgpt\qc_70ns"); dst.mkdir(parents=True,exist_ok=True)
for name in ('production.dcd','production.csv'):
 p=src/name
 for _ in range(6):
  a=(p.stat().st_size,p.stat().st_mtime_ns); time.sleep(2); b=(p.stat().st_size,p.stat().st_mtime_ns)
  if a==b: break
 else: raise RuntimeError(f'{name} did not stabilize')
 print(name,'stable',a,flush=True); shutil.copy2(p,dst/name)
print(json.dumps({n:(dst/n).stat().st_size for n in ('production.dcd','production.csv')}))
