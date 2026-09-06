from pathlib import Path
import time, shutil, json
src=Path(r"E:\chatgpt\ptger4_membrane_md_20260904")
dst=Path(r"E:\chatgpt\qc_50ns")
dst.mkdir(parents=True, exist_ok=True)
for name in ("production.dcd","production.csv","production.chk","final.pdb"):
    p=src/name
    if not p.exists():
        raise SystemExit(f"missing {p}")
    a=(p.stat().st_size,p.stat().st_mtime_ns)
    time.sleep(2)
    b=(p.stat().st_size,p.stat().st_mtime_ns)
    print(name, a, b, "stable", a==b, flush=True)
    shutil.copy2(p,dst/name)
print(json.dumps({"src":str(src),"dst":str(dst),"files":[x.name for x in dst.iterdir()]},ensure_ascii=False))
