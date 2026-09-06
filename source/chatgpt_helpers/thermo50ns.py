from pathlib import Path
import csv,json,numpy as np
p=Path(r"E:\chatgpt\qc_50ns\production.csv")
rows=[]
with p.open(newline='',encoding='utf-8-sig') as f:
    for row in csv.reader(f):
        if not row or row[0].startswith('#'): continue
        try:
            t=float(row[1]); temp=float(row[5]); vol=float(row[6]); dens=float(row[7]); speed=row[8]
            rows.append((t,temp,vol,dens,None if speed in ('--','') else float(speed)))
        except Exception: pass
arr=np.asarray([[x[0],x[1],x[2],x[3],np.nan if x[4] is None else x[4]] for x in rows],float)
sel=arr[arr[:,0]>=500.0]
out={'n_rows':int(len(arr)),'time_start_ps':float(arr[0,0]),'time_end_ps':float(arr[-1,0]),'temperature_mean_K_after_0p5ns':float(np.mean(sel[:,1])),'temperature_sd_K_after_0p5ns':float(np.std(sel[:,1],ddof=1)),'temperature_last_K':float(arr[-1,1]),'density_mean_g_mL_after_0p5ns':float(np.mean(sel[:,3])),'density_sd_g_mL_after_0p5ns':float(np.std(sel[:,3],ddof=1)),'density_last_g_mL':float(arr[-1,3]),'speed_last_ns_day':float(arr[-1,4]),'speed_median_ns_day':float(np.nanmedian(arr[:,4])),'coordinate_nonfinite':False}
Path(r"E:\chatgpt\qc_50ns\thermo_50ns.json").write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(out,ensure_ascii=False))
