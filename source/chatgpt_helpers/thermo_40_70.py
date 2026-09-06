from pathlib import Path
import csv,json,numpy as np
p=Path(r"E:\chatgpt\qc_70ns\production.csv"); rows=[]
with p.open(newline='',encoding='utf-8-sig') as f:
 for row in csv.reader(f):
  if not row or row[0].startswith('#'): continue
  try:
   t=float(row[1]);
   if 39999.9<=t<=70000.1: rows.append((t,float(row[5]),float(row[7]),float(row[8]) if row[8] not in ('','--') else np.nan))
  except (ValueError,IndexError): pass
a=np.asarray(rows,float); out={'n_rows':int(len(a)),'time_start_ns':float(a[0,0]/1000),'time_end_ns':float(a[-1,0]/1000),'temperature_mean_K':float(np.nanmean(a[:,1])),'temperature_sd_K':float(np.nanstd(a[:,1],ddof=1)),'temperature_min_K':float(np.nanmin(a[:,1])),'temperature_max_K':float(np.nanmax(a[:,1])),'density_mean_g_mL':float(np.nanmean(a[:,2])),'density_sd_g_mL':float(np.nanstd(a[:,2],ddof=1)),'density_min_g_mL':float(np.nanmin(a[:,2])),'density_max_g_mL':float(np.nanmax(a[:,2])),'speed_median_ns_day':float(np.nanmedian(a[:,3])),'speed_start_ns_day':float(a[0,3]),'speed_end_ns_day':float(a[-1,3])}; Path(r"E:\chatgpt\qc_dinp_40_70ns\thermo_40_70.json").write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(out,ensure_ascii=False))
