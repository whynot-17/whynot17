from pathlib import Path
import json,numpy as np
p=Path(r"E:\chatgpt\qc_dinp_40_50ns\status_40_50_local_reference.json"); s=json.loads(p.read_text(encoding='utf-8')); a=np.loadtxt(p.parent/'dinp_40_50_local_timeseries.csv',delimiter=',',skiprows=1)
s['two_ns_blocks']={}
for lo in (40,42,44,46,48):
 x=a[(a[:,0]>=lo*1000)&(a[:,0]<=(lo+2)*1000+0.1)]; s['two_ns_blocks'][f'{lo}-{lo+2}ns']={'core_rmsd_vs_40ns_mean_A':float(x[:,1].mean()),'branch1_rmsd_mean_A':float(x[:,3].mean()),'branch2_rmsd_mean_A':float(x[:,4].mean()),'com_drift_vs_40ns_mean_A':float(x[:,5].mean()),'pair_retention_vs_40ns_mean':float(x[:,7].mean()),'pocket_residue_fraction_vs_40ns_mean':float(x[:,8].mean()),'nearest_distance_to_40ns_pocket_mean_A':float(x[:,9].mean())}
for ns in (40,45,50):
 i=int(np.argmin(abs(a[:,0]-ns*1000))); s[f'nearest_{ns}ns']={'time_ps':float(a[i,0]),'core_rmsd_vs_40ns_A':float(a[i,1]),'core_internal_rmsd_A':float(a[i,2]),'branch1_rmsd_vs_40ns_A':float(a[i,3]),'branch2_rmsd_vs_40ns_A':float(a[i,4]),'com_drift_vs_40ns_A':float(a[i,5]),'com_displacement_vs_minimized_A':float(a[i,6]),'pair_retention_vs_40ns':float(a[i,7]),'pocket_residue_fraction_vs_40ns':float(a[i,8]),'nearest_distance_to_40ns_pocket_A':float(a[i,9])}
p.write_text(json.dumps(s,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps({'blocks':s['two_ns_blocks'],'exact50':s['nearest_50ns']},ensure_ascii=False))
