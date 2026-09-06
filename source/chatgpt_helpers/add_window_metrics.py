from pathlib import Path
import json,numpy as np
p=Path(r"E:\chatgpt\qc_dinp_minimal_0_50ns\minimal_test_summary.json"); s=json.loads(p.read_text(encoding='utf-8')); a=np.loadtxt(p.parent/'dinp_minimal_timeseries.csv',delimiter=',',skiprows=1); b=np.loadtxt(p.parent/'core_internal_and_pocket_residue_timeseries.csv',delimiter=',',skiprows=1)
def su(lo,hi):
 x=a[(a[:,0]>=lo)&(a[:,0]<=hi)]; y=b[(b[:,0]>=lo)&(b[:,0]<=hi)]
 return {'n_frames':int(len(x)),'core_in_pocket_mean_A':float(x[:,1].mean()),'core_in_pocket_start_A':float(x[0,1]),'core_in_pocket_end_A':float(x[-1,1]),'branch1_mean_A':float(x[:,2].mean()),'branch2_mean_A':float(x[:,3].mean()),'com_mean_A':float(x[:,4].mean()),'com_max_A':float(x[:,4].max()),'pair_retention_mean':float(y[:,2].mean()),'residue_contact_fraction_mean':float(y[:,3].mean()),'nearest_mean_A':float(y[:,4].mean()),'nearest_le_4p5_fraction':float(np.mean(y[:,4]<=4.5))}
s['windows']={'0-5ns':su(0,5000.1),'5-10ns':su(5000,10000.1),'10-50ns':su(10000,50000.1),'40-50ns':su(40000,50000.1)}
p.write_text(json.dumps(s,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(s['windows'],ensure_ascii=False))
