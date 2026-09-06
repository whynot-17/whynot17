import numpy as np,json
from pathlib import Path
p=Path(r"E:\chatgpt\qc_dinp_minimal_0_50ns"); a=np.loadtxt(p/'dinp_minimal_timeseries.csv',delimiter=',',skiprows=1); b=np.loadtxt(p/'core_internal_and_pocket_residue_timeseries.csv',delimiter=',',skiprows=1)
def su(x,y):
 return {'n_frames':len(x),'core_in_pocket_mean_A':float(x[:,1].mean()),'core_in_pocket_sd_A':float(x[:,1].std(ddof=1)),'core_in_pocket_start_A':float(x[0,1]),'core_in_pocket_end_A':float(x[-1,1]),'branch1_mean_A':float(x[:,2].mean()),'branch2_mean_A':float(x[:,3].mean()),'com_mean_A':float(x[:,4].mean()),'com_max_A':float(x[:,4].max()),'pair_retention_mean':float(y[:,2].mean()),'residue_contact_fraction_mean':float(y[:,3].mean()),'nearest_mean_A':float(y[:,4].mean()),'nearest_le_4p5_fraction':float(np.mean(y[:,4]<=4.5))}
for name,lo,hi in [('0-5ns',0,5000.1),('5-10ns',5000,10000.1),('10-50ns',10000,50000.1),('40-50ns',40000,50000.1)]:
 ma=(a[:,0]>=lo)&(a[:,0]<=hi); mb=(b[:,0]>=lo)&(b[:,0]<=hi); print(name,json.dumps(su(a[ma],b[mb]),ensure_ascii=False))
