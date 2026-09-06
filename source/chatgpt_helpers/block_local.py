import numpy as np,json
from pathlib import Path
p=Path(r"E:\chatgpt\qc_dinp_40_50ns\dinp_40_50_local_timeseries.csv");a=np.loadtxt(p,delimiter=',',skiprows=1)
for lo in (40,42,44,46,48):
 x=a[(a[:,0]>=lo*1000)&(a[:,0]<= (lo+2)*1000+0.1)]; print(f'{lo}-{lo+2}',json.dumps({'core':float(x[:,1].mean()),'b1':float(x[:,3].mean()),'b2':float(x[:,4].mean()),'com':float(x[:,5].mean()),'pair40':float(x[:,7].mean()),'res40':float(x[:,8].mean()),'near40':float(x[:,9].mean())},ensure_ascii=False))
