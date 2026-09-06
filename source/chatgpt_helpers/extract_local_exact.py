import numpy as np,json
from pathlib import Path
p=Path(r"E:\chatgpt\qc_dinp_40_50ns\dinp_40_50_local_timeseries.csv"); a=np.loadtxt(p,delimiter=',',skiprows=1)
for ns in (40,45,50):
 i=np.argmin(abs(a[:,0]-ns*1000)); print(ns,json.dumps({'time_ps':float(a[i,0]),'core_rmsd_vs40':float(a[i,1]),'core_internal':float(a[i,2]),'b1':float(a[i,3]),'b2':float(a[i,4]),'com_drift_vs40':float(a[i,5]),'com_vsmin':float(a[i,6]),'pair40':float(a[i,7]),'res40':float(a[i,8]),'near40':float(a[i,9])},ensure_ascii=False))
