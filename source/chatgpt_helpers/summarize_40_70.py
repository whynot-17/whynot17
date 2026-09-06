import numpy as np,json
from pathlib import Path
p=Path(r"E:\chatgpt\qc_dinp_40_70ns\dinp_40_70_local_timeseries.csv");a=np.loadtxt(p,delimiter=',',skiprows=1)
for lo in (40,45,50,55,60,65):
 x=a[(a[:,0]>=lo*1000)&(a[:,0]<=(lo+5)*1000+0.1)]; print(f'{lo}-{lo+5}',json.dumps({'core':float(x[:,1].mean()),'core_max':float(x[:,1].max()),'b1':float(x[:,3].mean()),'b2':float(x[:,4].mean()),'com40':float(x[:,5].mean()),'com40_max':float(x[:,5].max()),'pair40':float(x[:,7].mean()),'res40':float(x[:,8].mean()),'near40':float(x[:,9].mean()),'near40_max':float(x[:,9].max())},ensure_ascii=False))
for name,col,thr,op in [('core>2A',1,2,'gt'),('core>2.5A',1,2.5,'gt'),('com>1A',5,1,'gt'),('near>4A',9,4,'gt'),('residue<0.5',8,0.5,'lt')]:
 v=a[:,col]>thr if op=='gt' else a[:,col]<thr; print(name,'fraction',float(v.mean()),'n',int(v.sum()),'firstlast', (float(a[v,0].min()/1000) if v.any() else None,float(a[v,0].max()/1000) if v.any() else None))
for name,col in [('core',1),('com40',5),('near40',9),('pair40',7),('res40',8)]:
 i=int(np.argmax(a[:,col])) if name not in ('pair40','res40') else int(np.argmin(a[:,col])); print(name,'extreme',float(a[i,0]/1000),float(a[i,col]))
