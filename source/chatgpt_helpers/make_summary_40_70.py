from pathlib import Path
import json,numpy as np
p=Path(r"E:\chatgpt\qc_dinp_40_70ns\status_40_70_local_reference.json"); s=json.loads(p.read_text(encoding='utf-8')); a=np.loadtxt(p.parent/'dinp_40_70_local_timeseries.csv',delimiter=',',skiprows=1); thermo=json.loads((p.parent/'thermo_40_70.json').read_text(encoding='utf-8'))
s['five_ns_blocks']={}
for lo in (40,45,50,55,60,65):
 x=a[(a[:,0]>=lo*1000)&(a[:,0]<=(lo+5)*1000+0.1)]; s['five_ns_blocks'][f'{lo}-{lo+5}ns']={'core_rmsd_vs_40ns_mean_A':float(x[:,1].mean()),'branch1_rmsd_mean_A':float(x[:,3].mean()),'branch2_rmsd_mean_A':float(x[:,4].mean()),'com_drift_vs_40ns_mean_A':float(x[:,5].mean()),'pair_retention_vs_40ns_mean':float(x[:,7].mean()),'pocket_residue_fraction_vs_40ns_mean':float(x[:,8].mean()),'nearest_distance_to_40ns_pocket_mean_A':float(x[:,9].mean())}
# Event counts and locations.
s['transient_events']={'core_rmsd_gt_2A_fraction':float(np.mean(a[:,1]>2.0)),'core_rmsd_gt_2p5A_fraction':float(np.mean(a[:,1]>2.5)),'core_rmsd_max_time_ns':float(a[np.argmax(a[:,1]),0]/1000),'core_rmsd_max_A':float(a[:,1].max()),'com_drift_gt_1A_fraction':float(np.mean(a[:,5]>1.0)),'com_drift_max_time_ns':float(a[np.argmax(a[:,5]),0]/1000),'com_drift_max_A':float(a[:,5].max()),'residue_contact_fraction_lt_0p5_fraction':float(np.mean(a[:,8]<0.5)),'nearest_distance_gt_4A_fraction':float(np.mean(a[:,9]>4.0))}
s['thermodynamics_40_70ns']=thermo
s.setdefault('files',{}); s['files']['thermo_json']='E:/chatgpt/qc_dinp_40_70ns/thermo_40_70.json'; s['files']['timeseries_csv']='E:/chatgpt/qc_dinp_40_70ns/dinp_40_70_local_timeseries.csv'; s['interpretation']='With the 40 ns frame as the local reference, the 40-70 ns window remains in the same pocket microstate. Core RMSD and COM show modest fluctuations with a few brief excursions, but the nearest pocket distance stays within 4.5 A for every frame and pocket-residue contact fraction averages about 81%. There is no sustained dissociation signal; continuation to 100 ns remains reasonable.'
p.write_text(json.dumps(s,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps({'overall':s['overall_40_70ns'],'events':s['transient_events'],'thermo':thermo},ensure_ascii=False))

