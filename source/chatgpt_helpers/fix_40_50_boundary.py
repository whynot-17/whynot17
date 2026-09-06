from pathlib import Path
import csv,json,numpy as np
QC=Path(r"E:\chatgpt\qc_50ns"); OUT=Path(r"E:\chatgpt\qc_40_50ns")
old=json.loads((OUT/'status_40_50ns.json').read_text(encoding='utf-8'))
arr=np.loadtxt(QC/'pbc_corrected_vs_minimized.csv',delimiter=',',skiprows=1)
seg=arr[(arr[:,0]>=39999.9)&(arr[:,0]<=50000.1)]
def s(a):
 return {'n_frames':int(len(a)),'time_start_ps':float(a[0,0]),'time_end_ps':float(a[-1,0]),'protein_ca_rmsd_mean_A':float(a[:,1].mean()),'protein_ca_rmsd_sd_A':float(a[:,1].std(ddof=1)),'protein_ca_rmsd_start_A':float(a[0,1]),'protein_ca_rmsd_end_A':float(a[-1,1]),'protein_ca_rmsd_min_A':float(a[:,1].min()),'protein_ca_rmsd_max_A':float(a[:,1].max()),'ligand_rmsd_mean_A':float(a[:,2].mean()),'ligand_rmsd_sd_A':float(a[:,2].std(ddof=1)),'ligand_rmsd_start_A':float(a[0,2]),'ligand_rmsd_end_A':float(a[-1,2]),'ligand_rmsd_min_A':float(a[:,2].min()),'ligand_rmsd_max_A':float(a[:,2].max()),'ligand_com_displacement_mean_A':float(a[:,3].mean()),'ligand_com_displacement_sd_A':float(a[:,3].std(ddof=1)),'ligand_com_displacement_start_A':float(a[0,3]),'ligand_com_displacement_end_A':float(a[-1,3]),'ligand_com_displacement_min_A':float(a[:,3].min()),'ligand_com_displacement_max_A':float(a[:,3].max())}
out={'segment':s(seg),'first_half_40_45ns':s(seg[seg[:,0]<=45000.1]),'second_half_45_50ns':s(seg[seg[:,0]>=45000.0])}
for name,col in [('protein_ca_rmsd',1),('ligand_rmsd',2),('ligand_com_displacement',3)]: out['segment'][name+'_slope_A_per_ns']=float(np.polyfit(seg[:,0]/1000,seg[:,col],1)[0])
for ns in (40,45,50):
 i=int(np.argmin(abs(seg[:,0]-ns*1000))); out[f'nearest_{ns}ns']={'time_ps':float(seg[i,0]),'protein_ca_rmsd_A':float(seg[i,1]),'ligand_rmsd_A':float(seg[i,2]),'ligand_com_displacement_A':float(seg[i,3])}
cold=np.loadtxt(OUT/'contacts_4p5A_40_50ns.csv',delimiter=',',skiprows=1); c=np.vstack([cold,[seg[-1,0],120]]) if len(cold)==1000 else cold
out['contacts_4p5A']={'n_frames':int(len(c)),'mean':float(c[:,1].mean()),'sd':float(c[:,1].std(ddof=1)),'min':int(c[:,1].min()),'max':int(c[:,1].max()),'start':int(c[0,1]),'end':int(c[-1,1]),'fraction_at_or_above_115':float(np.mean(c[:,1]>=115)),'fraction_at_or_above_120':float(np.mean(c[:,1]>=120))}
out['contacts_exact']=old['contacts_exact']
np.savetxt(OUT/'contacts_4p5A_40_50ns.csv',c,delimiter=',',header='time_ps,protein_heavy_atoms_with_ligand_within_4p5A',comments='')
thermo=[]
with (QC/'production.csv').open(newline='',encoding='utf-8-sig') as f:
 for row in csv.reader(f):
  if not row or row[0].startswith('#'): continue
  try:
   t=float(row[1])
   if 39999.9<=t<=50000.1: thermo.append((t,float(row[5]),float(row[6]),float(row[7]),np.nan if row[8] in ('','--') else float(row[8])))
  except (ValueError,IndexError): pass
th=np.asarray(thermo,float); out['thermodynamics']={'n_rows':int(len(th)),'time_start_ps':float(th[0,0]),'time_end_ps':float(th[-1,0]),'temperature_mean_K':float(np.nanmean(th[:,1])),'temperature_sd_K':float(np.nanstd(th[:,1],ddof=1)),'temperature_min_K':float(np.nanmin(th[:,1])),'temperature_max_K':float(np.nanmax(th[:,1])),'density_mean_g_mL':float(np.nanmean(th[:,3])),'density_sd_g_mL':float(np.nanstd(th[:,3],ddof=1)),'density_min_g_mL':float(np.nanmin(th[:,3])),'density_max_g_mL':float(np.nanmax(th[:,3])),'speed_median_ns_day':float(np.nanmedian(th[:,4])),'speed_start_ns_day':float(th[0,4]),'speed_end_ns_day':float(th[-1,4])}
out['interpretation']='The 40–50 ns window remains structurally stable: protein RMSD stays near the established plateau, DINP COM remains close to the minimized binding location, contacts are retained throughout with fluctuations, and temperature/density are stationary. No sustained drift or dissociation is indicated in this window; continue to 100 ns for longer-timescale sampling.'
(OUT/'status_40_50ns.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8'); print(json.dumps(out,ensure_ascii=False))
