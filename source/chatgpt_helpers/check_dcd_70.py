import MDAnalysis as mda
u=mda.Universe(r'E:\mcop\mcop\analysis\dinp_crc_structural_pipeline\outputs\ptger4_membrane_smoke_1ns\system_built.pdb',r'E:\chatgpt\qc_70ns\production.dcd'); print('frames',len(u.trajectory),'start_ps',u.trajectory[0].time,'end_ps',u.trajectory[-1].time)
