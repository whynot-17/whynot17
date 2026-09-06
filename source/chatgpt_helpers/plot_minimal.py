from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
p=Path(r"E:\chatgpt\qc_dinp_minimal_0_50ns"); a=np.loadtxt(p/'dinp_minimal_timeseries.csv',delimiter=',',skiprows=1); b=np.loadtxt(p/'core_internal_and_pocket_residue_timeseries.csv',delimiter=',',skiprows=1); t=a[:,0]/1000
fig,ax=plt.subplots(2,2,figsize=(11,7),sharex=True)
ax[0,0].plot(t,b[:,1],lw=.7,label='core internal RMSD'); ax[0,0].plot(t,a[:,1],lw=.7,label='core in-pocket RMSD'); ax[0,0].set_ylabel('RMSD (A)'); ax[0,0].legend(frameon=False,fontsize=8); ax[0,0].set_title('DINP aromatic core')
ax[0,1].plot(t,a[:,2],lw=.7,label='branch 1'); ax[0,1].plot(t,a[:,3],lw=.7,label='branch 2'); ax[0,1].set_ylabel('RMSD after core fit (A)'); ax[0,1].legend(frameon=False,fontsize=8); ax[0,1].set_title('Flexible ester/alkyl arms')
ax[1,0].plot(t,a[:,4],lw=.7,color='tab:green'); ax[1,0].set_ylabel('COM displacement (A)'); ax[1,0].set_xlabel('Time (ns)'); ax[1,0].set_title('Whole-DINP COM')
ax[1,1].plot(t,b[:,2],lw=.7,label='pair retention'); ax[1,1].plot(t,b[:,3],lw=.7,label='residue contact fraction'); ax[1,1].set_ylim(0,1); ax[1,1].set_ylabel('Contact fraction'); ax[1,1].set_xlabel('Time (ns)'); ax[1,1].legend(frameon=False,fontsize=8); ax[1,1].set_title('Core–reference pocket contacts')
for x in ax.flat: x.grid(alpha=.2)
fig.tight_layout(); fig.savefig(p/'minimal_test_0_50ns.png',dpi=180); print(p/'minimal_test_0_50ns.png')
