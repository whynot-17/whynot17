# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Dose-response spline panels -> LineTrend/plot_trend.py -> param inherit
# (b) Sex overlay panels -> LineTrend/plot_sweep.py -> param inherit
# RULE: native scripts were not semantically compatible with survey-weighted
# spline inputs; the figure below inherits their restrained line/CI language.
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 8, "figure.titlesize": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "legend.frameon": False,
})
# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666", "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999"]
DIVERGING = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED = "#B2182B"; GREY = "#999999"; BLACK = "#222222"
# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "savefig.bbox": "tight", "savefig.dpi": 300})
def save_cns_figure(fig, filename):
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

import hashlib, importlib.util, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd
from scipy.special import expit
from scipy.stats import t, chi2

ROOT=Path(__file__).resolve().parents[3]; FW=ROOT/'analysis'/'disease_agnostic_environmental_framework'; OUT=Path(__file__).resolve().parent; DATA=ROOT/'work'/'nhanes_phase2a'/'data'
def load(p,n):
 s=importlib.util.spec_from_file_location(n,p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def num(x): return pd.to_numeric(x,errors='coerce')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rcs(x,k):
 k=np.asarray(k,float); out=[x]
 for j in range(len(k)-2): out.append(np.maximum(x-k[j],0)**3-np.maximum(x-k[-2],0)**3*(k[-1]-k[j])/(k[-1]-k[-2])+np.maximum(x-k[-1],0)**3*(k[-2]-k[j])/(k[-1]-k[j]))
 return np.column_stack(out)
def Xmat(d,sex=None,quart=False,spline=False,cr=False):
 x=pd.DataFrame({'I':1.},index=d.index); z=num(d.axis_log2)
 if quart: [x.__setitem__(f'Q{q}',d.q.eq(q).astype(float)) for q in [2,3,4]]; expcols=['Q2','Q3','Q4']
 elif spline:
  b=rcs(z,d.knots); expcols=[]
  for j in range(b.shape[1]): x[f'spline{j}']=b[:,j]; expcols.append(f'spline{j}')
 else: x['axis_log2']=z; expcols=['axis_log2']
 a=num(d.age)-50; x['age']=a; x['age2']=a*a; x['pir']=num(d.pir)
 if sex is None: x['female']=d.sex.eq('Female').astype(float); [x.__setitem__(f'{c}={v}',d[c].eq(v).astype(float)) for c,vs in [('race',['Mexican American','Other Hispanic','Non-Hispanic Black','Other/Multi']),('smoking',['Former','Current'])] for v in vs]; [x.__setitem__(f'cycle={v}',d.cycle.eq(v).astype(float)) for v in sorted(d.cycle.unique())[1:]]
 else: [x.__setitem__(f'{c}={v}',d[c].eq(v).astype(float)) for c,vs in [('race',['Mexican American','Other Hispanic','Non-Hispanic Black','Other/Multi']),('smoking',['Former','Current'])] for v in vs]; [x.__setitem__(f'cycle={v}',d.cycle.eq(v).astype(float)) for v in sorted(d.cycle.unique())[1:]]
 if sex is None and not quart and not spline: x['axis_log2:female']=x.axis_log2*x.female; expcols+=['axis_log2:female']
 if cr: x['creatinine_log2']=num(d.creatinine_log2); x['creatinine:female']=x.creatinine_log2*x.female
 return x,np.array(expcols)
def fit(d,sex=None,quart=False,spline=False,cr=False):
 req=['outcome','axis_log2','age','pir','race','smoking','pooled_weight','psu','strata','cycle']+(['creatinine_log2'] if cr else [])
 if sex is not None: d=d.loc[d.sex.eq(sex)]
 d=d.dropna(subset=req).loc[lambda q:q.pooled_weight.gt(0)].reset_index(drop=True); y=num(d.outcome).to_numpy(float); base={'N':len(d),'cases':int(y.sum()) if len(y) else 0,'controls':int(len(y)-y.sum()) if len(y) else 0}
 if len(d)==0 or y.sum() in (0,len(y)): return {**base,'status':'not_estimable'}
 x,terms=Xmat(d,sex,quart,spline,cr); w=num(d.pooled_weight).to_numpy(float); w=w/w.mean(); b=np.zeros(x.shape[1])
 def ll(q):
  p=expit(np.clip(x.to_numpy()@q,-35,35)); return float(np.sum(w*(y*np.log(p+1e-12)+(1-y)*np.log1p(-p+1e-12))))
 old=ll(b); ok=False
 for _ in range(160):
  p=expit(np.clip(x.to_numpy()@b,-35,35)); xx=x.to_numpy(); step=np.linalg.pinv(xx.T@((w*p*(1-p))[:,None]*xx))@(xx.T@(w*(y-p))); step=np.clip(step,-5,5); al=1
  while al>1e-8 and ll(b+al*step)<old-1e-10: al/=2
  if al<=1e-8: break
  b+=al*step; old=ll(b)
  if np.max(abs(al*step))<1e-8:ok=True;break
 xx=x.to_numpy(); p=expit(np.clip(xx@b,-35,35)); score=(w*(y-p))[:,None]*xx; bread=xx.T@((w*p*(1-p))[:,None]*xx); inv=np.linalg.pinv(bread); meat=np.zeros_like(bread)
 for _,g in d[['strata','psu']].groupby('strata',sort=False):
  ps=g.psu.unique()
  if len(ps)>1:
   z=np.vstack([score[g.index[g.psu.eq(q)],:].sum(0) for q in ps]); z-=z.mean(0); meat+=len(ps)/(len(ps)-1)*z.T@z
 cov=inv@meat@inv; cov=(cov+cov.T)/2; se=np.sqrt(np.maximum(np.diag(cov),0)); df=max(int(d.psu.nunique()-d.strata.nunique()),1); term='axis_log2:female' if sex is None and not quart and not spline else ('axis_log2' if not quart and not spline else None)
 if term is not None: ix=x.columns.tolist().index(term); stat=b[ix]/se[ix] if se[ix]>0 else np.nan; crit=t.ppf(.975,df); out={**base,'status':'ok' if ok else 'warning','beta':b[ix],'se':se[ix],'p':2*t.sf(abs(stat),df),'ci_low':b[ix]-crit*se[ix],'ci_high':b[ix]+crit*se[ix],'coef':b,'cov':cov,'names':x.columns.tolist(),'df':df}
 else: out={**base,'status':'ok' if ok else 'warning','coef':b,'cov':cov,'names':x.columns.tolist(),'df':df}
 return out
def main():
 model=load(ROOT/'work'/'scripts'/'mbzp_crc_phase2b.py','model'); reader=load(FW/'step05_crc_screen'/'run_step05_crc_screen.py','reader'); tests=pd.read_csv(FW/'step04_testset_freeze'/'unique_biomarker_test_set.csv',dtype=str,keep_default_na=False); reg=pd.read_csv(FW/'data_processed'/'detectability_registry_outcome_blinded.csv',low_memory=False); tr=tests.loc[tests.variable.eq('URXP25')].iloc[0]; exposure,src=reader.read_test_exposure(tr,reg); inv=load(FW/'outcome_inventory'/'run_outcome_sex_audit.py','inventory'); needed=src['cycles']; outs={o:load(FW/'sex_divergence_primary'/'run_primary_sex_divergence.py','prim').outcome_frame(o,needed,inv,model)[0] for o in ['obesity','myocardial_infarction','coronary_heart_disease','stroke','congestive_heart_failure','hypertension','T2D']}
 rows=[]; qrows=[]; spl=[]; spec=[]
 for o,d0 in outs.items():
  d=exposure.merge(d0,on=['SEQN','cycle'],how='inner',validate='one_to_one'); d['outcome']=d['outcome'].astype(int)
  if o in ['obesity','myocardial_infarction']:
   po=fit(d); ma=fit(d,'Male'); fe=fit(d,'Female'); rows.append({'outcome_id':o,'pooled_beta_interaction':po.get('beta',np.nan),'pooled_se':po.get('se',np.nan),'pooled_p':po.get('p',np.nan),'pooled_ci_low':po.get('ci_low',np.nan),'pooled_ci_high':po.get('ci_high',np.nan),'male_beta':ma.get('beta',np.nan),'male_se':ma.get('se',np.nan),'male_p':2*t.sf(abs(ma.get('beta',np.nan)/ma.get('se',np.nan)),ma.get('df',1)),'male_or':np.exp(ma.get('beta',np.nan)),'female_beta':fe.get('beta',np.nan),'female_se':fe.get('se',np.nan),'female_p':2*t.sf(abs(fe.get('beta',np.nan)/fe.get('se',np.nan)),fe.get('df',1)),'female_or':np.exp(fe.get('beta',np.nan)),'N':po.get('N',0),'cases':po.get('cases',0),'cycles':len(src['cycles'])})
   # Survey-weighted quartiles, locked to exposure-cycle pooled empirical cutpoints.
   d['q']=pd.qcut(d.axis_log2,4,labels=False,duplicates='drop')+1
   for sex in ['Male','Female']:
    for q in [2,3,4]:
     z=fit(d.loc[d.sex.eq(sex)],sex,quart=True); ix=z['names'].index(f'Q{q}') if f'Q{q}' in z['names'] else None; beta=z['coef'][ix] if ix is not None else np.nan; se=z['cov'][ix,ix]**.5 if ix is not None else np.nan; qrows.append({'outcome_id':o,'sex':sex,'contrast':f'Q{q}_vs_Q1','beta':beta,'OR':np.exp(beta),'SE':se,'P':2*t.sf(abs(beta/se),z['df']) if se>0 else np.nan,'N':z['N'],'cases':z['cases'],'quartile_cutpoints_log2':';'.join(map(str,sorted(d.axis_log2.quantile([.25,.5,.75]).round(8).tolist())))})
   # Spline prediction grid and joint nonlinear Wald diagnostic.
   for sex in ['Male','Female']:
    z=fit(d,sex,spline=True); grid=np.linspace(d.axis_log2.min(),d.axis_log2.max(),80); tmp=pd.DataFrame({'axis_log2':grid,'age':np.repeat(d.age.median(),80),'pir':np.repeat(d.pir.median(),80),'race':np.repeat('Non-Hispanic White',80),'smoking':np.repeat('Never',80),'cycle':np.repeat(sorted(d.cycle.unique())[0],80),'pooled_weight':1.,'psu':1,'strata':1,'outcome':0,'sex':sex,'knots':[np.quantile(d.axis_log2,x) for x in [.05,.35,.65,.95]]*80}); xx,_=Xmat(tmp,sex,spline=True); eta=xx.to_numpy()@z['coef']; var=np.einsum('ij,jk,ik->i',xx.to_numpy(),z['cov'],xx.to_numpy()); spl.extend([{'outcome_id':o,'sex':sex,'axis_log2':x,'OR':np.exp(y),'ci_low':np.exp(y-1.96*np.sqrt(max(v,0))),'ci_high':np.exp(y+1.96*np.sqrt(max(v,0))),'reference_axis_log2':float(np.median(d.axis_log2))} for x,y,v in zip(grid,eta,var)])
  else:
   po=fit(d); ma=fit(d,'Male'); fe=fit(d,'Female'); spec.append({'outcome_id':o,'male_beta':ma.get('beta',np.nan),'male_se':ma.get('se',np.nan),'female_beta':fe.get('beta',np.nan),'female_se':fe.get('se',np.nan),'interaction_beta':po.get('beta',np.nan),'interaction_p':po.get('p',np.nan)})
 # Existing robustness is a read-only recap.
 rb=pd.read_csv(FW/'sex_divergence_robustness'/'primary_robustness_uniform_results.csv'); hh=pd.read_csv(FW/'sex_divergence_robustness'/'primary_robustness_cycle_heterogeneity.csv'); rr=[]
 for o in ['obesity','myocardial_infarction']:
  z=rb[(rb.test_id=='NHANES_URXP25')&(rb.outcome_id==o)]; h=hh[(hh.test_id=='NHANES_URXP25')&(hh.outcome_id==o)].iloc[0]; rr.append({'outcome_id':o,**{f'{x}_beta':z.loc[z.sensitivity.eq(x),'beta_interaction'].iloc[0] if len(z.loc[z.sensitivity.eq(x)]) else np.nan for x in ['urinary_creatinine_sex_specific','winsorize_cycle_1_99pct','delete_cycle_upper_1pct','above_LOD_only']},'loco_n':int(z.sensitivity.str.startswith('LOCO:',na=False).sum()),'loco_min_beta':z.loc[z.sensitivity.str.startswith('LOCO:',na=False),'beta_interaction'].min(),'loco_max_beta':z.loc[z.sensitivity.str.startswith('LOCO:',na=False),'beta_interaction'].max(),'cycle_heterogeneity_p':h.p_cycle_heterogeneity})
 pd.DataFrame(rows).to_csv(OUT/'urxp25_continuous_results.csv',index=False); pd.DataFrame(qrows).to_csv(OUT/'urxp25_quartile_results.csv',index=False); pd.DataFrame(spl).to_csv(OUT/'urxp25_spline_results.csv',index=False); pd.DataFrame(spec).to_csv(OUT/'urxp25_specificity_results.csv',index=False); pd.DataFrame(rr).to_csv(OUT/'urxp25_robustness_summary.csv',index=False)
 # Two balanced double-column figures: each compares both sexes for one outcome.
 import matplotlib.pyplot as plt
 for outcome,title in [('obesity','URXP25 and obesity'),('myocardial_infarction','URXP25 and myocardial infarction')]:
  fig,axs=plt.subplots(1,2,figsize=(183/25.4,78/25.4),gridspec_kw={'width_ratios':[1.1,1]},constrained_layout=True); z=pd.DataFrame(spl).query('outcome_id==@outcome')
  for sex,c,label in [('Male',CATEGORICAL[0],'Male'),('Female',CATEGORICAL[1],'Female')]:
   a=z.query('sex==@sex'); axs[0].plot(a.axis_log2,a.OR,color=c,lw=1.5,label=label); axs[0].fill_between(a.axis_log2,a.ci_low,a.ci_high,color=c,alpha=.18); axs[1].plot(a.axis_log2,a.OR,color=c,lw=1.5,label=label)
  axs[0].axhline(1,color=GREY,ls='--',lw=.6); axs[1].axhline(1,color=GREY,ls='--',lw=.6); axs[0].set(xlabel='log2 URXP25 concentration',ylabel='Predicted odds ratio',title=title+' dose-response'); axs[1].set(xlabel='log2 URXP25 concentration',ylabel='Predicted odds ratio',title='Sex overlay'); axs[0].legend(); axs[1].legend(); save_cns_figure(fig,str(OUT/'figures'/f'urxp25_{outcome}_spline')); plt.close(fig)
 manifest={'analysis':'URXP25_FOCUSED','run_timestamp_utc':datetime.now(timezone.utc).isoformat(),'scope':'URXP25 only; obesity and MI primary; five specificity outcomes','no_external_data':True,'models_refit':'focused URXP25 only','input_hashes':{p.name:sha(p) for p in [FW/'sex_divergence_primary'/'sex_divergence_primary_406.csv',FW/'sex_divergence_robustness'/'primary_robustness_uniform_results.csv',FW/'sex_divergence_robustness'/'primary_robustness_cycle_heterogeneity.csv',FW/'sex_divergence_plan'/'SEX_DIVERGENCE_STATISTICAL_ANALYSIS_PLAN_v1.1.md']},'outputs':['urxp25_continuous_results.csv','urxp25_quartile_results.csv','urxp25_spline_results.csv','urxp25_specificity_results.csv','urxp25_robustness_summary.csv'],'figures':['figures/urxp25_obesity_spline.pdf','figures/urxp25_obesity_spline.png','figures/urxp25_myocardial_infarction_spline.pdf','figures/urxp25_myocardial_infarction_spline.png']}
 (OUT/'urxp25_focused_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');
if __name__=='__main__': main()
