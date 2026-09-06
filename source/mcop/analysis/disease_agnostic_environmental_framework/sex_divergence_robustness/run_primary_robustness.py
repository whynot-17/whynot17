"""Execute the locked, uniform Primary Robustness v1.0 package."""
from __future__ import annotations
import importlib.util, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import t, chi2

ROOT=Path(__file__).resolve().parents[3]; FW=ROOT/'analysis'/'disease_agnostic_environmental_framework'
PRIMARY=FW/'sex_divergence_primary'; OUT=Path(__file__).resolve().parent; DATA=ROOT/'work'/'nhanes_phase2a'/'data'
def load(p,n):
 s=importlib.util.spec_from_file_location(n,p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def n(s): return pd.to_numeric(s,errors='coerce')

def design(d, cr=False, hetero=False):
 x=pd.DataFrame({'Intercept':1.,'axis_log2':n(d.axis_log2)},index=d.index); a=n(d.age)-50
 x['age_centered']=a; x['age_centered_sq']=a*a; x['female']=d.sex.eq('Female').astype(float); x['axis_log2:female']=x.axis_log2*x.female; x['pir']=n(d.pir)
 for c,lv in [('race',['Mexican American','Other Hispanic','Non-Hispanic Black','Other/Multi']),('smoking',['Former','Current'])]:
  for z in lv:x[f'{c}={z}']=d[c].eq(z).astype(float)
 cycles=sorted(d.cycle.unique().tolist())
 for z in cycles[1:]: x[f'cycle={z}']=d.cycle.eq(z).astype(float)
 if cr:
  x['creatinine_log2']=n(d.creatinine_log2); x['creatinine_log2:female']=x.creatinine_log2*x.female
 if hetero:
  for z in cycles[1:]: x[f'axis_log2:female:cycle={z}']=x['axis_log2:female']*d.cycle.eq(z).astype(float)
 return x.to_numpy(float),x.columns.tolist()

def fit(d,cr=False,hetero=False):
 req=['outcome','axis_log2','age','pir','race','smoking','pooled_weight','psu','strata','sex','cycle']+(['creatinine_log2'] if cr else [])
 d=d.dropna(subset=req).loc[lambda z:z.pooled_weight.gt(0)].reset_index(drop=True); y=n(d.outcome).to_numpy(float); base={'N':len(d),'cases':int(y.sum()) if len(y) else 0,'controls':int(len(y)-y.sum()) if len(y) else 0}
 if len(d)==0 or y.sum() in (0,len(y)) or d.sex.nunique()!=2:return {**base,'status':'not_estimable','reason':'no complete-case outcome/sex variation'}
 x,names=design(d,cr,hetero); w=n(d.pooled_weight).to_numpy(float); w=w/w.mean(); b=np.zeros(x.shape[1])
 def ll(q):
  p=expit(np.clip(x@q,-35,35)); return float(np.sum(w*(y*np.log(p+1e-12)+(1-y)*np.log1p(-p+1e-12))))
 old=ll(b); ok=False
 for _ in range(200):
  p=expit(np.clip(x@b,-35,35)); step=np.linalg.pinv(x.T@((w*p*(1-p))[:,None]*x))@(x.T@(w*(y-p))); step=np.clip(step,-5,5); al=1.
  while al>=1e-8 and ll(b+al*step)<old-1e-10:al/=2
  if al<1e-8:break
  b+=al*step; old=ll(b)
  if np.max(abs(al*step))<1e-8:ok=True;break
 p=expit(np.clip(x@b,-35,35)); scores=(w*(y-p))[:,None]*x; bread=x.T@((w*p*(1-p))[:,None]*x); inv=np.linalg.pinv(bread); meat=np.zeros_like(bread)
 for _,g in d[['strata','psu']].groupby('strata',sort=False):
  ps=g.psu.unique()
  if len(ps)>1:
   z=np.vstack([scores[g.index[g.psu.eq(q)],:].sum(axis=0) for q in ps]); z-=z.mean(axis=0); meat+=len(ps)/(len(ps)-1)*z.T@z
 cov=inv@meat@inv; cov=(cov+cov.T)/2; se=np.sqrt(np.maximum(np.diag(cov),0)); df=max(int(d.psu.nunique()-d.strata.nunique()),1); ix=names.index('axis_log2:female')
 if not np.isfinite(se[ix]) or se[ix]<=0:return {**base,'status':'fit_failed','reason':'interaction variance unavailable'}
 z=b[ix]/se[ix]; out={**base,'status':'ok' if ok else 'converged_with_warning','reason':'','beta_interaction':float(b[ix]),'se_interaction':float(se[ix]),'p_interaction':float(2*t.sf(abs(z),df)),'design_df':df,'coef':dict(zip(names,b)),'cov':cov,'names':names}
 return out

def main():
 primary=load(PRIMARY/'run_primary_sex_divergence.py','primary'); reader=load(FW/'step05_crc_screen'/'run_step05_crc_screen.py','reader'); inv=load(FW/'outcome_inventory'/'run_outcome_sex_audit.py','inventory'); model=load(ROOT/'work'/'scripts'/'mbzp_crc_phase2b.py','model'); model.DATA_DIR=DATA
 tests=pd.read_csv(FW/'step04_testset_freeze'/'unique_biomarker_test_set.csv',dtype=str,keep_default_na=False); reg=pd.read_csv(FW/'data_processed'/'detectability_registry_outcome_blinded.csv',low_memory=False); frozen=pd.read_csv(FW/'outcome_inventory'/'frozen_outcome_set_v1.csv'); outcomes=frozen.loc[frozen.selection_for_followup.eq(True),'outcome_id'].tolist()
 cache={str(r.test_id):reader.read_test_exposure(r,reg) for _,r in tests.iterrows()}; cycles=sorted({c for _,s in cache.values() for c in s.get('cycles',[]) }); outs={o:primary.outcome_frame(o,cycles,inv,model)[0] for o in outcomes}; rows=[]; hetero=[]
 for _,test in tests.iterrows():
  e,src=cache[str(test.test_id)]
  for o in outcomes:
   base={'test_id':str(test.test_id),'outcome_id':o,'matrix':str(test.matrix),'cycles':';'.join(src.get('cycles',[]))}
   if e.empty: rows.append({**base,'sensitivity':'source','status':'not_estimable','reason':src.get('reason','empty exposure')}); continue
   d=e.merge(outs[o].loc[lambda q:q.cycle.isin(src['cycles'])],on=['SEQN','cycle'],how='inner',validate='one_to_one')
   def add(label,x):
    r=fit(x); rows.append({**base,'sensitivity':label,**{k:v for k,v in r.items() if k not in {'coef','cov','names'}}})
   if 'urine' in str(test.matrix).lower():
    parts=[]
    for c in src['cycles']:
     p=DATA/f'{c}_ALB_CR.XPT'; z=pd.read_sas(p,format='xport',encoding='latin1')[['SEQN','URXUCR']]; z['cycle']=c; z['creatinine_log2']=np.log2(n(z.URXUCR).where(n(z.URXUCR)>0)); parts.append(z[['SEQN','cycle','creatinine_log2']])
    add('urinary_creatinine_sex_specific',d.merge(pd.concat(parts),on=['SEQN','cycle'],how='left',validate='one_to_one'))
   else: rows.append({**base,'sensitivity':'urinary_creatinine_sex_specific','status':'not_applicable','reason':'non-urine exposure'})
   for c in src['cycles']: add(f'LOCO:{c}',d.loc[d.cycle.ne(c)])
   w=d.copy(); w['axis_log2']=w.groupby('cycle').axis_log2.transform(lambda q:q.clip(q.quantile(.01),q.quantile(.99))); add('winsorize_cycle_1_99pct',w)
   cut=d.groupby('cycle').axis_log2.transform(lambda q:q.quantile(.99)); add('delete_cycle_upper_1pct',d.loc[d.axis_log2.le(cut)])
   if 'above_lod' in d: add('above_LOD_only',d.loc[d.above_lod])
   else: rows.append({**base,'sensitivity':'above_LOD_only','status':'not_applicable','reason':'no usable LOD flag'})
   if len(src['cycles'])>=3:
    r=fit(d,hetero=True); terms=[i for i,z in enumerate(r.get('names',[])) if z.startswith('axis_log2:female:cycle=')]
    if r.get('status') in {'ok','converged_with_warning'} and terms:
     b=np.array([r['coef'][r['names'][i]] for i in terms]); v=r['cov'][np.ix_(terms,terms)]; stat=float(b@np.linalg.pinv(v)@b); p=float(chi2.sf(stat,len(terms))); hetero.append({**base,'status':r['status'],'n_cycle_deviations':len(terms),'wald_statistic':stat,'p_cycle_heterogeneity':p,'N':r['N'],'cases':r['cases']})
    else: hetero.append({**base,'status':r.get('status'),'reason':r.get('reason','')})
   else: hetero.append({**base,'status':'not_applicable','reason':'fewer than three cycles'})
 pd.DataFrame(rows).to_csv(OUT/'primary_robustness_uniform_results.csv',index=False); pd.DataFrame(hetero).to_csv(OUT/'primary_robustness_cycle_heterogeneity.csv',index=False)
 (OUT/'primary_robustness_manifest.json').write_text(json.dumps({'analysis':'PRIMARY_ROBUSTNESS_V1','run_timestamp_utc':datetime.now(timezone.utc).isoformat(),'pairs':406,'uniform_rows':len(rows),'heterogeneity_rows':len(hetero),'primary_results_used':'sex_divergence_primary_406.csv','interpretation_not_run':True},indent=2)+'\n')
if __name__=='__main__':main()
