"""Read-only post-primary search for transparent sex-divergent exemplars."""
from __future__ import annotations
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]; FW=ROOT/'analysis'/'disease_agnostic_environmental_framework'; OUT=Path(__file__).resolve().parent
P=FW/'sex_divergence_primary'/'sex_divergence_primary_406.csv'; E=FW/'sex_divergence_primary'/'sex_divergence_exposure_summary.csv'; S=FW/'sex_divergence_primary'/'sex_divergence_system_summary.csv'; R=FW/'sex_divergence_robustness'/'primary_robustness_uniform_results.csv'; H=FW/'sex_divergence_robustness'/'primary_robustness_cycle_heterogeneity.csv'
def sha(p):
 h=hashlib.sha256(); h.update(p.read_bytes()); return h.hexdigest()
def sg(x): return np.sign(pd.to_numeric(x,errors='coerce'))
def main():
 OUT.mkdir(parents=True,exist_ok=True); p=pd.read_csv(P); e=pd.read_csv(E); systems=pd.read_csv(S); r=pd.read_csv(R); h=pd.read_csv(H)
 p['male_z']=p.male_beta/p.male_se; p['female_z']=p.female_beta/p.female_se
 rows=[]; edges=[]; robust=[]
 for tid,g in p.groupby('test_id',sort=True):
  base=g.iloc[0]; record={'test_id':tid,'biomarker':base.biomarker,'exposure_axis':base.exposure_axis}
  for method,mc,fc in [('effect','male_beta','female_beta'),('z','male_z','female_z')]:
   mt=g.loc[g[mc].abs().idxmax()]; ft=g.loc[g[fc].abs().idxmax()]
   record.update({f'male_top_outcome_{method}':mt.outcome_id,f'male_top_system_{method}':mt.organ_system_id,f'male_top_beta_{method}':mt.male_beta,f'male_top_z_{method}':mt.male_z,f'female_top_outcome_{method}':ft.outcome_id,f'female_top_system_{method}':ft.organ_system_id,f'female_top_beta_{method}':ft.female_beta,f'female_top_z_{method}':ft.female_z,f'cross_disease_{method}':mt.outcome_id!=ft.outcome_id,f'cross_system_{method}':mt.organ_system_id!=ft.organ_system_id})
   for sex,label,x in [('male','male_top',mt),('female','female_top',ft)]:
    edges.append({'test_id':tid,'biomarker':base.biomarker,'exposure_axis':base.exposure_axis,'ranking_method':method,'destination_sex':sex,'outcome_id':x.outcome_id,'outcome_name':x.outcome_name,'organ_system_id':x.organ_system_id,'organ_system_name':x.organ_system_name,'sex_beta':x[f'{sex}_beta'],'sex_se':x[f'{sex}_se'],'sex_z':x[f'{sex}_z'],'beta_interaction':x.beta_interaction,'p_interaction':x.p_interaction,'bh_q_interaction_fixed406':x.bh_q_interaction_fixed406,'interaction_direction_supports_destination':(x.beta_interaction<0 if sex=='male' else x.beta_interaction>0)})
  rows.append(record)
 allx=pd.DataFrame(rows); edge=pd.DataFrame(edges)
 # Diagnostics are attached separately for every selected edge and ranking method.
 for _,x in edge.iterrows():
  rr=r[(r.test_id==x.test_id)&(r.outcome_id==x.outcome_id)]; primary=p[(p.test_id==x.test_id)&(p.outcome_id==x.outcome_id)].iloc[0]; pr=sg(primary.beta_interaction)
  def one(name):
   q=rr[rr.sensitivity.eq(name)]
   return q.iloc[0] if len(q) else None
  cr=one('urinary_creatinine_sex_specific'); wi=one('winsorize_cycle_1_99pct'); ta=one('delete_cycle_upper_1pct'); lo=one('above_LOD_only'); lc=rr[rr.sensitivity.str.startswith('LOCO:',na=False)]
  hh=h[(h.test_id==x.test_id)&(h.outcome_id==x.outcome_id)].iloc[0]
  good=lambda q: q is not None and q.status in ['ok','converged_with_warning']
  robust.append({**x.to_dict(),'primary_interaction_sign':pr,'creatinine_beta':cr.beta_interaction if good(cr) else np.nan,'creatinine_direction_agrees':bool(sg(cr.beta_interaction)==pr) if good(cr) else np.nan,'loco_successful_n':int(lc.status.isin(['ok','converged_with_warning']).sum()),'loco_min_beta':pd.to_numeric(lc.beta_interaction,errors='coerce').min(),'loco_max_beta':pd.to_numeric(lc.beta_interaction,errors='coerce').max(),'loco_sign_reversals':int((sg(lc.beta_interaction)*pr<0).sum()),'winsor_beta':wi.beta_interaction if good(wi) else np.nan,'tail_delete_beta':ta.beta_interaction if good(ta) else np.nan,'lod_beta':lo.beta_interaction if good(lo) else np.nan,'cycle_heterogeneity_p':hh.get('p_cycle_heterogeneity',np.nan),'cycle_heterogeneity_status':hh.status})
 rob=pd.DataFrame(robust)
 # Evidence labels use the predeclared rule.  Tiering is transparent and only
 # applies the requested cross-disease/system hierarchy to effect-ranking edges.
 shortlist=[]
 for _,a in allx.iterrows():
  ed=rob[(rob.test_id==a.test_id)&(rob.ranking_method=='effect')]; m=ed[ed.destination_sex=='male'].iloc[0]; f=ed[ed.destination_sex=='female'].iloc[0]
  formal=((m.bh_q_interaction_fixed406<.05) or (f.bh_q_interaction_fixed406<.05)) and bool(m.interaction_direction_supports_destination) and bool(f.interaction_direction_supports_destination)
  stable=lambda z: z.loco_sign_reversals==0 and (pd.isna(z.cycle_heterogeneity_p) or z.cycle_heterogeneity_p>=.05)
  if a.cross_disease_effect and a.cross_system_effect and formal and stable(m) and stable(f): tier='Tier A'; cls='Reciprocal cross-system exemplar'
  elif a.cross_disease_effect and formal: tier='Tier B'; cls='Reciprocal cross-disease candidate'
  elif formal: tier='Tier C'; cls='One-sex dominant disease landing'
  else: tier='Tier C'; cls='Unstable / ambiguous'
  shortlist.append({'test_id':a.test_id,'biomarker':a.biomarker,'exposure_axis':a.exposure_axis,'candidate_class':cls,'tier':tier,'male_top_outcome':m.outcome_name,'male_top_system':m.organ_system_name,'male_beta':m.sex_beta,'male_se':m.sex_se,'male_z':m.sex_z,'female_top_outcome':f.outcome_name,'female_top_system':f.organ_system_name,'female_beta':f.sex_beta,'female_se':f.sex_se,'female_z':f.sex_z,'interaction_for_male_top':m.beta_interaction,'q_for_male_top':m.bh_q_interaction_fixed406,'interaction_for_female_top':f.beta_interaction,'q_for_female_top':f.bh_q_interaction_fixed406,'cross_disease':a.cross_disease_effect,'cross_system':a.cross_system_effect,'main_strength':'formal interaction direction and separate robustness diagnostics shown in companion table','main_limitation':'descriptive destination ranking; no composite score or causal/mechanistic inference'})
 short=pd.DataFrame(shortlist).sort_values(['tier','test_id'])
 allx.to_csv(OUT/'exemplar_all_29_exposures.csv',index=False); edge.to_csv(OUT/'exemplar_candidate_edges.csv',index=False); rob.to_csv(OUT/'exemplar_robustness_summary.csv',index=False); short.to_csv(OUT/'exemplar_shortlist.csv',index=False)
 tier_counts=short.tier.value_counts().to_dict(); crossd=int(allx.cross_disease_effect.sum()); crosss=int(allx.cross_system_effect.sum())
 report=['# Sex-divergent disease exemplar search','', 'Read-only descriptive search of frozen primary and robustness results. No external source was accessed and no model was refit.','',f'- Exposures: 29; outcomes: 14; primary pairs: 406.',f'- Effect-ranked cross-disease exposures: {crossd}; cross-system exposures: {crosss}.',f'- Tier counts: {tier_counts}.','', 'Tier A requires an effect-ranked reciprocal cross-system split, formal interaction direction support with at least one fixed-family q<0.05, and no LOCO sign reversal/cycle-heterogeneity flag for either selected edge. This is a transparent descriptive selection rule, not a new significance score.','', '## Files','', '- `exemplar_all_29_exposures.csv`: effect- and z-ranked destination profiles.', '- `exemplar_candidate_edges.csv`: both rank-method edges.', '- `exemplar_robustness_summary.csv`: separate diagnostics for every candidate edge.', '- `exemplar_shortlist.csv`: transparent tiers; no composite ranking.', '']
 (OUT/'EXEMPLAR_SEARCH_REPORT.md').write_text('\n'.join(report),encoding='utf-8')
 (OUT/'exemplar_search_manifest.json').write_text(json.dumps({'run_timestamp_utc':datetime.now(timezone.utc).isoformat(),'inputs':{str(x.name):sha(x) for x in [P,E,S,R,H]},'script_sha256':sha(Path(__file__)),'external_data_accessed':False,'primary_or_robustness_models_refit':False,'pair_grid':406},indent=2)+'\n')
if __name__=='__main__':main()
