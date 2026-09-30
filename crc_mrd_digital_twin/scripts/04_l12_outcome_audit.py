from __future__ import annotations
from pathlib import Path
import csv, importlib.util, json, math

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'/'tables'; REPORTS=ROOT/'reports'; LOCAL=ROOT/'outputs'/'local_only'
OUT.mkdir(parents=True,exist_ok=True); REPORTS.mkdir(parents=True,exist_ok=True); LOCAL.mkdir(parents=True,exist_ok=True)
mod_path=ROOT/'scripts'/'03_multicohort_landmark_model.py'
spec=importlib.util.spec_from_file_location('mrd_model3',mod_path)
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
cohorts=[mod.read_chen(),mod.read_coloniaiq(),mod.read_cosmos()]
rows=[]; patient_rows=[]

def get_y24(p):
    et=p['event_time']; fu=p['followup_time']
    # A recurrence after month 12 establishes the patient reached the 12-month landmark event-free.
    # For event-free records, require observed follow-up through month 12 to establish landmark eligibility.
    if et is not None and et <= 12:
        return 'recurred_by_12', None
    if et is None and (fu is None or fu < 12):
        return 'unknown_before_12', None
    if et is not None and 12 < et <= 24:
        return 'at_risk_12', 1
    if et is not None and et > 24:
        return 'at_risk_12', 0
    if et is None and fu is not None and fu >= 24:
        return 'at_risk_12', 0
    return 'at_risk_12', None

def trajectory(hist):
    if not hist: return 'no_history'
    if len(hist)==1: return 'single_positive' if hist[0][1] else 'single_negative'
    first,last=hist[0][1],hist[-1][1]
    return { (0,0):'negative_to_negative',(0,1):'negative_to_positive',(1,0):'positive_to_negative',(1,1):'positive_to_positive'}[(first,last)]

for cohort in cohorts:
    name=cohort[0]['cohort'] if cohort else 'unknown'
    eligible=[]; all_at_risk=[]
    for p in cohort:
        status,y24=get_y24(p)
        hist=[(t,s) for t,s in p['history'] if t<=12]
        # Fail closed if any future MRD is accidentally passed to the state definition.
        assert max([t for t,s in hist],default=0.0)<=12.0
        at_risk=status=='at_risk_12'
        model_eligible=at_risk and y24 is not None and bool(hist)
        p_row={'cohort':name,'patient_id':p['patient_id'],'landmark12_status':status,'y12_24':y24,
               'event_time':p['event_time'],'followup_time':p['followup_time'],
               'n_mrd_le12':len(hist),'latest_mrd_le12':hist[-1][1] if hist else None,
               'prior_positive_before_current':int(any(s==1 for t,s in hist[:-1])) if hist else None,
               'transition_first_to_last':trajectory(hist)}
        patient_rows.append(p_row)
        if at_risk: all_at_risk.append(p_row)
        if model_eligible: eligible.append(p_row)
    statuses={x:sum(p['landmark12_status']==x for p in patient_rows if p['cohort']==name) for x in ['recurred_by_12','at_risk_12','unknown_before_12']}
    total=[p for p in patient_rows if p['cohort']==name]
    model=[p for p in eligible]
    row={
      'cohort':name,'raw_or_reconstructed_n':len(cohort),
      'recurrence_by_12m':statuses['recurred_by_12'],
      'confirmed_event_free_at_12m':statuses['at_risk_12'],
      '12m_status_unconfirmed_before_12m':statuses['unknown_before_12'],
      'recurrence_12_to_24m':sum(p['event_time'] is not None and 12<p['event_time']<=24 for p in cohort),
      'recurrence_after_24m':sum(p['event_time'] is not None and p['event_time']>24 for p in cohort),
      'confirmed_nonrecurrence_followup_ge24m':sum(p['event_time'] is None and p['followup_time'] is not None and p['followup_time']>=24 for p in cohort),
      'unknown_24m_outcome_among_12m_at_risk':sum(p['landmark12_status']=='at_risk_12' and p['y12_24'] is None for p in total),
      'mrd_history_le12_among_12m_at_risk':sum(p['landmark12_status']=='at_risk_12' and p['n_mrd_le12']>0 for p in total),
      'l12_model_n':len(model),'l12_model_events':sum(p['y12_24']==1 for p in model),
      'l12_latest_positive':sum(p['latest_mrd_le12']==1 for p in model),
      'l12_latest_negative':sum(p['latest_mrd_le12']==0 for p in model),
      'l12_prior_positive_yes':sum(p['prior_positive_before_current']==1 for p in model),
      'l12_prior_positive_no':sum(p['prior_positive_before_current']==0 for p in model),
      'l12_transition_single_positive':sum(p['transition_first_to_last']=='single_positive' for p in model),
      'l12_transition_single_negative':sum(p['transition_first_to_last']=='single_negative' for p in model),
      'l12_transition_negative_to_negative':sum(p['transition_first_to_last']=='negative_to_negative' for p in model),
      'l12_transition_negative_to_positive':sum(p['transition_first_to_last']=='negative_to_positive' for p in model),
      'l12_transition_positive_to_negative':sum(p['transition_first_to_last']=='positive_to_negative' for p in model),
      'l12_transition_positive_to_positive':sum(p['transition_first_to_last']=='positive_to_positive' for p in model),
    }
    row['l12_model_controls']=row['l12_model_n']-row['l12_model_events']
    row['l12_cohort_include']=int(row['l12_model_events']>=10 and row['l12_model_controls']>=10)
    row['l12_exclusion_reason']='' if row['l12_cohort_include'] else '可判定结局中事件或阴性对照少于10；不进入 L12 建模/LOCO'
    row['l12_primary_model_n']=row['l12_model_n'] if row['l12_cohort_include'] else 0
    row['l12_primary_model_events']=row['l12_model_events'] if row['l12_cohort_include'] else 0
    rows.append(row)

with (OUT/'l12_cohort_audit.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with (LOCAL/'l12_patient_outcome_audit.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=list(patient_rows[0]));w.writeheader();w.writerows(patient_rows)

aggregate={'recurrence_12_to_24m':sum(r['recurrence_12_to_24m'] for r in rows),
          'l12_model_n':sum(r['l12_primary_model_n'] for r in rows),
          'l12_model_events':sum(r['l12_primary_model_events'] for r in rows)}
cols=['cohort','raw_or_reconstructed_n','recurrence_by_12m','confirmed_event_free_at_12m','recurrence_12_to_24m','confirmed_nonrecurrence_followup_ge24m','recurrence_after_24m','unknown_24m_outcome_among_12m_at_risk','mrd_history_le12_among_12m_at_risk','l12_model_n','l12_model_events','l12_model_controls','l12_cohort_include','l12_exclusion_reason']
body='\n'.join('| '+ ' | '.join(str(r[k]) for k in cols)+' |' for r in rows)
metric='| 队列 | 原分析集N | ≤12月复发 | 确认到12月仍无复发 | (12,24]复发 | 确认未复发且观察≥24月 | >24月复发 | 24月结局未知（12月风险人群） | 12月前MRD历史（12月风险人群） | 候选L12 N | 事件 | 对照 | 纳入L12主模型 | 排除原因 |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n'+body

events_by_cohort={r['cohort']:(r['l12_primary_model_events'] if r['l12_cohort_include'] else '排除') for r in rows}
report=f'''# 12月 landmark → 24月复发结局可判定性审计

**日期：**2026-09-30  
**主分析队列：**Chen、ColonAiQ、COSMOS（Tie 已移出主分析）。

## 结局与风险集定义

- 先要求患者可被证明在术后12个月时仍未发生复发：临床复发时间必须大于12个月；若最终无复发记录，则需有至少12个月观察证据。
- L12 阳性定义为 `12 < recurrence_time ≤ 24`。
- L12 阴性定义为复发时间 `>24`，或最终未复发且逐患者观察达到至少24个月。
- 最终无复发但随访不足24个月者，24月结局未知，不编码为阴性。12个月状态本身不能确认者不进入风险集。
- 模型候选人还须有至少一条可判读且采样时间 `≤12` 月的术后 MRD 历史。
- 结局/随访证据沿用公开逐患者临床事件时间和可核实随访代理；Chen 与 ColonAiQ 的公开末次采血时间并不等同于影像学随访时间，阴性判定较为保守但仍有测量局限。

## 队列审计

{metric}

合计可建模人数 **{aggregate['l12_model_n']}**，12–24月事件 **{aggregate['l12_model_events']}**。L12纳入队列的事件数：{', '.join(f"{k} {v}" for k,v in events_by_cohort.items())}。

## MRD 特征核查

`l12_cohort_audit.csv` 报告 L12 建模人群截至12月的最新 MRD 阳性/阴性、当前样本之前曾阳性、以及首末 MRD 状态转换分布。特征提取强制过滤 `sample_time ≤ 12`，并逐患者运行上界断言。24月及更晚 MRD 未参与这些状态或模型候选集。

## 是否进入模型

按队列门槛纳入 Chen 与 COSMOS；ColonAiQ 的24月阴性对照为0（123名12月风险患者结局未知），故整个队列排除。纳入队列合计 L12 可建模 **{aggregate['l12_model_n']}人、{aggregate['l12_model_events']}例事件**。后续使用低维 ridge logistic，并仅作为 exploratory analysis；LOCO 的留出队列事件数和区分度将单独报告。若某一 LOCO 训练折没有足够事件、或结果不可识别，则保留审计表并明确不报告可推广的模型性能。该模型不会被表述为外部验证或临床工具。

## 输出

- 队列级审计：`outputs/tables/l12_cohort_audit.csv`
- 患者级审计明细（本地，不提交版本库）：`outputs/local_only/l12_patient_outcome_audit.csv`
'''
(REPORTS/'12m_landmark_24m_outcome_audit_zh.md').write_text(report,encoding='utf-8')
print(json.dumps({'audit':rows,'aggregate':aggregate,'report':str(REPORTS/'12m_landmark_24m_outcome_audit_zh.md')},ensure_ascii=False,indent=2))
