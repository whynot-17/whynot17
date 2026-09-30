from __future__ import annotations
from pathlib import Path
from collections import Counter, defaultdict
import csv, importlib.util, json, numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'/'tables'; FIG=ROOT/'outputs'/'figures'; LOCAL=ROOT/'outputs'/'local_only'; REPORTS=ROOT/'reports'
for d in (OUT,FIG,LOCAL,REPORTS): d.mkdir(parents=True,exist_ok=True)
HORIZON=18.0
LANDMARKS=(3,6,12)
BOOTSTRAP_REPS=1000
BOOTSTRAP_SEED=20260930

sp=importlib.util.spec_from_file_location('mrd_twin05',ROOT/'scripts'/'05_longitudinal_digital_twin.py')
base=importlib.util.module_from_spec(sp); sp.loader.exec_module(base)
mod=base.mod

def endpoint_at_l(p,landmark):
    et,fu=p['event_time'],p['followup_time']
    if et is not None and et<=landmark: return 'recurred_before_landmark',None
    if et is None and (fu is None or fu<landmark): return 'not_confirmed_at_landmark',None
    if et is not None and landmark<et<=HORIZON: return 'at_risk',1
    if et is not None and et>HORIZON: return 'at_risk',0
    if et is None and fu is not None and fu>=HORIZON: return 'at_risk',0
    return 'at_risk',None

def make_risk_rows(cohorts,landmark):
    grouped={c:[] for c in cohorts}
    for c,patients in cohorts.items():
        for p in patients:
            status,y=endpoint_at_l(p,landmark)
            hist=[(float(t),int(s)) for t,s in p['history'] if float(t)<=float(landmark) and s in (0,1)]
            max_time=max([t for t,_ in hist],default=0.0)
            assert max_time<=float(landmark)
            if status!='at_risk' or y is None or not hist: continue
            _,xc=mod.feature_row(p,landmark,False)
            _,xm=mod.feature_row(p,landmark,True)
            grouped[c].append({'patient':p,'y':int(y),'x_clin':np.asarray(xc,float),'x_mrd':np.asarray(xm,float),
                               'max_feature_mrd_time':max_time,'history_le':hist})
    return {c:v for c,v in grouped.items() if v}

def build_audit(cohorts):
    summary=[]; patient_rows=[]
    for c,patients in cohorts.items():
        recurrence_by18=sum(p['event_time'] is not None and p['event_time']<=HORIZON for p in patients)
        confirmed_nonrecurrence18=sum((p['event_time'] is not None and p['event_time']>HORIZON) or
                                      (p['event_time'] is None and p['followup_time'] is not None and p['followup_time']>=HORIZON)
                                      for p in patients)
        for L in LANDMARKS:
            status_counts=Counter(); event_count=control_count=unknown_count=hist_known_count=hist_event=hist_control=0
            for p in patients:
                status,y=endpoint_at_l(p,L)
                hist=[(t,s) for t,s in p['history'] if t<=L]
                status_counts[status]+=1
                if status=='at_risk':
                    if y is None: unknown_count+=1
                    else:
                        event_count+=y;control_count+=1-y
                        if hist:
                            hist_known_count+=1;hist_event+=y;hist_control+=1-y
                patient_rows.append({'cohort':c,'patient_id':p['patient_id'],'landmark_month':L,
                                     'status_at_landmark':status,'recurrence_by_18_from_landmark':y,
                                     'event_time':p['event_time'],'followup_proxy':p['followup_time'],
                                     'n_mrd_le_landmark':len(hist),'latest_mrd_le_landmark':hist[-1][1] if hist else None})
            model_n=hist_event+hist_control
            summary.append({'cohort':c,'raw_n':len(patients),'landmark_month':L,
                            'recurrence_by18_from_surgery':recurrence_by18,
                            'confirmed_nonrecurrence_through18_from_surgery':confirmed_nonrecurrence18,
                            'recurrence_before_or_at_landmark':status_counts['recurred_before_landmark'],
                            'not_confirmed_event_free_at_landmark':status_counts['not_confirmed_at_landmark'],
                            'at_risk_at_landmark':status_counts['at_risk'],
                            'events_after_landmark_through18':event_count,'confirmed_controls_through18':control_count,
                            'unknown_endpoint_among_landmark_risk':unknown_count,'known_outcome_with_mrd_history':hist_known_count,
                            'model_n_with_mrd':model_n,'model_events':hist_event,'model_controls':hist_control})
    return summary,patient_rows

def choose_examples_18(trajectories):
    patterns={}
    for r in trajectories:
        if len(r['history_le12'])<2: continue
        pat=r['mrd_pattern']
        if pat=='persistent_negative': patterns.setdefault('A sustained MRD-negative',[]).append(r)
        elif pat=='positive_to_negative_clearance': patterns.setdefault('B early-positive then cleared',[]).append(r)
        elif pat=='negative_to_positive_conversion': patterns.setdefault('C negative then converted positive',[]).append(r)
        elif pat=='persistent_positive': patterns.setdefault('D persistent MRD-positive',[]).append(r)
        if r['y12_18']==1 and pat=='negative_to_positive_conversion':
            patterns.setdefault('E new positive by 12m then recurrence by 18m',[]).append(r)
    selected=[]
    for label,cases in patterns.items():
        if label.startswith('A'): case=min(cases,key=lambda r:r['risk_12m'])
        elif label.startswith('B'): case=min(cases,key=lambda r:r['risk_12m']-r['risk_3m'])
        elif label.startswith('C'): case=max(cases,key=lambda r:r['risk_12m']-r['risk_3m'])
        else: case=max(cases,key=lambda r:r['risk_12m'])
        selected.append((label,case))
    return selected

def aggregate_trajectory_summary_18(trajectories):
    grouped={}
    for row in trajectories: grouped.setdefault(row['mrd_pattern'],[]).append(row)
    out=[]
    for pattern,rows in sorted(grouped.items()):
        item={'mrd_pattern':pattern,'n_with_three_landmarks':len(rows),
              'events_12_to_18':sum(r['y12_18']==1 for r in rows)}
        for col in ('risk_3m','risk_6m','risk_12m'):
            a=np.asarray([r[col] for r in rows],dtype=float)
            item[f'median_{col}']=float(np.median(a));item[f'mean_{col}']=float(np.mean(a))
        item['median_delta_risk_3_to_6']=float(np.median([r['risk_6m']-r['risk_3m'] for r in rows]))
        item['median_delta_risk_6_to_12']=float(np.median([r['risk_12m']-r['risk_6m'] for r in rows]))
        out.append(item)
    return out

def macro_metrics(fold_rows,pooled):
    def avg(field):
        vals=[float(r[field]) for r in fold_rows if np.isfinite(float(r[field]))]
        return float(np.mean(vals)) if vals else float('nan')
    out=dict(pooled)
    out['auc']=avg('auc')
    out['average_precision']=avg('average_precision')
    out['calibration_intercept']=float('nan')
    out['calibration_slope']=float('nan')
    return out

def bootstrap_loco_macro(groups,landmark,reps=BOOTSTRAP_REPS,seed=BOOTSTRAP_SEED):
    rng=np.random.default_rng(seed+int(landmark))
    store={k:defaultdict(list) for k in ('clinical','mrd','delta')}
    failures=Counter();cohorts=sorted(groups)
    for b in range(reps):
        sampled={}
        for c in cohorts:
            source=groups[c]
            draw=[]
            for y in (0,1):
                stratum=[r for r in source if r['y']==y]
                idx=rng.integers(0,len(stratum),size=len(stratum))
                draw.extend(stratum[int(i)] for i in idx)
            rng.shuffle(draw);sampled[c]=draw
        try:
            pc,fc,mc=base.loco(sampled,'x_clin')
            pm,fm,mm=base.loco(sampled,'x_mrd')
            if [r['outcome'] for r in pc] != [r['outcome'] for r in pm]:
                raise RuntimeError('paired bootstrap outcomes differ')
            vc=macro_metrics(fc,mc);vm=macro_metrics(fm,mm)
            for metric in ('auc','average_precision','brier'):
                store['clinical'][metric].append(vc[metric])
                store['mrd'][metric].append(vm[metric])
                store['delta'][metric].append(vm[metric]-vc[metric])
        except Exception as exc:
            failures[type(exc).__name__+': '+str(exc)[:100]]+=1
        if (b+1)%100==0:
            print(f'L{landmark} stratified bootstrap {b+1}/{reps}; valid={len(store["clinical"]["auc"])}; failed={sum(failures.values())}',flush=True)
    return store,failures

def main():
    cohort_patients={'Chen':mod.read_chen(),'ColonAiQ':mod.read_coloniaiq(),'COSMOS':mod.read_cosmos()}
    audit,patient_audit=build_audit(cohort_patients)
    base.csv_write(OUT/'l18_cohort_outcome_audit.csv',audit)
    base.csv_write(LOCAL/'l18_patient_outcome_audit.csv',patient_audit)

    expected={
      3:{'Chen':(155,25,130),'ColonAiQ':(70,40,30),'COSMOS':(326,23,303)},
      6:{'Chen':(153,23,130),'ColonAiQ':(56,26,30),'COSMOS':(322,19,303)},
      12:{'Chen':(140,10,130),'ColonAiQ':(55,16,39),'COSMOS':(307,4,303)}}
    groups_by_l={}
    for L in LANDMARKS:
        groups=make_risk_rows(cohort_patients,L)
        if set(groups)!={'Chen','ColonAiQ','COSMOS'}:
            raise RuntimeError(f'L{L}: not all three cohorts pass outcome/MRD ascertainment')
        for c,g in groups.items():
            got=(len(g),sum(r['y'] for r in g),len(g)-sum(r['y'] for r in g))
            if got!=expected[L][c]: raise RuntimeError(f'L{L} {c} counts changed: {got} != {expected[L][c]}')
            assert all(r['max_feature_mrd_time']<=L for r in g)
        train_events={held:sum(r['y'] for c,g in groups.items() if c!=held for r in g) for held in groups}
        if min(train_events.values())<10:
            raise RuntimeError(f'L{L}: a LOCO training fold has fewer than 10 events: {train_events}')
        groups_by_l[L]=groups

    metrics_rows=[];fold_rows=[];patient_predictions=[];predictions={};point_cache={};bootstrap_rows=[]
    for L,groups in groups_by_l.items():
        label=f'L{L}'
        pc,fc,mc_raw=base.loco(groups,'x_clin')
        pm,fm,mm_raw=base.loco(groups,'x_mrd')
        mc=macro_metrics(fc,mc_raw);mm=macro_metrics(fm,mm_raw)
        if [x['outcome'] for x in pc] != [x['outcome'] for x in pm]:
            raise RuntimeError(f'{label}: clinical and MRD models used different outcomes')
        for name,metric,folds,plist in [('clinical_only',mc,fc,pc),('clinical_plus_MRD',mm,fm,pm)]:
            metrics_rows.append({'landmark':label,'model':name,**metric,
                                 'features':'stage_III_IV;rectal_primary' if name=='clinical_only' else 'stage_III_IV;rectal_primary;current_mrd_positive;prior_mrd_positive_before_current',
                                 'cohorts':'Chen;ColonAiQ;COSMOS','endpoint':'recurrence from landmark through month 18',
                                 'ridge_lambda':base.RIDGE,'validation':'leave-one-cohort-out'})
            for r in folds: fold_rows.append({'landmark':label,'model':name,**r})
            for r in plist: patient_predictions.append({'landmark':label,'model':name,**r})
        predictions[L]={'clinical':pc,'mrd':pm}
        point_cache[L]={'clinical':mc,'mrd':mm}
    base.csv_write(OUT/'l18_loco_metrics.csv',metrics_rows)
    base.csv_write(OUT/'l18_loco_by_heldout_cohort.csv',fold_rows)
    base.csv_write(LOCAL/'l18_patient_predictions.csv',patient_predictions)

    for L,groups in groups_by_l.items():
        store,failures=bootstrap_loco_macro(groups,L,BOOTSTRAP_REPS,BOOTSTRAP_SEED)
        bootstrap_rows.extend(base.ci_rows_for_landmark(f'L{L}',point_cache[L]['clinical'],point_cache[L]['mrd'],
                                                         store,failures,BOOTSTRAP_REPS))
    base.csv_write(OUT/'l18_landmark_bootstrap_ci.csv',bootstrap_rows)

    pred_maps={L:{(r['cohort'],r['patient_id']):r for r in predictions[L]['mrd']} for L in LANDMARKS}
    trajectories=[]
    for c,records in groups_by_l[12].items():
        for rec in records:
            p=rec['patient'];key=(c,p['patient_id'])
            if not all(key in pred_maps[L] for L in LANDMARKS): continue
            hist=[(t,s) for t,s in p['history'] if t<=12]
            assert max([t for t,s in hist],default=0)<=12
            trajectories.append({'cohort':c,'patient_id':p['patient_id'],'event_time':p['event_time'],
                                 'followup_proxy':p['followup_time'],'y12_18':rec['y'],
                                 'mrd_pattern':base.state_pattern(hist),'history_le12':hist,
                                 'risk_3m':pred_maps[3][key]['predicted_risk'],
                                 'risk_6m':pred_maps[6][key]['predicted_risk'],
                                 'risk_12m':pred_maps[12][key]['predicted_risk']})
    traj_summary=aggregate_trajectory_summary_18(trajectories)
    base.csv_write(OUT/'l18_trajectory_summary.csv',traj_summary)
    details=[]
    for r in trajectories:
        details.append({**{k:v for k,v in r.items() if k!='history_le12'},
                        'mrd_history_le12':';'.join(f'{t:g}:{s}' for t,s in r['history_le12'])})
    base.csv_write(LOCAL/'l18_patient_risk_trajectories.csv',details)
    selected=choose_examples_18(trajectories)
    case_rows=[{'pattern':label,'cohort':r['cohort'],'patient_id':r['patient_id'],'event_time':r['event_time'],
                'outcome_12_18':r['y12_18'],'risk_3m':r['risk_3m'],'risk_6m':r['risk_6m'],'risk_12m':r['risk_12m'],
                'mrd_history_le12':';'.join(f'{t:g}:{s}' for t,s in r['history_le12'])} for label,r in selected]
    base.csv_write(LOCAL/'l18_typical_cases_local.csv',case_rows)
    base.svg_trajectories(selected,LOCAL/'l18_patient_trajectories.svg')
    base.svg_metric_plot(bootstrap_rows,'auc','18个月复发终点：LOCO AUC',FIG/'l18_auc_comparison.svg',0.4,1.0)
    base.svg_metric_plot(bootstrap_rows,'brier','18个月复发终点：LOCO Brier score',FIG/'l18_brier_comparison.svg',0.0,0.25)

    metric_lookup={(r['landmark'],r['model']):r for r in metrics_rows}
    ci={(r['landmark'],r['contrast'],r['metric']):r for r in bootstrap_rows}
    main_rows=[]
    for L in LANDMARKS:
        lm=f'L{L}';a=metric_lookup[(lm,'clinical_only')];b=metric_lookup[(lm,'clinical_plus_MRD')]
        def fm(kind,metric):
            r=ci[(lm,kind,metric)]
            return f'{r["estimate"]:.3f} [{r["ci_95_lower"]:.3f}, {r["ci_95_upper"]:.3f}]'
        main_rows.append(f'| L{L} | {a["n"]} | {a["events"]} | {a["auc"]:.3f} | {b["auc"]:.3f} [{ci[(lm,"clinical_plus_MRD","auc")]["ci_95_lower"]:.3f}, {ci[(lm,"clinical_plus_MRD","auc")]["ci_95_upper"]:.3f}] | {fm("MRD_minus_clinical","auc")} | {a["average_precision"]:.3f} | {b["average_precision"]:.3f} | {fm("MRD_minus_clinical","average_precision")} | {a["brier"]:.3f} | {b["brier"]:.3f} | {fm("MRD_minus_clinical","brier")} |')
    fold_md='\n'.join(f'| {r["landmark"]} | {r["model"]} | {r["held_out_cohort"]} | {r["n"]} | {r["events"]} | {r["auc"]:.3f} | {r["average_precision"]:.3f} | {r["brier"]:.3f} | {r["calibration_intercept"]:.3f} | {r["calibration_slope"]:.3f} |' for r in fold_rows)
    audit_md='\n'.join(f'| {r["landmark_month"]} | {r["cohort"]} | {r["model_n_with_mrd"]} | {r["model_events"]} | {r["model_controls"]} | {r["unknown_endpoint_among_landmark_risk"]} |' for r in audit)
    traj_md='\n'.join(f'| {r["mrd_pattern"]} | {r["n_with_three_landmarks"]} | {r["events_12_to_18"]} | {r["median_risk_3m"]:.3f} | {r["median_risk_6m"]:.3f} | {r["median_risk_12m"]:.3f} |' for r in traj_summary)
    report=f'''# 三队列 CRC MRD：18个月复发动态预测

**日期：**2026-09-30  
**队列：**Chen、ColonAiQ、COSMOS；Tie 不纳入。

## 18个月结局和 ColonAiQ 可用性

终点定义为每个 landmark 后、术后18个月内临床复发。阳性为 recurrence time ∈ (landmark, 18]；阴性为 recurrence time >18，或无复发记录且逐患者 follow-up proxy ≥18个月。仅用 landmark 当时及此前 MRD，缺失检测不作阴性填补。

| Landmark | 队列 | 有 MRD 历史且结局可判定 N | 事件 | 对照 | 风险集中18月结局未知 |
|---:|---|---:|---:|---:|---:|
{audit_md}

ColonAiQ 在三个 landmark 均有阳性、阴性病例可用：L3 为 70 人（40/30），L6 为 56 人（26/30），L12 为 55 人（16/39）。因此可进入三队列18个月模型。ColonAiQ 的阴性观察证据使用补充表中末次术后采血时间作为随访代理，并非逐患者影像学随访时间，解释时需考虑该限制。

三队列 L12 风险集的 MRD 可判定样本合计 502 人、30 例事件。Chen 为 140/10，ColonAiQ 为 55/16，COSMOS 为 307/4。每个 LOCO 训练折均至少有10例事件；但 COSMOS 留出折只有4例事件，L12 的 COSMOS 特异性能估计会很不稳定，故该部分为探索性。未因事件少而事后移除 COSMOS。

## 预测设置与无信息泄漏

- L3 使用 ≤3月 MRD，预测 (3,18] 月复发；L6 使用 ≤6月 MRD，预测 (6,18] 月复发；L12 使用 ≤12月 MRD，预测 (12,18] 月复发。
- 两个模型分别为 clinical-only（stage、site）与 clinical + MRD（再加入最近一次 MRD 状态、当前样本前曾阳性）。
- ridge logistic，λ=1.0；不加入 cohort 指示变量，以便留出队列预测。
- 所有患者逐一断言 MRD 特征时间不超过对应 landmark。

## LOCO 表现与完整 bootstrap 95% CI

在每个队列×结局类别层内以患者为单位有放回抽样，每次重新拟合所有 LOCO 折和两类模型。AUC/AP 为三个留出队列折内指标的宏平均，Brier 为全部 out-of-fold 预测的汇总均值；每个 landmark 1000次，随机种子 20260930。差值为 clinical + MRD 减 clinical-only。

| Landmark | N | 事件 | AUC clinical | AUC + MRD [95% CI] | ΔAUC [95% CI] | AP clinical | AP + MRD | ΔAP [95% CI] | Brier clinical | Brier + MRD | ΔBrier [95% CI] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(main_rows)}

## 各留出队列表现

| Landmark | 模型 | 留出队列 | N | 事件 | AUC | AP | Brier | 校准截距 | 校准斜率 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
{fold_md}

## 动态风险轨迹

同时具有 L3、L6、L12 cross-fitted 风险预测的 L12 患者数为 {len(trajectories)}。三个时点均预测到同一个术后18个月终点，风险随新 MRD 信息更新。

| MRD 状态模式 | 人数 | (12,18]月事件 | 中位 Risk3 | 中位 Risk6 | 中位 Risk12 |
|---|---:|---:|---:|---:|---:|
{traj_md}

患者级明细、预测、典型病例和患者轨迹图仅保存在 outputs/local_only/，未提交至 GitHub。

## 结果解释

该模型可在三个公开队列中计算，ColonAiQ 在18个月终点下可纳入。但 ColonAiQ 的随访判定来自末次采血时间代理，且 L12 总事件只有30例、COSMOS 留出折仅4例事件。结果属于多队列开发和 internal-external LOCO 评估，不是独立外部验证或临床效用证明。

## 输出

- outputs/tables/l18_cohort_outcome_audit.csv
- outputs/tables/l18_loco_metrics.csv
- outputs/tables/l18_loco_by_heldout_cohort.csv
- outputs/tables/l18_landmark_bootstrap_ci.csv
- outputs/tables/l18_trajectory_summary.csv
- outputs/figures/l18_auc_comparison.svg
- outputs/figures/l18_brier_comparison.svg
- 患者级输出保存在本地 outputs/local_only/
'''
    (REPORTS/'18m_landmark_multicohort_zh.md').write_text(report,encoding='utf-8')
    print(json.dumps({'landmark_n':{str(L):sum(len(g) for g in groups_by_l[L].values()) for L in LANDMARKS},
                      'landmark_events':{str(L):sum(r['y'] for g in groups_by_l[L].values() for r in g) for L in LANDMARKS},
                      'l12_by_cohort':{c:{'n':len(g),'events':sum(r['y'] for r in g)} for c,g in groups_by_l[12].items()},
                      'bootstrap_reps':BOOTSTRAP_REPS,'bootstrap_failures':{str(L):max(r['bootstrap_reps_failed'] for r in bootstrap_rows if r['landmark']==f'L{L}') for L in LANDMARKS},
                      'trajectory_n':len(trajectories),'typical_cases':len(selected),'report':str(REPORTS/'18m_landmark_multicohort_zh.md')},ensure_ascii=False,indent=2))

if __name__=='__main__': main()
