from __future__ import annotations
from pathlib import Path
from collections import Counter, defaultdict
import csv, importlib.util, html, json, math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'/'tables'; FIG=ROOT/'outputs'/'figures'; LOCAL=ROOT/'outputs'/'local_only'; REPORTS=ROOT/'reports'
for d in (OUT,FIG,LOCAL,REPORTS): d.mkdir(parents=True,exist_ok=True)
BOOTSTRAP_REPS=1000
BOOTSTRAP_SEED=20260930
RIDGE=1.0

spec=importlib.util.spec_from_file_location('mrd_model3',ROOT/'scripts'/'03_multicohort_landmark_model.py')
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)


def csv_write(path,rows):
    if not rows: return
    with Path(path).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def all_cohort_patients():
    return {'Chen':mod.read_chen(),'ColonAiQ':mod.read_coloniaiq(),'COSMOS':mod.read_cosmos()}


def make_risk_rows(cohorts, landmark):
    grouped={c:[] for c in cohorts}
    for c,patients in cohorts.items():
        for p in patients:
            et,fu=p['event_time'],p['followup_time']
            hist=[(float(t),int(s)) for t,s in p['history'] if float(t)<=float(landmark) and s in (0,1)]
            # Mandatory leakage guard: the feature history cannot cross its landmark.
            max_time=max([t for t,_ in hist],default=0.0)
            assert max_time<=float(landmark),f"{c}/{p['patient_id']}: MRD time {max_time} exceeds L{landmark}"
            if landmark<12:
                # Y12 must be known, the patient must be event-free at L, and have an observed MRD history.
                if p['y12'] is None or (et is not None and et<=landmark) or not hist: continue
                y=int(et is not None and landmark<et<=12)
            else:
                # Confirmed event-free through 12 months.
                if et is not None and et<=12: continue
                if et is None and (fu is None or fu<12): continue
                # Y24 must be known; otherwise keep out of the model.
                if et is not None and 12<et<=24: y=1
                elif et is not None and et>24: y=0
                elif et is None and fu is not None and fu>=24: y=0
                else: continue
                if not hist: continue
            n1,x1=mod.feature_row(p,landmark,False)
            n2,x2=mod.feature_row(p,landmark,True)
            grouped[c].append({'patient':p,'y':int(y),'x_clin':np.asarray(x1,float),'x_mrd':np.asarray(x2,float),
                               'max_feature_mrd_time':max_time,'history_le':hist})
    return {c:v for c,v in grouped.items() if v}


def logistic_predict(train_rows,test_rows,feature_key):
    ytr=np.array([r['y'] for r in train_rows],dtype=int)
    xtr=np.vstack([r[feature_key] for r in train_rows])
    xte=np.vstack([r[feature_key] for r in test_rows])
    beta=mod.fit_logit(np.column_stack([np.ones(len(xtr)),xtr]),ytr,penalty=RIDGE)
    pp=mod.sigmoid(np.column_stack([np.ones(len(xte)),xte])@beta)
    if not np.all(np.isfinite(pp)): raise FloatingPointError('nonfinite LOCO predictions')
    return pp


def loco(groups,feature_key):
    pred_rows=[]; fold_rows=[]
    cohort_order=sorted(groups)
    for held in cohort_order:
        test=groups[held]
        train=[r for c in cohort_order if c!=held for r in groups[c]]
        if not train: raise RuntimeError('LOCO training fold is empty')
        pp=logistic_predict(train,test,feature_key)
        yy=np.array([r['y'] for r in test],dtype=int)
        m=mod.score_metrics(yy,pp)
        fold_rows.append({'held_out_cohort':held,**m,'ridge_lambda':RIDGE})
        for r,pred in zip(test,pp):
            p=r['patient']
            pred_rows.append({'cohort':held,'patient_id':p['patient_id'],'outcome':r['y'],
                              'predicted_risk':float(pred),'event_time':p['event_time'],'followup_time':p['followup_time'],
                              'stage_III_IV':p['stage_high'],'rectal_primary':p['rectal_primary'],
                              'history_le':r['history_le'],'max_feature_mrd_time':r['max_feature_mrd_time']})
    y=np.asarray([r['outcome'] for r in pred_rows],dtype=int)
    p=np.asarray([r['predicted_risk'] for r in pred_rows],dtype=float)
    pooled=mod.score_metrics(y,p)
    return pred_rows,fold_rows,pooled


def quick_metrics(y,p):
    y=np.asarray(y,dtype=int);p=np.asarray(p,dtype=float)
    if len(np.unique(y))<2: raise ValueError('overall bootstrap sample has only one outcome class')
    return {'auc':mod.rank_auc(y,p),'average_precision':mod.average_precision(y,p),
            'brier':float(np.mean((y-p)**2)),'mean_predicted_risk':float(np.mean(p))}


def bootstrap_one_landmark(groups,landmark,reps=BOOTSTRAP_REPS,seed=BOOTSTRAP_SEED):
    rng=np.random.default_rng(seed+int(landmark))
    store={k:defaultdict(list) for k in ('clinical','mrd','delta')}
    failures=Counter()
    cohorts=sorted(groups)
    for b in range(reps):
        sampled={}
        for c in cohorts:
            source=groups[c]
            idx=rng.integers(0,len(source),size=len(source))
            sampled[c]=[source[int(i)] for i in idx]
        try:
            # Construct same bootstrap outcome vector/order to preserve paired comparison.
            # Both LOCO calls iterate cohorts and test rows in identical order.
            pred_c,_,_=loco(sampled,'x_clin')
            pred_m,_,_=loco(sampled,'x_mrd')
            y_c=np.asarray([r['outcome'] for r in pred_c],dtype=int)
            y_m=np.asarray([r['outcome'] for r in pred_m],dtype=int)
            if not np.array_equal(y_c,y_m): raise RuntimeError('paired bootstrap outcomes differ')
            vals_c=quick_metrics(y_c,np.asarray([r['predicted_risk'] for r in pred_c]))
            vals_m=quick_metrics(y_m,np.asarray([r['predicted_risk'] for r in pred_m]))
            for metric in ('auc','average_precision','brier'):
                store['clinical'][metric].append(vals_c[metric])
                store['mrd'][metric].append(vals_m[metric])
                store['delta'][metric].append(vals_m[metric]-vals_c[metric])
        except Exception as exc:
            failures[type(exc).__name__+': '+str(exc)[:100]]+=1
        if (b+1)%100==0:
            print(f"L{int(landmark)} bootstrap {b+1}/{reps}; valid={len(store['clinical']['auc'])}; failed={sum(failures.values())}",flush=True)
    return store,failures


def ci_rows_for_landmark(label,point_clin,point_mrd,store,failures,reps):
    out=[]
    point={'clinical':point_clin,'mrd':point_mrd,
           'delta':{k:point_mrd[k]-point_clin[k] for k in ('auc','average_precision','brier')}}
    names={'clinical':'clinical_only','mrd':'clinical_plus_MRD','delta':'MRD_minus_clinical'}
    for kind in ('clinical','mrd','delta'):
        for metric in ('auc','average_precision','brier'):
            vals=np.asarray(store[kind][metric],dtype=float)
            lo,hi=(np.quantile(vals,[0.025,0.975]).tolist() if len(vals) else (float('nan'),float('nan')))
            out.append({'landmark':label,'contrast':names[kind],'metric':metric,'estimate':point[kind][metric],
                        'ci_95_lower':lo,'ci_95_upper':hi,'bootstrap_reps_requested':reps,
                        'bootstrap_reps_valid':len(vals),'bootstrap_reps_failed':sum(failures.values()),
                        'failure_reasons':'; '.join(f'{k} ({v})' for k,v in failures.items())})
    return out


def state_pattern(hist):
    s=[x[1] for x in hist]
    if not s:return 'no_mrd'
    if len(s)==1:return 'single_positive' if s[0] else 'single_negative'
    if all(x==0 for x in s):return 'persistent_negative'
    if all(x==1 for x in s):return 'persistent_positive'
    if s[-1]==0 and any(x==1 for x in s[:-1]):return 'positive_to_negative_clearance'
    if s[-1]==1 and all(x==0 for x in s[:-1]):return 'negative_to_positive_conversion'
    return 'intermittent_mixed'


def aggregate_trajectory_summary(trajectories):
    groups=defaultdict(list)
    for row in trajectories: groups[row['mrd_pattern']].append(row)
    out=[]
    for pattern,rows in sorted(groups.items()):
        item={'mrd_pattern':pattern,'n_with_three_landmarks':len(rows),
              'events_12_to_24':sum(r['y12_24']==1 for r in rows)}
        for col in ('risk_3m','risk_6m','risk_12m'):
            a=np.asarray([r[col] for r in rows],dtype=float)
            item[f'median_{col}']=float(np.median(a));item[f'mean_{col}']=float(np.mean(a))
        item['median_delta_risk_3_to_6']=float(np.median([r['risk_6m']-r['risk_3m'] for r in rows]))
        item['median_delta_risk_6_to_12']=float(np.median([r['risk_12m']-r['risk_6m'] for r in rows]))
        out.append(item)
    return out


def choose_examples(trajectories):
    # Select observed public-data examples; absent patterns stay explicitly unavailable.
    patterns={}
    for r in trajectories:
        s=[x[1] for x in r['history_le12']]
        if len(s)<2: continue
        pat=r['mrd_pattern']
        if pat=='persistent_negative':
            patterns.setdefault('A sustained MRD-negative',[]).append(r)
        elif pat=='positive_to_negative_clearance':
            patterns.setdefault('B early-positive then cleared',[]).append(r)
        elif pat=='negative_to_positive_conversion':
            patterns.setdefault('C negative then converted positive',[]).append(r)
        elif pat=='persistent_positive':
            patterns.setdefault('D persistent MRD-positive',[]).append(r)
        if r['y12_24']==1 and pat=='negative_to_positive_conversion':
            patterns.setdefault('E new positive by 12m then recurrence by 24m',[]).append(r)
    selected=[]
    for label,cases in patterns.items():
        if label.startswith('A'): case=min(cases,key=lambda r:r['risk_12m'])
        elif label.startswith('B'): case=min(cases,key=lambda r:r['risk_12m']-r['risk_3m'])
        elif label.startswith('C'): case=max(cases,key=lambda r:r['risk_12m']-r['risk_3m'])
        else: case=max(cases,key=lambda r:r['risk_12m'])
        selected.append((label,case))
    return selected


def svg_trajectories(selected,path):
    W,H=1320,190+len(selected)*145
    colors={'neg':'#18864b','pos':'#c83b3b','landmark':'#4361a8','event':'#191919'}
    pieces=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
            '<rect width="100%" height="100%" fill="white"/>',
            '<style>text{font-family:Arial,"Microsoft YaHei",sans-serif;fill:#202124}.small{font-size:13px}.label{font-size:15px;font-weight:bold}.title{font-size:22px;font-weight:bold}</style>',
            '<text x="36" y="34" class="title">典型患者 MRD 与 cross-fitted 风险轨迹（患者级，仅本地）</text>',
            '<text x="36" y="62" class="small">风险表示从当前 landmark 到预设终点的复发风险：L3 (3,12]；L6 (6,12]；L12 (12,24]。实际复发时间仅作为结局标记。</text>']
    x0,x1=400,1170
    def xx(t): return x0+(x1-x0)*min(max(float(t),0),24)/24
    for m in range(0,25,3):
        x=xx(m);pieces.append(f'<line x1="{x}" y1="86" x2="{x}" y2="{H-38}" stroke="#d8dce2" stroke-width="1"/>')
        pieces.append(f'<text x="{x}" y="{H-15}" class="small" text-anchor="middle">{m}月</text>')
    for i,(label,r) in enumerate(selected):
        y=115+i*145
        pieces.append(f'<text x="34" y="{y-12}" class="label">{html.escape(label)}</text>')
        pieces.append(f'<text x="34" y="{y+13}" class="small">{html.escape(r["cohort"])} / {html.escape(r["patient_id"])}</text>')
        pieces.append(f'<line x1="{x0}" y1="{y+25}" x2="{x1}" y2="{y+25}" stroke="#555" stroke-width="2"/>')
        pieces.append(f'<polygon points="{xx(0)-7},{y+25} {xx(0)+7},{y+25} {xx(0)},{y+12}" fill="#222"/>')
        pieces.append(f'<text x="{xx(0)+10}" y="{y+50}" class="small">手术</text>')
        for l,risk,color in [(3,r['risk_3m'],'#7255a3'),(6,r['risk_6m'],'#3478a8'),(12,r['risk_12m'],'#b06d18')]:
            x=xx(l);pieces.append(f'<line x1="{x}" y1="{y-42}" x2="{x}" y2="{y+38}" stroke="{color}" stroke-width="2" stroke-dasharray="5,4"/>')
            pieces.append(f'<text x="{x+5}" y="{y-25}" class="small" fill="{color}">L{l} {risk:.0%}</text>')
        for t,s in r['history_le12']:
            x=xx(t);col=colors['pos'] if s else colors['neg']
            pieces.append(f'<circle cx="{x}" cy="{y+25}" r="7" fill="{col}" stroke="white" stroke-width="1.5"/>')
            pieces.append(f'<text x="{x}" y="{y+68}" class="small" text-anchor="middle">{t:g}月 {"+" if s else "−"}</text>')
        et=r['event_time']
        if et is not None:
            x=xx(et);pieces.append(f'<path d="M{x-8},{y+15} L{x+8},{y+35} M{x+8},{y+15} L{x-8},{y+35}" stroke="{colors["event"]}" stroke-width="3"/>')
            pieces.append(f'<text x="{x+10}" y="{y+5}" class="small">复发 {et:g}月</text>')
        else:
            pieces.append(f'<text x="{x1-4}" y="{y+56}" class="small" text-anchor="end">无复发记录；随访≥24月</text>')
    pieces.append('<circle cx="1010" cy="68" r="6" fill="#18864b"/><text x="1021" y="73" class="small">MRD−</text><circle cx="1080" cy="68" r="6" fill="#c83b3b"/><text x="1091" y="73" class="small">MRD+</text>')
    pieces.append('</svg>')
    Path(path).write_text('\n'.join(pieces),encoding='utf-8')


def svg_metric_plot(rows,metric,title,path,ymin,ymax):
    W,H=1000,610; left,right,top,bottom=95,930,75,525
    xs=[260,510,760]; labs=['3月：风险至12月','6月：风险至12月','12月：风险至24月']
    colors={'clinical_only':'#4472c4','clinical_plus_MRD':'#d97922'}
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">','<rect width="100%" height="100%" fill="white"/>',
         '<style>text{font-family:Arial,"Microsoft YaHei",sans-serif;fill:#202124}.axis{font-size:14px}.title{font-size:21px;font-weight:bold}</style>',
         f'<text x="{W/2}" y="32" text-anchor="middle" class="title">{html.escape(title)}</text>']
    def yy(v):return bottom-(float(v)-ymin)/(ymax-ymin)*(bottom-top)
    for k in range(6):
        v=ymin+(ymax-ymin)*k/5;y=yy(v)
        out.append(f'<line x1="{left}" y1="{y}" x2="{right}" y2="{y}" stroke="#e1e4e8"/>')
        out.append(f'<text x="{left-12}" y="{y+5}" text-anchor="end" class="axis">{v:.2f}</text>')
    for x,l in zip(xs,labs):
        out.append(f'<text x="{x}" y="{bottom+34}" text-anchor="middle" class="axis">{html.escape(l)}</text>')
    for j,lm in enumerate((3,6,12)):
        for kind,offset in [('clinical_only',-24),('clinical_plus_MRD',24)]:
            r=next(x for x in rows if x['landmark']==f'L{lm}' and x['contrast']==kind and x['metric']==metric)
            x=xs[j]+offset;lo=r['ci_95_lower'];hi=r['ci_95_upper'];est=r['estimate']
            y=yy(est);yl=yy(lo);yh=yy(hi);color=colors[kind]
            out.append(f'<line x1="{x}" y1="{yl}" x2="{x}" y2="{yh}" stroke="{color}" stroke-width="3"/><line x1="{x-7}" y1="{yl}" x2="{x+7}" y2="{yl}" stroke="{color}" stroke-width="2"/><line x1="{x-7}" y1="{yh}" x2="{x+7}" y2="{yh}" stroke="{color}" stroke-width="2"/><circle cx="{x}" cy="{y}" r="7" fill="{color}"/><text x="{x}" y="{y-12}" text-anchor="middle" class="axis">{est:.3f}</text>')
    out.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#333"/><line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#333"/>')
    out.append('<rect x="320" y="565" width="18" height="12" fill="#4472c4"/><text x="345" y="577" class="axis">Clinical only</text><rect x="505" y="565" width="18" height="12" fill="#d97922"/><text x="530" y="577" class="axis">Clinical + MRD</text>')
    out.append('<text x="500" y="602" text-anchor="middle" class="axis">误差线为完整 LOCO patient-stratified bootstrap percentile 95% CI；L12 仅 Chen/COSMOS。</text></svg>')
    Path(path).write_text('\n'.join(out),encoding='utf-8')


def main():
    all_cohorts=all_cohort_patients()
    # The L12 cohort gate comes from audited patient-level endpoint ascertainability.
    l12_groups_all=make_risk_rows(all_cohorts,12)
    l12_include=[]
    for c,records in l12_groups_all.items():
        ev=sum(r['y']==1 for r in records);ctrl=len(records)-ev
        if ev>=10 and ctrl>=10: l12_include.append(c)
    if set(l12_include)!={'Chen','COSMOS'}:
        raise RuntimeError(f"L12 cohort gate changed unexpectedly: {l12_include}; inspect l12_cohort_audit.csv")
    groups_by_l={3:make_risk_rows(all_cohorts,3),6:make_risk_rows(all_cohorts,6),12:{c:l12_groups_all[c] for c in l12_include}}
    # Cohort-level endpoint totals guard against accidental changes to the audited definitions.
    expected={3:{'Chen':(181,15),'ColonAiQ':(147,26),'COSMOS':(330,19)},
              6:{'Chen':(179,13),'ColonAiQ':(137,12),'COSMOS':(326,15)},
              12:{'Chen':(101,16),'COSMOS':(281,19)}}
    for L,groups in groups_by_l.items():
        for c,g in groups.items():
            got=(len(g),sum(r['y'] for r in g))
            if got!=expected[L][c]:raise RuntimeError(f"L{L} {c}: {got} != {expected[L][c]}")
            for r in g:
                assert r['max_feature_mrd_time']<=L

    metrics_rows=[];fold_rows=[];predictions=[];pred_by_landmark={}; boot_rows=[]
    point_cache={}
    for L,groups in groups_by_l.items():
        label=f'L{L}'
        pclin,fclin,mclin=loco(groups,'x_clin')
        pmrd,fmrd,mmrd=loco(groups,'x_mrd')
        # Same patients/labels for the two models at this landmark.
        yc=np.asarray([r['outcome'] for r in pclin],int);ym=np.asarray([r['outcome'] for r in pmrd],int)
        if not np.array_equal(yc,ym): raise RuntimeError(f'L{L}: clinical and MRD risks use different rows')
        for kind,met,folds,plist in [('clinical_only',mclin,fclin,pclin),('clinical_plus_MRD',mmrd,fmrd,pmrd)]:
            metrics_rows.append({'landmark':label,'model':kind,**met,'features':'stage_III_IV;rectal_primary' if kind=='clinical_only' else 'stage_III_IV;rectal_primary;current_mrd_positive;prior_mrd_positive_before_current',
                                 'cohorts':';'.join(sorted(groups)),'validation':'LOCO; no held-out cohort predictor'})
            for r in folds: fold_rows.append({'landmark':label,'model':kind,**r})
            for r in plist: predictions.append({'landmark':label,'model':kind,**r})
        pred_by_landmark[L]={'clinical':pclin,'mrd':pmrd}
        point_cache[L]={'clinical':mclin,'mrd':mmrd}
        print(f"Point LOCO L{L}: clinical AUC={mclin['auc']:.3f}; MRD AUC={mmrd['auc']:.3f}; n={mclin['n']}, events={mclin['events']}",flush=True)

    csv_write(OUT/'landmark_loco_metrics_all.csv',metrics_rows)
    csv_write(OUT/'landmark_loco_by_heldout_cohort_all.csv',fold_rows)
    csv_write(OUT/'l12_loco_metrics.csv',[r for r in metrics_rows if r['landmark']=='L12'])
    csv_write(OUT/'l12_loco_by_heldout_cohort.csv',[r for r in fold_rows if r['landmark']=='L12'])
    csv_write(LOCAL/'landmark_loco_patient_predictions.csv',predictions)

    # Stratified patient bootstrap; refit both models and every LOCO fold on each resample.
    rng_seed=BOOTSTRAP_SEED
    for L,groups in groups_by_l.items():
        store,failures=bootstrap_one_landmark(groups,L,BOOTSTRAP_REPS,rng_seed)
        boot_rows.extend(ci_rows_for_landmark(f'L{L}',point_cache[L]['clinical'],point_cache[L]['mrd'],store,failures,BOOTSTRAP_REPS))
    csv_write(OUT/'landmark_bootstrap_ci.csv',boot_rows)

    # Patient-specific trajectory among people with all three cross-fitted landmark predictions.
    pred_maps={}
    for L in (3,6,12):
        pred_maps[L]={(r['cohort'],r['patient_id']):r for r in pred_by_landmark[L]['mrd']}
    l12_candidates=groups_by_l[12]
    trajectories=[]
    for c,records in l12_candidates.items():
        for rec in records:
            p=rec['patient'];key=(c,p['patient_id'])
            if key not in pred_maps[3] or key not in pred_maps[6] or key not in pred_maps[12]: continue
            hist=[(t,s) for t,s in p['history'] if t<=12]
            assert max([t for t,s in hist],default=0)<=12
            trajectories.append({'cohort':c,'patient_id':p['patient_id'],'event_time':p['event_time'],
                                 'followup_time':p['followup_time'],'y12_24':rec['y'],'history_le12':hist,
                                 'mrd_pattern':state_pattern(hist),
                                 'risk_3m':pred_maps[3][key]['predicted_risk'],
                                 'risk_6m':pred_maps[6][key]['predicted_risk'],
                                 'risk_12m':pred_maps[12][key]['predicted_risk']})
    if len(trajectories)!=382: raise RuntimeError(f"Three-landmark trajectory table should cover L12 cohort n=382, got {len(trajectories)}")
    traj_summary=aggregate_trajectory_summary(trajectories)
    csv_write(OUT/'digital_twin_trajectory_summary.csv',traj_summary)
    trajectory_detail=[]
    for r in trajectories:
        hist=r['history_le12']
        trajectory_detail.append({**{k:v for k,v in r.items() if k!='history_le12'},
                                  'mrd_history_le12':';'.join(f'{t:g}:{s}' for t,s in hist)})
    csv_write(LOCAL/'digital_twin_patient_trajectories.csv',trajectory_detail)
    selected=choose_examples(trajectories)
    case_rows=[]
    for label,r in selected:
        case_rows.append({'pattern':label,'cohort':r['cohort'],'patient_id':r['patient_id'],'event_time':r['event_time'],
                          'y12_24':r['y12_24'],'risk_3m':r['risk_3m'],'risk_6m':r['risk_6m'],'risk_12m':r['risk_12m'],
                          'mrd_history_le12':';'.join(f'{t:g}:{s}' for t,s in r['history_le12'])})
    csv_write(LOCAL/'digital_twin_typical_cases_local.csv',case_rows)
    svg_trajectories(selected,FIG/'digital_twin_patient_trajectories.svg')
    svg_metric_plot(boot_rows,'auc','LOCO AUC：3、6、12月 landmark',FIG/'landmark_auc_comparison.svg',0.45,1.0)
    svg_metric_plot(boot_rows,'brier','LOCO Brier score：3、6、12月 landmark',FIG/'landmark_brier_comparison.svg',0.0,0.20)

    # Human-readable concise report. Detailed outcome ascertainment is in the audit report.
    by_l={}
    for r in metrics_rows: by_l[(r['landmark'],r['model'])]=r
    ci={(r['landmark'],r['contrast'],r['metric']):r for r in boot_rows}
    rows_md=[]
    for L in (3,6,12):
        lm=f'L{L}';a=by_l[(lm,'clinical_only')];b=by_l[(lm,'clinical_plus_MRD')]
        row=f"| L{L} | {a['n']} | {a['events']} | {a['auc']:.3f} [{ci[(lm,'clinical_only','auc')]['ci_95_lower']:.3f}, {ci[(lm,'clinical_only','auc')]['ci_95_upper']:.3f}] | {b['auc']:.3f} [{ci[(lm,'clinical_plus_MRD','auc')]['ci_95_lower']:.3f}, {ci[(lm,'clinical_plus_MRD','auc')]['ci_95_upper']:.3f}] | {ci[(lm,'MRD_minus_clinical','auc')]['estimate']:+.3f} [{ci[(lm,'MRD_minus_clinical','auc')]['ci_95_lower']:+.3f}, {ci[(lm,'MRD_minus_clinical','auc')]['ci_95_upper']:+.3f}] | {a['brier']:.3f} [{ci[(lm,'clinical_only','brier')]['ci_95_lower']:.3f}, {ci[(lm,'clinical_only','brier')]['ci_95_upper']:.3f}] | {b['brier']:.3f} [{ci[(lm,'clinical_plus_MRD','brier')]['ci_95_lower']:.3f}, {ci[(lm,'clinical_plus_MRD','brier')]['ci_95_upper']:.3f}] | {ci[(lm,'MRD_minus_clinical','brier')]['estimate']:+.3f} [{ci[(lm,'MRD_minus_clinical','brier')]['ci_95_lower']:+.3f}, {ci[(lm,'MRD_minus_clinical','brier')]['ci_95_upper']:+.3f}] |"
        rows_md.append(row)
    fold_md='\n'.join(f"| {r['landmark']} | {r['model']} | {r['held_out_cohort']} | {r['n']} | {r['events']} | {r['auc']:.3f} | {r['brier']:.3f} | {r['mean_predicted_risk']:.3f} | {r['calibration_intercept']:.3f} | {r['calibration_slope']:.3f} |" for r in fold_rows)
    summary_md='\n'.join(f"| {r['mrd_pattern']} | {r['n_with_three_landmarks']} | {r['events_12_to_24']} | {r['median_risk_3m']:.3f} | {r['median_risk_6m']:.3f} | {r['median_risk_12m']:.3f} |" for r in traj_summary)
    report = """# CRC MRD 纵向预测 digital twin：三队列纵向模型

**日期：**2026-09-30  
**项目：**E:/crc_mrd_digital_twin

## Digital twin 定义

预测型 digital twin 用患者基线临床特征初始化纵向疾病状态，并在获得新的术后 MRD 检测后更新。在每个 landmark，只使用该时点及此前可获得的信息预测后续复发风险。

形式化表示为 Z_t = f(X_baseline, MRD_{≤t})，Risk_t = P(T_recurrence ∈ (t,H_t] | Z_t)。本轮仅纳入可跨队列统一的 stage、site 和 assay-agnostic 二分类 MRD 状态；ACT 给药时间、CEA 及肿瘤分子基线无法在公开队列间可靠统一，未纳入模型。

## 风险集与时间窗

- L3：只使用采样时间≤3月的信息，预测术后 (3,12] 月复发。
- L6：只使用采样时间≤6月的信息，预测术后 (6,12] 月复发。
- L12：只使用采样时间≤12月的信息，预测术后 (12,24] 月复发。
- 三个风险对应不同的未来时间窗。MRD 特征包括当前状态及既往检测是否曾阳性。缺失检测不编码为阴性。
- 每位患者特征所用 MRD 时间均通过断言检查，必须≤landmark；L12 的任何>12月 MRD 检测均会过滤，若有泄漏脚本将中止。

## L12 结局可判定性与队列纳入

L12 结局审计见 reports/12m_landmark_24m_outcome_audit_zh.md。ColonAiQ 有19例在 (12,24] 月复发，但12月风险患者中没有可确认随访至24月的未复发对照；123人的24月结局未知，因此整个 ColonAiQ 队列从 L12 模型和 LOCO 排除。L12 主模型仅使用 Chen 与 COSMOS：N=382、事件=35（Chen 101人/16事件；COSMOS 281人/19事件）。L3/L6 仍使用 Chen、ColonAiQ、COSMOS 三队列。L12 按队列整体排除，不会从结局不完整的队列单独拼入事件病例。

## 模型与验证

低维 ridge logistic，λ=@@RIDGE@@。每个时间点分别拟合 clinical-only 和 clinical+MRD 模型；逐一留出整个队列。L3/L6 为三折 LOCO；L12 因 ColonAiQ 结局无法判定而只有两个队列，每轮用一个队列训练、另一个队列留出。这是有限的 internal-external transportability assessment，不是独立外部验证。

## 完整 LOCO 患者分层 bootstrap 95% CI

在每个队列内以患者为单位有放回抽样；每次对 clinical-only、clinical+MRD 模型以及每个 LOCO fold 重新拟合。同一 bootstrap replicate 用于配对差异。重复次数 @@REPS@@，随机种子 @@SEED@@。差值均为 MRD 模型减 clinical-only；ΔBrier 为负表示 MRD 模型 Brier 更低。区间为 percentile 95% CI，成功/失败次数见 outputs/tables/landmark_bootstrap_ci.csv。AUC、AP、Brier 差值未进行额外假设检验。

| Landmark | N | 事件 | AUC clinical | AUC + MRD | ΔAUC | Brier clinical | Brier + MRD | ΔBrier |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
@@ROWS@@

## LOCO 点估计和校准

| Landmark | 模型 | 留出队列 | N | 事件 | AUC | Brier | 平均预测风险 | 校准截距 | 校准斜率 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
@@FOLDS@@

校准截距和斜率在每个留出队列单独计算；若某折近似分离或不可识别则记为 NaN。汇总性能应与各留出队列表现一起解读，尤其是 L12 两队列的小样本结果。

## Digital twin 风险轨迹

对同时进入 L3、L6、L12 LOCO 的患者，保存 cross-fitted Risk_3m → Risk_6m → Risk_12m。三个风险各自对应不同的未来时间窗，图和数据表均按窗口标注。

| 截至12月 MRD 状态模式 | 三时点人数 | 12–24月事件 | 中位 Risk3 | 中位 Risk6 | 中位 Risk12 |
|---|---:|---:|---:|---:|---:|
@@SUMMARY@@

典型患者图展示实际 MRD 采样、landmark 预测和临床复发时间。患者级 ID 与轨迹仅保存在本地 outputs/local_only/，患者轨迹 SVG 也仅保存在 outputs/figures/digital_twin_patient_trajectories.svg，不会提交到 GitHub。若某种模式在可判定病例中不存在，图中不会虚构病例。

## 解释边界

本工作属于 multicohort development + internal-external validation。LOCO 评估跨公开队列的 transportability，不等于独立院外前瞻性验证；clinical utility 尚未通过临床决策实验或前瞻性研究验证。独立 external validation 仍预留给院内队列。L12 仅覆盖 Chen/COSMOS 两个来源，估计不确定性较大，应视为探索性分析。

## 输出文件

- outputs/tables/l12_cohort_audit.csv
- outputs/tables/l12_loco_metrics.csv
- outputs/tables/l12_loco_by_heldout_cohort.csv
- outputs/tables/landmark_bootstrap_ci.csv
- outputs/tables/digital_twin_trajectory_summary.csv
- outputs/figures/landmark_auc_comparison.svg
- outputs/figures/landmark_brier_comparison.svg
- 患者级轨迹数据和典型病例仅保存在 outputs/local_only/，患者轨迹 SVG 也不会提交到 GitHub。
"""
    report=report.replace("@@RIDGE@@",str(RIDGE)).replace("@@REPS@@",str(BOOTSTRAP_REPS)).replace("@@SEED@@",str(BOOTSTRAP_SEED))
    report=report.replace("@@ROWS@@",chr(10).join(rows_md)).replace("@@FOLDS@@",fold_md).replace("@@SUMMARY@@",summary_md)
    (REPORTS/'digital_twin_longitudinal_model_zh.md').write_text(report,encoding='utf-8')
    summary={'n_landmark':{str(L):sum(len(g) for g in groups_by_l[L].values()) for L in groups_by_l},
             'metrics':metrics_rows,'bootstrap_reps':BOOTSTRAP_REPS,'bootstrap_failures':{str(L):sum(1 for r in boot_rows if r['landmark']==f'L{L}' and r['contrast']=='clinical_only' and r['metric']=='auc' and r['bootstrap_reps_valid']<BOOTSTRAP_REPS) for L in groups_by_l},
             'trajectory_examples':[{'pattern':x[0],'cohort':x[1]['cohort'],'patient_id':x[1]['patient_id']} for x in selected],
             'report':str(REPORTS/'digital_twin_longitudinal_model_zh.md')}
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
