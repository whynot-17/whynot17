from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(r"D:\DINP\whynot17-new-project")
OUT=ROOT/"results"/"2026-10-12"/"RETSAT_CRC_epithelial_paired_pathways"

def bh(p):
    a=np.asarray(p,float); out=np.full(len(a),np.nan); ok=np.isfinite(a)
    if not ok.any(): return out
    v=a[ok]; order=np.argsort(v); r=v[order]
    q=r*len(r)/np.arange(1,len(r)+1); q=np.minimum.accumulate(q[::-1])[::-1]
    tmp=np.empty_like(q);tmp[order]=np.minimum(q,1)
    out[ok]=tmp;return out

def wil(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b);a=a[ok];b=b[ok]
    if len(a)<2:return np.nan
    if np.allclose(a-b,0):return 1.
    return float(wilcoxon(a,b,zero_method="wilcox",alternative="two-sided").pvalue)

def regenerate_weighted_localization():
    loc=pd.read_csv(OUT/"singlecell_RETSAT_localization_summary.csv")
    loc["weighted_detection_pct"]=100*loc["positive_cells"]/loc["tumor_cells"]
    loc.to_csv(OUT/"singlecell_RETSAT_localization_summary.csv",index=False)
    for dataset in loc.dataset.unique():
        loc[loc.dataset.eq(dataset)].to_csv(OUT/f"singlecell_RETSAT_localization_{dataset}.csv",index=False)
    primary=loc[loc.dataset.eq("GSE132465")].sort_values("weighted_detection_pct")
    fig,ax=plt.subplots(figsize=(8.8,4.8))
    colors=["#4575b4" if x=="Epithelial cells" else "#bdbdbd" for x in primary.cell_type]
    ax.barh(primary.cell_type,primary.weighted_detection_pct,color=colors)
    ax.set_xlabel("RETSAT-detected cells (%)");ax.set_ylabel("Author-annotated cell type (tumor samples)")
    ax.set_title("RETSAT localization across tumor cell types — GSE132465")
    for i,r in enumerate(primary.itertuples()):
        ax.text(r.weighted_detection_pct+0.2,i,f"{int(r.positive_cells)}/{int(r.tumor_cells)}",va="center",fontsize=8)
    fig.tight_layout();fig.savefig(OUT/"RETSAT_GSE132465_celltype_localization.png",dpi=220);plt.close(fig)
    return loc

def create_coverage(gt):
    records=[]
    for cohort in ["GSE132465","GSE200997","GSE44076","GSE32323"]:
        expr=pd.read_csv(OUT/f"{cohort}_paired_expression_and_scores.csv",index_col=0)
        for pathway in ["Retinol_metabolism","Ferroptosis"]:
            g=gt[(gt.cohort==cohort)&(gt.pathway==pathway)].gene.tolist()
            variable=[x for x in g if x in expr.columns and expr[x].std(ddof=0)>1e-12]
            records.append({"cohort":cohort,"pathway":pathway,"KEGG_members_after_target_exclusion":68 if pathway=="Retinol_metabolism" else 42,
                            "members_in_expression_matrix":len(g),"variable_members_used_in_score":len(variable),
                            "variable_member_symbols":";".join(variable)})
    df=pd.DataFrame(records)
    df.to_csv(OUT/"pathway_gene_coverage.csv",index=False)
    return df

def cellcount_sensitivity(pairs,gt):
    records=[];countrows=[]
    for cohort in ["GSE132465","GSE200997"]:
        d=pairs[pairs.cohort.eq(cohort)].copy()
        if "n_cells" not in d: continue
        counts=d.pivot(index="patient",columns="condition",values="n_cells")
        for patient,row in counts.iterrows():
            if {"Normal","Tumor"}.issubset(counts.columns):
                countrows.append({"cohort":cohort,"patient":patient,
                                  "normal_epithelial_cells":int(row["Normal"]),
                                  "tumor_epithelial_cells":int(row["Tumor"]),
                                  "min_pair_cells":int(min(row["Normal"],row["Tumor"])),
                                  "passes_10_cells":bool(min(row["Normal"],row["Tumor"])>=10),
                                  "passes_20_cells":bool(min(row["Normal"],row["Tumor"])>=20)})
        for threshold in [10,20]:
            keep=[pid for pid,row in counts.iterrows() if row.notna().all() and min(row["Normal"],row["Tumor"])>=threshold]
            sub=d[d.patient.isin(keep)].copy()
            for pathway in ["Retinol_metabolism","Ferroptosis"]:
                genes=gt[(gt.cohort==cohort)&(gt.pathway==pathway)].gene.tolist()
                genes=[g for g in genes if g in sub.columns and sub[g].std(ddof=0)>1e-12]
                if not genes:
                    sub[f"{pathway}_score"]=np.nan
                else:
                    x=sub[genes].astype(float)
                    sd=x.std(axis=0,ddof=0)
                    v=sd[sd>1e-12].index
                    sub[f"{pathway}_score"]=((x[v]-x[v].mean(axis=0))/sd[v]).mean(axis=1)
            for measure in ["RETSAT_log2expr","Retinol_metabolism_score","Ferroptosis_score"]:
                if measure=="RETSAT_log2expr": pass
                elif measure in sub.columns:
                    # Recompute scores in the current threshold-selected patients.
                    pathway=measure.replace("_score","")
                    genes=gt[(gt.cohort==cohort)&(gt.pathway==pathway)].gene.tolist()
                    genes=[g for g in genes if g in sub.columns]
                    if genes:
                        x=sub[genes].astype(float);sd=x.std(axis=0,ddof=0);v=sd[sd>1e-12].index
                        sub[measure]=((x[v]-x[v].mean(axis=0))/sd[v]).mean(axis=1)
                pv=[];nv=[];diff=[]
                for pid in keep:
                    g=sub[sub.patient.eq(pid)]
                    if (g.condition=="Tumor").sum()!=1 or (g.condition=="Normal").sum()!=1: continue
                    t=float(g.loc[g.condition.eq("Tumor"),measure].iloc[0]);n=float(g.loc[g.condition.eq("Normal"),measure].iloc[0])
                    pv.append(t);nv.append(n);diff.append(t-n)
                records.append({"cohort":cohort,"min_cells_per_sample":threshold,"measure":measure,
                                "n_pairs":len(diff),"median_paired_delta_tumor_minus_normal":float(np.median(diff)) if diff else np.nan,
                                "fraction_positive":float(np.mean(np.asarray(diff)>0)) if diff else np.nan,
                                "p_value":wil(pv,nv)})
    sens=pd.DataFrame(records)
    sens["BH_q_within_cohort_threshold_3measures"]=sens.groupby(["cohort","min_cells_per_sample"])["p_value"].transform(lambda x:bh(x.to_numpy()))
    sens.to_csv(OUT/"singlecell_cell_count_sensitivity.csv",index=False)
    pd.DataFrame(countrows).to_csv(OUT/"singlecell_epithelial_paired_cell_counts.csv",index=False)
    return sens,pd.DataFrame(countrows)

def report(loc,summary,corr,gt,coverage,sens,counts):
    pair_counts=pd.read_csv(OUT/"cohort_pair_counts.csv")
    # Localized RETSAT proportions in epithelial vs non-epithelial tumor cells
    loc_metrics={}
    for cohort in ["GSE132465","GSE200997"]:
        d=loc[loc.dataset.eq(cohort)].copy()
        e=d[d.cell_type.eq("Epithelial cells")] if cohort=="GSE132465" else d[d.cell_type.eq("epithelial_like")]
        ne=d[~d.cell_type.isin(e.cell_type)]
        epi_rate=float(e.positive_cells.sum()/e.tumor_cells.sum()*100)
        non_rate=float(ne.positive_cells.sum()/ne.tumor_cells.sum()*100)
        loc_metrics[cohort]=(e,epi_rate,non_rate,epi_rate/non_rate if non_rate else np.nan)
    # Cohort-specific pathway stats
    def st(co, measure):
        return summary[(summary.cohort==co)&(summary.measure==measure)].iloc[0]
    def cr(co,path):
        return corr[(corr.cohort==co)&corr.association.eq(f"paired_delta_RETSAT_vs_{path}_score")].iloc[0]
    overlap=set(gt[(gt.cohort=="GSE44076")&(gt.pathway=="Retinol_metabolism")&(gt.BH_q_within_pathway_cohort<.05)&(gt.median_paired_delta<0)].gene)
    overlap &= set(gt[(gt.cohort=="GSE32323")&(gt.pathway=="Retinol_metabolism")&(gt.BH_q_within_pathway_cohort<.05)&(gt.median_paired_delta<0)].gene)
    s132= sens[(sens.cohort=="GSE132465")&(sens.min_cells_per_sample==20)&(sens.measure=="RETSAT_log2expr")].iloc[0]
    lines=[
    "# RETSAT：CRC恶性上皮定位、患者级配对与独立通路线索","",
    "## 结论摘要","",
    f"- **细胞定位：** GSE132465中RETSAT主要落在肿瘤上皮细胞：{loc_metrics['GSE132465'][1]:.1f}%肿瘤上皮细胞检出RETSAT，非上皮细胞为{loc_metrics['GSE132465'][2]:.1f}%（检出率比{loc_metrics['GSE132465'][3]:.1f}倍）；83.7%的肿瘤样本RETSAT阳性细胞为上皮细胞。GSE200997的标记基因上皮样细胞检出率{loc_metrics['GSE200997'][1]:.1f}%，非上皮细胞{loc_metrics['GSE200997'][2]:.1f}%（{loc_metrics['GSE200997'][3]:.1f}倍）。",
    f"- **配对表达：** 四个队列中RETSAT的肿瘤−正常中位差均为负。GSE132465 10对时未达显著（Δ={st('GSE132465','RETSAT_log2expr').median_paired_delta:.2f}, P={st('GSE132465','RETSAT_log2expr').p_value:.3g}）；GSE200997 7对（Δ={st('GSE200997','RETSAT_log2expr').median_paired_delta:.2f}, P={st('GSE200997','RETSAT_log2expr').p_value:.3g}），GSE44076 98对和GSE32323 17对均显著下降。",
    f"- **视黄醇代谢：** 排除RETSAT自身后，两个bulk配对队列的视黄醇代谢分数均下降（GSE44076 Δ={st('GSE44076','Retinol_metabolism_score').median_paired_delta:.2f}, P={st('GSE44076','Retinol_metabolism_score').p_value:.2g}; GSE32323 Δ={st('GSE32323','Retinol_metabolism_score').median_paired_delta:.2f}, P={st('GSE32323','Retinol_metabolism_score').p_value:.2g}），且患者配对的RETSAT差值与通路差值正相关（ρ分别为{cr('GSE44076','Retinol_metabolism').spearman_rho:.2f}、{cr('GSE32323','Retinol_metabolism').spearman_rho:.2f}）。但是该通路总分在两套单细胞肿瘤上皮伪bulk中都未显著，故当前证据支持bulk层面的独立视黄醇代谢线索，尚未证明它是肿瘤上皮特异改变。",
    "- **铁死亡：** KEGG铁死亡转录分数在GSE132465肿瘤上皮伪bulk与GSE44076 bulk中上升；GSE200997及GSE32323不显著。通路成员方向有升有降，整体缺少跨队列一致性，不能据此推断铁死亡增强或减弱。",
    "",
    "## 研究设计与统计",
    "",
    "只读取本项目work目录下的原始GEO表达矩阵、原始细胞注释及平台探针注释；未读取历史分析结果。GSE132465采用作者给出的Cell_type标签；将肿瘤样本的上皮细胞称为肿瘤上皮候选群，未在本分析中进行inferCNV/拷贝数恶性判定。该队列含23位CRC患者肿瘤样本及10位配对正常黏膜样本。[GSE132465 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)。",
    "",
    "GSE200997的本地细胞注释没有通用细胞类型字段，因此以EPCAM、KRT8、KRT18、KRT19、KRT20中至少2个基因检出定义上皮样细胞；配对肿瘤样本细胞数均≥49，正常样本均≥161。该操作性注释用于敏感性定位，不等同于CNV恶性鉴定。[GSE200997 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997)。",
    "",
    "单细胞按患者/样本与细胞类型汇总原始UMI为伪bulk，RETSAT按总UMI计算log2(CPM+0.5)。配对患者为统计单位，不把单个细胞当作独立重复。所有单细胞样本上皮细胞数要求≥10；另报告≥20细胞敏感性分析。GSE44076与GSE32323使用GEO提供的log2 RMA处理表达值；GSE44076保留98对CRC与相邻正常，排除50名健康供者；GSE32323保留17对患者样本，排除细胞系。[GSE44076 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44076)，[GSE32323 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32323)。",
    "",
    "视黄醇代谢和铁死亡基因集分别来自KEGG人类通路[hsa00830](https://www.kegg.jp/pathway/hsa00830)与[hsa04216](https://www.kegg.jp/pathway/hsa04216)。视黄醇集合原有69个成员，其中RETSAT被剔除，余68个；通路分数是各队列每个基因跨配对样本标准化后Z分数的均值。所有配对差值为肿瘤−正常。配对检验使用双侧Wilcoxon signed-rank；三项主要指标（RETSAT、视黄醇分数、铁死亡分数）在队列内BH校正；通路成员在各队列/通路内BH校正。RETSAT与通路的关联用患者配对差值做Spearman相关。",
    "",
    "## 细胞定位","",
    "| 队列 | 上皮定义 | 肿瘤上皮细胞 | RETSAT阳性上皮细胞 | 上皮检出率 | 非上皮检出率 | 上皮/非上皮率比 |",
    "|---|---|---:|---:|---:|---:|---:|"]
    for co,label in [("GSE132465","作者注释Epithelial cells"),("GSE200997","≥2/5上皮标志物")]:
        e,er,nr,ratio=loc_metrics[co]
        ncell=int(e.tumor_cells.sum());npos=int(e.positive_cells.sum())
        lines.append(f"| {co} | {label} | {ncell} | {npos} | {er:.2f}% | {nr:.2f}% | {ratio:.2f} |")
    lines += ["","GSE132465各细胞类型的RETSAT检出比例见图；上皮细胞不是唯一RETSAT阳性类型，但在肿瘤样本中显著富集。GSE200997提示方向一致，但其标记分类是本次重建，证据等级较低。","",
    "## 患者配对结果","",
    "| 队列 | 指标 | 配对数 | 中位Δ | Δ为正比例 | Wilcoxon P | 队列内BH q |",
    "|---|---|---:|---:|---:|---:|---:|"]
    for r in summary.to_dict("records"):
        lines.append(f"| {r['cohort']} | {r['measure']} | {int(r['n_pairs'])} | {r['median_paired_delta']:.4g} | {r['fraction_delta_positive']:.1%} | {r['p_value']:.3g} | {r['BH_q_within_cohort_3measures']:.3g} |")
    lines += ["","RETSAT在四个队列的配对中位变化方向均为下降。GSE132465主单细胞配对分析10对，P=0.160；移除每侧上皮细胞<20的配对后为8对，RETSAT中位差仍为负且P={:.3g}（见敏感性表）。".format(s132.p_value),
    "",
    "### 单细胞每对上皮细胞数敏感性","",
    "| 队列 | 最低细胞数/样本 | 配对数 | RETSAT Δ | P | 视黄醇分数 Δ | P | 铁死亡分数 Δ | P |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for co in ["GSE132465","GSE200997"]:
        for threshold in [10,20]:
            d=sens[(sens.cohort==co)&(sens.min_cells_per_sample==threshold)].set_index("measure")
            rows=[d.loc[x] for x in ["RETSAT_log2expr","Retinol_metabolism_score","Ferroptosis_score"]]
            lines.append(f"| {co} | {threshold} | {int(rows[0].n_pairs)} | {rows[0].median_paired_delta_tumor_minus_normal:.3f} | {rows[0].p_value:.3g} | {rows[1].median_paired_delta_tumor_minus_normal:.3f} | {rows[1].p_value:.3g} | {rows[2].median_paired_delta_tumor_minus_normal:.3f} | {rows[2].p_value:.3g} |")
    lines += ["","### RETSAT变化与独立通路分数变化的患者级相关","",
    "| 队列 | 通路 | 配对数 | Spearman ρ | P | 队列内BH q |","|---|---|---:|---:|---:|---:|"]
    for r in corr.to_dict("records"):
        name="视黄醇代谢" if "Retinol" in r["association"] else "铁死亡"
        lines.append(f"| {r['cohort']} | {name} | {int(r['n_pairs'])} | {r['spearman_rho']:.3f} | {r['p_value']:.3g} | {r['BH_q_within_cohort_2correlations']:.3g} |")
    lines += ["","## 独立视黄醇代谢与铁死亡线索","",
    f"两个bulk队列共同有{len(overlap)}个视黄醇代谢通路成员在肿瘤中下降且各自队列内FDR<0.05：{', '.join(sorted(overlap))}。通路评分中没有RETSAT本身，因此这一结果是RETSAT之外成员提供的转录层面支持；不过两套单细胞上皮伪bulk均未显示视黄醇总分显著下降，提示需要将其解释为bulk组织层面的线索，不能直接归因于恶性上皮的内在代谢改变。","",
    "铁死亡通路评分在GSE132465上皮伪bulk与GSE44076中显著上升，但GSE200997和GSE32323未显著。成员既有上升也有下降，例如铁摄取/转运、脂质重塑、抗氧化与铁储存相关成员方向不一。KEGG转录成员集不直接测量脂质过氧化或细胞死亡活性，因此目前只能说有部分、非稳定的转录关联线索，不能声称铁死亡增强或受RETSAT调控。","",
    "## 文件","",
    "- RETSAT_GSE132465_celltype_localization.png：GSE132465肿瘤细胞类型定位及RETSAT阳性率。",
    "- RETSAT_paired_cohort_signals.png：四队列患者配对表达与通路分数。",
    "- cohort_paired_summary.csv：四队列配对检验及BH校正。",
    "- paired_delta_correlations.csv：患者级RETSAT差值与通路分数差值相关。",
    "- pathway_member_paired_tests.csv：单个通路成员配对效应及FDR。",
    "- pathway_gene_coverage.csv：每个队列进入通路评分的有效基因数。",
    "- singlecell_epithelial_paired_cell_counts.csv、singlecell_cell_count_sensitivity.csv：单细胞配对细胞数及敏感性分析。",
    "- singlecell_pseudobulk_all_groups_*.csv、*_paired_expression_and_scores.csv：样本级伪bulk表达数据。",
    "- 分析脚本：analysis/retsat_crc_epithelial_paired_pathways.py；本次输出整理脚本：analysis/finalize_retsat_crc_outputs.py。","",
    "## 复现","",
    "在项目根目录依次运行：python analysis/retsat_crc_epithelial_paired_pathways.py，然后运行 python analysis/finalize_retsat_crc_outputs.py。前者从work中的本地原始矩阵重算；后者重建细胞数敏感性、有效通路覆盖、图和报告。","",
    "## 边界","",
    "GSE44076/GSE32323为bulk组织，细胞组成和分化状态差异可能影响代谢通路分数；单细胞患者数较少，GSE200997上皮亚群标签依赖本次的标志基因定义。RETSAT与通路变化相关不代表因果。GSE21510未纳入主bulk分析，因为原始元数据显示其肿瘤样本为LCM、配对正常为homogenized，取材处理区室不匹配。",""]
    (OUT/"RETSAT_CRC_epithelial_paired_pathways_report.md").write_text("\n".join(lines),encoding="utf-8")

def main():
    loc=regenerate_weighted_localization()
    summary=pd.read_csv(OUT/"cohort_paired_summary.csv")
    corr=pd.read_csv(OUT/"paired_delta_correlations.csv")
    gt=pd.read_csv(OUT/"pathway_member_paired_tests.csv")
    pairs=pd.read_csv(OUT/"all_cohorts_paired_sample_expression_scores.csv")
    coverage=create_coverage(gt)
    sens,counts=cellcount_sensitivity(pairs,gt)
    report(loc,summary,corr,gt,coverage,sens,counts)
    print("Wrote finalized report, coverage, cell-count sensitivity, and corrected localization figure.")

if __name__=="__main__": main()
