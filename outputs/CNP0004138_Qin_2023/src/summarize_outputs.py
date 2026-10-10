from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

ROOT=Path(__file__).resolve().parents[1]
RES=ROOT/'results'; FIG=ROOT/'figures'; REPORT=ROOT/'report'
FIG.mkdir(exist_ok=True); REPORT.mkdir(exist_ok=True)
audit=json.loads((RES/'scRNA_qc_analysis_audit.json').read_text(encoding='utf-8'))
paired=pd.read_csv(RES/'paired_patient_tests.csv')
co=pd.read_csv(RES/'post_patient_state_cooccurrence.csv')
spatial=pd.read_csv(RES/'spatial_neighbor_permutation.csv')
samples=pd.read_csv(RES/'sample_level_state_summaries.csv')
cell_path=RES/'cell_metadata_qc_labels.csv.gz'

# Sensitivity: remove genes shared between epithelial EMT and the CAF program being compared.
genes=audit['feature_genes_found']; gi={g:i for i,g in enumerate(genes)}; n=audit['n_cells_raw']
X=np.memmap(RES/'cache'/'sc_target_counts.float32.mmap',mode='r',dtype=np.float32,shape=(n,len(genes)))
meta=pd.read_csv(cell_path,usecols=['patient_id','stage','sample_id','qc_pass','cell_type'])
paired_ids=set(audit['paired_patient_ids'])
programs={
 'epi_EMT_no_shared':['VIM','FN1','ZEB1','ZEB2','SNAI1','SNAI2','SERPINE1','ITGA5','TGFBI','MMP7'],
 'CAF_FAP_no_MIR4435':['FAP','COL11A1','INHBA','POSTN','COL5A1','TGFB1'],
 'CAF_MMP1_no_CXCL8':['MMP1','MMP3','MMP10','CXCL1','CXCL3','TGFB1'],
}
sc={k:np.asarray(X[:,[gi[g] for g in v]].mean(axis=1),dtype=np.float32) for k,v in programs.items()}
post=meta[(meta.stage=='POSTT')&meta.patient_id.isin(paired_ids)&meta.qc_pass]
postrows=[]
for sid,df in post.groupby('sample_id',sort=False):
    ix=df.index.to_numpy(); epi=ix[df.cell_type.to_numpy()=='Epithelial_like']; caf=ix[df.cell_type.to_numpy()=='CAF_like']
    postrows.append({'patient_id':df.patient_id.iloc[0],
        'epi_EMT_no_shared':float(sc['epi_EMT_no_shared'][epi].mean()) if len(epi) else np.nan,
        'CAF_FAP_no_MIR4435':float(sc['CAF_FAP_no_MIR4435'][caf].mean()) if len(caf) else np.nan,
        'CAF_MMP1_no_CXCL8':float(sc['CAF_MMP1_no_CXCL8'][caf].mean()) if len(caf) else np.nan})
postpatient=pd.DataFrame(postrows).set_index('patient_id')
sens=[]
for xx,yy in [('epi_EMT_no_shared','CAF_FAP_no_MIR4435'),('epi_EMT_no_shared','CAF_MMP1_no_CXCL8')]:
    d=postpatient[[xx,yy]].dropna()
    rho,p=stats.spearmanr(d[xx],d[yy]) if len(d)>=5 else (np.nan,np.nan)
    sens.append({'tumor_program':xx,'CAF_state':yy,'n_paired_patients':len(d),'spearman_rho':rho,'p_value':p})
sens=pd.DataFrame(sens)
p=sens.p_value.to_numpy(); order=np.argsort(p); q=np.empty_like(p); q[order]=np.minimum.accumulate((p[order]*len(p)/(np.arange(len(p))+1))[::-1])[::-1]; sens['q_BH']=np.minimum(q,1)
sens.to_csv(RES/'overlap_free_cooccurrence_sensitivity.csv',index=False)

# CR/NR directional summary across the two available sections per group; no group-level p-value.
spatial['group']=spatial.sample_id.str.extract(r'^ST-(CR|NR)')[0]
groupdir=spatial.groupby(['neighbor_state','group']).difference.agg(['mean','count']).unstack('group')
rows=[]
for state in sorted(spatial.neighbor_state.unique()):
    d=spatial[spatial.neighbor_state.eq(state)]
    cr=d[d.group.eq('CR')].difference.to_numpy(); nr=d[d.group.eq('NR')].difference.to_numpy()
    rows.append({'neighbor_state':state,'mean_CR_difference':float(np.mean(cr)),'n_CR_sections':len(cr),'n_CR_positive':int(np.sum(cr>0)),'mean_NR_difference':float(np.mean(nr)),'n_NR_sections':len(nr),'n_NR_positive':int(np.sum(nr>0))})
gdir=pd.DataFrame(rows)
gdir.to_csv(RES/'spatial_group_direction_summary.csv',index=False)

# Effect-size heatmap (within-section neighborhood difference).
order_states=['CAF_BMP4','CAF_FAP','CAF_MMP1','CAF_PI16','CAF_SLIT2','CAF_ADAM28','MAC_NLRP3_like','MAC_C1QA_like','MAC_SPP1_like','EPI_EMT']
order_samples=['ST-CR1','ST-CR2','ST-NR1','ST-NR2']
mat=spatial.pivot(index='neighbor_state',columns='sample_id',values='difference').reindex(index=order_states,columns=order_samples)
fig,ax=plt.subplots(figsize=(8.8,6.2),constrained_layout=True)
lim=float(np.nanmax(np.abs(mat.to_numpy())))
im=ax.imshow(mat.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim)
ax.set_xticks(range(len(order_samples)),order_samples)
ax.set_yticks(range(len(order_states)),order_states)
ax.set_title('Neighborhood scores around epithelial-marker-high bins')
ax.set_xlabel('Post-treatment spatial sections (2 CR, 2 NR)')
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        val=mat.iloc[i,j]
        if pd.notna(val): ax.text(j,i,f'{val:.2f}',ha='center',va='center',fontsize=8,color='black')
fig.colorbar(im,ax=ax,label='Mean score in 8-bin neighbors: epithelial-high minus other bins')
fig.savefig(FIG/'spatial_neighbor_effect_heatmap.png',dpi=220); plt.close(fig)

# Look up paired statistics for narrative.
def row(metric):
    z=paired.loc[paired.metric.eq(metric)]
    return z.iloc[0] if len(z) else None
def fmt(metric):
    z=row(metric)
    if z is None: return '缺少估计'
    return f"配对差中位数 {z.median_post_minus_pre:.3f}，n={int(z.n_patients)}，Wilcoxon p={z.p_wilcoxon:.3g}，BH q={z.q_BH:.3g}"
qc_count=audit['n_qc_pass']; nraw=audit['n_cells_raw']; paper_count=333285
stage_n=samples.groupby('stage').n_cells_raw.sum().to_dict()
co_top=co.loc[co.tumor_program.eq('epi_EMT_mean')&co.myeloid_or_CAF_state.eq('CAF_FAP_mean')].iloc[0]
c1q_hyp=co.loc[co.tumor_program.eq('epi_Hypoxia_stress_mean')&co.myeloid_or_CAF_state.eq('myeloid_C1QA_like_mean')].iloc[0]
epi_nlrp=co.loc[co.tumor_program.eq('epi_Epithelial_identity_mean')&co.myeloid_or_CAF_state.eq('myeloid_NLRP3_like_mean')].iloc[0]
fap_sens=sens.loc[sens.CAF_state.eq('CAF_FAP_no_MIR4435')].iloc[0]
crrow=gdirdict=None
lines=[
'# CNP0004138 队列复用分析报告',
'',
'## 结论摘要',
'- scRNA 患者内配对显示，治疗后上皮样细胞比例没有可检出的系统性改变；但上皮身份、WNT/干性、增殖、缺氧应激和 EMT 程序分数总体下降，其中 IFN/ISG 分数未见改变。由于公开条码没有 XFZ 患者编号到 TRG 应答组的映射，以上不能回答这些残留上皮状态是否在 NR 特异富集，也不能把上皮样群直接称为 CNV 确认的残留恶性细胞。',
'- CAF 样群在配对术后样本中增加，CAF_PI16/SLIT2/ADAM28 程序上升；CAF_BMP4/FAP/MMP1 程序总体下降。髓系总比例变化不显著；NLRP3-like 均分向下但不显著，C1QA-like 向上但经多重校正未达显著，SPP1-like 分数上升。所有 scRNA 配对统计以患者为单位（27 对）。',
'- 27 位配对患者的术后样本级共丰度中，上皮 EMT 与 CAF_FAP 的 Spearman 相关为 rho=%.2f（原定义模块）；去掉两边重叠的 MIR4435-2HG 后仍为 rho=%.2f，p=%.2g。该结果支持患者间共变化，不证明细胞直接互作。' % (co_top.spearman_rho,fap_sens.spearman_rho,fap_sens.p_value),
'- 四张术后空间切片呈现 CR 与 NR 的方向性差异：两张 NR 切片的上皮身份高分 bin 周围，CAF_BMP4/FAP/MMP1、PI16/SLIT2、NLRP3-like 与 C1QA-like 等程序均高于该切片其他有效 bin；CR 切片中相应方向多为零或负。空间结果是每张切片内的局部邻近支持，每组只有两位患者，不构成患者层面的应答组推断。',
'',
'## 数据与统计单位',
'- 原文队列为 29 位患者，27 份治疗前、29 份治疗后，27 位可配对患者；有 TRG 的患者为 CR 7、PR 15、NR 5。空间组学只覆盖 4 个术后样本（CR 两张、NR 两张；无 PR）。原文将恶性上皮细胞通过 CNV 与聚类识别。来源：[Qin et al., Cell Reports Medicine (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10591051/)。',
'- CNGBdb 公布矩阵含 378,344 个条码，按患者×阶段识别出 61 个样本组：PRET 27、POSTT 29、POSTN 5；POSTN 是邻近正常样本，不参加肿瘤配对比较。MatrixMarket 中 159 个技术文库已合并到 61 个患者阶段样本。',
'- 应用论文所述检出基因数（200–6000）和线粒体比例（≤20%）阈值后，矩阵中 378,344 个条码全部通过（每细胞检出 201–5,999 个基因）；这比原文报告的 333,285 个高质量细胞多 45,059 个。公开矩阵与原文最终细胞集合未能完全对齐，可能涉及额外筛选/版本差异；本报告保留此差异，不宣称复现作者最终 QC 集合。',
'- 细胞仅用于估计每个样本的细胞比例与程序分数；统计检验的 n 是患者数，配对前后比较均使用一位患者一个 PRET 与一个 POSTT 值。治疗后共丰度相关也用患者级汇总（n=27）。',
'',
'## 1. 化疗后保留的上皮样状态',
'- 过滤后 POSTT 中识别出 36,114 个上皮样标记细胞（29 个治疗后样本）。这是表达标志定义的 epithelial-like 计数，不是 CNV 确认的残留癌细胞数。',
'- 上皮样细胞占样本全部 QC 细胞比例：'+fmt('frac_Epithelial_like')+'。',
'- 上皮身份分数：'+fmt('epi_Epithelial_identity_mean')+'。',
'- EMT：'+fmt('epi_EMT_mean')+'；WNT/干性：'+fmt('epi_WNT_stemness_mean')+'；增殖：'+fmt('epi_Proliferation_mean')+'；缺氧/应激：'+fmt('epi_Hypoxia_stress_mean')+'；IFN/ISG：'+fmt('epi_IFN_ISG_mean')+'。BH 校正范围为对应配对指标集合。',
'- 因 XFZ—TRG 对照缺失，无法判断上述任何细胞状态是否在 NR 比 CR/PR 更突出。不能用 ST-CR/NR 标签给 scRNA 病人反向补组。',
'',
'## 2. 髓系/巨噬细胞与 CAF 状态',
'- 髓系细胞总体比例：'+fmt('frac_Myeloid')+'。巨噬细胞候选程序的每样本平均分数：NLRP3-like（'+fmt('myeloid_NLRP3_like_mean')+'）；C1QA-like（'+fmt('myeloid_C1QA_like_mean')+'）；SPP1-like（'+fmt('myeloid_SPP1_like_mean')+'）。这是连续 marker-module 分数，不等同于已经验证的巨噬细胞亚群比例。',
'- CAF 样群比例：'+fmt('frac_CAF_like')+'。候选状态程序的配对变化：PI16（'+fmt('CAF_PI16_mean')+'）；SLIT2（'+fmt('CAF_SLIT2_mean')+'）；ADAM28（'+fmt('CAF_ADAM28_mean')+'）；BMP4（'+fmt('CAF_BMP4_mean')+'）；FAP（'+fmt('CAF_FAP_mean')+'）；MMP1（'+fmt('CAF_MMP1_mean')+'）。',
'- 术后患者级共丰度中，EMT—CAF_FAP 关系最强（原签名 rho=%.2f，q=%.2g）；剔除两边共用的 MIR4435-2HG 后仍为 rho=%.2f，p=%.2g。髓系信号较弱：上皮缺氧/应激与 C1QA-like 为 rho=%.2f、q=%.3f；上皮身份与 NLRP3-like 为 rho=%.2f、q=%.3f。其余相关见 `post_patient_state_cooccurrence.csv`；样本级共丰度不代表同一空间位点或因果作用。' % (co_top.spearman_rho,co_top.q_BH,fap_sens.spearman_rho,fap_sens.p_value,c1q_hyp.spearman_rho,c1q_hyp.q_BH,epi_nlrp.spearman_rho,epi_nlrp.q_BH),
'',
'## 3. 空间邻近',
'- 对每张切片，在 panel-positive 的 bins 中将上皮身份分数最高四分位定义为 epithelial-high 区；对每个这样的 bin 取最近 8 个 bin，计算候选髓系/CAF/EMT 程序均分，与该切片其余 panel-positive bins 的邻居分数比较。计数矩阵只解析固定的标志基因 panel；为统一 bin 深度，采用该 panel 的总计数作切片内归一化，不是全转录组 UMI 归一化。',
'- “epithelial-high”仅是上皮标志富集区，不代表存活肿瘤或癌细胞定位；空间测量单位是约 25 μm 的 bin（非单细胞分割）。邻近统计使用切片内标签置换；置换 p 值是 bin-level 局部关联，不是患者层面的 CR 对 NR 显著性检验。',
'- 方向在两张 NR 切片一致为正的候选关系：CAF_BMP4（NR 均差 %.3f、%.3f；CR %.3f、%.3f）、CAF_PI16（NR %.3f、%.3f；CR %.3f、%.3f）、CAF_SLIT2（NR %.3f、%.3f；CR %.3f、%.3f）、MAC_NLRP3_like（NR %.3f、%.3f；CR %.3f、%.3f）、MAC_C1QA_like（NR %.3f、%.3f；CR %.3f、%.3f）。CAF_FAP 和 CAF_MMP1 也在两张 NR 切片为正；热图显示每张切片的完整效应值。'.replace('%s',''),
'- 这种方向一致性可作为后续靶向复核的空间候选证据，但病例数仅 CR=2、NR=2，且没有术前切片，不可外推到整个应答组或全体患者。',
'',
'## 分析定义与可复现文件',
'- 上皮状态：预设 EMT、IFN/ISG、增殖、缺氧应激、WNT/干性、上皮身份基因集；分数为每细胞 CP10K 后 log1p 表达的基因均分，再汇总至样本/患者。',
'- 巨噬/髓系与 CAF 状态：依原文描述的 marker 及常见标记建立预设程序；属于候选程序复核，不重跑作者 Seurat 聚类、inferCNV、NMF 或原文图像解卷积。',
'- 配对检验：双侧 Wilcoxon signed-rank，针对一组样本级指标作 BH 校正；术后程序共丰度以 Spearman 相关并按相关矩阵作 BH 校正。',
'- 复跑：先执行 `src/analyze_reuse.py`，再执行 `src/summarize_outputs.py`。',
'- 主要文件：`results/paired_patient_tests.csv`、`results/sample_level_state_summaries.csv`、`results/post_patient_state_cooccurrence.csv`、`results/overlap_free_cooccurrence_sensitivity.csv`、`results/spatial_neighbor_permutation.csv`、`results/spatial_group_direction_summary.csv`、`results/cell_counts_audit.csv`、`results/scRNA_qc_analysis_audit.json`；图件在 `figures/`。',
'',
'## 补做 NR 分层所缺的数据',
'需要 Table S1 中 XFZ 编号到 CR/PR/NR 的对照表。作者补充 PDF 的公开下载在本次会话遇到 reCAPTCHA/Cloudflare 验证，我没有尝试绕过验证；用户确认手头没有映射表。因此本次不输出 scRNA 的 NR 特异统计。如果之后把 Table S1 PDF 或一张 XFZ—TRG 对照表放进该 D 盘项目目录，可以在 27 位配对患者上补做 NR/PR/CR 分层模型及患者级效应量。',
]
# Replace the spatial direction placeholder with values directly from per-section effects.
z=spatial
get=lambda state, sample: float(z.loc[(z.neighbor_state.eq(state))&(z.sample_id.eq(sample)),'difference'].iloc[0])
# Fill the single direction sentence with a compact exact value list.
placeholder=lines.index(next(x for x in lines if x.startswith('- 方向在两张 NR 切片一致为正')))
lines[placeholder]='- 上皮高分区邻居差值（NR1、NR2；随后 CR1、CR2）分别为：CAF_BMP4 %.3f, %.3f; %.3f, %.3f；CAF_PI16 %.3f, %.3f; %.3f, %.3f；CAF_SLIT2 %.3f, %.3f; %.3f, %.3f；NLRP3-like %.3f, %.3f; %.3f, %.3f；C1QA-like %.3f, %.3f; %.3f, %.3f。FAP、MMP1 在两张 NR 切片均为正（FAP %.3f/%.3f；MMP1 %.3f/%.3f），而 CR 接近零或为负。' % (
 get('CAF_BMP4','ST-NR1'),get('CAF_BMP4','ST-NR2'),get('CAF_BMP4','ST-CR1'),get('CAF_BMP4','ST-CR2'),
 get('CAF_PI16','ST-NR1'),get('CAF_PI16','ST-NR2'),get('CAF_PI16','ST-CR1'),get('CAF_PI16','ST-CR2'),
 get('CAF_SLIT2','ST-NR1'),get('CAF_SLIT2','ST-NR2'),get('CAF_SLIT2','ST-CR1'),get('CAF_SLIT2','ST-CR2'),
 get('MAC_NLRP3_like','ST-NR1'),get('MAC_NLRP3_like','ST-NR2'),get('MAC_NLRP3_like','ST-CR1'),get('MAC_NLRP3_like','ST-CR2'),
 get('MAC_C1QA_like','ST-NR1'),get('MAC_C1QA_like','ST-NR2'),get('MAC_C1QA_like','ST-CR1'),get('MAC_C1QA_like','ST-CR2'),
 get('CAF_FAP','ST-NR1'),get('CAF_FAP','ST-NR2'),get('CAF_MMP1','ST-NR1'),get('CAF_MMP1','ST-NR2'))
(REPORT/'analysis_report_zh.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('sensitivity',sens.to_dict(orient='records'))
print('spatial direction summary',gdir.to_dict(orient='records'))
print('report',REPORT/'analysis_report_zh.md')
