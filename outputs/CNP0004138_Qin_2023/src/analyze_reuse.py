from __future__ import annotations
import csv, gzip, json, re, sys, time, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'raw'
RES = ROOT / 'results'
FIG = ROOT / 'figures'
CACHE = RES / 'cache'
REPORT = ROOT / 'report'
for p in (RES, FIG, CACHE, REPORT): p.mkdir(parents=True, exist_ok=True)
SC = RAW / 'scRNA'
SP = RAW / 'spatial'

# Curated exploratory marker programs. Scores are mean log1p(CP10K) across detected panel genes.
LINEAGE = {
 'Epithelial_like':['EPCAM','KRT8','KRT18','KRT19','KRT20','CDH1','CEACAM5','CEACAM6','TACSTD2','AGR2','REG4'],
 'T_NK':['CD3D','CD3E','TRAC','CD8A','CD4','NKG7','CCL5','GNLY','GZMB'],
 'B':['MS4A1','CD79A','CD79B','CD22','CD37','CD74'],
 'Plasma':['JCHAIN','MZB1','SDC1','IGKC','IGHG1','IGHG3','IGHA1'],
 'Myeloid':['LST1','TYROBP','FCER1G','LYZ','CTSS','FCGR3A','CD14','CSF1R','CD68','AIF1'],
 'CAF_like':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','FAP','PDPN'],
 'Endothelial':['PECAM1','VWF','EMCN','CLDN5','PLVAP','KDR'],
 'Pericyte':['RGS5','PDGFRB','CSPG4','MCAM','NOTCH3'],
 'SMC':['ACTA2','TAGLN','MYH11','CNN1'],
 'Mast':['CPA3','TPSAB1','TPSB2','KIT'],
 'Glial':['S100B','SOX10','PLP1','GFAP','SLC1A3'],
}
EPITHELIAL = {
 'EMT':['MIR4435-2HG','VIM','FN1','ZEB1','ZEB2','SNAI1','SNAI2','SERPINE1','ITGA5','TGFBI','MMP7','CXCL8'],
 'IFN_ISG':['IFIT1','IFIT2','IFIT3','ISG15','MX1','IFI6','IFI27','STAT1','EPSTI1','IRF7'],
 'Proliferation':['MKI67','TOP2A','UBE2C','BIRC5','CCNB1','CDK1','STMN1','CENPF'],
 'Hypoxia_stress':['CA9','SLC2A1','LDHA','VEGFA','ADM','P4HA1','HIF1A','BNIP3','DDIT4'],
 'WNT_stemness':['ASCL2','LGR5','OLFM4','AXIN2','MYC','SOX9','RNF43','WNT7B','PROM1'],
 'Epithelial_identity':['EPCAM','KRT8','KRT18','KRT19','KRT20','CEACAM5','CEACAM6','TACSTD2','CDH1','MUC2','TFF3'],
}
MACROPHAGE = {
 'NLRP3_like':['NLRP3','IL1B','CCL3','CCL4','SOD2','CXCL8','PTGS2','NFKBIA','TNFAIP3'],
 'C1QA_like':['C1QA','C1QB','C1QC','APOE','MRC1','CD163','FOLR2','MARCO','VSIG4'],
 'SPP1_like':['SPP1','LGALS3','CTSB','TREM2','APOC1','FN1','GPNMB'],
 'FCN1_monocyte':['FCN1','S100A8','S100A9','VCAN','LILRB1','LILRA1'],
}
CAF = {
 'CAF_PI16':['PI16','CFD','C7','CD34','DPT','COL14A1','C3','CLEC3B'],
 'CAF_SLIT2':['SLIT2','PI16','CFD','DCN','C3','C7','CXCL12'],
 'CAF_ADAM28':['ADAM28','CFH','HLA-DRA','CD74','C7','C3'],
 'CAF_BMP4':['BMP4','WNT2','FRZB','FGF7','RSPO3','COL1A1','COL1A2'],
 'CAF_FAP':['FAP','COL11A1','INHBA','MIR4435-2HG','POSTN','COL5A1','TGFB1'],
 'CAF_MMP1':['MMP1','MMP3','MMP10','CXCL1','CXCL3','CXCL8','TGFB1'],
}
MODULES = {}
MODULES.update({f'LIN_{k}':v for k,v in LINEAGE.items()})
MODULES.update({f'EPI_{k}':v for k,v in EPITHELIAL.items()})
MODULES.update({f'MAC_{k}':v for k,v in MACROPHAGE.items()})
MODULES.update(CAF)
TARGETS = sorted(set(g for genes in MODULES.values() for g in genes))


def log(msg):
    print(time.strftime('%H:%M:%S'), msg, flush=True)

def bh(pvals):
    p=np.asarray(pvals,dtype=float)
    out=np.full(p.shape,np.nan)
    good=np.isfinite(p)
    if not good.any(): return out
    v=p[good]; order=np.argsort(v); ranked=v[order]
    q=ranked*len(v)/(np.arange(len(v))+1)
    q=np.minimum.accumulate(q[::-1])[::-1]
    back=np.empty_like(q); back[order]=q
    out[good]=np.minimum(back,1.0)
    return out

def parse_sc():
    genes_path=SC/'MM_coordinate_matrix_genes.csv'
    barcode_path=SC/'MM_coordinate_matrix_barcodes.csv'
    genes=pd.read_csv(genes_path,header=None,dtype=str)[0].fillna('').tolist()
    bar=pd.read_csv(barcode_path,header=None,dtype=str,names=['barcode'])
    ncell=len(bar); ngene=len(genes)
    ids=bar.barcode.astype(str)
    extracted=ids.str.extract(r'^(?P<patient>XFZ-\d+)-(?P<stage>PRET|POSTT|POSTN)-lib-(?P<lib>\d+)')
    if extracted.patient.isna().any():
        raise ValueError('Could not parse patient/stage from all barcodes')
    bar['patient_id']=extracted.patient.values
    bar['stage']=extracted.stage.values
    bar['sample_id']=extracted.patient+'-'+extracted.stage
    log(f'barcodes={ncell:,}; genes={ngene:,}; patients={bar.patient_id.nunique()}; samples={bar.sample_id.nunique()}')
    # Gene name -> first matching feature row. The source feature list is one gene symbol per row.
    gene_to_row={}
    for i,g in enumerate(genes):
        gene_to_row.setdefault(g,i)
    found={g:gene_to_row[g] for g in TARGETS if g in gene_to_row}
    feature_names=list(found)
    row_to_feat=np.full(ngene,-1,dtype=np.int32)
    for j,g in enumerate(feature_names): row_to_feat[found[g]]=j
    mt_rows=np.array([i for i,g in enumerate(genes) if g.upper().startswith('MT-')],dtype=np.int32)
    mt_lookup=np.zeros(ngene,dtype=bool); mt_lookup[mt_rows]=True
    features=np.memmap(CACHE/'sc_target_counts.float32.mmap',mode='w+',dtype=np.float32,shape=(ncell,len(feature_names)))
    features[:]=0
    totals=np.zeros(ncell,dtype=np.float64)
    nfeat=np.zeros(ncell,dtype=np.int32)
    mtcounts=np.zeros(ncell,dtype=np.float64)
    mtx=SC/'MM_coordinate_matrix.mtx.gz'
    with gzip.open(mtx,'rb') as f:
        head=[f.readline() for _ in range(3)]
        dims=[int(x) for x in head[2].split()]
        if dims[0]!=ngene or dims[1]!=ncell: raise ValueError(f'Matrix dimensions mismatch {dims[:2]} vs {(ngene,ncell)}')
        nnz=dims[2]
    reader=pd.read_csv(mtx,sep=r'\s+',skiprows=3,header=None,dtype={0:np.int32,1:np.int32,2:np.float32},chunksize=3_000_000,engine='c',compression='gzip')
    processed=0; t0=time.time(); next_report=10_000_000
    for chunk in reader:
        rows=chunk.iloc[:,0].to_numpy(copy=False)-1
        cols=chunk.iloc[:,1].to_numpy(copy=False)-1
        vals=chunk.iloc[:,2].to_numpy(copy=False)
        totals += np.bincount(cols,weights=vals,minlength=ncell)
        nfeat += np.bincount(cols,minlength=ncell).astype(np.int32)
        if mt_rows.size:
            mask=mt_lookup[rows]
            if mask.any(): mtcounts += np.bincount(cols[mask],weights=vals[mask],minlength=ncell)
        fpos=row_to_feat[rows]
        keep=fpos>=0
        if keep.any(): features[cols[keep],fpos[keep]]=vals[keep]
        processed+=len(vals)
        if processed>=next_report:
            log(f'mtx entries {processed:,}/{nnz:,} ({100*processed/nnz:.1f}%), elapsed {time.time()-t0:.0f}s')
            next_report+=10_000_000
        del chunk,rows,cols,vals,fpos,keep
    if processed!=nnz: raise ValueError(f'MatrixMarket expected {nnz:,} entries; read {processed:,}')
    features.flush()
    log(f'count stream complete ({time.time()-t0:.0f}s); marker genes present {len(feature_names)}/{len(TARGETS)}')
    if len(mt_rows)==0: log('WARNING: no MT- genes found; QC mt% will be NA')
    # Gene-wise CP10K normalization in place, keeping memory bounded.
    denom=np.maximum(totals,1.0)
    for j in range(len(feature_names)):
        features[:,j]=np.log1p(features[:,j]/denom*10000).astype(np.float32)
    features.flush()
    np.savez_compressed(CACHE/'cell_qc.npz',total_umi=totals,n_genes=nfeat,mt_umi=mtcounts)
    genes_found=set(feature_names)
    def score(mod):
        use=[feature_names.index(g) for g in mod if g in genes_found]
        if not use: return np.zeros(ncell,dtype=np.float32),0
        return np.asarray(features[:,use].mean(axis=1),dtype=np.float32),len(use)
    module_scores={}; module_genecount={}
    for name,mod in MODULES.items():
        module_scores[name],module_genecount[name]=score(mod)
    lineage_names=list(LINEAGE)
    lmat=np.column_stack([module_scores['LIN_'+k] for k in lineage_names])
    best=np.argmax(lmat,axis=1); bestscore=lmat[np.arange(ncell),best]
    lineage=np.array([lineage_names[i] for i in best],dtype=object)
    lineage[bestscore<0.22]='Unassigned'
    qc=(nfeat>=200)&(nfeat<=6000)&((mtcounts/denom*100)<=20 if mt_rows.size else True)
    bar['qc_pass']=qc
    bar['n_genes']=nfeat
    bar['total_umi']=totals.astype(np.int64)
    bar['mt_percent']=np.where(denom>0,100*mtcounts/denom,np.nan)
    bar['cell_type']=lineage
    for name,vals in module_scores.items(): bar[name]=vals
    # Sample-level summary: each sample is the experimental observation; cells contribute to its estimates only.
    rows=[]
    indexes=bar.groupby('sample_id',sort=False).indices
    for sample,ix in indexes.items():
        ix=np.asarray(ix,dtype=np.int64)
        keep=ix[qc[ix]]
        rec={'sample_id':sample,'patient_id':bar.patient_id.iloc[ix[0]],'stage':bar.stage.iloc[ix[0]],'n_cells_raw':len(ix),'n_cells_qc':len(keep),'qc_fraction':len(keep)/max(1,len(ix))}
        subtypes=lineage[keep]
        for typ in list(LINEAGE)+['Unassigned']:
            rec['frac_'+typ]=float(np.mean(subtypes==typ)) if len(keep) else np.nan
        eix=keep[subtypes=='Epithelial_like']
        mix=keep[subtypes=='Myeloid']
        cix=keep[subtypes=='CAF_like']
        rec['n_epithelial_like']=len(eix); rec['n_myeloid']=len(mix); rec['n_CAF_like']=len(cix)
        for k in EPITHELIAL:
            rec['epi_'+k+'_mean']=float(module_scores['EPI_'+k][eix].mean()) if len(eix) else np.nan
        for k in MACROPHAGE:
            rec['myeloid_'+k+'_mean']=float(module_scores['MAC_'+k][mix].mean()) if len(mix) else np.nan
        for k in CAF:
            rec[k+'_mean']=float(module_scores[k][cix].mean()) if len(cix) else np.nan
        rows.append(rec)
    samples=pd.DataFrame(rows)
    samples.to_csv(RES/'sample_level_state_summaries.csv',index=False)
    bar[['barcode','patient_id','stage','sample_id','qc_pass','n_genes','total_umi','mt_percent','cell_type']].to_csv(RES/'cell_metadata_qc_labels.csv.gz',index=False,compression='gzip')
    # Cell type and stage audit counts.
    audit=bar.groupby(['stage','qc_pass','cell_type'],observed=True).size().rename('n_cells').reset_index()
    audit.to_csv(RES/'cell_counts_audit.csv',index=False)
    # Matched within-patient pre/post tests; patients, not cells, are the observations.
    paired=[]
    paired_ids=set(samples.loc[samples.stage.eq('PRET'),'patient_id'])&set(samples.loc[samples.stage.eq('POSTT'),'patient_id'])
    paired_s=samples[samples.patient_id.isin(paired_ids)&samples.stage.isin(['PRET','POSTT'])]
    wide=paired_s.set_index(['patient_id','stage'])
    metrics=[c for c in samples.columns if c.startswith('frac_') or c.startswith('epi_') or c.startswith('myeloid_') or c.startswith('CAF_')]
    for metric in metrics:
        pre=[]; post=[]
        for patient in sorted(paired_ids):
            try: a=wide.loc[(patient,'PRET'),metric]; b=wide.loc[(patient,'POSTT'),metric]
            except KeyError: continue
            if pd.notna(a) and pd.notna(b): pre.append(float(a)); post.append(float(b))
        pre=np.asarray(pre); post=np.asarray(post)
        if len(pre)<3: pval=np.nan
        elif np.allclose(post-pre,0): pval=1.0
        else:
            try: pval=float(stats.wilcoxon(post,pre,alternative='two-sided',zero_method='wilcox',method='auto').pvalue)
            except ValueError: pval=np.nan
        paired.append({'metric':metric,'n_patients':len(pre),'median_pre':np.nanmedian(pre) if len(pre) else np.nan,'median_post':np.nanmedian(post) if len(post) else np.nan,'median_post_minus_pre':np.nanmedian(post-pre) if len(pre) else np.nan,'n_increased':int(np.sum(post>pre)),'n_decreased':int(np.sum(post<pre)),'p_wilcoxon':pval})
    paired_df=pd.DataFrame(paired); paired_df['q_BH']=bh(paired_df.p_wilcoxon.to_numpy())
    paired_df.to_csv(RES/'paired_patient_tests.csv',index=False)
    # Post-treatment co-occurrence, one summary per paired patient.
    post=samples[samples.patient_id.isin(paired_ids)&samples.stage.eq('POSTT')].set_index('patient_id')
    co=[]
    tumor_metrics=['epi_'+k+'_mean' for k in EPITHELIAL]
    neighbor_metrics=['myeloid_'+k+'_mean' for k in MACROPHAGE]+[k+'_mean' for k in CAF]
    for tm in tumor_metrics:
        for nm in neighbor_metrics:
            xy=post[[tm,nm]].dropna()
            if len(xy)<5 or xy[tm].nunique()<2 or xy[nm].nunique()<2:
                rho=p=np.nan
            else: rho,p=stats.spearmanr(xy[tm],xy[nm])
            co.append({'tumor_program':tm,'myeloid_or_CAF_state':nm,'n_patients':len(xy),'spearman_rho':rho,'p_value':p})
    co=pd.DataFrame(co); co['q_BH']=bh(co.p_value.to_numpy()); co.to_csv(RES/'post_patient_state_cooccurrence.csv',index=False)
    # Output precise mapping-independent limitation.
    qc_summary={'n_cells_raw':int(ncell),'n_genes':int(ngene),'n_matrix_nonzero':int(nnz),'n_cell_barcodes':int(ncell),'n_patients':int(bar.patient_id.nunique()),'n_samples':int(bar.sample_id.nunique()),'n_qc_pass':int(qc.sum()),'n_qc_fail':int((~qc).sum()),'paired_patient_ids':sorted(paired_ids),'n_paired_patients':len(paired_ids),'feature_genes_found':feature_names,'feature_gene_missing':[x for x in TARGETS if x not in genes_found],'genes_total':ngene,'mt_gene_count':int(len(mt_rows))}
    (RES/'scRNA_qc_analysis_audit.json').write_text(json.dumps(qc_summary,indent=2,ensure_ascii=False),encoding='utf-8')
    log(f'QC pass {int(qc.sum()):,}/{ncell:,}; paired scRNA patients={len(paired_ids)}; sample metrics saved')
    return samples,bar,features,feature_names,module_scores,qc,paired_df,co


def parse_spatial(sample,feature_names):
    dat=SP/f'{sample}.exp.count.tsv'
    if not dat.exists(): dat=SP/f'{sample}.exp.count.tsv.gz'
    if not dat.exists() or dat.stat().st_size<1000: raise FileNotFoundError(dat)
    meta=pd.read_csv(SP/f'{sample}.meta.info.tsv.txt',sep='\t')
    ids=meta.iloc[:,0].astype(str).to_numpy()
    # Some distributed tables name the first bin column ''; count data has a wide gene x bin table.
    op=gzip.open if dat.suffix=='.gz' else open
    with op(dat,'rt',encoding='utf-8',newline='') as f:
        header=f.readline().rstrip('\r\n').split('\t')
        bins=header[1:]
        if len(bins)!=len(ids): raise ValueError(f'{sample}: count bins {len(bins)} != metadata {len(ids)}')
        order=pd.Index(bins).get_indexer(ids)
        if (order<0).any(): raise ValueError(f'{sample}: metadata bin IDs do not match count columns')
        use=set(feature_names)
        feat=np.zeros((len(bins),len(feature_names)),dtype=np.float32)
        total=np.zeros(len(bins),dtype=np.float64)
        fmap={g:i for i,g in enumerate(feature_names)}
        nrow=0; t0=time.time()
        for line in f:
            k=line.find('\t')
            if k<0: continue
            gene=line[:k]
            if gene in use:
                vals=np.fromstring(line[k+1:],sep='\t',dtype=np.float32)
                if len(vals)!=len(bins): raise ValueError(f'{sample}: row {gene} has {len(vals)} bins; expected {len(bins)}')
                feat[:,fmap[gene]]=vals
            nrow+=1
            if nrow%10000==0: log(f'{sample}: scanned gene labels {nrow:,}; elapsed {time.time()-t0:.0f}s')
    # Reorder from count matrix column order to metadata order.
    feat=feat[order]; total=total[order]
    # Within-section library factor from the fixed targeted marker panel; all bins use the same panel.
    total=feat.sum(axis=1,dtype=np.float64)
    denom=np.maximum(total,1.0)
    for j in range(feat.shape[1]): feat[:,j]=np.log1p(feat[:,j]/denom*10000).astype(np.float32)
    def smod(mod):
        js=[feature_names.index(g) for g in mod if g in feature_names]
        if not js: return np.zeros(len(ids),dtype=np.float32)
        return feat[:,js].mean(axis=1)
    scores={}
    for k,v in LINEAGE.items(): scores['LIN_'+k]=smod(v)
    for k,v in EPITHELIAL.items(): scores['EPI_'+k]=smod(v)
    for k,v in MACROPHAGE.items(): scores['MAC_'+k]=smod(v)
    for k,v in CAF.items(): scores[k]=smod(v)
    return meta,ids,scores,total


def spatial_analysis(feature_names):
    spatial_samples=['ST-CR1','ST-CR2','ST-NR1','ST-NR2']
    out=[]; mapdata={}
    state_keys=['MAC_NLRP3_like','MAC_C1QA_like','MAC_SPP1_like','CAF_PI16','CAF_SLIT2','CAF_ADAM28','CAF_BMP4','CAF_FAP','CAF_MMP1']
    rng=np.random.default_rng(4138)
    for sample in spatial_samples:
        log(f'begin {sample}')
        meta,ids,scores,total=parse_spatial(sample,feature_names)
        x=meta['x'].to_numpy(float); y=meta['y'].to_numpy(float)
        coords=np.column_stack([x,y]); valid=np.isfinite(total)&(total>0)
        # Epithelial-high spatial bins are a proxy for tumor-enriched regions, not a CNV call.
        epi=scores['EPI_Epithelial_identity']
        cut=np.quantile(epi[valid],0.75) if valid.any() else np.inf
        high=valid&(epi>=cut)
        k=min(9,len(ids)); _,nei=cKDTree(coords).query(coords,k=k)
        if nei.ndim==1: nei=nei[:,None]
        # Exclude the focal bin (first nearest neighbour); average the nearest eight bins.
        neigh_idx=nei[:,1:] if nei.shape[1]>1 else nei
        states={key:scores[key] for key in state_keys}
        local={key:v[neigh_idx].mean(axis=1) for key,v in states.items()}
        local_epi_emt=scores['EPI_EMT'][neigh_idx].mean(axis=1)
        # Permute tumor-bin labels within each tissue section; this assesses local association only.
        n=int(high.sum()); valid_ix=np.flatnonzero(valid); nvalid=len(valid_ix)
        for key,arr in list(local.items())+[('EPI_EMT',local_epi_emt)]:
            if n<3 or n>=nvalid-2:
                obs=np.nan; p=np.nan; mh=ml=np.nan
            else:
                mh=float(arr[high].mean()); ml=float(arr[valid&~high].mean()); obs=mh-ml
                vals=arr[valid_ix]
                ge=0; permN=999
                # batch independent random labels to limit transient allocation.
                for start in range(0,permN,100):
                    b=min(100,permN-start)
                    random_ix=np.argpartition(rng.random((b,nvalid)),n-1,axis=1)[:,:n]
                    means=vals[random_ix].mean(axis=1)
                    others=(vals.sum()-means*n)/(nvalid-n)
                    dif=means-others
                    ge+=int(np.sum(np.abs(dif)>=abs(obs)))
                p=(ge+1)/(permN+1)
            out.append({'sample_id':sample,'response_group_from_ST_label':sample.split('-')[1][:2],'index_region':'epithelial_identity_top_quartile','neighbor_state':key,'n_bins':len(ids),'n_valid_bins':nvalid,'n_epithelial_high_bins':n,'mean_neighbor_score_epi_high':mh,'mean_neighbor_score_other':ml,'difference':obs,'p_spatial_permutation':p})
        # Save compact coordinates/score data for maps.
        mapdata[sample]=(x,y,epi,local['CAF_FAP'],local['CAF_MMP1'],local['CAF_BMP4'],local_epi_emt,total)
        log(f'{sample} complete: bins={len(ids):,}, epi-high={n:,}')
    result=pd.DataFrame(out); result['q_BH']=bh(result.p_spatial_permutation.to_numpy())
    result.to_csv(RES/'spatial_neighbor_permutation.csv',index=False)
    # Figure: one row per section, epithelial identity and nearby CAF_FAP signal.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(4,2,figsize=(11,18),constrained_layout=True)
    for i,sample in enumerate(spatial_samples):
        x,y,epi,lfap,lmmp1,lbmp4,emt,total=mapdata[sample]
        good=total>0
        sc=axs[i,0].scatter(x[good],y[good],c=epi[good],s=2,cmap='viridis',rasterized=True)
        axs[i,0].set_title(f'{sample}: epithelial identity')
        fig.colorbar(sc,ax=axs[i,0],shrink=.75)
        sc2=axs[i,1].scatter(x[good],y[good],c=lfap[good],s=2,cmap='magma',rasterized=True)
        axs[i,1].set_title(f'{sample}: nearby CAF_FAP score')
        fig.colorbar(sc2,ax=axs[i,1],shrink=.75)
        for ax in axs[i,:]: ax.set_aspect('equal'); ax.set_xlabel('x'); ax.set_ylabel('y'); ax.invert_yaxis()
    fig.savefig(FIG/'spatial_epithelial_CAFFAP_neighbors.png',dpi=200); plt.close(fig)
    return result


def make_figures(samples,paired_df,co):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    paired=samples[samples.stage.isin(['PRET','POSTT'])].copy()
    pairs=paired.pivot(index='patient_id',columns='stage')
    panels=['frac_Epithelial_like','epi_EMT_mean','epi_IFN_ISG_mean','epi_Proliferation_mean','epi_Hypoxia_stress_mean','epi_WNT_stemness_mean','CAF_FAP_mean','CAF_BMP4_mean','myeloid_C1QA_like_mean','myeloid_NLRP3_like_mean']
    fig,axs=plt.subplots(2,5,figsize=(16,7),constrained_layout=True)
    for ax,m in zip(axs.flat,panels):
        x=paired.pivot(index='patient_id',columns='stage',values=m).dropna()
        for _,row in x.iterrows(): ax.plot([0,1],[row.get('PRET',np.nan),row.get('POSTT',np.nan)],color='#75808c',alpha=.32,lw=.8)
        if len(x):
            ax.scatter(np.zeros(len(x)),x['PRET'],s=13,color='#3b6fb6',label='PRET',zorder=2)
            ax.scatter(np.ones(len(x)),x['POSTT'],s=13,color='#d26a4a',label='POSTT',zorder=2)
        ax.set_xticks([0,1],['Pre','Post']); ax.set_title(m.replace('_',' '),fontsize=9); ax.grid(axis='y',alpha=.2)
        if m=='frac_Epithelial_like': ax.set_ylabel('fraction of QC cells')
    axs[0,0].legend(frameon=False,fontsize=8)
    fig.suptitle('CNP0004138: paired patient-level sample summaries (27 matched patients)')
    fig.savefig(FIG/'paired_patient_state_summaries.png',dpi=220); plt.close(fig)


def write_report(samples,bar,paired_df,co,spatial_df):
    paired_ids=sorted(set(samples.loc[samples.stage.eq('PRET'),'patient_id'])&set(samples.loc[samples.stage.eq('POSTT'),'patient_id']))
    qcpass=int(bar.qc_pass.sum()); n=int(len(bar)); nsp=len(spatial_df.sample_id.unique()) if spatial_df is not None and len(spatial_df) else 0
    lines=[
      '# CNP0004138 队列复用分析', '',
      '## 研究范围与数据',
      f'- scRNA 处理矩阵：{n:,} 个条码、{bar.patient_id.nunique()} 位患者；依据原文 QC 阈值（200–6000 检出基因、线粒体比例 ≤20%）筛得 {qcpass:,} 个细胞。',
      f'- 识别到 {len(paired_ids)} 位治疗前后配对患者；下游配对检验以患者为统计单位，每位患者每阶段先汇总成一个样本级数值。',
      '- 邻近正常组织样本（POSTN）不纳入肿瘤前后配对比较。',
      '- 状态标签由公开表达矩阵中的预定义标志基因程序评分得到，属于候选状态复核；未对全矩阵重跑原文的 Seurat 聚类、inferCNV 或 NMF。Epithelial_like 是上皮样群，不能直接等同于 CNV 确认的恶性细胞。',
      '', '## 主要限制：目前不能检验 scRNA 的 NR 特异性',
      '- 条码可解析患者编号与治疗阶段，但不含 CR/PR/NR。用户没有额外的编号对照表。原文完整患者—TRG 对照位于 Table S1；官方 PMC 下载链接目前被 reCAPTCHA/Cloudflare 阻断，因此本分析未绕过验证，也没有把无法确认的患者分组。',
      '- 因此，所有 scRNA 结果仅回答“配对治疗前后如何变化”和“患者级治疗后状态间是否共变”；不能把配对前后差异说成 NR 特异富集。',
      f'- 空间组学包含 {nsp} 张术后切片（ST-CR1/2、ST-NR1/2；CR、NR 各 2 张），没有 PR；样本内空间置换评估只是局部邻近支持，不是患者总体推断（每组 n=2）。',
      '', '## 分析问题与结果文件',
      '1. 上皮样残留状态：看 `paired_patient_tests.csv` 中 `frac_Epithelial_like` 与 `epi_*_mean`；前者是 QC 细胞中的上皮样比例，后者为上皮样细胞内状态程序均分。',
      '2. 髓系/巨噬细胞与 CAF 共现：`post_patient_state_cooccurrence.csv` 用配对患者术后汇总值做 Spearman 相关；是患者层面的共丰度，不代表同一细胞或直接互作。巨噬细胞标志程序与原有 NLRP3-like、C1QA-like 观察并列。',
      '3. 空间邻近：`spatial_neighbor_permutation.csv` 比较上皮身份高四分位空间 bin 周围最近 8 个 bin 的巨噬细胞/CAF 程序分数与该切片其他 bin，并在每张切片内打乱上皮高分 bin 标签。高分 bin 只是肿瘤富集代理指标，不能区分正常上皮与 CNV+肿瘤细胞。',
      '', '## 文件索引',
      '- `sample_level_state_summaries.csv`：每个样本一个样本级汇总。',
      '- `paired_patient_tests.csv`：配对 Wilcoxon、患者数、配对差中位数与 BH 校正。',
      '- `post_patient_state_cooccurrence.csv`：术后患者级状态共丰度相关。',
      '- `spatial_neighbor_permutation.csv`：4 张切片的空间邻近差异与切片内置换 p 值。',
      '- `cell_counts_audit.csv`、`scRNA_qc_analysis_audit.json`：QC 和状态分类审计。',
      '- `cell_metadata_qc_labels.csv.gz`：细胞级条码、QC 与粗略谱系标签；不含反应组标签。',
      '- `figures/paired_patient_state_summaries.png`、`figures/spatial_epithelial_CAFFAP_neighbors.png`：配对变化及空间分布图。',
      '- `src/analyze_reuse.py`：完整可复跑脚本。',
      '', '## 解释边界',
      '1. 细胞数仅用于估计样本比例与样本内程序分数；统计自由度是患者数。稀有细胞群会因样本内细胞数不同而不稳定。',
      '2. 人为定义的基因程序依赖 marker panel，不能替代原文完整聚类和恶性 CNV 注释。',
      '3. 空间置换保留单张切片中的基因表达和邻域结构，但不解决每组只有两位患者、切片选择偏倚或上皮 bin 肿瘤纯度问题。',
      '4. 获得 Table S1 中 XFZ 患者编号到 TRG 的映射后，才能补做 CR/PR/NR 分层的配对变化和 NR 特异性检验。',
    ]
    (ROOT/'report'/'analysis_report_zh.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

if __name__=='__main__':
    samples,bar,features,feature_names,module_scores,qc,paired_df,co=parse_sc()
    make_figures(samples,paired_df,co)
    log('Beginning spatial analysis (four sections, streaming gene rows)')
    spatial_df=spatial_analysis(feature_names)
    write_report(samples,bar,paired_df,co,spatial_df)
    log('Analysis complete; outputs in '+str(ROOT))




