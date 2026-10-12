from __future__ import annotations
import csv, gzip, re, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = Path(r"D:\DINP\whynot17-new-project")
WORK = PROJECT / "work"
SC_DIR = WORK / "single_cell_GSE132465_GSE200997"
GEO_DIR = WORK / "geo_cache"
OUT = PROJECT / "results" / "2026-10-12" / "RETSAT_CRC_epithelial_paired_pathways"
RETINOL_ID, FERROPTOSIS_ID = "hsa00830", "hsa04216"
EP_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20"]

def kegg_gene_set(pathway_id: str, outdir: Path) -> set[str]:
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    with urllib.request.urlopen(url, timeout=60) as response:
        txt = response.read().decode("utf-8")
    (outdir / f"{pathway_id}_KEGG_flatfile.txt").write_text(txt, encoding="utf-8")
    genes, in_genes = set(), False
    for line in txt.splitlines():
        if line.startswith("GENE"):
            in_genes, payload = True, line[4:]
        elif line.startswith("///"):
            break
        elif in_genes:
            payload = line
        else:
            continue
        match = re.match(r"\s*(\d+)\s+([A-Z0-9-]+);", payload)
        if match:
            genes.add(match.group(2).upper())
    return genes

def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any(): return out
    pv = p[ok]; order = np.argsort(pv); ranked = pv[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted = np.empty_like(q); adjusted[order] = np.minimum(q, 1.0)
    out[ok] = adjusted
    return out

def safe_wilcoxon(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    valid = np.isfinite(a) & np.isfinite(b); a, b = a[valid], b[valid]
    if len(a) < 2: return np.nan
    if np.allclose(a-b, 0): return 1.0
    try: return float(wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
    except Exception: return np.nan

def read_sc_annotation(dataset: str) -> pd.DataFrame:
    if dataset == "GSE132465":
        p = SC_DIR / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
        a = pd.read_csv(p, sep="\t", compression="gzip")
        return a.rename(columns={"Index":"cell_id","Patient":"patient","Class":"condition","Sample":"sample","Cell_type":"cell_type","Cell_subtype":"subtype"})
    p = SC_DIR / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"
    a = pd.read_csv(p, compression="gzip")
    a = a.rename(columns={a.columns[0]:"cell_id","samples":"sample","Condition":"condition"})
    a["patient"] = a["sample"].str.extract(r"cac(\d+)", expand=False).fillna(a["sample"])
    return a

def matrix_header(path: Path, delimiter: str):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        fields = next(csv.reader([f.readline()], delimiter=delimiter))
    return [x.strip().strip('"') for x in fields[1:]]

def read_selected_cell_rows(path: Path, delimiter: str, wanted: set[str]):
    cell_ids = matrix_header(path, delimiter); found = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        next(f)
        for line in f:
            gene, payload = line.rstrip("\r\n").split(delimiter, 1)
            gene = gene.strip().strip('"').upper()
            if gene in wanted:
                found[gene] = np.fromstring(payload, sep=delimiter, dtype=np.int32)
                if len(found) == len(wanted): break
    missing = wanted - set(found)
    if missing: print(f"WARNING absent genes in {path.name}: {sorted(missing)}")
    return cell_ids, found

def define_gse200997_epithelium(ann: pd.DataFrame):
    p = SC_DIR / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
    cell_ids, rows = read_selected_cell_rows(p, ",", set(EP_MARKERS))
    if ann["cell_id"].astype(str).tolist() != cell_ids:
        ann = ann.set_index("cell_id").loc[cell_ids].reset_index()
    mat = np.column_stack([rows.get(g, np.zeros(len(cell_ids), dtype=np.int32)) for g in EP_MARKERS])
    ann["epithelial_marker_count"] = (mat > 0).sum(axis=1)
    ann["cell_type"] = np.where(ann["epithelial_marker_count"] >= 2, "epithelial_like", "non_epithelial")
    return ann

def scan_single_cell(dataset: str, genesets: dict[str,set[str]], outdir: Path):
    ann = read_sc_annotation(dataset)
    if dataset == "GSE200997":
        ann = define_gse200997_epithelium(ann)
        matrix, delimiter = SC_DIR/"GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz", ","
    else:
        matrix, delimiter = SC_DIR/"GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz", "\t"
    cell_ids = matrix_header(matrix, delimiter)
    if ann["cell_id"].astype(str).tolist() != cell_ids:
        ann = ann.set_index("cell_id").loc[cell_ids].reset_index()
    group_cols = ["sample","patient","condition","cell_type"]
    ann[group_cols] = ann[group_cols].astype(str)
    codes, uniques = pd.factorize(pd.MultiIndex.from_frame(ann[group_cols]), sort=True)
    gm = pd.DataFrame(list(uniques), columns=group_cols)
    gm["n_cells"] = np.bincount(codes, minlength=len(gm))
    targets = {"RETSAT"}
    for gs in genesets.values(): targets.update(g.upper() for g in gs)
    gidx = {g:i for i,g in enumerate(sorted(targets))}
    sums = np.zeros((len(gidx),len(gm)),dtype=np.float64)
    totals = np.zeros(len(gm),dtype=np.float64)
    ret_detect = np.zeros(len(gm),dtype=np.float64)
    with gzip.open(matrix,"rt",encoding="utf-8",newline="") as f:
        next(f)
        for line_no,line in enumerate(f,start=2):
            gene,payload=line.rstrip("\r\n").split(delimiter,1)
            gene=gene.strip().strip('"').upper()
            values=np.fromstring(payload,sep=delimiter,dtype=np.int32)
            if len(values)!=len(cell_ids):
                raise ValueError(f"{dataset} malformed row {line_no}: {gene} {len(values)} != {len(cell_ids)}")
            totals += np.bincount(codes,weights=values,minlength=len(gm))
            if gene in gidx:
                sums[gidx[gene],:] += np.bincount(codes,weights=values,minlength=len(gm))
                if gene=="RETSAT":
                    ret_detect += np.bincount(codes,weights=(values>0),minlength=len(gm))
    ret_counts=sums[gidx["RETSAT"],:]
    ret_logcpm=np.log2(ret_counts/np.maximum(totals,1)*1e6+0.5)
    extra={"total_umis":totals,"RETSAT_counts":ret_counts,"RETSAT_positive_cells":ret_detect,
           "RETSAT_detection_pct":100*ret_detect/gm["n_cells"],"RETSAT_log2CPM":ret_logcpm}
    expr_cols={"RETSAT":ret_logcpm}
    for g,i in gidx.items():
        if g=="RETSAT": continue
        logcpm=np.log2(sums[i,:]/np.maximum(totals,1)*1e6+0.5)
        extra[f"count::{g}"]=sums[i,:]
        extra[f"log2CPM::{g}"]=logcpm
        expr_cols[g]=logcpm
    extra["dataset"]=np.repeat(dataset,len(gm))
    gm=pd.concat([gm,pd.DataFrame(extra)],axis=1)
    gm.to_csv(outdir/f"singlecell_pseudobulk_all_groups_{dataset}.csv",index=False)
    loc=gm[gm["condition"].eq("Tumor")].groupby("cell_type").agg(
        tumor_cells=("n_cells","sum"),patients=("patient","nunique"),
        positive_cells=("RETSAT_positive_cells","sum"),
        median_patient_log2CPM=("RETSAT_log2CPM","median"),
        mean_group_detection_pct=("RETSAT_detection_pct","mean"),
        groups_ge10=("n_cells",lambda s:int((s>=10).sum()))
    ).reset_index()
    if loc["positive_cells"].sum()>0:
        loc["pct_of_tumor_RETSAT_positive_cells"]=100*loc["positive_cells"]/loc["positive_cells"].sum()
    loc["dataset"]=dataset
    loc.to_csv(outdir/f"singlecell_RETSAT_localization_{dataset}.csv",index=False)
    if dataset=="GSE200997":
        ann.groupby(["sample","patient","condition","cell_type"]).size().rename("n_cells").reset_index().to_csv(outdir/"GSE200997_epithelial_marker_cell_counts.csv",index=False)
    epi_label="epithelial_like" if dataset=="GSE200997" else "Epithelial cells"
    epi=gm[gm["cell_type"].eq(epi_label)].copy()
    expr_data={"RETSAT":epi["RETSAT_log2CPM"].to_numpy()}
    for g in gidx:
        if g!="RETSAT": expr_data[g]=epi[f"log2CPM::{g}"].to_numpy()
    expr=pd.DataFrame(expr_data,index=epi["sample"].astype(str))
    expr.index.name="sample"
    meta=epi[["sample","patient","condition","n_cells"]].copy().set_index("sample")
    meta.index=meta.index.astype(str)
    return gm,loc,expr,meta

def parse_series_matrix_metadata(path: Path):
    fields={}
    with gzip.open(path,"rt",encoding="utf-8",errors="replace",newline="") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                header=next(csv.reader([next(f).rstrip("\r\n")],delimiter="\t")); break
            if line.startswith("!Sample_"):
                row=next(csv.reader([line.rstrip("\r\n")],delimiter="\t"))
                fields[row[0]]=row[1:]
    gsm=[s.strip().strip('"') for s in fields.get("!Sample_geo_accession",[])]
    titles=fields.get("!Sample_title",[])
    chars=[v for k,v in fields.items() if k.startswith("!Sample_characteristics_ch1")]
    samples=[]
    for i,sid in enumerate(gsm):
        title=titles[i] if i<len(titles) else ""
        text=" | ".join([title]+[x[i] for x in chars if i<len(x)])
        patient=condition=None
        if path.name.startswith("GSE44076"):
            m=re.search(r"(?:paired sample from|sample from)\s+([A-Z]\d{4})\s+patient",title,re.I)
            if m:
                patient=m.group(1).upper()
                if re.search(r"normal paired",title,re.I): condition="Normal"
                elif re.search(r"tumou?r|cancer",title,re.I): condition="Tumor"
        else:
            m=re.search(r"patient\s+([A-Za-z]*\d+[A-Za-z]?)",title,re.I)
            if m: patient=m.group(1).lower()
            if re.search(r"\bnormal\b",title,re.I): condition="Normal"
            elif re.search(r"tumou?r|cancer",title,re.I): condition="Tumor"
        samples.append({"sample":sid,"title":title,"patient":patient,"condition":condition,"metadata_text":text})
    return samples,header

def series_platform(path: Path):
    with gzip.open(path,"rt",encoding="utf-8",errors="replace") as f:
        for line in f:
            if line.startswith("!Series_platform_id"): return line.split("\t",1)[1].strip().strip('"')
            if line.startswith("!series_matrix_table_begin"): break
    raise ValueError(f"Platform missing: {path}")

def load_probe_mapping(platform: str):
    mapping={}
    if platform=="GPL570":
        p=GEO_DIR/"GPL570.annot.gz"
        with gzip.open(p,"rt",encoding="utf-8",errors="replace") as f:
            header=None
            for line in f:
                if line.startswith("ID\t"):
                    header=next(csv.reader([line.rstrip("\r\n")],delimiter="\t"))
                    break
            if header is None: raise ValueError("GPL570 symbol header missing")
            reader=csv.DictReader(f,fieldnames=header,delimiter="\t")
            for row in reader:
                probe=(row.get("ID") or "").strip(); symbol=(row.get("Gene symbol") or "").strip()
                if probe and symbol and symbol not in {"---","-"}:
                    ss=[x.strip().upper() for x in re.split(r"\s*///\s*|\s*\|\s*|\s*;\s*",symbol) if x.strip()]
                    if len(set(ss))==1: mapping[probe]=ss[0]
    elif platform=="GPL13667":
        p=GEO_DIR/"GPL13667_probe_to_symbol.tsv"
        with p.open("rt",encoding="utf-8",errors="replace") as f:
            for line in f:
                parts=line.rstrip("\r\n").split("\t")
                if len(parts)>=2 and parts[1].strip():
                    ss=[x.strip().upper() for x in re.split(r"\s*\|\s*|\s*///\s*",parts[1]) if x.strip()]
                    if len(set(ss))==1: mapping[parts[0].strip()]=ss[0]
    else: raise ValueError(f"Unsupported platform: {platform}")
    return mapping

def parse_bulk_pairs(acc: str,samples: list[dict]):
    df=pd.DataFrame(samples)
    df=df[df["patient"].notna() & df["condition"].isin(["Tumor","Normal"])].copy()
    keep=[]
    for pid,g in df.groupby("patient"):
        t=g[g["condition"].eq("Tumor")]; n=g[g["condition"].eq("Normal")]
        if len(t)==1 and len(n)==1: keep.extend([t.iloc[0].to_dict(),n.iloc[0].to_dict()])
    out=pd.DataFrame(keep).drop_duplicates("sample")
    if out.empty: raise ValueError(f"No unique patient pairs in {acc}")
    return out.sort_values(["patient","condition"]).reset_index(drop=True)

def load_bulk_cohort(acc: str,genesets: dict[str,set[str]],outdir: Path):
    path=GEO_DIR/f"{acc}_series_matrix.txt.gz"
    samples,header=parse_series_matrix_metadata(path)
    platform=series_platform(path); probe_to_gene=load_probe_mapping(platform)
    targets={"RETSAT"}
    for gs in genesets.values(): targets.update(gs)
    probe_to_gene={p:g for p,g in probe_to_gene.items() if g in targets}
    by_gene={}
    for p,g in probe_to_gene.items(): by_gene.setdefault(g,[]).append(p)
    sample_ids=[x.strip().strip('"') for x in header[1:]]
    values={g:[] for g in by_gene}
    with gzip.open(path,"rt",encoding="utf-8",errors="replace",newline="") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                next(f);break
        for line in f:
            if line.startswith("!series_matrix_table_end"): break
            probe,payload=line.rstrip("\r\n").split("\t",1)
            probe=probe.strip().strip('"')
            gene=probe_to_gene.get(probe)
            if gene:
                vals=np.fromstring(payload,sep="\t",dtype=float)
                if len(vals)!=len(sample_ids): raise ValueError(f"{acc}: sample count mismatch for {probe}")
                values[gene].append(vals)
    columns={g:np.median(np.vstack(a),axis=0) for g,a in values.items() if a}
    expr=pd.DataFrame(columns,index=sample_ids);expr.index.name="sample"
    meta=parse_bulk_pairs(acc,samples)
    meta.to_csv(outdir/f"{acc}_paired_sample_metadata.csv",index=False)
    return expr,meta,platform,{g:len(a) for g,a in values.items() if a}

def pair_table(expr: pd.DataFrame,meta: pd.DataFrame,genesets: dict[str,set[str]],cohort: str,outdir: Path):
    m=meta.copy();m["sample"]=m["sample"].astype(str);m["patient"]=m["patient"].astype(str);m["condition"]=m["condition"].astype(str)
    m=m[m["sample"].isin(expr.index)].copy()
    if "n_cells" in m.columns: m=m[m["n_cells"]>=10]
    pair_ids=[]
    for pid,g in m.groupby("patient"):
        if (g["condition"]=="Tumor").sum()==1 and (g["condition"]=="Normal").sum()==1: pair_ids.append(pid)
    m=m[m["patient"].isin(pair_ids)].copy()
    if len(pair_ids)<2: raise ValueError(f"{cohort}: <2 pairs after filtering")
    sm=m.set_index("sample")
    vals=expr.loc[m["sample"].values].copy();vals.index=m["sample"].values
    used={}
    for name,genes in genesets.items():
        gs=sorted((set(genes)-{"RETSAT"})&set(vals.columns))
        if len(gs)<5: raise ValueError(f"{cohort}: only {len(gs)} members measured for {name}")
        sub=vals[gs].astype(float);sd=sub.std(axis=0,ddof=0);keep=sd[sd>1e-12].index.tolist()
        vals[f"{name}_score"]=((sub[keep]-sub[keep].mean())/sd[keep]).mean(axis=1)
        used[name]=keep
    vals["RETSAT_log2expr"]=vals["RETSAT"] if "RETSAT" in vals.columns else np.nan
    extra=["n_cells"] if "n_cells" in sm.columns else []
    vals=vals.join(sm[["patient","condition"]+extra])
    measure_cols=["RETSAT_log2expr"]+[f"{x}_score" for x in genesets]
    deltas=[]
    for pid in sorted(pair_ids):
        g=vals[vals["patient"].eq(pid)];t=g[g["condition"].eq("Tumor")].iloc[0];n=g[g["condition"].eq("Normal")].iloc[0]
        d={"cohort":cohort,"patient":pid}
        for col in measure_cols: d[f"delta_{col}"]=float(t[col]-n[col]) if pd.notna(t[col]) and pd.notna(n[col]) else np.nan
        deltas.append(d)
    delta_df=pd.DataFrame(deltas);summ=[]
    for col in measure_cols:
        tv=[];nv=[]
        for pid in sorted(pair_ids):
            g=vals[vals["patient"].eq(pid)]
            tv.append(float(g.loc[g.condition.eq("Tumor"),col].iloc[0]));nv.append(float(g.loc[g.condition.eq("Normal"),col].iloc[0]))
        diff=np.asarray(tv)-np.asarray(nv)
        summ.append({"cohort":cohort,"measure":col,"n_pairs":len(pair_ids),"tumor_median":np.nanmedian(tv),"normal_median":np.nanmedian(nv),"median_paired_delta":np.nanmedian(diff),"mean_paired_delta":np.nanmean(diff),"fraction_delta_positive":np.nanmean(diff>0),"p_value":safe_wilcoxon(tv,nv)})
    summary=pd.DataFrame(summ);summary["BH_q_within_cohort_3measures"]=bh_adjust(summary["p_value"])
    corr=[]
    for name in genesets:
        x=delta_df["delta_RETSAT_log2expr"].to_numpy();y=delta_df[f"delta_{name}_score"].to_numpy()
        rho,p=spearmanr(x,y,nan_policy="omit") if len(x)>=4 and np.nanstd(x)>0 and np.nanstd(y)>0 else (np.nan,np.nan)
        corr.append({"cohort":cohort,"association":f"paired_delta_RETSAT_vs_{name}_score","n_pairs":len(delta_df),"spearman_rho":rho,"p_value":p})
    correlations=pd.DataFrame(corr);correlations["BH_q_within_cohort_2correlations"]=bh_adjust(correlations["p_value"])
    tests=[]
    for name,genes in genesets.items():
        measured=sorted((set(genes)-{"RETSAT"})&set(vals.columns))
        for gene in measured:
            tv=[];nv=[]
            for pid in sorted(pair_ids):
                g=vals[vals["patient"].eq(pid)];tv.append(float(g.loc[g.condition.eq("Tumor"),gene].iloc[0]));nv.append(float(g.loc[g.condition.eq("Normal"),gene].iloc[0]))
            diff=np.asarray(tv)-np.asarray(nv)
            tests.append({"cohort":cohort,"pathway":name,"gene":gene,"n_pairs":len(pair_ids),"median_paired_delta":np.median(diff),"fraction_delta_positive":np.mean(diff>0),"p_value":safe_wilcoxon(tv,nv)})
    gt=pd.DataFrame(tests)
    if len(gt): gt["BH_q_within_pathway_cohort"]=gt.groupby(["cohort","pathway"])["p_value"].transform(lambda x:bh_adjust(x.to_numpy()))
    vals.assign(cohort=cohort).to_csv(outdir/f"{cohort}_paired_expression_and_scores.csv")
    delta_df.to_csv(outdir/f"{cohort}_patient_paired_deltas.csv",index=False)
    return summary,correlations,gt,vals.reset_index(names="sample"),used

def make_figures(localization: pd.DataFrame,paired: pd.DataFrame,summary: pd.DataFrame,outdir: Path):
    primary=localization[localization.dataset.eq("GSE132465")].sort_values("mean_group_detection_pct")
    if len(primary):
        fig,ax=plt.subplots(figsize=(8.8,4.8))
        colors=["#4575b4" if x=="Epithelial cells" else "#bdbdbd" for x in primary.cell_type]
        ax.barh(primary.cell_type,primary.mean_group_detection_pct,color=colors)
        ax.set_xlabel("RETSAT-detected cells (%)");ax.set_ylabel("Author-annotated cell type (tumor samples)")
        ax.set_title("RETSAT localization across tumor cell types — GSE132465")
        for i,r in enumerate(primary.itertuples()): ax.text(r.mean_group_detection_pct+0.2,i,f"{int(r.positive_cells)}/{int(r.tumor_cells)}",va="center",fontsize=8)
        fig.tight_layout();fig.savefig(outdir/"RETSAT_GSE132465_celltype_localization.png",dpi=220);plt.close(fig)
    measures=["RETSAT_log2expr","Retinol_metabolism_score","Ferroptosis_score"]
    labels=["RETSAT log2 CPM/expression","Retinol pathway score","Ferroptosis pathway score"]
    cohorts=list(paired.cohort.drop_duplicates())
    fig,axes=plt.subplots(len(measures),len(cohorts),figsize=(4.1*len(cohorts),3.3*len(measures)),squeeze=False)
    for ri,(measure,label) in enumerate(zip(measures,labels)):
        for ci,cohort in enumerate(cohorts):
            ax=axes[ri,ci];d=paired[paired.cohort.eq(cohort)]
            p=d.pivot(index="patient",columns="condition",values=measure).dropna()
            if {"Normal","Tumor"}.issubset(p.columns):
                for _,r in p.iterrows(): ax.plot([0,1],[r["Normal"],r["Tumor"]],color="#888888",alpha=.28,lw=.65,zorder=1)
                ax.scatter(np.zeros(len(p)),p["Normal"],s=11,color="#4575b4",alpha=.65,zorder=2)
                ax.scatter(np.ones(len(p)),p["Tumor"],s=11,color="#d73027",alpha=.65,zorder=2)
                ax.set_xticks([0,1],["Normal","Tumor"])
                st=summary[(summary.cohort==cohort)&(summary.measure==measure)]
                pv=st.p_value.iloc[0] if len(st) else np.nan
                ax.set_title(f"{cohort} (n={len(p)}; p={pv:.3g})")
            else: ax.text(.5,.5,"paired data unavailable",ha="center",va="center")
            ax.set_ylabel(label);ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("Patient-level paired RETSAT and independent pathway signals",y=1.01,fontsize=13)
    fig.tight_layout();fig.savefig(outdir/"RETSAT_paired_cohort_signals.png",dpi=220,bbox_inches="tight");plt.close(fig)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    retinol=kegg_gene_set(RETINOL_ID,OUT);ferro=kegg_gene_set(FERROPTOSIS_ID,OUT)
    genesets={"Retinol_metabolism":retinol-{"RETSAT"},"Ferroptosis":ferro}
    pd.DataFrame([{"pathway":"Retinol_metabolism","gene":g,"excluded_target_gene":g=="RETSAT"} for g in sorted(retinol)]+[{"pathway":"Ferroptosis","gene":g,"excluded_target_gene":False} for g in sorted(ferro)]).to_csv(OUT/"KEGG_independent_pathway_gene_sets.tsv",sep="\t",index=False)
    print("KEGG sizes",{"Retinol_metabolism_excluding_RETSAT":len(genesets["Retinol_metabolism"]),"Ferroptosis":len(ferro)},"RETSAT originally in retinol pathway",("RETSAT" in retinol))
    summaries=[];correlations=[];gene_tests=[];all_pairs=[];localizations=[];pair_counts=[];coverage=[]
    for dataset in ["GSE132465","GSE200997"]:
        groups,loc,expr,meta=scan_single_cell(dataset,genesets,OUT);localizations.append(loc)
        meta["sample"]=meta.index.astype(str);meta["patient"]=meta["patient"].astype(str);meta["condition"]=meta["condition"].astype(str)
        meta=meta[meta.n_cells>=10].copy()
        pair_ids=[pid for pid,g in meta.groupby("patient") if (g.condition=="Tumor").sum()==1 and (g.condition=="Normal").sum()==1]
        pair_counts.append({"cohort":dataset,"n_pairs":len(pair_ids),"type":"single-cell epithelial pseudobulk","filter":"author epithelial label in GSE132465; ≥2/5 epithelial markers in GSE200997; ≥10 cells per sample"})
        for name,gs in genesets.items():
            measured=sorted(set(gs)&set(expr.columns))
            coverage.append({"cohort":dataset,"pathway":name,"gene_set_size":len(gs),"measured_genes":len(measured),"members":";".join(measured)})
        s,c,g,v,used=pair_table(expr,meta,genesets,dataset,OUT)
        summaries.append(s);correlations.append(c);gene_tests.append(g);all_pairs.append(v.assign(cohort=dataset))
    for acc in ["GSE44076","GSE32323"]:
        expr,meta,platform,probe_counts=load_bulk_cohort(acc,genesets,OUT)
        coverage.append({"cohort":acc,"pathway":"measured_gene_count","gene_set_size":len(set().union(*genesets.values())),"measured_genes":len(expr.columns),"members":";".join(sorted(expr.columns)),"platform":platform})
        pair_counts.append({"cohort":acc,"n_pairs":meta.patient.nunique(),"type":f"paired bulk microarray ({platform})","filter":"exactly one tumor and matched normal per patient; unrelated healthy donor normals excluded"})
        s,c,g,v,used=pair_table(expr,meta,genesets,acc,OUT)
        summaries.append(s);correlations.append(c);gene_tests.append(g);all_pairs.append(v.assign(cohort=acc))
    summary=pd.concat(summaries,ignore_index=True);corr=pd.concat(correlations,ignore_index=True)
    gt=pd.concat([x for x in gene_tests if len(x)],ignore_index=True);pairs=pd.concat(all_pairs,ignore_index=True);local=pd.concat(localizations,ignore_index=True)
    summary.to_csv(OUT/"cohort_paired_summary.csv",index=False);corr.to_csv(OUT/"paired_delta_correlations.csv",index=False)
    gt.to_csv(OUT/"pathway_member_paired_tests.csv",index=False);pairs.to_csv(OUT/"all_cohorts_paired_sample_expression_scores.csv",index=False)
    local.to_csv(OUT/"singlecell_RETSAT_localization_summary.csv",index=False)
    pd.DataFrame(pair_counts).to_csv(OUT/"cohort_pair_counts.csv",index=False);pd.DataFrame(coverage).to_csv(OUT/"pathway_gene_coverage.csv",index=False)
    make_figures(local,pairs,summary,OUT)
    # Report file
    lines=["# RETSAT：CRC恶性上皮定位、患者配对与独立通路线索","",
    "## 方法","",
    "主定位使用GSE132465作者提供的细胞类型注释；肿瘤样本中作者注释为上皮的细胞作为恶性上皮候选群。该数据没有在本次分析中新建拷贝数恶性判定。GSE200997没有通用细胞类型注释，以EPCAM、KRT8、KRT18、KRT19、KRT20中至少2个检出定义上皮样细胞，作为敏感性分析；不把该标记规则当作CNV恶性证据。","",
    "单细胞表达按患者/样本和上皮类型聚合原始UMI为伪bulk，再做患者配对比较；不将单个细胞当作独立生物学重复。每个样本至少10个上皮细胞。统计为双侧Wilcoxon signed-rank。单细胞RETSAT使用伪bulk CPM；bulk数据使用GEO系列矩阵处理后表达值。","",
    "视黄醇代谢与铁死亡分别使用KEGG人类hsa00830和hsa04216；视黄醇集合剔除RETSAT本身。通路评分是每个队列内基因标准化后的均值Z分数。RETSAT变化与通路分数变化的相关性按患者配对差值计算。通路成员单基因配对检验在各队列、各通路内BH校正。关联属于探索性表达支持，不代表因果调控。","",
    "## 配对样本数",""]
    for r in pair_counts: lines.append(f"- {r['cohort']}: {r['n_pairs']}对；{r['type']}；{r['filter']}")
    lines += ["","## RETSAT在肿瘤细胞中的定位",""]
    for r in local.to_dict("records"):
        lines.append(f"- {r['dataset']} / {r['cell_type']}: 肿瘤样本细胞 {int(r['tumor_cells'])} 个、患者 {int(r['patients'])} 位；RETSAT阳性 {int(r['positive_cells'])} 个（样本组检出率均值 {r['mean_group_detection_pct']:.2f}%），占全部肿瘤RETSAT阳性细胞 {r.get('pct_of_tumor_RETSAT_positive_cells',np.nan):.1f}%。")
    lines += ["","## 患者配对肿瘤–正常结果","",
    "| Cohort | 指标 | 配对数 | 中位配对差（肿瘤−正常） | 差值为正比例 | Wilcoxon P | 3项指标内BH q |","|---|---|---:|---:|---:|---:|---:|"]
    for r in summary.to_dict("records"):
        lines.append(f"| {r['cohort']} | {r['measure']} | {int(r['n_pairs'])} | {r['median_paired_delta']:.4g} | {r['fraction_delta_positive']:.1%} | {r['p_value']:.3g} | {r['BH_q_within_cohort_3measures']:.3g} |")
    lines += ["","## RETSAT差值与通路差值的患者级关联","",
    "| Cohort | 关联 | 配对数 | Spearman ρ | P | 两项关联内BH q |","|---|---|---:|---:|---:|---:|"]
    for r in corr.to_dict("records"):
        lines.append(f"| {r['cohort']} | {r['association']} | {int(r['n_pairs'])} | {r['spearman_rho']:.3f} | {r['p_value']:.3g} | {r['BH_q_within_cohort_2correlations']:.3g} |")
    sig=gt[gt.BH_q_within_pathway_cohort<0.05]
    lines += ["","## 通路成员线索",""]
    if len(sig):
        for (cohort,pathway),sub in sig.groupby(["cohort","pathway"]):
            members=", ".join(f"{r.gene} ({r.median_paired_delta:+.2f})" for r in sub.sort_values("BH_q_within_pathway_cohort").itertuples())
            lines.append(f"- {cohort} / {pathway}: BH-FDR<0.05成员 {len(sub)} 个：{members}")
    else: lines.append("- 所有通路成员均未达到队列内、通路内BH-FDR<0.05；完整表保留名义P值及方向。")
    lines += ["","## 解读边界","",
    "跨队列配对表达和通路评分只能支持RETSAT的定位与相关性，不能证明它驱动视黄醇代谢或铁死亡。单细胞通路有效基因覆盖数见pathway_gene_coverage.csv。GSE200997上皮样定义由标记基因建立；两个单细胞数据集均未在本分析中新建CNV判定。GSE21510未纳入主bulk比较，因为元数据标示其肿瘤样本为LCM而配对正常为homogenized，样本处理区室不匹配。","",
    "## 主要文件","",
    "- RETSAT_GSE132465_celltype_localization.png — GSE132465肿瘤细胞类型检出率。",
    "- RETSAT_paired_cohort_signals.png — 4个队列患者级配对表达与通路评分。",
    "- cohort_paired_summary.csv、paired_delta_correlations.csv、pathway_member_paired_tests.csv — 统计汇总。",
    "- singlecell_pseudobulk_all_groups_*.csv与*_paired_expression_and_scores.csv — 样本级伪bulk与表达数据。","",
    "## 输入与来源","",
    "仅读取本项目work目录下GSE132465/GSE200997原始UMI矩阵及原始注释、GSE44076/GSE32323 GEO系列矩阵和对应探针注释；没有读取旧分析结果。KEGG路径定义来自官方人类通路hsa00830（视黄醇代谢）及hsa04216（铁死亡）。",""]
    (OUT/"RETSAT_CRC_epithelial_paired_pathways_report.md").write_text("\n".join(lines),encoding="utf-8")
    print("\nPAIR SUMMARY\n",summary[["cohort","measure","n_pairs","median_paired_delta","p_value","BH_q_within_cohort_3measures"]].to_string(index=False))
    print("\nCORRELATIONS\n",corr.to_string(index=False))
    print("\nLOCALIZATION\n",local.to_string(index=False))
    print("\nOUTPUT",OUT)

if __name__=="__main__": main()
