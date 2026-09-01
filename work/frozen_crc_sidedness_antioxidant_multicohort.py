"""Frozen multi-cohort validation of CRC sidedness, PUFA incorporation and antioxidant buffering.

Frozen cohorts: GSE39582, GSE41258, GSE4554, GSE75316, TCGA-COAD.

Primary readouts:
- PUFA incorporation = mean z(ACSL4, LPCAT3, AGPAT3), require >=2/3 genes.
- Antioxidant buffering = mean z(SLC7A11, GPX4, AIFM2, GCH1, DHODH), require >=4/5 genes.
- SLC7A11 z-expression.
- GCH1 z-expression.

Primary sidedness rule:
Right = cecum/caecum/ascending/hepatic flexure/transverse/proximal/right-sided colon.
Left = descending/splenic flexure/sigmoid/rectosigmoid/distal/left-sided colon.
Rectum-only and generic/ambiguous colon are excluded.

Automatic sidedness uses SAMPLE-LEVEL title/source/characteristics only. Every decision is exported
for manual audit. Optional work/data/frozen_crc_sidedness_manual_overrides.csv can overwrite side
labels with columns: cohort,sample_id,side,reason. Allowed side values: right,left,exclude_rectum,
unknown,ambiguous.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm, ttest_ind

try:
    import GEOparse  # type: ignore
except ImportError as exc:
    raise SystemExit("Install GEOparse first: pip install GEOparse") from exc

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
DATA = WORK / "data"
OUT = ROOT / "outputs"
GEO_CACHE = DATA / "geo_frozen_sidedness"
OVERRIDE_FILE = DATA / "frozen_crc_sidedness_manual_overrides.csv"

FROZEN_GEO = ["GSE39582", "GSE41258", "GSE4554", "GSE75316"]
FROZEN_COHORTS = FROZEN_GEO + ["TCGA-COAD"]
PUFA_GENES = ["ACSL4", "LPCAT3", "AGPAT3"]
BUFFER_GENES = ["SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH"]
TARGET_GENES = sorted(set(PUFA_GENES + BUFFER_GENES))

TCGA_CASES = DATA / "tcga_coad_cases.json"
TCGA_EXPR_CANDIDATES = [
    DATA / "tcga_ferroptosis_three_layer_expression_uqfpkm.tsv",
    DATA / "tcga_aa_defense_candidate_expression_uqfpkm.tsv",
    DATA / "tcga_coad_gene_expression_uqfpkm.tsv",
]
TCGA_GENE_IDS = {
    "ACSL4": "ENSG00000068366", "LPCAT3": "ENSG00000111684", "AGPAT3": "ENSG00000160216",
    "SLC7A11": "ENSG00000151012", "GPX4": "ENSG00000167468", "AIFM2": "ENSG00000042286",
    "GCH1": "ENSG00000131979", "DHODH": "ENSG00000102967",
}

RIGHT = [r"\bcecum\b", r"\bcaecum\b", r"\bcecal\b", r"\bcaecal\b", r"\bascending\b",
         r"hepatic\s+flexure", r"\btransverse\b", r"\bproximal\b", r"right\s+(?:sided\s+)?colon",
         r"right[- ]sided"]
LEFT = [r"\bdescending\b", r"splenic\s+flexure", r"\bsigmoid\b", r"\brectosigmoid\b",
        r"\bdistal\b", r"left\s+(?:sided\s+)?colon", r"left[- ]sided"]
RECTUM = [r"\brectum\b", r"\brectal\b"]
NORMAL = [r"\bnormal\b", r"adjacent\s+normal", r"non[- ]tumou?r", r"healthy", r"control\s+mucosa"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=OUT)
    p.add_argument("--geo-cache", type=Path, default=GEO_CACHE)
    return p.parse_args()


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def bh_adjust(values: list[float]) -> list[float]:
    p = np.asarray(values, float); out = np.full(len(p), np.nan); ok = np.isfinite(p)
    if not ok.any(): return out.tolist()
    idx = np.where(ok)[0]; order = np.argsort(p[ok]); ranked = p[ok][order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked)+1); adj = np.minimum.accumulate(adj[::-1])[::-1]
    restored = np.empty_like(adj); restored[order] = np.minimum(adj, 1.0); out[idx] = restored
    return out.tolist()


def list_text(value: object) -> str:
    if isinstance(value, (list, tuple)): return " | ".join(map(str, value))
    return str(value)


def sample_level_text(metadata: dict[str, object]) -> str:
    keys = [k for k in metadata if k in {"title", "source_name_ch1"} or k.startswith("characteristics_ch1")]
    return " | ".join(f"{k}: {list_text(metadata[k])}" for k in keys)


def all_metadata_text(metadata: dict[str, object]) -> str:
    return " | ".join(f"{k}: {list_text(v)}" for k, v in metadata.items())


def classify_side(text: str) -> tuple[str, str]:
    low = text.lower(); rh = [p for p in RIGHT if re.search(p, low)]; lh = [p for p in LEFT if re.search(p, low)]
    rect = [p for p in RECTUM if re.search(p, low)]
    if rh and not lh: return "right", rh[0]
    if lh and not rh:
        if rect and not re.search(r"rectosigmoid", low): return "exclude_rectum", rect[0]
        return "left", lh[0]
    if rect and not rh and not lh: return "exclude_rectum", rect[0]
    if rh and lh: return "ambiguous", f"right={rh[0]};left={lh[0]}"
    return "unknown", ""


def is_tumor_sample(text: str) -> bool:
    return not any(re.search(p, text.lower()) for p in NORMAL)


def gene_symbol_column(gpl: pd.DataFrame) -> str:
    preferred = ["Gene Symbol", "GENE_SYMBOL", "Gene symbol", "Symbol", "SYMBOL", "gene_assignment"]
    for c in preferred:
        if c in gpl.columns: return c
    for c in gpl.columns:
        if "gene" in str(c).lower() and "symbol" in str(c).lower(): return c
    raise ValueError(f"No gene-symbol column found: {list(gpl.columns)[:30]}")


def symbol_tokens(value: object) -> set[str]:
    text = str(value).upper()
    return {t for t in re.split(r"\s*///\s*|\s*//\s*|\s*;\s*|\s*,\s*|\s+", text)
            if re.fullmatch(r"[A-Z0-9._-]+", t or "")}


def collapse_to_targets(expr_probe: pd.DataFrame, gpl: pd.DataFrame) -> pd.DataFrame:
    probes = gpl["ID"].astype(str) if "ID" in gpl.columns else pd.Series(gpl.index.astype(str), index=gpl.index)
    sym = gene_symbol_column(gpl); mapping = {g: [] for g in TARGET_GENES}
    for probe, raw in zip(probes, gpl[sym], strict=False):
        toks = symbol_tokens(raw)
        for gene in TARGET_GENES:
            if gene in toks: mapping[gene].append(str(probe))
    out = pd.DataFrame(index=expr_probe.columns)
    for gene, pids in mapping.items():
        avail = [p for p in pids if p in expr_probe.index]
        out[gene] = expr_probe.loc[avail].apply(pd.to_numeric, errors="coerce").median(axis=0) if avail else np.nan
    return out


def load_overrides() -> pd.DataFrame:
    if not OVERRIDE_FILE.exists(): return pd.DataFrame(columns=["cohort", "sample_id", "side", "reason"])
    x = pd.read_csv(OVERRIDE_FILE, dtype=str).fillna("")
    need = {"cohort", "sample_id", "side"}; missing = need - set(x.columns)
    if missing: raise ValueError(f"Override file missing columns: {sorted(missing)}")
    allowed = {"right", "left", "exclude_rectum", "unknown", "ambiguous"}
    bad = sorted(set(x["side"]) - allowed)
    if bad: raise ValueError(f"Invalid override side values: {bad}")
    if "reason" not in x: x["reason"] = ""
    return x


def apply_overrides(meta: pd.DataFrame, cohort: str, overrides: pd.DataFrame) -> pd.DataFrame:
    x = meta.copy(); ov = overrides[overrides["cohort"].eq(cohort)]
    for _, row in ov.iterrows():
        sid = row["sample_id"]
        if sid in x.index:
            x.loc[sid, "side_auto"] = x.loc[sid, "side"]
            x.loc[sid, "side"] = row["side"]
            x.loc[sid, "side_match"] = "manual_override"
            x.loc[sid, "override_reason"] = row.get("reason", "")
    return x


def local_series_matrix(accession: str) -> Path | None:
    candidates = [
        DATA / f"{accession}_series_matrix.clean.txt.gz",
        DATA / f"{accession}_series_matrix.txt.gz",
        GEO_CACHE / f"{accession}_series_matrix.txt.gz",
    ]
    return next((p for p in candidates if p.exists()), None)


def parse_series_matrix_metadata(path: Path) -> list[dict[str, str]]:
    samples: list[dict[str, str]] = []
    with gzip.open(path, "rt", newline="", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            fields = next(csv.reader([line.rstrip("\n\r")], delimiter="\t"))
            key = fields[0][len("!Sample_"):]
            while len(samples) < len(fields) - 1:
                samples.append({})
            for i, value in enumerate(fields[1:]):
                out_key = key
                suffix = 2
                while out_key in samples[i]:
                    out_key = f"{key}_{suffix}"
                    suffix += 1
                samples[i][out_key] = value
    return samples


def load_geo_matrix(accession: str, path: Path, overrides: pd.DataFrame):
    sample_rows = parse_series_matrix_metadata(path)
    sample_ids = [row.get("geo_accession", "") for row in sample_rows]
    if not sample_rows or not all(sample_ids):
        raise ValueError(f"{accession}: series matrix has no complete sample accession metadata")
    expr_probe = pd.read_csv(path, sep="\t", compression="gzip", comment="!", index_col=0, low_memory=False)
    expr_probe.columns = expr_probe.columns.astype(str)
    platform_ids = sorted({row.get("platform_id", "") for row in sample_rows if row.get("platform_id", "")})
    pieces = []
    for gpl_id in platform_ids:
        annotation_path = DATA / f"{gpl_id}.annot.gz"
        if not annotation_path.exists():
            raise FileNotFoundError(f"{accession}: missing local platform annotation {annotation_path}")
        gpl = GEOparse.get_GEO(filepath=str(annotation_path), geotype="GPL", silent=True)
        ids = [sid for sid, row in zip(sample_ids, sample_rows, strict=False) if row.get("platform_id", "") == gpl_id]
        ids = [sid for sid in ids if sid in expr_probe.columns]
        if ids:
            pieces.append(collapse_to_targets(expr_probe.loc[:, ids], gpl.table))
    if not pieces:
        raise ValueError(f"{accession}: no usable expression matrix/platform annotation")
    expr = pd.concat(pieces, axis=0)
    expr = expr[~expr.index.duplicated(keep="first")]
    rows = []
    for sid, row in zip(sample_ids, sample_rows, strict=False):
        core = sample_level_text(row)
        side, match = classify_side(core)
        rows.append({"sample_id": sid, "cohort": accession, "side": side, "side_auto": side, "side_match": match,
                     "override_reason": "", "is_tumor_like": is_tumor_sample(core), "sample_level_text": core,
                     "raw_metadata_text": all_metadata_text(row), "platform_id": row.get("platform_id", "")})
    meta = pd.DataFrame(rows).set_index("sample_id")
    meta = apply_overrides(meta, accession, overrides)
    common = expr.index.intersection(meta.index)
    expr = expr.loc[common]
    meta = meta.loc[common]
    return expr, meta, {"source": str(path), "platforms": platform_ids, "n": len(expr),
                        "side_counts": meta.side.value_counts().to_dict()}


def load_geo(accession: str, cache: Path, overrides: pd.DataFrame):
    matrix_path = local_series_matrix(accession)
    if matrix_path is not None:
        return load_geo_matrix(accession, matrix_path, overrides)
    cache.mkdir(parents=True, exist_ok=True)
    gse = GEOparse.get_GEO(geo=accession, destdir=str(cache), silent=True)
    pieces = []
    platform_ids = sorted({gsm.metadata.get("platform_id", [""])[0] for gsm in gse.gsms.values()})
    for gpl_id in platform_ids:
        if not gpl_id or gpl_id not in gse.gpls: continue
        ids = [n for n,gsm in gse.gsms.items() if gsm.metadata.get("platform_id", [""])[0] == gpl_id]
        matrix = {}
        for n in ids:
            tab = gse.gsms[n].table
            if {"ID_REF", "VALUE"}.issubset(tab.columns): matrix[n] = tab.set_index("ID_REF")["VALUE"]
        if matrix: pieces.append(collapse_to_targets(pd.DataFrame(matrix), gse.gpls[gpl_id].table))
    if not pieces: raise ValueError(f"{accession}: no usable expression matrix")
    expr = pd.concat(pieces, axis=0); expr = expr[~expr.index.duplicated(keep="first")]
    rows = []
    for sid, gsm in gse.gsms.items():
        core = sample_level_text(gsm.metadata); side, match = classify_side(core)
        rows.append({"sample_id":sid, "cohort":accession, "side":side, "side_auto":side, "side_match":match,
                     "override_reason":"", "is_tumor_like":is_tumor_sample(core), "sample_level_text":core,
                     "raw_metadata_text":all_metadata_text(gsm.metadata), "platform_id":list_text(gsm.metadata.get("platform_id", []))})
    meta = pd.DataFrame(rows).set_index("sample_id"); meta = apply_overrides(meta, accession, overrides)
    common = expr.index.intersection(meta.index); expr = expr.loc[common]; meta = meta.loc[common]
    return expr, meta, {"platforms":platform_ids, "n":len(expr), "side_counts":meta.side.value_counts().to_dict()}


def parse_tcga_expression(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", index_col=0)
    idx = raw.index.astype(str); cols = raw.columns.astype(str)
    # Case 1: symbol-labeled matrix.
    if sum(g in set(idx) for g in TARGET_GENES) >= 3:
        return raw.loc[[g for g in TARGET_GENES if g in raw.index]].T.apply(pd.to_numeric, errors="coerce")
    if sum(g in set(cols) for g in TARGET_GENES) >= 3:
        return raw.loc[:, [g for g in TARGET_GENES if g in raw.columns]].apply(pd.to_numeric, errors="coerce")
    # Case 2: Ensembl IDs x cases (the project GDC targeted caches use this form).
    id_to_gene = {v:k for k,v in TCGA_GENE_IDS.items()}
    stripped = pd.Index([x.split(".")[0] for x in idx])
    raw = raw.copy(); raw["__gene__"] = stripped.map(id_to_gene)
    sub = raw[raw["__gene__"].notna()].copy().set_index("__gene__")
    if sub.empty: raise ValueError(f"Could not map target genes in {path}")
    sub = sub[~sub.index.duplicated(keep="first")]
    return sub.T.apply(pd.to_numeric, errors="coerce")


def load_tcga(overrides: pd.DataFrame):
    if not TCGA_CASES.exists(): raise FileNotFoundError(TCGA_CASES)
    path = next((p for p in TCGA_EXPR_CANDIDATES if p.exists()), None)
    if path is None: raise FileNotFoundError("No TCGA expression cache found")
    expr = parse_tcga_expression(path); expr.index = expr.index.astype(str)
    cases = json.loads(TCGA_CASES.read_text(encoding="utf-8"))
    if isinstance(cases, dict):
        cases = cases.get("data", cases)
        if isinstance(cases, dict):
            cases = cases.get("hits", cases.get("data", []))
    rows = []
    for case in cases:
        sid = str(case.get("case_id") or case.get("submitter_id") or "")
        text = json.dumps(case, ensure_ascii=False); side, match = classify_side(text)
        rows.append({"sample_id":sid,"cohort":"TCGA-COAD","side":side,"side_auto":side,"side_match":match,
                     "override_reason":"","is_tumor_like":True,"sample_level_text":text,"raw_metadata_text":text,"platform_id":"RNA-seq"})
    meta = pd.DataFrame(rows).set_index("sample_id"); meta = apply_overrides(meta, "TCGA-COAD", overrides)
    common = expr.index.intersection(meta.index); expr = expr.loc[common]; meta = meta.loc[common]
    return expr, meta, {"expression_path":str(path),"n":len(expr),"side_counts":meta.side.value_counts().to_dict(),
                        "available_genes":sorted(set(expr.columns).intersection(TARGET_GENES))}


def add_scores(expr: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    x = expr.apply(pd.to_numeric, errors="coerce").copy(); finite = x.to_numpy(float); finite = finite[np.isfinite(finite)]
    if finite.size and np.nanpercentile(finite,95) > 50: x = np.log2(np.clip(x,0,None)+1)
    z = pd.DataFrame(index=x.index)
    for gene in TARGET_GENES: z[gene] = zscore(x[gene]) if gene in x else np.nan
    out = meta.copy()
    for gene in TARGET_GENES: out[gene] = z[gene]
    pufa = [g for g in PUFA_GENES if z[g].notna().any()]; buf = [g for g in BUFFER_GENES if z[g].notna().any()]
    out["pufa_available_genes"] = z[pufa].notna().sum(axis=1) if pufa else 0
    out["buffer_available_genes"] = z[buf].notna().sum(axis=1) if buf else 0
    out["pufa_incorporation_score"] = z[pufa].mean(axis=1,skipna=True) if pufa else np.nan
    out["antioxidant_buffering_score"] = z[buf].mean(axis=1,skipna=True) if buf else np.nan
    out.loc[out.pufa_available_genes < 2,"pufa_incorporation_score"] = np.nan
    out.loc[out.buffer_available_genes < 4,"antioxidant_buffering_score"] = np.nan
    return out


def compare_metric(data: pd.DataFrame, metric: str) -> dict[str, object]:
    d = data[data.is_tumor_like & data.side.isin(["right","left"])]
    r = pd.to_numeric(d.loc[d.side.eq("right"),metric],errors="coerce").dropna(); l = pd.to_numeric(d.loc[d.side.eq("left"),metric],errors="coerce").dropna()
    base = {"metric":metric,"n_right":len(r),"n_left":len(l)}
    if len(r)<2 or len(l)<2: return {**base,"right_mean":np.nan,"left_mean":np.nan,"right_minus_left":np.nan,"hedges_g":np.nan,"se_g":np.nan,"welch_p":np.nan,"mannwhitney_p":np.nan}
    vr,vl=r.var(ddof=1),l.var(ddof=1); n=len(r)+len(l); pooled=math.sqrt(((len(r)-1)*vr+(len(l)-1)*vl)/(n-2))
    d0=(r.mean()-l.mean())/pooled if pooled>0 else np.nan; j=1-3/(4*n-9); g=d0*j if np.isfinite(d0) else np.nan
    se=math.sqrt(n/(len(r)*len(l))+(g*g)/(2*(n-2))) if np.isfinite(g) else np.nan
    return {**base,"right_mean":r.mean(),"left_mean":l.mean(),"right_minus_left":r.mean()-l.mean(),"hedges_g":g,"se_g":se,
            "welch_p":ttest_ind(r,l,equal_var=False).pvalue,"mannwhitney_p":mannwhitneyu(r,l,alternative="two-sided").pvalue}


def random_effects_meta(stats: pd.DataFrame, metric: str) -> dict[str, object]:
    d=stats[(stats.metric.eq(metric))&stats.hedges_g.notna()&stats.se_g.notna()]
    if len(d)<2: return {"metric":metric,"k":len(d),"estimable":False}
    y=d.hedges_g.to_numpy(float); v=np.square(d.se_g.to_numpy(float)); w=1/v; fixed=np.sum(w*y)/np.sum(w)
    q=np.sum(w*np.square(y-fixed)); df=len(y)-1; c=np.sum(w)-np.sum(w*w)/np.sum(w); tau=max(0,(q-df)/c) if c>0 else 0
    wr=1/(v+tau); pooled=np.sum(wr*y)/np.sum(wr); se=math.sqrt(1/np.sum(wr)); z=pooled/se
    return {"metric":metric,"k":len(d),"estimable":True,"pooled_hedges_g":pooled,"se":se,"ci_low":pooled-1.96*se,
            "ci_high":pooled+1.96*se,"p":2*norm.sf(abs(z)),"tau2":tau,"i2_percent":max(0,(q-df)/q)*100 if q>0 else 0,
             "right_higher_cohorts":int(np.sum(y>0)),"right_lower_cohorts":int(np.sum(y<0))}


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(c) for c in frame.columns]
    rows = []
    for values in frame.itertuples(index=False, name=None):
        row = []
        for value in values:
            row.append("" if pd.isna(value) else str(value).replace("|", "\\|"))
        rows.append(row)
    lines = ["| " + " | ".join(columns) + " |",
             "| " + " | ".join("---" for _ in columns) + " |"]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def main() -> None:
    args=parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True); args.geo_cache.mkdir(parents=True,exist_ok=True)
    overrides=load_overrides(); scores=[]; audits=[]; provenance={}
    for cohort in FROZEN_GEO:
        expr,meta,info=load_geo(cohort,args.geo_cache,overrides); sc=add_scores(expr,meta); scores.append(sc); provenance[cohort]=info
        audits.append(sc[["cohort","side","side_auto","side_match","override_reason","is_tumor_like","sample_level_text","raw_metadata_text","platform_id"]])
    expr,meta,info=load_tcga(overrides); sc=add_scores(expr,meta); scores.append(sc); provenance["TCGA-COAD"]=info
    audits.append(sc[["cohort","side","side_auto","side_match","override_reason","is_tumor_like","sample_level_text","raw_metadata_text","platform_id"]])
    combined=pd.concat(scores,axis=0,sort=False); audit=pd.concat(audits,axis=0,sort=False)
    metrics=["pufa_incorporation_score","antioxidant_buffering_score","SLC7A11","GCH1"]
    rows=[]
    for cohort in FROZEN_COHORTS:
        for metric in metrics:
            row=compare_metric(combined[combined.cohort.eq(cohort)],metric); row["cohort"]=cohort; rows.append(row)
    stats=pd.DataFrame(rows)
    for metric in metrics:
        m=stats.metric.eq(metric); stats.loc[m,"welch_fdr_across_cohorts"] = bh_adjust(stats.loc[m,"welch_p"].tolist())
    meta=pd.DataFrame([random_effects_meta(stats,m) for m in metrics])
    prefix=args.out_dir/"frozen_crc_sidedness_antioxidant_multicohort"
    combined.to_csv(str(prefix)+"_scores.csv"); audit.to_csv(str(prefix)+"_side_audit.csv"); stats.to_csv(str(prefix)+"_stats.csv",index=False); meta.to_csv(str(prefix)+"_meta.csv",index=False)
    manifest={"frozen_cohorts":FROZEN_COHORTS,"metrics":metrics,"pufa_genes":PUFA_GENES,"buffer_genes":BUFFER_GENES,
              "minimum_genes":{"pufa":"2/3","buffer":"4/5"},"override_file":str(OVERRIDE_FILE),"provenance":provenance}
    Path(str(prefix)+"_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2,default=str),encoding="utf-8")
    lines=["# Frozen CRC sidedness antioxidant-buffering validation","","Cohorts: "+", ".join(FROZEN_COHORTS),"",
           "Primary metrics are transcriptional states, not direct biochemical measurements.","","## Cohort results","",markdown_table(stats),"",
           "## Random-effects meta-analysis","",markdown_table(meta),"","## Audit","",
           "Review `_side_audit.csv`; ambiguous/unknown samples are excluded unless manually overridden in work/data/frozen_crc_sidedness_manual_overrides.csv."]
    Path(str(prefix)+"_report.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print(json.dumps({"stats":stats.to_dict("records"),"meta":meta.to_dict("records"),"override_file":str(OVERRIDE_FILE)},ensure_ascii=False,indent=2,default=str))


if __name__ == "__main__":
    main()
