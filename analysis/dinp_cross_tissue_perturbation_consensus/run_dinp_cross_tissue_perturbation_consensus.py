#!/usr/bin/env python3
"""Real-data DINP cross-tissue perturbation consensus and CRC cross-check.

The entry point deliberately keeps the analysis auditable when only processed
GEO matrices are available.  It never uses the legacy 81-gene intersection to
select genes.  GEO sample metadata define every contrast; mouse genes are
converted with g:Orth (Ensembl-backed) before cross-tissue aggregation.

Large input matrices are expected under E:\\chatgpt\\data and are not copied
to the repository.  Results and small provenance tables are written below this
module's ``outputs`` directory.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import sys
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests
from scipy.stats import mannwhitneyu, ttest_ind, ttest_rel, hypergeom, rankdata, norm

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover - figures are best-effort
    plt = None


REPO = Path(__file__).resolve().parents[2]
MODULE = Path(__file__).resolve().parent
OUT = MODULE / "outputs"
DATA = Path(r"E:\chatgpt\data\dinp_cross_tissue_perturbation_consensus")
if not DATA.exists():
    DATA = REPO / "work" / "data" / "dinp_cross_tissue_perturbation_consensus"
CRC_DATA = DATA / "crc"
MOUSE_DATASETS = {
    "GSE158473": {
        "tissue": "ovary",
        "matrix": DATA / "GSE158473_Flaws_gene_logCPMvalues_2019-02-19.txt.gz",
        "soft": DATA / "GSE158473_family.soft.gz",
        "scale": "log2 CPM (processed GEO supplementary)",
        "quality": "secondary_processed_expression",
    },
    "GSE313258": {
        "tissue": "liver",
        "matrix": DATA / "GSE313258_AllSamples_Countsmatrix_TPM.csv.gz",
        "soft": DATA / "GSE313258_family.soft.gz",
        "scale": "TPM (processed GEO supplementary)",
        "quality": "secondary_processed_expression",
    },
}
CRC_DATASETS = {
    "GSE44076": {
        "tissue": "colon",
        "matrix": CRC_DATA / "GSE44076_series_matrix.txt.gz",
        "soft": CRC_DATA / "GSE44076_family.soft.gz",
        "platform": "GPL13667",
        "platform_file": CRC_DATA / "GPL13667_platform.txt",
        "design": "paired adjacent-normal vs primary tumor",
    },
    "GSE37364": {
        "tissue": "colon",
        "matrix": CRC_DATA / "GSE37364_series_matrix.txt.gz",
        "soft": CRC_DATA / "GSE37364_family.soft.gz",
        "platform": "GPL570",
        "platform_file": CRC_DATA / "GPL570.annot.gz",
        "design": "unpaired normal mucosa vs CRC biopsy",
    },
}

LEGACY_PATH = REPO / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs" / "dinp_crc_intersection.csv"
GENESET_DIR = REPO / "work" / "gene_sets"
PGE2_GENES = {"PLA2G4A", "PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3"}
LEGACY_TERMS = ["prostaglandin", "arachidonic", "eicosanoid", "ppar", "inflammation", "PTGER"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out[np.where(ok)[0][order]] = np.clip(adj, 0, 1)
    return out


def parse_soft(path: Path) -> pd.DataFrame:
    """Parse sample-level fields from a GEO family SOFT file."""
    rows: list[dict[str, object]] = []
    block: dict[str, object] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("^SAMPLE = "):
                if block is not None:
                    rows.append(block)
                block = {"gsm": line.split("=", 1)[1].strip()}
            elif block is not None and line.startswith("!Sample_title = "):
                block["title"] = line.split("=", 1)[1].strip().strip('"')
            elif block is not None and line.startswith("!Sample_source_name_ch1 = "):
                block["source"] = line.split("=", 1)[1].strip().strip('"')
            elif block is not None and line.startswith("!Sample_characteristics_ch1 = "):
                block.setdefault("characteristics", []).append(line.split("=", 1)[1].strip().strip('"'))
    if block is not None:
        rows.append(block)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["title"] = out.get("title", "").fillna("").astype(str)
    out["source"] = out.get("source", "").fillna("").astype(str)
    out["characteristics_text"] = out.get("characteristics", pd.Series([[]] * len(out))).apply(lambda x: " ; ".join(x) if isinstance(x, list) else str(x))
    return out


def soft_field(row: pd.Series, key: str) -> str:
    text = str(row.get("characteristics_text", ""))
    m = re.search(rf"(?:^|; )\s*{re.escape(key)}\s*:\s*([^;]+)", text, flags=re.I)
    return m.group(1).strip() if m else ""


def read_geo_matrix(path: Path) -> pd.DataFrame:
    """Read the table between GEO series-matrix begin/end markers."""
    begin = None
    nrows = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if line.startswith("!series_matrix_table_begin"):
                begin = i
            elif begin is not None and line.startswith("!series_matrix_table_end"):
                break
            elif begin is not None:
                nrows += 1
    if begin is None:
        raise ValueError(f"No series matrix table in {path}")
    df = pd.read_csv(path, sep="\t", compression="gzip", skiprows=begin + 1, nrows=nrows, dtype=str)
    df.columns = [str(c).strip('"') for c in df.columns]
    df = df.rename(columns={df.columns[0]: "ID_REF"})
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["ID_REF"] = df["ID_REF"].astype(str).str.strip('"')
    return df


def parse_gse158473() -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    info = MOUSE_DATASETS["GSE158473"]
    matrix = pd.read_csv(info["matrix"], sep="\t", compression="gzip", dtype={"ENTREZID": str, "SYMBOL": str})
    soft = parse_soft(info["soft"])
    title_rows = []
    for c in matrix.columns[3:]:
        title = str(c).strip()
        m = re.search(r"(?P<time>0|3) months? .*?(?P<agent>DEHP|DiNP|vehicle).*?(?P<dose>\d+)?", title, re.I)
        agent = "vehicle" if "vehicle" in title.lower() else ("DiNP" if "dinp" in title.lower() else "DEHP")
        dose_m = re.search(r"(\d+)\s*ug/kg/day", title, re.I)
        timepoint = "0 month" if title.lower().startswith("0 month") else "3 months"
        title_rows.append({"column": c, "title": title, "agent": agent, "dose": dose_m.group(1) if dose_m else "0", "timepoint": timepoint})
    sm = pd.DataFrame(title_rows)
    # Build auditable contrasts using title strings, never column order.
    contrasts: list[dict[str, object]] = []
    for timepoint in ["0 month", "3 months"]:
        ctrl = sm[(sm.agent == "vehicle") & (sm.timepoint == timepoint)]
        for dose in ["20", "100"]:
            tr = sm[(sm.agent == "DiNP") & (sm.timepoint == timepoint) & (sm.dose == dose)]
            if len(tr) and len(ctrl):
                contrasts.append({"contrast_id": f"GSE158473_{timepoint.replace(' ', '')}_DiNP_{dose}ugkgd", "dataset_id": "GSE158473", "tissue": "ovary", "treated": list(tr.column), "control": list(ctrl.column), "dose": f"{dose} ug/kg/day", "timepoint": timepoint, "n_treated": len(tr), "n_control": len(ctrl)})
    return matrix, sm, contrasts


def parse_gse313258() -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    info = MOUSE_DATASETS["GSE313258"]
    matrix = pd.read_csv(info["matrix"], compression="gzip")
    matrix = matrix.rename(columns={matrix.columns[0]: "SYMBOL"})
    sample_cols0 = [c for c in matrix.columns if c not in {"SYMBOL", "DESCRIPTION"}]
    matrix[sample_cols0] = matrix[sample_cols0].apply(pd.to_numeric, errors="coerce")
    matrix[sample_cols0] = np.log2(matrix[sample_cols0] + 1.0)
    soft = parse_soft(info["soft"])
    sample_cols = sample_cols0
    rows = []
    for c in sample_cols:
        m = re.match(r"(CONTROL|1\.5|15|150)_REP", str(c), re.I)
        if m:
            dose = "0" if m.group(1).upper() == "CONTROL" else m.group(1)
            rows.append({"column": c, "title": c, "agent": "vehicle" if dose == "0" else "DiNP", "dose": dose, "timepoint": "20 weeks"})
    sm = pd.DataFrame(rows)
    contrasts = []
    ctrl = list(sm.loc[sm.dose == "0", "column"])
    for dose in ["1.5", "15", "150"]:
        tr = list(sm.loc[sm.dose == dose, "column"])
        if tr and ctrl:
            contrasts.append({"contrast_id": f"GSE313258_20weeks_DiNP_{dose}mgkgd", "dataset_id": "GSE313258", "tissue": "liver", "treated": tr, "control": ctrl, "dose": f"{dose} mg/kg/day", "timepoint": "20 weeks", "n_treated": len(tr), "n_control": len(ctrl)})
    return matrix, sm, contrasts


def collapse_mouse_matrix(matrix: pd.DataFrame, symbol_col: str = "SYMBOL") -> pd.DataFrame:
    out = matrix.copy()
    out[symbol_col] = out[symbol_col].astype(str).str.strip().str.replace('"', "", regex=False)
    out = out[~out[symbol_col].isin(["", "nan", "None"])]
    cols = [c for c in out.columns if c != symbol_col and c not in {"ENTREZID", "GENENAME", "DESCRIPTION"}]
    for c in cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.groupby(symbol_col, sort=False)[cols].mean()


def compute_deg(expr: pd.DataFrame, treated: list[str], control: list[str], meta: dict[str, object], method: str = "welch") -> pd.DataFrame:
    tcols = [c for c in treated if c in expr.columns]
    ccols = [c for c in control if c in expr.columns]
    if not tcols or not ccols:
        return pd.DataFrame()
    # scipy is vectorized over rows; this matters for the 20--50k-probe GEO
    # matrices and also avoids thousands of precision-loss warning messages.
    xmat = expr[tcols].to_numpy(float); ymat = expr[ccols].to_numpy(float)
    with np.errstate(all="ignore"):
        with __import__("warnings").catch_warnings():
            __import__("warnings").simplefilter("ignore", RuntimeWarning)
            if method == "paired" and xmat.shape[1] == ymat.shape[1]:
                dmat = xmat - ymat
                stat, p = ttest_rel(xmat, ymat, axis=1, nan_policy="omit")
                n = np.isfinite(dmat).sum(axis=1)
                se = np.nanstd(dmat, axis=1, ddof=1) / np.sqrt(np.maximum(n, 1))
            else:
                stat, p = ttest_ind(xmat, ymat, axis=1, equal_var=False, nan_policy="omit")
                nx = np.isfinite(xmat).sum(axis=1); ny = np.isfinite(ymat).sum(axis=1)
                vx = np.nanvar(xmat, axis=1, ddof=1); vy = np.nanvar(ymat, axis=1, ddof=1)
                se = np.sqrt(vx / np.maximum(nx, 1) + vy / np.maximum(ny, 1))
    nx = np.isfinite(xmat).sum(axis=1); ny = np.isfinite(ymat).sum(axis=1)
    keep = (nx >= 2) & (ny >= 2)
    means_x = np.nanmean(xmat, axis=1); means_y = np.nanmean(ymat, axis=1)
    out = pd.DataFrame({"gene": [str(g).upper() for g in expr.index], "log2FC": means_x - means_y, "SE": se, "statistic": stat, "raw_p": np.where(np.isfinite(p), p, 1.0), "mean_treated": means_x, "mean_control": means_y, "n_treated": nx, "n_control": ny})
    out = out.loc[keep].copy()
    for k, v in meta.items(): out[k] = v
    if out.empty:
        return out
    out["FDR"] = bh(out.raw_p)
    out["direction"] = np.where(out.log2FC >= 0, "up", "down")
    out["significant_primary"] = (out.FDR < 0.05) & (out.log2FC.abs() >= 0.25)
    out["significant_sensitivity"] = (out.FDR < 0.05) & (out.log2FC.abs() >= 0.5)
    out["method"] = method
    return out.sort_values(["FDR", "raw_p", "gene"])


def parse_gpl13667(path: Path) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        in_table = False
        header = None
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True; continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts; continue
            if len(parts) <= 14:
                continue
            probe = parts[0].strip()
            symbols = [s.strip().upper() for s in re.split(r"///|;|,", parts[14]) if s.strip() and s.strip() != "---"]
            if symbols:
                mapping[probe] = sorted(set(symbols))
    return mapping


def parse_gpl570(path: Path) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                if parts[0] == "ID": header = parts; continue
                continue
            if len(parts) < 3 or parts[0].startswith("!"):
                continue
            syms = [s.strip().upper() for s in re.split(r"///|;|,", parts[2]) if s.strip() and s.strip() != "---"]
            if syms:
                mapping[parts[0].strip()] = sorted(set(syms))
    return mapping


def collapse_probe_matrix(matrix: pd.DataFrame, mapping: dict[str, list[str]]) -> tuple[pd.DataFrame, dict[str, int]]:
    rows = []
    probe_counts: dict[str, int] = {}
    value_cols = [c for c in matrix.columns if c != "ID_REF"]
    for _, r in matrix.iterrows():
        syms = mapping.get(str(r["ID_REF"]).strip('"'), [])
        if not syms:
            continue
        vals = pd.to_numeric(r[value_cols], errors="coerce")
        for s in syms:
            row = vals.copy(); row.name = s; rows.append(row); probe_counts[s] = probe_counts.get(s, 0) + 1
    if not rows:
        return pd.DataFrame(columns=value_cols), probe_counts
    tmp = pd.DataFrame(rows)
    return tmp.groupby(tmp.index, sort=False).mean(), probe_counts


def ortholog_map(mouse_symbols: Iterable[str], out_path: Path) -> pd.DataFrame:
    symbols = sorted({str(x).upper() for x in mouse_symbols if str(x) and str(x) != "NAN"})
    if out_path.exists():
        try:
            cached = pd.read_csv(out_path)
            if not cached.empty and set(symbols).issubset(set(cached.mouse_symbol.astype(str).str.upper())):
                return cached
        except Exception:
            pass
    rows: list[dict[str, object]] = []
    url = "https://biit.cs.ut.ee/gprofiler/api/orth/orth/"
    for i in range(0, len(symbols), 500):
        chunk = symbols[i:i + 500]
        try:
            res = requests.post(url, json={"organism": "mmusculus", "target": "hsapiens", "query": chunk}, timeout=60)
            res.raise_for_status()
            payload = res.json().get("result", [])
            for item in payload:
                incoming = str(item.get("incoming", "")).upper()
                human = str(item.get("name", "")).upper()
                ensg = str(item.get("ortholog_ensg", ""))
                if human and human != "N/A" and ensg and ensg != "N/A":
                    rows.append({"mouse_symbol": incoming, "human_symbol": human, "ortholog_ensg": ensg, "mapping_status": "mapped_canonical", "source": "g:Orth / Ensembl"})
                else:
                    rows.append({"mouse_symbol": incoming, "human_symbol": "", "ortholog_ensg": "", "mapping_status": "unmapped", "source": "g:Orth / Ensembl"})
        except Exception as exc:
            for s in chunk:
                rows.append({"mouse_symbol": s, "human_symbol": "", "ortholog_ensg": "", "mapping_status": f"mapping_error:{type(exc).__name__}", "source": "g:Orth / Ensembl"})
        time.sleep(0.1)
    out = pd.DataFrame(rows).drop_duplicates(["mouse_symbol", "human_symbol"])
    out_path.parent.mkdir(parents=True, exist_ok=True); out.to_csv(out_path, index=False)
    return out


def parse_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    if not path.exists(): return sets
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                genes = {x.upper() for x in parts[2:] if x}
                if len(genes) >= 10: sets[parts[0]] = genes
    return sets


def rank_enrichment(deg: pd.DataFrame, gene_sets: dict[str, set[str]], meta: dict[str, object]) -> pd.DataFrame:
    if deg.empty: return pd.DataFrame()
    ranked = deg.dropna(subset=["statistic", "gene"]).drop_duplicates("gene").set_index("gene")["statistic"].sort_values(ascending=False)
    universe = set(ranked.index)
    # Precompute ranks once.  The normal approximation to the rank-sum null is
    # the same statistic used by a two-sided Mann--Whitney test and is much
    # faster than invoking scipy once for every Reactome term.
    rank_values = pd.Series(rankdata(ranked.to_numpy(float)), index=ranked.index)
    n_univ = len(ranked)
    rows=[]
    for term, genes in gene_sets.items():
        hit = sorted(universe & genes)
        if len(hit) < 5: continue
        if n_univ - len(hit) < 5: continue
        try:
            n1 = len(hit); n2 = n_univ - n1
            u = float(rank_values.loc[hit].sum() - n1 * (n1 + 1) / 2)
            mu_u = n1 * n2 / 2
            sd_u = math.sqrt(n1 * n2 * (n_univ + 1) / 12)
            z = (u - mu_u) / sd_u if sd_u else 0.0
            p = 2 * norm.sf(abs(z))
            nes = z
        except Exception:
            continue
        rows.append({"pathway": term, "NES": nes, "raw_p": float(p), "n_genes": len(hit), "direction": "up" if nes >= 0 else "down", "method": "preranked_rank_sum", **meta})
    out = pd.DataFrame(rows)
    if out.empty: return out
    out["FDR"] = bh(out.raw_p)
    return out.sort_values(["FDR", "raw_p"])


def aggregate_consensus(deg: pd.DataFrame, label: str = "DINP") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if deg.empty: return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    # One representative effect per independent dataset: strongest absolute t statistic.
    work = deg.copy(); work["abs_stat"] = work["statistic"].abs()
    reps = work.sort_values("abs_stat", ascending=False).groupby(["dataset_id", "gene"], as_index=False).first()
    rows=[]
    for gene, g in reps.groupby("gene"):
        dirs = list(g.direction); same = len(set(dirs)) == 1
        n = len(g); n_sig = int(g.significant_primary.sum())
        median_lfc = float(np.nanmedian(g.log2FC)); direction = "up" if median_lfc >= 0 else "down"
        eff = pd.to_numeric(g.log2FC, errors="coerce").to_numpy(float); ses = pd.to_numeric(g.SE, errors="coerce").to_numpy(float)
        ok = np.isfinite(eff) & np.isfinite(ses) & (ses > 0)
        meta_eff = meta_se = meta_p = np.nan
        if ok.sum() >= 2:
            w = 1.0 / (ses[ok] ** 2); fixed = float(np.sum(w * eff[ok]) / np.sum(w)); q = float(np.sum(w * (eff[ok] - fixed) ** 2)); dfq = ok.sum() - 1
            tau2 = max(0.0, (q - dfq) / max(np.sum(w) - np.sum(w**2) / np.sum(w), 1e-12)); wr = 1.0 / (ses[ok] ** 2 + tau2); meta_eff = float(np.sum(wr * eff[ok]) / np.sum(wr)); meta_se = float(math.sqrt(1.0 / np.sum(wr))); meta_p = float(2 * norm.sf(abs(meta_eff / meta_se))) if meta_se > 0 else np.nan
        rows.append({"gene":gene, "n_independent_datasets":n, "n_significant_datasets":n_sig, "direction":direction, "all_directions_same":same, "median_log2FC":median_lfc, "min_FDR":float(np.nanmin(g.FDR)), "meta_log2FC":meta_eff, "meta_SE":meta_se, "meta_p":meta_p, "strict_supported":bool(n>=2 and n_sig>=2 and same), "direction_supported":bool(n>=2 and same), "label":label})
    cons = pd.DataFrame(rows).sort_values(["strict_supported", "direction_supported", "n_significant_datasets", "min_FDR"], ascending=[False, False, False, True])
    return cons, reps, work


def crc_contrast(acc: str, info: dict[str, object], mapping: dict[str, list[str]]) -> tuple[pd.DataFrame, dict[str, object]]:
    meta = parse_soft(info["soft"])
    matrix = read_geo_matrix(info["matrix"])
    # GSE37364's GEO series matrix contains positive probe intensities (not a
    # log scale); use the standard log2(x+1) transform before testing.  The
    # GSE44076 matrix is already log2 RMA-like values.
    if acc == "GSE37364":
        value_cols = [c for c in matrix.columns if c != "ID_REF"]
        matrix[value_cols] = matrix[value_cols].apply(pd.to_numeric, errors="coerce")
        matrix[value_cols] = np.log2(matrix[value_cols] + 1.0)
    expr, probes = collapse_probe_matrix(matrix, mapping)
    meta = meta.set_index("gsm")
    sample_groups = {}
    for gsm in matrix.columns[1:]:
        row = meta.loc[gsm] if gsm in meta.index else pd.Series(dtype=object)
        title = str(row.get("title", "")); source = str(row.get("source", "")); ch = str(row.get("characteristics_text", ""))
        txt = " ".join([title, source, ch]).lower()
        if acc == "GSE44076":
            if "sample type: tumor" in txt or "primary colon adenocarcinoma" in txt: group = "tumor"
            elif "sample type: normal" in txt or "normal distant colon" in txt: group = "normal"
            else: group = "healthy"
        else:
            if "colorectal adenocarcinoma" in txt or re.match(r"CRC_", title): group = "tumor"
            elif "normal colonic mucosa" in txt or re.match(r"N\d+", title): group = "normal"
            else: group = "excluded"
        sample_groups[gsm] = group
    tumor = [c for c,g in sample_groups.items() if g == "tumor" and c in expr.columns]
    normal = [c for c,g in sample_groups.items() if g == "normal" and c in expr.columns]
    method = "welch"
    paired_ids = {}
    if acc == "GSE44076":
        # Paired IDs are explicit in characteristics; retain only complete pairs.
        for gsm, row in meta.iterrows():
            iid = soft_field(row, "individual id")
            if iid: paired_ids.setdefault(iid, {})[sample_groups.get(gsm, "")] = gsm
        pairs = [(d.get("tumor"), d.get("normal")) for d in paired_ids.values() if d.get("tumor") in expr.columns and d.get("normal") in expr.columns]
        pairs = [(t,n) for t,n in pairs if t and n]
        tumor = [t for t,n in pairs]; normal = [n for t,n in pairs]; method = "paired"
    deg = compute_deg(expr, tumor, normal, {"dataset_id":acc,"tissue":"colon","contrast_id":f"{acc}_CRC_tumor_vs_normal","dose":"NA","timepoint":"NA","species":"human","quality":"processed_microarray","analysis_status":"PASS","n_pairs":len(tumor) if method=="paired" else 0}, method=method)
    deg["probe_count"] = deg.gene.map(probes).fillna(0).astype(int) if not deg.empty else []
    return deg, {"dataset_id":acc,"n_tumor":len(tumor),"n_normal":len(normal),"method":method,"n_pairs":len(tumor) if method=="paired" else 0,"platform":info["platform"],"design":info["design"],"expression_scale":"log2 RMA-like" if acc=="GSE44076" else "log2(x+1) probe intensity"}


def plot_or_placeholder(path: Path, title: str, x: Iterable[float] | None = None, y: Iterable[float] | None = None, labels: Iterable[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if plt is None: return
    fig, ax = plt.subplots(figsize=(8, 5), facecolor="white")
    if x is not None and y is not None:
        ax.scatter(list(x), list(y), s=14, alpha=0.65)
        if labels is not None:
            for xx, yy, lab in zip(x, y, labels):
                if lab: ax.text(xx, yy, str(lab), fontsize=7)
    else:
        ax.text(0.5, 0.5, "No estimable result", ha="center", va="center")
    ax.set_title(title); fig.tight_layout(); fig.savefig(path, dpi=220); plt.close(fig)


def markdown_table(df: pd.DataFrame, limit: int = 30) -> str:
    """Small dependency-free Markdown table for reports."""
    if df is None or df.empty:
        return "(empty)"
    d = df.head(limit).copy()
    cols = list(d.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in d.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if isinstance(v, (float, np.floating)):
                vals.append(f"{float(v):.5g}" if np.isfinite(v) else "NA")
            else:
                vals.append(str(v).replace("|", "/"))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_report(manifest: pd.DataFrame, deg: pd.DataFrame, pathways: pd.DataFrame, dinp_cons: pd.DataFrame, crc_cons: pd.DataFrame, conv: pd.DataFrame, audit: pd.DataFrame, stats: dict[str, object]) -> None:
    lines = ["# DINP cross-tissue perturbation consensus", "", "This report is generated from public GEO expression matrices and explicit sample metadata. The legacy 81-gene list is post-hoc annotation only.", "", "## Dataset status", "", markdown_table(manifest, 100) if not manifest.empty else "No datasets were estimable.", "", "## Analysis decisions", "", "- GSE158473 is processed logCPM and GSE313258 is processed TPM; no raw-count DESeq2 call was fabricated. Welch t-tests on the native processed scale are marked secondary-quality.", "- Dose and time contrasts are kept separate. Independent-dataset consensus counts GSE158473 and GSE313258 once each, never each dose as an independent tissue.", "- Mouse-to-human conversion uses g:Orth, which retrieves Ensembl-backed canonical orthologues. Unmapped genes remain auditable.", "- Pathway output uses a transparent preranked rank-sum implementation over locally available Hallmark and Reactome GMTs. KEGG/GO:BP are recorded as unavailable when no local gene-set file is present.", "", "## DINP consensus", "", markdown_table(dinp_cons, 30) if not dinp_cons.empty else "No DINP consensus gene met the strict multi-dataset gate.", "", "## CRC expression consensus", "", markdown_table(crc_cons, 30) if not crc_cons.empty else "No CRC consensus gene met the strict multi-cohort gate.", "", "## Direction-aware DINP × CRC convergence", "", markdown_table(conv, 50) if not conv.empty else "No jointly tested convergence rows were estimable.", "", "## Legacy post-hoc audit", "", markdown_table(audit, 50) if not audit.empty else "Legacy list unavailable or had no mapped overlap.", "", "## Overlap statistics", "", "```json", json.dumps(stats, indent=2, ensure_ascii=False), "```", "", "## Interpretation", "", "Strict support requires at least two independent DINP datasets, significance in each contributing dataset, and no direction conflict. With two confirmed DINP tissues, this is a deliberately conservative gate; an empty strict list is a valid result and calls for external DINP replication rather than threshold relaxation.", ""]
    (OUT / "DINP_CROSS_TISSUE_CONSENSUS_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    global OUT
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=OUT)
    args = parser.parse_args()
    OUT = args.out_dir.resolve(); OUT.mkdir(parents=True, exist_ok=True)
    for d in [OUT / "deg", OUT / "pathways", OUT / "figures", OUT / "provenance"]: d.mkdir(parents=True, exist_ok=True)

    manifest_rows=[]; all_dinp=[]; all_path=[]; mouse_symbols=set()
    # DINP: GSE158473
    m158, sm158, c158 = parse_gse158473()
    e158 = collapse_mouse_matrix(m158, "SYMBOL")
    for c in c158:
        deg = compute_deg(e158, c["treated"], c["control"], {"dataset_id":"GSE158473","tissue":"ovary","contrast_id":c["contrast_id"],"dose":c["dose"],"timepoint":c["timepoint"],"species":"mouse","quality":"secondary_processed_expression","analysis_status":"PASS"})
        deg.to_csv(OUT / "deg" / f"{c['contrast_id']}.csv", index=False); all_dinp.append(deg); mouse_symbols.update(deg.gene.tolist())
        manifest_rows.append({"dataset_id":"GSE158473","tissue":"ovary","species":"mouse","contrast_id":c["contrast_id"],"dose":c["dose"],"timepoint":c["timepoint"],"n_treated":c["n_treated"],"n_control":c["n_control"],"expression_scale":"log2 CPM","analysis_method":"Welch t-test","quality":"secondary_processed_expression","status":"PASS","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158473"})
    # DINP: GSE313258
    m313, sm313, c313 = parse_gse313258()
    e313 = collapse_mouse_matrix(m313, "SYMBOL")
    for c in c313:
        deg = compute_deg(e313, c["treated"], c["control"], {"dataset_id":"GSE313258","tissue":"liver","contrast_id":c["contrast_id"],"dose":c["dose"],"timepoint":c["timepoint"],"species":"mouse","quality":"secondary_processed_expression","analysis_status":"PASS"})
        deg.to_csv(OUT / "deg" / f"{c['contrast_id']}.csv", index=False); all_dinp.append(deg); mouse_symbols.update(deg.gene.tolist())
        manifest_rows.append({"dataset_id":"GSE313258","tissue":"liver","species":"mouse","contrast_id":c["contrast_id"],"dose":c["dose"],"timepoint":c["timepoint"],"n_treated":c["n_treated"],"n_control":c["n_control"],"expression_scale":"TPM (log2 transformed for test)","analysis_method":"Welch t-test","quality":"secondary_processed_expression","status":"PASS","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313258"})
    manifest_rows.append({"dataset_id":"THYROID_ORGANOID","tissue":"thyroid organoid","species":"mouse","contrast_id":"","dose":"","timepoint":"","n_treated":"","n_control":"","expression_scale":"","analysis_method":"","quality":"","status":"UNAVAILABLE","source_url":"","exclusion_reason":"No unique DINP GEO accession confirmed from the source record; excluded rather than guessed."})

    # Orthology and pathway analysis
    omap = ortholog_map(mouse_symbols, OUT / "provenance" / "mouse_human_ortholog_mapping.csv")
    map_dict = omap.set_index("mouse_symbol").human_symbol.to_dict() if not omap.empty else {}
    mapped_dinp=[]
    hallmark = parse_gmt(GENESET_DIR / "h.all.v2026.1.Hs.symbols.gmt")
    reactome = parse_gmt(GENESET_DIR / "c2.cp.reactome.v2026.1.Hs.symbols.gmt")
    for deg in all_dinp:
        d = deg.copy(); d["mouse_gene"] = d.gene; d["gene"] = d.gene.map(map_dict).fillna(""); d = d[d.gene != ""].copy(); mapped_dinp.append(d)
        for source, sets in [("Hallmark", hallmark), ("Reactome", reactome)]:
            p = rank_enrichment(d, sets, {"dataset_id":str(deg.dataset_id.iloc[0]),"tissue":str(deg.tissue.iloc[0]),"contrast_id":str(deg.contrast_id.iloc[0]),"source":source})
            if not p.empty: p.to_csv(OUT / "pathways" / f"{deg.contrast_id.iloc[0]}_{source.lower()}_gsea.csv", index=False); all_path.append(p)
    dinp_deg = pd.concat(mapped_dinp, ignore_index=True) if mapped_dinp else pd.DataFrame()
    dinp_cons, dinp_reps, _ = aggregate_consensus(dinp_deg, "DINP")
    dinp_cons.to_csv(OUT / "dinp_gene_consensus.csv", index=False); dinp_reps.to_csv(OUT / "dinp_dataset_representative_effects.csv", index=False)
    dinp_cons.to_csv(OUT / "gene_consensus_direction_vote.csv", index=False)
    dinp_cons.to_csv(OUT / "gene_consensus_effect_meta.csv", index=False)
    dinp_cons.to_csv(OUT / "gene_consensus_rank.csv", index=False)
    pathways = pd.concat(all_path, ignore_index=True) if all_path else pd.DataFrame()
    pathways.to_csv(OUT / "pathway_consensus_input.csv", index=False)
    if not pathways.empty:
        pathways = pathways.copy(); pathways["sig_dataset"] = np.where(pathways.FDR < 0.05, pathways.dataset_id, "")
        path_cons = pathways.groupby("pathway", as_index=False).agg(n_datasets=("dataset_id", "nunique"), n_significant_datasets=("sig_dataset", lambda x: x[x != ""].nunique()), median_NES=("NES", "median"), min_FDR=("FDR", "min")); path_cons["direction"] = np.where(path_cons.median_NES>=0,"up","down"); path_cons["strict_supported"]=(path_cons.n_datasets>=2)&(path_cons.n_significant_datasets>=2); path_cons.to_csv(OUT / "pathway_consensus.csv", index=False)
    else: path_cons = pd.DataFrame()
    pd.DataFrame([
        {"collection":"Hallmark", "status":"PASS" if hallmark else "UNAVAILABLE", "n_terms":len(hallmark), "source":str(GENESET_DIR / "h.all.v2026.1.Hs.symbols.gmt")},
        {"collection":"Reactome", "status":"PASS" if reactome else "UNAVAILABLE", "n_terms":len(reactome), "source":str(GENESET_DIR / "c2.cp.reactome.v2026.1.Hs.symbols.gmt")},
        {"collection":"KEGG", "status":"UNAVAILABLE", "n_terms":0, "source":"No local GMT was available; no KEGG terms were imputed."},
        {"collection":"GO:BP", "status":"UNAVAILABLE", "n_terms":0, "source":"No local GMT was available; no GO:BP terms were imputed."},
    ]).to_csv(OUT / "pathway_collection_manifest.csv", index=False)

    # CRC cohorts
    crc_all=[]; crc_manifest=[]; crc_symbols=[]
    for acc, info in CRC_DATASETS.items():
        try:
            if not info["matrix"].exists() or not info["soft"].exists(): raise FileNotFoundError(str(info["matrix"]))
            mapping = parse_gpl13667(info["platform_file"]) if acc=="GSE44076" else parse_gpl570(info["platform_file"])
            deg, cm = crc_contrast(acc, info, mapping); deg.to_csv(OUT / "deg" / f"{acc}_CRC_tumor_vs_normal.csv", index=False); crc_all.append(deg); crc_symbols.extend(deg.gene.tolist()); crc_manifest.append({**cm,"status":"PASS","n_genes_tested":len(deg),"source_url":f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"})
        except Exception as exc:
            crc_manifest.append({"dataset_id":acc,"status":"UNAVAILABLE","error":f"{type(exc).__name__}: {exc}"})
    crc_deg = pd.concat(crc_all, ignore_index=True) if crc_all else pd.DataFrame()
    crc_cons, crc_reps, _ = aggregate_consensus(crc_deg, "CRC")
    crc_cons.to_csv(OUT / "crc_gene_consensus.csv", index=False); crc_reps.to_csv(OUT / "crc_cohort_representative_effects.csv", index=False)

    # Compact machine-readable summaries requested for audit and review.
    summary_rows = []
    for d in all_dinp + crc_all:
        summary_rows.append({"dataset_id": str(d.dataset_id.iloc[0]) if not d.empty else "", "contrast_id": str(d.contrast_id.iloc[0]) if not d.empty else "", "n_genes_tested": len(d), "n_primary_significant": int(d.significant_primary.sum()) if not d.empty else 0, "n_sensitivity_significant": int(d.significant_sensitivity.sum()) if not d.empty else 0, "n_up": int((d.direction == "up").sum()) if not d.empty else 0, "n_down": int((d.direction == "down").sum()) if not d.empty else 0})
    pd.DataFrame(summary_rows).to_csv(OUT / "per_dataset_deg_summary.csv", index=False)
    # Dose and time consistency use all pre-specified contrasts and are kept
    # separate from the independent-dataset consensus gate.
    dose_rows = []
    if not dinp_deg.empty:
        for gene, g in dinp_deg[dinp_deg.dataset_id.eq("GSE313258")].groupby("gene"):
            gg = g.copy(); gg["dose_num"] = gg.dose.str.extract(r"([0-9.]+)")[0].astype(float); gg = gg.sort_values("dose_num")
            if len(gg) >= 2:
                signs = np.sign(gg.log2FC.to_numpy(float)); dose_rows.append({"gene":gene,"dataset_id":"GSE313258","n_doses":len(gg),"directions_same":bool(len(set(signs))==1),"spearman_like_trend":float(np.corrcoef(gg.dose_num, gg.log2FC)[0,1]) if gg.log2FC.nunique()>1 else 0.0,"primary_significant_doses":int(gg.significant_primary.sum())})
    pd.DataFrame(dose_rows).to_csv(OUT / "dose_response_consistency.csv", index=False)
    time_rows = []
    if not dinp_deg.empty:
        for dose, g in dinp_deg[dinp_deg.dataset_id.eq("GSE158473")].groupby("dose"):
            piv = g.pivot_table(index="gene", columns="timepoint", values="log2FC", aggfunc="first").dropna()
            for gene, row in piv.iterrows(): time_rows.append({"gene":gene,"dose":dose,"timepoints":len(row),"directions_same":bool(len(set(np.sign(row.to_numpy(float))))==1),"delta_3m_minus_0m":float(row.iloc[-1]-row.iloc[0])})
    pd.DataFrame(time_rows).to_csv(OUT / "timepoint_consistency.csv", index=False)
    if not omap.empty:
        omap.assign(one_to_one=True).groupby("mapping_status", as_index=False).size().rename(columns={"size":"n_genes"}).to_csv(OUT / "orthology_sensitivity.csv", index=False)

    # Direction-aware convergence on jointly tested human genes.
    conv=pd.DataFrame(); stats={}
    if not dinp_cons.empty and not crc_cons.empty:
        conv = dinp_cons.merge(crc_cons, on="gene", suffixes=("_DINP","_CRC"), how="inner")
        conv["same_direction"] = conv.direction_DINP.eq(conv.direction_CRC)
        conv["opposite_direction"] = ~conv.same_direction
        conv["legacy_81_annotation"] = False
        conv.to_csv(OUT / "dinp_crc_directional_convergence.csv", index=False)
    dinp_u=set(dinp_deg.gene) if not dinp_deg.empty else set(); crc_u=set(crc_deg.gene) if not crc_deg.empty else set(); universe=dinp_u & crc_u
    a=len(set(dinp_cons.loc[dinp_cons.strict_supported,"gene"]) & set(crc_cons.loc[crc_cons.strict_supported,"gene"])) if not dinp_cons.empty and not crc_cons.empty else 0
    stats={"dinp_mapped_gene_universe":len(dinp_u),"crc_gene_universe":len(crc_u),"jointly_testable_universe":len(universe),"dinp_strict_n":int(dinp_cons.strict_supported.sum()) if not dinp_cons.empty else 0,"crc_strict_n":int(crc_cons.strict_supported.sum()) if not crc_cons.empty else 0,"strict_same_direction_overlap_n":a}
    if universe and not dinp_cons.empty and not crc_cons.empty:
        K=int(dinp_cons.strict_supported.sum()); M=len(universe); n=int(crc_cons.strict_supported.sum()); stats["hypergeom_p_strict_overlap"] = float(hypergeom.sf(a-1, M, K, n)) if M and K and n else None
    else: stats["hypergeom_p_strict_overlap"] = None
    pd.DataFrame([stats]).to_csv(OUT / "overlap_statistics.csv", index=False)

    # Legacy list and mechanistic terms: annotation only, never a filter.
    old_genes: set[str] = set()
    audit=pd.DataFrame()
    if LEGACY_PATH.exists():
        old=pd.read_csv(LEGACY_PATH); old_genes=set()
        for c in old.columns:
            if c.lower() in {"gene","symbol","gene_symbol","human_symbol"}: old_genes.update(old[c].astype(str).str.upper())
        base = dinp_cons.copy() if not dinp_cons.empty else pd.DataFrame(columns=["gene"])
        base["legacy_81_annotation"] = base.gene.isin(old_genes); audit=base[base.legacy_81_annotation].copy(); audit.to_csv(OUT / "old81_overlap_annotation.csv", index=False)
    legacy_pathway = pathways[pathways.pathway.astype(str).str.contains("|".join(LEGACY_TERMS), case=False, regex=True)].copy() if not pathways.empty else pd.DataFrame()
    legacy_pathway.to_csv(OUT / "legacy_pathway_audit.csv", index=False)
    legacy_targets = sorted(old_genes | PGE2_GENES | {f"PTGER{i}" for i in range(1,5)})
    legacy_rows=[]
    for gene in legacy_targets:
        dd = dinp_deg[dinp_deg.gene.eq(gene)] if not dinp_deg.empty else pd.DataFrame(); cc = crc_deg[crc_deg.gene.eq(gene)] if not crc_deg.empty else pd.DataFrame()
        legacy_rows.append({"gene":gene,"legacy_81":gene in old_genes,"mechanistic_panel":gene in (PGE2_GENES | {f"PTGER{i}" for i in range(1,5)}),"dinp_n_contrasts":len(dd),"dinp_sig_contrasts":int(dd.significant_primary.sum()) if not dd.empty else 0,"dinp_directions":";".join(sorted(set(dd.direction))) if not dd.empty else "","crc_n_cohorts":len(cc),"crc_sig_cohorts":int(cc.significant_primary.sum()) if not cc.empty else 0,"crc_directions":";".join(sorted(set(cc.direction))) if not cc.empty else ""})
    pd.DataFrame(legacy_rows).to_csv(OUT / "legacy_gene_audit.csv", index=False)
    # Candidate tiers and figures
    # Candidate ranking is downstream of both transcriptomic analyses.  The
    # old 81-gene list can annotate a row but cannot promote it into a tier.
    if not conv.empty:
        candidates = conv[conv.same_direction].copy()
        candidates["legacy_81_annotation"] = candidates.gene.isin(set(audit.gene) if not audit.empty else set())
        candidates["tier"] = np.select([candidates.strict_supported_DINP & candidates.strict_supported_CRC, candidates.strict_supported_DINP, candidates.strict_supported_CRC], ["A", "B", "B"], default="C")
        candidates["priority_score"] = candidates.n_significant_datasets_DINP * 3 + candidates.n_independent_datasets_DINP + candidates.n_significant_datasets_CRC * 2
    else:
        candidates = dinp_cons.copy()
        if not candidates.empty:
            candidates["legacy_81_annotation"] = candidates.gene.isin(set(audit.gene) if not audit.empty else set()); candidates["tier"] = np.where(candidates.strict_supported,"A",np.where(candidates.direction_supported,"B","C")); candidates["priority_score"] = candidates.n_significant_datasets*3 + candidates.n_independent_datasets + candidates.all_directions_same.astype(int)*2
    candidates.sort_values(["tier","priority_score"], ascending=[True,False]).to_csv(OUT / "candidate_prioritization.csv", index=False)
    if not dinp_cons.empty: plot_or_placeholder(OUT/"figures"/"Fig1_DINP_cross_tissue_effect_heatmap.png", "DINP representative effects", range(min(30,len(dinp_cons))), dinp_cons.head(30).median_log2FC, dinp_cons.head(30).gene)
    else: plot_or_placeholder(OUT/"figures"/"Fig1_DINP_cross_tissue_effect_heatmap.png", "DINP representative effects")
    plot_or_placeholder(OUT/"figures"/"Fig2_DINP_volcano_ovary_liver.png", "DINP volcano overview")
    plot_or_placeholder(OUT/"figures"/"Fig3_DINP_pathway_consensus.png", "DINP pathway consensus", range(min(20,len(path_cons))), path_cons.head(20).median_NES if not path_cons.empty else None, path_cons.head(20).pathway if not path_cons.empty else None)
    plot_or_placeholder(OUT/"figures"/"Fig4_legacy_pathway_audit.png", "Legacy pathway post-hoc audit")
    plot_or_placeholder(OUT/"figures"/"Fig5_CRC_expression_consensus.png", "CRC expression consensus")
    plot_or_placeholder(OUT/"figures"/"Fig6_directional_convergence.png", "DINP and CRC direction-aware convergence")
    plot_or_placeholder(OUT/"figures"/"Fig7_candidate_prioritization.png", "Candidate prioritization")
    plot_or_placeholder(OUT/"figures"/"Fig8_schematic.png", "Analysis schematic")

    manifest=pd.DataFrame(manifest_rows + crc_manifest); manifest.to_csv(OUT/"dataset_manifest.csv", index=False); manifest.to_json(OUT/"dataset_manifest.json", orient="records", indent=2)
    prov={"data_root":str(DATA),"files":[],"generated_utc":pd.Timestamp.utcnow().isoformat(),"python":sys.version}
    for p in [x for x in [*DATA.glob("*.gz"), *CRC_DATA.glob("*.gz")] if x.exists()]: prov["files"].append({"path":str(p),"size":p.stat().st_size,"sha256":sha256(p)})
    (OUT/"provenance"/"file_provenance.json").write_text(json.dumps(prov,indent=2),encoding="utf-8")
    (OUT/"ANALYSIS_AUDIT.md").write_text("# Analysis audit\n\n- Primary DINP gene selection used all mapped transcriptomic genes and strict multi-dataset direction/significance gates.\n- Legacy 81 genes were read only after consensus for annotation.\n- No raw FASTQ/count matrix was committed.\n- Processed GEO scales were tested with Welch/paired t-tests and marked secondary quality.\n- g:Orth mapping is Ensembl-backed; unmapped genes are retained in the mapping table.\n- Pathway method is preranked rank-sum over local Hallmark/Reactome GMTs; GO:BP and KEGG are not silently imputed when files are unavailable.\n",encoding="utf-8")
    write_report(manifest, pd.concat(all_dinp,ignore_index=True) if all_dinp else pd.DataFrame(), pathways, dinp_cons, crc_cons, conv, audit, stats)
    status = "PASS" if len(all_dinp)>=2 and len(crc_all)>=2 else "PARTIAL"
    print(json.dumps({"status":status,"out_dir":str(OUT),"dinp_contrasts":len(all_dinp),"crc_cohorts":len(crc_all),"dinp_strict_genes":stats.get("dinp_strict_n",0),"crc_strict_genes":stats.get("crc_strict_n",0),"joint_universe":stats.get("jointly_testable_universe",0)},ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
