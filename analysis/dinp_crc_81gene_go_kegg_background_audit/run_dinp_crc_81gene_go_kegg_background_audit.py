#!/usr/bin/env python3
"""Audit GO/KEGG enrichment of the frozen 81-gene DINP–CRC intersection.

Primary inferential question:
    Within the genes that were eligible on the DINP exposure side, are the
    DINP∩CRC genes selectively enriched for particular functions?

Accordingly, the primary custom universe is the complete frozen DINP exposure
set.  Three additional universes are retained for interpretation/sensitivity:
CRC disease union, legacy exposure∪CRC union, and g:Profiler's annotated human
genome background.  The genome analysis is descriptive and must not be used as
an intersection-specific null model.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs"
OUT_DIR = Path(__file__).resolve().parent / "outputs"
API_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG"]
UA = "whynot17-dinp-crc-background-audit/1.0"


def unique_symbols(values: Iterable[object]) -> List[str]:
    return sorted({str(v).strip().upper() for v in values if pd.notna(v) and str(v).strip()})


def sha256_lines(values: List[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def bh_fixed_family(p: pd.Series, family_size: int) -> pd.Series:
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna() & p.between(0, 1)
    if not valid.any() or family_size <= 0:
        return out
    vals = p.loc[valid].to_numpy(float)
    order = np.argsort(vals, kind="mergesort")
    ranked = vals[order] * float(family_size) / np.arange(1, len(vals) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.clip(ranked, 0, 1)
    tmp = np.empty_like(ranked)
    tmp[order] = ranked
    out.loc[valid] = tmp
    return out


def request_gprofiler(query: List[str], label: str, background: Optional[List[str]]) -> Dict:
    payload = {
        "organism": "hsapiens",
        "query": query,
        "sources": SOURCES,
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "no_evidences": True,
        "all_results": True,
        "ordered": False,
        "combined": False,
        "measure_underrepresentation": False,
    }
    if background is None:
        payload["domain_scope"] = "annotated"
    else:
        payload["domain_scope"] = "custom"
        payload["background"] = background
    r = requests.post(API_URL, json=payload, headers={"User-Agent": UA}, timeout=300)
    r.raise_for_status()
    ans = r.json()
    ans["_payload"] = payload
    ans["_status"] = r.status_code
    ans["_sha256"] = hashlib.sha256(r.content).hexdigest()
    ans["_label"] = label
    return ans


def flatten(resp: Dict) -> pd.DataFrame:
    meta = resp.get("meta", {})
    rm = meta.get("result_metadata", {})
    label = resp["_label"]
    rows = []
    for x in resp.get("result", []):
        source = x.get("source")
        N = int(x.get("effective_domain_size", 0) or 0)
        K = int(x.get("term_size", 0) or 0)
        n = int(x.get("query_size", 0) or 0)
        k = int(x.get("intersection_size", 0) or 0)
        raw = float(hypergeom.sf(k - 1, N, K, n)) if N > 0 and k > 0 else np.nan
        rows.append({
            "background": label,
            "source": source,
            "term_id": x.get("native"),
            "term_name": x.get("name"),
            "effective_domain_size": N,
            "term_size": K,
            "query_size": n,
            "intersection_size": k,
            "raw_p_hypergeom": raw,
            "gprofiler_fdr": x.get("p_value"),
            "tested_term_family_size": int(rm.get(source, {}).get("number_of_terms", 0) or 0),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    total_family = sum(int(v.get("number_of_terms", 0) or 0) for v in rm.values())
    df["BH_FDR_all_GO_KEGG"] = bh_fixed_family(df["raw_p_hypergeom"], total_family)
    df["BH_FDR_within_source"] = np.nan
    for source, idx in df.groupby("source").groups.items():
        fam = int(rm.get(source, {}).get("number_of_terms", 0) or 0)
        df.loc[idx, "BH_FDR_within_source"] = bh_fixed_family(df.loc[idx, "raw_p_hypergeom"], fam)
    return df.sort_values(["BH_FDR_all_GO_KEGG", "raw_p_hypergeom"]).reset_index(drop=True)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    required = {
        "intersection": INPUT_DIR / "dinp_crc_intersection.csv",
        "exposure": INPUT_DIR / "dinp_exposure_gene_matrix.csv",
        "crc": INPUT_DIR / "crc_gene_matrix.csv",
    }
    missing = [str(p) for p in required.values() if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing convergence outputs: " + "; ".join(missing))

    q = unique_symbols(pd.read_csv(required["intersection"])["gene_symbol"])
    exposure = unique_symbols(pd.read_csv(required["exposure"])["gene_symbol"])
    crc = unique_symbols(pd.read_csv(required["crc"])["gene_symbol"])
    legacy_union = sorted(set(exposure) | set(crc))

    if len(q) != 81:
        raise ValueError(f"Expected 81 query genes, found {len(q)}")
    if not set(q).issubset(exposure):
        raise ValueError("Query is not a subset of the DINP exposure universe")
    if not set(q).issubset(crc):
        raise ValueError("Query is not a subset of the CRC disease universe")

    configs = [
        ("primary_dinp_exposure_universe", exposure),
        ("sensitivity_crc_union", crc),
        ("legacy_exposure_crc_union", legacy_union),
        ("descriptive_annotated_genome", None),
    ]

    manifest = {
        "analysis": "DINP–CRC 81-gene GO/KEGG background audit",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "query_n": len(q),
        "dinp_exposure_n": len(exposure),
        "crc_union_n": len(crc),
        "legacy_union_n": len(legacy_union),
        "query_fraction_of_dinp": len(q) / len(exposure),
        "primary_null": "DINP exposure universe; tests whether CRC-overlapping DINP genes are selectively enriched relative to all DINP-eligible genes",
        "descriptive_genome_warning": "Annotated-genome background is descriptive only and is not an intersection-specific null model",
        "hashes": {"query": sha256_lines(q), "exposure": sha256_lines(exposure), "crc": sha256_lines(crc)},
        "runs": {},
    }

    frames = []
    for label, bg in configs:
        resp = request_gprofiler(q, label, bg)
        public_resp = {k: v for k, v in resp.items() if not k.startswith("_")}
        (OUT_DIR / f"api_response_{label}.json").write_text(json.dumps(public_resp, indent=2), encoding="utf-8")
        (OUT_DIR / f"request_payload_{label}.json").write_text(json.dumps(resp["_payload"], indent=2), encoding="utf-8")
        df = flatten(resp)
        df.to_csv(OUT_DIR / f"enrichment_{label}.csv", index=False)
        frames.append(df)
        manifest["runs"][label] = {
            "background_n": None if bg is None else len(bg),
            "response_sha256": resp["_sha256"],
            "returned_terms": int(len(df)),
            "significant_global_fdr": int((df["BH_FDR_all_GO_KEGG"] < 0.05).sum()) if len(df) else 0,
            "significant_within_source_fdr": int((df["BH_FDR_within_source"] < 0.05).sum()) if len(df) else 0,
        }

    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_csv(OUT_DIR / "enrichment_all_backgrounds.csv", index=False)

    key_pattern = r"prostaglandin|arachidonic|eicosanoid|inflamm|lipid|PPAR|receptor"
    key = all_df[all_df["term_name"].fillna("").str.contains(key_pattern, case=False, regex=True)].copy()
    key.to_csv(OUT_DIR / "key_pathway_background_comparison.csv", index=False)

    summary = []
    for label, _ in configs:
        sub = all_df[all_df.background == label]
        summary.append({
            "background": label,
            "n_terms": len(sub),
            "n_global_fdr_lt_0_05": int((sub.BH_FDR_all_GO_KEGG < 0.05).sum()) if len(sub) else 0,
            "best_term": sub.iloc[0].term_name if len(sub) else None,
            "best_global_fdr": float(sub.iloc[0].BH_FDR_all_GO_KEGG) if len(sub) else np.nan,
        })
    pd.DataFrame(summary).to_csv(OUT_DIR / "background_audit_summary.csv", index=False)
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    report = f"""# DINP–CRC enrichment background audit\n\n## Why this audit was required\nThe frozen query contains {len(q)} genes, while the entire DINP exposure universe contains only {len(exposure)} genes. Thus {len(q)}/{len(exposure)} ({100*len(q)/len(exposure):.1f}%) of DINP-eligible genes already overlap the CRC disease union. A very large CRC/union/genome background can therefore make pathway structure already present on the DINP side appear extremely significant without demonstrating intersection-specific enrichment.\n\n## Primary analysis\nThe primary inferential universe is now the full DINP exposure set ({len(exposure)} genes). This asks whether the CRC-overlapping subset is functionally selective relative to genes that could have entered the intersection from the DINP side.\n\n## Sensitivity analyses\n- CRC disease union: {len(crc)} genes.\n- Legacy DINP∪CRC union: {len(legacy_union)} genes.\n- g:Profiler annotated human genome: descriptive functional annotation only.\n\n## Interpretation rule\nA term that is significant only against CRC/union/genome but not against the DINP exposure universe should be described as a feature of the 81-gene set or of DINP-related biology, not as a CRC-convergence-specific enrichment.\n\nSee `key_pathway_background_comparison.csv` for prostaglandin/arachidonic/eicosanoid/inflammatory/lipid terms across all four backgrounds.\n"""
    (OUT_DIR / "BACKGROUND_AUDIT_REPORT.md").write_text(report, encoding="utf-8")

    print("DINP–CRC GO/KEGG BACKGROUND AUDIT: PASS")
    print(f"Query: {len(q)}; DINP universe: {len(exposure)}; CRC universe: {len(crc)}; legacy union: {len(legacy_union)}")
    print(f"Intersection fraction of DINP universe: {len(q)/len(exposure):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
