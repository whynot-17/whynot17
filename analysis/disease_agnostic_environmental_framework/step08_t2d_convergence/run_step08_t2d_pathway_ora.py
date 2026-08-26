"""Run the frozen Tier A T2D pathway-convergence ORA.

The script uses the Step 7 CTD gene unions for the four Tier A exposure axes,
an 11-cluster CTD background, and g:Profiler's custom-background API. It saves
raw responses and recomputes raw hypergeometric P values plus one global
Benjamini–Hochberg correction across all returned axis x source x term tests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.request
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


DEFAULT_STEP7 = Path("analysis/disease_agnostic_environmental_framework/step07_genecard_convergence")
DEFAULT_OUT = Path("analysis/disease_agnostic_environmental_framework/step08_t2d_convergence")
GP_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "REAC", "KEGG"]
TIER_A = {"cluster_6", "cluster_5", "cluster_8", "cluster_11"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_frozen_inputs(step7: Path) -> tuple[pd.DataFrame, dict[str, set[str]], set[str]]:
    joint = pd.read_csv(step7 / "t2d_step7_joint_prioritization.csv")
    if set(joint.loc[joint["final_tier"].eq("Tier_A"), "cluster_id"]) != TIER_A:
        raise ValueError("Step 7 Tier A set is not the frozen four-cluster set")
    membership = pd.read_csv(step7 / "t2d_cluster_ctd_gene_membership.csv")
    if not {"cluster_id", "gene_symbol"}.issubset(membership.columns):
        raise ValueError("Step 7 membership file lacks cluster_id/gene_symbol")
    genes = {
        str(cid): set(group["gene_symbol"].astype(str).str.upper().str.strip())
        for cid, group in membership.groupby("cluster_id")
    }
    genes = {cid: vals - {"", "NAN", "NONE"} for cid, vals in genes.items()}
    background = set().union(*genes.values())
    if len(genes) != 11 or not all(genes[cid] for cid in genes):
        raise ValueError("Step 7 CTD gene membership is incomplete")
    return joint, genes, background


def query_gprofiler(
    cluster_id: str,
    query_genes: set[str],
    background: set[str],
    sources: list[str],
) -> dict[str, object]:
    payload = {
        "organism": "hsapiens",
        "query": sorted(query_genes),
        "sources": sources,
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": sorted(background),
        "all_results": True,
        "no_evidences": True,
        "no_iea": False,
        "measure_underrepresentation": False,
        "highlight": False,
    }
    req = urllib.request.Request(
        GP_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": "Step8-T2D-Convergence/1.0"},
        method="POST",
    )
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(req, timeout=180) as response:
                body = response.read()
            result = json.loads(body.decode("utf-8"))
            if "result" not in result or "meta" not in result:
                raise ValueError(f"g:Profiler response for {cluster_id} lacks result/meta")
            return {"payload": payload, "response": result}
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            last_error = exc
            if attempt < 3:
                time.sleep(15 * attempt)
    raise RuntimeError(f"g:Profiler request failed for {cluster_id} after 3 attempts: {last_error}")


def parse_results(
    cluster_id: str,
    response: dict[str, object],
    source_filter: str | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    rows = response.get("result", [])
    parsed: list[dict[str, object]] = []
    for item in rows:
        if not isinstance(item, dict):
            continue
        if source_filter is not None and item.get("source") != source_filter:
            continue
        try:
            m = int(item["effective_domain_size"])
            n = int(item["term_size"])
            q = int(item["query_size"])
            k = int(item["intersection_size"])
            raw_p = float(hypergeom.sf(k - 1, m, n, q))
        except (KeyError, TypeError, ValueError):
            continue
        intersections = item.get("intersections", [])
        if isinstance(intersections, list):
            intersection_text = ";".join(str(x) for x in intersections)
        else:
            intersection_text = str(intersections) if intersections else ""
        parsed.append({
            "cluster_id": cluster_id,
            "source": item.get("source", ""),
            "native": item.get("native", ""),
            "term": item.get("name", ""),
            "description": item.get("description", ""),
            "term_size": n,
            "query_size": q,
            "intersection_size": k,
            "effective_domain_size": m,
            "raw_hypergeom_p": raw_p,
            "gprofiler_adjusted_p": item.get("p_value", np.nan),
            "gprofiler_significant_at_threshold_1": item.get("significant", False),
            "intersection_genes": intersection_text,
        })
    meta = response.get("meta", {})
    genes_meta = meta.get("genes_metadata", {}) if isinstance(meta, dict) else {}
    result_meta = meta.get("result_metadata", {}) if isinstance(meta, dict) else {}
    query_meta = genes_meta.get("query", {}) if isinstance(genes_meta, dict) else {}
    query_record = query_meta.get("query_1", {}) if isinstance(query_meta, dict) else {}
    mapped_ids = query_record.get("ensgs", []) if isinstance(query_record, dict) else []
    failed_symbols = genes_meta.get("failed", []) if isinstance(genes_meta, dict) else []
    ambiguous = genes_meta.get("ambiguous", {}) if isinstance(genes_meta, dict) else {}
    duplicates = genes_meta.get("duplicates", []) if isinstance(genes_meta, dict) else []
    audit = {
        "cluster_id": cluster_id,
        "source_filter": source_filter,
        "n_returned_terms": len(parsed),
        "n_query_mapped": len(mapped_ids),
        "n_query_failed": len(failed_symbols),
        "failed_symbols": failed_symbols,
        "n_query_ambiguous": len(ambiguous),
        "n_query_duplicates": len(duplicates),
        "result_metadata": result_meta,
        "timestamp": meta.get("timestamp") if isinstance(meta, dict) else None,
        "gprofiler_version": meta.get("version") if isinstance(meta, dict) else None,
    }
    return pd.DataFrame(parsed), audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step7-dir", type=Path, default=DEFAULT_STEP7)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.output_dir / "raw_gprofiler"
    raw_dir.mkdir(parents=True, exist_ok=True)

    joint, cluster_genes, background = load_frozen_inputs(args.step7_dir)
    axis_audit = []
    all_results = []
    api_audits = []
    for cluster_id in sorted(TIER_A):
        query_genes = cluster_genes[cluster_id]
        axis_audit.append({
            "cluster_id": cluster_id,
            "final_tier": "Tier_A",
            "n_query_ctd_genes": len(query_genes),
            "n_background_ctd_genes": len(background),
            "query_source": "Step 7 frozen CTD cluster gene union",
            "background_source": "Step 7 frozen union of all 11 CTD cluster gene sets",
        })
        legacy_cache_path = raw_dir / f"gprofiler_{cluster_id}.json"
        for source in SOURCES:
            source_tag = source.replace(":", "_")
            cache_path = raw_dir / f"gprofiler_{cluster_id}_{source_tag}.json"
            if cache_path.exists() and cache_path.stat().st_size > 0:
                result = json.loads(cache_path.read_text(encoding="utf-8"))
                if "response" not in result or "result" not in result["response"]:
                    raise ValueError(f"Invalid cached g:Profiler response for {cluster_id}/{source}")
            elif legacy_cache_path.exists() and legacy_cache_path.stat().st_size > 0:
                result = json.loads(legacy_cache_path.read_text(encoding="utf-8"))
                if "response" not in result or "result" not in result["response"]:
                    raise ValueError(f"Invalid legacy g:Profiler response for {cluster_id}")
            else:
                result = query_gprofiler(cluster_id, query_genes, background, [source])
                cache_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
            parsed, audit = parse_results(cluster_id, result["response"], source)
            all_results.append(parsed)
            api_audits.append(audit)
            time.sleep(1.0)

    results = pd.concat(all_results, ignore_index=True)
    if results.empty:
        raise RuntimeError("g:Profiler returned no pathway terms for any Tier A axis")
    results["global_bh_fdr"] = multipletests(results["raw_hypergeom_p"].to_numpy(), method="fdr_bh")[1]
    results = results.sort_values(["global_bh_fdr", "raw_hypergeom_p", "cluster_id", "source"]).reset_index(drop=True)
    results.to_csv(args.output_dir / "t2d_step8_pathway_ora_all.csv", index=False)
    results[results["global_bh_fdr"] < 0.05].to_csv(args.output_dir / "t2d_step8_pathway_ora_significant.csv", index=False)
    pd.DataFrame(axis_audit).to_csv(args.output_dir / "t2d_step8_axis_input_audit.csv", index=False)
    (args.output_dir / "t2d_step8_gprofiler_query_audit.json").write_text(json.dumps(api_audits, indent=2), encoding="utf-8")

    sig = results[results["global_bh_fdr"] < 0.05]
    lines = [
        "# Step 8 — T2D pathway convergence",
        "",
        "- Status: **complete_pathway_ora**",
        f"- Tier A axes analyzed: **{len(TIER_A)}** — {', '.join(sorted(TIER_A))}",
        f"- Frozen CTD background genes: **{len(background):,}**",
        f"- Returned pathway tests: **{len(results):,}**",
        f"- Global BH-FDR < 0.05 terms: **{len(sig):,}**",
        "- Sources: **GO:BP, Reactome, KEGG**; custom background; human symbols.",
        "- Evidence annotations were not requested to control response size; term statistics and intersection sizes are unchanged.",
        "",
        "## Axis-level significant-term counts",
        "",
        "| Cluster | Significant terms | Top terms |",
        "|---|---:|---|",
    ]
    for cluster_id in sorted(TIER_A):
        x = sig[sig["cluster_id"].eq(cluster_id)].sort_values("global_bh_fdr").head(5)
        top = "; ".join(f"{r.term} ({r.source})" for r in x.itertuples()) or "none"
        lines.append(f"| {cluster_id} | {len(sig[sig['cluster_id'].eq(cluster_id)])} | {top} |")
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "This is direction-agnostic pathway over-representation analysis of CTD-associated genes. It does not establish pathway activation, exposure causality, or mediation of T2D. Transcriptomic directionality and interaction-network context are separate analyses.",
        "",
    ]
    (args.output_dir / "STEP8_T2D_PATHWAY_REPORT.md").write_text("\n".join(lines), encoding="utf-8")

    manifest = {
        "analysis": "Step 8 T2D pathway convergence",
        "status": "complete_pathway_ora",
        "tier_a_clusters": sorted(TIER_A),
        "n_tier_a": len(TIER_A),
        "n_background_ctd_genes": len(background),
        "n_returned_pathway_tests": len(results),
        "n_global_bh_fdr_lt_0_05": len(sig),
        "gene_set_sources": SOURCES,
        "organism": "hsapiens",
        "background_rule": "union of all 11 frozen Step 7 CTD cluster gene sets",
        "query_rule": "frozen CTD gene union for each Step 7 Tier A cluster",
        "ora": "one-sided hypergeometric; raw P reconstructed from API domain/term/query/intersection sizes",
        "multiple_testing": "single BH-FDR family across all returned Tier A x source x term tests",
        "gprofiler_api": GP_URL,
        "api_options": {"all_results": True, "no_evidences": True, "domain_scope": "custom", "significance_threshold_method": "fdr"},
        "gprofiler_audit": api_audits,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "canonical_outputs": [
            "t2d_step8_axis_input_audit.csv",
            "t2d_step8_pathway_ora_all.csv",
            "t2d_step8_pathway_ora_significant.csv",
            "t2d_step8_gprofiler_query_audit.json",
            "STEP8_T2D_PATHWAY_REPORT.md",
            "STEP8_MANIFEST.json",
        ],
    }
    (args.output_dir / "STEP8_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
