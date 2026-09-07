#!/usr/bin/env python3
"""GO/KEGG ORA for the 16-gene DINP–CRC GeneCards-threshold intersection.

The query is the CTD DINP human gene set intersected with GeneCards CRC genes
having RelevanceScore >= 10.  The ORA background is the complete 86-gene CTD
DINP human gene set.  Open Targets and all outcome statistics are excluded.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUT = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs" / "dinp_crc_genecards_relevance10_intersection.csv"
EXPOSURE = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs" / "dinp_exposure_gene_matrix.csv"
OUT = HERE / "outputs"
API_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG"]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def unique(values: list[str]) -> list[str]:
    return sorted({v.strip().upper() for v in values if v and v.strip()})


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(M: int, n: int, N: int, k: int) -> float:
    """P[X >= k] for X~Hypergeom(M, n, N), computed stably in log-space."""
    if min(M, n, N) < 0 or n > M or N > M or k <= 0:
        return 1.0 if k <= 0 else float("nan")
    lo = max(k, N - (M - n))
    hi = min(n, N)
    if lo > hi:
        return 0.0
    denom = log_comb(M, N)
    terms = [log_comb(n, x) + log_comb(M - n, N - x) - denom for x in range(lo, hi + 1)]
    peak = max(terms)
    return min(1.0, math.exp(peak) * sum(math.exp(v - peak) for v in terms))


def bh_fixed(p_values: list[float], family_size: int) -> list[float]:
    out = [float("nan")] * len(p_values)
    valid = [(i, p) for i, p in enumerate(p_values) if math.isfinite(p) and 0 <= p <= 1]
    if not valid or family_size <= 0:
        return out
    ordered = sorted(valid, key=lambda x: (x[1], x[0]))
    adjusted = [0.0] * len(ordered)
    running = 1.0
    for j in range(len(ordered) - 1, -1, -1):
        rank = j + 1
        running = min(running, ordered[j][1] * family_size / rank)
        adjusted[j] = min(1.0, running)
    for (idx, _), value in zip(ordered, adjusted):
        out[idx] = value
    return out


def request(query: list[str], background: list[str]) -> tuple[dict, bytes]:
    payload = {
        "organism": "hsapiens",
        "query": query,
        "background": background,
        "sources": SOURCES,
        "domain_scope": "custom",
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "no_evidences": True,
        "all_results": True,
        "ordered": False,
        "combined": False,
        "measure_underrepresentation": False,
    }
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        API_URL,
        data=body,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "whynot17-dinp-crc-genecards-relevance10-ora/1.0",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=300) as resp:
        raw = resp.read()
        return json.loads(raw.decode("utf-8")), raw


def main() -> None:
    query_rows = read_csv(INPUT)
    exposure_rows = read_csv(EXPOSURE)
    query = unique([r.get("gene_symbol", "") for r in query_rows])
    background = unique([r.get("gene_symbol", "") for r in exposure_rows if r.get("CTD", "0") == "1"])
    if len(query) != 16:
        raise ValueError(f"Expected 16 GeneCards-threshold genes; found {len(query)}")
    if len(background) != 86:
        raise ValueError(f"Expected 86 CTD DINP human genes; found {len(background)}")
    if not set(query).issubset(background):
        raise ValueError("Every query gene must be present in the CTD background")

    result, raw = request(query, background)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "api_response.json").write_bytes(raw)
    payload = {
        "organism": "hsapiens",
        "query": query,
        "background": background,
        "sources": SOURCES,
        "domain_scope": "custom",
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "all_results": True,
    }
    (OUT / "request_payload.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    metadata = result.get("meta", {}).get("result_metadata", {})
    rows: list[dict[str, object]] = []
    for rec in result.get("result", []):
        source = rec.get("source", "")
        M = int(rec.get("effective_domain_size", 0) or 0)
        n = int(rec.get("term_size", 0) or 0)
        N = int(rec.get("query_size", 0) or 0)
        k = int(rec.get("intersection_size", 0) or 0)
        raw_p = hypergeom_sf(M, n, N, k) if M and n and N and k else float("nan")
        rows.append(
            {
                "source": source,
                "term_id": rec.get("native", ""),
                "term_name": rec.get("name", ""),
                "description": rec.get("description", ""),
                "term_size": n,
                "effective_domain_size": M,
                "query_size": N,
                "intersection_size": k,
                "precision": rec.get("precision", ""),
                "recall": rec.get("recall", ""),
                "raw_p_hypergeom": raw_p,
                "gprofiler_fdr": rec.get("p_value", ""),
                "gprofiler_significant": rec.get("significant", ""),
                "parents": ";".join(rec.get("parents", []) or []),
                "tested_term_family_size": int(metadata.get(source, {}).get("number_of_terms", 0) or 0),
            }
        )
    family_all = sum(int(v.get("number_of_terms", 0) or 0) for v in metadata.values())
    raw_p_values = [float(r["raw_p_hypergeom"]) for r in rows]
    all_fdr = bh_fixed(raw_p_values, family_all)
    for row, value in zip(rows, all_fdr):
        row["BH_FDR_all_GO_KEGG"] = value
        row["BH_FDR_all_GO_KEGG_significant"] = bool(math.isfinite(value) and value < 0.05)
    for source in SOURCES:
        idx = [i for i, row in enumerate(rows) if row["source"] == source]
        values = bh_fixed([raw_p_values[i] for i in idx], int(metadata.get(source, {}).get("number_of_terms", 0) or 0))
        for i, value in zip(idx, values):
            rows[i]["BH_FDR_within_source"] = value
            rows[i]["BH_FDR_within_source_significant"] = bool(math.isfinite(value) and value < 0.05)

    fields = list(rows[0]) if rows else ["source"]
    with (OUT / "enrichment_all_GO_KEGG.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for source, name in [("GO:BP", "go_bp"), ("GO:MF", "go_mf"), ("GO:CC", "go_cc"), ("KEGG", "kegg")]:
        subset = [row for row in rows if row["source"] == source]
        with (OUT / f"{name}.csv").open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(subset)

    significant = [r for r in rows if r.get("BH_FDR_all_GO_KEGG_significant")]
    sig_by_source = {source: sum(r["source"] == source for r in significant) for source in SOURCES}
    source_specific_sig_by_source = {
        source: sum(
            r["source"] == source
            and isinstance(r.get("BH_FDR_within_source"), float)
            and math.isfinite(r["BH_FDR_within_source"])
            and r["BH_FDR_within_source"] < 0.05
            for r in rows
        )
        for source in SOURCES
    }
    lines = [
        "# GO/KEGG enrichment of the 16-gene DINP–CRC GeneCards-threshold intersection",
        "",
        "- Method: ORA with g:Profiler; custom CTD-DINP background",
        "- Query: CTD DINP human genes ∩ GeneCards CRC `RelevanceScore >= 10`",
        "- Background: all 86 CTD DINP human genes",
        "- Open Targets: not used",
        "- Multiple testing: BH-FDR across the GO:BP, GO:MF, GO:CC and KEGG family; source-specific FDR also retained",
        "",
        f"- Query genes: **{len(query)}**",
        f"- Background genes: **{len(background)}**",
        f"- Returned overlapping terms: **{len(rows)}**",
        f"- Global BH-FDR <0.05 terms: **{len(significant)}**",
        "- Significant terms by source: " + ", ".join(f"{k}={v}" for k, v in sig_by_source.items()),
        "- Source-specific BH-FDR <0.05 terms: " + ", ".join(f"{k}={v}" for k, v in source_specific_sig_by_source.items()),
        "",
        "## Interpretation boundary",
        "",
        "This is functional over-representation, not evidence of DINP causality, direction of regulation, or a validated mechanism. The GeneCards input is the available archived ordinary CRC top-2000 reference; the result is not presented as a full-ranking analysis.",
        "",
        "## Files",
        "",
        "- `enrichment_all_GO_KEGG.csv`",
        "- `go_bp.csv`, `go_mf.csv`, `go_cc.csv`, `kegg.csv`",
        "- `api_response.json`, `request_payload.json`, `manifest.json`",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest = {
        "analysis": "DINP-CRC GeneCards relevance-threshold intersection GO/KEGG ORA",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "query_count": len(query),
        "background_count": len(background),
        "query_sha256": hashlib.sha256("\n".join(query).encode()).hexdigest(),
        "background_sha256": hashlib.sha256("\n".join(background).encode()).hexdigest(),
        "sources": SOURCES,
        "threshold": {"field": "GeneCards RelevanceScore", "operator": ">=", "value": 10},
        "method": "ORA; hypergeometric enrichment with custom CTD-DINP background",
        "multiple_testing": {"family": "GO:BP + GO:MF + GO:CC + KEGG", "family_size": family_all, "adjustment": "BH"},
        "inputs": {"query": str(INPUT), "query_sha256": sha256(INPUT), "background": str(EXPOSURE), "background_sha256": sha256(EXPOSURE)},
        "firewall": {"open_targets_used": False, "crc_outcome_statistics_used": False, "enrichment_used_for_gene_selection": False},
        "api": {"endpoint": API_URL, "response_sha256": hashlib.sha256(raw).hexdigest(), "version": result.get("meta", {}).get("version"), "timestamp": result.get("meta", {}).get("timestamp")},
        "counts": {
            "returned_terms": len(rows),
            "global_fdr_lt_0_05": len(significant),
            "significant_by_source": sig_by_source,
            "source_specific_fdr_lt_0_05_by_source": source_specific_sig_by_source,
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest["counts"], ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
