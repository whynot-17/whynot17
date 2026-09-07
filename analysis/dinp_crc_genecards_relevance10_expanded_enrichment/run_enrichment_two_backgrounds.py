#!/usr/bin/env python3
"""Run source-preserving GO/KEGG/Reactome ORA with two frozen backgrounds.

Primary background: all 881 archived ordinary GeneCards CRC genes with
RelevanceScore >= 10. Sensitivity background: all 93 DINP multi-source genes.
The query is the 18-gene expanded DINP-GeneCards intersection. Open Targets is
not used. Each background is analyzed as its own multiplicity family.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
INPUT_DIR = REPO_ROOT / "analysis" / "dinp_crc_multi_database_target_convergence"
OUT = HERE / "outputs"
QUERY_FILE = INPUT_DIR / "outputs" / "dinp_crc_genecards_relevance10_expanded_intersection.csv"
EXPOSURE_FILE = INPUT_DIR / "outputs" / "dinp_exposure_gene_matrix_ctd_toxcast_tox21_t3db.csv"
GENECARDS_FILE = INPUT_DIR / "outputs" / "source_records" / "genecards_crc_top2000_reference.csv"
API_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "KEGG", "REAC"]
USER_AGENT = "whynot17-dinp-crc-expanded-enrichment/1.0"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_genes(genes: list[str]) -> str:
    return hashlib.sha256("\n".join(genes).encode()).hexdigest()


def repo_relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def unique_symbols(values: list[str]) -> list[str]:
    return sorted({str(value).strip().upper() for value in values if str(value).strip()})


def log_combination(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(population: int, successes: int, draws: int, observed: int) -> float:
    """P[X >= observed] for a hypergeometric random variable."""
    if population <= 0 or successes < 0 or draws < 0 or observed < 0:
        return float("nan")
    if successes > population or draws > population:
        return float("nan")
    if observed == 0:
        return 1.0
    lo = max(observed, draws - (population - successes))
    hi = min(successes, draws)
    if lo > hi:
        return 0.0
    denominator = log_combination(population, draws)
    log_terms = [
        log_combination(successes, x)
        + log_combination(population - successes, draws - x)
        - denominator
        for x in range(lo, hi + 1)
    ]
    peak = max(log_terms)
    return min(1.0, math.exp(peak) * sum(math.exp(value - peak) for value in log_terms))


def bh_fixed(values: list[float], family_size: int) -> list[float]:
    output = [float("nan")] * len(values)
    valid = [(index, value) for index, value in enumerate(values) if math.isfinite(value) and 0 <= value <= 1]
    if not valid or family_size <= 0:
        return output
    ordered = sorted(valid, key=lambda item: (item[1], item[0]))
    running = 1.0
    adjusted = [0.0] * len(ordered)
    for position in range(len(ordered) - 1, -1, -1):
        rank = position + 1
        running = min(running, ordered[position][1] * family_size / rank)
        adjusted[position] = min(1.0, running)
    for (index, _), value in zip(ordered, adjusted):
        output[index] = value
    return output


def request_gprofiler(query: list[str], background: list[str]) -> tuple[dict[str, Any], bytes, dict[str, Any]]:
    payload = {
        "organism": "hsapiens",
        "query": query,
        "background": background,
        "sources": SOURCES,
        "domain_scope": "custom",
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "all_results": True,
        "no_evidences": True,
        "ordered": False,
        "combined": False,
        "measure_underrepresentation": False,
    }
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        API_URL,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        method="POST",
    )
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                raw = response.read()
            result = json.loads(raw.decode("utf-8"))
            if "result" not in result or "meta" not in result:
                raise ValueError("g:Profiler response lacks result/meta")
            return result, raw, payload
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ValueError) as exc:
            last_error = exc
            if attempt < 3:
                time.sleep(10 * attempt)
    raise RuntimeError(f"g:Profiler request failed after 3 attempts: {last_error}")


def parse_response(result: dict[str, Any], background_name: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    metadata = result.get("meta", {}).get("result_metadata", {}) or {}
    rows: list[dict[str, Any]] = []
    for item in result.get("result", []):
        if not isinstance(item, dict) or item.get("source") not in SOURCES:
            continue
        try:
            domain = int(item.get("effective_domain_size", 0) or 0)
            term_size = int(item.get("term_size", 0) or 0)
            query_size = int(item.get("query_size", 0) or 0)
            overlap = int(item.get("intersection_size", 0) or 0)
        except (TypeError, ValueError):
            continue
        raw_p = hypergeom_sf(domain, term_size, query_size, overlap)
        intersections = item.get("intersections", [])
        rows.append(
            {
                "background": background_name,
                "source": item.get("source", ""),
                "term_id": item.get("native", ""),
                "term_name": item.get("name", ""),
                "description": item.get("description", ""),
                "term_size": term_size,
                "effective_domain_size": domain,
                "query_size": query_size,
                "intersection_size": overlap,
                "raw_hypergeom_p": raw_p,
                "gprofiler_p_value": item.get("p_value", ""),
                "gprofiler_significant": item.get("significant", ""),
                "intersection_genes": ";".join(str(value) for value in intersections) if isinstance(intersections, list) else str(intersections),
                "parents": ";".join(str(value) for value in item.get("parents", []) or []),
            }
        )

    source_family_sizes = {
        source: int((metadata.get(source) or {}).get("number_of_terms", 0) or 0)
        for source in SOURCES
    }
    audit = {
        "background": background_name,
        "n_returned_terms": len(rows),
        "source_family_sizes": source_family_sizes,
        "total_source_family_size": sum(source_family_sizes.values()),
        "timestamp": result.get("meta", {}).get("timestamp"),
        "gprofiler_version": result.get("meta", {}).get("version"),
        "genes_metadata": result.get("meta", {}).get("genes_metadata", {}),
    }
    return rows, audit


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw_dir = OUT / "raw_gprofiler"
    raw_dir.mkdir(parents=True, exist_ok=True)
    query_rows = read_csv(QUERY_FILE)
    exposure_rows = read_csv(EXPOSURE_FILE)
    cards_rows = read_csv(GENECARDS_FILE)

    query = unique_symbols([row.get("gene_symbol", "") for row in query_rows])
    exposure_background = unique_symbols([row.get("gene_symbol", "") for row in exposure_rows])
    cards_background = unique_symbols(
        [
            row.get("gene_symbol") or row.get("GeneSymbol", "")
            for row in cards_rows
            if _numeric(row.get("RelevanceScore", "")) is not None and float(row["RelevanceScore"]) >= 10
        ]
    )
    if len(query) != 18:
        raise ValueError(f"Expected 18 query genes, found {len(query)}")
    if len(exposure_background) != 93:
        raise ValueError(f"Expected 93 exposure background genes, found {len(exposure_background)}")
    if len(cards_background) != 881:
        raise ValueError(f"Expected 881 GeneCards background genes, found {len(cards_background)}")
    if not set(query).issubset(exposure_background) or not set(query).issubset(cards_background):
        raise ValueError("Query is not a subset of both frozen backgrounds")

    backgrounds = {
        "primary_genecards881": cards_background,
        "sensitivity_dinp93": exposure_background,
    }
    all_rows: list[dict[str, Any]] = []
    api_audit: list[dict[str, Any]] = []
    source_counts: dict[str, dict[str, int]] = {}
    for background_name, background in backgrounds.items():
        raw_path = raw_dir / f"gprofiler_{background_name}.json"
        payload_path = OUT / f"request_payload_{background_name}.json"
        if raw_path.exists() and raw_path.stat().st_size > 0 and payload_path.exists() and payload_path.stat().st_size > 0:
            raw = raw_path.read_bytes()
            result = json.loads(raw.decode("utf-8"))
            payload = json.loads(payload_path.read_text(encoding="utf-8"))
        else:
            result, raw, payload = request_gprofiler(query, background)
            raw_path.write_bytes(raw)
            payload_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        rows, audit = parse_response(result, background_name)
        raw_p = [float(row["raw_hypergeom_p"]) for row in rows]
        family_size = int(audit["total_source_family_size"] or len(rows))
        global_fdr = bh_fixed(raw_p, family_size)
        for row, value in zip(rows, global_fdr):
            row["BH_FDR_global_source_family"] = value
            row["BH_FDR_global_significant"] = bool(math.isfinite(value) and value < 0.05)
        for source in SOURCES:
            indices = [index for index, row in enumerate(rows) if row["source"] == source]
            source_family_size = int(audit["source_family_sizes"].get(source, 0) or len(indices))
            values = bh_fixed([raw_p[index] for index in indices], source_family_size)
            for index, value in zip(indices, values):
                rows[index]["BH_FDR_within_source"] = value
                rows[index]["BH_FDR_within_source_significant"] = bool(math.isfinite(value) and value < 0.05)
        all_rows.extend(rows)
        api_audit.append(
            {
                **audit,
                "raw_response": repo_relative(raw_path),
                "raw_response_sha256": sha256_file(raw_path),
                "request_payload": repo_relative(payload_path),
                "query_count": len(query),
                "background_count": len(background),
                "query_mapped_count": len((audit.get("genes_metadata") or {}).get("query", {}).get("query_1", {}).get("ensgs", []) or []),
                "failed_symbols": (audit.get("genes_metadata") or {}).get("failed", []) or [],
            }
        )
        source_counts[background_name] = {
            source: sum(1 for row in rows if row["source"] == source and row["BH_FDR_global_significant"])
            for source in SOURCES
        }

    fields = list(all_rows[0]) if all_rows else ["background"]
    write_csv(OUT / "dinp_crc_enrichment_all_backgrounds.csv", all_rows, fields)
    for background_name in backgrounds:
        subset = [row for row in all_rows if row["background"] == background_name]
        write_csv(OUT / f"dinp_crc_enrichment_{background_name}.csv", subset, fields)
        for source, tag in [("GO:BP", "go_bp"), ("KEGG", "kegg"), ("REAC", "reactome")]:
            write_csv(OUT / f"{tag}_{background_name}.csv", [row for row in subset if row["source"] == source], fields)

    gene_list_path = OUT / "dinp_crc_18_gene_list.csv"
    write_csv(gene_list_path, [{"gene_symbol": gene} for gene in query], ["gene_symbol"])
    significant = [row for row in all_rows if row["BH_FDR_global_significant"]]
    counts = {
        background_name: {
            "returned_terms": sum(1 for row in all_rows if row["background"] == background_name),
            "global_fdr_lt_0_05": sum(1 for row in significant if row["background"] == background_name),
            "global_significant_by_source": source_counts[background_name],
            "within_source_fdr_lt_0_05": {
                source: sum(
                    1
                    for row in all_rows
                    if row["background"] == background_name
                    and row["source"] == source
                    and row["BH_FDR_within_source_significant"]
                )
                for source in SOURCES
            },
        }
        for background_name in backgrounds
    }
    manifest = {
        "analysis": "Expanded 18-gene DINP-CRC GeneCards-threshold GO/KEGG/Reactome ORA",
        "run_utc": utc_now(),
        "query": {
            "path": repo_relative(QUERY_FILE),
            "sha256": sha256_file(QUERY_FILE),
            "count": len(query),
            "gene_list_sha256": sha256_genes(query),
        },
        "backgrounds": {
            "primary_genecards881": {
                "path": repo_relative(GENECARDS_FILE),
                "sha256": sha256_file(GENECARDS_FILE),
                "rule": "ordinary archived CRC GeneCards reference, RelevanceScore >= 10",
                "count": len(cards_background),
                "gene_list_sha256": sha256_genes(cards_background),
            },
            "sensitivity_dinp93": {
                "path": repo_relative(EXPOSURE_FILE),
                "sha256": sha256_file(EXPOSURE_FILE),
                "rule": "all genes in source-preserving CTD + ToxCast + Tox21 + T3DB DINP matrix",
                "count": len(exposure_background),
                "gene_list_sha256": sha256_genes(exposure_background),
            },
        },
        "sources": SOURCES,
        "organism": "hsapiens",
        "method": "one-sided hypergeometric ORA with custom background",
        "multiple_testing": "each background is a separate family; global BH uses the sum of g:Profiler source term-family sizes; within-source BH retained",
        "api": {"endpoint": API_URL, "gprofiler_audit": api_audit},
        "counts": counts,
        "firewall": {"open_targets_used": False, "crc_outcome_statistics_used": False, "enrichment_used_for_gene_selection": False},
        "interpretation_boundary": "Enrichment is direction-agnostic and does not establish DINP causality, pathway activation, or mediation. The primary GeneCards background is limited to the archived ordinary top-2000 reference.",
        "canonical_outputs": [
            "dinp_crc_18_gene_list.csv",
            "dinp_crc_enrichment_all_backgrounds.csv",
            "dinp_crc_enrichment_primary_genecards881.csv",
            "dinp_crc_enrichment_sensitivity_dinp93.csv",
            "dinp_crc_enrichment_manifest.json",
            "dinp_crc_enrichment_summary.md",
        ],
    }
    (OUT / "dinp_crc_enrichment_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = [
        "# Expanded 18-gene DINP–CRC enrichment",
        "",
        "- Method: one-sided hypergeometric ORA with custom backgrounds",
        "- Query: the expanded 18-gene DINP exposure ∩ GeneCards CRC threshold set",
        "- Sources: GO:BP, KEGG, Reactome",
        "- Open Targets: not used",
        "- Each background is a separate BH-FDR family",
        "",
        "## Frozen inputs",
        "",
        "| Analysis | Query | Background |",
        "|---|---:|---:|",
        f"| Primary — GeneCards high-relevance background | {len(query)} | {len(cards_background)} |",
        f"| Sensitivity — DINP multi-source background | {len(query)} | {len(exposure_background)} |",
        "",
        "## Results",
        "",
        "| Background | Returned terms | Global BH-FDR <0.05 | GO:BP | KEGG | Reactome |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for background_name, label in [("primary_genecards881", "Primary GeneCards 881"), ("sensitivity_dinp93", "Sensitivity DINP 93")]:
        item = counts[background_name]
        by_source = item["global_significant_by_source"]
        lines.append(
            f"| {label} | {item['returned_terms']} | {item['global_fdr_lt_0_05']} | {by_source['GO:BP']} | {by_source['KEGG']} | {by_source['REAC']} |"
        )
    lines += [
        "",
        "## Top terms by background",
        "",
    ]
    for background_name, label in [("primary_genecards881", "Primary GeneCards 881"), ("sensitivity_dinp93", "Sensitivity DINP 93")]:
        lines.append(f"### {label}")
        rows = sorted(
            [row for row in all_rows if row["background"] == background_name],
            key=lambda row: (float(row["BH_FDR_global_source_family"]), float(row["raw_hypergeom_p"])),
        )[:10]
        if not rows:
            lines.append("No terms returned.")
        else:
            for row in rows:
                lines.append(
                    f"- `{row['source']}` {row['term_name']} — overlap {row['intersection_size']}/{row['query_size']}; global BH-FDR={float(row['BH_FDR_global_source_family']):.4g}"
                )
        lines.append("")
    lines += [
        "## Interpretation boundary",
        "",
        "This analysis tests over-representation only. It is direction-agnostic and does not show pathway activation, DINP causality, or mediation. The GeneCards primary background is the available archived ordinary CRC top-2000 reference rather than a full GeneCards export.",
        "",
        "## Files",
        "",
        "See the manifest for input hashes, g:Profiler version/timestamps, raw responses, and source-specific result files.",
    ]
    (OUT / "dinp_crc_enrichment_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("DINP-CRC ENRICHMENT: PASS")
    print(json.dumps({"query": len(query), "primary_background": len(cards_background), "sensitivity_background": len(exposure_background), "counts": counts}, indent=2, ensure_ascii=False))


def _numeric(value: Any) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


if __name__ == "__main__":
    main()
