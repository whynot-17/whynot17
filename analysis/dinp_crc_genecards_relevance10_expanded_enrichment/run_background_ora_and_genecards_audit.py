#!/usr/bin/env python3
"""Independent audit of the expanded DINP--CRC ORA backgrounds.

This audit does not overwrite the historical enrichment outputs.  It rebuilds
the term counts from a frozen g:Profiler Ensembl GMT snapshot, recomputes the
one-sided hypergeometric tail locally, and compares the result to the stored
g:Profiler responses.  It also audits the archived GeneCards CRC top-2000
reference used to define the RelevanceScore >= 10 background.

The historical API requests used ``domain_scope=custom``.  Therefore the
g:Profiler-equivalent local calculation uses the mapped Ensembl domain
(93 symbols can expand to 100 Ensembl IDs).  A second calculation uses only
background IDs that occur in the relevant source's fixed gene-set mapping;
this is reported as an annotated-domain sensitivity and is not substituted
into the historical result.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
ENRICHMENT_OUT = HERE / "outputs"
AUDIT_OUT = ENRICHMENT_OUT / "audit"
RAW_DIR = ENRICHMENT_OUT / "raw_gprofiler"
MULTI_DB = REPO_ROOT / "analysis" / "dinp_crc_multi_database_target_convergence"
QUERY_FILE = MULTI_DB / "outputs" / "dinp_crc_genecards_relevance10_expanded_intersection.csv"
EXPOSURE_FILE = MULTI_DB / "outputs" / "dinp_exposure_gene_matrix_ctd_toxcast_tox21_t3db.csv"
GENECARDS_FILE = MULTI_DB / "outputs" / "source_records" / "genecards_crc_top2000_reference.csv"
GMT_FILE = REPO_ROOT / "work" / "gprofiler_audit_cache" / "gprofiler_full_hsapiens.ENSG.gmt"
GMT_URL = "https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_hsapiens.ENSG.gmt"
CONVERT_URL = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
SOURCES = ("GO:BP", "REAC")
PREFIXES = {"GO:BP": "GO:", "KEGG": "KEGG:", "REAC": "REAC:"}
USER_AGENT = "whynot17-expanded-ora-audit/1.0"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    """Return P[X >= observed] for Hypergeometric(population, successes, draws)."""
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
    terms = [
        log_combination(successes, x)
        + log_combination(population - successes, draws - x)
        - denominator
        for x in range(lo, hi + 1)
    ]
    peak = max(terms)
    return min(1.0, math.exp(peak) * sum(math.exp(value - peak) for value in terms))


def bh(values: list[float]) -> list[float]:
    output = [float("nan")] * len(values)
    valid = [(i, value) for i, value in enumerate(values) if math.isfinite(value) and 0 <= value <= 1]
    if not valid:
        return output
    ordered = sorted(valid, key=lambda item: (item[1], item[0]))
    running = 1.0
    adjusted = [0.0] * len(ordered)
    for position in range(len(ordered) - 1, -1, -1):
        rank = position + 1
        running = min(running, ordered[position][1] * len(ordered) / rank)
        adjusted[position] = min(1.0, running)
    for (index, _), value in zip(ordered, adjusted):
        output[index] = value
    return output


def convert_symbols(symbols: list[str]) -> dict[str, Any]:
    payload = {"organism": "hsapiens", "target": "ENSG", "query": symbols}
    request = urllib.request.Request(
        CONVERT_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        raw = response.read()
    result = json.loads(raw.decode("utf-8"))
    by_symbol: dict[str, list[str]] = defaultdict(list)
    for row in result.get("result", []):
        converted = row.get("converted")
        if converted not in (None, "None", ""):
            by_symbol[str(row.get("incoming", "")).strip().upper()].append(str(converted))
    return {
        "payload": payload,
        "response": result,
        "by_symbol": {key: sorted(set(value)) for key, value in sorted(by_symbol.items())},
        "raw_sha256": hashlib.sha256(raw).hexdigest(),
    }


def ensg_set(mapping: dict[str, list[str]]) -> set[str]:
    return {ensg for values in mapping.values() for ensg in values}


def load_gmt(term_ids_by_source: dict[str, set[str]]) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, int]]:
    if not GMT_FILE.exists() or GMT_FILE.stat().st_size == 0:
        raise FileNotFoundError(f"Missing frozen GMT snapshot: {GMT_FILE}")
    term_genes: dict[str, set[str]] = {}
    source_union: dict[str, set[str]] = {source: set() for source in SOURCES}
    source_term_counts: dict[str, int] = {source: 0 for source in SOURCES}
    wanted = set().union(*(term_ids_by_source[source] for source in SOURCES))
    with GMT_FILE.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            term_id = fields[0]
            genes = {value.strip() for value in fields[2:] if value.strip()}
            source = next((name for name, prefix in PREFIXES.items() if term_id.startswith(prefix)), None)
            if source is None:
                continue
            source_term_counts[source] += 1
            source_union[source].update(genes)
            if term_id in wanted:
                term_genes[term_id] = genes
    return term_genes, source_union, source_term_counts


def parse_gene_cards_audit() -> dict[str, Any]:
    rows = read_csv(GENECARDS_FILE)
    ranked = [
        {
            "rank": int(row["GeneCards_Rank"]),
            "symbol": row["GeneSymbol"].strip().upper(),
            "relevance_score": float(row["RelevanceScore"]),
        }
        for row in rows
    ]
    by_rank = {row["rank"]: row for row in ranked}
    threshold = 10.0
    count_ge_threshold = sum(row["relevance_score"] >= threshold for row in ranked)
    monotonic = all(
        ranked[index]["relevance_score"] >= ranked[index + 1]["relevance_score"]
        for index in range(len(ranked) - 1)
    )
    rank_881 = by_rank.get(881)
    rank_882 = by_rank.get(882)
    rank_2000 = by_rank.get(2000)
    conditions = {
        "exactly_2000_rows": len(ranked) == 2000,
        "complete_rank_sequence_1_to_2000": sorted(by_rank) == list(range(1, 2001)),
        "unique_symbols": len({row["symbol"] for row in ranked}) == len(ranked),
        "score_monotonic_nonincreasing": monotonic,
        "count_relevance_ge_10_is_881": count_ge_threshold == 881,
        "rank_881_score_ge_10": bool(rank_881 and rank_881["relevance_score"] >= threshold),
        "rank_882_score_lt_10": bool(rank_882 and rank_882["relevance_score"] < threshold),
        "rank_2000_score_lt_10": bool(rank_2000 and rank_2000["relevance_score"] < threshold),
    }
    return {
        "source_file": repo_relative(GENECARDS_FILE),
        "source_sha256": sha256_file(GENECARDS_FILE),
        "source_scope": "archived ordinary GeneCards CRC top-2000 reference",
        "full_crc_genecards_export_present_in_this_source": False,
        "threshold_rule": "RelevanceScore >= 10.0",
        "rows": len(ranked),
        "unique_symbols": len({row["symbol"] for row in ranked}),
        "count_relevance_ge_10": count_ge_threshold,
        "boundary": {"rank_881": rank_881, "rank_882": rank_882, "rank_2000": rank_2000},
        "conditions": conditions,
        "audit_status": "PASS_WITH_ARCHIVED_TOP2000_SCOPE" if all(conditions.values()) else "FAIL",
        "interpretation": (
            "Within the archived ranked top-2000 source, the threshold crosses from 10.0 at rank 881 to 9.9 at rank 882, "
            "and the scores are nonincreasing through rank 2000 (score 3.1). This supports 881 as the complete >=10 set "
            "under the source's globally sorted-ranking assumption; it is not evidence that a full login-restricted export was obtained."
        ),
    }


def read_api_rows(background_name: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    raw_path = RAW_DIR / f"gprofiler_{background_name}.json"
    result = json.loads(raw_path.read_text(encoding="utf-8"))
    rows = [row for row in result.get("result", []) if row.get("source") in SOURCES]
    return result, rows


def read_historical_rows(background_name: str) -> dict[tuple[str, str], dict[str, str]]:
    path = ENRICHMENT_OUT / f"dinp_crc_enrichment_{background_name}.csv"
    return {(row["source"], row["term_id"]): row for row in read_csv(path)}


def number(value: Any) -> float:
    return float(value)


def main() -> None:
    AUDIT_OUT.mkdir(parents=True, exist_ok=True)
    query_rows = read_csv(QUERY_FILE)
    exposure_rows = read_csv(EXPOSURE_FILE)
    cards_rows = read_csv(GENECARDS_FILE)
    query_symbols = unique_symbols([row.get("gene_symbol", "") for row in query_rows])
    sensitivity_symbols = unique_symbols([row.get("gene_symbol", "") for row in exposure_rows])
    primary_symbols = unique_symbols(
        [
            row.get("gene_symbol") or row.get("GeneSymbol", "")
            for row in cards_rows
            if float(row.get("RelevanceScore", "nan")) >= 10.0
        ]
    )
    input_sets = {
        "query18": query_symbols,
        "sensitivity_dinp93": sensitivity_symbols,
        "primary_genecards881": primary_symbols,
    }
    if len(query_symbols) != 18 or len(sensitivity_symbols) != 93 or len(primary_symbols) != 881:
        raise ValueError(f"Unexpected frozen counts: { {key: len(value) for key, value in input_sets.items()} }")

    conversions: dict[str, dict[str, Any]] = {}
    for label, symbols in input_sets.items():
        conversions[label] = convert_symbols(symbols)
        (AUDIT_OUT / f"gprofiler_convert_{label}.json").write_text(
            json.dumps(
                {
                    "retrieved_at": utc_now(),
                    "endpoint": CONVERT_URL,
                    "target": "ENSG",
                    "organism": "hsapiens",
                    "input_count": len(symbols),
                    "conversion": conversions[label],
                },
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

    query_ensg = ensg_set(conversions["query18"]["by_symbol"])
    backgrounds = {
        "primary_genecards881": ensg_set(conversions["primary_genecards881"]["by_symbol"]),
        "sensitivity_dinp93": ensg_set(conversions["sensitivity_dinp93"]["by_symbol"]),
    }
    term_ids_by_source: dict[str, set[str]] = {source: set() for source in SOURCES}
    api_cache: dict[str, tuple[dict[str, Any], list[dict[str, Any]]]] = {}
    for background_name in backgrounds:
        api_cache[background_name] = read_api_rows(background_name)
        for row in api_cache[background_name][1]:
            term_ids_by_source[row["source"]].add(row["native"])
    term_genes, source_union, source_term_counts = load_gmt(term_ids_by_source)

    mapping_rows: list[dict[str, Any]] = []
    for label, conversion in conversions.items():
        for symbol, ensgs in conversion["by_symbol"].items():
            for ensg in ensgs:
                mapping_rows.append({"input_set": label, "input_symbol": symbol, "ensembl_id": ensg})
    with (AUDIT_OUT / "canonical_id_mapping.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["input_set", "input_symbol", "ensembl_id"])
        writer.writeheader()
        writer.writerows(mapping_rows)

    comparison_rows: list[dict[str, Any]] = []
    background_summaries: dict[str, Any] = {}
    for background_name, background_ensg in backgrounds.items():
        _, api_rows = api_cache[background_name]
        historical = read_historical_rows(background_name)
        by_source: dict[str, list[dict[str, Any]]] = {source: [] for source in SOURCES}
        for api_row in api_rows:
            by_source[api_row["source"]].append(api_row)
        source_summaries: dict[str, Any] = {}
        for source in SOURCES:
            rows = by_source[source]
            custom_raw: list[float] = []
            annotated_raw: list[float] = []
            source_background_annotated = background_ensg & source_union[source]
            source_query_annotated = query_ensg & source_union[source]
            metric_match = True
            historical_raw_match = True
            term_map_complete = True
            source_comparisons: list[dict[str, Any]] = []
            for api_row in rows:
                term_id = api_row["native"]
                genes = term_genes.get(term_id)
                if genes is None:
                    term_map_complete = False
                    continue
                n_custom = len(query_ensg)
                k_custom = len(query_ensg & genes)
                K_custom = len(background_ensg & genes)
                n_annotated = len(source_query_annotated)
                k_annotated = len(source_query_annotated & genes)
                K_annotated = len(source_background_annotated & genes)
                local_custom = hypergeom_sf(len(background_ensg), K_custom, n_custom, k_custom)
                local_annotated = hypergeom_sf(len(source_background_annotated), K_annotated, n_annotated, k_annotated)
                custom_raw.append(local_custom)
                annotated_raw.append(local_annotated)
                api_metric = (
                    int(api_row.get("term_size", 0)),
                    int(api_row.get("intersection_size", 0)),
                    int(api_row.get("query_size", 0)),
                    int(api_row.get("effective_domain_size", 0)),
                )
                local_metric = (K_custom, k_custom, n_custom, len(background_ensg))
                if api_metric != local_metric:
                    metric_match = False
                historical_row = historical.get((source, term_id))
                historical_value = number(historical_row["raw_hypergeom_p"]) if historical_row else float("nan")
                if historical_row is None or not math.isclose(historical_value, local_custom, rel_tol=1e-10, abs_tol=1e-12):
                    historical_raw_match = False
                source_comparisons.append(
                    {
                        "background": background_name,
                        "source": source,
                        "term_id": term_id,
                        "term_name": api_row.get("name", ""),
                        "term_size_api": int(api_row.get("term_size", 0)),
                        "term_size_local": K_custom,
                        "intersection_size_api": int(api_row.get("intersection_size", 0)),
                        "intersection_size_local": k_custom,
                        "query_size_api": int(api_row.get("query_size", 0)),
                        "query_size_local": n_custom,
                        "effective_domain_api": int(api_row.get("effective_domain_size", 0)),
                        "effective_domain_local": len(background_ensg),
                        "local_raw_p_custom_domain": local_custom,
                        "local_raw_p_annotated_domain": local_annotated,
                        "gprofiler_p_value_adjusted": number(api_row.get("p_value")),
                        "historical_runner_raw_p": historical_value,
                    }
                )
            reconstructed = bh(custom_raw)
            annotated_bh = bh(annotated_raw)
            for item, reconstructed_p, annotated_p in zip(source_comparisons, reconstructed, annotated_bh):
                item["gprofiler_p_value_reconstructed_bh"] = reconstructed_p
                item["local_annotated_domain_bh_fdr"] = annotated_p
                item["api_adjusted_p_abs_diff"] = abs(item["gprofiler_p_value_adjusted"] - reconstructed_p)
                item["audit_focus"] = item["gprofiler_p_value_adjusted"] < 0.05 or item["gprofiler_p_value_reconstructed_bh"] < 0.05
            comparison_rows.extend(source_comparisons)
            api_diffs = [item["api_adjusted_p_abs_diff"] for item in source_comparisons]
            source_summaries[source] = {
                "api_rows": len(rows),
                "fixed_gmt_term_map_complete": term_map_complete and len(source_comparisons) == len(rows),
                "api_term_metrics_match_fixed_gmt": metric_match,
                "historical_runner_raw_p_matches_fixed_gmt_local_p": historical_raw_match,
                "gprofiler_adjusted_p_reproduced_by_sourcewise_bh": bool(api_diffs) and max(api_diffs) <= 1e-10,
                "api_adjusted_p_max_abs_difference": max(api_diffs) if api_diffs else None,
                "mapped_background_N_custom_domain": len(background_ensg),
                "annotated_background_N_alternative_domain": len(source_background_annotated),
                "mapped_query_n_custom_domain": len(query_ensg),
                "annotated_query_n_alternative_domain": len(source_query_annotated),
                "gprofiler_effective_domain_sizes": sorted({int(row.get("effective_domain_size", 0)) for row in rows}),
                "gprofiler_adjusted_fdr_lt_0_05": sum(item["gprofiler_p_value_adjusted"] < 0.05 for item in source_comparisons),
                "local_annotated_domain_bh_fdr_lt_0_05": sum(item["local_annotated_domain_bh_fdr"] < 0.05 for item in source_comparisons),
            }
        background_summaries[background_name] = {
            "input_symbol_count": len(input_sets[background_name]),
            "canonical_ensembl_count": len(background_ensg),
            "one_to_many_symbols": {
                symbol: values
                for symbol, values in conversions[background_name]["by_symbol"].items()
                if len(values) > 1
            },
            "unmapped_symbols": sorted(set(input_sets[background_name]) - set(conversions[background_name]["by_symbol"])),
            "source_summaries": source_summaries,
        }

    comparison_fields = [
        "background", "source", "term_id", "term_name", "term_size_api", "term_size_local",
        "intersection_size_api", "intersection_size_local", "query_size_api", "query_size_local",
        "effective_domain_api", "effective_domain_local", "local_raw_p_custom_domain",
        "local_raw_p_annotated_domain", "gprofiler_p_value_adjusted", "historical_runner_raw_p",
        "gprofiler_p_value_reconstructed_bh", "local_annotated_domain_bh_fdr", "api_adjusted_p_abs_diff",
        "audit_focus",
    ]
    with (AUDIT_OUT / "local_vs_gprofiler_ora.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=comparison_fields)
        writer.writeheader()
        writer.writerows(comparison_rows)

    conversion_summary = {}
    for label, conversion in conversions.items():
        symbols = input_sets[label]
        conversion_summary[label] = {
            "input_symbols": len(symbols),
            "conversion_rows": len(conversion["response"].get("result", [])),
            "canonical_ensembl_ids": len(ensg_set(conversion["by_symbol"])),
            "unmapped_symbols": sorted(set(symbols) - set(conversion["by_symbol"])),
            "one_to_many_symbols": {
                symbol: values for symbol, values in conversion["by_symbol"].items() if len(values) > 1
            },
            "response_sha256": conversion["raw_sha256"],
        }
    all_source_checks = [
        check
        for summary in background_summaries.values()
        for check in summary["source_summaries"].values()
        for check in [
            check["fixed_gmt_term_map_complete"],
            check["api_term_metrics_match_fixed_gmt"],
            check["historical_runner_raw_p_matches_fixed_gmt_local_p"],
            check["gprofiler_adjusted_p_reproduced_by_sourcewise_bh"],
        ]
    ]
    overall = all(all_source_checks) and conversion_summary["query18"]["canonical_ensembl_ids"] == 18 and conversion_summary["sensitivity_dinp93"]["canonical_ensembl_ids"] == 100
    audit = {
        "audit": "expanded_dinp_crc_background_ora_and_genecards_audit",
        "run_utc": utc_now(),
        "inputs": {
            "query": {"path": repo_relative(QUERY_FILE), "sha256": sha256_file(QUERY_FILE), "count": len(query_symbols)},
            "sensitivity_background": {"path": repo_relative(EXPOSURE_FILE), "sha256": sha256_file(EXPOSURE_FILE), "count": len(sensitivity_symbols)},
            "primary_background": {"path": repo_relative(GENECARDS_FILE), "sha256": sha256_file(GENECARDS_FILE), "count": len(primary_symbols)},
        },
        "canonical_id_namespace": "Ensembl gene IDs via g:Profiler g:Convert, target=ENSG, organism=hsapiens",
        "gprofiler_version": json.loads((RAW_DIR / "gprofiler_sensitivity_dinp93.json").read_text(encoding="utf-8"))["meta"].get("version"),
        "fixed_gene_set_mapping": {
            "url": GMT_URL,
            "path": repo_relative(GMT_FILE),
            "sha256": sha256_file(GMT_FILE),
            "bytes": GMT_FILE.stat().st_size,
            "source_term_counts_all_gmt": source_term_counts,
            "terms_rechecked_from_api": {source: len(term_ids_by_source[source]) for source in SOURCES},
        },
        "frozen_hypergeometric_definition": {
            "formula": "P(X>=k)=sum_{i=k}^{min(K,n)} choose(K,i)*choose(N-K,n-i)/choose(N,n)",
            "custom_domain_N": "all unique canonical Ensembl IDs mapped from the frozen input background; matches domain_scope=custom",
            "annotated_domain_N": "background IDs occurring in the relevant fixed GO/Reactome mapping; reported as an alternative sensitivity only",
            "bh": "Benjamini-Hochberg applied source-wise to the returned terms with nonzero query intersection, reproducing g:Profiler p_value",
        },
        "conversion_summary": conversion_summary,
        "background_summaries": background_summaries,
        "genecards_threshold_audit": parse_gene_cards_audit(),
        "overall_status": "PASS" if overall and parse_gene_cards_audit()["audit_status"] == "PASS_WITH_ARCHIVED_TOP2000_SCOPE" else "REVIEW",
        "interpretation_boundary": (
            "The 93-symbol sensitivity background maps to 100 canonical Ensembl IDs because SNORA72 expands to eight IDs. "
            "This is identifier expansion, not an additional biological input gene list. The historical g:Profiler run used N=100 "
            "under domain_scope=custom. The annotated-domain N values are not substituted into historical results."
        ),
    }
    (AUDIT_OUT / "step_gene_background_ora_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (AUDIT_OUT / "genecards_881_threshold_audit.json").write_text(
        json.dumps(audit["genecards_threshold_audit"], indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# Expanded DINP--CRC ORA and GeneCards background audit",
        "",
        f"Audit status: **{audit['overall_status']}**",
        "",
        "## Executive findings",
        "",
        "- The sensitivity input is 93 gene symbols, but g:Profiler maps it to 100 unique Ensembl IDs.",
        "- The complete one-to-many expansion is `SNORA72` to 8 Ensembl IDs; the other sensitivity symbols map one-to-one and none fail.",
        "- Fixed GMT term membership, N/K/n/k metrics, historical local raw P values, and g:Profiler adjusted P values were checked source-by-source.",
        "- g:Profiler `p_value` is a multiple-testing-adjusted value; the local raw hypergeometric P is compared to a source-wise BH reconstruction rather than to that adjusted field directly.",
        "",
        "## Domain accounting",
        "",
        "| Background | Input symbols | Canonical Ensembl IDs | GO:BP annotated N | Reactome annotated N |",
        "|---|---:|---:|---:|---:|",
    ]
    for background_name, summary in background_summaries.items():
        source_summaries = summary["source_summaries"]
        lines.append(
            f"| {background_name} | {summary['input_symbol_count']} | {summary['canonical_ensembl_count']} | "
            f"{source_summaries['GO:BP']['annotated_background_N_alternative_domain']} | "
            f"{source_summaries['REAC']['annotated_background_N_alternative_domain']} |"
        )
    lines += [
        "",
        "## Source-wise recheck",
        "",
        "| Background | Source | API rows | API N | Fixed GMT metrics | Historical raw P | API adjusted P reproduced |",
        "|---|---|---:|---:|---|---|---|",
    ]
    for background_name, summary in background_summaries.items():
        for source in SOURCES:
            item = summary["source_summaries"][source]
            lines.append(
                f"| {background_name} | {source} | {item['api_rows']} | {','.join(map(str, item['gprofiler_effective_domain_sizes']))} | "
                f"{'PASS' if item['api_term_metrics_match_fixed_gmt'] else 'FAIL'} | "
                f"{'PASS' if item['historical_runner_raw_p_matches_fixed_gmt_local_p'] else 'FAIL'} | "
                f"{'PASS' if item['gprofiler_adjusted_p_reproduced_by_sourcewise_bh'] else 'FAIL'} |"
            )
    gc = audit["genecards_threshold_audit"]
    lines += [
        "",
        "## GeneCards RelevanceScore >= 10 audit",
        "",
        f"- Source: `{gc['source_scope']}`",
        f"- Rows/ranks: {gc['rows']} / ranks 1--2000; unique symbols: {gc['unique_symbols']}",
        f"- Score >=10: **{gc['count_relevance_ge_10']}**",
        f"- Rank 881: `{gc['boundary']['rank_881']['symbol']}` score {gc['boundary']['rank_881']['relevance_score']}",
        f"- Rank 882: `{gc['boundary']['rank_882']['symbol']}` score {gc['boundary']['rank_882']['relevance_score']}",
        f"- Rank 2000: `{gc['boundary']['rank_2000']['symbol']}` score {gc['boundary']['rank_2000']['relevance_score']}",
        f"- Within-file score monotonicity: **{'PASS' if gc['conditions']['score_monotonic_nonincreasing'] else 'FAIL'}**",
        "",
        gc["interpretation"],
        "",
        "## Files",
        "",
        "- `local_vs_gprofiler_ora.csv`: one row per API-returned GO:BP/Reactome term, with custom-domain and annotated-domain local calculations.",
        "- `canonical_id_mapping.csv`: symbol-to-Ensembl mappings used in the local calculation.",
        "- `step_gene_background_ora_audit.json`: complete machine-readable audit.",
        "- `genecards_881_threshold_audit.json`: GeneCards threshold proof.",
        "",
        "## Boundary",
        "",
        "The audit validates the historical `domain_scope=custom` implementation and separately reports an annotated-domain alternative. The official static g:Profiler GMT snapshot used here does not include KEGG; KEGG is therefore outside this fixed-mapping audit and its historical output is untouched. The audit does not silently rewrite the canonical enrichment results or claim that the archived top-2000 file is a full GeneCards export.",
    ]
    (AUDIT_OUT / "STEP_GENE_BACKGROUND_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("ORA + GENECARDS AUDIT: " + audit["overall_status"])
    print(json.dumps({
        "sensitivity_input_symbols": len(sensitivity_symbols),
        "sensitivity_canonical_ensg": len(backgrounds["sensitivity_dinp93"]),
        "sensitivity_one_to_many": background_summaries["sensitivity_dinp93"]["one_to_many_symbols"],
        "genecards_ge10": audit["genecards_threshold_audit"]["count_relevance_ge_10"],
        "audit_dir": str(AUDIT_OUT.resolve()),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
