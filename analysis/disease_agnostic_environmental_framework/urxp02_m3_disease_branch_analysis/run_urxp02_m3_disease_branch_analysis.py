"""URXP02 M3: disease-branch enrichment and PPI/module analysis.

This script starts from the frozen M1b expanded molecular universe and M2
disease-branch classification.  It performs no NHANES modelling and makes no
sex-specific molecular claim.  The three disease branches are analysed with
custom-background ORA and STRING functional PPI networks; all evidence-layer
fields from M1b are carried into the results.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import time
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
CLASSIFICATION = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_m2_disease_branch" / "03_branch_gene_classification.csv"
M2_SUMMARY = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_m2_disease_branch" / "04_branch_evidence_summary.csv"
GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
STRING_BASE = "https://string-db.org/api/json"
SPECIES = 9606
STRING_REQUIRED_SCORE = 700
STRING_NETWORK_TYPE = "functional"
ENRICHMENT_SOURCES = ["GO:BP", "REAC", "KEGG"]
BRANCHES = ["thyroid-specific", "hypertension-specific", "shared"]
NETWORK_BRANCHES = ["thyroid-branch", "hypertension-branch", "shared-core"]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: str | bool | None) -> bool:
    return str(value or "").strip().lower() == "true"


def as_int(value: str | int | float | None) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def http_json(url: str, *, method: str = "GET", body: dict | None = None, timeout: int = 180) -> dict | list:
    encoded = None if body is None else json.dumps(body, ensure_ascii=False).encode("utf-8")
    headers = {"User-Agent": "URXP02-M3/1.0", "Accept": "application/json"}
    if body is not None:
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=encoded, headers=headers, method=method)
    last_error = None
    for attempt in range(2):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:  # retry once for transient public-API failures
            last_error = exc
            if attempt == 0:
                time.sleep(2)
    raise RuntimeError(f"API request failed: {url}: {last_error}")


def branch_evidence(rows: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {}
    for branch in BRANCHES:
        subset = [r for r in rows if r["branch_class"] == branch]
        out[branch] = {
            "n_genes": len(subset),
            "exact_2NAP_human": sum(as_bool(r.get("exact_2NAP_human_support")) for r in subset),
            "exact_2NAP_experimental": sum(as_bool(r.get("exact_2NAP_experimental_support")) for r in subset),
            "parent_naphthalene": sum(as_bool(r.get("parent_naphthalene_support")) for r in subset),
            "multi_source_ge2": sum(as_int(r.get("number_of_sources")) >= 2 for r in subset),
            "CTD": sum(as_bool(r.get("CTD_support")) for r in subset),
            "toxicogenomic": sum(as_bool(r.get("toxicogenomic_support")) for r in subset),
            "bioassay_target": sum(as_bool(r.get("bioassay_target_support")) for r in subset),
        }
    return out


def evidence_for_genes(rows: list[dict[str, str]], genes: list[str]) -> dict[str, int]:
    """Summarize M1b evidence for an arbitrary network gene set."""
    wanted = set(genes)
    subset = [r for r in rows if r.get("gene_symbol") in wanted]
    return {
        "n_genes": len(subset),
        "exact_2NAP_human": sum(as_bool(r.get("exact_2NAP_human_support")) for r in subset),
        "exact_2NAP_experimental": sum(as_bool(r.get("exact_2NAP_experimental_support")) for r in subset),
        "parent_naphthalene": sum(as_bool(r.get("parent_naphthalene_support")) for r in subset),
        "multi_source_ge2": sum(as_int(r.get("number_of_sources")) >= 2 for r in subset),
        "CTD": sum(as_bool(r.get("CTD_support")) for r in subset),
        "toxicogenomic": sum(as_bool(r.get("toxicogenomic_support")) for r in subset),
        "bioassay_target": sum(as_bool(r.get("bioassay_target_support")) for r in subset),
    }


def enrichment(branch: str, genes: list[str], background: list[str], evidence: dict[str, int]) -> tuple[list[dict], dict]:
    accessed = now_utc()
    payload = {
        "organism": "hsapiens",
        "query": genes,
        "background": background,
        "domain_scope": "custom",
        "sources": ENRICHMENT_SOURCES,
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "all_results": True,
        "no_evidences": True,
    }
    try:
        response = http_json(GPROFILER_URL, method="POST", body=payload)
        if not isinstance(response, dict):
            raise RuntimeError("g:Profiler returned a non-object response")
        if response.get("error"):
            raise RuntimeError(str(response["error"]))
        records = response.get("result") or []
        rows = []
        for rec in records:
            rows.append({
                "branch_class": branch,
                "status": "OK",
                "source": rec.get("source", ""),
                "native_id": rec.get("native", ""),
                "term_name": rec.get("name", ""),
                "description": rec.get("description", ""),
                "intersection_size": rec.get("intersection_size", ""),
                "term_size": rec.get("term_size", ""),
                "effective_domain_size": rec.get("effective_domain_size", ""),
                "query_size_mappable": rec.get("query_size", ""),
                "p_value_corrected": rec.get("p_value", ""),
                "significant_g_scs": rec.get("significant", False),
                "precision": rec.get("precision", ""),
                "recall": rec.get("recall", ""),
                "source_order": rec.get("source_order", ""),
                "group_id": rec.get("group_id", ""),
                "branch_n_genes": evidence["n_genes"],
                "branch_exact_2NAP_human_n": evidence["exact_2NAP_human"],
                "branch_exact_2NAP_experimental_n": evidence["exact_2NAP_experimental"],
                "branch_parent_naphthalene_n": evidence["parent_naphthalene"],
                "branch_multi_source_ge2_n": evidence["multi_source_ge2"],
                "branch_CTD_n": evidence["CTD"],
                "branch_toxicogenomic_n": evidence["toxicogenomic"],
                "branch_bioassay_target_n": evidence["bioassay_target"],
                "background_n_genes": len(background),
                "background_definition": "all 828 M1b expanded-universe gene symbols",
                "multiple_testing": "g:SCS",
                "accessed_at_utc": accessed,
                "source_database": "g:Profiler g:GOSt",
            })
        if not rows:
            rows = [{
                "branch_class": branch,
                "status": "NO_TERMS_RETURNED",
                "branch_n_genes": evidence["n_genes"],
                "background_n_genes": len(background),
                "background_definition": "all 828 M1b expanded-universe gene symbols",
                "multiple_testing": "g:SCS",
                "accessed_at_utc": accessed,
                "source_database": "g:Profiler g:GOSt",
            }]
        audit = {
            "branch_class": branch,
            "status": "OK",
            "query_n_genes": len(genes),
            "retrieved_n_terms": len(records),
            "significant_n_terms": sum(bool(r.get("significant")) for r in records),
            "sources": ";".join(ENRICHMENT_SOURCES),
            "background_n_genes": len(background),
            "background_definition": "all 828 M1b expanded-universe gene symbols",
            "domain_scope": "custom",
            "multiple_testing": "g:SCS",
            "user_threshold": 0.05,
            "accessed_at_utc": accessed,
            "endpoint": GPROFILER_URL,
        }
        return rows, audit
    except Exception as exc:
        return [{
            "branch_class": branch,
            "status": "API_ERROR",
            "error": str(exc),
            "branch_n_genes": evidence["n_genes"],
            "background_n_genes": len(background),
            "background_definition": "all 828 M1b expanded-universe gene symbols",
            "multiple_testing": "g:SCS",
            "accessed_at_utc": accessed,
            "source_database": "g:Profiler g:GOSt",
        }], {
            "branch_class": branch,
            "status": "API_ERROR",
            "query_n_genes": len(genes),
            "retrieved_n_terms": 0,
            "significant_n_terms": 0,
            "sources": ";".join(ENRICHMENT_SOURCES),
            "background_n_genes": len(background),
            "background_definition": "all 828 M1b expanded-universe gene symbols",
            "domain_scope": "custom",
            "multiple_testing": "g:SCS",
            "user_threshold": 0.05,
            "accessed_at_utc": accessed,
            "endpoint": GPROFILER_URL,
            "error": str(exc),
        }


def string_call(endpoint: str, genes: list[str], **extra) -> list[dict] | dict:
    params = {"identifiers": "\r".join(genes), "species": SPECIES, **extra}
    url = f"{STRING_BASE}/{endpoint}?{urllib.parse.urlencode(params)}"
    return http_json(url, timeout=180)


def network_analysis(branch: str, genes: list[str], evidence: dict[str, int]) -> tuple[list[dict], list[dict], list[dict], dict]:
    accessed = now_utc()
    mapping = string_call("get_string_ids", genes)
    mapping_rows = mapping if isinstance(mapping, list) else []
    query_to_pref = {str(r.get("queryItem")): str(r.get("preferredName")) for r in mapping_rows if r.get("queryItem") and r.get("preferredName")}
    query_to_id = {str(r.get("queryItem")): str(r.get("stringId")) for r in mapping_rows if r.get("queryItem") and r.get("stringId")}
    pref_to_query = {pref: query for query, pref in query_to_pref.items()}

    graph = nx.Graph()
    graph.add_nodes_from(genes)
    edge_response = string_call(
        "network",
        genes,
        required_score=STRING_REQUIRED_SCORE,
        network_type=STRING_NETWORK_TYPE,
        add_nodes=0,
    )
    edge_rows = []
    for rec in edge_response if isinstance(edge_response, list) else []:
        pref_a = str(rec.get("preferredName_A") or "")
        pref_b = str(rec.get("preferredName_B") or "")
        gene_a = pref_to_query.get(pref_a, pref_a if pref_a in genes else "")
        gene_b = pref_to_query.get(pref_b, pref_b if pref_b in genes else "")
        if not gene_a or not gene_b or gene_a == gene_b:
            continue
        score = float(rec.get("score") or 0)
        graph.add_edge(gene_a, gene_b, score=score)
        edge_rows.append({
            "branch_class": branch,
            "gene_a": gene_a,
            "gene_b": gene_b,
            "string_id_a": rec.get("stringId_A", ""),
            "string_id_b": rec.get("stringId_B", ""),
            "preferred_name_a": pref_a,
            "preferred_name_b": pref_b,
            "score": score,
            "required_score": STRING_REQUIRED_SCORE,
            "network_type": STRING_NETWORK_TYPE,
            "species": SPECIES,
            "accessed_at_utc": accessed,
            "source_database": "STRING",
        })

    ppi_response = string_call(
        "ppi_enrichment",
        genes,
        required_score=STRING_REQUIRED_SCORE,
        network_type=STRING_NETWORK_TYPE,
    )
    ppi = (ppi_response[0] if isinstance(ppi_response, list) and ppi_response else {}) if isinstance(ppi_response, (list, dict)) else {}
    components = list(nx.connected_components(graph))
    largest = max((len(c) for c in components), default=0)
    n = graph.number_of_nodes()
    e = graph.number_of_edges()
    density = nx.density(graph) if n > 1 else 0.0
    clustering = nx.average_clustering(graph, weight="score") if n else 0.0
    try:
        betweenness = nx.betweenness_centrality(graph, normalized=True)
    except Exception:
        betweenness = {node: 0.0 for node in graph.nodes}
    try:
        pagerank = nx.pagerank(graph, weight="score", max_iter=500)
    except Exception:
        pagerank = {node: 0.0 for node in graph.nodes}
    closeness = nx.closeness_centrality(graph)
    degree = dict(graph.degree())
    weighted_degree = {node: sum(float(data.get("score", 0.0)) for _, _, data in graph.edges(node, data=True)) for node in graph.nodes}
    core = nx.core_number(graph) if graph.number_of_nodes() else {}

    def rank_metric(values: dict[str, float], reverse: bool = True) -> dict[str, int]:
        ordered = sorted(values, key=lambda x: (-values[x] if reverse else values[x], x))
        return {node: index + 1 for index, node in enumerate(ordered)}

    degree_rank = rank_metric(degree)
    betweenness_rank = rank_metric(betweenness)
    pagerank_rank = rank_metric(pagerank)
    top_n = max(1, math.ceil(max(1, n) * 0.10))
    node_rows = []
    for gene in sorted(graph.nodes):
        node_rows.append({
            "branch_class": branch,
            "gene_symbol": gene,
            "string_id": query_to_id.get(gene, ""),
            "string_preferred_name": query_to_pref.get(gene, gene),
            "mapped_to_STRING": gene in query_to_id,
            "degree": degree.get(gene, 0),
            "weighted_degree": weighted_degree.get(gene, 0.0),
            "degree_centrality": (degree.get(gene, 0) / (n - 1)) if n > 1 else 0.0,
            "degree_rank": degree_rank.get(gene, ""),
            "betweenness_centrality": betweenness.get(gene, 0.0),
            "betweenness_rank": betweenness_rank.get(gene, ""),
            "closeness_centrality": closeness.get(gene, 0.0),
            "pagerank": pagerank.get(gene, 0.0),
            "pagerank_rank": pagerank_rank.get(gene, ""),
            "core_number": core.get(gene, 0),
            "network_central_top10pct_any": (degree_rank.get(gene, n + 1) <= top_n or betweenness_rank.get(gene, n + 1) <= top_n or pagerank_rank.get(gene, n + 1) <= top_n),
            "network_top10pct_n": top_n,
            "network_nodes_total": n,
            "accessed_at_utc": accessed,
            "source_database": "STRING + NetworkX",
        })

    communities = list(nx.algorithms.community.greedy_modularity_communities(graph, weight="score")) if n else []
    communities = sorted((sorted(c) for c in communities), key=lambda c: (-len(c), c[0] if c else ""))
    module_rows = []
    for module_index, members in enumerate(communities, start=1):
        sub = graph.subgraph(members)
        internal_edges = sub.number_of_edges()
        module_density = nx.density(sub) if len(members) > 1 else 0.0
        module_genes = ";".join(members)
        for gene in members:
            module_rows.append({
                "branch_class": branch,
                "module_id": f"{branch[:3].upper()}_M{module_index:03d}",
                "gene_symbol": gene,
                "module_size": len(members),
                "module_internal_edges": internal_edges,
                "module_density": module_density,
                "module_genes": module_genes,
                "algorithm": "greedy_modularity_communities",
                "edge_weight": "STRING combined score",
                "accessed_at_utc": accessed,
                "source_database": "STRING + NetworkX",
            })

    metrics = {
        "branch_class": branch,
        "input_n_genes": len(genes),
        "string_mapped_n_genes": len(query_to_id),
        "network_n_nodes": n,
        "network_n_edges": e,
        "network_n_components": len(components),
        "largest_component_n": largest,
        "network_density": density,
        "average_node_degree": (2 * e / n) if n else 0.0,
        "average_clustering_weighted": clustering,
        "module_n": len(communities),
        "module_largest_n": max((len(c) for c in communities), default=0),
        "ppi_enrichment_n_nodes": ppi.get("number_of_nodes", ""),
        "ppi_enrichment_n_edges": ppi.get("number_of_edges", ""),
        "ppi_enrichment_average_node_degree": ppi.get("average_node_degree", ""),
        "ppi_enrichment_expected_edges": ppi.get("expected_number_of_edges", ""),
        "ppi_enrichment_p_value": ppi.get("p_value", ""),
        "branch_exact_2NAP_human_n": evidence["exact_2NAP_human"],
        "branch_exact_2NAP_experimental_n": evidence["exact_2NAP_experimental"],
        "branch_parent_naphthalene_n": evidence["parent_naphthalene"],
        "branch_multi_source_ge2_n": evidence["multi_source_ge2"],
        "required_score": STRING_REQUIRED_SCORE,
        "network_type": STRING_NETWORK_TYPE,
        "species": SPECIES,
        "accessed_at_utc": accessed,
        "source_database": "STRING",
    }
    return edge_rows, node_rows, module_rows, metrics


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    classification = read_csv(CLASSIFICATION)
    _m2_summary = read_csv(M2_SUMMARY)  # read to keep the M2 evidence audit as an explicit input
    if len(classification) != 828:
        raise RuntimeError(f"Expected 828 M1b genes, found {len(classification)}")
    background = sorted({r["gene_symbol"] for r in classification if r.get("gene_symbol")})
    by_branch = {branch: sorted(r["gene_symbol"] for r in classification if r["branch_class"] == branch) for branch in BRANCHES}
    evidence = branch_evidence(classification)
    network_gene_sets = {
        "thyroid-branch": sorted(r["gene_symbol"] for r in classification if as_bool(r.get("in_thyroid_disease_set"))),
        "hypertension-branch": sorted(r["gene_symbol"] for r in classification if as_bool(r.get("in_hypertension_disease_set"))),
        "shared-core": by_branch["shared"],
    }
    network_evidence = {label: evidence_for_genes(classification, genes) for label, genes in network_gene_sets.items()}

    enrichment_rows: list[dict] = []
    enrichment_audit: list[dict] = []
    for branch in BRANCHES:
        rows, audit = enrichment(branch, by_branch[branch], background, evidence[branch])
        enrichment_rows.extend(rows)
        enrichment_audit.append(audit)
        time.sleep(1)

    enrichment_fields = [
        "branch_class", "status", "source", "native_id", "term_name", "description", "intersection_size", "term_size", "effective_domain_size", "query_size_mappable", "p_value_corrected", "significant_g_scs", "precision", "recall", "source_order", "group_id", "branch_n_genes", "branch_exact_2NAP_human_n", "branch_exact_2NAP_experimental_n", "branch_parent_naphthalene_n", "branch_multi_source_ge2_n", "branch_CTD_n", "branch_toxicogenomic_n", "branch_bioassay_target_n", "background_n_genes", "background_definition", "multiple_testing", "accessed_at_utc", "source_database", "error",
    ]
    write_csv(OUT / "01_enrichment_results.csv", enrichment_rows, enrichment_fields)
    write_csv(OUT / "02_enrichment_audit.csv", enrichment_audit, list(enrichment_audit[0].keys()))

    edge_rows: list[dict] = []
    node_rows: list[dict] = []
    module_rows: list[dict] = []
    network_metrics: list[dict] = []
    for branch in NETWORK_BRANCHES:
        edges, nodes, modules, metrics = network_analysis(branch, network_gene_sets[branch], network_evidence[branch])
        edge_rows.extend(edges)
        node_rows.extend(nodes)
        module_rows.extend(modules)
        network_metrics.append(metrics)
        time.sleep(1)

    write_csv(OUT / "03_ppi_edges.csv", edge_rows, [
        "branch_class", "gene_a", "gene_b", "string_id_a", "string_id_b", "preferred_name_a", "preferred_name_b", "score", "required_score", "network_type", "species", "accessed_at_utc", "source_database",
    ])
    write_csv(OUT / "04_ppi_node_centrality.csv", node_rows, [
        "branch_class", "gene_symbol", "string_id", "string_preferred_name", "mapped_to_STRING", "degree", "weighted_degree", "degree_centrality", "degree_rank", "betweenness_centrality", "betweenness_rank", "closeness_centrality", "pagerank", "pagerank_rank", "core_number", "network_central_top10pct_any", "network_top10pct_n", "network_nodes_total", "accessed_at_utc", "source_database",
    ])
    write_csv(OUT / "05_ppi_modules.csv", module_rows, [
        "branch_class", "module_id", "gene_symbol", "module_size", "module_internal_edges", "module_density", "module_genes", "algorithm", "edge_weight", "accessed_at_utc", "source_database",
    ])
    write_csv(OUT / "06_ppi_network_metrics.csv", network_metrics, list(network_metrics[0].keys()))

    evidence_by_symbol = {r["gene_symbol"]: r for r in classification}
    shared_centrality = [r for r in node_rows if r["branch_class"] == "shared-core"]
    hub_rows = []
    for node in sorted(shared_centrality, key=lambda r: (not r["network_central_top10pct_any"], int(r["degree_rank"] or 10**9), r["gene_symbol"]),):
        source = evidence_by_symbol[node["gene_symbol"]]
        exact_any = as_bool(source.get("exact_2NAP_human_support")) or as_bool(source.get("exact_2NAP_experimental_support"))
        multi_source = as_int(source.get("number_of_sources")) >= 2
        central = as_bool(node.get("network_central_top10pct_any"))
        hub_rows.append({
            "gene_symbol": node["gene_symbol"],
            "exact_2NAP_human_support": source.get("exact_2NAP_human_support", ""),
            "exact_2NAP_experimental_support": source.get("exact_2NAP_experimental_support", ""),
            "exact_2NAP_any": exact_any,
            "number_of_sources": source.get("number_of_sources", ""),
            "multi_source_ge2": multi_source,
            "parent_naphthalene_support": source.get("parent_naphthalene_support", ""),
            "CTD_support": source.get("CTD_support", ""),
            "toxicogenomic_support": source.get("toxicogenomic_support", ""),
            "bioassay_target_support": source.get("bioassay_target_support", ""),
            "degree": node.get("degree", ""),
            "degree_rank": node.get("degree_rank", ""),
            "betweenness_centrality": node.get("betweenness_centrality", ""),
            "betweenness_rank": node.get("betweenness_rank", ""),
            "pagerank": node.get("pagerank", ""),
            "pagerank_rank": node.get("pagerank_rank", ""),
            "network_central_top10pct_any": central,
            "priority_hub_candidate": exact_any and multi_source and central,
            "candidate_rule": "exact 2-NAP evidence (human or experimental) AND number_of_sources >=2 AND top 10% by degree, betweenness, or PageRank in shared STRING network",
        })
    write_csv(OUT / "07_shared_hub_candidates.csv", hub_rows, list(hub_rows[0].keys()) if hub_rows else ["gene_symbol"])

    generated_files = sorted(OUT.glob("*.csv"))
    report_lines = [
        "# URXP02 M3 disease-branch analysis",
        "",
        f"Generated {now_utc()}. This package analyzes the frozen M1b 828-gene universe against the M2 thyroid and hypertension disease branches.",
        "",
        "## Scope and safeguards",
        "",
        "- Custom-background ORA was run separately for thyroid-specific (30), hypertension-specific (251), and shared (189) genes.",
        "- Enrichment sources were GO Biological Process, Reactome, and KEGG using g:Profiler g:GOSt with the complete 828-gene expanded universe as background and g:SCS correction.",
        "- STRING functional networks were retrieved for the full 219-gene thyroid branch, 440-gene hypertension branch, and 189-gene shared core at combined score >=700; no added interactors were requested.",
        "- Modules use weighted greedy modularity communities; centrality is descriptive and calculated on the retrieved network.",
        "- Shared hub candidates require exact 2-NAP evidence (human or experimental), at least two evidence sources, and top-10% network centrality by any of degree, betweenness, or PageRank. This is a transparent gate, not a composite score.",
        "- No NHANES model, pathway-to-disease causal claim, sex-specific molecular claim, tissue/cell mapping, or figures were produced.",
        "",
        "## Branch evidence composition",
        "",
    ]
    for branch in BRANCHES:
        e = evidence[branch]
        report_lines.append(f"- **{branch}**: {e['n_genes']} genes; exact 2-NAP human {e['exact_2NAP_human']}; exact 2-NAP experimental {e['exact_2NAP_experimental']}; parent naphthalene {e['parent_naphthalene']}; multi-source (>=2) {e['multi_source_ge2']}.")
    report_lines.extend(["", "## Enrichment audit", ""])
    for a in enrichment_audit:
        report_lines.append(f"- **{a['branch_class']}**: status {a['status']}; retrieved {a['retrieved_n_terms']} terms; g:SCS-significant terms {a['significant_n_terms']}; query genes {a['query_n_genes']}; effective background is API-reported per term.")
    report_lines.extend(["", "## PPI/network audit", ""])
    for m in network_metrics:
        report_lines.append(f"- **{m['branch_class']}**: {m['network_n_nodes']} network nodes, {m['network_n_edges']} edges, {m['network_n_components']} components, largest component {m['largest_component_n']}, {m['module_n']} modules; STRING PPI-enrichment p={m['ppi_enrichment_p_value']}.")
    priority = [r["gene_symbol"] for r in hub_rows if as_bool(r.get("priority_hub_candidate"))]
    report_lines.extend(["", "## Shared-core hub candidates", "", f"The shared-core evidence/network gate identifies **{len(priority)}** priority candidates: {', '.join(priority) if priority else 'none'}.", "", "The complete ranked shared centrality table is in `07_shared_hub_candidates.csv`; it includes genes that fail one or more evidence/centrality gates so that exclusions remain auditable.", ""])
    (OUT / "URXP02_M3_DISEASE_BRANCH_REPORT.md").write_text("\n".join(report_lines), encoding="utf-8")

    manifest = {
        "analysis": "URXP02 M3 disease-branch enrichment and PPI/module analysis",
        "generated_at_utc": now_utc(),
        "inputs": {str(CLASSIFICATION): sha256(CLASSIFICATION), str(M2_SUMMARY): sha256(M2_SUMMARY)},
        "enrichment_branch_counts": {branch: len(by_branch[branch]) for branch in BRANCHES},
        "network_branch_counts": {branch: len(network_gene_sets[branch]) for branch in NETWORK_BRANCHES},
        "intersection_counts": {"thyroid_branch_total": len(network_gene_sets["thyroid-branch"]), "hypertension_branch_total": len(network_gene_sets["hypertension-branch"]), "shared": len(network_gene_sets["shared-core"])},
        "evidence_composition": evidence,
        "enrichment": {
            "service": "g:Profiler g:GOSt",
            "endpoint": GPROFILER_URL,
            "organism": "hsapiens",
            "sources": ENRICHMENT_SOURCES,
            "background": "all 828 M1b expanded-universe gene symbols",
            "domain_scope": "custom",
            "user_threshold": 0.05,
            "multiple_testing": "g:SCS",
            "all_results": True,
            "no_evidences": True,
            "audit": enrichment_audit,
        },
        "ppi": {
            "database": "STRING",
            "base": STRING_BASE,
            "species": SPECIES,
            "required_score": STRING_REQUIRED_SCORE,
            "network_type": STRING_NETWORK_TYPE,
            "add_nodes": 0,
            "module_algorithm": "NetworkX greedy_modularity_communities(weight=score)",
            "centrality_rule": "top 10% by degree, betweenness, or PageRank",
            "metrics": network_metrics,
        },
        "constraints": ["No figures", "No NHANES rerun", "No single-cell/tissue analysis", "No sex-specific molecular claim", "No pathway-to-disease causality claim", "No composite hub score"],
        "files": {},
    }
    for path in generated_files + [OUT / "URXP02_M3_DISEASE_BRANCH_REPORT.md"]:
        manifest["files"][path.name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
