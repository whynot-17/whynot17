from __future__ import annotations

import io
import json
import time
from collections import defaultdict
from datetime import date, datetime, timezone
from pathlib import Path

import networkx as nx
import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "2026-10-10" / "geo_replication_113"
OUT = ROOT / "results" / "2026-10-10" / "go_kegg_ppi_23_vs_113"
QUERY_FILE = RESULTS / "stable_genes_replicated_in_3_or_more_cohorts.tsv"
BACKGROUND_FILE = ROOT / "results" / "2026-10-10" / "DINP_286_x_TCGA_GDC_FDR_7748.tsv"

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
STRING_BASE = "https://version-12-5.string-db.org/api"
STRING_SCORE = 700
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG"]


def read_symbols(path: Path, column: str) -> list[str]:
    frame = pd.read_csv(path, sep="\t")
    values = frame[column].astype(str).str.strip().str.upper()
    return sorted(set(v for v in values if v and v != "NAN"))


def post_tsv(url: str, params: dict) -> pd.DataFrame:
    response = None
    for attempt in range(1, 4):
        try:
            response = requests.post(
                url,
                data=params,
                headers={"User-Agent": "DINP-CRC-enrichment/1.0"},
                timeout=180,
            )
            break
        except requests.exceptions.RequestException:
            if attempt == 3:
                raise
            time.sleep(attempt * 2)
    assert response is not None
    response.raise_for_status()
    if not response.text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(response.text), sep="\t")


def query_string_endpoint(method: str, params: dict) -> pd.DataFrame:
    return post_tsv(f"{STRING_BASE}/tsv/{method}", params)


def run_gprofiler(query: list[str], background: list[str], correction: str) -> tuple[dict, pd.DataFrame]:
    payload = {
        "organism": "hsapiens",
        "query": query,
        "sources": SOURCES,
        "user_threshold": 1.0,
        "significance_threshold_method": correction,
        "domain_scope": "custom",
        "background": background,
        "no_iea": False,
        "no_evidences": False,
    }
    request_path = OUT / f"gprofiler_request_{correction}.json"
    response_path = OUT / f"gprofiler_response_{correction}.json"
    if request_path.exists() and response_path.exists() and json.loads(request_path.read_text(encoding="utf-8")) == payload:
        data = json.loads(response_path.read_text(encoding="utf-8"))
    else:
        response = None
        for attempt in range(1, 4):
            try:
                response = requests.post(
                    GPROFILER_URL,
                    json=payload,
                    headers={"User-Agent": "DINP-CRC-enrichment/1.0"},
                    timeout=180,
                )
                break
            except requests.exceptions.RequestException:
                if attempt == 3:
                    raise
                time.sleep(attempt * 2)
        assert response is not None
        response.raise_for_status()
        data = response.json()
    meta = data.get("meta", {})
    query_order = meta.get("query_metadata", {}).get("queries", {}).get("query_1", query)
    flat_rows = []
    for term in data.get("result", []):
        if term.get("source") == "KEGG" and term.get("native") == "KEGG:00000":
            continue  # database root aggregator, not a biological pathway
        evidence = term.get("intersections") or []
        hit_genes = [gene for gene, codes in zip(query_order, evidence) if codes]
        flat_rows.append({
            "source": term.get("source"),
            "term_id": term.get("native"),
            "term_name": term.get("name"),
            f"p_value_{correction}": term.get("p_value"),
            f"significant_{correction}_lt_0.05": bool(term.get("p_value", 1) < 0.05),
            "query_size_mapped": term.get("query_size"),
            "background_size_effective": term.get("effective_domain_size"),
            "overlap_size": term.get("intersection_size"),
            "term_size": term.get("term_size"),
            "overlap_genes": ";".join(hit_genes),
            "description": term.get("description", ""),
        })
    frame = pd.DataFrame(flat_rows)
    if not frame.empty:
        frame = frame.sort_values(["source", f"p_value_{correction}", "term_id"], kind="stable")
    request_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    response_path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    return data, frame


def ppi_edges(query_ids: list[str], network_type: str, score: int) -> pd.DataFrame:
    params = {
        "identifiers": "\r".join(query_ids),
        "species": 9606,
        "required_score": score,
        "add_nodes": 0,
        "network_type": network_type,
        "caller_identity": "DINP_CRC_functional_analysis",
    }
    return query_string_endpoint("network", params)


def get_network_edges(query_ids: list[str], network_type: str, score: int, id_to_symbol: dict[str, str]) -> pd.DataFrame:
    suffix = f"{score / 1000:.3f}"
    path = OUT / f"PPI_edges_{network_type}_score_ge_{suffix}.tsv"
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return map_edges(ppi_edges(query_ids, network_type, score), id_to_symbol)


def map_edges(edges: pd.DataFrame, id_to_symbol: dict[str, str]) -> pd.DataFrame:
    if edges.empty:
        return pd.DataFrame(columns=["gene_A", "gene_B", "combined_score", "network_type"])
    out = edges.copy()
    out["gene_A"] = out["stringId_A"].map(id_to_symbol).fillna(out["preferredName_A"])
    out["gene_B"] = out["stringId_B"].map(id_to_symbol).fillna(out["preferredName_B"])
    out["combined_score"] = pd.to_numeric(out["score"], errors="coerce")
    out = out[out["gene_A"] != out["gene_B"]].copy()
    out["_pair"] = out.apply(lambda r: "|".join(sorted([r["gene_A"], r["gene_B"]])), axis=1)
    out = out.drop_duplicates("_pair").drop(columns="_pair")
    return out


def graph_metrics(edges: pd.DataFrame, genes: list[str], label: str, score: int) -> tuple[dict, pd.DataFrame, nx.Graph]:
    graph = nx.Graph(network_type=label, confidence_threshold=score / 1000)
    graph.add_nodes_from(genes)
    for _, row in edges.iterrows():
        graph.add_edge(
            row["gene_A"], row["gene_B"],
            combined_score=float(row["combined_score"]),
        )
    components = sorted(nx.connected_components(graph), key=lambda c: (-len(c), sorted(c)[0] if c else ""))
    comp_id = {gene: f"C{i + 1}" for i, comp in enumerate(components) for gene in comp}
    communities = []
    if graph.number_of_edges():
        communities = list(nx.community.greedy_modularity_communities(graph, weight="combined_score"))
    community_id = {gene: f"M{i + 1}" for i, group in enumerate(communities) for gene in group}
    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph, normalized=True, weight=None)
    degree_column = f"degree_at_score_ge_{score / 1000:.3f}"
    degree_centrality = nx.degree_centrality(graph)
    nodes = pd.DataFrame([
        {
            "gene_symbol": gene,
            degree_column: degree.get(gene, 0),
            "degree_centrality": degree_centrality.get(gene, 0.0),
            "betweenness_centrality": betweenness.get(gene, 0.0),
            "connected_component": comp_id.get(gene, ""),
            "modularity_community": community_id.get(gene, "singleton"),
        }
        for gene in genes
    ]).sort_values([degree_column, "gene_symbol"], ascending=[False, True])
    stats = {
        "network_type": label,
        "confidence_threshold": score / 1000,
        "n_query_genes": len(genes),
        "n_nodes": graph.number_of_nodes(),
        f"n_edges_score_ge_{score / 1000:.3f}": graph.number_of_edges(),
        "mean_degree_over_all_query_nodes": (2 * graph.number_of_edges() / graph.number_of_nodes()) if graph.number_of_nodes() else 0,
        "connected_components": len(components),
        "component_sizes": [len(x) for x in components],
        "largest_component_genes": sorted(components[0]) if components else [],
        "greedy_modularity_communities": [sorted(x) for x in communities],
    }
    return stats, nodes, graph


def reduce_go_redundancy(significant: pd.DataFrame, jaccard_cutoff: float = 0.5) -> pd.DataFrame:
    rows = []
    go = significant[significant["source"].str.startswith("GO:")].copy()
    for source, block in go.groupby("source", sort=True):
        block = block.reset_index(drop=True)
        graph = nx.Graph()
        graph.add_nodes_from(range(len(block)))
        gene_sets = [set(str(x).split(";")) - {""} for x in block["overlap_genes"]]
        for i in range(len(block)):
            for j in range(i + 1, len(block)):
                union = gene_sets[i] | gene_sets[j]
                jaccard = len(gene_sets[i] & gene_sets[j]) / len(union) if union else 0.0
                if jaccard >= jaccard_cutoff:
                    graph.add_edge(i, j, jaccard=jaccard)
        for cluster_n, cluster in enumerate(nx.connected_components(graph), start=1):
            members = block.iloc[sorted(cluster)].copy()
            representative = members.sort_values(["p_value_g_SCS", "overlap_size"], ascending=[True, False]).iloc[0].to_dict()
            representative["redundancy_cluster"] = f"{source}_R{cluster_n}"
            representative["cluster_term_count"] = len(members)
            representative["cluster_terms"] = ";".join(members["term_name"].astype(str).tolist())
            rows.append(representative)
    if not rows:
        return pd.DataFrame(columns=list(significant.columns) + ["redundancy_cluster", "cluster_term_count", "cluster_terms"])
    return pd.DataFrame(rows).sort_values(["source", "p_value_gSCS"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    query = read_symbols(QUERY_FILE, "gene_symbol")
    background = read_symbols(BACKGROUND_FILE, "gene_symbol")
    if len(query) != 23 or len(background) != 113:
        raise ValueError(f"Expected 23 query and 113 background genes; found {len(query)} and {len(background)}")
    if set(query) - set(background):
        raise ValueError(f"Query genes outside the 113-gene background: {sorted(set(query) - set(background))}")
    (OUT / "query_23_genes.txt").write_text("\n".join(query) + "\n", encoding="utf-8")
    (OUT / "background_113_genes.txt").write_text("\n".join(background) + "\n", encoding="utf-8")

    gprofiler_data, all_terms = run_gprofiler(query, background, "g_SCS")
    time.sleep(1)
    gprofiler_fdr_data, fdr_terms = run_gprofiler(query, background, "fdr")
    significant = all_terms[all_terms["significant_g_SCS_lt_0.05"]].copy()
    significant_fdr = fdr_terms[fdr_terms["significant_fdr_lt_0.05"]].copy()
    significant.to_csv(OUT / "GO_KEGG_significant_gSCS.tsv", sep="\t", index=False)
    all_terms.to_csv(OUT / "GO_KEGG_all_results.tsv", sep="\t", index=False)
    significant_fdr.to_csv(OUT / "GO_KEGG_significant_FDR_sensitivity.tsv", sep="\t", index=False)
    fdr_terms.to_csv(OUT / "GO_KEGG_all_results_FDR_sensitivity.tsv", sep="\t", index=False)
    if not fdr_terms.empty:
        top_non_sig = fdr_terms[~fdr_terms["significant_fdr_lt_0.05"]].groupby("source", group_keys=False).head(5)
    else:
        top_non_sig = fdr_terms.copy()
    top_non_sig.to_csv(OUT / "GO_KEGG_top_ranked_non_significant.tsv", sep="\t", index=False)
    go_representatives = reduce_go_redundancy(significant)
    go_representatives.to_csv(OUT / "GO_significant_representative_terms.tsv", sep="\t", index=False)
    significant[significant["source"] == "KEGG"].to_csv(OUT / "KEGG_significant_gSCS.tsv", sep="\t", index=False)

    cached_metadata_path = OUT / "analysis_metadata.json"
    try:
        string_version = requests.get(f"{STRING_BASE}/tsv/version", timeout=60)
        string_version.raise_for_status()
        version_frame = pd.read_csv(io.StringIO(string_version.text), sep="\t")
    except requests.exceptions.RequestException:
        cached_meta = json.loads(cached_metadata_path.read_text(encoding="utf-8")) if cached_metadata_path.exists() else {}
        version_frame = pd.DataFrame(cached_meta.get("string_version", []))
        if version_frame.empty:
            raise
    map_params = {
        "identifiers": "\r".join(background),
        "species": 9606,
        "echo_query": 1,
        "caller_identity": "DINP_CRC_functional_analysis",
    }
    mapping_path = OUT / "STRING_mapping_background_113.tsv"
    if mapping_path.exists():
        background_map = pd.read_csv(mapping_path, sep="\t")
    else:
        time.sleep(1)
        background_map = query_string_endpoint("get_string_ids", map_params)
    mapping_columns = ["queryItem", "queryIndex", "stringId", "ncbiTaxonId", "taxonName", "preferredName"]
    background_map[mapping_columns].to_csv(OUT / "STRING_mapping_background_113.tsv", sep="\t", index=False)
    background_map[background_map["queryItem"].isin(query)][mapping_columns].to_csv(OUT / "STRING_mapping_query_23.tsv", sep="\t", index=False)
    query_map = background_map[background_map["queryItem"].isin(query)].drop_duplicates("queryItem")
    if set(query) - set(query_map["queryItem"].astype(str)):
        raise ValueError(f"Unmapped STRING query genes: {sorted(set(query) - set(query_map['queryItem'].astype(str)))}")
    if set(background) - set(background_map["queryItem"].astype(str)):
        raise ValueError(f"Unmapped STRING background genes: {sorted(set(background) - set(background_map['queryItem'].astype(str)))}")
    query_ids = query_map["stringId"].astype(str).tolist()
    background_ids = background_map["stringId"].drop_duplicates().astype(str).tolist()
    id_to_symbol = dict(zip(query_map["stringId"].astype(str), query_map["queryItem"].astype(str)))
    string_requests = {
        "species": 9606,
        "string_version": version_frame.to_dict(orient="records"),
        "confidence_thresholds": [400, 700],
        "query_23_string_ids": query_ids,
        "background_113_string_ids": background_ids,
        "network_parameters": {"add_nodes": 0, "network_types": ["functional", "physical"]},
        "ppi_enrichment_uses_background_string_identifiers": True,
    }
    (OUT / "STRING_request_parameters.json").write_text(json.dumps(string_requests, indent=2, ensure_ascii=False), encoding="utf-8")

    time.sleep(1)
    functional_edges = get_network_edges(query_ids, "functional", 700, id_to_symbol)
    functional_edges.to_csv(OUT / "PPI_edges_functional_score_ge_0.700.tsv", sep="\t", index=False)
    time.sleep(1)
    physical_edges = get_network_edges(query_ids, "physical", 700, id_to_symbol)
    physical_edges.to_csv(OUT / "PPI_edges_physical_score_ge_0.700.tsv", sep="\t", index=False)
    time.sleep(1)
    ppi_params = {
        "identifiers": "\r".join(query_ids),
        "species": 9606,
        "required_score": 700,
        "background_string_identifiers": "\r".join(background_ids),
        "caller_identity": "DINP_CRC_functional_analysis",
    }
    ppi_path = OUT / "STRING_PPI_enrichment_custom_background.tsv"
    ppi_enrichment = pd.read_csv(ppi_path, sep="\t") if ppi_path.exists() else query_string_endpoint("ppi_enrichment", ppi_params)
    ppi_enrichment.to_csv(OUT / "STRING_PPI_enrichment_custom_background.tsv", sep="\t", index=False)

    time.sleep(1)
    functional_medium = get_network_edges(query_ids, "functional", 400, id_to_symbol)
    functional_medium.to_csv(OUT / "PPI_edges_functional_score_ge_0.400.tsv", sep="\t", index=False)
    time.sleep(1)
    physical_medium = get_network_edges(query_ids, "physical", 400, id_to_symbol)
    physical_medium.to_csv(OUT / "PPI_edges_physical_score_ge_0.400.tsv", sep="\t", index=False)
    time.sleep(1)
    ppi_medium_params = dict(ppi_params)
    ppi_medium_params["required_score"] = 400
    ppi_medium_path = OUT / "STRING_PPI_enrichment_custom_background_score_ge_0.400.tsv"
    ppi_medium = pd.read_csv(ppi_medium_path, sep="\t") if ppi_medium_path.exists() else query_string_endpoint("ppi_enrichment", ppi_medium_params)
    ppi_medium.to_csv(OUT / "STRING_PPI_enrichment_custom_background_score_ge_0.400.tsv", sep="\t", index=False)

    functional_stats, node_metrics, functional_graph = graph_metrics(functional_edges, query, "functional", 700)
    physical_stats, physical_nodes, physical_graph = graph_metrics(physical_edges, query, "physical", 700)
    functional_medium_stats, functional_medium_nodes, functional_medium_graph = graph_metrics(functional_medium, query, "functional", 400)
    physical_medium_stats, physical_medium_nodes, physical_medium_graph = graph_metrics(physical_medium, query, "physical", 400)
    node_metrics.to_csv(OUT / "PPI_functional_node_metrics.tsv", sep="\t", index=False)
    physical_nodes.to_csv(OUT / "PPI_physical_node_metrics.tsv", sep="\t", index=False)
    functional_medium_nodes.to_csv(OUT / "PPI_functional_node_metrics_score_ge_0.400.tsv", sep="\t", index=False)
    physical_medium_nodes.to_csv(OUT / "PPI_physical_node_metrics_score_ge_0.400.tsv", sep="\t", index=False)
    nx.write_graphml(functional_graph, OUT / "PPI_functional_network.graphml")
    nx.write_graphml(physical_graph, OUT / "PPI_physical_network.graphml")
    nx.write_graphml(functional_medium_graph, OUT / "PPI_functional_network_score_ge_0.400.graphml")
    nx.write_graphml(physical_medium_graph, OUT / "PPI_physical_network_score_ge_0.400.graphml")

    ppi_enrichment_rows = ppi_enrichment.to_dict(orient="records") if not ppi_enrichment.empty else []
    ppi_medium_rows = ppi_medium.to_dict(orient="records") if not ppi_medium.empty else []
    gprof_meta = gprofiler_data.get("meta", {})
    gprof_fdr_meta = gprofiler_fdr_data.get("meta", {})
    source_meta = {
        "analysis_date": date.today().isoformat(),
        "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "organism": "Homo sapiens (NCBI taxon 9606)",
        "query_size": len(query),
        "background_size": len(background),
        "query_all_within_background": True,
        "gprofiler_database_metadata": gprof_meta.get("result_metadata", {}),
        "gprofiler_query_metadata": gprof_meta.get("query_metadata", {}),
        "gprofiler_mapping_metadata": gprof_meta.get("genes_metadata", {}),
        "gprofiler_fdr_metadata": gprof_fdr_meta.get("query_metadata", {}),
        "gprofiler_data_version": gprof_meta.get("version"),
        "gprofiler_timestamp": gprof_meta.get("timestamp"),
        "gprofiler_fdr_data_version": gprof_fdr_meta.get("version"),
        "gprofiler_fdr_timestamp": gprof_fdr_meta.get("timestamp"),
        "string_version": version_frame.to_dict(orient="records"),
        "string_api_base": STRING_BASE,
        "string_confidence_cutoff": STRING_SCORE,
        "string_functional_network_stats": functional_stats,
        "string_physical_network_stats": physical_stats,
        "string_functional_network_stats_score_ge_0.400": functional_medium_stats,
        "string_physical_network_stats_score_ge_0.400": physical_medium_stats,
        "string_ppi_enrichment_custom_background": ppi_enrichment_rows,
        "string_ppi_enrichment_custom_background_score_ge_0.400": ppi_medium_rows,
        "gprofiler_significant_terms_by_source": significant.groupby("source").size().to_dict() if not significant.empty else {},
        "gprofiler_significant_terms_by_source_FDR_sensitivity": significant_fdr.groupby("source").size().to_dict() if not significant_fdr.empty else {},
        "gprofiler_min_corrected_p_by_source_FDR_sensitivity": fdr_terms.groupby("source")["p_value_fdr"].min().to_dict() if not fdr_terms.empty else {},
        "gprofiler_result_effective_domain_sizes_by_source": all_terms.groupby("source")["background_size_effective"].agg(["min", "max"]).to_dict(orient="index") if not all_terms.empty else {},
        "gprofiler_source_urls": {
            "API": GPROFILER_URL,
            "API_documentation": "https://biit.cs.ut.ee/gprofiler/page/apis",
        },
        "string_source_urls": {
            "API_documentation": "https://string-db.org/help/api/",
            "version": f"{STRING_BASE}/tsv/version",
        },
    }
    (OUT / "analysis_metadata.json").write_text(json.dumps(source_meta, indent=2, ensure_ascii=False, default=str), encoding="utf-8")

    report = [
        "# 23 个稳定复现基因的 GO/KEGG 富集与 PPI 网络分析",
        "",
        f"- 日期：{date.today().isoformat()}",
        "- 物种：Homo sapiens（taxon 9606）；基因 ID：HGNC gene symbol。",
        f"- 查询集：23 个在至少 3/4 GEO accession 达到此前复现标准的基因；背景集：本项目 113 个 DINP×CRC 交集基因。查询集全部包含于背景集。",
        "- GO/KEGG：使用 g:Profiler g:GOSt，来源 GO:BP、GO:MF、GO:CC、KEGG；custom background（113 symbols）；g:SCS 多重检验校正；显著阈值为校正后 p < 0.05。ORA 不区分肿瘤上调/下调方向。",
        "- STRING PPI：STRING v12.5，人类；功能关联网络与物理相互作用网络分开导出；置信分数阈值 700/1000；STRING PPI enrichment 使用 113 基因背景中成功映射到 STRING 的蛋白 ID。功能关联包含直接物理证据以外的已知/预测功能联系，不能将每条边都解释为直接结合。",
        "- PPI 网络节点仅为 23 个查询基因；add_nodes=0，不添加外部互作蛋白。网络指标在指定 score 阈值下计算；degree/betweenness 是拓扑描述，不代表因果或调控方向。",
        "",
        "## 显著富集项数量",
        "",
        "| 数据源 | g:SCS < 0.05 项数 |",
        "|---|---:|",
    ]
    sig_counts = significant.groupby("source").size().to_dict() if not significant.empty else {}
    fdr_sig_counts = significant_fdr.groupby("source").size().to_dict() if not significant_fdr.empty else {}
    fdr_min = fdr_terms.groupby("source")["p_value_fdr"].min().to_dict() if not fdr_terms.empty else {}
    for source in SOURCES:
        report.append(f"| {source} | {sig_counts.get(source, 0)} |")
    report.extend(["", "BH-FDR 结果作为敏感性分析单独报告；与 g:SCS 分开，不混作主分析。", "", "| 数据源 | BH-FDR < 0.05 项数 | 最小 BH-FDR |", "|---|---:|---:|"])
    for source in SOURCES:
        min_value = fdr_min.get(source)
        min_text = f"{min_value:.3g}" if min_value is not None else "NA"
        report.append(f"| {source} | {fdr_sig_counts.get(source, 0)} | {min_text} |")
    report.extend(["", "## GO 代表性条目", ""])
    if go_representatives.empty:
        report.append("无 GO 条目达到 g:SCS < 0.05。")
    else:
        report.extend([
            "下表将显著 GO 条目按重叠基因 Jaccard ≥0.5 合并相近项，每簇列出校正 p 值最小的代表条目。完整显著项见 GO_KEGG_significant_gSCS.tsv。",
            "",
            "| GO 类别 | 代表条目 | g:SCS p | 重叠基因 | 簇内条目数 |",
            "|---|---|---:|---|---:|",
        ])
        for _, row in go_representatives.head(30).iterrows():
            report.append(f"| {row['source']} | {row['term_name']} ({row['term_id']}) | {row['p_value_g_SCS']:.3g} | {row['overlap_genes']} | {row['cluster_term_count']} |")
    report.extend(["", "## KEGG 显著通路", ""])
    kegg_sig = significant[significant["source"] == "KEGG"]
    if kegg_sig.empty:
        report.append("无 KEGG 通路达到 g:SCS < 0.05。")
    else:
        report.extend(["| 通路 | g:SCS p | 重叠基因 | 命中/通路基因数 |", "|---|---:|---|---:|"])
        for _, row in kegg_sig.iterrows():
            report.append(f"| {row['term_name']} ({row['term_id']}) | {row['p_value_g_SCS']:.3g} | {row['overlap_genes']} | {row['overlap_size']}/{row['term_size']} |")
    if not significant_fdr.empty:
        report.extend(["", "### BH-FDR 敏感性分析中显著的条目", "", "| 数据源 | 条目 | BH-FDR p | 重叠基因 |", "|---|---|---:|---|"])
        for _, row in significant_fdr.sort_values(["source", "p_value_fdr"]).iterrows():
            report.append(f"| {row['source']} | {row['term_name']} ({row['term_id']}) | {row['p_value_fdr']:.3g} | {row['overlap_genes']} |")
    else:
        report.extend(["", "### 最靠前但未达显著的条目（仅作探索）", "", "以下按 BH-FDR 排序；全部仍未达到 0.05，不作为显著富集结论。", "", "| 数据源 | 条目 | BH-FDR p | 重叠基因 |", "|---|---|---:|---|"])
        for _, row in top_non_sig.iterrows():
            report.append(f"| {row['source']} | {row['term_name']} ({row['term_id']}) | {row['p_value_fdr']:.3g} | {row['overlap_genes']} |")
    report.extend(["", "## PPI 网络摘要", ""])
    for stats in (functional_stats, physical_stats, functional_medium_stats, physical_medium_stats):
        edge_count = stats[f"n_edges_score_ge_{stats['confidence_threshold']:.3f}"]
        report.append(
            f"- {stats['network_type']}（score≥{stats['confidence_threshold']:.3f}）：{stats['n_nodes']} 个查询节点、{edge_count} 条边，"
            f"平均度 {stats['mean_degree_over_all_query_nodes']:.2f}，连通分量 {stats['connected_components']} 个；"
            f"分量规模 {stats['component_sizes']}。"
        )
    hubs = functional_medium_nodes.sort_values(["degree_at_score_ge_0.400", "gene_symbol"], ascending=[False, True]).head(5)
    report.append("- score≥0.400 功能网络的度最高节点：" + ", ".join(f"{row['gene_symbol']} (degree={int(row['degree_at_score_ge_0.400'])})" for _, row in hubs.iterrows()) + ".")
    report.append("- score≥0.400 功能网络最大连通分量：" + ", ".join(functional_medium_stats["largest_component_genes"]) + ".")
    report.append("- score≥0.700 功能网络边：" + "; ".join(f"{r.gene_A}–{r.gene_B} ({r.combined_score:.3f})" for r in functional_edges.sort_values("combined_score", ascending=False).itertuples()) + ".")
    report.append("- score≥0.400 物理网络边：" + ("; ".join(f"{r.gene_A}–{r.gene_B} ({r.combined_score:.3f})" for r in physical_medium.sort_values("combined_score", ascending=False).itertuples()) if not physical_medium.empty else "无") + ".")
    if ppi_enrichment_rows:
        report.append(f"- STRING PPI enrichment（自定义背景，score≥0.700）：{json.dumps(ppi_enrichment_rows[0], ensure_ascii=False, default=str)}")
    if ppi_medium_rows:
        report.append(f"- STRING PPI enrichment（自定义背景，score≥0.400）：{json.dumps(ppi_medium_rows[0], ensure_ascii=False, default=str)}")
    report.extend([
        "",
        "## 文件清单",
        "",
        "- GO_KEGG_all_results.tsv / GO_KEGG_significant_gSCS.tsv：完整结果与显著项。",
        "- GO_KEGG_all_results_FDR_sensitivity.tsv / GO_KEGG_significant_FDR_sensitivity.tsv：BH-FDR 敏感性分析。",
        "- GO_KEGG_top_ranked_non_significant.tsv：每个数据源 BH-FDR 最靠前的未显著条目。",
        "- GO_significant_representative_terms.tsv：GO 重叠条目聚类后的代表项。",
        "- KEGG_significant_gSCS.tsv：KEGG 显著通路单独表。",
        "- PPI_edges_*_score_ge_0.700.tsv 与 PPI_edges_*_score_ge_0.400.tsv：STRING 边表，可导入 Cytoscape。",
        "- PPI_functional_node_metrics.tsv、PPI_physical_node_metrics.tsv：节点度、中心性、分量与模块。",
        "- PPI_*_network*.graphml：图网络文件，可导入 Cytoscape/Gephi。",
        "- STRING_PPI_enrichment_custom_background*.tsv：STRING PPI enrichment 数值，使用 113 基因自定义背景。",
        "- STRING_mapping_background_113.tsv、STRING_mapping_query_23.tsv：STRING ID 映射。",
        "- STRING_request_parameters.json：查询及背景的 STRING ID 与参数。",
        "- query_23_genes.txt、background_113_genes.txt：本次实际使用的输入集。",
        "- gprofiler_request_g_SCS.json、gprofiler_response_g_SCS.json、gprofiler_request_fdr.json、gprofiler_response_fdr.json、analysis_metadata.json：API 参数、原始响应与版本/统计元数据。",
        "- 复现脚本：analysis/enrichment_ppi_23_vs_113.py。",
        "",
        "## 解释边界",
        "",
        "这是基于阈值化基因集的过度代表分析（ORA），结果表示注释条目的基因在 23 个基因中相对 113 个背景基因富集，不表示通路被激活或抑制。GO 各项存在层级和重叠，代表表按基因交集进行简化。STRING 功能关联网络汇总多种证据来源；应区分功能关联与物理相互作用。",
    ])
    (OUT / "README.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"query={len(query)} background={len(background)}")
    print("GO/KEGG significant by source:", sig_counts)
    print("STRING mapped background/query:", len(background_map), len(query_map))
    print("Functional PPI stats:", json.dumps(functional_stats, ensure_ascii=False))
    print("Physical PPI stats:", json.dumps(physical_stats, ensure_ascii=False))
    print("STRING PPI enrichment:", json.dumps(ppi_enrichment_rows, ensure_ascii=False, default=str))
    print("Output:", OUT)


if __name__ == "__main__":
    main()
