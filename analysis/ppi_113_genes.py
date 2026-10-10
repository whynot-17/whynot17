from __future__ import annotations

import io
import json
import time
from datetime import date, datetime, timezone
from pathlib import Path

import networkx as nx
import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
QUERY_FILE = ROOT / "results" / "2026-10-10" / "geo_replication_113_replacement" / "query_113_genes.txt"
OUT = ROOT / "results" / "2026-10-10" / "ppi_113_genes"
STRING_BASE = "https://version-12-5.string-db.org/api"
SPECIES = 9606
CALLER = "DINP_CRC_PPI_113"


def post_tsv(method: str, params: dict) -> pd.DataFrame:
    response = requests.post(
        f"{STRING_BASE}/tsv/{method}",
        data=params,
        headers={"User-Agent": "DINP-CRC-PPI/1.0"},
        timeout=180,
    )
    response.raise_for_status()
    if not response.text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(response.text), sep="\t")


def read_query() -> list[str]:
    genes = [line.strip().upper() for line in QUERY_FILE.read_text(encoding="utf-8").splitlines()]
    genes = sorted(set(gene for gene in genes if gene))
    if len(genes) != 113:
        raise ValueError(f"Expected 113 unique query genes; found {len(genes)}")
    return genes


def map_query(genes: list[str]) -> pd.DataFrame:
    params = {
        "identifiers": "\r".join(genes),
        "species": SPECIES,
        "echo_query": 1,
        "caller_identity": CALLER,
    }
    mapped = post_tsv("get_string_ids", params)
    required = ["queryItem", "queryIndex", "stringId", "ncbiTaxonId", "taxonName", "preferredName"]
    missing_cols = [column for column in required if column not in mapped.columns]
    if missing_cols:
        raise RuntimeError(f"STRING ID mapping response missing columns: {missing_cols}")
    return mapped[required].drop_duplicates("queryItem")


def get_edges(string_ids: list[str], id_to_symbol: dict[str, str], network_type: str, score: int) -> pd.DataFrame:
    params = {
        "identifiers": "\r".join(string_ids),
        "species": SPECIES,
        "required_score": score,
        "add_nodes": 0,
        "network_type": network_type,
        "caller_identity": CALLER,
    }
    raw = post_tsv("network", params)
    columns = ["gene_A", "gene_B", "combined_score", "network_type"]
    if raw.empty:
        return pd.DataFrame(columns=columns)
    edges = raw.copy()
    edges["gene_A"] = edges["stringId_A"].map(id_to_symbol).fillna(edges["preferredName_A"])
    edges["gene_B"] = edges["stringId_B"].map(id_to_symbol).fillna(edges["preferredName_B"])
    edges["combined_score"] = pd.to_numeric(edges["score"], errors="coerce") / 1000
    edges = edges[edges["gene_A"] != edges["gene_B"]].copy()
    edges["_pair"] = edges.apply(lambda row: "|".join(sorted([row["gene_A"], row["gene_B"]])), axis=1)
    edges = edges.drop_duplicates("_pair").drop(columns="_pair")
    edges["network_type"] = network_type
    return edges[columns].sort_values(["combined_score", "gene_A", "gene_B"], ascending=[False, True, True])


def summarize_network(edges: pd.DataFrame, genes: list[str], network_type: str, threshold: float) -> tuple[dict, pd.DataFrame, nx.Graph]:
    graph = nx.Graph(network_type=network_type, confidence_threshold=threshold)
    graph.add_nodes_from(genes)
    for row in edges.itertuples(index=False):
        graph.add_edge(row.gene_A, row.gene_B, combined_score=float(row.combined_score))
    components = sorted(nx.connected_components(graph), key=lambda group: (-len(group), sorted(group)[0]))
    component_ids = {gene: f"C{index}" for index, group in enumerate(components, start=1) for gene in group}
    communities = list(nx.community.greedy_modularity_communities(graph, weight="combined_score")) if graph.number_of_edges() else []
    community_ids = {gene: f"M{index}" for index, group in enumerate(communities, start=1) for gene in group}
    degree = dict(graph.degree())
    degree_centrality = nx.degree_centrality(graph)
    betweenness = nx.betweenness_centrality(graph, normalized=True)
    nodes = pd.DataFrame([
        {
            "gene_symbol": gene,
            "degree": degree.get(gene, 0),
            "degree_centrality": degree_centrality.get(gene, 0.0),
            "betweenness_centrality": betweenness.get(gene, 0.0),
            "connected_component": component_ids.get(gene, ""),
            "modularity_community": community_ids.get(gene, "singleton"),
        }
        for gene in genes
    ]).sort_values(["degree", "gene_symbol"], ascending=[False, True])
    summary = {
        "network_type": network_type,
        "confidence_threshold": threshold,
        "query_gene_count": len(genes),
        "mapped_query_gene_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "mean_degree_over_query_genes": (2 * graph.number_of_edges() / graph.number_of_nodes()) if graph.number_of_nodes() else 0,
        "connected_component_count": len(components),
        "component_sizes": [len(group) for group in components],
        "largest_component_genes": sorted(components[0]) if components else [],
        "modularity_communities": [sorted(group) for group in communities],
    }
    return summary, nodes, graph


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = read_query()
    (OUT / "query_113_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    mapping = map_query(genes)
    mapping.to_csv(OUT / "STRING_mapping_query_113.tsv", sep="\t", index=False)
    mapped_symbols = set(mapping["queryItem"].astype(str).str.upper())
    unmapped = sorted(set(genes) - mapped_symbols)
    if unmapped:
        raise RuntimeError(f"Unmapped query genes in STRING: {unmapped}")
    id_to_symbol = dict(zip(mapping["stringId"].astype(str), mapping["queryItem"].astype(str)))
    ids = mapping["stringId"].astype(str).tolist()

    network_summaries = {}
    for network_type in ("functional", "physical"):
        for score in (400, 700):
            time.sleep(1)
            edges = get_edges(ids, id_to_symbol, network_type, score)
            suffix = f"{score / 1000:.3f}"
            edges.to_csv(OUT / f"PPI_edges_{network_type}_score_ge_{suffix}.tsv", sep="\t", index=False)
            summary, nodes, graph = summarize_network(edges, genes, network_type, score / 1000)
            nodes.to_csv(OUT / f"PPI_{network_type}_node_metrics_score_ge_{suffix}.tsv", sep="\t", index=False)
            nx.write_graphml(graph, OUT / f"PPI_{network_type}_network_score_ge_{suffix}.graphml")
            network_summaries[f"{network_type}_score_ge_{suffix}"] = summary

    try:
        version_response = requests.get(f"{STRING_BASE}/tsv/version", timeout=60)
        version_response.raise_for_status()
        string_version = pd.read_csv(io.StringIO(version_response.text), sep="\t").to_dict(orient="records")
    except requests.RequestException:
        string_version = []
    metadata = {
        "analysis_date": date.today().isoformat(),
        "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "organism": "Homo sapiens (NCBI taxon 9606)",
        "query_gene_count": len(genes),
        "query_file": str(QUERY_FILE.relative_to(ROOT)),
        "unmapped_query_genes": unmapped,
        "string_api_base": STRING_BASE,
        "string_version": string_version,
        "network_thresholds": [0.4, 0.7],
        "network_types": ["functional", "physical"],
        "network_parameters": {"add_nodes": 0},
        "network_summaries": network_summaries,
        "note": "This network-only analysis does not require a DEG background. GO/KEGG enrichment will use the user's specified 7,748-DEG custom background when its complete gene-symbol list is available.",
    }
    (OUT / "analysis_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = [
        "# 113 个交集基因的 STRING PPI 网络",
        "",
        f"- 分析日期：{metadata['analysis_date']}",
        "- 查询集：DINP 与 TCGA/GDC CRC DEG 的 113 个交集基因。",
        "- 物种：Homo sapiens（NCBI taxon 9606）；STRING v12.5 API。",
        "- STRING 网络分为 functional（功能关联）和 physical（物理互作）两类，分别导出；置信分数阈值为 0.400 和 0.700。",
        "- `add_nodes=0`，网络只保留查询集内的节点，不增加外部互作蛋白。",
        "- Functional edges combine multiple evidence types and do not all represent direct physical binding.",
        "- GO/KEGG 富集不在本次 PPI 网络运行中计算；后续按用户指定的 7,748 个 DEG 自定义背景进行。",
        "",
        "## 网络摘要",
        "",
        "| 网络 | 阈值 | 节点 | 边 | 连通分量 | 分量规模 |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for key, summary in network_summaries.items():
        lines.append(
            f"| {summary['network_type']} | {summary['confidence_threshold']:.3f} | {summary['query_gene_count']} | {summary['edge_count']} | {summary['connected_component_count']} | {summary['component_sizes']} |"
        )
    lines.extend(["", "## 文件", "", "- `query_113_genes.txt`：本次查询集。", "- `STRING_mapping_query_113.tsv`：STRING ID 映射。", "- `PPI_edges_*.tsv`：边表，可导入 Cytoscape。", "- `PPI_*_node_metrics_*.tsv`：节点度、中心性、连通分量和模块。", "- `PPI_*_network_*.graphml`：GraphML 网络，可导入 Cytoscape/Gephi。", "- `analysis_metadata.json`：API 版本、请求参数和网络统计。", ""])
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"query_gene_count": len(genes), "mapped": len(mapping), "summaries": network_summaries, "output": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
