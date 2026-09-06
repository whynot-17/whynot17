#!/usr/bin/env python3
"""Step 8C: STRING high-confidence network convergence.

The four Tier A networks are analyzed separately.  Primary inputs are the
frozen Step 7 overlap genes; the frozen union of all Step 7 cluster genes is
used only as an empirical randomization background.  STRING functional
associations at combined score >= 700 are used, with no added interactors.

The acquisition script uses pairwise blocks for large gene sets so every pair
of input blocks is queried despite STRING's request-size limit.  Edges are
deduplicated before analysis.  Community detection uses NetworkX Louvain with
a fixed seed.  Network enrichment is tested by degree-stratified permutation
against the frozen background graph, not by treating individual edges as
independent observations.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx


DEFAULT_DIR = Path(__file__).resolve().parent
TIER_A = ("cluster_11", "cluster_5", "cluster_6", "cluster_8")
SCORE_THRESHOLD = 0.700
N_PERMUTATIONS = 1000
RANDOM_SEED = 20260827


def read_gene_sets(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            sets[row["cluster_id"]].add(row["gene_symbol"])
    return dict(sets)


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def edge_key(row: dict[str, object]) -> tuple[str, str]:
    left = str(row.get("stringId_A", ""))
    right = str(row.get("stringId_B", ""))
    return tuple(sorted((left, right)))


def load_edges(raw_dir: Path, prefix: str) -> dict[tuple[str, str], float]:
    edges: dict[tuple[str, str], float] = {}
    for path in sorted(raw_dir.glob(f"{prefix}*.json")):
        payload = read_json(path)
        if not isinstance(payload, list):
            continue
        for row in payload:
            if not isinstance(row, dict):
                continue
            key = edge_key(row)
            if not key[0] or key[0] == key[1]:
                continue
            score = float(row.get("score", 0.0))
            if score >= SCORE_THRESHOLD:
                edges[key] = max(score, edges.get(key, 0.0))
    return edges


def load_mapping(raw_dir: Path, cluster: str) -> tuple[dict[str, str], dict[str, list[str]], dict[str, object]]:
    payload = read_json(raw_dir / f"string_{cluster}_mapping.json")
    records = payload if isinstance(payload, list) else []
    by_query: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in records:
        if isinstance(row, dict) and row.get("queryItem"):
            by_query[str(row["queryItem"])].append(row)
    chosen: dict[str, str] = {}
    all_ids: dict[str, list[str]] = {}
    for query, rows in by_query.items():
        ids = sorted({str(row.get("stringId")) for row in rows if row.get("stringId")})
        all_ids[query] = ids
        if ids:
            chosen[query] = ids[0]
    audit = {
        "n_query_genes_in_mapping": len(by_query),
        "n_mapped_genes_in_mapping": len(chosen),
        "n_ambiguous_queries": sum(len(ids) > 1 for ids in all_ids.values()),
    }
    return chosen, all_ids, audit


def make_graph(edges: dict[tuple[str, str], float]) -> nx.Graph:
    graph = nx.Graph()
    for (left, right), score in edges.items():
        graph.add_edge(left, right, score=score)
    return graph


def degree_bins(degrees: dict[str, int], n_bins: int = 10) -> dict[str, int]:
    ordered = sorted(degrees, key=lambda node: (degrees[node], node))
    bins: dict[str, int] = {}
    for index, node in enumerate(ordered):
        bins[node] = min(n_bins - 1, index * n_bins // max(len(ordered), 1))
    return bins


def stratified_sample(
    candidates: list[str],
    target_nodes: set[str],
    bin_map: dict[str, int],
    rng: random.Random,
) -> set[str]:
    target_counts = Counter(bin_map[node] for node in target_nodes if node in bin_map)
    by_bin: dict[int, list[str]] = defaultdict(list)
    for node in candidates:
        if node not in target_nodes and node in bin_map:
            by_bin[bin_map[node]].append(node)
    sampled: set[str] = set()
    for bin_id, count in target_counts.items():
        pool = by_bin.get(bin_id, [])
        if len(pool) < count:
            # Preserve the requested size if a very sparse degree bin occurs.
            pool = [node for node in candidates if node not in target_nodes and node not in sampled]
        if len(pool) >= count:
            sampled.update(rng.sample(pool, count))
    if len(sampled) < len(target_nodes):
        pool = [node for node in candidates if node not in target_nodes and node not in sampled]
        sampled.update(rng.sample(pool, min(len(target_nodes) - len(sampled), len(pool))))
    return sampled


def parse_module_enrichment(raw_dir: Path, cluster: str, module_id: str, name_to_id: dict[str, str]) -> tuple[list[dict[str, object]], dict[str, int]]:
    path = raw_dir / f"string_enrichment_{cluster}_{module_id}.json"
    if not path.exists() or path.stat().st_size <= 2:
        return [], {}
    payload = read_json(path)
    if not isinstance(payload, list):
        return [], {}
    allowed = {"Process", "KEGG", "REACTOME", "Reactome", "Process"}
    rows = []
    recurrence: Counter[str] = Counter()
    for item in payload:
        if not isinstance(item, dict) or item.get("category") not in allowed:
            continue
        try:
            fdr = float(item.get("fdr", 1.0))
        except (TypeError, ValueError):
            fdr = 1.0
        if fdr >= 0.05:
            continue
        genes = item.get("preferredNames", item.get("inputGenes", []))
        if isinstance(genes, str):
            genes = [gene.strip() for gene in genes.split(",") if gene.strip()]
        if not isinstance(genes, list):
            genes = []
        string_ids = []
        for gene in genes:
            gene = str(gene)
            recurrence[name_to_id.get(gene, gene)] += 1
            string_ids.append(name_to_id.get(gene, gene))
        rows.append({
            "cluster_id": cluster,
            "module_id": module_id,
            "category": item.get("category", ""),
            "term": item.get("term", ""),
            "description": item.get("description", ""),
            "p_value": item.get("p_value", ""),
            "fdr": fdr,
            "number_of_genes": item.get("number_of_genes", ""),
            "input_string_ids": ";".join(sorted(set(string_ids))),
        })
    rows.sort(key=lambda row: (float(row["fdr"]), str(row["category"]), str(row["term"])))
    return rows, dict(recurrence)


def centrality_table(
    graph: nx.Graph,
    axis_nodes: set[str],
    community_by_node: dict[str, set[str]],
    pathway_recurrence: dict[str, int],
    name_by_id: dict[str, str],
) -> list[dict[str, object]]:
    sub = graph.subgraph(axis_nodes).copy()
    degree = dict(sub.degree())
    betweenness = nx.betweenness_centrality(sub, normalized=True) if sub.number_of_edges() else {node: 0.0 for node in sub}
    eigen = {node: 0.0 for node in sub}
    for component in nx.connected_components(sub):
        component_graph = sub.subgraph(component).copy()
        if len(component_graph) > 1 and component_graph.number_of_edges() > 0:
            eigen.update(nx.eigenvector_centrality(component_graph, max_iter=1000, tol=1e-08))
    max_degree = max(degree.values(), default=1)
    max_betweenness = max(betweenness.values(), default=1.0) or 1.0
    max_eigen = max(eigen.values(), default=1.0) or 1.0
    max_recurrence = max(pathway_recurrence.values(), default=1) or 1
    rows = []
    for node in sorted(axis_nodes):
        d = degree.get(node, 0)
        b = betweenness.get(node, 0.0)
        e = eigen.get(node, 0.0)
        community = community_by_node.get(node, {node})
        internal_degree = sub.subgraph(community).degree(node) if node in sub else 0
        module_connectivity = internal_degree / max(len(community) - 1, 1)
        recurrence = pathway_recurrence.get(node, 0)
        score = (
            d / max_degree
            + b / max_betweenness
            + e / max_eigen
            + module_connectivity
            + recurrence / max_recurrence
        ) / 5.0
        rows.append({
            "node": node,
            "preferred_name": name_by_id.get(node, node),
            "degree": d,
            "betweenness": b,
            "eigenvector": e,
            "module_size": len(community),
            "within_module_degree": internal_degree,
            "within_module_connectivity": module_connectivity,
            "pathway_recurrence_count": recurrence,
            "network_priority_score": score,
        })
    return sorted(rows, key=lambda row: (-float(row["network_priority_score"]), -int(row["degree"]), str(row["node"])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_DIR)
    parser.add_argument("--step7-membership", type=Path, default=DEFAULT_DIR.parent / "step07_genecard_convergence" / "t2d_cluster_ctd_gene_membership.csv")
    args = parser.parse_args()
    raw_dir = args.output_dir / "raw_string"
    if not raw_dir.exists():
        raise FileNotFoundError("STRING raw cache not found; run run_step08c_string_acquisition.ps1 first")
    gene_sets = read_gene_sets(args.step7_membership)
    background_genes = set().union(*(gene_sets.values()))
    background_edges = load_edges(raw_dir, "string_background_11clusters_network_")
    background_graph = make_graph(background_edges)
    background_chosen, _, background_mapping_audit = load_mapping(raw_dir, "background_11clusters")
    background_node_set = {background_chosen[gene] for gene in background_genes if gene in background_chosen} & set(background_graph.nodes())
    degree_map = dict(background_graph.degree())
    bins = degree_bins(degree_map)
    rng = random.Random(RANDOM_SEED)

    mapping_audit = [{
        "cluster_id": "background_11clusters",
        **background_mapping_audit,
        "n_step7_overlap_genes": len(background_genes),
        "n_mapped_string_ids": len(background_node_set),
        "n_unmapped_input_genes": len(background_genes - set(background_chosen)),
        "input_mapped_fraction": len(background_chosen) / len(background_genes) if background_genes else 0.0,
    }]
    all_edge_rows = []
    all_node_rows = []
    module_rows = []
    module_annotations = []
    network_summary = []
    for cluster in TIER_A:
        chosen, all_ids, audit = load_mapping(raw_dir, cluster)
        name_by_id = {}
        for query, string_ids in all_ids.items():
            for string_id in string_ids:
                name_by_id[string_id] = query
        # Prefer STRING's display name where available.
        mapping_records = read_json(raw_dir / f"string_{cluster}_mapping.json")
        for record in mapping_records if isinstance(mapping_records, list) else []:
            if isinstance(record, dict) and record.get("stringId") and record.get("preferredName"):
                name_by_id[str(record["stringId"])] = str(record["preferredName"])
        name_to_id = {name: string_id for string_id, name in name_by_id.items()}
        input_genes = gene_sets[cluster]
        mapped_ids = {chosen[gene] for gene in input_genes if gene in chosen}
        audit.update({
            "cluster_id": cluster,
            "n_step7_overlap_genes": len(input_genes),
            "n_mapped_string_ids": len(mapped_ids),
            "n_unmapped_input_genes": len(input_genes - set(chosen)),
            "input_mapped_fraction": len(chosen) / len(input_genes) if input_genes else 0.0,
        })
        mapping_audit.append(audit)
        edges = load_edges(raw_dir, f"string_{cluster}_network_")
        graph = make_graph(edges)
        axis_nodes = mapped_ids & set(graph.nodes())
        sub = graph.subgraph(axis_nodes).copy()
        components = list(nx.connected_components(sub))
        communities = nx.community.louvain_communities(sub, seed=RANDOM_SEED, weight=None) if sub.number_of_edges() else [{node} for node in sorted(axis_nodes)]
        module_map: dict[str, set[str]] = {}
        pathway_recurrence: Counter[str] = Counter()
        for index, community in enumerate(sorted(communities, key=lambda nodes: (-len(nodes), sorted(nodes)[0] if nodes else "")), start=1):
            module_id = f"{cluster}_M{index:03d}"
            for node in community:
                module_map[node] = set(community)
            enrichment_rows, recurrence = parse_module_enrichment(raw_dir, cluster, module_id, name_to_id)
            pathway_recurrence.update(recurrence)
            module_annotations.extend(enrichment_rows)
            module_rows.append({"cluster_id": cluster, "module_id": module_id, "n_nodes": len(community), "n_internal_edges": sub.subgraph(community).number_of_edges(), "module_density": nx.density(sub.subgraph(community)) if len(community) > 1 else 0.0, "nodes": ";".join(sorted(community)), "n_significant_annotations": len(enrichment_rows), "top_annotation": enrichment_rows[0]["description"] if enrichment_rows else ""})
        permutation_edges = []
        candidate_background = sorted(background_node_set)
        for _ in range(N_PERMUTATIONS):
            sampled = stratified_sample(candidate_background, axis_nodes, bins, rng)
            permutation_edges.append(background_graph.subgraph(sampled).number_of_edges())
        observed_edges = sub.number_of_edges()
        expected_edges = sum(permutation_edges) / len(permutation_edges) if permutation_edges else math.nan
        p_empirical = (1 + sum(value >= observed_edges for value in permutation_edges)) / (1 + len(permutation_edges)) if permutation_edges else math.nan
        centrality = centrality_table(graph, axis_nodes, module_map, pathway_recurrence, name_by_id)
        for row in centrality:
            row.update({"cluster_id": cluster, "n_axis_nodes": len(axis_nodes), "n_axis_edges": observed_edges, "top_network_prioritized": centrality.index(row) < 10})
            all_node_rows.append(row)
        for (left, right), score in sorted(edges.items()):
            all_edge_rows.append({"cluster_id": cluster, "string_id_a": left, "string_id_b": right, "combined_score": score, "in_mapped_axis": str(left in mapped_ids and right in mapped_ids)})
        network_summary.append({
            "cluster_id": cluster,
            "n_step7_overlap_genes": len(input_genes),
            "n_mapped_string_ids": len(mapped_ids),
            "n_network_nodes": sub.number_of_nodes(),
            "n_network_edges_score_ge_700": observed_edges,
            "n_connected_components": len(components),
            "n_louvain_modules": len(communities),
            "largest_module_nodes": max((len(nodes) for nodes in communities), default=0),
            "background_nodes": len(background_node_set),
            "background_edges": background_graph.number_of_edges(),
            "degree_stratified_perm_n": len(permutation_edges),
            "perm_expected_edges": expected_edges,
            "observed_to_expected_edge_ratio": observed_edges / expected_edges if expected_edges else math.nan,
            "empirical_p_ge_observed": p_empirical,
            "random_seed": RANDOM_SEED,
            "network_type": "functional",
            "required_score": SCORE_THRESHOLD,
        })

    write_csv(args.output_dir / "t2d_step8c_string_mapping_audit.csv", mapping_audit)
    write_csv(args.output_dir / "t2d_step8c_network_summary.csv", network_summary)
    write_csv(args.output_dir / "t2d_step8c_network_edges.csv", all_edge_rows)
    write_csv(args.output_dir / "t2d_step8c_network_nodes.csv", all_node_rows)
    write_csv(args.output_dir / "t2d_step8c_network_modules.csv", module_rows)
    write_csv(args.output_dir / "t2d_step8c_module_annotations.csv", module_annotations)
    manifest = {
        "status": "complete_network_convergence",
        "tier_a_axes": list(TIER_A),
        "input_rule": "Step 7 overlap genes only for each Tier A axis",
        "species": 9606,
        "string_network_type": "functional",
        "required_score": SCORE_THRESHOLD,
        "add_nodes": 0,
        "background_rule": "frozen union of all 11 Step 7 CTD cluster gene sets used only for degree-stratified empirical randomization",
        "permutation_n": N_PERMUTATIONS,
        "random_seed": RANDOM_SEED,
        "community_method": "NetworkX louvain_communities; weight=None; fixed seed",
        "hub_rule": "network-priority composite averages normalized degree, betweenness, eigenvector, within-module connectivity, and pathway recurrence; module membership reported separately; no degree-only hub claims",
        "outputs": ["t2d_step8c_string_mapping_audit.csv", "t2d_step8c_network_summary.csv", "t2d_step8c_network_edges.csv", "t2d_step8c_network_nodes.csv", "t2d_step8c_network_modules.csv", "t2d_step8c_module_annotations.csv"],
        "interpretation": "network connectivity and module organization are descriptive; no exposure causality, activation, or mediation inference",
    }
    (args.output_dir / "STEP8C_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    report = [
        "# Step 8C — STRING network convergence",
        "",
        "- Status: **complete_network_convergence**",
        "- STRING species: Homo sapiens (9606)",
        "- Network: functional associations, combined score >= 700, no added interactors",
        "- Input: frozen Step 7 overlap genes per Tier A axis",
        "- Randomization background: frozen union of all 11 Step 7 cluster genes",
        "- Permutations: 1,000 degree-stratified samples per axis",
        "- Network-priority score: normalized degree, betweenness, eigenvector, within-module connectivity, and pathway recurrence; descriptive ranking only",
        "",
        "## Axis summary",
        "",
        "| Axis | Input genes | Mapped IDs | Network nodes | Edges | Components | Louvain modules | Observed/expected | Empirical P |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in network_summary:
        report.append(f"| {row['cluster_id']} | {row['n_step7_overlap_genes']} | {row['n_mapped_string_ids']} | {row['n_network_nodes']} | {row['n_network_edges_score_ge_700']} | {row['n_connected_components']} | {row['n_louvain_modules']} | {row['observed_to_expected_edge_ratio']:.3f} | {row['empirical_p_ge_observed']:.4g} |")
    report.extend([
        "",
        "## Interpretation boundary",
        "",
        "Network-prioritized genes are not causal targets.  The network-priority score is a descriptive ranking aid that averages normalized degree, betweenness, eigenvector centrality, within-module connectivity, and pathway recurrence; module membership and annotations remain separate audit fields.  STRING functional associations may include evidence beyond direct physical binding, so results are described as high-confidence functional network convergence rather than definitive direct PPI.",
        "",
        "All axes were analyzed separately.  No transcriptomic data, T2D expression data, or flagship selection entered this stage.",
    ])
    (args.output_dir / "STEP8C_T2D_NETWORK_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
