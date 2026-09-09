from __future__ import annotations

from collections import deque
from pathlib import Path
import json
import math

import networkx as nx
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT
EPC_RUNS = 1000
EPC_THRESHOLD = 0.5
EPC_SEED = 20260909
TOP_K = 25
CONSENSUS_MIN = 7

ALGORITHMS = [
    "MCC",
    "MNC",
    "EPC",
    "Degree",
    "Closeness",
    "Betweenness",
    "Radiality",
    "Stress",
]


def load_network() -> tuple[pd.DataFrame, pd.DataFrame, nx.Graph, list[str]]:
    overlap = pd.read_csv(OUTPUTS / "DINP_CRC_overlap.csv", dtype=str).fillna("")
    hub = pd.read_csv(OUTPUTS / "DINP_CRC_STRING_hub_ranking.csv", dtype=str).fillna("")
    edges = pd.read_csv(OUTPUTS / "DINP_CRC_STRING_network_edges.csv", dtype=str).fillna("")
    overlap["gene_symbol"] = overlap["gene_symbol"].str.strip().str.upper()
    hub["gene_symbol"] = hub["gene_symbol"].str.strip().str.upper()
    hub = hub.drop_duplicates("gene_symbol")
    string_to_gene = dict(zip(hub["string_id"], hub["gene_symbol"]))
    edges["score_numeric"] = pd.to_numeric(edges["score"], errors="coerce")
    edges = edges[edges["score_numeric"].ge(0.7)].copy()
    G = nx.Graph()
    unresolved_edges: list[str] = []
    for _, row in edges.iterrows():
        gene_a = string_to_gene.get(row["stringId_A"], "")
        gene_b = string_to_gene.get(row["stringId_B"], "")
        if not gene_a or not gene_b or gene_a == gene_b:
            unresolved_edges.append(f"{row['stringId_A']}--{row['stringId_B']}")
            continue
        G.add_edge(gene_a, gene_b, weight=float(row["score_numeric"]))
    return overlap, hub, G, unresolved_edges


def component_map(G: nx.Graph) -> dict[str, set[str]]:
    output: dict[str, set[str]] = {}
    for component in nx.connected_components(G):
        component_set = set(component)
        for node in component_set:
            output[node] = component_set
    return output


def mnc_scores(G: nx.Graph) -> dict[str, float]:
    scores: dict[str, float] = {}
    for node in G:
        neighborhood = list(G.neighbors(node))
        if not neighborhood:
            scores[node] = 0.0
            continue
        induced = G.subgraph(neighborhood)
        scores[node] = float(max(len(component) for component in nx.connected_components(induced)))
    return scores


def mcc_scores(G: nx.Graph) -> dict[str, float]:
    scores = {node: 0.0 for node in G}
    for clique in nx.find_cliques(G):
        contribution = float(math.factorial(max(len(clique) - 1, 0)))
        for node in clique:
            scores[node] += contribution
    return scores


def stress_scores(G: nx.Graph) -> dict[str, float]:
    """Count shortest paths through each node, excluding path endpoints."""
    scores = {node: 0.0 for node in G}
    for component in nx.connected_components(G):
        nodes = sorted(component)
        distances: dict[str, dict[str, int]] = {}
        path_counts: dict[str, dict[str, int]] = {}
        for source in nodes:
            distance = {source: 0}
            count = {source: 1}
            queue: deque[str] = deque([source])
            while queue:
                current = queue.popleft()
                for neighbor in G.neighbors(current):
                    if neighbor not in component:
                        continue
                    if neighbor not in distance:
                        distance[neighbor] = distance[current] + 1
                        count[neighbor] = count[current]
                        queue.append(neighbor)
                    elif distance[neighbor] == distance[current] + 1:
                        count[neighbor] += count[current]
            distances[source] = distance
            path_counts[source] = count
        for source_index, source in enumerate(nodes):
            for target in nodes[source_index + 1 :]:
                shortest_distance = distances[source][target]
                for node in nodes:
                    if node in {source, target}:
                        continue
                    if distances[source][node] + distances[node][target] == shortest_distance:
                        scores[node] += path_counts[source][node] * path_counts[node][target]
    return scores


def radiality_scores(G: nx.Graph, components: dict[str, set[str]]) -> dict[str, float]:
    scores: dict[str, float] = {}
    total_nodes = len(G)
    for node in G:
        component = components[node]
        subgraph = G.subgraph(component)
        distances = nx.single_source_shortest_path_length(subgraph, node)
        diameter = nx.diameter(subgraph) if len(component) > 1 else 0
        if diameter == 0:
            scores[node] = 0.0
        else:
            total = sum((diameter + 1 - distance) for distance in distances.values())
            scores[node] = (len(component) / total_nodes) * total / diameter
    return scores


def closeness_scores(G: nx.Graph, components: dict[str, set[str]]) -> dict[str, float]:
    # CytoHubba's connected-component formulation is proportional to the
    # reciprocal-distance sum; this preserves the component-aware ranking.
    scores: dict[str, float] = {}
    for node in G:
        subgraph = G.subgraph(components[node])
        distances = nx.single_source_shortest_path_length(subgraph, node)
        scores[node] = float(sum(1.0 / distance for target, distance in distances.items() if target != node))
    return scores


def epc_scores(G: nx.Graph) -> dict[str, float]:
    """Monte Carlo Edge Percolated Component with CytoHubba's default p=0.5."""
    nodes = list(G.nodes())
    edge_list = list(G.edges())
    scores = {node: 0.0 for node in nodes}
    rng = np.random.default_rng(EPC_SEED)
    for _ in range(EPC_RUNS):
        retained = [edge for edge in edge_list if rng.random() >= EPC_THRESHOLD]
        percolated = nx.Graph()
        percolated.add_nodes_from(nodes)
        percolated.add_edges_from(retained)
        component_sizes: dict[str, int] = {}
        for component in nx.connected_components(percolated):
            size = len(component)
            for node in component:
                component_sizes[node] = size
        for node in nodes:
            scores[node] += component_sizes[node]
    return {node: value / EPC_RUNS for node, value in scores.items()}


def rank_scores(scores: dict[str, float], algorithm: str) -> pd.DataFrame:
    ranked = pd.DataFrame({"gene_symbol": list(scores), "score": list(scores.values())})
    ranked = ranked.sort_values(["score", "gene_symbol"], ascending=[False, True]).reset_index(drop=True)
    ranked.insert(0, "rank", np.arange(1, len(ranked) + 1))
    ranked["algorithm"] = algorithm
    ranked["top25_flag"] = ranked["rank"].le(TOP_K)
    return ranked[["algorithm", "rank", "gene_symbol", "score", "top25_flag"]]


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def main() -> None:
    overlap, hub, G, unresolved_edges = load_network()
    universe = sorted(set(overlap["gene_symbol"]) - {""})
    if len(universe) != 97:
        raise AssertionError(f"Expected 97 overlap genes, got {len(universe)}")
    if G.number_of_edges() != 385:
        raise AssertionError(f"Expected 385 STRING edges at score >=0.7, got {G.number_of_edges()}")
    components = component_map(G)
    largest_component_size = max((len(component) for component in nx.connected_components(G)), default=0)

    scores: dict[str, dict[str, float]] = {}
    components_scores = {
        "MCC": mcc_scores(G),
        "MNC": mnc_scores(G),
        "EPC": epc_scores(G),
        "Degree": {node: float(G.degree(node)) for node in G},
        "Closeness": closeness_scores(G, components),
        "Betweenness": nx.betweenness_centrality(G, normalized=False),
        "Radiality": radiality_scores(G, components),
        "Stress": stress_scores(G),
    }
    scores.update(components_scores)
    long_rows = [rank_scores(scores[algorithm], algorithm) for algorithm in ALGORITHMS]
    rankings_long = pd.concat(long_rows, ignore_index=True)
    rankings_long.to_csv(OUTPUTS / "DINP_CRC_STRING_cytohubba_8_top25_by_algorithm.csv", index=False)

    ranking_wide = pd.DataFrame({"gene_symbol": universe})
    for algorithm in ALGORITHMS:
        current = rankings_long[rankings_long["algorithm"].eq(algorithm)].copy()
        current = current.rename(columns={"rank": f"{algorithm}_rank", "score": f"{algorithm}_score", "top25_flag": f"{algorithm}_top25"})
        ranking_wide = ranking_wide.merge(current[["gene_symbol", f"{algorithm}_rank", f"{algorithm}_score", f"{algorithm}_top25"]], on="gene_symbol", how="left")
    top25_columns = [f"{algorithm}_top25" for algorithm in ALGORITHMS]
    ranking_wide["top25_algorithm_count"] = pd.concat(
        [ranking_wide[column].eq(True) for column in top25_columns], axis=1
    ).sum(axis=1)
    ranking_wide["consensus_hub"] = ranking_wide["top25_algorithm_count"].ge(CONSENSUS_MIN)
    ranking_wide["top25_algorithm_fraction"] = ranking_wide["top25_algorithm_count"] / len(ALGORITHMS)
    ranking_wide["top25_algorithms"] = ranking_wide.apply(
        lambda row: ";".join(algorithm for algorithm in ALGORITHMS if bool(row[f"{algorithm}_top25"])), axis=1
    )

    hub["hub_rank"] = pd.to_numeric(hub["hub_rank"], errors="coerce")
    hub["degree"] = pd.to_numeric(hub["degree"], errors="coerce")
    hub = hub.drop_duplicates("gene_symbol")
    ranking_wide = ranking_wide.merge(
        hub[["gene_symbol", "string_id", "mapping_status", "in_network", "hub_rank", "degree"]],
        on="gene_symbol",
        how="left",
    )
    ranking_wide["original_degree_top20"] = ranking_wide["hub_rank"].le(20)
    ranking_wide["algorithm_panel"] = ",".join(ALGORITHMS)
    ranking_wide = ranking_wide.sort_values(["consensus_hub", "top25_algorithm_count", "gene_symbol"], ascending=[False, False, True]).reset_index(drop=True)
    ranking_wide.insert(0, "consensus_rank", np.arange(1, len(ranking_wide) + 1))

    stable_path = OUTPUTS / "DINP_CRC_101ML_gene_stability.csv"
    cross_path = OUTPUTS / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv"
    stable = pd.read_csv(stable_path, dtype=str).fillna("")
    cross = pd.read_csv(cross_path, dtype=str).fillna("")
    stable["gene_symbol"] = stable["gene_symbol"].str.strip().str.upper()
    cross["gene_symbol"] = cross["gene_symbol"].str.strip().str.upper()
    stable = stable.drop_duplicates("gene_symbol")
    cross = cross.drop_duplicates("gene_symbol")
    ranking_wide = ranking_wide.merge(
        stable[["gene_symbol", "ml_rank", "stable_ml_flag", "ml_priority_score"]], on="gene_symbol", how="left"
    ).merge(
        cross[["gene_symbol", "cross_rank", "candidate_tier", "cross_support_flag", "cross_rank_score"]],
        on="gene_symbol",
        how="left",
    )
    ranking_wide["stable_ml"] = bool_series(ranking_wide["stable_ml_flag"])
    ranking_wide["tier1"] = bool_series(ranking_wide["cross_support_flag"])
    ranking_wide.to_csv(OUTPUTS / "DINP_CRC_STRING_cytohubba_8_algorithm_rankings.csv", index=False)
    consensus = ranking_wide[ranking_wide["consensus_hub"]].copy()
    consensus.to_csv(OUTPUTS / "DINP_CRC_STRING_cytohubba_8_consensus_hubs.csv", index=False)

    comparison_genes = sorted({"CEBPB", "CD36"} | set(cross.loc[cross["cross_support_flag"].astype(str).str.lower().eq("true"), "gene_symbol"]) | set(stable.loc[stable["stable_ml_flag"].astype(str).str.lower().eq("true"), "gene_symbol"]))
    comparison = ranking_wide[ranking_wide["gene_symbol"].isin(comparison_genes)].copy()
    comparison = comparison.sort_values(["tier1", "stable_ml", "consensus_hub", "top25_algorithm_count", "gene_symbol"], ascending=[False, False, False, False, True])
    comparison.to_csv(OUTPUTS / "DINP_CRC_STRING_cytohubba_8_sensitivity_comparison.csv", index=False)

    audit_sets = {
        "consensus_hub_ge7of8": set(consensus["gene_symbol"]),
        "original_degree_top20": set(ranking_wide.loc[ranking_wide["original_degree_top20"], "gene_symbol"]),
        "tier1": set(ranking_wide.loc[ranking_wide["tier1"], "gene_symbol"]),
        "stable_ml": set(ranking_wide.loc[ranking_wide["stable_ml"], "gene_symbol"]),
    }
    audit_rows: list[dict] = []
    for set_name, genes in audit_sets.items():
        audit_rows.append(
            {
                "set_name": set_name,
                "background_n": len(universe),
                "set_n": len(genes),
                "genes": ";".join(sorted(genes)),
            }
        )
    for left_name, right_name in [("consensus_hub_ge7of8", "original_degree_top20"), ("consensus_hub_ge7of8", "tier1"), ("consensus_hub_ge7of8", "stable_ml"), ("tier1", "stable_ml")]:
        left = audit_sets[left_name]
        right = audit_sets[right_name]
        audit_rows.append(
            {
                "set_name": f"{left_name} ∩ {right_name}",
                "background_n": len(universe),
                "set_n": len(left & right),
                "genes": ";".join(sorted(left & right)),
            }
        )
    pd.DataFrame(audit_rows).to_csv(OUTPUTS / "DINP_CRC_STRING_cytohubba_8_overlap_audit.csv", index=False)

    summary_lines = [
        "# DINP–CRC CytoHubba 8-algorithm PPI sensitivity analysis",
        "",
        "## Prespecified sensitivity rule",
        "",
        "The original CRC report explicitly states that genes ranked in the top 25 in at least 7 of 8 CytoHubba algorithms were considered hub genes. Its text does not enumerate the eight algorithms and elsewhere refers to all 12 CytoHubba methods; therefore this analysis records an operationally fixed canonical eight-method panel rather than claiming an exact undocumented panel replication.",
        "",
        f"Panel: {', '.join(ALGORITHMS)}.",
        "",
        f"A gene is a sensitivity consensus hub only when it is in the Top {TOP_K} for at least {CONSENSUS_MIN}/8 algorithms. The threshold and Top25 cutoff were not tuned to the results.",
        "",
        "## Network",
        "",
        f"- STRING network edge threshold: combined score ≥0.700.",
        f"- Edges: {G.number_of_edges()}; mapped/in-network nodes: {G.number_of_nodes()}; largest connected component: {largest_component_size} nodes.",
        f"- Unresolved/self-loop edges skipped: {len(unresolved_edges)}.",
        "- Degree, MNC, MCC, shortest-path centralities, and Stress use the unweighted topology; STRING combined score is retained only as the network-construction filter/edge metadata.",
        f"- EPC: {EPC_RUNS} reproducible Monte Carlo edge-percolation runs, edge-removal threshold {EPC_THRESHOLD}, seed {EPC_SEED}.",
        "",
        "## Results",
        "",
        f"- Consensus hubs (≥7/8 Top25): {len(consensus)}",
        f"- Original degree Top20: {int(ranking_wide['original_degree_top20'].sum())}",
        f"- Consensus hubs retained from original degree Top20: {int((ranking_wide['consensus_hub'] & ranking_wide['original_degree_top20']).sum())}",
        f"- Original Tier 1 genes retained as consensus hubs: {int((ranking_wide['consensus_hub'] & ranking_wide['tier1']).sum())}/{int(ranking_wide['tier1'].sum())}",
        f"- Stable-ML genes that are consensus hubs: {int((ranking_wide['consensus_hub'] & ranking_wide['stable_ml']).sum())}/{int(ranking_wide['stable_ml'].sum())}",
        f"- CEBPB consensus status: {'yes' if bool(ranking_wide.loc[ranking_wide['gene_symbol'].eq('CEBPB'), 'consensus_hub'].iloc[0]) else 'no'}",
        f"- CD36 consensus status: {'yes' if bool(ranking_wide.loc[ranking_wide['gene_symbol'].eq('CD36'), 'consensus_hub'].iloc[0]) else 'no'}",
        "",
        "## Consensus hubs",
        "",
        "| Consensus rank | Gene | Algorithms in Top25 | Count | Original degree rank | Tier 1 | Stable ML |",
        "|---:|---|---|---:|---:|---|---|",
    ]
    for _, row in consensus.sort_values(["top25_algorithm_count", "gene_symbol"], ascending=[False, True]).iterrows():
        degree_rank = "" if pd.isna(row["hub_rank"]) else str(int(row["hub_rank"]))
        summary_lines.append(
            f"| {int(row['consensus_rank'])} | {row['gene_symbol']} | {row['top25_algorithms']} | {int(row['top25_algorithm_count'])} | {degree_rank} | {'yes' if row['tier1'] else 'no'} | {'yes' if row['stable_ml'] else 'no'} |"
        )
    summary_lines.extend(
        [
            "",
            "## Audit interpretation",
            "",
            "The degree Top20 and 8-algorithm consensus are separate hub definitions. The sensitivity analysis does not overwrite the original degree-based PPI result or the original Tier 1 label.",
            "",
            "The individual algorithm rankings and the gene-level comparison table should be reviewed before using consensus status as a downstream prioritization criterion. Consensus hub status is network-topology evidence, not direct DINP binding or causal evidence.",
        ]
    )
    (OUTPUTS / "DINP_CRC_STRING_cytohubba_8_sensitivity_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    log = {
        "network_input": "DINP_CRC_STRING_network_edges.csv filtered to score >= 0.700",
        "network_edge_count": int(G.number_of_edges()),
        "network_node_count": int(G.number_of_nodes()),
        "largest_connected_component": int(largest_component_size),
        "input_background_n": len(universe),
        "algorithm_panel": ALGORITHMS,
        "top_k": TOP_K,
        "consensus_min_algorithms": CONSENSUS_MIN,
        "epc": {"runs": EPC_RUNS, "edge_removal_threshold": EPC_THRESHOLD, "seed": EPC_SEED},
        "unresolved_edges_skipped": unresolved_edges,
        "consensus_hub_n": int(len(consensus)),
        "original_degree_top20_n": int(ranking_wide["original_degree_top20"].sum()),
        "tier1_n": int(ranking_wide["tier1"].sum()),
        "stable_ml_n": int(ranking_wide["stable_ml"].sum()),
        "cebpb_consensus_hub": bool(ranking_wide.loc[ranking_wide["gene_symbol"].eq("CEBPB"), "consensus_hub"].iloc[0]),
        "cd36_consensus_hub": bool(ranking_wide.loc[ranking_wide["gene_symbol"].eq("CD36"), "consensus_hub"].iloc[0]),
        "tier1_consensus_overlap": sorted(ranking_wide.loc[ranking_wide["consensus_hub"] & ranking_wide["tier1"], "gene_symbol"]),
        "stable_ml_consensus_overlap": sorted(ranking_wide.loc[ranking_wide["consensus_hub"] & ranking_wide["stable_ml"], "gene_symbol"]),
        "qc": {
            "ranking_rows": int(len(ranking_wide)),
            "unique_gene_symbols": int(ranking_wide["gene_symbol"].nunique()),
            "long_ranking_rows": int(len(rankings_long)),
            "long_rows_expected": int(len(G) * len(ALGORITHMS)),
            "one_top25_set_per_algorithm": {algorithm: int(rankings_long.loc[rankings_long["algorithm"].eq(algorithm), "top25_flag"].sum()) for algorithm in ALGORITHMS},
        },
    }
    (OUTPUTS / "DINP_CRC_STRING_cytohubba_8_sensitivity_log.md").write_text(
        "# DINP–CRC CytoHubba 8-algorithm sensitivity audit log\n\n```json\n"
        + json.dumps(log, indent=2, ensure_ascii=False)
        + "\n```\n",
        encoding="utf-8",
    )

    print(f"network_nodes={G.number_of_nodes()} network_edges={G.number_of_edges()}")
    print(f"consensus_hubs_ge7of8={len(consensus)}")
    print("consensus genes:", ", ".join(consensus.sort_values(["top25_algorithm_count", "gene_symbol"], ascending=[False, True])["gene_symbol"]))
    print("CEBPB/CD36:")
    print(comparison[comparison["gene_symbol"].isin(["CEBPB", "CD36"])][["gene_symbol", "top25_algorithm_count", "top25_algorithms", "consensus_hub", "original_degree_top20", "tier1", "stable_ml"]].to_string(index=False))


if __name__ == "__main__":
    main()
