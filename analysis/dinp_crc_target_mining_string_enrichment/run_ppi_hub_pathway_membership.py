from __future__ import annotations

import ast
import math
from pathlib import Path

import networkx as nx
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"


def parse_list(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    try:
        parsed = ast.literal_eval(text)
    except (SyntaxError, ValueError):
        return [x.strip() for x in text.split(";") if x.strip()]
    if isinstance(parsed, (list, tuple)):
        return [str(x).strip() for x in parsed if str(x).strip()]
    return [str(parsed).strip()]


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def main() -> None:
    mapping = pd.read_csv(OUT / "DINP_CRC_STRING_mapping_selected.csv", dtype=str).fillna("")
    edges = pd.read_csv(OUT / "DINP_CRC_STRING_network_edges.csv", dtype=str).fillna("")
    overlap = pd.read_csv(OUT / "DINP_CRC_overlap.csv", dtype=str).fillna("")
    reps = pd.read_csv(OUT / "DINP_CRC_GO_KEGG_representatives.csv", dtype=str).fillna("")

    mapping["string_id"] = mapping["string_id"].astype(str).str.strip()
    mapping["gene_symbol"] = mapping["gene_symbol"].astype(str).str.strip().str.upper()
    mapping = mapping.drop_duplicates("gene_symbol", keep="first")
    mapped = mapping[mapping["string_id"].ne("")].drop_duplicates("string_id", keep="first").copy()
    id_to_gene = dict(zip(mapped["string_id"], mapped["gene_symbol"]))
    id_to_preferred = dict(zip(mapped["string_id"], mapped["preferred_name"]))

    graph = nx.Graph()
    graph.add_nodes_from(mapped["string_id"].tolist())
    for row in edges.itertuples(index=False):
        a = str(row.stringId_A).strip()
        b = str(row.stringId_B).strip()
        if not a or not b or a == b:
            continue
        score = safe_float(row.score)
        if graph.has_edge(a, b):
            graph[a][b]["score"] = max(score, safe_float(graph[a][b].get("score", 0)))
        else:
            graph.add_edge(a, b, score=score)

    weighted_degree = dict(graph.degree(weight="score"))
    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph, normalized=True, weight=None)
    closeness = nx.closeness_centrality(graph)
    try:
        eigenvector = nx.eigenvector_centrality(graph, max_iter=5000, weight="score")
    except nx.NetworkXException:
        eigenvector = {node: 0.0 for node in graph.nodes}
    pagerank = nx.pagerank(graph, weight="score") if graph.number_of_edges() else {node: 0.0 for node in graph.nodes}

    hub = mapping[["gene_symbol", "string_id", "preferred_name", "mapping_status"]].copy()
    hub["in_network"] = hub["string_id"].isin(graph.nodes).map({True: "yes", False: "no"})
    hub["degree"] = hub["string_id"].map(degree).fillna(0).astype(int)
    hub["weighted_degree"] = hub["string_id"].map(weighted_degree).fillna(0.0)
    hub["betweenness_centrality"] = hub["string_id"].map(betweenness).fillna(0.0)
    hub["closeness_centrality"] = hub["string_id"].map(closeness).fillna(0.0)
    hub["eigenvector_centrality"] = hub["string_id"].map(eigenvector).fillna(0.0)
    hub["pagerank"] = hub["string_id"].map(pagerank).fillna(0.0)
    hub = hub.sort_values(
        ["degree", "weighted_degree", "betweenness_centrality", "eigenvector_centrality", "gene_symbol"],
        ascending=[False, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    hub["hub_rank"] = pd.Series(pd.NA, index=hub.index, dtype="Int64")
    mapped_rows = hub["string_id"].ne("")
    hub.loc[mapped_rows, "hub_rank"] = range(1, int(mapped_rows.sum()) + 1)
    hub["hub_tier"] = hub["hub_rank"].map(lambda x: "top20" if pd.notna(x) and x <= 20 else ("ranked_all" if pd.notna(x) else "unmapped"))
    hub_columns = [
        "hub_rank",
        "hub_tier",
        "gene_symbol",
        "string_id",
        "preferred_name",
        "mapping_status",
        "in_network",
        "degree",
        "weighted_degree",
        "betweenness_centrality",
        "closeness_centrality",
        "eigenvector_centrality",
        "pagerank",
    ]
    hub[hub_columns].to_csv(OUT / "DINP_CRC_STRING_hub_ranking.csv", index=False, encoding="utf-8-sig")
    hub.loc[hub["hub_rank"] <= 20, hub_columns].to_csv(
        OUT / "DINP_CRC_STRING_hubs_top20.csv", index=False, encoding="utf-8-sig"
    )

    overlap_genes = set(overlap["gene_symbol"].astype(str).str.strip().str.upper())
    reps = reps.copy()
    if "is_representative" in reps.columns:
        reps = reps[reps["is_representative"].astype(str).isin({"1", "1.0", "True", "true"})].copy()
    if "fdr_numeric" in reps.columns:
        reps["fdr_numeric"] = pd.to_numeric(reps["fdr_numeric"], errors="coerce")
        reps = reps[reps["fdr_numeric"].lt(0.05)].copy()

    hub_lookup = hub.set_index("gene_symbol").to_dict("index")
    membership_rows: list[dict[str, object]] = []
    for row in reps.itertuples(index=False):
        ids = parse_list(getattr(row, "inputGenes", ""))
        names = parse_list(getattr(row, "preferredNames", ""))
        if len(ids) != len(names):
            names = [id_to_gene.get(item, id_to_preferred.get(item, "")) for item in ids]
        for string_id, gene_symbol in zip(ids, names):
            gene_symbol = str(gene_symbol).strip().upper()
            if not gene_symbol or gene_symbol not in overlap_genes:
                continue
            metrics = hub_lookup.get(gene_symbol, {})
            membership_rows.append(
                {
                    "analysis_category": getattr(row, "analysis_category", ""),
                    "category": getattr(row, "category", ""),
                    "term": getattr(row, "term", ""),
                    "description": getattr(row, "description", ""),
                    "reduced_rank": getattr(row, "reduced_rank", ""),
                    "fdr": getattr(row, "fdr", ""),
                    "fdr_numeric": getattr(row, "fdr_numeric", ""),
                    "p_value": getattr(row, "p_value", ""),
                    "gene_symbol": gene_symbol,
                    "string_id": string_id,
                    "hub_rank": metrics.get("hub_rank", ""),
                    "degree": metrics.get("degree", 0),
                    "membership_basis": "significant nonredundant STRING representative terms",
                }
            )

    membership = pd.DataFrame(membership_rows)
    membership = membership.sort_values(
        ["analysis_category", "fdr_numeric", "hub_rank", "gene_symbol", "term"],
        ascending=[True, True, True, True, True],
        kind="mergesort",
    )
    membership.to_csv(OUT / "DINP_CRC_pathway_membership.csv", index=False, encoding="utf-8-sig")

    top20_genes = set(hub.loc[hub["hub_rank"] <= 20, "gene_symbol"])
    hub_membership = membership[membership["gene_symbol"].isin(top20_genes)].copy()
    hub_membership.to_csv(OUT / "DINP_CRC_hub_pathway_membership.csv", index=False, encoding="utf-8-sig")

    summary_rows: list[dict[str, object]] = []
    for row in hub.loc[hub["hub_rank"] <= 20].itertuples(index=False):
        current = membership[membership["gene_symbol"].eq(row.gene_symbol)].copy()
        summary_rows.append(
            {
                "hub_rank": row.hub_rank,
                "gene_symbol": row.gene_symbol,
                "degree": row.degree,
                "weighted_degree": row.weighted_degree,
                "betweenness_centrality": row.betweenness_centrality,
                "pathway_or_term_membership_count": len(current),
                "GO_BP_membership_count": int(current["analysis_category"].eq("GO_BP").sum()),
                "GO_MF_membership_count": int(current["analysis_category"].eq("GO_MF").sum()),
                "GO_CC_membership_count": int(current["analysis_category"].eq("GO_CC").sum()),
                "KEGG_membership_count": int(current["analysis_category"].eq("KEGG").sum()),
                "member_terms": "; ".join(
                    f"{term} | {desc}" for term, desc in zip(current["term"], current["description"])
                ),
            }
        )
    pd.DataFrame(summary_rows).to_csv(OUT / "DINP_CRC_hub_pathway_summary.csv", index=False, encoding="utf-8-sig")

    nonzero = int((hub["degree"] > 0).sum())
    unmapped_count = int((hub["mapping_status"].str.lower() != "mapped").sum())
    generated_at = pd.Timestamp.now(tz="Asia/Shanghai").isoformat()
    log_lines = [
        "# DINP–CRC PPI hub and pathway membership audit log",
        "",
        f"- Generated locally: {generated_at}",
        "- Input overlap: `DINP_CRC_overlap.csv` (97 genes).",
        "- PPI network: `DINP_CRC_STRING_network_edges.csv`, STRING functional edges at combined score >=700.",
        "- Hub ranking is continuous and deterministic; it is sorted by degree, weighted degree, betweenness, eigenvector centrality, then gene symbol.",
        "- `top20` is a reporting convenience only; it is not a statistical cutoff.",
        "- Pathway membership uses significant (`FDR < 0.05`) nonredundant STRING representative terms from `DINP_CRC_GO_KEGG_representatives.csv`.",
        "- Membership is reported as gene-to-term evidence; it does not imply DINP-specific causal regulation.",
        "",
        "## Count reconciliation",
        "",
        f"- Input genes represented in hub table: {len(hub)}",
        f"- Mapped STRING proteins ranked: {len(hub) - unmapped_count}",
        f"- Input genes not mapped to STRING: {unmapped_count}",
        f"- Nodes with >=1 high-confidence edge: {nonzero}",
        f"- Network edges used: {graph.number_of_edges()}",
        f"- Significant representative terms: {len(reps)}",
        f"- Gene–term membership rows: {len(membership)}",
        f"- Top-20 hub–term membership rows: {len(hub_membership)}",
        "",
        "## Output files",
        "",
        "- `DINP_CRC_STRING_hub_ranking.csv`: all mapped STRING proteins with degree and centrality metrics.",
        "- `DINP_CRC_STRING_hubs_top20.csv`: top 20 ranked hubs for focused review.",
        "- `DINP_CRC_pathway_membership.csv`: all input-gene memberships in significant representative GO/KEGG terms.",
        "- `DINP_CRC_hub_pathway_membership.csv`: the same membership table restricted to top-20 hubs.",
        "- `DINP_CRC_hub_pathway_summary.csv`: one row per top-20 hub with term counts and member terms.",
    ]
    (OUT / "DINP_CRC_PPI_hub_pathway_log.md").write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
