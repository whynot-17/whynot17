#!/usr/bin/env python
"""PPI/topology and target-actionability prioritization for seven macrophage genes.

The seven inputs are fixed by the preceding donor-aware macrophage driver
decomposition.  STRING functional and physical networks are queried without
adding interactors.  UniProt, PDBe and ChEMBL are used only for protein,
structure and measured-activity context; their evidence is kept separate from
the STRING network evidence.

The output is a transparent shortlist, not a causal target claim.  NEAT1 is
kept as a non-protein state/regulatory node and is excluded from docking-ready
protein candidates.  The proposed docking shortlist is role-based: two
high-confidence network bridges plus the most actionable direct prostaglandin
node.  PTGES3 remains an explicit pathway-direct reserve candidate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "analysis" / "dinp_crc_81gene_macrophage_driver_decomposition" / "outputs"
OUT_DIR = Path(__file__).resolve().parent / "outputs"
CANDIDATE_INPUT = INPUT_DIR / "macrophage_driver_candidates.csv"
PATHWAY_INPUT = INPUT_DIR / "pathway_membership_audit.csv"
DEFAULT_CANDIDATES = ["NEAT1", "MMP9", "TIMP1", "STAT3", "PTGER4", "PTGES3", "CXCR4"]
SPECIES = 9606
STRING_SCORE = 0.700
STRING_BASE = "https://string-db.org/api/tsv"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
PDBE_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/best_structures"
CHEMBL_TARGET_URL = "https://www.ebi.ac.uk/chembl/api/data/target.json"
CHEMBL_ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
GP_RELEASE = "e114_eg62_p19_27110d83"
RANDOM_SEED = 20260903
# Only these exact pathway records constitute direct prostaglandin/AA
# evidence.  The broad GO inflammatory-response annotation is retained as
# context, but must not promote every inflammatory gene to a direct pathway
# node.
DIRECT_PROSTAGLANDIN_TERM_IDS = {"GO:0006693", "GO:0001516", "KEGG:00590"}
EXPLICIT_DIRECT_PROSTAGLANDIN_GENES = {"PTGER4", "PTGES3"}


def request_get(session: requests.Session, url: str, **kwargs) -> requests.Response:
    kwargs.setdefault("timeout", 60)
    response = session.get(url, **kwargs)
    response.raise_for_status()
    return response


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_candidates() -> list[str]:
    if CANDIDATE_INPUT.exists():
        frame = pd.read_csv(CANDIDATE_INPUT)
        if "gene_symbol" in frame.columns:
            values = frame["gene_symbol"].dropna().astype(str).str.upper().tolist()
            if values:
                return list(dict.fromkeys(values))
    return DEFAULT_CANDIDATES


def read_pathway_membership(genes: list[str]) -> pd.DataFrame:
    if not PATHWAY_INPUT.exists():
        return pd.DataFrame(columns=["gene_symbol", "term_id", "source", "term_name", "pathway_category"])
    frame = pd.read_csv(PATHWAY_INPUT)
    return frame.loc[frame["gene_symbol"].astype(str).str.upper().isin(genes)].copy()


def query_string(session: requests.Session, genes: list[str], network_type: str) -> pd.DataFrame:
    identifiers = "\r".join(genes)
    url = f"{STRING_BASE}/network"
    response = request_get(
        session,
        url,
        params={
            "identifiers": identifiers,
            "species": SPECIES,
            "required_score": int(STRING_SCORE * 1000),
            "network_type": network_type,
        },
    )
    if not response.text.strip() or response.text.startswith("stringId_A") is False:
        return pd.DataFrame()
    return pd.read_csv(StringIO(response.text), sep="\t")


def query_string_mapping(session: requests.Session, genes: list[str]) -> pd.DataFrame:
    response = request_get(
        session,
        f"{STRING_BASE}/get_string_ids",
        params={
            "identifiers": "\r".join(genes),
            "species": SPECIES,
            "limit": 5,
        },
    )
    if not response.text.strip() or not response.text.startswith("queryIndex"):
        return pd.DataFrame()
    return pd.read_csv(StringIO(response.text), sep="\t")


def query_uniprot(session: requests.Session, gene: str) -> dict:
    response = request_get(
        session,
        UNIPROT_URL,
        params={
            "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
            "format": "tsv",
            "fields": "accession,id,protein_name,length,gene_names",
            "size": 5,
        },
    )
    lines = [line for line in response.text.splitlines() if line.strip()]
    if len(lines) < 2:
        return {"gene_symbol": gene, "uniprot_accession": "", "uniprot_entry": "", "protein_name": "", "protein_length": np.nan, "gene_names": "", "uniprot_reviewed": False}
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    record = dict(zip(header, row))
    return {
        "gene_symbol": gene,
        "uniprot_accession": record.get("Entry", ""),
        "uniprot_entry": record.get("Entry Name", ""),
        "protein_name": record.get("Protein names", ""),
        "protein_length": int(record["Length"]) if record.get("Length", "").isdigit() else np.nan,
        "gene_names": record.get("Gene Names", ""),
        "uniprot_reviewed": True,
    }


def query_pdbe(session: requests.Session, accession: str) -> list[dict]:
    if not accession:
        return []
    response = request_get(session, f"{PDBE_URL}/{accession}")
    payload = response.json()
    return payload.get(accession, []) if isinstance(payload, dict) else []


def query_chembl(session: requests.Session, accession: str) -> tuple[list[dict], int, int]:
    if not accession:
        return [], 0, 0
    response = request_get(
        session,
        CHEMBL_TARGET_URL,
        params={"target_components__accession": accession, "limit": 100},
    )
    payload = response.json()
    targets = payload.get("targets", []) if isinstance(payload, dict) else []
    total_count = int(payload.get("page_meta", {}).get("total_count", len(targets)))
    single = [row for row in targets if str(row.get("target_type", "")).upper() == "SINGLE PROTEIN"]
    activity_counts: list[int] = []
    for row in single[:3]:
        target_id = row.get("target_chembl_id")
        if not target_id:
            continue
        activity_response = request_get(
            session,
            CHEMBL_ACTIVITY_URL,
            params={"target_chembl_id": target_id, "limit": 1},
        )
        activity_payload = activity_response.json()
        count = int(activity_payload.get("page_meta", {}).get("total_count", 0))
        row["activity_count"] = count
        activity_counts.append(count)
    return targets, total_count, max(activity_counts) if activity_counts else 0


def as_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def parse_edges(frame: pd.DataFrame, genes: set[str]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=["gene_a", "gene_b", "score", "string_id_a", "string_id_b"])
    rows: list[dict] = []
    for _, row in frame.iterrows():
        a = str(row.get("preferredName_A", "")).upper()
        b = str(row.get("preferredName_B", "")).upper()
        if a not in genes or b not in genes or a == b:
            continue
        left, right = sorted((a, b))
        rows.append(
            {
                "gene_a": left,
                "gene_b": right,
                "score": as_float(row.get("score")),
                "string_id_a": str(row.get("stringId_A", "")),
                "string_id_b": str(row.get("stringId_B", "")),
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=["gene_a", "gene_b", "score", "string_id_a", "string_id_b"])
    return result.groupby(["gene_a", "gene_b"], as_index=False).agg(
        score=("score", "max"), string_id_a=("string_id_a", "first"), string_id_b=("string_id_b", "first")
    )


def minmax(series: pd.Series) -> pd.Series:
    values = series.astype(float)
    if values.max() == values.min():
        return pd.Series(0.0, index=series.index)
    return (values - values.min()) / (values.max() - values.min())


def graph_topology(genes: list[str], edges: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    graph = nx.Graph()
    graph.add_nodes_from(genes)
    for _, row in edges.iterrows():
        graph.add_edge(row["gene_a"], row["gene_b"], score=float(row["score"]))
    degree = dict(graph.degree())
    weighted_degree = dict(graph.degree(weight="score"))
    betweenness = nx.betweenness_centrality(graph, normalized=True, weight=None)
    closeness = nx.closeness_centrality(graph)
    try:
        eigenvector = nx.eigenvector_centrality(graph, max_iter=2000, weight="score")
    except nx.NetworkXException:
        eigenvector = {gene: 0.0 for gene in genes}
    if len(graph) and hasattr(nx.community, "louvain_communities"):
        communities = nx.community.louvain_communities(graph, weight="score", seed=RANDOM_SEED)
    else:
        communities = list(nx.community.greedy_modularity_communities(graph, weight="score")) if graph.number_of_edges() else [{gene} for gene in genes]
    module_by_gene = {}
    for module_id, members in enumerate(sorted(communities, key=lambda x: (min(x), len(x))), start=1):
        for gene in members:
            module_by_gene[gene] = f"M{module_id:02d}"
    node = pd.DataFrame({
        "gene_symbol": genes,
        "functional_degree": [degree.get(g, 0) for g in genes],
        "functional_weighted_degree": [weighted_degree.get(g, 0.0) for g in genes],
        "functional_betweenness": [betweenness.get(g, 0.0) for g in genes],
        "functional_closeness": [closeness.get(g, 0.0) for g in genes],
        "functional_eigenvector": [eigenvector.get(g, 0.0) for g in genes],
        "functional_module": [module_by_gene.get(g, "") for g in genes],
    })
    node["topology_score"] = (
        minmax(node["functional_degree"])
        + minmax(node["functional_weighted_degree"])
        + minmax(node["functional_betweenness"])
        + minmax(node["functional_eigenvector"])
    ) / 4.0
    modules = pd.DataFrame(
        [{"functional_module": module_id, "genes": ";".join(sorted(members)), "n_genes": len(members)} for module_id, members in sorted(((module_by_gene[g], {x for x, m in module_by_gene.items() if m == module_by_gene[g]}) for g in genes), key=lambda x: x[0])]
    ).drop_duplicates("functional_module")
    summary = {
        "n_input_nodes": len(genes),
        "n_edges": int(edges.shape[0]),
        "possible_edges": int(len(genes) * (len(genes) - 1) / 2),
        "edge_density": float(nx.density(graph)) if len(graph) > 1 else 0.0,
        "n_connected_components": int(nx.number_connected_components(graph)) if len(graph) else 0,
        "largest_component_nodes": int(max((len(c) for c in nx.connected_components(graph)), default=0)),
        "n_louvain_modules": int(len(communities)),
    }
    return node, modules, summary


def pathway_features(genes: list[str], pathway: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in genes:
        subset = pathway.loc[pathway["gene_symbol"].astype(str).str.upper().eq(gene)]
        direct_terms = subset.loc[subset["term_id"].astype(str).isin(DIRECT_PROSTAGLANDIN_TERM_IDS)]
        inflammation_terms = subset.loc[subset["term_id"].astype(str).eq("GO:0006954")]
        cats = sorted(set(subset["pathway_category"].dropna().astype(str)))
        terms = sorted(set(subset["term_id"].dropna().astype(str)))
        if gene in EXPLICIT_DIRECT_PROSTAGLANDIN_GENES:
            relation_basis = "explicit direct prostaglandin node"
        elif not direct_terms.empty:
            relation_basis = "exact prostaglandin/AA pathway term"
        elif not inflammation_terms.empty:
            relation_basis = "inflammatory-context only"
        else:
            relation_basis = "no prior pathway relation"
        rows.append(
            {
                "gene_symbol": gene,
                "prior_pathway_categories": ";".join(cats),
                "prior_pathway_term_ids": ";".join(terms),
                "prior_pathway_term_count": int(len(terms)),
                "core_prostaglandin_AA_inflammatory_term_count": int(direct_terms["term_id"].nunique()),
                "inflammatory_context_term_count": int(inflammation_terms["term_id"].nunique()),
                "direct_prostaglandin_pathway_node": bool(not direct_terms.empty or gene in EXPLICIT_DIRECT_PROSTAGLANDIN_GENES),
                "pathway_relation_basis": relation_basis,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", nargs="*", default=None)
    args = parser.parse_args()
    genes = [str(g).upper() for g in (args.candidates or read_candidates())]
    genes = list(dict.fromkeys(genes))
    if genes != DEFAULT_CANDIDATES:
        # The preceding output is the authoritative candidate family; fail
        # loudly if a later run silently changes it.
        if set(genes) != set(DEFAULT_CANDIDATES):
            raise ValueError(f"Expected the frozen seven-gene family, found {genes}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pathway = read_pathway_membership(genes)
    pathway.to_csv(OUT_DIR / "pathway_relation_audit.csv", index=False)

    session = requests.Session()
    session.headers.update({"User-Agent": "whynot17-dinp-crc-network-prioritization/1.0"})
    mapping = query_string_mapping(session, genes)
    functional_raw = query_string(session, genes, "functional")
    physical_raw = query_string(session, genes, "physical")
    mapping.to_csv(OUT_DIR / "string_mapping.tsv", sep="\t", index=False)
    functional_raw.to_csv(OUT_DIR / "string_functional_network.tsv", sep="\t", index=False)
    physical_raw.to_csv(OUT_DIR / "string_physical_network.tsv", sep="\t", index=False)
    functional_edges = parse_edges(functional_raw, set(genes))
    physical_edges = parse_edges(physical_raw, set(genes))
    functional_edges.to_csv(OUT_DIR / "functional_edges.csv", index=False)
    physical_edges.to_csv(OUT_DIR / "physical_edges.csv", index=False)

    uniprot_rows = []
    pdb_rows = []
    chembl_rows = []
    for gene in genes:
        info = query_uniprot(session, gene)
        info["is_protein_coding"] = bool(info.get("uniprot_accession"))
        uniprot_rows.append(info)
        structures = query_pdbe(session, info.get("uniprot_accession", ""))
        for structure in structures:
            pdb_rows.append(
                {
                    "gene_symbol": gene,
                    "uniprot_accession": info.get("uniprot_accession", ""),
                    "pdb_id": structure.get("pdb_id", ""),
                    "chain_id": structure.get("chain_id", ""),
                    "experimental_method": structure.get("experimental_method", ""),
                    "resolution": structure.get("resolution"),
                    "coverage": structure.get("coverage"),
                    "unp_start": structure.get("unp_start"),
                    "unp_end": structure.get("unp_end"),
                }
            )
        targets, target_count, activity_count = query_chembl(session, info.get("uniprot_accession", ""))
        for target in targets:
            chembl_rows.append(
                {
                    "gene_symbol": gene,
                    "uniprot_accession": info.get("uniprot_accession", ""),
                    "target_chembl_id": target.get("target_chembl_id", ""),
                    "pref_name": target.get("pref_name", ""),
                    "target_type": target.get("target_type", ""),
                    "activity_count": target.get("activity_count", np.nan),
                }
            )
        info["chembl_target_record_count"] = target_count
        info["chembl_single_protein_target_count"] = int(sum(str(t.get("target_type", "")).upper() == "SINGLE PROTEIN" for t in targets))
        info["chembl_max_single_protein_activity_count"] = activity_count

    uniprot_df = pd.DataFrame(uniprot_rows)
    pdb_df = pd.DataFrame(pdb_rows)
    chembl_df = pd.DataFrame(chembl_rows)
    uniprot_df.to_csv(OUT_DIR / "uniprot_target_annotations.csv", index=False)
    pdb_df.to_csv(OUT_DIR / "pdbe_best_structures.csv", index=False)
    chembl_df.to_csv(OUT_DIR / "chembl_target_annotations.csv", index=False)

    if pdb_df.empty:
        pdb_summary = pd.DataFrame({"gene_symbol": genes, "pdb_structure_count": 0, "best_pdb_id": "", "best_structure_method": "", "best_resolution": np.nan, "best_coverage": np.nan})
    else:
        rows = []
        for gene in genes:
            subset = pdb_df.loc[pdb_df["gene_symbol"].eq(gene)].copy()
            if subset.empty:
                rows.append({"gene_symbol": gene, "pdb_structure_count": 0, "best_pdb_id": "", "best_structure_method": "", "best_resolution": np.nan, "best_coverage": np.nan})
                continue
            subset["coverage_num"] = pd.to_numeric(subset["coverage"], errors="coerce").fillna(0.0)
            subset["resolution_num"] = pd.to_numeric(subset["resolution"], errors="coerce").fillna(np.inf)
            best = subset.sort_values(["coverage_num", "resolution_num"], ascending=[False, True]).iloc[0]
            rows.append({"gene_symbol": gene, "pdb_structure_count": int(len(subset)), "best_pdb_id": best["pdb_id"], "best_structure_method": best["experimental_method"], "best_resolution": best["resolution"], "best_coverage": best["coverage"]})
        pdb_summary = pd.DataFrame(rows)

    topology, modules, network_summary = graph_topology(genes, functional_edges)
    physical_degree = physical_edges.groupby("gene_a").size().add(physical_edges.groupby("gene_b").size(), fill_value=0) if not physical_edges.empty else pd.Series(dtype=float)
    topology["physical_degree"] = topology["gene_symbol"].map(physical_degree).fillna(0).astype(int)
    topology["physical_edge_present"] = topology["physical_degree"].gt(0)
    features = pathway_features(genes, pathway)
    node = topology.merge(features, on="gene_symbol", how="left").merge(uniprot_df, on="gene_symbol", how="left").merge(pdb_summary, on="gene_symbol", how="left")
    node["docking_protein_eligible"] = node["is_protein_coding"].fillna(False) & node["pdb_structure_count"].fillna(0).gt(0)
    node["known_measured_target_evidence"] = node["chembl_max_single_protein_activity_count"].fillna(0).gt(0)
    node["candidate_role"] = "supporting/context"
    node.loc[node["topology_score"].ge(node["topology_score"].nlargest(2).min()), "candidate_role"] = "network bridge candidate"
    node.loc[node["direct_prostaglandin_pathway_node"].fillna(False), "candidate_role"] = "direct prostaglandin node"
    node.loc[node["gene_symbol"].eq("NEAT1"), "candidate_role"] = "non-protein state/regulatory node"
    node["network_priority_rank"] = node["topology_score"].rank(method="min", ascending=False).astype(int)
    node = node.sort_values(["candidate_role", "network_priority_rank", "gene_symbol"]).reset_index(drop=True)
    node.to_csv(OUT_DIR / "network_target_evidence_matrix.csv", index=False)
    modules.to_csv(OUT_DIR / "functional_network_modules.csv", index=False)

    # Role-based shortlist: the two highest-topology mapped proteins outside
    # the direct prostaglandin class plus the most actionable direct
    # prostaglandin node. PTGES3 is retained as reserve.  Importantly,
    # inflammatory-context membership alone does not remove a node from the
    # network-bridge pool.
    protein_nodes = node.loc[node["docking_protein_eligible"]].copy()
    bridge = protein_nodes.loc[~protein_nodes["direct_prostaglandin_pathway_node"].fillna(False)].sort_values(["topology_score", "functional_weighted_degree", "gene_symbol"], ascending=[False, False, True])
    bridge_genes = bridge["gene_symbol"].head(2).tolist()
    direct = protein_nodes.loc[protein_nodes["direct_prostaglandin_pathway_node"].fillna(False)].copy()
    direct["actionability_proxy"] = minmax(direct["chembl_max_single_protein_activity_count"].fillna(0)) + minmax(direct["best_coverage"].fillna(0))
    direct_genes = direct.sort_values(["actionability_proxy", "chembl_max_single_protein_activity_count", "gene_symbol"], ascending=[False, False, True])["gene_symbol"].head(1).tolist()
    shortlist_genes = list(dict.fromkeys(bridge_genes + direct_genes))
    shortlist = node.loc[node["gene_symbol"].isin(shortlist_genes)].copy()
    shortlist["shortlist_role"] = np.where(shortlist["gene_symbol"].isin(bridge_genes), "network bridge", "direct prostaglandin node")
    shortlist["shortlist_rank"] = shortlist["gene_symbol"].map({gene: idx for idx, gene in enumerate(shortlist_genes, start=1)})
    shortlist = shortlist.sort_values("shortlist_rank")
    shortlist.to_csv(OUT_DIR / "docking_md_shortlist.csv", index=False)
    reserve = node.loc[node["gene_symbol"].eq("PTGES3")].copy()
    reserve["reserve_reason"] = "direct prostaglandin synthase with experimental structure; isolated in seven-node STRING subnetwork and limited ChEMBL activity count"
    reserve.to_csv(OUT_DIR / "docking_md_reserve_candidates.csv", index=False)

    physical_summary = {
        "physical_n_edges": int(physical_edges.shape[0]),
        "physical_edge_density": float(physical_edges.shape[0] / network_summary["possible_edges"]) if network_summary["possible_edges"] else 0.0,
        "physical_edges": [f"{r.gene_a}--{r.gene_b}" for r in physical_edges.itertuples()] if not physical_edges.empty else [],
    }
    network_summary.update({
        "string_species": SPECIES,
        "string_required_score": STRING_SCORE,
        "string_functional_edges": int(functional_edges.shape[0]),
        "string_physical_edges": int(physical_edges.shape[0]),
        "mapped_string_input_genes": int(mapping["preferredName"].astype(str).str.upper().isin(genes).sum()) if not mapping.empty and "preferredName" in mapping.columns else 0,
        "protein_coding_nodes": int(node["is_protein_coding"].fillna(False).sum()),
        "structure_eligible_nodes": int(node["docking_protein_eligible"].sum()),
        "docking_md_shortlist": ";".join(shortlist_genes),
        **physical_summary,
    })
    pd.DataFrame([network_summary]).to_csv(OUT_DIR / "network_summary.csv", index=False)

    request_audit = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "genes": genes,
        "string": {"mapping_endpoint": f"{STRING_BASE}/get_string_ids", "network_endpoint": f"{STRING_BASE}/network", "species": SPECIES, "required_score": STRING_SCORE, "functional_network_type": "functional", "physical_network_type": "physical", "added_interactors": 0},
        "uniprot": {"endpoint": UNIPROT_URL, "reviewed_only": True, "organism_id": 9606},
        "pdbe": {"endpoint": PDBE_URL, "selection": "all best_structures records retained; best structure selected by maximum coverage then minimum resolution"},
        "chembl": {"target_endpoint": CHEMBL_TARGET_URL, "activity_endpoint": CHEMBL_ACTIVITY_URL, "target_rule": "target_components accession exact; SINGLE PROTEIN records summarized separately"},
        "gprofiler_pathway_membership": {"release": GP_RELEASE, "input": str(PATHWAY_INPUT), "term_scope": "previously observed exact terms only"},
    }
    (OUT_DIR / "api_request_audit.json").write_text(json.dumps(request_audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    manifest = {
        "analysis": "Macrophage seven-gene PPI, topology and target-actionability prioritization",
        "run_timestamp_utc": request_audit["run_timestamp_utc"],
        "candidate_input": str(CANDIDATE_INPUT),
        "candidate_genes": genes,
        "candidate_family_frozen": True,
        "string_required_score": STRING_SCORE,
        "string_functional_network": "functional associations; no added interactors",
        "string_physical_network": "physical associations reported separately; no added interactors",
        "network_inference": "descriptive induced subgraph topology; no causal direction and no network-level P-value claimed",
        "pathway_input": str(PATHWAY_INPUT),
        "pathway_membership_release": GP_RELEASE,
        "pathway_membership_rule": "exact memberships from the previously observed g:Profiler terms; no term selected using macrophage results",
        "target_context_sources": ["UniProt reviewed human entries", "PDBe UniProt best structures", "ChEMBL target/activity records"],
        "docking_md_shortlist_rule": "two highest-topology protein candidates outside the direct prostaglandin class plus the most actionable direct prostaglandin node; PTGES3 retained as reserve",
        "docking_md_shortlist": shortlist_genes,
        "interpretation_boundary": "shortlist candidates are not causal targets; structure and measured ligand records do not establish DINP binding",
        "network_summary": network_summary,
        "outputs": [str(p) for p in sorted(OUT_DIR.glob("*.csv"))],
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    report_lines = [
        "# Macrophage seven-gene PPI, topology and target-actionability prioritization",
        "",
        f"Generated: {request_audit['run_timestamp_utc']}",
        "",
        "## Frozen input and evidence separation",
        "",
        f"- Frozen candidates: `{', '.join(genes)}`.",
        f"- STRING functional network: `{functional_edges.shape[0]}` high-confidence edges at combined score ≥0.700; no added interactors.",
        f"- STRING physical network: `{physical_edges.shape[0]}` high-confidence edge(s), reported separately from functional associations.",
        "- STRING topology is descriptive. It does not establish direction, causality, or DINP binding.",
        "- UniProt/PDBe/ChEMBL evidence is not merged into the STRING edge set; it is reported as target context.",
        "",
        "## Network result",
        "",
        f"The high-confidence functional subnetwork contains `{network_summary['n_input_nodes']}` input nodes, `{network_summary['n_edges']}` edges, `{network_summary['n_connected_components']}` connected components, and a largest component of `{network_summary['largest_component_nodes']}` nodes.",
        "",
        "Observed functional edges:",
        "",
    ]
    if functional_edges.empty:
        report_lines.append("- No high-confidence functional edges were returned among the seven inputs.")
    else:
        report_lines.extend([f"- `{row.gene_a} — {row.gene_b}` (score `{row.score:.3f}`)" for row in functional_edges.sort_values("score", ascending=False).itertuples()])
    report_lines += [
        "",
        "## Role-based interpretation",
        "",
        "- **Network bridge candidates:** highest topology among mapped proteins outside the direct prostaglandin class.",
        "- **Direct prostaglandin nodes:** pathway/protein-context candidates; they need not have a STRING edge inside this seven-node induced subgraph.",
        "- **Supporting/context nodes:** relevant expression or network context without enough evidence for the primary shortlist.",
        "- **NEAT1:** non-protein state/regulatory node; not docking-eligible.",
        "",
        "## Proposed docking/MD shortlist",
        "",
        "| Rank | Gene | Role | Functional degree | Topology score | PDB structures | Best coverage | ChEMBL measured activity |",
        "|---:|---|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in shortlist.iterrows():
        report_lines.append(
            f"| {int(row['shortlist_rank'])} | **{row['gene_symbol']}** | {row['shortlist_role']} | {int(row['functional_degree'])} | {row['topology_score']:.3f} | {int(row['pdb_structure_count'])} | {row['best_coverage']:.2f} | {int(row['chembl_max_single_protein_activity_count'])} |"
        )
    report_lines += [
        "",
        "### Reserve",
        "",
        "`PTGES3` remains a direct prostaglandin-synthase reserve candidate because it has an experimental structure and direct pathway relevance, but it is isolated in the seven-node STRING functional network and has a much smaller ChEMBL activity record count than the main shortlist.",
        "",
        "## Interpretation boundary",
        "",
        "The shortlist prioritizes candidates for a later structural workflow. It does not show that DINP binds any target, that a STRING edge is a direct physical interaction, or that any gene mediates the epidemiologic association. Docking/MD should begin only after ligand identity, protein construct, binding site, and assay evidence are independently frozen.",
        "",
        "Full STRING tables, pathway relation audit, node topology, UniProt/PDBe/ChEMBL context, request audit, and manifest are retained in `outputs/`.",
    ]
    (OUT_DIR / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps({"functional_edges": int(functional_edges.shape[0]), "physical_edges": int(physical_edges.shape[0]), "shortlist": shortlist_genes, "reserve": "PTGES3"}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
