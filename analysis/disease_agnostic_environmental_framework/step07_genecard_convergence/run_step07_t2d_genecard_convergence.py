"""Run the frozen Step 7 T2D CTD × GeneCards convergence analysis.

The script can be run in ``--ctd-only`` mode to audit the exposure-cluster and
CTD side before a disease-specific GeneCards export is supplied. It deliberately
does not use the old CRC GeneCards file as a fallback.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom


DEFAULT_CLUSTER = Path("analysis/disease_agnostic_environmental_framework/step06_t2d_robustness/t2d_exposure_clusters.csv")
DEFAULT_ROBUST = Path("analysis/disease_agnostic_environmental_framework/step06_t2d_robustness/t2d_robustness_results.csv")
DEFAULT_MEMBERSHIP = Path("analysis/disease_agnostic_environmental_framework/hypothesis_unit_audit/step4_test_chemical_membership.csv")
DEFAULT_OUT = Path("analysis/disease_agnostic_environmental_framework/step07_genecard_convergence")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical_id(value: object) -> str:
    text = str(value or "").strip()
    return text[5:] if text.startswith("MESH:") else text


def split_semicolon(value: object) -> list[str]:
    return [x.strip() for x in str(value or "").split(";") if x.strip()]


def split_pipe(value: object) -> list[str]:
    return [x.strip() for x in str(value or "").split("|") if x.strip()]


def read_ctd(path: Path) -> pd.DataFrame:
    """Read CTD's commented field header without treating it as a data row."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("# Fields:"):
                header = next(handle).lstrip("# ").rstrip("\r\n").split("\t")
                return pd.read_csv(
                    handle,
                    sep="\t",
                    names=header,
                    dtype=str,
                    keep_default_na=False,
                    comment="#",
                    low_memory=False,
                )
    raise ValueError(f"Could not find CTD '# Fields:' header in {path}")


def pmids(value: object) -> set[str]:
    return {x.strip() for x in re.split(r"[|;,]", str(value or "")) if x.strip()}


def load_cluster_inputs(cluster_path: Path, robust_path: Path, membership_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clusters = pd.read_csv(cluster_path, dtype=str, keep_default_na=False)
    robust = pd.read_csv(robust_path, dtype=str, keep_default_na=False)
    membership = pd.read_csv(membership_path, dtype=str, keep_default_na=False)
    if "test_id" not in clusters.columns and "variable" in clusters.columns:
        clusters = clusters.rename(columns={"variable": "test_id"})
    required_cluster = {"test_id", "biomarker", "cluster_id"}
    required_membership = {"variable", "test_id", "chemical_id", "chemical_name"}
    missing = sorted(required_cluster - set(clusters.columns))
    if missing:
        raise ValueError(f"Cluster mapping missing columns: {missing}")
    missing = sorted(required_membership - set(membership.columns))
    if missing:
        raise ValueError(f"Chemical membership missing columns: {missing}")
    clusters = clusters.copy()
    clusters["test_id"] = clusters["test_id"].astype(str)
    clusters["cluster_id"] = clusters["cluster_id"].astype(str)
    if robust.empty:
        raise ValueError("Step 6 robustness table is empty")
    return clusters, robust, membership


def build_cluster_chemicals(clusters: pd.DataFrame, membership: pd.DataFrame) -> pd.DataFrame:
    left = clusters[["test_id", "biomarker", "cluster_id"]].drop_duplicates()
    m = membership.copy()
    m["test_id"] = m["test_id"].astype(str)
    m["chemical_id"] = m["chemical_id"].map(canonical_id)
    m = m[m["chemical_id"].ne("")].copy()
    out = left.merge(m, on="test_id", how="left", suffixes=("", "_membership"))
    out["cluster_id"] = out["cluster_id"].astype(str)
    out["mapping_inherited_from_step4"] = out["chemical_id"].notna()
    return out.drop_duplicates(["cluster_id", "chemical_id", "test_id"])


def load_ctd_side(path: Path, cluster_chemicals: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    raw = read_ctd(path)
    required = {"ChemicalID", "ChemicalName", "GeneID", "GeneSymbol", "Organism", "PubMedIDs"}
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"CTD interaction file missing columns: {missing}")
    human = raw[(raw["Organism"].eq("Homo sapiens")) & raw["ChemicalID"].ne("") & raw["GeneID"].ne("")].copy()
    human["chemical_id"] = human["ChemicalID"].map(canonical_id)
    human["gene_id"] = human["GeneID"].astype(str).str.strip()
    human["gene_symbol"] = human["GeneSymbol"].astype(str).str.strip().str.upper()
    human["pmid_set"] = human["PubMedIDs"].map(pmids)
    pair = human[["chemical_id", "ChemicalName", "gene_id", "gene_symbol"]].drop_duplicates(["chemical_id", "gene_id"])
    ids = set(cluster_chemicals["chemical_id"].dropna().astype(str))
    pair = pair[pair["chemical_id"].isin(ids)].copy()
    evidence = human[human["chemical_id"].isin(ids)].groupby("chemical_id").agg(
        n_raw_interaction_rows=("chemical_id", "size"),
        n_unique_pmids=("pmid_set", lambda x: len(set().union(*x) if len(x) else set())),
    ).reset_index()
    evidence = evidence.merge(
        pair.groupby("chemical_id").size().rename("n_unique_chemical_gene_pairs").reset_index(),
        on="chemical_id",
        how="outer",
    ).fillna(0)
    gene_sets = pair.groupby("chemical_id")["gene_symbol"].agg(lambda x: {g for g in x if g}).to_dict()
    cluster_genes = {}
    for cluster_id, group in cluster_chemicals.groupby("cluster_id"):
        chem_ids = set(group["chemical_id"].dropna().astype(str))
        cluster_genes[str(cluster_id)] = set().union(*(gene_sets.get(cid, set()) for cid in chem_ids))
    return pair, evidence, cluster_genes


def find_column(columns: list[str], candidates: list[str]) -> str | None:
    normal = {re.sub(r"[^a-z0-9]", "", c.lower()): c for c in columns}
    for candidate in candidates:
        key = re.sub(r"[^a-z0-9]", "", candidate.lower())
        if key in normal:
            return normal[key]
    return None


def load_genecards(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False, engine="python")
    symbol_col = find_column(list(df.columns), ["GeneSymbol", "Symbol", "Gene Symbol"])
    rank_col = find_column(list(df.columns), ["GeneCards_Rank", "Rank", "#", "GeneCards Rank"])
    if symbol_col is None or rank_col is None:
        raise ValueError("GeneCards export must contain symbol and rank columns")
    out = pd.DataFrame({"gene_symbol": df[symbol_col].astype(str).str.strip().str.upper(), "rank": pd.to_numeric(df[rank_col], errors="coerce")})
    rel_col = find_column(list(df.columns), ["RelevanceScore", "Relevance Score", "Score"])
    know_col = find_column(list(df.columns), ["KnowledgeScore", "Knowledge Score"])
    if rel_col:
        out["relevance_score"] = pd.to_numeric(df[rel_col], errors="coerce")
    if know_col:
        out["knowledge_score"] = pd.to_numeric(df[know_col], errors="coerce")
    out = out.dropna(subset=["gene_symbol", "rank"])
    out = out[out["gene_symbol"].ne("") & (out["rank"] > 0)].sort_values("rank").drop_duplicates("gene_symbol")
    if out.empty:
        raise ValueError("GeneCards export yielded no ranked genes")
    return out.reset_index(drop=True)


def bh(values: list[float]) -> list[float]:
    arr = np.asarray(values, dtype=float)
    order = np.argsort(arr)
    adj = np.empty(len(arr), dtype=float)
    running = 1.0
    for i in range(len(arr) - 1, -1, -1):
        rank = i + 1
        running = min(running, arr[order[i]] * len(arr) / rank)
        adj[order[i]] = running
    return adj.tolist()


def enrichment(cluster_genes: dict[str, set[str]], cards: pd.DataFrame, ks: list[int], out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    all_background = set().union(*cluster_genes.values()) if cluster_genes else set()
    for k in ks:
        gc = cards[cards["rank"] <= k].copy()
        gc_genes = set(gc["gene_symbol"])
        universe = all_background
        if not universe or not gc_genes:
            continue
        n = len(universe)
        m = len(universe & gc_genes)
        pvals = []
        staged: list[dict[str, object]] = []
        for cluster_id, genes in sorted(cluster_genes.items()):
            hits = sorted(genes & gc_genes)
            x = len(hits)
            k_cluster = len(genes)
            p = float(hypergeom.sf(x - 1, n, m, k_cluster)) if k_cluster else 1.0
            a = x
            b = max(k_cluster - x, 0)
            c = max(m - x, 0)
            d = max(n - m - b, 0)
            try:
                odds = float(fisher_exact([[a, b], [c, d]], alternative="greater").statistic)
            except Exception:
                odds = math.inf if a and not c else float("nan")
            weighted = float(sum(1.0 / math.log2(float(r) + 1.0) for r in gc.loc[gc["gene_symbol"].isin(hits), "rank"]))
            record = {
                "gene_cards_k": k,
                "cluster_id": cluster_id,
                "n_cluster_ctd_genes": k_cluster,
                "n_background_ctd_genes": n,
                "n_t2d_genecard_genes": m,
                "n_overlap": x,
                "odds_ratio": odds,
                "hypergeom_p": p,
                "rank_weighted_overlap": weighted,
                "overlap_genes": ";".join(hits),
            }
            staged.append(record)
            pvals.append(p)
            for gene in hits:
                overlap_rows.append({"gene_cards_k": k, "cluster_id": cluster_id, "gene_symbol": gene})
        adjusted = bh(pvals)
        for record, q in zip(staged, adjusted):
            record["bh_fdr"] = q
            rows.append(record)
    return pd.DataFrame(rows), pd.DataFrame(overlap_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-map", type=Path, default=DEFAULT_CLUSTER)
    parser.add_argument("--robustness", type=Path, default=DEFAULT_ROBUST)
    parser.add_argument("--membership", type=Path, default=DEFAULT_MEMBERSHIP)
    parser.add_argument("--ctd", type=Path, required=True)
    parser.add_argument("--genecards", type=Path)
    parser.add_argument("--genecards-query", default="")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--k", type=int, nargs="+", default=[500, 1000, 2000])
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    clusters, robust, membership = load_cluster_inputs(args.cluster_map, args.robustness, args.membership)
    cluster_chemicals = build_cluster_chemicals(clusters, membership)
    pair, evidence, cluster_genes = load_ctd_side(args.ctd, cluster_chemicals)
    evidence = evidence.merge(
        cluster_chemicals[["cluster_id", "chemical_id", "chemical_name", "test_id", "biomarker", "step2_mapping_type", "step2_mapping_confidence", "recorded_parent_compound", "candidate_hypothesis_unit"]].drop_duplicates(),
        on="chemical_id",
        how="right",
    )
    evidence.to_csv(args.output_dir / "t2d_step7_ctd_chemical_evidence.csv", index=False)
    cluster_chemicals.to_csv(args.output_dir / "t2d_step7_cluster_chemical_map.csv", index=False)
    pair.merge(
        cluster_chemicals[["cluster_id", "chemical_id"]].drop_duplicates(),
        on="chemical_id",
        how="left",
    ).to_csv(args.output_dir / "t2d_step7_ctd_interaction_pairs.csv", index=False)

    gene_rows = []
    for cluster_id, genes in sorted(cluster_genes.items()):
        for symbol in sorted(genes):
            gene_rows.append({"cluster_id": cluster_id, "gene_symbol": symbol})
    pd.DataFrame(gene_rows).to_csv(args.output_dir / "t2d_step7_cluster_ctd_genes.csv", index=False)

    enrichment_table = pd.DataFrame()
    overlap_table = pd.DataFrame()
    genecard_meta: dict[str, object] = {"status": "missing", "query": args.genecards_query}
    if args.genecards:
        cards = load_genecards(args.genecards)
        genecard_meta = {
            "status": "loaded",
            "query": args.genecards_query,
            "path": str(args.genecards),
            "sha256": sha256(args.genecards),
            "n_ranked_rows": int(len(cards)),
            "max_rank": float(cards["rank"].max()),
        }
        available_ks = [k for k in args.k if (cards["rank"] <= k).sum() > 0]
        enrichment_table, overlap_table = enrichment(cluster_genes, cards, available_ks, args.output_dir)
        enrichment_table.to_csv(args.output_dir / "t2d_step7_cluster_genecard_enrichment.csv", index=False)
        overlap_table.to_csv(args.output_dir / "t2d_step7_genecard_overlap_genes.csv", index=False)

    cluster_summary = []
    robust_lookup = robust.set_index("test_id", drop=False).to_dict("index")
    for cluster_id, group in clusters.groupby("cluster_id"):
        test_ids = sorted(group["test_id"].unique())
        rs = [robust_lookup[x] for x in test_ids if x in robust_lookup]
        cluster_summary.append({
            "cluster_id": cluster_id,
            "biomarkers": ";".join(sorted(group["biomarker"].unique())),
            "n_tests": len(test_ids),
            "n_chemicals": int(cluster_chemicals.loc[cluster_chemicals["cluster_id"].eq(cluster_id), "chemical_id"].nunique()),
            "n_ctd_genes": len(cluster_genes.get(str(cluster_id), set())),
            "max_primary_bh_fdr": min(float(x.get("BH_FDR", "nan")) for x in rs) if rs else float("nan"),
            "robust_tests": ";".join(x["biomarker"] for x in rs if x.get("priority_tier") == "robust_fdr_candidate"),
        })
    pd.DataFrame(cluster_summary).to_csv(args.output_dir / "t2d_step7_cluster_summary.csv", index=False)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "Step 7 T2D-specific CTD x GeneCards convergence",
        "n_step6_clusters": int(clusters["cluster_id"].nunique()),
        "cluster_ids": sorted(clusters["cluster_id"].unique().tolist()),
        "n_cluster_chemical_rows": int(len(cluster_chemicals)),
        "n_ctd_interaction_pairs": int(len(pair)),
        "n_ctd_cluster_genes": int(sum(len(x) for x in cluster_genes.values())),
        "genecards": genecard_meta,
        "outcome_firewall": "T2D outcome is not used to construct the 29 frozen tests or 11 clusters; Step 7 is post-firewall biological prioritization.",
        "parent_mapping": "Inherited from Step 4; no parent relationships inferred from names.",
        "ctd_deduplication": "unique ChemicalID x GeneID among Homo sapiens interactions",
        "status": "complete_with_genecards" if args.genecards else "ctd_preflight_complete_genecards_input_missing",
    }
    (args.output_dir / "t2d_step7_analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    report = [
        "# Step 7 T2D-specific CTD × GeneCards convergence",
        "",
        f"- Status: **{manifest['status']}**",
        f"- Frozen exposure clusters analyzed: **{manifest['n_step6_clusters']}**",
        f"- CTD human chemical–gene pairs represented: **{manifest['n_ctd_interaction_pairs']:,}**",
        f"- Cluster CTD gene memberships (cluster-level counts summed): **{manifest['n_ctd_cluster_genes']:,}**",
        "- Parent/proxy mapping: inherited from Step 4; no name-based parent inference.",
        "",
    ]
    if args.genecards and not enrichment_table.empty:
        primary = enrichment_table[enrichment_table["gene_cards_k"].eq(1000)].sort_values(["bh_fdr", "hypergeom_p"])
        report += ["## Primary T2D GeneCards enrichment (K=1000)", "", primary.to_markdown(index=False), ""]
    else:
        report += ["## GeneCards dependency", "", "No T2D GeneCards export was supplied. The CTD-side preflight is complete; overlap/enrichment was intentionally not fabricated or substituted with the CRC GeneCards set.", ""]
    (args.output_dir / "t2d_step7_convergence_report.md").write_text("\n".join(report), encoding="utf-8")


if __name__ == "__main__":
    main()
