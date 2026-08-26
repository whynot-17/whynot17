"""Run the formal two-input Step 7 T2D CTD x GeneCards analysis.

The primary input is the complete ranked result of the ordinary GeneCards
query ``type 2 diabetes mellitus``.  The 111-row Disorders-scoped exact-phrase
capture is retained as a high-specificity sensitivity analysis.  Both inputs
are applied to all 11 frozen exposure clusters with the same CTD background
and the same one-sided hypergeometric test.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_step07_t2d_genecard_convergence as core


DEFAULT_CLUSTER = Path("analysis/disease_agnostic_environmental_framework/step06_t2d_robustness/t2d_exposure_clusters.csv")
DEFAULT_ROBUST = Path("analysis/disease_agnostic_environmental_framework/step06_t2d_robustness/t2d_robustness_results.csv")
DEFAULT_MEMBERSHIP = Path("analysis/disease_agnostic_environmental_framework/hypothesis_unit_audit/step4_test_chemical_membership.csv")
DEFAULT_OUT = Path("analysis/disease_agnostic_environmental_framework/step07_genecard_convergence")


def raw_gene_audit(path: Path, label: str, ctd_universe: set[str]) -> pd.DataFrame:
    """Preserve the public/official GeneCards row fields and mapping audit."""
    raw = pd.read_csv(path, dtype=str, keep_default_na=False, engine="python")
    symbol_col = core.find_column(list(raw.columns), ["GeneSymbol", "Symbol", "Gene Symbol"])
    rank_col = core.find_column(list(raw.columns), ["GeneCards_Rank", "Rank", "#", "GeneCards Rank"])
    if symbol_col is None or rank_col is None:
        raise ValueError(f"GeneCards {label} input must contain symbol and rank columns")
    rel_col = core.find_column(list(raw.columns), ["RelevanceScore", "Relevance Score", "Score"])
    know_col = core.find_column(list(raw.columns), ["KnowledgeScore", "Knowledge Score"])
    name_col = core.find_column(list(raw.columns), ["GeneName", "Name", "Gene Name"])
    type_col = core.find_column(list(raw.columns), ["GeneType", "Type", "Gene Type"])
    out = pd.DataFrame({
        "set_label": label,
        "rank": pd.to_numeric(raw[rank_col], errors="coerce"),
        "gene_symbol": raw[symbol_col].astype(str).str.strip().str.upper(),
        "gene_name": raw[name_col].astype(str).str.strip() if name_col else "",
        "gene_type": raw[type_col].astype(str).str.strip() if type_col else "",
        "relevance_score": pd.to_numeric(raw[rel_col], errors="coerce") if rel_col else np.nan,
        "knowledge_score": pd.to_numeric(raw[know_col], errors="coerce") if know_col else np.nan,
    })
    out = out.dropna(subset=["rank"])
    out = out[out["gene_symbol"].ne("") & (out["rank"] > 0)].sort_values("rank")
    out["duplicate_symbol"] = out["gene_symbol"].duplicated(keep=False)
    ranks = out["rank"].to_numpy()
    integer_ranks = ranks.astype(int)
    if not np.all(ranks == integer_ranks) or not np.array_equal(integer_ranks, np.arange(1, len(out) + 1)):
        raise ValueError(f"GeneCards {label} ranks are not continuous from 1 to N")
    if bool(out["duplicate_symbol"].any()):
        dupes = sorted(out.loc[out["duplicate_symbol"], "gene_symbol"].unique())
        raise ValueError(f"GeneCards {label} contains duplicate symbols: {dupes[:10]}")
    out["mapped_to_ctd_background"] = out["gene_symbol"].isin(ctd_universe)
    out["gene_type"] = out["gene_type"].replace("", "Unknown")
    return out.reset_index(drop=True)


def rank_cutoffs(cards: pd.DataFrame) -> list[int]:
    """Primary full list plus prespecified descriptive rank cutoffs."""
    max_rank = int(cards["rank"].max())
    return sorted(set([100, 500, 1000, 2000, max_rank]))


def one_set(
    label: str,
    cards_path: Path,
    cluster_genes: dict[str, set[str]],
    out_dir: Path,
    output_stem: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    cards = core.load_genecards(cards_path)
    cards["rank"] = pd.to_numeric(cards["rank"], errors="coerce")
    cards = cards.dropna(subset=["rank"]).copy()
    ks = rank_cutoffs(cards)
    enrichment, overlaps = core.enrichment(cluster_genes, cards, ks, out_dir)
    enrichment["set_label"] = label
    overlaps["set_label"] = label
    full_k = int(cards["rank"].max())
    full = enrichment[enrichment["gene_cards_k"].eq(full_k)].copy()
    full_overlap = overlaps[overlaps["gene_cards_k"].eq(full_k)].copy()
    full.to_csv(out_dir / f"t2d_cluster_enrichment_{output_stem}.csv", index=False)
    full_overlap.to_csv(out_dir / f"t2d_cluster_genecards_overlap_{output_stem}.csv", index=False)
    rank_table = enrichment[[
        "set_label", "gene_cards_k", "cluster_id", "n_cluster_ctd_genes",
        "n_background_ctd_genes", "n_t2d_genecard_genes", "n_overlap",
        "rank_weighted_overlap", "odds_ratio", "hypergeom_p", "bh_fdr",
    ]].copy()
    metadata = {
        "set_label": label,
        "path": str(cards_path),
        "n_ranked_rows": int(len(cards)),
        "max_rank": int(cards["rank"].max()),
        "rank_cutoffs_reported": ks,
        "primary_full_k": full_k,
        "n_primary_clusters": int(len(full)),
        "n_q_lt_0_05": int((full["bh_fdr"] < 0.05).sum()),
        "minimum_q": float(full["bh_fdr"].min()) if not full.empty else None,
    }
    return full, full_overlap, rank_table, cards, metadata


def make_joint(
    clusters: pd.DataFrame,
    cluster_chemicals: pd.DataFrame,
    robust: pd.DataFrame,
    primary: pd.DataFrame,
    strict: pd.DataFrame,
    out_dir: Path,
) -> pd.DataFrame:
    base = []
    robust_lookup = robust.set_index("test_id", drop=False).to_dict("index")
    for cluster_id, group in clusters.groupby("cluster_id"):
        tests = sorted(group["test_id"].astype(str).unique())
        rs = [robust_lookup[x] for x in tests if x in robust_lookup]
        p = primary[primary["cluster_id"].eq(cluster_id)].iloc[0]
        s = strict[strict["cluster_id"].eq(cluster_id)].iloc[0]
        n_robust = sum(x.get("priority_tier") == "robust_fdr_candidate" for x in rs)
        if float(p["bh_fdr"]) < 0.05 and n_robust:
            tier = "Tier_A"
        elif n_robust:
            tier = "Tier_B"
        else:
            tier = "Tier_C"
        base.append({
            "cluster_id": cluster_id,
            "biomarkers": ";".join(sorted(group["biomarker"].unique())),
            "n_tests": len(tests),
            "n_chemicals": int(cluster_chemicals.loc[cluster_chemicals["cluster_id"].eq(cluster_id), "chemical_id"].nunique()),
            "n_ctd_genes": int(p["n_cluster_ctd_genes"]),
            "robust_test_count": n_robust,
            "robust_tests": ";".join(x["biomarker"] for x in rs if x.get("priority_tier") == "robust_fdr_candidate"),
            "primary_full_n_overlap": int(p["n_overlap"]),
            "primary_full_odds_ratio": float(p["odds_ratio"]),
            "primary_full_hypergeom_p": float(p["hypergeom_p"]),
            "primary_full_bh_fdr": float(p["bh_fdr"]),
            "primary_full_rank_weighted_overlap": float(p["rank_weighted_overlap"]),
            "strict_111_n_overlap": int(s["n_overlap"]),
            "strict_111_odds_ratio": float(s["odds_ratio"]),
            "strict_111_hypergeom_p": float(s["hypergeom_p"]),
            "strict_111_bh_fdr": float(s["bh_fdr"]),
            "strict_111_rank_weighted_overlap": float(s["rank_weighted_overlap"]),
            "final_tier": tier,
        })
    out = pd.DataFrame(base).sort_values(["final_tier", "primary_full_bh_fdr", "primary_full_hypergeom_p"])
    out.to_csv(out_dir / "t2d_step7_joint_prioritization.csv", index=False)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-map", type=Path, default=DEFAULT_CLUSTER)
    parser.add_argument("--robustness", type=Path, default=DEFAULT_ROBUST)
    parser.add_argument("--membership", type=Path, default=DEFAULT_MEMBERSHIP)
    parser.add_argument("--ctd", type=Path, required=True)
    parser.add_argument("--genecards-primary", type=Path, required=True)
    parser.add_argument("--genecards-strict", type=Path, required=True)
    parser.add_argument("--primary-query", default="type 2 diabetes mellitus")
    parser.add_argument("--strict-query", default='[Disorders] "type 2 diabetes mellitus"')
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    clusters, robust, membership = core.load_cluster_inputs(args.cluster_map, args.robustness, args.membership)
    cluster_chemicals = core.build_cluster_chemicals(clusters, membership)
    pair, evidence, cluster_genes = core.load_ctd_side(args.ctd, cluster_chemicals)
    evidence = evidence.merge(
        cluster_chemicals[["cluster_id", "chemical_id", "chemical_name", "test_id", "biomarker", "step2_mapping_type", "step2_mapping_confidence", "recorded_parent_compound", "candidate_hypothesis_unit"]].drop_duplicates(),
        on="chemical_id", how="right",
    )
    evidence.to_csv(args.output_dir / "t2d_step7_ctd_chemical_evidence.csv", index=False)
    cluster_chemicals.to_csv(args.output_dir / "t2d_step7_cluster_chemical_map.csv", index=False)
    pair.merge(cluster_chemicals[["cluster_id", "chemical_id"]].drop_duplicates(), on="chemical_id", how="left").to_csv(
        args.output_dir / "t2d_step7_ctd_interaction_pairs.csv", index=False
    )

    gene_rows = [{"cluster_id": cid, "gene_symbol": gene} for cid, genes in sorted(cluster_genes.items()) for gene in sorted(genes)]
    ctd_universe = set().union(*cluster_genes.values()) if cluster_genes else set()
    pd.DataFrame(gene_rows).to_csv(args.output_dir / "t2d_cluster_ctd_gene_membership.csv", index=False)

    primary_audit = raw_gene_audit(args.genecards_primary, "primary_anywhere", ctd_universe)
    strict_audit = raw_gene_audit(args.genecards_strict, "strict_disorders_exact", ctd_universe)
    primary_audit.to_csv(args.output_dir / "t2d_genecards_primary_gene_audit.csv", index=False)
    strict_audit.to_csv(args.output_dir / "t2d_genecards_strict_gene_audit.csv", index=False)
    pd.concat([primary_audit, strict_audit], ignore_index=True).groupby(["set_label", "gene_type"], dropna=False).size().rename("n_rows").reset_index().to_csv(
        args.output_dir / "t2d_genecards_gene_type_audit.csv", index=False
    )

    primary_full, primary_overlap, primary_rank, primary_cards, primary_meta = one_set(
        "primary_anywhere", args.genecards_primary, cluster_genes, args.output_dir, "primary"
    )
    strict_full, strict_overlap, strict_rank, strict_cards, strict_meta = one_set(
        "strict_disorders_exact", args.genecards_strict, cluster_genes, args.output_dir, "strict"
    )
    pd.concat([primary_rank, strict_rank], ignore_index=True).to_csv(args.output_dir / "t2d_rank_weighted_convergence.csv", index=False)

    overlap = pd.concat([primary_overlap, strict_overlap], ignore_index=True)
    overlap.to_csv(args.output_dir / "t2d_overlap_gene_details.csv", index=False)

    joint = make_joint(clusters, cluster_chemicals, robust, primary_full, strict_full, args.output_dir)

    manifest = {
        "analysis": "Step 7 T2D-specific CTD x GeneCards biological convergence",
        "status": "complete_two_gene_sets",
        "n_step6_clusters": int(clusters["cluster_id"].nunique()),
        "n_cluster_chemical_rows": int(len(cluster_chemicals)),
        "n_ctd_interaction_pairs": int(len(pair)),
        "n_ctd_cluster_genes": int(sum(len(x) for x in cluster_genes.values())),
        "ctd_background": "union of all 11 frozen cluster human CTD genes",
        "primary": {**primary_meta, "query": args.primary_query, "acquisition_mode": "public pagination table capture", "sha256": core.sha256(args.genecards_primary)},
        "strict_sensitivity": {**strict_meta, "query": args.strict_query, "acquisition_mode": "public pagination table capture", "sha256": core.sha256(args.genecards_strict)},
        "outcome_firewall": "T2D outcome was not used to construct the 29 tests or 11 clusters; GeneCards enters only at post-firewall Step 7.",
        "parent_mapping": "Inherited from Step 4; no parent relationships inferred from names.",
        "ctd_deduplication": "unique ChemicalID x GeneID among Homo sapiens interactions",
        "multiple_testing": "BH-FDR across 11 clusters separately for primary and strict sets",
        "primary_q_lt_0_05": int((primary_full["bh_fdr"] < 0.05).sum()),
        "strict_q_lt_0_05": int((strict_full["bh_fdr"] < 0.05).sum()),
        "canonical_outputs": [
            "t2d_cluster_ctd_gene_membership.csv",
            "t2d_cluster_enrichment_primary.csv",
            "t2d_cluster_enrichment_strict.csv",
            "t2d_cluster_genecards_overlap_primary.csv",
            "t2d_cluster_genecards_overlap_strict.csv",
            "t2d_genecards_primary_gene_audit.csv",
            "t2d_genecards_strict_gene_audit.csv",
            "t2d_genecards_gene_type_audit.csv",
            "t2d_rank_weighted_convergence.csv",
            "t2d_overlap_gene_details.csv",
            "t2d_step7_joint_prioritization.csv",
            "STEP7_MANIFEST.json",
            "STEP7_T2D_GENE_CONVERGENCE_REPORT.md",
        ],
    }
    (args.output_dir / "STEP7_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    primary_sig = primary_full[primary_full["bh_fdr"] < 0.05]["cluster_id"].tolist()
    strict_sig = strict_full[strict_full["bh_fdr"] < 0.05]["cluster_id"].tolist()
    report = [
        "# Step 7 — T2D-specific CTD × GeneCards biological convergence",
        "",
        "- Status: **complete_two_gene_sets**",
        f"- Frozen exposure clusters: **{manifest['n_step6_clusters']}**",
        f"- CTD human chemical–gene pairs represented: **{manifest['n_ctd_interaction_pairs']:,}**",
        f"- Cluster CTD gene memberships summed over clusters: **{manifest['n_ctd_cluster_genes']:,}**",
        "",
        "## GeneCards input sets",
        "",
        f"- Primary ordinary query: `{args.primary_query}`; complete public result: **{primary_meta['n_ranked_rows']:,} rows**.",
        f"- High-specificity sensitivity: `{args.strict_query}`; complete public result: **{strict_meta['n_ranked_rows']:,} rows**.",
        "- Primary analysis uses the complete ordinary-query list; the 111-row exact Disorders result is not used as the primary set.",
        "- Primary full-list and strict-set enrichment are each corrected across the same 11-cluster family.",
        "",
        "## Primary full-list result",
        "",
        f"- Primary q < 0.05 clusters: **{len(primary_sig)}** — {', '.join(primary_sig) if primary_sig else 'none'}.",
        f"- Minimum primary q: **{primary_meta['minimum_q']:.4g}**.",
        "",
        "## Strict exact-phrase sensitivity",
        "",
        f"- Strict q < 0.05 clusters: **{len(strict_sig)}** — {', '.join(strict_sig) if strict_sig else 'none'}.",
        f"- Minimum strict q: **{strict_meta['minimum_q']:.4g}**.",
        "",
        "## Cluster-level primary summary",
        "",
        "| Cluster | Biomarker(s) | CTD genes | GeneCards overlap | OR | q |",
        "|---|---|---:|---:|---:|---:|",
    ]
    biomarker_lookup = clusters.groupby("cluster_id")["biomarker"].apply(lambda x: ";".join(sorted(x.astype(str).unique())))
    for _, row in primary_full.sort_values("bh_fdr").iterrows():
        report.append(
            f"| {row['cluster_id']} | {biomarker_lookup.get(row['cluster_id'], '')} | "
            f"{int(row['n_cluster_ctd_genes'])} | {int(row['n_overlap'])} | "
            f"{float(row['odds_ratio']):.3g} | {float(row['bh_fdr']):.4g} |"
        )
    report += [
        "",
        "## Interpretation boundary",
        "",
        "The GeneCards analysis is post-firewall biological prioritization. It does not modify the 29-test T2D screen, the 14 FDR-positive tests, the 13 robustness-supported tests, or the 11 exposure clusters. CTD chemical–gene associations and GeneCards disease associations are not causal evidence.",
        "",
    ]
    (args.output_dir / "STEP7_T2D_GENE_CONVERGENCE_REPORT.md").write_text("\n".join(report), encoding="utf-8")


if __name__ == "__main__":
    main()
