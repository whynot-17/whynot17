from __future__ import annotations

from itertools import combinations
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def bh_adjust(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    q = np.full(len(p), np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(q, index=values.index)
    valid_idx = np.flatnonzero(valid)
    order = valid_idx[np.argsort(p[valid_idx])]
    m = len(order)
    running = 1.0
    for rank in range(m, 0, -1):
        idx = order[rank - 1]
        running = min(running, p[idx] * m / rank)
        q[idx] = running
    return pd.Series(q, index=values.index)


def require(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")


def compare_sets(
    set_a_name: str,
    set_a_label: str,
    set_a: set[str],
    set_b_name: str,
    set_b_label: str,
    set_b: set[str],
    universe: set[str],
) -> dict:
    if not set_a.issubset(universe) or not set_b.issubset(universe):
        raise ValueError(f"Set outside the 97-gene universe: {set_a_name} or {set_b_name}")
    intersection = sorted(set_a & set_b)
    a_only = sorted(set_a - set_b)
    b_only = sorted(set_b - set_a)
    neither = sorted(universe - set_a - set_b)
    table = [[len(intersection), len(a_only)], [len(b_only), len(neither)]]
    odds_ratio, p_two_sided = fisher_exact(table, alternative="two-sided")
    _, p_greater = fisher_exact(table, alternative="greater")
    _, p_less = fisher_exact(table, alternative="less")
    expected_overlap = len(set_a) * len(set_b) / len(universe)
    return {
        "set_a": set_a_name,
        "set_a_label": set_a_label,
        "set_b": set_b_name,
        "set_b_label": set_b_label,
        "background_n": len(universe),
        "set_a_n": len(set_a),
        "set_b_n": len(set_b),
        "overlap_n": len(intersection),
        "expected_overlap": expected_overlap,
        "fold_enrichment_observed_over_expected": len(intersection) / expected_overlap if expected_overlap else np.nan,
        "odds_ratio": float(odds_ratio),
        "fisher_p_two_sided": float(p_two_sided),
        "fisher_p_overlap_greater": float(p_greater),
        "fisher_p_overlap_less": float(p_less),
        "jaccard_index": len(intersection) / len(set_a | set_b) if (set_a | set_b) else np.nan,
        "a_only_n": len(a_only),
        "b_only_n": len(b_only),
        "neither_n": len(neither),
        "overlap_genes": ";".join(intersection),
        "a_only_genes": ";".join(a_only),
        "b_only_genes": ";".join(b_only),
    }


def main() -> None:
    overlap = pd.read_csv(OUTPUTS / "DINP_CRC_overlap.csv", dtype=str).fillna("")
    stability = pd.read_csv(OUTPUTS / "DINP_CRC_101ML_gene_stability.csv", dtype=str).fillna("")
    cross_rank = pd.read_csv(OUTPUTS / "DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv", dtype=str).fillna("")
    require(overlap, ["gene_symbol"], "DINP-CRC overlap")
    require(stability, ["gene_symbol", "ml_rank", "stable_ml_flag", "ml_priority_score"], "101ML gene stability")
    require(
        cross_rank,
        [
            "gene_symbol",
            "cross_rank",
            "cross_support_flag",
            "ppi_pathway_candidate",
            "independent_datasets_significant_fdr_lt_0_05",
            "cross_rank_score",
        ],
        "PPI/pathway/transcriptomic cross-ranking",
    )

    universe = set(overlap["gene_symbol"].str.strip().str.upper()) - {""}
    if len(universe) != 97:
        raise AssertionError(f"Expected 97 background genes, got {len(universe)}")
    for frame in [stability, cross_rank]:
        frame["gene_symbol"] = frame["gene_symbol"].str.strip().str.upper()
    stability = stability.drop_duplicates("gene_symbol")
    cross_rank = cross_rank.drop_duplicates("gene_symbol")
    if len(stability) != 97 or len(cross_rank) != 97:
        raise AssertionError("Input audit tables must each contain 97 unique genes")

    membership = pd.DataFrame({"gene_symbol": sorted(universe)})
    membership = membership.merge(
        stability[["gene_symbol", "ml_rank", "stable_ml_flag", "ml_priority_score", "expression_measured"]],
        on="gene_symbol",
        how="left",
    ).merge(
        cross_rank[
            [
                "gene_symbol",
                "cross_rank",
                "cross_support_flag",
                "ppi_pathway_candidate",
                "independent_datasets_significant_fdr_lt_0_05",
                "representative_term_count",
                "hub_rank",
                "cross_rank_score",
                "candidate_tier",
                "consensus_direction",
            ]
        ],
        on="gene_symbol",
        how="left",
    )
    membership["stable_ml"] = as_bool(membership["stable_ml_flag"])
    membership["tier1"] = as_bool(membership["cross_support_flag"])
    membership["ppi_pathway_candidate_flag"] = as_bool(membership["ppi_pathway_candidate"])
    membership["transcriptomic_replicated_ge2"] = pd.to_numeric(
        membership["independent_datasets_significant_fdr_lt_0_05"], errors="coerce"
    ).ge(2)
    membership["crossrank_top20"] = pd.to_numeric(membership["cross_rank"], errors="coerce").le(20)
    membership["overlap_group_stable_vs_tier1"] = np.select(
        [membership["stable_ml"] & membership["tier1"], membership["stable_ml"], membership["tier1"]],
        ["both_stable_ml_and_tier1", "stable_ml_only", "tier1_only"],
        default="neither",
    )
    membership = membership.sort_values(["stable_ml", "tier1", "gene_symbol"], ascending=[False, False, True]).reset_index(drop=True)

    sets = {
        "stable_ml": ("17 stable-ML genes", set(membership.loc[membership["stable_ml"], "gene_symbol"])),
        "tier1": ("16 Tier 1 cross-supported genes", set(membership.loc[membership["tier1"], "gene_symbol"])),
        "ppi_pathway_candidate": (
            "PPI top-20 hub + pathway membership",
            set(membership.loc[membership["ppi_pathway_candidate_flag"], "gene_symbol"]),
        ),
        "transcriptomic_replicated_ge2": (
            "Transcriptomic replicated in ≥2/3 primary datasets",
            set(membership.loc[membership["transcriptomic_replicated_ge2"], "gene_symbol"]),
        ),
        "crossrank_top20": ("Top 20 cross-ranked genes", set(membership.loc[membership["crossrank_top20"], "gene_symbol"])),
    }
    expected_sizes = {"stable_ml": 17, "tier1": 16}
    for key, expected in expected_sizes.items():
        if len(sets[key][1]) != expected:
            raise AssertionError(f"Expected {expected} genes in {key}, got {len(sets[key][1])}")

    pair_rows: list[dict] = []
    for (key_a, (label_a, set_a)), (key_b, (label_b, set_b)) in combinations(sets.items(), 2):
        pair_rows.append(compare_sets(key_a, label_a, set_a, key_b, label_b, set_b, universe))
    pairwise = pd.DataFrame(pair_rows)
    pairwise["fisher_q_two_sided_bh"] = bh_adjust(pairwise["fisher_p_two_sided"])
    pairwise["primary_comparison"] = (pairwise["set_a"].eq("stable_ml") & pairwise["set_b"].eq("tier1"))
    pairwise = pairwise.sort_values(["primary_comparison", "fisher_p_two_sided", "set_a", "set_b"], ascending=[False, True, True, True]).reset_index(drop=True)
    primary = pairwise[pairwise["primary_comparison"]].copy()
    primary.to_csv(OUTPUTS / "DINP_CRC_101ML_Tier1_fisher_exact.csv", index=False)
    pairwise.to_csv(OUTPUTS / "DINP_CRC_101ML_overlap_pairwise_audit.csv", index=False)
    membership.to_csv(OUTPUTS / "DINP_CRC_101ML_Tier1_overlap_membership.csv", index=False)

    row = primary.iloc[0]
    summary_lines = [
        "# DINP–CRC 101ML × Tier 1 overlap and enrichment audit",
        "",
        "## Primary comparison",
        "",
        "The primary test compares the 17 default stable-ML genes with the 16 pre-existing Tier 1 cross-supported genes inside the frozen 97-gene DINP–CRC overlap universe.",
        "",
        "| Quantity | Value |",
        "|---|---:|",
        f"| Background universe | {int(row['background_n'])} |",
        f"| Stable-ML genes | {int(row['set_a_n'])} |",
        f"| Tier 1 genes | {int(row['set_b_n'])} |",
        f"| Observed overlap | {int(row['overlap_n'])} |",
        f"| Expected overlap under random overlap | {row['expected_overlap']:.3f} |",
        f"| Observed / expected | {row['fold_enrichment_observed_over_expected']:.3f} |",
        f"| Odds ratio | {row['odds_ratio']:.4f} |",
        f"| Fisher exact p, two-sided | {row['fisher_p_two_sided']:.6g} |",
        f"| Fisher exact p, overlap greater | {row['fisher_p_overlap_greater']:.6g} |",
        f"| Fisher exact p, overlap less | {row['fisher_p_overlap_less']:.6g} |",
        f"| BH q, two-sided pairwise audit | {row['fisher_q_two_sided_bh']:.6g} |",
        f"| Jaccard index | {row['jaccard_index']:.4f} |",
        f"| Overlap genes | {row['overlap_genes']} |",
        "",
        "## 2×2 table",
        "",
        "Rows are stable-ML status; columns are Tier 1 status.",
        "",
        "| | Tier 1 yes | Tier 1 no |",
        "|---|---:|---:|",
        f"| Stable ML yes | {int(row['overlap_n'])} | {int(row['a_only_n'])} |",
        f"| Stable ML no | {int(row['b_only_n'])} | {int(row['neither_n'])} |",
        "",
        "## Pairwise overlap audit",
        "",
        "The full pairwise table tests the stable-ML and Tier 1 sets against PPI/pathway candidates, transcriptomic replication ≥2/3, and the top-20 cross-ranked set, always using the same 97-gene background. BH q-values adjust the two-sided Fisher tests across all pairwise comparisons.",
        "",
        "| Set A | Set B | A | B | Overlap | Expected | Fold enrichment | Odds ratio | Fisher p | BH q |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, current in pairwise.iterrows():
        summary_lines.append(
            f"| {current['set_a_label']} | {current['set_b_label']} | {int(current['set_a_n'])} | {int(current['set_b_n'])} | {int(current['overlap_n'])} | {current['expected_overlap']:.3f} | {current['fold_enrichment_observed_over_expected']:.3f} | {current['odds_ratio']:.3f} | {current['fisher_p_two_sided']:.4g} | {current['fisher_q_two_sided_bh']:.4g} |"
        )
    summary_lines.extend(
        [
            "",
            "## Gene-level membership",
            "",
            "The 97-row membership table preserves stable-ML status, Tier 1 status, PPI/pathway candidate status, transcriptomic replication status, and both ranking positions for every background gene.",
            "",
            "## Interpretation",
            "",
            "The observed stable-ML × Tier 1 overlap is interpreted against the 97-gene background only. Fisher's exact test evaluates set overlap, not biological causality, model performance, or MR validity. The two sets were constructed independently, so limited overlap is informative rather than a failure by itself.",
        ]
    )
    (OUTPUTS / "DINP_CRC_101ML_Tier1_overlap_summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    log = {
        "background_definition": "97 unique DINP-CRC overlap gene symbols from DINP_CRC_overlap.csv",
        "set_definitions": {
            "stable_ml": "stable_ml_flag == True in DINP_CRC_101ML_gene_stability.csv",
            "tier1": "cross_support_flag == True in DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv",
            "ppi_pathway_candidate": "ppi_pathway_candidate == True",
            "transcriptomic_replicated_ge2": "independent_datasets_significant_fdr_lt_0_05 >= 2",
            "crossrank_top20": "cross_rank <= 20",
        },
        "stable_ml_n": int(len(sets["stable_ml"][1])),
        "tier1_n": int(len(sets["tier1"][1])),
        "primary_overlap_n": int(row["overlap_n"]),
        "primary_overlap_genes": row["overlap_genes"].split(";") if row["overlap_genes"] else [],
        "fisher_table": [[int(row["overlap_n"]), int(row["a_only_n"])], [int(row["b_only_n"]), int(row["neither_n"])]],
        "multiple_testing": "Benjamini-Hochberg adjustment across the 10 two-sided pairwise Fisher tests",
        "qc": {
            "background_n": len(universe),
            "membership_rows": len(membership),
            "membership_gene_symbols_unique": int(membership["gene_symbol"].nunique()) == 97,
            "pairwise_comparisons": len(pairwise),
        },
    }
    (OUTPUTS / "DINP_CRC_101ML_Tier1_overlap_log.md").write_text(
        "# DINP–CRC 101ML × Tier 1 overlap audit log\n\n```json\n"
        + json.dumps(log, indent=2, ensure_ascii=False)
        + "\n```\n",
        encoding="utf-8",
    )

    print("Primary stable-ML x Tier 1 comparison:")
    print(primary[["background_n", "set_a_n", "set_b_n", "overlap_n", "expected_overlap", "fold_enrichment_observed_over_expected", "odds_ratio", "fisher_p_two_sided", "fisher_p_overlap_greater", "fisher_p_overlap_less", "overlap_genes"]].to_string(index=False))
    print("Pairwise comparisons:", len(pairwise))


if __name__ == "__main__":
    main()
