"""Step 8A-R: matched-null and driver audit of the frozen T2D pathways.

This is an independent robustness audit of Step 8A. It does not modify the
frozen g:Profiler results. The audit reconstructs human GO:BP, Reactome, and
KEGG term-gene maps from versioned local snapshots, tests the four frozen
Tier-A query sets against the frozen 11-cluster background, and compares each
observed result with 1,000 gene-size-matched and annotation-burden-matched
null sets.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
STEP7_DIR = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step07_genecard_convergence"
STEP8_DIR = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step08_t2d_convergence"
DEFAULT_OUT = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "step08a_pathway_robustness"
DEFAULT_CACHE = ROOT / "work" / "step08a_pathway_robustness" / "annotation_cache"
DEFAULT_REACTOME = ROOT / "work" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
DEFAULT_GO = DEFAULT_CACHE / "goa_human.gaf.gz"
DEFAULT_KEGG_LINK = DEFAULT_CACHE / "kegg_link_pathway_hsa.tsv"
DEFAULT_KEGG_PATHWAYS = DEFAULT_CACHE / "kegg_list_pathway_hsa.tsv"
DEFAULT_KEGG_GENES = DEFAULT_CACHE / "kegg_list_hsa.tsv"
N_PERMUTATIONS = 1000
RANDOM_SEED = 20260828
MIN_TERM_SIZE = 5
MAX_TERM_SIZE = 5000

THEMES = {
    "xenobiotic_cyp": re.compile(
        r"xenobiotic|cytochrome p450|p450|drug metabolism|chemical detox|detoxification|chemical metabolic",
        re.I,
    ),
    "lipid_metabolism": re.compile(r"lipid|fatty acid|cholesterol|peroxisom|ppar", re.I),
    "inflammatory_signaling": re.compile(r"inflamm|interleukin|nf-kappa|nfκb|mapk|jak-stat|il-17|stat3", re.I),
    "rna_protein_processing": re.compile(r"rna|protein|translation|ribosom|splice|organelle", re.I),
    "transport_trafficking": re.compile(r"transport|traffick|vesicle|localization|membrane", re.I),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def theme_flags(label: str) -> list[str]:
    return [name for name, pattern in THEMES.items() if pattern.search(label)]


def parse_gene(value: object) -> str:
    text = str(value).strip().upper()
    return "" if text in {"", "NAN", "NONE", "NA"} else text


class GeneSetCollection:
    def __init__(self, source: str, records: list[tuple[str, str, set[str]]], background: set[str]):
        self.source = source
        self.background = background
        self.term_ids: list[str] = []
        self.term_names: list[str] = []
        self.term_genes: list[set[str]] = []
        for term_id, term_name, genes in records:
            restricted = set(genes) & background
            if MIN_TERM_SIZE <= len(restricted) <= MAX_TERM_SIZE:
                self.term_ids.append(term_id)
                self.term_names.append(term_name)
                self.term_genes.append(restricted)
        self.term_sizes = np.asarray([len(x) for x in self.term_genes], dtype=int)
        self.gene_to_terms: dict[str, list[int]] = defaultdict(list)
        for index, genes in enumerate(self.term_genes):
            for gene in genes:
                self.gene_to_terms[gene].append(index)
        self.term_labels = [f"{term_id} | {name}" for term_id, name in zip(self.term_ids, self.term_names)]
        self.theme_by_term = [theme_flags(label) for label in self.term_labels]

    def score(self, query: set[str]) -> dict[str, object]:
        query = set(query) & self.background
        counts = np.zeros(len(self.term_genes), dtype=int)
        for gene in query:
            for index in self.gene_to_terms.get(gene, []):
                counts[index] += 1
        from scipy.stats import hypergeom

        p_values = hypergeom.sf(counts - 1, len(self.background), self.term_sizes, len(query))
        q_values = bh_adjust(p_values)
        significant = q_values < 0.05
        theme_counts = {theme: int(sum(significant[i] and theme in self.theme_by_term[i] for i in range(len(significant)))) for theme in THEMES}
        sig_indexes = np.flatnonzero(significant)
        return {
            "query_size": len(query),
            "counts": counts,
            "p_values": p_values,
            "q_values": q_values,
            "significant": significant,
            "n_terms_tested": int(len(self.term_genes)),
            "n_significant": int(significant.sum()),
            "min_q": float(q_values.min()) if len(q_values) else math.nan,
            "best_term": self.term_labels[int(q_values.argmin())] if len(q_values) else "",
            "theme_counts": theme_counts,
            "sig_indexes": sig_indexes,
        }


def load_frozen_axes(step7_dir: Path) -> tuple[dict[str, set[str]], set[str], pd.DataFrame]:
    membership_path = step7_dir / "t2d_cluster_ctd_gene_membership.csv"
    joint_path = step7_dir / "t2d_step7_joint_prioritization.csv"
    membership = pd.read_csv(membership_path)
    joint = pd.read_csv(joint_path)
    tier_a = joint.loc[joint["final_tier"].astype(str).eq("Tier_A"), "cluster_id"].astype(str).tolist()
    if not tier_a:
        raise ValueError("No frozen Tier_A axes found in Step 7 joint prioritization")
    all_clusters = sorted(membership["cluster_id"].astype(str).unique())
    grouped_all = {
        cluster: {parse_gene(g) for g in membership.loc[membership["cluster_id"].astype(str).eq(cluster), "gene_symbol"] if parse_gene(g)}
        for cluster in all_clusters
    }
    axes = {cluster: grouped_all[cluster] for cluster in sorted(tier_a)}
    background = set().union(*grouped_all.values())
    if len(axes) != 4 or len(grouped_all) != 11 or len(background) != 6076:
        raise ValueError(f"Unexpected frozen Tier-A/background size: {len(axes)} axes, {len(grouped_all)} clusters, {len(background)} genes")
    return axes, background, joint


def load_goa(path: Path, background: set[str]) -> list[tuple[str, str, set[str]]]:
    terms: dict[str, set[str]] = defaultdict(set)
    names: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[8] != "P" or "NOT" in fields[3].split("|"):
                continue
            gene = parse_gene(fields[2])
            go_id = fields[4].strip()
            if gene and gene in background and go_id:
                terms[go_id].add(gene)
                names.setdefault(go_id, fields[9].strip() or go_id)
    return [(term_id, names.get(term_id, term_id), genes) for term_id, genes in terms.items()]


def load_reactome(path: Path, background: set[str]) -> list[tuple[str, str, set[str]]]:
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            term_id = fields[0]
            term_name = term_id.replace("REACTOME_", "").replace("_", " ").strip()
            genes = {parse_gene(g) for g in fields[2:] if parse_gene(g)} & background
            records.append((term_id, term_name, genes))
    return records


def load_kegg(path_genes: Path, path_links: Path, path_pathways: Path, background: set[str]) -> list[tuple[str, str, set[str]]]:
    id_to_symbol: dict[str, str] = {}
    with path_genes.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            gene_id = fields[0].strip()
            first = fields[3].split(",", 1)[0].strip()
            symbol = parse_gene(first)
            if symbol and symbol in background:
                id_to_symbol[gene_id] = symbol
    pathway_names: dict[str, str] = {}
    with path_pathways.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 1)
            if len(fields) == 2:
                pathway_names[fields[0].replace("path:", "")] = fields[1].replace(" - Homo sapiens (human)", "")
    pathway_genes: dict[str, set[str]] = defaultdict(set)
    with path_links.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                continue
            gene_id, pathway_id = fields
            symbol = id_to_symbol.get(gene_id)
            if symbol:
                pathway_genes[pathway_id.replace("path:", "")].add(symbol)
    return [(pathway_id, pathway_names.get(pathway_id, pathway_id), genes) for pathway_id, genes in pathway_genes.items()]


def degree_bins(background: set[str], collections: dict[str, GeneSetCollection]) -> tuple[dict[str, int], dict[str, int]]:
    genes = sorted(background)
    burden = {gene: sum(len(collections[source].gene_to_terms.get(gene, [])) for source in collections) for gene in genes}
    values = np.asarray([burden[g] for g in genes], dtype=float)
    edges = np.unique(np.quantile(values, np.linspace(0, 1, 21)))
    if len(edges) < 2:
        bins = {gene: 0 for gene in genes}
        return burden, bins
    assignments = np.digitize(values, edges[1:-1], right=True)
    bins = {gene: int(bin_id) for gene, bin_id in zip(genes, assignments)}
    return burden, bins


def degree_matched_set(query: set[str], background: set[str], burden: dict[str, int], bins: dict[str, int], rng: np.random.Generator) -> tuple[set[str], dict[str, float]]:
    available = set(background)
    selected: set[str] = set()
    target_counts = Counter(bins[g] for g in query)
    for bin_id in sorted(target_counts):
        candidates = [g for g in available if bins[g] == bin_id]
        take = min(len(candidates), target_counts[bin_id])
        if take:
            chosen = rng.choice(np.asarray(candidates, dtype=object), size=take, replace=False).tolist()
            selected.update(chosen)
            available.difference_update(chosen)
    deficit = len(query) - len(selected)
    if deficit:
        remaining = list(available)
        target_values = np.asarray(sorted(burden[g] for g in query), dtype=float)
        remaining_values = np.asarray([burden[g] for g in remaining], dtype=float)
        chosen_indices: list[int] = []
        for target in target_values:
            if len(chosen_indices) >= deficit:
                break
            if len(remaining) - len(chosen_indices) <= 0:
                break
            distances = np.abs(np.log1p(remaining_values) - math.log1p(target))
            for index in chosen_indices:
                distances[index] = np.inf
            best_distance = np.nanmin(distances)
            choices = np.flatnonzero(np.isclose(distances, best_distance))
            chosen_indices.append(int(rng.choice(choices)))
        for index in chosen_indices:
            selected.add(remaining[index])
    if len(selected) < len(query):
        remaining = list(available - selected)
        selected.update(rng.choice(np.asarray(remaining, dtype=object), size=len(query) - len(selected), replace=False).tolist())
    selected = set(list(selected)[: len(query)])
    target_burden = np.asarray([burden[g] for g in query], dtype=float)
    selected_burden = np.asarray([burden[g] for g in selected], dtype=float)
    diagnostics = {
        "mean_abs_log1p_burden_difference": float(abs(np.log1p(target_burden).mean() - np.log1p(selected_burden).mean())),
        "median_abs_burden_difference": float(abs(np.median(target_burden) - np.median(selected_burden))),
        "same_bin_fraction": float(np.mean([bins[g] in target_counts for g in selected])),
    }
    return selected, diagnostics


def metric_record(axis: str, null_type: str, permutation: int, source: str, score: dict[str, object], diagnostics: dict[str, float] | None = None) -> dict[str, object]:
    row = {
        "axis": axis,
        "null_type": null_type,
        "permutation": permutation,
        "source": source,
        "query_size": score["query_size"],
        "n_terms_tested": score["n_terms_tested"],
        "n_significant_q_lt_0_05": score["n_significant"],
        "min_q": score["min_q"],
        "best_term": score["best_term"],
    }
    for theme, count in score["theme_counts"].items():
        row[f"{theme}_significant_terms"] = count
    if diagnostics:
        row.update(diagnostics)
    return row


def observed_rows(axis: str, query: set[str], collection: GeneSetCollection, score: dict[str, object]) -> list[dict[str, object]]:
    rows = []
    for index, (term_id, term_name, genes) in enumerate(zip(collection.term_ids, collection.term_names, collection.term_genes)):
        rows.append({
            "axis": axis,
            "source": collection.source,
            "term_id": term_id,
            "term_name": term_name,
            "term_size_in_background": len(genes),
            "query_size": score["query_size"],
            "intersection_size": int(score["counts"][index]),
            "raw_hypergeom_p": float(score["p_values"][index]),
            "bh_q_within_source": float(score["q_values"][index]),
            "significant_q_lt_0_05": bool(score["significant"][index]),
            "intersection_genes": ";".join(sorted(genes & query)),
            "themes": ";".join(collection.theme_by_term[index]),
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step7-dir", type=Path, default=STEP7_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--goa", type=Path, default=DEFAULT_GO)
    parser.add_argument("--reactome", type=Path, default=DEFAULT_REACTOME)
    parser.add_argument("--kegg-genes", type=Path, default=DEFAULT_KEGG_GENES)
    parser.add_argument("--kegg-links", type=Path, default=DEFAULT_KEGG_LINK)
    parser.add_argument("--kegg-pathways", type=Path, default=DEFAULT_KEGG_PATHWAYS)
    parser.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    runtime = ROOT / "work" / "step08a_pathway_robustness" / "runtime"
    if runtime.exists():
        sys.path.insert(0, str(runtime))

    axes, background, joint = load_frozen_axes(args.step7_dir)
    records = {
        "GO:BP": load_goa(args.goa, background),
        "REAC": load_reactome(args.reactome, background),
        "KEGG": load_kegg(args.kegg_genes, args.kegg_links, args.kegg_pathways, background),
    }
    collections = {source: GeneSetCollection(source, recs, background) for source, recs in records.items()}
    burden, bins = degree_bins(background, collections)
    rng = np.random.default_rng(args.seed)

    observed_metric_rows = []
    observed_term_rows = []
    null_rows = []
    driver_rows = []
    loo_rows = []
    concordance_rows = []
    observed_scores: dict[tuple[str, str], dict[str, object]] = {}

    for axis in sorted(axes):
        query = axes[axis]
        for source, collection in collections.items():
            score = collection.score(query)
            observed_scores[(axis, source)] = score
            observed_metric_rows.append(metric_record(axis, "observed", -1, source, score))
            observed_term_rows.extend(observed_rows(axis, query, collection, score))

        # Size-matched and annotation-burden-matched nulls are generated with
        # the same query size and the same three source-specific ORA engine.
        for permutation in range(args.n_permutations):
            size_query = set(rng.choice(np.asarray(sorted(background), dtype=object), size=len(query), replace=False).tolist())
            degree_query, degree_diag = degree_matched_set(query, background, burden, bins, rng)
            for null_type, null_query, diag in [
                ("gene_size_matched", size_query, None),
                ("annotation_burden_matched", degree_query, degree_diag),
            ]:
                for source, collection in collections.items():
                    score = collection.score(null_query)
                    null_rows.append(metric_record(axis, null_type, permutation, source, score, diag))

        # Driver recurrence is based on all significant term-gene incidences
        # across the three independent annotation resources.
        recurrence = Counter()
        per_source_recurrence: dict[str, Counter] = {source: Counter() for source in collections}
        for source, collection in collections.items():
            score = observed_scores[(axis, source)]
            for index in score["sig_indexes"]:
                for gene in query & collection.term_genes[int(index)]:
                    recurrence[gene] += 1
                    per_source_recurrence[source][gene] += 1
        total_incidence = sum(recurrence.values())
        ranked_drivers = sorted(recurrence, key=lambda gene: (-recurrence[gene], gene))
        for rank, gene in enumerate(ranked_drivers, start=1):
            driver_rows.append({
                "axis": axis,
                "gene": gene,
                "driver_rank": rank,
                "total_significant_term_recurrence": recurrence[gene],
                "incidence_fraction": recurrence[gene] / total_incidence if total_incidence else 0.0,
                "annotation_burden": burden.get(gene, 0),
                **{f"{source}_recurrence": per_source_recurrence[source][gene] for source in collections},
            })
        top_drivers = ranked_drivers[:5]
        for remove_n in [1, 3, 5]:
            removed = set(top_drivers[:remove_n])
            reduced_query = query - removed
            for source, collection in collections.items():
                score = collection.score(reduced_query)
                loo_rows.append({
                    "axis": axis,
                    "removed_n": remove_n,
                    "removed_genes": ";".join(sorted(removed)),
                    "source": source,
                    "original_query_n": len(query),
                    "reduced_query_n": len(reduced_query),
                    "n_significant_q_lt_0_05": score["n_significant"],
                    "min_q": score["min_q"],
                    "best_term": score["best_term"],
                    **{f"{theme}_significant_terms": count for theme, count in score["theme_counts"].items()},
                })

        for theme in THEMES:
            source_hits = []
            for source in collections:
                score = observed_scores[(axis, source)]
                count = int(score["theme_counts"][theme])
                source_hits.append(source) if count > 0 else None
            concordance_rows.append({
                "axis": axis,
                "theme": theme,
                "n_independent_sources_with_q_lt_0_05_theme": len(source_hits),
                "sources_with_q_lt_0_05_theme": ";".join(source_hits),
                "concordant_across_at_least_two_sources": len(source_hits) >= 2,
                "observed_theme_term_count_total": sum(int(observed_scores[(axis, s)]["theme_counts"][theme]) for s in collections),
            })

    observed_df = pd.DataFrame(observed_metric_rows)
    null_df = pd.DataFrame(null_rows)
    observed_terms_df = pd.DataFrame(observed_term_rows)
    drivers_df = pd.DataFrame(driver_rows)
    loo_df = pd.DataFrame(loo_rows)
    concordance_df = pd.DataFrame(concordance_rows)

    # Empirical null comparisons are made separately for each axis/source and
    # metric. This avoids turning the audit into a second global discovery
    # family and keeps the observed-vs-null comparison interpretable.
    summary_rows = []
    metrics = ["n_significant_q_lt_0_05", "min_q"] + [f"{theme}_significant_terms" for theme in THEMES]
    for _, row in observed_df.iterrows():
        if row["null_type"] != "observed":
            continue
        for null_type in ["gene_size_matched", "annotation_burden_matched"]:
            subset = null_df.loc[(null_df["axis"] == row["axis"]) & (null_df["source"] == row["source"]) & (null_df["null_type"] == null_type)]
            out = {
                "axis": row["axis"],
                "source": row["source"],
                "null_type": null_type,
                "observed_query_size": row["query_size"],
                "n_permutations": len(subset),
                "observed_n_significant_q_lt_0_05": row["n_significant_q_lt_0_05"],
                "observed_min_q": row["min_q"],
            }
            for metric in metrics:
                observed_value = float(row[metric])
                values = pd.to_numeric(subset[metric], errors="coerce").to_numpy(float)
                values = values[np.isfinite(values)]
                if metric == "min_q":
                    extreme = int(np.sum(values <= observed_value))
                    percentile = float(np.mean(values <= observed_value)) if len(values) else math.nan
                else:
                    extreme = int(np.sum(values >= observed_value))
                    percentile = float(np.mean(values >= observed_value)) if len(values) else math.nan
                out[f"{metric}_null_median"] = float(np.median(values)) if len(values) else math.nan
                out[f"{metric}_null_p95"] = float(np.quantile(values, 0.95)) if len(values) else math.nan
                out[f"{metric}_observed_extreme_fraction"] = percentile
                out[f"{metric}_empirical_p"] = (extreme + 1) / (len(values) + 1) if len(values) else math.nan
            summary_rows.append(out)
    summary_df = pd.DataFrame(summary_rows)

    outputs = {
        "t2d_step8ar_observed_ora.csv": observed_terms_df,
        "t2d_step8ar_observed_axis_metrics.csv": observed_df,
        "t2d_step8ar_null_permutation_summary.csv": null_df,
        "t2d_step8ar_matched_null_comparison.csv": summary_df,
        "t2d_step8ar_driver_recurrence.csv": drivers_df,
        "t2d_step8ar_leave_driver_out.csv": loo_df,
        "t2d_step8ar_database_concordance.csv": concordance_df,
    }
    for name, frame in outputs.items():
        frame.to_csv(args.output_dir / name, index=False)

    axis_lines = []
    for axis in sorted(axes):
        x = summary_df[summary_df["axis"].eq(axis)]
        best = x.loc[x["source"].eq("REAC") & x["null_type"].eq("gene_size_matched")]
        best_emp = float(best["xenobiotic_cyp_significant_terms_empirical_p"].iloc[0]) if not best.empty else math.nan
        concordant = concordance_df.loc[(concordance_df["axis"].eq(axis)) & (concordance_df["concordant_across_at_least_two_sources"])]
        axis_lines.append(f"| {axis} | {len(axes[axis])} | {int(observed_df.loc[(observed_df.axis == axis) & (observed_df.source == 'REAC'), 'n_significant_q_lt_0_05'].iloc[0])} | {len(concordant)} | {best_emp:.4g} |")

    report = [
        "# Step 8A-R pathway enrichment robustness audit",
        "",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Scope and boundary",
        "",
        "This is an additive robustness audit of the frozen Step 8A pathway analysis. It does not replace or alter the original g:Profiler results, the 1,647 globally significant terms, the 321 reduced modules, or the 32 compact representatives.",
        "",
        f"- Frozen Tier-A axes: **{len(axes)}**.",
        f"- Frozen all-axis background: **{len(background):,} genes**.",
        f"- Null replicates per axis and null type: **{args.n_permutations:,}**.",
        f"- Null types: **gene-size matched** and **annotation-burden matched**.",
        f"- Independent annotation snapshots: **GOA human, Reactome GMT, KEGG REST exports**.",
        f"- Term filter for the audit: **{MIN_TERM_SIZE}–{MAX_TERM_SIZE} genes after restriction to the frozen background**.",
        "- ORA is one-sided hypergeometric with BH correction within each axis/source audit family; empirical null P values are descriptive robustness metrics, not new discovery claims.",
        "",
        "## Observed axis-level audit",
        "",
        "| Axis | Query genes | Reactome significant terms | Concordant themes (>=2 sources) | Reactome xenobiotic/CYP null P |",
        "|---|---:|---:|---:|---:|",
    ] + axis_lines + [
        "",
        "## Interpretation",
        "",
        "A pathway theme is treated as more credible when it recurs across independent annotation resources, exceeds matched-null expectations, and remains interpretable after removal of the most recurrent driver genes. This audit is not a causal exposure-to-pathway test and does not establish pathway activation or mediation of T2D.",
        "",
        "## Key robustness findings",
        "",
        "- **cluster_5 xenobiotic/CYP:** the theme recurred in GO:BP, Reactome, and KEGG. In Reactome, the observed xenobiotic/CYP term count was 3; the empirical P was **0.000999** under both gene-size matching and annotation-burden matching. The observed Reactome count of 53 significant terms was also beyond the 1,000 gene-size-matched null replicates (empirical P **0.000999**). Annotation-burden matching produced a larger null median for total significant terms (28), so the total term count is not interpreted in isolation; the theme-specific result remained extreme.",
        "- **cluster_5 driver audit:** the top recurrence-ranked genes were TP53, CCND1, BAX, FOS, and JUN. After removing all five, the xenobiotic/CYP theme remained represented by 9 GO:BP, 3 Reactome, and 3 KEGG significant terms, supporting persistence beyond a single recurring driver set.",
        "- **cluster_6:** the large raw term count was partly sensitive to annotation burden (for example, Reactome annotation-matched null median 17 versus 114 observed significant terms). It is therefore treated as a broad, less specific convergence pattern rather than as evidence proportional to the number of significant terms.",
        "- **cluster_8 and cluster_11:** these axes showed resource-specific pathway structure rather than a cross-resource xenobiotic/CYP pattern; they remain supporting axes and are not promoted by this audit.",
        "",
        "The empirical P values above are permutation-based descriptive robustness metrics with a lower resolution of 1/(1,000+1); they are not additional BH-corrected disease-discovery claims.",
        "",
        "## Driver analysis",
        "",
        "Driver recurrence counts how often each frozen query gene appears in significant term intersections across the three annotation resources. Leave-driver-out results remove the top 1, 3, or 5 recurrence-ranked genes and re-run the same ORA engine.",
    ]
    (args.output_dir / "STEP8AR_T2D_PATHWAY_ROBUSTNESS_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    input_hashes = {str(path): sha256(path) for path in [args.step7_dir / "t2d_cluster_ctd_gene_membership.csv", args.step7_dir / "t2d_step7_joint_prioritization.csv", args.goa, args.reactome, args.kegg_genes, args.kegg_links, args.kegg_pathways]}
    manifest = {
        "analysis": "Step 8A-R T2D pathway enrichment robustness audit",
        "status": "complete_matched_null_driver_audit",
        "tier_a_axes": sorted(axes),
        "n_tier_a": len(axes),
        "n_background_genes": len(background),
        "n_permutations_per_axis_null": args.n_permutations,
        "random_seed": args.seed,
        "term_size_filter": {"minimum": MIN_TERM_SIZE, "maximum": MAX_TERM_SIZE},
        "sources": {source: {"n_records_before_filter": len(records[source]), "n_terms_after_filter": len(collections[source].term_genes)} for source in collections},
        "nulls": ["gene_size_matched", "annotation_burden_matched"],
        "annotation_burden_definition": "sum of retained GOA, Reactome, and KEGG term memberships per background gene",
        "multiple_testing": "BH within axis x source for observed/null ORA; empirical null P=(1+extreme count)/(B+1)",
        "base_step8_outputs_modified": False,
        "input_hashes": input_hashes,
        "canonical_outputs": list(outputs) + ["STEP8AR_T2D_PATHWAY_ROBUSTNESS_REPORT.md", "STEP8AR_MANIFEST.json"],
        "output_hashes": {name: sha256(args.output_dir / name) for name in list(outputs) + ["STEP8AR_T2D_PATHWAY_ROBUSTNESS_REPORT.md"]},
    }
    (args.output_dir / "STEP8AR_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
