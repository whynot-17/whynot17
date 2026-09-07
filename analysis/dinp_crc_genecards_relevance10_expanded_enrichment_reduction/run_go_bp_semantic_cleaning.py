#!/usr/bin/env python3
"""Perform a GO:BP-specific semantic redundancy audit.

The existing primary enrichment result contains 99 significant GO:BP terms.
This additive script uses the official go-basic ontology DAG plus the
18-gene term overlaps already reconstructed in the preceding audit.  More
specific terms are retained over enriched ancestors when their query-gene
sets are sufficiently redundant.  Non-ancestor terms are never folded solely
because they share genes.  No P values, FDR values, gene list, or background
are recomputed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from run_term_reduction_and_driver_contribution import (
    API_URL,
    ENRICHMENT_DIR,
    INPUT_QUERY_FILE,
    OUT_DIR,
    REPO_ROOT,
    THEME_RULES,
    finite_float,
    read_csv,
    sha256_file,
    theme_for_term,
    unique_symbols,
    write_csv,
)


GO_BASIC_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
GO_BASIC_CACHE = REPO_ROOT / "work" / "gene_sets" / "go-basic.obo"
INPUT_TERM_FILE = OUT_DIR / "primary_106_terms_with_genes.csv"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def download_ontology(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        GO_BASIC_URL,
        headers={"User-Agent": "whynot17-go-bp-cleaning/1.0"},
        method="GET",
    )
    with urllib.request.urlopen(request, timeout=300) as response:
        path.write_bytes(response.read())


def parse_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    parents: dict[str, set[str]] = {}
    names: dict[str, str] = {}
    current: dict[str, Any] | None = None

    def flush() -> None:
        nonlocal current
        if not current or current.get("id", "").startswith("GO:") is False:
            current = None
            return
        go_id = str(current["id"])
        parents[go_id] = set(current.get("parents", set()))
        names[go_id] = str(current.get("name", ""))
        current = None

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "[Term]":
                flush()
                current = {"id": "", "name": "", "parents": set()}
                continue
            if line.startswith("["):
                flush()
                current = None
                continue
            if current is None:
                continue
            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("is_a: "):
                parent = line[6:].split("!", 1)[0].strip()
                if parent.startswith("GO:"):
                    current["parents"].add(parent)
    flush()
    return parents, names


def ancestor_closure(go_id: str, parents: dict[str, set[str]], memo: dict[str, set[str]]) -> set[str]:
    if go_id in memo:
        return memo[go_id]
    ancestors: set[str] = set()
    for parent in parents.get(go_id, set()):
        ancestors.add(parent)
        ancestors.update(ancestor_closure(parent, parents, memo))
    memo[go_id] = ancestors
    return ancestors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ancestor-jaccard", type=float, default=0.50)
    args = parser.parse_args()
    if not 0 < args.ancestor_jaccard <= 1:
        raise ValueError("--ancestor-jaccard must be between 0 and 1")

    out = OUT_DIR / "go_bp_cleaning"
    out.mkdir(parents=True, exist_ok=True)
    if not GO_BASIC_CACHE.exists() or GO_BASIC_CACHE.stat().st_size < 1000:
        download_ontology(GO_BASIC_CACHE)
    ontology_parents, ontology_names = parse_obo(GO_BASIC_CACHE)

    rows = [row for row in read_csv(INPUT_TERM_FILE) if row["source"] == "GO:BP"]
    if len(rows) != 99:
        raise ValueError(f"Expected 99 GO:BP terms, found {len(rows)}")
    query_genes = unique_symbols([row["gene_symbol"] for row in read_csv(INPUT_QUERY_FILE)])
    if len(query_genes) != 18:
        raise ValueError(f"Expected 18 query genes, found {len(query_genes)}")

    terms: dict[str, dict[str, Any]] = {}
    memo: dict[str, set[str]] = {}
    for row in rows:
        term_id = row["term_id"]
        genes = set(filter(None, row["intersection_genes"].split(";")))
        if not genes:
            raise ValueError(f"Missing reconstructed query genes for {term_id}")
        if term_id not in ontology_parents:
            raise ValueError(f"GO term missing from go-basic ontology: {term_id}")
        terms[term_id] = {
            **row,
            "genes": genes,
            "term_size_int": int(row["term_size"]),
            "fdr_float": finite_float(row["global_fdr"], 1.0),
            "ancestors": ancestor_closure(term_id, ontology_parents, memo),
            "depth": len(ancestor_closure(term_id, ontology_parents, memo)),
            "theme_id": theme_for_term(row["term_name"])[0],
            "theme_label": theme_for_term(row["term_name"])[1],
        }

    # Specificity-first greedy pruning.  A term can be folded into an already
    # retained, more-specific term when it is a DAG ancestor with sufficiently
    # redundant query-gene membership.  Non-ancestor terms are never folded
    # solely because they share query genes.
    order = sorted(
        terms,
        key=lambda term_id: (
            -terms[term_id]["depth"],
            terms[term_id]["term_size_int"],
            terms[term_id]["fdr_float"],
            term_id,
        ),
    )
    retained: list[str] = []
    folded_into: dict[str, str] = {}
    fold_reason: dict[str, str] = {}
    fold_similarity: dict[str, float] = {}
    for term_id in order:
        term = terms[term_id]
        candidates: list[tuple[float, str, str]] = []
        for kept_id in retained:
            kept = terms[kept_id]
            similarity = jaccard(term["genes"], kept["genes"])
            if kept_id in term["ancestors"] and similarity >= args.ancestor_jaccard:
                candidates.append((similarity, kept_id, "ancestor-descendant redundancy"))
            elif term_id in kept["ancestors"] and similarity >= args.ancestor_jaccard:
                # This can occur when depth ties are resolved by term size/FDR.
                candidates.append((similarity, kept_id, "ancestor-descendant redundancy"))
        if candidates:
            similarity, kept_id, reason = max(candidates, key=lambda item: (item[0], -terms[item[1]]["term_size_int"], item[1]))
            folded_into[term_id] = kept_id
            fold_reason[term_id] = reason
            fold_similarity[term_id] = similarity
        else:
            retained.append(term_id)

    # If a term was initially folded into a term that is itself folded, point
    # the record to the final retained representative.
    def final_rep(term_id: str) -> str:
        seen: set[str] = set()
        while term_id in folded_into:
            if term_id in seen:
                raise ValueError(f"Cycle in representative mapping at {term_id}")
            seen.add(term_id)
            term_id = folded_into[term_id]
        return term_id

    for term_id in list(folded_into):
        folded_into[term_id] = final_rep(term_id)

    representative_rows: list[dict[str, Any]] = []
    cleaned_rows: list[dict[str, Any]] = []
    for term_id in sorted(terms):
        term = terms[term_id]
        representative = final_rep(term_id)
        kept = term_id == representative
        cleaned_rows.append(
            {
                "term_id": term_id,
                "term_name": term["term_name"],
                "global_fdr": term["global_fdr"],
                "intersection_size": term["intersection_size"],
                "intersection_genes": ";".join(sorted(term["genes"])),
                "go_depth": term["depth"],
                "term_theme_id": term["theme_id"],
                "term_theme_label": term["theme_label"],
                "status": "retained_representative" if kept else "folded",
                "representative_term_id": representative,
                "fold_reason": "" if kept else fold_reason[term_id],
                "gene_jaccard_to_representative": "" if kept else f"{fold_similarity[term_id]:.6f}",
            }
        )
        if kept:
            representative_rows.append(
                {
                    "representative_term_id": term_id,
                    "term_name": term["term_name"],
                    "global_fdr": term["global_fdr"],
                    "intersection_size": term["intersection_size"],
                    "intersection_genes": ";".join(sorted(term["genes"])),
                    "go_depth": term["depth"],
                    "term_theme_id": term["theme_id"],
                    "term_theme_label": term["theme_label"],
                    "folded_term_count": sum(1 for value in folded_into.values() if value == term_id),
                }
            )

    # Re-run the driver contribution on the cleaned representative set only.
    family_order = {theme_id: order for order, (theme_id, _label, _patterns) in enumerate(THEME_RULES)}
    family_order["other"] = len(THEME_RULES)
    family_ids = {
        theme_id: f"family_{position:02d}"
        for position, theme_id in enumerate(
            sorted({terms[term_id]["theme_id"] for term_id in retained}, key=lambda value: (family_order.get(value, 999), value)),
            start=1,
        )
    }
    family_groups: dict[str, list[str]] = defaultdict(list)
    for term_id in retained:
        family_groups[terms[term_id]["theme_id"]].append(term_id)
    family_rows: list[dict[str, Any]] = []
    family_driver_rows: list[dict[str, Any]] = []
    global_driver_counts: Counter[str] = Counter()
    for term_id in retained:
        for gene in terms[term_id]["genes"]:
            global_driver_counts[gene] += 1
    global_driver_rank = {gene: rank for rank, (gene, _count) in enumerate(global_driver_counts.most_common(), start=1)}
    for theme_id in sorted(family_groups, key=lambda value: (family_order.get(value, 999), value)):
        term_ids = family_groups[theme_id]
        family_id = family_ids[theme_id]
        label = terms[term_ids[0]]["theme_label"]
        union_genes = set().union(*(terms[term_id]["genes"] for term_id in term_ids))
        source_count = len({terms[term_id]["source"] for term_id in term_ids})
        representative = min(term_ids, key=lambda term_id: (terms[term_id]["fdr_float"], terms[term_id]["term_size_int"], term_id))
        family_rows.append(
            {
                "family_id": family_id,
                "theme_id": theme_id,
                "theme_label": label,
                "representative_term_count": len(term_ids),
                "source_count": source_count,
                "representative_term": terms[representative]["term_name"],
                "union_gene_count": len(union_genes),
            }
        )
        local_counts = Counter(gene for term_id in term_ids for gene in terms[term_id]["genes"])
        for gene, count in sorted(local_counts.items(), key=lambda item: (-item[1], global_driver_rank.get(item[0], 999), item[0])):
            family_driver_rows.append(
                {
                    "family_id": family_id,
                    "theme_id": theme_id,
                    "theme_label": label,
                    "gene_symbol": gene,
                    "representative_term_count": count,
                    "fraction_of_family_representatives": count / len(term_ids),
                    "family_top_driver": count == max(local_counts.values()),
                    "global_cleaned_term_count": global_driver_counts[gene],
                    "global_driver_rank": global_driver_rank[gene],
                }
            )

    out_clean = out / "go_bp_99_terms_cleaning.csv"
    out_reps = out / "go_bp_clean_representatives.csv"
    out_families = out / "go_bp_clean_theme_families.csv"
    out_family_drivers = out / "go_bp_clean_theme_driver_contribution.csv"
    out_summary = out / "go_bp_cleaning_summary.md"
    out_manifest = out / "go_bp_cleaning_manifest.json"
    write_csv(out_clean, cleaned_rows, list(cleaned_rows[0].keys()))
    write_csv(out_reps, representative_rows, list(representative_rows[0].keys()))
    write_csv(out_families, family_rows, list(family_rows[0].keys()))
    write_csv(out_family_drivers, family_driver_rows, list(family_driver_rows[0].keys()))

    summary_lines = [
        "# GO:BP semantic cleaning of the 99 primary significant terms",
        "",
        f"Generated (UTC): {utc_now()}",
        "",
        "## Scope",
        "",
        "This is an additive GO-specific redundancy audit. It does not change the frozen 18-gene query, 881-gene GeneCards background, original P values/FDR values, or the original 106-term output. Reactome terms were not included in this GO cleaning.",
        "",
        "## Method",
        "",
        f"The official `go-basic.obo` is used to recover the GO `is_a` DAG. Terms were ordered from more specific to broader (ontology depth, then term size, then original global FDR). A term was folded into a retained representative only when it was an enriched ancestor/descendant pair with query-gene Jaccard >= **{args.ancestor_jaccard:.2f}**. Non-ancestor terms were never folded based only on shared query genes, because a small query set can make unrelated terms look similar. All mappings are recorded term by term.",
        "",
        f"- GO:BP significant terms: **99**",
        f"- Retained GO representatives: **{len(representative_rows)}**",
        f"- Folded GO terms: **{len(folded_into)}**",
        f"- GO theme families among representatives: **{len(family_rows)}**",
        "",
        "## Clean theme families",
        "",
        "| Family | Theme | Retained representatives | Representative sources | Example representative | Union genes |",
        "|---|---|---:|---:|---|---:|",
    ]
    for row in family_rows:
        summary_lines.append(
            f"| {row['family_id']} | {row['theme_label']} | {row['representative_term_count']} | {row['source_count']} | {row['representative_term']} | {row['union_gene_count']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Cleaned driver contribution",
            "",
            "Driver counts below use retained GO representatives only. They are descriptive repeated-membership counts, not gene-level tests.",
            "",
            "| Rank | Gene | Retained GO terms |",
            "|---:|---|---:|",
        ]
    )
    for rank, (gene, count) in enumerate(global_driver_counts.most_common(10), start=1):
        summary_lines.append(f"| {rank} | {gene} | {count} |")
    summary_lines.extend(
        [
            "",
            "## Boundary",
            "",
            "A cleaned representative set is not a new enrichment family and does not make the underlying 99 GO terms independent. It is intended to support readable figures/tables and prevent broad GO ancestors from being counted as separate biological stories. The original 99-term results remain the statistical record.",
            "",
        ]
    )
    out_summary.write_text("\n".join(summary_lines), encoding="utf-8")

    manifest = {
        "generated_utc": utc_now(),
        "script": str(Path(__file__).resolve()),
        "input_term_file": str(INPUT_TERM_FILE.resolve()),
        "input_term_file_sha256": sha256_file(INPUT_TERM_FILE),
        "go_ontology_url": GO_BASIC_URL,
        "go_ontology_cache": str(GO_BASIC_CACHE.resolve()),
        "go_ontology_cache_sha256": sha256_file(GO_BASIC_CACHE),
        "go_ontology_term_count": len(ontology_parents),
        "primary_go_bp_term_count": 99,
        "retained_representative_count": len(representative_rows),
        "folded_term_count": len(folded_into),
        "theme_family_count": len(family_rows),
        "ancestor_jaccard_threshold": args.ancestor_jaccard,
        "non_ancestor_jaccard_threshold": None,
        "original_enrichment_unchanged": True,
        "reactome_mixed_in": False,
        "query_genes": query_genes,
        "outputs": {
            name: {
                "path": str((out / name).resolve()),
                "sha256": sha256_file(out / name),
            }
            for name in [
                "go_bp_99_terms_cleaning.csv",
                "go_bp_clean_representatives.csv",
                "go_bp_clean_theme_families.csv",
                "go_bp_clean_theme_driver_contribution.csv",
                "go_bp_cleaning_summary.md",
            ]
        },
    }
    out_manifest.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"GO:BP terms: 99")
    print(f"Retained GO representatives: {len(representative_rows)}")
    print(f"Folded GO terms: {len(folded_into)}")
    print(f"GO theme families: {len(family_rows)}")
    print(f"Summary: {out_summary}")
    print(f"Manifest: {out_manifest}")


if __name__ == "__main__":
    main()
