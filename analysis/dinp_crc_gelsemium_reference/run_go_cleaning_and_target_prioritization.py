#!/usr/bin/env python3
"""Clean the fresh 1,893-gene CRC enrichment and prioritize targets.

This script is deliberately scoped to the published Gelsemium-derived CRC
reference set and the fresh 41-gene DINP--CRC intersection.  It performs an
additive interpretation layer only:

* GO:BP significant terms are reduced using the official GO DAG and
  ancestor/descendant query-gene Jaccard redundancy.
* Cleaned GO representatives are aligned descriptively with significant
  Reactome and KEGG terms using frozen name-based biological themes.
* Genes are ordered by transparent cross-resource recurrence and same-theme
  support.  No new enrichment P values, FDR values, or hidden weighted score
  are introduced.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"
ENRICHMENT = OUTPUTS / "enrichment_all_backgrounds.csv"
QUERY = OUTPUTS / "query_41_genes.csv"
# ROOT = <repo>/analysis/dinp_crc_gelsemium_reference; its second parent is
# the repository root.  Keep the ontology cache in the repository-level work
# directory rather than creating a nested <repo>/work/work directory.
GO_CACHE = ROOT.parents[1] / "work" / "gene_sets" / "go-basic.obo"
GO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
BACKGROUND = "gelsemium_crc1893"
FDR_CUTOFF = 0.05
ANCESTOR_JACCARD = 0.50


THEMES: list[tuple[str, str, list[str]]] = [
    (
        "lipid_eicosanoid",
        "Lipid / fatty-acid / eicosanoid metabolism",
        [
            r"fatty acid",
            r"arachidon",
            r"eicosanoid",
            r"prostagland",
            r"prostanoid",
            r"thrombox",
            r"lipid",
            r"steroid",
            r"cholesterol",
            r"sphingolipid",
            r"glycerophospholipid",
        ],
    ),
    (
        "chemical_xenobiotic_stress",
        "Chemical / xenobiotic / oxygen-stress response",
        [r"xenobiotic", r"chemical", r"oxygen", r"oxidative", r"detox", r"cytochrome p450"],
    ),
    (
        "inflammation_immune",
        "Inflammation / immune / stimulus response",
        [r"inflamm", r"immune", r"interleukin", r"cytokine", r"defense", r"toll-like", r"stimulus"],
    ),
    (
        "signaling_receptor",
        "Signaling / receptor / hormone response",
        [r"signal", r"receptor", r"hormone", r"gpcr", r"nuclear receptor", r"transcription factor"],
    ),
    (
        "matrix_migration",
        "Matrix remodeling / adhesion / migration",
        [r"matrix", r"metalloproteinase", r"collagen", r"adhesion", r"migration", r"extracellular"],
    ),
    (
        "cell_death_cycle",
        "Cell death / cell cycle / genome maintenance",
        [r"apopt", r"cell death", r"cell cycle", r"dna repair", r"p53", r"chromosome"],
    ),
    (
        "transport_localization",
        "Transport / localization / vesicle trafficking",
        [r"transport", r"localization", r"traffick", r"vesicle", r"endocyt", r"membrane organization"],
    ),
    (
        "general_metabolism",
        "General metabolism / biosynthesis",
        [r"metabolic", r"biosynthetic", r"catabolic", r"metabolism", r"biosynthesis"],
    ),
]
THEME_LABEL = {theme_id: label for theme_id, label, _ in THEMES}
THEME_ORDER = {theme_id: index for index, (theme_id, _, _) in enumerate(THEMES)}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fdr(row: dict[str, str]) -> float:
    try:
        value = float(row["gprofiler_adjusted_p"])
    except (KeyError, TypeError, ValueError):
        return math.inf
    return value if math.isfinite(value) else math.inf


def int_value(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def symbols(text: str) -> set[str]:
    return {part.strip().upper() for part in (text or "").split(";") if part.strip()}


def theme_for(name: str) -> tuple[str, str, int]:
    lowered = (name or "").lower()
    scored: list[tuple[int, int, str, str]] = []
    for order, (theme_id, label, patterns) in enumerate(THEMES):
        score = sum(1 for pattern in patterns if re.search(pattern, lowered))
        if score:
            scored.append((score, -order, theme_id, label))
    if not scored:
        return "other", "Other / mixed biological processes", 0
    score, _neg_order, theme_id, label = max(scored)
    return theme_id, label, score


def download_go_obo(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(GO_URL, headers={"User-Agent": "whynot17-gelsemium-go-cleaning/1.0"})
    with urllib.request.urlopen(request, timeout=300) as response:
        path.write_bytes(response.read())


def parse_go_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    parents: dict[str, set[str]] = {}
    names: dict[str, str] = {}
    current: dict[str, Any] | None = None

    def flush() -> None:
        nonlocal current
        if current and current.get("id", "").startswith("GO:"):
            go_id = str(current["id"])
            # Preserve obsolete terms returned by the enrichment service.  If
            # the ontology supplies a replacement term, treat that replacement
            # as a DAG parent for redundancy auditing while retaining the
            # original term ID/name in the output and audit trail.
            term_parents = set(current.get("parents", set()))
            replacement = current.get("replacement", "")
            if current.get("obsolete", False) and replacement:
                term_parents.add(replacement)
            parents[go_id] = term_parents
            names[go_id] = str(current.get("name", ""))
        current = None

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                flush()
                current = {"id": "", "name": "", "parents": set(), "obsolete": False, "replacement": ""}
            elif line.startswith("["):
                flush()
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:].strip()
                elif line.startswith("name: "):
                    current["name"] = line[6:].strip()
                elif line.startswith("is_a: "):
                    parent = line[6:].split("!", 1)[0].strip()
                    if parent.startswith("GO:"):
                        current["parents"].add(parent)
                elif line.startswith("is_obsolete: true"):
                    current["obsolete"] = True
                elif line.startswith("replaced_by: GO:"):
                    current["replacement"] = line.split(": ", 1)[1].strip()
    flush()
    return parents, names


def ancestor_closure(go_id: str, parents: dict[str, set[str]], memo: dict[str, set[str]]) -> set[str]:
    if go_id in memo:
        return memo[go_id]
    if go_id not in parents:
        memo[go_id] = set()
        return memo[go_id]
    ancestors: set[str] = set()
    for parent in parents[go_id]:
        ancestors.add(parent)
        ancestors.update(ancestor_closure(parent, parents, memo))
    memo[go_id] = ancestors
    return ancestors


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def final_representative(term_id: str, folded_into: dict[str, str]) -> str:
    seen: set[str] = set()
    while term_id in folded_into:
        if term_id in seen:
            raise RuntimeError(f"cycle in representative mapping at {term_id}")
        seen.add(term_id)
        term_id = folded_into[term_id]
    return term_id


def clean_go(go_rows: list[dict[str, str]], go_parents: dict[str, set[str]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, str]]:
    memo: dict[str, set[str]] = {}
    terms: dict[str, dict[str, Any]] = {}
    for row in go_rows:
        term_id = row["term_id"]
        genes = symbols(row.get("intersection_genes", ""))
        if not genes:
            raise RuntimeError(f"GO term has no query-gene membership: {term_id}")
        if term_id not in go_parents:
            raise RuntimeError(f"GO term missing from go-basic ontology: {term_id}")
        theme_id, theme_label, _score = theme_for(row["term_name"])
        terms[term_id] = {
            **row,
            "genes": genes,
            "ancestors": ancestor_closure(term_id, go_parents, memo),
            "depth": len(ancestor_closure(term_id, go_parents, memo)),
            "term_size_int": int_value(row.get("term_size")),
            "theme_id": theme_id,
            "theme_label": theme_label,
        }

    order = sorted(
        terms,
        key=lambda term_id: (-terms[term_id]["depth"], terms[term_id]["term_size_int"], fdr(terms[term_id]), term_id),
    )
    retained: list[str] = []
    folded_into: dict[str, str] = {}
    fold_similarity: dict[str, float] = {}
    for term_id in order:
        current = terms[term_id]
        candidates: list[tuple[float, str]] = []
        for kept_id in retained:
            kept = terms[kept_id]
            is_related = kept_id in current["ancestors"] or term_id in kept["ancestors"]
            similarity = jaccard(current["genes"], kept["genes"])
            if is_related and similarity >= ANCESTOR_JACCARD:
                candidates.append((similarity, kept_id))
        if candidates:
            similarity, kept_id = max(candidates, key=lambda item: (item[0], -terms[item[1]]["term_size_int"], item[1]))
            folded_into[term_id] = kept_id
            fold_similarity[term_id] = similarity
        else:
            retained.append(term_id)

    for term_id in list(folded_into):
        folded_into[term_id] = final_representative(folded_into[term_id], folded_into)

    cleaned_rows: list[dict[str, Any]] = []
    representatives: list[dict[str, Any]] = []
    for term_id in sorted(terms):
        term = terms[term_id]
        representative = final_representative(term_id, folded_into)
        retained_here = representative == term_id
        row = {
            "source": "GO:BP",
            "background": BACKGROUND,
            "term_id": term_id,
            "term_name": term["term_name"],
            "adjusted_p": term["gprofiler_adjusted_p"],
            "intersection_size": term["intersection_size"],
            "intersection_genes": ";".join(sorted(term["genes"])),
            "go_ancestor_count": term["depth"],
            "theme_id": term["theme_id"],
            "theme_label": term["theme_label"],
            "status": "retained_representative" if retained_here else "folded",
            "representative_term_id": representative,
            "fold_similarity": "" if retained_here else f"{fold_similarity[term_id]:.6f}",
        }
        cleaned_rows.append(row)
        if retained_here:
            representatives.append(
                {
                    **row,
                    "folded_term_count": sum(1 for value in folded_into.values() if value == term_id),
                }
            )
    return cleaned_rows, representatives, folded_into


def build_aligned_terms(representatives: list[dict[str, Any]], rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    aligned: list[dict[str, Any]] = []
    for row in representatives:
        aligned.append(
            {
                "source": "GO:BP",
                "background": BACKGROUND,
                "theme_id": row["theme_id"],
                "theme_label": row["theme_label"],
                "term_id": row["term_id"],
                "term_name": row["term_name"],
                "adjusted_p": row["adjusted_p"],
                "intersection_size": row["intersection_size"],
                "intersection_genes": row["intersection_genes"],
                "alignment_basis": "GO-DAG-cleaned representative",
            }
        )
    for row in rows:
        theme_id, theme_label, _score = theme_for(row["term_name"])
        aligned.append(
            {
                "source": row["source"],
                "background": row["background"],
                "theme_id": theme_id,
                "theme_label": theme_label,
                "term_id": row["term_id"],
                "term_name": row["term_name"],
                "adjusted_p": row["gprofiler_adjusted_p"],
                "intersection_size": row["intersection_size"],
                "intersection_genes": row["intersection_genes"],
                "alignment_basis": "frozen term-name theme rule",
            }
        )
    return aligned


def build_theme_alignment(aligned: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, dict[str, dict[str, set[str]]]]]:
    by_theme: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for row in aligned:
        by_theme[row["theme_id"]][row["source"]].append(row)
    source_genes: dict[str, dict[str, dict[str, set[str]]]] = defaultdict(lambda: defaultdict(dict))
    summary: list[dict[str, Any]] = []
    theme_ids = [theme_id for theme_id, _, _ in THEMES] + ["other"]
    for theme_id in theme_ids:
        sources = by_theme.get(theme_id, {})
        gene_sets: dict[str, set[str]] = {}
        for source in ("GO:BP", "REAC", "KEGG"):
            gene_sets[source] = set().union(*(symbols(row["intersection_genes"]) for row in sources.get(source, []))) if sources.get(source) else set()
        source_genes[theme_id] = {"genes": gene_sets}
        union = set().union(*gene_sets.values()) if gene_sets else set()
        shared_go_reac = gene_sets["GO:BP"] & gene_sets["REAC"]
        shared_go_kegg = gene_sets["GO:BP"] & gene_sets["KEGG"]
        shared_reac_kegg = gene_sets["REAC"] & gene_sets["KEGG"]
        shared_all = set.intersection(*[values for values in gene_sets.values() if values]) if sum(bool(values) for values in gene_sets.values()) == 3 else set()
        summary.append(
            {
                "theme_id": theme_id,
                "theme_label": THEME_LABEL.get(theme_id, "Other / mixed biological processes"),
                "GO_representative_count": len(sources.get("GO:BP", [])),
                "Reactome_term_count": len(sources.get("REAC", [])),
                "KEGG_term_count": len(sources.get("KEGG", [])),
                "GO_gene_count": len(gene_sets["GO:BP"]),
                "Reactome_gene_count": len(gene_sets["REAC"]),
                "KEGG_gene_count": len(gene_sets["KEGG"]),
                "GO_Reactome_shared_gene_count": len(shared_go_reac),
                "GO_KEGG_shared_gene_count": len(shared_go_kegg),
                "Reactome_KEGG_shared_gene_count": len(shared_reac_kegg),
                "three_source_shared_gene_count": len(shared_all),
                "union_gene_count": len(union),
                "GO_Reactome_union_Jaccard": f"{len(shared_go_reac) / len(gene_sets['GO:BP'] | gene_sets['REAC']):.4f}" if gene_sets["GO:BP"] and gene_sets["REAC"] else "",
            }
        )
    return summary, source_genes


def build_driver_rows(aligned: list[dict[str, Any]]) -> list[dict[str, Any]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    term_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in aligned:
        theme_id = row["theme_id"]
        genes = symbols(row["intersection_genes"])
        counts[theme_id].update(genes)
        term_counts[theme_id].update({gene: 1 for gene in genes})
    rows: list[dict[str, Any]] = []
    for theme_id, counter in counts.items():
        source_term_count = sum(1 for row in aligned if row["theme_id"] == theme_id)
        for rank, (gene, count) in enumerate(counter.most_common(), start=1):
            rows.append(
                {
                    "theme_id": theme_id,
                    "theme_label": THEME_LABEL.get(theme_id, "Other / mixed biological processes"),
                    "gene_symbol": gene,
                    "term_membership_count": count,
                    "fraction_of_theme_terms": f"{count / source_term_count:.4f}" if source_term_count else "",
                    "theme_driver_rank": rank,
                }
            )
    return sorted(rows, key=lambda row: (THEME_ORDER.get(row["theme_id"], 999), int(row["theme_driver_rank"]), row["gene_symbol"]))


def build_target_rows(query_genes: set[str], aligned: list[dict[str, Any]]) -> list[dict[str, Any]]:
    per_gene: dict[str, dict[str, dict[str, set[str]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in aligned:
        source = row["source"]
        theme_id = row["theme_id"]
        genes = symbols(row["intersection_genes"])
        for gene in genes:
            if gene in query_genes:
                per_gene[gene][source]["themes"].add(theme_id)
                per_gene[gene][source]["terms"].add(row["term_id"])
                counts[gene][source] += 1

    rows: list[dict[str, Any]] = []
    informative_theme_ids = set(THEME_LABEL)
    for gene in sorted(query_genes):
        data = per_gene[gene]
        # Keep the unassigned/other bucket for reporting, but never let it
        # count as positive cross-resource biological convergence.
        source_theme_sets = {
            source: {theme for theme in data[source]["themes"] if theme in informative_theme_ids}
            for source in ("GO:BP", "REAC", "KEGG")
        }
        present_sources = [source for source, values in source_theme_sets.items() if values]
        two_source_themes = sorted(
            [theme_id for theme_id in set().union(*source_theme_sets.values()) if sum(theme_id in values for values in source_theme_sets.values()) >= 2],
            key=lambda value: (THEME_ORDER.get(value, 999), value),
        )
        three_source_themes = sorted(
            [theme_id for theme_id in set().union(*source_theme_sets.values()) if all(theme_id in source_theme_sets[source] for source in ("GO:BP", "REAC", "KEGG") if source_theme_sets[source]) and len(present_sources) == 3],
            key=lambda value: (THEME_ORDER.get(value, 999), value),
        )
        same_theme_union = set().union(*source_theme_sets.values()) if source_theme_sets else set()
        biological_key = (
            -len(three_source_themes),
            -len(two_source_themes),
            -len(present_sources),
            -counts[gene]["REAC"],
            -counts[gene]["GO:BP"],
            -counts[gene]["KEGG"],
            -len(same_theme_union),
            gene,
        )
        rows.append(
            {
                "gene_symbol": gene,
                "GO_representative_count": counts[gene]["GO:BP"],
                "Reactome_term_count": counts[gene]["REAC"],
                "KEGG_term_count": counts[gene]["KEGG"],
                "total_aligned_term_count": sum(counts[gene].values()),
                "source_count": len(present_sources),
                "same_theme_two_source_count": len(two_source_themes),
                "same_theme_three_source_count": len(three_source_themes),
                "same_theme_labels": "; ".join(THEME_LABEL.get(theme, "Other / mixed biological processes") for theme in two_source_themes),
                "source_status": "+".join(present_sources) if present_sources else "none",
                "GO_theme_count": len(source_theme_sets["GO:BP"]),
                "Reactome_theme_count": len(source_theme_sets["REAC"]),
                "KEGG_theme_count": len(source_theme_sets["KEGG"]),
                "theme_breadth": len(same_theme_union),
                "_biological_key": biological_key,
            }
        )
    rows.sort(key=lambda row: row["_biological_key"])
    for rank, row in enumerate(rows, start=1):
        row["biological_priority_rank"] = rank
        row.pop("_biological_key", None)
    return rows


def main() -> None:
    if not ENRICHMENT.exists() or not QUERY.exists():
        raise FileNotFoundError("fresh Gelsemium reference enrichment outputs are missing")
    if not GO_CACHE.exists() or GO_CACHE.stat().st_size < 1000:
        download_go_obo(GO_CACHE)
    go_parents, go_names = parse_go_obo(GO_CACHE)
    rows = [row for row in read_csv(ENRICHMENT) if row["background"] == BACKGROUND and fdr(row) < FDR_CUTOFF]
    go_rows = [row for row in rows if row["source"] == "GO:BP"]
    reactome_rows = [row for row in rows if row["source"] == "REAC"]
    kegg_rows = [row for row in rows if row["source"] == "KEGG"]
    expected = {"GO:BP": 202, "REAC": 50, "KEGG": 13}
    actual = {"GO:BP": len(go_rows), "REAC": len(reactome_rows), "KEGG": len(kegg_rows)}
    if actual != expected:
        raise RuntimeError(f"unexpected significant-term counts: {actual}; expected {expected}")
    query_genes = {row["gene_symbol"].strip().upper() for row in read_csv(QUERY) if row.get("gene_symbol", "").strip()}
    if len(query_genes) != 41:
        raise RuntimeError(f"expected 41 fresh query genes, found {len(query_genes)}")
    if not all(term_id in go_names for term_id in {row["term_id"] for row in go_rows}):
        raise RuntimeError("one or more current GO terms are absent from go-basic.obo")

    cleaned_rows, representatives, folded_into = clean_go(go_rows, go_parents)
    aligned = build_aligned_terms(representatives, reactome_rows + kegg_rows)
    theme_rows, _theme_gene_sets = build_theme_alignment(aligned)
    go_driver_rows = build_driver_rows(
        [row for row in aligned if row["source"] == "GO:BP"]
    )
    driver_rows = build_driver_rows(aligned)
    target_rows = build_target_rows(query_genes, aligned)

    clean_dir = OUTPUTS / "go_bp_cleaning_1893"
    align_dir = OUTPUTS / "go_reactome_kegg_aligned_1893"
    target_dir = OUTPUTS / "target_prioritization_1893"
    clean_dir.mkdir(parents=True, exist_ok=True)
    align_dir.mkdir(parents=True, exist_ok=True)
    target_dir.mkdir(parents=True, exist_ok=True)

    cleaned_fields = list(cleaned_rows[0].keys())
    representative_fields = list(representatives[0].keys())
    write_csv(clean_dir / "go_bp_202_terms_cleaning.csv", cleaned_rows, cleaned_fields)
    write_csv(clean_dir / "go_bp_clean_representatives.csv", representatives, representative_fields)
    write_csv(clean_dir / "go_bp_clean_theme_driver_contribution.csv", go_driver_rows, list(go_driver_rows[0].keys()))
    write_csv(align_dir / "go_reactome_kegg_aligned_terms.csv", aligned, list(aligned[0].keys()))
    write_csv(align_dir / "go_reactome_kegg_theme_alignment.csv", theme_rows, list(theme_rows[0].keys()))
    write_csv(align_dir / "go_reactome_kegg_driver_contribution.csv", driver_rows, list(driver_rows[0].keys()))

    target_fields = [
        "biological_priority_rank",
        "gene_symbol",
        "GO_representative_count",
        "Reactome_term_count",
        "KEGG_term_count",
        "total_aligned_term_count",
        "source_count",
        "same_theme_two_source_count",
        "same_theme_three_source_count",
        "same_theme_labels",
        "source_status",
        "GO_theme_count",
        "Reactome_theme_count",
        "KEGG_theme_count",
        "theme_breadth",
    ]
    write_csv(target_dir / "go_reactome_kegg_target_prioritization.csv", target_rows, target_fields)
    shortlist = [row for row in target_rows if int(row["same_theme_two_source_count"]) > 0]
    write_csv(target_dir / "cross_resource_target_shortlist.csv", shortlist, target_fields)

    top_drivers = driver_rows[:]
    summary_lines = [
        "# Fresh DINP--CRC GO cleaning and target prioritization",
        "",
        f"Generated (UTC): {now_utc()}",
        "",
        "## Frozen scope",
        "",
        "This analysis uses only the published Supplementary Table S2 CRC reference set, the frozen 93-gene DINP multi-source set, and their fresh 41-gene symbol-level intersection. It does not reuse the old 881-gene background or 18-gene query.",
        "",
        f"- Fresh query genes: **{len(query_genes)}**",
        f"- Significant GO:BP terms before cleaning: **{len(go_rows)}**",
        f"- Significant Reactome terms: **{len(reactome_rows)}**",
        f"- Significant KEGG terms: **{len(kegg_rows)}**",
        f"- Retained GO representatives: **{len(representatives)}**",
        f"- Folded GO terms: **{len(folded_into)}**",
        "",
        "## GO cleaning",
        "",
        f"GO:BP terms were ordered by GO-DAG depth and reduced only for ancestor/descendant pairs with query-gene Jaccard >= **{ANCESTOR_JACCARD:.2f}**. Non-ancestor terms were not collapsed solely because they shared genes.",
        "",
        "The original enrichment P values and FDR values remain unchanged. The cleaned representatives are an interpretation/readability layer, not a new statistical family.",
        "",
        "## Cross-resource themes",
        "",
        "GO representatives, Reactome, and KEGG were assigned to the same frozen term-name themes. Theme overlap is descriptive and does not treat GO, Reactome, and KEGG as independent tests. The unassigned/other bucket is reported for completeness but is excluded from positive target prioritization.",
        "",
        "| Theme | GO reps | Reactome | KEGG | GO–Reactome shared genes | GO–KEGG shared genes | 3-source shared genes | Union genes |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in theme_rows:
        if int(row["GO_representative_count"]) + int(row["Reactome_term_count"]) + int(row["KEGG_term_count"]) == 0:
            continue
        summary_lines.append(
            f"| {row['theme_label']} | {row['GO_representative_count']} | {row['Reactome_term_count']} | {row['KEGG_term_count']} | {row['GO_Reactome_shared_gene_count']} | {row['GO_KEGG_shared_gene_count']} | {row['three_source_shared_gene_count']} | {row['union_gene_count']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Target prioritization",
            "",
            "The ranking is lexicographic, not a fitted composite score: three-source same-theme support, two-source same-theme support, number of supporting resources, Reactome recurrence, GO recurrence, KEGG recurrence, then theme breadth and gene symbol. This prioritizes reproducible cross-resource anchors without claiming causality.",
            "",
            "| Rank | Gene | GO | Reactome | KEGG | Resources | Same-theme 2-source | Same-theme labels |",
            "|---:|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in target_rows[:15]:
        summary_lines.append(
            f"| {row['biological_priority_rank']} | **{row['gene_symbol']}** | {row['GO_representative_count']} | {row['Reactome_term_count']} | {row['KEGG_term_count']} | {row['source_count']} | {row['same_theme_two_source_count']} | {row['same_theme_labels'] or '—'} |"
        )
    summary_lines.extend(
        [
            "",
            "## Boundary",
            "",
            "A recurrent gene is a pathway-context candidate, not a proven DINP target. Database overlap, pathway recurrence, and target prioritization do not establish direct binding, directionality, causality, or in-vivo mediation. Any docking or structural follow-up must remain a separate hypothesis-generating layer.",
            "",
            "## Files",
            "",
            "- `go_bp_cleaning_1893/`: GO:BP term-level cleaning and representatives.",
            "- `go_reactome_kegg_aligned_1893/`: descriptive cross-resource theme alignment and drivers.",
            "- `target_prioritization_1893/`: current 41-gene target ranking and cross-resource shortlist.",
        ]
    )
    summary_path = OUTPUTS / "go_cleaning_target_prioritization_summary.md"
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest = {
        "generated_utc": now_utc(),
        "script": str(Path(__file__).resolve()),
        "input_enrichment": str(ENRICHMENT.resolve()),
        "input_enrichment_sha256": sha256_file(ENRICHMENT),
        "input_query": str(QUERY.resolve()),
        "input_query_sha256": sha256_file(QUERY),
        "background": BACKGROUND,
        "query_gene_count": len(query_genes),
        "significant_term_counts": actual,
        "go_ontology_url": GO_URL,
        "go_ontology_cache": str(GO_CACHE.resolve()),
        "go_ontology_sha256": sha256_file(GO_CACHE),
        "go_ontology_term_count": len(go_parents),
        "ancestor_jaccard_threshold": ANCESTOR_JACCARD,
        "go_representative_count": len(representatives),
        "go_folded_term_count": len(folded_into),
        "theme_count": len(theme_rows),
        "target_gene_count": len(target_rows),
        "cross_resource_shortlist_count": len(shortlist),
        "statistical_boundary": "No enrichment P/FDR recalculation; cleaning and prioritization are descriptive post-enrichment layers.",
        "outputs": {
            "summary": str(summary_path.resolve()),
            "go_cleaning": str((clean_dir / "go_bp_202_terms_cleaning.csv").resolve()),
            "go_representatives": str((clean_dir / "go_bp_clean_representatives.csv").resolve()),
            "aligned_terms": str((align_dir / "go_reactome_kegg_aligned_terms.csv").resolve()),
            "theme_alignment": str((align_dir / "go_reactome_kegg_theme_alignment.csv").resolve()),
            "target_prioritization": str((target_dir / "go_reactome_kegg_target_prioritization.csv").resolve()),
            "target_shortlist": str((target_dir / "cross_resource_target_shortlist.csv").resolve()),
        },
    }
    (OUTPUTS / "go_cleaning_target_prioritization_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print("FRESH DINP--CRC GO CLEANING + TARGET PRIORITIZATION: PASS")
    print(f"Query genes: {len(query_genes)}")
    print(f"GO:BP terms: {len(go_rows)} -> representatives: {len(representatives)}")
    print(f"Reactome terms: {len(reactome_rows)}; KEGG terms: {len(kegg_rows)}")
    print(f"Cross-resource target shortlist: {len(shortlist)}")
    print(f"Summary: {summary_path.resolve()}")
    print(f"Target table: {(target_dir / 'go_reactome_kegg_target_prioritization.csv').resolve()}")


if __name__ == "__main__":
    main()
