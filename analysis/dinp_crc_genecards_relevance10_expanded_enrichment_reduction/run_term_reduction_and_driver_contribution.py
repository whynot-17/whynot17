#!/usr/bin/env python3
"""Reduce the 106 primary enrichment terms and quantify 18-gene drivers.

This is an additive post-enrichment audit.  It does not change the frozen
18-gene query, the 881-gene primary background, or the primary enrichment
threshold.  The input term universe is exactly the 106 terms with global
BH-FDR < 0.05 in the existing primary GO:BP/KEGG/Reactome output.

Term-to-gene membership is obtained from a second, provenance-preserving
g:Profiler request with ``no_evidences=false``.  The API returns one evidence
list per query gene; a non-empty evidence list means that the query gene is
in that term.  Reduction is deliberately descriptive: terms are grouped by
average-linkage similarity using query-gene Jaccard overlap, normalized term
tokens, and GO/Reactome parent linkage.  It is not a new enrichment test.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
ENRICHMENT_DIR = REPO_ROOT / "analysis" / "dinp_crc_genecards_relevance10_expanded_enrichment"
INPUT_TERM_FILE = ENRICHMENT_DIR / "outputs" / "dinp_crc_enrichment_primary_genecards881.csv"
INPUT_QUERY_FILE = ENRICHMENT_DIR / "outputs" / "dinp_crc_18_gene_list.csv"
INPUT_INTERSECTION_FILE = (
    REPO_ROOT
    / "analysis"
    / "dinp_crc_multi_database_target_convergence"
    / "outputs"
    / "dinp_crc_genecards_relevance10_expanded_intersection.csv"
)
OUT_DIR = HERE / "outputs"
RAW_DIR = OUT_DIR / "raw_gprofiler"
API_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "KEGG", "REAC"]
DEFAULT_THRESHOLD = 0.58
USER_AGENT = "whynot17-dinp-crc-term-reduction/1.0"


THEME_RULES: list[tuple[str, str, list[str]]] = [
    (
        "lipid_eicosanoid",
        "Lipid / fatty-acid / eicosanoid metabolism",
        [
            r"fatty acid",
            r"lipid",
            r"prostagland",
            r"prostanoid",
            r"icosanoid",
            r"steroid",
            r"cholesterol",
            r"spm",
            r"dha",
            r"monocarboxylic acid",
            r"carboxylic acid",
            r"oxoacid",
            r"organic acid",
            r"olefinic",
        ],
    ),
    (
        "inflammation_defense",
        "Inflammation / defense / stimulus response",
        [
            r"inflamm",
            r"defense",
            r"response to stimulus",
            r"response to stress",
            r"biotic",
            r"interspecies interaction",
            r"other organism",
            r"external stimulus",
        ],
    ),
    (
        "chemical_oxygen_stress",
        "Chemical / oxygen / abiotic stress response",
        [
            r"chemical stimulus",
            r"chemical stress",
            r"response to chemical",
            r"oxygen",
            r"hypoxia",
            r"uv-a",
            r"abiotic",
            r"decreased oxygen",
        ],
    ),
    (
        "apoptosis_cell_death",
        "Apoptosis / programmed cell death",
        [r"apopt", r"cell death", r"programmed cell death"],
    ),
    (
        "mirna_regulation",
        "miRNA transcription / metabolism",
        [r"mirna"],
    ),
    (
        "signaling_receptor",
        "Signaling / receptor / hormone response",
        [r"signal", r"signaling", r"receptor", r"hormone", r"cell communication"],
    ),
    (
        "matrix_migration",
        "Matrix remodeling / cell migration",
        [r"matrix metalloproteinase", r"cell migration", r"extracellular matrix"],
    ),
    (
        "development_reproduction",
        "Development / reproduction / circulation",
        [
            r"development",
            r"embryo",
            r"pregnancy",
            r"reproductive",
            r"organ development",
            r"circulation",
            r"circulatory",
            r"blood",
            r"proliferation",
            r"multicellular organism",
            r"biological regulation",
            r"positive regulation of cellular process",
            r"positive regulation of biological process",
        ],
    ),
    (
        "general_metabolism",
        "General metabolism / biosynthesis / storage",
        [r"metabolic", r"biosynthetic", r"catabolic", r"storage", r"biological quality"],
    ),
    (
        "ontology_umbrella",
        "Ontology umbrella / root term",
        [r"reactome root term"],
    ),
]

STOPWORDS = {
    "a",
    "an",
    "and",
    "cell",
    "cellular",
    "of",
    "organism",
    "process",
    "regulation",
    "response",
    "the",
    "to",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def finite_float(value: Any, default: float = float("nan")) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def unique_symbols(values: list[str]) -> list[str]:
    return sorted({str(value).strip().upper() for value in values if str(value).strip()})


def request_gprofiler(payload: dict[str, Any]) -> tuple[dict[str, Any], bytes]:
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        API_URL,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        method="POST",
    )
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                raw = response.read()
            result = json.loads(raw.decode("utf-8"))
            if "result" not in result or "meta" not in result:
                raise ValueError("g:Profiler response lacks result/meta")
            return result, raw
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ValueError) as exc:
            last_error = exc
            if attempt < 3:
                time.sleep(10 * attempt)
    raise RuntimeError(f"g:Profiler request failed after 3 attempts: {last_error}")


def normalized_tokens(text: str) -> set[str]:
    tokens = set(re.findall(r"[a-z0-9]+", text.lower()))
    return {token for token in tokens if token not in STOPWORDS and len(token) > 1}


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def theme_for_term(name: str) -> tuple[str, str, int]:
    lowered = name.lower()
    scored: list[tuple[int, int, str, str]] = []
    for order, (theme_id, label, patterns) in enumerate(THEME_RULES):
        score = sum(1 for pattern in patterns if re.search(pattern, lowered))
        if score:
            scored.append((score, -order, theme_id, label))
    if not scored:
        return "other", "Other / mixed biological processes", 0
    score, neg_order, theme_id, label = max(scored)
    return theme_id, label, score


def term_similarity(left: dict[str, Any], right: dict[str, Any]) -> float:
    gene_score = jaccard(left["genes"], right["genes"])
    token_score = jaccard(left["tokens"], right["tokens"])
    parent_score = jaccard(left["parents"], right["parents"])
    direct_parent = bool(
        left["term_id"] in right["parents"]
        or right["term_id"] in left["parents"]
    )
    # Gene overlap carries the most information because the input query is only
    # 18 genes. Parent linkage prevents GO child/parent terms from being split
    # merely because their names differ.
    return 0.65 * gene_score + 0.20 * token_score + 0.10 * parent_score + 0.05 * float(direct_parent)


def average_linkage_clusters(terms: list[dict[str, Any]], threshold: float) -> list[list[int]]:
    clusters: list[list[int]] = [[index] for index in range(len(terms))]
    while True:
        best: tuple[float, int, int] | None = None
        for i in range(len(clusters)):
            for j in range(i):
                scores = [
                    term_similarity(terms[left], terms[right])
                    for left in clusters[i]
                    for right in clusters[j]
                ]
                score = sum(scores) / len(scores) if scores else 0.0
                candidate = (score, -i, -j)
                if best is None or candidate > best:
                    best = candidate
                    best_pair = (i, j)
        if best is None or best[0] < threshold:
            break
        i, j = best_pair
        merged = sorted(clusters[i] + clusters[j])
        clusters[i] = merged
        clusters.pop(j)
    return sorted(clusters, key=lambda cluster: (cluster[0], len(cluster)))


def build_term_gene_map(result: dict[str, Any], expected_genes: list[str]) -> dict[str, set[str]]:
    query_meta = result.get("meta", {}).get("query_metadata", {}).get("queries", {})
    if not query_meta:
        raise ValueError("g:Profiler response lacks query_metadata.queries")
    query_key = sorted(query_meta)[0]
    query_genes = [str(value).upper() for value in query_meta[query_key]]
    if query_genes != expected_genes:
        raise ValueError(f"Query gene order changed unexpectedly: {query_genes}")
    membership: dict[str, set[str]] = {}
    for item in result.get("result", []):
        native = str(item.get("native", ""))
        evidences = item.get("intersections", [])
        if not native or not isinstance(evidences, list) or len(evidences) != len(query_genes):
            continue
        genes = {gene for gene, evidence in zip(query_genes, evidences) if evidence}
        membership[native] = genes
    return membership


def make_payload(query: list[str], background: list[str]) -> dict[str, Any]:
    return {
        "organism": "hsapiens",
        "query": query,
        "background": background,
        "sources": SOURCES,
        "domain_scope": "custom",
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "all_results": True,
        "no_evidences": False,
        "ordered": False,
        "combined": False,
        "measure_underrepresentation": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cluster-threshold", type=float, default=DEFAULT_THRESHOLD)
    args = parser.parse_args()
    if not 0 < args.cluster_threshold < 1:
        raise ValueError("--cluster-threshold must be between 0 and 1")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    primary_rows = [
        row
        for row in read_csv(INPUT_TERM_FILE)
        if parse_bool(row.get("BH_FDR_global_significant", ""))
    ]
    if len(primary_rows) != 106:
        raise ValueError(f"Expected exactly 106 primary significant terms, found {len(primary_rows)}")

    query_rows = read_csv(INPUT_QUERY_FILE)
    query_genes = unique_symbols([row.get("gene_symbol", "") for row in query_rows])
    if len(query_genes) != 18:
        raise ValueError(f"Expected 18 query genes, found {len(query_genes)}")

    intersection_rows = read_csv(INPUT_INTERSECTION_FILE)
    genecard_lookup = {
        str(row.get("gene_symbol", "")).upper(): row for row in intersection_rows
    }
    background = unique_symbols(
        [
            row.get("gene_symbol", "")
            for row in read_csv(
                REPO_ROOT
                / "analysis"
                / "dinp_crc_multi_database_target_convergence"
                / "outputs"
                / "dinp_crc_genecards_relevance10_intersection.csv"
            )
        ]
    )
    # The 881-gene background is the source registry used by the existing
    # enrichment script. Reconstruct it directly from the archived top-2000
    # GeneCards export to make the request provenance self-contained.
    cards_path = (
        REPO_ROOT
        / "analysis"
        / "dinp_crc_multi_database_target_convergence"
        / "outputs"
        / "source_records"
        / "genecards_crc_top2000_reference.csv"
    )
    if cards_path.exists():
        background = unique_symbols(
            [
                row.get("gene_symbol") or row.get("GeneSymbol", "")
                for row in read_csv(cards_path)
                if finite_float(row.get("RelevanceScore")) >= 10
            ]
        )
    if len(background) != 881:
        raise ValueError(f"Expected 881 primary background genes, found {len(background)}")
    if not set(query_genes).issubset(background):
        raise ValueError("Query genes are not a subset of the 881-gene background")

    payload = make_payload(query_genes, background)
    payload_path = OUT_DIR / "request_payload_primary_106_intersections.json"
    raw_path = RAW_DIR / "gprofiler_primary_106_intersections.json"
    if raw_path.exists() and raw_path.stat().st_size > 0 and payload_path.exists():
        raw = raw_path.read_bytes()
        result = json.loads(raw.decode("utf-8"))
        cached_payload = json.loads(payload_path.read_text(encoding="utf-8"))
        if cached_payload != payload:
            raise ValueError("Cached payload does not match the frozen primary request")
    else:
        result, raw = request_gprofiler(payload)
        raw_path.write_bytes(raw)
        payload_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    membership = build_term_gene_map(result, query_genes)
    if len(membership) < len(primary_rows):
        missing = sorted({row["term_id"] for row in primary_rows} - set(membership))
        raise ValueError(f"Missing term-to-gene evidence for {len(missing)} primary terms: {missing[:5]}")

    terms: list[dict[str, Any]] = []
    for row in primary_rows:
        term_id = row["term_id"]
        genes = membership[term_id]
        expected_overlap = int(row["intersection_size"])
        if len(genes) != expected_overlap:
            raise ValueError(
                f"Intersection mismatch for {term_id}: API evidence={len(genes)} vs stored={expected_overlap}"
            )
        theme_id, theme_label, theme_score = theme_for_term(row["term_name"])
        terms.append(
            {
                **row,
                "genes": genes,
                "tokens": normalized_tokens(row["term_name"]),
                "parents_set": set(filter(None, row.get("parents", "").split(";"))),
                "parents": set(filter(None, row.get("parents", "").split(";"))),
                "term_theme_id": theme_id,
                "term_theme_label": theme_label,
                "term_theme_score": theme_score,
            }
        )

    clusters = average_linkage_clusters(terms, args.cluster_threshold)
    term_to_cluster: dict[int, str] = {}
    cluster_rows: list[dict[str, Any]] = []
    cluster_term_rows: list[dict[str, Any]] = []
    for cluster_index, indices in enumerate(clusters, start=1):
        cluster_id = f"theme_{cluster_index:02d}"
        for index in indices:
            term_to_cluster[index] = cluster_id
        theme_counts = Counter(terms[index]["term_theme_id"] for index in indices)
        theme_id = sorted(theme_counts, key=lambda key: (-theme_counts[key], key))[0]
        theme_label = next(
            terms[index]["term_theme_label"]
            for index in indices
            if terms[index]["term_theme_id"] == theme_id
        )
        union_genes = set().union(*(terms[index]["genes"] for index in indices))
        source_counts = Counter(terms[index]["source"] for index in indices)
        representative_index = min(
            indices,
            key=lambda index: (
                finite_float(terms[index]["BH_FDR_global_source_family"], 1.0),
                int(terms[index]["term_size"]),
                terms[index]["term_id"],
            ),
        )
        cluster_rows.append(
            {
                "cluster_id": cluster_id,
                "theme_label": theme_label,
                "term_count": len(indices),
                "source_count": len(source_counts),
                "sources": ";".join(sorted(source_counts)),
                "union_gene_count": len(union_genes),
                "union_genes": ";".join(sorted(union_genes)),
                "representative_term_id": terms[representative_index]["term_id"],
                "representative_term_name": terms[representative_index]["term_name"],
                "theme_term_counts": ";".join(
                    f"{key}={theme_counts[key]}" for key in sorted(theme_counts)
                ),
            }
        )
        for index in sorted(indices, key=lambda value: (finite_float(terms[value]["BH_FDR_global_source_family"], 1.0), terms[value]["term_id"])):
            term = terms[index]
            cluster_term_rows.append(
                {
                    "cluster_id": cluster_id,
                    "theme_label": theme_label,
                    "term_id": term["term_id"],
                    "source": term["source"],
                    "term_name": term["term_name"],
                    "global_fdr": term["BH_FDR_global_source_family"],
                    "intersection_size": term["intersection_size"],
                    "term_theme_id": term["term_theme_id"],
                    "term_theme_label": term["term_theme_label"],
                    "intersection_genes": ";".join(sorted(term["genes"])),
                }
            )

    cluster_lookup = {row["cluster_id"]: row for row in cluster_rows}
    # The fine-grained similarity clusters above preserve local redundancy
    # structure.  For interpretation, also provide a disjoint set of broad
    # theme families using the fixed term-name rule set.  This avoids calling
    # dozens of GO descendants separate biological stories while retaining the
    # underlying term-level assignment in the output files.
    family_order = {theme_id: order for order, (theme_id, _label, _patterns) in enumerate(THEME_RULES)}
    family_order["other"] = len(THEME_RULES)
    family_groups: dict[str, list[int]] = defaultdict(list)
    for index, term in enumerate(terms):
        family_groups[term["term_theme_id"]].append(index)
    family_ids = {
        theme_id: f"family_{position:02d}"
        for position, theme_id in enumerate(
            sorted(family_groups, key=lambda value: (family_order.get(value, 999), value)),
            start=1,
        )
    }
    term_to_family = {
        index: family_ids[term["term_theme_id"]]
        for index, term in enumerate(terms)
    }
    family_rows: list[dict[str, Any]] = []
    family_term_rows: list[dict[str, Any]] = []
    for theme_id in sorted(family_groups, key=lambda value: (family_order.get(value, 999), value)):
        indices = family_groups[theme_id]
        label = terms[indices[0]]["term_theme_label"]
        union_genes = set().union(*(terms[index]["genes"] for index in indices))
        source_counts = Counter(terms[index]["source"] for index in indices)
        representative_index = min(
            indices,
            key=lambda index: (
                finite_float(terms[index]["BH_FDR_global_source_family"], 1.0),
                int(terms[index]["term_size"]),
                terms[index]["term_id"],
            ),
        )
        family_rows.append(
            {
                "family_id": family_ids[theme_id],
                "theme_id": theme_id,
                "theme_label": label,
                "term_count": len(indices),
                "source_count": len(source_counts),
                "sources": ";".join(sorted(source_counts)),
                "union_gene_count": len(union_genes),
                "union_genes": ";".join(sorted(union_genes)),
                "representative_term_id": terms[representative_index]["term_id"],
                "representative_term_name": terms[representative_index]["term_name"],
            }
        )
        for index in sorted(
            indices,
            key=lambda value: (
                finite_float(terms[value]["BH_FDR_global_source_family"], 1.0),
                terms[value]["term_id"],
            ),
        ):
            term = terms[index]
            family_term_rows.append(
                {
                    "family_id": family_ids[theme_id],
                    "theme_id": theme_id,
                    "theme_label": label,
                    "term_id": term["term_id"],
                    "source": term["source"],
                    "term_name": term["term_name"],
                    "global_fdr": term["BH_FDR_global_source_family"],
                    "intersection_size": term["intersection_size"],
                    "intersection_genes": ";".join(sorted(term["genes"])),
                }
            )

    gene_term_counts: Counter[str] = Counter()
    gene_cluster_sets: dict[str, set[str]] = defaultdict(set)
    gene_family_sets: dict[str, set[str]] = defaultdict(set)
    gene_family_term_counts: dict[str, Counter[str]] = defaultdict(Counter)
    gene_source_counts: dict[str, Counter[str]] = defaultdict(Counter)
    cluster_gene_term_counts: Counter[tuple[str, str]] = Counter()
    for index, term in enumerate(terms):
        cluster_id = term_to_cluster[index]
        family_id = term_to_family[index]
        for gene in term["genes"]:
            gene_term_counts[gene] += 1
            gene_cluster_sets[gene].add(cluster_id)
            gene_family_sets[gene].add(family_id)
            gene_family_term_counts[gene][family_id] += 1
            gene_source_counts[gene][term["source"]] += 1
            cluster_gene_term_counts[(cluster_id, gene)] += 1

    driver_rows: list[dict[str, Any]] = []
    for gene in query_genes:
        card = genecard_lookup.get(gene, {})
        family_counts = gene_family_term_counts[gene]
        driver_rows.append(
            {
                "gene_symbol": gene,
                "significant_term_count": gene_term_counts[gene],
                "fraction_of_106_terms": gene_term_counts[gene] / len(terms),
                "redundancy_cluster_count": len(gene_cluster_sets[gene]),
                "fraction_of_redundancy_clusters": len(gene_cluster_sets[gene]) / len(clusters),
                "theme_family_count": len(gene_family_sets[gene]),
                "fraction_of_theme_families": len(gene_family_sets[gene]) / len(family_rows),
                "GO_BP_term_count": gene_source_counts[gene]["GO:BP"],
                "KEGG_term_count": gene_source_counts[gene]["KEGG"],
                "Reactome_term_count": gene_source_counts[gene]["REAC"],
                "theme_family_contributions": ";".join(
                    f"{family_id}={family_counts[family_id]}"
                    for family_id in sorted(family_counts, key=lambda key: (-family_counts[key], key))
                ),
                "theme_family_membership": ";".join(sorted(gene_family_sets[gene])),
                "redundancy_cluster_membership": ";".join(sorted(gene_cluster_sets[gene])),
                "gene_cards_rank": card.get("gene_cards_rank", ""),
                "gene_cards_relevance_score": card.get("gene_cards_relevance_score", ""),
            }
        )
    driver_rows.sort(key=lambda row: (-int(row["significant_term_count"]), -int(row["theme_family_count"]), row["gene_symbol"]))
    overall_rank = {row["gene_symbol"]: index + 1 for index, row in enumerate(driver_rows)}

    cluster_driver_rows: list[dict[str, Any]] = []
    for cluster in cluster_rows:
        cluster_id = cluster["cluster_id"]
        term_count = int(cluster["term_count"])
        members = [gene for gene in query_genes if cluster_gene_term_counts[(cluster_id, gene)] > 0]
        local_counts = {gene: cluster_gene_term_counts[(cluster_id, gene)] for gene in members}
        max_count = max(local_counts.values(), default=0)
        for gene in sorted(members, key=lambda value: (-local_counts[value], overall_rank[value], value)):
            cluster_driver_rows.append(
                {
                    "cluster_id": cluster_id,
                    "theme_label": cluster["theme_label"],
                    "gene_symbol": gene,
                    "cluster_term_count": local_counts[gene],
                    "fraction_of_cluster_terms": local_counts[gene] / term_count,
                    "cluster_top_driver": local_counts[gene] == max_count,
                    "global_significant_term_count": gene_term_counts[gene],
                    "global_theme_cluster_count": len(gene_cluster_sets[gene]),
                    "overall_driver_rank": overall_rank[gene],
                }
            )

    family_driver_rows: list[dict[str, Any]] = []
    for family in family_rows:
        family_id = family["family_id"]
        family_term_count = int(family["term_count"])
        members = [gene for gene in query_genes if gene_family_term_counts[gene][family_id] > 0]
        local_counts = {gene: gene_family_term_counts[gene][family_id] for gene in members}
        max_count = max(local_counts.values(), default=0)
        for gene in sorted(members, key=lambda value: (-local_counts[value], overall_rank[value], value)):
            family_driver_rows.append(
                {
                    "family_id": family_id,
                    "theme_id": family["theme_id"],
                    "theme_label": family["theme_label"],
                    "gene_symbol": gene,
                    "family_term_count": local_counts[gene],
                    "fraction_of_family_terms": local_counts[gene] / family_term_count,
                    "family_top_driver": local_counts[gene] == max_count,
                    "global_significant_term_count": gene_term_counts[gene],
                    "theme_family_count": len(gene_family_sets[gene]),
                    "overall_driver_rank": overall_rank[gene],
                }
            )

    term_output_rows: list[dict[str, Any]] = []
    for index, term in enumerate(terms):
        term_output_rows.append(
            {
                "cluster_id": term_to_cluster[index],
                "theme_label": cluster_lookup[term_to_cluster[index]]["theme_label"],
                "family_id": term_to_family[index],
                "family_label": next(row["theme_label"] for row in family_rows if row["family_id"] == term_to_family[index]),
                "source": term["source"],
                "term_id": term["term_id"],
                "term_name": term["term_name"],
                "global_fdr": term["BH_FDR_global_source_family"],
                "term_size": term["term_size"],
                "effective_domain_size": term["effective_domain_size"],
                "intersection_size": term["intersection_size"],
                "term_theme_id": term["term_theme_id"],
                "term_theme_label": term["term_theme_label"],
                "parents": term["parents"],
                "intersection_genes": ";".join(sorted(term["genes"])),
            }
        )

    write_csv(
        OUT_DIR / "primary_106_terms_with_genes.csv",
        term_output_rows,
        list(term_output_rows[0].keys()),
    )
    write_csv(
        OUT_DIR / "primary_theme_clusters.csv",
        cluster_rows,
        list(cluster_rows[0].keys()),
    )
    write_csv(
        OUT_DIR / "primary_cluster_term_membership.csv",
        cluster_term_rows,
        list(cluster_term_rows[0].keys()),
    )
    write_csv(
        OUT_DIR / "primary_18_gene_driver_contribution.csv",
        driver_rows,
        list(driver_rows[0].keys()),
    )
    write_csv(
        OUT_DIR / "primary_cluster_driver_contribution.csv",
        cluster_driver_rows,
        list(cluster_driver_rows[0].keys()),
    )
    write_csv(
        OUT_DIR / "primary_theme_family_term_membership.csv",
        family_term_rows,
        list(family_term_rows[0].keys()),
    )
    write_csv(
        OUT_DIR / "primary_theme_family_driver_contribution.csv",
        family_driver_rows,
        list(family_driver_rows[0].keys()),
    )
    write_csv(
        OUT_DIR / "primary_theme_families.csv",
        family_rows,
        list(family_rows[0].keys()),
    )

    global_top = driver_rows[:8]
    family_driver_lookup: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in family_driver_rows:
        family_driver_lookup[row["family_id"]].append(row)
    summary_lines = [
        "# Primary 106-term redundancy reduction and 18-gene driver contribution",
        "",
        f"Generated (UTC): {utc_now()}",
        "",
        "## Frozen scope",
        "",
        "- Input: the existing primary GO:BP/KEGG/Reactome enrichment output only.",
        "- Query: the frozen 18-gene DINP–CRC intersection; no genes were added or removed.",
        "- Background: the archived GeneCards CRC high-relevance set with RelevanceScore >= 10 (881 genes).",
        "- Primary term universe: exactly 106 terms with global source-family BH-FDR < 0.05.",
        "- The sensitivity background (93 genes; 0 globally significant terms) was not mixed into this reduction.",
        "",
        "## Term-to-gene evidence",
        "",
        "Term membership was reconstructed from a separate g:Profiler request with `no_evidences=false`. The response provides one evidence list per query gene; non-empty evidence denotes membership. All 106 terms passed the validation that the evidence-derived overlap count equals the stored `intersection_size`.",
        "",
        "## Reduction method",
        "",
        f"Average-linkage agglomeration used a fixed similarity threshold of **{args.cluster_threshold:.2f}**. Pairwise similarity was 0.65 × query-gene Jaccard + 0.20 × normalized term-token Jaccard + 0.10 × parent-set Jaccard + 0.05 × direct parent/child linkage. The result is a descriptive redundancy reduction, not a new enrichment test and not a causal model.",
        "",
        f"- 106 significant terms -> **{len(clusters)} fine-grained redundancy clusters** -> **{len(family_rows)} broad theme families**.",
        f"- Sources represented: {', '.join(sorted(set(row['source'] for row in term_output_rows)))}.",
        "",
        "## Broad theme families",
        "",
        "The family layer is the main interpretation view. It is a transparent, disjoint grouping of terms using the fixed term-name rules in the script; the fine-grained similarity clusters remain available in `primary_theme_clusters.csv` and `primary_cluster_term_membership.csv`.",
        "",
        "| Family | Theme | Terms | Sources | Representative term | Union genes |",
        "|---|---|---:|---:|---|---:|",
    ]
    for row in family_rows:
        summary_lines.append(
            f"| {row['family_id']} | {row['theme_label']} | {row['term_count']} | {row['source_count']} | {row['representative_term_name']} | {row['union_gene_count']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Fine-grained redundancy clusters",
            "",
            "| Cluster | Theme | Terms | Sources | Representative term | Union genes |",
            "|---|---|---:|---:|---|---:|",
        ]
    )
    for row in cluster_rows:
        summary_lines.append(
            f"| {row['cluster_id']} | {row['theme_label']} | {row['term_count']} | {row['source_count']} | {row['representative_term_name']} | {row['union_gene_count']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Theme-family driver contribution",
            "",
            "For each family, genes are ranked by the number of significant terms in that family containing the gene. A tie is resolved by the overall contribution rank. `family_top_driver` is descriptive and does not denote a gene-level test.",
            "",
            "| Family | Theme | Top contributing genes (gene:term count) |",
            "|---|---|---|",
        ]
    )
    for family in family_rows:
        family_drivers = family_driver_lookup[family["family_id"]]
        top_text = "; ".join(
            f"{row['gene_symbol']}:{row['family_term_count']}"
            for row in family_drivers[:6]
        )
        summary_lines.append(
            f"| {family['family_id']} | {family['theme_label']} | {top_text} |"
        )
    summary_lines.extend(
        [
            "",
            "## Broadest recurrent drivers",
            "",
            "These are ranked by the number of the 106 significant terms containing each gene; this is a contribution count, not an independent gene-level significance claim.",
            "",
            "| Rank | Gene | Significant terms | Fine clusters | Theme families | GO:BP | Reactome |",
            "|---:|---|---:|---:|---:|---:|---:|",
        ]
    )
    for rank, row in enumerate(global_top, start=1):
        summary_lines.append(
            f"| {rank} | {row['gene_symbol']} | {row['significant_term_count']} | {row['redundancy_cluster_count']} | {row['theme_family_count']} | {row['GO_BP_term_count']} | {row['Reactome_term_count']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "The 106 terms are not 106 independent biological findings. Most are related GO descendants/ancestors or share the same 18 query genes. The reduced clusters, theme families, and driver counts summarize redundancy and repeated support; they do not add evidence beyond the frozen enrichment analysis. Cluster/family labels are transparent descriptive labels and should not be read as pathway activation or exposure causality.",
            "",
        ]
    )
    (OUT_DIR / "primary_term_reduction_summary.md").write_text("\n".join(summary_lines), encoding="utf-8")

    manifest = {
        "generated_utc": utc_now(),
        "script": str(Path(__file__).resolve()),
        "input_term_file": str(INPUT_TERM_FILE.resolve()),
        "input_term_file_sha256": sha256_file(INPUT_TERM_FILE),
        "query_gene_file": str(INPUT_QUERY_FILE.resolve()),
        "query_gene_file_sha256": sha256_file(INPUT_QUERY_FILE),
        "query_genes": query_genes,
        "background_gene_count": len(background),
        "background_rule": "archived GeneCards CRC export, RelevanceScore >= 10",
        "primary_significant_term_count": len(terms),
        "term_to_gene_evidence": "g:Profiler no_evidences=false; non-empty evidence list per query gene",
        "api_url": API_URL,
        "api_payload_file": str(payload_path.resolve()),
        "api_payload_sha256": sha256_file(payload_path),
        "api_raw_response_file": str(raw_path.resolve()),
        "api_raw_response_sha256": sha256_file(raw_path),
        "cluster_method": "average-linkage agglomeration",
        "cluster_threshold": args.cluster_threshold,
        "similarity_formula": "0.65 gene Jaccard + 0.20 normalized-term-token Jaccard + 0.10 parent-set Jaccard + 0.05 direct parent linkage",
        "fine_redundancy_cluster_count": len(clusters),
        "cluster_term_counts": {row["cluster_id"]: int(row["term_count"]) for row in cluster_rows},
        "theme_family_count": len(family_rows),
        "theme_family_term_counts": {row["family_id"]: int(row["term_count"]) for row in family_rows},
        "validation": {
            "all_106_terms_have_evidence": True,
            "all_api_overlap_counts_match_stored_intersection_size": True,
            "sensitivity_terms_mixed": False,
            "query_modified": False,
            "new_enrichment_test_run": False,
        },
        "output_sha256": {
            name: sha256_file(OUT_DIR / name)
            for name in [
                "primary_106_terms_with_genes.csv",
                "primary_theme_clusters.csv",
                "primary_cluster_term_membership.csv",
                "primary_18_gene_driver_contribution.csv",
                "primary_cluster_driver_contribution.csv",
                "primary_theme_families.csv",
                "primary_theme_family_term_membership.csv",
                "primary_theme_family_driver_contribution.csv",
                "primary_term_reduction_summary.md",
            ]
        },
    }
    (OUT_DIR / "primary_term_reduction_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print(f"Primary significant terms: {len(terms)}")
    print(f"Fine-grained redundancy clusters: {len(clusters)}")
    print(f"Broad theme families: {len(family_rows)}")
    print("Theme family term counts:", "; ".join(f"{row['family_id']}={row['term_count']}" for row in family_rows))
    print("Top drivers:", ", ".join(row["gene_symbol"] for row in global_top))
    print(f"Summary: {OUT_DIR / 'primary_term_reduction_summary.md'}")
    print(f"Manifest: {OUT_DIR / 'primary_term_reduction_manifest.json'}")


if __name__ == "__main__":
    main()
