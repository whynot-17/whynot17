#!/usr/bin/env python3
"""Step 8B: pathway redundancy reduction and module-level prioritization.

This stage is intentionally descriptive.  It does not re-test the 1,647
global-BH-significant pathway terms and it does not select a biological
mechanism.  It keeps every significant term, assigns it to a reproducible
module, and creates a compact representative table for interpretation.

The hierarchy is taken from the g:Profiler ``parents`` field captured in the
local raw response cache.  The resulting term hierarchy is written to a
small canonical CSV so later reruns do not require the large raw responses.
For GO/Reactome/KEGG, module edges use parent/ancestor structure; a guarded
ancestor-set similarity rule is used for terms without a direct significant
parent.  Generic terms covering >25% of the frozen background are not used
as graph bridges.  Cross-source consolidation is lexical and is explicitly
labelled as such; it is not treated as semantic proof.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from itertools import combinations
from pathlib import Path


DEFAULT_DIR = Path(__file__).resolve().parent
DEFAULT_RAW = DEFAULT_DIR / "raw_gprofiler"
ALL_RESULTS = DEFAULT_DIR / "t2d_step8_pathway_ora_all.csv"
HIERARCHY = DEFAULT_DIR / "t2d_step8_term_hierarchy.csv"

ROOT_TERMS = {
    "GO:0008150",  # biological_process
    "GO:0009987",  # cellular process
    "GO:0003674",  # molecular_function; defensive for mixed responses
    "GO:0005575",  # cellular_component; defensive for mixed responses
    "KEGG:00000",
    "REAC:0000000",
}

SOURCES = ("GO:BP", "REAC", "KEGG")
TIER_A = ("cluster_11", "cluster_5", "cluster_6", "cluster_8")
FDR_THRESHOLD = 0.05
MAX_BRIDGE_FRACTION = 0.25
ANCESTOR_JACCARD_THRESHOLD = 0.35
MIN_COMMON_ANCESTORS = 2
MAX_SELECTED_MODULES_PER_AXIS = 8

STOPWORDS = {
    "a", "an", "and", "as", "at", "by", "cell", "cellular", "of", "or",
    "process", "processes", "regulation", "response", "the", "to", "via",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def raw_paths(raw_dir: Path, cluster_id: str) -> list[Path]:
    legacy = raw_dir / f"gprofiler_{cluster_id}.json"
    if legacy.exists() and legacy.stat().st_size > 0:
        return [legacy]
    paths = []
    for source in SOURCES:
        path = raw_dir / f"gprofiler_{cluster_id}_{source.replace(':', '_')}.json"
        if path.exists() and path.stat().st_size > 0:
            paths.append(path)
    return paths


def build_hierarchy_from_raw(raw_dir: Path) -> list[dict[str, str]]:
    records: dict[tuple[str, str], dict[str, str]] = {}
    for cluster_id in TIER_A:
        paths = raw_paths(raw_dir, cluster_id)
        if not paths:
            raise FileNotFoundError(
                f"No g:Profiler cache found for {cluster_id}; provide {HIERARCHY} "
                "or rerun Step 8 pathway ORA first."
            )
        for path in paths:
            payload = json.loads(path.read_text(encoding="utf-8"))
            for item in payload.get("response", {}).get("result", []):
                source = str(item.get("source", ""))
                native = str(item.get("native", ""))
                if source not in SOURCES or not native:
                    continue
                key = (source, native)
                parents = item.get("parents", [])
                if not isinstance(parents, list):
                    parents = []
                parent_text = "|".join(sorted({str(parent) for parent in parents if parent}))
                records[key] = {
                    "source": source,
                    "native": native,
                    "term": str(item.get("name", "")),
                    "parents": parent_text,
                }
    rows = sorted(records.values(), key=lambda row: (row["source"], row["native"]))
    with HIERARCHY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "native", "term", "parents"])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def load_hierarchy(path: Path, raw_dir: Path) -> dict[tuple[str, str], set[str]]:
    if not path.exists() or path.stat().st_size == 0:
        build_hierarchy_from_raw(raw_dir)
    parents: dict[tuple[str, str], set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            source = row["source"]
            native = row["native"]
            for parent in row.get("parents", "").split("|"):
                if parent:
                    parents[(source, native)].add(parent)
    return parents


def ancestor_closure(
    source: str,
    native: str,
    parents: dict[tuple[str, str], set[str]],
    memo: dict[tuple[str, str], set[str]],
) -> set[str]:
    key = (source, native)
    if key in memo:
        return memo[key]
    seen: set[str] = set()
    stack = list(parents.get(key, set()))
    while stack:
        node = stack.pop()
        if node in seen or node in ROOT_TERMS:
            continue
        seen.add(node)
        stack.extend(parents.get((source, node), set()))
    memo[key] = seen
    return seen


class UnionFind:
    def __init__(self, items: list[tuple[str, str]]) -> None:
        self.parent = {item: item for item in items}

    def find(self, item: tuple[str, str]) -> tuple[str, str]:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            nxt = self.parent[item]
            self.parent[item] = root
            item = nxt
        return root

    def union(self, left: tuple[str, str], right: tuple[str, str]) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def tokens(text: str) -> set[str]:
    values = set(re.findall(r"[a-z0-9]+", text.lower())) - STOPWORDS
    normalized = set()
    for value in values:
        if value.endswith("ies") and len(value) > 5:
            value = value[:-3] + "y"
        elif value.endswith("ing") and len(value) > 6:
            value = value[:-3]
        elif value.endswith("ed") and len(value) > 5:
            value = value[:-2]
        elif value.endswith("es") and len(value) > 5:
            value = value[:-2]
        elif value.endswith("s") and len(value) > 4:
            value = value[:-1]
        normalized.add(value)
    return normalized


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def score_term(row: dict[str, str], background_size: int) -> float:
    q = max(float(row["global_bh_fdr"]), 1e-300)
    k = max(int(row["intersection_size"]), 1)
    term_size = max(int(row["term_size"]), 1)
    specificity = max(math.log(background_size / term_size), 0.0)
    return -math.log10(q) * math.log1p(k) * specificity


def choose_representative(members: list[dict[str, str]], background_size: int) -> dict[str, str] | None:
    eligible = [
        row for row in members
        if int(row["intersection_size"]) >= 3
        and int(row["term_size"]) / background_size <= MAX_BRIDGE_FRACTION
    ]
    if not eligible:
        return None
    candidates = eligible
    return max(
        candidates,
        key=lambda row: (
            score_term(row, background_size),
            int(row["intersection_size"]),
            -int(row["term_size"]),
            row["native"],
        ),
    )


def are_related(
    left: dict[str, str],
    right: dict[str, str],
    closures: dict[tuple[str, str], set[str]],
    term_tokens: dict[tuple[str, str], set[str]],
) -> bool:
    left_key = (left["source"], left["native"])
    right_key = (right["source"], right["native"])
    left_anc = closures[left_key]
    right_anc = closures[right_key]
    left_bridge = int(left["term_size"]) / int(left["effective_domain_size"]) > MAX_BRIDGE_FRACTION
    right_bridge = int(right["term_size"]) / int(right["effective_domain_size"]) > MAX_BRIDGE_FRACTION
    direct_or_ancestral = left["native"] in right_anc or right["native"] in left_anc
    if direct_or_ancestral and not (left_bridge or right_bridge):
        return True
    shared = (left_anc & right_anc) - ROOT_TERMS
    union = (left_anc | right_anc) - ROOT_TERMS
    same_source_similarity = len(shared) >= MIN_COMMON_ANCESTORS and (
        len(shared) / len(union) if union else 0.0
    ) >= ANCESTOR_JACCARD_THRESHOLD
    if same_source_similarity and not (left_bridge or right_bridge):
        return True
    # Cross-source consolidation is only used for the compact summary, and
    # is deliberately lexical rather than called semantic.
    return jaccard(term_tokens[left_key], term_tokens[right_key]) >= 0.60


def reduce_axis_source(
    rows: list[dict[str, str]],
    parents: dict[tuple[str, str], set[str]],
    background_size: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows = sorted(rows, key=lambda row: (row["native"], row["term"]))
    keys = [(row["source"], row["native"]) for row in rows]
    uf = UnionFind(keys)
    memo: dict[tuple[str, str], set[str]] = {}
    closures = {
        key: ancestor_closure(key[0], key[1], parents, memo)
        for key in keys
    }
    term_tokens = {key: tokens(row["term"]) for key, row in zip(keys, rows)}
    for left, right in combinations(rows, 2):
        left_key = (left["source"], left["native"])
        right_key = (right["source"], right["native"])
        left_bridge = int(left["term_size"]) / background_size > MAX_BRIDGE_FRACTION
        right_bridge = int(right["term_size"]) / background_size > MAX_BRIDGE_FRACTION
        if left_bridge or right_bridge:
            continue
        left_anc = closures[left_key]
        right_anc = closures[right_key]
        if left["native"] in right_anc or right["native"] in left_anc:
            uf.union(left_key, right_key)
            continue
        shared = (left_anc & right_anc) - ROOT_TERMS
        union = (left_anc | right_anc) - ROOT_TERMS
        similarity = len(shared) / len(union) if union else 0.0
        if len(shared) >= MIN_COMMON_ANCESTORS and similarity >= ANCESTOR_JACCARD_THRESHOLD:
            uf.union(left_key, right_key)
    components: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row["source"], row["native"])
        components[uf.find(key)].append(row)

    summaries: list[dict[str, str]] = []
    annotated: list[dict[str, str]] = []
    ordered_components = sorted(
        components.values(),
        key=lambda members: (
            -score_term(choose_representative(members, background_size), background_size)
            if choose_representative(members, background_size) else 0,
            choose_representative(members, background_size)["native"]
            if choose_representative(members, background_size) else "",
        ),
    )
    for module_number, members in enumerate(ordered_components, start=1):
        representative = choose_representative(members, background_size)
        module_id = f"{members[0]['cluster_id']}_{members[0]['source'].replace(':', '_')}_M{module_number:03d}"
        if representative is None:
            module_score = 0.0
            representative_eligible = "False"
        else:
            module_score = score_term(representative, background_size)
            representative_eligible = "True"
        for row in members:
            annotated.append({
                **row,
                "module_id": module_id,
                "module_n_terms": str(len(members)),
                "module_representative_native": representative["native"] if representative else "",
                "module_representative_term": representative["term"] if representative else "",
                "module_representative_score": f"{module_score:.8g}",
                "is_module_representative": str(representative is not None and row["native"] == representative["native"]),
                "module_representative_eligible": representative_eligible,
            })
        summaries.append({
            "cluster_id": members[0]["cluster_id"],
            "source": members[0]["source"],
            "module_id": module_id,
            "n_terms": str(len(members)),
            "representative_native": representative["native"] if representative else "",
            "representative_term": representative["term"] if representative else "",
            "representative_description": representative["description"] if representative else "",
            "representative_intersection_size": representative["intersection_size"] if representative else "",
            "representative_term_size": representative["term_size"] if representative else "",
            "representative_global_bh_fdr": representative["global_bh_fdr"] if representative else "",
            "module_score": f"{module_score:.8g}",
            "representative_eligible": representative_eligible,
            "reduction_method": "parent/ancestor graph; broad-term bridges excluded",
        })
    return annotated, summaries


def select_axis_modules(summaries: list[dict[str, str]]) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    for row in sorted(
        summaries,
        key=lambda item: (-float(item["module_score"]), item["source"], item["module_id"]),
    ):
        if row["representative_eligible"] != "True":
            continue
        if len([item for item in selected if item["cluster_id"] == row["cluster_id"]]) >= MAX_SELECTED_MODULES_PER_AXIS:
            continue
        row_tokens = tokens(row["representative_term"])
        too_close = False
        for prior in selected:
            if prior["cluster_id"] != row["cluster_id"]:
                continue
            if jaccard(row_tokens, tokens(prior["representative_term"])) >= 0.60:
                too_close = True
                break
            if prior["source"] == row["source"] and prior["representative_native"] == row["representative_native"]:
                too_close = True
                break
        if not too_close:
            selected.append(row)
    selected_by_id: dict[tuple[str, str], int] = {}
    axis_rank: dict[str, int] = defaultdict(int)
    for row in selected:
        axis_rank[row["cluster_id"]] += 1
        selected_by_id[(row["cluster_id"], row["module_id"])] = axis_rank[row["cluster_id"]]
    output = []
    for row in summaries:
        key = (row["cluster_id"], row["module_id"])
        output.append({**row, "selected_for_compact_summary": str(key in selected_by_id), "selected_rank": str(selected_by_id.get(key, ""))})
    return output


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_DIR)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_rows = [row for row in read_rows(args.output_dir / "t2d_step8_pathway_ora_all.csv") if float(row["global_bh_fdr"]) < FDR_THRESHOLD]
    if not all_rows:
        raise ValueError("No global-BH-significant pathway terms found")
    background_size = max(int(row["effective_domain_size"]) for row in all_rows)
    parents = load_hierarchy(args.output_dir / "t2d_step8_term_hierarchy.csv", args.raw_dir)

    annotated: list[dict[str, str]] = []
    summaries: list[dict[str, str]] = []
    for cluster_id in TIER_A:
        for source in SOURCES:
            subset = [row for row in all_rows if row["cluster_id"] == cluster_id and row["source"] == source]
            if not subset:
                continue
            reduced, module_rows = reduce_axis_source(subset, parents, background_size)
            annotated.extend(reduced)
            summaries.extend(module_rows)
    summaries = select_axis_modules(summaries)
    selected_keys = {
        (row["cluster_id"], row["module_id"])
        for row in summaries if row["selected_for_compact_summary"] == "True"
    }
    annotated = [
        {**row, "selected_for_compact_summary": str((row["cluster_id"], row["module_id"]) in selected_keys)}
        for row in annotated
    ]

    term_fields = list(all_rows[0].keys()) + [
        "module_id", "module_n_terms", "module_representative_native",
        "module_representative_term", "module_representative_score",
        "is_module_representative", "module_representative_eligible",
        "selected_for_compact_summary",
    ]
    term_fields = list(dict.fromkeys(term_fields))
    summary_fields = list(summaries[0].keys())
    write_csv(args.output_dir / "t2d_step8_pathway_modules.csv", annotated, term_fields)
    write_csv(args.output_dir / "t2d_step8_module_summary.csv", summaries, summary_fields)
    selected_rows = [row for row in summaries if row["selected_for_compact_summary"] == "True"]
    write_csv(args.output_dir / "t2d_step8_module_representatives.csv", selected_rows, summary_fields)

    counts = defaultdict(int)
    for row in annotated:
        counts[row["cluster_id"]] += 1
    report_lines = [
        "# Step 8B — T2D pathway redundancy reduction",
        "",
        "- Status: **complete_redundancy_reduction**",
        f"- Input: **{len(all_rows):,}** pathway terms with global BH-FDR < 0.05",
        f"- Frozen effective background: **{background_size:,} genes**",
        f"- Modules after parent/ancestor reduction: **{len(summaries):,}**",
        f"- Compact representatives retained: **{len(selected_rows):,}** (maximum 8 per axis)",
        "",
        "## Axis audit",
        "",
        "| Axis | Significant terms | Modules | Compact representatives |",
        "|---|---:|---:|---:|",
    ]
    for cluster_id in TIER_A:
        axis_terms = sum(1 for row in annotated if row["cluster_id"] == cluster_id)
        axis_modules = sum(1 for row in summaries if row["cluster_id"] == cluster_id)
        axis_selected = sum(1 for row in selected_rows if row["cluster_id"] == cluster_id)
        report_lines.append(f"| {cluster_id} | {axis_terms:,} | {axis_modules:,} | {axis_selected} |")
    report_lines.extend([
        "",
        "## Interpretation boundary",
        "",
        "All significant terms remain in `t2d_step8_pathway_modules.csv`; module representatives are a deterministic reduction layer, not a new statistical test.",
        "GO/Reactome/KEGG parent structure was used where available. Broad terms covering >25% of the effective background were not allowed to bridge modules. Cross-source similarity is lexical only and is not described as semantic evidence.",
        "A representative term is chosen using frozen global q value, intersection size, and term specificity. This stage does not infer pathway direction, activation, exposure causality, or mediation of T2D.",
        "",
        "## Compact representatives",
        "",
        "| Axis | Source | Representative | Overlap | Term size | Global q |",
        "|---|---|---|---:|---:|---:|",
    ])
    for row in sorted(selected_rows, key=lambda item: (item["cluster_id"], int(item["selected_rank"]))):
        report_lines.append(
            f"| {row['cluster_id']} | {row['source']} | {row['representative_term']} | {row['representative_intersection_size']} | {row['representative_term_size']} | {row['representative_global_bh_fdr']} |"
        )
    (args.output_dir / "STEP8B_T2D_REDUNDANCY_REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest = {
        "status": "complete_redundancy_reduction",
        "input": "t2d_step8_pathway_ora_all.csv filtered at global_bh_fdr < 0.05",
        "n_significant_terms": len(all_rows),
        "effective_background_size": background_size,
        "n_modules": len(summaries),
        "n_compact_representatives": len(selected_rows),
        "tier_a_axes": list(TIER_A),
        "parameters": {
            "max_bridge_fraction": MAX_BRIDGE_FRACTION,
            "ancestor_jaccard_threshold": ANCESTOR_JACCARD_THRESHOLD,
            "min_common_ancestors": MIN_COMMON_ANCESTORS,
            "max_selected_modules_per_axis": MAX_SELECTED_MODULES_PER_AXIS,
            "representative_min_intersection": 3,
        },
        "interpretation": "descriptive redundancy reduction; no new hypothesis test; no direction or causal inference",
        "raw_cache": "local-only raw_gprofiler cache; not committed",
    }
    (args.output_dir / "STEP8B_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
