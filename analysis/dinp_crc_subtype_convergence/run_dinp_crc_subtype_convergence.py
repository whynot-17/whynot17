#!/usr/bin/env python3
"""Localize the accessible DINP–CRC intersection across CRC anatomic anchors.

This is deliberately not a patient-level sidedness analysis. Open Targets
does not expose one exact right-sided versus left-sided CRC concept in the
search used here, so the analysis freezes ascending colon cancer as the right
anchor and sigmoid/rectosigmoid concepts as left-sided anchors. The output is
source-preserving and exploratory.
"""
from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
UPSTREAM = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs"
OUT = HERE / "outputs"
SRC = OUT / "source_records"
OT_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
UA = "whynot17-dinp-crc-subtype-convergence/1.0"

CONCEPTS = {
    "right_ascending": {
        "search_term": "ascending colon cancer",
        "expected_id": "MONDO_0002238",
        "expected_name": "ascending colon cancer",
        "side": "right_anchor",
    },
    "left_sigmoid_strict": {
        "search_term": "sigmoid colon cancer",
        "expected_id": "MONDO_0001464",
        "expected_name": "sigmoid colon cancer",
        "side": "left_anchor_strict",
    },
    "left_rectosigmoid_sensitivity": {
        "search_term": "rectosigmoid carcinoma",
        "expected_id": "MONDO_0002424",
        "expected_name": "rectosigmoid carcinoma",
        "side": "left_anchor_sensitivity",
    },
}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class HTTPClient:
    def __init__(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": UA, "Accept": "application/json"})
        self.calls: list[dict[str, Any]] = []

    def post(self, body: dict[str, Any]) -> tuple[Any, dict[str, Any]]:
        last_error = ""
        for attempt in range(1, 4):
            try:
                response = self.session.post(OT_ENDPOINT, json=body, timeout=(15, 120))
                raw = response.content
                meta: dict[str, Any] = {
                    "method": "POST",
                    "url": response.url,
                    "status_code": response.status_code,
                    "attempt": attempt,
                    "response_sha256": hashlib.sha256(raw).hexdigest(),
                    "content_length": len(raw),
                }
                if response.ok:
                    payload = response.json()
                    self.calls.append(meta)
                    return payload, meta
                last_error = f"HTTP {response.status_code}: {response.text[:300]}"
                meta["error"] = last_error
                if response.status_code not in {408, 425, 429, 500, 502, 503, 504}:
                    self.calls.append(meta)
                    return None, meta
            except Exception as exc:  # source availability is part of provenance
                last_error = f"{type(exc).__name__}: {exc}"
            if attempt < 3:
                time.sleep(1.5 * attempt)
        meta = {"method": "POST", "url": OT_ENDPOINT, "status_code": None, "attempt": 3, "error": last_error}
        self.calls.append(meta)
        return None, meta


def clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def dump_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def search_concept(client: HTTPClient, term: str) -> tuple[dict[str, Any], dict[str, Any]]:
    query = """query Search($queryString: String!, $entityNames: [String!], $page: Pagination!) {
      search(queryString: $queryString, entityNames: $entityNames, page: $page) {
        total hits { id entity name score }
      }
    }"""
    body = {"query": query, "variables": {"queryString": term, "entityNames": ["disease"], "page": {"index": 0, "size": 50}}}
    payload, meta = client.post(body)
    hits = payload.get("data", {}).get("search", {}).get("hits", []) if isinstance(payload, dict) else []
    exact = next((h for h in hits if clean(h.get("name")).lower() == term.lower()), None)
    chosen = exact or (hits[0] if hits else {})
    return {"request": body, "response": payload, "hits": hits, "chosen": chosen}, meta


def fetch_targets(client: HTTPClient, disease_id: str) -> tuple[list[dict[str, Any]], int, list[dict[str, Any]]]:
    query = """query DiseaseAssociations($efoId: String!, $page: Pagination!) {
      disease(efoId: $efoId) {
        name
        associatedTargets(page: $page) {
          count
          rows { target { id approvedSymbol } score datasourceScores { id score } }
        }
      }
    }"""
    rows: list[dict[str, Any]] = []
    calls: list[dict[str, Any]] = []
    total = 0
    for page in range(100):
        body = {"query": query, "variables": {"efoId": disease_id, "page": {"index": page, "size": 500}}}
        payload, meta = client.post(body)
        calls.append(meta)
        block = payload.get("data", {}).get("disease", {}) if isinstance(payload, dict) else {}
        assoc = block.get("associatedTargets", {}) if isinstance(block, dict) else {}
        if page == 0:
            total = int(assoc.get("count", 0) or 0)
        batch = assoc.get("rows", []) if isinstance(assoc, dict) else []
        if not batch:
            break
        rows.extend(batch)
        if len(rows) >= total:
            break
    return rows, total, calls


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    n = len(p_values)
    order = sorted(range(n), key=lambda i: p_values[i])
    adjusted = [1.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * n / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)
    intersection = pd.read_csv(UPSTREAM / "dinp_crc_intersection.csv", dtype=str).fillna("")
    exposure = pd.read_csv(UPSTREAM / "dinp_exposure_gene_matrix.csv", dtype=str).fillna("")
    crc = pd.read_csv(UPSTREAM / "crc_gene_matrix.csv", dtype=str).fillna("")
    intersection["gene_symbol"] = intersection["gene_symbol"].str.upper().str.strip()
    exposure["gene_symbol"] = exposure["gene_symbol"].str.upper().str.strip()
    crc["gene_symbol"] = crc["gene_symbol"].str.upper().str.strip()
    genes_81 = set(intersection["gene_symbol"])
    exposure_genes = set(exposure["gene_symbol"])
    crc_genes = set(crc["gene_symbol"])
    client = HTTPClient()
    resolved: dict[str, dict[str, Any]] = {}
    target_sets: dict[str, set[str]] = {}
    target_scores: dict[str, dict[str, float]] = {}
    query_archive: dict[str, Any] = {"generated_at": utc(), "concepts": {}}
    for label, config in CONCEPTS.items():
        search, search_meta = search_concept(client, config["search_term"])
        chosen = search.get("chosen", {})
        concept_id = clean(chosen.get("id"))
        exact_ok = concept_id == config["expected_id"] and clean(chosen.get("name")).lower() == config["expected_name"]
        rows: list[dict[str, Any]] = []
        total = 0
        assoc_calls: list[dict[str, Any]] = []
        if exact_ok:
            rows, total, assoc_calls = fetch_targets(client, concept_id)
        genes: set[str] = set()
        scores: dict[str, float] = {}
        flat_rows: list[dict[str, Any]] = []
        for row in rows:
            target = row.get("target") or {}
            symbol = clean(target.get("approvedSymbol")).upper()
            if not symbol:
                continue
            genes.add(symbol)
            score = row.get("score")
            if isinstance(score, (int, float)):
                scores[symbol] = float(score)
            flat_rows.append({
                "concept_label": label,
                "disease_id": concept_id,
                "disease_name": clean(chosen.get("name")),
                "target_id": clean(target.get("id")),
                "gene_symbol": symbol,
                "overall_score": score,
                "datasource_scores_json": json.dumps(row.get("datasourceScores", []), ensure_ascii=False, sort_keys=True),
            })
        target_sets[label] = genes
        target_scores[label] = scores
        resolved[label] = {
            **config,
            "resolved_id": concept_id,
            "resolved_name": clean(chosen.get("name")),
            "exact_expected_concept": exact_ok,
            "reported_target_count": total,
            "returned_target_rows": len(flat_rows),
            "unique_approved_symbols": len(genes),
            "search_request": search.get("request"),
            "search_request_meta": search_meta,
            "association_request_meta": assoc_calls,
        }
        query_archive["concepts"][label] = {"config": config, "search": search, "search_meta": search_meta, "association_request_meta": assoc_calls, "rows": rows}
        pd.DataFrame(flat_rows).to_csv(SRC / f"opentargets_{label}_targets.csv", index=False)
    dump_json(SRC / "opentargets_subtype_queries.json", query_archive)

    # Preserve the 81-gene input and add source-specific anatomical support.
    out = intersection[intersection["gene_symbol"].isin(genes_81)].copy()
    for label, genes in target_sets.items():
        out[f"OpenTargets_{label}"] = out["gene_symbol"].isin(genes).astype(int)
        out[f"OpenTargets_{label}_score"] = out["gene_symbol"].map(target_scores[label]).fillna("")
    out["left_expanded_support"] = ((out["OpenTargets_left_sigmoid_strict"] == 1) | (out["OpenTargets_left_rectosigmoid_sensitivity"] == 1)).astype(int)
    def pattern(row: pd.Series) -> str:
        right = int(row["OpenTargets_right_ascending"])
        left = int(row["left_expanded_support"])
        return "both" if right and left else "right_only" if right else "left_only" if left else "neither"
    out["anatomic_pattern"] = out.apply(pattern, axis=1)
    out.to_csv(OUT / "dinp_crc_subtype_gene_support.csv", index=False)

    # Exploratory enrichment uses the frozen general CRC union as the universe.
    comparisons = [
        ("right_ascending", target_sets["right_ascending"]),
        ("left_sigmoid_strict", target_sets["left_sigmoid_strict"]),
        ("left_expanded_sigmoid_or_rectosigmoid", target_sets["left_sigmoid_strict"] | target_sets["left_rectosigmoid_sensitivity"]),
    ]
    enrichment_rows: list[dict[str, Any]] = []
    for label, targets in comparisons:
        subtype = targets & crc_genes
        a = len(genes_81 & subtype)
        b = len(genes_81 - subtype)
        c = len(subtype - genes_81)
        d = len(crc_genes - genes_81 - subtype)
        odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        enrichment_rows.append({
            "comparison": label,
            "background_crc_union_n": len(crc_genes),
            "dinp_crc_intersection_n": len(genes_81),
            "subtype_target_n_in_background": len(subtype),
            "intersection_genes_in_subtype": a,
            "intersection_genes_not_in_subtype": b,
            "subtype_genes_not_in_intersection": c,
            "background_genes_neither": d,
            "fisher_greater_odds_ratio": odds,
            "fisher_greater_p": p_value,
        })
    enrichment = pd.DataFrame(enrichment_rows)
    enrichment["BH_FDR_3_subtype_comparisons"] = benjamini_hochberg(enrichment["fisher_greater_p"].tolist())
    enrichment.to_csv(OUT / "dinp_crc_subtype_enrichment.csv", index=False)

    pattern_counts = out["anatomic_pattern"].value_counts().to_dict()
    summary = {
        "generated_at": utc(),
        "input": {
            "dinp_crc_intersection_n": len(genes_81),
            "dinp_exposure_union_n": len(exposure_genes),
            "general_crc_union_n": len(crc_genes),
            "upstream_intersection_sha256": sha256_file(UPSTREAM / "dinp_crc_intersection.csv"),
        },
        "concepts": resolved,
        "pattern_counts": pattern_counts,
        "limitations": [
            "Open Targets has no single exact right-sided-versus-left-sided CRC concept in this run.",
            "Ascending colon cancer is used as a right-sided anatomical anchor; sigmoid colon cancer is the strict left anchor; rectosigmoid carcinoma is an expanded left sensitivity anchor.",
            "This is not patient-level sidedness, expression, survival, or epidemiologic effect modification.",
            "The Fisher tests are exploratory knowledge-base localization tests against the frozen general CRC union; they do not establish DINP causality or subtype-specific exposure effects.",
        ],
    }
    dump_json(OUT / "subtype_manifest.json", {"analysis": summary, "http_call_log": client.calls})
    lines = [
        "# DINP–CRC anatomical subtype convergence summary",
        "",
        f"Generated: `{summary['generated_at']}`",
        "",
        "## Frozen interpretation",
        "",
        "The 81-gene accessible DINP–CRC intersection was localized using source-native Open Targets anatomical CRC concepts. This is not a patient-level right-versus-left CRC analysis.",
        "",
        "## Input counts",
        "",
        f"- DINP exposure union: **{len(exposure_genes)} genes**",
        f"- General CRC union: **{len(crc_genes)} genes**",
        f"- Frozen DINP–CRC intersection: **{len(genes_81)} genes**",
        "",
        "## Open Targets concepts",
        "",
        "| Label | Concept | ID | Unique targets | Exact expected concept |",
        "|---|---|---|---:|---|",
    ]
    for label, info in resolved.items():
        lines.append(f"| `{label}` | {info['resolved_name']} | `{info['resolved_id']}` | {info['unique_approved_symbols']} | {info['exact_expected_concept']} |")
    lines += [
        "",
        "## 81-gene localization",
        "",
        "| Pattern | Genes |",
        "|---|---:|",
    ]
    for label in ["both", "right_only", "left_only", "neither"]:
        lines.append(f"| `{label}` | {pattern_counts.get(label, 0)} |")
    lines += [
        "",
        "## Exploratory enrichment",
        "",
        "| Comparison | Target genes in background | 81-gene overlap | Fisher OR | P | BH-FDR |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in enrichment.iterrows():
        lines.append(f"| `{row['comparison']}` | {int(row['subtype_target_n_in_background'])} | {int(row['intersection_genes_in_subtype'])} | {row['fisher_greater_odds_ratio']:.4g} | {row['fisher_greater_p']:.4g} | {row['BH_FDR_3_subtype_comparisons']:.4g} |")
    lines += [
        "",
        "## Boundaries",
        "",
        "- Anatomical anchor concepts are not equivalent to a clinical right/left-sided CRC phenotype label.",
        "- Source support is descriptive; no heterogeneous source evidence was collapsed into a biological truth score.",
        "- The upstream 81-gene intersection was not modified; these outputs are a subtype-localization follow-up only.",
        "",
        "## Files",
        "",
        "- `dinp_crc_subtype_gene_support.csv`: 81-gene source-preserving anatomical support table.",
        "- `dinp_crc_subtype_enrichment.csv`: exploratory enrichment calculations.",
        "- `subtype_manifest.json` and `source_records/`: concept resolution, target rows, requests, and hashes.",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"pattern_counts": pattern_counts, "enrichment": enrichment_rows}, indent=2, default=str))


if __name__ == "__main__":
    main()
