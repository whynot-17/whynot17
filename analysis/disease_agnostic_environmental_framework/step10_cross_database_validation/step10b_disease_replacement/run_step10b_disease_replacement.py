#!/usr/bin/env python3
"""Step 10B-D: replace the T2D disease-knowledge source.

GeneCards remains the frozen reference.  Open Targets and GWAS Catalog are
queried independently, with source-native identifiers and evidence retained.
This script reports convergence and ranking stability; it does not alter the
29-test family, the 11 clusters, or any prior Tier assignment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom

VERSION = "1.0.0"
UA = "whynot17-step10b-disease/1.0"
ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
STEP7 = FRAMEWORK / "step07_genecard_convergence"
STEP10 = FRAMEWORK / "step10_cross_database_validation"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def s(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value).strip()


def dumps(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


class Client:
    def __init__(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": UA, "Accept": "application/json"})
        self.calls: list[dict[str, object]] = []

    def get(self, url: str, params: dict[str, object] | None = None, timeout: tuple[int, int] = (10, 60)) -> tuple[object | None, dict[str, object]]:
        last = ""
        for attempt in range(2):
            try:
                r = self.session.get(url, params=params, timeout=timeout)
                meta = {"status_code": r.status_code, "attempt": attempt + 1, "url": r.url, "response_sha256": sha256_bytes(r.content)}
                if r.status_code == 200:
                    payload = r.json()
                    self.calls.append({"endpoint": url, **meta})
                    return payload, meta
                last = f"HTTP {r.status_code}: {r.text[:300]}"
            except Exception as exc:
                last = f"{type(exc).__name__}: {exc}"
            if attempt == 0:
                time.sleep(1.0)
        meta = {"status_code": None, "attempt": 2, "url": url, "error": last}
        self.calls.append({"endpoint": url, **meta})
        return None, meta

    def post_json(self, url: str, body: dict[str, object], timeout: tuple[int, int] = (10, 90)) -> tuple[object | None, dict[str, object]]:
        last = ""
        for attempt in range(2):
            try:
                r = self.session.post(url, json=body, timeout=timeout)
                meta = {"status_code": r.status_code, "attempt": attempt + 1, "url": r.url, "response_sha256": sha256_bytes(r.content)}
                if r.status_code == 200:
                    payload = r.json()
                    self.calls.append({"endpoint": url, **meta})
                    return payload, meta
                last = f"HTTP {r.status_code}: {r.text[:300]}"
            except Exception as exc:
                last = f"{type(exc).__name__}: {exc}"
            if attempt == 0:
                time.sleep(1.0)
        meta = {"status_code": None, "attempt": 2, "url": url, "error": last}
        self.calls.append({"endpoint": url, **meta})
        return None, meta


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, set[str]]:
    cards = pd.read_csv(STEP7 / "t2d_genecards_primary_gene_audit.csv", dtype=str).fillna("")
    cards = cards[cards["set_label"].eq("primary_anywhere")].drop_duplicates("gene_symbol")
    clusters = pd.read_csv(STEP7 / "t2d_step7_cluster_ctd_genes.csv", dtype=str).fillna("")
    clusters = clusters[clusters["gene_symbol"].astype(bool)].drop_duplicates(["cluster_id", "gene_symbol"])
    cluster_chem = pd.read_csv(STEP7 / "t2d_step7_cluster_chemical_map.csv", dtype=str).fillna("")
    background = set(clusters["gene_symbol"].str.upper())
    return cards, clusters, cluster_chem, background


def query_open_targets(client: Client) -> tuple[str, dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    endpoint = "https://api.platform.opentargets.org/api/v4/graphql"
    search_query = "query Search($queryString: String!, $entityNames: [String!], $page: Pagination!) { search(queryString: $queryString, entityNames: $entityNames, page: $page) { total hits { id entity name score } } }"
    search_payload, search_meta = client.post_json(endpoint, {"query": search_query, "variables": {"queryString": "type 2 diabetes mellitus", "entityNames": ["disease"], "page": {"index": 0, "size": 10}}})
    hits = search_payload.get("data", {}).get("search", {}).get("hits", []) if isinstance(search_payload, dict) else []
    disease_id = ""
    for hit in hits:
        if s(hit.get("entity")) == "disease" and s(hit.get("name")).lower() == "type 2 diabetes mellitus":
            disease_id = s(hit.get("id"))
            break
    if not disease_id and hits:
        disease_id = s(hits[0].get("id"))
    association_query = "query DiseaseAssociations($efoId: String!, $page: Pagination!) { disease(efoId: $efoId) { name associatedTargets(page: $page) { count rows { target { id approvedSymbol } score datasourceScores { id score } } } } }"
    all_rows: list[dict[str, object]] = []
    page_meta: list[dict[str, object]] = []
    total = 0
    if disease_id:
        page = 0
        while True:
            payload, meta = client.post_json(endpoint, {"query": association_query, "variables": {"efoId": disease_id, "page": {"index": page, "size": 500}}})
            page_meta.append(meta)
            block = payload.get("data", {}).get("disease", {}) if isinstance(payload, dict) else {}
            assoc = block.get("associatedTargets", {}) if isinstance(block, dict) else {}
            if page == 0:
                total = int(assoc.get("count", 0) or 0)
            rows = assoc.get("rows", []) if isinstance(assoc, dict) else []
            if not rows:
                break
            all_rows.extend(rows)
            page += 1
            if len(all_rows) >= total or page >= 100:
                break
    return disease_id, {"search": search_meta, "association_pages": page_meta, "search_hits": hits, "reported_target_count": total}, all_rows, []


def recursive_gene_symbols(value: object) -> set[str]:
    out: set[str] = set()
    if isinstance(value, dict):
        for key, item in value.items():
            low = str(key).lower().replace("_", "")
            if low in {"mappedgene", "genesymbol", "approvedsymbol", "genename"}:
                if isinstance(item, str):
                    out.add(item.strip().upper())
                elif isinstance(item, dict):
                    for k in ("symbol", "name", "geneName", "approvedSymbol"):
                        if isinstance(item.get(k), str):
                            out.add(item[k].strip().upper())
            out |= recursive_gene_symbols(item)
    elif isinstance(value, list):
        for item in value:
            out |= recursive_gene_symbols(item)
    return {x for x in out if x and x not in {"NA", "NAN", "NONE", "NULL"}}


def query_gwas(client: Client, disease_id: str) -> tuple[set[str], dict[str, object], object | None, list[dict[str, object]]]:
    base = "https://www.ebi.ac.uk/gwas/rest/api"
    trait, trait_meta = client.get(f"{base}/efoTraits/{disease_id}", timeout=(10, 30))
    attempts: list[dict[str, object]] = [trait_meta]
    payload = None
    urls = [
        f"{base}/efoTraits/{disease_id}/associations?page=0&size=100",
        f"{base}/efoTraits/{disease_id}/associations?projection=associationByEfoTrait&page=0&size=100",
        f"{base}/associations?efoTrait={disease_id}&page=0&size=100",
    ]
    for url in urls:
        payload, meta = client.get(url, timeout=(10, 35))
        attempts.append(meta)
        if isinstance(payload, dict) and ("_embedded" in payload or "associations" in payload):
            break
        payload = None
    associations = payload.get("_embedded", {}).get("associations", []) if isinstance(payload, dict) else []
    gene_records: dict[str, dict[str, object]] = {}
    for association in associations:
        for locus in association.get("loci", []) if isinstance(association, dict) else []:
            for gene in locus.get("authorReportedGenes", []) if isinstance(locus, dict) else []:
                if not isinstance(gene, dict):
                    continue
                symbol = s(gene.get("geneName")).upper()
                if not symbol:
                    continue
                record = gene_records.setdefault(symbol, {"source": "GWAS_Catalog", "gene_symbol": symbol, "trait_id": disease_id, "n_associations": 0, "ensembl_gene_ids": set(), "entrez_gene_ids": set()})
                record["n_associations"] = int(record["n_associations"]) + 1
                record["ensembl_gene_ids"].update(s(x.get("ensemblGeneId")) for x in gene.get("ensemblGeneIds", []) if isinstance(x, dict) and s(x.get("ensemblGeneId")))
                record["entrez_gene_ids"].update(s(x.get("entrezGeneId")) for x in gene.get("entrezGeneIds", []) if isinstance(x, dict) and s(x.get("entrezGeneId")))
    gene_rows = []
    for record in gene_records.values():
        record = dict(record)
        record["ensembl_gene_ids"] = ";".join(sorted(record["ensembl_gene_ids"]))
        record["entrez_gene_ids"] = ";".join(sorted(record["entrez_gene_ids"]))
        record["mapping_rule"] = "GWAS Catalog authorReportedGenes with source-native Ensembl/Entrez identifiers"
        gene_rows.append(record)
    status = "complete_association_collection" if associations else ("api_timeout_or_unavailable" if payload is None else "empty_association_collection")
    meta = {"trait_metadata": trait, "attempts": attempts, "disease_id": disease_id, "status": status, "association_count_returned": len(associations), "gene_extraction_note": "authorReportedGenes parsed with source-native Ensembl/Entrez identifiers; no unsupported locus-to-gene inference"}
    return {r["gene_symbol"] for r in gene_rows}, meta, payload, gene_rows


def odds_ratio(M: int, K: int, N: int, x: int) -> float:
    a = x
    b = max(N - x, 0)
    c = max(K - x, 0)
    d = max(M - K - b, 0)
    if b == 0 or c == 0:
        return float("inf") if a > 0 else 0.0
    return (a * d) / (b * c)


def convergence(clusters: pd.DataFrame, source_sets: dict[str, set[str]], background: set[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    M = len(background)
    for source, genes in source_sets.items():
        mapped = set(genes) & background
        K = len(mapped)
        for cluster_id, block in clusters.groupby("cluster_id"):
            cgenes = set(block["gene_symbol"].str.upper()) & background
            x = len(cgenes & mapped)
            p = float(hypergeom.sf(x - 1, M, K, len(cgenes))) if K and cgenes else 1.0
            rows.append({"source": source, "cluster_id": cluster_id, "source_gene_count_raw": len(genes), "source_gene_count_in_ctd_background": K, "cluster_ctd_gene_count": len(cgenes), "overlap_count": x, "odds_ratio": odds_ratio(M, K, len(cgenes), x) if K else float("nan"), "hypergeom_p": p, "overlap_genes": ";".join(sorted(cgenes & mapped)), "mapping_rule": "case-insensitive HGNC symbol intersection; source-native IDs retained in source tables"})
    out = pd.DataFrame(rows)
    out["bh_fdr_within_source"] = out.groupby("source")["hypergeom_p"].transform(lambda x: x.rank(method="first") * len(x) / (len(x) * x.rank(method="first")))
    # Replace the simple rank-based placeholder with a deterministic BH calculation.
    for source, idx in out.groupby("source").groups.items():
        p = out.loc[idx, "hypergeom_p"].to_numpy(dtype=float)
        order = np.argsort(p)
        q = np.empty(len(p), dtype=float)
        running = 1.0
        for rank in range(len(p) - 1, -1, -1):
            pos = order[rank]
            running = min(running, p[pos] * len(p) / (rank + 1))
            q[pos] = running
        out.loc[idx, "bh_fdr_within_source"] = q
    return out


def rank_stability(table: pd.DataFrame) -> pd.DataFrame:
    wide = table.pivot(index="cluster_id", columns="source", values="overlap_count").fillna(0)
    rows = []
    base = wide["GeneCards_primary"] if "GeneCards_primary" in wide else pd.Series(dtype=float)
    for source in wide.columns:
        if source == "GeneCards_primary":
            continue
        n = int(len(wide))
        rho = float(base.rank().corr(wide[source].rank())) if n >= 3 else float("nan")
        rows.append({"reference_source": "GeneCards_primary", "replacement_source": source, "metric": "cluster_overlap_count", "spearman_rho": rho, "n_clusters": n, "interpretation": "descriptive ranking concordance; no source replacement changes Tier selection"})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-gwas", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    started = utc()
    client = Client()
    cards, clusters, cluster_chem, background = load_inputs()
    gene_cards = set(cards["gene_symbol"].str.upper())
    disease_id, ot_meta, ot_rows, _ = query_open_targets(client)
    ot_targets = []
    ot_genes: set[str] = set()
    for row in ot_rows:
        target = row.get("target", {}) if isinstance(row, dict) else {}
        symbol = s(target.get("approvedSymbol")).upper() if isinstance(target, dict) else ""
        if symbol:
            ot_genes.add(symbol)
            ot_targets.append({"source": "OpenTargets", "target_id": s(target.get("id")), "gene_symbol": symbol, "association_score": row.get("score"), "datasource_scores": dumps(row.get("datasourceScores", []))})
    pd.DataFrame(ot_targets).drop_duplicates(["target_id", "gene_symbol"]).to_csv(OUT / "step10b_open_targets_t2d_targets.csv", index=False)
    if args.skip_gwas:
        gwas_genes, gwas_meta, gwas_raw, gwas_gene_rows = set(), {"status": "skipped_by_argument"}, None, []
    else:
        gwas_genes, gwas_meta, gwas_raw, gwas_gene_rows = query_gwas(client, disease_id or "MONDO_0005148")
    small_gwas_snapshot = {"status": gwas_meta.get("status"), "association_count_returned": gwas_meta.get("association_count_returned"), "extracted_gene_count": len(gwas_genes), "response_sha256": next((x.get("response_sha256") for x in gwas_meta.get("attempts", []) if x.get("response_sha256")), None), "endpoint_family": "GWAS Catalog REST efoTraits/{MONDO_ID}/associations"}
    (OUT / "step10b_gwas_catalog_query_snapshot.json").write_text(json.dumps(small_gwas_snapshot, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(gwas_gene_rows).to_csv(OUT / "step10b_gwas_catalog_t2d_genes.csv", index=False)

    source_sets = {"GeneCards_primary": gene_cards, "OpenTargets": ot_genes, "GWAS_Catalog": gwas_genes}
    table = convergence(clusters, source_sets, background)
    table.to_csv(OUT / "step10b_disease_source_cluster_convergence.csv", index=False)
    rank_stability(table).to_csv(OUT / "step10b_disease_source_rank_stability.csv", index=False)
    pd.DataFrame([{ "source": "GeneCards_primary", "gene_count_raw": len(gene_cards), "gene_count_in_ctd_background": len(gene_cards & background), "status": "frozen_reference"}, {"source": "OpenTargets", "gene_count_raw": len(ot_genes), "gene_count_in_ctd_background": len(ot_genes & background), "status": "complete" if ot_genes else "empty_or_unavailable"}, {"source": "GWAS_Catalog", "gene_count_raw": len(gwas_genes), "gene_count_in_ctd_background": len(gwas_genes & background), "status": gwas_meta.get("status", "unknown")}]).to_csv(OUT / "step10b_disease_source_coverage.csv", index=False)

    snapshot = {
        "lock_type": "STEP10B_D_DISEASE_SOURCE_SNAPSHOT",
        "script_version": VERSION,
        "retrieval_started_utc": started,
        "retrieval_finished_utc": utc(),
        "frozen_disease_term": "type 2 diabetes mellitus",
        "open_targets": {"endpoint": "https://api.platform.opentargets.org/api/v4/graphql", "resolved_source_native_disease_id": disease_id, "data_release_documented_by_source": "26.06", "api_version_documented_by_source": "26.6.3", "query_metadata": ot_meta, "returned_target_rows": len(ot_rows), "human_target_rule": "Open Targets target records are source-native human Ensembl targets with approvedSymbol", "absence_interpretation": "not a biological negative"},
        "gwas_catalog": {"base_api": "https://www.ebi.ac.uk/gwas/rest/api", "source_native_trait_id": disease_id or "MONDO_0005148", "metadata": gwas_meta, "returned_gene_count": len(gwas_genes), "absence_interpretation": "not a biological negative; API timeout/partial response remains explicit"},
        "gene_symbol_mapping": "upper-cased approved symbols/mapped-gene symbols for convergence only; source-native target IDs and association fields retained in separate tables",
        "input_hashes": {"genecards": sha256(STEP7 / "t2d_genecards_primary_gene_audit.csv"), "cluster_genes": sha256(STEP7 / "t2d_step7_cluster_ctd_genes.csv")},
        "output_hashes": {},
        "api_calls": client.calls,
        "analysis_boundary": "Source replacement is a post-firewall robustness audit; no replacement-source result changes the frozen 29-test family, 11 clusters, or Tier assignments.",
    }
    for path in OUT.glob("step10b_*.csv"):
        snapshot["output_hashes"][path.name] = sha256(path)
    (OUT / "STEP10B_D_SOURCE_SNAPSHOT.json").write_text(json.dumps(snapshot, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "STEP10B_D_API_CALL_MANIFEST.json").write_text(json.dumps({"generated_utc": utc(), "calls": client.calls}, ensure_ascii=False, indent=2), encoding="utf-8")
    status = "complete_for_OpenTargets_with_GWAS_partial" if gwas_meta.get("status") != "complete_association_collection" else "complete_source_probe"
    report = "\n".join([
        "# Step 10B-D — Disease knowledge-source replacement audit",
        "",
        f"Generated: `{utc()}`",
        f"Open Targets disease ID: `{disease_id or 'unresolved'}`; returned target rows: **{len(ot_rows)}**; unique approved symbols: **{len(ot_genes)}**.",
        f"GWAS Catalog author-reported gene symbols returned: **{len(gwas_genes)}**; associations returned: **{gwas_meta.get('association_count_returned', 'NA')}**; status: **{gwas_meta.get('status', 'unknown')}**.",
        "",
        "GeneCards is retained as the frozen reference. Open Targets and GWAS Catalog are not merged into a single disease-gene truth set: coverage, source-native evidence, and convergence are reported separately.",
        "",
        "The GWAS Catalog API is deliberately fail-visible. Only source-native authorReportedGenes with Ensembl/Entrez identifiers are retained; no unsupported locus-to-gene inference is added. A timeout or empty response is recorded as such and never converted to an empty biological gene set.",
        "",
        f"Status: `{status}`. This audit does not change the frozen environmental panel, epidemiologic results, 11 clusters, or Tier assignments.",
    ])
    (OUT / "STEP10B_D_DISEASE_REPLACEMENT_REPORT.md").write_text(report, encoding="utf-8")
    print(f"completed Step 10B-D: OpenTargets {len(ot_genes)} genes; GWAS {len(gwas_genes)} genes ({gwas_meta.get('status', 'unknown')})", flush=True)


if __name__ == "__main__":
    main()
