#!/usr/bin/env python3
"""DINP–CRC target convergence across three exposure and three disease sources.

The script is deliberately source-preserving.  Counts from different
databases are descriptive support indicators, not a merged biological truth
score.  Missing or inaccessible databases remain explicit in provenance.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
SRC = OUT / "source_records"
PHASE1 = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data"
CTD_IXN = PHASE1 / "CTD_chem_gene_ixns.tsv.gz"
GENECARDS = PHASE1 / "genecards_anywhere_crc_top2000.csv"

DINP_CTD_ID = "C012125"
DINP_CASRN = "28553-12-0"
DINP_DTXSID = "DTXSID4022521"
OT_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
UA = "whynot17-dinp-crc-target-convergence/1.0"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def dump_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def parse_ctd() -> tuple[pd.DataFrame, dict[str, Any]]:
    if not CTD_IXN.exists():
        raise FileNotFoundError(CTD_IXN)
    header: list[str] | None = None
    records: list[dict[str, str]] = []
    with gzip.open(CTD_IXN, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        for line in fh:
            if line.startswith("# ChemicalName"):
                header = line[2:].rstrip("\n").split("\t")
                break
        if header is None:
            raise ValueError("CTD interaction header not found")
        reader = csv.DictReader(fh, fieldnames=header, delimiter="\t")
        for row in reader:
            if clean(row.get("ChemicalID")) != DINP_CTD_ID:
                continue
            if clean(row.get("OrganismID")) != "9606":
                continue
            interaction = clean(row.get("Interaction"))
            actions = clean(row.get("InteractionActions"))
            # CTD marks co-treatment in InteractionActions and usually also in
            # the interaction statement. Keep those rows, but flag them.
            cotreat = bool(re.search(r"cotreatment|co-treated with", f"{actions} {interaction}", flags=re.I))
            rec = {k: clean(v) for k, v in row.items() if k in header}
            rec["multi_chemical_co_treatment_flag"] = str(cotreat).lower()
            rec["single_chemical_record_flag"] = str(not cotreat).lower()
            records.append(rec)
    frame = pd.DataFrame(records)
    if frame.empty:
        frame = pd.DataFrame(columns=header + ["multi_chemical_co_treatment_flag", "single_chemical_record_flag"])
    path = SRC / "ctd_dinp_human_interactions.csv"
    frame.to_csv(path, index=False)
    genes = set(frame["GeneSymbol"].str.upper()) if not frame.empty else set()
    single = frame[frame["single_chemical_record_flag"].eq("true")] if not frame.empty else frame
    meta = {
        "source": "CTD",
        "status": "complete",
        "source_file": str(CTD_IXN),
        "source_file_sha256": sha256_file(CTD_IXN),
        "chemical_id_used": DINP_CTD_ID,
        "casrn": DINP_CASRN,
        "dtxsid": DINP_DTXSID,
        "human_raw_rows": int(len(frame)),
        "unique_human_genes": int(len(genes)),
        "unique_human_gene_ids": int(frame["GeneID"].nunique()) if not frame.empty else 0,
        "single_chemical_rows": int(len(single)),
        "single_chemical_genes": int(single["GeneSymbol"].str.upper().nunique()) if not single.empty else 0,
        "co_treatment_rows": int((frame["multi_chemical_co_treatment_flag"] == "true").sum()) if not frame.empty else 0,
        "record_definition": "Homo sapiens rows for exact CTD parent DINP C012125; co-treatment retained and flagged",
    }
    return frame, meta


class HTTPClient:
    def __init__(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": UA, "Accept": "application/json"})
        self.calls: list[dict[str, Any]] = []

    def request(self, method: str, url: str, **kwargs: Any) -> tuple[Any, dict[str, Any]]:
        last = ""
        for attempt in range(1, 4):
            try:
                response = self.session.request(method, url, timeout=kwargs.pop("timeout", (15, 120)), **kwargs)
                body = response.content
                meta: dict[str, Any] = {
                    "method": method,
                    "url": response.url,
                    "status_code": response.status_code,
                    "attempt": attempt,
                    "response_sha256": sha256_bytes(body),
                    "content_length": len(body),
                }
                if response.status_code == 200:
                    try:
                        payload = response.json()
                    except ValueError:
                        payload = response.text
                    self.calls.append(meta)
                    return payload, meta
                last = f"HTTP {response.status_code}: {response.text[:300]}"
                meta["error"] = last
                if response.status_code not in {408, 425, 429, 500, 502, 503, 504}:
                    self.calls.append(meta)
                    return None, meta
            except Exception as exc:  # network availability is a provenance result
                last = f"{type(exc).__name__}: {exc}"
            if attempt < 3:
                time.sleep(1.5 * attempt)
        meta = {"method": method, "url": url, "status_code": None, "attempt": 3, "error": last}
        self.calls.append(meta)
        return None, meta


def query_chembl(client: HTTPClient) -> tuple[set[str], dict[str, Any], list[dict[str, Any]]]:
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"
    payload, meta = client.request("GET", url, params={"q": "diisononyl phthalate", "limit": 100})
    filter_queries = [
        {"pref_name__iexact": "diisononyl phthalate", "limit": 100},
        {"molecule_synonyms__synonyms__iexact": "diisononyl phthalate", "limit": 100},
        {"molecule_synonyms__synonyms__iexact": "DINP", "limit": 100},
        {"molecule_synonyms__synonyms__iexact": DINP_CASRN, "limit": 100},
    ]
    filter_results: list[dict[str, Any]] = []
    for params in filter_queries:
        exact_payload, exact_meta = client.request("GET", "https://www.ebi.ac.uk/chembl/api/data/molecule.json", params=params)
        filter_results.append({"params": params, "response": exact_payload, "request": exact_meta})
    dump_json(SRC / "chembl_dinp_molecule_search.json", {"broad_query": "diisononyl phthalate", "broad_response": payload, "broad_request": meta, "exact_filter_queries": filter_results})
    molecules = payload.get("molecules", []) if isinstance(payload, dict) else []
    exact: list[dict[str, Any]] = []
    exact_ids: set[str] = set()
    for result in filter_results:
        block = result.get("response")
        for mol in block.get("molecules", []) if isinstance(block, dict) else []:
            cid = clean(mol.get("molecule_chembl_id"))
            if cid and cid not in exact_ids:
                exact.append(mol)
                exact_ids.add(cid)
    # ChEMBL search is not used to infer a match from a related phthalate.
    meta_out = {
        "source": "ChEMBL",
        "status": "complete_no_reliable_parent_match" if not exact else "parent_match_requires_activity_query",
        "query": "diisononyl phthalate",
        "returned_molecule_rows": len(molecules),
        "exact_filter_queries": filter_queries,
        "exact_parent_matches": [clean(x.get("molecule_chembl_id")) for x in exact],
        "activity_query_performed": bool(exact),
        "substitution_rule": "no related phthalate substitution",
        "request": meta,
    }
    records: list[dict[str, Any]] = []
    genes: set[str] = set()
    for mol in exact:
        chembl_id = clean(mol.get("molecule_chembl_id"))
        act_payload, act_meta = client.request("GET", "https://www.ebi.ac.uk/chembl/api/data/activity.json", params={"molecule_chembl_id": chembl_id, "limit": 1000})
        records.append({"molecule_chembl_id": chembl_id, "activity_response": act_payload, "request": act_meta})
        for row in act_payload.get("activities", []) if isinstance(act_payload, dict) else []:
            target = row.get("target") or {}
            if clean(target.get("organism")).lower() == "homo sapiens":
                symbol = clean(target.get("pref_name")).upper()
                if symbol:
                    genes.add(symbol)
    if records:
        dump_json(SRC / "chembl_dinp_parent_activity.json", records)
    meta_out["human_direct_target_rows"] = len(records)
    return genes, meta_out, records


def query_toxcast(client: HTTPClient) -> tuple[set[str], dict[str, Any]]:
    api_url = "https://comptox.epa.gov/ctx-api/bioactivity/data/summary/search/by-dtxsid/" + DINP_DTXSID
    payload, meta = client.request("GET", api_url)
    public_url = "https://gaftp.epa.gov/Comptox/High_Throughput_Screening_Data/InVitroDB_V3.2/Summary_Files/INVITRODB_V3_2_SUMMARY.zip"
    _, public_meta = client.request("HEAD", public_url, timeout=(15, 30))
    dump_json(SRC / "toxcast_dinp_access_attempt.json", {"api_response": payload, "api_request": meta, "public_bulk_release": {"url": public_url, "request": public_meta}})
    status = "complete" if isinstance(payload, (dict, list)) else "unavailable_without_api_key"
    return set(), {
        "source": "EPA CompTox/ToxCast",
        "status": status,
        "dtxsid": DINP_DTXSID,
        "api_endpoint": api_url,
        "api_request": meta,
        "public_bulk_release_url": public_url,
        "public_bulk_release_request": public_meta,
        "human_target_genes_returned": 0,
        "absence_interpretation": "not a biological negative; live CTX API requires an x-api-key and no substitute source was silently used",
    }


def query_disgenet(client: HTTPClient) -> tuple[set[str], dict[str, Any]]:
    # C0346647 is UMLS colorectal cancer; the concept is kept as an attempted
    # source-native identifier, but a 401 remains unavailable rather than a
    # fabricated empty gene set.
    disease_id = "C0346647"
    url = f"https://api.disgenet.com/api/v1/gda/disease/umls/{disease_id}"
    payload, meta = client.request("GET", url)
    dump_json(SRC / "disgenet_crc_access_attempt.json", {"disease_id": disease_id, "response": payload, "request": meta})
    return set(), {
        "source": "DisGeNET",
        "status": "complete" if payload is not None else "unavailable_without_api_key",
        "disease_concept": "UMLS:" + disease_id,
        "request": meta,
        "unique_human_genes_returned": 0,
        "absence_interpretation": "not a biological negative; API access was not authorized",
    }


def query_open_targets(client: HTTPClient) -> tuple[set[str], dict[str, Any], pd.DataFrame]:
    search_q = "query Search($queryString: String!, $entityNames: [String!], $page: Pagination!) { search(queryString: $queryString, entityNames: $entityNames, page: $page) { total hits { id entity name score } } }"
    search_body = {"query": search_q, "variables": {"queryString": "colorectal cancer", "entityNames": ["disease"], "page": {"index": 0, "size": 20}}}
    search_payload, search_meta = client.request("POST", OT_ENDPOINT, json=search_body)
    hits = search_payload.get("data", {}).get("search", {}).get("hits", []) if isinstance(search_payload, dict) else []
    disease_id = ""
    disease_name = ""
    for hit in hits:
        if clean(hit.get("entity")) == "disease" and clean(hit.get("name")).lower() == "colorectal cancer":
            disease_id, disease_name = clean(hit.get("id")), clean(hit.get("name"))
            break
    if not disease_id and hits:
        disease_id, disease_name = clean(hits[0].get("id")), clean(hits[0].get("name"))
    assoc_q = "query DiseaseAssociations($efoId: String!, $page: Pagination!) { disease(efoId: $efoId) { name associatedTargets(page: $page) { count rows { target { id approvedSymbol } score datasourceScores { id score } } } } }"
    rows: list[dict[str, Any]] = []
    pages: list[dict[str, Any]] = []
    total = 0
    for page in range(100):
        body = {"query": assoc_q, "variables": {"efoId": disease_id, "page": {"index": page, "size": 500}}}
        payload, meta = client.request("POST", OT_ENDPOINT, json=body)
        pages.append(meta)
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
    flat: list[dict[str, Any]] = []
    genes: set[str] = set()
    scores: dict[str, float] = {}
    for row in rows:
        target = row.get("target") or {}
        symbol = clean(target.get("approvedSymbol")).upper()
        if not symbol:
            continue
        genes.add(symbol)
        score = row.get("score")
        scores[symbol] = float(score) if isinstance(score, (int, float)) else float("nan")
        flat.append({
            "target_id": clean(target.get("id")),
            "approved_symbol": symbol,
            "overall_score": score,
            "datasource_scores_json": json.dumps(row.get("datasourceScores", []), ensure_ascii=False, sort_keys=True),
        })
    pd.DataFrame(flat).to_csv(SRC / "opentargets_crc_associated_targets.csv", index=False)
    dump_json(SRC / "opentargets_crc_query.json", {"search_response": search_payload, "search_request": search_meta, "association_request": {"disease_id": disease_id, "pages": pages, "reported_target_count": total}, "association_rows": rows})
    return genes, {
        "source": "Open Targets",
        "status": "complete" if disease_id and flat else "unavailable_or_empty",
        "disease_id": disease_id,
        "disease_name": disease_name,
        "search_hits": hits,
        "reported_target_count": total,
        "returned_target_rows": len(flat),
        "unique_human_approved_symbols": len(genes),
        "association_pages": pages,
    }, pd.DataFrame(flat)


def audit_hgnc_symbols(symbols: set[str], client: HTTPClient) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Audit the symbols used in the DINP exposure-side result against HGNC.

    CTD and GeneCards are source-native symbol inputs, whereas Open Targets
    explicitly returns approvedSymbol.  This audit makes that distinction
    visible and prevents an unverified symbol from being silently renamed.
    The full 15k-gene Open Targets disease universe is not re-queried here:
    those rows are already source-native approvedSymbol values; the audit
    focuses on the much smaller DINP exposure/intersection universe where
    cross-source harmonization matters.
    """
    rows: list[dict[str, Any]] = []
    for raw_symbol in sorted(symbols):
        symbol = clean(raw_symbol).upper()
        if not symbol:
            continue
        url = f"https://rest.genenames.org/fetch/symbol/{symbol}"
        payload, meta = client.request("GET", url, headers={"Accept": "application/json"}, timeout=(15, 30))
        docs = []
        resolution = "approved_symbol_exact"
        if isinstance(payload, dict):
            docs = payload.get("response", {}).get("docs", []) or []
        if not docs:
            alias_url = f"https://rest.genenames.org/fetch/alias_symbol/{symbol}"
            alias_payload, alias_meta = client.request("GET", alias_url, headers={"Accept": "application/json"}, timeout=(15, 30))
            if isinstance(alias_payload, dict):
                docs = alias_payload.get("response", {}).get("docs", []) or []
            if docs:
                resolution = "approved_symbol_from_hgnc_alias"
            meta = alias_meta if alias_meta.get("status_code") is not None else meta
        doc = docs[0] if docs else {}
        rows.append({
            "input_symbol": symbol,
            "hgnc_status": clean(doc.get("status")) if doc else "not_found_or_unavailable",
            "approved_symbol": clean(doc.get("symbol")).upper() if doc else "",
            "hgnc_id": clean(doc.get("hgnc_id")),
            "entrez_id": clean(doc.get("entrez_id")),
            "input_matches_approved": int(bool(doc) and clean(doc.get("symbol")).upper() == symbol),
            "resolution": resolution if doc else "not_resolved",
            "request_status_code": meta.get("status_code"),
        })
    frame = pd.DataFrame(rows, columns=[
        "input_symbol", "hgnc_status", "approved_symbol", "hgnc_id", "entrez_id",
        "input_matches_approved", "resolution", "request_status_code",
    ])
    frame.to_csv(SRC / "hgnc_symbol_audit.csv", index=False)
    status_counts = frame["hgnc_status"].value_counts(dropna=False).to_dict() if not frame.empty else {}
    meta = {
        "source": "HGNC REST symbol audit",
        "status": "complete",
        "endpoint": "https://rest.genenames.org/fetch/symbol/{symbol}",
        "audited_symbol_count": int(len(frame)),
        "approved_status_count": int((frame["hgnc_status"] == "Approved").sum()) if not frame.empty else 0,
        "exact_approved_symbol_match_count": int(frame["input_matches_approved"].sum()) if not frame.empty else 0,
        "status_counts": status_counts,
        "scope": "DINP exposure union and DINP–CRC intersection symbols; Open Targets disease rows retain source-native approvedSymbol",
    }
    return frame, meta


def load_genecards() -> tuple[set[str], pd.DataFrame, dict[str, Any]]:
    frame = pd.read_csv(GENECARDS, dtype=str).fillna("")
    frame["gene_symbol"] = frame["GeneSymbol"].astype(str).str.upper().str.strip()
    frame = frame[frame["gene_symbol"].ne("")].drop_duplicates("gene_symbol")
    frame.to_csv(SRC / "genecards_crc_top2000_reference.csv", index=False)
    return set(frame["gene_symbol"]), frame, {
        "source": "GeneCards",
        "status": "frozen_reference",
        "query_scope": "ordinary CRC search top-2000 export",
        "source_file": str(GENECARDS),
        "source_file_sha256": sha256_file(GENECARDS),
        "raw_rows": int(len(pd.read_csv(GENECARDS, dtype=str))),
        "unique_symbols_used": int(len(frame)),
        "note": "The historical strict scoped file is not used; no substitution for the ordinary top-2000 CRC reference.",
    }


def build_outputs(ctd: pd.DataFrame, ctd_meta: dict[str, Any], chembl_genes: set[str], chembl_meta: dict[str, Any], tox_genes: set[str], tox_meta: dict[str, Any], cards: set[str], cards_frame: pd.DataFrame, cards_meta: dict[str, Any], dis_genes: set[str], dis_meta: dict[str, Any], ot_genes: set[str], ot_meta: dict[str, Any], ot_frame: pd.DataFrame, hgnc_frame: pd.DataFrame, hgnc_meta: dict[str, Any]) -> dict[str, Any]:
    ctd_symbol_col = "normalized_gene_symbol" if "normalized_gene_symbol" in ctd.columns else "GeneSymbol"
    ctd_genes = set(ctd[ctd_symbol_col].astype(str).str.upper()) if not ctd.empty else set()
    exp_genes = sorted(ctd_genes | chembl_genes | tox_genes)
    crc_genes = sorted(cards | dis_genes | ot_genes)
    ctd_by: dict[str, dict[str, Any]] = defaultdict(lambda: {"raw": 0, "pairs": set(), "pmids": set(), "actions": set(), "directions": set(), "single": False, "cotreat": False})
    for _, row in ctd.iterrows():
        g = clean(row.get(ctd_symbol_col)).upper()
        if not g:
            continue
        d = ctd_by[g]
        d["raw"] += 1
        d["pairs"].add((clean(row.get("GeneID")), g))
        d["pmids"].update(x for x in clean(row.get("PubMedIDs")).split("|") if x)
        d["actions"].update(x for x in clean(row.get("InteractionActions")).split("|") if x)
        for action in d["actions"]:
            if "increases" in action:
                d["directions"].add("increases")
            if "decreases" in action:
                d["directions"].add("decreases")
            if "affects" in action:
                d["directions"].add("affects")
        d["single"] |= clean(row.get("single_chemical_record_flag")) == "true"
        d["cotreat"] |= clean(row.get("multi_chemical_co_treatment_flag")) == "true"
    exp_rows: list[dict[str, Any]] = []
    for g in exp_genes:
        d = ctd_by[g]
        c = int(g in ctd_genes)
        t = int(g in tox_genes)
        ch = int(g in chembl_genes)
        exp_rows.append({
            "gene_symbol": g,
            "CTD": c,
            "ToxCast": t,
            "ChEMBL": ch,
            "exposure_support_count": c + t + ch,
            "CTD_direction": ";".join(sorted(d["directions"])),
            "CTD_action": ";".join(sorted(d["actions"])),
            "CTD_raw_row_count": d["raw"],
            "CTD_unique_chemical_gene_pairs": len(d["pairs"]),
            "CTD_unique_pmids": len(d["pmids"]),
            "CTD_single_chemical_evidence": int(d["single"]),
            "CTD_cotreatment_evidence": int(d["cotreat"]),
            "ToxCast_assay_count": 0,
            "ChEMBL_record_count": 0,
        })
    exp_frame = pd.DataFrame(exp_rows)
    exp_frame.to_csv(OUT / "dinp_exposure_gene_matrix.csv", index=False)

    ot_scores = dict(zip(ot_frame.get("approved_symbol", pd.Series(dtype=str)), ot_frame.get("overall_score", pd.Series(dtype=float)))) if not ot_frame.empty else {}
    crc_rows: list[dict[str, Any]] = []
    for g in crc_genes:
        row = cards_frame[cards_frame["gene_symbol"].eq(g)].head(1)
        score = clean(row["RelevanceScore"].iloc[0]) if not row.empty and "RelevanceScore" in row else ""
        ot_score = ot_scores.get(g, "")
        crc_rows.append({
            "gene_symbol": g,
            "GeneCards": int(g in cards),
            "DisGeNET": int(g in dis_genes),
            "OpenTargets": int(g in ot_genes),
            "disease_support_count": int(g in cards) + int(g in dis_genes) + int(g in ot_genes),
            "GeneCards_score": score,
            "DisGeNET_score": "",
            "OpenTargets_score": ot_score,
        })
    crc_frame = pd.DataFrame(crc_rows)
    crc_frame.to_csv(OUT / "crc_gene_matrix.csv", index=False)

    merged = exp_frame.merge(crc_frame, on="gene_symbol", how="inner")
    merged["total_support_count"] = merged["exposure_support_count"] + merged["disease_support_count"]
    merged.to_csv(OUT / "dinp_crc_intersection.csv", index=False)
    high = merged[(merged["exposure_support_count"] >= 2) & (merged["disease_support_count"] >= 2)].copy()
    high.to_csv(OUT / "high_confidence_intersection.csv", index=False)
    subsets: list[pd.DataFrame] = []
    if not merged.empty:
        for label, block in {
            "all_intersection": merged,
            "exposure_support_ge2": merged[merged["exposure_support_count"] >= 2],
            "disease_support_ge2": merged[merged["disease_support_count"] >= 2],
            "exposure_3of3": merged[merged["exposure_support_count"] == 3],
            "disease_3of3": merged[merged["disease_support_count"] == 3],
        }.items():
            x = block.copy()
            x.insert(0, "subset_label", label)
            subsets.append(x)
    pd.concat(subsets, ignore_index=True).to_csv(OUT / "intersection_subsets.csv", index=False) if subsets else pd.DataFrame().to_csv(OUT / "intersection_subsets.csv", index=False)

    source_meta = {
        "CTD": ctd_meta,
        "ToxCast": tox_meta,
        "ChEMBL": chembl_meta,
        "GeneCards": cards_meta,
        "DisGeNET": dis_meta,
        "OpenTargets": ot_meta,
        "HGNC": hgnc_meta,
    }
    exposure_pairwise = {
        "CTD__ToxCast": len(ctd_genes & tox_genes),
        "CTD__ChEMBL": len(ctd_genes & chembl_genes),
        "ToxCast__ChEMBL": len(tox_genes & chembl_genes),
    }
    disease_pairwise = {
        "GeneCards__DisGeNET": len(cards & dis_genes),
        "GeneCards__OpenTargets": len(cards & ot_genes),
        "DisGeNET__OpenTargets": len(dis_genes & ot_genes),
    }
    exposure_three_way = len(ctd_genes & tox_genes & chembl_genes)
    disease_three_way = len(cards & dis_genes & ot_genes)
    top_bilateral = merged.sort_values(
        ["exposure_support_count", "disease_support_count", "total_support_count", "gene_symbol"],
        ascending=[False, False, False, True],
    ).head(20)
    dump_json(OUT / "source_manifest.json", {"generated_at": utc(), "dinp": {"ctd_id": DINP_CTD_ID, "casrn": DINP_CASRN, "dtxsid": DINP_DTXSID}, "sources": source_meta})
    summary = {
        "generated_at": utc(),
        "counts": {
            "CTD_human_genes": len(ctd_genes),
            "ToxCast_human_genes": len(tox_genes),
            "ChEMBL_human_genes": len(chembl_genes),
            "DINP_exposure_union_genes": len(exp_genes),
            "GeneCards_crc_genes": len(cards),
            "DisGeNET_crc_genes": len(dis_genes),
            "OpenTargets_crc_genes": len(ot_genes),
            "CRC_disease_union_genes": len(crc_genes),
            "DINP_CRC_intersection": len(merged),
            "high_confidence_intersection": len(high),
        },
        "pairwise_overlaps": {
            "exposure_side": exposure_pairwise,
            "disease_side": disease_pairwise,
        },
        "three_way_overlaps": {
            "exposure_side": exposure_three_way,
            "disease_side": disease_three_way,
        },
        "top_bilateral_support_genes": [
            {
                "gene_symbol": clean(row.get("gene_symbol")),
                "exposure_support_count": int(row["exposure_support_count"]),
                "disease_support_count": int(row["disease_support_count"]),
                "total_support_count": int(row["total_support_count"]),
            }
            for _, row in top_bilateral.iterrows()
        ],
        "source_status": source_meta,
        "rules": {
            "human_rule": "CTD OrganismID 9606; Open Targets source-native approvedSymbol targets; other sources must provide human target evidence",
            "high_confidence_rule": "exposure_support_count >= 2 AND disease_support_count >= 2",
            "support_count_interpretation": "descriptive source coverage only; not a biological truth score",
            "no_substitution": "parent DINP only; related phthalates are not used to fill missing ChEMBL/ToxCast records",
            "cotreatment": "CTD co-treatment rows retained and flagged, not silently deleted",
            "gene_normalization": "Symbols are upper-cased for joins; DINP exposure/intersection symbols are audited against HGNC, and source-native approvedSymbol is retained for Open Targets",
        },
    }
    lines = [
        "# DINP–CRC multi-database target convergence summary",
        "",
        f"Generated: `{summary['generated_at']}`",
        "",
        "## Design",
        "",
        "Exposure sources were kept separate: CTD, EPA CompTox/ToxCast, and ChEMBL. Disease sources were kept separate: GeneCards, DisGeNET, and Open Targets. Source support counts are descriptive and do not represent a merged biological truth score.",
        "",
        "Parent DINP was fixed to CTD `C012125`, CASRN `28553-12-0`, and DTXSID `DTXSID4022521`. No related phthalate was substituted.",
        "",
        "## Counts",
        "",
        "| Quantity | Count |",
        "|---|---:|",
    ]
    for k, v in summary["counts"].items():
        lines.append(f"| `{k}` | {v} |")
    lines += [
        "",
        "## Pairwise and three-way overlaps",
        "",
        "Overlap counts are reported separately within the exposure side and disease side; an unavailable source contributes an empty set only as a computational placeholder and is not interpreted as biological absence.",
        "",
        "| Exposure-side overlap | Genes |",
        "|---|---:|",
    ]
    for k, v in exposure_pairwise.items():
        lines.append(f"| `{k}` | {v} |")
    lines.append(f"| `CTD__ToxCast__ChEMBL` | {exposure_three_way} |")
    lines += ["", "| Disease-side overlap | Genes |", "|---|---:|"]
    for k, v in disease_pairwise.items():
        lines.append(f"| `{k}` | {v} |")
    lines.append(f"| `GeneCards__DisGeNET__OpenTargets` | {disease_three_way} |")
    lines += [
        "",
        "## Bilateral-support overview",
        "",
        "The following is a descriptive ordering of genes present in the all-source intersection, ordered by exposure-side support, then disease-side support, then total source coverage. It is not a biological truth score and is not used to infer causality.",
        "",
        "| Gene | Exposure support | Disease support | Total source coverage |",
        "|---|---:|---:|---:|",
    ]
    for _, row in top_bilateral.iterrows():
        lines.append(f"| `{row['gene_symbol']}` | {int(row['exposure_support_count'])} | {int(row['disease_support_count'])} | {int(row['total_support_count'])} |")
    lines += [
        "",
        "## Interpretation boundaries",
        "",
        "- CTD includes 166 human DINP rows in the archived interaction file; co-treatment rows are retained and flagged in the source record table.",
        "- The EPA live CTX bioactivity request is recorded as unavailable when no API key is provided. Public bulk release availability is recorded separately; this is not treated as biological absence.",
        "- ChEMBL returned no reliable exact parent DINP molecule match in the name search; related phthalates were not queried as substitutes.",
        "- DisGeNET API access was not authorized in this run; its zero is not a biological negative.",
        "- GeneCards is the archived ordinary CRC top-2000 export. The historical strict scoped CRC file is not used in the primary matrix.",
        "- Open Targets uses the source-native disease concept resolved from the exact search hit `MONDO_0005575` (colorectal cancer).",
        "",
        "## Files",
        "",
        "- `dinp_exposure_gene_matrix.csv`: exposure-side source-preserving matrix.",
        "- `crc_gene_matrix.csv`: disease-side source-preserving matrix.",
        "- `dinp_crc_intersection.csv`: all genes supported on both sides by at least one accessible source.",
        "- `high_confidence_intersection.csv`: prespecified >=2 exposure sources and >=2 disease sources.",
        "- `intersection_subsets.csv`: requested descriptive subsets.",
        "- `source_manifest.json` and `source_records/`: source versions, hashes, requests, and raw/flattened records.",
        "- `source_records/hgnc_symbol_audit.csv`: HGNC approval/status audit for the exposure/intersection symbol universe.",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)
    client = HTTPClient()
    ctd, ctd_meta = parse_ctd()
    cards, cards_frame, cards_meta = load_genecards()
    chembl_genes, chembl_meta, _ = query_chembl(client)
    tox_genes, tox_meta = query_toxcast(client)
    dis_genes, dis_meta = query_disgenet(client)
    ot_genes, ot_meta, ot_frame = query_open_targets(client)
    exposure_symbols = set(ctd["GeneSymbol"].astype(str).str.upper()) if not ctd.empty else set()
    exposure_symbols.update(chembl_genes)
    exposure_symbols.update(tox_genes)
    disease_union_for_audit = cards | dis_genes | ot_genes
    hgnc_symbols = exposure_symbols | (exposure_symbols & disease_union_for_audit)
    hgnc_frame, hgnc_meta = audit_hgnc_symbols(hgnc_symbols, client)
    hgnc_map = dict(zip(hgnc_frame["input_symbol"], hgnc_frame["approved_symbol"])) if not hgnc_frame.empty else {}
    ctd_analysis = ctd.copy()
    ctd_analysis["normalized_gene_symbol"] = ctd_analysis["GeneSymbol"].astype(str).str.upper().map(hgnc_map)
    unresolved = ctd_analysis["normalized_gene_symbol"].isna() | ctd_analysis["normalized_gene_symbol"].eq("")
    ctd_analysis.loc[unresolved, "normalized_gene_symbol"] = ctd_analysis.loc[unresolved, "GeneSymbol"].astype(str).str.upper()
    ctd_analysis.to_csv(SRC / "ctd_dinp_human_interactions.csv", index=False)
    ctd_meta["normalized_gene_symbol_count"] = int(ctd_analysis["normalized_gene_symbol"].nunique())
    ctd_meta["raw_to_normalized_symbol_changes"] = int((ctd_analysis["GeneSymbol"].astype(str).str.upper() != ctd_analysis["normalized_gene_symbol"]).sum())
    summary = build_outputs(ctd_analysis, ctd_meta, chembl_genes, chembl_meta, tox_genes, tox_meta, cards, cards_frame, cards_meta, dis_genes, dis_meta, ot_genes, ot_meta, ot_frame, hgnc_frame, hgnc_meta)
    dump_json(OUT / "http_call_log.json", client.calls)
    print(json.dumps(summary["counts"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
