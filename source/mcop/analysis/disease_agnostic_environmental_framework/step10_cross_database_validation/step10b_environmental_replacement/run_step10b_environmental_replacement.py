#!/usr/bin/env python3
"""Step 10B-E: source-specific environmental replacement audit.

ChEMBL, BindingDB, and PubChem BioAssay remain separate evidence layers.
Absence or unresolved identity is never treated as biological absence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests

VERSION = "1.0.0"
UA = "whynot17-step10b-environmental/1.0"
ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
FRAMEWORK = ROOT / "analysis" / "disease_agnostic_environmental_framework"
STEP7 = FRAMEWORK / "step07_genecard_convergence"
STEP10 = FRAMEWORK / "step10_cross_database_validation"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


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


def cid(value: object) -> str:
    m = re.search(r"(?:CID:)?\s*(\d+)", s(value), flags=re.I)
    return m.group(1) if m else ""


class Client:
    def __init__(self) -> None:
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": UA, "Accept": "application/json"})
        self.calls: list[dict[str, object]] = []

    def get(self, url: str, params: dict[str, object] | None = None, delay: float = 0.0, timeout: tuple[int, int] = (15, 90), attempts: int = 3) -> tuple[object | None, dict[str, object]]:
        if delay:
            time.sleep(delay)
        last = ""
        attempts = max(1, attempts)
        for attempt in range(attempts):
            try:
                r = self.session.get(url, params=params, timeout=timeout)
                meta = {"status_code": r.status_code, "attempt": attempt + 1, "url": r.url, "response_sha256": hashlib.sha256(r.content).hexdigest()}
                if r.status_code == 200:
                    payload = r.json()
                    self.calls.append({"endpoint": url, **meta})
                    return payload, meta
                last = f"HTTP {r.status_code}: {r.text[:300]}"
            except Exception as exc:
                last = f"{type(exc).__name__}: {exc}"
            if attempt < attempts - 1:
                time.sleep(1.5 * (attempt + 1))
        meta = {"status_code": None, "attempt": attempts, "url": url, "error": last}
        self.calls.append({"endpoint": url, **meta})
        return None, meta


def load_inputs() -> pd.DataFrame:
    candidates = pd.read_csv(FRAMEWORK / "t2d_exposure_opportunity/01_candidate_master/unique_candidate_chemicals.csv", dtype=str).fillna("")
    ctd = pd.read_csv(ROOT / "outputs/environmental_toxicology_crc_phase1_chemical_classification.csv", dtype=str, usecols=["ChemicalID", "ChemicalName", "CasRN", "PubChemCID", "DTXSID", "InChIKey", "n_unique_chemical_gene_pairs", "n_unique_pmids"]).fillna("").drop_duplicates("ChemicalID")
    cmap = pd.read_csv(STEP7 / "t2d_step7_cluster_chemical_map.csv", dtype=str).fillna("")
    agg = cmap.groupby("chemical_id", as_index=False).agg(
        step7_cluster_ids=("cluster_id", lambda x: ";".join(sorted(set(v for v in x if v)))),
        step7_biomarkers=("biomarker", lambda x: ";".join(sorted(set(v for v in x if v)))),
    )
    out = candidates.merge(ctd, left_on="chemical_id", right_on="ChemicalID", how="left").merge(agg, on="chemical_id", how="left")
    out["pubchem_cid_ctd"] = out["PubChemCID"].map(cid)
    out["candidate_source_row"] = np.arange(1, len(out) + 1)
    cols = ["candidate_source_row", "chemical_id", "chemical_name", "chemical_class", "positive_biomarkers", "exposure_axes", "step7_cluster_ids", "step7_biomarkers", "ChemicalName", "CasRN", "pubchem_cid_ctd", "DTXSID", "InChIKey", "n_unique_chemical_gene_pairs", "n_unique_pmids", "mapping_gate_disposition", "mapping_gate_reasons", "manual_review_required", "cycle_lists", "matrices"]
    return out[cols].rename(columns={"ChemicalName": "ctd_chemical_name", "CasRN": "ctd_casrn", "DTXSID": "ctd_dtxsid", "InChIKey": "ctd_inchikey"})


def pubchem_resolve(client: Client, row: pd.Series) -> tuple[str, dict[str, object], dict[str, object]]:
    resolved = s(row.get("pubchem_cid_ctd"))
    method = "ctd_pubchem_cid" if resolved else ""
    meta: dict[str, object] = {}
    for field, rule in (("ctd_inchikey", "ctd_inchikey"), ("ctd_casrn", "ctd_casrn")):
        if resolved or not s(row.get(field)):
            continue
        namespace = "inchikey" if field == "ctd_inchikey" else "name"
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{namespace}/{quote(s(row[field]))}/cids/JSON"
        payload, call_meta = client.get(url, delay=0.22, timeout=(10, 30), attempts=2)
        meta[rule] = call_meta
        values = payload.get("IdentifierList", {}).get("CID", []) if isinstance(payload, dict) else []
        if values:
            resolved, method = str(values[0]), rule
    props: dict[str, object] = {}
    prop_meta: dict[str, object] = {}
    if resolved:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{quote(resolved)}/property/InChIKey,CanonicalSMILES,IsomericSMILES/JSON"
        payload, prop_meta = client.get(url, delay=0.22, timeout=(10, 30), attempts=2)
        values = payload.get("PropertyTable", {}).get("Properties", []) if isinstance(payload, dict) else []
        if values:
            props = values[0]
    return resolved, {"resolution_method": method, **props}, {"resolve": meta, "property": prop_meta}


def chembl_molecules(client: Client, inchikey: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    if not inchikey:
        return [], {"status": "not_queried", "reason": "no_structure_identifier"}
    return_payload, meta = client.get("https://www.ebi.ac.uk/chembl/api/data/molecule.json", {"molecule_structures__standard_inchi_key": inchikey, "limit": 100}, delay=0.12)
    return (return_payload.get("molecules", []) if isinstance(return_payload, dict) else []), meta


def chembl_activity(client: Client, chembl_id: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    fields = "molecule_chembl_id,target_chembl_id,target_organism,target_type,assay_chembl_id,pchembl_value,standard_type,standard_value,standard_units"
    payload, meta = client.get("https://www.ebi.ac.uk/chembl/api/data/activity.json", {"molecule_chembl_id": chembl_id, "limit": 1000, "only": fields}, delay=0.12)
    return (payload.get("activities", []) if isinstance(payload, dict) else []), meta


def bindingdb(client: Client, smiles: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    if not smiles:
        return [], {"status": "not_queried", "reason": "no_canonical_smiles"}
    payload, meta = client.get("https://bindingdb.org/rest/getTargetByCompound", {"smiles": smiles, "cutoff": "1.0", "response": "application/json"}, delay=0.75)
    if not isinstance(payload, dict):
        return [], meta
    values = payload.get("getLindsByUniprotResponse", {}).get("bdb.affinities", [])
    return (values if isinstance(values, list) else []), meta


def pubchem_aids(client: Client, resolved_cid: str) -> tuple[list[int], dict[str, object]]:
    if not resolved_cid:
        return [], {"status": "not_queried", "reason": "no_pubchem_cid"}
    payload, meta = client.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{quote(resolved_cid)}/aids/JSON", delay=0.22, timeout=(10, 20), attempts=1)
    info = payload.get("InformationList", {}).get("Information", []) if isinstance(payload, dict) else []
    values = info[0].get("AID", []) if info else []
    return [int(v) for v in values], meta


def spearman(a: pd.Series, b: pd.Series) -> tuple[float, int]:
    mask = a.notna() & b.notna()
    n = int(mask.sum())
    return (float(a[mask].rank().corr(b[mask].rank())) if n >= 3 else float("nan"), n)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-candidates", type=int, default=0)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    started = utc()
    client = Client()
    crosswalk = load_inputs()
    if args.max_candidates:
        crosswalk = crosswalk.head(args.max_candidates).copy()
    crosswalk.to_csv(OUT / "step10b_environmental_candidate_crosswalk.csv", index=False)

    chembl_status, chembl_status_meta = client.get("https://www.ebi.ac.uk/chembl/api/data/status.json")
    bindingdb_page, bindingdb_page_meta = client.get("https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp")
    if bindingdb_page is not None:
        (OUT / "bindingdb_download_page_snapshot.json").write_text(json.dumps(bindingdb_page, ensure_ascii=False, indent=2), encoding="utf-8")

    rows: list[dict[str, object]] = []
    for i, row in crosswalk.iterrows():
        resolved_cid, props, pub_meta = pubchem_resolve(client, row)
        inchikey = s(props.get("InChIKey")) or s(row.get("ctd_inchikey"))
        smiles = s(props.get("ConnectivitySMILES")) or s(props.get("CanonicalSMILES")) or s(props.get("IsomericSMILES"))
        molecules, molecule_meta = chembl_molecules(client, inchikey)
        chembl_ids = sorted(set(s(m.get("molecule_chembl_id")) for m in molecules if s(m.get("molecule_chembl_id"))))
        activities: list[dict[str, object]] = []
        activity_meta: list[dict[str, object]] = []
        for chembl_id in chembl_ids:
            values, meta = chembl_activity(client, chembl_id)
            activities.extend(values)
            activity_meta.append(meta)
        human_activity = [a for a in activities if s(a.get("target_organism")).lower() in {"human", "homo sapiens"} or "homo sapiens" in s(a.get("target_organism")).lower()]
        binding_hits, binding_meta = bindingdb(client, smiles)
        human_binding = [h for h in binding_hits if s(h.get("bdb.species")).lower() == "human"]
        binding_targets = sorted(set(s(h.get("bdb.target")) for h in human_binding if s(h.get("bdb.target"))))
        aids, aids_meta = pubchem_aids(client, resolved_cid)
        rows.append({
            **row.to_dict(),
            "pubchem_cid_resolved": resolved_cid,
            "pubchem_resolution_method": s(props.get("resolution_method")),
            "pubchem_inchikey": inchikey,
            "pubchem_canonical_smiles": smiles,
            "pubchem_resolution_status": "resolved" if resolved_cid else "unresolved",
            "pubchem_aid_count": len(aids),
            "pubchem_assay_evidence_status": "observed" if aids else ("not_observed" if resolved_cid else "unresolved"),
            "pubchem_target_annotation_status": "not_derived_from_CID_AID_membership",
            "chembl_molecule_count": len(chembl_ids),
            "chembl_ids": ";".join(chembl_ids),
            "chembl_activity_count": len(activities),
            "chembl_human_activity_count": len(human_activity),
            "chembl_human_target_count": len(set(s(a.get("target_chembl_id")) for a in human_activity if s(a.get("target_chembl_id")))),
            "chembl_evidence_status": "observed" if human_activity else ("resolved_no_human_activity" if chembl_ids else "not_observed_or_unresolved"),
            "bindingdb_match_rule": "Tanimoto_cutoff_1.0_via_getTargetByCompound",
            "bindingdb_affinity_row_count": len(binding_hits),
            "bindingdb_human_affinity_row_count": len(human_binding),
            "bindingdb_human_target_count": len(binding_targets),
            "bindingdb_human_targets": ";".join(binding_targets),
            "bindingdb_evidence_status": "observed" if human_binding else ("observed_nonhuman_only" if binding_hits else "not_observed_or_unresolved"),
            "api_meta_pubchem": dumps(pub_meta),
            "api_meta_chembl": dumps({"molecule": molecule_meta, "activity": activity_meta}),
            "api_meta_bindingdb": dumps(binding_meta),
            "api_meta_pubchem_aids": dumps(aids_meta),
        })
        if (i + 1) % 10 == 0:
            print(f"processed {i + 1}/{len(crosswalk)}", flush=True)
    records = pd.DataFrame(rows)
    scores = {"E1_ChEMBL": "chembl_human_activity_count", "E2_BindingDB": "bindingdb_human_affinity_row_count", "E3_PubChem_BioAssay": "pubchem_aid_count"}
    topk: list[dict[str, object]] = []
    for source, score_col in scores.items():
        ranked = records.sort_values([score_col, "chemical_id"], ascending=[False, True]).reset_index(drop=True)
        ranked["source_rank"] = np.arange(1, len(ranked) + 1)
        positive = ranked[ranked[score_col].fillna(0) > 0]
        for k in (5, 10, 15):
            top = positive.head(k)
            topk.append({"source_id": source, "top_k": k, "positive_candidates": len(positive), "observed_top_k": len(top), "chemical_ids": ";".join(top["chemical_id"].astype(str)), "interpretation": "within-source descriptive rank; not cross-source potency ranking"})
        records[f"{source}_rank"] = ranked.set_index("chemical_id")["source_rank"].reindex(records["chemical_id"]).to_numpy()
    records.to_csv(OUT / "step10b_environmental_source_records.csv", index=False)
    pd.DataFrame(topk).to_csv(OUT / "step10b_environmental_fixed_topk_retention.csv", index=False)
    support = list(scores.values())
    records["n_sources_with_observed_evidence"] = (records[support].fillna(0) > 0).sum(axis=1)
    records["observed_evidence_sources"] = records.apply(lambda r: ";".join(source for source, col in zip(scores, support) if float(r[col] or 0) > 0), axis=1)
    records[["chemical_id", "chemical_name", "step7_cluster_ids", "n_sources_with_observed_evidence", "observed_evidence_sources"] + support].to_csv(OUT / "step10b_environmental_source_dropout.csv", index=False)
    concordance = []
    pairs = [("E1_ChEMBL", support[0], "E2_BindingDB", support[1]), ("E1_ChEMBL", support[0], "E3_PubChem_BioAssay", support[2]), ("E2_BindingDB", support[1], "E3_PubChem_BioAssay", support[2])]
    for s1, c1, s2, c2 in pairs:
        rho, n = spearman(records[c1].astype(float), records[c2].astype(float))
        concordance.append({"source_1": s1, "source_2": s2, "score_1": c1, "score_2": c2, "spearman_rho": rho, "n_candidates": n, "interpretation": "descriptive count concordance; source semantics remain non-equivalent"})
    for source, score_col in scores.items():
        for ctd_col, label in [("n_unique_chemical_gene_pairs", "CTD_human_chemical_gene_pair_count"), ("n_unique_pmids", "CTD_unique_PMID_count")]:
            rho, n = spearman(records[score_col].astype(float), pd.to_numeric(records[ctd_col], errors="coerce"))
            concordance.append({"source_1": source, "source_2": label, "score_1": score_col, "score_2": ctd_col, "spearman_rho": rho, "n_candidates": n, "interpretation": "descriptive source-record versus CTD annotation-burden concordance; not causal and not used for promotion"})
    pd.DataFrame(concordance).to_csv(OUT / "step10b_environmental_rank_concordance.csv", index=False)

    summary = []
    for source, col in scores.items():
        summary.append({"source_id": source, "candidate_count": len(records), "positive_record_count": int((records[col].fillna(0) > 0).sum()), "positive_record_fraction": float((records[col].fillna(0) > 0).mean()), "record_count_field": col, "target_filter_status": "human applied" if source != "E3_PubChem_BioAssay" else "not inferable from CID_AID membership", "absence_interpretation": "not biological negative"})
    pd.DataFrame(summary).to_csv(OUT / "step10b_environmental_source_summary.csv", index=False)

    snapshot = {
        "lock_type": "STEP10B_E_ENVIRONMENTAL_SOURCE_SNAPSHOT",
        "script_version": VERSION,
        "retrieval_started_utc": started,
        "retrieval_finished_utc": utc(),
        "candidate_count_analyzed": len(records),
        "frozen_input": str(STEP10 / "STEP10_FROZEN_SOURCE_SET.json"),
        "sources": {
            "E1": {"name": "ChEMBL", "reference": "https://www.ebi.ac.uk/chembl/", "status_endpoint": "https://www.ebi.ac.uk/chembl/api/data/status.json", "status_payload": chembl_status, "status_response_meta": chembl_status_meta, "human_rule": "target_organism contains Homo sapiens or Human", "activity_rule": "no potency threshold; all returned activity rows retained", "missingness": "absence is not biological negative"},
            "E2": {"name": "BindingDB", "reference": "https://www.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp", "download_page": "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp", "download_page_response_meta": bindingdb_page_meta, "download_page_sha256": sha256(OUT / "bindingdb_download_page_snapshot.json") if (OUT / "bindingdb_download_page_snapshot.json").exists() else None, "human_rule": "bdb.species equals Human", "query_rule": "Tanimoto cutoff 1.0 through getTargetByCompound; not represented as a stable ligand-ID join", "missingness": "absence is not biological negative"},
            "E3": {"name": "PubChem BioAssay", "reference": "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial", "query_rule": "CID to AID membership; assay count is coverage, not target count", "human_rule": "not inferred from CID to AID membership", "missingness": "absence is not biological negative"},
        },
        "api_call_count": len(client.calls),
        "api_call_manifest": client.calls,
        "input_hashes": {"crosswalk": sha256(OUT / "step10b_environmental_candidate_crosswalk.csv")},
        "output_hashes": {},
        "analysis_boundary": "Source replacement outputs do not change the frozen candidate universe, epidemiologic results, or Tier assignments.",
    }
    for path in OUT.glob("step10b_environmental_*.csv"):
        snapshot["output_hashes"][path.name] = sha256(path)
    (OUT / "STEP10B_E_SOURCE_SNAPSHOT.json").write_text(json.dumps(snapshot, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "STEP10B_E_API_CALL_MANIFEST.json").write_text(json.dumps({"generated_utc": utc(), "calls": client.calls}, ensure_ascii=False, indent=2), encoding="utf-8")
    report = "\n".join([
        "# Step 10B-E — Environmental knowledge-source replacement audit",
        "",
        f"Generated: `{utc()}`",
        f"Candidates analyzed: **{len(records)}** (inherited frozen Step 4/7 universe)",
        "",
        "ChEMBL measured bioactivity, BindingDB affinity, and PubChem BioAssay compound-to-assay membership were queried as separate evidence layers. They were not merged into a chemical-gene edge list.",
        "",
        "Missing or unresolved records are reported as not observed/unresolved and are not interpreted as biological negatives. PubChem CID-to-AID membership does not establish a human protein target, so a PubChem target count is not fabricated.",
        "",
        "Exact source metadata, query rules, API response hashes, input/output hashes, and call records are in `STEP10B_E_SOURCE_SNAPSHOT.json` and `STEP10B_E_API_CALL_MANIFEST.json`.",
        "",
        "Status: complete as a source-specific coverage/replacement audit; mechanistic target inference requires additional source-native assay parsing.",
    ])
    (OUT / "STEP10B_E_ENVIRONMENTAL_REPLACEMENT_REPORT.md").write_text(report, encoding="utf-8")
    print(f"completed Step 10B-E: {len(records)} candidates; {len(client.calls)} API calls", flush=True)


if __name__ == "__main__":
    main()
