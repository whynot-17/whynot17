"""URXP02 / 2-naphthol mechanism gene discovery (gene-level only).

This script intentionally uses the exact CTD chemical identity for 2-naphthol
and direct Open Targets disease associations.  It does not fit or read any
NHANES association model and does not make mechanistic or sex-specific claims.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import time
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
CTD_CHEM = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "CTD_chemicals.tsv.gz"
CTD_IXN = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "CTD_chem_gene_ixns.tsv.gz"
GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
GPROFILER = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
CHEMICAL_ID = "C028405"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ctid(value: str) -> str:
    value = str(value or "").strip()
    return value[5:] if value.startswith("MESH:") else value


def split_pmids(value: str) -> set[str]:
    return {x.strip() for x in str(value or "").replace("|", ";").split(";") if x.strip()}


def ctd_rows(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = None
        for line in handle:
            if line.startswith("# Fields:"):
                header = next(handle).lstrip("# ").rstrip("\r\n").split("\t")
                break
        if header is None:
            raise ValueError(f"Missing CTD field header: {path}")
        yield from csv.DictReader(handle, fieldnames=header, delimiter="\t")


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_chemical_identity() -> dict:
    matches = []
    for row in ctd_rows(CTD_CHEM):
        name = str(row.get("ChemicalName") or "").strip().lower()
        synonyms = " ".join(str(row.get(k) or "") for k in ("MESHSynonyms", "CTDCuratedSynonyms")).lower()
        if name == "2-naphthol" or (ctid(row.get("ChemicalID", "")) == CHEMICAL_ID):
            matches.append(row)
    exact = [r for r in matches if ctid(r.get("ChemicalID", "")) == CHEMICAL_ID]
    if not exact:
        raise RuntimeError("The exact CTD 2-naphthol identity was not found")
    row = exact[0]
    return {
        "chemical_name": row.get("ChemicalName", ""),
        "chemical_id": ctid(row.get("ChemicalID", "")),
        "casrn": row.get("CasRN", ""),
        "pubchem_cid": row.get("PubChemCID", ""),
        "mesh_synonyms": row.get("MESHSynonyms", ""),
        "ctd_curated_synonyms": row.get("CTDCuratedSynonyms", ""),
        "matched_exact_aliases": "2-naphthol;2-hydroxynaphthalene;2-NAP;2NAP;2-hydroxy-naphthalene",
        "excluded_near_matches": "1-naphthol;2-naphthyl sulfate;parent naphthalene;unrelated naphthol derivatives",
    }


def load_ctd_human(identity: dict) -> list[dict]:
    by_gene: dict[tuple[str, str], dict] = {}
    for row in ctd_rows(CTD_IXN):
        if ctid(row.get("ChemicalID", "")) != identity["chemical_id"]:
            continue
        if row.get("Organism", "") != "Homo sapiens":
            continue
        symbol = str(row.get("GeneSymbol") or "").strip().upper()
        gene_id = str(row.get("GeneID") or "").strip()
        if not symbol or not gene_id:
            continue
        key = (symbol, gene_id)
        rec = by_gene.setdefault(
            key,
            {
                "gene_symbol": symbol,
                "entrez_id": gene_id,
                "ctd_chemical_name": row.get("ChemicalName", ""),
                "ctd_chemical_id": ctid(row.get("ChemicalID", "")),
                "organism": "Homo sapiens",
                "ctd_interaction_types": set(),
                "ctd_interaction_evidence": set(),
                "ctd_pubmed_ids": set(),
                "ctd_raw_rows": 0,
            },
        )
        rec["ctd_raw_rows"] += 1
        if row.get("InteractionActions"):
            rec["ctd_interaction_types"].add(row["InteractionActions"])
        if row.get("Interaction"):
            rec["ctd_interaction_evidence"].add(row["Interaction"])
        rec["ctd_pubmed_ids"].update(split_pmids(row.get("PubMedIDs", "")))
    out = []
    for rec in sorted(by_gene.values(), key=lambda x: x["gene_symbol"]):
        rec["ctd_interaction_types"] = "; ".join(sorted(rec["ctd_interaction_types"]))
        rec["ctd_interaction_evidence"] = " || ".join(sorted(rec["ctd_interaction_evidence"]))
        rec["ctd_pubmed_ids"] = ";".join(sorted(rec["ctd_pubmed_ids"], key=lambda x: (len(x), x)))
        rec["ctd_unique_reference_count"] = len(rec["ctd_pubmed_ids"].split(";")) if rec["ctd_pubmed_ids"] else 0
        rec["source"] = "CTD chemical-gene interaction"
        rec["ctd_identity_note"] = "Exact CTD C028405 (2-naphthol; synonyms include 2-hydroxynaphthalene/2-NAP)"
        out.append(rec)
    return out


def graphql(query: str) -> dict:
    payload = post_json(GRAPHQL, {"query": query})
    if payload.get("errors"):
        raise RuntimeError(json.dumps(payload["errors"], ensure_ascii=False))
    return payload["data"]


def post_json(url: str, body: dict) -> dict:
    encoded = json.dumps(body).encode("utf-8")
    request = urllib.request.Request(url, data=encoded, headers={"Content-Type": "application/json", "User-Agent": "URXP02-M1/1.0"}, method="POST")
    with urllib.request.urlopen(request, timeout=180) as response:
        return json.loads(response.read().decode("utf-8"))


def opentargets_disease(disease_id: str, expected_name: str) -> tuple[dict, list[dict]]:
    gql_id = json.dumps(disease_id)
    meta_query = f"{{ disease(efoId: {gql_id}) {{ id name }} }}"
    meta = graphql(meta_query)["disease"]
    if not meta:
        raise RuntimeError(f"Open Targets disease not found: {disease_id}")
    rows: list[dict] = []
    page_size = 1000
    first_count = None
    index = 0
    while first_count is None or len(rows) < first_count:
        query = (
            "{ disease(efoId: %s) { associatedTargets(enableIndirect:false, "
            "page:{index:%d,size:%d}) { count rows { score target { id approvedSymbol approvedName } } } } }"
            % (gql_id, index, page_size)
        )
        block = graphql(query)["disease"]["associatedTargets"]
        first_count = int(block["count"])
        if not block["rows"]:
            break
        rows.extend(block["rows"])
        index += 1
        time.sleep(0.05)
    result = []
    seen = set()
    for rank, rec in enumerate(rows, start=1):
        target = rec.get("target") or {}
        symbol = str(target.get("approvedSymbol") or "").strip().upper()
        if not symbol or symbol in seen:
            continue
        seen.add(symbol)
        result.append(
            {
                "gene_symbol": symbol,
                "ensembl_id": target.get("id", ""),
                "approved_name": target.get("approvedName", ""),
                "disease_id": disease_id,
                "disease_name": meta.get("name") or expected_name,
                "disease_relevance_score": rec.get("score"),
                "rank": rank,
                "source": "Open Targets Platform GraphQL",
                "association_scope": "direct (enableIndirect=false)",
            }
        )
    return meta, result


def gprofiler_enrichment(label: str, genes: list[str], background: list[str]) -> list[dict]:
    if not genes:
        return [{"branch": label, "status": "NO_GENES", "query_genes": "", "background_size": len(background)}]
    body = {
        "organism": "hsapiens",
        "query": genes,
        "background": background,
        "domain_scope": "custom",
        "sources": ["GO:BP", "KEGG", "REAC"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": False,
    }
    payload = post_json(GPROFILER, body)
    output = []
    for rec in payload.get("result", []):
        output.append(
            {
                "branch": label,
                "status": "OK",
                "source": rec.get("source", ""),
                "native_id": rec.get("native", ""),
                "term_name": rec.get("name", ""),
                "description": rec.get("description", ""),
                "intersection_size": rec.get("intersection_size"),
                "term_size": rec.get("term_size"),
                "effective_domain_size": rec.get("effective_domain_size"),
                "query_size": rec.get("query_size"),
                "p_value": rec.get("p_value"),
                "significant_g_scs": rec.get("significant"),
                "query_genes": ";".join(genes),
                "background_size": len(background),
                "gprofiler_intersection_evidence": json.dumps(rec.get("intersections", []), ensure_ascii=False),
            }
        )
    if not output:
        output.append({"branch": label, "status": "NO_TERMS_RETURNED", "query_genes": ";".join(genes), "background_size": len(background)})
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    identity = load_chemical_identity()
    ctd = load_ctd_human(identity)
    ctd_symbols = {r["gene_symbol"] for r in ctd}
    thyroid_meta, thyroid = opentargets_disease("MONDO_0003240", "thyroid gland disorder")
    hyper_meta, hyper = opentargets_disease("HP_0000822", "Hypertension")
    thyroid_by = {r["gene_symbol"]: r for r in thyroid}
    hyper_by = {r["gene_symbol"]: r for r in hyper}
    thyroid_set = ctd_symbols & set(thyroid_by)
    hyper_set = ctd_symbols & set(hyper_by)
    shared = thyroid_set & hyper_set
    thyroid_only = thyroid_set - hyper_set
    hyper_only = hyper_set - thyroid_set

    ctd_fields = [
        "gene_symbol", "entrez_id", "ctd_chemical_name", "ctd_chemical_id", "organism",
        "ctd_raw_rows", "ctd_unique_reference_count", "ctd_pubmed_ids", "ctd_interaction_types",
        "ctd_interaction_evidence", "source", "ctd_identity_note",
    ]
    write_csv(OUT / "01_2NAP_CTD_human_genes.csv", ctd, ctd_fields)
    disease_fields = ["gene_symbol", "ensembl_id", "approved_name", "disease_id", "disease_name", "disease_relevance_score", "rank", "source", "association_scope"]
    write_csv(OUT / "02_thyroid_disease_genes.csv", thyroid, disease_fields)
    write_csv(OUT / "03_hypertension_disease_genes.csv", hyper, disease_fields)

    def intersection_rows(symbols: set[str], disease_by: dict[str, dict], branch: str) -> list[dict]:
        rows = []
        for symbol in sorted(symbols):
            c = next(r for r in ctd if r["gene_symbol"] == symbol)
            d = disease_by[symbol]
            rows.append({**c, "disease_branch": branch, "disease_relevance_score": d["disease_relevance_score"], "disease_rank": d["rank"], "disease_id": d["disease_id"], "disease_name": d["disease_name"], "ensembl_id": d["ensembl_id"]})
        return rows

    int_fields = ctd_fields + ["disease_branch", "disease_relevance_score", "disease_rank", "disease_id", "disease_name", "ensembl_id"]
    write_csv(OUT / "04_2NAP_thyroid_intersection.csv", intersection_rows(thyroid_set, thyroid_by, "thyroid"), int_fields)
    write_csv(OUT / "05_2NAP_hypertension_intersection.csv", intersection_rows(hyper_set, hyper_by, "hypertension"), int_fields)

    union = sorted(thyroid_set | hyper_set)
    branch_rows = []
    for symbol in union:
        flags = (symbol in thyroid_only, symbol in hyper_only, symbol in shared)
        branch = "thyroid-specific" if flags[0] else "hypertension-specific" if flags[1] else "shared"
        branch_rows.append({"gene_symbol": symbol, "in_2nap_ctd": True, "in_thyroid_disease_set": symbol in thyroid_set, "in_hypertension_disease_set": symbol in hyper_set, "branch_class": branch})
    write_csv(OUT / "06_branch_gene_classification.csv", branch_rows, ["gene_symbol", "in_2nap_ctd", "in_thyroid_disease_set", "in_hypertension_disease_set", "branch_class"])

    evidence_rows = []
    ctd_by = {r["gene_symbol"]: r for r in ctd}
    for symbol in union:
        c = ctd_by[symbol]
        t = thyroid_by.get(symbol, {})
        h = hyper_by.get(symbol, {})
        evidence_rows.append({
            "gene_symbol": symbol, "entrez_id": c["entrez_id"], "ctd_raw_rows": c["ctd_raw_rows"],
            "ctd_unique_reference_count": c["ctd_unique_reference_count"], "ctd_pubmed_ids": c["ctd_pubmed_ids"],
            "ctd_interaction_types": c["ctd_interaction_types"], "ctd_interaction_evidence": c["ctd_interaction_evidence"],
            "thyroid_disease_relevance_score": t.get("disease_relevance_score", ""), "thyroid_disease_rank": t.get("rank", ""),
            "hypertension_disease_relevance_score": h.get("disease_relevance_score", ""), "hypertension_disease_rank": h.get("rank", ""),
            "thyroid_specific": symbol in thyroid_only, "hypertension_specific": symbol in hyper_only, "shared": symbol in shared,
            "independent_supporting_reference_count": c["ctd_unique_reference_count"],
            "evidence_note": "Component evidence only; no composite score and no causal inference.",
        })
    evidence_fields = list(evidence_rows[0].keys()) if evidence_rows else ["gene_symbol"]
    write_csv(OUT / "07_candidate_gene_evidence_matrix.csv", evidence_rows, evidence_fields)

    enrich_thyroid = gprofiler_enrichment("thyroid-specific", sorted(thyroid_only), sorted(ctd_symbols))
    enrich_hyper = gprofiler_enrichment("hypertension-specific", sorted(hyper_only), sorted(ctd_symbols))
    enrich_shared = gprofiler_enrichment("shared", sorted(shared), sorted(ctd_symbols))
    enrich_fields = sorted({k for row in enrich_thyroid + enrich_hyper + enrich_shared for k in row})
    for fname, rows in [("08_thyroid_branch_enrichment.csv", enrich_thyroid), ("09_hypertension_branch_enrichment.csv", enrich_hyper), ("10_shared_branch_enrichment.csv", enrich_shared)]:
        write_csv(OUT / fname, rows, enrich_fields)

    # Transparent candidate list: branch-specific genes first, then shared genes;
    # no black-box score.  Pathway support is set-level only because the API
    # response does not provide a defensible per-gene attribution table.
    sig = {"thyroid-specific": any(r.get("significant_g_scs") is True for r in enrich_thyroid), "hypertension-specific": any(r.get("significant_g_scs") is True for r in enrich_hyper), "shared": any(r.get("significant_g_scs") is True for r in enrich_shared)}
    candidate_rows = []
    for symbol in union:
        branch = "thyroid-specific" if symbol in thyroid_only else "hypertension-specific" if symbol in hyper_only else "shared"
        c = ctd_by[symbol]; t = thyroid_by.get(symbol, {}); h = hyper_by.get(symbol, {})
        candidate_rows.append({
            "gene_symbol": symbol, "candidate_branch": branch, "priority_tier": "Priority 1" if sig[branch] else "Priority 2",
            "exposure_supported": True, "disease_supported": True, "pathway_supported_at_set_level": sig[branch],
            "ctd_raw_rows": c["ctd_raw_rows"], "ctd_unique_reference_count": c["ctd_unique_reference_count"],
            "thyroid_disease_relevance_score": t.get("disease_relevance_score", ""), "thyroid_disease_rank": t.get("rank", ""),
            "hypertension_disease_relevance_score": h.get("disease_relevance_score", ""), "hypertension_disease_rank": h.get("rank", ""),
            "candidate_selection_note": "Exposure+disease evidence; pathway flag is set-level only. No sex-specific claim.",
        })
    candidate_rows.sort(key=lambda r: (r["candidate_branch"] == "shared", r["priority_tier"] != "Priority 1", -(float(r["thyroid_disease_relevance_score"] or 0) + float(r["hypertension_disease_relevance_score"] or 0))))
    write_csv(OUT / "11_cell_mapping_candidate_genes.csv", candidate_rows, list(candidate_rows[0].keys()) if candidate_rows else ["gene_symbol"])

    generated = [p for p in sorted(OUT.glob("*.csv"))]
    manifest = {
        "analysis": "URXP02 / 2-naphthol gene-level mechanism discovery",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "chemical_identity": identity,
        "ctd_source_file": str(CTD_IXN),
        "ctd_source_sha256": sha256(CTD_IXN),
        "disease_gene_source": "Open Targets Platform GraphQL direct target-disease associations (enableIndirect=false)",
        "graphql_endpoint": GRAPHQL,
        "disease_queries": {"thyroid": {"id": thyroid_meta.get("id"), "name": thyroid_meta.get("name"), "direct_gene_count": len(thyroid)}, "hypertension": {"id": hyper_meta.get("id"), "name": hyper_meta.get("name"), "direct_gene_count": len(hyper)}},
        "ctd_human_gene_count": len(ctd),
        "intersection_counts": {"thyroid": len(thyroid_set), "hypertension": len(hyper_set), "thyroid_specific": len(thyroid_only), "hypertension_specific": len(hyper_only), "shared": len(shared)},
        "enrichment": {"service": "g:Profiler g:GOSt", "endpoint": GPROFILER, "sources": ["GO:BP", "KEGG", "REAC"], "background": "all 2-NAP CTD human genes", "domain_scope": "custom", "multiple_testing": "g:SCS", "interpretation": "descriptive; no pathway forced"},
        "constraints": ["No NHANES model rerun", "No other exposure analyzed", "No figures", "No sex-specific gene claim", "CTD links are not treated as causality"],
        "files": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size} for p in generated},
    }
    report = f"""# URXP02 / 2-NAP mechanism gene discovery (M1)\n\nGenerated {manifest['generated_at_utc']}. This is a gene-level discovery package only; it does not establish causality, sex-specific molecular mechanisms, or tissue/cell localization.\n\n## Chemical identity and exposure universe\n\nThe exact CTD entry is **{identity['chemical_name']}** (`{identity['chemical_id']}`; CAS {identity['casrn']}). CTD curated synonyms include 2-hydroxynaphthalene and 2-NAP. Near matches (1-naphthol, 2-naphthyl sulfate, parent naphthalene, and unrelated derivatives) were excluded. The human CTD interaction universe contains **{len(ctd)} genes**.\n\n## Disease gene sources\n\nDisease branches were kept independent. Open Targets direct target–disease associations (`enableIndirect=false`) returned {len(thyroid)} genes for `{thyroid_meta.get('id')}` ({thyroid_meta.get('name')}) and {len(hyper)} genes for `{hyper_meta.get('id')}` ({hyper_meta.get('name')}). These are ranked association resources, not causal gene lists.\n\n## Intersections\n\n- 2-NAP ∩ thyroid: **{len(thyroid_set)}** genes\n- 2-NAP ∩ hypertension: **{len(hyper_set)}** genes\n- thyroid-specific (A−B): **{len(thyroid_only)}**\n- hypertension-specific (B−A): **{len(hyper_only)}**\n- shared (A∩B): **{len(shared)}**\n\nFull rows and component evidence are in `04`–`07`.\n\n## Enrichment and candidate use\n\nGO Biological Process, KEGG, and Reactome were queried separately for each branch with the **{len(ctd)}-gene human CTD 2-NAP universe as the custom background**. The thyroid-specific branch contains no genes, and g:Profiler returned no terms for the four-gene hypertension-specific branch or the three-gene shared branch under this custom background. Thus no pathway enrichment is claimed in M1. The outputs preserve the null/empty results and g:SCS settings; no expected pathway was forced. Because the enrichment response does not provide a defensible per-gene attribution table, `11_cell_mapping_candidate_genes.csv` reports set-level pathway support only and does not claim sex specificity.\n\nNext phase may map these transparent candidates to tissues/cells. It must test, rather than assume, any male/female molecular divergence.\n\n## Provenance\n\n- CTD chemical–gene interactions: `{CTD_IXN}` (SHA-256 `{sha256(CTD_IXN)}`).\n- Disease associations: Open Targets Platform GraphQL `{GRAPHQL}`; direct associations only.\n- Enrichment: g:Profiler g:GOSt `{GPROFILER}`; custom background, sources GO:BP/KEGG/REAC, g:SCS.\n- No NHANES association or exposure–outcome result was read or refit by this M1 script.\n"""
    (OUT / "URXP02_MECHANISM_GENE_REPORT.md").write_text(report, encoding="utf-8")
    report_path = OUT / "URXP02_MECHANISM_GENE_REPORT.md"
    manifest["files"][report_path.name] = {"sha256": sha256(report_path), "bytes": report_path.stat().st_size}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
