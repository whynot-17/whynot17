"""M1b: expanded molecular universe for 2-naphthol (URXP02).

The script broadens the evidence tiers without making causal or sex-specific
claims.  It keeps exact 2-naphthol separate from parent naphthalene and keeps
assay-level records long until the explicit gene summary step.
"""
from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
CTD_CHEM = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "CTD_chemicals.tsv.gz"
CTD_IXN = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data" / "CTD_chem_gene_ixns.tsv.gz"
OT_CHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={ids}&retmode=json"
CHEMBL_ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={}&limit=1000"
CHEMBL_TARGET = "https://www.ebi.ac.uk/chembl/api/data/target/{}.json"

EXACT = {"chemical": "2-naphthol", "chemical_identity": "CTD:C028405; PubChem CID:8663; ChEMBL:CHEMBL14126", "ctd_id": "C028405", "pubchem_cid": 8663, "chembl_id": "CHEMBL14126"}
PARENT = {"chemical": "naphthalene", "chemical_identity": "CTD:C031721; PubChem CID:931; ChEMBL:CHEMBL16293", "ctd_id": "C031721", "pubchem_cid": 931, "chembl_id": "CHEMBL16293"}

LONG_FIELDS = ["gene_symbol", "gene_id", "chemical", "chemical_identity", "evidence_tier", "evidence_type", "species", "tissue_or_cell", "direction", "effect_size", "p_value", "fdr", "dose", "duration", "dataset_or_reference", "source_database", "notes"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def ctid(value: object) -> str:
    s = str(value or "").strip()
    return s[5:] if s.startswith("MESH:") else s


def split_pmids(value: object) -> list[str]:
    return [x.strip() for x in re.split(r"[|;,]", str(value or "")) if x.strip()]


def ctd_rows(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = None
        for line in handle:
            if line.startswith("# Fields:"):
                header = next(handle).lstrip("# ").rstrip("\r\n").split("\t")
                break
        if header is None:
            raise RuntimeError(f"Missing CTD header in {path}")
        yield from csv.DictReader(handle, fieldnames=header, delimiter="\t")


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def ctd_action_direction(actions: str) -> str:
    vals = []
    for item in str(actions or "").split("|"):
        item = item.strip()
        if not item:
            continue
        vals.append(item.split("^", 1)[0])
    return ";".join(sorted(set(vals)))


def exact_ctd_rows(chemical: dict) -> list[dict]:
    out = []
    for row in ctd_rows(CTD_IXN):
        if ctid(row.get("ChemicalID")) != chemical["ctd_id"]:
            continue
        symbol = str(row.get("GeneSymbol") or "").strip()
        gid = str(row.get("GeneID") or "").strip()
        if not symbol or not gid:
            continue
        species = str(row.get("Organism") or "Unknown").strip() or "Unknown"
        tier = "Tier A" if species == "Homo sapiens" else "Tier B"
        pmids = split_pmids(row.get("PubMedIDs"))
        out.append({
            "gene_symbol": symbol,
            "gene_id": gid,
            "chemical": chemical["chemical"],
            "chemical_identity": chemical["chemical_identity"],
            "evidence_tier": tier,
            "evidence_type": "curated chemical-gene interaction",
            "species": species,
            "tissue_or_cell": "",
            "direction": ctd_action_direction(row.get("InteractionActions")),
            "effect_size": "",
            "p_value": "",
            "fdr": "",
            "dose": "",
            "duration": "",
            "dataset_or_reference": ";".join("CTD-PMID:" + p for p in pmids),
            "source_database": "CTD",
            "notes": str(row.get("Interaction") or ""),
        })
    return out


def post_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "URXP02-M1b/1.0"})
    with urllib.request.urlopen(req, timeout=180) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_pubchem(chemical: dict) -> list[dict]:
    payload = post_json(OT_CHEM.format(cid=chemical["pubchem_cid"]))
    table = payload.get("Table", {})
    cols = table.get("Columns", {}).get("Column", [])
    out = []
    for row in table.get("Row", []):
        cells = row.get("Cell", [])
        rec = dict(zip(cols, cells))
        gid = str(rec.get("Target GeneID") or "").strip()
        if not gid:
            continue
        out.append({"chemical": chemical, "gene_id": gid, "activity": rec})
    return out


def ncbi_gene_map(ids: set[str]) -> dict[str, dict]:
    result = {}
    ids = sorted(x for x in ids if x)
    for start in range(0, len(ids), 100):
        chunk = ids[start:start + 100]
        payload = post_json(NCBI_ESUMMARY.format(ids=urllib.parse.quote(",".join(chunk))))
        for uid in payload.get("result", {}).get("uids", []):
            rec = payload["result"].get(str(uid), {})
            result[str(uid)] = rec
        time.sleep(0.05)
    return result


def pubchem_rows(records: list[dict], gene_map: dict[str, dict], tier_parent: bool = False) -> list[dict]:
    out = []
    for record in records:
        chemical = record["chemical"]
        rec = record["activity"]
        gid = record["gene_id"]
        gene = gene_map.get(gid, {})
        species = ((gene.get("organism") or {}).get("scientificname") or "Unknown")
        symbol = str(gene.get("name") or "").strip() or "UNKNOWN_GENE_" + gid
        tier = "Tier C" if tier_parent else ("Tier A" if species == "Homo sapiens" else "Tier B")
        activity_value = str(rec.get("Activity Value [uM]") or "").strip()
        effect = str(rec.get("Activity Outcome") or "").strip()
        notes = "assay_name=" + str(rec.get("Assay Name") or "") + "; assay_type=" + str(rec.get("Assay Type") or "") + "; activity_name=" + str(rec.get("Activity Name") or "") + "; target_accession=" + str(rec.get("Target Accession") or "")
        ref = "PubChem-AID:" + str(rec.get("AID") or "")
        if rec.get("PubMed ID"):
            ref += "; PMID:" + str(rec.get("PubMed ID"))
        out.append({
            "gene_symbol": symbol,
            "gene_id": gid,
            "chemical": chemical["chemical"],
            "chemical_identity": chemical["chemical_identity"],
            "evidence_tier": tier,
            "evidence_type": "bioassay activity",
            "species": species,
            "tissue_or_cell": "",
            "direction": effect,
            "effect_size": activity_value,
            "p_value": "",
            "fdr": "",
            "dose": "",
            "duration": "",
            "dataset_or_reference": ref,
            "source_database": "PubChem BioAssay",
            "notes": notes,
        })
    return out


def chembl_target_details(target_ids: set[str]) -> dict[str, dict]:
    details = {}
    valid = sorted(x for x in target_ids if x)
    def fetch(target_id: str) -> tuple[str, dict]:
        try:
            return target_id, post_json(CHEMBL_TARGET.format(target_id))
        except Exception:
            return target_id, {}
    with ThreadPoolExecutor(max_workers=12) as pool:
        futures = [pool.submit(fetch, target_id) for target_id in valid]
        for future in as_completed(futures):
            target_id, detail = future.result()
            details[target_id] = detail
    return details


def chembl_records(chemical: dict) -> list[dict]:
    payload = post_json(CHEMBL_ACTIVITY.format(chemical["chembl_id"]))
    activities = payload.get("activities", [])
    target_ids = {str(a.get("target_chembl_id") or "") for a in activities}
    targets = chembl_target_details(target_ids)
    out = []
    for a in activities:
        target_id = str(a.get("target_chembl_id") or "")
        target = targets.get(target_id, {})
        organism = str(target.get("organism") or "Unknown")
        components = target.get("target_components") or []
        if not components:
            components = [{"component_description": target.get("pref_name") or "", "target_component_synonyms": []}]
        for component in components:
            symbols = [str(s.get("component_synonym") or "").strip() for s in (component.get("target_component_synonyms") or []) if s.get("syn_type") in {"GENE_SYMBOL", "GENE_SYMBOL_OTHER"}]
            symbol = (symbols[0] if symbols else str(component.get("component_description") or "").strip()) or "UNKNOWN_TARGET_" + target_id
            tier = "Tier C" if chemical["chemical"] == "naphthalene" else ("Tier A" if organism == "Homo sapiens" else "Tier B")
            ref = "ChEMBL-assay:" + str(a.get("assay_chembl_id") or "")
            if a.get("document_chembl_id"):
                ref += "; ChEMBL-document:" + str(a.get("document_chembl_id"))
            out.append({
                "gene_symbol": symbol,
                "gene_id": str(component.get("accession") or component.get("component_id") or ""),
                "chemical": chemical["chemical"],
                "chemical_identity": chemical["chemical_identity"],
                "evidence_tier": tier,
                "evidence_type": "target activity / physicochemical assay",
                "species": organism,
                "tissue_or_cell": "",
                # ChEMBL does not provide a common signed direction for these
                # heterogeneous activity/physicochemical records; retain the
                # measurement type and units in notes instead of mislabelling
                # IC50/LogP/etc. as a biological direction.
                "direction": "",
                "effect_size": str(a.get("standard_value") or ""),
                "p_value": "",
                "fdr": "",
                "dose": "",
                "duration": "",
                "dataset_or_reference": ref,
                "source_database": "ChEMBL",
                "notes": "target=" + target_id + "; target_type=" + str(target.get("target_type") or "") + "; standard_type=" + str(a.get("standard_type") or "") + "; standard_units=" + str(a.get("standard_units") or "") + "; relation=" + str(a.get("standard_relation") or "") + "; assay_description=" + str(a.get("assay_description") or ""),
            })
    return out


def resource_audit() -> dict:
    audit = {}
    for db, base in [("GEO", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term="), ("BioStudies", "https://www.ebi.ac.uk/biostudies/api/v1/search?query=")]:
        counts = {}
        for term in ["2-naphthol", "2-hydroxynaphthalene", "2-NAP"]:
            try:
                if db == "GEO":
                    payload = post_json(base + urllib.parse.quote(term) + "&retmode=json&retmax=0")
                    counts[term] = int(payload.get("esearchresult", {}).get("count", 0))
                else:
                    payload = post_json(base + urllib.parse.quote(term) + "&pageSize=1")
                    counts[term] = int(payload.get("totalHits", payload.get("hits", 0) if isinstance(payload.get("hits", 0), int) else len(payload.get("hits", []))))
            except Exception as exc:
                counts[term] = "unavailable:" + type(exc).__name__
        audit[db] = counts
    audit["Tox21_ToxCast"] = "No separate bulk source ingested; PubChem BioAssay records retained as the assay-level proxy and source-labelled PubChem rows."
    audit["LINCS"] = "No exact 2-naphthol perturbation dataset ingested."
    return audit


def summarize(long_rows: list[dict]) -> list[dict]:
    by_gene = defaultdict(list)
    for row in long_rows:
        by_gene[row["gene_symbol"]].append(row)
    out = []
    for symbol in sorted(by_gene):
        rows = by_gene[symbol]
        sources = {r["source_database"] for r in rows}
        refs = {r["dataset_or_reference"] for r in rows if r["dataset_or_reference"]}
        species = {r["species"] for r in rows if r["species"] and r["species"] != "Unknown"}
        contexts = {r["tissue_or_cell"] for r in rows if r["tissue_or_cell"]}
        directions = {r["direction"].lower() for r in rows if r["direction"]}
        direction_class = "not_assessable"
        if directions:
            if len(directions) == 1:
                direction_class = "consistent:" + next(iter(directions))
            else:
                direction_class = "mixed"
        exact_human = any(r["chemical"] == "2-naphthol" and r["species"] == "Homo sapiens" for r in rows)
        exact_exp = any(r["chemical"] == "2-naphthol" and r["species"] != "Homo sapiens" for r in rows)
        parent = any(r["chemical"] == "naphthalene" for r in rows)
        ctd = any(r["source_database"] == "CTD" for r in rows)
        tox = any(r["source_database"] in {"PubChem BioAssay", "ChEMBL"} or "transcriptomic" in r["evidence_type"] for r in rows)
        bio = any(r["source_database"] == "PubChem BioAssay" for r in rows)
        target = any(r["source_database"] == "ChEMBL" or "binding" in r["direction"].lower() or "activity" in r["evidence_type"] for r in rows)
        metabolism = any(any(x in (r["evidence_type"] + " " + r["direction"] + " " + r["notes"]).lower() for x in ["metabol", "sulfat", "glucuron", "abundance"]) for r in rows)
        binding = any(any(x in (r["evidence_type"] + " " + r["direction"] + " " + r["notes"]).lower() for x in ["binding", "activity", "ic50", "inhibition"]) for r in rows)
        transcript = any("transcriptomic" in r["evidence_type"].lower() for r in rows)
        gene_ids = [r["gene_id"] for r in rows if r["gene_id"]]
        out.append({
            "gene_symbol": symbol,
            "gene_id": gene_ids[0] if gene_ids else "",
            "exact_2NAP_human_support": exact_human,
            "exact_2NAP_experimental_support": exact_exp,
            "parent_naphthalene_support": parent,
            "number_of_sources": len(sources),
            "number_of_independent_references": len(refs),
            "number_of_species": len(species),
            "number_of_tissues_or_cells": len(contexts),
            "perturbation_direction_consistency": direction_class,
            "CTD_support": ctd,
            "toxicogenomic_support": tox,
            "bioassay_target_support": bio,
            "direct_target_or_binding_support": target,
            "exposure_metabolism_flag": metabolism,
            "direct_target_binding_flag": binding,
            "transcriptomic_response_flag": transcript,
            "evidence_rows": len(rows),
        })
    return out


def audit_rows(long_rows: list[dict], summary: list[dict], resource: dict) -> list[dict]:
    genes = {r["gene_symbol"] for r in long_rows}
    exact_genes = {r["gene_symbol"] for r in long_rows if r["chemical"] == "2-naphthol"}
    parent_only = {r["gene_symbol"] for r in long_rows if r["chemical"] == "naphthalene"} - exact_genes
    multi_source = {r["gene_symbol"] for r in summary if int(r["number_of_sources"]) >= 2}
    human = {r["gene_symbol"] for r in long_rows if r["species"] == "Homo sapiens"}
    perturb = {r["gene_symbol"] for r in long_rows if r["source_database"] in {"PubChem BioAssay", "ChEMBL"} or "transcriptomic" in r["evidence_type"]}
    refs = {r["dataset_or_reference"] for r in long_rows if r["dataset_or_reference"]}
    species = {r["species"] for r in long_rows if r["species"] and r["species"] != "Unknown"}
    contexts = {r["tissue_or_cell"] for r in long_rows if r["tissue_or_cell"]}
    metabolism = {r["gene_symbol"] for r in long_rows if any(x in (r["evidence_type"] + " " + r["direction"] + " " + r["notes"]).lower() for x in ["metabol", "sulfat", "glucuron", "abundance"])}
    binding = {r["gene_symbol"] for r in long_rows if any(x in (r["evidence_type"] + " " + r["direction"] + " " + r["notes"]).lower() for x in ["binding", "activity", "ic50", "inhibition"])}
    transcript = {r["gene_symbol"] for r in long_rows if "transcriptomic" in r["evidence_type"].lower()}
    metrics = [
        ("total_unique_genes", len(genes), "Unique gene_symbol values across all retained tiers"),
        ("genes_with_exact_2NAP_evidence", len(exact_genes), "Genes with exact 2-naphthol CTD/PubChem/ChEMBL evidence"),
        ("genes_with_parent_naphthalene_only_support", len(parent_only), "Genes present for naphthalene but absent from exact 2-naphthol evidence"),
        ("genes_supported_by_at_least_2_sources", len(multi_source), "Gene symbols with >=2 source_database values"),
        ("genes_supported_in_human_data", len(human), "Gene symbols with at least one Homo sapiens row"),
        ("genes_supported_by_perturbation_response", len(perturb), "Gene symbols with PubChem BioAssay/ChEMBL or transcriptomic evidence"),
        ("number_of_datasets_or_references", len(refs), "Unique dataset_or_reference values"),
        ("number_of_species", len(species), "Distinct species values excluding Unknown"),
        ("number_of_tissue_or_cell_contexts", len(contexts), "Distinct non-empty tissue_or_cell values"),
        ("exposure_metabolism_flagged_genes", len(metabolism), "Rows mentioning metabolism, sulfation, glucuronidation, or abundance"),
        ("direct_target_binding_flagged_genes", len(binding), "Rows mentioning binding, activity, IC50, or inhibition"),
        ("transcriptomic_response_flagged_genes", len(transcript), "Rows explicitly labelled transcriptomic response"),
        ("exact_resource_audit", json.dumps(resource, ensure_ascii=False), "Search audit for exact terms; not treated as evidence rows"),
    ]
    return [{"metric": m, "value": v, "definition": d} for m, v, d in metrics]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    exact_ctd = exact_ctd_rows(EXACT)
    parent_ctd = exact_ctd_rows(PARENT)
    exact_pub = fetch_pubchem(EXACT)
    parent_pub = fetch_pubchem(PARENT)
    gene_ids = {r["gene_id"] for r in exact_pub + parent_pub}
    gene_map = ncbi_gene_map(gene_ids)
    exact_pub_rows = pubchem_rows(exact_pub, gene_map)
    parent_pub_rows = pubchem_rows(parent_pub, gene_map, tier_parent=True)
    exact_chembl = chembl_records(EXACT)
    parent_chembl = chembl_records(PARENT)
    # CTD curated interaction rows are retained in file 01.  File 02 is the
    # distinct assay/target perturbation layer so the unified table does not
    # double-count the same CTD evidence record.
    exact_perturb = exact_pub_rows + exact_chembl
    parent_long = parent_ctd + parent_pub_rows + parent_chembl
    unified = exact_ctd + exact_perturb + parent_long
    # Retain the exact CTD table as a curated-evidence table; 02 is intentionally
    # assay/perturbation-long and may contain repeated CTD records by design.
    write_csv(OUT / "01_exact_2NAP_evidence.csv", exact_ctd, LONG_FIELDS)
    write_csv(OUT / "02_2NAP_perturbation_evidence.csv", exact_perturb, LONG_FIELDS)
    write_csv(OUT / "03_parent_naphthalene_evidence.csv", parent_long, LONG_FIELDS)
    write_csv(OUT / "04_unified_gene_evidence_long.csv", unified, LONG_FIELDS)
    summary = summarize(unified)
    summary_fields = list(summary[0]) if summary else ["gene_symbol"]
    write_csv(OUT / "05_gene_evidence_summary.csv", summary, summary_fields)
    resource = resource_audit()
    audit = audit_rows(unified, summary, resource)
    write_csv(OUT / "06_molecular_universe_audit.csv", audit, ["metric", "value", "definition"])
    exact_genes = sorted({r["gene_symbol"] for r in unified if r["chemical"] == "2-naphthol"})
    parent_genes = sorted({r["gene_symbol"] for r in unified if r["chemical"] == "naphthalene"})
    report = f"""# URXP02 / 2-NAP expanded molecular universe (M1b)\n\nGenerated {datetime.now(timezone.utc).isoformat()}. This package maximizes evidence coverage while keeping exact 2-naphthol separate from parent naphthalene. It does not re-prove the NHANES association, infer causality, or make sex-specific molecular claims.\n\n## Evidence tiers\n\n- **Tier A:** exact 2-naphthol human evidence.\n- **Tier B:** exact 2-naphthol experimental/non-human evidence.\n- **Tier C:** parent naphthalene support; never labelled exact 2-NAP.\n\nThe exact CTD identity is `C028405` (2-naphthol; synonyms 2-hydroxynaphthalene/2-NAP). Parent naphthalene is `C031721`. PubChem and ChEMBL identifiers are recorded in every long-table row.\n\n## Size audit\n\n- Exact 2-NAP genes: **{len(exact_genes)}**\n- Parent naphthalene genes: **{len(parent_genes)}**\n- Unified evidence rows: **{len(unified)}**\n- CTD exact rows: **{len(exact_ctd)}**; PubChem exact assay rows: **{len(exact_pub_rows)}**; ChEMBL exact activity rows: **{len(exact_chembl)}**\n- Parent CTD rows: **{len(parent_ctd)}**; parent PubChem assay rows: **{len(parent_pub_rows)}**; parent ChEMBL activity rows: **{len(parent_chembl)}**\n\nSee `06_molecular_universe_audit.csv` for the full metric table.\n\n## Source and coverage notes\n\nCTD records retain interaction actions, interaction sentences, species, and PubMed identifiers. PubChem BioAssay records retain activity outcome/value, assay identifiers, target GeneID, and NCBI-mapped species/symbol when available. ChEMBL records retain target organism, assay type, standard measurement, target identifier, and component accession/symbol. GEO/BioStudies exact-term search results and unintegrated resource status are recorded in the audit/manifest; no exact expression dataset was silently converted into a gene list.\n\nThe current expanded universe is therefore suitable for downstream pathway/network/tissue/cell work, but the small exact-2-NAP human evidence base and the assay-heavy Tier C records should be kept visibly separate in later analyses.\n\n## Restrictions honored\n\nNo disease intersection, pathway enrichment, tissue/cell analysis, figures, or new NHANES association model was run.\n"""
    report_path = OUT / "URXP02_EXPANDED_MOLECULAR_UNIVERSE_REPORT.md"
    report_path.write_text(report, encoding="utf-8")
    files = sorted(list(OUT.glob("*.csv")) + [report_path])
    manifest = {
        "analysis": "URXP02 / 2-NAP expanded molecular universe (M1b)",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "chemical_identities": {"exact_2NAP": EXACT, "parent_naphthalene": PARENT},
        "evidence_tiers": {"Tier A": "exact 2-naphthol human", "Tier B": "exact 2-naphthol experimental/non-human", "Tier C": "parent naphthalene; not exact 2-NAP"},
        "source_files": {"ctd_chem_gene_ixns": {"path": str(CTD_IXN), "sha256": sha256(CTD_IXN)}, "ctd_chemicals": {"path": str(CTD_CHEM), "sha256": sha256(CTD_CHEM)}},
        "external_endpoints": {"pubchem_assaysummary": OT_CHEM, "ncbi_gene_esummary": NCBI_ESUMMARY, "chembl_activity": CHEMBL_ACTIVITY, "chembl_target": CHEMBL_TARGET, "geo_search": "NCBI E-utilities gds esearch", "biostudies_search": "EMBL-EBI BioStudies API"},
        "resource_audit": resource,
        "row_counts": {"exact_ctd": len(exact_ctd), "exact_perturbation": len(exact_perturb), "parent": len(parent_long), "unified": len(unified), "unique_genes": len({r["gene_symbol"] for r in unified})},
        "constraints": ["No figures", "No disease intersection", "No pathway enrichment", "No tissue/cell analysis", "No new NHANES model", "Parent naphthalene is Tier C and is not exact 2-NAP", "No causal inference", "No sex-specific gene claim"],
        "files": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size} for p in files},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
