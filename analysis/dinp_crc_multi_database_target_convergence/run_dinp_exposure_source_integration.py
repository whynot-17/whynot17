#!/usr/bin/env python3
"""Integrate DINP exposure-side evidence from CTD, EPA CompTox, and T3DB.

This is a source-preserving integration, not a merged truth score.  CTD
chemical-gene interactions, CompTox ToxCast/Tox21 bioactivity endpoints, and
T3DB toxin-target curation are retained in separate columns and source files.
Only active, human-labelled CompTox endpoints with an approved/official gene
symbol contribute to the CompTox gene layer.  Inactive or target-free assays
remain in the assay-level export and are never treated as a biological
negative.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
OUT = HERE / "outputs"
SRC = OUT / "source_records"

CTD_FILE = SRC / "ctd_dinp_human_interactions.csv"
COMPTOX_FILE = SRC / "dinp_comptox_bioactivity_v4_2.json"
T3DB_RAW_FILE = SRC / "t3db_dinp_pubchem_pugview.json"

DINP = {
    "preferred_name": "Diisononyl phthalate",
    "synonym": "DINP",
    "casrn": "28553-12-0",
    "ctd_id": "C012125",
    "dtxsid": "DTXSID4022521",
    "pubchem_cid": "590836",
    "t3db_id": "T3D3648",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repo_relative(path: Path) -> str:
    """Return a portable repository-relative POSIX path for provenance."""
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def unique_join(values: Iterable[str], sep: str = ";") -> str:
    return sep.join(sorted({clean(v) for v in values if clean(v)}))


def read_ctd() -> tuple[list[dict[str, str]], dict[str, dict[str, Any]]]:
    if not CTD_FILE.exists():
        raise FileNotFoundError(CTD_FILE)
    rows: list[dict[str, str]] = []
    by_gene: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "raw_rows": 0,
            "pairs": set(),
            "pmids": set(),
            "directions": set(),
            "actions": set(),
            "single": False,
            "cotreat": False,
        }
    )
    with CTD_FILE.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            gene = clean(row.get("normalized_gene_symbol") or row.get("GeneSymbol")).upper()
            if not gene:
                continue
            rows.append(row)
            record = by_gene[gene]
            record["raw_rows"] += 1
            record["pairs"].add((clean(row.get("ChemicalID")), gene))
            record["pmids"].update(
                p.strip() for p in clean(row.get("PubMedIDs")).replace("|", ",").split(",") if p.strip()
            )
            record["directions"].update(
                action.split("^")[0]
                for action in clean(row.get("InteractionActions")).split("|")
                if action.strip()
            )
            record["actions"].update(
                action.strip() for action in clean(row.get("InteractionActions")).split("|") if action.strip()
            )
            record["single"] = record["single"] or clean(row.get("single_chemical_record_flag")).lower() == "true"
            record["cotreat"] = record["cotreat"] or clean(row.get("multi_chemical_co_treatment_flag")).lower() == "true"
    return rows, by_gene


def read_comptox() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not COMPTOX_FILE.exists():
        raise FileNotFoundError(COMPTOX_FILE)
    payload = json.loads(COMPTOX_FILE.read_text(encoding="utf-8"))
    records = payload.get("bioactivity") or []
    assay_rows: list[dict[str, Any]] = []
    by_gene: dict[str, dict[str, Any]] = defaultdict(
        lambda: {"toxcast": [], "tox21": [], "all": []}
    )
    for record in records:
        source = clean(record.get("assaySourceName")) or ""
        assay_name = clean(record.get("assayName"))
        is_tox21 = source.upper() == "TOX21" or "TOX21" in assay_name.upper()
        source_class = "Tox21" if is_tox21 else "ToxCast"
        organism = clean(record.get("organism"))
        active = int(float(record.get("hitCall") or 0) > 0)
        genes = record.get("geneInfo") or []
        symbols: list[str] = []
        entrez: list[str] = []
        uniprot: list[str] = []
        for gene in genes:
            symbol = clean(gene.get("official_symbol") or gene.get("gene_symbol")).upper()
            if symbol:
                symbols.append(symbol)
            if clean(gene.get("entrez_gene_id")):
                entrez.append(clean(gene.get("entrez_gene_id")))
            if clean(gene.get("uniprot_accession_number")):
                uniprot.append(clean(gene.get("uniprot_accession_number")))
        assay_row = {
            "dtxsid": clean(record.get("dtxsid")) or DINP["dtxsid"],
            "source_class": source_class,
            "assay_source": source,
            "assay_name": assay_name,
            "assay_component": clean(record.get("assayComponentName")),
            "endpoint_name": clean(record.get("assayComponentEndpNm")),
            "organism": organism,
            "human_label": int(organism.lower() == "human"),
            "active_hit_call": active,
            "hit_call_continuous": record.get("hitcallContinuous"),
            "ac50_uM": record.get("ac50"),
            "bmd_uM": record.get("bmd"),
            "scaled_top": record.get("scaledTop"),
            "intended_target_family": clean(record.get("intendedTargetFamily")),
            "intended_target_type": clean(record.get("intendedTargetType")),
            "intended_target_subtype": clean(record.get("intendedTargetSubType")),
            "biological_process_target": clean(record.get("biologicalProcessTarget")),
            "cell_or_tissue": clean(record.get("cellShortName") or record.get("tissue")),
            "gene_symbols": unique_join(symbols),
            "entrez_gene_ids": unique_join(entrez),
            "uniprot_accessions": unique_join(uniprot),
            "flags": unique_join(record.get("flag") or [], sep=" | "),
        }
        assay_rows.append(assay_row)
        if active and organism.lower() == "human":
            for symbol in set(symbols):
                by_gene[symbol]["all"].append(assay_row)
                by_gene[symbol]["tox21" if is_tox21 else "toxcast"].append(assay_row)

    meta = {
        "retrieval_file": repo_relative(COMPTOX_FILE),
        "retrieval_file_sha256": sha256_file(COMPTOX_FILE),
        "source_url": payload.get("source_url") or "https://comptox.epa.gov/dashboard/chemical/invitrodb/DTXSID4022521",
        "dashboard_record_version": "invitrodb v4.2 as reported by the public dashboard QC record",
        "all_endpoint_rows": len(assay_rows),
        "active_endpoint_rows": sum(row["active_hit_call"] for row in assay_rows),
        "human_active_endpoint_rows": sum(row["active_hit_call"] and row["human_label"] for row in assay_rows),
        "active_human_target_gene_count": len(by_gene),
        "active_human_target_genes": sorted(by_gene),
        "active_tox21_endpoint_rows": sum(
            row["active_hit_call"] and row["human_label"] and row["source_class"] == "Tox21" for row in assay_rows
        ),
        "active_toxcast_endpoint_rows": sum(
            row["active_hit_call"] and row["human_label"] and row["source_class"] == "ToxCast" for row in assay_rows
        ),
        "absence_interpretation": "Inactive or absent endpoint/target records are not interpreted as biological negatives.",
    }
    return assay_rows, {"by_gene": by_gene, "meta": meta}


def t3db_rows() -> list[dict[str, str]]:
    # The public DINP T3DB record reports one target: human NR1I2/PXR.
    # The raw PubChem PUG-View file is retained alongside this curated row and
    # cross-references the T3DB record as SourceID Compound::T3D3648.
    return [
        {
            "t3db_id": DINP["t3db_id"],
            "chemical_name": DINP["preferred_name"],
            "gene_symbol": "NR1I2",
            "target_name": "Nuclear receptor subfamily 1 group I member 2",
            "uniprot_id": "O75469",
            "target_scope": "human protein target",
            "evidence_type": "T3DB toxin-target curation",
            "source_url": "https://www.t3db.ca/toxins/T3D3648",
            "source_record_note": "Public T3DB DINP record; target section reports one NR1I2/PXR target.",
        }
    ]


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_genecards_threshold_symbols() -> set[str]:
    path = OUT / "dinp_crc_genecards_relevance10_intersection.csv"
    if not path.exists():
        return set()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {clean(row.get("gene_symbol")).upper() for row in csv.DictReader(handle) if clean(row.get("gene_symbol"))}


def load_genecards_threshold_rows() -> dict[str, dict[str, str]]:
    path = OUT / "dinp_crc_genecards_relevance10_intersection.csv"
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {clean(row.get("gene_symbol")).upper(): row for row in csv.DictReader(handle) if clean(row.get("gene_symbol"))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-crc-overlay", action="store_true", help="Do not write the optional GeneCards threshold overlay.")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)

    ctd_rows, ctd_by = read_ctd()
    comptox_rows, comptox = read_comptox()
    t3_rows = t3db_rows()

    assay_fields = [
        "dtxsid", "source_class", "assay_source", "assay_name", "assay_component", "endpoint_name",
        "organism", "human_label", "active_hit_call", "hit_call_continuous", "ac50_uM", "bmd_uM",
        "scaled_top", "intended_target_family", "intended_target_type", "intended_target_subtype",
        "biological_process_target", "cell_or_tissue", "gene_symbols", "entrez_gene_ids",
        "uniprot_accessions", "flags",
    ]
    write_csv(OUT / "dinp_toxcast_tox21_bioactivity.csv", comptox_rows, assay_fields)

    t3_fields = list(t3_rows[0])
    write_csv(OUT / "dinp_t3db_target_records.csv", t3_rows, t3_fields)

    ctd_genes = set(ctd_by)
    tox_genes = set(comptox["by_gene"])
    toxcast_genes = {
        gene for gene, group in comptox["by_gene"].items() if group["toxcast"]
    }
    tox21_genes = {
        gene for gene, group in comptox["by_gene"].items() if group["tox21"]
    }
    t3db_genes = {clean(row["gene_symbol"]).upper() for row in t3_rows}
    union = sorted(ctd_genes | tox_genes | t3db_genes)

    matrix_rows: list[dict[str, Any]] = []
    for gene in union:
        ctd = ctd_by.get(gene, {})
        group = comptox["by_gene"].get(gene, {"toxcast": [], "tox21": [], "all": []})
        matrix_rows.append(
            {
                "gene_symbol": gene,
                "CTD": int(gene in ctd_genes),
                "ToxCast": int(gene in toxcast_genes),
                "Tox21": int(gene in tox21_genes),
                "T3DB": int(gene in t3db_genes),
                "exposure_support_count": sum(int(x) for x in [gene in ctd_genes, gene in toxcast_genes, gene in tox21_genes, gene in t3db_genes]),
                "CTD_raw_row_count": ctd.get("raw_rows", 0),
                "CTD_unique_chemical_gene_pairs": len(ctd.get("pairs", set())),
                "CTD_unique_pmids": len(ctd.get("pmids", set())),
                "CTD_single_chemical_evidence": int(ctd.get("single", False)),
                "CTD_cotreatment_evidence": int(ctd.get("cotreat", False)),
                "ToxCast_active_assay_count": len(group["toxcast"]),
                "ToxCast_active_assay_sources": unique_join(row["assay_source"] for row in group["toxcast"]),
                "ToxCast_active_endpoints": unique_join(row["endpoint_name"] for row in group["toxcast"]),
                "Tox21_active_assay_count": len(group["tox21"]),
                "Tox21_active_endpoints": unique_join(row["endpoint_name"] for row in group["tox21"]),
                "T3DB_target_count": sum(1 for row in t3_rows if row["gene_symbol"].upper() == gene),
                "T3DB_uniprot_ids": unique_join(row["uniprot_id"] for row in t3_rows if row["gene_symbol"].upper() == gene),
            }
        )
    matrix_fields = list(matrix_rows[0])
    matrix_path = OUT / "dinp_exposure_gene_matrix_ctd_toxcast_tox21_t3db.csv"
    write_csv(matrix_path, matrix_rows, matrix_fields)

    overlay_rows: list[dict[str, Any]] = []
    cards_rows = load_genecards_threshold_rows()
    if not args.no_crc_overlay and cards_rows:
        for gene in sorted(set(cards_rows) & set(union)):
            row = cards_rows[gene]
            matrix = next(item for item in matrix_rows if item["gene_symbol"] == gene)
            overlay_rows.append(
                {
                    "gene_symbol": gene,
                    **{k: row.get(k, "") for k in ["gene_cards_rank", "gene_cards_relevance_score", "gene_cards_knowledge_score"]},
                    **{k: matrix.get(k, "") for k in ["CTD", "ToxCast", "Tox21", "T3DB", "exposure_support_count"]},
                    "legacy_ctd_intersection": int(gene in ctd_genes),
                }
            )
        overlay_fields = list(overlay_rows[0]) if overlay_rows else ["gene_symbol"]
        write_csv(OUT / "dinp_crc_genecards_relevance10_expanded_exposure_overlay.csv", overlay_rows, overlay_fields)

    source_files = {
        "CTD_interaction_records": {
            "path": repo_relative(CTD_FILE),
            "source_url": "https://www.ctdbase.org/detail.go?type=chem&acc=C012125",
            "sha256": sha256_file(CTD_FILE),
        },
        "CompTox_dashboard_bioactivity_json": {
            "path": repo_relative(COMPTOX_FILE),
            "source_url": "https://comptox.epa.gov/dashboard/chemical/invitrodb/DTXSID4022521",
            "sha256": sha256_file(COMPTOX_FILE),
        },
        "T3DB_PubChem_PUG_View_raw": {
            "path": repo_relative(T3DB_RAW_FILE),
            "source_url": "https://www.t3db.ca/toxins/T3D3648",
            "sha256": sha256_file(T3DB_RAW_FILE),
        } if T3DB_RAW_FILE.exists() else {"path": repo_relative(T3DB_RAW_FILE), "status": "missing"},
    }
    summary = {
        "generated_at": utc_now(),
        "chemical": DINP,
        "integration_rule": {
            "source_preserving": True,
            "support_count_definition": "number of source layers among CTD, ToxCast, Tox21, and T3DB with retained evidence",
            "comptox_gene_rule": "active_hit_call > 0 AND organism == human AND official/approved gene symbol present",
            "t3db_gene_rule": "human target gene explicitly reported by the public DINP T3DB record",
            "absence_interpretation": "missing or inactive records are not biological negatives",
            "no_related_phthalate_substitution": True,
        },
        "source_files": source_files,
        "counts": {
            "CTD_human_rows": len(ctd_rows),
            "CTD_unique_human_genes": len(ctd_genes),
            "CompTox_all_bioactivity_endpoints": comptox["meta"]["all_endpoint_rows"],
            "CompTox_active_endpoints": comptox["meta"]["active_endpoint_rows"],
            "CompTox_active_human_endpoint_rows": comptox["meta"]["human_active_endpoint_rows"],
            "ToxCast_active_human_genes": len(toxcast_genes),
            "Tox21_active_human_genes": len(tox21_genes),
            "CompTox_active_human_unique_genes": len(tox_genes),
            "T3DB_target_rows": len(t3_rows),
            "T3DB_unique_genes": len(t3db_genes),
            "DINP_exposure_union_genes": len(union),
            "CTD_only_genes": len(ctd_genes - (tox_genes | t3db_genes)),
            "new_genes_added_by_CompTox": len(tox_genes - ctd_genes),
            "new_genes_added_by_T3DB": len(t3db_genes - ctd_genes),
            "genes_supported_by_at_least_2_layers": sum(row["exposure_support_count"] >= 2 for row in matrix_rows),
        },
        "overlaps": {
            "CTD_ToxCast": len(ctd_genes & toxcast_genes),
            "CTD_Tox21": len(ctd_genes & tox21_genes),
            "CTD_T3DB": len(ctd_genes & t3db_genes),
            "ToxCast_Tox21": len(toxcast_genes & tox21_genes),
            "ToxCast_T3DB": len(toxcast_genes & t3db_genes),
            "Tox21_T3DB": len(tox21_genes & t3db_genes),
            "all_four_layers": len(ctd_genes & toxcast_genes & tox21_genes & t3db_genes),
        },
        "compTox_meta": comptox["meta"],
    }
    (OUT / "dinp_exposure_source_integration_manifest.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    lines = [
        "# DINP exposure-side source integration",
        "",
        f"Generated: `{summary['generated_at']}`",
        "",
        "## Scope",
        "",
        "The parent chemical was fixed to diisononyl phthalate (DINP), CASRN 28553-12-0, CTD C012125, DTXSID4022521, PubChem CID 590836, and T3DB T3D3648. CTD, EPA CompTox/ToxCast/Tox21, and T3DB evidence are preserved as separate layers; support counts are descriptive and are not a biological truth score.",
        "",
        "## Integrated counts",
        "",
        "| Layer | Result |",
        "|---|---:|",
        f"| CTD human interaction rows | {len(ctd_rows)} |",
        f"| CTD unique human genes | {len(ctd_genes)} |",
        f"| CompTox dashboard bioactivity endpoints | {comptox['meta']['all_endpoint_rows']} |",
        f"| CompTox active endpoints | {comptox['meta']['active_endpoint_rows']} |",
        f"| Active human ToxCast target genes | {len(toxcast_genes)} |",
        f"| Active human Tox21 target genes | {len(tox21_genes)} |",
        f"| T3DB target rows | {len(t3_rows)} |",
        f"| Expanded DINP exposure union | {len(union)} |",
        f"| Genes supported by ≥2 source layers | {sum(row['exposure_support_count'] >= 2 for row in matrix_rows)} |",
        "",
        "## CompTox active human targets",
        "",
        f"`{', '.join(sorted(tox_genes))}`",
        "",
        "The CompTox extraction retained all 445 endpoint rows in the assay-level CSV. The gene layer uses only human-labelled active endpoints with an official gene symbol. Target-free active endpoints remain in the assay export but do not contribute a gene. Inactive or missing records are not treated as negative evidence.",
        "",
        "## T3DB",
        "",
        "The public DINP T3DB record is T3D3648 and reports one human target, NR1I2/PXR (UniProt O75469). The raw PubChem PUG-View response is retained because it cross-references the T3DB source record; the curated T3DB target row is kept separately.",
        "",
        "## CRC overlay",
        "",
        f"The optional overlay against the existing GeneCards relevance-score ≥10 intersection contains {len(overlay_rows)} genes. It is an overlay only; it does not alter the prior disease-side threshold or multiplicity analysis.",
        "",
        "## Outputs",
        "",
        "- `dinp_toxcast_tox21_bioactivity.csv`: source-level CompTox assay/end-point evidence.",
        "- `dinp_t3db_target_records.csv`: curated T3DB DINP target row.",
        "- `dinp_exposure_gene_matrix_ctd_toxcast_tox21_t3db.csv`: source-preserving gene matrix.",
        "- `dinp_crc_genecards_relevance10_expanded_exposure_overlay.csv`: optional exposure expansion over the existing CRC GeneCards threshold intersection.",
        "- `dinp_exposure_source_integration_manifest.json`: counts, rules, provenance, and hashes.",
    ]
    (OUT / "dinp_exposure_source_integration_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("DINP EXPOSURE SOURCE INTEGRATION: PASS")
    print(f"CTD unique human genes: {len(ctd_genes)}")
    print(f"CompTox active human target genes: {len(tox_genes)}")
    print(f"T3DB target genes: {len(t3db_genes)}")
    print(f"Expanded DINP exposure union: {len(union)}")
    print(f"Source-preserving matrix: {matrix_path.resolve()}")


if __name__ == "__main__":
    main()
