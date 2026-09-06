#!/usr/bin/env python3
"""Run a reproducible first-pass PubMed collision audit for Phase 2.

The input candidate universe is the Phase 1 mapping-gate output.  This script
does not infer that a PubMed hit is an eligible study: it records retrieval
counts and the exact E-utilities query so that title/abstract screening can be
performed transparently in the next pass.  It deliberately does not run
docking, target discovery, mechanism expansion, or experimental feasibility.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


SCRIPT_VERSION = "phase2_pubmed_collision_audit_v1.0"
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
SEARCH_DATE = datetime.now(timezone.utc).date().isoformat()

SEARCH_GROUPS = OrderedDict(
    [
        (
            "uranium",
            {
                "aliases": '"uranium"[Title/Abstract] OR "urinary uranium"[Title/Abstract] OR "uranium exposure"[Title/Abstract]',
                "chemical_ids": ["D014501"],
            },
        ),
        (
            "molybdenum",
            {
                "aliases": '"molybdenum"[Title/Abstract] OR "urinary molybdenum"[Title/Abstract] OR "molybdenum exposure"[Title/Abstract]',
                "chemical_ids": ["D008982"],
            },
        ),
        (
            "tungsten",
            {
                "aliases": '"tungsten"[Title/Abstract] OR "tungstate"[Title/Abstract] OR "tungsten exposure"[Title/Abstract]',
                "chemical_ids": ["D014414"],
            },
        ),
        (
            "lead",
            {
                "aliases": '"lead exposure"[Title/Abstract] OR "blood lead"[Title/Abstract] OR "urinary lead"[Title/Abstract] OR "lead"[MeSH Terms]',
                "chemical_ids": ["D007854"],
            },
        ),
        (
            "barium",
            {
                "aliases": '"barium"[Title/Abstract] OR "barium exposure"[Title/Abstract] OR "barium"[MeSH Terms]',
                "chemical_ids": ["D001464"],
            },
        ),
        (
            "tin",
            {
                "aliases": '"tin exposure"[Title/Abstract] OR "organotin"[Title/Abstract] OR "stannous"[Title/Abstract] OR "tin"[MeSH Terms]',
                "chemical_ids": ["D014001"],
            },
        ),
        (
            "silver",
            {
                "aliases": '"silver exposure"[Title/Abstract] OR "silver nanoparticles"[Title/Abstract] OR "silver"[MeSH Terms]',
                "chemical_ids": ["D012834"],
            },
        ),
        (
            "mibp",
            {
                "aliases": '"mono-isobutyl phthalate"[Title/Abstract] OR "MiBP"[Title/Abstract] OR "MIBP"[Title/Abstract]',
                "chemical_ids": ["C575690"],
            },
        ),
        (
            "dehp",
            {
                "aliases": '"di(2-ethylhexyl) phthalate"[Title/Abstract] OR "diethylhexyl phthalate"[Title/Abstract] OR "DEHP"[Title/Abstract]',
                "chemical_ids": ["D004051"],
            },
        ),
        (
            "dinp_parent",
            {
                "aliases": '"diisononyl phthalate"[Title/Abstract] OR "di-n-isononyl phthalate"[Title/Abstract] OR "DINP"[Title/Abstract] OR "dinonyl phthalate"[Title/Abstract]',
                "chemical_ids": ["C012125", "C019174"],
            },
        ),
        (
            "meohp",
            {
                "aliases": '"mono(2-ethyl-5-oxohexyl) phthalate"[Title/Abstract] OR "MEOHP"[Title/Abstract] OR "mono-2-ethyl-5-oxohexyl phthalate"[Title/Abstract]',
                "chemical_ids": ["C080276"],
            },
        ),
        (
            "mcop",
            {
                "aliases": '"mono(carboxy-isooctyl) phthalate"[Title/Abstract] OR "MCOP"[Title/Abstract] OR "mono-carboxy-isooctyl phthalate"[Title/Abstract]',
                "chemical_ids": ["C573544"],
            },
        ),
        (
            "mecpp",
            {
                "aliases": '"mono(2-ethyl-5-carboxypentyl) phthalate"[Title/Abstract] OR "MECPP"[Title/Abstract] OR "2-ethyl-5-carboxypentyl phthalate"[Title/Abstract]',
                "chemical_ids": ["C051450"],
            },
        ),
        (
            "pfhxs",
            {
                "aliases": '"perfluorohexanesulfonic acid"[Title/Abstract] OR "perfluorohexanesulfonate"[Title/Abstract] OR "PFHxS"[Title/Abstract]',
                "chemical_ids": ["C471071"],
            },
        ),
    ]
)

SEARCH_CATEGORIES = OrderedDict(
    [
        (
            "diabetes_total",
            '(diabet*[Title/Abstract] OR "insulin resistance"[Title/Abstract] OR hyperglyc*[Title/Abstract] OR "glucose metabolism"[Title/Abstract])',
        ),
        (
            "human_epidemiology",
            '(diabet*[Title/Abstract] OR "insulin resistance"[Title/Abstract]) AND (human*[Title/Abstract] OR cohort[Title/Abstract] OR population[Title/Abstract] OR epidemiol*[Title/Abstract] OR NHANES[Title/Abstract] OR cross-sectional[Title/Abstract] OR prospective[Title/Abstract])',
        ),
        (
            "prospective",
            '(diabet*[Title/Abstract] OR "insulin resistance"[Title/Abstract]) AND (prospective[Title/Abstract] OR longitudinal[Title/Abstract] OR incidence[Title/Abstract] OR incident[Title/Abstract] OR follow-up[Title/Abstract])',
        ),
        (
            "mechanism",
            '(diabet*[Title/Abstract] OR "insulin resistance"[Title/Abstract]) AND (mechanism*[Title/Abstract] OR beta-cell[Title/Abstract] OR "insulin secretion"[Title/Abstract] OR hepatocyte*[Title/Abstract] OR adipocyte*[Title/Abstract] OR macrophage*[Title/Abstract])',
        ),
        (
            "animal_cell",
            '(diabet*[Title/Abstract] OR "insulin resistance"[Title/Abstract]) AND (animal*[Title/Abstract] OR mouse[Title/Abstract] OR mice[Title/Abstract] OR rat[Title/Abstract] OR cell*[Title/Abstract] OR in vitro[Title/Abstract])',
        ),
        (
            "target",
            '(diabet*[Title/Abstract] OR "insulin resistance"[Title/Abstract]) AND (target[Title/Abstract] OR receptor[Title/Abstract] OR pathway[Title/Abstract] OR signaling[Title/Abstract])',
        ),
        (
            "network_bioinformatics",
            '("network toxicology"[Title/Abstract] OR bioinformatics[Title/Abstract] OR docking[Title/Abstract] OR "molecular docking"[Title/Abstract])',
        ),
    ]
)


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def query_pubmed(term: str, retmax: int = 5) -> Dict[str, object]:
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "retmax": str(retmax),
        "sort": "relevance",
        "tool": "whynot17_t2d_opportunity_audit",
        "email": "2163421056@qq.com",
    }
    url = NCBI_BASE + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": "whynot17-t2d-opportunity-audit/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=45) as response:
            payload = json.loads(response.read().decode("utf-8"))
        result = payload.get("esearchresult", {})
        return {
            "status": "ok",
            "count": int(result.get("count", 0)),
            "ids": ";".join(result.get("idlist", [])),
            "url": url,
            "error": "",
        }
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        return {"status": "error", "count": "", "ids": "", "url": url, "error": repr(exc)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    repo = args.repo.resolve()
    framework = repo / "analysis" / "disease_agnostic_environmental_framework"
    root = framework / "t2d_exposure_opportunity"
    candidate_path = root / "01_candidate_master" / "unique_candidate_chemicals.csv"
    candidate_rows = read_csv(candidate_path)
    candidates = {str(row.get("chemical_id", "")): row for row in candidate_rows if str(row.get("mapping_gate_disposition", "")).startswith("advance")}
    selected_ids = set().union(*(set(spec["chemical_ids"]) for spec in SEARCH_GROUPS.values()))
    missing_ids = sorted(selected_ids - set(candidates))
    if missing_ids:
        raise RuntimeError(f"Candidate IDs missing from Phase 1 advance pool: {missing_ids}")

    count_rows: List[Dict[str, object]] = []
    errors: List[Dict[str, object]] = []
    group_results: Dict[str, Dict[str, Dict[str, object]]] = {}
    for group_id, spec in SEARCH_GROUPS.items():
        group_results[group_id] = {}
        base_alias = f"({spec['aliases']})"
        for category, category_terms in SEARCH_CATEGORIES.items():
            term = f"{base_alias} AND ({category_terms})"
            result = query_pubmed(term)
            row = {
                "search_date": SEARCH_DATE,
                "search_group": group_id,
                "chemical_ids": ";".join(spec["chemical_ids"]),
                "category": category,
                "query": term,
                "result_count": result["count"],
                "top_pmids": result["ids"],
                "status": result["status"],
                "error": result["error"],
                "source_url": result["url"],
            }
            count_rows.append(row)
            group_results[group_id][category] = result
            if result["status"] != "ok":
                errors.append(row)
            time.sleep(0.36)

    count_fields = ["search_date", "search_group", "chemical_ids", "category", "query", "result_count", "top_pmids", "status", "error", "source_url"]
    counts_path = root / "03_literature_collision" / "pubmed_collision_counts.csv"
    write_csv(counts_path, count_fields, count_rows)

    candidate_count_rows: List[Dict[str, object]] = []
    for group_id, spec in SEARCH_GROUPS.items():
        values = group_results[group_id]
        for chemical_id in spec["chemical_ids"]:
            candidate = candidates[chemical_id]
            all_count = values["diabetes_total"]["count"]
            human_count = values["human_epidemiology"]["count"]
            prospective_count = values["prospective"]["count"]
            mechanism_count = values["mechanism"]["count"]
            animal_cell_count = values["animal_cell"]["count"]
            target_count = values["target"]["count"]
            network_count = values["network_bioinformatics"]["count"]
            if all_count == "":
                epi_level = "not_estimable_search_error"
            elif int(human_count or 0) == 0:
                epi_level = "0_no_human_epidemiology_hit"
            elif int(prospective_count or 0) == 0:
                epi_level = "2_human_hit_no_prospective_hit"
            elif int(prospective_count or 0) == 1:
                epi_level = "3_one_prospective_screening_hit"
            else:
                epi_level = "4_multiple_prospective_screening_hits"
            total_numeric = int(all_count or 0) if all_count != "" else -1
            if total_numeric == 0:
                crowding = "no_diabetes_related_hit"
            elif total_numeric <= 10:
                crowding = "sparse_screening_signal"
            elif total_numeric <= 50:
                crowding = "intermediate_screening_signal"
            else:
                crowding = "crowded_screening_signal"
            candidate_count_rows.append(
                {
                    "chemical_id": chemical_id,
                    "chemical_name": candidate["chemical_name"],
                    "chemical_class": candidate["chemical_class"],
                    "search_group": group_id,
                    "positive_biomarkers": candidate["positive_biomarkers"],
                    "mapping_grades": candidate["mapping_grades"],
                    "mapping_gate_disposition": candidate["mapping_gate_disposition"],
                    "diabetes_total_pubmed_hits": all_count,
                    "human_epidemiology_pubmed_hits": human_count,
                    "prospective_pubmed_hits": prospective_count,
                    "mechanism_pubmed_hits": mechanism_count,
                    "animal_cell_pubmed_hits": animal_cell_count,
                    "target_pubmed_hits": target_count,
                    "network_bioinformatics_pubmed_hits": network_count,
                    "screening_epidemiology_level": epi_level,
                    "crowding_screen_flag": crowding,
                    "manual_title_abstract_screen_required": "True",
                    "search_date": SEARCH_DATE,
                }
            )
    literature_fields = list(candidate_count_rows[0].keys()) if candidate_count_rows else []
    literature_path = root / "03_literature_collision" / "literature_counts.csv"
    write_csv(literature_path, literature_fields, candidate_count_rows)

    vocab_rows = []
    for group_id, spec in SEARCH_GROUPS.items():
        for chemical_id in spec["chemical_ids"]:
            candidate = candidates[chemical_id]
            vocab_rows.append(
                {
                    "chemical_id": chemical_id,
                    "chemical_name": candidate["chemical_name"],
                    "search_group": group_id,
                    "aliases": spec["aliases"],
                    "mapping_gate_disposition": candidate["mapping_gate_disposition"],
                    "source": "Phase 1 unique_candidate_chemicals.csv",
                }
            )
    vocab_path = root / "03_literature_collision" / "candidate_search_vocabulary.csv"
    write_csv(vocab_path, ["chemical_id", "chemical_name", "search_group", "aliases", "mapping_gate_disposition", "source"], vocab_rows)

    report = f"""# Phase 2 — PubMed literature collision audit

## Scope

- **Search date (UTC):** `{SEARCH_DATE}`
- **Candidates:** 15 chemical IDs in the Phase 1 advance pool, grouped into 14 search groups because the two DINP parent IDs share the same alias set.
- **Database:** PubMed via NCBI E-utilities `esearch.fcgi`.
- **Search categories:** diabetes-related total, human epidemiology, prospective/longitudinal, mechanism, animal/cell, target/pathway, and network/bioinformatics.
- **Interpretation boundary:** retrieval counts are screening signals, not counts of eligible studies. Every candidate remains subject to title/abstract and, where needed, full-text adjudication.

## Search design

Each query was constructed as `(candidate alias set) AND (category terms)`. The complete query, result count, top relevance-sorted PMIDs, status, and source URL are retained in `pubmed_collision_counts.csv`. DINP parent IDs `C012125` and `C019174` deliberately share the `dinp_parent` search group to avoid duplicate literature counts while preserving both chemical identities.

## Outputs

- `03_literature_collision/pubmed_collision_counts.csv`: one row per search group × category, 98 queries in total.
- `03_literature_collision/literature_counts.csv`: one row per chemical ID with category-level retrieval counts and conservative screening flags.
- `03_literature_collision/candidate_search_vocabulary.csv`: chemical IDs, names, search groups, and aliases.

## Important limitation

The category filters are intentionally sensitive and can retrieve papers that mention an analyte, a mixture, a therapeutic formulation, or a non-T2D metabolic endpoint. Therefore `human_epidemiology_pubmed_hits`, `prospective_pubmed_hits`, and `mechanism_pubmed_hits` are not yet adjudicated evidence counts. A later screening pass must verify whether the exposure is the candidate chemical, whether the outcome is T2D/diabetes, study design, and whether a complete chemical→T2D mechanism exists.

## Phase 2 decision status

This audit establishes the reproducible collision-search layer. It does **not** yet assign final novelty grades, docking feasibility, target status, or a definitive Top 5. Those require review of the retrieved records rather than ranking candidates from raw search counts alone.
"""
    report_path = root / "03_literature_collision" / "PHASE2_PUBMED_COLLISION_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    generated = [counts_path, literature_path, vocab_path, report_path]
    manifest = {
        "material_id": "T2D-EXPOSURE-OPPORTUNITY-PHASE2-PUBMED",
        "script_version": SCRIPT_VERSION,
        "status": "complete_first_pass_pubmed_collision_audit",
        "search_date_utc": SEARCH_DATE,
        "candidate_chemical_count": len(candidate_count_rows),
        "search_group_count": len(SEARCH_GROUPS),
        "category_count": len(SEARCH_CATEGORIES),
        "query_count": len(count_rows),
        "successful_query_count": sum(row["status"] == "ok" for row in count_rows),
        "failed_query_count": len(errors),
        "search_groups": {key: spec["chemical_ids"] for key, spec in SEARCH_GROUPS.items()},
        "categories": list(SEARCH_CATEGORIES),
        "input": {"path": str(candidate_path.relative_to(repo)), "sha256": sha256_file(candidate_path)},
        "outputs": {str(path.relative_to(repo)): {"sha256": sha256_file(path), "bytes": path.stat().st_size} for path in generated},
        "not_performed": ["title_abstract_adjudication", "full_text_screening", "mechanism_validation", "target_discovery", "docking", "experimental_feasibility", "final_opportunity_score"],
        "errors": errors,
    }
    (root / "03_literature_collision" / "PHASE2_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": manifest["status"],
        "candidate_chemical_count": len(candidate_count_rows),
        "search_group_count": len(SEARCH_GROUPS),
        "query_count": len(count_rows),
        "successful_query_count": manifest["successful_query_count"],
        "failed_query_count": manifest["failed_query_count"],
        "output_dir": str(root / "03_literature_collision"),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
