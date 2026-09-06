#!/usr/bin/env python3
"""Retrieve titles for the relevance-ranked PubMed records in Phase 2."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def esummary(ids):
    params = {
        "db": "pubmed",
        "id": ",".join(ids),
        "retmode": "json",
        "tool": "whynot17_t2d_opportunity_audit",
        "email": "2163421056@qq.com",
    }
    url = ESUMMARY_URL + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": "whynot17-t2d-opportunity-audit/1.0"})
    with urllib.request.urlopen(request, timeout=45) as response:
        payload = json.loads(response.read().decode("utf-8"))
    return payload.get("result", {}), url


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    repo = args.repo.resolve()
    root = repo / "analysis" / "disease_agnostic_environmental_framework" / "t2d_exposure_opportunity"
    counts_path = root / "03_literature_collision" / "pubmed_collision_counts.csv"
    counts = read_csv(counts_path)
    wanted_categories = {"diabetes_total", "human_epidemiology", "prospective", "mechanism", "network_bioinformatics"}
    by_group = defaultdict(list)
    for row in counts:
        if row["category"] not in wanted_categories or row["status"] != "ok":
            continue
        for pmid in row["top_pmids"].split(";"):
            if pmid:
                by_group[(row["search_group"], row["category"])].append((pmid, row))

    rows = []
    failed = []
    for (group, category), pairs in sorted(by_group.items()):
        ids = [pmid for pmid, _ in pairs]
        try:
            summary, source_url = esummary(ids)
            for pmid, query_row in pairs:
                info = summary.get(pmid, {})
                rows.append(
                    {
                        "retrieval_date": datetime.now(timezone.utc).date().isoformat(),
                        "search_group": group,
                        "category": category,
                        "pmid": pmid,
                        "title": info.get("title", ""),
                        "pubdate": info.get("pubdate", ""),
                        "sort_first_author": info.get("sortfirstauthor", ""),
                        "source_query": query_row["query"],
                        "source_url": source_url,
                        "status": "ok" if info else "missing_summary",
                    }
                )
        except Exception as exc:  # record the failure without losing prior groups
            failed.append({"search_group": group, "category": category, "error": repr(exc)})
        time.sleep(0.36)

    fields = ["retrieval_date", "search_group", "category", "pmid", "title", "pubdate", "sort_first_author", "source_query", "source_url", "status"]
    out_path = root / "03_literature_collision" / "top_pubmed_records.csv"
    write_csv(out_path, fields, rows)
    failures_path = root / "03_literature_collision" / "top_pubmed_retrieval_errors.csv"
    write_csv(failures_path, ["search_group", "category", "error"], failed)

    report_path = root / "03_literature_collision" / "TOP_PUBMED_RECORDS_REPORT.md"
    report_text = (
        f"# Top PubMed record retrieval\n\n"
        f"- Retrieval date (UTC): {datetime.now(timezone.utc).date().isoformat()}\n"
        f"- Search groups: {len(by_group)} group/category batches\n"
        f"- Records requested: {sum(len(v) for v in by_group.values())}\n"
        f"- Records with summaries: {len(rows)}\n"
        f"- Retrieval errors: {len(failed)}\n\n"
        "The records are the relevance-ranked top PMIDs from the Phase 2 PubMed search log. They are screening aids, not manually adjudicated eligible studies. Titles must be reviewed for candidate identity, T2D relevance, exposure type, and study design before assigning final literature grades.\n"
    )
    report_path.write_text(report_text, encoding="utf-8")
    generated = [out_path, failures_path, report_path]
    manifest = {
        "script": "retrieve_top_pubmed_records.py",
        "status": "complete_top_record_retrieval",
        "retrieval_date_utc": datetime.now(timezone.utc).date().isoformat(),
        "input": {"path": str(counts_path.relative_to(repo)), "sha256": sha256_file(counts_path)},
        "batch_count": len(by_group),
        "records_requested": sum(len(v) for v in by_group.values()),
        "records_with_summaries": len(rows),
        "retrieval_error_count": len(failed),
        "outputs": {str(path.relative_to(repo)): {"sha256": sha256_file(path), "bytes": path.stat().st_size} for path in generated},
        "not_performed": ["title_abstract_adjudication", "full_text_screening", "final_novelty_grade"],
    }
    (root / "03_literature_collision" / "TOP_PUBMED_RECORDS_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "batch_count": len(by_group), "records_requested": manifest["records_requested"], "records_with_summaries": manifest["records_with_summaries"], "retrieval_error_count": len(failed)}, indent=2))


if __name__ == "__main__":
    main()
