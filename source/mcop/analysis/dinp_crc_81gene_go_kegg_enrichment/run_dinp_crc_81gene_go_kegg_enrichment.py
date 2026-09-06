#!/usr/bin/env python
"""ORA of the frozen 81-gene DINP--CRC intersection.

The analysis uses g:Profiler's human GO and KEGG annotations with a custom
gene universe.  The primary universe is the union of the frozen DINP exposure
gene set and the frozen CRC disease-gene union used in the three-source
convergence analysis.  A CRC-union-only background is retained as a
pre-specified sensitivity analysis.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs"
OUT_DIR = Path(__file__).resolve().parent / "outputs"
API_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG"]
USER_AGENT = "whynot17-dinp-crc-81gene-ora/1.0"


def unique_symbols(values: Iterable[object]) -> List[str]:
    values = [str(v).strip().upper() for v in values if pd.notna(v) and str(v).strip()]
    return sorted(set(values))


def bh_with_fixed_family(p_values: pd.Series, family_size: int) -> pd.Series:
    """BH adjustment using the full pre-specified family size.

    g:Profiler returns terms with a non-zero overlap.  The API metadata gives
    the number of terms tested in each source, so the adjustment below keeps
    the full tested family as the denominator instead of silently dropping
    unreturned/no-overlap terms.
    """
    out = pd.Series(float("nan"), index=p_values.index, dtype=float)
    valid = p_values.notna() & (p_values >= 0) & (p_values <= 1)
    if not valid.any() or family_size <= 0:
        return out
    vals = p_values.loc[valid].to_numpy(dtype=float)
    order = vals.argsort(kind="mergesort")
    ranked = vals[order] * float(family_size) / (pd.Series(range(1, len(vals) + 1)).to_numpy())
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.clip(ranked, 0.0, 1.0)
    adjusted = ranked[order.argsort()]
    out.loc[valid] = adjusted
    return out


def request_gprofiler(query: List[str], background: List[str]) -> Dict:
    payload = {
        "organism": "hsapiens",
        "query": query,
        "background": background,
        "sources": SOURCES,
        "domain_scope": "custom",
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "no_evidences": True,
        "all_results": True,
        "ordered": False,
        "combined": False,
        "measure_underrepresentation": False,
    }
    response = requests.post(
        API_URL,
        json=payload,
        headers={"User-Agent": USER_AGENT},
        timeout=300,
    )
    response.raise_for_status()
    result = response.json()
    result["_request_payload"] = payload
    result["_http_status"] = response.status_code
    result["_response_sha256"] = hashlib.sha256(response.content).hexdigest()
    return result


def flatten_results(response: Dict, background_label: str) -> pd.DataFrame:
    meta = response.get("meta", {})
    result_metadata = meta.get("result_metadata", {})
    rows = []
    for record in response.get("result", []):
        source = record.get("source")
        family_size = int(result_metadata.get(source, {}).get("number_of_terms", 0))
        effective_domain = int(record.get("effective_domain_size", 0))
        term_size = int(record.get("term_size", 0))
        query_size = int(record.get("query_size", 0))
        overlap = int(record.get("intersection_size", 0))
        raw_p = float("nan")
        if effective_domain > 0 and term_size >= 0 and query_size >= 0 and overlap > 0:
            raw_p = float(hypergeom.sf(overlap - 1, effective_domain, term_size, query_size))
        rows.append(
            {
                "background": background_label,
                "source": source,
                "term_id": record.get("native"),
                "term_name": record.get("name"),
                "description": record.get("description"),
                "term_size": term_size,
                "effective_domain_size": effective_domain,
                "query_size": query_size,
                "intersection_size": overlap,
                "precision": record.get("precision"),
                "recall": record.get("recall"),
                "raw_p_hypergeom": raw_p,
                "gprofiler_fdr": record.get("p_value"),
                "gprofiler_significant": record.get("significant"),
                "parents": ";".join(record.get("parents", []) or []),
                "source_order": record.get("source_order"),
                "tested_term_family_size": family_size,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    total_family = int(sum(int(v.get("number_of_terms", 0)) for v in result_metadata.values()))
    df["BH_FDR_all_GO_KEGG"] = bh_with_fixed_family(df["raw_p_hypergeom"], total_family)
    df["BH_FDR_within_source"] = float("nan")
    for source, idx in df.groupby("source").groups.items():
        family_size = int(result_metadata.get(source, {}).get("number_of_terms", 0))
        df.loc[idx, "BH_FDR_within_source"] = bh_with_fixed_family(
            df.loc[idx, "raw_p_hypergeom"], family_size
        )
    df["BH_FDR_all_GO_KEGG_significant"] = df["BH_FDR_all_GO_KEGG"] < 0.05
    df["BH_FDR_within_source_significant"] = df["BH_FDR_within_source"] < 0.05
    return df.sort_values(["BH_FDR_all_GO_KEGG", "raw_p_hypergeom", "source", "term_id"]).reset_index(drop=True)


def write_response(response: Dict, path: Path) -> None:
    serializable = {k: v for k, v in response.items() if not k.startswith("_")}
    path.write_text(json.dumps(serializable, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    intersection = pd.read_csv(INPUT_DIR / "dinp_crc_intersection.csv")
    exposure = pd.read_csv(INPUT_DIR / "dinp_exposure_gene_matrix.csv")
    disease = pd.read_csv(INPUT_DIR / "crc_gene_matrix.csv")

    query_genes = unique_symbols(intersection["gene_symbol"])
    exposure_genes = unique_symbols(exposure["gene_symbol"])
    crc_genes = unique_symbols(disease["gene_symbol"])
    combined_background = sorted(set(exposure_genes) | set(crc_genes))
    if len(query_genes) != 81:
        raise ValueError(f"Expected 81 input genes, found {len(query_genes)}")
    if not set(query_genes).issubset(set(combined_background)):
        raise ValueError("Some query genes are absent from the combined background")
    if not set(query_genes).issubset(set(crc_genes)):
        raise ValueError("Some query genes are absent from the CRC background")

    pd.DataFrame({"gene_symbol": query_genes}).to_csv(OUT_DIR / "input_81_genes.csv", index=False)
    pd.DataFrame({"gene_symbol": combined_background}).to_csv(
        OUT_DIR / "background_combined_exposure_crc.csv", index=False
    )
    pd.DataFrame({"gene_symbol": crc_genes}).to_csv(OUT_DIR / "background_crc_union.csv", index=False)

    analyses = {
        "primary_combined_exposure_crc": combined_background,
        "sensitivity_crc_union": crc_genes,
    }
    all_frames = []
    manifest = {
        "analysis": "81-gene DINP-CRC intersection GO and KEGG over-representation analysis",
        "method": "ORA; hypergeometric enrichment with custom background",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "endpoint": API_URL,
        "organism": "hsapiens",
        "input_gene_count": len(query_genes),
        "input_gene_sha256": hashlib.sha256("\n".join(query_genes).encode()).hexdigest(),
        "sources": SOURCES,
        "multiple_testing": {
            "api_method": "Benjamini-Hochberg FDR via g:Profiler significance_threshold_method=fdr",
            "reported_primary_columns": "BH_FDR_all_GO_KEGG and BH_FDR_within_source",
            "custom_family_sizes": "g:Profiler result_metadata number_of_terms; no-overlap terms retained in denominator",
        },
        "analyses": {},
    }

    for label, background in analyses.items():
        response = request_gprofiler(query_genes, background)
        write_response(response, OUT_DIR / f"api_response_{label}.json")
        (OUT_DIR / f"request_payload_{label}.json").write_text(
            json.dumps(response["_request_payload"], indent=2, ensure_ascii=False), encoding="utf-8"
        )
        frame = flatten_results(response, label)
        frame.to_csv(OUT_DIR / f"enrichment_{label}.csv", index=False)
        for source, name in [("GO:BP", "go_bp"), ("GO:MF", "go_mf"), ("GO:CC", "go_cc"), ("KEGG", "kegg")]:
            frame.loc[frame["source"] == source].to_csv(
                OUT_DIR / f"{name}_{label}.csv", index=False
            )
        all_frames.append(frame)
        meta = response.get("meta", {})
        result_metadata = meta.get("result_metadata", {})
        manifest["analyses"][label] = {
            "background_gene_count": len(background),
            "background_gene_sha256": hashlib.sha256("\n".join(background).encode()).hexdigest(),
            "query_mapped_failed": meta.get("genes_metadata", {}).get("failed", []),
            "query_duplicates": meta.get("genes_metadata", {}).get("duplicates", []),
            "response_sha256": response.get("_response_sha256"),
            "request_payload_file": f"request_payload_{label}.json",
            "gprofiler_version": meta.get("version"),
            "gprofiler_timestamp": meta.get("timestamp"),
            "result_metadata": result_metadata,
            "returned_term_count": int(len(frame)),
            "returned_by_source": frame["source"].value_counts().to_dict() if not frame.empty else {},
            "significant_global_count": int(frame["BH_FDR_all_GO_KEGG_significant"].sum()) if not frame.empty else 0,
            "significant_within_source_count": int(frame["BH_FDR_within_source_significant"].sum()) if not frame.empty else 0,
        }

    combined = pd.concat(all_frames, ignore_index=True)
    combined.to_csv(OUT_DIR / "enrichment_all_backgrounds.csv", index=False)
    primary = all_frames[0].copy()
    primary = primary[primary["term_id"] != "KEGG:00000"]
    headline = (
        primary.sort_values(["source", "BH_FDR_all_GO_KEGG", "raw_p_hypergeom"])
        .groupby("source", group_keys=False)
        .head(10)
        .reset_index(drop=True)
    )
    headline.to_csv(OUT_DIR / "headline_terms_primary.csv", index=False)
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    summary_lines = [
        "# 81-gene DINP–CRC intersection: GO and KEGG enrichment",
        "",
        "## Scope",
        "",
        "ORA was run on the frozen 81-gene intersection from the three-source DINP–CRC convergence analysis.",
        "The primary custom universe is the union of the frozen DINP exposure-gene union (86 genes) and the frozen CRC disease-gene union (15,885 genes), yielding 15,890 unique symbols (the two source unions overlap substantially).",
        "A CRC-union-only universe (15,885 genes) is reported as sensitivity analysis.",
        "",
        "## Method",
        "",
        "Human GO Biological Process (GO:BP), Molecular Function (GO:MF), Cellular Component (GO:CC), and KEGG annotations were queried through g:Profiler with `domain_scope=custom`.",
        "The API was asked for all returned terms (`all_results=true`, threshold 1.0); hypergeometric raw P values were reconstructed from the returned effective domain, term size, query size, and overlap.",
        "BH-FDR was then recomputed across the full pre-specified tested term family using the source-specific term counts in the API metadata; no-overlap terms remain in that denominator.",
        "",
        "## Primary readout",
        "",
        "Use `enrichment_primary_combined_exposure_crc.csv` for the primary analysis. The principal columns are `raw_p_hypergeom`, `BH_FDR_all_GO_KEGG`, and `BH_FDR_within_source`; `intersection_size` is the number of query genes overlapping each term.",
        "GO branches and KEGG are also emitted as separate CSV files. Term counts and significant-term counts are recorded in `manifest.json`.",
        "",
        "## Reproducibility",
        "",
        "The run was performed against g:Profiler's live API; the exact request payload, API response, timestamp, and response SHA-256 are retained in `manifest.json` and the two `api_response_*.json` files. The g:Profiler source release is recorded in each analysis manifest.",
        "",
        "This analysis is functional enrichment of the 81-gene intersection; it does not establish DINP causality or direction of regulation.",
    ]
    (OUT_DIR / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "query_genes": len(query_genes),
        "combined_background": len(combined_background),
        "crc_background": len(crc_genes),
        "primary_terms": int(len(all_frames[0])),
        "primary_global_fdr_05": int(all_frames[0]["BH_FDR_all_GO_KEGG_significant"].sum()),
        "primary_source_counts": all_frames[0]["source"].value_counts().to_dict(),
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
