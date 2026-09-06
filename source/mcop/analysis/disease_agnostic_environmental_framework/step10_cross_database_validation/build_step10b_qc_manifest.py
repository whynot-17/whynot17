#!/usr/bin/env python3
"""Build the small, version-control-friendly Step 10B QC manifest."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE = Path(__file__).resolve().parent
E = BASE / "step10b_environmental_replacement"
D = BASE / "step10b_disease_replacement"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    env = pd.read_csv(E / "step10b_environmental_source_records.csv", dtype=str).fillna("")
    env_summary = pd.read_csv(E / "step10b_environmental_source_summary.csv")
    ot = pd.read_csv(D / "step10b_open_targets_t2d_targets.csv", dtype=str).fillna("")
    gw = pd.read_csv(D / "step10b_gwas_catalog_t2d_genes.csv", dtype=str).fillna("")
    conv = pd.read_csv(D / "step10b_disease_source_cluster_convergence.csv", dtype=str).fillna("")
    e_snapshot = read_json(E / "STEP10B_E_SOURCE_SNAPSHOT.json")
    d_snapshot = read_json(D / "STEP10B_D_SOURCE_SNAPSHOT.json")
    gwas_snapshot = read_json(D / "step10b_gwas_catalog_query_snapshot.json")

    checks = {
        "environmental_candidate_rows": {"observed": int(len(env)), "expected": 134, "pass": len(env) == 134},
        "environmental_source_rows": {"observed": int(len(env_summary)), "expected": 3, "pass": len(env_summary) == 3},
        "open_targets_rows": {"observed": int(len(ot)), "expected_min": 9000, "pass": len(ot) >= 9000},
        "gwas_author_reported_gene_rows": {"observed": int(len(gw)), "expected_min": 1, "pass": len(gw) >= 1},
        "disease_source_cluster_rows": {"observed": int(len(conv)), "expected": 33, "pass": len(conv) == 33},
        "gwas_association_count_returned": {"observed": gwas_snapshot.get("association_count_returned"), "expected_min": 1, "pass": int(gwas_snapshot.get("association_count_returned") or 0) >= 1},
    }
    errors = [name for name, check in checks.items() if not check["pass"]]
    output_files = [
        E / "run_step10b_environmental_replacement.py",
        E / "capture_step10b_bindingdb_metadata.py",
        E / "step10b_environmental_candidate_crosswalk.csv",
        E / "step10b_environmental_source_records.csv",
        E / "step10b_environmental_source_summary.csv",
        E / "step10b_environmental_fixed_topk_retention.csv",
        E / "step10b_environmental_source_dropout.csv",
        E / "step10b_environmental_rank_concordance.csv",
        E / "STEP10B_E_SOURCE_SNAPSHOT.json",
        E / "STEP10B_E_API_CALL_MANIFEST.json",
        E / "STEP10B_E_BINDINGDB_METADATA.json",
        E / "bindingdb_download_page_snapshot.html",
        D / "run_step10b_disease_replacement.py",
        D / "step10b_open_targets_t2d_targets.csv",
        D / "step10b_gwas_catalog_t2d_genes.csv",
        D / "step10b_gwas_catalog_query_snapshot.json",
        D / "step10b_disease_source_coverage.csv",
        D / "step10b_disease_source_cluster_convergence.csv",
        D / "step10b_disease_source_rank_stability.csv",
        D / "STEP10B_D_SOURCE_SNAPSHOT.json",
        D / "STEP10B_D_API_CALL_MANIFEST.json",
        BASE / "STEP10_FROZEN_SOURCE_SET.json",
        BASE / "STEP10_SOURCE_REPLACEMENT_REGISTRY.csv",
        BASE / "build_step10b_qc_manifest.py",
    ]
    output_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in output_files if p.exists()}
    source_summary = {
        "E1_ChEMBL": {"candidate_count": 134, "positive_human_activity_candidates": int((pd.to_numeric(env["chembl_human_activity_count"], errors="coerce").fillna(0) > 0).sum()), "release": e_snapshot.get("sources", {}).get("E1", {}).get("status_payload", {})},
        "E2_BindingDB": {"candidate_count": 134, "positive_human_affinity_candidates": int((pd.to_numeric(env["bindingdb_human_affinity_row_count"], errors="coerce").fillna(0) > 0).sum()), "release_basis": "official download-page snapshot and REST query metadata; date-stamped files are not collapsed into one database-wide release", "absence_is_not_negative": True},
        "E3_PubChem_BioAssay": {"candidate_count": 134, "positive_cid_aid_candidates": int((pd.to_numeric(env["pubchem_aid_count"], errors="coerce").fillna(0) > 0).sum()), "release_basis": "rolling PUG-REST retrieval snapshot", "target_count_not_inferred": True},
        "D1_OpenTargets": {"unique_approved_symbols": int(ot["gene_symbol"].nunique()), "data_release": d_snapshot.get("open_targets", {}).get("data_release_documented_by_source"), "api_version": d_snapshot.get("open_targets", {}).get("api_version_documented_by_source")},
        "D2_GWAS_Catalog": {"author_reported_gene_symbols": int(gw["gene_symbol"].nunique()), "association_count_returned": gwas_snapshot.get("association_count_returned"), "release_basis": "source-native trait association API snapshot; no file release tag was inferred", "gene_rule": "authorReportedGenes with source-native Ensembl/Entrez IDs"},
    }
    manifest = {
        "lock_type": "STEP10B_CROSS_DATABASE_REPLACEMENT_AUDIT",
        "generated_utc": utc(),
        "status": "complete_with_explicit_source_boundaries" if not errors else "qc_failed",
        "checks": checks,
        "errors": errors,
        "frozen_inputs": {
            "source_set": str((BASE / "STEP10_FROZEN_SOURCE_SET.json").relative_to(ROOT)),
            "step7_cluster_gene_input": str((ROOT / "analysis/disease_agnostic_environmental_framework/step07_genecard_convergence/t2d_step7_cluster_ctd_genes.csv").relative_to(ROOT)),
            "genecards_reference": str((ROOT / "analysis/disease_agnostic_environmental_framework/step07_genecard_convergence/t2d_genecards_primary_gene_audit.csv").relative_to(ROOT)),
        },
        "source_summary": source_summary,
        "output_hashes": output_hashes,
        "boundary": "Environmental sources remain semantically separate; disease sources remain separate; replacement results do not change the frozen 29-test family, 11 clusters, or any Tier assignment.",
    }
    (BASE / "STEP10B_MANIFEST.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    report = "\n".join([
        "# Step 10B — Cross-database robustness and annotation-bias stress test",
        "",
        f"Generated: `{utc()}`",
        "",
        "## Environmental replacement (E)",
        "",
        f"The frozen 134-candidate set was queried separately against ChEMBL, BindingDB, and PubChem BioAssay. Observed source-layer coverage: ChEMBL human activity in {source_summary['E1_ChEMBL']['positive_human_activity_candidates']}/134 candidates; BindingDB human affinity in {source_summary['E2_BindingDB']['positive_human_affinity_candidates']}/134; PubChem CID-to-AID membership in {source_summary['E3_PubChem_BioAssay']['positive_cid_aid_candidates']}/134.",
        "",
        "These counts are evidence-source coverage indicators, not comparable biological scores. BindingDB absence is not interpreted as a negative; PubChem CID-to-AID membership is not converted into a human target count. ChEMBL/BindingDB/PubChem are not merged into a single chemical-gene edge list.",
        "",
        "## Disease replacement (D)",
        "",
        f"The T2D concept resolved to `MONDO_0005148` in both source-native routes. Open Targets returned {source_summary['D1_OpenTargets']['unique_approved_symbols']} unique approved symbols (data release {source_summary['D1_OpenTargets']['data_release']}, API {source_summary['D1_OpenTargets']['api_version']}). GWAS Catalog returned {source_summary['D2_GWAS_Catalog']['association_count_returned']} trait associations and {source_summary['D2_GWAS_Catalog']['author_reported_gene_symbols']} author-reported gene symbols in the source-native association collection.",
        "",
        "GeneCards remains the frozen reference. Open Targets and GWAS Catalog are reported as independent disease-knowledge layers with source-native fields retained; their results are not merged into a single disease-gene truth set.",
        "",
        "## QC and interpretation boundary",
        "",
        f"QC checks passed: **{len(checks) - len(errors)}/{len(checks)}**. Exact API response hashes, query metadata, source snapshots, input hashes, and output hashes are in `STEP10B_MANIFEST.json`, `STEP10B_E_SOURCE_SNAPSHOT.json`, and `STEP10B_D_SOURCE_SNAPSHOT.json`.",
        "",
        "This is a post-firewall source-replacement audit. It does not promote/demote candidates, change the 29-test family, recompute epidemiologic FDR, or select a flagship axis.",
    ])
    (BASE / "STEP10B_CROSS_DATABASE_REPLACEMENT_REPORT.md").write_text(report, encoding="utf-8")
    (BASE / "STEP10B_QC_SUMMARY.json").write_text(json.dumps({"generated_utc": utc(), "status": manifest["status"], "checks": checks, "source_summary": source_summary}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Step 10B QC: {manifest['status']}; {len(checks) - len(errors)}/{len(checks)} checks passed", flush=True)


if __name__ == "__main__":
    main()
