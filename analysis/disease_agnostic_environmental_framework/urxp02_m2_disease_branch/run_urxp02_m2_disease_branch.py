"""M2 first step: disease branching of the frozen expanded URXP02 universe."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
UNIVERSE = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_expanded_molecular_universe" / "05_gene_evidence_summary.csv"
THYROID = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_mechanism_gene_discovery" / "02_thyroid_disease_genes.csv"
HYPERTENSION = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_mechanism_gene_discovery" / "03_hypertension_disease_genes.csv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def as_bool(v: str) -> bool:
    return str(v).strip().lower() == "true"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    universe = read_csv(UNIVERSE)
    thyroid_rows = read_csv(THYROID)
    hyper_rows = read_csv(HYPERTENSION)
    thyroid = {r["gene_symbol"]: r for r in thyroid_rows if r.get("gene_symbol")}
    hyper = {r["gene_symbol"]: r for r in hyper_rows if r.get("gene_symbol")}

    branch_rows = []
    thyroid_intersection = []
    hyper_intersection = []
    for u in universe:
        symbol = u["gene_symbol"]
        t = thyroid.get(symbol)
        h = hyper.get(symbol)
        in_t = t is not None
        in_h = h is not None
        if in_t and not in_h:
            branch = "thyroid-specific"
        elif in_h and not in_t:
            branch = "hypertension-specific"
        elif in_t and in_h:
            branch = "shared"
        else:
            branch = "neither"
        row = dict(u)
        row.update({
            "in_thyroid_disease_set": in_t,
            "in_hypertension_disease_set": in_h,
            "branch_class": branch,
            "thyroid_disease_relevance_score": t.get("disease_relevance_score", "") if t else "",
            "thyroid_disease_rank": t.get("rank", "") if t else "",
            "hypertension_disease_relevance_score": h.get("disease_relevance_score", "") if h else "",
            "hypertension_disease_rank": h.get("rank", "") if h else "",
            "thyroid_disease_id": t.get("disease_id", "") if t else "",
            "hypertension_disease_id": h.get("disease_id", "") if h else "",
            "branch_evidence_note": "Branch membership is based on gene symbol membership in the independent disease lists; evidence-layer fields are inherited from the expanded molecular universe.",
        })
        branch_rows.append(row)
        if in_t:
            thyroid_intersection.append(row)
        if in_h:
            hyper_intersection.append(row)

    extra = ["in_thyroid_disease_set", "in_hypertension_disease_set", "branch_class", "thyroid_disease_relevance_score", "thyroid_disease_rank", "hypertension_disease_relevance_score", "hypertension_disease_rank", "thyroid_disease_id", "hypertension_disease_id", "branch_evidence_note"]
    fields = list(universe[0].keys()) + extra
    write_csv(OUT / "01_828_thyroid_intersection.csv", thyroid_intersection, fields)
    write_csv(OUT / "02_828_hypertension_intersection.csv", hyper_intersection, fields)
    write_csv(OUT / "03_branch_gene_classification.csv", branch_rows, fields)

    counts = Counter(r["branch_class"] for r in branch_rows)
    summary = []
    for branch in ["thyroid-specific", "hypertension-specific", "shared", "neither"]:
        rows = [r for r in branch_rows if r["branch_class"] == branch]
        summary.append({
            "branch_class": branch,
            "n_genes": len(rows),
            "n_exact_2NAP_human": sum(as_bool(r["exact_2NAP_human_support"]) for r in rows),
            "n_exact_2NAP_experimental": sum(as_bool(r["exact_2NAP_experimental_support"]) for r in rows),
            "n_parent_naphthalene_support": sum(as_bool(r["parent_naphthalene_support"]) for r in rows),
            "n_multi_source_ge2": sum(int(float(r["number_of_sources"] or 0)) >= 2 for r in rows),
            "n_human_any": sum(as_bool(r["exact_2NAP_human_support"]) for r in rows),
            "n_bioassay_target_support": sum(as_bool(r["bioassay_target_support"]) for r in rows),
            "n_toxicogenomic_support": sum(as_bool(r["toxicogenomic_support"]) for r in rows),
            "n_CTD_support": sum(as_bool(r["CTD_support"]) for r in rows),
        })
    write_csv(OUT / "04_branch_evidence_summary.csv", summary, list(summary[0].keys()))

    report = f"""# URXP02 M2 disease branching\n\nGenerated {datetime.now(timezone.utc).isoformat()}. This step only maps the frozen expanded URXP02 molecular universe to the independent thyroid and hypertension gene lists. No pathway, PPI, module, hub, tissue, cell, or new NHANES analysis was run.\n\n## Branch counts\n\n- Expanded universe: **{len(universe)} genes**\n- 828 ∩ thyroid: **{len(thyroid_intersection)}**\n- 828 ∩ hypertension: **{len(hyper_intersection)}**\n- Thyroid-specific (A−B): **{counts['thyroid-specific']}**\n- Hypertension-specific (B−A): **{counts['hypertension-specific']}**\n- Shared (A∩B): **{counts['shared']}**\n- Neither disease list: **{counts['neither']}**\n\nThe classification table carries the three requested evidence layers from M1b: exact 2-NAP human support, exact 2-NAP experimental support, and multi-source support. Disease membership is symbol-based and does not imply causality or sex specificity.\n\n## Files\n\n`01` and `02` contain the two intersections; `03` contains all 828 genes including `neither`; `04` summarizes evidence-layer counts by branch.\n"""
    report_path = OUT / "URXP02_M2_DISEASE_BRANCH_REPORT.md"
    report_path.write_text(report, encoding="utf-8")
    files = sorted(list(OUT.glob("*.csv")) + [report_path])
    manifest = {
        "analysis": "URXP02 M2 disease branching",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {str(p): sha256(p) for p in [UNIVERSE, THYROID, HYPERTENSION]},
        "disease_sources": {"thyroid": "Open Targets direct associations inherited from M1", "hypertension": "Open Targets direct associations inherited from M1"},
        "counts": {"expanded_universe": len(universe), "thyroid_intersection": len(thyroid_intersection), "hypertension_intersection": len(hyper_intersection), "thyroid_specific": counts["thyroid-specific"], "hypertension_specific": counts["hypertension-specific"], "shared": counts["shared"], "neither": counts["neither"]},
        "evidence_layers": ["exact_2NAP_human_support", "exact_2NAP_experimental_support", "number_of_sources >= 2"],
        "constraints": ["No pathway enrichment", "No PPI/module/hub analysis", "No tissue/cell analysis", "No NHANES model", "No sex-specific molecular claim"],
        "files": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size} for p in files},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
