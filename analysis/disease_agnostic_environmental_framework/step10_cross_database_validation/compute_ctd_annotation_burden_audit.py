#!/usr/bin/env python
"""Audit whether CTD annotation density could bias candidate prioritisation.

This is a descriptive post-selection audit.  It does not change any frozen
screen, FDR family, robustness tier, or candidate ranking.  The comparison
population is the outcome-free actionability ledger (all CTD chemicals that
entered the formal environmental workflow); the T2D candidate flag is carried
through only to quantify possible annotation-density imbalance.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CTD_IXNS = ROOT / "work/environmental_toxicology_crc_phase1/data/CTD_chem_gene_ixns.tsv.gz"
DEFAULT_CTD_CHEMICALS = ROOT / "work/environmental_toxicology_crc_phase1/data/CTD_chemicals.tsv.gz"
DEFAULT_LEDGER = ROOT / "analysis/disease_agnostic_environmental_framework/step03_actionability/actionability_ledger.csv"
DEFAULT_CLASSIFICATION = ROOT / "analysis/disease_agnostic_environmental_framework/clean_room_reconstruction/outputs/clean_room_classification_ledger.csv"
DEFAULT_CANDIDATES = ROOT / "analysis/disease_agnostic_environmental_framework/t2d_exposure_opportunity/01_candidate_master/unique_candidate_chemicals.csv"

PMID_RE = re.compile(r"\d+")


def ctd_header(path: Path) -> list[str]:
    """Read the commented CTD header without treating it as data."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "ChemicalID" in line:
                return line.lstrip("# \t").rstrip("\r\n").split("\t")
    raise ValueError(f"No CTD header found in {path}")


def split_pipe(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text or text in {"-", "NA", "NaN"}:
        return []
    return [token.strip() for token in text.split("|") if token.strip() and token.strip() != "-"]


def parse_pmids(value: object) -> set[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return set()
    return set(PMID_RE.findall(str(value)))


def clean_text(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    text = str(value).strip()
    return "" if text in {"", "-", "NA", "NaN"} else text


def first_nonempty(series: pd.Series) -> str:
    for value in series:
        text = clean_text(value)
        if text:
            return text
    return ""


def join_unique(series: pd.Series, limit: int = 20) -> str:
    values: list[str] = []
    seen: set[str] = set()
    for value in series:
        text = clean_text(value)
        if text and text not in seen:
            values.append(text)
            seen.add(text)
    return "|".join(values[:limit])


def load_ledger(path: Path) -> pd.DataFrame:
    ledger = pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    ledger["actionable_mapping_bool"] = ledger["actionable_mapping"].str.lower().eq("true")
    aggregations: dict[str, tuple[str, object]] = {
        "chemical_name": ("chemical_name", first_nonempty),
        "actionable_mapping": ("actionable_mapping_bool", "max"),
        "n_mapping_rows": ("chemical_id", "size"),
        "n_actionable_mapping_rows": ("actionable_mapping_bool", "sum"),
        "human_biomarker": ("human_biomarker", first_nonempty),
        "mapping_status": ("mapping_status", join_unique),
        "weight_variable": ("weight_variable", join_unique),
        "matrix": ("human_biomarker", lambda s: ""),
    }
    grouped = ledger.groupby("chemical_id", sort=True).agg(**aggregations).reset_index()
    grouped["actionable_mapping"] = grouped["actionable_mapping"].astype(bool)
    return grouped


def load_classification(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["chemical_id", "chemical_class"])
    columns = ["chemical_id", "chemical_class"]
    frame = pd.read_csv(path, usecols=lambda c: c in columns, encoding="utf-8-sig", dtype=str)
    return frame.drop_duplicates("chemical_id")


def load_candidates(path: Path) -> pd.DataFrame:
    candidates = pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    keep = [
        "chemical_id", "chemical_name", "chemical_class", "positive_biomarkers",
        "exposure_axes", "NHANES_variables", "T2D_OR", "T2D_P", "T2D_FDR",
        "T2D_analytic_N", "T2D_case_N", "T2D_robustness_statuses", "T2D_robustness_tiers",
    ]
    keep = [column for column in keep if column in candidates.columns]
    result = candidates[keep].drop_duplicates("chemical_id").copy()
    result["is_t2d_candidate"] = True
    return result


def empty_stat() -> dict[str, object]:
    return {
        "n_raw_human_interaction_rows": 0,
        "gene_pairs": set(),
        "gene_ids": set(),
        "gene_symbols": set(),
        "pmids": set(),
        "actions": set(),
        "interaction_terms": set(),
    }


def read_ctd_interaction_stats(path: Path, chemical_ids: set[str], chunksize: int = 200_000) -> dict[str, dict[str, object]]:
    header = ctd_header(path)
    required = ["ChemicalID", "GeneID", "GeneSymbol", "Organism", "OrganismID", "Interaction", "InteractionActions", "PubMedIDs"]
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"CTD interaction file is missing columns: {missing}")

    stats: dict[str, dict[str, object]] = defaultdict(empty_stat)
    usecols = [column for column in required if column in header]
    for chunk in pd.read_csv(
        path,
        sep="\t",
        names=header,
        header=None,
        comment="#",
        usecols=usecols,
        dtype=str,
        keep_default_na=False,
        chunksize=chunksize,
        compression="gzip",
    ):
        chunk = chunk[chunk["ChemicalID"].isin(chemical_ids)]
        if chunk.empty:
            continue
        organism_id = chunk["OrganismID"].astype(str).str.strip()
        organism = chunk["Organism"].astype(str).str.lower().str.strip()
        chunk = chunk[organism_id.eq("9606") | organism.eq("homo sapiens")]
        for row in chunk.itertuples(index=False):
            record = row._asdict()
            chemical_id = clean_text(record.get("ChemicalID"))
            gene_id = clean_text(record.get("GeneID"))
            gene_symbol = clean_text(record.get("GeneSymbol"))
            if not chemical_id or not (gene_id or gene_symbol):
                continue
            target = stats[chemical_id]
            target["n_raw_human_interaction_rows"] = int(target["n_raw_human_interaction_rows"]) + 1
            pair_gene = gene_id or gene_symbol
            target["gene_pairs"].add(f"{chemical_id}::{pair_gene}")
            if gene_id:
                target["gene_ids"].add(gene_id)
            if gene_symbol:
                target["gene_symbols"].add(gene_symbol)
            target["pmids"].update(parse_pmids(record.get("PubMedIDs")))
            target["actions"].update(split_pipe(record.get("InteractionActions")))
            interaction = clean_text(record.get("Interaction"))
            if interaction:
                target["interaction_terms"].add(interaction)
    return stats


def read_ctd_vocabulary(path: Path, chemical_ids: set[str]) -> pd.DataFrame:
    header = ctd_header(path)
    required = [
        "ChemicalID", "ChemicalName", "Definition", "ParentIDs", "TreeNumbers",
        "ParentTreeNumbers", "MESHSynonyms", "CTDCuratedSynonyms",
    ]
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"CTD chemical vocabulary is missing columns: {missing}")
    pieces: list[pd.DataFrame] = []
    usecols = required
    for chunk in pd.read_csv(
        path,
        sep="\t",
        names=header,
        header=None,
        comment="#",
        usecols=usecols,
        dtype=str,
        keep_default_na=False,
        chunksize=100_000,
        compression="gzip",
    ):
        selected = chunk[chunk["ChemicalID"].isin(chemical_ids)]
        if not selected.empty:
            pieces.append(selected)
    if not pieces:
        return pd.DataFrame(columns=required)
    return pd.concat(pieces, ignore_index=True).drop_duplicates("ChemicalID")


def digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_float(value: object) -> float | None:
    try:
        number = float(value)
        return number if np.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ctd-interactions", type=Path, default=DEFAULT_CTD_IXNS)
    parser.add_argument("--ctd-chemicals", type=Path, default=DEFAULT_CTD_CHEMICALS)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--classification", type=Path, default=DEFAULT_CLASSIFICATION)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ledger = load_ledger(args.ledger)
    classification = load_classification(args.classification)
    candidates = load_candidates(args.candidates)
    candidate_ids = set(candidates["chemical_id"])
    chemical_ids = set(ledger["chemical_id"])

    print(f"Actionability universe: {len(chemical_ids)} chemicals; T2D candidate flag: {len(candidate_ids)}", flush=True)
    print("Reading human CTD chemical–gene interaction rows…", flush=True)
    interaction_stats = read_ctd_interaction_stats(args.ctd_interactions, chemical_ids)
    print("Reading CTD chemical vocabulary…", flush=True)
    vocabulary = read_ctd_vocabulary(args.ctd_chemicals, chemical_ids)

    rows: list[dict[str, object]] = []
    vocabulary = vocabulary.set_index("ChemicalID") if not vocabulary.empty else vocabulary
    for chemical_id in sorted(chemical_ids):
        record = interaction_stats.get(chemical_id, empty_stat())
        vocab = vocabulary.loc[chemical_id] if not vocabulary.empty and chemical_id in vocabulary.index else {}
        row: dict[str, object] = {
            "chemical_id": chemical_id,
            "n_raw_human_interaction_rows": int(record["n_raw_human_interaction_rows"]),
            "n_unique_human_gene_pairs": len(record["gene_pairs"]),
            "n_unique_human_genes": len(record["gene_ids"] or record["gene_symbols"]),
            "n_unique_gene_symbols": len(record["gene_symbols"]),
            "n_unique_pmids": len(record["pmids"]),
            "n_unique_interaction_actions": len(record["actions"]),
            "n_unique_interaction_terms": len(record["interaction_terms"]),
            "n_parent_ids": len(split_pipe(vocab.get("ParentIDs", "") if isinstance(vocab, dict) else vocab["ParentIDs"])),
            "n_tree_numbers": len(split_pipe(vocab.get("TreeNumbers", "") if isinstance(vocab, dict) else vocab["TreeNumbers"])),
            "n_parent_tree_numbers": len(split_pipe(vocab.get("ParentTreeNumbers", "") if isinstance(vocab, dict) else vocab["ParentTreeNumbers"])),
            "n_mesh_synonyms": len(split_pipe(vocab.get("MESHSynonyms", "") if isinstance(vocab, dict) else vocab["MESHSynonyms"])),
            "n_ctd_curated_synonyms": len(split_pipe(vocab.get("CTDCuratedSynonyms", "") if isinstance(vocab, dict) else vocab["CTDCuratedSynonyms"])),
            "ctd_definition_characters": len(clean_text(vocab.get("Definition", "") if isinstance(vocab, dict) else vocab["Definition"])),
            "is_t2d_candidate": chemical_id in candidate_ids,
        }
        rows.append(row)
    audit = pd.DataFrame(rows)
    audit = audit.merge(ledger, on="chemical_id", how="left", validate="one_to_one")
    audit = audit.merge(classification, on="chemical_id", how="left", validate="one_to_one")
    audit = audit.merge(candidates, on="chemical_id", how="left", suffixes=("", "_candidate"), validate="one_to_one")
    audit["is_t2d_candidate"] = audit["is_t2d_candidate"].fillna(False).astype(bool)
    if "chemical_name_candidate" in audit:
        audit["chemical_name"] = audit["chemical_name"].where(audit["chemical_name"].astype(str).str.len().gt(0), audit["chemical_name_candidate"])
    if "chemical_class_candidate" in audit:
        audit["chemical_class"] = audit["chemical_class"].where(audit["chemical_class"].astype(str).str.len().gt(0), audit["chemical_class_candidate"])

    numeric = [
        "n_raw_human_interaction_rows", "n_unique_human_gene_pairs", "n_unique_human_genes",
        "n_unique_gene_symbols", "n_unique_pmids", "n_unique_interaction_actions",
        "n_unique_interaction_terms", "n_parent_ids", "n_tree_numbers", "n_parent_tree_numbers",
        "n_mesh_synonyms", "n_ctd_curated_synonyms", "ctd_definition_characters",
    ]
    for column in numeric:
        audit[column] = pd.to_numeric(audit[column], errors="coerce").fillna(0).astype(int)

    actionable = audit["actionable_mapping"].astype(bool)
    for column in ["n_raw_human_interaction_rows", "n_unique_human_gene_pairs", "n_unique_human_genes", "n_unique_pmids"]:
        audit[f"{column}_percentile_all_2042"] = audit[column].rank(method="average", pct=True)
        audit[f"{column}_rank_desc_all_2042"] = audit[column].rank(method="min", ascending=False).astype(int)
        audit[f"{column}_percentile_actionable_409"] = np.nan
        audit.loc[actionable, f"{column}_percentile_actionable_409"] = audit.loc[actionable, column].rank(method="average", pct=True)
        audit[f"{column}_rank_desc_actionable_409"] = np.nan
        audit.loc[actionable, f"{column}_rank_desc_actionable_409"] = audit.loc[actionable, column].rank(method="min", ascending=False).astype(int)

    if audit["chemical_class"].notna().any():
        for column in ["n_unique_human_genes", "n_unique_pmids"]:
            audit[f"{column}_percentile_within_class"] = np.nan
            nonempty = audit["chemical_class"].fillna("").astype(str).str.len().gt(0)
            audit.loc[nonempty, f"{column}_percentile_within_class"] = audit.loc[nonempty].groupby("chemical_class")[column].rank(method="average", pct=True)
    actionable_gene_q90 = audit.loc[actionable, "n_unique_human_genes"].quantile(0.90)
    actionable_pmid_q90 = audit.loc[actionable, "n_unique_pmids"].quantile(0.90)
    audit["high_gene_annotation_burden_top_decile_actionable"] = actionable & audit["n_unique_human_genes"].ge(actionable_gene_q90)
    audit["high_pmid_burden_top_decile_actionable"] = actionable & audit["n_unique_pmids"].ge(actionable_pmid_q90)

    output_full = args.output_dir / "ctd_annotation_burden_by_chemical.csv"
    output_candidates = args.output_dir / "ctd_annotation_burden_candidate_summary.csv"
    audit.sort_values(["is_t2d_candidate", "n_unique_human_genes", "n_unique_pmids", "chemical_id"], ascending=[False, False, False, True]).to_csv(output_full, index=False, encoding="utf-8-sig")
    candidate_columns = [
        "chemical_id", "chemical_name", "chemical_class", "is_t2d_candidate", "actionable_mapping",
        "human_biomarker", "positive_biomarkers", "exposure_axes", "T2D_OR", "T2D_P", "T2D_FDR",
        "T2D_robustness_statuses", "T2D_robustness_tiers", "n_raw_human_interaction_rows",
        "n_unique_human_gene_pairs", "n_unique_human_genes", "n_unique_pmids",
        "n_unique_interaction_actions", "n_unique_interaction_terms",
        "n_unique_human_genes_percentile_all_2042", "n_unique_human_genes_percentile_actionable_409",
        "n_unique_pmids_percentile_all_2042", "n_unique_pmids_percentile_actionable_409",
        "high_gene_annotation_burden_top_decile_actionable", "high_pmid_burden_top_decile_actionable",
    ]
    candidate_columns = [column for column in candidate_columns if column in audit.columns]
    audit.loc[audit["is_t2d_candidate"], candidate_columns].sort_values(["n_unique_human_genes", "n_unique_pmids", "chemical_id"], ascending=[False, False, True]).to_csv(output_candidates, index=False, encoding="utf-8-sig")

    summary: dict[str, object] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "descriptive CTD annotation-burden audit; no candidate reranking",
        "ctd_interaction_scope": "human-only (OrganismID=9606 or Organism=Homo sapiens), unique chemical-gene pairs; raw rows and PMIDs retained separately",
        "actionability_universe_n": int(len(audit)),
        "actionable_chemical_n": int(actionable.sum()),
        "t2d_candidate_n": int(audit["is_t2d_candidate"].sum()),
        "actionable_t2d_candidate_n": int((actionable & audit["is_t2d_candidate"]).sum()),
        "actionable_gene_annotation_q90": safe_float(actionable_gene_q90),
        "actionable_pmid_q90": safe_float(actionable_pmid_q90),
        "candidate_high_gene_burden_n": int((audit["is_t2d_candidate"] & audit["high_gene_annotation_burden_top_decile_actionable"]).sum()),
        "candidate_high_pmid_burden_n": int((audit["is_t2d_candidate"] & audit["high_pmid_burden_top_decile_actionable"]).sum()),
        "candidate_medians": {},
        "actionable_background_medians": {},
        "source_files": {
            "ctd_interactions": str(args.ctd_interactions),
            "ctd_chemicals": str(args.ctd_chemicals),
            "actionability_ledger": str(args.ledger),
            "classification_ledger": str(args.classification),
            "t2d_candidate_chemicals": str(args.candidates),
        },
        "source_sha256": {
            "ctd_interactions": digest_file(args.ctd_interactions),
            "ctd_chemicals": digest_file(args.ctd_chemicals),
            "actionability_ledger": digest_file(args.ledger),
            "classification_ledger": digest_file(args.classification),
            "t2d_candidate_chemicals": digest_file(args.candidates),
        },
    }
    burden_metrics = ["n_raw_human_interaction_rows", "n_unique_human_gene_pairs", "n_unique_human_genes", "n_unique_pmids"]
    for label, mask in [("candidate", audit["is_t2d_candidate"]), ("actionable_background", actionable)]:
        summary[f"{label}_medians"] = {column: safe_float(audit.loc[mask, column].median()) for column in burden_metrics}

    try:
        from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

        noncandidate = actionable & ~audit["is_t2d_candidate"]
        tests: dict[str, object] = {}
        for column in burden_metrics:
            if audit.loc[audit["is_t2d_candidate"], column].notna().any() and audit.loc[noncandidate, column].notna().any():
                result = mannwhitneyu(audit.loc[audit["is_t2d_candidate"], column], audit.loc[noncandidate, column], alternative="two-sided")
                tests[f"candidate_vs_actionable_noncandidate_{column}_mann_whitney_p"] = safe_float(result.pvalue)
        top_gene = actionable & audit["n_unique_human_genes"].ge(actionable_gene_q90)
        table = [
            [int((audit["is_t2d_candidate"] & top_gene).sum()), int((audit["is_t2d_candidate"] & ~top_gene).sum())],
            [int((noncandidate & top_gene).sum()), int((noncandidate & ~top_gene).sum())],
        ]
        odds, pvalue = fisher_exact(table, alternative="greater")
        tests["candidate_enrichment_in_actionable_top_decile_gene_burden_OR"] = safe_float(odds)
        tests["candidate_enrichment_in_actionable_top_decile_gene_burden_p"] = safe_float(pvalue)
        rho, pvalue = spearmanr(audit["n_unique_human_genes"], audit["n_unique_pmids"])
        tests["all_chemical_gene_degree_vs_pmid_spearman_rho"] = safe_float(rho)
        tests["all_chemical_gene_degree_vs_pmid_spearman_p"] = safe_float(pvalue)
        summary["descriptive_tests"] = tests
    except ImportError:
        summary["descriptive_tests"] = {"status": "scipy_unavailable"}

    summary_path = args.output_dir / "CTD_ANNOTATION_BURDEN_AUDIT_SUMMARY.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    manifest = {
        "audit": "CTD annotation burden audit",
        "generated_utc": summary["generated_utc"],
        "script": str(Path(__file__)),
        "inputs": summary["source_files"],
        "input_sha256": summary["source_sha256"],
        "outputs": {
            output_full.name: digest_file(output_full),
            output_candidates.name: digest_file(output_candidates),
            summary_path.name: digest_file(summary_path),
        },
        "rules": [
            "Human CTD interactions only: OrganismID=9606 or Homo sapiens.",
            "Interaction evidence is deduplicated at chemical ID x GeneID (GeneSymbol fallback only if GeneID is absent).",
            "Raw interaction rows, unique PMIDs, actions, and terms are descriptive audit fields, not weights in enrichment.",
            "Reference population is the formal 2042-chemical actionability ledger; actionable status is summarized at chemical level.",
            "T2D candidate flag is post-selection metadata and cannot change any frozen screen, FDR, tier, or rank.",
        ],
        "counts": {
            "actionability_universe_n": summary["actionability_universe_n"],
            "actionable_chemical_n": summary["actionable_chemical_n"],
            "t2d_candidate_n": summary["t2d_candidate_n"],
        },
    }
    (args.output_dir / "CTD_ANNOTATION_BURDEN_AUDIT_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
