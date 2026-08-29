#!/usr/bin/env python3
"""Build the candidate-level Phase 2 collision and opportunity audit.

The PubMed search counts are intentionally retained as screening signals rather
than treated as eligible-study counts.  This script adds a transparent,
title/abstract-level first-pass adjudication for all 14 search groups (15
chemical IDs because the two DINP records share one parent group).  It does
not perform a systematic-review-grade full-text eligibility review, docking,
or new target discovery.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


SCRIPT_VERSION = "phase2_candidate_adjudication_v1.1_internal_collision"

# These labels are a first-pass title/abstract screen of the retrieved records.
# They are deliberately qualitative and each row carries a concrete gap.  The
# counts in literature_counts.csv are not changed by this dictionary.
ADJUDICATION = {
    "mcop": {
        "direct_t2d_relevance": "sparse_candidate_specific_signal",
        "human_epidemiology_status": "one_candidate_specific_diabetes_record_plus_non_target_hit",
        "prospective_status": "no_candidate_specific_prospective_record_retrieved",
        "mechanism_maturity": "no_candidate_specific_t2d_mechanism_retrieved",
        "existing_target_status": "none_confirmed_in_first_pass",
        "network_toxicology_status": "none_confirmed_in_first_pass",
        "novelty_grade": "A_high_novelty_high_evidence_gap",
        "opportunity_bucket": "Top5_1_high_novelty_validation_gap",
        "representative_pmids": "33152652",
        "largest_gap": "Parent-specific MCOP exposure-to-T2D evidence is sparse; no prospective or candidate-specific mechanistic study was retrieved.",
        "decision_reason": "Advance as a high-novelty opportunity, but require parent-specific full-text adjudication before any mechanistic claim.",
    },
    "dinp_parent": {
        "direct_t2d_relevance": "mostly_parent_or_phthalate_family_level",
        "human_epidemiology_status": "human_phthalate_replacement_and_family_level_records; direct_DINP_specificity_unclear",
        "prospective_status": "screening_hits_are_family_level_or_not_DINP_specific",
        "mechanism_maturity": "one_phthalate_substitute_adipogenesis_record; direct_DINP_T2D_mechanism_not_confirmed",
        "existing_target_status": "non_T2D_phthalate_target_literature_present",
        "network_toxicology_status": "DINP_network_record_in_non_T2D_cancer_retrieved",
        "novelty_grade": "A_high_novelty_specificity_uncertain",
        "opportunity_bucket": "Top5_2_high_novelty_specificity_gap",
        "representative_pmids": "25993640;35567983;41601071",
        "largest_gap": "The relevant human and mechanistic records often concern phthalate replacements or the phthalate family, not analytically resolved DINP exposure.",
        "decision_reason": "Advance as a parent-specificity opportunity; do not count family-level phthalate records as direct DINP evidence.",
    },
    "mecpp": {
        "direct_t2d_relevance": "sparse_to_intermediate_metabolite_level_signal",
        "human_epidemiology_status": "human_metabolic_and_gestational_diabetes_records_retrieved;_specificity_requires_full_text",
        "prospective_status": "prospective_or_nested_records_retrieved_but_not_all_T2D_or_candidate_specific",
        "mechanism_maturity": "emerging_candidate_specific_metabolic_network_record",
        "existing_target_status": "metabolic_mechanism_and_metabolite_target_records_retrieved",
        "network_toxicology_status": "candidate_specific_network_toxicology_record_retrieved;_independent_validation_needed",
        "novelty_grade": "B_emerging_literature_not_saturated",
        "opportunity_bucket": "Top5_3_emerging_mechanistic_gap",
        "representative_pmids": "28095285;36480145;42074192",
        "largest_gap": "The exposure-to-T2D chain remains fragmented across pilot, gestational, and metabolic-disease studies; candidate-specific prospective T2D evidence is not yet established.",
        "decision_reason": "Advance as an emerging opportunity with a promising but not yet complete T2D-specific mechanism.",
    },
    "tin": {
        "direct_t2d_relevance": "limited_elemental_tin_specificity",
        "human_epidemiology_status": "human_organotin_and_urinary_tin_records;_elemental_tin_vs_organotin_not_resolved",
        "prospective_status": "no_clear_environmental_elemental_tin_T2D_prospective_record_retrieved",
        "mechanism_maturity": "organotin_metabolic_mechanisms_retrieved_but_not_transferable_to_elemental_tin_assay",
        "existing_target_status": "organotin_target_records_retrieved",
        "network_toxicology_status": "organotin_network_records_retrieved;_not_elemental_tin_T2D_specific",
        "novelty_grade": "B_high_species_resolution_gap",
        "opportunity_bucket": "Top5_4_speciation_gap",
        "representative_pmids": "15993011;30458368;31306684;31057082",
        "largest_gap": "NHANES urinary tin does not by itself establish the organotin species responsible for the biological records; chemical-speciation linkage is the central limitation.",
        "decision_reason": "Advance only as a clearly labeled speciation-limited opportunity; do not merge elemental tin and organotin evidence.",
    },
    "mibp": {
        "direct_t2d_relevance": "phthalate_metabolite_human_signal_with_family_level_overlap",
        "human_epidemiology_status": "multiple_human_phthalate_metabolite_and_diabetes_or_insulin_resistance_records",
        "prospective_status": "early_life_or_longitudinal_phthalate_records_retrieved; MiBP_specificity_varies",
        "mechanism_maturity": "limited_candidate_specific_mechanism; broader_phthalate_beta_cell_and_receptor_records",
        "existing_target_status": "phthalate_monoester_CAR_binding_record_retrieved",
        "network_toxicology_status": "mostly_family_level_or_non_target_network_records",
        "novelty_grade": "B_moderate_novelty_family_collision",
        "opportunity_bucket": "Top5_5_actionable_but_family_crowded",
        "representative_pmids": "22498808;28898934;25938866;36527833",
        "largest_gap": "MiBP is measurable and has human metabolic literature, but much of the evidence is phthalate-family or mixture-level rather than a complete MiBP-specific T2D mechanism.",
        "decision_reason": "Advance as a tractable secondary opportunity, below the sparse-collision candidates because phthalate-family literature is crowded.",
    },
    "pfhxs": {
        "direct_t2d_relevance": "PFAS_class_level_human_signal",
        "human_epidemiology_status": "multiple_PFAS_and_T2D_or_gestational_diabetes_records",
        "prospective_status": "prospective_PFAS_metabolic_records_retrieved",
        "mechanism_maturity": "candidate_specific_non_T2D_or_animal_mechanistic_records",
        "existing_target_status": "PPAR_and_hepatotoxicity_records_retrieved",
        "network_toxicology_status": "PFAS_network_toxicology_records_retrieved",
        "novelty_grade": "C_moderate_to_low_novelty_crowded_class",
        "opportunity_bucket": "monitor_not_Top5_class_crowded",
        "representative_pmids": "33984575;39680074;39286118;40865220",
        "largest_gap": "PFHxS-specific causal and T2D-specific mechanistic evidence is less developed than the broader PFAS literature.",
        "decision_reason": "Keep as a comparator/monitor candidate; prioritize sparse-collision candidates first.",
    },
    "uranium": {
        "direct_t2d_relevance": "direct_metal_epidemiology_with_growing_prospective_literature",
        "human_epidemiology_status": "direct_urinary_uranium_and_diabetes_records_retrieved",
        "prospective_status": "multiple_prospective_screening_records_including_incident_diabetes_record",
        "mechanism_maturity": "limited_direct_T2D_mechanism; broader_metal_or_toxicity_mechanisms",
        "existing_target_status": "no_candidate_specific_T2D_target_confirmed",
        "network_toxicology_status": "mostly_non_T2D_or_broad_toxicity_network_records",
        "novelty_grade": "C_existing_human_and_prospective_collision",
        "opportunity_bucket": "exclude_from_Top5_mature_epidemiology",
        "representative_pmids": "26542316;42486290;41734627",
        "largest_gap": "Mechanistic specificity remains incomplete, but the human/prospective collision burden lowers its novelty as a new opportunity.",
        "decision_reason": "Do not prioritize as a novelty opportunity; retain as a positive-control-like comparator for the audit.",
    },
    "molybdenum": {
        "direct_t2d_relevance": "established_metal_and_metabolic_literature",
        "human_epidemiology_status": "multiple_human_urinary_metal_and_glucose_records",
        "prospective_status": "multiple_prospective_or_longitudinal_screening_records",
        "mechanism_maturity": "broad_and_mixed_essential_biology_and_metabolic_mechanisms",
        "existing_target_status": "multiple_biochemical_or_metabolic_target_records",
        "network_toxicology_status": "high_collision_and_many_non_T2D_records",
        "novelty_grade": "C_crowded_and_biologically_mixed",
        "opportunity_bucket": "exclude_from_Top5_crowded",
        "representative_pmids": "26542316;32927284;36182150;39159974",
        "largest_gap": "The major limitation is exposure interpretation across essential molybdenum biology, mixtures, and diverse metabolic endpoints rather than a simple unexplored gap.",
        "decision_reason": "Do not prioritize over the sparse-collision candidates.",
    },
    "tungsten": {
        "direct_t2d_relevance": "direct_metal_signal_but_tungstate_therapeutic_collision",
        "human_epidemiology_status": "human_urinary_tungsten_and_diabetes_records_retrieved",
        "prospective_status": "multiple_prospective_screening_hits; environmental_vs_therapeutic_context_mixed",
        "mechanism_maturity": "many_tungstate_therapeutic_or_animal_records; environmental_mechanism_unclear",
        "existing_target_status": "tungstate_metabolic_records_retrieved",
        "network_toxicology_status": "mixed_non_T2D_network_records",
        "novelty_grade": "C_therapeutic_collision_and_mixed_context",
        "opportunity_bucket": "exclude_from_Top5_context_mixed",
        "representative_pmids": "26542316;34909553;30132852;24253047",
        "largest_gap": "Separating environmental tungsten exposure from sodium tungstate pharmacology is essential and currently limits a clean T2D opportunity claim.",
        "decision_reason": "Do not prioritize until exposure-form specificity is resolved.",
    },
    "barium": {
        "direct_t2d_relevance": "metal_mixture_and_barium_epidemiology_collision",
        "human_epidemiology_status": "multiple_human_metal_and_diabetes_records_barium_specificity_varies",
        "prospective_status": "multiple_prospective_or_follow_up_records_retrieved",
        "mechanism_maturity": "mostly_therapeutic_or_non_environmental_barium_records",
        "existing_target_status": "no_clean_candidate_specific_T2D_target_confirmed",
        "network_toxicology_status": "low_specificity_non_T2D_records",
        "novelty_grade": "C_crowded_mixture_context",
        "opportunity_bucket": "exclude_from_Top5_crowded",
        "representative_pmids": "26542316;32927284;38772737;41161236",
        "largest_gap": "Barium-specific exposure and mechanism remain difficult to separate from metal mixtures and therapeutic barium formulations.",
        "decision_reason": "Retain as an audit comparator rather than a leading opportunity.",
    },
    "lead": {
        "direct_t2d_relevance": "well_studied_environmental_metal_signal",
        "human_epidemiology_status": "substantial_human_environmental_lead_and_metabolic_literature",
        "prospective_status": "multiple_prospective_or_longitudinal_screening_records",
        "mechanism_maturity": "broad_environmental_toxicology_and_metabolic_mechanisms",
        "existing_target_status": "many_existing_toxicology_targets_and_pathways",
        "network_toxicology_status": "high_collision",
        "novelty_grade": "C_low_novelty_crowded",
        "opportunity_bucket": "exclude_from_Top5_crowded",
        "representative_pmids": "30200608;29460222;28936425",
        "largest_gap": "Lead is an established environmental toxicology topic; the remaining gap is refinement, not a clean unexplored exposure opportunity.",
        "decision_reason": "Exclude from the opportunity shortlist because literature collision is high.",
    },
    "silver": {
        "direct_t2d_relevance": "nanoparticle_and_therapeutic_literature_dominated",
        "human_epidemiology_status": "records_are_mostly_silver_nanoparticle_or_diabetic_wound_context",
        "prospective_status": "prospective_records_are_primarily_therapeutic_nanoparticle_studies",
        "mechanism_maturity": "nanoparticle_biosafety_and_wound_therapy_mechanisms",
        "existing_target_status": "many_nanoparticle_targets",
        "network_toxicology_status": "high_nanoparticle_network_collision",
        "novelty_grade": "C_identity_context_mismatch",
        "opportunity_bucket": "exclude_assay_context_mismatch",
        "representative_pmids": "32802176;32409075;35138974",
        "largest_gap": "The retrieved literature does not cleanly represent environmental elemental silver exposure measured by the NHANES analyte.",
        "decision_reason": "Exclude from the primary opportunity shortlist because the literature context does not match the exposure identity.",
    },
    "dehp": {
        "direct_t2d_relevance": "well_studied_phthalate_diabetes_axis",
        "human_epidemiology_status": "substantial_phthalate_and_diabetes_or_metabolic_literature",
        "prospective_status": "multiple_prospective_or_review_records",
        "mechanism_maturity": "mature_DEHP_insulin_and_metabolic_mechanisms",
        "existing_target_status": "multiple_existing_targets_and_pathways",
        "network_toxicology_status": "very_high_collision",
        "novelty_grade": "C_low_novelty_crowded",
        "opportunity_bucket": "exclude_from_Top5_crowded",
        "representative_pmids": "24130215;37367903;32549860;41349316",
        "largest_gap": "The main opportunity gap is not discovery of a new association but resolution of mixtures, metabolites, and causal specificity.",
        "decision_reason": "Exclude from a novelty-led Top 5; retain as a crowded-class benchmark.",
    },
    "meohp": {
        "direct_t2d_relevance": "DEHP_metabolite_and_phthalate_family_signal",
        "human_epidemiology_status": "human_metabolic_records_retrieved; metabolite_specificity_varies",
        "prospective_status": "screening_hits_are_mixed_and_not_all_T2D_specific",
        "mechanism_maturity": "PPAR_and_oxidative_stress_records_exist_but_family_context_dominates",
        "existing_target_status": "PPAR_gamma_and_other_phthalate_receptor_records",
        "network_toxicology_status": "family_level_network_collision",
        "novelty_grade": "C_phthalate_family_crowded",
        "opportunity_bucket": "exclude_from_Top5_family_crowded",
        "representative_pmids": "23977034;30085373;28095285",
        "largest_gap": "A clean MEOHP-specific exposure-to-T2D mechanism is difficult to distinguish from the larger DEHP metabolite family.",
        "decision_reason": "Do not prioritize over MECPP or sparse-collision candidates.",
    },
}

# Internal-collision reference: the existing project is the DINP-related
# exposure axis measured with urinary MCOP/URXCOP and developed around a
# phthalate/PPAR-style mechanism. This is a design constraint for a new T2D
# paper, not a new biological result. Penalties are prioritization aids, not
# inferential statistics.
INTERNAL_COLLISION = {
    "mcop": {
        "internal_exposure_overlap": "exact_DINP_related_axis",
        "internal_biomarker_overlap": "exact_URXCOP_MCOP",
        "internal_mechanism_overlap": "high_shared_phthalate_PPAR_architecture",
        "internal_collision_grade": "high_exact_reuse",
        "internal_collision_penalty": 6.0,
        "self_overlap_reason": "Same DINP-related exposure axis and same urinary MCOP/URXCOP biomarker as the existing CRC project; the downstream phthalate/PPAR framing would also substantially overlap.",
        "new_paper_position": "holdout_or_secondary_only",
    },
    "dinp_parent": {
        "internal_exposure_overlap": "exact_DINP_related_axis",
        "internal_biomarker_overlap": "exact_URXCOP_MCOP_axis",
        "internal_mechanism_overlap": "high_shared_phthalate_PPAR_architecture",
        "internal_collision_grade": "high_exact_axis_reuse",
        "internal_collision_penalty": 6.0,
        "self_overlap_reason": "Same parent exposure axis and URXCOP/MCOP biomarker family as the existing CRC project, with closely related phthalate mechanism framing.",
        "new_paper_position": "holdout_or_secondary_only",
    },
    "mecpp": {
        "internal_exposure_overlap": "related_phthalate_class",
        "internal_biomarker_overlap": "different_URXECP_biomarker",
        "internal_mechanism_overlap": "moderate_shared_phthalate_PPAR_architecture",
        "internal_collision_grade": "moderate_related_class",
        "internal_collision_penalty": 2.5,
        "self_overlap_reason": "Different metabolite, but still a phthalate exposure and likely to reuse the same broad PPAR/obesogen mechanistic language.",
        "new_paper_position": "eligible_with_explicit_differentiation",
    },
    "meohp": {
        "internal_exposure_overlap": "related_DEHP_phthalate_class",
        "internal_biomarker_overlap": "different_URXMOH_biomarker",
        "internal_mechanism_overlap": "moderate_shared_phthalate_PPAR_architecture",
        "internal_collision_grade": "moderate_related_class",
        "internal_collision_penalty": 2.5,
        "self_overlap_reason": "Different metabolite but same broader phthalate/PPAR narrative space; DEHP-family collision is already substantial.",
        "new_paper_position": "eligible_with_explicit_differentiation",
    },
    "mibp": {
        "internal_exposure_overlap": "related_phthalate_class",
        "internal_biomarker_overlap": "different_URXMIB_biomarker",
        "internal_mechanism_overlap": "moderate_shared_phthalate_PPAR_architecture",
        "internal_collision_grade": "moderate_related_class",
        "internal_collision_penalty": 2.0,
        "self_overlap_reason": "Different phthalate parent/metabolite and biomarker, but family-level exposure and receptor biology overlap with the prior project.",
        "new_paper_position": "eligible_with_explicit_differentiation",
    },
    "dehp": {
        "internal_exposure_overlap": "related_phthalate_class",
        "internal_biomarker_overlap": "related_DEHP_metabolite_panel",
        "internal_mechanism_overlap": "high_shared_phthalate_PPAR_architecture",
        "internal_collision_grade": "high_related_class_and_crowded",
        "internal_collision_penalty": 4.0,
        "self_overlap_reason": "Same phthalate family and overlapping metabolite/PPAR architecture, with an already crowded external literature.",
        "new_paper_position": "not_suitable_for_novelty_led_paper",
    },
    "pfhxs": {
        "internal_exposure_overlap": "distinct_PFAS_exposure",
        "internal_biomarker_overlap": "distinct_LBXPFHS_biomarker",
        "internal_mechanism_overlap": "low_partial_PPAR_overlap_only",
        "internal_collision_grade": "low_distinct_exposure",
        "internal_collision_penalty": 0.5,
        "self_overlap_reason": "Different chemical class and biomarker; only a limited conceptual overlap exists because PPAR signaling appears in some PFHxS toxicology records.",
        "new_paper_position": "eligible_distinct_axis",
    },
}

# External-opportunity baseline used only to make the re-ranking reproducible.
# It encodes the first-pass audit's qualitative novelty/actionability judgment,
# not a p-value or a biological score.
EXTERNAL_OPPORTUNITY_SCORE = {
    "mcop": 10.0,
    "dinp_parent": 9.0,
    "mecpp": 8.0,
    "tin": 7.5,
    "mibp": 7.0,
    "meohp": 6.0,
    "pfhxs": 6.5,
    "uranium": 5.5,
    "barium": 4.5,
    "molybdenum": 4.0,
    "tungsten": 4.0,
    "lead": 2.0,
    "silver": 1.0,
    "dehp": 2.0,
}

DEFAULT_INTERNAL_COLLISION = {
    "internal_exposure_overlap": "distinct_exposure",
    "internal_biomarker_overlap": "distinct_biomarker",
    "internal_mechanism_overlap": "distinct_or_not_comparable",
    "internal_collision_grade": "none_distinct_axis",
    "internal_collision_penalty": 0.0,
    "self_overlap_reason": "No material exposure, biomarker, or mechanism-architecture overlap with the existing DINP/MCOP–CRC project was identified in this first-pass design audit.",
    "new_paper_position": "eligible_distinct_axis",
}


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[4])
    args = parser.parse_args()
    repo = args.repo.resolve()
    root = repo / "analysis" / "disease_agnostic_environmental_framework" / "t2d_exposure_opportunity"
    lit_dir = root / "03_literature_collision"
    candidates = read_csv(root / "01_candidate_master" / "unique_candidate_chemicals.csv")
    counts = {row["chemical_id"]: row for row in read_csv(lit_dir / "literature_counts.csv")}
    if set(ADJUDICATION) != {row["search_group"] for row in counts.values()}:
        raise RuntimeError("Adjudication dictionary and literature search groups do not match exactly")

    rows: List[Dict[str, object]] = []
    for candidate in candidates:
        if not str(candidate.get("mapping_gate_disposition", "")).startswith("advance"):
            continue
        chemical_id = candidate["chemical_id"]
        count = counts[chemical_id]
        group = count["search_group"]
        manual = ADJUDICATION[group]
        internal = INTERNAL_COLLISION.get(group, DEFAULT_INTERNAL_COLLISION)
        external_score = EXTERNAL_OPPORTUNITY_SCORE[group]
        adjusted_score = external_score - float(internal["internal_collision_penalty"])
        rows.append(
            {
                "chemical_id": chemical_id,
                "chemical_name": candidate["chemical_name"],
                "chemical_class": candidate["chemical_class"],
                "positive_biomarkers": candidate["positive_biomarkers"],
                "exposure_axes": candidate.get("exposure_axes", ""),
                "search_group": group,
                "mapping_grades": candidate["mapping_grades"],
                "mapping_gate_disposition": candidate["mapping_gate_disposition"],
                "diabetes_total_pubmed_hits": count["diabetes_total_pubmed_hits"],
                "human_epidemiology_pubmed_hits": count["human_epidemiology_pubmed_hits"],
                "prospective_pubmed_hits": count["prospective_pubmed_hits"],
                "mechanism_pubmed_hits": count["mechanism_pubmed_hits"],
                "target_pubmed_hits": count["target_pubmed_hits"],
                "network_bioinformatics_pubmed_hits": count["network_bioinformatics_pubmed_hits"],
                "direct_t2d_relevance": manual["direct_t2d_relevance"],
                "human_epidemiology_status": manual["human_epidemiology_status"],
                "prospective_status": manual["prospective_status"],
                "mechanism_maturity": manual["mechanism_maturity"],
                "existing_target_status": manual["existing_target_status"],
                "network_toxicology_status": manual["network_toxicology_status"],
                "novelty_grade": manual["novelty_grade"],
                "opportunity_bucket": manual["opportunity_bucket"],
                "external_opportunity_score": f"{external_score:.1f}",
                "internal_exposure_overlap": internal["internal_exposure_overlap"],
                "internal_biomarker_overlap": internal["internal_biomarker_overlap"],
                "internal_mechanism_overlap": internal["internal_mechanism_overlap"],
                "internal_collision_grade": internal["internal_collision_grade"],
                "internal_collision_penalty": f"{float(internal['internal_collision_penalty']):.1f}",
                "revised_opportunity_score": f"{adjusted_score:.1f}",
                "new_paper_position": internal["new_paper_position"],
                "self_overlap_reason": internal["self_overlap_reason"],
                "representative_pmids": manual["representative_pmids"],
                "largest_gap": manual["largest_gap"],
                "decision_reason": manual["decision_reason"],
                "screening_basis": "top_relevance_sorted_PubMed_records_plus_category_counts; not_full_text_adjudicated",
            }
        )

    ranked_groups = sorted(
        ADJUDICATION,
        key=lambda group: (-float(next(row["revised_opportunity_score"] for row in rows if row["search_group"] == group)), group),
    )
    revised_top5_groups = ranked_groups[:5]
    group_rank = {group: str(rank) for rank, group in enumerate(revised_top5_groups, start=1)}
    for row in rows:
        row["revised_opportunity_rank"] = group_rank.get(str(row["search_group"]), "")
    rows.sort(key=lambda row: (0 if row["revised_opportunity_rank"] else 1, row["revised_opportunity_rank"] or "99", row["chemical_id"]))
    fields = list(rows[0].keys()) if rows else []
    adjudication_path = lit_dir / "candidate_evidence_adjudication.csv"
    write_csv(adjudication_path, fields, rows)

    top5_chemical_rows = [row for row in rows if row["search_group"] in revised_top5_groups]
    grouped_rows: Dict[str, List[Dict[str, object]]] = {}
    for row in top5_chemical_rows:
        grouped_rows.setdefault(str(row["search_group"]), []).append(row)
    top5 = []
    for group, group_rows in sorted(grouped_rows.items(), key=lambda item: int(str(item[1][0]["revised_opportunity_rank"]))):
        first = group_rows[0]
        top5.append(
            {
                "revised_opportunity_rank": first["revised_opportunity_rank"],
                "chemical_ids": ";".join(str(row["chemical_id"]) for row in group_rows),
                "chemical_names": "; ".join(str(row["chemical_name"]) for row in group_rows),
                "positive_biomarkers": ";".join(sorted({str(row["positive_biomarkers"]) for row in group_rows})),
                "exposure_axes": ";".join(sorted({str(row["exposure_axes"]) for row in group_rows})),
                "search_group": group,
                "mapping_grades": ";".join(sorted({str(row["mapping_grades"]) for row in group_rows})),
                "diabetes_total_pubmed_hits": first["diabetes_total_pubmed_hits"],
                "human_epidemiology_pubmed_hits": first["human_epidemiology_pubmed_hits"],
                "prospective_pubmed_hits": first["prospective_pubmed_hits"],
                "mechanism_maturity": first["mechanism_maturity"],
                "novelty_grade": first["novelty_grade"],
                "opportunity_bucket": first["opportunity_bucket"],
                "external_opportunity_score": first["external_opportunity_score"],
                "internal_exposure_overlap": first["internal_exposure_overlap"],
                "internal_biomarker_overlap": first["internal_biomarker_overlap"],
                "internal_mechanism_overlap": first["internal_mechanism_overlap"],
                "internal_collision_grade": first["internal_collision_grade"],
                "internal_collision_penalty": first["internal_collision_penalty"],
                "revised_opportunity_score": first["revised_opportunity_score"],
                "revised_opportunity_rank": first["revised_opportunity_rank"],
                "new_paper_position": first["new_paper_position"],
                "self_overlap_reason": first["self_overlap_reason"],
                "representative_pmids": first["representative_pmids"],
                "largest_gap": first["largest_gap"],
                "decision_reason": first["decision_reason"],
            }
        )
    top5_fields = [
        "revised_opportunity_rank", "chemical_ids", "chemical_names", "positive_biomarkers", "exposure_axes",
        "search_group", "mapping_grades", "diabetes_total_pubmed_hits", "human_epidemiology_pubmed_hits",
        "prospective_pubmed_hits", "mechanism_maturity", "novelty_grade", "opportunity_bucket",
        "external_opportunity_score", "internal_exposure_overlap", "internal_biomarker_overlap",
        "internal_mechanism_overlap", "internal_collision_grade", "internal_collision_penalty",
        "revised_opportunity_score", "new_paper_position", "self_overlap_reason",
        "representative_pmids", "largest_gap", "decision_reason",
    ]
    shortlist_path = lit_dir / "PHASE2_PROVISIONAL_TOP5_OPPORTUNITIES.csv"
    write_csv(shortlist_path, top5_fields, top5)

    internal_fields = [
        "revised_opportunity_rank", "chemical_id", "chemical_name", "chemical_class", "positive_biomarkers",
        "exposure_axes", "search_group", "mapping_grades", "external_opportunity_score",
        "internal_exposure_overlap", "internal_biomarker_overlap", "internal_mechanism_overlap",
        "internal_collision_grade", "internal_collision_penalty", "revised_opportunity_score",
        "new_paper_position", "self_overlap_reason",
    ]
    internal_audit_path = lit_dir / "PHASE2_INTERNAL_COLLISION_AUDIT.csv"
    write_csv(internal_audit_path, internal_fields, rows)

    search_date = max(row.get("search_date", "") for row in counts.values())
    report_lines = [
        "# Phase 2 — Literature collision, novelty, and opportunity audit",
        "",
        "## Status",
        "",
        "- **Status:** `complete_first_pass_internal_collision_reranking`",
        f"- **Audit date (UTC):** `{datetime.now(timezone.utc).date().isoformat()}`",
        f"- **PubMed search date:** `{search_date}`",
        "- **Scope:** all 15 chemical IDs in the Phase 1 advance pool, represented by 14 search groups because the two DINP parent IDs share one parent-specific vocabulary.",
        "- **Methods:** reproducible PubMed category searches, retrieval of top relevance-sorted records, and a transparent first-pass title/abstract-level candidate adjudication.",
        "",
        "## Interpretation boundary",
        "",
        "The category counts are collision-screening signals, not eligible-study counts. The first-pass adjudication checks whether retrieved titles support candidate identity, human/T2D relevance, prospective design, mechanism maturity, target status, and network-toxicology collision. It is not a full-text systematic review and does not establish causality.",
        "",
        "## Internal project-collision audit",
        "",
        "The existing DINP/MCOP–CRC project is treated as an internal-collision reference. Each candidate is annotated for overlap in exposure axis, urinary biomarker, and mechanism architecture. Exact reuse receives the largest penalty; related phthalate/PPAR framing receives an intermediate penalty; distinct chemical classes receive zero or minimal penalty. This is a transparent paper-design constraint, not a biological or inferential statistic.",
        "",
        "The external-opportunity baseline and internal-collision penalty are embedded in `build_phase2_candidate_adjudication.py`; they are ordinal prioritization aids and must not be interpreted as effect sizes, probabilities, or evidence scores.",
        "",
        "## Revised Top 5 opportunity pool",
        "",
        "The revised Top 5 is an opportunity shortlist after applying the internal-collision constraint. It prioritizes a combination of mapping actionability, external literature headroom, and low self-overlap. It is not a claim that these candidates have the strongest causal evidence, and ranking may change after full-text eligibility review.",
        "",
        "| Rank | Exposure opportunity | Mapping | External score | Self-overlap penalty | Revised score | Position |",
        "|---:|---|---|---:|---:|---:|---|",
    ]
    for row in top5:
        report_lines.append(
            f"| {row['revised_opportunity_rank']} | {row['chemical_names']} ({row['positive_biomarkers']}) | {row['mapping_grades']} | {row['external_opportunity_score']} | {row['internal_collision_penalty']} | {row['revised_opportunity_score']} | {row['new_paper_position']} |"
        )
    report_lines.extend(
        [
            "",
            "## Candidate-level disposition",
            "",
            "- `candidate_evidence_adjudication.csv` contains one row per Phase 1 advance-pool chemical ID, including all 15 IDs and the shared DINP parent search group.",
            "- `PHASE2_PROVISIONAL_TOP5_OPPORTUNITIES.csv` is the compact self-overlap-adjusted shortlist for the next full-text collision pass.",
            "- `PHASE2_INTERNAL_COLLISION_AUDIT.csv` exposes the overlap dimensions, penalty, revised score, and new-paper position for all 15 candidate IDs.",
            "- `candidate_evidence_adjudication.csv` contains the same internal overlap fields for every candidate ID; MCOP/DINP are retained as explicit high-overlap holdouts rather than silently deleted.",
            "- All non-Top-5 candidates remain in the adjudication table with an explicit reason for monitoring or exclusion from the novelty-led shortlist.",
            "",
            "## Main opportunity gaps",
            "",
            "1. MCOP and DINP remain externally sparse but are intentionally downranked because they overlap with the existing exposure, biomarker, and mechanistic architecture.",
            "2. Several phthalate metabolites remain related-class opportunities only if the new paper explicitly differentiates exposure, biomarker, and mechanism from the existing project.",
            "3. Tin requires chemical-speciation resolution: elemental urinary tin cannot be treated as interchangeable with organotin studies.",
            "4. Candidates with high external collision counts (lead, DEHP, silver nanoparticle literature, and established metal/PFAS axes) are useful comparators but are not prioritized as novelty opportunities.",
            "",
            "## Next required step",
            "",
            "Perform title/abstract adjudication followed by targeted full-text retrieval for the revised Top 5 and any candidate whose mapping identity or exposure form is ambiguous. Freeze eligibility criteria before reviewing outcome-specific details; do not run docking or expand target discovery in this phase.",
        ]
    )
    report_path = lit_dir / "PHASE2_CANDIDATE_ADJUDICATION_REPORT.md"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    generated = [adjudication_path, shortlist_path, internal_audit_path, report_path]
    manifest = {
        "material_id": "T2D-EXPOSURE-OPPORTUNITY-PHASE2-ADJUDICATION",
        "script_version": SCRIPT_VERSION,
        "status": "complete_first_pass_internal_collision_reranking",
        "audit_date_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "phase1_candidate_path": str((root / "01_candidate_master" / "unique_candidate_chemicals.csv").relative_to(repo)),
            "phase1_candidate_sha256": sha256_file(root / "01_candidate_master" / "unique_candidate_chemicals.csv"),
            "literature_counts_path": str((lit_dir / "literature_counts.csv").relative_to(repo)),
            "literature_counts_sha256": sha256_file(lit_dir / "literature_counts.csv"),
            "top_pubmed_records_path": str((lit_dir / "top_pubmed_records.csv").relative_to(repo)),
            "top_pubmed_records_sha256": sha256_file(lit_dir / "top_pubmed_records.csv"),
        },
        "candidate_chemical_count": len(rows),
        "search_group_count": len(ADJUDICATION),
        "revised_top5_group_count": len(top5),
        "revised_top5_search_groups": [row["search_group"] for row in top5],
        "revised_top5_chemical_ids": [row["chemical_ids"] for row in top5],
        "internal_collision_reference": "existing DINP/MCOP–CRC project: DINP-related exposure axis, urinary MCOP/URXCOP, phthalate/PPAR-style mechanism framing",
        "scoring_guardrail": "External baseline minus internal-collision penalty is ordinal prioritization only; it is not an inferential statistic.",
        "outputs": {str(path.relative_to(repo)): {"sha256": sha256_file(path), "bytes": path.stat().st_size} for path in generated},
        "not_performed": ["full_text_systematic_review", "formal_eligible_study_count", "causal_inference", "docking", "new_target_discovery", "experimental_feasibility"],
        "manual_adjudication_dictionary_version": "embedded_in_script_v1.0",
    }
    manifest_path = lit_dir / "PHASE2_ADJUDICATION_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": manifest["status"],
        "candidate_chemical_count": len(rows),
        "search_group_count": len(ADJUDICATION),
        "revised_top5_group_count": len(top5),
        "output_dir": str(lit_dir),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
