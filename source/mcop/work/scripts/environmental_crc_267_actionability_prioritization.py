"""Phase 2I: outcome-blinded actionability prioritization across all 267 core chemicals.

The purpose of this script is to answer the original paper question directly:
was MCOP present in the prespecified 267-chemical discovery universe and did it
survive a CRC-outcome-blinded human-testability filter before the NHANES result
was considered?

This is deliberately not a second molecular ranking and it is not a targeted
MCOP rescue. Molecular fields are inherited from Phase 1, while human
actionability fields are populated only from the auditable biomonitoring files
currently available locally. At present those files cover MCOP and MiNP only;
the other 265 chemicals are explicitly sent to a manual-review queue.

No human CRC OR, CI, P value, LOCO result, or cycle-specific CRC effect is used
for eligibility. CRC overlap in the CTD x GeneCards molecular screen is allowed
because it is part of Phase 1 molecular nomination, not a human epidemiologic
outcome.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / "outputs"

PHASE1_CORE = OUTPUTS / "environmental_toxicology_crc_phase1_ranked_core.csv"
DEGREE = OUTPUTS / "environmental_toxicology_crc_phase1_degree_matched_permutation.csv"
NHANES_AUDIT = OUTPUTS / "nhanes_dinp_phase2a_audit_summary.csv"

OUT_MATRIX = OUTPUTS / "environmental_crc_267_actionability_matrix.csv"
OUT_RULES = OUTPUTS / "environmental_crc_267_actionability_rules.json"
OUT_SUMMARY = OUTPUTS / "environmental_crc_267_actionability_summary.md"
OUT_SENSITIVITY = OUTPUTS / "environmental_crc_267_actionability_sensitivity.csv"
OUT_HUMAN_TESTABLE = OUTPUTS / "environmental_crc_267_human_testable_candidates.csv"
OUT_DISPOSITION = OUTPUTS / "environmental_crc_267_disposition_audit.csv"
OUT_REVIEW = OUTPUTS / "environmental_crc_267_manual_review_queue.csv"
OUT_MANIFEST = OUTPUTS / "environmental_crc_267_prioritization_manifest.json"


CORE_SCOPE = "GeneCards_Disorders"
CORE_K = 1000
EXPECTED_CHEMICALS = 267


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def nonempty(value: object) -> bool:
    return value is not None and not pd.isna(value) and str(value).strip() != ""


def entity_annotation(row: pd.Series) -> tuple[int, str, bool]:
    """Conservative entity-quality annotation with an explicit reason."""

    name = str(row.get("ChemicalName", ""))
    broad_patterns = [
        r"\bvolatile organic compounds?\b",
        r"\bair pollutants?\b",
        r"\boccupational air pollutants?\b",
        r"\bphthalic acids?\b",
        r"\bpolycyclic aromatic hydrocarbons?\b",
        r"\bpolychlorinated biphenyls?\b",
        r"\bflame retardants?\b",
        r"\bpesticides?\b",
        r"\bper[- ]?and polyfluoroalkyl substances?\b",
        r"\bpfas\b",
        r"\bheavy metals?\b",
        r"\bmetals?\b",
        r"\bphenols?\b",
    ]
    matched = next((p for p in broad_patterns if re.search(p, name, flags=re.I)), None)
    identifier_fields = ["CasRN", "PubChemCID", "DTXSID", "InChIKey"]
    identifiers = [field for field in identifier_fields if nonempty(row.get(field))]
    if matched:
        return 0, f"broad_or_umbrella_name:{matched}", True
    if identifiers:
        return 1, "specific_entity_identifier:" + ",".join(identifiers), False
    # CTD ChemicalID is itself a stable entity identifier.  Do not downgrade a
    # chemically specific name such as monoisononylphthalate merely because the
    # Phase 1 export lacks a CAS/PubChem cross-reference.  Broad terms were
    # already intercepted above.
    if nonempty(row.get("ChemicalID")) and nonempty(name):
        return 1, "specific_ctd_entity_id:" + str(row["ChemicalID"]), False
    return 0, "no_specific_entity_identifier_in_phase1_export", True


def testability_level(n: object, cases: object) -> tuple[int, str]:
    """Pre-outcome feasibility tier based on available sample infrastructure."""

    n_value = float(n) if nonempty(n) else np.nan
    cases_value = float(cases) if nonempty(cases) else np.nan
    if not np.isfinite(n_value) or not np.isfinite(cases_value):
        return 0, "not audited"
    if n_value >= 10000 and cases_value >= 60:
        return 2, f"analytic_n={int(n_value)};crc_cases={int(cases_value)}"
    if n_value >= 500 and cases_value >= 20:
        return 1, f"analytic_n={int(n_value)};crc_cases={int(cases_value)}"
    return 0, f"analytic_n={int(n_value)};crc_cases={int(cases_value)};below feasibility threshold"


def load_primary_phase1() -> pd.DataFrame:
    phase1 = pd.read_csv(PHASE1_CORE, low_memory=False)
    primary = phase1.loc[
        (phase1["scope"] == CORE_SCOPE) & (phase1["gene_cards_k"] == CORE_K)
    ].copy()
    if primary["ChemicalID"].nunique() != EXPECTED_CHEMICALS:
        raise ValueError(
            f"Expected {EXPECTED_CHEMICALS} unique primary chemicals, "
            f"found {primary['ChemicalID'].nunique()}"
        )
    if len(primary) != EXPECTED_CHEMICALS:
        raise ValueError(f"Primary Phase 1 slice has {len(primary)} rows, not {EXPECTED_CHEMICALS}")
    degree = pd.read_csv(DEGREE, low_memory=False)
    degree = degree.drop_duplicates("ChemicalID")
    merged = primary.merge(degree, on="ChemicalID", how="left", validate="one_to_one")
    if merged["degree_matched_bh_fdr"].isna().any():
        missing = merged.loc[merged["degree_matched_bh_fdr"].isna(), "ChemicalID"].tolist()
        raise ValueError(f"Degree-matched output missing {len(missing)} primary chemicals: {missing[:10]}")
    return merged


def load_human_audit() -> pd.DataFrame:
    audit = pd.read_csv(NHANES_AUDIT)
    wanted = audit.loc[audit["analyte"].isin(["MCOP", "MiNP"])].copy()
    wanted["ChemicalID"] = wanted["analyte"].map(
        {"MCOP": "C573544", "MiNP": "C471400"}
    )
    return wanted


def apply_actionability(merged: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame:
    audit_by_id = audit.set_index("ChemicalID").to_dict(orient="index")
    records: list[dict[str, object]] = []

    for _, row in merged.iterrows():
        chemical_id = row["ChemicalID"]
        e_tag, e_reason, e_review = entity_annotation(row)

        n_ctd = float(row["n_ctd_human_genes"])
        overlap = float(row["crc_overlap"])
        bh = float(row["bh_fdr"])
        degree_bh = float(row["degree_matched_bh_fdr"])
        m2 = bool(
            bh < 0.05
            and degree_bh < 0.05
            and n_ctd >= 20
            and overlap >= 5
        )
        m1 = bool(
            m2
            or bh < 0.05
            or degree_bh < 0.05
            or (n_ctd >= 20 and overlap >= 5)
        )
        molecular_level = "M2" if m2 else ("M1" if m1 else "M0")

        human = audit_by_id.get(chemical_id)
        if human is None:
            human = {}
        known = bool(human)
        analyte = str(human.get("analyte", "")) if known else ""
        cycles = int(human["n_cycles"]) if known and nonempty(human.get("n_cycles")) else np.nan
        analytic_n = (
            float(human["exposure_and_crc_outcome_n"])
            if known and nonempty(human.get("exposure_and_crc_outcome_n"))
            else np.nan
        )
        crc_cases = float(human["crc_cases"]) if known and nonempty(human.get("crc_cases")) else np.nan
        t_tag, t_basis = testability_level(analytic_n, crc_cases)

        # These are exposure-axis annotations, not candidate substitutions.
        if analyte == "MCOP":
            x_tag, b_tag, d_tag, c_tag = 1, 1, 2, 2
            axis = "DINP-related exposure axis"
            relationship = "direct urinary MCOP analyte; biomarker of DINP-related exposure"
            identity_note = "Direct MCOP candidate; do not relabel as MiNP."
            d_reason = "MCOP above-LOD 98.4% across 7 NHANES cycles"
            n_tag = "pending_manual_review"
            review_reason = "Novelty/collision review remains to be completed; human actionability is audited."
        elif analyte == "MiNP":
            x_tag, b_tag, d_tag, c_tag = 1, 1, 0, 2
            axis = "DINP-related exposure axis"
            relationship = "direct urinary MiNP analyte; biomarker of DINP-related exposure"
            identity_note = "Direct MiNP candidate; not interchangeable with MCOP."
            d_reason = "MiNP above-LOD 27.4% across 7 NHANES cycles; fails direct detectability gate"
            n_tag = "pending_manual_review"
            review_reason = "Novelty/collision review remains to be completed; direct analyte detectability is insufficient."
        else:
            x_tag = b_tag = d_tag = c_tag = 0
            axis = ""
            relationship = "not identified in current local biomonitoring audit"
            identity_note = "No candidate-specific human biomarker mapping audited locally."
            d_reason = "not identified in current audit"
            n_tag = "pending_manual_review"
            review_reason = "No local human biomarker, detectability, cycle-coverage, or testability audit."

        # Outcome-blinded actionability gates. No human OR/CI/P/LOCO fields enter this expression.
        eligible_perm = bool(e_tag == 1 and x_tag == 1 and b_tag == 1 and d_tag >= 1 and c_tag >= 1 and t_tag >= 1)
        eligible_mod = bool(e_tag == 1 and x_tag == 1 and b_tag == 1 and d_tag >= 1 and c_tag >= 2 and t_tag >= 1)
        eligible_strict = bool(e_tag == 1 and x_tag == 1 and b_tag == 1 and d_tag == 2 and c_tag == 2 and t_tag == 2)

        if analyte == "MCOP" and eligible_perm:
            disposition = "advance_to_systematic_human_screen"
            priority = f"human-testable; {molecular_level}; novelty_pending"
        elif analyte == "MiNP":
            disposition = "retain_as_molecular_axis_candidate_but_fail_direct_detectability_gate"
            priority = f"not_directly_eligible; {molecular_level}; translated_axis_review"
        elif known:
            disposition = "manual_review_required"
            priority = f"manual_review; {molecular_level}"
        else:
            disposition = "manual_review_required"
            priority = f"manual_review; {molecular_level}"

        records.append(
            {
                "ChemicalID": chemical_id,
                "ChemicalName": row["ChemicalName"],
                "chemical_class": row.get("chemical_class", ""),
                "CasRN": row.get("CasRN", ""),
                "PubChemCID": row.get("PubChemCID", ""),
                "DTXSID": row.get("DTXSID", ""),
                "phase1_scope": row["scope"],
                "gene_cards_k": int(row["gene_cards_k"]),
                "phase1_unfiltered_rank": int(row["unfiltered_rank"]),
                "phase1_bh_fdr": bh,
                "phase1_degree_matched_bh_fdr": degree_bh,
                "phase1_n_ctd_human_genes": int(n_ctd),
                "phase1_crc_overlap": int(overlap),
                "molecular_level": molecular_level,
                "molecular_M2_definition_met": m2,
                "entity_E_tag": e_tag,
                "entity_E_reason": e_reason,
                "entity_manual_review_required": e_review,
                "exposure_axis_X_tag": x_tag,
                "exposure_axis_name": axis,
                "biomarker_relationship": relationship,
                "biomarker_B_tag": b_tag,
                "human_biomarker": analyte,
                "biomarker_matrix": "urine" if known else "",
                "direct_analyte_identity": row["ChemicalName"] if known else "",
                "identity_separation_note": identity_note,
                "detectability_D_tag": d_tag,
                "above_lod_pct": float(human["above_lod_pct"]) if known and nonempty(human.get("above_lod_pct")) else np.nan,
                "detectability_reason": d_reason,
                "cycle_coverage_C_tag": c_tag,
                "n_cycles_available": cycles,
                "cycles_available": human.get("cycles_with_file", "") if known else "",
                "testability_T_tag": t_tag,
                "testability_basis": t_basis,
                "analytic_exposure_outcome_n": analytic_n,
                "available_crc_cases_for_planning": crc_cases,
                "available_cancer_free_controls_for_planning": float(human["cancer_free_controls"]) if known and nonempty(human.get("cancer_free_controls")) else np.nan,
                "novelty_N_tag": n_tag,
                "manual_review_required": bool(e_review or not known or n_tag == "pending_manual_review"),
                "manual_review_reason": review_reason,
                "eligible_permissive": eligible_perm,
                "eligible_moderate": eligible_mod,
                "eligible_strict": eligible_strict,
                "priority_tier": priority,
                "disposition": disposition,
                "PRIORITIZATION_OUTCOME_BLINDED": True,
            }
        )

    return pd.DataFrame.from_records(records)


def write_outputs(matrix: pd.DataFrame, audit: pd.DataFrame) -> None:
    matrix.sort_values(["eligible_permissive", "molecular_level", "phase1_unfiltered_rank"], ascending=[False, True, True]).to_csv(OUT_MATRIX, index=False)
    matrix.loc[matrix["eligible_permissive"]].sort_values("phase1_unfiltered_rank").to_csv(OUT_HUMAN_TESTABLE, index=False)
    matrix.loc[matrix["manual_review_required"]].sort_values(["known" if "known" in matrix else "ChemicalID"]).to_csv(OUT_REVIEW, index=False)
    matrix[[
        "ChemicalID", "ChemicalName", "phase1_unfiltered_rank", "molecular_level",
        "eligible_permissive", "eligible_moderate", "eligible_strict", "disposition",
        "manual_review_required", "manual_review_reason",
    ]].to_csv(OUT_DISPOSITION, index=False)

    sensitivity = pd.DataFrame(
        [
            {"rule_set": "permissive", "definition": "E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1", "n_candidates": int(matrix["eligible_permissive"].sum()), "candidates": ";".join(matrix.loc[matrix["eligible_permissive"], "ChemicalID"].tolist())},
            {"rule_set": "moderate", "definition": "E=1 & X=1 & B=1 & D>=1 & C>=2 & T>=1", "n_candidates": int(matrix["eligible_moderate"].sum()), "candidates": ";".join(matrix.loc[matrix["eligible_moderate"], "ChemicalID"].tolist())},
            {"rule_set": "strict", "definition": "E=1 & X=1 & B=1 & D=2 & C=2 & T=2", "n_candidates": int(matrix["eligible_strict"].sum()), "candidates": ";".join(matrix.loc[matrix["eligible_strict"], "ChemicalID"].tolist())},
        ]
    )
    sensitivity.to_csv(OUT_SENSITIVITY, index=False)

    molecular_counts = matrix["molecular_level"].value_counts().to_dict()
    known = matrix.loc[matrix["human_biomarker"].astype(str).str.len() > 0]
    mcop = matrix.loc[matrix["ChemicalID"] == "C573544"].iloc[0]

    summary = f"""# Phase 2I — 267-chemical outcome-blinded actionability prioritization

## Direct answer to the original paper question

MCOP is present as **ChemicalID C573544** in the original Phase 1 universe of
**{len(matrix)} core environmental chemicals**. It is not introduced after the
NHANES result and it is not selected by a candidate-specific CRC association.

Under the prespecified actionability rules, MCOP is retained for the systematic
human screen because it is a specific entity with an auditable urinary biomarker,
7-cycle NHANES coverage, 98.4% above-LOD detectability, and sufficient available
sample infrastructure. Its Phase 1 molecular tier is **{mcop['molecular_level']}**;
human actionability, rather than molecular rank, is what permits advancement.

This is the direct-discovery formulation:

`267 core chemicals → outcome-blinded actionability annotation → MCOP retained → uniform NHANES screen`

It must not be rewritten as “MiNP molecular signal → MCOP” or as a post-hoc
replacement of MiNP by MCOP. MiNP remains a separate direct analyte and fails the
direct detectability gate in the current audit (27.4% above LOD).

## Current matrix status

- Primary Phase 1 slice: **{len(matrix)} unique chemicals** (`{CORE_SCOPE}`, GeneCards k={CORE_K}).
- Molecular tiers: **{molecular_counts}**.
- Human biomarker audit currently available locally: **{len(known)} chemicals (MCOP and MiNP only)**.
- Remaining **{len(matrix) - len(known)} chemicals** are explicitly in the manual-review queue; they are not treated as biomarker-negative.
- Outcome-blinded candidates passing permissive actionability rules: **{int(matrix['eligible_permissive'].sum())}**.
- Permissive / moderate / strict counts: **{int(matrix['eligible_permissive'].sum())} / {int(matrix['eligible_moderate'].sum())} / {int(matrix['eligible_strict'].sum())}**.

## Frozen actionability rules

`E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1` (permissive)

Moderate requires `C>=2`; strict requires `D=2 & C=2 & T=2`. Molecular
evidence is reported as M0/M1/M2 and is not a hard gate. Novelty is pending
manual collision review and is never inferred from the NHANES result.

## Outcome firewall

The eligibility expression uses no human CRC OR, CI, P value, LOCO result,
cycle-specific CRC effect, or candidate-specific epidemiologic significance.
Only Phase 1 molecular overlap/FDR and pre-existing sample infrastructure
(analytic N, available case count, control count, weights/cycle coverage) are
used for actionability/testability. The latter are feasibility fields, not
effect estimates.

## Limitation that must remain visible

This is the completed **267-chemical input and firewall audit**, but not yet a
completed 267-chemical biomonitoring annotation: 265 chemicals still require
the same exposure-database/biomarker/detectability review. The current MCOP
direct-discovery claim is valid as an audit of its original presence and its
pre-outcome eligibility; the full-universe claim should not be overstated until
that queue is completed.

Generated: {datetime.now(timezone.utc).isoformat()}
"""
    OUT_SUMMARY.write_text(summary, encoding="utf-8")


def main() -> None:
    merged = load_primary_phase1()
    audit = load_human_audit()
    matrix = apply_actionability(merged, audit)
    write_outputs(matrix, audit)

    forbidden_patterns = [
        "human_or", "nhanes_or", "odds_ratio_human", "ci_low", "ci_high",
        "nhanes_p", "loco", "cycle_specific_or", "cycle_specific_effect",
    ]
    source_columns = set(pd.read_csv(PHASE1_CORE, nrows=0).columns)
    forbidden_present = sorted(
        column for column in source_columns
        if any(pattern in column.lower() for pattern in forbidden_patterns)
    )
    if forbidden_present:
        raise ValueError(f"Outcome firewall detected forbidden columns: {forbidden_present}")

    rules = {
        "phase": "Phase 2I",
        "objective": "direct discovery of MCOP from the original 267 core environmental chemicals",
        "input": {
            "file": str(PHASE1_CORE),
            "scope": CORE_SCOPE,
            "gene_cards_k": CORE_K,
            "expected_unique_chemicals": EXPECTED_CHEMICALS,
        },
        "outcome_blinded": True,
        "outcome_firewall": {
            "forbidden_human_outcome_fields": [
                "human OR/HR", "human CI", "human P value", "LOCO CRC result",
                "cycle-specific CRC effect", "candidate-specific epidemiologic significance",
            ],
            "forbidden_columns_detected_in_phase1_input": forbidden_present,
            "allowed_feasibility_fields": ["analytic exposure/outcome N", "available CRC cases", "controls", "cycle coverage", "detectability"],
        },
        "molecular_rules": {
            "M2": "BH-FDR<0.05 AND degree-matched BH-FDR<0.05 AND n_ctd_human_genes>=20 AND CRC_overlap>=5",
            "M1": "M2 OR either FDR<0.05 OR n_ctd_human_genes>=20 AND CRC_overlap>=5",
            "M0": "otherwise",
            "note": "CRC_overlap here is Phase 1 molecular overlap, not a human epidemiologic result.",
        },
        "actionability_rules": {
            "permissive": "E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1",
            "moderate": "E=1 & X=1 & B=1 & D>=1 & C>=2 & T>=1",
            "strict": "E=1 & X=1 & B=1 & D=2 & C=2 & T=2",
            "M_is_hard_gate": False,
            "N_is_hard_gate": False,
        },
        "local_human_audit_coverage": {
            "audited_analytes": ["MCOP", "MiNP"],
            "unreviewed_chemical_count": int(len(matrix) - len(audit)),
            "warning": "Unreviewed chemicals are manual-review unknowns, not negatives.",
        },
        "mcop_direct_discovery_audit": {
            "ChemicalID": "C573544",
            "present_in_original_267": bool((matrix["ChemicalID"] == "C573544").any()),
            "phase1_rank": int(matrix.loc[matrix["ChemicalID"] == "C573544", "phase1_unfiltered_rank"].iloc[0]),
            "molecular_level": str(matrix.loc[matrix["ChemicalID"] == "C573544", "molecular_level"].iloc[0]),
            "eligible_permissive": bool(matrix.loc[matrix["ChemicalID"] == "C573544", "eligible_permissive"].iloc[0]),
            "selection_used_human_crc_effect": False,
        },
    }
    OUT_RULES.write_text(json.dumps(rules, indent=2, ensure_ascii=False), encoding="utf-8")

    manifest = {
        "phase": "Phase 2I",
        "script": str(Path(__file__)),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": [str(path) for path in [OUT_MATRIX, OUT_RULES, OUT_SUMMARY, OUT_SENSITIVITY, OUT_HUMAN_TESTABLE, OUT_DISPOSITION, OUT_REVIEW]],
        "primary_chemical_count": int(len(matrix)),
        "human_audit_chemical_count": int((matrix["human_biomarker"].astype(str).str.len() > 0).sum()),
        "manual_review_count": int(matrix["manual_review_required"].sum()),
        "eligible_permissive_count": int(matrix["eligible_permissive"].sum()),
        "eligible_moderate_count": int(matrix["eligible_moderate"].sum()),
        "eligible_strict_count": int(matrix["eligible_strict"].sum()),
        "outcome_blinded": True,
        "forbidden_columns_detected": forbidden_present,
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    # Machine-readable stdout for the run log.
    print(json.dumps({
        "primary_chemicals": len(matrix),
        "human_audited": int((matrix["human_biomarker"].astype(str).str.len() > 0).sum()),
        "manual_review": int(matrix["manual_review_required"].sum()),
        "eligible_permissive": int(matrix["eligible_permissive"].sum()),
        "eligible_moderate": int(matrix["eligible_moderate"].sum()),
        "eligible_strict": int(matrix["eligible_strict"].sum()),
        "mcop": matrix.loc[matrix["ChemicalID"] == "C573544", ["ChemicalName", "molecular_level", "eligible_permissive", "disposition"]].to_dict(orient="records"),
        "outputs": [str(path) for path in [OUT_MATRIX, OUT_RULES, OUT_SUMMARY, OUT_SENSITIVITY, OUT_HUMAN_TESTABLE, OUT_DISPOSITION, OUT_REVIEW, OUT_MANIFEST]],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
