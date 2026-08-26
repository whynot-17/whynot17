#!/usr/bin/env python3
"""Step 9D: synthesize the CRC negative branch into a diagnostic output.

This script only integrates the locked 9A/9B/9C outputs.  It does not refit
models, recompute multiplicity, or add new biological/exposure information.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

READINESS_OVERVIEW = OUT / "step9a_readiness_overview.csv"
READINESS_PER_TEST = OUT / "step9a_per_test_readiness.csv"
ATTRIBUTION_SUMMARY = OUT / "step9b_crc_failure_attribution_summary.csv"
COUNTERFACTUAL_SUMMARY = OUT / "step9c_counterfactual_summary.csv"
PRIORITY_MATRIX = OUT / "step9c_data_improvement_priority_matrix.csv"
STRUCTURAL = OUT / "step9c_structural_outcome_improvements.csv"

SYNTHESIS_OUT = OUT / "step9d_negative_branch_synthesis.csv"
SYMMETRY_OUT = OUT / "step9d_positive_negative_branch_symmetry.csv"
MANIFEST_OUT = OUT / "STEP9D_MANIFEST.json"
REPORT_OUT = OUT / "STEP9D_NEGATIVE_BRANCH_SYNTHESIS_REPORT.md"


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def metric(summary: pd.DataFrame, name: str):
    hit = summary.loc[summary["metric"].eq(name), "value"]
    if len(hit) != 1:
        raise ValueError(f"Missing or duplicated metric: {name}")
    return hit.iloc[0]


def build_synthesis(overview: pd.DataFrame, per_test: pd.DataFrame, attribution: pd.DataFrame, counterfactual: pd.DataFrame, priority: pd.DataFrame, structural: pd.DataFrame) -> pd.DataFrame:
    crc = overview.loc[overview["outcome"].eq("CRC")].iloc[0]
    t2d = overview.loc[overview["outcome"].eq("T2D")].iloc[0]
    attr = {row["metric"]: row["value"] for _, row in attribution.iterrows()}
    cf = {row["metric"]: row["value"] for _, row in counterfactual.iterrows()}
    crc_retention = float(per_test.loc[per_test["outcome"].eq("CRC"), "analytic_retention_from_merge"].median())
    t2d_retention = float(per_test.loc[per_test["outcome"].eq("T2D"), "analytic_retention_from_merge"].median())
    rows = [
        {
            "domain": "case_event_density",
            "question": "Why did the CRC application fail to cross the discovery threshold?",
            "evidence": f"CRC assay-specific outcome QC: {int(crc['outcome_case_n_total']):,} cases; T2D: {int(t2d['outcome_case_n_total']):,} cases; median analytic cases {int(crc['analytic_case_n_median']):,} vs {int(t2d['analytic_case_n_median']):,}.",
            "classification": "statistical_readiness_limitation",
            "observed_finding": "CRC event density is much lower despite the same frozen 29-test panel.",
            "improves_precision": "High",
            "improves_temporality": "No",
            "improves_phenotype_resolution": "No",
            "recommended_data_addition": "Increase CRC event count, ideally with incident outcome ascertainment.",
            "boundary": "Does not imply that any future association will be positive.",
        },
        {
            "domain": "analytic_retention",
            "question": "Why did the CRC application lose effective information?",
            "evidence": f"Median analytic retention from merged exposure frame: CRC {100*crc_retention:.1f}% vs T2D {100*t2d_retention:.1f}%.",
            "classification": "statistical_readiness_limitation",
            "observed_finding": "Assay-specific exposure/outcome/covariate overlap and complete-case loss are materially worse for CRC.",
            "improves_precision": "Medium–high",
            "improves_temporality": "No",
            "improves_phenotype_resolution": "No",
            "recommended_data_addition": "Improve assay continuity, covariate completeness, and exposure–outcome overlap.",
            "boundary": "Retention improvement is modeled under non-differential missingness; it does not repair timing.",
        },
        {
            "domain": "observed_signal_architecture",
            "question": "Are all CRC nulls the same kind of null?",
            "evidence": f"Step 9B: {int(attr['signal_limited_observed_near_null'])}/29 near-null; {int(attr['power_limited_nonnull_observed_signal'])}/29 non-near-null but precision-limited under the OR=1.20 reference.",
            "classification": "mixed_statistical_architecture",
            "observed_finding": "The negative screen combines observed near-null signals with directional but imprecise signals.",
            "improves_precision": "Test-dependent",
            "improves_temporality": "No",
            "improves_phenotype_resolution": "No",
            "recommended_data_addition": "Use per-test event/SE diagnostics before deciding whether expansion is informative.",
            "boundary": "Power-limited is a detectability label, not evidence of a true CRC association.",
        },
        {
            "domain": "event_expansion_counterfactual",
            "question": "What statistical improvement matters most?",
            "evidence": f"Approximate MDE<=OR 1.20: current {int(cf['current_tests_meeting_mde_1_20'])}/29; cases×2 {int(cf['event_expansion_tests_meeting_mde_1_20_at_2x'])}/29; ×4 {int(cf['event_expansion_tests_meeting_mde_1_20_at_4x'])}/29; ×8 {int(cf['event_expansion_tests_meeting_mde_1_20_at_8x'])}/29.",
            "classification": "counterfactual_statistical_guidance",
            "observed_finding": "Event expansion improves detectability more strongly than retention improvement in the frozen approximation.",
            "improves_precision": "High",
            "improves_temporality": "No",
            "improves_phenotype_resolution": "No",
            "recommended_data_addition": "Prioritize additional CRC events while preserving assay coverage.",
            "boundary": "SE scaling is an approximation; no new P values or FDR values were generated.",
        },
        {
            "domain": "retention_counterfactual",
            "question": "Would better complete-case retention alone solve the problem?",
            "evidence": f"Approximate MDE<=OR 1.20: retention 70% {int(cf['retention_70pct_tests_meeting_mde_1_20'])}/29; retention 85% {int(cf['retention_85pct_tests_meeting_mde_1_20'])}/29.",
            "classification": "counterfactual_statistical_guidance",
            "observed_finding": "Retention improvement helps but is less influential than large event expansion under the frozen assumptions.",
            "improves_precision": "Medium–high",
            "improves_temporality": "No",
            "improves_phenotype_resolution": "No",
            "recommended_data_addition": "Reduce assay-specific missingness as a complementary, not standalone, intervention.",
            "boundary": "Expected cases and SEs assume non-differential missingness and stable case fraction.",
        },
        {
            "domain": "temporality",
            "question": "Which limitations cannot be fixed by adding N alone?",
            "evidence": "CRC is prevalent and cross-sectional; biomarker collection is contemporaneous with outcome ascertainment and not prediagnostic.",
            "classification": "structural_design_limitation",
            "observed_finding": "The current design cannot establish exposure-before-CRC ordering.",
            "improves_precision": "Low/indirect",
            "improves_temporality": "Very high",
            "improves_phenotype_resolution": "Low",
            "recommended_data_addition": "Prediagnostic biospecimens plus reliable diagnosis date.",
            "boundary": "Temporal ordering reduces a major interpretive limitation but does not guarantee association.",
        },
        {
            "domain": "phenotype_granularity",
            "question": "Which outcome details are missing?",
            "evidence": "No CRC stage, detailed anatomical site beyond colon/rectum, treatment, recurrence/progression, or prospective follow-up fields.",
            "classification": "structural_design_limitation",
            "observed_finding": "The screen cannot resolve clinically meaningful CRC heterogeneity or longitudinal outcomes.",
            "improves_precision": "None",
            "improves_temporality": "Low–high depending on field",
            "improves_phenotype_resolution": "High",
            "recommended_data_addition": "Stage/site plus treatment and follow-up/recurrence ascertainment.",
            "boundary": "Improves answerability and interpretation, not necessarily the pooled effect size.",
        },
        {
            "domain": "exposure_history",
            "question": "How can exposure misclassification and transient measurement be addressed?",
            "evidence": "The current screen has a single survey biomarker measurement per participant.",
            "classification": "structural_exposure_design_limitation",
            "observed_finding": "A single biomarker measurement may not represent persistent exposure before CRC.",
            "improves_precision": "Low–medium",
            "improves_temporality": "High",
            "improves_phenotype_resolution": "Medium",
            "recommended_data_addition": "Repeated urine/serum measurements before outcome ascertainment.",
            "boundary": "Repeated measures improve exposure characterization but do not remove all confounding.",
        },
    ]
    # Keep the supplied priority matrix and structural table as provenance
    # inputs even though this synthesis uses a compact narrative matrix.
    if len(priority) == 0 or len(structural) == 0:
        raise ValueError("Step 9C priority/structural inputs must not be empty")
    return pd.DataFrame(rows)


def build_symmetry(overview: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([
        {
            "branch": "positive",
            "outcome": "T2D",
            "role": "Outcome-firewalled environmental panel followed by disease-specific biological prioritization",
            "terminal_output": "FDR-supported biomarkers/clusters with T2D pathway, network, and transcriptomic convergence",
            "interpretation_boundary": "Biological prioritization is not causal proof.",
        },
        {
            "branch": "negative",
            "outcome": "CRC",
            "role": "Same outcome-firewalled environmental panel followed by question-specific failure diagnosis",
            "terminal_output": "Readiness profile → failure attribution → actionable data-improvement guidance",
            "interpretation_boundary": "No FDR-supported CRC test does not prove absence of environmental relevance.",
        },
    ])


def report_text(synthesis: pd.DataFrame, overview: pd.DataFrame) -> str:
    crc = overview.loc[overview.outcome.eq("CRC")].iloc[0]
    t2d = overview.loc[overview.outcome.eq("T2D")].iloc[0]
    return f"""# Step 9D — Question-Specific Data Readiness Diagnostic

## Purpose

Step 9D formally closes the CRC negative branch by integrating Steps 9A–9C.
It does not rerun models, revise the 29-test family, change BH-FDR, or add
outcome-informed candidate selection. The purpose is to distinguish why the CRC
application did not cross the discovery threshold from what the current design
cannot interpret even if precision improves.

## Final diagnostic statement

The CRC negative screen is best described as a **mixed failure architecture**:

1. **Statistical/readiness limitations:** CRC has {int(crc['outcome_case_n_total']):,} pooled assay-specific QC cases versus {int(t2d['outcome_case_n_total']):,} for T2D, and a median of {int(crc['analytic_case_n_median']):,} analytic cases per test versus {int(t2d['analytic_case_n_median']):,}. Assay-specific complete-case retention is also lower for CRC.
2. **Observed signal heterogeneity:** some tests are near-null, whereas others have directional estimates whose current uncertainty is too large for the prespecified OR=1.20 detectability reference.
3. **Structural design limitations:** prevalent cross-sectional CRC ascertainment, no prediagnostic biospecimen, and no diagnosis date/stage/site/treatment/recurrence/follow-up prevent a strong exposure-before-disease interpretation.

This means **CRC failure does not equal biological absence**, but neither does it
mean that additional cases would guarantee a positive result. The framework
diagnoses answerability rather than manufacturing an explanation for every null.

## Statistical versus structural limitations

| Layer | What the current data show | What it affects | What it does not establish |
|---|---|---|---|
| Statistical/readiness | Low CRC event density and lower assay-specific complete-case retention | Precision, event support, detectability | Does not prove a latent association |
| Observed signal | Near-null and directional/imprecise tests coexist | Which tests merit any future expansion | Does not convert power-limited labels into true effects |
| Structural design | Prevalent cross-sectional outcome and limited phenotype/timing fields | Temporality, reverse-causation protection, phenotype resolution | Cannot be repaired by N alone |

## Counterfactual guidance

Under the frozen SE-scaling approximation, current MDE<=OR 1.20 was met by only
1/29 tests, versus 9/29 at 2× cases, 20/29 at 4×, and 24/29 at 8×. Retention
improvement to 70% and 85% reached 2/29 and 7/29, respectively. These are
detectability scenarios, not new inferential results: no P values, FDR values,
or associations were simulated.

The practical conclusion is that additional CRC events are the strongest
precision intervention in this approximation, while better retention is a
useful complementary intervention. Prediagnostic samples, diagnosis date,
stage/site, longitudinal follow-up, and repeated exposure measurements address
different structural limitations and should not be treated as substitutes for
event support.

## Symmetric framework output

The positive T2D branch ends in disease-specific biological prioritization after
the exposure panel is frozen. The CRC branch ends in a diagnosed failure
architecture and actionable data-improvement guidance after the same firewall.
This symmetry is the methodological result; it is not a claim that positive and
negative disease applications have identical data conditions.

## Scope and provenance

- Source inputs are the locked Step 9A readiness, Step 9B attribution, and Step 9C counterfactual outputs.
- No new model fit, P value, FDR, candidate ranking, or GeneCards/CTD search was performed.
- The full synthesis matrix is in `step9d_negative_branch_synthesis.csv`.
- The branch symmetry table is in `step9d_positive_negative_branch_symmetry.csv`.
"""


def main() -> int:
    overview = read_csv(READINESS_OVERVIEW)
    per_test = read_csv(READINESS_PER_TEST)
    attribution = read_csv(ATTRIBUTION_SUMMARY)
    counterfactual = read_csv(COUNTERFACTUAL_SUMMARY)
    priority = read_csv(PRIORITY_MATRIX)
    structural = read_csv(STRUCTURAL)
    if set(overview["outcome"]) != {"T2D", "CRC"} or len(overview) != 2:
        raise ValueError("Step 9D requires a two-row T2D/CRC readiness overview")
    if len(per_test[per_test["outcome"].eq("CRC")]) != 29:
        raise ValueError("Step 9D requires all 29 CRC per-test readiness rows")

    synthesis = build_synthesis(overview, per_test, attribution, counterfactual, priority, structural)
    symmetry = build_symmetry(overview)
    synthesis.to_csv(SYNTHESIS_OUT, index=False)
    symmetry.to_csv(SYMMETRY_OUT, index=False)
    REPORT_OUT.write_text(report_text(synthesis, overview), encoding="utf-8")

    inputs = [READINESS_OVERVIEW, READINESS_PER_TEST, ATTRIBUTION_SUMMARY, COUNTERFACTUAL_SUMMARY, PRIORITY_MATRIX, STRUCTURAL]
    manifest = {
        "step": "9D",
        "title": "Question-Specific Data Readiness Diagnostic",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in inputs},
        "outputs": [
            str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
            str(SYNTHESIS_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(SYMMETRY_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(REPORT_OUT.relative_to(ROOT)).replace("\\", "/"),
        ],
        "scope_boundaries": {
            "refit_models": False,
            "new_p_values": False,
            "new_fdr_values": False,
            "change_frozen_29_test_family": False,
            "new_gene_cards_or_ctd_search": False,
            "new_candidate_selection": False,
            "purpose": "integrate Steps 9A-9C into a formal negative-branch diagnostic",
        },
        "qc": {
            "overview_rows": int(len(overview)),
            "per_test_crc_rows": int(len(per_test[per_test["outcome"].eq("CRC")])),
            "synthesis_rows": int(len(synthesis)),
            "symmetry_rows": int(len(symmetry)),
            "priority_rows_consumed": int(len(priority)),
            "structural_rows_consumed": int(len(structural)),
        },
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "step": "9D", "synthesis_rows": len(synthesis), "symmetry_rows": len(symmetry)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
