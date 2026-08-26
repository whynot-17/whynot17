#!/usr/bin/env python3
"""Step 9C: counterfactual data-improvement sensitivity analysis for CRC.

This script is intentionally not a power-analysis replacement.  It rescales
the observed Step 9B standard errors under explicit inverse-square-root
assumptions and models improved complete-case retention as an expected-data
scenario.  It never creates new P values, FDR values, or causal claims.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
INPUT = OUT / "step9b_crc_failure_attribution_29_tests.csv"

EVENT_OUT = OUT / "step9c_event_expansion_simulation.csv"
RETENTION_OUT = OUT / "step9c_retention_improvement_simulation.csv"
STRUCTURAL_OUT = OUT / "step9c_structural_outcome_improvements.csv"
PRIORITY_OUT = OUT / "step9c_data_improvement_priority_matrix.csv"
SUMMARY_OUT = OUT / "step9c_counterfactual_summary.csv"
MANIFEST_OUT = OUT / "STEP9C_MANIFEST.json"
REPORT_OUT = OUT / "STEP9C_COUNTERFACTUAL_DATA_IMPROVEMENT_REPORT.md"

Z_ALPHA = 1.959963984540054
Z_POWER = 0.8416212335729143
Z_SUM = Z_ALPHA + Z_POWER
MEANINGFUL_OR = 1.20
TARGET_RETENTIONS = [0.70, 0.85]
EVENT_MULTIPLIERS = [1, 2, 4, 8]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_input() -> pd.DataFrame:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)
    df = pd.read_csv(INPUT)
    if len(df) != 29 or df["test_id"].nunique() != 29:
        raise ValueError("Step 9C requires exactly 29 Step 9B CRC attribution rows")
    numeric_cols = [
        "analytic_n", "analytic_crc_cases", "analytic_controls", "merge_n",
        "analytic_case_fraction", "analytic_retention_from_merge", "beta_logOR",
        "SE_logOR", "OR_per_doubling", "CI_low", "CI_high",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def approximation_fields(beta: float, se: float) -> dict[str, float | bool]:
    if pd.isna(se) or se <= 0:
        return {
            "ci_low_or": np.nan,
            "ci_high_or": np.nan,
            "ci_width_or": np.nan,
            "ci_width_log_or": np.nan,
            "mde_abs_logOR": np.nan,
            "mde_or_upper": np.nan,
            "mde_or_lower": np.nan,
            "detects_reference_or_1_20": False,
        }
    ci_low = math.exp(beta - Z_ALPHA * se)
    ci_high = math.exp(beta + Z_ALPHA * se)
    mde = Z_SUM * se
    mde_upper = math.exp(mde)
    mde_lower = math.exp(-mde)
    return {
        "ci_low_or": ci_low,
        "ci_high_or": ci_high,
        "ci_width_or": ci_high - ci_low,
        "ci_width_log_or": 2 * Z_ALPHA * se,
        "mde_abs_logOR": mde,
        "mde_or_upper": mde_upper,
        "mde_or_lower": mde_lower,
        "detects_reference_or_1_20": bool(mde_upper <= MEANINGFUL_OR),
    }


def event_expansion(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        for multiplier in EVENT_MULTIPLIERS:
            current_se = row["SE_logOR"]
            new_se = current_se / math.sqrt(multiplier) if not pd.isna(current_se) else np.nan
            approx = approximation_fields(row["beta_logOR"], new_se)
            current_cases = row["analytic_crc_cases"]
            rows.append({
                "scenario": "event_expansion",
                "test_id": row["test_id"],
                "biomarker": row["biomarker"],
                "current_analytic_crc_cases": current_cases,
                "case_multiplier": multiplier,
                "expected_crc_cases": current_cases * multiplier if not pd.isna(current_cases) else np.nan,
                "observed_OR_held_constant": row["OR_per_doubling"],
                "current_SE_logOR": current_se,
                "counterfactual_SE_logOR": new_se,
                **approx,
                "reference_effect_OR": MEANINGFUL_OR,
                "assumption": "Effect estimate and nuisance structure held constant; SE scales as 1/sqrt(case multiplier).",
                "new_P_or_FDR_computed": False,
            })
    return pd.DataFrame(rows)


def retention_scenarios(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        current_retention = row["analytic_retention_from_merge"]
        current_n = row["analytic_n"]
        current_cases = row["analytic_crc_cases"]
        current_case_fraction = row["analytic_case_fraction"]
        current_se = row["SE_logOR"]
        for target in [current_retention] + TARGET_RETENTIONS:
            scenario = "current_observed_retention" if target == current_retention else f"retention_{int(target * 100)}pct"
            effective_target = max(current_retention, target) if not pd.isna(current_retention) else target
            expected_n = row["merge_n"] * effective_target if not pd.isna(row["merge_n"]) else np.nan
            if not pd.isna(row["merge_n"]):
                expected_n = min(expected_n, row["merge_n"])
            expected_cases = expected_n * current_case_fraction if not pd.isna(expected_n) and not pd.isna(current_case_fraction) else np.nan
            if not pd.isna(current_se) and not pd.isna(expected_cases) and expected_cases > 0 and current_cases > 0:
                new_se = current_se * math.sqrt(current_cases / expected_cases)
            else:
                new_se = np.nan
            approx = approximation_fields(row["beta_logOR"], new_se)
            rows.append({
                "scenario": scenario,
                "test_id": row["test_id"],
                "biomarker": row["biomarker"],
                "current_analytic_n": current_n,
                "current_analytic_crc_cases": current_cases,
                "merge_n_frame": row["merge_n"],
                "current_retention": current_retention,
                "target_retention": target,
                "effective_retention_used": effective_target,
                "expected_analytic_n": expected_n,
                "expected_crc_cases": expected_cases,
                "observed_OR_held_constant": row["OR_per_doubling"],
                "current_SE_logOR": current_se,
                "counterfactual_SE_logOR": new_se,
                **approx,
                "reference_effect_OR": MEANINGFUL_OR,
                "assumption": "Improved complete-case retention is non-differential with respect to outcome; observed case fraction and effect are held constant.",
                "new_P_or_FDR_computed": False,
            })
    return pd.DataFrame(rows)


def structural_improvements() -> pd.DataFrame:
    rows = [
        {
            "improvement": "Increase CRC event count",
            "current_state": "Low effective CRC event count per biomarker model",
            "target_state": "More incident/prevalent CRC events with the same frozen exposure panel",
            "precision_gain": "High",
            "temporality_gain": "None by itself",
            "phenotype_resolution_gain": "None by itself",
            "priority": "High",
            "interpretive_role": "Improves detectability; does not establish temporality or causality.",
        },
        {
            "improvement": "Improve assay-specific analytic retention",
            "current_state": "CRC median complete-case retention is lower than T2D",
            "target_state": "More exposure/outcome/covariate overlap and fewer complete-case losses",
            "precision_gain": "Medium–high",
            "temporality_gain": "None",
            "phenotype_resolution_gain": "None",
            "priority": "High",
            "interpretive_role": "Improves effective sample size under a non-differential missingness assumption.",
        },
        {
            "improvement": "Prediagnostic biospecimen",
            "current_state": "Current survey biomarker measured with prevalent CRC status",
            "target_state": "Biomarker collected before CRC diagnosis",
            "precision_gain": "Low/indirect",
            "temporality_gain": "Very high",
            "phenotype_resolution_gain": "Low",
            "priority": "High",
            "interpretive_role": "Directly reduces the simplest reverse-causation concern; does not guarantee association.",
        },
        {
            "improvement": "Diagnosis date",
            "current_state": "No diagnosis date in the primary CRC screen",
            "target_state": "Reliable diagnosis date and exposure-to-diagnosis interval",
            "precision_gain": "None",
            "temporality_gain": "High",
            "phenotype_resolution_gain": "Medium",
            "priority": "High",
            "interpretive_role": "Enables lag and timing analyses without changing the exposure panel.",
        },
        {
            "improvement": "Stage and anatomical site",
            "current_state": "Only colon/rectum classification; no stage",
            "target_state": "Stage, subsite, and clinically meaningful tumor phenotype",
            "precision_gain": "None",
            "temporality_gain": "Low",
            "phenotype_resolution_gain": "High",
            "priority": "Medium–high",
            "interpretive_role": "Improves phenotype resolution and etiologic heterogeneity assessment.",
        },
        {
            "improvement": "Treatment, follow-up, and recurrence/progression",
            "current_state": "No prospective follow-up or recurrence timeline",
            "target_state": "Longitudinal outcomes after biomarker collection",
            "precision_gain": "None/indirect",
            "temporality_gain": "High",
            "phenotype_resolution_gain": "High",
            "priority": "High",
            "interpretive_role": "Supports incidence, prognosis, and recurrence analyses rather than a prevalent screen.",
        },
        {
            "improvement": "Repeated exposure measurements",
            "current_state": "Single survey biomarker measurement",
            "target_state": "Repeated urine/serum measurements before outcome ascertainment",
            "precision_gain": "Low–medium",
            "temporality_gain": "High",
            "phenotype_resolution_gain": "Medium",
            "priority": "High",
            "interpretive_role": "Reduces exposure misclassification and helps distinguish persistent from transient exposure.",
        },
    ]
    return pd.DataFrame(rows)


def priority_matrix() -> pd.DataFrame:
    rows = [
        ("Increase CRC events", "High", "None", "None", "High"),
        ("Improve assay retention", "Medium–high", "None", "None", "High"),
        ("Prediagnostic biospecimen", "Low/indirect", "Very high", "Low", "High"),
        ("Diagnosis date", "None", "High", "Medium", "High"),
        ("Stage/site", "None", "Low", "High", "Medium–high"),
        ("Follow-up/recurrence", "None", "High", "High", "High"),
        ("Repeated exposure measurements", "Low–medium", "High", "Medium", "High"),
    ]
    return pd.DataFrame(rows, columns=["improvement", "precision", "temporality", "phenotype_resolution", "priority"])


def summary_table(event: pd.DataFrame, retention: pd.DataFrame, df: pd.DataFrame) -> pd.DataFrame:
    current_event = event[event["case_multiplier"] == 1]
    rows = [
        {"metric": "frozen_crc_tests", "value": len(df), "unit": "tests"},
        {"metric": "event_expansion_tests_meeting_mde_1_20_at_2x", "value": int(event.loc[event.case_multiplier.eq(2), "detects_reference_or_1_20"].sum()), "unit": "tests"},
        {"metric": "event_expansion_tests_meeting_mde_1_20_at_4x", "value": int(event.loc[event.case_multiplier.eq(4), "detects_reference_or_1_20"].sum()), "unit": "tests"},
        {"metric": "event_expansion_tests_meeting_mde_1_20_at_8x", "value": int(event.loc[event.case_multiplier.eq(8), "detects_reference_or_1_20"].sum()), "unit": "tests"},
        {"metric": "event_expansion_tests_not_meeting_mde_1_20_at_8x", "value": int((~event.loc[event.case_multiplier.eq(8), "detects_reference_or_1_20"]).sum()), "unit": "tests"},
        {"metric": "retention_70pct_tests_meeting_mde_1_20", "value": int(retention.loc[retention.scenario.eq("retention_70pct"), "detects_reference_or_1_20"].sum()), "unit": "tests"},
        {"metric": "retention_85pct_tests_meeting_mde_1_20", "value": int(retention.loc[retention.scenario.eq("retention_85pct"), "detects_reference_or_1_20"].sum()), "unit": "tests"},
        {"metric": "current_tests_meeting_mde_1_20", "value": int(current_event["detects_reference_or_1_20"].sum()), "unit": "tests"},
    ]
    return pd.DataFrame(rows)


def report_text(event: pd.DataFrame, retention: pd.DataFrame, structural: pd.DataFrame, summary: pd.DataFrame) -> str:
    def value(metric):
        return int(summary.loc[summary.metric.eq(metric), "value"].iloc[0])
    current = event[event.case_multiplier.eq(1)]
    e8 = event[event.case_multiplier.eq(8)]
    return f"""# Step 9C — Counterfactual data-improvement sensitivity analysis

## Scope and firewall

This is a planning sensitivity analysis for the CRC negative branch. It does
not re-fit the CRC models, create new P values or BH-FDR values, change the
29-test family, or claim that an enlarged study would produce an association.
The event and retention sections only rescale the recorded Step 9B estimates
under explicit assumptions.

## 1. Event expansion

For each test, cases were scaled by ×1, ×2, ×4, and ×8. The observed log(OR)
was held constant and `SE_new = SE_current / sqrt(case multiplier)`. The output
reports expected case count, CI width, and the approximate 80%-power MDE. The
reference effect is OR=1.20 per exposure doubling; it is a planning reference,
not a clinical threshold.

| Scenario | Tests with approximate MDE <= OR 1.20 |
|---|---:|
| Current | {value('current_tests_meeting_mde_1_20')} / 29 |
| Cases ×2 | {value('event_expansion_tests_meeting_mde_1_20_at_2x')} / 29 |
| Cases ×4 | {value('event_expansion_tests_meeting_mde_1_20_at_4x')} / 29 |
| Cases ×8 | {value('event_expansion_tests_meeting_mde_1_20_at_8x')} / 29 |
| Not meeting the MDE reference even at ×8 | {value('event_expansion_tests_not_meeting_mde_1_20_at_8x')} / 29 |

This simulation answers detectability under the stated scaling assumption. It
does not identify which tests are biologically real, and it does not simulate
selection, confounding, exposure correlation, survey-design changes, or FDR.
Tests whose observed OR is near the null remain near-null in the scenario even
if their precision improves.

## 2. Analytic-retention improvement

For targets of 70% and 85% retention from the exposure-merged frame, expected
analytic N was set to `merge N × target retention` (capped at merge N), expected
case fraction was held at the observed analytic case fraction, and SE was
scaled by the square root of the current/expected case count. This is a
non-differential missingness approximation. It is not a reanalysis of any
participant-level data.

| Scenario | Tests with approximate MDE <= OR 1.20 |
|---|---:|
| Retention 70% | {value('retention_70pct_tests_meeting_mde_1_20')} / 29 |
| Retention 85% | {value('retention_85pct_tests_meeting_mde_1_20')} / 29 |

The retention table preserves each biomarker's current merge N, current case
fraction, expected analytic N, expected cases, and counterfactual MDE. It does
not assume that better retention repairs outcome timing or phenotype definition.

## 3. Structural outcome improvements

Structural changes cannot be assigned a new P value from the present NHANES
data. Their value is inferential rather than purely statistical:

- prediagnostic biospecimens primarily improve temporal ordering and protection against reverse causation;
- diagnosis date enables lag/timing analyses;
- stage and site improve phenotype resolution;
- treatment/follow-up/recurrence support incidence, prognosis, and longitudinal outcome analyses;
- repeated exposure measurements reduce exposure misclassification and identify persistent exposure.

The full structural matrix and priority matrix are machine-readable outputs.

## Overall interpretation

The CRC negative branch is therefore actionable but not self-exonerating. Event
expansion and better retention can improve detectability under transparent
assumptions, while prediagnostic and longitudinal outcome data address a
different limitation—temporality and interpretation. No row in this report
means that a future study is expected to be positive; it means that a specified
data addition would improve the ability to answer the CRC question.

## Outputs

- `step9c_event_expansion_simulation.csv`
- `step9c_retention_improvement_simulation.csv`
- `step9c_structural_outcome_improvements.csv`
- `step9c_data_improvement_priority_matrix.csv`
- `step9c_counterfactual_summary.csv`
- `STEP9C_MANIFEST.json`
"""


def main() -> int:
    df = read_input()
    event = event_expansion(df)
    retention = retention_scenarios(df)
    structural = structural_improvements()
    priority = priority_matrix()
    summary = summary_table(event, retention, df)

    event.to_csv(EVENT_OUT, index=False)
    retention.to_csv(RETENTION_OUT, index=False)
    structural.to_csv(STRUCTURAL_OUT, index=False)
    priority.to_csv(PRIORITY_OUT, index=False)
    summary.to_csv(SUMMARY_OUT, index=False)
    REPORT_OUT.write_text(report_text(event, retention, structural, summary), encoding="utf-8")

    inputs = [INPUT]
    manifest = {
        "step": "9C",
        "title": "Counterfactual data-improvement sensitivity analysis",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in inputs},
        "outputs": [
            str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
            str(EVENT_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(RETENTION_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(STRUCTURAL_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(PRIORITY_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(SUMMARY_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(REPORT_OUT.relative_to(ROOT)).replace("\\", "/"),
        ],
        "frozen_assumptions": {
            "event_multipliers": EVENT_MULTIPLIERS,
            "event_se_scaling": "SE_new = SE_current / sqrt(case multiplier)",
            "retention_targets": TARGET_RETENTIONS,
            "retention_se_scaling": "SE_new = SE_current * sqrt(current cases / expected cases)",
            "retention_missingness": "non-differential with respect to outcome",
            "reference_effect_or": MEANINGFUL_OR,
            "mde": "two-sided alpha=0.05, 80% power, normal approximation",
            "refit_models": False,
            "new_p_values": False,
            "new_fdr_values": False,
            "structural_changes_get_no_fabricated_p_values": True,
        },
        "qc": {
            "input_tests": int(len(df)),
            "event_rows": int(len(event)),
            "retention_rows": int(len(retention)),
            "structural_rows": int(len(structural)),
            "priority_rows": int(len(priority)),
            "summary_rows": int(len(summary)),
        },
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "step": "9C", "summary": summary.to_dict(orient="records")}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
