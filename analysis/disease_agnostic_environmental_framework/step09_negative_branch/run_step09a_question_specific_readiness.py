#!/usr/bin/env python3
"""Step 9A: question-specific data-readiness audit for the frozen 29-test panel.

This module compares the assay-specific rebuilt CRC screen with the frozen T2D
screen using the same descriptive metrics.  It is deliberately descriptive:
it does not re-fit models, alter the 29-test family, attribute individual nulls,
or run counterfactual power simulations.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
CRC_DIR = ROOT / "analysis/disease_agnostic_environmental_framework/step05_crc_screen_rebuilt"
T2D_DIR = ROOT / "analysis/disease_agnostic_environmental_framework/step05_t2d_screen"
CRC_LEDGER = ROOT / "analysis/disease_agnostic_environmental_framework/step05_crc_screen/CRC_case_control_ledger.csv"

CRC_RESULTS = CRC_DIR / "full_29_test_crc_screen_rebuilt.csv"
CRC_MERGE = CRC_DIR / "assay_specific_merge_audit.csv"
CRC_OUTCOME_QC = CRC_DIR / "assay_specific_outcome_frame_qc.csv"
T2D_RESULTS = T2D_DIR / "t2d_primary_29_tests.csv"
T2D_MERGE = T2D_DIR / "t2d_merge_audit.csv"
T2D_OUTCOME_QC = T2D_DIR / "t2d_outcome_qc.csv"

PER_TEST_OUT = OUT / "step9a_per_test_readiness.csv"
OVERVIEW_OUT = OUT / "step9a_readiness_overview.csv"
GRANULARITY_OUT = OUT / "step9a_outcome_granularity.csv"
MANIFEST_OUT = OUT / "STEP9A_MANIFEST.json"
REPORT_OUT = OUT / "STEP9A_QUESTION_SPECIFIC_READINESS_REPORT.md"


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def clean_num(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def parse_cycles(value) -> list[str]:
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def fmt_num(value, digits=2):
    if pd.isna(value):
        return "NA"
    return f"{float(value):,.{digits}f}"


def fmt_int(value):
    if pd.isna(value):
        return "NA"
    return f"{int(round(float(value))):,}"


def fmt_pct(value, digits=1):
    if pd.isna(value):
        return "NA"
    return f"{100 * float(value):.{digits}f}%"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")


def outcome_granularity(crc_ledger: pd.DataFrame, crc_qc: pd.DataFrame, t2d_qc: pd.DataFrame) -> pd.DataFrame:
    require_columns(crc_ledger, ["crc_case", "cancer_free", "crc_diagnosis_age"], "CRC ledger")
    # The assay-specific rebuilt screen has its own universal outcome-frame QC.
    # Use that frame for the primary readiness counts.  The older case/control
    # ledger is retained only for diagnosis-age availability and reconciliation;
    # its totals are not silently mixed with the rebuilt screen totals.
    crc_qc_case_n = int(pd.to_numeric(crc_qc["crc_cases"], errors="coerce").sum())
    crc_qc_control_n = int(pd.to_numeric(crc_qc["cancer_free_controls"], errors="coerce").sum())
    crc_qc_frame_n = int(pd.to_numeric(crc_qc["universal_frame_rows"], errors="coerce").sum())
    crc_qc_eligible_n = int(pd.to_numeric(crc_qc["crc_vs_cancer_free_rows"], errors="coerce").sum())
    crc_ledger_case_n = int(pd.to_numeric(crc_ledger["crc_case"], errors="coerce").fillna(0).sum())
    crc_ledger_control_n = int(pd.to_numeric(crc_ledger["cancer_free"], errors="coerce").fillna(0).sum())
    crc_diag_age_n = int(pd.to_numeric(crc_ledger.loc[pd.to_numeric(crc_ledger["crc_case"], errors="coerce") == 1, "crc_diagnosis_age"], errors="coerce").notna().sum())

    t2d_case_n = int(pd.to_numeric(t2d_qc["t2d_case_n"], errors="coerce").sum())
    t2d_control_n = int(pd.to_numeric(t2d_qc["non_diabetic_control_n"], errors="coerce").sum())
    t2d_eligible_n = int(pd.to_numeric(t2d_qc["t2d_eligible_n"], errors="coerce").sum())
    t2d_adult_frame_n = int(pd.to_numeric(t2d_qc["adult_outcome_frame_rows"], errors="coerce").sum())

    rows = [
        {
            "outcome": "CRC",
            "outcome_type": "prevalent CRC case-control screen",
            "primary_case_definition": "Adult (>=20 y): MCQ220=1 and colon/rectum cancer-type code 16 or 31",
            "control_definition": "MCQ220=2; non-CRC cancers excluded",
            "adult_or_outcome_frame_n": crc_qc_frame_n,
            "pooled_eligible_outcome_n": crc_qc_eligible_n,
            "pooled_case_n": crc_qc_case_n,
            "pooled_control_n": crc_qc_control_n,
            "outcome_case_fraction": crc_qc_case_n / (crc_qc_case_n + crc_qc_control_n),
            "diagnosis_age_available": f"{crc_diag_age_n}/{crc_ledger_case_n} CRC cases in the separate outcome ledger",
            "diagnosis_date_available": "No",
            "stage_available": "No",
            "site_available": "No beyond colon/rectum classification",
            "treatment_available": "No",
            "recurrence_progression_available": "No",
            "prediagnostic_biospecimen": "No; current survey biomarker with prevalent outcome",
            "repeated_exposure_measurements": "No participant-level repeated urine/serum measurements in this screen",
            "follow_up": "No incident follow-up",
            "temporality_assessment": "Limited: cross-sectional/prevalent outcome",
            "interpretation_risk": "Reverse causation and survivor bias cannot be excluded; diagnosis age supports timing sensitivity only.",
            "secondary_ledger_rows_n": int(len(crc_ledger)),
            "secondary_ledger_case_n": crc_ledger_case_n,
            "secondary_ledger_control_n": crc_ledger_control_n,
            "secondary_ledger_reconciliation": "Not directly comparable to assay-specific rebuilt outcome QC; retained only as provenance for diagnosis-age availability.",
        },
        {
            "outcome": "T2D",
            "outcome_type": "cross-sectional T2D case-control screen",
            "primary_case_definition": "Adult T2D: diagnosed diabetes except likely early-onset insulin-dependent cases, or undiagnosed HbA1c >=6.5%",
            "control_definition": "DIQ010=No with available HbA1c <6.5%; ambiguous/missing categories indeterminate",
            "adult_or_outcome_frame_n": t2d_adult_frame_n,
            "pooled_eligible_outcome_n": t2d_eligible_n,
            "pooled_case_n": t2d_case_n,
            "pooled_control_n": t2d_control_n,
            "outcome_case_fraction": t2d_case_n / t2d_eligible_n,
            "diagnosis_age_available": "No diagnosis-age variable used for this screen",
            "diagnosis_date_available": "No",
            "stage_available": "Not applicable",
            "site_available": "Not applicable",
            "treatment_available": "No incident-treatment timeline",
            "recurrence_progression_available": "Not applicable",
            "prediagnostic_biospecimen": "No; same-visit biomarker and disease classification",
            "repeated_exposure_measurements": "No participant-level repeated urine/serum measurements in this screen",
            "follow_up": "No incident follow-up",
            "temporality_assessment": "Limited: cross-sectional disease classification, with objective HbA1c component",
            "interpretation_risk": "Temporal ordering and long-term exposure history are not established, but outcome case density is substantially higher than CRC.",
            "secondary_ledger_rows_n": np.nan,
            "secondary_ledger_case_n": np.nan,
            "secondary_ledger_control_n": np.nan,
            "secondary_ledger_reconciliation": "Not applicable",
        },
    ]
    return pd.DataFrame(rows)


def build_per_test(outcome: str, results: pd.DataFrame, merge: pd.DataFrame, outcome_qc: pd.DataFrame) -> pd.DataFrame:
    if outcome == "CRC":
        case_col = "analytic_crc_cases"
        control_col = "analytic_controls"
        p_col = "P"
        q_col = "BH_FDR"
        outcome_label = "CRC"
        merge_outcome_col = "crc_vs_cancer_free_rows"
        eligible_retention_col = "crc_vs_cancer_free_rows"
        fit_status_col = "status"
    else:
        case_col = "analytic_t2d_cases"
        control_col = "analytic_controls"
        p_col = "P"
        q_col = "BH_FDR"
        outcome_label = "T2D"
        merge_outcome_col = "t2d_eligible_rows"
        eligible_retention_col = "t2d_eligible_rows"
        fit_status_col = "status"

    audit = merge.copy()
    numeric_audit_cols = [
        "exposure_rows", "exposure_nonmissing", "universal_outcome_covariate_rows",
        "outcome_frame_rows", "merge_rows", "crc_vs_cancer_free_rows",
        "t2d_eligible_rows", "complete_case_rows", "complete_case_crc_cases",
        "complete_case_t2d_cases",
    ]
    for col in numeric_audit_cols:
        if col in audit:
            audit[col] = pd.to_numeric(audit[col], errors="coerce")

    outcome_frame_col = "universal_outcome_covariate_rows" if outcome == "CRC" else "outcome_frame_rows"
    audit["complete_case_cases"] = audit["complete_case_crc_cases"] if outcome == "CRC" else audit["complete_case_t2d_cases"]
    audit["eligible_rows"] = audit[eligible_retention_col] if eligible_retention_col in audit else np.nan

    records = []
    for _, row in results.iterrows():
        test_id = row["test_id"]
        sub = audit[audit["test_id"] == test_id].copy()
        cycles = parse_cycles(row.get("frozen_cycle_list", ""))
        source_cycle_list = parse_cycles(row.get("source_cycles_read", ""))
        source_rows = pd.to_numeric(sub.get("exposure_rows", pd.Series(dtype=float)), errors="coerce").sum(min_count=1)
        source_nonmissing = pd.to_numeric(sub.get("exposure_nonmissing", pd.Series(dtype=float)), errors="coerce").sum(min_count=1)
        merge_rows = pd.to_numeric(sub.get("merge_rows", pd.Series(dtype=float)), errors="coerce").sum(min_count=1)
        eligible_rows = pd.to_numeric(sub.get("eligible_rows", pd.Series(dtype=float)), errors="coerce").sum(min_count=1)
        complete_rows = pd.to_numeric(sub.get("complete_case_rows", pd.Series(dtype=float)), errors="coerce").sum(min_count=1)
        complete_cases = pd.to_numeric(sub.get("complete_case_cases", pd.Series(dtype=float)), errors="coerce").sum(min_count=1)
        cycle_complete = pd.to_numeric(sub.get("complete_case_rows", pd.Series(dtype=float)), errors="coerce")
        cycle_cases = pd.to_numeric(sub.get("complete_case_cases", pd.Series(dtype=float)), errors="coerce")
        result_n = clean_num(row.get("analytic_n", row.get("N")))
        result_cases = clean_num(row.get(case_col, row.get("CRC_N")))
        result_controls = clean_num(row.get(control_col, row.get("Control_N")))
        result_merge = clean_num(row.get("merge_N"))
        p_value = clean_num(row.get(p_col))
        q_value = clean_num(row.get(q_col))
        ci_low = clean_num(row.get("CI_low"))
        ci_high = clean_num(row.get("CI_high"))

        records.append({
            "outcome": outcome_label,
            "test_id": test_id,
            "biomarker": row.get("biomarker", row.get("variable")),
            "variable": row.get("variable"),
            "exposure_axis": row.get("exposure_axis"),
            "matrix": row.get("matrix"),
            "cycles_n_frozen": len(cycles),
            "cycle_list_frozen": ";".join(cycles),
            "cycles_n_source_read": len(source_cycle_list),
            "cycle_list_source_read": ";".join(source_cycle_list),
            "source_registry_n": clean_num(row.get("source_registry_n")),
            "source_rows_total_from_cycle_audit": source_rows,
            "source_nonmissing_total_from_cycle_audit": source_nonmissing,
            "source_missing_fraction": (1 - source_nonmissing / source_rows) if source_rows and not pd.isna(source_rows) else np.nan,
            "outcome_frame_rows_total_from_cycle_audit": pd.to_numeric(sub.get(outcome_frame_col, pd.Series(dtype=float)), errors="coerce").sum(min_count=1),
            "eligible_rows_total_from_cycle_audit": eligible_rows,
            "merge_n_result": result_merge,
            "merge_rows_total_from_cycle_audit": merge_rows,
            "analytic_n_result": result_n,
            "complete_case_rows_total_from_cycle_audit": complete_rows,
            "analytic_cases_result": result_cases,
            "complete_case_cases_total_from_cycle_audit": complete_cases,
            "analytic_controls_result": result_controls,
            "case_fraction": (result_cases / result_n) if result_n and not pd.isna(result_n) else np.nan,
            "analytic_retention_from_merge": (result_n / result_merge) if result_merge and not pd.isna(result_merge) else np.nan,
            "analytic_retention_from_eligible": (result_n / eligible_rows) if eligible_rows and not pd.isna(eligible_rows) else np.nan,
            "cycles_with_merge_rows": int((pd.to_numeric(sub.get("merge_rows", pd.Series(dtype=float)), errors="coerce") > 0).sum()),
            "cycles_with_complete_cases": int((cycle_complete > 0).sum()),
            "min_cycle_complete_case_n": cycle_complete[cycle_complete > 0].min() if (cycle_complete > 0).any() else np.nan,
            "max_cycle_complete_case_n": cycle_complete.max() if len(cycle_complete) else np.nan,
            "min_cycle_case_n": cycle_cases[cycle_cases > 0].min() if (cycle_cases > 0).any() else np.nan,
            "max_cycle_case_n": cycle_cases.max() if len(cycle_cases) else np.nan,
            "OR_per_doubling": clean_num(row.get("OR_per_doubling", row.get("OR"))),
            "CI_low": ci_low,
            "CI_high": ci_high,
            "CI_width": ci_high - ci_low if not pd.isna(ci_low) and not pd.isna(ci_high) else np.nan,
            "nominal_P": p_value,
            "BH_FDR": q_value,
            "fit_status": row.get(fit_status_col),
            "design_df": clean_num(row.get("design_df")),
            "outcome_definition_class": "Prevalent CRC case-control" if outcome == "CRC" else "Cross-sectional T2D classification",
            "temporality_class": "Current biomarker vs prevalent outcome" if outcome == "CRC" else "Same-visit biomarker and disease classification",
            "diagnosis_date_available": "No",
            "stage_site_treatment_recurrence_available": "No" if outcome == "CRC" else "Not applicable",
            "follow_up_available": "No",
        })

    result = pd.DataFrame(records)
    result["audit_result_n_difference"] = result["analytic_n_result"] - result["complete_case_rows_total_from_cycle_audit"]
    result["audit_case_n_difference"] = result["analytic_cases_result"] - result["complete_case_cases_total_from_cycle_audit"]
    return result


def build_overview(per_test: pd.DataFrame, granularity: pd.DataFrame, outcome_qc_map: dict[str, pd.DataFrame]) -> pd.DataFrame:
    records = []
    for outcome in ["T2D", "CRC"]:
        sub = per_test[per_test["outcome"] == outcome].copy()
        gran = granularity[granularity["outcome"] == outcome].iloc[0]
        qc = outcome_qc_map[outcome]
        if outcome == "CRC":
            outcome_frame_rows = pd.to_numeric(qc["universal_frame_rows"], errors="coerce").sum()
            eligible_rows = pd.to_numeric(qc["crc_vs_cancer_free_rows"], errors="coerce").sum()
            case_n = pd.to_numeric(qc["crc_cases"], errors="coerce").sum()
            control_n = pd.to_numeric(qc["cancer_free_controls"], errors="coerce").sum()
        else:
            outcome_frame_rows = pd.to_numeric(qc["adult_outcome_frame_rows"], errors="coerce").sum()
            eligible_rows = pd.to_numeric(qc["t2d_eligible_n"], errors="coerce").sum()
            case_n = pd.to_numeric(qc["t2d_case_n"], errors="coerce").sum()
            control_n = pd.to_numeric(qc["non_diabetic_control_n"], errors="coerce").sum()

        p = pd.to_numeric(sub["nominal_P"], errors="coerce")
        q = pd.to_numeric(sub["BH_FDR"], errors="coerce")
        analytic_n = pd.to_numeric(sub["analytic_n_result"], errors="coerce")
        analytic_cases = pd.to_numeric(sub["analytic_cases_result"], errors="coerce")
        model_estimable = pd.to_numeric(sub["nominal_P"], errors="coerce").notna() & pd.to_numeric(sub["analytic_n_result"], errors="coerce").gt(0)
        records.append({
            "outcome": outcome,
            "outcome_type": gran["outcome_type"],
            "nhanes_cycle_n_in_outcome_qc": int(len(qc)),
            "outcome_frame_rows_total": int(outcome_frame_rows),
            "eligible_outcome_rows_total": int(eligible_rows),
            "outcome_case_n_total": int(case_n),
            "outcome_control_n_total": int(control_n),
            "outcome_case_fraction": float(case_n / eligible_rows),
            "frozen_tests": int(len(sub)),
            "estimable_tests": int(model_estimable.sum()),
            "fit_warning_n": int((sub["fit_status"].astype(str).str.lower() != "ok").sum()),
            "analytic_n_min": float(analytic_n.min()),
            "analytic_n_median": float(analytic_n.median()),
            "analytic_n_max": float(analytic_n.max()),
            "analytic_case_n_min": float(analytic_cases.min()),
            "analytic_case_n_median": float(analytic_cases.median()),
            "analytic_case_n_max": float(analytic_cases.max()),
            "analytic_case_fraction_median": float(pd.to_numeric(sub["case_fraction"], errors="coerce").median()),
            "cycle_coverage_min": int(sub["cycles_n_frozen"].min()),
            "cycle_coverage_median": float(sub["cycles_n_frozen"].median()),
            "cycle_coverage_max": int(sub["cycles_n_frozen"].max()),
            "cycles_with_complete_cases_min": int(sub["cycles_with_complete_cases"].min()),
            "cycles_with_complete_cases_median": float(sub["cycles_with_complete_cases"].median()),
            "cycles_with_complete_cases_max": int(sub["cycles_with_complete_cases"].max()),
            "nominal_p_lt_0_05_n": int((p < 0.05).sum()),
            "bh_fdr_lt_0_05_n": int((q < 0.05).sum()),
            "bh_fdr_denominator": 29,
            "temporality_assessment": gran["temporality_assessment"],
            "diagnosis_date_available": gran["diagnosis_date_available"],
            "stage_site_treatment_recurrence_available": "; ".join([gran["stage_available"], gran["site_available"], gran["treatment_available"], gran["recurrence_progression_available"]]),
            "prediagnostic_biospecimen": gran["prediagnostic_biospecimen"],
            "follow_up": gran["follow_up"],
        })
    return pd.DataFrame(records)


def report_text(overview: pd.DataFrame, per_test: pd.DataFrame, granularity: pd.DataFrame) -> str:
    t2d = overview[overview["outcome"] == "T2D"].iloc[0]
    crc = overview[overview["outcome"] == "CRC"].iloc[0]
    t2d_tests = per_test[per_test["outcome"] == "T2D"]
    crc_tests = per_test[per_test["outcome"] == "CRC"]
    crc_n_loss = 1 - crc_tests["analytic_n_result"].median() / crc_tests["merge_n_result"].median()
    t2d_n_loss = 1 - t2d_tests["analytic_n_result"].median() / t2d_tests["merge_n_result"].median()
    return f"""# Step 9A — Question-specific data-readiness audit

## Scope and status

This is a descriptive audit of the same frozen 29-test environmental biomarker
panel in the T2D and assay-specific rebuilt CRC outcome screens. It is the first
layer of the CRC negative branch. It does **not** re-fit models, change the
29-test multiplicity family, perform failure attribution, or simulate larger
studies. The purpose is to document whether the two disease questions are being
asked under comparable data conditions.

The analysis deliberately preserves a multidimensional readiness profile rather
than collapsing it into a single score.

## Executive comparison

| Metric | T2D | CRC |
|---|---:|---:|
| NHANES cycles in outcome QC | {int(t2d['nhanes_cycle_n_in_outcome_qc'])} | {int(crc['nhanes_cycle_n_in_outcome_qc'])} |
| Outcome frame rows | {fmt_int(t2d['outcome_frame_rows_total'])} | {fmt_int(crc['outcome_frame_rows_total'])} |
| Eligible case/control rows | {fmt_int(t2d['eligible_outcome_rows_total'])} | {fmt_int(crc['eligible_outcome_rows_total'])} |
| Outcome cases in pooled outcome QC | {fmt_int(t2d['outcome_case_n_total'])} | {fmt_int(crc['outcome_case_n_total'])} |
| Outcome case fraction among eligible rows | {fmt_pct(t2d['outcome_case_fraction'])} | {fmt_pct(crc['outcome_case_fraction'])} |
| Frozen tests | {int(t2d['frozen_tests'])} | {int(crc['frozen_tests'])} |
| Estimable tests | {int(t2d['estimable_tests'])} | {int(crc['estimable_tests'])} |
| Fit warnings retained in status | {int(t2d['fit_warning_n'])} | {int(crc['fit_warning_n'])} |
| Analytic N across tests (min / median / max) | {fmt_int(t2d['analytic_n_min'])} / {fmt_int(t2d['analytic_n_median'])} / {fmt_int(t2d['analytic_n_max'])} | {fmt_int(crc['analytic_n_min'])} / {fmt_int(crc['analytic_n_median'])} / {fmt_int(crc['analytic_n_max'])} |
| Analytic cases across tests (min / median / max) | {fmt_int(t2d['analytic_case_n_min'])} / {fmt_int(t2d['analytic_case_n_median'])} / {fmt_int(t2d['analytic_case_n_max'])} | {fmt_int(crc['analytic_case_n_min'])} / {fmt_int(crc['analytic_case_n_median'])} / {fmt_int(crc['analytic_case_n_max'])} |
| Median analytic case fraction | {fmt_pct(t2d['analytic_case_fraction_median'])} | {fmt_pct(crc['analytic_case_fraction_median'])} |
| Frozen cycle coverage (min / median / max) | {int(t2d['cycle_coverage_min'])} / {fmt_num(t2d['cycle_coverage_median'], 1)} / {int(t2d['cycle_coverage_max'])} | {int(crc['cycle_coverage_min'])} / {fmt_num(crc['cycle_coverage_median'], 1)} / {int(crc['cycle_coverage_max'])} |
| Cycles with complete cases (min / median / max) | {int(t2d['cycles_with_complete_cases_min'])} / {fmt_num(t2d['cycles_with_complete_cases_median'], 1)} / {int(t2d['cycles_with_complete_cases_max'])} | {int(crc['cycles_with_complete_cases_min'])} / {fmt_num(crc['cycles_with_complete_cases_median'], 1)} / {int(crc['cycles_with_complete_cases_max'])} |
| Median analytic retention from merged exposure frame | {fmt_pct(t2d_tests['analytic_retention_from_merge'].median())} | {fmt_pct(crc_tests['analytic_retention_from_merge'].median())} |
| Approx. median loss after merge to analytic complete cases | {fmt_pct(t2d_n_loss)} | {fmt_pct(crc_n_loss)} |
| Nominal P<0.05 (context only) | {int(t2d['nominal_p_lt_0_05_n'])}/29 | {int(crc['nominal_p_lt_0_05_n'])}/29 |
| BH-FDR<0.05 (context only) | {int(t2d['bh_fdr_lt_0_05_n'])}/29 | {int(crc['bh_fdr_lt_0_05_n'])}/29 |

The main readiness contrast is not simply “significant versus null.” CRC has a
much lower case density and fewer outcome events per biomarker-specific analytic
sample, while its primary outcome is prevalent and cross-sectional. T2D has a
larger outcome case pool and all 29 tests are estimable, but it is also a
same-visit cross-sectional classification rather than prospective incidence.

The CRC rebuilt output contains one `converged_with_warning` fit. It remains
estimable because it has a finite model P value and a non-zero analytic sample;
the warning is retained as a technical-quality flag, not silently converted to
`ok`.

## CRC readiness profile

1. **Case density and precision:** the assay-specific rebuilt outcome QC contains {fmt_int(crc['outcome_case_n_total'])} CRC cases across the pooled CRC-versus-cancer-free frame, with a median of {fmt_int(crc['analytic_case_n_median'])} complete-case CRC events per biomarker model. The median analytic CRC case fraction is {fmt_pct(crc['analytic_case_fraction_median'])}.
2. **Assay-specific analytic overlap:** the rebuilt screen uses the correct laboratory file and subsample weight for each assay family. Even after that correction, biomarker-specific merge and complete-case retention vary across tests; this is a measurable source of effective sample-size loss rather than a generic “NHANES N.”
3. **Cycle coverage:** the frozen panel has CRC assay coverage ranging from {int(crc['cycle_coverage_min'])} to {int(crc['cycle_coverage_max'])} cycles. The per-test table records the exact cycles and the number of cycles retaining complete cases.
4. **Outcome granularity:** CRC is defined as prevalent CRC versus cancer-free control. Diagnosis age is retained for {granularity.loc[granularity['outcome']=='CRC', 'diagnosis_age_available'].iloc[0]}, but diagnosis date, stage, detailed site, treatment, recurrence/progression, and prospective follow-up are unavailable in this screen. The diagnosis-age ledger is a separate provenance artifact and is not used to replace the assay-specific rebuilt outcome totals.
5. **Temporality:** the urine/serum measurement is contemporaneous with the survey and is not a prediagnostic biospecimen. Therefore the screen can identify an association under the frozen outcome definition, but cannot establish that exposure preceded CRC.

## T2D comparator profile

T2D has {fmt_int(t2d['outcome_case_n_total'])} pooled eligible outcome cases and a median of {fmt_int(t2d['analytic_case_n_median'])} cases per biomarker model. All 29 tests are estimable under their assay-specific frames. Its case definition uses diagnosed diabetes plus an objective HbA1c rule for probable undiagnosed disease, with indeterminate categories excluded from the primary case/control comparison. T2D remains cross-sectional and does not provide prospective exposure-to-onset timing.

## What this audit does and does not establish

This audit establishes that the CRC and T2D questions have materially different
readiness profiles: CRC is particularly constrained by case density, effective
case/sample overlap, and outcome temporality/granularity. It does **not** show
that any CRC association is absent, that a larger sample would guarantee an
association, or that the T2D findings are causal. Those questions belong to the
next negative-branch layer (failure attribution) and to external prospective
validation, respectively.

## Reproducibility and QC

- CRC results are from the assay-specific rebuilt 29-test screen, not the superseded phthalate-shaped frame.
- T2D results are from the frozen assay-specific 29-test screen.
- No CRC or T2D model was re-fit by this script.
- The 29-test BH-FDR results are included only as context; no threshold was used to define readiness.
- Per-test audit totals are compared with the recorded model-level analytic N and case counts. Any non-zero audit discrepancy is retained in the output for inspection rather than silently corrected.
- The outcome-granularity table explicitly records the timing and missing follow-up fields.
- The rebuilt assay-specific outcome QC and the legacy CRC case/control ledger have different row/case totals; they are not pooled. The rebuilt QC supplies primary readiness counts, while the legacy ledger is used only to document diagnosis-age availability.

Output files:

- `step9a_readiness_overview.csv`
- `step9a_per_test_readiness.csv`
- `step9a_outcome_granularity.csv`
- `STEP9A_MANIFEST.json`
"""


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    crc_results = read_csv(CRC_RESULTS)
    crc_merge = read_csv(CRC_MERGE)
    crc_qc = read_csv(CRC_OUTCOME_QC)
    crc_ledger = read_csv(CRC_LEDGER)
    t2d_results = read_csv(T2D_RESULTS)
    t2d_merge = read_csv(T2D_MERGE)
    t2d_qc = read_csv(T2D_OUTCOME_QC)

    require_columns(crc_results, ["test_id", "biomarker", "matrix", "frozen_cycle_list", "analytic_n", "analytic_crc_cases", "analytic_controls", "merge_N", "P", "BH_FDR"], "CRC results")
    require_columns(t2d_results, ["test_id", "biomarker", "matrix", "frozen_cycle_list", "analytic_n", "analytic_t2d_cases", "analytic_controls", "merge_N", "P", "BH_FDR"], "T2D results")

    crc_per_test = build_per_test("CRC", crc_results, crc_merge, crc_qc)
    t2d_per_test = build_per_test("T2D", t2d_results, t2d_merge, t2d_qc)
    per_test = pd.concat([t2d_per_test, crc_per_test], ignore_index=True)
    granularity = outcome_granularity(crc_ledger, crc_qc, t2d_qc)
    overview = build_overview(per_test, granularity, {"CRC": crc_qc, "T2D": t2d_qc})

    # Machine-check the frozen panel and model/audit reconciliation.
    if set(t2d_results["test_id"]) != set(crc_results["test_id"]):
        raise ValueError("CRC and T2D test universes do not match")
    if len(t2d_results) != 29 or len(crc_results) != 29:
        raise ValueError("Expected exactly 29 frozen tests in each outcome")

    per_test.to_csv(PER_TEST_OUT, index=False)
    overview.to_csv(OVERVIEW_OUT, index=False)
    granularity.to_csv(GRANULARITY_OUT, index=False)
    REPORT_OUT.write_text(report_text(overview, per_test, granularity), encoding="utf-8")

    inputs = [CRC_RESULTS, CRC_MERGE, CRC_OUTCOME_QC, CRC_LEDGER, T2D_RESULTS, T2D_MERGE, T2D_OUTCOME_QC]
    manifest = {
        "step": "9A",
        "title": "Question-specific data-readiness audit",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository_relative_inputs": [str(p.relative_to(ROOT)).replace("\\", "/") for p in inputs],
        "input_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in inputs},
        "outputs": [
            str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
            str(PER_TEST_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(OVERVIEW_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(GRANULARITY_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(REPORT_OUT.relative_to(ROOT)).replace("\\", "/"),
        ],
        "frozen_panel_n": 29,
        "scope_boundaries": {
            "refit_models": False,
            "change_fdr_family": False,
            "failure_attribution": False,
            "counterfactual_simulation": False,
            "readiness_score_collapsed": False,
        },
        "qc": {
            "crc_test_n": int(len(crc_results)),
            "t2d_test_n": int(len(t2d_results)),
            "crc_estimable_n": int((pd.to_numeric(crc_results["P"], errors="coerce").notna() & pd.to_numeric(crc_results["analytic_n"], errors="coerce").gt(0)).sum()),
            "t2d_estimable_n": int((pd.to_numeric(t2d_results["P"], errors="coerce").notna() & pd.to_numeric(t2d_results["analytic_n"], errors="coerce").gt(0)).sum()),
            "crc_audit_n_discrepancy_abs_max": float(pd.to_numeric(crc_per_test["audit_result_n_difference"], errors="coerce").abs().max()),
            "t2d_audit_n_discrepancy_abs_max": float(pd.to_numeric(t2d_per_test["audit_result_n_difference"], errors="coerce").abs().max()),
            "crc_audit_case_discrepancy_abs_max": float(pd.to_numeric(crc_per_test["audit_case_n_difference"], errors="coerce").abs().max()),
            "t2d_audit_case_discrepancy_abs_max": float(pd.to_numeric(t2d_per_test["audit_case_n_difference"], errors="coerce").abs().max()),
            "crc_rebuilt_outcome_qc_case_n": int(pd.to_numeric(crc_qc["crc_cases"], errors="coerce").sum()),
            "crc_legacy_ledger_case_n_secondary": int(pd.to_numeric(crc_ledger["crc_case"], errors="coerce").fillna(0).sum()),
            "crc_outcome_qc_and_legacy_ledger_not_pooled": True,
        },
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    crc_estimable = int((pd.to_numeric(crc_results["P"], errors="coerce").notna() & pd.to_numeric(crc_results["analytic_n"], errors="coerce").gt(0)).sum())
    t2d_estimable = int((pd.to_numeric(t2d_results["P"], errors="coerce").notna() & pd.to_numeric(t2d_results["analytic_n"], errors="coerce").gt(0)).sum())
    print(json.dumps({
        "status": "complete",
        "step": "9A",
        "outputs": [str(PER_TEST_OUT), str(OVERVIEW_OUT), str(GRANULARITY_OUT), str(REPORT_OUT), str(MANIFEST_OUT)],
        "crc_tests": len(crc_results),
        "t2d_tests": len(t2d_results),
        "crc_estimable": crc_estimable,
        "t2d_estimable": t2d_estimable,
        "crc_max_abs_audit_n_difference": manifest["qc"]["crc_audit_n_discrepancy_abs_max"],
        "t2d_max_abs_audit_n_difference": manifest["qc"]["t2d_audit_n_discrepancy_abs_max"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
