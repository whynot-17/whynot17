#!/usr/bin/env python3
"""Step 9B: descriptive failure attribution for the frozen CRC 29-test screen.

The script does not refit models.  It combines the rebuilt Step 5 primary
results, the Step 9A readiness table, and the already-frozen Step 6 robustness
artifacts to classify observable limitations.  Labels are intentionally
multi-valued: a test may be near-null, precision-limited, technically warned,
and interpretation/design-limited at the same time.
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
CRC_DIR = ROOT / "analysis/disease_agnostic_environmental_framework/step05_crc_screen_rebuilt"
ROBUST_DIR = ROOT / "analysis/disease_agnostic_environmental_framework/step06_robustness"

PRIMARY = CRC_DIR / "full_29_test_crc_screen_rebuilt.csv"
READINESS = OUT / "step9a_per_test_readiness.csv"
ROBUSTNESS = ROBUST_DIR / "robustness_results.csv"
LOCO = ROBUST_DIR / "loco_results.csv"
HETEROGENEITY = ROBUST_DIR / "cycle_heterogeneity.csv"
RUBRIC = ROBUST_DIR / "ROBUSTNESS_RUBRIC_LOCK.md"

ATTRIBUTION_OUT = OUT / "step9b_crc_failure_attribution_29_tests.csv"
SUMMARY_OUT = OUT / "step9b_crc_failure_attribution_summary.csv"
MANIFEST_OUT = OUT / "STEP9B_MANIFEST.json"
REPORT_OUT = OUT / "STEP9B_CRC_FAILURE_ATTRIBUTION_REPORT.md"

ALPHA = 0.05
TARGET_POWER = 0.80
MEANINGFUL_OR = 1.20
NEAR_NULL_OR_LOW = 1 / 1.10
NEAR_NULL_OR_HIGH = 1.10
Z_ALPHA = 1.959963984540054
Z_POWER = 0.8416212335729143


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def num(frame: pd.DataFrame, columns: list[str]) -> None:
    for col in columns:
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt_int(v) -> str:
    return "NA" if pd.isna(v) else f"{int(round(float(v))):,}"


def fmt_num(v, digits=3) -> str:
    return "NA" if pd.isna(v) else f"{float(v):.{digits}f}"


def fmt_pct(v, digits=1) -> str:
    return "NA" if pd.isna(v) else f"{100 * float(v):.{digits}f}%"


def join_labels(row: pd.Series) -> str:
    labels = []
    if bool(row["signal_limited_flag"]):
        labels.append("signal-limited")
    if bool(row["power_limited_nonnull_flag"]):
        labels.append("power-limited")
    if bool(row["stability_limited_flag"]):
        labels.append("stability-limited")
    if bool(row["technical_warning_flag"]):
        labels.append("technical-warning")
    if bool(row["design_limited_interpretation_flag"]):
        labels.append("design-limited-interpretation")
    return ";".join(labels) if labels else "none_from_current_descriptors"


def build_attribution(primary: pd.DataFrame, readiness: pd.DataFrame, robustness: pd.DataFrame, loco: pd.DataFrame, heterogeneity: pd.DataFrame) -> pd.DataFrame:
    primary = primary.copy()
    readiness = readiness.copy()
    robustness = robustness.copy()
    loco = loco.copy()
    heterogeneity = heterogeneity.copy()
    num(primary, ["P", "BH_FDR", "OR", "SE", "CI_low", "CI_high", "analytic_n", "analytic_crc_cases", "analytic_controls", "merge_N"])
    num(readiness, ["source_missing_fraction", "analytic_retention_from_merge", "analytic_n_result", "analytic_cases_result", "CI_width", "nominal_P", "BH_FDR"])
    num(robustness, ["F", "L", "C", "H", "D", "T", "A", "E", "primary_P", "primary_SE", "primary_fit_crc_cases", "primary_analytic_n", "primary_BH_FDR_29tests"])
    num(heterogeneity, ["interaction_P", "P", "fit_crc_cases"])
    num(loco, ["P", "fit_crc_cases", "OR"])

    # Step 6 was frozen to the two Step 5 FDR-supported signals.  A missing
    # join therefore means "not assessed under that frozen audit scope", not
    # "stable" or "unstable".
    robustness_key = "primary_biomarker"
    rob = robustness.drop_duplicates(robustness_key).set_index(robustness_key) if robustness_key in robustness else pd.DataFrame()
    het = heterogeneity[heterogeneity.get("analysis", "") == "cycle_interaction"].drop_duplicates("primary_biomarker").set_index("primary_biomarker") if len(heterogeneity) else pd.DataFrame()
    primary_beta_by_biomarker = primary.set_index("biomarker")["beta"].to_dict()
    loco_counts = {}
    for biomarker, group in loco.groupby("primary_biomarker", dropna=False):
        pooled_beta = pd.to_numeric(pd.Series([primary_beta_by_biomarker.get(biomarker, np.nan)]), errors="coerce").iloc[0]
        loco_beta = pd.to_numeric(group["beta"], errors="coerce")
        same_direction = int((np.sign(loco_beta) == np.sign(pooled_beta)).sum()) if not pd.isna(pooled_beta) else np.nan
        loco_counts[biomarker] = {"loco_n": int(len(group)), "loco_same_direction_n": same_direction}

    recs = []
    for _, row in primary.iterrows():
        biomarker = str(row["biomarker"])
        rr = readiness[readiness["biomarker"].astype(str) == biomarker]
        rr = rr.iloc[0] if len(rr) else pd.Series(dtype=object)
        rb = rob.loc[biomarker] if len(rob) and biomarker in rob.index else pd.Series(dtype=object)
        hh = het.loc[biomarker] if len(het) and biomarker in het.index else pd.Series(dtype=object)
        lc = pd.Series(loco_counts.get(biomarker, {}), dtype=object)

        beta = float(row["beta"]) if not pd.isna(row["beta"]) else np.nan
        se = float(row["SE"]) if not pd.isna(row["SE"]) else np.nan
        abs_log_or = abs(beta) if not pd.isna(beta) else np.nan
        mde_abs_log_or = (Z_ALPHA + Z_POWER) * se if not pd.isna(se) else np.nan
        mde_or_upper = math.exp(mde_abs_log_or) if not pd.isna(mde_abs_log_or) else np.nan
        mde_or_lower = math.exp(-mde_abs_log_or) if not pd.isna(mde_abs_log_or) else np.nan
        observed_or = float(row["OR"]) if not pd.isna(row["OR"]) else np.nan
        case_n = float(row["analytic_crc_cases"]) if not pd.isna(row["analytic_crc_cases"]) else np.nan
        fit_status = str(row.get("status", ""))

        signal_limited = not pd.isna(observed_or) and NEAR_NULL_OR_LOW <= observed_or <= NEAR_NULL_OR_HIGH
        precision_limited = not pd.isna(mde_or_upper) and (mde_or_upper > MEANINGFUL_OR or mde_or_lower < 1 / MEANINGFUL_OR)
        technical_warning = fit_status.lower() != "ok"
        stability_assessed = len(rb) > 0
        H = rb.get("H", np.nan) if len(rb) else np.nan
        L = rb.get("L", np.nan) if len(rb) else np.nan
        C = rb.get("C", np.nan) if len(rb) else np.nan
        stability_limited = stability_assessed and ((not pd.isna(H) and H <= 0) or (not pd.isna(L) and L <= 0) or (not pd.isna(C) and C <= 0))
        if pd.isna(case_n):
            event_support = "E0"
        elif case_n >= 60:
            event_support = "E2"
        elif case_n >= 30:
            event_support = "E1"
        else:
            event_support = "E0"

        out = {
            "test_id": row.get("test_id"),
            "biomarker": biomarker,
            "variable": row.get("variable"),
            "matrix": row.get("matrix"),
            "exposure_axis": row.get("exposure_axis"),
            "cycles_n_frozen": row.get("source_cycles_read", row.get("frozen_cycle_list", "")).count(";") + 1 if str(row.get("source_cycles_read", row.get("frozen_cycle_list", ""))) else 0,
            "analytic_n": float(row["analytic_n"]),
            "analytic_crc_cases": case_n,
            "analytic_controls": float(row["analytic_controls"]),
            "merge_n": float(row["merge_N"]),
            "analytic_case_fraction": case_n / float(row["analytic_n"]) if not pd.isna(case_n) and float(row["analytic_n"]) > 0 else np.nan,
            "analytic_retention_from_merge": rr.get("analytic_retention_from_merge", np.nan),
            "source_missing_fraction": rr.get("source_missing_fraction", np.nan),
            "OR_per_doubling": observed_or,
            "beta_logOR": beta,
            "SE_logOR": se,
            "abs_logOR": abs_log_or,
            "CI_low": row.get("CI_low"),
            "CI_high": row.get("CI_high"),
            "CI_width": row.get("CI_high") - row.get("CI_low") if not pd.isna(row.get("CI_high")) and not pd.isna(row.get("CI_low")) else np.nan,
            "nominal_P": row.get("P"),
            "BH_FDR_29": row.get("BH_FDR"),
            "observed_direction": "positive" if beta > 0 else ("inverse" if beta < 0 else "null"),
            "event_support_class": event_support,
            "mde_abs_logOR_approx_80pct": mde_abs_log_or,
            "mde_OR_upper_approx_80pct": mde_or_upper,
            "mde_OR_lower_approx_80pct": mde_or_lower,
            "mde_reference": "two-sided alpha=0.05, 80% power, normal approximation using current SE",
            "meaningful_effect_reference_OR": MEANINGFUL_OR,
            "signal_limited_flag": signal_limited,
            "precision_limited_20pct_flag": precision_limited,
            "power_limited_nonnull_flag": bool(precision_limited and not signal_limited),
            "technical_warning_flag": technical_warning,
            "fit_status": fit_status,
            "stability_assessment_status": "assessed_in_frozen_step06" if stability_assessed else "not_assessed_step06_scope_was_two_fdr_supported_signals",
            "loco_n": lc.get("loco_n", np.nan),
            "loco_same_direction_n": lc.get("loco_same_direction_n", np.nan),
            "loco_same_direction_fraction": (lc.get("loco_same_direction_n", np.nan) / lc.get("loco_n", np.nan)) if lc.get("loco_n", np.nan) else np.nan,
            "robustness_L": L,
            "robustness_C": C,
            "heterogeneity_H": H,
            "cycle_interaction_P": hh.get("interaction_P", np.nan),
            "heterogeneity_status": hh.get("heterogeneity_status", np.nan),
            "stability_limited_flag": stability_limited,
            "design_limited_interpretation_flag": True,
            "design_limitation_labels": "prevalent_cross_sectional;no_prediagnostic_biospecimen;no_stage_site_treatment_recurrence_followup",
            "failure_labels": "",
        }
        recs.append(out)

    result = pd.DataFrame(recs)
    result["failure_labels"] = result.apply(join_labels, axis=1)
    result["primary_statistical_pattern"] = np.select(
        [result["signal_limited_flag"], result["power_limited_nonnull_flag"]],
        ["observed_near_null", "non-null_direction_but_20pct_MDE_large"],
        default="not_classified_by_signal_precision_rule",
    )
    return result


def summary_table(attr: pd.DataFrame) -> pd.DataFrame:
    def n(col):
        return int(attr[col].fillna(False).astype(bool).sum())
    rows = [
        {"metric": "frozen_crc_tests", "value": len(attr), "unit": "tests"},
        {"metric": "nominal_p_lt_0_05", "value": int((pd.to_numeric(attr["nominal_P"], errors="coerce") < 0.05).sum()), "unit": "tests"},
        {"metric": "bh_fdr_lt_0_05", "value": int((pd.to_numeric(attr["BH_FDR_29"], errors="coerce") < 0.05).sum()), "unit": "tests"},
        {"metric": "signal_limited_observed_near_null", "value": n("signal_limited_flag"), "unit": "tests"},
        {"metric": "mde_exceeds_20pct_reference_all_tests", "value": n("precision_limited_20pct_flag"), "unit": "tests"},
        {"metric": "power_limited_nonnull_observed_signal", "value": n("power_limited_nonnull_flag"), "unit": "tests"},
        {"metric": "technical_warning", "value": n("technical_warning_flag"), "unit": "tests"},
        {"metric": "stability_assessed_under_frozen_step06", "value": int((attr["stability_assessment_status"] == "assessed_in_frozen_step06").sum()), "unit": "tests"},
        {"metric": "stability_limited_among_assessed", "value": n("stability_limited_flag"), "unit": "tests"},
        {"metric": "event_support_E2_ge_60_cases", "value": int((attr["event_support_class"] == "E2").sum()), "unit": "tests"},
        {"metric": "event_support_E1_30_to_59_cases", "value": int((attr["event_support_class"] == "E1").sum()), "unit": "tests"},
        {"metric": "event_support_E0_lt_30_cases", "value": int((attr["event_support_class"] == "E0").sum()), "unit": "tests"},
        {"metric": "design_limited_interpretation", "value": n("design_limited_interpretation_flag"), "unit": "tests; interpretation limitation not a significance cause"},
    ]
    return pd.DataFrame(rows)


def report_text(attr: pd.DataFrame, summary: pd.DataFrame) -> str:
    def sval(metric):
        return int(summary.loc[summary["metric"] == metric, "value"].iloc[0])
    return f"""# Step 9B — CRC negative-branch failure attribution

## Scope

This analysis diagnoses the frozen 29-test CRC screen after the assay-specific
rebuild. It does not re-fit any model, alter the 29-test BH-FDR family, or
retroactively select candidates. It uses the recorded primary estimate and
standard error, the Step 9A merge/readiness fields, and the Step 6 robustness
artifacts. Step 6 was frozen to the two CRC FDR-supported signals, so stability
is **not** imputed for the other 27 tests.

## Frozen attribution rules

The labels are intentionally multi-valued.

- **Signal-limited:** observed OR is within 0.90–1.10, a descriptive near-null rule.
- **Power-limited:** the observed OR is outside 0.90–1.10 and the approximate two-sided alpha=0.05, 80%-power MDE derived from the current SE exceeds OR=1.20 in either direction. Near-null tests still retain their continuous MDE, but are not called power-limited merely because their MDE is wide.
- **Stability-limited:** among tests actually audited by the frozen Step 6 scope, the recorded L/C score is 0 or the H heterogeneity tag is 0. Missing Step 6 coverage is reported as not assessed.
- **Technical-warning:** the recorded primary status is not `ok`.
- **Design-limited-interpretation:** all CRC tests carry the same limitation: prevalent cross-sectional outcome, no prediagnostic biospecimen, and no stage/site/treatment/recurrence/follow-up fields. This label limits interpretation and temporality; it is not claimed to explain a P value.
- **Event support:** E2 = >=60 analytic CRC cases, E1 = 30–59, E0 = <30, following the frozen Step 6 event-support rubric.

## Summary

| Metric | Count |
|---|---:|
| Frozen CRC tests | {len(attr)} |
| Nominal P<0.05 | {sval('nominal_p_lt_0_05')} |
| BH-FDR<0.05 | {sval('bh_fdr_lt_0_05')} |
| Observed near-null / signal-limited | {sval('signal_limited_observed_near_null')} |
| All tests with 20% MDE above reference | {sval('mde_exceeds_20pct_reference_all_tests')} |
| Power-limited among non-near-null tests | {sval('power_limited_nonnull_observed_signal')} |
| Primary technical warnings | {sval('technical_warning')} |
| Stability assessed under frozen Step 6 scope | {sval('stability_assessed_under_frozen_step06')} |
| Stability-limited among assessed tests | {sval('stability_limited_among_assessed')} |
| E2 (>=60 cases) | {sval('event_support_E2_ge_60_cases')} |
| E1 (30–59 cases) | {sval('event_support_E1_30_to_59_cases')} |
| E0 (<30 cases) | {sval('event_support_E0_lt_30_cases')} |

## Interpretation

The CRC screen has complete model-level estimability for all 29 tests, but the
event-support profile is much weaker than T2D. The attribution table separates
three different statements that should not be collapsed: (i) an observed effect
can be near the null, (ii) the current SE may be too large to detect a prespecified
20% effect reliably, and (iii) even a statistically precise cross-sectional
association would still have limited temporality and outcome granularity.

The MDE column is therefore a sensitivity descriptor. It does not prove that the
CRC negative screen is caused by low power, and it does not imply that adding
cases would guarantee a discovery. Conversely, a test that is not flagged
power-limited is not thereby proven biologically null.

Step 6 stability results are not generalized to all 29 tests: only the two
FDR-supported CRC signals were in that frozen audit. This prevents the negative
branch from manufacturing stability claims by absence of analysis.

## Machine-readable outputs

- `step9b_crc_failure_attribution_29_tests.csv`
- `step9b_crc_failure_attribution_summary.csv`
- `STEP9B_MANIFEST.json`

The complete per-test table retains the observed OR, SE, CI, analytic cases,
retention, approximate MDE, event-support class, available stability fields,
and all applicable attribution labels.
"""


def main() -> int:
    primary = read_csv(PRIMARY)
    readiness = read_csv(READINESS)
    robustness = read_csv(ROBUSTNESS)
    loco = read_csv(LOCO)
    heterogeneity = read_csv(HETEROGENEITY)

    if len(primary) != 29 or primary["test_id"].nunique() != 29:
        raise ValueError("CRC primary input must contain exactly 29 unique frozen tests")
    if len(readiness[readiness["outcome"] == "CRC"]) != 29:
        raise ValueError("Step 9A readiness input must contain 29 CRC rows")

    attr = build_attribution(primary, readiness[readiness["outcome"] == "CRC"], robustness, loco, heterogeneity)
    summary = summary_table(attr)
    ATTRIBUTION_OUT.parent.mkdir(parents=True, exist_ok=True)
    attr.to_csv(ATTRIBUTION_OUT, index=False)
    summary.to_csv(SUMMARY_OUT, index=False)
    REPORT_OUT.write_text(report_text(attr, summary), encoding="utf-8")

    inputs = [PRIMARY, READINESS, ROBUSTNESS, LOCO, HETEROGENEITY, RUBRIC]
    manifest = {
        "step": "9B",
        "title": "CRC negative-branch failure attribution",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in inputs},
        "outputs": [
            str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
            str(ATTRIBUTION_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(SUMMARY_OUT.relative_to(ROOT)).replace("\\", "/"),
            str(REPORT_OUT.relative_to(ROOT)).replace("\\", "/"),
        ],
        "frozen_rules": {
            "tests": 29,
            "near_null_or_interval": [NEAR_NULL_OR_LOW, NEAR_NULL_OR_HIGH],
            "meaningful_or_reference": MEANINGFUL_OR,
            "alpha": ALPHA,
            "target_power": TARGET_POWER,
            "mde_formula": "(z_(1-alpha/2)+z_power)*SE on log(OR), exponentiated",
            "power_limited_requires_non_near_null_observed_or": True,
            "stability_scope": "Step 6 frozen two FDR-supported CRC signals only",
            "design_flag_is_interpretation_limitation_not_statistical_failure_cause": True,
            "refit_models": False,
            "change_fdr_family": False,
            "counterfactual_simulation": False,
        },
        "qc": {
            "primary_tests": int(len(primary)),
            "attribution_rows": int(len(attr)),
            "all_primary_tests_joined": bool(set(primary["test_id"]) == set(attr["test_id"])),
            "step9a_crc_rows": int(len(readiness[readiness["outcome"] == "CRC"])),
            "step6_stability_rows": int((attr["stability_assessment_status"] == "assessed_in_frozen_step06").sum()),
        },
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "step": "9B", "rows": len(attr), "summary": summary.to_dict(orient="records")}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
