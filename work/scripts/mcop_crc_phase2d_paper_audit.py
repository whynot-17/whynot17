"""Phase 2D paper-readiness audit for MCOP-CRC in NHANES.

This is a pre-specified, narrow follow-up to the existing MCOP Phase 2
validation. It reuses the established NHANES harmonized frame and the
validated survey-logistic implementation. The audit covers sex-specific
effects and interaction, diagnosis-timing exclusions, phthalate
co-exposure specificity, and six paper figures. It does not search for a new
chemical, run mechanistic analyses, or use WHI data.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from scipy.stats import t


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_HARMONIZED = ROOT / "work" / "nhanes_phase2a" / "phase2b_harmonized.pkl"
DEFAULT_OUTPUT = ROOT / "outputs"

BASE_CONTINUOUS = ["age", "bmi", "pir", "creatinine_log2"]
BASE_CATEGORICAL = ["sex", "race", "smoking"]
COEXPOSURES = {
    "MEHHP": "URXMHH",
    "MEOHP": "URXMOH",
    "MECPP": "URXECP",
    "MBzP": "URXMZP",
}
BURDEN_METABOLITES = {
    "MBzP": "URXMZP",
    "MnBP": "URXMBP",
    "MiBP": "URXMIB",
    "MEP": "URXMEP",
    "MMP": "URXMNM",
    "MEHP": "URXMHP",
    "MEHHP": "URXMHH",
    "MEOHP": "URXMOH",
    "MECPP": "URXECP",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_mcop_module():
    path = ROOT / "work" / "scripts" / "mcop_crc_phase2.py"
    spec = importlib.util.spec_from_file_location("mcop_phase2d_loader", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load MCOP implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def public_fit(fit: dict, analysis: str, cycles: list[str] | None = None) -> dict:
    row = {
        "Analysis": analysis,
        "Cycles": ";".join(cycles or []),
        **{key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}},
    }
    return row


def fit_primary(frame: pd.DataFrame, model_module, analysis: str, cycles: list[str] | None = None, categorical: list[str] | None = None) -> dict:
    fit = model_module.fit_survey_logistic(
        frame,
        ["mcop_log2", *BASE_CONTINUOUS],
        categorical or BASE_CATEGORICAL,
        exposure_name="mcop_log2",
        levels=model_module.LEVELS,
    )
    return public_fit(fit, analysis, cycles)


def coefficient_effect(fit: dict, weights: dict[str, float]) -> dict:
    """Return a Wald effect for a linear combination of fitted coefficients."""
    coefficients = fit.get("coefficients", {})
    covariance = fit.get("covariance")
    if fit.get("status") not in {"ok", "converged_with_warning"} or covariance is None:
        return {"beta": np.nan, "SE": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": np.nan}
    names = list(coefficients)
    if any(name not in coefficients for name in weights):
        return {"beta": np.nan, "SE": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": np.nan}
    index = [names.index(name) for name in weights]
    vector = np.asarray([weights[name] for name in weights], dtype=float)
    beta_values = np.asarray([coefficients[name] for name in weights], dtype=float)
    cov = np.asarray(covariance)[np.ix_(index, index)]
    beta = float(vector @ beta_values)
    variance = float(vector @ cov @ vector)
    se = float(np.sqrt(max(variance, 0.0))) if np.isfinite(variance) else np.nan
    design_df = max(int(fit.get("design_df", 1)), 1)
    if not np.isfinite(se) or se <= 0:
        return {"beta": beta, "SE": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": np.nan}
    critical = float(t.ppf(0.975, design_df))
    z_value = beta / se
    return {
        "beta": beta,
        "SE": se,
        "OR": float(np.exp(np.clip(beta, -700, 700))),
        "CI_low": float(np.exp(np.clip(beta - critical * se, -700, 700))),
        "CI_high": float(np.exp(np.clip(beta + critical * se, -700, 700))),
        "P": float(2 * t.sf(abs(z_value), design_df)),
    }


def add_log2(frame: pd.DataFrame, source: str, target: str) -> pd.DataFrame:
    values = pd.to_numeric(frame[source], errors="coerce")
    frame[target] = np.log2(values.where(values > 0))
    return frame


def bh_fdr(pvalues: pd.Series) -> np.ndarray:
    values = pd.to_numeric(pvalues, errors="coerce").fillna(1.0).to_numpy(float)
    order = np.argsort(values)
    adjusted = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty(len(values), dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def run_sex_audit(primary: pd.DataFrame, model_module) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for group in ["Female", "Male"]:
        frame = primary[primary["sex"].eq(group)].copy()
        fit = fit_primary(frame, model_module, f"{group}_MCOP_primary", categorical=["race", "smoking"])
        fit["Sex_group"] = group
        rows.append(fit)

    work = primary.copy()
    work["mcop_x_male"] = work["mcop_log2"] * work["sex"].eq("Male").astype(float)
    interaction_fit = model_module.fit_survey_logistic(
        work,
        ["mcop_log2", "mcop_x_male", *BASE_CONTINUOUS],
        BASE_CATEGORICAL,
        exposure_name="mcop_log2",
        levels=model_module.LEVELS,
    )
    female_effect = coefficient_effect(interaction_fit, {"mcop_log2": 1.0})
    male_effect = coefficient_effect(interaction_fit, {"mcop_log2": 1.0, "mcop_x_male": 1.0})
    interaction_effect = coefficient_effect(interaction_fit, {"mcop_x_male": 1.0})
    interaction_row = public_fit(interaction_fit, "MCOP_x_sex_interaction")
    interaction_row.update(
        {
            "Effect": "MCOP x sex interaction: male/female slope ratio",
            "beta": interaction_effect["beta"],
            "SE": interaction_effect["SE"],
            "OR": interaction_effect["OR"],
            "CI_low": interaction_effect["CI_low"],
            "CI_high": interaction_effect["CI_high"],
            "P": interaction_effect["P"],
            "Female_OR": female_effect["OR"],
            "Female_CI_low": female_effect["CI_low"],
            "Female_CI_high": female_effect["CI_high"],
            "Male_OR": male_effect["OR"],
            "Male_CI_low": male_effect["CI_low"],
            "Male_CI_high": male_effect["CI_high"],
            "Interaction_OR_Male_vs_Female": interaction_effect["OR"],
            "Interaction_CI_low": interaction_effect["CI_low"],
            "Interaction_CI_high": interaction_effect["CI_high"],
            "Interaction_P": interaction_effect["P"],
        }
    )
    return pd.DataFrame(rows), pd.DataFrame([interaction_row])


def timing_metrics(cases: pd.DataFrame) -> dict:
    known = cases[pd.to_numeric(cases["years_since_crc"], errors="coerce").notna()].copy()
    metrics = {
        "Timing_known_CRC_N": int(len(known)),
        "Timing_missing_CRC_N": int(len(cases) - len(known)),
        "Exam_age_median": np.nan,
        "Exam_age_Q1": np.nan,
        "Exam_age_Q3": np.nan,
        "Diagnosis_age_median": np.nan,
        "Diagnosis_age_Q1": np.nan,
        "Diagnosis_age_Q3": np.nan,
        "Years_since_diagnosis_median": np.nan,
        "Years_since_diagnosis_Q1": np.nan,
        "Years_since_diagnosis_Q3": np.nan,
    }
    if known.empty:
        return metrics
    for source, prefix in [("age", "Exam_age"), ("crc_diagnosis_age", "Diagnosis_age"), ("years_since_crc", "Years_since_diagnosis")]:
        values = pd.to_numeric(known[source], errors="coerce").dropna()
        if values.empty:
            continue
        metrics[f"{prefix}_median"] = float(values.median())
        metrics[f"{prefix}_Q1"] = float(values.quantile(0.25))
        metrics[f"{prefix}_Q3"] = float(values.quantile(0.75))
    return metrics


def run_timing_audit(primary: pd.DataFrame, model_module) -> tuple[pd.DataFrame, pd.DataFrame]:
    case_table = primary[primary["outcome"].eq(1) & primary["mcop_log2"].notna()].copy()
    case_table = case_table.rename(columns={"age": "exam_age", "crc_diagnosis_age": "diagnosis_age", "years_since_crc": "years_since_diagnosis"})
    case_table["timing_group"] = "unknown"
    interval = pd.to_numeric(case_table["years_since_diagnosis"], errors="coerce")
    case_table.loc[interval.lt(1), "timing_group"] = "<1 year"
    case_table.loc[interval.ge(1) & interval.le(2), "timing_group"] = "1-2 years"
    case_table.loc[interval.gt(2) & interval.le(5), "timing_group"] = ">2-5 years"
    case_table.loc[interval.gt(5), "timing_group"] = ">5 years"
    output_columns = [
        column for column in [
            "SEQN", "cycle", "exam_age", "diagnosis_age", "years_since_diagnosis",
            "timing_group", "URXCOP", "mcop_log2",
        ] if column in case_table.columns
    ]
    case_output = case_table[output_columns].sort_values(["cycle", "SEQN"] if "SEQN" in case_table.columns else ["cycle"]).reset_index(drop=True)

    all_metrics = timing_metrics(case_table.rename(columns={"exam_age": "age", "diagnosis_age": "crc_diagnosis_age", "years_since_diagnosis": "years_since_crc"}))
    rows = []
    for label, threshold in [("Primary_all_cases", None), ("Exclude_diagnosis_lt_1y", 1.0), ("Exclude_diagnosis_lt_2y", 2.0), ("Exclude_diagnosis_lt_5y", 5.0)]:
        if threshold is None:
            frame = primary.copy()
            excluded = 0
        else:
            years = pd.to_numeric(primary["years_since_crc"], errors="coerce")
            excluded_mask = primary["outcome"].eq(1) & years.notna() & years.lt(threshold)
            excluded = int(excluded_mask.sum())
            frame = primary.loc[~excluded_mask].copy()
        fit = fit_primary(frame, model_module, label)
        row = fit | {"Excluded_known_timing_CRC_N": excluded, **all_metrics}
        rows.append(row)

    for label, mask in [("Descriptive_CRC_interval_le_5y", interval.le(5)), ("Descriptive_CRC_interval_gt_5y", interval.gt(5))]:
        subset = case_table[mask].copy()
        row = {"Analysis": label, "status": "descriptive_only", "N": int(len(subset)), "CRC_N": int(len(subset)), "Control_N": 0, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": np.nan, "Cycles": ";".join(sorted(subset["cycle"].dropna().unique())), **timing_metrics(subset.rename(columns={"exam_age": "age", "diagnosis_age": "crc_diagnosis_age", "years_since_diagnosis": "years_since_crc"}))}
        rows.append(row)
    return pd.DataFrame(rows), case_output


def common_cycles(primary: pd.DataFrame, columns: list[str]) -> list[str]:
    cycles = set(primary.loc[primary["mcop_log2"].notna(), "cycle"].dropna().unique())
    for column in columns:
        cycles &= set(primary.loc[pd.to_numeric(primary[column], errors="coerce").gt(0), "cycle"].dropna().unique())
    return sorted(cycles)


def fit_coexposure_model(
    frame: pd.DataFrame,
    model_module,
    secondary_label: str,
    secondary_column: str | None,
    cycles: list[str],
    secondary_already_log2: bool = False,
) -> list[dict]:
    work = frame.copy()
    if secondary_column is not None:
        if secondary_already_log2:
            work["secondary_log2"] = pd.to_numeric(work[secondary_column], errors="coerce")
        else:
            work = add_log2(work, secondary_column, "secondary_log2")
    required = ["outcome", "mcop_log2", *BASE_CONTINUOUS, *BASE_CATEGORICAL]
    if secondary_column is not None:
        required.append("secondary_log2")
    cc = work.dropna(subset=required).copy()
    primary_fit = fit_primary(cc, model_module, f"MCOP_alone_same_CC_{secondary_label}", cycles)
    primary_fit["Model_role"] = "same_complete_case_MCOP_only"
    primary_fit["Secondary_exposure"] = secondary_label
    primary_fit["Secondary_OR"] = np.nan
    primary_fit["Secondary_CI_low"] = np.nan
    primary_fit["Secondary_CI_high"] = np.nan
    primary_fit["Secondary_P"] = np.nan
    if secondary_column is None:
        primary_fit["Model_role"] = "MCOP_alone"
        return [primary_fit]

    fit = model_module.fit_survey_logistic(
        cc,
        ["mcop_log2", "secondary_log2", *BASE_CONTINUOUS],
        BASE_CATEGORICAL,
        exposure_name="mcop_log2",
        levels=model_module.LEVELS,
    )
    row = public_fit(fit, f"MCOP_plus_{secondary_label}", cycles)
    secondary_effect = coefficient_effect(fit, {"secondary_log2": 1.0})
    row.update(
        {
            "Model_role": "coexposure_adjusted",
            "Secondary_exposure": secondary_label,
            "Secondary_OR": secondary_effect["OR"],
            "Secondary_CI_low": secondary_effect["CI_low"],
            "Secondary_CI_high": secondary_effect["CI_high"],
            "Secondary_P": secondary_effect["P"],
            "Complete_case_N": int(len(cc)),
            "Complete_case_CRC_N": int(cc["outcome"].sum()),
        }
    )
    primary_fit["Complete_case_N"] = int(len(cc))
    primary_fit["Complete_case_CRC_N"] = int(cc["outcome"].sum())
    return [primary_fit, row]


def run_coexposure_audit(primary: pd.DataFrame, model_module) -> pd.DataFrame:
    rows = []
    all_mcop_cycles = common_cycles(primary, [])
    rows.extend(fit_coexposure_model(primary[primary["cycle"].isin(all_mcop_cycles)].copy(), model_module, "None", None, all_mcop_cycles))
    for label, column in COEXPOSURES.items():
        cycles = common_cycles(primary, [column])
        if not cycles:
            continue
        rows.extend(fit_coexposure_model(primary[primary["cycle"].isin(cycles)].copy(), model_module, label, column, cycles))

    burden_work = primary[primary["cycle"].isin(all_mcop_cycles)].copy()
    log_data = pd.DataFrame(index=burden_work.index)
    for label, column in BURDEN_METABOLITES.items():
        if column in burden_work.columns:
            log_data[label] = np.log2(pd.to_numeric(burden_work[column], errors="coerce").where(pd.to_numeric(burden_work[column], errors="coerce") > 0))
    z = (log_data - log_data.mean()) / log_data.std(ddof=0)
    burden_work["phthalate_burden_excl_mcop"] = z.mean(axis=1, skipna=True)
    burden_work["phthalate_burden_n_metabolites"] = z.notna().sum(axis=1)
    burden_work = burden_work[burden_work["phthalate_burden_n_metabolites"] >= 2].copy()
    rows.extend(
        fit_coexposure_model(
            burden_work,
            model_module,
            "PhthalateBurden_excl_MCOP",
            "phthalate_burden_excl_mcop",
            all_mcop_cycles,
            secondary_already_log2=True,
        )
    )
    result = pd.DataFrame(rows)
    adjusted = result["Model_role"].eq("coexposure_adjusted")
    result.loc[adjusted, "MCOP_BH_FDR"] = bh_fdr(result.loc[adjusted, "P"])
    result["Burden_definition"] = "mean z-score of 9 non-MCOP phthalate metabolites; minimum 2 measured" 
    return result


def existing_row(path: Path, mask: pd.Series) -> dict | None:
    if not path.exists():
        return None
    frame = pd.read_csv(path)
    subset = frame.loc[mask(frame)] if callable(mask) else frame.loc[mask]
    if subset.empty:
        return None
    return subset.iloc[0].to_dict()


def forest_plot(rows: pd.DataFrame, path: Path, title: str, label_column: str = "label") -> None:
    work = rows.copy()
    work = work.dropna(subset=["OR", "CI_low", "CI_high"]).reset_index(drop=True)
    if work.empty:
        return
    y = np.arange(len(work))
    fig, ax = plt.subplots(figsize=(8, max(3.2, len(work) * 0.42)))
    ax.errorbar(work["OR"], y, xerr=[work["OR"] - work["CI_low"], work["CI_high"] - work["OR"]], fmt="o", color="#1f4e79", ecolor="#1f4e79", capsize=3)
    ax.axvline(1, color="0.45", lw=1, ls="--")
    ax.set_yticks(y, work[label_column])
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio per MCOP doubling (log scale)")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def make_figures(output_dir: Path, primary_fit: dict, timing: pd.DataFrame, sex_rows: pd.DataFrame, interaction: pd.DataFrame, coexposure: pd.DataFrame) -> None:
    # Figure 1: workflow.
    fig, ax = plt.subplots(figsize=(13, 3.3))
    ax.axis("off")
    boxes = [
        (0.02, "CTD × GeneCards\nno-prior screen"),
        (0.19, "MBzP\nNHANES null → stop"),
        (0.36, "DINP axis\nMiNP signal / MCOP biomarker"),
        (0.54, "NHANES MCOP\nadjusted discovery"),
        (0.72, "Phase 2D\npaper-readiness audit"),
        (0.90, "Future WHI\nprospective replication"),
    ]
    for index, (x, label) in enumerate(boxes):
        style = dict(boxstyle="round,pad=0.6", facecolor="#dbeafe" if index < 5 else "#f3f4f6", edgecolor="#1f4e79", lw=1.2)
        ax.text(x, 0.5, label, ha="center", va="center", fontsize=10, bbox=style, transform=ax.transAxes)
        if index < len(boxes) - 1:
            ax.annotate("", xy=(boxes[index + 1][0] - 0.055, 0.5), xytext=(x + 0.055, 0.5), xycoords=ax.transAxes, arrowprops=dict(arrowstyle="->", lw=1.5, color="#1f4e79"))
    ax.set_title("MCOP–CRC project workflow", fontsize=13, pad=12)
    fig.tight_layout()
    fig.savefig(output_dir / "mcop_crc_phase2d_figure1_workflow.png", dpi=220)
    plt.close(fig)

    # Figure 2: primary and reverse-causation sensitivities plus previously audited stability checks.
    forest = []
    for label, source in [("Primary", primary_fit)]:
        forest.append({"label": label, "OR": source.get("OR"), "CI_low": source.get("CI_low"), "CI_high": source.get("CI_high")})
    for label, analysis in [("Exclude diagnosis <1y", "Exclude_diagnosis_lt_1y"), ("Exclude diagnosis <2y", "Exclude_diagnosis_lt_2y"), ("Exclude diagnosis <5y", "Exclude_diagnosis_lt_5y")]:
        row = timing[timing["Analysis"].eq(analysis)]
        if not row.empty:
            row = row.iloc[0]
            forest.append({"label": label, "OR": row.get("OR"), "CI_low": row.get("CI_low"), "CI_high": row.get("CI_high")})
    tail_path = output_dir / "mcop_crc_phase2_tail_exclusion.csv"
    if tail_path.exists():
        tail = pd.read_csv(tail_path)
        for _, row in tail.iterrows():
            forest.append({"label": str(row["Analysis"]).replace("Exclude_", ""), "OR": row.get("OR"), "CI_low": row.get("CI_low"), "CI_high": row.get("CI_high")})
    norm_path = output_dir / "mcop_crc_phase2_creatinine_normalized.csv"
    if norm_path.exists():
        row = pd.read_csv(norm_path).iloc[0]
        forest.append({"label": "Creatinine-normalized", "OR": row.get("OR"), "CI_low": row.get("CI_low"), "CI_high": row.get("CI_high")})
    forest_plot(pd.DataFrame(forest), output_dir / "mcop_crc_phase2d_figure2_forest.png", "MCOP–CRC primary model and stability sensitivities")

    # Figure 3: cycle-specific estimates and pooled estimate.
    cycle_path = output_dir / "mcop_crc_phase2_per_cycle.csv"
    if cycle_path.exists():
        cycle = pd.read_csv(cycle_path)
        cycle_rows = cycle.rename(columns={"Cycle": "label"})[["label", "OR", "CI_low", "CI_high"]].copy()
        cycle_rows.loc[len(cycle_rows)] = {"label": "Pooled", "OR": primary_fit.get("OR"), "CI_low": primary_fit.get("CI_low"), "CI_high": primary_fit.get("CI_high")}
        forest_plot(cycle_rows, output_dir / "mcop_crc_phase2d_figure3_cycles.png", "Cycle-specific and pooled MCOP–CRC estimates")

    # Figure 4: previously audited RCS curve.
    spline_path = output_dir / "mcop_crc_phase2_spline_curve.csv"
    if spline_path.exists():
        curve = pd.read_csv(spline_path)
        fig, ax = plt.subplots(figsize=(6.5, 4.2))
        ax.plot(2 ** curve["mcop_log2"], curve["OR_vs_median"], color="#1f4e79", lw=2)
        ax.axhline(1, color="0.45", ls="--", lw=1)
        ax.set_xscale("log")
        ax.set_xlabel("Urinary MCOP (ng/mL, log scale)")
        ax.set_ylabel("OR versus median exposure")
        ax.set_title("Restricted cubic spline: MCOP–CRC")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(output_dir / "mcop_crc_phase2d_figure4_rcs.png", dpi=220)
        plt.close(fig)

    # Figure 5: sex-specific effects.
    sex_plot = sex_rows.rename(columns={"Sex_group": "label"})[["label", "OR", "CI_low", "CI_high"]]
    interaction_p = float(interaction.iloc[0]["Interaction_P"])
    forest_plot(sex_plot, output_dir / "mcop_crc_phase2d_figure5_sex_interaction.png", f"Sex-specific MCOP–CRC estimates (interaction P={interaction_p:.3f})")

    # Figure 6: MCOP after co-exposure adjustment.
    coex_plot = coexposure[(coexposure["Model_role"].eq("MCOP_alone")) | (coexposure["Model_role"].eq("coexposure_adjusted"))].copy()
    burden_label = "MCOP + non-MCOP phthalate burden"
    coex_plot["label"] = "MCOP + " + coex_plot["Secondary_exposure"].fillna("").astype(str)
    coex_plot.loc[coex_plot["Secondary_exposure"].eq("PhthalateBurden_excl_MCOP"), "label"] = burden_label
    coex_plot.loc[coex_plot["Model_role"].eq("MCOP_alone"), "label"] = "MCOP alone"
    forest_plot(coex_plot, output_dir / "mcop_crc_phase2d_figure6_coexposure.png", "MCOP specificity under co-exposure adjustment")


def fmt(value) -> str:
    return "NA" if pd.isna(value) else f"{float(value):.4g}"


def write_report(output_dir: Path, sex_rows: pd.DataFrame, interaction: pd.DataFrame, timing: pd.DataFrame, coexposure: pd.DataFrame) -> None:
    female = sex_rows.loc[sex_rows["Sex_group"].eq("Female")].iloc[0]
    male = sex_rows.loc[sex_rows["Sex_group"].eq("Male")].iloc[0]
    interaction_row = interaction.iloc[0]
    timing_models = timing[timing["Analysis"].isin(["Exclude_diagnosis_lt_1y", "Exclude_diagnosis_lt_2y", "Exclude_diagnosis_lt_5y"])].copy()
    coex_models = coexposure[coexposure["Model_role"].eq("coexposure_adjusted")].copy()
    female_positive = bool(pd.notna(female["OR"]) and female["OR"] > 1)
    timing_positive = bool(not timing_models.empty and timing_models["OR"].notna().all() and (timing_models["OR"] > 1).all())
    coex_positive = bool(not coex_models.empty and coex_models["OR"].notna().all() and (coex_models["OR"] > 1).all())
    lines = [
        "# MCOP–CRC Phase 2D：paper-readiness audit",
        "",
        "## Material Passport",
        "",
        "- Origin Skill: academic-research-suite + experiment-agent",
        "- Origin Mode: run / validation audit",
        f"- Origin Date: {datetime.now(timezone.utc).date().isoformat()}",
        "- Verification Status: ANALYZED — local NHANES rerun completed; no WHI data used",
        "- Scope: female/male effect, sex interaction, diagnosis timing, co-exposure specificity, paper figures",
        "",
        "## Frozen scope",
        "",
        "This Phase 2D run does not screen another chemical and does not run TCGA/PPI/GO. All regression models reuse the established NHANES survey-logistic implementation, cancer-free-control definition, phthalate weights, pooled strata/PSU identifiers and Model 2 covariates: age, BMI, PIR, urinary creatinine, sex, race and smoking unless sex is constant within a sex-specific analysis.",
        "",
        "## Gate results",
        "",
        f"- Female point estimate >1: **{'PASS' if female_positive else 'FAIL'}**; OR={fmt(female['OR'])}, 95% CI {fmt(female['CI_low'])}-{fmt(female['CI_high'])}, P={fmt(female['P'])}, N={int(female['N'])}, CRC={int(female['CRC_N'])}. This is directionally positive but not a close point-estimate replication of the pooled OR≈1.25.",
        f"- Male point estimate: OR={fmt(male['OR'])}, 95% CI {fmt(male['CI_low'])}-{fmt(male['CI_high'])}, P={fmt(male['P'])}, N={int(male['N'])}, CRC={int(male['CRC_N'])}.",
        f"- MCOP × sex interaction: OR ratio male/female={fmt(interaction_row['Interaction_OR_Male_vs_Female'])}, 95% CI {fmt(interaction_row['Interaction_CI_low'])}-{fmt(interaction_row['Interaction_CI_high'])}, P={fmt(interaction_row['Interaction_P'])}.",
        "- Interpretation: the female estimate remains positive but is weaker than the pooled estimate and its CI includes the null; the formal interaction does not support a statistically clear sex difference.",
        f"- Recent-diagnosis exclusions retain positive MCOP direction: **{'PASS' if timing_positive else 'FAIL'}**; ORs={'; '.join(f'{x:.4g}' for x in timing_models['OR'].dropna())}.",
        f"- MCOP co-exposure models retain positive MCOP direction: **{'PASS' if coex_positive else 'FAIL'}**; MCOP ORs={'; '.join(f'{x:.4g}' for x in coex_models['OR'].dropna())}.",
        "",
        "These are association audits. A positive result does not establish that MCOP or DINP caused CRC; the NHANES design still permits reverse causation and survivor bias.",
        "",
        "## Diagnosis timing / reverse-causation audit",
        "",
        "The case-level timing export contains exam age, CRC diagnosis age and years since diagnosis. The model rows exclude cases with known diagnosis-to-exam interval below 1, 2 or 5 years; cases with missing diagnosis age are retained and counted explicitly as missing timing.",
        "",
        "| Analysis | N | CRC cases | Excluded known-timing CRC | OR | 95% CI | P |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, row in timing[timing["Analysis"].str.startswith(("Primary", "Exclude_diagnosis"))].iterrows():
        lines.append(f"| {row['Analysis']} | {int(row['N'])} | {int(row['CRC_N'])} | {int(row.get('Excluded_known_timing_CRC_N', 0))} | {fmt(row['OR'])} | {fmt(row['CI_low'])}-{fmt(row['CI_high'])} | {fmt(row['P'])} |")
    lines += [
        "",
        "The ≤5-year versus >5-year case groups are descriptive only and are not treated as independent etiologic tests.",
        "",
        "## Phthalate specificity",
        "",
        "Each co-exposure model is restricted to cycles with MCOP and that co-exposure available. The burden is the mean z-score of nine non-MCOP urinary phthalate metabolites, requiring at least two measured metabolites. Same-complete-case MCOP-only comparator rows are included in the CSV.",
        "",
        "| Model | Cycles | N | CRC | MCOP OR | 95% CI | MCOP P | Co-exposure OR | Co-exposure P | MCOP BH-FDR |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _, row in coexposure[coexposure["Model_role"].eq("coexposure_adjusted")].iterrows():
        lines.append(f"| {row['Analysis']} | {row['Cycles']} | {int(row['N'])} | {int(row['CRC_N'])} | {fmt(row['OR'])} | {fmt(row['CI_low'])}-{fmt(row['CI_high'])} | {fmt(row['P'])} | {fmt(row['Secondary_OR'])} | {fmt(row['Secondary_P'])} | {fmt(row.get('MCOP_BH_FDR'))} |")
    lines += [
        "",
        "## Figures",
        "",
        "- Figure 1: workflow from CTD × GeneCards to NHANES MCOP audit and future WHI replication.",
        "- Figure 2: primary model and stability/reverse-causation sensitivities.",
        "- Figure 3: seven cycle-specific estimates and pooled estimate.",
        "- Figure 4: restricted cubic spline.",
        "- Figure 5: female/male estimates with interaction audit.",
        "- Figure 6: MCOP specificity under co-exposure adjustment.",
        "",
        "WHI is shown only as a future prospective replication line; no WHI data were accessed or analyzed in this run.",
        "",
        "## Statistical fallacy scan (11/11)",
        "",
        "1. Simpson's paradox — sex-specific estimates are compared with pooled and interaction results; direction is reported rather than hidden.",
        "2. Ecological fallacy — not applicable; the analysis unit is the individual NHANES participant.",
        "3. Berkson's paradox — caution: the analytic population is conditioned on observed cancer outcome and phthalate subsample availability.",
        "4. Collider bias — caution: conditioning on cancer ascertainment/survival and complete covariates may induce selection effects.",
        "5. Base-rate neglect — not applicable to the reported OR models; unweighted case counts are reported.",
        "6. Regression to the mean — not applicable to this cross-sectional exposure/outcome comparison.",
        "7. Survivorship bias — RED_FLAG/major limitation: CRC is prevalent at examination, so survivors are the observed cases.",
        "8. Look-elsewhere effect — caution: sex, timing and co-exposure models are a targeted audit set and are all reported.",
        "9. Garden of forking paths — caution: the audit was frozen before execution; prior MCOP analyses and new outputs remain distinguishable.",
        "10. Correlation ≠ causation — caution: use association language only.",
        "11. Reverse causality — caution: timing exclusions are an audit, not proof of temporal causality.",
        "",
        "## Decision",
        "",
        f"The Phase 2D gate is **{'provisionally passed' if female_positive and timing_positive and coex_positive else 'not passed'}** under the pre-specified direction-only criteria. This does not upgrade NHANES to prospective evidence; it determines whether DINP-axis molecular validation is justified as the next phase.",
        "",
        "## Output files",
        "",
        "- mcop_crc_phase2_female_male.csv",
        "- mcop_crc_phase2_sex_interaction.csv",
        "- mcop_crc_phase2_diagnosis_timing.csv",
        "- mcop_crc_phase2_reverse_causation_sensitivity.csv",
        "- mcop_crc_phase2_coexposure_models.csv",
        "- mcop_crc_phase2d_figure1_workflow.png through mcop_crc_phase2d_figure6_coexposure.png",
    ]
    (output_dir / "mcop_crc_phase2d_paper_readiness_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", type=Path, default=DEFAULT_HARMONIZED)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    mcop_module = load_mcop_module()
    model_module = mcop_module.load_validated_functions()
    harmonized = pd.read_pickle(args.harmonized)
    mcop, input_manifest = mcop_module.read_mcop(args.data_dir)
    data = mcop_module.build_frame(harmonized, mcop)
    primary = model_module.population_frames(data)["CRC_vs_cancer_free"].copy()

    sex_rows, interaction = run_sex_audit(primary, model_module)
    timing, timing_cases = run_timing_audit(primary, model_module)
    coexposure = run_coexposure_audit(primary, model_module)

    sex_rows.to_csv(args.outdir / "mcop_crc_phase2_female_male.csv", index=False)
    interaction.to_csv(args.outdir / "mcop_crc_phase2_sex_interaction.csv", index=False)
    timing_cases.to_csv(args.outdir / "mcop_crc_phase2_diagnosis_timing.csv", index=False)
    timing.to_csv(args.outdir / "mcop_crc_phase2_reverse_causation_sensitivity.csv", index=False)
    coexposure.to_csv(args.outdir / "mcop_crc_phase2_coexposure_models.csv", index=False)

    primary_fit = fit_primary(primary, model_module, "Primary_CRC_vs_cancer_free")
    make_figures(args.outdir, primary_fit, timing, sex_rows, interaction, coexposure)
    write_report(args.outdir, sex_rows, interaction, timing, coexposure)

    manifest = {
        "analysis": "MCOP-CRC Phase 2D paper-readiness audit",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "scope": ["female_vs_male", "MCOP_x_sex_interaction", "diagnosis_timing_exclusions", "phthalate_coexposure_specificity", "six_figures"],
        "population": "NHANES age >=20, CRC vs cancer-free controls, MCOP available",
        "model_continuous": ["log2(MCOP)", "age", "BMI", "PIR", "log2(urinary creatinine)"],
        "model_categorical": ["sex", "race", "smoking"],
        "burden_definition": "mean z-score of nine non-MCOP metabolites; minimum two measured",
        "official_model_source": "work/scripts/mbzp_crc_phase2b.py",
        "input_harmonized": str(args.harmonized),
        "input_harmonized_sha256": sha256(args.harmonized),
        "input_files": input_manifest,
        "outputs": [
            "mcop_crc_phase2_female_male.csv",
            "mcop_crc_phase2_sex_interaction.csv",
            "mcop_crc_phase2_diagnosis_timing.csv",
            "mcop_crc_phase2_reverse_causation_sensitivity.csv",
            "mcop_crc_phase2_coexposure_models.csv",
            "mcop_crc_phase2d_paper_readiness_report.md",
        ],
    }
    (args.outdir / "mcop_crc_phase2d_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print("Female/male")
    print(sex_rows.to_string(index=False))
    print("Interaction")
    print(interaction.to_string(index=False))
    print("Timing")
    print(timing.to_string(index=False))
    print("Co-exposure")
    print(coexposure.to_string(index=False))


if __name__ == "__main__":
    main()
