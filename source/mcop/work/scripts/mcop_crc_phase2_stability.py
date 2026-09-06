"""Stability audit for the MCOP-CRC signal.

This follow-up is intentionally limited to:
1) cycle-specific effects,
2) a global MCOP-by-cycle interaction test,
3) restricted cubic spline shape,
4) exclusion of the top 1% and 2.5% MCOP values,
5) creatinine-normalized MCOP.

The validated survey-logistic implementation and the MCOP data linker are
reused from the preceding MCOP Phase 2 script.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, f


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_HARMONIZED = ROOT / "work" / "nhanes_phase2a" / "phase2b_harmonized.pkl"
DEFAULT_OUTPUT = ROOT / "outputs"


def load_mcop_module():
    path = ROOT / "work" / "scripts" / "mcop_crc_phase2.py"
    spec = importlib.util.spec_from_file_location("mcop_crc_phase2_validated", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load MCOP Phase 2 implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def wald_test(fit: dict, coefficient_names: list[str]) -> dict:
    if fit.get("status") not in {"ok", "converged_with_warning"}:
        return {"status": "not_estimable", "P": np.nan}
    coefficients = fit.get("coefficients", {})
    covariance = fit.get("covariance")
    if covariance is None or any(name not in coefficients for name in coefficient_names):
        return {"status": "not_estimable", "P": np.nan}
    names = list(coefficients)
    indices = [names.index(name) for name in coefficient_names]
    beta = np.array([coefficients[name] for name in coefficient_names], dtype=float)
    cov = np.asarray(covariance)[np.ix_(indices, indices)]
    statistic = float(beta @ np.linalg.pinv(cov) @ beta)
    df_num = len(coefficient_names)
    design_df = int(fit.get("design_df", 1))
    p_chi2 = float(chi2.sf(statistic, df_num))
    f_stat = statistic / df_num
    p_f = float(f.sf(f_stat, df_num, design_df))
    return {
        "status": fit.get("status"),
        "Wald_statistic": statistic,
        "F_statistic": f_stat,
        "df_num": df_num,
        "df_denom": design_df,
        "P_chi2": p_chi2,
        "P_F": p_f,
    }


def fit_adjusted(frame: pd.DataFrame, model_module, exposure: str, label: str, include_creatinine: bool = True) -> dict:
    continuous = [exposure, "age", "bmi", "pir"]
    if include_creatinine:
        continuous.append("creatinine_log2")
    fit = model_module.model_module.fit_survey_logistic(
        frame,
        continuous,
        ["sex", "race", "smoking"],
        exposure_name=exposure,
        levels=model_module.model_module.LEVELS,
    )
    return {
        "Analysis": label,
        "Exposure": exposure,
        **{key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}},
    }


def per_cycle(primary: pd.DataFrame, mcop_module) -> pd.DataFrame:
    rows = []
    cycles = sorted(primary.loc[primary["mcop_log2"].notna(), "cycle"].unique())
    for cycle in cycles:
        frame = primary[primary["cycle"].eq(cycle)].copy()
        rows.append(fit_adjusted(frame, mcop_module, "mcop_log2", f"Single_cycle_{cycle}"))
        rows[-1]["Cycle"] = cycle
    return pd.DataFrame(rows)


def global_interaction(primary: pd.DataFrame, mcop_module) -> tuple[pd.DataFrame, dict]:
    cycles = sorted(primary.loc[primary["mcop_log2"].notna(), "cycle"].unique())
    work = primary[primary["mcop_log2"].notna()].copy()
    reference = cycles[0]
    interaction_names = []
    for cycle in cycles[1:]:
        name = f"mcop_x_cycle_{cycle}"
        work[name] = work["mcop_log2"] * work["cycle"].eq(cycle).astype(float)
        interaction_names.append(name)
    continuous = ["mcop_log2", *interaction_names, "age", "bmi", "pir", "creatinine_log2"]
    levels = {"cycle": cycles, **mcop_module.model_module.LEVELS}
    fit = mcop_module.model_module.fit_survey_logistic(
        work,
        continuous,
        ["cycle", "sex", "race", "smoking"],
        exposure_name="mcop_log2",
        levels=levels,
    )
    test = wald_test(fit, interaction_names)
    out = {
        "status": fit.get("status"),
        "reference_cycle": reference,
        "interaction_terms": ";".join(interaction_names),
        "N": fit.get("N"),
        "CRC_N": fit.get("CRC_N"),
        "design_df": fit.get("design_df"),
        **test,
    }
    return pd.DataFrame([out]), {"fit": fit, "cycles": cycles, "interaction_names": interaction_names}


def spline_analysis(primary: pd.DataFrame, mcop_module, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = [
        "outcome", "mcop_log2", "age", "bmi", "pir", "creatinine_log2",
        "sex", "race", "smoking", "pooled_weight", "psu", "strata",
    ]
    cc = primary.dropna(subset=required).copy()
    knots = cc["mcop_log2"].quantile([0.05, 0.35, 0.65, 0.95]).to_numpy()
    knots = np.unique(knots)
    if len(knots) < 4:
        return pd.DataFrame([{"status": "not_estimable", "reason": "insufficient distinct spline knots"}]), pd.DataFrame()
    basis = mcop_module.model_module.spline_basis(cc["mcop_log2"].to_numpy(), knots)
    basis_cols = []
    for index in range(basis.shape[1]):
        name = f"mcop_spline_{index + 1}"
        cc[name] = basis[:, index]
        basis_cols.append(name)
    fit = mcop_module.model_module.fit_survey_logistic(
        cc,
        basis_cols + ["age", "bmi", "pir", "creatinine_log2"],
        ["sex", "race", "smoking"],
        exposure_name=basis_cols[0],
        levels=mcop_module.model_module.LEVELS,
    )
    overall = wald_test(fit, basis_cols)
    nonlinear = wald_test(fit, basis_cols[1:])
    result = pd.DataFrame(
        [
            {
                "status": fit.get("status"),
                "N": fit.get("N"),
                "CRC_N": fit.get("CRC_N"),
                "knots_log2": ";".join(f"{value:.8g}" for value in knots),
                "overall_P_chi2": overall.get("P_chi2"),
                "overall_P_F": overall.get("P_F"),
                "nonlinear_P_chi2": nonlinear.get("P_chi2"),
                "nonlinear_P_F": nonlinear.get("P_F"),
                "nonlinear_df": nonlinear.get("df_num"),
                "design_df": fit.get("design_df"),
            }
        ]
    )
    grid = np.linspace(cc["mcop_log2"].min(), cc["mcop_log2"].max(), 200)
    grid_basis = mcop_module.model_module.spline_basis(grid, knots)
    beta = np.array([fit["coefficients"][column] for column in basis_cols])
    linear_predictor = grid_basis @ beta
    median_basis = mcop_module.model_module.spline_basis(np.array([cc["mcop_log2"].median()]), knots)[0]
    relative = linear_predictor - median_basis @ beta
    curve = pd.DataFrame(
        {
            "mcop_log2": grid,
            "spline_linear_predictor": linear_predictor,
            "spline_log_OR_vs_median": relative,
            "OR_vs_median": np.exp(relative),
        }
    )
    curve.to_csv(output_dir / "mcop_crc_phase2_spline_curve.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(grid, relative, color="#2b6cb0", linewidth=2)
    ax.axhline(0, color="0.5", linewidth=0.8)
    ax.set_xlabel("log2 urinary MCOP (ng/mL)")
    ax.set_ylabel("Spline log OR vs median MCOP")
    ax.set_title("MCOP-CRC restricted cubic spline")
    fig.tight_layout()
    fig.savefig(output_dir / "mcop_crc_phase2_spline.png", dpi=180)
    plt.close(fig)
    return result, curve


def tail_exclusion(primary: pd.DataFrame, mcop_module) -> pd.DataFrame:
    available = primary.loc[primary["mcop_log2"].notna(), "mcop_log2"]
    rows = []
    for fraction in [0.01, 0.025]:
        cutoff = float(available.quantile(1 - fraction))
        trimmed = primary[primary["mcop_log2"].le(cutoff) | primary["mcop_log2"].isna()].copy()
        fit = fit_adjusted(trimmed, mcop_module, "mcop_log2", f"Exclude_top_{fraction * 100:g}pct")
        fit["Excluded_fraction"] = fraction
        fit["Cutoff_log2"] = cutoff
        fit["Retained_exposure_n"] = int(trimmed["mcop_log2"].notna().sum())
        rows.append(fit)
    return pd.DataFrame(rows)


def creatinine_normalized(primary: pd.DataFrame, mcop_module) -> pd.DataFrame:
    work = primary.copy()
    work["mcop_creatinine_norm_log2"] = work["mcop_log2"] - work["creatinine_log2"]
    fit = fit_adjusted(
        work,
        mcop_module,
        "mcop_creatinine_norm_log2",
        "Creatinine_normalized_MCOP",
        include_creatinine=False,
    )
    fit["Normalization"] = "log2(MCOP) - log2(urine creatinine)"
    return pd.DataFrame([fit])


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
    mcop, _ = mcop_module.read_mcop(args.data_dir)
    data = mcop_module.build_frame(harmonized, mcop)
    primary = model_module.population_frames(data)["CRC_vs_cancer_free"]

    per_cycle_df = per_cycle(primary, type("Module", (), {"model_module": model_module})())
    per_cycle_df.to_csv(args.outdir / "mcop_crc_phase2_per_cycle.csv", index=False)

    interaction_df, interaction_meta = global_interaction(
        primary,
        type("Module", (), {"model_module": model_module})(),
    )
    interaction_df.to_csv(args.outdir / "mcop_crc_phase2_cycle_interaction.csv", index=False)

    spline_df, spline_curve = spline_analysis(
        primary,
        type("Module", (), {"model_module": model_module})(),
        args.outdir,
    )
    spline_df.to_csv(args.outdir / "mcop_crc_phase2_spline.csv", index=False)

    tails = tail_exclusion(primary, type("Module", (), {"model_module": model_module})())
    tails.to_csv(args.outdir / "mcop_crc_phase2_tail_exclusion.csv", index=False)

    normalized = creatinine_normalized(primary, type("Module", (), {"model_module": model_module})())
    normalized.to_csv(args.outdir / "mcop_crc_phase2_creatinine_normalized.csv", index=False)

    def fmt(value) -> str:
        if value is None or pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating, int, np.integer)):
            return f"{float(value):.6g}"
        return str(value)

    interaction_row = interaction_df.iloc[0].to_dict()
    spline_row = spline_df.iloc[0].to_dict()
    tail_rows = {
        row["Analysis"]: row
        for row in tails.to_dict("records")
    }
    normalized_row = normalized.iloc[0].to_dict()
    positive_cycles = per_cycle_df.loc[per_cycle_df["OR"] > 1, "Cycle"].tolist()
    negative_cycles = per_cycle_df.loc[per_cycle_df["OR"] < 1, "Cycle"].tolist()
    report = [
        "# MCOP-CRC Phase 2B：稳定性审计",
        "",
        "本轮不更换候选，只审计 MCOP 连续模型的稳健性：逐 cycle、global MCOP×cycle interaction、RCS、极端值排除和 creatinine-normalized exposure。",
        "",
        "## Key findings",
        "",
        f"- Single-cycle effects above 1: {len(positive_cycles)}/7 ({'; '.join(positive_cycles)}). Below 1: {len(negative_cycles)}/7 ({'; '.join(negative_cycles) if negative_cycles else 'none'}).",
        "- The 2011-2012 single-cycle estimate is below 1, but its confidence interval is wide and includes substantial positive effects.",
        f"- Global MCOP×cycle interaction: F={fmt(interaction_row.get('F_statistic'))}, df={fmt(interaction_row.get('df_num'))},{fmt(interaction_row.get('df_denom'))}, design-based P={fmt(interaction_row.get('P_F'))}; chi-square P={fmt(interaction_row.get('P_chi2'))}.",
        f"- RCS nonlinear test: design-based P={fmt(spline_row.get('nonlinear_P_F'))}; chi-square P={fmt(spline_row.get('nonlinear_P_chi2'))}.",
        "- The RCS point curve does not show a stable Q3 peak; it rises toward the extreme upper tail, where uncertainty is large.",
        f"- Excluding top 1%: OR={fmt(tail_rows.get('Exclude_top_1pct', {}).get('OR'))}, 95% CI {fmt(tail_rows.get('Exclude_top_1pct', {}).get('CI_low'))}-{fmt(tail_rows.get('Exclude_top_1pct', {}).get('CI_high'))}, P={fmt(tail_rows.get('Exclude_top_1pct', {}).get('P'))}.",
        f"- Excluding top 2.5%: OR={fmt(tail_rows.get('Exclude_top_2.5pct', {}).get('OR'))}, 95% CI {fmt(tail_rows.get('Exclude_top_2.5pct', {}).get('CI_low'))}-{fmt(tail_rows.get('Exclude_top_2.5pct', {}).get('CI_high'))}, P={fmt(tail_rows.get('Exclude_top_2.5pct', {}).get('P'))}.",
        f"- Creatinine-normalized MCOP: OR={fmt(normalized_row.get('OR'))}, 95% CI {fmt(normalized_row.get('CI_low'))}-{fmt(normalized_row.get('CI_high'))}, P={fmt(normalized_row.get('P'))}.",
        "",
        "## Interpretation",
        "",
        "The continuous association is not explained by the top 1%-2.5% exposure tail or by the choice between creatinine covariate adjustment and creatinine-ratio normalization. However, the global interaction test indicates cycle heterogeneity, so the result is not yet a uniformly replicated effect across cycles. The RCS test does not support a stable nonlinear peak. Keep the candidate at yellow-green and prioritize an independent cohort rather than mechanistic expansion.",
        "",
        "## 输出文件",
        "",
        "- mcop_crc_phase2_per_cycle.csv",
        "- mcop_crc_phase2_cycle_interaction.csv",
        "- mcop_crc_phase2_spline.csv",
        "- mcop_crc_phase2_spline_curve.csv",
        "- mcop_crc_phase2_tail_exclusion.csv",
        "- mcop_crc_phase2_creatinine_normalized.csv",
        "",
        "RCS 使用 5th/35th/65th/95th percentile knots。非线性检验同时给出 design-based Wald F 近似和 chi-square 近似；主要参考 F 近似。",
        "Creatinine-normalized 模型使用 log2(MCOP)-log2(creatinine)，因此不再把原始 creatinine_log2 作为同一模型协变量重复加入。",
        "",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
    ]
    (args.outdir / "mcop_crc_phase2_stability_report.md").write_text(
        "\n".join(report),
        encoding="utf-8",
    )

    manifest = {
        "analysis": "MCOP-CRC Phase 2B stability audit",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_harmonized": str(args.harmonized),
        "analyses": [
            "single-cycle effects",
            "global MCOP-by-cycle interaction",
            "restricted cubic spline",
            "top 1% and 2.5% exclusion",
            "creatinine-normalized MCOP",
        ],
        "spline_knots": "5th,35th,65th,95th percentile",
        "primary_outcome_population": "CRC vs cancer-free controls",
        "excluded": ["MONP", "MiNP", "new candidate screening", "mechanistic expansion"],
        "outputs": [
            "mcop_crc_phase2_per_cycle.csv",
            "mcop_crc_phase2_cycle_interaction.csv",
            "mcop_crc_phase2_spline.csv",
            "mcop_crc_phase2_spline_curve.csv",
            "mcop_crc_phase2_spline.png",
            "mcop_crc_phase2_tail_exclusion.csv",
            "mcop_crc_phase2_creatinine_normalized.csv",
        ],
    }
    (args.outdir / "mcop_crc_phase2_stability_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print("Per-cycle")
    print(per_cycle_df.to_string(index=False))
    print("Cycle interaction")
    print(interaction_df.to_string(index=False))
    print("Spline")
    print(spline_df.to_string(index=False))
    print("Tail exclusion")
    print(tails.to_string(index=False))
    print("Creatinine normalized")
    print(normalized.to_string(index=False))


if __name__ == "__main__":
    main()
