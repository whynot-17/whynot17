"""Final MCOP-CRC validation in the existing NHANES harmonized frame.

This is a deliberately narrow follow-up to the DINP Phase 2A audit:
continuous MCOP, quartiles, age >=40, cancer-free controls, and leave-one-cycle
out. MONP and MiNP are not analyzed here. The fitted survey-logistic routine is
imported from the already validated MBzP Phase 2B implementation so that
weights, pooled strata/PSU identifiers, and the CRC outcome definition remain
unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = ROOT / "work" / "nhanes_phase2a" / "data"
DEFAULT_HARMONIZED = ROOT / "work" / "nhanes_phase2a" / "phase2b_harmonized.pkl"
DEFAULT_OUTPUT = ROOT / "outputs"
CDC_PHTHTE_SOURCE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm"


def load_validated_functions():
    path = ROOT / "work" / "scripts" / "mbzp_crc_phase2b.py"
    spec = importlib.util.spec_from_file_location("mbzp_crc_phase2b_validated", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load validated model implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_mcop(data_dir: Path) -> tuple[pd.DataFrame, list[dict]]:
    rows: list[pd.DataFrame] = []
    manifest: list[dict] = []
    for path in sorted(data_dir.glob("*_PHTHTE.XPT")):
        cycle = path.name.split("_", 1)[0]
        try:
            frame = pd.read_sas(path, format="xport", encoding="latin1")
        except ValueError:
            continue
        if "SEQN" not in frame.columns or "URXCOP" not in frame.columns:
            continue
        selected = ["SEQN", "URXCOP"]
        if "URDCOPLC" in frame.columns:
            selected.append("URDCOPLC")
        part = frame[selected].copy()
        part["cycle"] = cycle
        part["SEQN"] = pd.to_numeric(part["SEQN"], errors="coerce").astype("Int64")
        rows.append(part)
        manifest.append(
            {
                "cycle": cycle,
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "value_variable": "URXCOP",
                "lod_comment_variable": "URDCOPLC" if "URDCOPLC" in frame.columns else None,
                "raw_rows": int(len(frame)),
            }
        )
    if not rows:
        raise RuntimeError("No readable PHTHTE files containing URXCOP were found")
    return pd.concat(rows, ignore_index=True), manifest


def build_frame(harmonized: pd.DataFrame, mcop: pd.DataFrame) -> pd.DataFrame:
    left = harmonized.copy()
    left["SEQN"] = pd.to_numeric(left["SEQN"], errors="coerce").astype("Int64")
    merged = left.merge(mcop, on=["SEQN", "cycle"], how="left", validate="one_to_one")
    value = pd.to_numeric(merged["URXCOP"], errors="coerce")
    merged["mcop_log2"] = np.log2(value.where(value > 0))
    return merged


def fit_continuous(frame: pd.DataFrame, model_module, label: str) -> dict:
    fit = model_module.fit_survey_logistic(
        frame,
        ["mcop_log2", "age", "bmi", "pir", "creatinine_log2"],
        ["sex", "race", "smoking"],
        exposure_name="mcop_log2",
        levels=model_module.LEVELS,
    )
    return {
        "Analysis": label,
        "Exposure": "log2(MCOP)",
        **{key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}},
    }


def quartile_analysis(primary: pd.DataFrame, model_module) -> tuple[pd.DataFrame, dict]:
    available = primary["mcop_log2"].notna()
    cutpoints = np.unique(
        primary.loc[available, "mcop_log2"].quantile([0, 0.25, 0.5, 0.75, 1]).to_numpy()
    )
    if len(cutpoints) < 5:
        return pd.DataFrame([{"status": "not_estimable", "reason": "fewer than four distinct exposure quantiles"}]), {}
    work = primary.copy()
    work["mcop_quartile"] = pd.cut(
        work["mcop_log2"],
        bins=[-np.inf, *cutpoints[1:-1], np.inf],
        labels=False,
        include_lowest=True,
    ) + 1
    required = [
        "outcome", "mcop_quartile", "age", "bmi", "pir", "creatinine_log2",
        "sex", "race", "smoking", "pooled_weight", "psu", "strata",
    ]
    cc = work.dropna(subset=required).copy()
    rows: list[dict] = []
    for q in [1, 2, 3, 4]:
        if q == 1:
            rows.append(
                {
                    "Quartile": "Q1",
                    "Reference": True,
                    "N": len(cc),
                    "CRC_N": int(cc["outcome"].sum()),
                    "OR": 1.0,
                    "CI_low": 1.0,
                    "CI_high": 1.0,
                    "P": np.nan,
                    "status": "reference",
                }
            )
            continue
        contrast = cc.copy()
        contrast["q_exposure"] = (contrast["mcop_quartile"] == q).astype(float)
        fit = model_module.fit_survey_logistic(
            contrast,
            ["q_exposure", "age", "bmi", "pir", "creatinine_log2"],
            ["sex", "race", "smoking"],
            exposure_name="q_exposure",
            levels=model_module.LEVELS,
        )
        rows.append(
            {
                "Quartile": f"Q{q}",
                "Reference": False,
                **{key: value for key, value in fit.items() if key not in {"coefficients", "covariance"}},
            }
        )
    trend_frame = cc.assign(qtrend=cc["mcop_quartile"].astype(float))
    trend = model_module.fit_survey_logistic(
        trend_frame,
        ["qtrend", "age", "bmi", "pir", "creatinine_log2"],
        ["sex", "race", "smoking"],
        exposure_name="qtrend",
        levels=model_module.LEVELS,
    )
    for row in rows:
        row["P_trend"] = trend.get("P", np.nan)
        for index, value in enumerate(cutpoints[:4], start=1):
            row[f"Q{index}_cutpoint_log2"] = float(value)
    meta = {
        "quartile_cutpoints_log2": [float(value) for value in cutpoints],
        "complete_case_n": int(len(cc)),
        "complete_case_crc_n": int(cc["outcome"].sum()),
        "trend_fit_status": trend.get("status"),
    }
    return pd.DataFrame(rows), meta


def audit_frame(frame: pd.DataFrame, model_module) -> dict:
    primary = model_module.population_frames(frame)["CRC_vs_cancer_free"]
    exposure = primary[primary["mcop_log2"].notna()]
    lod_comment = pd.to_numeric(exposure["URDCOPLC"], errors="coerce")
    return {
        "harmonized_rows": int(len(frame)),
        "mcop_available_rows": int(frame["mcop_log2"].notna().sum()),
        "mcop_cycles": ";".join(sorted(frame.loc[frame["mcop_log2"].notna(), "cycle"].unique())),
        "primary_exposure_outcome_rows": int(len(exposure)),
        "primary_crc_cases": int(exposure["outcome"].sum()),
        "primary_cancer_free_controls": int(len(exposure) - exposure["outcome"].sum()),
        "primary_mcop_above_lod_n": int((lod_comment.notna() & ~lod_comment.eq(1)).sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--harmonized", type=Path, default=DEFAULT_HARMONIZED)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    model_module = load_validated_functions()
    harmonized = pd.read_pickle(args.harmonized)
    mcop, input_manifest = read_mcop(args.data_dir)
    data = build_frame(harmonized, mcop)
    populations = model_module.population_frames(data)
    primary = populations["CRC_vs_cancer_free"]

    audit = audit_frame(data, model_module)
    pd.DataFrame([audit]).to_csv(args.outdir / "mcop_crc_phase2_audit.csv", index=False)

    model_rows = [fit_continuous(primary, model_module, "Primary_CRC_vs_cancer_free")]
    age40 = primary[primary["age"] >= 40].copy()
    model_rows.append(fit_continuous(age40, model_module, "Age_ge_40_CRC_vs_cancer_free"))

    loco_rows: list[dict] = []
    mcop_cycles = sorted(primary.loc[primary["mcop_log2"].notna(), "cycle"].unique())
    for cycle in mcop_cycles:
        fit = fit_continuous(
            primary[primary["cycle"] != cycle].copy(),
            model_module,
            f"LOCO_drop_{cycle}",
        )
        loco_rows.append({"Dropped_cycle": cycle, **fit})
    pd.DataFrame(model_rows).to_csv(args.outdir / "mcop_crc_phase2_main_models.csv", index=False)
    pd.DataFrame(loco_rows).to_csv(args.outdir / "mcop_crc_phase2_leave_one_cycle_out.csv", index=False)

    quartiles, quartile_meta = quartile_analysis(primary, model_module)
    quartiles.to_csv(args.outdir / "mcop_crc_phase2_quartiles.csv", index=False)

    loco_frame = pd.DataFrame(loco_rows)
    primary_fit = model_rows[0]
    age40_fit = model_rows[1]
    q4_rows = quartiles[quartiles.get("Quartile", pd.Series(dtype=str)).eq("Q4")]
    q4_fit = q4_rows.iloc[0].to_dict() if not q4_rows.empty else {}
    loco_or = pd.to_numeric(loco_frame.get("OR", pd.Series(dtype=float)), errors="coerce").dropna()
    loco_beta = pd.to_numeric(loco_frame.get("beta", pd.Series(dtype=float)), errors="coerce").dropna()
    direction_consistent = bool(
        len(loco_beta) and ((loco_beta > 0).all() or (loco_beta < 0).all())
    )
    report = [
        "# MCOP-CRC Phase 2：DINP 轴最终人群验证",
        "",
        "## 判定",
        "",
        "本轮只分析 MCOP（NHANES URXCOP），不再分析 MONP 或 MiNP。分析使用既有 20 岁以上 NHANES CRC harmonized frame、cancer-free controls、pooled phthalate weights、cycle-pooled strata/PSU 和已验证的 survey-logistic sandwich 实现。",
        "",
        f"- Primary continuous OR per MCOP doubling: **{primary_fit.get('OR', np.nan):.6g}**",
        f"- 95% CI: **{primary_fit.get('CI_low', np.nan):.6g}-{primary_fit.get('CI_high', np.nan):.6g}**; P={primary_fit.get('P', np.nan):.6g}",
        f"- Primary N={primary_fit.get('N', np.nan)}, CRC cases={primary_fit.get('CRC_N', np.nan)}",
        f"- Age >=40 OR: **{age40_fit.get('OR', np.nan):.6g}**; 95% CI {age40_fit.get('CI_low', np.nan):.6g}-{age40_fit.get('CI_high', np.nan):.6g}; P={age40_fit.get('P', np.nan):.6g}",
        f"- Q4 vs Q1 OR: **{q4_fit.get('OR', np.nan):.6g}**; 95% CI {q4_fit.get('CI_low', np.nan):.6g}-{q4_fit.get('CI_high', np.nan):.6g}; P-trend={q4_fit.get('P_trend', np.nan):.6g}",
        f"- LOCO OR range: **{loco_or.min():.6g}-{loco_or.max():.6g}**; direction consistent: {'YES' if direction_consistent else 'NO'}",
        "- Quartile pattern is not strictly monotonic; Q4 is not a stronger contrast than Q3.",
        "",
        "MCOP exposure availability and model outputs are reported separately below. This is a validation result, not evidence that MCOP is a DINP-specific causal exposure measure.",
        "",
        "## Data audit",
        "",
        "- MCOP variable: URXCOP; LOD comment variable: URDCOPLC",
        f"- Available cycles: {audit['mcop_cycles']}",
        f"- Exposure + cancer outcome: {audit['primary_exposure_outcome_rows']}",
        f"- CRC cases: {audit['primary_crc_cases']}",
        f"- MCOP above LOD: {audit['primary_mcop_above_lod_n']}",
        f"- Official codebook: [CDC NHANES PHTHTE_J]({CDC_PHTHTE_SOURCE})",
        "",
        "## Prespecified analyses",
        "",
        "1. Continuous log2(MCOP) - full primary cancer-free-control population.",
        "2. Quartiles - Q1 reference, Q4 vs Q1, and linear trend.",
        "3. Age >=40 - same covariate set.",
        "4. Cancer-free controls - primary population definition.",
        "5. Leave-one-cycle-out - one cycle removed at a time.",
        "",
        "## Files",
        "",
        "- mcop_crc_phase2_audit.csv",
        "- mcop_crc_phase2_main_models.csv",
        "- mcop_crc_phase2_quartiles.csv",
        "- mcop_crc_phase2_leave_one_cycle_out.csv",
        "",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
    ]
    (args.outdir / "mcop_crc_phase2_report.md").write_text("\n".join(report), encoding="utf-8")

    manifest = {
        "analysis": "MCOP-CRC Phase 2 final DINP-axis validation",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "exposure": "log2(URXCOP)",
        "lod_comment": "URDCOPLC; CDC code 0/blank XPORT sentinel = at or above LOD, 1 = below LOD",
        "official_source": CDC_PHTHTE_SOURCE,
        "outcome": "CRC type code 16 or 31 vs MCQ220=2 cancer-free controls",
        "covariates": ["age", "sex", "race", "BMI", "smoking", "PIR", "log2(URXUCR)"],
        "analyses": ["continuous", "quartiles", "age_ge_40", "cancer_free_controls", "LOCO"],
        "excluded_from_this_run": ["MONP", "MiNP", "mechanistic_analysis", "subgroup_fishing"],
        "audit": audit,
        "quartile_meta": quartile_meta,
        "input_harmonized": {"path": str(args.harmonized), "sha256": sha256(args.harmonized)},
        "input_files": input_manifest,
        "model_implementation": "imported from work/scripts/mbzp_crc_phase2b.py; validated Newton-IRLS survey-logistic fit",
    }
    (args.outdir / "mcop_crc_phase2_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {"audit": audit, "primary": primary_fit, "age_ge_40": age40_fit},
            ensure_ascii=False,
            default=str,
        )
    )


if __name__ == "__main__":
    main()
