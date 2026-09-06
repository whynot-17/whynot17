"""Cycle heterogeneity audit for MCOP-CRC NHANES analyses.

The audit describes the seven MCOP-observed cycles and does not fit new
etiologic models. It quantifies exposure distributions, LOD coverage,
creatinine, CRC age structure, demographic composition, weighted CRC prevalence,
assay/LLOD metadata, and pooled quartile event counts.
"""

from __future__ import annotations

import argparse
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


ASSAY_METADATA = {
    "2005-2006": {
        "codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PHTHTE_D.htm",
        "llod_ng_mL": 0.7,
        "platform": "HPLC-ESI-MS/MS",
        "laboratory": "CDC Division of Laboratory Sciences / NCEH",
        "method_note": "Online SPE, reversed-phase HPLC-ESI-MS/MS, isotopically labeled internal standards.",
        "change_note": "No cycle-specific no-change statement located in the codebook; MCOP included in the listed phthalate panel.",
    },
    "2007-2008": {
        "codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/PHTHTE_E.htm",
        "llod_ng_mL": 0.7,
        "platform": "HPLC-ESI-MS/MS",
        "laboratory": "CDC Division of Laboratory Sciences / NCEH",
        "method_note": "Online SPE, reversed-phase HPLC-ESI-MS/MS, isotopically labeled internal standards.",
        "change_note": "No cycle-specific no-change statement located in the codebook; MCOP included in the listed phthalate panel.",
    },
    "2009-2010": {
        "codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/PHTHTE_F.htm",
        "llod_ng_mL": 0.2,
        "platform": "HPLC-ESI-MS/MS",
        "laboratory": "CDC Division of Laboratory Sciences / NCEH",
        "method_note": "Online SPE, reversed-phase HPLC-ESI-MS/MS, isotopically labeled internal standards.",
        "change_note": "No cycle-specific no-change statement located in the codebook; MCOP included in the listed phthalate panel.",
    },
    "2011-2012": {
        "codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PHTHTE_G.htm",
        "method_manual": "https://wwwn.cdc.gov/nchs/data/nhanes/public/2011/labmethods/phthte_g_met.pdf",
        "llod_ng_mL": 0.2,
        "platform": "HPLC/ESI-MS/MS (Method 6306.04)",
        "laboratory": "Personal Care Products Laboratory, Organic Analytical Toxicology Branch, CDC NCEH",
        "method_note": "Online SPE, reversed-phase HPLC-ESI-MS/MS, isotopically labeled standards; manual explicitly maps MCOP to the di-isononyl phthalate axis.",
        "change_note": "Method manual revised April 2013; no explicit cycle-level no-change statement in the codebook.",
    },
    "2013-2014": {
        "codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/PHTHTE_H.htm",
        "llod_ng_mL": 0.3,
        "platform": "HPLC-ESI-MS/MS",
        "laboratory": "CDC Division of Laboratory Sciences / NCEH",
        "method_note": "Online SPE, reversed-phase HPLC-ESI-MS/MS, isotopically labeled internal standards.",
        "change_note": "CDC codebook states no changes to lab method, equipment, or site in 2013-2014.",
    },
    "2015-2016": {
        "codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PHTHTE_I.htm",
        "llod_ng_mL": 0.3,
        "platform": "HPLC-ESI-MS/MS",
        "laboratory": "CDC Division of Laboratory Sciences / NCEH",
        "method_note": "Online SPE, reversed-phase HPLC-ESI-MS/MS, isotopically labeled internal standards.",
        "change_note": "CDC codebook states no changes to lab method, equipment, or site in 2015-2016.",
    },
    "2017-2018": {
        "codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm",
        "llod_ng_mL": 0.3,
        "platform": "HPLC-ESI-MS/MS",
        "laboratory": "CDC Division of Laboratory Sciences / NCEH",
        "method_note": "Online SPE, reversed-phase HPLC-ESI-MS/MS, isotopically labeled internal standards; MONP was additionally measured with a separate isotope-dilution method.",
        "change_note": "CDC codebook states no changes to lab method, equipment, or site in 2017-2018.",
    },
}


def load_mcop_module():
    path = ROOT / "work" / "scripts" / "mcop_crc_phase2.py"
    spec = importlib.util.spec_from_file_location("mcop_crc_phase2_validated", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load MCOP Phase 2 implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def weighted_quantile(values: pd.Series, weights: pd.Series, quantiles: list[float]) -> list[float]:
    frame = pd.DataFrame(
        {
            "value": pd.to_numeric(values, errors="coerce"),
            "weight": pd.to_numeric(weights, errors="coerce"),
        }
    ).dropna()
    frame = frame[frame["weight"] > 0].sort_values("value")
    if frame.empty:
        return [np.nan for _ in quantiles]
    cumulative = frame["weight"].cumsum().to_numpy()
    cumulative = cumulative / cumulative[-1]
    return [
        float(np.interp(q, cumulative, frame["value"].to_numpy()))
        for q in quantiles
    ]


def pct(series: pd.Series, predicate) -> float:
    values = series.notna()
    denominator = int(values.sum())
    return float(100.0 * (predicate(series) & values).sum() / denominator) if denominator else np.nan


def iqr_values(series: pd.Series) -> tuple[float, float, float]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return np.nan, np.nan, np.nan
    q = values.quantile([0.25, 0.5, 0.75])
    return float(q.loc[0.25]), float(q.loc[0.5]), float(q.loc[0.75])


def cycle_summary(primary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cycles = sorted(primary.loc[primary["mcop_log2"].notna(), "cycle"].unique())
    model_required = [
        "outcome", "mcop_log2", "age", "bmi", "pir", "creatinine_log2",
        "sex", "race", "smoking", "pooled_weight", "psu", "strata",
    ]
    for cycle in cycles:
        frame = primary[primary["cycle"].eq(cycle) & primary["mcop_log2"].notna()].copy()
        mcop = pd.to_numeric(frame["URXCOP"], errors="coerce")
        creat = pd.to_numeric(frame["URXUCR"], errors="coerce")
        q25, median, q75 = iqr_values(mcop)
        weighted_q25, weighted_median, weighted_q75, weighted_p95, weighted_p99 = weighted_quantile(
            mcop,
            frame["pooled_weight"],
            [0.25, 0.5, 0.75, 0.95, 0.99],
        )
        crc = frame[frame["outcome"].eq(1)]
        controls = frame[frame["outcome"].eq(0)]
        model_cc = frame.dropna(subset=model_required).copy()
        model_crc = model_cc[model_cc["outcome"].eq(1)]
        model_controls = model_cc[model_cc["outcome"].eq(0)]
        model_crc_median = float(model_crc["URXCOP"].median()) if not model_crc.empty else np.nan
        model_control_median = float(model_controls["URXCOP"].median()) if not model_controls.empty else np.nan
        crc_age_q25, crc_age_median, crc_age_q75 = iqr_values(crc["age"])
        lod_comment = pd.to_numeric(frame["URDCOPLC"], errors="coerce")
        race = frame["race"].astype("string")
        weight = pd.to_numeric(frame["pooled_weight"], errors="coerce")
        valid_weight = weight.notna() & weight.gt(0)
        weighted_crc_prevalence = (
            float(
                (weight[valid_weight] * frame.loc[valid_weight, "outcome"]).sum()
                / weight[valid_weight].sum()
            )
            if valid_weight.any()
            else np.nan
        )
        rows.append(
            {
                "cycle": cycle,
                "N_primary_exposure_outcome": int(len(frame)),
                "CRC_cases": int(frame["outcome"].sum()),
                "CRC_case_pct_unweighted": float(100.0 * frame["outcome"].mean()),
                "CRC_prevalence_weighted": weighted_crc_prevalence,
                "Model_complete_case_N": int(len(model_cc)),
                "Model_complete_case_CRC_cases": int(model_cc["outcome"].sum()),
                "Model_CRC_MCOP_median_ng_mL": model_crc_median,
                "Model_control_MCOP_median_ng_mL": model_control_median,
                "Model_case_control_MCOP_median_ratio": (
                    model_crc_median / model_control_median
                    if pd.notna(model_crc_median) and pd.notna(model_control_median) and model_control_median > 0
                    else np.nan
                ),
                "MCOP_min_ng_mL": float(mcop.min()),
                "MCOP_Q1_ng_mL": q25,
                "MCOP_median_ng_mL": median,
                "MCOP_Q3_ng_mL": q75,
                "MCOP_P95_ng_mL": float(mcop.quantile(0.95)),
                "MCOP_P99_ng_mL": float(mcop.quantile(0.99)),
                "MCOP_weighted_Q1_ng_mL": weighted_q25,
                "MCOP_weighted_median_ng_mL": weighted_median,
                "MCOP_weighted_Q3_ng_mL": weighted_q75,
                "MCOP_weighted_P95_ng_mL": weighted_p95,
                "MCOP_weighted_P99_ng_mL": weighted_p99,
                "MCOP_above_LOD_pct": float(100.0 * (lod_comment.notna() & ~lod_comment.eq(1)).mean()),
                "MCOP_below_LOD_n": int(lod_comment.eq(1).sum()),
                "LLOD_ng_mL_codebook": ASSAY_METADATA[cycle]["llod_ng_mL"],
                "Creatinine_median_mg_dL": float(creat.median()),
                "Creatinine_Q1_mg_dL": float(creat.quantile(0.25)),
                "Creatinine_Q3_mg_dL": float(creat.quantile(0.75)),
                "Creatinine_P95_mg_dL": float(creat.quantile(0.95)),
                "Age_median_all": float(frame["age"].median()),
                "Age_median_CRC": crc_age_median,
                "CRC_age_Q1": crc_age_q25,
                "CRC_age_Q3": crc_age_q75,
                "Female_pct": pct(frame["sex"], lambda x: x.eq("Female")),
                "NonHispanicWhite_pct": pct(race, lambda x: x.eq("Non-Hispanic White")),
                "NonHispanicBlack_pct": pct(race, lambda x: x.eq("Non-Hispanic Black")),
                "Hispanic_pct": pct(race, lambda x: x.isin(["Mexican American", "Other Hispanic"])),
                "OtherRace_pct": pct(race, lambda x: x.eq("Other/Multi")),
                "Smoking_current_pct": pct(frame["smoking"], lambda x: x.eq("Current")),
                "Smoking_former_pct": pct(frame["smoking"], lambda x: x.eq("Former")),
                "Smoking_never_pct": pct(frame["smoking"], lambda x: x.eq("Never")),
                "assay_platform": ASSAY_METADATA[cycle]["platform"],
            }
        )
    return pd.DataFrame(rows)


def quartile_case_counts(primary: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    available = primary["mcop_log2"].notna()
    cutpoints = np.unique(
        primary.loc[available, "mcop_log2"].quantile([0, 0.25, 0.5, 0.75, 1]).to_numpy()
    )
    required = [
        "outcome", "mcop_log2", "age", "bmi", "pir", "creatinine_log2",
        "sex", "race", "smoking", "pooled_weight", "psu", "strata",
    ]
    cc = primary.dropna(subset=required).copy()
    cc["mcop_quartile"] = pd.cut(
        cc["mcop_log2"],
        bins=[-np.inf, *cutpoints[1:-1], np.inf],
        labels=False,
        include_lowest=True,
    ) + 1
    rows = []
    for q in [1, 2, 3, 4]:
        frame = cc[cc["mcop_quartile"].eq(q)]
        weights = pd.to_numeric(frame["pooled_weight"], errors="coerce")
        valid = weights.notna() & weights.gt(0)
        rows.append(
            {
                "Quartile": f"Q{q}",
                "N": int(len(frame)),
                "CRC_cases": int(frame["outcome"].sum()),
                "Controls": int(len(frame) - frame["outcome"].sum()),
                "CRC_pct_unweighted": float(100.0 * frame["outcome"].mean()) if len(frame) else np.nan,
                "CRC_prevalence_weighted": float(
                    (weights[valid] * frame.loc[valid, "outcome"]).sum() / weights[valid].sum()
                ) if valid.any() else np.nan,
                "MCOP_log2_min": float(frame["mcop_log2"].min()),
                "MCOP_log2_max": float(frame["mcop_log2"].max()),
                "MCOP_ng_mL_median": float(frame["URXCOP"].median()),
            }
        )
    meta = {
        "cutpoints_log2": [float(value) for value in cutpoints],
        "complete_case_N": int(len(cc)),
        "complete_case_CRC_cases": int(cc["outcome"].sum()),
    }
    return pd.DataFrame(rows), meta


def assay_table() -> pd.DataFrame:
    return pd.DataFrame(
        [{"cycle": cycle, **meta} for cycle, meta in ASSAY_METADATA.items()]
    )


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

    summary = cycle_summary(primary)
    summary.to_csv(args.outdir / "mcop_crc_phase2_cycle_heterogeneity_summary.csv", index=False)
    quartiles, quartile_meta = quartile_case_counts(primary)
    quartiles.to_csv(args.outdir / "mcop_crc_phase2_quartile_case_counts.csv", index=False)
    assay = assay_table()
    assay.to_csv(args.outdir / "mcop_crc_phase2_assay_lod_audit.csv", index=False)

    row_2011 = summary.loc[summary["cycle"].eq("2011-2012")].iloc[0].to_dict()
    pooled_median = float(primary.loc[primary["mcop_log2"].notna(), "URXCOP"].median())
    report = [
        "# MCOP-CRC Phase 2C：cycle heterogeneity audit",
        "",
        "本轮只解释 MCOP cycle heterogeneity，不更换候选、不增加机制分析。所有人口学和 CRC 结构描述使用 MCOP 有效且 CRC outcome 为 CRC 或 cancer-free 的 primary population；模型型 quartile case counts 使用与前一轮相同的完整协变量 complete-case frame。",
        "",
        "## 主要发现",
        "",
        f"- 2011-2012 的 primary MCOP median={row_2011['MCOP_median_ng_mL']:.4g} ng/mL，Q1-Q3={row_2011['MCOP_Q1_ng_mL']:.4g}-{row_2011['MCOP_Q3_ng_mL']:.4g}；CRC cases={int(row_2011['CRC_cases'])}，CRC median age={row_2011['Age_median_CRC']:.4g}.",
        f"- 2011-2012 MCOP above LOD={row_2011['MCOP_above_LOD_pct']:.2f}%，codebook LLOD={row_2011['LLOD_ng_mL_codebook']:.3g} ng/mL；因此 2011-2012 的反向 OR 不是由最高 LOD/censoring 直接解释。",
        f"- 在与主模型相同的 complete-case frame 中，2011-2012 CRC cases={int(row_2011['Model_complete_case_CRC_cases'])}，病例 MCOP median={row_2011['Model_CRC_MCOP_median_ng_mL']:.4g} ng/mL，对照 median={row_2011['Model_control_MCOP_median_ng_mL']:.4g} ng/mL，病例/对照 median ratio={row_2011['Model_case_control_MCOP_median_ratio']:.3f}；这与该周期未调整 OR<1 的方向一致。",
        f"- Pooled primary MCOP median={pooled_median:.4g} ng/mL.",
        "- 公开代码本的方法描述均属于 HPLC-ESI-MS/MS 类平台，但这不能证明跨周期完全没有校准/批次尺度差异。可见的 documented assay-scale change 是 MCOP LLOD：2005-2008 为 0.7 ng/mL，2009-2012 为 0.2，2013-2018 为 0.3。",
        "- CDC explicitly states no lab method, equipment, or site changes for 2013-2014, 2015-2016, and 2017-2018; the 2011-2012 manual is Method 6306.04 and explicitly maps MCOP to the di-isononyl phthalate axis.",
        "- 因此当前最合理的解释不是单一 LOD 故障，而是 2011-2012 暴露分布整体偏高、该周期仅 10 个 complete-case CRC cases，且病例与对照的暴露排序反向；年龄、种族、吸烟和肌酐构成差异可能进一步改变调整后效应。仅凭当前数据不能把其中任何一个因素定为唯一原因。",
        "",
        "## Quartile case counts",
        "",
        "| Quartile | N | CRC cases | Controls | Unweighted CRC % | Weighted CRC prevalence |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in quartiles.itertuples(index=False):
        report.append(
            f"| {row.Quartile} | {row.N} | {row.CRC_cases} | {row.Controls} | {row.CRC_pct_unweighted:.3f}% | {row.CRC_prevalence_weighted:.6g} |"
        )
    report += [
        "",
        f"Quartile complete-case frame: N={quartile_meta['complete_case_N']}, CRC cases={quartile_meta['complete_case_CRC_cases']}; log2 cutpoints={';'.join(f'{x:.6g}' for x in quartile_meta['cutpoints_log2'])}.",
        "",
        "## Files",
        "",
        "- mcop_crc_phase2_cycle_heterogeneity_summary.csv",
        "- mcop_crc_phase2_assay_lod_audit.csv",
        "- mcop_crc_phase2_quartile_case_counts.csv",
        "",
        "官方依据：[NHANES 2005-2006 PHTHTE_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PHTHTE_D.htm)、[2011-2012 PHTHTE_G](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PHTHTE_G.htm)、[2011-2012 laboratory manual](https://wwwn.cdc.gov/nchs/data/nhanes/public/2011/labmethods/phthte_g_met.pdf)、[2013-2014 PHTHTE_H](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/PHTHTE_H.htm)、[2015-2016 PHTHTE_I](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PHTHTE_I.htm)、[2017-2018 PHTHTE_J](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm)。",
        "",
        f"Run timestamp (UTC): {datetime.now(timezone.utc).isoformat()}",
    ]
    (args.outdir / "mcop_crc_phase2_heterogeneity_report.md").write_text("\n".join(report), encoding="utf-8")

    manifest = {
        "analysis": "MCOP-CRC Phase 2C cycle heterogeneity audit",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_harmonized": str(args.harmonized),
        "population": "CRC vs cancer-free controls, age >=20, MCOP available",
        "outputs": [
            "mcop_crc_phase2_cycle_heterogeneity_summary.csv",
            "mcop_crc_phase2_assay_lod_audit.csv",
            "mcop_crc_phase2_quartile_case_counts.csv",
        ],
        "quartile_meta": quartile_meta,
        "assay_sources": sorted(meta["codebook"] for meta in ASSAY_METADATA.values()),
    }
    (args.outdir / "mcop_crc_phase2_heterogeneity_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(quartiles.to_string(index=False))


if __name__ == "__main__":
    main()
