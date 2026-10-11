from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/2026-10-11/single_cell_16genes_GSE132465_GSE200997"
INPUT = RESULTS / "expression_by_sample_cell_type.tsv"
WORK = ROOT / "work/single_cell_GSE132465_GSE200997"
OUTPUT = RESULTS / "sample_level_paired_tests.tsv"
REPORT = RESULTS / "sample_level_paired_test_report.md"


def bh_adjust(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted p values, preserving input index."""
    p = values.astype(float)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    n = len(ranked)
    adjusted_sorted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    adjusted = np.empty(n, dtype=float)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return pd.Series(adjusted, index=values.index)


def exact_paired_wilcoxon(differences: np.ndarray) -> float:
    """Two-sided exact signed-rank permutation p, enumerating signs."""
    d = np.asarray(differences, dtype=float)
    d = d[np.isfinite(d)]
    d = d[d != 0]
    if len(d) == 0:
        return 1.0
    ranks = pd.Series(np.abs(d)).rank(method="average").to_numpy()
    observed = abs(float(np.sum(np.sign(d) * ranks)))
    extreme = 0
    total = 1 << len(d)
    # At most 10 paired patients here; exhaustive sign permutations are exact and cheap.
    for signs in itertools.product((-1.0, 1.0), repeat=len(d)):
        statistic = abs(float(np.dot(np.asarray(signs), ranks)))
        if statistic >= observed - 1e-12:
            extreme += 1
    return extreme / total


def paired_samples(cohort: str, sample_conditions: pd.DataFrame) -> dict[str, dict[str, str]]:
    if cohort == "GSE132465":
        ann = pd.read_csv(
            WORK / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
            sep="\t", usecols=["Patient", "Class", "Sample"],
        ).drop_duplicates()
        patient_samples = {}
        for patient, group in ann.groupby("Patient", sort=True):
            by_condition = group.groupby("Class")["Sample"].agg(lambda x: sorted(set(x)))
            if {"Normal", "Tumor"}.issubset(by_condition.index):
                normals, tumors = by_condition["Normal"], by_condition["Tumor"]
                if len(normals) == 1 and len(tumors) == 1:
                    patient_samples[str(patient)] = {"Normal": normals[0], "Tumor": tumors[0]}
    elif cohort == "GSE200997":
        ann = pd.read_csv(
            WORK / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz",
            usecols=["samples", "Condition"],
        ).drop_duplicates()
        ann["patient"] = ann["samples"].str.extract(r"^[BT]_cac(\d+)$", expand=False)
        patient_samples = {}
        for patient, group in ann.dropna(subset=["patient"]).groupby("patient", sort=True):
            by_condition = group.groupby("Condition")["samples"].agg(lambda x: sorted(set(x)))
            if {"Normal", "Tumor"}.issubset(by_condition.index):
                normals, tumors = by_condition["Normal"], by_condition["Tumor"]
                if len(normals) == 1 and len(tumors) == 1:
                    patient_samples[str(patient)] = {"Normal": normals[0], "Tumor": tumors[0]}
    else:
        raise ValueError(f"Unexpected cohort: {cohort}")

    available = sample_conditions.groupby("sample")["condition"].first().to_dict()
    return {
        patient: samples for patient, samples in patient_samples.items()
        if all(available.get(samples[condition]) == condition for condition in ("Normal", "Tumor"))
    }


def main() -> None:
    data = pd.read_csv(INPUT, sep="\t")
    required = {"cohort", "sample", "condition", "cell_type", "gene", "n_cells", "mean_log1p_CP10K"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing required input columns: {sorted(missing)}")

    rows = []
    pairing_summary = []
    for cohort, cohort_data in data.groupby("cohort", sort=True):
        sample_conditions = cohort_data[["sample", "condition"]].drop_duplicates()
        pairs = paired_samples(cohort, sample_conditions)
        pairing_summary.append((cohort, len(pairs), ", ".join(sorted(pairs))))
        for (cell_type, gene), group in cohort_data.groupby(["cell_type", "gene"], sort=True):
            lookup = group.set_index(["sample", "condition"])["mean_log1p_CP10K"].to_dict()
            paired_values = []
            for patient, sample_map in pairs.items():
                normal_key = (sample_map["Normal"], "Normal")
                tumor_key = (sample_map["Tumor"], "Tumor")
                if normal_key in lookup and tumor_key in lookup:
                    normal = float(lookup[normal_key])
                    tumor = float(lookup[tumor_key])
                    paired_values.append((patient, normal, tumor, tumor - normal))
            if not paired_values:
                continue

            vals = np.asarray([x[3] for x in paired_values], dtype=float)
            rows.append({
                "cohort": cohort,
                "gene": gene,
                "cell_type": cell_type,
                "n_complete_pairs": len(vals),
                "n_nonzero_pair_differences": int(np.count_nonzero(vals)),
                "n_tumor_higher_pairs": int(np.sum(vals > 0)),
                "n_normal_higher_pairs": int(np.sum(vals < 0)),
                "n_tied_pairs": int(np.sum(vals == 0)),
                "test": "two-sided exact paired Wilcoxon signed-rank (sample-level means)",
                "normal_mean_log1p_CP10K": float(np.mean([x[1] for x in paired_values])),
                "tumor_mean_log1p_CP10K": float(np.mean([x[2] for x in paired_values])),
                "mean_paired_difference_tumor_minus_normal": float(np.mean(vals)),
                "median_paired_difference_tumor_minus_normal": float(np.median(vals)),
                "p_value": exact_paired_wilcoxon(vals),
                "paired_patient_ids": ";".join(x[0] for x in paired_values),
            })

    out = pd.DataFrame(rows)
    if out.empty:
        raise RuntimeError("No complete paired comparisons were available")
    out["q_value_BH_global"] = bh_adjust(out["p_value"])
    out["q_value_BH_within_cohort"] = out.groupby("cohort", group_keys=False)["p_value"].apply(bh_adjust)
    out["p_lt_0_05"] = out["p_value"] < 0.05
    out["q_global_lt_0_05"] = out["q_value_BH_global"] < 0.05
    out["q_within_cohort_lt_0_05"] = out["q_value_BH_within_cohort"] < 0.05
    out = out.sort_values(["q_value_BH_global", "p_value", "cohort", "cell_type", "gene"]).reset_index(drop=True)
    out.to_csv(OUTPUT, sep="\t", index=False, float_format="%.10g")

    global_sig = out[out["q_global_lt_0_05"]].sort_values(["cohort", "cell_type", "gene"])
    primary_sig = out[(out["cohort"] == "GSE132465") & out["q_within_cohort_lt_0_05"]].sort_values(["cell_type", "gene"])
    primary_genes = primary_sig.groupby("gene", sort=True)
    cohort_counts = out.groupby("cohort")[["p_lt_0_05", "q_global_lt_0_05", "q_within_cohort_lt_0_05"]].sum()
    lines = [
        "# 16个基因单细胞队列的样本级配对检验",
        "",
        "## 统计单位与配对",
        "",
        "以样本/患者为独立统计单位；每个样本、细胞类型、基因先汇总为该细胞类型所有细胞的平均 `log1p(CP10K)`，随后进行双侧精确配对 Wilcoxon 符号秩检验。未把细胞当作独立重复。仅使用该患者该细胞类型在肿瘤和正常样本中均有细胞的完整配对；缺失细胞类型不当作零表达。",
        "",
        "- GSE132465：从 GEO 注释中的 `Patient`、`Class`、`Sample` 字段识别10对患者样本。",
        "- GSE200997：从注释中的 `B_cacN`/`T_cacN`样本编号识别7对正常/肿瘤患者样本；此对应也与GEO样本记录中的患者编号一致。",
        f"- 共检验 {len(out)} 个可配对的基因×细胞类型×队列组合；全局 BH-FDR 在全部组合上校正。另附每队列内BH校正作为敏感性结果。",
        "",
        "## 结果",
        "",
        f"名义 `p < 0.05` 的比较有 {int(out['p_lt_0_05'].sum())} 个；全局 BH 校正后 `q < 0.05` 有 {int(out['q_global_lt_0_05'].sum())} 个；分队列 BH 校正后 `q < 0.05` 有 {int(out['q_within_cohort_lt_0_05'].sum())} 个。",
        "",
        "全局 BH 校正后显著的结果：",
        "",
        "| 队列 | 基因 | 细胞类型 | 完整配对数 | 肿瘤−正常均值差（log1p CP10K） | 精确 p | 全局 BH q |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in global_sig.itertuples(index=False):
        lines.append(
            f"| {row.cohort} | {row.gene} | {row.cell_type} | {row.n_complete_pairs} | "
            f"{row.mean_paired_difference_tumor_minus_normal:.4f} | {row.p_value:.6g} | {row.q_value_BH_global:.6g} |"
        )
    lines += [
        "", "## 主分析：GSE132465",
        "",
        f"若将GSE132465作为主分析、GSE200997作为敏感性探索，GSE132465队列内BH校正后有 {len(primary_sig)} 个显著的基因×细胞类型比较；去重后是 {primary_sig['gene'].nunique()} 个不同基因（共16个算法基因）。",
        "",
        "| 基因 | 显著细胞类型 | 显著比较数 | 最小 p | 最小队列内 BH q |",
        "|---|---|---:|---:|---:|",
    ]
    for gene, group in primary_genes:
        lines.append(
            f"| {gene} | {', '.join(sorted(group['cell_type'].unique()))} | {len(group)} | "
            f"{group['p_value'].min():.6g} | {group['q_value_BH_within_cohort'].min():.6g} |"
        )
    lines += ["", "队列内检验计数（数字是基因×细胞类型比较，不是唯一基因数）：", "", "| 队列 | 名义 p<0.05 | 全局BH q<0.05 | 队列内BH q<0.05 |", "|---|---:|---:|---:|"]
    for cohort, row in cohort_counts.iterrows():
        lines.append(f"| {cohort} | {int(row['p_lt_0_05'])} | {int(row['q_global_lt_0_05'])} | {int(row['q_within_cohort_lt_0_05'])} |")
    lines += [
        "",
        f"全局BH的13个显著比较去重后为 {global_sig['gene'].nunique()} 个基因，全部来自GSE132465。GSE200997作为敏感性探索，有13个名义 `p<0.05` 比较，但没有通过队列内或全局BH校正的比较；部分细胞类型完整配对数仅2–4，检验效能有限。",
        "",
        "## 解读",
        "",
        "统计显著表示样本层面表达差异与零差异不一致，不等同于DINP机制证据。GSE200997细胞类型为标志基因推断标签；两队列的细胞类型标签并非完全等价。效应值为平均 `log1p(CP10K)` 差值，适合描述方向和尺度，不应直接解释为倍数变化。",
        "",
        "## 文件",
        "",
        "- `sample_level_paired_tests.tsv`：所有组合的配对数、方向、精确p值与BH校正q值。",
        "- `analysis/single_cell_sample_level_pvalues.py`：可复现分析脚本。",
        "",
        "队列来源：[GSE132465 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)；[GSE200997 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997)。",
        "",
    ]
    REPORT.write_text("\n".join(lines), encoding="utf-8")

    print(f"Input rows: {len(data):,}")
    print(f"Paired tests: {len(out):,}")
    print("Patient-level pairs verified from annotation sample IDs:")
    for cohort, n, ids in pairing_summary:
        print(f"  {cohort}: {n} pairs ({ids})")
    print("Complete pair counts by cohort/cell type:")
    print(out.groupby(["cohort", "cell_type"])["n_complete_pairs"].max().to_string())
    print("Nominal p < 0.05:", int(out["p_lt_0_05"].sum()))
    print("Global BH q < 0.05:", int(out["q_global_lt_0_05"].sum()))
    print("Within-cohort BH q < 0.05 comparisons:", int(out["q_within_cohort_lt_0_05"].sum()))
    print("Primary GSE132465 unique genes passing within-cohort BH:", primary_sig["gene"].nunique(), sorted(primary_sig["gene"].unique()))
    print("Wrote:", OUTPUT)
    print("Wrote:", REPORT)


if __name__ == "__main__":
    main()

