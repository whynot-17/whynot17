from __future__ import annotations

import json
import pathlib
import sys
from datetime import date, datetime, timezone

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
import validate_113_across_geo_replacement as base  # noqa: E402

BASE = ROOT / "results" / "2026-10-10" / "geo_replication_113_replacement"
OUT = ROOT / "results" / "2026-10-10" / "GSE41258_paired_44_recheck"
GSE = "GSE41258"


def direction(fc):
    if not np.isfinite(fc):
        return "not represented"
    return "Tumor_up" if fc > 0 else "Tumor_down" if fc < 0 else "No_change"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base.CACHE = ROOT / "work" / "geo_cache_replacement_2026-10-10"
    input_df = pd.read_csv(base.INPUT_113, sep="\t")
    targets = sorted(set(input_df.gene_symbol.astype(str).str.strip().str.upper()) - {""})
    if len(targets) != 113:
        raise ValueError(f"Expected 113 candidates, found {len(targets)}")

    metadata, table, matrix_hash = base.load_series_matrix(GSE)
    records = base.assign_groups(GSE, base.sample_records(metadata))
    sample_columns = list(table.columns[1:])
    records_by_gsm = {r["gsm"]: r for r in records}
    if set(sample_columns) - set(records_by_gsm):
        raise ValueError("Expression columns are missing GEO sample metadata")
    records = [records_by_gsm[gsm] for gsm in sample_columns]
    selected = [i for i, r in enumerate(records) if r["group"] in {"Tumor", "Normal"}]
    pair_ids = base.complete_pair_ids(records, selected)
    if len(pair_ids) != 44:
        raise ValueError(f"Expected 44 complete pairs, found {len(pair_ids)}")

    platform = base.PLATFORM_BY_GSE[GSE]
    mapping, annotation_hash = base.platform_mapping(platform)
    genes, matrix, probes_per_gene, transformed = base.map_probe_matrix(table, mapping, GSE)
    fc, p, n_pairs, used_pairs = base.paired_test(matrix, records, selected)
    q = base.adjust_bh(p)
    if n_pairs != 44:
        raise ValueError(f"Paired t-test used {n_pairs} pairs")
    idx = {g: i for i, g in enumerate(genes)}
    all_gene = pd.DataFrame({
        "gene_symbol": genes,
        "n_complete_pairs": n_pairs,
        "probes_collapsed": [probes_per_gene[g] for g in genes],
        "log2FC_tumor_minus_normal": fc,
        "P_value": p,
        "FDR_BH_all_mapped_genes": q,
    }).sort_values(["FDR_BH_all_mapped_genes", "gene_symbol"])
    all_gene.to_csv(OUT / "GSE41258_all_mapped_genes_paired_DE.tsv", sep="\t", index=False)

    rows = []
    for gene in targets:
        if gene not in idx:
            rows.append({
                "gene_symbol": gene, "accession": GSE, "platform": platform,
                "n_tumor": 182, "n_normal": 53, "n_pairs_used": n_pairs,
                "statistical_method": "paired t-test; replicate arrays averaged within patient and tissue",
                "probes_collapsed": 0, "log2FC_tumor_minus_normal": np.nan,
                "P_value": np.nan, "FDR_BH_all_mapped_genes": np.nan,
                "passes_FDR_and_abs_log2FC": False, "direction": "not represented",
            })
            continue
        j = idx[gene]
        value_p = float(p[j]) if np.isfinite(p[j]) else np.nan
        value_q = float(q[j]) if np.isfinite(q[j]) else np.nan
        value_fc = float(fc[j])
        rows.append({
            "gene_symbol": gene, "accession": GSE, "platform": platform,
            "n_tumor": 182, "n_normal": 53, "n_pairs_used": n_pairs,
            "statistical_method": "paired t-test; replicate arrays averaged within patient and tissue",
            "probes_collapsed": probes_per_gene[gene],
            "log2FC_tumor_minus_normal": value_fc,
            "P_value": value_p, "FDR_BH_all_mapped_genes": value_q,
            "passes_FDR_and_abs_log2FC": bool(
                np.isfinite(value_q) and value_q < base.STABLE_FDR
                and abs(value_fc) > base.STABLE_ABS_LOG2FC
            ),
            "direction": direction(value_fc),
        })
    paired_candidates = pd.DataFrame(rows).sort_values("gene_symbol")
    paired_candidates.to_csv(OUT / "GSE41258_113_candidates_paired_DE.tsv", sep="\t", index=False)

    previous = pd.read_csv(BASE / "113_gene_per_cohort_DE_results.tsv", sep="\t")
    previous_summary = pd.read_csv(BASE / "113_gene_stability_summary.tsv", sep="\t")
    old_gse = previous[previous.accession.eq(GSE)].copy()
    if len(old_gse) != 113:
        raise ValueError(f"Expected 113 previous GSE41258 results, found {len(old_gse)}")
    old29 = old_gse[old_gse.passes_FDR_and_abs_log2FC].copy()
    recheck29 = old29.merge(
        paired_candidates[[
            "gene_symbol", "log2FC_tumor_minus_normal", "P_value",
            "FDR_BH_all_mapped_genes", "passes_FDR_and_abs_log2FC", "direction"
        ]],
        on="gene_symbol", suffixes=("_original_model", "_paired_44"), how="left",
    )
    recheck29["direction_preserved"] = (
        recheck29.direction_original_model.eq(recheck29.direction_paired_44)
    )
    recheck29["stable_same_threshold_and_direction"] = (
        recheck29.passes_FDR_and_abs_log2FC_paired_44 & recheck29.direction_preserved
    )
    recheck29.to_csv(OUT / "original_model_29_candidates_rechecked.tsv", sep="\t", index=False)

    # Replace only GSE41258 rows; retain the other three cohort estimates.
    updated = pd.concat([
        previous[~previous.accession.eq(GSE)],
        paired_candidates[previous.columns],
    ], ignore_index=True).sort_values(["gene_symbol", "accession"])
    summary_rows = []
    for gene in targets:
        block = updated[updated.gene_symbol.eq(gene)]
        hits = block[block.passes_FDR_and_abs_log2FC]
        directions = set(hits.direction) - {"not represented", "No_change"}
        consistent = len(directions) <= 1
        row = {
            "gene_symbol": gene,
            "significant_cohorts_FDR_lt_0.05_abs_log2FC_gt_1": len(hits),
            "cohorts_meeting_threshold": ";".join(hits.accession.tolist()),
            "significant_up_cohorts": int(hits.direction.eq("Tumor_up").sum()),
            "significant_down_cohorts": int(hits.direction.eq("Tumor_down").sum()),
            "direction_consistent_among_significant": consistent,
            "replicated_in_at_least_3_of_4": len(hits) >= 3 and consistent,
            "replicated_in_all_4": len(hits) == 4 and consistent,
        }
        for accession in base.ACCESSIONS:
            item = block[block.accession.eq(accession)].iloc[0]
            row[f"{accession}_log2FC"] = item.log2FC_tumor_minus_normal
            row[f"{accession}_FDR"] = item.FDR_BH_all_mapped_genes
            row[f"{accession}_passes"] = bool(item.passes_FDR_and_abs_log2FC)
        summary_rows.append(row)
    updated_summary = pd.DataFrame(summary_rows)
    stable3 = updated_summary[updated_summary.replicated_in_at_least_3_of_4].copy()
    stable4 = updated_summary[updated_summary.replicated_in_all_4].copy()
    updated.to_csv(OUT / "updated_113_gene_per_cohort_results.tsv", sep="\t", index=False)
    updated_summary.to_csv(OUT / "updated_113_gene_stability_summary.tsv", sep="\t", index=False)
    stable3.to_csv(OUT / "updated_stable_genes_3_of_4.tsv", sep="\t", index=False)
    stable4.to_csv(OUT / "updated_stable_genes_4_of_4.tsv", sep="\t", index=False)

    old19 = set(previous_summary.loc[previous_summary.replicated_in_at_least_3_of_4, "gene_symbol"])
    old9 = set(previous_summary.loc[previous_summary.replicated_in_all_4, "gene_symbol"])
    new19 = set(stable3.gene_symbol)
    new9 = set(stable4.gene_symbol)
    status19 = previous_summary[previous_summary.gene_symbol.isin(old19)].merge(
        updated_summary[["gene_symbol", "significant_cohorts_FDR_lt_0.05_abs_log2FC_gt_1",
                        "cohorts_meeting_threshold", "replicated_in_at_least_3_of_4",
                        "replicated_in_all_4", "GSE41258_log2FC", "GSE41258_FDR", "GSE41258_passes"]],
        on="gene_symbol", suffixes=("_original", "_paired_44"),
    )
    status19.to_csv(OUT / "original_19_status_after_paired_reanalysis.tsv", sep="\t", index=False)
    status9 = previous_summary[previous_summary.gene_symbol.isin(old9)].merge(
        updated_summary[["gene_symbol", "significant_cohorts_FDR_lt_0.05_abs_log2FC_gt_1",
                        "cohorts_meeting_threshold", "replicated_in_all_4",
                        "GSE41258_log2FC", "GSE41258_FDR", "GSE41258_passes"]],
        on="gene_symbol", suffixes=("_original", "_paired_44"),
    )
    status9.to_csv(OUT / "original_9_four_of_four_status.tsv", sep="\t", index=False)

    pair_count_rows = []
    for patient in used_pairs:
        group_records = [r for r in records if r["pair_id"] == patient and r["group"] in {"Tumor", "Normal"}]
        pair_count_rows.append({
            "patient_id": patient,
            "tumor_samples": sum(r["group"] == "Tumor" for r in group_records),
            "normal_samples": sum(r["group"] == "Normal" for r in group_records),
        })
    pd.DataFrame(pair_count_rows).to_csv(OUT / "GSE41258_44_pair_sample_counts.tsv", sep="\t", index=False)

    old29_names = set(old29.gene_symbol)
    new_gse_hits = set(paired_candidates.loc[paired_candidates.passes_FDR_and_abs_log2FC, "gene_symbol"])
    n29_sig = int(recheck29.passes_FDR_and_abs_log2FC_paired_44.sum())
    n29_stable = int(recheck29.stable_same_threshold_and_direction.sum())
    old_fc = recheck29["log2FC_tumor_minus_normal_original_model"]
    paired_fc = recheck29["log2FC_tumor_minus_normal_paired_44"]
    fc_pearson = float(old_fc.corr(paired_fc, method="pearson"))
    fc_spearman = float(old_fc.corr(paired_fc, method="spearman"))
    fc_median_abs_difference = float((old_fc - paired_fc).abs().median())
    summary = {
        "analysis_date": date.today().isoformat(),
        "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "accession": GSE,
        "pairs_tested": n_pairs,
        "candidate_genes": len(targets),
        "mapped_genes_used_for_BH_FDR": len(genes),
        "candidate_platform_coverage": int(paired_candidates.probes_collapsed.gt(0).sum()),
        "threshold": {"BH_FDR_lt": base.STABLE_FDR, "abs_log2FC_gt": base.STABLE_ABS_LOG2FC},
        "previous_model_gse41258_hits": len(old29),
        "previous_29_still_significant": n29_sig,
        "previous_29_still_significant_same_direction": n29_stable,
        "previous_29_log2fc_pearson_correlation": fc_pearson,
        "previous_29_log2fc_spearman_correlation": fc_spearman,
        "previous_29_log2fc_median_absolute_difference": fc_median_abs_difference,
        "old_29_lost_threshold": sorted(old29_names - new_gse_hits),
        "new_gse41258_hits_not_in_old_29": sorted(new_gse_hits - old29_names),
        "old_19_list": sorted(old19),
        "updated_3_of_4_list": sorted(new19),
        "old_19_no_longer_in_3_of_4": sorted(old19 - new19),
        "newly_in_3_of_4": sorted(new19 - old19),
        "old_9_four_of_four": sorted(old9),
        "updated_4_of_4_list": sorted(new9),
        "old_9_no_longer_4_of_4": sorted(old9 - new9),
        "newly_in_4_of_4": sorted(new9 - old9),
        "checksums": {
            "GSE41258_series_matrix_sha256": matrix_hash,
            "GPL96_annotation_sha256": annotation_hash,
            "input_113_sha256": base.sha256(base.INPUT_113),
        },
    }
    (OUT / "analysis_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    pd.DataFrame([{
        "total_series_matrix_samples": len(records),
        "selected_tumor_samples": sum(r["group"] == "Tumor" for r in records),
        "selected_normal_samples": sum(r["group"] == "Normal" for r in records),
        "complete_pairs_tested": n_pairs,
        "group_means_tested": 2 * n_pairs,
        "mapped_genes_for_BH_FDR": len(genes),
        "candidate_platform_coverage": int(paired_candidates.probes_collapsed.gt(0).sum()),
        "test": "paired t-test after averaging replicate arrays within patient and tissue",
        "expression_transform": "log2(PLIER signal)" if transformed else "none",
    }]).to_csv(OUT / "GSE41258_paired_model_counts.tsv", sep="\t", index=False)

    lines = [
        "# GSE41258 44 对患者配对分析与复现名单复核",
        "",
        f"- 分析日期：{date.today().isoformat()}",
        "- 本次只重算 GSE41258；其他三个队列沿用本项目上一轮的结果。",
        "- GSE41258/GPL96 Series Matrix 和平台注释来自 NCBI GEO；按 included in analysis=Yes 选 Primary Tumor 与 Normal Colon。",
        "- 对同一患者、同一组织的重复芯片先取表达均值，再对 44 个完整配对做 paired t-test；PLIER 信号按 log2 转换。",
        "- 多 probe 按 gene symbol 取表达中位数；BH-FDR 在全部可映射基因上重新计算。阈值保持 FDR<0.05 且 |log2FC|>1。",
        "",
        "## 原 GSE41258 模型的 29 个显著候选",
        "",
        f"- 新配对分析仍达到阈值：{n29_sig}/29。",
        f"- 同时达到阈值且方向与原模型一致：{n29_stable}/29。",
        f"- 原模型与配对模型的 29 个 log2FC 高度相关：Pearson r={fc_pearson:.3f}，Spearman ρ={fc_spearman:.3f}；log2FC 绝对差异中位数={fc_median_abs_difference:.3f}。",
        f"- 未达到阈值：{', '.join(summary['old_29_lost_threshold']) or '无'}。",
        f"- 新配对模型达到阈值、但不在原 29 个中的候选：{', '.join(summary['new_gse41258_hits_not_in_old_29']) or '无'}。",
        "",
        "## 原 19 个 3/4 基因",
        "",
        f"- 原名单仍保留：{len(old19 & new19)}/19。",
        f"- 移出：{', '.join(sorted(old19 - new19)) or '无'}。",
        f"- 新增：{', '.join(sorted(new19 - old19)) or '无'}。",
        f"- 更新名单（{len(new19)} 个）：{', '.join(sorted(new19)) or '无'}。",
        "",
        "## 原 9 个 4/4 基因",
        "",
        f"- 仍满足 4/4：{len(old9 & new9)}/9。",
        f"- 不再满足：{', '.join(sorted(old9 - new9)) or '无'}。",
        f"- 更新 4/4 名单（{len(new9)} 个）：{', '.join(sorted(new9)) or '无'}。",
        "",
        "## 文件",
        "",
        "- GSE41258_113_candidates_paired_DE.tsv：候选基因配对检验结果。",
        "- original_model_29_candidates_rechecked.tsv：原 29 个候选逐项比较。",
        "- original_19_status_after_paired_reanalysis.tsv、original_9_four_of_four_status.tsv：原名单复核。",
        "- updated_stable_genes_3_of_4.tsv、updated_stable_genes_4_of_4.tsv：更新后的名单。",
        "- updated_113_gene_per_cohort_results.tsv、updated_113_gene_stability_summary.tsv：替换单队列后完整汇总。",
        "- GSE41258_all_mapped_genes_paired_DE.tsv：全平台基因的 P 值与 BH-FDR。",
        "- GSE41258_44_pair_sample_counts.tsv：44 对的患者内重复芯片数。",
        "- analysis_summary.json：关键数字和输入数据校验和。",
        "- 可复现脚本：analysis/recheck_gse41258_44_paired.py。",
        "",
        "数据记录：[NCBI GEO GSE41258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41258)。",
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Output: {OUT}")


if __name__ == "__main__":
    main()
