"""Phase 2F-A: targeted CRC disease-state validation of the DINP-axis program.

This is deliberately narrow. It tests only the frozen Phase 2E genes:

    PPAR/nuclear-receptor core: PPARA, PPARD, PPARG, NR1I2, NR1I3,
                                NR1H2, NR1H3
    inflammatory complement:    RELA, STAT3

The script uses the public UCSC Toil Xena hub, with the same harmonized
TCGA/GTEx RSEM gene-TPM dataset and its phenotype table. It compares:

1. TCGA COAD+READ primary tumors vs TCGA CRC solid-tissue normals;
2. TCGA COAD+READ primary tumors vs GTEx transverse/sigmoid colon normals.

It reports gene-level median shifts and two-sided Mann–Whitney tests, plus
sample-level pathway scores formed from within-contrast z-scored expression.
This is disease-state validation, not causal exposure mediation. Single-cell
validation is intentionally not silently substituted here: the companion
CELLxGENE Census step is recorded as pending when the local runtime cannot
install TileDB-SOMA.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

try:
    import xenaPython as xena
except ImportError as exc:  # pragma: no cover - exercised only without optional dependency
    raise SystemExit(
        "xenaPython is required. Install with: pip install xenaPython==1.0.14"
    ) from exc


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs"

XENA_HUB = "https://toil.xenahubs.net"
EXPRESSION_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"

GENES = ["PPARA", "PPARD", "PPARG", "NR1I2", "NR1I3", "NR1H2", "NR1H3", "RELA", "STAT3"]
NR_GENES = ["PPARA", "PPARD", "PPARG", "NR1I2", "NR1I3", "NR1H2", "NR1H3"]
INFLAMMATORY_GENES = ["RELA", "STAT3"]


def bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").fillna(1.0).to_numpy(float)
    order = np.argsort(p)
    adjusted = np.ones(len(p), dtype=float)
    running = 1.0
    m = max(len(p), 1)
    for rank in range(len(p) - 1, -1, -1):
        idx = order[rank]
        running = min(running, p[idx] * m / (rank + 1))
        adjusted[idx] = running
    return pd.Series(adjusted, index=values.index)


def decode_codes(field: str) -> list[str]:
    raw = xena.field_codes(XENA_HUB, PHENOTYPE_DATASET, [field])[0]["code"]
    return str(raw).split("\t")


def decode_value(value: object, codes: list[str]) -> str:
    if value is None or str(value).lower() in {"nan", "none"}:
        return "NaN"
    try:
        return codes[int(float(value))]
    except (ValueError, TypeError, IndexError):
        return str(value)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def load_sample_metadata() -> pd.DataFrame:
    samples = xena.dataset_samples(XENA_HUB, EXPRESSION_DATASET, None)
    fields = ["_study", "primary disease or tissue", "_sample_type"]
    values = [
        xena.dataset_probe_values(XENA_HUB, PHENOTYPE_DATASET, samples, [field])[1][0]
        for field in fields
    ]
    code_map = {field: decode_codes(field) for field in fields}
    metadata = pd.DataFrame({"sample_id": samples})
    for field, field_values in zip(fields, values):
        metadata[field] = [decode_value(value, code_map[field]) for value in field_values]

    metadata["group"] = "outside_scope"
    tcga_crc = (
        metadata["_study"].eq("TCGA")
        & metadata["primary disease or tissue"].isin(["Colon Adenocarcinoma", "Rectum Adenocarcinoma"])
    )
    metadata.loc[tcga_crc & metadata["_sample_type"].eq("Primary Tumor"), "group"] = "TCGA_CRC_primary_tumor"
    metadata.loc[tcga_crc & metadata["_sample_type"].eq("Solid Tissue Normal"), "group"] = "TCGA_CRC_solid_normal"
    gtex_colon = (
        metadata["_study"].eq("GTEX")
        & metadata["primary disease or tissue"].isin(["Colon - Transverse", "Colon - Sigmoid"])
        & metadata["_sample_type"].eq("Normal Tissue")
    )
    metadata.loc[gtex_colon, "group"] = "GTEx_colon_normal"
    metadata["include"] = metadata["group"].ne("outside_scope")
    metadata["dataset"] = EXPRESSION_DATASET
    metadata["patient_id"] = metadata["sample_id"].where(
        metadata["_study"].eq("TCGA"),
        np.nan,
    ).astype("string").str.split("-").str[:3].str.join("-")
    return metadata


def load_expression(sample_ids: list[str]) -> pd.DataFrame:
    values = xena.dataset_gene_probe_avg(XENA_HUB, EXPRESSION_DATASET, sample_ids, GENES)
    expression = pd.DataFrame(index=sample_ids)
    for record in values:
        gene = record["gene"]
        scores = record.get("scores", [[]])
        expression[gene] = scores[0] if scores and scores[0] else np.nan
    expression.index.name = "sample_id"
    return expression.apply(pd.to_numeric, errors="coerce")


def compare_gene_sets(expression: pd.DataFrame, metadata: pd.DataFrame, contrast: str, normal_group: str) -> pd.DataFrame:
    tumor_ids = metadata.loc[metadata["group"].eq("TCGA_CRC_primary_tumor"), "sample_id"]
    normal_ids = metadata.loc[metadata["group"].eq(normal_group), "sample_id"]
    rows = []
    for gene in GENES:
        tumor = expression.loc[expression.index.intersection(tumor_ids), gene].dropna().to_numpy(float)
        normal = expression.loc[expression.index.intersection(normal_ids), gene].dropna().to_numpy(float)
        if len(tumor) == 0 or len(normal) == 0:
            p_value = 1.0
            statistic = np.nan
        else:
            statistic, p_value = mannwhitneyu(tumor, normal, alternative="two-sided")
        rows.append({
            "contrast": contrast,
            "gene": gene,
            "tumor_n": len(tumor),
            "normal_n": len(normal),
            "tumor_median": float(np.median(tumor)) if len(tumor) else np.nan,
            "normal_median": float(np.median(normal)) if len(normal) else np.nan,
            "tumor_mean": float(np.mean(tumor)) if len(tumor) else np.nan,
            "normal_mean": float(np.mean(normal)) if len(normal) else np.nan,
            "median_delta_tumor_minus_normal": float(np.median(tumor) - np.median(normal)) if len(tumor) and len(normal) else np.nan,
            "mann_whitney_U": float(statistic) if np.isfinite(statistic) else np.nan,
            "p_value": float(p_value),
        })
    result = pd.DataFrame(rows)
    result["BH_FDR_within_contrast"] = bh(result["p_value"])
    return result


def score_samples(expression: pd.DataFrame, metadata: pd.DataFrame, contrast: str, normal_group: str) -> pd.DataFrame:
    ids = metadata.loc[metadata["group"].isin(["TCGA_CRC_primary_tumor", normal_group]), "sample_id"]
    ids = expression.index.intersection(ids)
    values = expression.loc[ids, GENES].copy()
    z = (values - values.mean(axis=0)) / values.std(axis=0, ddof=1).replace(0, np.nan)
    scored = pd.DataFrame(index=ids)
    scored["contrast"] = contrast
    scored["sample_id"] = ids
    scored["group"] = metadata.set_index("sample_id").loc[ids, "group"].to_numpy()
    scored["PPAR_nuclear_receptor_score"] = z[NR_GENES].mean(axis=1)
    scored["inflammatory_RELA_STAT3_score"] = z[INFLAMMATORY_GENES].mean(axis=1)
    scored["DINP_axis_9_gene_score"] = z[GENES].mean(axis=1)
    return scored.reset_index(drop=True)


def compare_scores(score_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for contrast, subset in score_df.groupby("contrast", sort=False):
        for score_name in ["PPAR_nuclear_receptor_score", "inflammatory_RELA_STAT3_score", "DINP_axis_9_gene_score"]:
            tumor = subset.loc[subset["group"].eq("TCGA_CRC_primary_tumor"), score_name].dropna().to_numpy(float)
            normal = subset.loc[~subset["group"].eq("TCGA_CRC_primary_tumor"), score_name].dropna().to_numpy(float)
            statistic, p_value = mannwhitneyu(tumor, normal, alternative="two-sided")
            rows.append({
                "contrast": contrast,
                "score": score_name,
                "tumor_n": len(tumor),
                "normal_n": len(normal),
                "tumor_median": float(np.median(tumor)),
                "normal_median": float(np.median(normal)),
                "median_delta_tumor_minus_normal": float(np.median(tumor) - np.median(normal)),
                "mann_whitney_U": float(statistic),
                "p_value": float(p_value),
            })
    result = pd.DataFrame(rows)
    result["BH_FDR_within_contrast"] = result.groupby("contrast", group_keys=False)["p_value"].apply(bh).reset_index(drop=True)
    return result


def paired_tcga_stats(expression: pd.DataFrame, metadata: pd.DataFrame, score_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Use patient-matched TCGA tumor/solid-normal pairs when available."""
    tcga = metadata[metadata["group"].isin(["TCGA_CRC_primary_tumor", "TCGA_CRC_solid_normal"])].copy()
    tumor = tcga[tcga["group"].eq("TCGA_CRC_primary_tumor")].drop_duplicates("patient_id")
    normal = tcga[tcga["group"].eq("TCGA_CRC_solid_normal")].drop_duplicates("patient_id")
    pair_ids = sorted(set(tumor["patient_id"].dropna()) & set(normal["patient_id"].dropna()))
    tumor_map = tumor.set_index("patient_id")["sample_id"]
    normal_map = normal.set_index("patient_id")["sample_id"]
    rows = []
    for gene in GENES:
        tumor_values = expression.loc[[tumor_map[p] for p in pair_ids], gene].to_numpy(float)
        normal_values = expression.loc[[normal_map[p] for p in pair_ids], gene].to_numpy(float)
        keep = np.isfinite(tumor_values) & np.isfinite(normal_values)
        differences = tumor_values[keep] - normal_values[keep]
        if len(differences) == 0 or np.allclose(differences, 0):
            statistic, p_value = np.nan, 1.0
        else:
            statistic, p_value = wilcoxon(differences, alternative="two-sided", method="auto")
        rows.append({
            "comparison": "TCGA_paired_primary_vs_solid_normal",
            "gene": gene,
            "paired_n": int(len(differences)),
            "tumor_median": float(np.median(tumor_values[keep])) if keep.any() else np.nan,
            "normal_median": float(np.median(normal_values[keep])) if keep.any() else np.nan,
            "median_delta_tumor_minus_normal": float(np.median(differences)) if len(differences) else np.nan,
            "wilcoxon_statistic": float(statistic) if np.isfinite(statistic) else np.nan,
            "p_value": float(p_value),
        })
    gene_result = pd.DataFrame(rows)
    gene_result["BH_FDR_within_comparison"] = bh(gene_result["p_value"])

    score_rows = []
    score_subset = score_df[score_df["contrast"].eq("TCGA_primary_vs_TCGA_solid_normal")].copy()
    score_subset = score_subset.merge(metadata[["sample_id", "patient_id"]], on="sample_id", how="left")
    tumor_scores = score_subset[score_subset["group"].eq("TCGA_CRC_primary_tumor")].drop_duplicates("patient_id").set_index("patient_id")
    normal_scores = score_subset[score_subset["group"].eq("TCGA_CRC_solid_normal")].drop_duplicates("patient_id").set_index("patient_id")
    paired_score_ids = sorted(set(tumor_scores.index.dropna()) & set(normal_scores.index.dropna()))
    for score_name in ["PPAR_nuclear_receptor_score", "inflammatory_RELA_STAT3_score", "DINP_axis_9_gene_score"]:
        tumor_values = tumor_scores.loc[paired_score_ids, score_name].to_numpy(float)
        normal_values = normal_scores.loc[paired_score_ids, score_name].to_numpy(float)
        keep = np.isfinite(tumor_values) & np.isfinite(normal_values)
        differences = tumor_values[keep] - normal_values[keep]
        if len(differences) == 0 or np.allclose(differences, 0):
            statistic, p_value = np.nan, 1.0
        else:
            statistic, p_value = wilcoxon(differences, alternative="two-sided", method="auto")
        score_rows.append({
            "comparison": "TCGA_paired_primary_vs_solid_normal",
            "score": score_name,
            "paired_n": int(len(differences)),
            "tumor_median": float(np.median(tumor_values[keep])) if keep.any() else np.nan,
            "normal_median": float(np.median(normal_values[keep])) if keep.any() else np.nan,
            "median_delta_tumor_minus_normal": float(np.median(differences)) if len(differences) else np.nan,
            "wilcoxon_statistic": float(statistic) if np.isfinite(statistic) else np.nan,
            "p_value": float(p_value),
        })
    score_result = pd.DataFrame(score_rows)
    score_result["BH_FDR_within_comparison"] = bh(score_result["p_value"])
    return gene_result, score_result


def make_figure(gene_stats: pd.DataFrame, score_stats: pd.DataFrame, output_path: Path) -> None:
    contrasts = list(gene_stats["contrast"].drop_duplicates())
    pivot = gene_stats.pivot(index="gene", columns="contrast", values="median_delta_tumor_minus_normal").reindex(GENES)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), gridspec_kw={"height_ratios": [1, 1.05]})
    ax = axes[0, 0]
    im = ax.imshow(pivot.to_numpy(float), aspect="auto", cmap="coolwarm", vmin=-3, vmax=3)
    ax.set_yticks(np.arange(len(GENES)), GENES)
    ax.set_xticks(np.arange(len(contrasts)), [x.replace("_", "\n") for x in contrasts])
    ax.set_title("Tumor − normal median expression")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.iloc[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Xena expression units")

    ax = axes[0, 1]
    fdr = gene_stats.pivot(index="gene", columns="contrast", values="BH_FDR_within_contrast").reindex(GENES)
    image = -np.log10(np.clip(fdr.to_numpy(float), 1e-300, 1))
    im2 = ax.imshow(image, aspect="auto", cmap="viridis", vmin=0, vmax=max(5, np.nanpercentile(image, 95)))
    ax.set_yticks(np.arange(len(GENES)), GENES)
    ax.set_xticks(np.arange(len(contrasts)), [x.replace("_", "\n") for x in contrasts])
    ax.set_title("Gene-level evidence: −log10(BH-FDR)")
    fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04, label="−log10 FDR")

    ax = axes[1, 0]
    score_names = ["PPAR_nuclear_receptor_score", "inflammatory_RELA_STAT3_score", "DINP_axis_9_gene_score"]
    labels = ["PPAR/NR", "RELA+STAT3", "9-gene axis"]
    xpos = np.arange(len(score_names))
    width = 0.36
    for k, contrast in enumerate(contrasts):
        subset = score_stats[score_stats["contrast"].eq(contrast)].set_index("score").reindex(score_names)
        ax.bar(xpos + (k - (len(contrasts) - 1) / 2) * width, subset["median_delta_tumor_minus_normal"], width=width, label=contrast.replace("_", " "))
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(xpos, labels)
    ax.set_ylabel("Tumor − normal median z-score")
    ax.set_title("Targeted pathway scores")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    show = gene_stats[gene_stats["contrast"].eq(contrasts[0])].copy().sort_values("median_delta_tumor_minus_normal")
    ax.barh(show["gene"], show["median_delta_tumor_minus_normal"], color=["#4C78A8" if x >= 0 else "#E45756" for x in show["median_delta_tumor_minus_normal"]])
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Tumor − TCGA solid normal median")
    ax.set_title("Primary TCGA contrast")
    fig.suptitle("Phase 2F-A: DINP-axis CRC disease-state validation", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    metadata = load_sample_metadata()
    included = metadata.loc[metadata["include"], "sample_id"].tolist()
    expression = load_expression(included)
    metadata = metadata[metadata["sample_id"].isin(expression.index)].copy()

    stats_tcga = compare_gene_sets(expression, metadata, "TCGA_primary_vs_TCGA_solid_normal", "TCGA_CRC_solid_normal")
    stats_gtex = compare_gene_sets(expression, metadata, "TCGA_primary_vs_GTEx_colon_normal", "GTEx_colon_normal")
    gene_stats = pd.concat([stats_tcga, stats_gtex], ignore_index=True)

    scores_tcga = score_samples(expression, metadata, "TCGA_primary_vs_TCGA_solid_normal", "TCGA_CRC_solid_normal")
    scores_gtex = score_samples(expression, metadata, "TCGA_primary_vs_GTEx_colon_normal", "GTEx_colon_normal")
    score_samples_df = pd.concat([scores_tcga, scores_gtex], ignore_index=True)
    score_stats = compare_scores(score_samples_df)
    paired_gene_stats, paired_score_stats = paired_tcga_stats(expression, metadata, score_samples_df)

    direction_pivot = gene_stats.pivot(index="gene", columns="contrast", values="median_delta_tumor_minus_normal").reindex(GENES)
    direction_rows = []
    for gene, row in direction_pivot.iterrows():
        signs = np.sign(row.to_numpy(float))
        direction_rows.append({
            "gene": gene,
            "TCGA_normal_direction": "up" if signs[0] > 0 else "down" if signs[0] < 0 else "zero_or_missing",
            "GTEx_normal_direction": "up" if signs[1] > 0 else "down" if signs[1] < 0 else "zero_or_missing",
            "direction_concordant": bool(signs[0] == signs[1] and signs[0] != 0),
        })
    direction_df = pd.DataFrame(direction_rows)

    # A simple COAD vs READ descriptive check is included, but it is not a
    # primary gate and is not interpreted as a replication contrast.
    tumor_meta = metadata[metadata["group"].eq("TCGA_CRC_primary_tumor")].copy()
    tumor_meta["tcga_subtype"] = tumor_meta["primary disease or tissue"].map({"Colon Adenocarcinoma": "COAD", "Rectum Adenocarcinoma": "READ"})
    subtype_rows = []
    for gene in GENES:
        coad = expression.loc[expression.index.intersection(tumor_meta.loc[tumor_meta["tcga_subtype"].eq("COAD"), "sample_id"]), gene].dropna().to_numpy(float)
        read = expression.loc[expression.index.intersection(tumor_meta.loc[tumor_meta["tcga_subtype"].eq("READ"), "sample_id"]), gene].dropna().to_numpy(float)
        u, p = mannwhitneyu(coad, read, alternative="two-sided")
        subtype_rows.append({"comparison": "TCGA_COAD_vs_READ_primary_tumor", "gene": gene, "COAD_n": len(coad), "READ_n": len(read), "COAD_median": float(np.median(coad)), "READ_median": float(np.median(read)), "median_delta_COAD_minus_READ": float(np.median(coad) - np.median(read)), "p_value": float(p), "mann_whitney_U": float(u)})
    subtype_stats = pd.DataFrame(subtype_rows)
    subtype_stats["BH_FDR_within_comparison"] = bh(subtype_stats["p_value"])

    metadata.to_csv(OUTPUT / "mcop_phase2f_bulk_sample_manifest.csv", index=False)
    gene_stats.to_csv(OUTPUT / "mcop_phase2f_bulk_gene_stats.csv", index=False)
    score_samples_df.to_csv(OUTPUT / "mcop_phase2f_bulk_pathway_scores_by_sample.csv", index=False)
    score_stats.to_csv(OUTPUT / "mcop_phase2f_bulk_pathway_score_stats.csv", index=False)
    paired_gene_stats.to_csv(OUTPUT / "mcop_phase2f_tcga_paired_gene_stats.csv", index=False)
    paired_score_stats.to_csv(OUTPUT / "mcop_phase2f_tcga_paired_pathway_score_stats.csv", index=False)
    direction_df.to_csv(OUTPUT / "mcop_phase2f_bulk_direction_concordance.csv", index=False)
    subtype_stats.to_csv(OUTPUT / "mcop_phase2f_bulk_coad_read_descriptive.csv", index=False)
    make_figure(gene_stats, score_stats, OUTPUT / "mcop_phase2f_figure_bulk_expression.png")

    group_counts = metadata.loc[metadata["include"], "group"].value_counts().to_dict()
    primary_gene = gene_stats[gene_stats["contrast"].eq("TCGA_primary_vs_TCGA_solid_normal")].sort_values("BH_FDR_within_contrast")
    primary_score = score_stats[score_stats["contrast"].eq("TCGA_primary_vs_TCGA_solid_normal")].sort_values("BH_FDR_within_contrast")
    nr_primary = primary_score[primary_score["score"].eq("PPAR_nuclear_receptor_score")].iloc[0]
    nr_paired = paired_score_stats[paired_score_stats["score"].eq("PPAR_nuclear_receptor_score")].iloc[0]
    gtex_nr = score_stats[(score_stats["contrast"].eq("TCGA_primary_vs_GTEx_colon_normal")) & score_stats["score"].eq("PPAR_nuclear_receptor_score")].iloc[0]
    concordant_n = int(direction_df["direction_concordant"].sum())
    report = [
        "# MCOP–CRC Phase 2F-A：CRC PPAR/nuclear-receptor disease-state validation",
        "",
        "## 判定",
        "",
        "本轮只验证 Phase 2E 冻结的 DINP-axis molecular program，不把 MCOP 当作 CTD-specific 分子发现物，也不进行泛 GO/PPI 扩张。",
        "",
        f"- PPAR/nuclear-receptor score 在 TCGA 内部对照中下降：tumor median={float(nr_primary['tumor_median']):.3g}，normal median={float(nr_primary['normal_median']):.3g}，delta={float(nr_primary['median_delta_tumor_minus_normal']):.3g}，P={float(nr_primary['p_value']):.3g}，BH-FDR={float(nr_primary['BH_FDR_within_contrast']):.3g}；病人级配对分析同方向，delta={float(nr_paired['median_delta_tumor_minus_normal']):.3g}，P={float(nr_paired['p_value']):.3g}。",
        f"- 但换成 GTEx colon normal 后 PPAR/NR score 方向变为上升：delta={float(gtex_nr['median_delta_tumor_minus_normal']):.3g}，P={float(gtex_nr['p_value']):.3g}，BH-FDR={float(gtex_nr['BH_FDR_within_contrast']):.3g}。9 个基因中只有 {concordant_n}/9 个基因在两个 normal reference 下方向一致。",
        "- 因此 Phase 2F-A 当前判定是 **reference-dependent / provisional，不升级为机制验证通过**。TCGA 内部配对结果支持疾病状态相关性，但 TCGA–GTEx 对照不支持一个稳定的统一方向。",
        f"- 样本量：TCGA primary tumor={int(group_counts.get('TCGA_CRC_primary_tumor', 0))}；TCGA CRC solid normal={int(group_counts.get('TCGA_CRC_solid_normal', 0))}；GTEx transverse/sigmoid colon normal={int(group_counts.get('GTEx_colon_normal', 0))}。",
        "- 这一步可以回答“该 program 是否在 CRC tumor state 中呈现疾病相关表达状态”，但不能回答“MCOP 是否导致该 program 改变”。",
        "",
        "## 1. 数据与对照",
        "",
        "| Contrast | Tumor | Normal | 用途 |",
        "|---|---:|---:|---|",
        f"| TCGA primary CRC | {int(group_counts.get('TCGA_CRC_primary_tumor', 0))} | {int(group_counts.get('TCGA_CRC_solid_normal', 0))} | primary tumor vs within-TCGA solid normal |",
        f"| TCGA vs GTEx | {int(group_counts.get('TCGA_CRC_primary_tumor', 0))} | {int(group_counts.get('GTEx_colon_normal', 0))} | external normal-tissue reference |",
        "",
        "Expression source: UCSC Toil Xena hub `TcgaTargetGtex_rsem_gene_tpm`; phenotype source: `TcgaTargetGTEX_phenotype.txt`. Values are analyzed on the Xena-delivered scale; no cross-platform re-normalization is applied beyond within-contrast z-scoring for pathway scores.",
        "",
        "## 2. Frozen gene sets",
        "",
        f"- PPAR/nuclear-receptor core: `{'; '.join(NR_GENES)}`",
        f"- inflammatory complement: `{'; '.join(INFLAMMATORY_GENES)}`",
        f"- combined DINP-axis score: all `{'; '.join(GENES)}`",
        "",
        "## 3. How to read the outputs",
        "",
        "`mcop_phase2f_bulk_gene_stats.csv` reports gene-level medians, tumor-minus-normal shifts, Mann–Whitney P values and BH-FDR separately for the TCGA-normal and GTEx-normal contrasts.",
        "",
        "`mcop_phase2f_bulk_pathway_score_stats.csv` reports sample-level mean z-scores for the PPAR/NR, RELA+STAT3 and 9-gene programs. The score is a disease-state readout, not a causal mediation score.",
        "",
        "`mcop_phase2f_tcga_paired_gene_stats.csv` and `mcop_phase2f_tcga_paired_pathway_score_stats.csv` use available patient-matched TCGA primary tumor/solid-normal pairs; these are the most internally comparable public checks in this run.",
        "",
        "## 4. Single-cell status",
        "",
        "CELLxGENE Census was selected for Phase 2F-B because its versioned metadata and primary-data filter are appropriate for cross-dataset single-cell queries. The current Windows runtime could not install `cellxgene-census` because the required TileDB-SOMA package had no compatible wheel and fell back to a local CMake build without a compiler. Therefore **no single-cell result is claimed in this commit**. The next single-cell run must pin the Census release and include `is_primary_data == True` before comparing malignant epithelial, myeloid and stromal compartments.",
        "",
        "## 5. Interpretation boundary",
        "",
        "The TCGA-internal paired shift supports CRC disease-state relevance, but the reference discordance prevents a clean mechanism upgrade. The next step should resolve normal-reference and tissue-composition effects before any TCGA mechanistic narrative is written. This result still does not prove exposure-specific direction or replace an exposure-linked tissue/perturbation experiment.",
        "",
        "## Files",
        "",
        "- `mcop_phase2f_bulk_gene_stats.csv`",
        "- `mcop_phase2f_bulk_pathway_score_stats.csv`",
        "- `mcop_phase2f_bulk_pathway_scores_by_sample.csv`",
        "- `mcop_phase2f_bulk_sample_manifest.csv`",
        "- `mcop_phase2f_bulk_coad_read_descriptive.csv`",
        "- `mcop_phase2f_tcga_paired_gene_stats.csv`",
        "- `mcop_phase2f_tcga_paired_pathway_score_stats.csv`",
        "- `mcop_phase2f_bulk_direction_concordance.csv`",
        "- `mcop_phase2f_figure_bulk_expression.png`",
        "",
        "## Reproducibility",
        "",
        f"- Run UTC: `{datetime.now(timezone.utc).isoformat()}`",
        f"- Script: `{Path(__file__).relative_to(ROOT)}`",
        f"- Xena hub: `{XENA_HUB}`",
        f"- Expression dataset: `{EXPRESSION_DATASET}`",
        f"- Phenotype dataset: `{PHENOTYPE_DATASET}`",
        f"- Gene-set hash: `{sha256_text('|'.join(GENES))}`",
        f"- Python: `{platform.python_version()}`",
        "",
        "**Phase 2F-A status: bulk tumor-normal analysis completed; result is reference-dependent/provisional. Single-cell validation remains a separately labeled Phase 2F-B task.**",
    ]
    (OUTPUT / "mcop_phase2f_crc_expression_validation_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "MCOP-CRC Phase 2F-A targeted CRC expression validation",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(ROOT)),
        "xena_hub": XENA_HUB,
        "expression_dataset": EXPRESSION_DATASET,
        "phenotype_dataset": PHENOTYPE_DATASET,
        "genes": GENES,
        "PPAR_nuclear_receptor_genes": NR_GENES,
        "inflammatory_genes": INFLAMMATORY_GENES,
        "group_counts": {str(k): int(v) for k, v in group_counts.items()},
        "single_cell_status": "pending_cellxgene_census_runtime",
        "single_cell_blocker": "cellxgene-census installation requires TileDB-SOMA build unavailable in the current Windows runtime",
        "analysis_boundary": "disease-state expression validation, not exposure causality or mediation",
    }
    (OUTPUT / "mcop_phase2f_crc_expression_validation_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
