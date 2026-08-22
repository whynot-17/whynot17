"""Phase 2E molecular validation for the DINP/MiNP -> MCOP CRC axis.

This is intentionally a bridge-validation analysis, not a new chemical screen.
It reuses the frozen CTD human interaction export and GeneCards Disorders CRC
export from Phase 1, and separates:

* MiNP (monoisononylphthalate), the Phase 1 DINP-axis molecular discovery;
* MCOP (mono(carboxy-isooctyl)phthalate), the NHANES urinary biomarker;
* diisononyl phthalate, the parent DINP comparator; and
* MBzP, the historical positive/negative human-validation comparator.

The script audits evidence concentration and co-treatment dependence, reports
the exact Phase 1 primary enrichment rows, and runs deliberately exploratory
ORA against local Hallmark and Reactome GMT files using U_core as the
background.  No pathway result is treated as confirmatory when the input list
is too small or dominated by co-treatment evidence.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "work" / "environmental_toxicology_crc_phase1" / "data"
OUTPUT = ROOT / "outputs"

CTD_CHEM = DATA / "CTD_chemicals.tsv.gz"
CTD_IXN = DATA / "CTD_chem_gene_ixns.tsv.gz"
GENECARDS = DATA / "genecards_disorders_crc.csv"
CLASSIFICATION = OUTPUT / "environmental_toxicology_crc_phase1_chemical_classification.csv"
PHASE1_RANKED = OUTPUT / "environmental_toxicology_crc_phase1_ranked_core.csv"
HALLMARK_GMT = ROOT / "work" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
REACTOME_GMT = ROOT / "work" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

CHEMICALS = {
    "MiNP": {"id": "C471400", "name": "monoisononylphthalate", "role": "DINP-axis molecular discovery"},
    "MCOP": {"id": "C573544", "name": "mono(carboxy-isooctyl)phthalate", "role": "NHANES urinary biomarker"},
    "DINP_parent": {"id": "C012125", "name": "diisononyl phthalate", "role": "parent DINP comparator"},
    "MBzP": {"id": "C103325", "name": "mono-benzyl phthalate", "role": "historical human-validation comparator"},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_ctd_chemicals() -> pd.DataFrame:
    names = [
        "ChemicalName", "ChemicalID", "CasRN", "PubChemCID", "PubChemSID",
        "DTXSID", "InChIKey", "Definition", "ParentIDs", "TreeNumbers",
        "ParentTreeNumbers", "MESHSynonyms", "CTDCuratedSynonyms",
    ]
    out = pd.read_csv(CTD_CHEM, sep="\t", comment="#", header=None, names=names, dtype=str)
    out["ChemicalID"] = out["ChemicalID"].astype(str).str.replace(r"^MESH:", "", regex=True)
    return out


def read_ctd_interactions() -> pd.DataFrame:
    names = [
        "ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
        "GeneForms", "Organism", "OrganismID", "Interaction", "InteractionActions",
        "PubMedIDs",
    ]
    return pd.read_csv(CTD_IXN, sep="\t", comment="#", header=None, names=names, dtype=str)


def read_genecards() -> pd.DataFrame:
    gc = pd.read_csv(GENECARDS, dtype=str)
    gc["GeneCards_Rank"] = pd.to_numeric(gc["GeneCards_Rank"], errors="coerce")
    gc["RelevanceScore"] = pd.to_numeric(gc["RelevanceScore"], errors="coerce")
    gc["GeneSymbol"] = gc["GeneSymbol"].astype(str).str.strip().str.upper()
    return gc.sort_values("GeneCards_Rank").drop_duplicates("GeneSymbol")


def split_pmids(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {x for x in re.split(r"[|;,\s]+", str(value)) if x and x.isdigit()}


def flag_cotreatment(row: pd.Series) -> bool:
    actions = str(row.get("InteractionActions", "") or "").lower()
    interaction = str(row.get("Interaction", "") or "").lower()
    return "cotreatment" in actions or "co-treated" in interaction or "co-treated" in actions


def flag_single_chemical(row: pd.Series) -> bool:
    if bool(row["is_co_treatment"]):
        return False
    name = str(row["ChemicalName"] or "").strip().lower()
    interaction = str(row["Interaction"] or "").strip().lower()
    return bool(name and interaction.startswith(name))


def infer_bridge_class(gene_symbol: str, gene_type: object, gene_forms: object) -> str:
    gt = str(gene_type or "").lower()
    sym = str(gene_symbol or "").upper()
    forms = str(gene_forms or "").lower()
    if "protein coding" in gt or "protein" in forms:
        return "protein_coding"
    if "mir" in sym or "microrna" in gt:
        return "miRNA"
    if "antisense" in gt or sym.endswith("-AS1") or sym.endswith("AS1"):
        return "lncRNA_or_antisense"
    if "rna" in gt or "gene" in forms or "mrna" in forms:
        return "RNA_or_unresolved"
    return "unresolved"


def enrichment_stats(genes: set[str], crc_genes: set[str], universe: set[str]) -> dict[str, float | int]:
    genes = set(genes) & set(universe)
    crc_genes = set(crc_genes) & set(universe)
    n_universe = len(universe)
    n_crc = len(crc_genes)
    n_genes = len(genes)
    overlap = len(genes & crc_genes)
    table = [
        [overlap, n_genes - overlap],
        [n_crc - overlap, n_universe - n_crc - n_genes + overlap],
    ]
    odds_ratio, fisher_p = fisher_exact(table, alternative="greater")
    hypergeom_p = float(hypergeom.sf(overlap - 1, n_universe, n_crc, n_genes)) if overlap else 1.0
    return {
        "n_universe": n_universe,
        "n_crc_genes": n_crc,
        "n_interacting_genes": n_genes,
        "crc_overlap": overlap,
        "odds_ratio": float(odds_ratio) if np.isfinite(odds_ratio) else float("inf"),
        "fisher_p": float(fisher_p),
        "hypergeom_p": hypergeom_p,
        "overlap_genes": ";".join(sorted(genes & crc_genes)),
    }


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


def read_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                sets[fields[0]] = {x.upper() for x in fields[2:] if x}
    return sets


def pathway_ora(query_name: str, query_genes: set[str], library_name: str, library: dict[str, set[str]], universe: set[str]) -> pd.DataFrame:
    query = set(query_genes) & set(universe)
    rows: list[dict[str, object]] = []
    for term, term_genes in library.items():
        term_in_universe = term_genes & universe
        if len(term_in_universe) < 5:
            continue
        overlap = query & term_in_universe
        p = float(hypergeom.sf(len(overlap) - 1, len(universe), len(term_in_universe), len(query))) if overlap else 1.0
        rows.append({
            "query": query_name,
            "library": library_name,
            "term": term,
            "query_n_mappable": len(query),
            "term_n_in_universe": len(term_in_universe),
            "overlap_n": len(overlap),
            "overlap_genes": ";".join(sorted(overlap)),
            "p_value": p,
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["library_bh_fdr"] = bh(out["p_value"])
    return out


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    chemicals = read_ctd_chemicals()
    interactions = read_ctd_interactions()
    genecards = read_genecards()
    classification = pd.read_csv(CLASSIFICATION, dtype=str)
    classification["is_core"] = classification["is_core"].astype(str).str.lower().eq("true")

    target_ids = {v["id"] for v in CHEMICALS.values()}
    target_names = {v["name"] for v in CHEMICALS.values()}
    target_vocab = chemicals[chemicals["ChemicalID"].isin(target_ids)].copy()
    if set(target_vocab["ChemicalID"]) != target_ids:
        raise ValueError("One or more frozen target chemical IDs are missing from the CTD vocabulary")

    human = interactions[interactions["OrganismID"].eq("9606")].copy()
    human["GeneSymbol"] = human["GeneSymbol"].astype(str).str.strip().str.upper()
    human["is_co_treatment"] = human.apply(flag_cotreatment, axis=1)
    human["is_single_chemical"] = human.apply(flag_single_chemical, axis=1)
    human["evidence_context"] = np.select(
        [human["is_co_treatment"], human["is_single_chemical"]],
        ["co-treatment", "single-chemical statement"],
        default="indirect/other",
    )

    core_ids = set(classification.loc[classification["is_core"], "ChemicalID"].astype(str))
    core_human = human[human["ChemicalID"].isin(core_ids)]
    core_universe = set(core_human["GeneSymbol"].dropna())
    if len(core_universe) != 22786:
        raise ValueError(f"U_core mismatch: expected 22786 genes, found {len(core_universe)}")

    # Phase 1 primary set: GeneCards Disorders, k=1000, U_core. The export has
    # 706 rows but only the 585 genes in U_core enter the primary test.
    crc_genes = set(genecards.loc[genecards["GeneSymbol"].isin(core_universe), "GeneSymbol"])
    if len(crc_genes) != 585:
        raise ValueError(f"GeneCards primary set mismatch: expected 585 genes in U_core, found {len(crc_genes)}")

    phase1 = pd.read_csv(PHASE1_RANKED)
    phase1_primary = phase1[
        phase1["scope"].eq("GeneCards_Disorders")
        & phase1["gene_cards_k"].eq(1000)
        & phase1["background"].eq("U_core")
        & phase1["ChemicalID"].isin(target_ids)
    ].copy()
    if len(phase1_primary) != len(CHEMICALS):
        raise ValueError("Could not recover exactly one Phase 1 primary row per target chemical")

    bridge_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    target_gene_sets: dict[str, set[str]] = {}

    for label, meta in CHEMICALS.items():
        subset = human[human["ChemicalID"].eq(meta["id"])].copy()
        unique_pairs = subset.drop_duplicates(["ChemicalID", "GeneID"])
        genes = set(unique_pairs["GeneSymbol"].dropna())
        target_gene_sets[label] = genes
        bridge = genes & crc_genes
        no_cotreatment = set(unique_pairs.loc[~unique_pairs["is_co_treatment"], "GeneSymbol"].dropna())
        phase_row = phase1_primary[phase1_primary["ChemicalID"].eq(meta["id"])].iloc[0]
        no_cot_stats = enrichment_stats(no_cotreatment, crc_genes, core_universe)
        for gene in sorted(bridge):
            gene_rows = subset[subset["GeneSymbol"].eq(gene)]
            pmids = set().union(*(split_pmids(x) for x in gene_rows["PubMedIDs"]))
            co_pmids = set().union(*(split_pmids(x) for x in gene_rows.loc[gene_rows["is_co_treatment"], "PubMedIDs"]))
            direct_pmids = set().union(*(split_pmids(x) for x in gene_rows.loc[gene_rows["is_single_chemical"], "PubMedIDs"]))
            gc_row = genecards[genecards["GeneSymbol"].eq(gene)]
            gc_row = gc_row.iloc[0] if not gc_row.empty else None
            gene_type = gc_row["GeneType"] if gc_row is not None else ""
            rank = gc_row["GeneCards_Rank"] if gc_row is not None else np.nan
            relevance = gc_row["RelevanceScore"] if gc_row is not None else np.nan
            bridge_rows.append({
                "chemical_label": label,
                "ChemicalID": meta["id"],
                "ChemicalName": meta["name"],
                "role": meta["role"],
                "GeneID": ";".join(sorted(set(gene_rows["GeneID"].dropna().astype(str)))),
                "GeneSymbol": gene,
                "GeneCards_rank": rank,
                "GeneCards_relevance_score": relevance,
                "GeneCards_gene_type": gene_type,
                "bridge_class": infer_bridge_class(gene, gene_type, ";".join(gene_rows["GeneForms"].dropna().astype(str))),
                "n_raw_interaction_rows": len(gene_rows),
                "n_unique_pmids": len(pmids),
                "n_co_treatment_rows": int(gene_rows["is_co_treatment"].sum()),
                "n_co_treatment_pmids": len(co_pmids),
                "n_single_chemical_rows": int(gene_rows["is_single_chemical"].sum()),
                "n_single_chemical_pmids": len(direct_pmids),
                "has_single_chemical_evidence": bool(len(direct_pmids)),
                "PMIDs": ";".join(sorted(pmids)),
                "co_treatment_PMIDs": ";".join(sorted(co_pmids)),
                "single_chemical_PMIDs": ";".join(sorted(direct_pmids)),
            })

        summary_rows.append({
            "chemical_label": label,
            "ChemicalID": meta["id"],
            "ChemicalName": meta["name"],
            "role": meta["role"],
            "n_raw_interaction_rows": len(subset),
            "n_unique_chemical_gene_pairs": len(unique_pairs),
            "n_unique_human_genes": len(genes),
            "n_unique_pmids": len(set().union(*(split_pmids(x) for x in subset["PubMedIDs"]))),
            "n_co_treatment_rows": int(subset["is_co_treatment"].sum()),
            "n_co_treatment_genes": subset.loc[subset["is_co_treatment"], "GeneSymbol"].nunique(),
            "n_co_treatment_pmids": len(set().union(*(split_pmids(x) for x in subset.loc[subset["is_co_treatment"], "PubMedIDs"]))),
            "n_single_chemical_rows": int(subset["is_single_chemical"].sum()),
            "n_single_chemical_genes": subset.loc[subset["is_single_chemical"], "GeneSymbol"].nunique(),
            "n_crc_overlap": len(bridge),
            "overlap_genes": ";".join(sorted(bridge)),
            "primary_OR": phase_row["odds_ratio"],
            "primary_Fisher_P": phase_row["fisher_p"],
            "primary_Hypergeom_P": phase_row["hypergeom_p"],
            "primary_BH_FDR": phase_row["bh_fdr"],
            "primary_rank_weighted_overlap": phase_row["rank_weighted_overlap"],
            "primary_top_overlap_genes": phase_row["top_overlap_genes"],
            "no_cotreatment_n_genes": no_cot_stats["n_interacting_genes"],
            "no_cotreatment_crc_overlap": no_cot_stats["crc_overlap"],
            "no_cotreatment_OR": no_cot_stats["odds_ratio"],
            "no_cotreatment_Fisher_P_unadjusted": no_cot_stats["fisher_p"],
            "no_cotreatment_overlap_genes": no_cot_stats["overlap_genes"],
        })

    bridge_df = pd.DataFrame(bridge_rows).sort_values(["chemical_label", "GeneCards_rank"])
    summary_df = pd.DataFrame(summary_rows)
    evidence_df = human[human["ChemicalID"].isin(target_ids)].copy()
    evidence_df = evidence_df[[
        "ChemicalName", "ChemicalID", "GeneSymbol", "GeneID", "GeneForms",
        "Interaction", "InteractionActions", "PubMedIDs", "is_co_treatment",
        "is_single_chemical", "evidence_context",
    ]].sort_values(["ChemicalID", "GeneSymbol", "PubMedIDs"])

    # Exploratory local ORA. The query sets are intentionally explicit so the
    # report cannot silently switch to all genes or all CTD chemicals.
    hallmark = read_gmt(HALLMARK_GMT)
    reactome = read_gmt(REACTOME_GMT)
    pathway_frames = []
    for label in ["MiNP", "MCOP", "DINP_parent", "MBzP"]:
        pathway_frames.append(pathway_ora(label + "_all_human_interacting_genes", target_gene_sets[label], "Hallmark_2026.1", hallmark, core_universe))
        pathway_frames.append(pathway_ora(label + "_all_human_interacting_genes", target_gene_sets[label], "Reactome_2026.1", reactome, core_universe))
        subset = evidence_df[(evidence_df["ChemicalID"].eq(CHEMICALS[label]["id"])) & (~evidence_df["is_co_treatment"])]
        pathway_frames.append(pathway_ora(label + "_no_cotreatment_genes", set(subset["GeneSymbol"]), "Hallmark_2026.1", hallmark, core_universe))
        pathway_frames.append(pathway_ora(label + "_no_cotreatment_genes", set(subset["GeneSymbol"]), "Reactome_2026.1", reactome, core_universe))
    pathway_df = pd.concat([x for x in pathway_frames if not x.empty], ignore_index=True) if any(not x.empty for x in pathway_frames) else pd.DataFrame()
    if not pathway_df.empty:
        pathway_df["global_bh_fdr"] = bh(pathway_df["p_value"])
        pathway_df = pathway_df.sort_values(["global_bh_fdr", "p_value", "overlap_n"])

    # A compact publication-oriented evidence figure.
    plot_df = summary_df.copy().set_index("chemical_label").loc[["MiNP", "MCOP", "DINP_parent", "MBzP"]].reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1, 1.3]})
    x = np.arange(len(plot_df))
    axes[0].bar(x - 0.18, plot_df["n_unique_human_genes"], width=0.36, label="CTD human genes", color="#4C78A8")
    axes[0].bar(x + 0.18, plot_df["n_crc_overlap"], width=0.36, label="CRC overlap", color="#E45756")
    axes[0].set_xticks(x, plot_df["chemical_label"])
    axes[0].set_ylabel("Gene count")
    axes[0].set_title("CTD bridge size")
    axes[0].legend(frameon=False, fontsize=8)
    for i, row in plot_df.iterrows():
        axes[0].text(i + 0.18, row["n_crc_overlap"] + 0.6, str(int(row["n_crc_overlap"])), ha="center", fontsize=9)

    evidence_cols = ["n_unique_pmids", "n_co_treatment_genes", "n_single_chemical_genes"]
    labels = ["unique PMIDs", "co-treatment genes", "single-chemical genes"]
    colors = ["#72B7B2", "#F58518", "#54A24B"]
    bottom = np.zeros(len(plot_df))
    for col, label, color in zip(evidence_cols, labels, colors):
        vals = plot_df[col].astype(float).to_numpy()
        axes[1].bar(x, vals, bottom=bottom, label=label, color=color)
        bottom += vals
    axes[1].set_xticks(x, plot_df["chemical_label"])
    axes[1].set_ylabel("Count")
    axes[1].set_title("Evidence composition")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].text(0.02, 0.96, "MiNP CRC bridge is mostly co-treatment evidence;\nMCOP has no matching CTD bridge strength.", transform=axes[1].transAxes, va="top", fontsize=8.5)
    fig.suptitle("Phase 2E DINP-axis molecular bridge audit", fontsize=13)
    fig.tight_layout()
    figure_path = OUTPUT / "mcop_phase2e_figure_bridge_evidence.png"
    fig.savefig(figure_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    bridge_path = OUTPUT / "mcop_phase2e_molecular_bridge_qc.csv"
    evidence_path = OUTPUT / "mcop_phase2e_molecular_evidence_long.csv"
    summary_path = OUTPUT / "mcop_phase2e_molecular_candidate_summary.csv"
    pathway_path = OUTPUT / "mcop_phase2e_pathway_ora.csv"
    bridge_df.to_csv(bridge_path, index=False)
    evidence_df.to_csv(evidence_path, index=False)
    summary_df.to_csv(summary_path, index=False)
    pathway_df.to_csv(pathway_path, index=False)

    minp = summary_df.loc[summary_df["chemical_label"].eq("MiNP")].iloc[0]
    mcop = summary_df.loc[summary_df["chemical_label"].eq("MCOP")].iloc[0]
    minp_bridge = bridge_df[bridge_df["chemical_label"].eq("MiNP")]
    minp_overlap_co = int((minp_bridge["n_co_treatment_rows"] > 0).sum())
    minp_overlap_direct = int(minp_bridge["has_single_chemical_evidence"].sum())
    pathway_sig = pathway_df[(pathway_df["global_bh_fdr"] < 0.05) & (pathway_df["overlap_n"] >= 2)] if not pathway_df.empty else pathway_df
    minp_pathway_sig = pathway_sig[pathway_sig["query"].str.startswith("MiNP_")] if not pathway_sig.empty else pathway_sig
    pathway_sig_distinct_terms = int(pathway_sig["term"].nunique()) if not pathway_sig.empty else 0
    minp_pathway_sig_distinct_terms = int(minp_pathway_sig["term"].nunique()) if not minp_pathway_sig.empty else 0

    report_lines = [
        "# MCOP–CRC Phase 2E：DINP-axis molecular validation",
        "",
        "## 结论先行",
        "",
        "本轮不是重新筛化学物，也不是把 MCOP 事后包装成 CTD 的直接发现物。冻结的逻辑是：",
        "",
        "> **CTD 发现 DINP/MiNP 轴的分子桥接 → NHANES 以 MCOP 作为 DINP 轴尿液 biomarker 出现 CRC association → 机制验证先审计桥接证据的可重复性。**",
        "",
        f"- **MiNP**：{int(minp['n_unique_human_genes'])} 个 human interacting genes，{int(minp['n_crc_overlap'])} 个 GeneCards Disorders CRC overlap；Phase 1 primary OR={float(minp['primary_OR']):.3g}，BH-FDR={float(minp['primary_BH_FDR']):.3g}。",
        f"- **MCOP**：{int(mcop['n_unique_human_genes'])} 个 human interacting genes，{int(mcop['n_crc_overlap'])} 个 overlap；Phase 1 primary BH-FDR={float(mcop['primary_BH_FDR']):.3g}，因此不能写成 CTD 已直接发现 MCOP–CRC 分子桥。",
        f"- MiNP 的 {int(minp['n_crc_overlap'])} 个 CRC overlap gene 中，{minp_overlap_co} 个有 co-treatment evidence，{minp_overlap_direct} 个有 single-chemical statement；共处理依赖性是本桥接目前最大的限制。",
        f"- 本轮本地 Hallmark/Reactome ORA 使用 **U_core={len(core_universe):,} 个 CTD 可检测 human genes** 作背景。全局 BH-FDR<0.05 且 overlap≥2 的 exploratory rows：**{len(pathway_sig)}**（{pathway_sig_distinct_terms} 个不同 term；MiNP query 为 {len(minp_pathway_sig)} rows/{minp_pathway_sig_distinct_terms} terms）；Reactome 层级 term 高度冗余，因此不把这个计数当作独立机制数量。",
        "",
        "## 1. 分析边界与输入冻结",
        "",
        "| 元素 | 冻结定义 |",
        "|---|---|",
        "| CTD interaction | Homo sapiens；统计按 unique ChemicalID × GeneID；PMID 只用于证据审计 |",
        "| CRC gene set | GeneCards Disorders-scoped export；与 U_core 相交后 585 genes；Phase 1 primary row 为 k=1000 |",
        "| U_core | 最终 core environmental toxicants 的所有 CTD human interacting genes；本轮复核为 22,786 genes |",
        "| 分子比较 | MiNP、MCOP、DINP parent、MBzP；不重新筛选候选 |",
        "| pathway | local MSigDB Hallmark 2026.1 与 Reactome 2026.1 GMT；ORA exploratory；不使用全基因组背景 |",
        "",
        "## 2. MiNP 与 MCOP 的桥接对照",
        "",
        "| 化学物 | 角色 | CTD human genes | CRC overlap | unique PMIDs | co-treatment genes | single-chemical genes | Phase 1 BH-FDR |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in summary_df.set_index("chemical_label").loc[["MiNP", "MCOP", "DINP_parent", "MBzP"]].reset_index().iterrows():
        report_lines.append(
            f"| {row['chemical_label']} | {row['role']} | {int(row['n_unique_human_genes'])} | {int(row['n_crc_overlap'])} | {int(row['n_unique_pmids'])} | {int(row['n_co_treatment_genes'])} | {int(row['n_single_chemical_genes'])} | {float(row['primary_BH_FDR']):.3g} |"
        )
    report_lines += [
        "",
        "### MiNP CRC overlap genes",
        "",
        f"`{minp['overlap_genes']}`",
        "",
        "逐基因证据见 `mcop_phase2e_molecular_bridge_qc.csv`。关键限制是：BAX、CASP8、MIR141、CDKN2A 的 MiNP CTD 记录来自同一个 phthalate co-treatment study/PMID，而 PPARG 才有单化学物的 binding/activity 记录。因此，这个桥接支持“DINP/MiNP 轴存在 CRC-relevant molecular evidence”，但目前不支持“MiNP 单独通过 5 个 CRC genes 直接驱动该桥接”。",
        "",
        "## 3. Co-treatment sensitivity",
        "",
        "去掉所有带 co-treatment 标记的 CTD rows 后，重新使用同一个 U_core 与 GeneCards CRC background 计算的是**未进行跨化学物重排的敏感性 P 值**，不是新的 confirmatory FDR。",
        "",
        "| 化学物 | no-co-treatment genes | CRC overlap | OR | Fisher P（未调整） | overlap genes |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for _, row in summary_df.iterrows():
        report_lines.append(
            f"| {row['chemical_label']} | {int(row['no_cotreatment_n_genes'])} | {int(row['no_cotreatment_crc_overlap'])} | {float(row['no_cotreatment_OR']):.3g} | {float(row['no_cotreatment_Fisher_P_unadjusted']):.3g} | {row['no_cotreatment_overlap_genes'] or '—'} |"
        )
    report_lines += [
        "",
        "MiNP 的结果在去除共处理证据后明显变弱，说明 Phase 1 的强富集不能被解释为一个完全独立、纯 MiNP 单化学物的 CRC gene signature。它仍可作为轴级别机制假设，但需要单化学物实验或独立暴露组学来确认。",
        "",
        "## 4. Pathway ORA",
        "",
        "本轮 ORA 只用于定位可检验的方向，不用于把 24 个 CTD genes 变成机制定论。对于 MiNP，4/5 CRC overlap gene 受共处理证据影响，且 gene list 很短；因此 pathway 结果必须按 exploratory 处理。",
        "",
        "结果文件：`mcop_phase2e_pathway_ora.csv`。所有可测试通路（包括 0 overlap）均纳入 BH；全局 BH-FDR 作为跨 query 与 Hallmark+Reactome 的保守参考。",
        "",
        "MiNP 的领先方向（只作候选机制）：",
        "",
        "最稳定的方向是核受体/PPAR 相关转录调控；这个方向在去掉 co-treatment 后仍由 PPARA、PPARG、PPARD、NR1H2/NR1H3、NR1I2/NR1I3 等基因驱动，但它是 **DINP-axis level** 的候选机制，不等于 CRC-specific pathway 已被验证。Reactome 的 SUMOylation、脂质代谢和凋亡 terms 之间存在明显基因集重叠。",
    ]
    if not minp_pathway_sig.empty:
        for _, row in minp_pathway_sig.sort_values(["global_bh_fdr", "p_value"]).head(8).iterrows():
            report_lines.append(f"- `{row['term']}`：query={row['query']}，overlap={int(row['overlap_n'])}，genes={row['overlap_genes']}，global BH-FDR={float(row['global_bh_fdr']):.3g}。")
    else:
        report_lines.append("- 经完整 BH 后没有达到 global BH-FDR<0.05 且 overlap≥2 的 MiNP pathway row；不预设机制方向。")
    report_lines += [
        "",
        "## 5. 当前机制判定",
        "",
        "### 支持",
        "",
        "- MiNP 的 Phase 1 富集在冻结的 CTD×GeneCards 流程中可复核；",
        "- MiNP overlap 中包含 PPARG 的 single-chemical binding/activity evidence；",
        "- MCOP 人群信号可以合理描述为 DINP-axis biomarker signal，而不是 MCOP 已被 CTD 分子桥接直接发现；",
        "- 机制验证的下一步应围绕 PPAR/核受体、氧化应激/凋亡和表观遗传调控提出可检验假设，而不是继续扩大无先验通路清单。",
        "",
        "### 限制",
        "",
        "- MiNP CRC overlap 的主要证据集中在共处理研究；",
        "- MCOP 本身 CTD overlap 很弱（2 genes，FDR≈0.284），不能作为 MCOP-specific molecular proof；",
        "- ORA 的输入列表很短且来自数据库交集，显著性不等于暴露因果；",
        "- NHANES 仍是横断面人群发现，WHI 尚未获得或分析真实 biospecimen。",
        "",
        "## 文件",
        "",
        "- `mcop_phase2e_molecular_candidate_summary.csv`：四个化学物的桥接与证据汇总",
        "- `mcop_phase2e_molecular_bridge_qc.csv`：CRC overlap gene-level evidence QC",
        "- `mcop_phase2e_molecular_evidence_long.csv`：CTD interaction 原始证据长表",
        "- `mcop_phase2e_pathway_ora.csv`：custom-background Hallmark/Reactome exploratory ORA",
        "- `mcop_phase2e_figure_bridge_evidence.png`：桥接大小与证据构成图",
        "",
        "## 可复现性",
        "",
        f"- 运行时间（UTC）：{datetime.now(timezone.utc).isoformat()}",
        f"- 脚本：`{Path(__file__).relative_to(ROOT)}`",
        f"- CTD chemicals SHA256：`{sha256(CTD_CHEM)}`",
        f"- CTD interactions SHA256：`{sha256(CTD_IXN)}`",
        f"- GeneCards Disorders SHA256：`{sha256(GENECARDS)}`",
        f"- Phase 1 ranked output SHA256：`{sha256(PHASE1_RANKED)}`",
        "",
        "**最终判断：Phase 2E 通过的是“继续做 DINP-axis molecular validation”的门，不是“MCOP 已有一对一 CTD 机制证明”的门。**",
    ]
    report_path = OUTPUT / "mcop_phase2e_molecular_validation_report.md"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    manifest = {
        "analysis": "MCOP-CRC Phase 2E DINP-axis molecular validation",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(ROOT)),
        "chemical_targets": CHEMICALS,
        "ctd_species": "Homo sapiens",
        "primary_gene_cards_scope": "Disorders",
        "primary_gene_cards_crc_n_in_U_core": len(crc_genes),
        "U_core_n": len(core_universe),
        "pathway_libraries": {"Hallmark_2026.1": str(HALLMARK_GMT.relative_to(ROOT)), "Reactome_2026.1": str(REACTOME_GMT.relative_to(ROOT))},
        "input_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in [CTD_CHEM, CTD_IXN, GENECARDS, CLASSIFICATION, PHASE1_RANKED, HALLMARK_GMT, REACTOME_GMT]},
        "interpretation": "MiNP axis bridge supports a hypothesis but is partly co-treatment driven; MCOP remains a human biomarker, not a CTD-specific molecular discovery.",
    }
    (OUTPUT / "mcop_phase2e_molecular_validation_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
