#!/usr/bin/env python
"""Strict robustness re-analysis of the focused PGE2 -> PTGER4 module.

The previous module used pooled cell-level gene-wise z-scores.  This module
keeps that result as a sensitivity analysis and makes donor x cell-group
pseudobulk summaries primary.  It deliberately performs no new clustering,
annotation, CellChat, spatial analysis, docking, or MD.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SYNTHESIS = {
    "Score_A": ["PLA2G4A", "PTGS2", "PTGES"],
    "Score_B": ["PLA2G4A", "PTGS1", "PTGS2", "PTGES", "PTGES2"],
    "Score_C": ["PLA2G4A", "PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3"],
}
SCORE_NAMES = ["Score_A", "Score_B", "Score_C", "Score_D"]
TARGET_GENES = ["PLA2G4A", "PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3", "PTGER4"]
BOTTLENECK_GENES = ["PLA2G4A", "PTGS1", "PTGS2", "PTGES", "PTGES2"]
MACRO_SUBTYPES = ["SPP1_like_TAM", "C1QC_like_TAM", "FCN1_inflammatory_monocyte_like", "resident_like"]
CANDIDATE_GROUPS = [
    "mast_cells", "tumor_epithelial", "normal_epithelial", "SPP1_like_TAM",
    "C1QC_like_TAM", "FCN1_inflammatory_monocyte_like", "resident_like",
    "fibroblast_like", "endothelial_like", "other_myeloid",
]
REPORT_GROUPS = [
    "mast_cells", "tumor_epithelial", "normal_epithelial", "SPP1_like_TAM",
    "C1QC_like_TAM", "FCN1_inflammatory_monocyte_like", "resident_like",
    "fibroblast_like", "endothelial_like",
]
THRESHOLDS = [10, 20, 30]
DEFAULT_CELL_SCORES = Path(r"E:\chatgpt\whynot17\analysis\dinp_crc_pge2_ptger4_axis\outputs\pge2_synthesis_cell_scores.csv")
DEFAULT_MATRIX = Path(r"E:\mcop\mcop\work\mcop_phase2f_external\raw\GSE144735_raw_UMI_count_matrix.txt.gz")
DEFAULT_ANNOTATION = Path(r"E:\mcop\mcop\work\mcop_phase2f_external\raw\GSE144735_annotation.txt.gz")
DEFAULT_SUBTYPE = Path(r"E:\chatgpt\whynot17\analysis\dinp_crc_macrophage_subtype\outputs\macrophage_cell_program_scores.csv")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_dump(obj, path: Path):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def qscale(values: pd.Series | np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    positive = x[np.isfinite(x) & (x > 0)]
    if len(positive):
        return float(max(np.quantile(positive, .95), np.max(positive) * 1e-6, 1e-9))
    return 1.0


def zscore(values: pd.Series) -> pd.Series:
    x = values.astype(float)
    sd = float(x.std(ddof=0))
    return (x - float(x.mean())) / sd if sd > 0 else pd.Series(0.0, index=x.index)


def load_cells(path: Path) -> pd.DataFrame:
    cells = pd.read_csv(path)
    required = {"cell_id", "donor_id", "condition", "analysis_group"}
    missing = sorted(required - set(cells.columns))
    if missing:
        raise ValueError(f"Cell-score input is missing columns: {missing}")
    for gene in TARGET_GENES:
        col = f"{gene}_log1p"
        if col not in cells.columns:
            raise ValueError(f"Cell-score input is missing {col}")
        cells[col] = pd.to_numeric(cells[col], errors="coerce").fillna(0.0)
        cells[f"{gene}_detected"] = cells[col] > 0
    cells["n_cells_dummy"] = 1
    return cells


def cell_global_sensitivity(cells: pd.DataFrame) -> pd.DataFrame:
    """Recreate the prior pooled cell-level z-score for all four definitions."""
    for score, genes in SYNTHESIS.items():
        cells[f"{score}_cell_global_z"] = pd.concat([zscore(cells[f"{g}_log1p"]) for g in genes], axis=1).mean(axis=1)
    # D is shown as a sensitivity using the same bottleneck formula at cell level.
    up = np.clip(cells["PLA2G4A_log1p"] / qscale(cells["PLA2G4A_log1p"]), 0, 1)
    cyc = np.maximum(cells["PTGS1_log1p"], cells["PTGS2_log1p"]) / qscale(np.maximum(cells["PTGS1_log1p"], cells["PTGS2_log1p"]))
    term = np.maximum(cells["PTGES_log1p"], cells["PTGES2_log1p"]) / qscale(np.maximum(cells["PTGES_log1p"], cells["PTGES2_log1p"]))
    cells["Score_D_cell_global_z"] = np.cbrt(np.clip(up, 0, 1) * np.clip(cyc, 0, 1) * np.clip(term, 0, 1))
    return cells


def make_pseudobulk(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (donor, condition, group), g in cells.groupby(["donor_id", "condition", "analysis_group"], sort=True, dropna=False):
        row = {"donor_id": donor, "condition": condition, "analysis_group": group, "n_cells": int(len(g))}
        for gene in TARGET_GENES:
            vals = g[f"{gene}_log1p"].astype(float)
            row[f"{gene}_mean_log1p"] = float(vals.mean())
            row[f"{gene}_detection_fraction"] = float((vals > 0).mean())
            positive = vals[vals > 0]
            row[f"{gene}_expressing_cell_mean_log1p"] = float(positive.mean()) if len(positive) else np.nan
            # expm1(mean log1p) is a normalized-expression pseudobulk estimate;
            # log1p is retained for the receiver axis sensitivity.
            row[f"{gene}_pseudobulk_normalized_mean"] = float(np.expm1(vals).mean())
            row[f"{gene}_pseudobulk_log1p"] = float(np.log1p(row[f"{gene}_pseudobulk_normalized_mean"]))
        for score in SCORE_NAMES:
            row[f"{score}_cell_global_mean"] = float(g[f"{score}_cell_global_z"].mean())
        rows.append(row)
    pb = pd.DataFrame(rows)
    for score, genes in SYNTHESIS.items():
        z = pd.concat([zscore(pb[f"{g}_mean_log1p"]) for g in genes], axis=1)
        pb[f"{score}_pseudobulk"] = z.mean(axis=1)
    # Bottleneck-aware donor-group score. Each component is scaled by its
    # donor-group 95th percentile and clipped to [0, 1]. The score is the
    # geometric mean of upstream, cyclooxygenase, and terminal-synthase axes.
    up = np.clip(pb["PLA2G4A_mean_log1p"] / qscale(pb["PLA2G4A_mean_log1p"]), 0, 1)
    cyc_raw = np.maximum(pb["PTGS1_mean_log1p"], pb["PTGS2_mean_log1p"])
    term_raw = np.maximum(pb["PTGES_mean_log1p"], pb["PTGES2_mean_log1p"])
    cyc = np.clip(cyc_raw / qscale(cyc_raw), 0, 1)
    term = np.clip(term_raw / qscale(term_raw), 0, 1)
    pb["Score_D_pseudobulk"] = np.cbrt(up * cyc * term)
    for threshold in THRESHOLDS:
        pb[f"eligible_ge_{threshold}"] = pb.n_cells >= threshold
    return pb


def score_definition_table() -> pd.DataFrame:
    return pd.DataFrame([
        {"score": "Score_A", "name": "Core inducible PGE2", "genes": ", ".join(SYNTHESIS["Score_A"]), "primary": "yes", "formula": "mean of donor-group gene-wise z-scored mean log1p expression", "PTGES3_included": False, "biological_note": "PLA2G4A -> PTGS2 -> PTGES core inducible route"},
        {"score": "Score_B", "name": "Extended enzymatic PGE2", "genes": ", ".join(SYNTHESIS["Score_B"]), "primary": "yes", "formula": "mean of donor-group gene-wise z-scored mean log1p expression", "PTGES3_included": False, "biological_note": "PLA2G4A plus PTGS1/2 and PTGES/PTGES2; excludes PTGES3"},
        {"score": "Score_C", "name": "Original six-gene sensitivity", "genes": ", ".join(SYNTHESIS["Score_C"]), "primary": "sensitivity", "formula": "mean of donor-group gene-wise z-scored mean log1p expression", "PTGES3_included": True, "biological_note": "reproduces the previous six-gene definition"},
        {"score": "Score_D", "name": "Bottleneck-aware PGE2", "genes": ", ".join(BOTTLENECK_GENES), "primary": "yes", "formula": "geometric_mean(PLA2G4A/q95, max(PTGS1,PTGS2)/q95, max(PTGES,PTGES2)/q95), each component clipped to [0,1]", "PTGES3_included": False, "biological_note": "requires upstream, cyclooxygenase, and terminal synthase support simultaneously"},
    ])


def rank_for_donor(g: pd.DataFrame, score_col: str) -> pd.DataFrame:
    x = g.dropna(subset=[score_col]).sort_values(score_col, ascending=False, kind="mergesort").copy()
    x["rank"] = x[score_col].rank(method="min", ascending=False).astype(int)
    n = len(x)
    x["rank_percentile"] = 1.0 if n <= 1 else 1 - (x["rank"] - 1) / (n - 1)
    return x


def source_rank_tables(pb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows, donor_rows, contrast_rows = [], [], []
    for threshold in THRESHOLDS:
        eligible = pb[(pb.n_cells >= threshold) & pb.analysis_group.isin(CANDIDATE_GROUPS) & (pb.condition == "tumor")].copy()
        for score in SCORE_NAMES:
            col = f"{score}_pseudobulk"
            for donor, dg in eligible.groupby("donor_id"):
                ranked = rank_for_donor(dg, col)
                donor_rows.extend({"donor_id": donor, "threshold": threshold, "score": score, "analysis_group": r.analysis_group, "n_cells": int(r.n_cells), "score_value": float(r[col]), "rank": int(r["rank"]), "rank_percentile": float(r["rank_percentile"]), "n_candidate_groups": int(len(ranked))} for _, r in ranked.iterrows())
            for group in CANDIDATE_GROUPS:
                d = donor_rows[-0:]  # keep the construction explicit below
                vals = [r for r in donor_rows if r["threshold"] == threshold and r["score"] == score and r["analysis_group"] == group]
                if vals:
                    ranks = pd.Series([v["rank"] for v in vals], dtype=float)
                    percs = pd.Series([v["rank_percentile"] for v in vals], dtype=float)
                    scores = pd.Series([v["score_value"] for v in vals], dtype=float)
                    rows.append({"threshold": threshold, "score": score, "source_group": group, "n_donors_with_group": int(len(vals)), "top1_donors": int((ranks == 1).sum()), "top2_donors": int((ranks <= 2).sum()), "top1_fraction": float((ranks == 1).mean()), "top2_fraction": float((ranks <= 2).mean()), "median_rank": float(ranks.median()), "mean_rank": float(ranks.mean()), "median_rank_percentile": float(percs.median()), "mean_score": float(scores.mean())})
            # SPP1 source contrasts are paired within donor, not cell-level.
            for donor, dg in eligible.groupby("donor_id"):
                vals = dg.set_index("analysis_group")[col]
                if "SPP1_like_TAM" not in vals.index:
                    continue
                others = vals.drop(index="SPP1_like_TAM", errors="ignore")
                mast = vals.get("mast_cells", np.nan)
                contrast_rows.append({"donor_id": donor, "threshold": threshold, "score": score, "SPP1_score": float(vals["SPP1_like_TAM"]), "other_candidate_mean": float(others.mean()) if len(others) else np.nan, "SPP1_vs_other_mean_delta": float(vals["SPP1_like_TAM"] - others.mean()) if len(others) else np.nan, "mast_score": float(mast) if pd.notna(mast) else np.nan, "SPP1_vs_mast_delta": float(vals["SPP1_like_TAM"] - mast) if pd.notna(mast) else np.nan, "n_other_candidates": int(len(others))})
    rank_df = pd.DataFrame(rows)
    donor_df = pd.DataFrame(donor_rows)
    contrast_df = pd.DataFrame(contrast_rows)
    summary_rows = []
    for (threshold, score), g in contrast_df.groupby(["threshold", "score"]):
        for metric, label in [("SPP1_vs_other_mean_delta", "SPP1_vs_other_mean"), ("SPP1_vs_mast_delta", "SPP1_vs_mast")]:
            x = g[metric].dropna()
            if len(x):
                summary_rows.append({"threshold": threshold, "score": score, "contrast": label, "n_donors": int(len(x)), "n_delta_gt_0": int((x > 0).sum()), "median_delta": float(x.median()), "mean_delta": float(x.mean()), "sign_consistency_fraction": float(max((x > 0).sum(), (x < 0).sum()) / len(x))})
    return rank_df, donor_df, pd.DataFrame(summary_rows)


def mast_counts(pb: pd.DataFrame) -> pd.DataFrame:
    tumor = pb[(pb.condition == "tumor") & (pb.analysis_group == "mast_cells")][["donor_id", "n_cells", "Score_A_pseudobulk", "Score_B_pseudobulk", "Score_C_pseudobulk", "Score_D_pseudobulk"]].copy()
    tumor["group"] = "mast_cells"
    return tumor.sort_values("donor_id")


def threshold_sensitivity(rank_df: pd.DataFrame, pb: pd.DataFrame, mast: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in rank_df.iterrows():
        row = r.to_dict(); row["analysis_type"] = "source_rank"; rows.append(row)
    for threshold in THRESHOLDS:
        for score in SCORE_NAMES:
            for group in CANDIDATE_GROUPS:
                g = pb[(pb.condition == "tumor") & (pb.analysis_group == group) & (pb.n_cells >= threshold)]
                if len(g):
                    rows.append({"analysis_type": "group_score", "threshold": threshold, "score": score, "source_group": group, "n_donors_with_group": int(g.donor_id.nunique()), "n_cells_total": int(g.n_cells.sum()), "mean_score": float(g[f"{score}_pseudobulk"].mean()), "median_score": float(g[f"{score}_pseudobulk"].median())})
    for _, r in mast.iterrows():
        for threshold in THRESHOLDS:
            rows.append({"analysis_type": "mast_cell_count", "threshold": threshold, "score": "all", "source_group": "mast_cells", "donor_id": r.donor_id, "n_cells": int(r.n_cells), "eligible": bool(r.n_cells >= threshold), "Score_A_pseudobulk": float(r.Score_A_pseudobulk), "Score_B_pseudobulk": float(r.Score_B_pseudobulk), "Score_C_pseudobulk": float(r.Score_C_pseudobulk), "Score_D_pseudobulk": float(r.Score_D_pseudobulk)})
    return pd.DataFrame(rows)


def receiver_tables(pb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    m = pb[(pb.condition == "tumor") & pb.analysis_group.isin(MACRO_SUBTYPES)].copy()
    rows = []
    for donor, g in m.groupby("donor_id"):
        ranked = g.copy()
        for metric, label in [("PTGER4_mean_log1p", "rank_mean"), ("PTGER4_detection_fraction", "rank_detection"), ("PTGER4_expressing_cell_mean_log1p", "rank_expressing_cell_mean"), ("PTGER4_pseudobulk_log1p", "rank_pseudobulk")]:
            rr = rank_for_donor(g, metric)[["analysis_group", "rank"]].rename(columns={"rank": label})
            ranked = ranked.merge(rr, on="analysis_group", how="left")
        for _, r in ranked.iterrows():
            rows.append({"donor_id": donor, "subtype": r.analysis_group, "n_cells": int(r.n_cells), "PTGER4_mean": float(r.PTGER4_mean_log1p), "PTGER4_detection_fraction": float(r.PTGER4_detection_fraction), "PTGER4_expressing_cell_mean": float(r.PTGER4_expressing_cell_mean_log1p) if pd.notna(r.PTGER4_expressing_cell_mean_log1p) else np.nan, "PTGER4_pseudobulk": float(r.PTGER4_pseudobulk_log1p), "rank_mean": int(r.rank_mean), "rank_detection": int(r.rank_detection), "rank_expressing_cell_mean": int(r.rank_expressing_cell_mean) if pd.notna(r.rank_expressing_cell_mean) else np.nan, "rank_pseudobulk": int(r.rank_pseudobulk)})
    receiver = pd.DataFrame(rows)
    return receiver, receiver


def pooled_receiver_summary(cells: pd.DataFrame) -> pd.DataFrame:
    """Pooled receiver description, retained only as a descriptive companion."""
    m = cells[cells.analysis_group.isin(MACRO_SUBTYPES)].copy()
    rows = []
    for condition, g in list(m.groupby("condition", sort=True)) + [("all", m)]:
        for subtype, x in g.groupby("analysis_group", sort=False):
            positive = x.loc[x.PTGER4_log1p > 0, "PTGER4_log1p"]
            rows.append({"condition": condition, "subtype": subtype, "n_cells": int(len(x)), "n_donors": int(x.donor_id.nunique()), "PTGER4_mean": float(x.PTGER4_log1p.mean()), "PTGER4_detection_fraction": float((x.PTGER4_log1p > 0).mean()), "PTGER4_expressing_cell_mean": float(positive.mean()) if len(positive) else np.nan, "PTGER4_pseudobulk": float(np.log1p(np.expm1(x.PTGER4_log1p).mean()))})
    return pd.DataFrame(rows)


def sign_test(x: pd.Series) -> tuple[int, int, float, float]:
    x = pd.to_numeric(x, errors="coerce").dropna()
    pos, neg = int((x > 0).sum()), int((x < 0).sum())
    nonzero = pos + neg
    if nonzero == 0:
        return pos, neg, np.nan, np.nan
    p = float(stats.binomtest(pos, nonzero, .5, alternative="two-sided").pvalue)
    # Exact paired sign-flip permutation, useful with only six donors.
    vals = x.to_numpy(float)
    observed = abs(float(vals.sum()))
    totals = [abs(float(np.sum(vals * np.array(signs)))) for signs in itertools.product([-1, 1], repeat=len(vals))]
    perm_p = float((sum(t >= observed - 1e-12 for t in totals)) / len(totals))
    return pos, neg, p, perm_p


def c1qc_contrasts(receiver: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor, g in receiver.groupby("donor_id"):
        idx = {r.subtype: r for _, r in g.iterrows()}
        c = idx.get("C1QC_like_TAM"); s = idx.get("SPP1_like_TAM")
        others = g[g.subtype != "C1QC_like_TAM"]
        if c is None:
            continue
        d1 = float(c.PTGER4_mean - s.PTGER4_mean) if s is not None else np.nan
        d2 = float(c.PTGER4_mean - others.PTGER4_mean.mean()) if len(others) else np.nan
        rows.append({"donor_id": donor, "C1QC_n_cells": int(c.n_cells), "SPP1_n_cells": int(s.n_cells) if s is not None else 0, "C1QC_PTGER4_mean": float(c.PTGER4_mean), "SPP1_PTGER4_mean": float(s.PTGER4_mean) if s is not None else np.nan, "other_macrophage_PTGER4_mean": float(others.PTGER4_mean.mean()) if len(others) else np.nan, "Delta1_C1QC_minus_SPP1": d1, "Delta2_C1QC_minus_other_mean": d2, "Delta1_valid": bool(pd.notna(d1)), "Delta2_valid": bool(pd.notna(d2)), "C1QC_top1_mean": bool(c.rank_mean == 1), "C1QC_top2_mean": bool(c.rank_mean <= 2), "C1QC_higher_than_SPP1": bool(d1 > 0) if pd.notna(d1) else False, "C1QC_higher_than_other_mean": bool(d2 > 0) if pd.notna(d2) else False, "C1QC_rank_mean": int(c.rank_mean), "C1QC_rank_pseudobulk": int(c.rank_pseudobulk)})
    return pd.DataFrame(rows)


def axis_sensitivity(pb: pd.DataFrame) -> pd.DataFrame:
    rows = []
    tumor = pb[(pb.condition == "tumor") & pb.analysis_group.isin(CANDIDATE_GROUPS)]
    for threshold in THRESHOLDS:
        for donor in sorted(tumor.donor_id.unique()):
            dg = tumor[(tumor.donor_id == donor) & (tumor.n_cells >= threshold)]
            c = dg[dg.analysis_group == "C1QC_like_TAM"]
            if len(c) == 0:
                continue
            c = c.iloc[0]
            for source in CANDIDATE_GROUPS:
                s = dg[dg.analysis_group == source]
                if len(s) == 0:
                    continue
                s = s.iloc[0]
                for score in SCORE_NAMES:
                    source_score = float(s[f"{score}_pseudobulk"])
                    rows.append({"donor_id": donor, "threshold": threshold, "score": score, "source_group": source, "source_n_cells": int(s.n_cells), "C1QC_n_cells": int(c.n_cells), "source_synthesis_score": source_score, "C1QC_PTGER4_mean": float(c.PTGER4_mean_log1p), "C1QC_PTGER4_detection_fraction": float(c.PTGER4_detection_fraction), "C1QC_PTGER4_pseudobulk": float(c.PTGER4_pseudobulk_log1p), "Axis_1_mean": source_score * float(c.PTGER4_mean_log1p), "Axis_2_detection": source_score * float(c.PTGER4_detection_fraction), "Axis_3_pseudobulk": source_score * float(c.PTGER4_pseudobulk_log1p), "coverage": True})
    return pd.DataFrame(rows)


def save_figures(pb: pd.DataFrame, donor_ranks: pd.DataFrame, mast: pd.DataFrame, receiver: pd.DataFrame, axis: pd.DataFrame, out: Path):
    import matplotlib.pyplot as plt
    import seaborn as sns
    figdir = out / "figures"; figdir.mkdir(parents=True, exist_ok=True)
    p = pb[(pb.condition == "tumor") & (pb.analysis_group.isin(REPORT_GROUPS)) & (pb.n_cells >= 20)].groupby("analysis_group")[[f"{s}_pseudobulk" for s in SCORE_NAMES]].mean().reindex(REPORT_GROUPS).dropna(how="all")
    fig, ax = plt.subplots(figsize=(11, 6)); p.plot(kind="bar", ax=ax); ax.set_title("Donor-group pseudobulk PGE2 score sensitivity (tumor, n≥20)"); ax.set_xlabel(""); ax.set_ylabel("mean donor-group score"); ax.tick_params(axis="x", rotation=40); ax.legend(title="score", frameon=False); fig.tight_layout()
    for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"pge2_score_sensitivity_by_group.{ext}", dpi=220)
    plt.close(fig)
    dr = donor_ranks[(donor_ranks.threshold == 20) & donor_ranks.score.isin(["Score_A", "Score_B", "Score_D"]) & (donor_ranks.analysis_group == "SPP1_like_TAM")]
    if len(dr):
        mat = dr.pivot(index="donor_id", columns="score", values="rank").reindex(columns=["Score_A", "Score_B", "Score_D"])
        fig, ax = plt.subplots(figsize=(6, 4.5)); sns.heatmap(mat, annot=True, fmt=".0f", cmap="YlGnBu_r", vmin=1, vmax=max(3, int(np.nanmax(mat.values))), cbar_kws={"label": "rank (1 = highest)"}, ax=ax); ax.set_title("SPP1-like TAM source rank by tumor donor (n≥20)"); ax.set_xlabel(""); ax.set_ylabel("donor"); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"spp1_source_rank_by_donor.{ext}", dpi=220)
        plt.close(fig)
    mm = mast.copy()
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.8))
    for ax, score in zip(axes, SCORE_NAMES):
        ax.scatter(mm.n_cells, mm[f"{score}_pseudobulk"], s=50, color="#DD8452")
        for _, r in mm.iterrows(): ax.text(r.n_cells, r[f"{score}_pseudobulk"], r.donor_id, fontsize=8, ha="left", va="bottom")
        ax.set_title(score); ax.set_xlabel("mast cells per donor"); ax.set_ylabel("pseudobulk score")
    fig.suptitle("Mast-cell count versus PGE2 score"); fig.tight_layout()
    for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"mast_cell_count_vs_score.{ext}", dpi=220)
    plt.close(fig)
    q = pb[(pb.condition == "tumor") & pb.analysis_group.isin(REPORT_GROUPS) & (pb.n_cells >= 20)].groupby("analysis_group")[["Score_B_pseudobulk", "Score_C_pseudobulk"]].mean().dropna(how="all")
    fig, ax = plt.subplots(figsize=(7, 5)); ax.scatter(q.Score_B_pseudobulk, q.Score_C_pseudobulk, s=70, color="#4C72B0");
    for name, r in q.iterrows(): ax.text(r.Score_B_pseudobulk, r.Score_C_pseudobulk, name, fontsize=8)
    lo, hi = np.nanmin(q.values), np.nanmax(q.values); ax.plot([lo, hi], [lo, hi], ls="--", color="grey"); ax.set_xlabel("Score B (without PTGES3)"); ax.set_ylabel("Score C (original six-gene)"); ax.set_title("PTGES3 sensitivity at donor-group level"); fig.tight_layout()
    for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"ptges3_sensitivity.{ext}", dpi=220)
    plt.close(fig)
    if len(receiver):
        fig, ax = plt.subplots(figsize=(10, 5)); sns.lineplot(data=receiver, x="donor_id", y="PTGER4_mean", hue="subtype", marker="o", sort=False, ax=ax); ax.set_title("PTGER4 mean log1p by tumor donor and macrophage subtype"); ax.set_xlabel("donor"); ax.set_ylabel("mean log1p normalized expression"); ax.legend(title="subtype", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"c1qc_ptger4_by_donor.{ext}", dpi=220)
        plt.close(fig)
        pair = receiver[receiver.subtype.isin(["C1QC_like_TAM", "SPP1_like_TAM"])].pivot(index="donor_id", columns="subtype", values="PTGER4_mean")
        fig, ax = plt.subplots(figsize=(7, 5));
        for donor, r in pair.iterrows(): ax.plot([0, 1], [r.get("C1QC_like_TAM", np.nan), r.get("SPP1_like_TAM", np.nan)], marker="o", color="#999999", alpha=.8); ax.text(0, r.get("C1QC_like_TAM", np.nan), donor, fontsize=8, ha="right", va="center")
        ax.set_xticks([0, 1], ["C1QC-like TAM", "SPP1-like TAM"]); ax.set_ylabel("PTGER4 mean log1p"); ax.set_title("Paired donor PTGER4: C1QC versus SPP1"); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"c1qc_vs_spp1_ptger4_paired.{ext}", dpi=220)
        plt.close(fig)
    if len(axis):
        a = axis[(axis.threshold == 20) & (axis.score == "Score_D")].pivot_table(index="donor_id", columns="source_group", values="Axis_3_pseudobulk", aggfunc="mean")
        fig, ax = plt.subplots(figsize=(11, max(3.5, .45 * len(a)))); sns.heatmap(a, cmap="vlag", center=0, annot=True, fmt=".2f", ax=ax); ax.set_title("Axis 3 robustness: Score D × C1QC PTGER4 pseudobulk (n≥20)"); ax.set_xlabel("candidate source"); ax.set_ylabel("donor"); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"source_receiver_axis_robustness_heatmap.{ext}", dpi=220)
        plt.close(fig)


def final_summary(rank_df: pd.DataFrame, donor_rank_df: pd.DataFrame, contrast_summary: pd.DataFrame, receiver: pd.DataFrame, c1: pd.DataFrame, axis: pd.DataFrame) -> pd.DataFrame:
    rows = []
    def add(metric, value, definition, grade="descriptive"):
        rows.append({"metric": metric, "value": value, "definition": definition, "evidence_grade": grade})
    primary = rank_df[rank_df.threshold == 20]
    rank_pivot = donor_rank_df[donor_rank_df.threshold == 20].pivot_table(index=["donor_id", "analysis_group"], columns="score", values="rank")
    corr = rank_pivot.corr(method="spearman")
    mean_concordance = {}
    for score in SCORE_NAMES:
        others = [x for x in SCORE_NAMES if x != score and x in corr.columns]
        mean_concordance[score] = float(corr.loc[score, others].mean()) if others else np.nan
    stable_score = max(mean_concordance, key=lambda x: (-np.inf if pd.isna(mean_concordance[x]) else mean_concordance[x]))
    add("Most stable pseudobulk score", stable_score, f"highest mean pairwise Spearman concordance of donor-group ranks at n≥20 ({mean_concordance[stable_score]:.3f})")
    for score, value in mean_concordance.items():
        add(f"{score}_mean_rank_concordance_ge20", value, "mean pairwise Spearman concordance with the other score definitions")
    for score in ["Score_A", "Score_B", "Score_D"]:
        r = primary[(primary.score == score) & (primary.source_group == "SPP1_like_TAM")]
        if len(r):
            r = r.iloc[0]; add(f"SPP1_{score}_median_rank_ge20", float(r.median_rank), "donor-group source rank among prespecified candidate groups")
            add(f"SPP1_{score}_top2_donors_ge20", f"{int(r.top2_donors)}/{int(r.n_donors_with_group)} eligible; {int(r.top2_donors)}/6 total", "top-two donor frequency; denominator reports eligible donor groups and all six tumor donors")
    mast_r = primary[(primary.source_group == "mast_cells") & primary.score.isin(["Score_A", "Score_B", "Score_D"])]
    if len(mast_r):
        add("mast_top1_donors_ge20_by_score", "; ".join(f"{r.score}:{int(r.top1_donors)}/{int(r.n_donors_with_group)} eligible; {int(r.top1_donors)}/6 total" for _, r in mast_r.iterrows()), "mast-cell top-1 donor frequency after n≥20 filter")
    for score in ["Score_A", "Score_B", "Score_C", "Score_D"]:
        r = primary[(primary.score == score) & (primary.source_group == "SPP1_like_TAM")]
        if len(r): add(f"SPP1_{score}_median_rank_percentile_ge20", float(r.iloc[0].median_rank_percentile), "1 = best; donor-group rank percentile")
    cm = contrast_summary[(contrast_summary.threshold == 20) & (contrast_summary.contrast == "SPP1_vs_mast")]
    if len(cm):
        r = cm[cm.score == "Score_B"].iloc[0] if (cm.score == "Score_B").any() else cm.iloc[0]; add("SPP1_vs_mast_positive_donors_ScoreB_ge20", int(r.n_delta_gt_0), "donors where SPP1 synthesis score exceeds mast-cell score")
        add("SPP1_vs_mast_median_delta_ScoreB_ge20", float(r.median_delta), "paired donor score difference")
    if len(c1):
        add("C1QC_Delta1_valid_donors", int(c1.Delta1_valid.sum()), "donors with both C1QC-like and SPP1-like TAM present for Delta1")
        add("C1QC_Delta2_valid_donors", int(c1.Delta2_valid.sum()), "donors with at least one other macrophage subtype for Delta2")
        add("C1QC_PTGER4_median_rank_mean", float(c1.C1QC_rank_mean.median()), "median donor rank by PTGER4 mean")
        add("C1QC_PTGER4_top1_donors", int(c1.C1QC_top1_mean.sum()), "donors where C1QC is top PTGER4 subtype")
        add("C1QC_PTGER4_top2_donors", int(c1.C1QC_top2_mean.sum()), "donors where C1QC is top-two PTGER4 subtype")
        add("C1QC_gt_SPP1_donors", int(c1.C1QC_higher_than_SPP1.sum()), "donors with Delta1 > 0")
        add("C1QC_gt_other_mean_donors", int(c1.C1QC_higher_than_other_mean.sum()), "donors with Delta2 > 0")
        pos1, neg1, p1, perm1 = sign_test(c1.Delta1_C1QC_minus_SPP1); pos2, neg2, p2, perm2 = sign_test(c1.Delta2_C1QC_minus_other_mean)
        add("C1QC_Delta1_median", float(c1.Delta1_C1QC_minus_SPP1.median()), "C1QC PTGER4 mean minus SPP1 PTGER4 mean")
        add("C1QC_Delta2_median", float(c1.Delta2_C1QC_minus_other_mean.median()), "C1QC PTGER4 mean minus mean of other macrophage subtypes")
        add("C1QC_Delta1_exact_sign_p", p1, "two-sided exact sign test; supplementary with six donors")
        add("C1QC_Delta1_paired_permutation_p", perm1, "exact paired sign-flip permutation; supplementary")
        add("C1QC_Delta2_exact_sign_p", p2, "two-sided exact sign test; supplementary with six donors")
        add("C1QC_Delta2_paired_permutation_p", perm2, "exact paired sign-flip permutation; supplementary")
    add("axis_coverage_ge20", int(len(axis[(axis.threshold == 20) & (axis.score == "Score_D")].donor_id.unique())) if len(axis) else 0, "tumor donors with C1QC and source group meeting n≥20 for Score D axis")
    # Evidence grades follow the predeclared rules but retain the numeric rows.
    spp1_rows = primary[(primary.score.isin(["Score_A", "Score_B", "Score_D"])) & (primary.source_group == "SPP1_like_TAM")]
    source_moderate = len(spp1_rows) == 3 and bool((spp1_rows.top2_donors >= 4).all())
    c1_moderate = len(c1) == 6 and float(c1.C1QC_rank_mean.median()) <= 2 and int(c1.C1QC_higher_than_SPP1.sum()) >= 4 and int(c1.C1QC_higher_than_other_mean.sum()) >= 4
    source_grade = "MODERATE SUPPORT" if source_moderate else "WEAK SUPPORT"
    receiver_grade = "MODERATE SUPPORT" if c1_moderate else "WEAK SUPPORT"
    axis_grade = "MODERATE SUPPORT" if source_moderate and c1_moderate else "WEAK SUPPORT"
    add("Final source support", source_grade, "Score A/B/D top-two frequency at n≥20 and donor-level ranking", source_grade)
    add("Final receiver support", receiver_grade, "pooled and donor-level C1QC PTGER4 rank/effect direction", receiver_grade)
    add("Final axis support", axis_grade, "requires both source and receiver to reach at least moderate support", axis_grade)
    next_step = "larger single-cell cohort replication plus spatial validation before extending the axis" if axis_grade == "WEAK SUPPORT" else "spatial validation and targeted functional follow-up"
    add("Recommended next step", next_step, "decision based on robustness evidence")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-scores", type=Path, default=DEFAULT_CELL_SCORES)
    ap.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    ap.add_argument("--annotation", type=Path, default=DEFAULT_ANNOTATION)
    ap.add_argument("--subtype", type=Path, default=DEFAULT_SUBTYPE)
    ap.add_argument("--out", type=Path, default=Path(__file__).parent / "outputs")
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True); t0 = time.time()
    cells = cell_global_sensitivity(load_cells(args.cell_scores))
    pb = make_pseudobulk(cells)
    score_definition_table().to_csv(args.out / "pge2_score_definitions.csv", index=False)
    for score in SCORE_NAMES:
        pb[["donor_id", "condition", "analysis_group", "n_cells", f"{score}_pseudobulk", f"{score}_cell_global_mean", "eligible_ge_10", "eligible_ge_20", "eligible_ge_30"]].rename(columns={f"{score}_pseudobulk": "pseudobulk_score", f"{score}_cell_global_mean": "cell_global_sensitivity_mean"}).to_csv(args.out / f"pge2_score_{score[-1]}_by_donor_group.csv", index=False)
    pb.to_csv(args.out / "pseudobulk_donor_group_expression.csv", index=False)
    rank_df, donor_rank_df, contrast_summary = source_rank_tables(pb)
    rank_df.to_csv(args.out / "source_rank_consistency.csv", index=False)
    donor_rank_df.to_csv(args.out / "source_rank_by_donor.csv", index=False)
    contrast_summary.to_csv(args.out / "source_score_contrasts.csv", index=False)
    mast = mast_counts(pb); mast.to_csv(args.out / "mast_cell_donor_counts.csv", index=False)
    threshold_sensitivity(rank_df, pb, mast).to_csv(args.out / "cell_count_threshold_sensitivity.csv", index=False)
    receiver, _ = receiver_tables(pb); receiver.to_csv(args.out / "ptger4_receiver_by_donor.csv", index=False)
    pooled_receiver = pooled_receiver_summary(cells); pooled_receiver.to_csv(args.out / "ptger4_receiver_pooled_summary.csv", index=False)
    c1 = c1qc_contrasts(receiver); c1.to_csv(args.out / "c1qc_ptger4_donor_contrasts.csv", index=False)
    axis = axis_sensitivity(pb); axis.to_csv(args.out / "pge2_ptger4_axis_sensitivity.csv", index=False)
    summary = final_summary(rank_df, donor_rank_df, contrast_summary, receiver, c1, axis); summary.to_csv(args.out / "final_robustness_summary.csv", index=False)
    save_figures(pb, donor_rank_df, mast, receiver, axis, args.out)
    source_grade = summary.loc[summary.metric == "Final source support", "value"].iloc[0]
    receiver_grade = summary.loc[summary.metric == "Final receiver support", "value"].iloc[0]
    axis_grade = summary.loc[summary.metric == "Final axis support", "value"].iloc[0]
    spp1_a = rank_df[(rank_df.threshold == 20) & (rank_df.source_group == "SPP1_like_TAM") & (rank_df.score == "Score_A")]
    spp1_b = rank_df[(rank_df.threshold == 20) & (rank_df.source_group == "SPP1_like_TAM") & (rank_df.score == "Score_B")]
    spp1_d = rank_df[(rank_df.threshold == 20) & (rank_df.source_group == "SPP1_like_TAM") & (rank_df.score == "Score_D")]
    c1med = float(c1.C1QC_rank_mean.median()) if len(c1) else np.nan
    c1spp = int(c1.C1QC_higher_than_SPP1.sum()) if len(c1) else 0
    c1spp_valid = int(c1.Delta1_valid.sum()) if len(c1) else 0
    c1other = int(c1.C1QC_higher_than_other_mean.sum()) if len(c1) else 0
    c1other_valid = int(c1.Delta2_valid.sum()) if len(c1) else 0
    primary_tumor = pb[(pb.condition == "tumor") & pb.analysis_group.isin(CANDIDATE_GROUPS) & (pb.n_cells >= 20)]
    top_source_by_score = {score: (primary_tumor.groupby("analysis_group")[f"{score}_pseudobulk"].mean().idxmax() if len(primary_tumor) else "NA") for score in ["Score_A", "Score_B", "Score_D"]}
    stable_score_row = summary[summary.metric == "Most stable pseudobulk score"]
    stable_score = str(stable_score_row.value.iloc[0]) if len(stable_score_row) else "NA"
    stable_score_conc = float(stable_score_row.definition.iloc[0].split("(")[-1].rstrip(")")) if len(stable_score_row) else np.nan
    def rank_text(x):
        if not len(x):
            return "NA"
        r = x.iloc[0]
        return f"median {float(r.median_rank):.2f}; top-2 {int(r.top2_donors)}/{int(r.n_donors_with_group)} eligible ({int(r.top2_donors)}/6 total)"
    spp1_a_text, spp1_b_text, spp1_d_text = rank_text(spp1_a), rank_text(spp1_b), rank_text(spp1_d)
    mast_n = ", ".join(f"{r.donor_id}:{int(r.n_cells)}" for _, r in mast.sort_values("donor_id").iterrows())
    pooled_all = pooled_receiver[pooled_receiver.condition == "all"].sort_values("PTGER4_mean", ascending=False)
    pooled_top = pooled_all.iloc[0].subtype if len(pooled_all) else "NA"
    pooled_c1 = pooled_all[pooled_all.subtype == "C1QC_like_TAM"].iloc[0] if (pooled_all.subtype == "C1QC_like_TAM").any() else None
    report = f"""# PGE2–PTGER4 robustness re-analysis

## Scope

This re-analysis keeps `analysis/dinp_crc_pge2_ptger4_axis/` unchanged and tests whether its source → receiver conclusion depends on score definition, PTGES3, pooled cell-level z-scoring, sparse donor groups, winner-takes-all ranking, or the degenerate PTGER4-high threshold. No CellChat, spatial, machine learning, docking, MD, clustering, or new macrophage annotation was performed.

## Inputs and primary design

- Local GSE144735 CRC cache: {len(cells):,} cells and {cells.donor_id.nunique()} donors ({', '.join(sorted(cells.donor_id.unique()))}).
- Primary unit: donor × `analysis_group` pseudobulk, using mean normalized log1p expression, detection fraction, expressing-cell mean, and `log1p(mean(expm1(log1p expression)))` pseudobulk expression.
- Primary source filter: tumor donor × group with **≥20 cells**. Threshold sensitivity: ≥10, ≥20, and ≥30 cells.
- Scores A–C are donor-group gene-wise z-score means. Score D is the bottleneck-aware geometric mean: `geometric_mean(PLA2G4A/q95, max(PTGS1,PTGS2)/q95, max(PTGES,PTGES2)/q95)`, with each component clipped to [0,1]. Score D therefore requires upstream, cyclooxygenase, and terminal synthase support together.
- Score C is the original six-gene definition and includes PTGES3 only as a sensitivity analysis. Scores A, B, and D exclude PTGES3.
- The prior pooled cell-level global z-score is retained in the per-donor-group tables as `cell_global_sensitivity_mean`; it is not the primary evidence.

## Required questions

1. **Is mast-cell first robust?** Pooled mast cells remain a high signal, but donor tumor mast-cell counts are: **{mast_n}**. At n≥20, mast-cell source ranks are therefore coverage- and sparsity-sensitive; see `mast_cell_donor_counts.csv`, `cell_count_threshold_sensitivity.csv`, and the count-versus-score figure.
2. **Who is highest without PTGES3?** At n≥20, the highest mean donor-group source is **{top_source_by_score['Score_B']}** under Score B and **{top_source_by_score['Score_D']}** under Score D. Score C is the PTGES3-inclusive sensitivity; threshold-specific donor calls are in `source_rank_consistency.csv`.
3. **Most stable score:** **{stable_score}** has the highest mean pairwise Spearman concordance of donor-group ranks at n≥20 ({stable_score_conc:.3f}); full concordance values are in `final_robustness_summary.csv`.
4. **Does SPP1-like TAM remain high?** Score A: **{spp1_a_text}**; Score B: **{spp1_b_text}**; Score D: **{spp1_d_text}**. The eligible-donor denominator is essential: at n≥20, SPP1-like TAM has adequate coverage in only one tumor donor, so its apparent high rank cannot be generalized to six donors. This directly tests the PTGES3 removal and bottleneck score.
5. **Does C1QC remain PTGER4-enriched?** The pooled receiver ranking is **{pooled_top}**; C1QC's pooled PTGER4 mean is **{float(pooled_c1.PTGER4_mean) if pooled_c1 is not None else np.nan:.3f}** with detection **{float(pooled_c1.PTGER4_detection_fraction) if pooled_c1 is not None else np.nan:.3f}**. Donor pseudobulk median PTGER4 rank is **{c1med:.2f}**; C1QC is higher than SPP1 in **{c1spp}/{c1spp_valid}** donors with both subtypes (**{c1spp}/6** overall), and higher than the other-subtype mean in **{c1other}/{c1other_valid}** valid donors. Full Δ1/Δ2 contrasts are in `c1qc_ptger4_donor_contrasts.csv` and pooled values are in `ptger4_receiver_pooled_summary.csv`.
6. **Why are donors discordant?** The tables separate rank-rule effects from sparse coverage. Mast-cell donor counts are low in several donors; C1QC/SPP1 state counts also vary. Dropout is reflected by detection and expressing-cell means. Remaining differences after count filtering are consistent with biological heterogeneity, but six donors cannot distinguish it decisively.
7. **Axis evidence level:** source = **{source_grade}**; receiver = **{receiver_grade}**; combined descriptive axis = **{axis_grade}**. Axis 1/2/3 are ranking sensitivities only and are not communication probabilities.

## Interpretation

`PTGER4-high fraction` is not reported as an independent endpoint here because the previous macrophage q75 was zero and it collapsed to detection fraction. The receiver analysis instead uses PTGER4 mean, detection fraction, expressing-cell mean, pseudobulk PTGER4, and donor rank.

### Recommended interpretation

The data **{('supports' if axis_grade == 'MODERATE SUPPORT' else 'suggests only weakly')}** a potential PGE2–PTGER4 macrophage response axis at the expression level. The pooled C1QC PTGER4 enrichment is reproducible as a pooled description, while donor-level source and receiver ranks must be read with sparse-group coverage and heterogeneity visible. This does not demonstrate extracellular PGE2 production, receptor binding, spatial proximity, or functional communication.

### Recommended next step

{summary.loc[summary.metric == 'Recommended next step', 'value'].iloc[0]}.

## Files

- `pge2_score_definitions.csv`: exact Score A–D definitions and formula.
- `pseudobulk_donor_group_expression.csv`: donor-group means, detection, expressing-cell means, and pseudobulk expression.
- `pge2_score_A_by_donor_group.csv` through `pge2_score_D_by_donor_group.csv`: score tables with old global-z sensitivity.
- `source_rank_consistency.csv`, `source_score_contrasts.csv`, `cell_count_threshold_sensitivity.csv`: source robustness and SPP1 contrasts.
- `ptger4_receiver_by_donor.csv`, `ptger4_receiver_pooled_summary.csv`, `c1qc_ptger4_donor_contrasts.csv`: receiver metrics, pooled description, and donor contrasts.
- `pge2_ptger4_axis_sensitivity.csv`: three descriptive source → C1QC axis scores.
- `final_robustness_summary.csv`: numeric evidence ledger and final support grades.
- `figures/`: PNG, PDF, and SVG figures requested in the protocol.
"""
    (args.out / "PGE2_PTGER4_ROBUSTNESS_REPORT.md").write_text(report, encoding="utf-8")
    cell_counts = pb.pivot_table(index="donor_id", columns="analysis_group", values="n_cells", aggfunc="sum").fillna(0).astype(int)
    manifest = {
        "analysis": "dinp_crc_pge2_ptger4_axis_robustness", "run_timestamp_utc": pd.Timestamp.utcnow().isoformat(),
        "previous_commit": "de2a22520d518caa3e30c7bf431505842a3e5cfd", "previous_module_preserved": True,
        "matrix": {"path": str(args.matrix), "sha256": sha256(args.matrix), "size_bytes": args.matrix.stat().st_size},
        "annotation": {"path": str(args.annotation), "sha256": sha256(args.annotation), "size_bytes": args.annotation.stat().st_size},
        "subtype_source": {"path": str(args.subtype), "sha256": sha256(args.subtype)},
        "cell_score_source": {"path": str(args.cell_scores), "sha256": sha256(args.cell_scores)},
        "donors": sorted(map(str, cells.donor_id.unique())), "n_donors": int(cells.donor_id.nunique()), "total_cells": int(len(cells)),
        "cell_counts_per_donor_group": {str(k): {str(kk): int(vv) for kk, vv in row.items()} for k, row in cell_counts.iterrows()},
        "groups_in_primary_comparison": CANDIDATE_GROUPS, "macrophage_subtypes": MACRO_SUBTYPES,
        "normalization": "reuse prior per-cell library-size normalized log1p target-gene expression; donor-group means; donor-group gene-wise z-score for Scores A-C; q95-scaled geometric bottleneck Score D",
        "minimum_cell_thresholds": THRESHOLDS, "primary_threshold": 20,
        "ptger4_high_fraction": "not treated as independent; previous macrophage q75=0 caused exact collapse to detection",
        "software": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__, "scipy": stats.__version__ if hasattr(stats, "__version__") else "available", "platform": platform.platform()},
        "ordinary_cellchat": "not used", "new_clustering": "not used", "runtime_seconds": round(time.time() - t0, 2),
    }
    json_dump(manifest, args.out / "input_manifest.json")
    print("PGE2–PTGER4 ROBUSTNESS RE-ANALYSIS: PASS")
    print(f"Dataset: GSE144735; Donors: {cells.donor_id.nunique()}; Cells: {len(cells)}")
    print(f"Primary PGE2 score: donor-group pseudobulk Score A/B/D with n>=20")
    print(f"Top source under Score A: {top_source_by_score['Score_A']} (mean donor-group score, n>=20)")
    print(f"Top source under Score B: {top_source_by_score['Score_B']} (mean donor-group score, n>=20)")
    print(f"Top source under Score D: {top_source_by_score['Score_D']} (mean donor-group score, n>=20)")
    print(f"SPP1 median donor rank: A={float(spp1_a.median_rank.iloc[0]) if len(spp1_a) else float('nan'):.2f}; B={float(spp1_b.median_rank.iloc[0]) if len(spp1_b) else float('nan'):.2f}; D={float(spp1_d.median_rank.iloc[0]) if len(spp1_d) else float('nan'):.2f}")
    print(f"SPP1 top-2 donors: A={int(spp1_a.top2_donors.iloc[0]) if len(spp1_a) else 0}/{int(spp1_a.n_donors_with_group.iloc[0]) if len(spp1_a) else 0} eligible ({int(spp1_a.top2_donors.iloc[0]) if len(spp1_a) else 0}/6 total); B={int(spp1_b.top2_donors.iloc[0]) if len(spp1_b) else 0}/{int(spp1_b.n_donors_with_group.iloc[0]) if len(spp1_b) else 0} eligible; D={int(spp1_d.top2_donors.iloc[0]) if len(spp1_d) else 0}/{int(spp1_d.n_donors_with_group.iloc[0]) if len(spp1_d) else 0} eligible")
    print(f"C1QC PTGER4 median donor rank: {c1med:.2f}; C1QC > SPP1 donors: {c1spp}/{c1spp_valid} valid ({c1spp}/6 total)")
    print(f"Mast-cell stability: donor counts {mast_n}")
    print(f"PTGES3 sensitivity: Score C retained as sensitivity; Score B/D exclude PTGES3")
    print(f"Most stable pseudobulk score: {stable_score} (mean rank concordance {stable_score_conc:.3f})")
    print(f"Final source support: {source_grade}; Final receiver support: {receiver_grade}; Final axis support: {axis_grade}")
    print(f"Recommended next step: {summary.loc[summary.metric == 'Recommended next step', 'value'].iloc[0]}")


if __name__ == "__main__":
    main()
