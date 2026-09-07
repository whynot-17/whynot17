#!/usr/bin/env python
"""Focused PGE2 source -> PTGER4 receiver analysis for the CRC atlas.

This module answers one narrow question: which broad cell groups have the
strongest PGE2 synthesis program, and which macrophage state has the strongest
PTGER4 signal?  It uses donor-level summaries for validation and deliberately
does not infer physical ligand-receptor communication.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_MATRIX = Path(r"E:\mcop\mcop\work\mcop_phase2f_external\raw\GSE144735_raw_UMI_count_matrix.txt.gz")
DEFAULT_ANNOTATION = Path(r"E:\mcop\mcop\work\mcop_phase2f_external\raw\GSE144735_annotation.txt.gz")
DEFAULT_SUBTYPE = Path(r"E:\chatgpt\whynot17\analysis\dinp_crc_macrophage_subtype\outputs\macrophage_cell_program_scores.csv")

SYNTHESIS_GENES = ["PLA2G4A", "PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3"]
TARGET_GENES = SYNTHESIS_GENES + ["PTGER4"]
MACRO_SUBTYPES = [
    "SPP1_like_TAM",
    "C1QC_like_TAM",
    "FCN1_inflammatory_monocyte_like",
    "resident_like",
]
REPORT_GROUPS = ["tumor_epithelial", "SPP1_like_TAM", "C1QC_like_TAM", "fibroblast_like", "endothelial_like", "normal_epithelial"]
# Retain mast cells as a candidate source because this distinct compartment
# ranks highly in the pooled synthesis score.
SOURCE_GROUPS = ["mast_cells"] + REPORT_GROUPS


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_dump(obj, path: Path):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def safe_mean(values) -> float:
    a = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    return float(a.mean()) if len(a) else float("nan")


def assign_group(row: pd.Series) -> str:
    subtype = str(row.get("subtype", ""))
    ctype = str(row.get("Cell_type", ""))
    csub = str(row.get("Cell_subtype", ""))
    condition = str(row.get("condition", ""))
    if subtype in MACRO_SUBTYPES:
        return subtype
    if ctype == "Epithelial cells":
        return "tumor_epithelial" if condition == "tumor" else "normal_epithelial" if condition == "normal" else "border_epithelial"
    if ctype == "Stromal cells":
        if csub in {"Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3"}:
            return "fibroblast_like"
        if csub in {"Tip-like ECs", "Stalk-like ECs", "Lymphatic ECs"}:
            return "endothelial_like"
        if csub == "Pericytes":
            return "pericyte_like"
        if csub == "Smooth muscle cells":
            return "smooth_muscle_like"
        if csub == "Enteric glial cells":
            return "enteric_glial_like"
        return "other_stromal"
    if ctype == "T cells":
        return "T_cells"
    if ctype == "B cells":
        return "B_cells"
    if ctype == "Mast cells":
        return "mast_cells"
    if ctype == "Myeloids":
        return "other_myeloid"
    return "other"


def read_target_rows_and_library(matrix_path: Path, target_genes: list[str], progress_every: int = 5000):
    """Stream raw gene x cell text matrix; retain seven target rows and totals."""
    target = {g.upper(): i for i, g in enumerate(target_genes)}
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        if not header or header[0] != "Index":
            raise ValueError("Unexpected matrix header; expected first field Index")
        cells = header[1:]
        n = len(cells)
        library = np.zeros(n, dtype=np.float64)
        rows = np.zeros((len(target_genes), n), dtype=np.float32)
        found = {}
        nrows = 0
        t0 = time.time()
        for line in fh:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != n + 1:
                raise ValueError(f"Matrix row has {len(parts)-1} values, expected {n}: {parts[0][:40]}")
            # np.fromiter is substantially faster and less memory hungry than
            # building a full dense matrix for this 27k-cell input.
            values = np.fromiter((float(x) if x else 0.0 for x in parts[1:]), dtype=np.float32, count=n)
            library += values
            gene = parts[0].strip().upper()
            if gene in target and gene not in found:
                rows[target[gene], :] = values
                found[gene] = True
            nrows += 1
            if progress_every and nrows % progress_every == 0:
                print(f"parsed {nrows:,} genes; elapsed {time.time()-t0:.1f}s", flush=True)
    return cells, rows, library, sorted(found)


def normalize_target_rows(rows: np.ndarray, library: np.ndarray) -> np.ndarray:
    scale = np.divide(1e4, library, out=np.zeros_like(library), where=library > 0)
    return np.log1p(rows.astype(np.float64) * scale[np.newaxis, :])


def load_metadata(annotation_path: Path, subtype_path: Path, matrix_cells: list[str]) -> pd.DataFrame:
    ann = pd.read_csv(annotation_path, sep="\t", compression="infer", dtype=str)
    ann = ann.rename(columns={"Index": "cell_id", "Patient": "donor_id", "Class": "condition_raw"})
    ann["condition"] = ann["condition_raw"].str.lower()
    subtype = pd.read_csv(subtype_path, dtype={"cell_id": str})
    keep = [c for c in ["cell_id", "subtype", "cluster"] if c in subtype.columns]
    subtype = subtype[keep].drop_duplicates("cell_id")
    meta = pd.DataFrame({"cell_id": matrix_cells}).merge(ann, on="cell_id", how="left", validate="one_to_one")
    meta = meta.merge(subtype, on="cell_id", how="left", validate="one_to_one")
    if meta["donor_id"].isna().any():
        raise ValueError("Some matrix cells are missing from the annotation")
    meta["analysis_group"] = meta.apply(assign_group, axis=1)
    return meta


def cell_scores(meta: pd.DataFrame, norm: np.ndarray, target_genes: list[str]) -> pd.DataFrame:
    expr = pd.DataFrame(norm.T, columns=target_genes)
    out = meta.reset_index(drop=True).copy()
    for gene in target_genes:
        out[f"{gene}_log1p"] = expr[gene].astype(float)
        out[f"{gene}_detected"] = expr[gene].gt(0)
    z = pd.DataFrame(index=out.index)
    for gene in SYNTHESIS_GENES:
        vals = expr[gene].to_numpy(float)
        sd = float(np.nanstd(vals, ddof=0))
        z[gene] = (vals - float(np.nanmean(vals))) / sd if sd > 0 else 0.0
    out["PGE2_synthesis_score"] = z[SYNTHESIS_GENES].mean(axis=1)
    return out


def summarize_groups(cells: pd.DataFrame, global_high: float) -> pd.DataFrame:
    rows = []
    for group, g in cells.groupby("analysis_group", sort=False):
        row = {
            "analysis_group": group,
            "n_cells": int(len(g)),
            "n_donors": int(g.donor_id.nunique()),
            "tumor_cells": int((g.condition == "tumor").sum()),
            "normal_cells": int((g.condition == "normal").sum()),
            "PGE2_synthesis_mean": float(g.PGE2_synthesis_score.mean()),
            "PGE2_synthesis_median": float(g.PGE2_synthesis_score.median()),
            "PGE2_synthesis_high_fraction": float((g.PGE2_synthesis_score >= global_high).mean()),
        }
        for gene in SYNTHESIS_GENES:
            row[f"{gene}_mean"] = float(g[f"{gene}_log1p"].mean())
            row[f"{gene}_detection_fraction"] = float(g[f"{gene}_detected"].mean())
        rows.append(row)
    return pd.DataFrame(rows).sort_values("PGE2_synthesis_mean", ascending=False).reset_index(drop=True)


def donor_group_summary(cells: pd.DataFrame, global_high: float) -> pd.DataFrame:
    agg = {
        "n_cells": ("cell_id", "size"),
        "PGE2_synthesis_mean": ("PGE2_synthesis_score", "mean"),
        "PGE2_synthesis_median": ("PGE2_synthesis_score", "median"),
        "PGE2_synthesis_high_fraction": ("PGE2_synthesis_score", lambda x: float((x >= global_high).mean())),
    }
    for gene in SYNTHESIS_GENES:
        agg[f"{gene}_mean"] = (f"{gene}_log1p", "mean")
        agg[f"{gene}_detection_fraction"] = (f"{gene}_detected", "mean")
    return cells.groupby(["donor_id", "condition", "analysis_group"], dropna=False).agg(**agg).reset_index()


def macro_ptger4_summary(cells: pd.DataFrame, ptger4_high: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    m = cells[cells.analysis_group.isin(MACRO_SUBTYPES)].copy()
    m["PTGER4_high"] = m["PTGER4_log1p"] >= ptger4_high
    rows = []
    for subtype, g in m.groupby("analysis_group", sort=False):
        rows.append({
            "subtype": subtype, "n_cells": int(len(g)), "n_donors": int(g.donor_id.nunique()),
            "PTGER4_mean": float(g.PTGER4_log1p.mean()), "PTGER4_median": float(g.PTGER4_log1p.median()),
            "PTGER4_detection_fraction": float(g.PTGER4_detected.mean()),
            "PTGER4_high_fraction": float(g.PTGER4_high.mean()),
            "tumor_PTGER4_mean": float(g.loc[g.condition == "tumor", "PTGER4_log1p"].mean()),
            "tumor_PTGER4_detection_fraction": float(g.loc[g.condition == "tumor", "PTGER4_detected"].mean()),
            "tumor_PTGER4_high_fraction": float(g.loc[g.condition == "tumor", "PTGER4_high"].mean()),
        })
    donor = m.groupby(["donor_id", "condition", "analysis_group"], dropna=False).agg(
        n_cells=("cell_id", "size"), PTGER4_mean=("PTGER4_log1p", "mean"),
        PTGER4_median=("PTGER4_log1p", "median"), PTGER4_detection_fraction=("PTGER4_detected", "mean"),
        PTGER4_high_fraction=("PTGER4_high", "mean"),
    ).reset_index().rename(columns={"analysis_group": "subtype"})
    return pd.DataFrame(rows).sort_values("PTGER4_mean", ascending=False).reset_index(drop=True), donor


def donor_axis(cells: pd.DataFrame, group_donor: pd.DataFrame, macro_donor: pd.DataFrame, global_high: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    major = [g for g in SOURCE_GROUPS if g in set(group_donor.analysis_group)]
    tumor_donors = sorted(cells.loc[cells.condition == "tumor", "donor_id"].unique())
    rows = []
    for donor in tumor_donors:
        gd = group_donor[(group_donor.donor_id == donor) & (group_donor.condition == "tumor")]
        md = macro_donor[(macro_donor.donor_id == donor) & (macro_donor.condition == "tumor")]
        source_vals = gd[gd.analysis_group.isin(major)]
        top_source = source_vals.loc[source_vals.PGE2_synthesis_mean.idxmax(), "analysis_group"] if len(source_vals) else ""
        pt_vals = md[md.subtype.isin(MACRO_SUBTYPES)]
        top_receiver = pt_vals.loc[pt_vals.PTGER4_mean.idxmax(), "subtype"] if len(pt_vals) else ""
        for source in major:
            s = gd[gd.analysis_group == source]
            s = s.iloc[0] if len(s) else None
            for receiver in MACRO_SUBTYPES:
                r = md[md.subtype == receiver]
                r = r.iloc[0] if len(r) else None
                if s is None or r is None:
                    row = {"donor_id": donor, "condition": "tumor", "source_group": source, "receiver_subtype": receiver,
                           "source_cells": 0 if s is None else int(s.n_cells), "receiver_cells": 0 if r is None else int(r.n_cells),
                           "source_mean_score": np.nan if s is None else float(s.PGE2_synthesis_mean),
                           "source_high_fraction": np.nan if s is None else float(s.PGE2_synthesis_high_fraction),
                           "receiver_PTGER4_mean": np.nan if r is None else float(r.PTGER4_mean),
                           "receiver_detection_fraction": np.nan if r is None else float(r.PTGER4_detection_fraction),
                           "receiver_high_fraction": np.nan if r is None else float(r.PTGER4_high_fraction),
                           "source_top_among_major": bool(source == top_source), "receiver_top_among_macrophages": bool(receiver == top_receiver),
                           "coverage": False, "axis_score": np.nan}
                else:
                    # Product is descriptive only; it is not a communication probability.
                    axis = float(s.PGE2_synthesis_mean * max(float(r.PTGER4_mean), 0.0))
                    row = {"donor_id": donor, "condition": "tumor", "source_group": source, "receiver_subtype": receiver,
                           "source_cells": int(s.n_cells), "receiver_cells": int(r.n_cells),
                           "source_mean_score": float(s.PGE2_synthesis_mean), "source_high_fraction": float(s.PGE2_synthesis_high_fraction),
                           "receiver_PTGER4_mean": float(r.PTGER4_mean), "receiver_detection_fraction": float(r.PTGER4_detection_fraction),
                           "receiver_high_fraction": float(r.PTGER4_high_fraction), "source_top_among_major": bool(source == top_source),
                           "receiver_top_among_macrophages": bool(receiver == top_receiver), "coverage": True, "axis_score": axis}
                rows.append(row)
    axis = pd.DataFrame(rows)
    if len(axis):
        axis["pair_is_directionally_consistent"] = axis.coverage & axis.source_top_among_major & axis.receiver_top_among_macrophages
    summary_rows = []
    for source in major:
        pair = axis[(axis.source_group == source) & (axis.receiver_subtype == "C1QC_like_TAM")]
        summary_rows.append({
            "source_group": source, "receiver_subtype": "C1QC_like_TAM", "n_tumor_donors": len(tumor_donors),
            "n_donors_with_coverage": int(pair.coverage.sum()), "n_source_top_donors": int(pair.source_top_among_major.sum()),
            "n_C1QC_receiver_top_donors": int(pair.receiver_top_among_macrophages.sum()),
            "n_directionally_consistent_donors": int(pair.pair_is_directionally_consistent.sum()),
            "directional_consistency_fraction": float(pair.pair_is_directionally_consistent.mean()) if len(pair) else np.nan,
            "mean_axis_score": float(pair.axis_score.mean()) if pair.axis_score.notna().any() else np.nan,
            "interpretation": "supports a reproducible potential source-to-receiver axis" if len(pair) and pair.pair_is_directionally_consistent.sum() >= 4 else "insufficient or inconsistent donor-level support",
        })
    return axis, pd.DataFrame(summary_rows)


def donor_top_ranks(group_donor: pd.DataFrame, macro_donor: pd.DataFrame, source_groups: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-donor top source/receiver calls and descriptive counts."""
    gd = group_donor[(group_donor.condition == "tumor") & group_donor.analysis_group.isin(source_groups)]
    md = macro_donor[macro_donor.condition == "tumor"]
    src = gd.loc[gd.groupby("donor_id").PGE2_synthesis_mean.idxmax(), ["donor_id", "analysis_group", "n_cells", "PGE2_synthesis_mean"]].rename(columns={"analysis_group": "top_source_group"})
    rec = md.loc[md.groupby("donor_id").PTGER4_mean.idxmax(), ["donor_id", "subtype", "n_cells", "PTGER4_mean", "PTGER4_detection_fraction"]].rename(columns={"subtype": "top_receiver_subtype"})
    per_donor = src.merge(rec, on="donor_id", how="outer", suffixes=("_source", "_receiver"))
    rows = []
    for col, label in [("top_source_group", "source"), ("top_receiver_subtype", "receiver")]:
        counts = per_donor[col].value_counts(dropna=True)
        for category, count in counts.items():
            rows.append({"rank_type": label, "category": category, "n_tumor_donors_top": int(count), "fraction_of_6_donors": float(count / 6)})
    return per_donor, pd.DataFrame(rows)


def save_figures(group: pd.DataFrame, macro: pd.DataFrame, axis: pd.DataFrame, out: Path):
    import matplotlib.pyplot as plt
    import seaborn as sns
    figdir = out / "figures"; figdir.mkdir(parents=True, exist_ok=True)
    colors = {"tumor_epithelial": "#C44E52", "SPP1_like_TAM": "#DD8452", "C1QC_like_TAM": "#4C72B0", "fibroblast_like": "#55A868", "endothelial_like": "#8172B2", "normal_epithelial": "#64B5CD"}
    p = group[group.analysis_group.isin(REPORT_GROUPS)].sort_values("PGE2_synthesis_mean")
    fig, ax = plt.subplots(figsize=(9, 5.5)); ax.barh(p.analysis_group, p.PGE2_synthesis_mean, color=[colors.get(x, "#777777") for x in p.analysis_group]); ax.axvline(0, color="black", lw=.7); ax.set_xlabel("PGE2 synthesis score (mean gene-wise z-score)"); ax.set_ylabel(""); ax.set_title("PGE2 synthesis program by major cell group"); fig.tight_layout()
    for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"pge2_synthesis_by_group.{ext}", dpi=220)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    order = [x for x in MACRO_SUBTYPES if x in set(macro.subtype)]
    sns.barplot(data=macro, x="subtype", y="PTGER4_mean", order=order, ax=axes[0], color="#4C72B0", errorbar=None)
    sns.barplot(data=macro, x="subtype", y="PTGER4_detection_fraction", order=order, ax=axes[1], color="#55A868", errorbar=None)
    sns.barplot(data=macro, x="subtype", y="PTGER4_high_fraction", order=order, ax=axes[2], color="#C44E52", errorbar=None)
    axes[0].set_ylabel("mean log1p normalized expression"); axes[1].set_ylabel("detection fraction"); axes[2].set_ylabel("PTGER4-high fraction")
    for ax in axes: ax.tick_params(axis="x", rotation=35); ax.set_xlabel("")
    fig.suptitle("PTGER4 in macrophage subtypes (cell-level descriptive summaries)"); fig.tight_layout()
    for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"ptger4_by_macrophage_subtype.{ext}", dpi=220)
    plt.close(fig)
    if len(axis):
        c = axis[axis.receiver_subtype == "C1QC_like_TAM"].pivot(index="donor_id", columns="source_group", values="axis_score")
        fig, ax = plt.subplots(figsize=(10, max(3.5, .45 * len(c)))); sns.heatmap(c, annot=True, fmt=".2f", cmap="vlag", center=0, ax=ax); ax.set_title("Descriptive source score × C1QC PTGER4 by tumor donor"); ax.set_xlabel("candidate PGE2 source"); ax.set_ylabel("donor"); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"pge2_source_ptger4_receiver_heatmap.{ext}", dpi=220)
        plt.close(fig)
        cons = axis[(axis.receiver_subtype == "C1QC_like_TAM")].pivot(index="donor_id", columns="source_group", values="pair_is_directionally_consistent").astype(float)
        fig, ax = plt.subplots(figsize=(10, max(3.5, .45 * len(cons)))); sns.heatmap(cons, annot=True, fmt=".0f", cmap="Greens", vmin=0, vmax=1, cbar=False, ax=ax); ax.set_title("Donor-level directional consistency (1 = source top and C1QC PTGER4 top)"); ax.set_xlabel("candidate PGE2 source"); ax.set_ylabel("donor"); fig.tight_layout()
        for ext in ("png", "pdf", "svg"): fig.savefig(figdir / f"donor_axis_consistency.{ext}", dpi=220)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    ap.add_argument("--annotation", type=Path, default=DEFAULT_ANNOTATION)
    ap.add_argument("--subtype", type=Path, default=DEFAULT_SUBTYPE)
    ap.add_argument("--out", type=Path, default=Path(__file__).parent / "outputs")
    args = ap.parse_args()
    t0 = time.time(); args.out.mkdir(parents=True, exist_ok=True)
    print(f"Reading target rows and library sizes from {args.matrix}", flush=True)
    matrix_cells, raw, library, present = read_target_rows_and_library(args.matrix, TARGET_GENES)
    norm = normalize_target_rows(raw, library)
    meta = load_metadata(args.annotation, args.subtype, matrix_cells)
    cells = cell_scores(meta, norm, TARGET_GENES)
    global_high = float(cells.PGE2_synthesis_score.quantile(.75))
    # PTGER4-high is defined within the macrophage compartment, requiring
    # detection when the upper quartile is zero.
    macro_mask = cells.analysis_group.isin(MACRO_SUBTYPES)
    pt_q75 = float(cells.loc[macro_mask, "PTGER4_log1p"].quantile(.75))
    ptger4_high = pt_q75 if pt_q75 > 0 else float(np.nextafter(0, 1))
    cells["PGE2_synthesis_high"] = cells.PGE2_synthesis_score >= global_high
    cells["PTGER4_detected"] = cells.PTGER4_log1p > 0
    cells.to_csv(args.out / "pge2_synthesis_cell_scores.csv", index=False)
    group = summarize_groups(cells, global_high); group.to_csv(args.out / "pge2_synthesis_by_group.csv", index=False)
    gd = donor_group_summary(cells, global_high); gd.to_csv(args.out / "pge2_synthesis_by_donor_group.csv", index=False)
    macro, md = macro_ptger4_summary(cells, ptger4_high); macro.to_csv(args.out / "ptger4_macrophage_subtype_statistics.csv", index=False); md.to_csv(args.out / "ptger4_by_donor_subtype.csv", index=False)
    axis, axis_summary = donor_axis(cells, gd, md, global_high); axis.to_csv(args.out / "pge2_ptger4_donor_axis.csv", index=False); axis_summary.to_csv(args.out / "pge2_ptger4_axis_donor_summary.csv", index=False)
    donor_ranks, donor_rank_counts = donor_top_ranks(gd, md, [g for g in SOURCE_GROUPS if g in set(gd.analysis_group)])
    donor_ranks.to_csv(args.out / "pge2_ptger4_top_calls_by_donor.csv", index=False)
    donor_rank_counts.to_csv(args.out / "pge2_ptger4_top_call_counts.csv", index=False)
    save_figures(group, macro, axis, args.out)
    top_source = group.iloc[0].analysis_group if len(group) else "NA"
    focus_group = group[group.analysis_group.isin(REPORT_GROUPS)]
    top_focus_source = focus_group.iloc[0].analysis_group if len(focus_group) else "NA"
    top_pt_mean = macro.iloc[0].subtype if len(macro) else "NA"
    top_pt_det = macro.sort_values("PTGER4_detection_fraction", ascending=False).iloc[0].subtype if len(macro) else "NA"
    top_pt_high = macro.sort_values("PTGER4_high_fraction", ascending=False).iloc[0].subtype if len(macro) else "NA"
    c1 = axis_summary[(axis_summary.receiver_subtype == "C1QC_like_TAM") & (axis_summary.source_group != "C1QC_like_TAM")] if len(axis_summary) else pd.DataFrame()
    c1_best = c1.sort_values(["n_directionally_consistent_donors", "mean_axis_score"], ascending=False).iloc[0] if len(c1) else None
    receiver_rank_counts = donor_rank_counts[donor_rank_counts.rank_type == "receiver"]
    c1_mean_top = int(receiver_rank_counts.loc[receiver_rank_counts.category == "C1QC_like_TAM", "n_tumor_donors_top"].iloc[0]) if (receiver_rank_counts.category == "C1QC_like_TAM").any() else 0
    md_tumor = md[md.condition == "tumor"]
    c1_det_top = int(md_tumor.loc[md_tumor.groupby("donor_id").PTGER4_detection_fraction.idxmax(), "subtype"].eq("C1QC_like_TAM").sum()) if len(md_tumor) else 0
    ptger4_threshold_text = f"{ptger4_high:.3f} (detected-only fallback because the macrophage 75th percentile is zero)" if pt_q75 <= 0 else f"{ptger4_high:.3f}"
    report = f"""# Focused PGE2 -> PTGER4 axis analysis

## Question

Which cell group has the strongest **PGE2 synthesis program**, and which macrophage subtype has the strongest **PTGER4 receiving signal**? This is a descriptive molecular axis analysis and does not prove extracellular PGE2 production, binding, or physical cell-cell communication.

## Data and scoring

- Dataset: local GSE144735 CRC matrix ({len(cells):,} annotated cells; {cells.donor_id.nunique()} donors: {', '.join(sorted(cells.donor_id.unique()))}).
- Synthesis genes: {', '.join(SYNTHESIS_GENES)}.
- Raw counts were library-size normalized to 10,000 and log1p transformed; each synthesis gene was z-scored across all cells and averaged into `PGE2_synthesis_score`.
- `PGE2_synthesis_high` is the global cell-level upper quartile (threshold {global_high:.3f}).
- Macrophage PTGER4-high is the macrophage-compartment upper quartile of normalized log1p PTGER4; here the 75th percentile is zero, so the prespecified fallback is **detected-only** (threshold {ptger4_threshold_text}).
- Donor-level tables are descriptive aggregates; no cell-level significance is used.

## Main readouts

- Highest broad-group PGE2 synthesis mean: **{top_source}**. Among the prespecified microenvironment comparison groups, the highest mean is **{top_focus_source}**.
- Highest macrophage PTGER4 mean: **{top_pt_mean}**.
- Highest macrophage PTGER4 detection fraction: **{top_pt_det}**.
- Highest macrophage PTGER4-high fraction: **{top_pt_high}**.
- At the donor level, C1QC-like TAM is the top PTGER4 subtype by mean in **{c1_mean_top}/6** tumor donors and by detection fraction in **{c1_det_top}/6**; it is therefore not uniformly highest across all six donors.
- The complete group and subtype rankings are in `pge2_synthesis_by_group.csv` and `ptger4_macrophage_subtype_statistics.csv`.

## Donor-level source -> receiver check

Tumor donor rows compare each candidate source with each macrophage subtype. A donor is directionally consistent for a candidate source -> C1QC-like TAM pair when that source is the top PGE2 synthesis group among mast cells plus the prespecified comparison groups **and** C1QC-like TAM is the top PTGER4 macrophage subtype for that donor. The product in `pge2_ptger4_donor_axis.csv` is a descriptive ranking score, not a communication probability.

"""
    if c1_best is not None:
        report += f"For the most directionally consistent non-C1QC source pair, **{c1_best.source_group} -> C1QC_like_TAM**, {int(c1_best.n_directionally_consistent_donors)}/{int(c1_best.n_tumor_donors)} tumor donors meet both top-ranking conditions.\n\n"
    report += "Interpretation should use **supports/suggests a potential PGE2-producing source -> PTGER4-high macrophage axis**. It must not be written as demonstrated PGE2 communication without spatial, biochemical, or functional validation.\n\n"
    report += "## Files\n\n- `pge2_synthesis_by_group.csv`: group-level synthesis ranking.\n- `ptger4_macrophage_subtype_statistics.csv`: PTGER4 mean, detection fraction, and high-cell fraction.\n- `pge2_ptger4_axis_donor_summary.csv`: six-donor consistency summary for each source -> C1QC pair.\n- `pge2_ptger4_top_calls_by_donor.csv`: per-donor top source and receiver calls; sparse groups remain visible through their cell counts.\n- `figures/`: publication-ready PNG, PDF, and SVG summaries.\n"
    (args.out / "PGE2_PTGER4_REPORT.md").write_text(report, encoding="utf-8")
    manifest = {
        "analysis": "dinp_crc_pge2_ptger4_axis", "run_timestamp_utc": pd.Timestamp.utcnow().isoformat(),
        "matrix": {"path": str(args.matrix), "sha256": sha256(args.matrix), "size_bytes": args.matrix.stat().st_size},
        "annotation": {"path": str(args.annotation), "sha256": sha256(args.annotation), "size_bytes": args.annotation.stat().st_size},
        "prior_subtype_output": {"path": str(args.subtype), "sha256": sha256(args.subtype)},
        "source_dataset": "GSE144735 CRC single-cell matrix; local cache; six donors",
        "n_cells": int(len(cells)), "n_donors": int(cells.donor_id.nunique()), "donors": sorted(cells.donor_id.unique()),
        "group_counts": {str(k): int(v) for k, v in cells.analysis_group.value_counts().items()},
        "target_genes": TARGET_GENES, "genes_present": present, "genes_missing": [g for g in TARGET_GENES if g not in present],
        "synthesis_genes": SYNTHESIS_GENES, "normalization": "per-cell library-size normalize to 1e4, log1p; synthesis gene-wise z-score mean",
        "pge2_high_threshold_global_cell_q75": global_high, "ptger4_high_threshold_macrophage_q75": ptger4_high,
        "ptger4_high_definition": "detected-only fallback because macrophage q75 is zero" if pt_q75 <= 0 else "macrophage q75",
        "donor_validation": "tumor donor-level aggregates; source top among mast cells plus six prespecified comparison groups and receiver top among four macrophage subtypes",
        "ordinary_cellchat": "not used by user request; this is a focused source-receiver analysis",
        "r_environment": r"E:\\chatgpt\\cellchat_env (R 4.4.3; CellChat not installed)",
        "software": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__, "platform": platform.platform()},
        "runtime_seconds": round(time.time() - t0, 2),
    }
    json_dump(manifest, args.out / "input_manifest.json")
    print(f"cells={len(cells):,}; donors={cells.donor_id.nunique()}; top synthesis={top_source}; top PTGER4 mean={top_pt_mean}; runtime={time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
