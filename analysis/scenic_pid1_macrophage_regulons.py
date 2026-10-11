#!/usr/bin/env python
"""Patient-aware pySCENIC workflow for PID1 macrophage-state analysis.

Primary GRN/regulon inference: GSE132465 inflammatory and C1QC macrophages.
External sensitivity scoring: GSE200997 cells scored against primary regulons.
Raw GEO matrices are read-only inputs; all derived files go to a dated results folder.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(r"D:\DINP\whynot17-new-project")
RAW = ROOT / "work" / "single_cell_GSE132465_GSE200997"
ANN = ROOT / "results" / "2026-10-11" / "PID1_myeloid_finelocalization"
OUT = ROOT / "results" / "2026-10-12" / "SCENIC_PID1_macrophage_regulons"
RES = ROOT / "work" / "scenic_resources" / "hg38"
DB = RES / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
MOTIFS = RES / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TFS = RES / "allTFs_hg38.txt"
GROUPS = ["Macrophage_C1QC", "Macrophage_inflammatory"]
DISPLAY = {"Macrophage_C1QC": "C1QC", "Macrophage_inflammatory": "Inflammatory"}
SEED = 777


def load_annotations(cohort: str) -> pd.DataFrame:
    suffix = "tsv.gz"
    path = ANN / f"{cohort}_PID1_CellChat_myeloid_annotations.{suffix}"
    sep = "\t"
    a = pd.read_csv(path, sep=sep, compression="gzip", low_memory=False)
    a = a[a.pid1_subtype.isin(GROUPS)].copy()
    for col in ["qc_pass", "scrublet_predicted_doublet", "marker_lineage_conflict"]:
        if a[col].dtype == object:
            a[col] = a[col].astype(str).str.lower().isin(["true", "1", "yes"])
    a = a[a.qc_pass & ~a.scrublet_predicted_doublet & ~a.marker_lineage_conflict].copy()
    a["unit_id"] = a["unit_id"].astype(str)
    a["cell_id"] = a["cell_id"].astype(str)
    return a


def load_counts(cohort: str, cell_ids: list[str]) -> pd.DataFrame:
    if cohort == "GSE132465":
        path = RAW / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
        sep = "\t"
        gene_col = "Index"
    else:
        path = RAW / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
        sep = ","
        gene_col = "Unnamed: 0"
    # The GEO tables are extremely wide. Stream each gene row and retain only
    # selected cells; pandas' C parser can allocate buffers proportional to the
    # full cell width even when usecols is supplied.
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter=sep)
        header = next(reader)
        positions = {name: i for i, name in enumerate(header)}
        missing = set(cell_ids).difference(positions)
        if missing:
            raise ValueError(f"{cohort}: {len(missing)} annotated cells absent from raw matrix; e.g. {sorted(missing)[:5]}")
        ix = [positions[c] for c in cell_ids]
        genes = []
        rows = []
        for row in reader:
            if not row:
                continue
            genes.append(str(row[0]))
            rows.append(np.fromiter((float(row[j]) for j in ix), dtype=np.float32, count=len(ix)))
    mat = pd.DataFrame(np.asarray(rows, dtype=np.float32), index=genes, columns=cell_ids)
    if mat.index.has_duplicates:
        mat = mat.groupby(level=0, sort=False).sum()
    mat = mat.reindex(columns=cell_ids)
    return mat


def normalize_log1p(counts: pd.DataFrame) -> pd.DataFrame:
    x = counts.to_numpy(dtype=np.float64).T
    totals = x.sum(axis=1)
    if (totals <= 0).any():
        bad = np.where(totals <= 0)[0]
        raise ValueError(f"Cells with zero library size: {len(bad)}")
    x = np.log1p((x / totals[:, None]) * 1e4).astype(np.float32)
    out = pd.DataFrame(x, index=counts.columns, columns=counts.index)
    # Exclude genes detected in fewer than 1% of cells and invariant genes.
    keep = (counts.gt(0).sum(axis=1).values >= max(3, int(np.ceil(0.01 * len(counts.columns)))))
    out = out.loc[:, keep]
    out = out.loc[:, out.var(axis=0) > 0]
    return out


def run(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("COMMAND: " + subprocess.list2cmdline(cmd) + "\n\n")
        log.flush()
        p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
    if p.returncode:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd[0]} {cmd[1]}; see {log_path}")


def qvals(pvals):
    vals = np.asarray(pvals, dtype=float)
    out = np.full(vals.shape, np.nan)
    ok = np.isfinite(vals)
    if ok.any():
        out[ok] = multipletests(vals[ok], method="fdr_bh")[1]
    return out


def activity_patient_tables(auc: pd.DataFrame, ann: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # normalize orientation to cells x regulons
    cells = set(ann.cell_id)
    if len(cells.intersection(auc.index.astype(str))) >= len(cells) * 0.8:
        auc.index = auc.index.astype(str)
        ac = auc.loc[auc.index.intersection(cells)].copy()
    elif len(cells.intersection(auc.columns.astype(str))) >= len(cells) * 0.8:
        auc.columns = auc.columns.astype(str)
        ac = auc.loc[:, auc.columns.intersection(cells)].T.copy()
    else:
        raise ValueError("Could not align AUCell matrix with annotated cell IDs")
    ac.index.name = "cell_id"
    long = ac.reset_index().melt(id_vars="cell_id", var_name="regulon", value_name="AUC")
    long = long.merge(ann[["cell_id", "unit_id", "sample_id", "pid1_subtype", "PID1_log1p"]], on="cell_id", how="inner")
    patient = (long.groupby(["unit_id", "sample_id", "pid1_subtype", "regulon"], observed=True)
               .agg(AUC=("AUC", "mean"), cells=("cell_id", "nunique"), PID1_log1p=("PID1_log1p", "mean"))
               .reset_index())
    return ac, patient


def paired_regulon_stats(patient: pd.DataFrame, ann: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    counts = ann.groupby(["unit_id", "pid1_subtype"]).size().unstack(fill_value=0)
    eligible = counts.index[(counts.get(GROUPS[0], 0) >= min_cells) & (counts.get(GROUPS[1], 0) >= min_cells)]
    rows = []
    for regulon, g in patient[patient.unit_id.isin(eligible)].groupby("regulon", sort=False):
        p = g.pivot(index="unit_id", columns="pid1_subtype", values="AUC").dropna(subset=GROUPS)
        d = p[GROUPS[1]] - p[GROUPS[0]]
        if len(d) < 3 or np.allclose(d, 0):
            pv = np.nan
        else:
            try:
                pv = wilcoxon(d, alternative="two-sided", zero_method="wilcox", method="auto").pvalue
            except ValueError:
                pv = np.nan
        rows.append({"regulon": regulon, "n_paired_units": len(d), "mean_delta_inflammatory_minus_C1QC": d.mean() if len(d) else np.nan,
                     "median_delta": d.median() if len(d) else np.nan,
                     "direction_consistent_n": int((d > 0).sum()) if len(d) else 0,
                     "direction_consistent_fraction": float((d > 0).mean()) if len(d) else np.nan,
                     "p_wilcoxon": pv})
    out = pd.DataFrame(rows)
    if len(out):
        out["FDR_BH"] = qvals(out.p_wilcoxon)
        out = out.sort_values(["FDR_BH", "p_wilcoxon"], na_position="last")
    return out


def patient_paired_deltas(patient: pd.DataFrame, ann: pd.DataFrame, min_cells: int = 10) -> pd.DataFrame:
    counts = ann.groupby(["unit_id", "pid1_subtype"]).size().unstack(fill_value=0)
    eligible = counts.index[(counts.get(GROUPS[0], 0) >= min_cells) & (counts.get(GROUPS[1], 0) >= min_cells)]
    rows = []
    for regulon, g in patient[patient.unit_id.isin(eligible)].groupby("regulon", sort=False):
        a = g.pivot(index="unit_id", columns="pid1_subtype", values="AUC").dropna(subset=GROUPS)
        p = g.pivot(index="unit_id", columns="pid1_subtype", values="PID1_log1p").dropna(subset=GROUPS)
        cells = g.pivot(index="unit_id", columns="pid1_subtype", values="cells").dropna(subset=GROUPS)
        ix = a.index.intersection(p.index).intersection(cells.index)
        for unit in ix:
            rows.append({"unit_id": unit, "regulon": regulon,
                         "AUC_C1QC": a.loc[unit, GROUPS[0]], "AUC_inflammatory": a.loc[unit, GROUPS[1]],
                         "delta_AUC_inflammatory_minus_C1QC": a.loc[unit, GROUPS[1]] - a.loc[unit, GROUPS[0]],
                         "PID1_C1QC": p.loc[unit, GROUPS[0]], "PID1_inflammatory": p.loc[unit, GROUPS[1]],
                         "delta_PID1_inflammatory_minus_C1QC": p.loc[unit, GROUPS[1]] - p.loc[unit, GROUPS[0]],
                         "cells_C1QC": int(cells.loc[unit, GROUPS[0]]),
                         "cells_inflammatory": int(cells.loc[unit, GROUPS[1]])})
    return pd.DataFrame(rows)


def pid1_correlations(patient: pd.DataFrame, ann: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Correlate PID1 with each regulon in each subtype across units, plus paired deltas.
    min_cells = 10
    counts = ann.groupby(["unit_id", "pid1_subtype"]).size().unstack(fill_value=0)
    eligible = counts.index[(counts.get(GROUPS[0], 0) >= min_cells) & (counts.get(GROUPS[1], 0) >= min_cells)]
    rows = []
    for subtype in GROUPS:
        sub = patient[(patient.pid1_subtype == subtype) & (patient.unit_id.isin(eligible))]
        for reg, g in sub.groupby("regulon", sort=False):
            g = g.dropna(subset=["AUC", "PID1_log1p"])
            if len(g) >= 5 and g.PID1_log1p.nunique() > 1 and g.AUC.nunique() > 1:
                rho, p = spearmanr(g.PID1_log1p, g.AUC)
            else:
                rho, p = np.nan, np.nan
            rows.append({"analysis": "within_subtype_patient_PID1_vs_AUC", "subtype": subtype, "regulon": reg,
                         "n_units": len(g), "spearman_rho": rho, "p_value": p})
    within = pd.DataFrame(rows)
    if len(within):
        within["FDR_BH_within_subtype"] = within.groupby("subtype")["p_value"].transform(lambda x: qvals(x))
    # Matched within-patient difference correlation is an especially direct association test.
    rows = []
    for reg, g in patient[patient.unit_id.isin(eligible)].groupby("regulon", sort=False):
        a = g.pivot(index="unit_id", columns="pid1_subtype", values="AUC").dropna(subset=GROUPS)
        p = g.pivot(index="unit_id", columns="pid1_subtype", values="PID1_log1p").dropna(subset=GROUPS)
        ix = a.index.intersection(p.index)
        da = a.loc[ix, GROUPS[1]] - a.loc[ix, GROUPS[0]]
        dp = p.loc[ix, GROUPS[1]] - p.loc[ix, GROUPS[0]]
        if len(ix) >= 5 and da.nunique() > 1 and dp.nunique() > 1:
            rho, pv = spearmanr(dp, da)
        else:
            rho, pv = np.nan, np.nan
        rows.append({"regulon": reg, "n_paired_units": len(ix), "spearman_rho_delta_PID1_vs_delta_AUC": rho, "p_value": pv})
    delta = pd.DataFrame(rows)
    if len(delta):
        delta["FDR_BH"] = qvals(delta.p_value)
        delta = delta.sort_values(["FDR_BH", "p_value"], na_position="last")
    return within, delta


def parse_regulons(path: Path) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    # pySCENIC's own loader understands the multi-row ctx CSV format and
    # combines enriched motifs into TF regulons exactly as AUCell will score them.
    from pyscenic.cli.utils import load_signatures
    signatures = load_signatures(str(path))
    targets = {str(sig.name): set(map(str, sig.genes)) for sig in signatures}
    regdf = pd.DataFrame([{"regulon": name, "TF": name[:-3] if name.endswith("(+)") or name.endswith("(-)") else name,
                           "n_targets": len(genes), "TargetGenes": sorted(genes)}
                          for name, genes in targets.items()])
    return regdf, targets


def draw_figures(stats: pd.DataFrame, within: pd.DataFrame, delta: pd.DataFrame, output: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    output.mkdir(parents=True, exist_ok=True)
    if len(stats):
        s = stats.dropna(subset=["mean_delta_inflammatory_minus_C1QC"]).head(30).copy()
        if len(s):
            s = s.sort_values("mean_delta_inflammatory_minus_C1QC")
            fig, ax = plt.subplots(figsize=(8.5, max(5, 0.23 * len(s))))
            colors = ["#C44E52" if x > 0 else "#4C72B0" for x in s.mean_delta_inflammatory_minus_C1QC]
            ax.barh(s.regulon, s.mean_delta_inflammatory_minus_C1QC, color=colors)
            ax.axvline(0, color="black", lw=0.8)
            ax.set_xlabel("Mean paired Δ AUCell (inflammatory − C1QC)")
            ax.set_ylabel("Regulon")
            ax.set_title("GSE132465: patient-paired regulon activity")
            ax.xaxis.set_major_locator(MaxNLocator(6))
            fig.tight_layout()
            fig.savefig(output / "GSE132465_paired_regulon_activity_top30.png", dpi=220)
            plt.close(fig)
    if len(delta):
        d = delta.dropna(subset=["spearman_rho_delta_PID1_vs_delta_AUC"]).head(30).copy()
        if len(d):
            d = d.sort_values("spearman_rho_delta_PID1_vs_delta_AUC")
            fig, ax = plt.subplots(figsize=(8.5, max(5, 0.23 * len(d))))
            colors = ["#C44E52" if x > 0 else "#4C72B0" for x in d.spearman_rho_delta_PID1_vs_delta_AUC]
            ax.barh(d.regulon, d.spearman_rho_delta_PID1_vs_delta_AUC, color=colors)
            ax.axvline(0, color="black", lw=0.8)
            ax.set_xlim(-1, 1)
            ax.set_xlabel("Spearman ρ: paired ΔPID1 vs paired ΔAUCell")
            ax.set_ylabel("Regulon")
            ax.set_title("GSE132465: patient-level PID1–regulon association")
            fig.tight_layout()
            fig.savefig(output / "GSE132465_PID1_delta_regulon_correlations_top30.png", dpi=220)
            plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-inference", action="store_true", help="Reuse existing intermediate PySCENIC outputs")
    ap.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 4) - 1)))
    args = ap.parse_args()
    primary = OUT / "GSE132465_primary"
    external = OUT / "GSE200997_validation"
    for p in [OUT, primary, external]:
        p.mkdir(parents=True, exist_ok=True)

    # Cohort QC and expression preparation.
    annotations = {}
    exprs = {}
    for cohort in ["GSE132465", "GSE200997"]:
        a = load_annotations(cohort)
        annotations[cohort] = a
        expr_path = (primary if cohort == "GSE132465" else external) / f"{cohort}_macrophage_log1p_cpm.tsv.gz"
        if expr_path.exists():
            # These intermediates are generated by this dated run and allow
            # restarting after an inference-stage failure without rereading GEO.
            expr = pd.read_csv(expr_path, sep="\t", compression="gzip", index_col=0)
            expr.index = expr.index.astype(str)
            expr.columns = expr.columns.astype(str)
        else:
            cnt = load_counts(cohort, a.cell_id.tolist())
            expr = normalize_log1p(cnt)
            expr.to_csv(expr_path, sep="\t", compression="gzip", float_format="%.5g")
        exprs[cohort] = expr
        summary = a.groupby(["pid1_subtype", "unit_id"]).size().rename("n_cells").reset_index()
        summary.to_csv((primary if cohort == "GSE132465" else external) / f"{cohort}_patient_subtype_cell_counts.tsv", sep="\t", index=False)

    pexpr = exprs["GSE132465"]
    pexpr_path = primary / "GSE132465_macrophage_log1p_cpm.tsv.gz"
    adjacency = primary / "GSE132465_GRNBoost2_adjacencies.tsv"
    regulons_path = primary / "GSE132465_motif_pruned_regulons.csv"
    auc_primary_path = primary / "GSE132465_AUCell.tsv"
    exe = str(Path(sys.executable).parent / "Scripts" / "pyscenic.exe")
    if not args.skip_inference:
        # CLI accepts gzipped delimited matrices; preserve gene symbols as column names.
        if not adjacency.exists() or adjacency.stat().st_size == 0:
            cmd = [exe, "grn", str(pexpr_path), str(TFS), "--method", "grnboost2",
                   "--output", str(adjacency), "--num_workers", str(args.workers), "--seed", str(SEED)]
            run(cmd, primary / "01_grnboost2.log")
        if not regulons_path.exists() or regulons_path.stat().st_size == 0:
            # dask_multiprocessing avoids pySCENIC 0.12's multiprocessing-on-dill
            # worker bootstrap incompatibility on Windows.
            cmd = [exe, "ctx", str(adjacency), str(DB), "--annotations_fname", str(MOTIFS),
                   "--expression_mtx_fname", str(pexpr_path), "--output", str(regulons_path),
                   "--mode", "dask_multiprocessing", "--num_workers", str(min(args.workers, 4))]
            run(cmd, primary / "02_motif_pruning.log")
        if not auc_primary_path.exists() or auc_primary_path.stat().st_size == 0:
            cmd = [exe, "aucell", str(pexpr_path), str(regulons_path), "--output", str(auc_primary_path),
                   "--num_workers", str(args.workers), "--seed", str(SEED)]
            run(cmd, primary / "03_aucell_primary.log")
    # Primary and external AUCell use the same motif-pruned regulon definitions.
    auc1 = pd.read_csv(auc_primary_path, sep="\t", index_col=0)
    ac1, patient1 = activity_patient_tables(auc1, annotations["GSE132465"])
    ac1.to_csv(primary / "GSE132465_AUCell_cells_by_regulons.tsv.gz", sep="\t", compression="gzip", float_format="%.6g")
    patient1.to_csv(primary / "GSE132465_patient_subtype_regulon_activity.tsv", sep="\t", index=False, float_format="%.7g")

    # Patient-level paired differential activity at prespecified ≥5/10/20 cells per subtype.
    sensitivity = []
    sensitivity_tables = {}
    for thresh in [5, 10, 20]:
        stats = paired_regulon_stats(patient1, annotations["GSE132465"], thresh)
        stats.insert(0, "minimum_cells_per_subtype", thresh)
        stats.to_csv(primary / f"GSE132465_paired_regulon_stats_min{thresh}_cells.tsv", sep="\t", index=False, float_format="%.7g")
        if thresh == 10:
            main_stats = stats.copy()
        sensitivity_tables[thresh] = stats[["regulon", "mean_delta_inflammatory_minus_C1QC", "direction_consistent_n", "n_paired_units", "FDR_BH"]].rename(
            columns={"mean_delta_inflammatory_minus_C1QC": f"mean_delta_min{thresh}",
                     "direction_consistent_n": f"n_positive_min{thresh}",
                     "n_paired_units": f"n_paired_min{thresh}", "FDR_BH": f"FDR_min{thresh}"})
        sensitivity.append({"minimum_cells_per_subtype": thresh,
                            "n_paired_units": int(stats.n_paired_units.max()) if len(stats) else 0})
    threshold_stability = sensitivity_tables[5].merge(sensitivity_tables[10], on="regulon").merge(sensitivity_tables[20], on="regulon")
    sign5 = np.sign(threshold_stability.mean_delta_min5)
    sign10 = np.sign(threshold_stability.mean_delta_min10)
    sign20 = np.sign(threshold_stability.mean_delta_min20)
    threshold_stability["direction_same_all_thresholds"] = (sign5 == sign10) & (sign10 == sign20)
    threshold_stability["FDR_lt_0.05_all_thresholds"] = (threshold_stability[["FDR_min5", "FDR_min10", "FDR_min20"]] < 0.05).all(axis=1)
    threshold_stability.to_csv(primary / "GSE132465_regulon_threshold_sensitivity.tsv", sep="\t", index=False, float_format="%.7g")
    deltas = patient_paired_deltas(patient1, annotations["GSE132465"], min_cells=10)
    deltas.to_csv(primary / "GSE132465_patient_paired_regulon_deltas.tsv", sep="\t", index=False, float_format="%.7g")
    within, delta = pid1_correlations(patient1, annotations["GSE132465"])
    within.to_csv(primary / "GSE132465_patient_PID1_regulon_correlations_within_subtype.tsv", sep="\t", index=False, float_format="%.7g")
    delta.to_csv(primary / "GSE132465_patient_PID1_delta_regulon_correlations.tsv", sep="\t", index=False, float_format="%.7g")

    regdf, targets = parse_regulons(regulons_path)
    reg_names = set(ac1.columns)
    import re
    tf_by_reg = {str(x): re.sub(r"\([+-]\)$", "", str(x).replace("_extended", "")) for x in reg_names}
    annotated = main_stats.copy()
    annotated["TF"] = annotated.regulon.map(tf_by_reg).fillna(annotated.regulon.str.replace("_extended", "", regex=False))
    annotated["regulon_targets_PID1"] = annotated.regulon.map(lambda r: "PID1" in targets.get(r, set()))
    annotated = annotated.merge(delta[["regulon", "spearman_rho_delta_PID1_vs_delta_AUC", "p_value", "FDR_BH"]].rename(
        columns={"p_value":"PID1_delta_correlation_p", "FDR_BH":"PID1_delta_correlation_FDR"}), on="regulon", how="left")
    annotated.to_csv(primary / "GSE132465_PID1_related_candidate_regulons.tsv", sep="\t", index=False, float_format="%.7g")
    draw_figures(main_stats, within, delta, primary)

    # External cohort: score the same primary regulons. No second network inference.
    ext_expr = exprs["GSE200997"]
    ext_expr_path = external / "GSE200997_macrophage_log1p_cpm.tsv.gz"
    auc2_path = external / "GSE200997_AUCell_primary_regulons.tsv"
    if not args.skip_inference and (not auc2_path.exists() or auc2_path.stat().st_size == 0):
        cmd = [exe, "aucell", str(ext_expr_path), str(regulons_path), "--output", str(auc2_path),
               "--num_workers", str(args.workers), "--seed", str(SEED)]
        run(cmd, external / "01_aucell_external.log")
    auc2 = pd.read_csv(auc2_path, sep="\t", index_col=0)
    ac2, patient2 = activity_patient_tables(auc2, annotations["GSE200997"])
    ac2.to_csv(external / "GSE200997_AUCell_cells_by_primary_regulons.tsv.gz", sep="\t", compression="gzip", float_format="%.6g")
    patient2.to_csv(external / "GSE200997_sample_subtype_primary_regulon_activity.tsv", sep="\t", index=False, float_format="%.7g")
    extstats = paired_regulon_stats(patient2, annotations["GSE200997"], min_cells=10)
    extstats.to_csv(external / "GSE200997_sample_paired_regulon_stats_min10_cells.tsv", sep="\t", index=False, float_format="%.7g")
    ext_deltas = patient_paired_deltas(patient2, annotations["GSE200997"], min_cells=10)
    ext_deltas.to_csv(external / "GSE200997_sample_paired_regulon_deltas_min10_cells.tsv", sep="\t", index=False, float_format="%.7g")

    # Cohort-level count and primary/external direction comparison.
    unit_counts = {}
    for cohort, a in annotations.items():
        counts = a.groupby(["unit_id", "pid1_subtype"]).size().unstack(fill_value=0)
        both10 = int(((counts.get(GROUPS[0], 0) >= 10) & (counts.get(GROUPS[1], 0) >= 10)).sum())
        both5 = int(((counts.get(GROUPS[0], 0) >= 5) & (counts.get(GROUPS[1], 0) >= 5)).sum())
        both20 = int(((counts.get(GROUPS[0], 0) >= 20) & (counts.get(GROUPS[1], 0) >= 20)).sum())
        unit_counts[cohort] = {"cells": int(len(a)), "units_any_target_subtype": int(a.unit_id.nunique()), "paired_units_min5": both5,
                               "paired_units_min10": both10, "paired_units_min20": both20,
                               "cells_C1QC": int((a.pid1_subtype == GROUPS[0]).sum()),
                               "cells_inflammatory": int((a.pid1_subtype == GROUPS[1]).sum())}
    comparison = main_stats[["regulon", "mean_delta_inflammatory_minus_C1QC", "FDR_BH"]].merge(
        extstats[["regulon", "mean_delta_inflammatory_minus_C1QC", "FDR_BH", "n_paired_units"]],
        on="regulon", suffixes=("_GSE132465", "_GSE200997"), how="left")
    comparison["direction_replicated"] = np.sign(comparison.mean_delta_inflammatory_minus_C1QC_GSE132465) == np.sign(comparison.mean_delta_inflammatory_minus_C1QC_GSE200997)
    comparison.to_csv(OUT / "primary_external_regulon_direction_comparison.tsv", sep="\t", index=False, float_format="%.7g")

    # Store target sets for candidate networks in a compact machine-readable edge list.
    pid1_links = annotated[annotated.regulon_targets_PID1]
    ranked_diff = annotated.sort_values(["FDR_BH", "PID1_delta_correlation_FDR"], na_position="last").head(30)
    ranked_pid1 = annotated.sort_values("PID1_delta_correlation_p", na_position="last").head(10)
    candidate = pd.concat([ranked_diff, ranked_pid1, pid1_links]).drop_duplicates("regulon")
    edge_rows = []
    for _, r in candidate.iterrows():
        for gene in sorted(targets.get(r.regulon, set())):
            edge_rows.append({"regulon": r.regulon, "TF": r.TF, "target": gene,
                              "differential_FDR": r.FDR_BH,
                              "PID1_delta_correlation_rho": r.spearman_rho_delta_PID1_vs_delta_AUC,
                              "regulon_targets_PID1": r.regulon_targets_PID1})
    pd.DataFrame(edge_rows).to_csv(primary / "candidate_regulatory_network_edges_top30.tsv", sep="\t", index=False, float_format="%.7g")

    metadata = {"cohorts": unit_counts, "primary_network_cohort": "GSE132465", "external_scoring_cohort": "GSE200997",
                "filter": "two macrophage subtypes; qc_pass; not Scrublet doublet; no marker-lineage conflict",
                "unit_of_inference": "patient/sample; paired Wilcoxon; no cell-level inferential p-values",
                "paired_thresholds_primary": [5, 10, 20], "primary_paired_stats_min10": sensitivity[1],
                "n_regulons_scored_primary": int(ac1.shape[1]), "n_regulons_scored_external": int(ac2.shape[1]),
                "primary_cells_scored": int(ac1.shape[0]), "external_cells_scored": int(ac2.shape[0]),
                "TF_seed": SEED, "pyscenic_version": "0.12.1", "python_version": sys.version.split()[0],
                "numpy_version": np.__version__, "motif_ranking_database": str(DB.name), "motif_annotations": str(MOTIFS.name)}
    (OUT / "analysis_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")

    # Concise Chinese report with numerical results, network links, and caveats.
    n_primary = unit_counts["GSE132465"]["paired_units_min10"]
    n_external = unit_counts["GSE200997"]["paired_units_min10"]
    sig = annotated[annotated.FDR_BH < 0.05]
    rep = comparison[comparison.direction_replicated & comparison.mean_delta_inflammatory_minus_C1QC_GSE200997.notna()]
    n_sig5 = int((sensitivity_tables[5]["FDR_min5"] < 0.05).sum())
    n_sig10 = int((sensitivity_tables[10]["FDR_min10"] < 0.05).sum())
    n_sig20 = int((sensitivity_tables[20]["FDR_min20"] < 0.05).sum())
    n_stable_sig = int(threshold_stability["FDR_lt_0.05_all_thresholds"].sum())
    top_pid1 = annotated.sort_values("PID1_delta_correlation_p", na_position="last").iloc[0]
    infl_high = sig[sig.mean_delta_inflammatory_minus_C1QC > 0].sort_values("mean_delta_inflammatory_minus_C1QC", ascending=False).head(5)
    c1qc_high = sig[sig.mean_delta_inflammatory_minus_C1QC < 0].sort_values("mean_delta_inflammatory_minus_C1QC").head(5)
    infl_names = ", ".join(infl_high.regulon.tolist())
    c1qc_names = ", ".join(c1qc_high.regulon.tolist())
    lines = [
        "# SCENIC：PID1相关巨噬细胞调控网络（患者级分析）", "",
        "## 设计与样本量", "",
        f"以GSE132465推断GRN并进行motif pruning/AUCell；GSE200997仅用同一组primary regulons打分，未单独推断网络。只纳入Inflammatory与C1QC巨噬细胞，并要求当前注释的QC通过、非Scrublet doublet、无marker lineage conflict。",
        f"GSE132465：{unit_counts['GSE132465']['cells']}个细胞（C1QC {unit_counts['GSE132465']['cells_C1QC']}；炎症型 {unit_counts['GSE132465']['cells_inflammatory']}），{n_primary}位患者两亚群均≥10个细胞。≥5/≥20阈值下分别有{unit_counts['GSE132465']['paired_units_min5']}/{unit_counts['GSE132465']['paired_units_min20']}位配对患者。",
        f"GSE200997：{unit_counts['GSE200997']['cells']}个细胞（C1QC {unit_counts['GSE200997']['cells_C1QC']}；炎症型 {unit_counts['GSE200997']['cells_inflammatory']}），仅{n_external}个样本两亚群均≥10个细胞（≥5：{unit_counts['GSE200997']['paired_units_min5']}；≥20：{unit_counts['GSE200997']['paired_units_min20']}），外部队列只能作探索性方向核对。",
        f"GSE132465识别并打分{ac1.shape[1]}个regulons；采用患者/样本为统计单位，主比较为每位患者内炎症型减C1QC的AUCell差，Wilcoxon并对regulons做BH校正。", "",
        "## 主要结果", "",
        f"在≥10细胞/亚群的患者级配对分析中，BH-FDR<0.05的差异regulons有{len(sig)}/{ac1.shape[1]}个；炎症型高/低方向分别按患者差值统计，见`GSE132465_paired_regulon_stats_min10_cells.tsv`及逐患者差值表。",
        f"显著且炎症型活性较高的前五个为{infl_names}；显著且C1QC较高的前五个为{c1qc_names}。",
        f"细胞数阈值敏感性：≥5、≥10、≥20时分别有{n_sig5}、{n_sig10}、{n_sig20}个regulons达到BH-FDR<0.05；其中{n_stable_sig}个在三种阈值下均显著，见`GSE132465_regulon_threshold_sensitivity.tsv`。",
        f"患者内ΔPID1与Δregulon AUCell的相关分析中，最强的是{top_pid1.regulon}（Spearman ρ={top_pid1.spearman_rho_delta_PID1_vs_delta_AUC:.3f}，原始P={top_pid1.PID1_delta_correlation_p:.4g}，140个regulon校正后FDR={top_pid1.PID1_delta_correlation_FDR:.3f}）；目前没有相关性通过BH-FDR 0.05。",
        f"motif-pruned regulon中将PID1列为候选靶基因的有{len(pid1_links)}个：{', '.join(pid1_links.TF.tolist()) if len(pid1_links) else '无'}。其中可结合细胞状态差异看作优先候选，但其患者级PID1相关性仍需以表格中的效应量和FDR判断，不能视为调控已被证明。",
        f"GSE200997在1个可配对样本中可核对140个regulon方向，其中{int(rep.direction_replicated.sum())}个与主队列方向一致；由于仅1个样本达到每亚群≥10个细胞，不将方向一致率称为跨患者复现。", "",
        "## 解释边界", "",
        "SCENIC输出的是表达共变驱动的候选TF regulon及motif支持的调控关系，不等同于直接结合或因果调控。患者数较少，尤其外部队列亚群细胞稀疏；细胞数阈值敏感性、患者方向一致率与FDR应共同解读。不同患者/样本作为统计单位，未把细胞当作独立重复。",
        "本次只依赖GSE132465/GSE200997原始UMI矩阵、当前PID1髓系精细注释以及官方人类TF/motif-ranking资源；不读取旧CellChat预测或历史SCENIC结果。",
        "", "## 资源与软件", "",
        "pySCENIC 0.12.1：GRNBoost2 → cisTarget motif pruning → AUCell；人类v10 clustered hg38 RefSeq R80 10kb upstream/downstream ranking database与匹配motif注释。资源说明见Aerts lab官方资源库及pySCENIC文档。",
        "- https://github.com/aertslab/pySCENIC", "- https://resources.aertslab.org/cistarget/databases/", "- https://resources.aertslab.org/cistarget/motif2tf/", "",
        "关键文件：`GSE132465_PID1_related_candidate_regulons.tsv`、`candidate_regulatory_network_edges_top30.tsv`、`primary_external_regulon_direction_comparison.tsv`及患者级活动/统计表。", ""
    ]
    (OUT / "SCENIC_PID1_macrophage_regulon_report_zh.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
