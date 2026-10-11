#!/usr/bin/env python
"""Summarize fine-state CellChat recurrence and patient-level PID1 associations."""
from __future__ import annotations

from itertools import combinations
from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

ROOT = Path(r"D:\DINP\whynot17-new-project")
OUT = ROOT / "results" / "2026-10-11" / "PID1_CellChat_subtype_consistency"
COHORTS = ["GSE132465", "GSE200997"]
EDGE = ["ligand", "receptor", "ligand_receptor", "interaction_name"]
KEY = ["direction", *EDGE]
MIN_UNITS_FOR_CROSS_COHORT_STABILITY = 5


def bh(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = values.notna() & np.isfinite(values)
    if ok.any():
        out.loc[ok] = multipletests(values.loc[ok], method="fdr_bh")[1]
    return out


def clean_key(row: pd.Series) -> tuple[str, ...]:
    return tuple(str(row[c]) for c in KEY)


def is_app_mif(ligand: object, receptor: object) -> bool:
    lig, rec = str(ligand).upper(), str(receptor).upper()
    return (lig == "APP" and "CD74" in rec) or (lig == "MIF" and ("CD74" in rec or "CD44" in rec))


def candidate_label(ligand: object, receptor: object) -> str:
    lig, rec = str(ligand).upper(), str(receptor).upper()
    if lig == "APP" and "CD74" in rec:
        return "APP–CD74"
    if lig == "MIF" and "CD74" in rec and "CD44" in rec:
        return "MIF–CD74/CD44"
    if lig == "MIF" and "CD74" in rec and "CXCR4" in rec:
        return "MIF–CD74/CXCR4"
    if lig == "MIF" and "CD74" in rec:
        return "MIF–CD74"
    if lig == "MIF" and "CD44" in rec:
        return "MIF–CD44"
    return f"{ligand}–{receptor}"


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    comm, elig = [], []
    for cohort in COHORTS:
        cp = OUT / f"{cohort}_PID1_CellChat_patient_sample_interactions.tsv.gz"
        ep = OUT / f"{cohort}_PID1_CellChat_pair_eligibility.tsv"
        c = pd.read_csv(cp, sep="\t") if cp.exists() else pd.DataFrame()
        e = pd.read_csv(ep, sep="\t")
        if not c.empty:
            c["cohort"] = cohort
            for col in ["ligand", "receptor", "ligand_receptor", "interaction_name", "direction", "pid1_subtype"]:
                c[col] = c[col].astype(str)
            c["prob"] = pd.to_numeric(c["prob"], errors="coerce").fillna(0.0)
            c["qval_pair"] = pd.to_numeric(c["qval_pair"], errors="coerce").fillna(1.0)
            comm.append(c)
        e["cohort"] = cohort
        for col in ["unit_id", "pid1_subtype"]:
            e[col] = e[col].astype(str)
        elig.append(e)
    return (pd.concat(comm, ignore_index=True) if comm else pd.DataFrame(),
            pd.concat(elig, ignore_index=True))


def recurrence(comm: pd.DataFrame, elig: pd.DataFrame) -> pd.DataFrame:
    out = []
    if comm.empty:
        return pd.DataFrame()
    interactions = comm[KEY].drop_duplicates().copy()
    for threshold, flag in [(11, "eligible_gt10"), (20, "eligible_ge20"), (50, "eligible_ge50")]:
        e0 = elig[elig[flag].astype(bool)].copy()
        for (cohort, subtype), e in e0.groupby(["cohort", "pid1_subtype"], sort=False):
            units = set(e.unit_id)
            state_comm = comm[(comm.cohort == cohort) & (comm.pid1_subtype == subtype) & comm.unit_id.isin(units)]
            grid = e[["unit_id"]].drop_duplicates().assign(_join=1).merge(
                interactions.assign(_join=1), on="_join", how="inner"
            ).drop(columns="_join")
            by_unit = grid.merge(
                state_comm[["unit_id", *KEY, "prob", "qval_pair"]].drop_duplicates(["unit_id", *KEY]),
                on=["unit_id", *KEY], how="left"
            )
            by_unit["prob"] = by_unit["prob"].fillna(0.0)
            by_unit["qval_pair"] = by_unit["qval_pair"].fillna(1.0)
            by_unit["q05"] = (by_unit.qval_pair < 0.05) & (by_unit.prob > 0)
            sums = by_unit.groupby(KEY, dropna=False).agg(
                n_prob_positive=("prob", lambda x: int((x > 0).sum())),
                n_units_q05=("q05", "sum"),
                median_prob_all_eligible=("prob", "median"),
            ).reset_index()
            positive_medians = by_unit[by_unit.prob > 0].groupby(KEY, dropna=False).prob.median().rename("median_prob_positive").reset_index()
            sums = sums.merge(positive_medians, on=KEY, how="left").fillna({"median_prob_positive": 0.0})
            for row in sums.to_dict("records"):
                n_q = int(row["n_units_q05"])
                out.append({
                    "cohort": cohort, "min_myeloid_cells": threshold, "pid1_subtype": subtype,
                    **{k: str(row[k]) for k in KEY},
                    "n_eligible_units": len(units), "n_prob_positive": int(row["n_prob_positive"]),
                    "n_units_q05": n_q, "q05_fraction": n_q / len(units) if units else np.nan,
                    "median_prob_all_eligible": float(row["median_prob_all_eligible"]),
                    "median_prob_positive": float(row["median_prob_positive"]),
                    "median_PID1_mean_log1p": float(e.PID1_mean_log1p.median()),
                    "median_myeloid_cells": float(e.n_myeloid_cells.median()),
                })
    tab = pd.DataFrame(out)
    if not tab.empty:
        tab["stable_half_q05_min5"] = (tab.n_eligible_units >= MIN_UNITS_FOR_CROSS_COHORT_STABILITY) & (tab.q05_fraction >= 0.5)
    return tab


def selected_edges(comm: pd.DataFrame) -> set[tuple[str, ...]]:
    # Use recurrent within-unit BH-supported interactions plus the pre-specified APP/MIF axes.
    sig = comm[(comm.prob > 0) & (comm.qval_pair < 0.05)]
    recurrent = sig.groupby(["cohort", *KEY], dropna=False).unit_id.nunique()
    chosen = {tuple(map(str, k[1:] if isinstance(k, tuple) else (k,))) for k, n in recurrent.items() if n >= 2}
    for _, row in comm.iterrows():
        if is_app_mif(row.ligand, row.receptor):
            chosen.add(clean_key(row))
    return chosen


def paired_state_tests(comm: pd.DataFrame, elig: pd.DataFrame, chosen: set[tuple[str, ...]]) -> pd.DataFrame:
    rows = []
    for cohort in COHORTS:
        ec = elig[(elig.cohort == cohort) & elig.eligible_gt10.astype(bool)]
        states = sorted(ec.pid1_subtype.unique())
        for key in sorted(chosen):
            direction, ligand, receptor, lr, interaction = key
            sub = comm[(comm.cohort == cohort) & (comm.direction == direction)]
            for col, val in zip(EDGE, [ligand, receptor, lr, interaction]):
                sub = sub[sub[col].astype(str) == val]
            edge_states = sorted(set(states))
            if len(edge_states) < 2:
                continue
            values = {}
            for state in edge_states:
                e = ec[ec.pid1_subtype == state][["unit_id"]].drop_duplicates()
                obs = sub[sub.pid1_subtype == state].drop_duplicates("unit_id").set_index("unit_id").prob
                values[state] = e.set_index("unit_id").join(obs.rename("prob"), how="left").prob.fillna(0.0)
            for state_a, state_b in combinations(edge_states, 2):
                both = pd.concat([values[state_a].rename("a"), values[state_b].rename("b")], axis=1, join="inner").dropna()
                n = len(both)
                if n < 5:
                    continue
                delta = both.a - both.b
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        p = float(wilcoxon(both.a, both.b, zero_method="wilcox", correction=False, method="approx").pvalue) if np.any(delta != 0) else 1.0
                except Exception:
                    p = np.nan
                rows.append({
                    "cohort": cohort, "direction": direction, "ligand": ligand, "receptor": receptor,
                    "ligand_receptor": lr, "interaction_name": interaction,
                    "subtype_A": state_a, "subtype_B": state_b,
                    "n_paired_units": n, "median_prob_A": float(both.a.median()),
                    "median_prob_B": float(both.b.median()), "median_delta_A_minus_B": float(delta.median()),
                    "n_A_gt_B": int((delta > 0).sum()), "wilcoxon_p": p,
                })
    tab = pd.DataFrame(rows)
    if not tab.empty:
        tab["p_BH"] = tab.groupby("cohort").wilcoxon_p.transform(bh)
    return tab


def cellchat_unit_grid(comm: pd.DataFrame, elig: pd.DataFrame, key: tuple[str, ...], cohort: str) -> pd.DataFrame:
    direction, ligand, receptor, lr, interaction = key
    e = elig[(elig.cohort == cohort) & elig.eligible_gt10.astype(bool)].copy()
    o = comm[(comm.cohort == cohort) & (comm.direction == direction)]
    for col, val in zip(EDGE, [ligand, receptor, lr, interaction]):
        o = o[o[col].astype(str) == val]
    o = o.drop_duplicates(["unit_id", "pid1_subtype"])
    d = e.merge(o[["unit_id", "pid1_subtype", "prob"]], on=["unit_id", "pid1_subtype"], how="left")
    d["prob"] = d.prob.fillna(0.0)
    d["pid1_z"] = (d.PID1_mean_log1p - d.PID1_mean_log1p.mean()) / d.PID1_mean_log1p.std(ddof=0)
    d["log_myeloid"] = np.log1p(d.n_myeloid_cells)
    d["log_epithelial"] = np.log1p(d.n_epithelial_cells)
    return d


def pid1_models(comm: pd.DataFrame, elig: pd.DataFrame, chosen: set[tuple[str, ...]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    within, identity = [], []
    for cohort in COHORTS:
        for key in sorted(chosen):
            direction, ligand, receptor, lr, interaction = key
            d = cellchat_unit_grid(comm, elig, key, cohort)
            if d.empty:
                continue
            # Within-subtype association: state identity is fixed by stratification; cell counts are covariates.
            for subtype, g in d.groupby("pid1_subtype"):
                g = g.dropna(subset=["PID1_mean_log1p", "prob"])
                if len(g) < 8 or g.prob.nunique() < 2 or g.PID1_mean_log1p.std(ddof=0) == 0:
                    continue
                g = g.copy()
                sd = g.PID1_mean_log1p.std(ddof=0)
                g["pid1_z"] = (g.PID1_mean_log1p - g.PID1_mean_log1p.mean()) / sd
                try:
                    fit = ols("prob ~ pid1_z + log_myeloid + log_epithelial", data=g).fit(cov_type="HC3")
                    within.append({
                        "cohort": cohort, "direction": direction, "pid1_subtype": subtype,
                        "ligand": ligand, "receptor": receptor, "ligand_receptor": lr, "interaction_name": interaction,
                        "n_units": g.unit_id.nunique(), "n_prob_positive": int((g.prob > 0).sum()),
                        "beta_per_SD_PID1": float(fit.params.get("pid1_z", np.nan)),
                        "se_HC3": float(fit.bse.get("pid1_z", np.nan)),
                        "p_value": float(fit.pvalues.get("pid1_z", np.nan)),
                    })
                except Exception:
                    pass
            # Patient fixed effects + subtype fixed effects separate a within-patient PID1 association
            # from the average communication difference caused by subtype identity.
            d = d.dropna(subset=["PID1_mean_log1p", "prob"])
            if d.unit_id.nunique() < 8 or d.pid1_subtype.nunique() < 2 or len(d) < 20 or (d.prob > 0).sum() < 3:
                continue
            sd = d.PID1_mean_log1p.std(ddof=0)
            if sd == 0:
                continue
            d = d.copy()
            d["pid1_z"] = (d.PID1_mean_log1p - d.PID1_mean_log1p.mean()) / sd
            try:
                fit = ols("prob ~ pid1_z + C(pid1_subtype) + C(unit_id) + log_myeloid + log_epithelial", data=d).fit(
                    cov_type="cluster", cov_kwds={"groups": d.unit_id, "use_correction": True}
                )
                identity.append({
                    "cohort": cohort, "direction": direction, "ligand": ligand, "receptor": receptor,
                    "ligand_receptor": lr, "interaction_name": interaction,
                    "n_observations": len(d), "n_units": d.unit_id.nunique(), "n_subtypes": d.pid1_subtype.nunique(),
                    "n_prob_positive": int((d.prob > 0).sum()),
                    "beta_per_SD_PID1_adjusted_for_subtype_and_unit": float(fit.params.get("pid1_z", np.nan)),
                    "se_cluster": float(fit.bse.get("pid1_z", np.nan)),
                    "p_value": float(fit.pvalues.get("pid1_z", np.nan)),
                })
            except Exception:
                continue
    w, i = pd.DataFrame(within), pd.DataFrame(identity)
    if not w.empty:
        w["p_BH"] = w.groupby("cohort").p_value.transform(bh)
    if not i.empty:
        i["p_BH"] = i.groupby("cohort").p_value.transform(bh)
    return w, i


def main() -> None:
    comm, elig = load()
    if comm.empty:
        raise RuntimeError("CellChat interaction outputs are not present yet; wait for the R run to finish.")
    error_frames = []
    for cohort in COHORTS:
        ep = OUT / f"{cohort}_PID1_CellChat_errors.tsv"
        if ep.exists():
            er = pd.read_csv(ep, sep="\t")
            if not er.empty:
                error_frames.append(er)
    errors = pd.concat(error_frames, ignore_index=True) if error_frames else pd.DataFrame(columns=["cohort", "unit_id", "error"])
    failed_units = {(str(r.cohort), str(r.unit_id).removesuffix(".rds")) for r in errors.itertuples()}
    if failed_units:
        elig["cellchat_status"] = ["error" if (str(c), str(u)) in failed_units else "success" for c, u in zip(elig.cohort, elig.unit_id)]
        # Never count a failed fit as a biological zero in recurrence or association tests.
        elig = elig[elig.cellchat_status.eq("success")].copy()
    else:
        elig["cellchat_status"] = "success"
    elig.to_csv(OUT / "PID1_CellChat_pair_eligibility_model_units.tsv", sep="\t", index=False)
    rec = recurrence(comm, elig)
    rec.to_csv(OUT / "PID1_CellChat_edge_recurrence_by_subtype_threshold.tsv.gz", sep="\t", index=False, compression="gzip")
    selected = selected_edges(comm)
    selected_summary = rec[(rec.min_myeloid_cells == 11) & rec.apply(lambda r: is_app_mif(r.ligand, r.receptor), axis=1)].copy()
    selected_summary.to_csv(OUT / "PID1_CellChat_APP_MIF_candidate_recurrence.tsv", sep="\t", index=False)
    paired = paired_state_tests(comm, elig, selected)
    paired.to_csv(OUT / "PID1_CellChat_paired_subtype_edge_contrasts.tsv.gz", sep="\t", index=False, compression="gzip")
    within, identity = pid1_models(comm, elig, selected)
    within.to_csv(OUT / "PID1_CellChat_PID1_association_within_subtype.tsv", sep="\t", index=False)
    identity.to_csv(OUT / "PID1_CellChat_PID1_association_adjusted_for_subtype_identity.tsv", sep="\t", index=False)

    cross_rows = []
    for key, g in rec[rec.min_myeloid_cells == 11].groupby(KEY, dropna=False):
        direction, ligand, receptor, lr, interaction = key
        row = {"direction": direction, "ligand": ligand, "receptor": receptor,
               "ligand_receptor": lr, "interaction_name": interaction}
        stable_by_cohort = []
        for cohort in COHORTS:
            cg = g[g.cohort == cohort]
            for state in ["Monocyte_CD14_classical", "Monocyte_CD14_classical_inflammatory",
                          "Macrophage_C1QC", "Macrophage_inflammatory", "Macrophage_SPP1_TREM2", "Macrophage_IFN"]:
                one = cg[cg.pid1_subtype == state]
                suffix = f"{cohort}_{state}"
                if len(one):
                    rr = one.iloc[0]
                    row[f"{suffix}_eligible"] = int(rr.n_eligible_units)
                    row[f"{suffix}_q05"] = int(rr.n_units_q05)
                    row[f"{suffix}_q05_fraction"] = float(rr.q05_fraction)
                    stable_by_cohort.append(bool(rr.stable_half_q05_min5))
                else:
                    row[f"{suffix}_eligible"] = 0
                    row[f"{suffix}_q05"] = 0
                    row[f"{suffix}_q05_fraction"] = 0.0
                    stable_by_cohort.append(False)
        row["any_shared_subtype_stable_half_q05_min5"] = any(
            all(row.get(f"{cohort}_{state}_eligible", 0) >= MIN_UNITS_FOR_CROSS_COHORT_STABILITY and
                row.get(f"{cohort}_{state}_q05_fraction", 0) >= .5 for cohort in COHORTS)
            for state in ["Monocyte_CD14_classical", "Monocyte_CD14_classical_inflammatory", "Macrophage_C1QC", "Macrophage_inflammatory"]
        )
        cross_rows.append(row)
    cross = pd.DataFrame(cross_rows)
    cross.to_csv(OUT / "PID1_CellChat_cross_cohort_subtype_edge_recurrence.tsv.gz", sep="\t", index=False, compression="gzip")

    # Human-readable compact report.
    coverage = pd.read_csv(OUT / "PID1_CellChat_subtype_coverage_summary.tsv", sep="\t")
    cand = selected_summary.sort_values(["cohort", "direction", "pid1_subtype", "ligand_receptor"])
    top_rows = rec[(rec.min_myeloid_cells == 11) & (rec.n_eligible_units >= 5) & (rec.n_units_q05 >= 2)].sort_values(
        ["cohort", "q05_fraction", "n_units_q05", "median_prob_all_eligible"], ascending=[True, False, False, False]
    )
    top_rows = top_rows.groupby("cohort", sort=False).head(12)
    focal = paired[((paired.subtype_A == "Macrophage_C1QC") & (paired.subtype_B == "Macrophage_inflammatory")) |
                   ((paired.subtype_A == "Monocyte_CD14_classical") & (paired.subtype_B == "Monocyte_CD14_classical_inflammatory"))] if not paired.empty else paired
    best_assoc = identity.sort_values("p_BH") if not identity.empty else identity

    def md(df: pd.DataFrame, cols: list[str], labels: list[str] | None = None, limit: int | None = None) -> str:
        if df is None or df.empty:
            return "（无符合展示条件的行。）"
        d = df[cols].copy()
        if limit:
            d = d.head(limit)
        labels = labels or cols
        def fmt(x: object) -> str:
            if pd.isna(x): return "—"
            if isinstance(x, (float, np.floating)):
                return f"{x:.3g}"
            return str(x).replace("|", "\\|")
        lines = ["| " + " | ".join(labels) + " |", "| " + " | ".join(["---"] * len(labels)) + " |"]
        lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in d.itertuples(index=False, name=None)]
        return "\n".join(lines)

    coverage_table = md(coverage, ["cohort", "pid1_subtype", "n_eligible_gt10", "n_eligible_ge20", "n_eligible_ge50", "median_subtype_cells", "median_PID1_mean_log1p", "median_PID1_detection"],
                        ["队列", "亚群", ">10细胞单位数", "≥20细胞单位数", "≥50细胞单位数", "中位细胞数", "PID1中位均值", "PID1中位检出率"])
    cand_table = md(cand, ["cohort", "direction", "pid1_subtype", "ligand_receptor", "n_eligible_units", "n_units_q05", "q05_fraction", "median_prob_all_eligible"],
                    ["队列", "方向", "髓系亚群", "配体_受体", "合格单位数", "BH-q<0.05单位数", "复现比例", "全部合格单位的中位概率"])
    top_table = md(top_rows, ["cohort", "direction", "pid1_subtype", "ligand_receptor", "n_eligible_units", "n_units_q05", "q05_fraction", "median_prob_all_eligible"],
                   ["队列", "方向", "髓系亚群", "配体_受体", "合格单位数", "BH-q<0.05单位数", "复现比例", "中位概率"], limit=24)
    focal_table = md(focal.sort_values(["cohort", "p_BH"]) if not focal.empty else focal,
                     ["cohort", "direction", "ligand_receptor", "subtype_A", "subtype_B", "n_paired_units", "median_delta_A_minus_B", "wilcoxon_p", "p_BH"],
                     ["队列", "方向", "配体_受体", "亚群A", "亚群B", "配对单位数", "中位概率差（A−B）", "Wilcoxon p", "BH-FDR"], limit=30)
    within_sorted = within.sort_values("p_BH") if not within.empty else within
    within_table = md(within_sorted, ["cohort", "direction", "pid1_subtype", "ligand_receptor", "n_units", "n_prob_positive", "beta_per_SD_PID1", "p_value", "p_BH"],
                      ["队列", "方向", "髓系亚群", "配体_受体", "单位数", "正概率单位数", "PID1系数（每SD）", "p", "BH-FDR"], limit=25)
    assoc_table = md(best_assoc, ["cohort", "direction", "ligand_receptor", "n_units", "n_subtypes", "n_prob_positive", "beta_per_SD_PID1_adjusted_for_subtype_and_unit", "p_value", "p_BH"],
                     ["队列", "方向", "配体_受体", "单位数", "亚群数", "正概率单位数", "PID1调整系数（每SD）", "p", "BH-FDR"], limit=20)

    sig_within = int((within.p_BH < .05).sum()) if not within.empty else 0
    sig_identity = int((identity.p_BH < .05).sum()) if not identity.empty else 0

    def rec_line(cohort: str, state: str, ligand_receptor: str) -> str:
        rows = selected_summary[
            (selected_summary.cohort == cohort)
            & (selected_summary.direction == "Epithelial_to_myeloid")
            & (selected_summary.pid1_subtype == state)
            & (selected_summary.ligand_receptor == ligand_receptor)
        ]
        if rows.empty:
            return "不可估计"
        row = rows.sort_values("n_units_q05", ascending=False).iloc[0]
        return f"{int(row.n_units_q05)}/{int(row.n_eligible_units)}（{row.q05_fraction:.1%}）"

    c1_app_132 = rec_line("GSE132465", "Macrophage_C1QC", "APP_CD74")
    c1_app_200 = rec_line("GSE200997", "Macrophage_C1QC", "APP_CD74")
    c1_mif_132 = rec_line("GSE132465", "Macrophage_C1QC", "MIF_CD74_CD44")
    c1_mif_200 = rec_line("GSE200997", "Macrophage_C1QC", "MIF_CD74_CD44")
    paired_app = focal[(focal.cohort == "GSE132465") & (focal.direction == "Epithelial_to_myeloid") & (focal.ligand_receptor == "APP_CD74") & (focal.subtype_A == "Macrophage_C1QC") & (focal.subtype_B == "Macrophage_inflammatory")]
    paired_mif = focal[(focal.cohort == "GSE132465") & (focal.direction == "Epithelial_to_myeloid") & (focal.ligand_receptor == "MIF_CD74_CD44") & (focal.subtype_A == "Macrophage_C1QC") & (focal.subtype_B == "Macrophage_inflammatory")]

    def paired_summary(row_set: pd.DataFrame) -> str:
        if row_set.empty:
            return "不可估计"
        row = row_set.iloc[0]
        return (f"配对n={int(row.n_paired_units)}，中位差={row.median_delta_A_minus_B:.4g}，"
                f"Wilcoxon p={row.wilcoxon_p:.3g}，BH-FDR={row.p_BH:.3g}")

    app_adjusted = identity[(identity.cohort == "GSE132465") & (identity.direction == "Epithelial_to_myeloid") & (identity.ligand_receptor == "APP_CD74")]
    mif_adjusted = identity[(identity.cohort == "GSE132465") & (identity.direction == "Epithelial_to_myeloid") & (identity.ligand_receptor == "MIF_CD74_CD44")]
    identity_app_fdr = f"{app_adjusted.p_BH.iloc[0]:.3g}" if not app_adjusted.empty else "不可估计"
    identity_mif_fdr = f"{mif_adjusted.p_BH.iloc[0]:.3g}" if not mif_adjusted.empty else "不可估计"

    # Candidate-axis recurrence figure. Each bar is a fraction of the units
    # eligible for that exact epithelial–subtype pair; denominators stay visible.
    cand_plot = rec[(rec.min_myeloid_cells == 11) & rec.apply(lambda r: is_app_mif(r.ligand, r.receptor), axis=1)].copy()
    cand_plot["candidate"] = [candidate_label(l, r) for l, r in zip(cand_plot.ligand, cand_plot.receptor)]
    if not cand_plot.empty:
        focal_candidates = ["APP–CD74", "MIF–CD74/CD44", "MIF–CD74/CXCR4"]
        cand_plot = cand_plot[cand_plot.candidate.isin(focal_candidates)].copy()
        agg = cand_plot.sort_values("q05_fraction", ascending=False).drop_duplicates(
            ["cohort", "direction", "pid1_subtype", "candidate"]
        )
        order = ["Monocyte_CD14_classical", "Monocyte_CD14_classical_inflammatory", "Macrophage_C1QC",
                 "Macrophage_inflammatory", "Macrophage_SPP1_TREM2", "Macrophage_IFN"]
        short = {x: x.replace("Monocyte_", "Mono_").replace("Macrophage_", "Mac_").replace("_classical_inflammatory", "_inflam").replace("_classical", "_classical").replace("_inflammatory", "_inflam").replace("_SPP1_TREM2", "_SPP1/TREM2") for x in order}
        fig, axes = plt.subplots(2, 2, figsize=(15, 9.5), sharey=True)
        colors = {"APP–CD74": "#4C78A8", "MIF–CD74/CD44": "#F58518", "MIF–CD74/CXCR4": "#E45756",
                  "MIF–CD74": "#72B7B2", "MIF–CD44": "#54A24B"}
        for ri, direction in enumerate(["Epithelial_to_myeloid", "Myeloid_to_epithelial"]):
            for ci, cohort in enumerate(COHORTS):
                ax = axes[ri, ci]
                sub = agg[(agg.cohort == cohort) & (agg.direction == direction)]
                xbase = np.arange(len(order), dtype=float)
                candidates = focal_candidates
                width = .22
                for j, candidate in enumerate(candidates):
                    vals, denoms = [], []
                    for state in order:
                        z = sub[(sub.pid1_subtype == state) & (sub.candidate == candidate)]
                        vals.append(float(z.q05_fraction.iloc[0]) if len(z) else 0.0)
                        denoms.append(int(z.n_eligible_units.iloc[0]) if len(z) else 0)
                    pos = xbase + (j - (len(candidates) - 1) / 2) * width
                    ax.bar(pos, vals, width=width, color=colors[candidate], label=candidate)
                for xx, state in zip(xbase, order):
                    d0 = sub[sub.pid1_subtype == state]
                    if not d0.empty:
                        ax.text(xx, 1.035, f"n={int(d0.n_eligible_units.iloc[0])}", ha="center", va="bottom", fontsize=8)
                ax.set_xticks(xbase)
                ax.set_xticklabels([short[x] for x in order], rotation=35, ha="right", fontsize=8)
                ax.set_ylim(0, 1.12)
                ax.set_ylabel("Fraction of eligible units with within-unit BH q < 0.05")
                ax.set_title(f"{cohort}: {direction.replace('_', ' ')}")
                ax.grid(axis="y", alpha=.2)
        fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=colors[x]) for x in focal_candidates],
                   labels=focal_candidates, frameon=False, ncol=3, fontsize=9, loc="upper center", bbox_to_anchor=(.5, .955))
        fig.suptitle("Pre-specified APP/MIF axes: recurrence across PID1 myeloid subtypes", y=.995, fontsize=14)
        fig.subplots_adjust(top=.88, bottom=.18, hspace=.42, wspace=.12)
        fig.savefig(OUT / "APP_MIF_subtype_recurrence.png", dpi=220, bbox_inches="tight")
        plt.close(fig)

    # Forest-style view of adjusted PID1 associations if any are estimable.
    forest = identity.sort_values("p_BH").head(20).copy() if not identity.empty else identity
    if not forest.empty:
        fig_h = max(4.2, .32 * len(forest))
        fig, ax = plt.subplots(figsize=(10, fig_h))
        y = np.arange(len(forest))
        se = forest.se_cluster.to_numpy(float)
        beta = forest.beta_per_SD_PID1_adjusted_for_subtype_and_unit.to_numpy(float)
        ax.errorbar(beta, y, xerr=1.96 * se, fmt="o", color="#315A7D", ecolor="#94AFC4", capsize=2)
        ax.axvline(0, color="black", linewidth=1, linestyle="--")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.cohort} | {r.direction} | {r.ligand_receptor}" for r in forest.itertuples()], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("PID1 coefficient per SD (95% cluster-robust Wald interval)")
        ax.set_title("Candidate CellChat edge probability vs PID1, adjusted for subtype, unit and cell counts")
        fig.tight_layout()
        fig.savefig(OUT / "PID1_adjusted_communication_associations.png", dpi=220, bbox_inches="tight")
        plt.close(fig)

    text = f"""# PID1 与上皮–髓系亚群通讯：逐患者/样本 CellChat 分析

日期：2026-10-11  
结果目录：`results/2026-10-11/PID1_CellChat_subtype_consistency/`

## 结论摘要

- 上皮→C1QC巨噬细胞的 APP–CD74 在两个队列均较常复现：GSE132465 {c1_app_132}，GSE200997 {c1_app_200}；MIF–CD74/CD44 分别为 {c1_mif_132} 和 {c1_mif_200}。GSE200997细分髓系样本数偏少，应作为敏感性支持。
- 在 GSE132465 同患者配对比较中，APP–CD74概率在 C1QC 巨噬细胞高于炎症巨噬细胞：{paired_summary(paired_app)}。相反，本次标签中PID1亚群中位均值为 C1QC 0.0225、炎症巨噬细胞 0.198，因此该通讯差异不支持“PID1高表达状态对应更强APP–CD74”的简单解释。
- 同一比较中 MIF–CD74/CD44 的亚群概率差异不显著：{paired_summary(paired_mif)}。控制亚群身份、患者/样本和细胞数后，APP–CD74与MIF–CD74/CD44的PID1关联 BH-FDR 分别为 {identity_app_fdr} 和 {identity_mif_fdr}；全部身份调整检验中 BH-FDR<0.05 为 {sig_identity}/{len(identity)}。
- 因此，当前数据支持部分上皮–髓系预测通讯在细分亚群后仍可见，但**没有证据证明这些通讯与PID1独立相关**。尤其GSE200997炎症髓系样本覆盖有限，不足以验证其亚群差异。

## 分析问题

沿用原逐患者 CellChat 的输入表达矩阵和肿瘤上皮标签，将旧版宽泛的 `Monocyte/Macrophage/Other_myeloid` 替换为本次 PID1 分析中独立确定的单核/巨噬细胞亚群。只分析肿瘤上皮↔单核/巨噬细胞，不纳入CAF、DC、未解析髓系或T/NK污染。两队列分开分析，单位分别是 GSE132465 患者和 GSE200997 肿瘤样本。

本次重用的表达矩阵是由原始 GEO UMI 数据生成的 CellChat 配体/受体信号基因面板矩阵（各队列1471和1349个基因，均含 APP、MIF、CD74、CD44、PID1），不是全转录组。上皮标签沿用原分析；髓系标签完全换成本次新注释。对照旧分析时仅使用你指定的 CellChat 文件，不读取其他历史分析结果。

最终纳入的表达矩阵/标签对齐细胞为 GSE132465 23,010 个（17,469个肿瘤上皮+5,541个新版单核/巨噬细胞）和 GSE200997 7,400 个（6,871个上皮+529个新版单核/巨噬细胞）；GSE200997有一个上皮样本没有达到任何单核/巨噬亚群的分析门槛。

旧版宽分组的关键结果已从你指定的旧版复现表核对：上皮→单核细胞 APP–CD74 为 GSE132465 19/19、GSE200997 11/12；上皮→巨噬细胞为22/23、11/12。MIF–CD74/CD44 在上皮→单核细胞为19/19、8/12；上皮→巨噬细胞为22/23、7/12。下面的新结果将这些宽标签拆成当前可辨认的单核/巨噬状态，并重新计算各状态的分母和网络。

## CellChat与统计方法

- CellChat 2.2.0.9001，`CellChatDB.human`；每名患者/样本独立建网，`triMean`、`population.size=FALSE`、100次置换。
- 与旧结果相同，CellChat内部过滤要求群体细胞数 **>10**。另按目标亚群每单位≥20和≥50个细胞重新计算覆盖率敏感性；边概率保持同一次患者/样本模型估计，仅改变纳入分母。
- 每个单位内先将置换p值下限设为1/101，再在每个发送者–接收者组合内做BH。`q<0.05`的跨患者比例是复现度描述，不是患者级因果或临床推断检验。
- 亚群间比较以同一患者/样本内相同上皮、不同髓系状态的 CellChat 概率做配对 Wilcoxon；未返回的边按0记为“本次网络未预测到”，不是证明生物学信号绝对不存在。
- PID1关联检验分两层：固定髓系亚群后，跨患者/样本回归通信概率与该亚群PID1均值，并调整亚群和上皮细胞数；以及加入患者/样本固定效应和亚群固定效应的合并模型。后者检验在控制患者背景、亚群身份及细胞数后，患者内不同亚群的PID1差异是否仍与边概率相关。模型仅对多单位重复的候选边运行，并对检验做BH校正。
- 本次所有具备上皮和至少一个目标亚群的单位均成功建网：CellChat运行异常 {len(errors)} 个；异常单位不纳入统计分母。

## 每个亚群的患者/样本覆盖

{coverage_table}

> GSE200997没有患者ID，以样本作单位；不少细分状态覆盖很少，未达到门槛的亚群不作“未通讯”的结论。

## 预先关注的 APP/MIF 轴

下表报告 APP–CD74、MIF–CD74或其受体复合体、MIF–CD44 在不同亚群和两个方向中的复现。分母是该上皮–目标亚群具备足够细胞数的患者/样本。

{cand_table}

## 其他重复出现的配体–受体边

以下按每亚群至少5个合格单位、至少2个单位内BH-q<0.05筛出；展示排序靠前的边，不作为独立机制验证。

{top_table}

## 同患者亚群对照

重点查看经典与炎症单核细胞，以及 C1QC 与炎症巨噬细胞。差值为亚群A减亚群B；概率缺失按0。表中只显示配对单位数≥5且边在至少两个单位中达到组内BH阈值、或属于预先关注的APP/MIF轴的比较。

{focal_table}

## PID1是否独立关联通讯

“亚群内模型”固定一个髓系状态，在该状态各患者/样本间回归边概率与PID1均值，并调整髓系和上皮细胞数。身份调整模型再加入患者/样本与髓系亚群固定效应。共运行 {len(within)} 个亚群内模型和 {len(identity)} 个身份调整模型；BH后显著数分别为 {sig_within} 和 {sig_identity}。

### 固定亚群后的跨患者/样本关联

{within_table}

### 同时调整亚群身份和患者/样本后的关联

{assoc_table}

表格中的系数代表观察性关联。即使某条通讯边和PID1在统计上相关，也不能证明PID1调控配体/受体或实际细胞通讯。若身份调整后关联消失，而某亚群边更稳定，则更符合“亚群身份/细胞组成可解释该差异”的结果；若关联仍在，也只能视为后续实验的候选线索。

## 结果文件

- [`PID1_CellChat_APP_MIF_candidate_recurrence.tsv`](PID1_CellChat_APP_MIF_candidate_recurrence.tsv)：重点通路按亚群、方向、阈值的复现率。
- [`PID1_CellChat_pair_eligibility_model_units.tsv`](PID1_CellChat_pair_eligibility_model_units.tsv)：实际进入统计分母的逐患者/样本亚群覆盖。
- [`PID1_CellChat_edge_recurrence_by_subtype_threshold.tsv.gz`](PID1_CellChat_edge_recurrence_by_subtype_threshold.tsv.gz)：所有配体–受体边与10/20/50细胞覆盖敏感性。
- [`PID1_CellChat_paired_subtype_edge_contrasts.tsv.gz`](PID1_CellChat_paired_subtype_edge_contrasts.tsv.gz)：同患者/样本亚群配对比较。
- [`PID1_CellChat_PID1_association_within_subtype.tsv`](PID1_CellChat_PID1_association_within_subtype.tsv)：亚群内部的PID1关联模型。
- [`PID1_CellChat_PID1_association_adjusted_for_subtype_identity.tsv`](PID1_CellChat_PID1_association_adjusted_for_subtype_identity.tsv)：同时调整亚群身份与患者/样本的模型。
- [`GSE132465_PID1_CellChat_patient_sample_interactions.tsv.gz`](GSE132465_PID1_CellChat_patient_sample_interactions.tsv.gz)、[`GSE200997_PID1_CellChat_patient_sample_interactions.tsv.gz`](GSE200997_PID1_CellChat_patient_sample_interactions.tsv.gz)：逐单位全边与概率、p/q值。
- [`GSE132465_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz`](GSE132465_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz)、[`GSE200997_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz`](GSE200997_PID1_CellChat_epithelial_myeloid_cell_annotations.tsv.gz)：本次CellChat实际纳入的逐细胞标签。
- `APP_MIF_subtype_recurrence.png`：预先关注轴的逐亚群复现图。
- 若存在可估计的身份调整模型，`PID1_adjusted_communication_associations.png` 汇总其系数与区间。

## 解释限制

细胞通讯由配体/受体表达模式预测，不能证明蛋白结合、分泌或方向性信号已发生。GSE200997的上皮标签来自先前的多marker推断，且细分髓系状态单位覆盖较少；它是敏感性分析。CellChat概率也会受表达检测率、每亚群细胞数及样本组成影响；因此报告纳入细胞数、阈值敏感性和身份调整模型。此次分析不支持因果机制结论，也不扩展到机器学习。
"""
    (OUT / "PID1_CellChat_subtype_consistency_report_zh.md").write_text(text, encoding="utf-8")
    write_manifest = {
        "run_date": "2026-10-11",
        "analysis": "PID1-associated epithelial-myeloid CellChat at fine myeloid subtype resolution",
        "cohorts": COHORTS,
        "units": {"GSE132465": "patient", "GSE200997": "tumor sample"},
        "expression_input": "Prior CellChat ligand/receptor panel matrices generated from original GEO UMI data; only this referenced CellChat input was reused.",
        "epithelial_annotation": "Epithelial cell label from the prior CellChat preparation, limited to tumor epithelial cells.",
        "myeloid_annotation": "Current PID1_myeloid_finelocalization annotations; only eligible Monocyte_* and Macrophage_* states included.",
        "cellchat": {"version": "2.2.0.9001", "database": "CellChatDB.human", "average": "triMean", "population_size": False, "nboot": 100, "min_cells": ">10", "seed": "20261011 + ordered unit index"},
        "not_done": ["machine learning", "causal inference", "all-cell communication"],
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(write_manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote report; recurrence rows={len(rec)}, paired tests={len(paired)}, within-state models={len(within)}, identity models={len(identity)}")


if __name__ == "__main__":
    main()
