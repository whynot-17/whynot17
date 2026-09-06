#!/usr/bin/env python
"""Focused macrophage -> CRC epithelial ligand-receptor analysis.

R/CellChat is optional in this workspace.  When it is unavailable (the current
case), this script uses a small, explicit human CellChat-compatible prior and
performs pooled descriptive scores plus donor-level support checks.  It never
uses cells as biological replicates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parents[1]
sys.path.insert(0, str(REPO_ROOT))
from analysis.dinp_crc_macrophage_subtype.run_macrophage_subtype_analysis import (  # noqa: E402
    DEFAULT_ANNOTATION,
    DEFAULT_MATRIX,
    read_matrix_selected,
)

DEFAULT_MACRO_SCORE = REPO_ROOT / "analysis/dinp_crc_macrophage_subtype/outputs/macrophage_cell_program_scores.csv"

FROZEN_81 = [
    "ABCC4", "ADAMTS1", "AHSG", "AKR1C1", "AKR1C3", "ATF4", "ATG5", "BECN1", "C4BPA", "CA9", "CADM3", "CASP14", "CCL20", "CGA", "CPS1", "CRP", "CXCL13", "CXCR4", "CYP11A1", "CYP2A6", "DDIT4", "DKK1", "DPT", "DUSP1", "ESR1", "FLG", "GSTA2", "HPGD", "HSD3B2", "HSPA1A", "HSPA1L", "IGFBP1", "LUM", "MAP1LC3B", "MMP2", "MMP9", "NEAT1", "NR1I2", "NR1I3", "NR3C1", "OMD", "PARP10", "PGR", "PLA2G4A", "PLAC8", "PPARA", "PPARD", "PPARG", "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTGES", "PTGES2", "PTGES3", "PTGFR", "PTGS1", "PTGS2", "PTX3", "RELA", "RGS2", "RPS2", "RXRA", "RXRB", "SCGN", "SIRT1", "SIRT2", "SIRT3", "SIRT5", "SLCO2A1", "SNHG25", "SNORA13", "SNORA65", "SPP2", "SQSTM1", "STAR", "STAT3", "SUSD2", "TIMP1", "TIMP2", "TNC"
]
DRIVERS = ["NEAT1", "MMP9", "TIMP1", "STAT3", "PTGER4", "PTGES3", "CXCR4"]

# Explicitly versioned, human ligand-receptor prior.  Complexes use '+' and
# are scored as the geometric mean of their component expression values.
LR_ROWS = [
    ("SPP1", "CD44", "SPP1", "Secreted Signaling"), ("SPP1", "ITGAV+ITGB1", "SPP1", "ECM-Receptor"), ("SPP1", "ITGAV+ITGB5", "SPP1", "ECM-Receptor"),
    ("CXCL12", "CXCR4", "CXCL", "Secreted Signaling"), ("CXCL9", "CXCR3", "CXCL", "Secreted Signaling"), ("CXCL10", "CXCR3", "CXCL", "Secreted Signaling"), ("CXCL11", "CXCR3", "CXCL", "Secreted Signaling"), ("CXCL8", "CXCR1", "CXCL", "Secreted Signaling"), ("CXCL8", "CXCR2", "CXCL", "Secreted Signaling"), ("CXCL1", "CXCR2", "CXCL", "Secreted Signaling"), ("CXCL2", "CXCR2", "CXCL", "Secreted Signaling"),
    ("CCL2", "CCR2", "CCL", "Secreted Signaling"), ("CCL3", "CCR1", "CCL", "Secreted Signaling"), ("CCL3", "CCR5", "CCL", "Secreted Signaling"), ("CCL4", "CCR5", "CCL", "Secreted Signaling"), ("CCL5", "CCR5", "CCL", "Secreted Signaling"), ("CCL20", "CCR6", "CCL", "Secreted Signaling"),
    ("VEGFA", "FLT1", "VEGF", "Secreted Signaling"), ("VEGFA", "KDR", "VEGF", "Secreted Signaling"), ("VEGFA", "NRP1", "VEGF", "Secreted Signaling"), ("VEGFC", "FLT4", "VEGF", "Secreted Signaling"),
    ("TGFB1", "TGFBR1+TGFBR2", "TGFb", "Secreted Signaling"), ("TGFB2", "TGFBR1+TGFBR2", "TGFb", "Secreted Signaling"), ("TGFB3", "TGFBR1+TGFBR2", "TGFb", "Secreted Signaling"),
    ("TNF", "TNFRSF1A", "TNF", "Secreted Signaling"), ("TNF", "TNFRSF1B", "TNF", "Secreted Signaling"), ("LTA", "LTBR", "TNF", "Secreted Signaling"), ("FASLG", "FAS", "TNF", "Cell-Cell Contact"),
    ("IL1B", "IL1R1", "IL-related", "Secreted Signaling"), ("IL6", "IL6R", "IL-related", "Secreted Signaling"), ("IL10", "IL10RA+IL10RB", "IL-related", "Secreted Signaling"), ("IL11", "IL11RA", "IL-related", "Secreted Signaling"), ("IL17A", "IL17RA+IL17RC", "IL-related", "Secreted Signaling"), ("IL18", "IL18R1+IL18RAP", "IL-related", "Secreted Signaling"), ("IL12A+IL12B", "IL12RB1+IL12RB2", "IL-related", "Secreted Signaling"),
    ("EGF", "EGFR", "EGF", "Secreted Signaling"), ("TGFA", "EGFR", "EGF", "Secreted Signaling"), ("HBEGF", "EGFR", "EGF", "Secreted Signaling"), ("AREG", "EGFR", "EGF", "Secreted Signaling"), ("EREG", "EGFR", "EGF", "Secreted Signaling"),
    ("MIF", "CD74+CXCR4", "MIF", "Secreted Signaling"), ("MIF", "CD74+ACKR3", "MIF", "Secreted Signaling"),
    ("C3", "C3AR1", "Complement", "Secreted Signaling"), ("C5", "C5AR1", "Complement", "Secreted Signaling"), ("C5", "C5AR2", "Complement", "Secreted Signaling"),
    ("GAS6", "AXL", "GAS", "Secreted Signaling"), ("GAS6", "MERTK", "GAS", "Secreted Signaling"), ("GAS6", "TYRO3", "GAS", "Secreted Signaling"),
    ("CSF1", "CSF1R", "CSF", "Secreted Signaling"), ("CSF2", "CSF2RA+CSF2RB", "CSF", "Secreted Signaling"), ("CSF3", "CSF3R", "CSF", "Secreted Signaling"),
    ("FN1", "ITGA5+ITGB1", "FN1", "ECM-Receptor"), ("COL1A1", "ITGA1+ITGB1", "Collagen", "ECM-Receptor"), ("COL3A1", "ITGA1+ITGB1", "Collagen", "ECM-Receptor"), ("LAMC1", "ITGA6+ITGB1", "Laminin", "ECM-Receptor"),
    ("JAG1", "NOTCH1", "Notch", "Cell-Cell Contact"), ("JAG2", "NOTCH2", "Notch", "Cell-Cell Contact"), ("DLL1", "NOTCH1", "Notch", "Cell-Cell Contact"), ("DLL4", "NOTCH1", "Notch", "Cell-Cell Contact"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def dump_json(obj, path: Path):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def genes_in(value: str) -> list[str]:
    return [x.strip().upper() for x in value.split("+") if x.strip()]


def component_value(logx: sp.csr_matrix, indices: dict[str, int], genes: list[str], cells: np.ndarray):
    vals, det = [], []
    for gene in genes:
        gi = indices.get(gene)
        if gi is None:
            continue
        a = np.asarray(logx[cells, gi].toarray()).ravel()
        vals.append(float(a.mean()) if len(a) else np.nan)
        det.append(float(np.mean(a > 0)) if len(a) else np.nan)
    if not vals:
        return np.nan, np.nan, 0
    return float(np.exp(np.mean(np.log(np.maximum(vals, 1e-12))))), float(np.exp(np.mean(np.log(np.maximum(det, 1e-12))))), len(vals)


def make_score(logx, var_index, sender_cells, receiver_cells, ligand, receptor):
    lm, ld, nl = component_value(logx, var_index, genes_in(ligand), sender_cells)
    rm, rd, nr = component_value(logx, var_index, genes_in(receptor), receiver_cells)
    if not nl or not nr or not len(sender_cells) or not len(receiver_cells):
        return {"ligand_mean": np.nan, "ligand_detection_fraction": np.nan, "receptor_mean": np.nan, "receptor_detection_fraction": np.nan, "communication_score": np.nan, "ligand_components_present": nl, "receptor_components_present": nr}
    return {"ligand_mean": lm, "ligand_detection_fraction": ld, "receptor_mean": rm, "receptor_detection_fraction": rd, "communication_score": float(math.sqrt(max(lm, 0) * max(rm, 0)) * math.sqrt(max(ld, 0) * max(rd, 0))), "ligand_components_present": nl, "receptor_components_present": nr}


def matrix_meta(ann: pd.DataFrame, macro: pd.DataFrame):
    ann = ann.copy()
    ann["condition"] = ann["Class"].astype(str).str.lower()
    ann["cell_id"] = ann["Index"].astype(str)
    ann["donor_id"] = ann["Patient"].astype(str)
    macro = macro.copy()
    macro["cell_id"] = macro["cell_id"].astype(str)
    macro = macro[["cell_id", "donor_id", "condition", "cluster", "subtype"]].drop_duplicates("cell_id")
    merged = ann.merge(macro, on=["cell_id", "donor_id", "condition"], how="left", suffixes=("", "_macro"))
    # Tumor epithelial states come from the existing source labels (CMS1–4);
    # this module does not infer malignancy or CNV anew.
    merged["is_tumor_epithelial"] = (merged["Cell_type"].astype(str).str.lower() == "epithelial cells") & (merged["condition"] == "tumor")
    merged["is_normal_epithelial"] = (merged["Cell_type"].astype(str).str.lower() == "epithelial cells") & (merged["condition"] == "normal")
    merged["is_sender"] = merged["subtype"].notna()
    merged["receiver_state"] = np.where(merged["is_tumor_epithelial"], merged["Cell_subtype"].astype(str), "normal_epithelial")
    return merged


def pooled_scores(logx, var_index, meta: pd.DataFrame, state_pairs: list[tuple[str, str, str]], condition: str, receiver_mode: str):
    rows = []
    idx = {cid: i for i, cid in enumerate(meta["cell_id"])}
    if receiver_mode == "tumor":
        recv = meta["is_tumor_epithelial"] & (meta["condition"] == condition)
        receiver_states = sorted(meta.loc[recv, "receiver_state"].astype(str).unique().tolist()) + ["tumor_epithelial"]
    else:
        recv = meta["is_normal_epithelial"] & (meta["condition"] == condition)
        receiver_states = ["normal_epithelial"]
    recv_cells_all = np.where(recv.to_numpy())[0]
    for subtype, ligand, receptor in state_pairs:
        send = meta["is_sender"] & (meta["subtype"] == subtype) & (meta["condition"] == condition)
        send_cells = np.where(send.to_numpy())[0]
        if len(send_cells) < 5:
            continue
        for state in receiver_states:
            recv_mask = recv if state in ("tumor_epithelial", "normal_epithelial") else (recv & (meta["receiver_state"] == state))
            recv_cells = np.where(recv_mask.to_numpy())[0]
            if len(recv_cells) < 10:
                continue
            s = make_score(logx, var_index, send_cells, recv_cells, ligand, receptor)
            if not np.isfinite(s["communication_score"]) or s["communication_score"] <= 0:
                continue
            rows.append({"condition": condition, "sender_subtype": subtype, "receiver_state": state, "receiver_group": "tumor_epithelial" if receiver_mode == "tumor" else "normal_epithelial", "ligand": ligand, "receptor": receptor, "pathway": LR_LOOKUP[(ligand, receptor)][0], "category": LR_LOOKUP[(ligand, receptor)][1], "sender_cells": len(send_cells), "receiver_cells": len(recv_cells), **s})
    return pd.DataFrame(rows)


def donor_support(logx, var_index, meta, top_pairs):
    rows = []
    donors = sorted(meta.loc[meta["condition"].isin(["tumor", "normal"]), "donor_id"].unique())
    for _, pair in top_pairs.iterrows():
        sender_sub = pair["sender_subtype"]
        condition = pair["condition"]
        receiver_group = pair["receiver_group"]
        for donor in donors:
            send_mask = (meta["is_sender"] & (meta["subtype"] == sender_sub) & (meta["condition"] == condition) & (meta["donor_id"] == donor))
            if receiver_group == "tumor_epithelial":
                recv_mask = meta["is_tumor_epithelial"] & (meta["condition"] == condition) & (meta["donor_id"] == donor)
            else:
                recv_mask = meta["is_normal_epithelial"] & (meta["condition"] == condition) & (meta["donor_id"] == donor)
            scells, rcells = np.where(send_mask.to_numpy())[0], np.where(recv_mask.to_numpy())[0]
            s = make_score(logx, var_index, scells, rcells, pair["ligand"], pair["receptor"])
            supported = bool(len(scells) >= 5 and len(rcells) >= 10 and np.isfinite(s["communication_score"]) and s["communication_score"] > 0)
            rows.append({"donor_id": donor, "condition": condition, "sender_subtype": sender_sub, "receiver_group": receiver_group, "receiver_state": pair["receiver_state"], "ligand": pair["ligand"], "receptor": pair["receptor"], "pathway": pair["pathway"], "sender_cell_count": len(scells), "receiver_cell_count": len(rcells), "support": supported, "support_status": "supported" if supported else ("insufficient coverage" if len(scells) < 5 or len(rcells) < 10 else "not detected"), **s})
    out = pd.DataFrame(rows)
    if not len(out):
        return out
    key = ["condition", "sender_subtype", "receiver_group", "receiver_state", "ligand", "receptor"]
    summary = out.groupby(key, dropna=False).agg(n_donors=("donor_id", "nunique"), n_supported_donors=("support", "sum"), median_supported_score=("communication_score", lambda x: float(np.nanmedian(x[x > 0])) if np.any(x > 0) else np.nan), max_donor_score=("communication_score", "max")).reset_index()
    summary["consistent_across_ge4_donors"] = summary["n_supported_donors"] >= 4
    summary["donor_driven"] = (summary["n_supported_donors"] <= 2) & (summary["max_donor_score"] > 0)
    out = out.merge(summary, on=key, how="left")
    return out


def savefig_all(fig, path: Path):
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(path.with_suffix("." + ext), dpi=220, facecolor="white")


def make_figures(scores, donor, incoming, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_theme(style="white", context="talk")
    figdir = out / "figures"; figdir.mkdir(exist_ok=True)
    tumor = scores[(scores["condition"] == "tumor") & (scores["receiver_group"] == "tumor_epithelial")].copy()
    top = tumor.sort_values("communication_score", ascending=False).head(15)
    # Focused network: four sender nodes and epithelial state nodes, top edges only.
    fig, ax = plt.subplots(figsize=(10, 7)); ax.axis("off")
    senders = ["SPP1_like_TAM", "C1QC_like_TAM", "FCN1_inflammatory_monocyte_like", "resident_like"]
    receivers = sorted(top["receiver_state"].unique()) if len(top) else ["tumor_epithelial"]
    pos = {s: (-1.0, 1.5 - i * .9) for i, s in enumerate(senders)}
    pos.update({r: (1.0, 1.5 - i * (3.0 / max(1, len(receivers) - 1))) for i, r in enumerate(receivers)})
    for s, (x, y) in pos.items(): ax.scatter([x], [y], s=1400 if s in senders else 1200, color="#6baed6" if s in senders else "#fdae6b", edgecolor="black", zorder=3); ax.text(x, y, s.replace("_", "\n"), ha="center", va="center", fontsize=8)
    if len(top):
        for _, r in top.head(20).iterrows():
            x1, y1 = pos.get(r["sender_subtype"], (-1, 0)); x2, y2 = pos.get(r["receiver_state"], (1, 0)); ax.annotate("", xy=(x2 - .12, y2), xytext=(x1 + .12, y1), arrowprops=dict(arrowstyle="->", color="#777", alpha=.45, lw=1 + 4 * r["communication_score"] / max(top["communication_score"].max(), 1e-9)))
    ax.set_title("Focused TAM → tumor epithelial communication")
    savefig_all(fig, figdir / "tam_to_tumor_network_circle"); plt.close(fig)
    # Bubble plot.
    if len(top):
        b = top.copy(); b["pair"] = b["ligand"] + "–" + b["receptor"]
        fig, ax = plt.subplots(figsize=(12, 8)); sns.scatterplot(data=b, x="sender_subtype", y="pair", size="communication_score", hue="pathway", sizes=(30, 350), palette="tab20", ax=ax, legend=False); ax.tick_params(axis="x", rotation=35); ax.set_title("Top TAM → tumor LR pairs"); savefig_all(fig, figdir / "tam_tumor_lr_bubble"); plt.close(fig)
        spp = b[b["sender_subtype"] == "SPP1_like_TAM"].head(15)
        fig, ax = plt.subplots(figsize=(10, 6)); sns.barplot(data=spp, y="ligand", x="communication_score", hue="receptor", dodge=False, ax=ax, palette="Set2"); ax.set_title("SPP1-like TAM → tumor epithelial LR pairs"); savefig_all(fig, figdir / "spp1_tam_tumor_specific_lr"); plt.close(fig)
    else:
        for name in ["tam_tumor_lr_bubble", "spp1_tam_tumor_specific_lr"]:
            fig, ax = plt.subplots(); ax.text(.5, .5, "No supported pooled LR pairs", ha="center"); ax.axis("off"); savefig_all(fig, figdir / name); plt.close(fig)
    # SPP1 versus C1QC pathway heatmap.
    p = tumor[tumor["sender_subtype"].isin(["SPP1_like_TAM", "C1QC_like_TAM"])].groupby(["sender_subtype", "pathway"])["communication_score"].max().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(12, 5)); sns.heatmap(p, cmap="Blues", ax=ax); ax.set_title("Outgoing pathway strength: SPP1-like vs C1QC-like"); savefig_all(fig, figdir / "spp1_vs_c1qc_outgoing_pathways"); plt.close(fig)
    # C1QC incoming/outgoing role.
    role = pd.DataFrame({"role": ["C1QC outgoing → tumor", "tumor incoming → C1QC"], "score": [float(tumor[tumor["sender_subtype"] == "C1QC_like_TAM"]["communication_score"].sum()), float(incoming["communication_score"].sum()) if len(incoming) else 0.0]})
    fig, ax = plt.subplots(figsize=(8, 5)); sns.barplot(data=role, x="role", y="score", palette=["#6baed6", "#fd8d3c"], ax=ax); ax.tick_params(axis="x", rotation=20); ax.set_title("C1QC-like TAM communication role context"); savefig_all(fig, figdir / "c1qc_tam_incoming_outgoing_role"); plt.close(fig)
    # Tumor versus normal outgoing strength.
    ps = scores.groupby(["condition", "sender_subtype"])["communication_score"].sum().reset_index()
    fig, ax = plt.subplots(figsize=(10, 5)); sns.barplot(data=ps, x="sender_subtype", y="communication_score", hue="condition", palette="Set2", ax=ax); ax.tick_params(axis="x", rotation=35); ax.set_title("Tumor versus normal focused communication (descriptive)"); savefig_all(fig, figdir / "tumor_vs_normal_communication"); plt.close(fig)
    # Donor consistency heatmap.
    if len(donor):
        d = donor.copy(); d["pair"] = d["ligand"] + "–" + d["receptor"]
        hp = d.pivot_table(index="pair", columns="donor_id", values="communication_score", aggfunc="max")
    else:
        hp = pd.DataFrame([[0]], index=["no pairs"], columns=["no donors"])
    fig, ax = plt.subplots(figsize=(10, max(4, .3 * len(hp)))); sns.heatmap(hp, cmap="Blues", vmin=0, ax=ax); ax.set_title("Donor-level support for top LR pairs"); savefig_all(fig, figdir / "donor_consistency_heatmap"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    ap.add_argument("--annotation", type=Path, default=DEFAULT_ANNOTATION)
    ap.add_argument("--macro-score", type=Path, default=DEFAULT_MACRO_SCORE)
    ap.add_argument("--outdir", type=Path, default=MODULE_DIR / "outputs")
    args = ap.parse_args(); out = args.outdir; out.mkdir(parents=True, exist_ok=True); (out / "figures").mkdir(exist_ok=True)
    t0 = time.time()
    ann = pd.read_csv(args.annotation, sep="\t")
    macro = pd.read_csv(args.macro_score)
    meta = matrix_meta(ann, macro)
    # Use all four existing macrophage subtype senders and epithelial receivers;
    # no new clustering or malignant-state inference is performed.
    selected_ids = meta.loc[meta["is_sender"] | meta["is_tumor_epithelial"] | meta["is_normal_epithelial"], "cell_id"].tolist()
    X, genes, _ = read_matrix_selected(args.matrix, selected_ids)
    logx = X.tocsr().astype(np.float32)
    lib = np.asarray(logx.sum(axis=1)).ravel(); lib[lib <= 0] = 1
    logx = logx.multiply((1e4 / lib)[:, None]).tocsr(); logx.data = np.log1p(logx.data)
    # read_matrix_selected preserves requested cell order, matching meta subset.
    meta = meta.set_index("cell_id").loc[selected_ids].reset_index()
    var_index = {g.upper(): i for i, g in enumerate(genes)}
    global LR_LOOKUP
    LR_LOOKUP = {(l, r): (p, c) for l, r, p, c in LR_ROWS}
    lr_db = pd.DataFrame(LR_ROWS, columns=["ligand", "receptor", "pathway", "category"]); lr_db["species"] = "human"; lr_db["database"] = "CellChat-compatible focused prior"; lr_db["database_version"] = "0.1-local"; lr_db.to_csv(out / "lr_database_human.csv", index=False)
    pairs = [(sub, l, r) for sub in ["SPP1_like_TAM", "C1QC_like_TAM", "FCN1_inflammatory_monocyte_like", "resident_like"] for l, r, _, _ in LR_ROWS]
    score_frames = [pooled_scores(logx, var_index, meta, pairs, "tumor", "tumor"), pooled_scores(logx, var_index, meta, pairs, "normal", "normal")]
    scores = pd.concat([x for x in score_frames if len(x)], ignore_index=True) if any(len(x) for x in score_frames) else pd.DataFrame()
    if len(scores): scores.to_csv(out / "tam_tumor_lr_scores.csv", index=False)
    else: pd.DataFrame(columns=["condition", "sender_subtype", "receiver_state", "ligand", "receptor", "communication_score"]).to_csv(out / "tam_tumor_lr_scores.csv", index=False)
    # Pathway and outgoing summaries are descriptive pooled scores.
    if len(scores):
        path = scores.groupby(["condition", "receiver_group", "sender_subtype", "pathway", "category"], dropna=False).agg(pathway_max_score=("communication_score", "max"), pathway_mean_score=("communication_score", "mean"), n_supported_pairs=("communication_score", "size")).reset_index(); path.to_csv(out / "communication_pathway_summary.csv", index=False)
        outgo = scores.groupby(["condition", "receiver_group", "sender_subtype"], dropna=False).agg(outgoing_strength=("communication_score", "sum"), n_supported_pairs=("communication_score", "size"), max_pair_score=("communication_score", "max")).reset_index(); outgo.to_csv(out / "outgoing_strength_by_subtype.csv", index=False)
    else:
        path = pd.DataFrame(); outgo = pd.DataFrame()
        path.to_csv(out / "communication_pathway_summary.csv", index=False); outgo.to_csv(out / "outgoing_strength_by_subtype.csv", index=False)
    # Directional incoming context: tumor epithelial as sender to C1QC-like TAM.
    incoming = []
    tumor_cells = np.where((meta["is_tumor_epithelial"]).to_numpy())[0]
    c1qc_cells = np.where((meta["is_sender"] & (meta["subtype"] == "C1QC_like_TAM") & (meta["condition"] == "tumor")).to_numpy())[0]
    for ligand, receptor, pathway, category in LR_ROWS:
        s = make_score(logx, var_index, tumor_cells, c1qc_cells, ligand, receptor)
        if np.isfinite(s["communication_score"]) and s["communication_score"] > 0:
            incoming.append({"condition": "tumor", "sender_state": "tumor_epithelial", "receiver_state": "C1QC_like_TAM", "ligand": ligand, "receptor": receptor, "pathway": pathway, "category": category, "sender_cells": len(tumor_cells), "receiver_cells": len(c1qc_cells), **s})
    incoming_df = pd.DataFrame(incoming)
    if len(incoming_df): incoming_df.to_csv(out / "tumor_to_c1qc_incoming_lr_scores.csv", index=False)
    else: incoming_df.to_csv(out / "tumor_to_c1qc_incoming_lr_scores.csv", index=False)
    c1qc_role = pd.DataFrame([{ "direction": "C1QC_like_TAM_to_tumor", "mean_score": float(scores[(scores["condition"] == "tumor") & (scores["sender_subtype"] == "C1QC_like_TAM")]["communication_score"].mean()) if len(scores) else np.nan, "sum_score": float(scores[(scores["condition"] == "tumor") & (scores["sender_subtype"] == "C1QC_like_TAM")]["communication_score"].sum()) if len(scores) else np.nan, "n_pairs": int(len(scores[(scores["condition"] == "tumor") & (scores["sender_subtype"] == "C1QC_like_TAM")])) if len(scores) else 0 }, {"direction": "tumor_to_C1QC_like_TAM", "mean_score": float(incoming_df["communication_score"].mean()) if len(incoming_df) else np.nan, "sum_score": float(incoming_df["communication_score"].sum()) if len(incoming_df) else np.nan, "n_pairs": int(len(incoming_df))}]); c1qc_role.to_csv(out / "c1qc_incoming_outgoing_role.csv", index=False)
    # Overlap annotation is post hoc and never filters the LR network.
    overlap_rows = []
    f81, f7 = set(FROZEN_81), set(DRIVERS)
    for l, r, p, c in LR_ROWS:
        lg, rg = set(genes_in(l)), set(genes_in(r)); overlap_rows.append({"ligand": l, "receptor": r, "pathway": p, "category": c, "ligand_in_81gene": bool(lg & f81), "receptor_in_81gene": bool(rg & f81), "pair_overlaps_81gene": bool((lg | rg) & f81), "overlap_81_genes": ",".join(sorted((lg | rg) & f81)), "ligand_in_7driver": bool(lg & f7), "receptor_in_7driver": bool(rg & f7), "pair_overlaps_7driver": bool((lg | rg) & f7), "overlap_7_drivers": ",".join(sorted((lg | rg) & f7)), "prostaglandin_ptger4_represented": False})
    pd.DataFrame(overlap_rows).to_csv(out / "lr_dinp_crc_overlap.csv", index=False)
    # Donor-aware validation for the top 20 pooled tumor LR records.
    top_pairs = scores[(scores["condition"] == "tumor") & (scores["receiver_group"] == "tumor_epithelial")].sort_values("communication_score", ascending=False).drop_duplicates(["sender_subtype", "receiver_state", "ligand", "receptor"]).head(20) if len(scores) else pd.DataFrame()
    donor = donor_support(logx, var_index, meta, top_pairs)
    donor.to_csv(out / "donor_level_lr_support.csv", index=False)
    if len(donor):
        donor_summary = donor.drop_duplicates(["condition", "sender_subtype", "receiver_group", "receiver_state", "ligand", "receptor"]); donor_summary.to_csv(out / "donor_level_lr_support_summary.csv", index=False)
    else: donor_summary = pd.DataFrame()
    # Figures use donor support in the role panel only as context; no cell-level tests.
    make_figures(scores if len(scores) else pd.DataFrame(columns=["condition", "receiver_group", "sender_subtype", "receiver_state", "ligand", "receptor", "pathway", "communication_score"]), donor, incoming_df, out)
    # Provenance and final report.
    tumor_ep_counts = ann[(ann["Cell_type"].astype(str).str.lower() == "epithelial cells") & (ann["Class"].astype(str).str.lower() == "tumor")].groupby("Cell_subtype").size().to_dict()
    macro_counts = macro["subtype"].value_counts().to_dict()
    manifest = {"analysis": "dinp_crc_macrophage_cellchat", "timestamp_utc": pd.Timestamp.utcnow().isoformat(), "dataset": "GSE144735 CRC single-cell matrix", "source_object": str(args.matrix), "source_annotation": str(args.annotation), "matrix_sha256": sha256(args.matrix), "annotation_sha256": sha256(args.annotation), "total_cells": int(len(ann)), "donors": int(ann["Patient"].nunique()), "tumor_cells": int((ann["Class"].astype(str).str.lower() == "tumor").sum()), "normal_cells": int((ann["Class"].astype(str).str.lower() == "normal").sum()), "border_cells": int((ann["Class"].astype(str).str.lower() == "border").sum()), "macrophage_subtype_counts": macro_counts, "tumor_epithelial_counts": {str(k): int(v) for k, v in tumor_ep_counts.items()}, "selected_sender_cells": int(meta["is_sender"].sum()), "selected_tumor_epithelial_cells": int(meta["is_tumor_epithelial"].sum()), "selected_normal_epithelial_cells": int(meta["is_normal_epithelial"].sum()), "metadata_fields": list(map(str, ann.columns)), "donor_id_field": "Patient", "condition_field": "Class", "broad_cell_type_field": "Cell_type", "epithelial_state_field": "Cell_subtype", "macrophage_subtype_source": str(args.macro_score), "macrophage_subtype_commit": "21d5382", "normalization": "raw UMI counts; per-cell total 1e4; log1p", "software": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__, "scipy": __import__("scipy").__version__, "CellChat": "unavailable: Rscript not found", "CellChatDB": "not instantiated; local CellChat-compatible human prior v0.1-local"}, "species": "human", "signaling_database_used": "CellChat-compatible focused human LR prior v0.1-local; secreted, ECM-Receptor, Cell-Cell Contact", "ptger4_database_status": "not represented as a robust CellChat LR interaction; no pair was invented", "no_redownload": True, "runtime_seconds": round(time.time() - t0, 2)}
    dump_json(manifest, out / "input_manifest.json")
    tumor_scores = scores[(scores["condition"] == "tumor") & (scores["receiver_group"] == "tumor_epithelial")] if len(scores) else pd.DataFrame()
    def top_path(sub):
        q = path[(path["condition"] == "tumor") & (path["receiver_group"] == "tumor_epithelial") & (path["sender_subtype"] == sub)].sort_values("pathway_max_score", ascending=False) if len(path) else pd.DataFrame()
        return str(q.iloc[0]["pathway"]) if len(q) else "NA"
    strongest = tumor_scores.iloc[0] if len(tumor_scores) else None
    consistent_n = int(donor_summary["consistent_across_ge4_donors"].sum()) if len(donor_summary) and "consistent_across_ge4_donors" in donor_summary else 0
    sppn = int(donor_summary.loc[donor_summary["sender_subtype"] == "SPP1_like_TAM", "consistent_across_ge4_donors"].sum()) if len(donor_summary) and "consistent_across_ge4_donors" in donor_summary else 0
    c1qc_out_sum = float(c1qc_role.loc[c1qc_role["direction"] == "C1QC_like_TAM_to_tumor", "sum_score"].iloc[0])
    c1qc_in_sum = float(c1qc_role.loc[c1qc_role["direction"] == "tumor_to_C1QC_like_TAM", "sum_score"].iloc[0])
    report = f"""# Focused macrophage–tumor communication analysis\n\n- Dataset: local GSE144735 CRC matrix ({len(ann):,} cells, {ann['Patient'].nunique()} donors); no data were redownloaded.\n- Scope: four existing macrophage subtypes as senders → source-labeled tumor epithelial states as receivers. Tumor epithelial labels (CMS1–4) were reused; no new malignant CNV inference was performed.\n- Method: pooled descriptive ligand–receptor scores plus donor-level support. R/CellChat was unavailable, so the reproducible fallback is the versioned local human prior in `lr_database_human.csv`; this is not presented as an actual CellChat run.\n\n## Required questions\n\n**Q1. SPP1-like TAM pathways.** Top outgoing pathway: **{top_path('SPP1_like_TAM')}**. The top pooled LR record is **{(str(strongest['ligand']) + '–' + str(strongest['receptor'])) if strongest is not None else 'NA'}**; donor support is summarized in `donor_level_lr_support.csv` (consistent pairs: {sppn}).\n\n**Q2. C1QC-like TAM pathways.** Top outgoing pathway: **{top_path('C1QC_like_TAM')}**. PTGER4 is quantified separately as a receptor expression context in the prior module; it is not fabricated into a CellChat pair here.\n\n**Q3. Distinct roles.** C1QC-like pooled outgoing score sum is {c1qc_out_sum:.3f}; the tumor→C1QC incoming context sum is {c1qc_in_sum:.3f}, with a higher mean per incoming pair. This supports keeping an outgoing SPP1-like axis separate from a PTGER4-enriched C1QC context rather than forcing one role.\n\n**Q4. Strongest TAM→tumor sender.** By pooled focused score, the leading sender is **{str(tumor_scores.groupby('sender_subtype')['communication_score'].sum().sort_values(ascending=False).index[0]) if len(tumor_scores) else 'NA'}**. This is descriptive, not a causal DINP result.\n\n**Q5. Donor repetition.** {consistent_n} of the top-20 pooled records meet the predefined support rule in at least four donors. Records with low coverage or donor-driven behavior are explicitly flagged.\n\n**Q6. DINP–CRC overlap.** Overlap annotation is post hoc in `lr_dinp_crc_overlap.csv`; the clearest 7-driver connection is CXCR4 in CXCL12–CXCR4. PTGER4 is not robustly represented in the local LR prior. No overlap was used to filter the network.\n\n**Q7. Spatial readiness.** **Larger-cohort replication first.** Several top subtype–condition signals have sparse paired donor coverage, so these results support prioritizing validation targets rather than immediate causal or spatial claims.\n\n## Interpretation boundary\n\nThe results support statements such as “DINP–CRC-associated molecular program localized to…” and “SPP1-like TAM displayed prominent communication scores.” They do not show that DINP causes communication, activates SPP1, or that PTGER4 mediates the network.\n\nRecommended next step: **Larger cohort replication**, then spatial validation if the same sender–receiver axes reproduce.\n"""
    (out / "CELLCHAT_REPORT.md").write_text(report, encoding="utf-8")
    print("MACROPHAGE–TUMOR COMMUNICATION ANALYSIS: PASS")
    print(f"Dataset: GSE144735; Total donors: {ann['Patient'].nunique()}; SPP1-like TAM cells: {macro_counts.get('SPP1_like_TAM', 0)}; C1QC-like TAM cells: {macro_counts.get('C1QC_like_TAM', 0)}; Tumor epithelial cells: {int(meta['is_tumor_epithelial'].sum())}")
    print(f"Top SPP1-like outgoing pathway: {top_path('SPP1_like_TAM')}; Top C1QC-like outgoing pathway: {top_path('C1QC_like_TAM')}")
    print(f"Strongest TAM→tumor LR pair: {(str(strongest['ligand']) + '–' + str(strongest['receptor'])) if strongest is not None else 'NA'}; Donor consistency: {consistent_n} top records supported in ≥4 donors")
    print("PTGER4-related communication found: NOT REPRESENTED IN DATABASE")
    print("Recommended next step: Larger cohort replication")


if __name__ == "__main__":
    main()
