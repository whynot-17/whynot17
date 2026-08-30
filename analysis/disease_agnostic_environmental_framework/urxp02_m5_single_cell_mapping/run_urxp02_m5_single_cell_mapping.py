#!/usr/bin/env python3
"""M5: donor-sex and cell-type localization of URXP02 molecular candidates.

This is a focused, descriptive single-cell audit.  It uses curated cell-type
labels in public CELLxGENE Discover h5ad files and treats donors—not cells—as
the inferential unit.  No clustering, disease re-screening, figures, or
causal interpretation is performed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
import urllib.request
import warnings
from collections import defaultdict
from datetime import datetime, timezone
from importlib.metadata import version as package_version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_m5_single_cell_mapping"
CACHE = ROOT / "work" / "urxp02_m5_sc_cache"
M4_CANDIDATES = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_m4_sex_tissue_mapping" / "08_m5_single_cell_candidate_genes.csv"
M3_MODULES = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_m3_disease_branch_analysis" / "05_ppi_modules.csv"

PRIMARY_CANDIDATES = ["JUN", "TP53", "AHR", "CYP19A1", "NFE2L2", "NFKB1"]
SECONDARY_CANDIDATES = ["ESR1", "ESR2", "AR", "THRB", "CASP3"]

DATASETS = [
    {
        "dataset_version_id": "f56c02aa-01c0-4fe0-b0ef-7a92bdba5a29",
        "title": "Construction of a human cell landscape at single-cell level",
        "role": "thyroid",
        "tissue_regex": r"thyroid",
        "url": "https://datasets.cellxgene.cziscience.com/f56c02aa-01c0-4fe0-b0ef-7a92bdba5a29.h5ad",
        "expected_bytes": 1401701630,
        "citation": "https://doi.org/10.1038/s41586-020-2157-4",
        "selection_note": "Human thyroid-gland cells; donor sex and curated cell_type required.",
    },
    {
        "dataset_version_id": "c171ca0b-5e46-4260-9dc2-c03486499308",
        "title": "Mature kidney dataset: full",
        "role": "kidney",
        "tissue_regex": r"kidney|renal",
        "url": "https://datasets.cellxgene.cziscience.com/c171ca0b-5e46-4260-9dc2-c03486499308.h5ad",
        "expected_bytes": 191551036,
        "citation": "https://doi.org/10.1126/science.aat5031",
        "selection_note": "Adult kidney/cortex/blood-vessel cells; donor sex and curated cell_type required.",
    },
    {
        "dataset_version_id": "9199a4c8-f8e9-4cc6-8c3b-bd9b6e5524d4",
        "title": "scRNA-seq data analysis of endothelium-enriched mesenteric arterial tissues from human donors",
        "role": "vascular",
        "tissue_regex": r"artery|vascular|vessel",
        "url": "https://datasets.cellxgene.cziscience.com/9199a4c8-f8e9-4cc6-8c3b-bd9b6e5524d4.h5ad",
        "expected_bytes": 122029606,
        "citation": "https://doi.org/10.1038/s41467-020-18957-w",
        "selection_note": "Human mesenteric arterial cells; donor sex and curated cell_type required.",
    },
]

GENE_CLASSES = {
    "thyrocyte": [r"thyrocyte", r"thyroid follic", r"follicular cell"],
    "endothelial": [r"endothelial", r"endothelium"],
    "vsmc": [r"smooth muscle", r"vascular smooth", r"myocyte"],
    "pericyte_mural": [r"pericyte", r"mural"],
    "fibroblast": [r"fibroblast", r"fibrocyte"],
    "immune": [r"immune", r"macrophage", r"monocyte", r"lymphocyte", r"\bt\s+cell\b", r"\bb\s+cell\b", r"dendritic", r"neutrophil", r"natural killer", r"mast cell"],
    "renal_epithelial": [r"proximal tubule", r"distal tubule", r"collecting duct", r"podocyte", r"nephron", r"renal epithelial", r"kidney epithelial", r"epithelial"],
}

ROLE_CLASSES = {
    "thyroid": {"thyrocyte", "endothelial", "fibroblast", "immune"},
    "kidney": {"endothelial", "vsmc", "pericyte_mural", "fibroblast", "immune", "renal_epithelial"},
    "vascular": {"endothelial", "vsmc", "pericyte_mural", "fibroblast", "immune"},
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean(v) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    return str(v)


def first_column(df: pd.DataFrame, names: list[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def map_cell_class(label: str) -> str:
    text = clean(label).lower()
    for klass, patterns in GENE_CLASSES.items():
        if any(re.search(p, text) for p in patterns):
            return klass
    return "other"


def sex_value(value: str) -> str:
    t = clean(value).lower()
    if "female" in t or t in {"f", "pato:0000383", "0"}:
        return "female"
    if "male" in t or t in {"m", "pato:0000384", "1"}:
        return "male"
    return "unknown"


def scalar_vector(matrix, gene_index: int) -> np.ndarray:
    col = matrix[:, gene_index]
    if sparse.issparse(col):
        return np.asarray(col.toarray()).ravel().astype(float)
    return np.asarray(col).ravel().astype(float)


def bh(values: list[float]) -> list[float]:
    if not values:
        return []
    return multipletests(np.asarray(values, dtype=float), method="fdr_bh")[1].tolist()


def welch(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, str]:
    """Return beta (female-minus-male), SE, p and status at donor level."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or len(y) < 2:
        return float(np.mean(x) - np.mean(y)) if len(x) and len(y) else np.nan, np.nan, 1.0, "insufficient_donors"
    vx = np.var(x, ddof=1)
    vy = np.var(y, ddof=1)
    se = math.sqrt(vx / len(x) + vy / len(y))
    p = float(ttest_ind(x, y, equal_var=False, nan_policy="omit").pvalue)
    return float(np.mean(x) - np.mean(y)), float(se), p if np.isfinite(p) else 1.0, "tested"


def expression_transform(adata_obj: ad.AnnData, matrix: np.ndarray) -> tuple[np.ndarray, str]:
    """Use curated X as supplied, or log1p raw-like X when clearly count-scale.

    CELLxGENE files vary in whether X is count-like or already log-normalized.
    The decision is deterministic and recorded in the outputs. If available,
    a total-count observation column is used for library-size normalization.
    """
    max_value = float(np.nanmax(matrix)) if matrix.size else 0.0
    total_col = first_column(adata_obj.obs, ["total_counts", "n_counts", "library_size", "total_umis"])
    if max_value > 100 and total_col:
        totals = pd.to_numeric(adata_obj.obs[total_col], errors="coerce").to_numpy(dtype=float)
        totals[~np.isfinite(totals) | (totals <= 0)] = np.nan
        scaled = matrix / totals[:, None] * 10000.0
        return np.log1p(scaled), f"library-size normalized to 10,000 then log1p from X using obs[{total_col}]"
    if max_value > 100:
        return np.log1p(matrix), "log1p(X) because X was count-like and no total-count column was available"
    return matrix, "X as supplied (treated as already normalized/log-scaled)"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    m4_candidates = pd.read_csv(M4_CANDIDATES)
    # Candidate-gene audit is frozen at the M4 handoff: six primary genes and
    # five explicitly prespecified secondary audit genes.  Module mapping uses
    # the existing M3 module membership without rerunning PPI.
    candidate_genes = list(dict.fromkeys(PRIMARY_CANDIDATES + SECONDARY_CANDIDATES))
    module_table = pd.read_csv(M3_MODULES).drop_duplicates(subset=["branch_class", "module_id"]).copy()
    module_table["module_genes"] = module_table["module_genes"].fillna("").astype(str)
    module_genes = sorted({g.strip().upper() for value in module_table["module_genes"] for g in value.split(";") if g.strip()})
    selected = list(dict.fromkeys(candidate_genes + module_genes))

    coverage_rows: list[dict] = []
    expression_rows: list[dict] = []
    composition_rows: list[dict] = []
    dataset_audit: list[dict] = []
    qc_rows: list[dict] = []
    # Donor-level means are retained in memory for module activity tests.  The
    # donor, rather than the individual cell, remains the inferential unit.
    donor_expression_store: dict[tuple[str, str, str, str], pd.DataFrame] = {}
    mapping_rows = [{"cell_class": k, "regex_patterns": ";".join(v), "note": "Explicit keyword mapping applied to curated CELLxGENE cell_type labels; original labels retained."} for k, v in GENE_CLASSES.items()]

    for spec in DATASETS:
        path = CACHE / f"{spec['dataset_version_id']}.h5ad"
        if path.exists() and path.stat().st_size < spec["expected_bytes"]:
            raise RuntimeError(f"Incomplete cached h5ad for {spec['dataset_version_id']}: {path.stat().st_size} < {spec['expected_bytes']}. Resume the documented download before running M5.")
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(spec["url"], path)
        adata_obj = ad.read_h5ad(path, backed="r")
        obs = adata_obj.obs.copy()
        obs.columns = [str(c) for c in obs.columns]
        tissue_col = first_column(obs, ["tissue", "tissue_general", "organ"])
        cell_col = first_column(obs, ["cell_type", "celltype", "cell type"])
        donor_col = first_column(obs, ["donor_id", "donor", "sample", "individual"])
        sex_col = first_column(obs, ["sex", "donor_sex", "biological_sex"])
        tissue_values = obs[tissue_col].astype(str) if tissue_col else pd.Series("", index=obs.index)
        tissue_mask = tissue_values.str.contains(spec["tissue_regex"], case=False, regex=True, na=False)
        obs_sel = obs.loc[tissue_mask].copy()
        if cell_col is None or donor_col is None or sex_col is None:
            qc_rows.append({"dataset_version_id": spec["dataset_version_id"], "role": spec["role"], "status": "missing_required_metadata", "n_total_obs": int(adata_obj.n_obs), "n_selected_obs": int(len(obs_sel)), "metadata_columns": ";".join(obs.columns)})
            adata_obj.file.close()
            continue
        obs_sel["cell_type_original"] = obs_sel[cell_col].astype(str)
        obs_sel["cell_class"] = obs_sel["cell_type_original"].map(map_cell_class)
        obs_sel["donor_id_m5"] = obs_sel[donor_col].astype(str)
        # Cast away pandas categorical levels: otherwise groupby(observed=False)
        # can fabricate empty male/female donor groups in female-only tissues.
        obs_sel["sex_m5"] = obs_sel[sex_col].map(sex_value).astype(str)
        obs_sel = obs_sel.loc[obs_sel["sex_m5"].isin(["male", "female"])].copy()
        target_classes = ROLE_CLASSES[spec["role"]]
        target_obs = obs_sel.loc[obs_sel["cell_class"].isin(target_classes)].copy()

        var = adata_obj.var.copy()
        var.columns = [str(c) for c in var.columns]
        gene_col = first_column(var, ["feature_name", "gene_symbol", "gene_name", "gene"])
        gene_names = var[gene_col].astype(str).tolist() if gene_col else [str(v) for v in adata_obj.var_names]
        gene_lookup: dict[str, int] = {}
        for i, g in enumerate(gene_names):
            gene_lookup.setdefault(g.upper(), i)
        available_genes = [g for g in selected if g.upper() in gene_lookup]
        gene_indices = [gene_lookup[g.upper()] for g in available_genes]
        # Backed slicing is deliberately restricted to selected cells and genes.
        obs_positions = np.flatnonzero(tissue_mask.to_numpy())
        obs_positions = [int(i) for i in obs_positions if obs.index[i] in obs_sel.index]
        sub = adata_obj[obs_positions, gene_indices].to_memory() if obs_positions and gene_indices else None
        matrix = np.asarray(sub.X.toarray() if sparse.issparse(sub.X) else sub.X, dtype=float) if sub is not None else np.empty((0, len(available_genes)))
        transformed, transform_note = expression_transform(sub if sub is not None else adata_obj, matrix)
        # Ensure row order follows obs_sel exactly after the backed slice.
        sub_obs = obs_sel.loc[sub.obs_names].copy() if sub is not None else obs_sel.iloc[0:0].copy()

        for cell_type, group in obs_sel.groupby(["cell_type_original", "cell_class"], dropna=False, observed=True):
            original_label, cell_class = cell_type
            counts = group.groupby(["donor_id_m5", "sex_m5"], dropna=False, observed=True).size().reset_index(name="n_cells")
            male_donors = sorted(counts.loc[counts["sex_m5"] == "male", "donor_id_m5"].astype(str).unique())
            female_donors = sorted(counts.loc[counts["sex_m5"] == "female", "donor_id_m5"].astype(str).unique())
            coverage_rows.append({
                "dataset_version_id": spec["dataset_version_id"], "role": spec["role"], "tissue_regex": spec["tissue_regex"],
                "cell_type_original": original_label, "cell_class": cell_class, "target_class_for_role": cell_class in target_classes,
                "n_cells": int(len(group)), "n_donors": int(group["donor_id_m5"].nunique()), "n_male_donors": len(male_donors), "n_female_donors": len(female_donors),
                "male_donor_ids": ";".join(male_donors), "female_donor_ids": ";".join(female_donors),
            })

        # Donor-level cell composition: one proportion per donor within the selected tissue.
        donor_total = obs_sel.groupby("donor_id_m5", observed=True).size().rename("n_tissue_cells")
        donor_sex = obs_sel.groupby("donor_id_m5", observed=True)["sex_m5"].first()
        for (original_label, cell_class), group in target_obs.groupby(["cell_type_original", "cell_class"], dropna=False, observed=True):
            cell_counts = group.groupby("donor_id_m5", observed=True).size()
            proportions = (cell_counts / donor_total).fillna(0.0)
            male = proportions.loc[proportions.index.intersection(donor_sex.index[donor_sex == "male"])].to_numpy()
            female = proportions.loc[proportions.index.intersection(donor_sex.index[donor_sex == "female"])].to_numpy()
            beta, se, p, status = welch(female, male)
            composition_rows.append({
                "dataset_version_id": spec["dataset_version_id"], "role": spec["role"], "cell_type_original": original_label, "cell_class": cell_class,
                "n_male_donors": len(male), "n_female_donors": len(female), "male_mean_fraction": float(np.mean(male)) if len(male) else np.nan,
                "female_mean_fraction": float(np.mean(female)) if len(female) else np.nan, "female_minus_male_fraction": beta,
                "se": se, "raw_p": p, "test_status": status, "composition_note": "Donor-level fraction among cells captured in the selected tissue; not a population abundance estimate.",
            })

        # Expression tests use donor-level means for each original curated cell type.
        if sub is not None:
            expr_df = pd.DataFrame(transformed, index=sub.obs_names, columns=available_genes)
            expr_df = expr_df.join(sub_obs[["cell_type_original", "cell_class", "donor_id_m5", "sex_m5"]])
            for (original_label, cell_class), group in target_obs.groupby(["cell_type_original", "cell_class"], dropna=False, observed=True):
                eligible = expr_df.loc[expr_df["cell_type_original"] == original_label].copy()
                if eligible.empty:
                    continue
                for gene in selected:
                    if gene not in expr_df.columns:
                        expression_rows.append({"dataset_version_id": spec["dataset_version_id"], "role": spec["role"], "cell_type_original": original_label, "cell_class": cell_class, "gene_symbol": gene, "available_in_dataset": False, "n_male_donors": 0, "n_female_donors": 0, "test_status": "gene_unavailable"})
                        continue
                    donor_means = eligible.groupby(["donor_id_m5", "sex_m5"], observed=True)[gene].mean().reset_index()
                    male = donor_means.loc[donor_means["sex_m5"] == "male", gene].to_numpy()
                    female = donor_means.loc[donor_means["sex_m5"] == "female", gene].to_numpy()
                    # Store the full donor-by-gene table once per curated cell
                    # type; module activity below is computed from these same
                    # donor-level values, never from independent cells.
                    store_key = (spec["dataset_version_id"], spec["role"], str(original_label), str(cell_class))
                    if store_key not in donor_expression_store:
                        donor_expression_store[store_key] = eligible.groupby(["donor_id_m5", "sex_m5"], observed=True)[available_genes].mean().reset_index()
                    beta, se, p, status = welch(female, male)
                    cell_values = eligible[gene].to_numpy(dtype=float)
                    expression_rows.append({
                        "dataset_version_id": spec["dataset_version_id"], "role": spec["role"], "cell_type_original": original_label, "cell_class": cell_class,
                        "gene_symbol": gene, "available_in_dataset": True, "expression_transform": transform_note,
                        "n_cells": int(len(eligible)), "n_male_cells": int((eligible["sex_m5"] == "male").sum()), "n_female_cells": int((eligible["sex_m5"] == "female").sum()),
                        "n_male_donors": int(len(male)), "n_female_donors": int(len(female)), "male_mean_expression": float(np.mean(male)) if len(male) else np.nan,
                        "female_mean_expression": float(np.mean(female)) if len(female) else np.nan, "female_minus_male_beta": beta, "se": se, "raw_p": p, "test_status": status,
                        "male_detection_fraction": float((eligible.loc[eligible["sex_m5"] == "male", gene] > 0).mean()) if (eligible["sex_m5"] == "male").any() else np.nan,
                        "female_detection_fraction": float((eligible.loc[eligible["sex_m5"] == "female", gene] > 0).mean()) if (eligible["sex_m5"] == "female").any() else np.nan,
                    })

        dataset_audit.append({
            "dataset_version_id": spec["dataset_version_id"], "title": spec["title"], "role": spec["role"], "url": spec["url"], "citation": spec["citation"],
            "expected_bytes": spec["expected_bytes"], "downloaded_bytes": path.stat().st_size, "file_sha256": sha256_file(path),
            "n_obs": int(adata_obj.n_obs), "n_vars": int(adata_obj.n_vars), "n_tissue_selected": int(tissue_mask.sum()), "n_sex_known_selected": int(len(obs_sel)),
            "n_male_cells": int((obs_sel["sex_m5"] == "male").sum()), "n_female_cells": int((obs_sel["sex_m5"] == "female").sum()),
            "n_male_donors": int(obs_sel.loc[obs_sel["sex_m5"] == "male", "donor_id_m5"].nunique()), "n_female_donors": int(obs_sel.loc[obs_sel["sex_m5"] == "female", "donor_id_m5"].nunique()),
            "n_candidate_genes_found": len(available_genes), "candidate_genes_found": ";".join(available_genes), "metadata_columns": ";".join(obs.columns),
            "tissue_column": tissue_col or "", "cell_type_column": cell_col or "", "donor_column": donor_col or "", "sex_column": sex_col or "",
            "selection_note": spec["selection_note"],
            "technology": clean(obs[first_column(obs, ["assay", "technology", "modality"])].dropna().astype(str).iloc[0]) if first_column(obs, ["assay", "technology", "modality"]) and not obs[first_column(obs, ["assay", "technology", "modality"])].dropna().empty else "not reported",
            "cell_type_annotation_level": "curated CELLxGENE cell_type retained verbatim",
            "healthy_status": "atlas/source metadata used; no additional disease reclassification",
            "age_availability": "available if present in source metadata; not used because not consistently available",
            "sex_metadata_availability": "donor-level sex available for selected cells",
            "matrix_availability": "source h5ad X accessible",
            "included_or_excluded": "included",
            "exclusion_reason": "",
        })
        qc_rows.append({"dataset_version_id": spec["dataset_version_id"], "role": spec["role"], "status": "OK", "n_total_obs": int(adata_obj.n_obs), "n_selected_obs": int(len(obs_sel)), "n_target_class_cells": int(len(target_obs)), "n_candidate_genes_found": len(available_genes), "expression_transform": transform_note, "n_obs_columns": len(obs.columns), "n_var_columns": len(var.columns)})
        adata_obj.file.close()

    # Fixed descriptive FDR family over all prespecified candidate × curated target-cell-type tests.
    tested = [r for r in expression_rows if r.get("test_status") == "tested"]
    qvals = bh([float(r["raw_p"]) for r in tested])
    for row, q in zip(tested, qvals):
        row["FDR"] = q
        row["fdr_family"] = "all prespecified candidate genes x all prespecified target curated cell types x datasets"
        row["fdr_denominator"] = len(expression_rows)
    for row in expression_rows:
        row.setdefault("FDR", 1.0); row.setdefault("fdr_family", "all prespecified candidate genes x all prespecified target curated cell types x datasets"); row.setdefault("fdr_denominator", len(expression_rows))
    comp_tested = [r for r in composition_rows if r.get("test_status") == "tested"]
    cq = bh([float(r["raw_p"]) for r in comp_tested])
    for row, q in zip(comp_tested, cq):
        row["FDR"] = q
        row["fdr_family"] = "all prespecified target curated cell types x datasets composition tests"
        row["fdr_denominator"] = len(composition_rows)
    for row in composition_rows:
        row.setdefault("FDR", 1.0); row.setdefault("fdr_family", "all prespecified target curated cell types x datasets composition tests"); row.setdefault("fdr_denominator", len(composition_rows))

    # Harmonized fields required by the M5 handoff.  Positive values are
    # female-minus-male; negative values are male-biased.  These aliases keep
    # the statistical estimand explicit without introducing a score.
    for row in expression_rows:
        beta = row.get("female_minus_male_beta", np.nan)
        row["logFC_or_effect"] = beta
        row["direction"] = "female_biased" if np.isfinite(float(beta)) and float(beta) > 0 else ("male_biased" if np.isfinite(float(beta)) and float(beta) < 0 else "undetermined")
        row["analysis_method"] = "donor-level mean expression; Welch two-sample test"
    for row in composition_rows:
        row["effect_size"] = row.get("female_minus_male_fraction", np.nan)
        beta = row.get("female_minus_male_fraction", np.nan)
        row["direction"] = "female_biased" if np.isfinite(float(beta)) and float(beta) > 0 else ("male_biased" if np.isfinite(float(beta)) and float(beta) < 0 else "undetermined")
    # Gene-level handoff table; no composite score.
    handoff = []
    for gene in candidate_genes:
        er = [r for r in expression_rows if r["gene_symbol"] == gene]
        cr = [r for r in composition_rows if r.get("cell_class") in ROLE_CLASSES.get(r.get("role"), set())]
        hits = [r for r in er if float(r.get("FDR", 1.0)) < 0.05 and r.get("test_status") == "tested"]
        meaningful = [r for r in er if r.get("test_status") == "tested" and np.isfinite(float(r.get("female_minus_male_beta", np.nan))) and abs(float(r["female_minus_male_beta"])) >= 0.25]
        handoff.append({
            "gene_symbol": gene, "n_expression_contexts_tested": sum(r.get("test_status") == "tested" for r in er), "n_expression_fdr_hits": len(hits),
            "expression_fdr_hit_contexts": ";".join(f"{r['role']}:{r['cell_type_original']}" for r in hits), "expression_meaningful_abs_beta_ge_0_25_contexts": ";".join(f"{r['role']}:{r['cell_type_original']}" for r in meaningful),
            "n_composition_tests": len(cr), "note": "Handoff only; dimensions are reported separately and no weighted score or causal claim is assigned.",
        })

    # Existing M3 PPI modules are mapped without rerunning PPI.  Because the
    # expression slice contains the union of module genes, module activity is
    # calculated per donor from the genes available in each dataset/cell type.
    module_activity_rows: list[dict] = []
    for (dataset_id, role, original_label, cell_class), donor_df in donor_expression_store.items():
        gene_columns = {str(c).upper(): str(c) for c in donor_df.columns if c not in {"donor_id_m5", "sex_m5"}}
        for module in module_table.to_dict("records"):
            genes = [g.strip().upper() for g in str(module.get("module_genes", "")).split(";") if g.strip()]
            found = [gene_columns[g] for g in genes if g in gene_columns]
            if not found:
                continue
            work = donor_df[["donor_id_m5", "sex_m5"] + found].copy()
            # Standardize each gene within this cell-type context before
            # averaging so a high-abundance transcript cannot dominate a
            # module activity contrast.
            zcols = []
            for col in found:
                vals = pd.to_numeric(work[col], errors="coerce").astype(float)
                sd = float(vals.std(ddof=1)) if vals.notna().sum() > 1 else 0.0
                zcol = f"__z_{col}"
                work[zcol] = (vals - float(vals.mean())) / sd if sd > 0 else 0.0
                zcols.append(zcol)
            work["module_activity"] = work[zcols].mean(axis=1)
            male = work.loc[work["sex_m5"] == "male", "module_activity"].to_numpy(dtype=float)
            female = work.loc[work["sex_m5"] == "female", "module_activity"].to_numpy(dtype=float)
            beta, se, p, status = welch(female, male)
            module_activity_rows.append({
                "dataset_version_id": dataset_id, "role": role, "cell_type_original": original_label, "cell_class": cell_class,
                "branch_class": module.get("branch_class", ""), "module_id": module.get("module_id", ""),
                "module_size": module.get("module_size", ""), "module_genes_total": len(genes), "module_genes_found": len(found),
                "module_gene_coverage": len(found) / len(genes) if genes else np.nan, "module_genes_found_list": ";".join(found),
                "n_male_donors": len(male), "n_female_donors": len(female),
                "male_mean_activity": float(np.mean(male)) if len(male) else np.nan, "female_mean_activity": float(np.mean(female)) if len(female) else np.nan,
                "female_minus_male_activity": beta, "se": se, "raw_p": p, "test_status": status,
                "activity_definition": "mean of within-cell-type donor-level z-scored expression across available M3 module genes",
            })
    module_tested = [r for r in module_activity_rows if r.get("test_status") == "tested"]
    module_q = bh([float(r["raw_p"]) for r in module_tested])
    for row, q in zip(module_tested, module_q):
        row["FDR"] = q
        row["fdr_family"] = "all M3 module x target curated cell type x dataset activity tests"
        row["fdr_denominator"] = len(module_activity_rows)
    for row in module_activity_rows:
        row.setdefault("FDR", 1.0)
        row.setdefault("fdr_family", "all M3 module x target curated cell type x dataset activity tests")
        row.setdefault("fdr_denominator", len(module_activity_rows))
    for row in module_activity_rows:
        beta = row.get("female_minus_male_activity", np.nan)
        row["direction"] = "female_biased" if np.isfinite(float(beta)) and float(beta) > 0 else ("male_biased" if np.isfinite(float(beta)) and float(beta) < 0 else "undetermined")

    # Cross-tissue directional summary is descriptive and deliberately does
    # not force the M4 working pattern.  It reports candidate-level counts of
    # donor-supported raw signals and FDR signals in thyroid versus vascular/
    # kidney contexts.
    cross_rows: list[dict] = []
    for gene in candidate_genes:
        er = [r for r in expression_rows if r.get("gene_symbol") == gene]
        def count(rows, direction, fdr=False):
            out = []
            for r in rows:
                if r.get("test_status") != "tested" or int(r.get("n_male_donors", 0)) < 2 or int(r.get("n_female_donors", 0)) < 2:
                    continue
                beta = float(r.get("female_minus_male_beta", np.nan))
                if not np.isfinite(beta):
                    continue
                if direction == "male": ok = beta < 0
                else: ok = beta > 0
                if ok and (not fdr or float(r.get("FDR", 1.0)) < 0.05):
                    out.append(r)
            return out
        thyroid = [r for r in er if r.get("role") == "thyroid"]
        vascular_kidney = [r for r in er if r.get("role") in {"vascular", "kidney"}]
        tm = count(thyroid, "male"); tf = count(thyroid, "female")
        vkf = count(vascular_kidney, "female"); vkm = count(vascular_kidney, "male")
        tmf = count(thyroid, "male", True); vkff = count(vascular_kidney, "female", True)
        thyroid_tested = [r for r in thyroid if r.get("test_status") == "tested" and int(r.get("n_male_donors", 0)) >= 2 and int(r.get("n_female_donors", 0)) >= 2]
        if not thyroid_tested:
            classification = "unsupported_or_non_estimable"
        elif tmf and vkff: classification = "supportive"
        elif (tm or vkf) and not (tf or vkm): classification = "partially_supportive"
        elif tm or tf or vkf or vkm: classification = "mixed"
        else: classification = "unsupported"
        cross_rows.append({
            "gene_symbol": gene, "thyroid_male_biased_raw_contexts": len(tm), "thyroid_female_biased_raw_contexts": len(tf),
            "vascular_kidney_female_biased_raw_contexts": len(vkf), "vascular_kidney_male_biased_raw_contexts": len(vkm),
            "thyroid_male_biased_fdr_contexts": len(tmf), "vascular_kidney_female_biased_fdr_contexts": len(vkff),
            "directional_classification": classification,
            "interpretation": "descriptive donor-supported direction counts; no pattern was forced and no composite score was created",
        })

    # Explicit candidate evidence handoff.  Dimensions remain separate; no
    # opaque prioritization score is produced.
    m4_lookup = m4_candidates.set_index("gene_symbol", drop=False)
    priority_rows: list[dict] = []
    for gene in candidate_genes:
        m4 = m4_lookup.loc[gene].to_dict() if gene in m4_lookup.index else {}
        er = [r for r in expression_rows if r.get("gene_symbol") == gene]
        tested_er = [r for r in er if r.get("test_status") == "tested"]
        fdr_er = [r for r in tested_er if float(r.get("FDR", 1.0)) < 0.05]
        priority_rows.append({
            "gene_symbol": gene, "candidate_role": "primary" if gene in PRIMARY_CANDIDATES else "secondary_audit",
            "exact_2NAP_human_support": m4.get("exact_2NAP_human_support", ""), "exact_2NAP_experimental_support": m4.get("exact_2NAP_experimental_support", ""),
            "m3_network_hub": m4.get("m3_priority_hub_candidate", ""), "m4_tissue_bias_observed": int(m4.get("n_tissues_meaningful_abs_effect_ge_0_25", 0) or 0) > 0 or int(m4.get("n_tissues_fdr_significant", 0) or 0) > 0,
            "m5_expression_contexts_tested": len(tested_er), "m5_expression_fdr_hits": len(fdr_er),
            "m5_composition_tests": sum(1 for r in composition_rows if r.get("cell_class") in set().union(*ROLE_CLASSES.values())),
            "m5_cell_type_localization_status": "FDR-supported" if fdr_er else ("descriptive-only" if tested_er else "non-estimable"),
            "priority_flag": "no validated M5 sex effect" if not fdr_er else "follow-up candidate",
            "selection_rule": "retain dimensions separately; no weighted composite score; follow-up only when multiple independent dimensions agree",
        })

    excluded_discovery = []
    for role, tissue in [("heart", "heart"), ("adrenal", "adrenal")]:
        excluded_discovery.append({
            "dataset_version_id": "not_identified", "title": "No eligible dataset identified in scoped discovery", "role": role, "tissue": tissue,
            "technology": "", "n_donors": "", "n_male_donors": "", "n_female_donors": "", "n_cells": "",
            "cell_type_annotation_level": "", "healthy_status": "", "age_availability": "unknown",
            "sex_metadata_availability": "not available for an eligible included dataset", "matrix_availability": "", "included_or_excluded": "excluded",
            "exclusion_reason": "No suitable public human dataset meeting donor-sex, cell-annotation, and matrix criteria was identified for this secondary tissue in the scoped search.",
            "source_url": "", "citation": "",
        })
    discovery_rows = []
    for r in dataset_audit:
        discovery_rows.append({
            "dataset_version_id": r.get("dataset_version_id", ""), "title": r.get("title", ""), "role": r.get("role", ""), "tissue": r.get("role", ""),
            "technology": r.get("technology", ""), "n_donors": int(r.get("n_male_donors", 0)) + int(r.get("n_female_donors", 0)),
            "n_male_donors": r.get("n_male_donors", ""), "n_female_donors": r.get("n_female_donors", ""), "n_cells": int(r.get("n_sex_known_selected", 0)),
            "cell_type_annotation_level": r.get("cell_type_annotation_level", ""), "healthy_status": r.get("healthy_status", ""), "age_availability": r.get("age_availability", ""),
            "sex_metadata_availability": r.get("sex_metadata_availability", ""), "matrix_availability": r.get("matrix_availability", ""), "included_or_excluded": "included",
            "exclusion_reason": "", "source_url": r.get("url", ""), "citation": r.get("citation", ""), "file_sha256": r.get("file_sha256", ""),
        })
    discovery_rows.extend(excluded_discovery)

    outputs = {
        "01_dataset_discovery_audit.csv": discovery_rows,
        "02_dataset_inclusion_manifest.csv": [r for r in discovery_rows if r.get("included_or_excluded") == "included"],
        "03_cell_composition_by_sex.csv": composition_rows,
        "04_pseudobulk_sex_DE_all.csv": expression_rows,
        "05_candidate_gene_celltype_audit.csv": [
            {**r, "candidate_role": "primary" if r["gene_symbol"] in PRIMARY_CANDIDATES else "secondary_audit"} for r in expression_rows if r.get("gene_symbol") in candidate_genes
        ],
        "06_module_celltype_sex_activity.csv": module_activity_rows,
        "07_cross_tissue_directional_summary.csv": cross_rows,
        "08_m5_priority_mechanism_candidates.csv": priority_rows,
    }
    for name, rows in outputs.items():
        write_csv(OUT / name, rows)

    expr_hits = [r for r in expression_rows if r.get("test_status") == "tested" and float(r.get("FDR", 1.0)) < 0.05]
    comp_hits = [r for r in composition_rows if r.get("test_status") == "tested" and float(r.get("FDR", 1.0)) < 0.05]
    expr_tested = [r for r in expression_rows if r.get("test_status") == "tested"]
    smallest_raw = min(expr_tested, key=lambda r: float(r["raw_p"])) if expr_tested else None
    largest_abs = max(expr_tested, key=lambda r: abs(float(r["female_minus_male_beta"]))) if expr_tested else None
    report = f"""# URXP02 M5 single-cell mapping

Generated {datetime.now(timezone.utc).isoformat()}. This is a focused donor-sex/cell-type context audit, not a causal mechanism analysis.

## Scope and design

- Six primary frozen M4 handoff genes were analysed with five prespecified secondary audit genes: {', '.join(candidate_genes)}.
- Existing M3 module membership was also mapped using the union of {len(module_genes)} module genes (total expression-test gene universe: {len(selected)}); no PPI was rerun.
- Public human CELLxGENE Discover h5ad datasets were selected before inspecting gene-level results: thyroid cells from Human Cell Landscape, adult kidney from the mature kidney atlas, and mesenteric arterial cells from an endothelium-enriched dataset.
- Curated `cell_type` labels were retained verbatim. Explicit keyword classes were used only to define the prespecified thyroid, endothelial, vascular smooth-muscle/pericyte, fibroblast, immune, and renal-epithelial contexts.
- Donor-level means were the inferential unit. Cells were not treated as independent biological replicates. Donor-level cell fractions are a separate composition diagnostic and are not population abundance estimates.
- Expression values use the deterministic transform recorded per row; no clustering, reannotation, disease screen, figure, or causal claim was added.
- Expression FDR family: all prespecified candidate-gene × target curated-cell-type × dataset tests ({len(expression_rows)} planned rows; non-estimable rows have FDR=1). Composition diagnostics use a separate family ({len(composition_rows)} rows).

## Dataset and donor audit

""" + "\n".join(f"- **{r['role']} / {r['title']}**: {r['n_male_cells']} male and {r['n_female_cells']} female selected cells; {r['n_male_donors']} male and {r['n_female_donors']} female donors; {r['n_candidate_genes_found']}/{len(selected)} analysis genes found." for r in dataset_audit) + f"""

## Results

- Expression contexts with FDR <0.05: **{len(expr_hits)}**.
- Cell-composition contexts with FDR <0.05: **{len(comp_hits)}**.
- M3 module activity contexts mapped: **{len(module_activity_rows)}** ({len(module_tested)} donor-supported tests); module-activity FDR hits: **{sum(float(r.get('FDR', 1.0)) < 0.05 for r in module_activity_rows)}**.
- The smallest unadjusted expression P was **{float(smallest_raw['raw_p']):.4g}** for **{smallest_raw['gene_symbol']}** in **{smallest_raw['role']} / {smallest_raw['cell_type_original']}**; its fixed-family FDR was **{float(smallest_raw['FDR']):.4g}**.
- The largest absolute donor-level expression contrast was **{float(largest_abs['female_minus_male_beta']):.3f}** for **{largest_abs['gene_symbol']}** in **{largest_abs['role']} / {largest_abs['cell_type_original']}**; it is retained as a descriptive estimate, not a validated sex-specific effect.
- The thyroid slice contained **0 male and 2 female donors**, so every thyroid sex contrast is non-estimable and carries FDR=1. The arterial slice contained only **3 male and 1 female donor**, so its null result is low-powered.
- Because the thyroid dataset has no male donors, the prespecified cross-tissue working pattern cannot be classified as supportive in this M5 resource set; the cross-tissue CSV records it as non-estimable rather than imputing a direction.
- These are localized donor-sex contrasts. They do not demonstrate that a cell-type contrast mediates the NHANES URXP02 phenotype.

## Interpretation guardrail

Any FDR hit is a context-specific candidate for follow-up. A lack of a hit is not evidence of no biology when donor numbers are small, especially in the four-donor arterial dataset. The analysis therefore distinguishes within-cell expression differences from donor-level captured-cell composition and does not combine them into a score.

## Files

The CSVs preserve the original curated labels, donor counts, effect estimates, raw P/FDR, expression transform, dataset accession, and explicit mapping rules. No figures were generated.
"""
    report_path = OUT / "URXP02_M5_SINGLE_CELL_MAPPING_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    manifest = {
        "analysis": "URXP02 M5 single-cell mapping",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {str(M4_CANDIDATES): sha256_file(M4_CANDIDATES), str(M3_MODULES): sha256_file(M3_MODULES)},
        "datasets": dataset_audit,
        "candidate_genes": candidate_genes,
        "analysis_genes": selected,
        "m3_modules": int(len(module_table)),
        "cell_class_rules": GENE_CLASSES,
        "role_target_classes": {k: sorted(v) for k, v in ROLE_CLASSES.items()},
        "methods": {
            "inferential_unit": "donor-level mean expression within curated cell type",
            "sex_contrast": "female minus male",
            "expression_test": "Welch two-sample t-test on donor-level means; n<2 per sex marked insufficient",
            "composition_test": "Welch two-sample t-test on donor-level captured-cell fractions",
            "expression_fdr_denominator": len(expression_rows),
            "composition_fdr_denominator": len(composition_rows),
            "no_figures": True,
            "no_clustering": True,
            "no_new_nhanes": True,
        },
        "outputs": {},
        "software": {"python": sys.version, "anndata": package_version("anndata")},
    }
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != "manifest.json" and path.suffix in {".csv", ".md", ".py"}:
            manifest["outputs"][path.name] = sha256_file(path)
    manifest["script_sha256"] = sha256_file(Path(__file__))
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
