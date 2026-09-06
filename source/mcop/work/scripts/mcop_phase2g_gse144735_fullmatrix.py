"""Phase 2G alternative primary analysis using the complete GSE144735 matrix.

This path is deliberately independent of the live CELLxGENE Census/TileDB
raw-expression query.  It uses the GEO Series processed natural-log TPM
matrix and the matched patient annotations from GSE144735.  The frozen
7-gene PPAR/NR core, epithelial state universe, donor-aware inference, and
DoRothEA A-C regulator network are retained.  The dataset is small (six
paired patients), so results are reported as an alternative primary analysis
with explicit power and multi-dataset limitations.
"""

from __future__ import annotations

import gzip
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as t_dist, wilcoxon
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

import mcop_phase2g_epithelial_state_analysis as phase2g


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs"
RAW = ROOT / "work" / "mcop_phase2g" / "raw"
MATRIX = RAW / "GSE144735_processed_KUL3_CRC_10X_natural_log_TPM_matrix.full.gz"
RAW_CORE_MATRIX = ROOT / "work" / "mcop_phase2f_external" / "raw" / "GSE144735_raw_UMI_count_matrix.txt.gz"
ANNOTATION = ROOT / "work" / "mcop_phase2f_external" / "raw" / "GSE144735_annotation.txt.gz"
DOROTHEA_RAW = ROOT / "work" / "mcop_phase2g" / "dorothea_raw.tsv"

PREFIX = "mcop_phase2g_gse144735"
OUT_CELL = OUTPUT / f"{PREFIX}_cell_state_scores.csv"
OUT_DE = OUTPUT / f"{PREFIX}_ppar_low_high_de.csv"
OUT_PATHWAY = OUTPUT / f"{PREFIX}_pathway_state_scores.csv"
OUT_CORR = OUTPUT / f"{PREFIX}_state_correlations.csv"
OUT_DONOR = OUTPUT / f"{PREFIX}_donor_level_validation.csv"
OUT_SUBTYPE = OUTPUT / f"{PREFIX}_subtype_localization.csv"
OUT_INTERACTION = OUTPUT / f"{PREFIX}_tumor_ppar_interaction.csv"
OUT_REGULATOR = OUTPUT / f"{PREFIX}_regulator_activity.csv"
OUT_ANCHOR = OUTPUT / f"{PREFIX}_regulatory_anchor_ranking.csv"
OUT_EXTERNAL = OUTPUT / f"{PREFIX}_external_validation.csv"
OUT_BRIDGE = OUTPUT / f"{PREFIX}_bridge_evidence_table.csv"
OUT_PSEUDOBULK = OUTPUT / f"{PREFIX}_donor_state_pseudobulk.csv"
OUT_REPORT = OUTPUT / f"{PREFIX}_report.md"
OUT_MANIFEST = OUTPUT / f"{PREFIX}_manifest.json"


PPAR_NR_GENES = phase2g.PPAR_NR_GENES
ANCHOR_CANDIDATES = phase2g.ANCHOR_CANDIDATES


def bh(values: pd.Series | np.ndarray) -> np.ndarray:
    values = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    out = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if valid.any():
        out[valid] = multipletests(values[valid], method="fdr_bh")[1]
    return out


def zscore_columns(frame: pd.DataFrame, reference: pd.DataFrame | None = None) -> pd.DataFrame:
    reference = frame if reference is None else reference
    mean = reference.mean(axis=0)
    sd = reference.std(axis=0, ddof=1).replace(0, np.nan)
    return frame.subtract(mean, axis=1).divide(sd, axis=1)


def load_annotation() -> pd.DataFrame:
    annotation = pd.read_csv(ANNOTATION, sep="\t", index_col=0)
    annotation.index = annotation.index.astype(str)
    annotation["group"] = annotation["Class"].map(
        {"Tumor": "tumor", "Normal": "normal", "Border": "border"}
    )
    annotation["cell_type"] = annotation["Cell_type"].astype(str)
    annotation["cell_subtype"] = annotation["Cell_subtype"].astype(str)
    keep = annotation["cell_type"].eq("Epithelial cells") & annotation["group"].isin(
        ["tumor", "normal", "border"]
    )
    annotation = annotation.loc[keep, ["Patient", "group", "cell_type", "cell_subtype"]].copy()
    annotation["donor_id"] = annotation["Patient"].astype(str)
    return annotation


def read_target_matrix(cell_ids: list[str], wanted_genes: set[str], matrix_path: Path = MATRIX) -> tuple[pd.DataFrame, list[str]]:
    """Stream only selected gene rows from the wide GEO matrix."""
    wanted = {str(g).upper() for g in wanted_genes}
    rows: dict[str, list[np.ndarray]] = {}
    with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        matrix_cells_all = header[1:]
        matrix_index = {cell_id: i for i, cell_id in enumerate(matrix_cells_all)}
        matrix_cells = [cell_id for cell_id in cell_ids if cell_id in matrix_index]
        positions = [matrix_index[cell_id] for cell_id in matrix_cells]
        if not matrix_cells:
            raise ValueError("No annotation cell IDs were found in the expression matrix header.")
        print(
            f"[GSE144735 matrix] annotation cells={len(cell_ids):,}; matrix cells={len(matrix_cells_all):,}; aligned={len(matrix_cells):,}; missing_annotation={len(set(cell_ids) - set(matrix_cells)):,}",
            flush=True,
        )
        for line_number, line in enumerate(handle, start=1):
            first, sep, rest = line.rstrip("\r\n").partition("\t")
            if not sep or first.upper() not in wanted:
                continue
            values = np.fromstring(rest, sep="\t", dtype=float)
            if len(values) != len(matrix_cells_all):
                raise ValueError(f"Unexpected matrix width at row {line_number}: {first}")
            rows.setdefault(first.upper(), []).append(values[positions])
            if line_number % 5000 == 0:
                print(f"[GSE144735 matrix] scanned rows={line_number:,}; selected={len(rows):,}", flush=True)
    missing = sorted(wanted - set(rows))
    selected = {
        gene: np.nanmean(np.vstack(values), axis=0)
        for gene, values in rows.items()
    }
    expression = pd.DataFrame(selected, index=matrix_cells)
    expression.index.name = "cell_id"
    return expression, missing


def make_state_scores(expression: pd.DataFrame, annotation: pd.DataFrame, state_sets: dict[str, set[str]]) -> tuple[pd.DataFrame, dict[str, int]]:
    primary = annotation["group"].isin(["tumor", "normal"])
    reference = expression.loc[annotation.index[primary]]
    state_scores = pd.DataFrame(index=expression.index)
    sizes: dict[str, int] = {}
    for state, genes in state_sets.items():
        present = sorted(set(genes) & set(expression.columns))
        sizes[state] = len(present)
        if not present:
            state_scores[state] = np.nan
            continue
        state_scores[state] = zscore_columns(expression[present], reference[present]).mean(axis=1)
    return state_scores, sizes


def paired_summary(frame: pd.DataFrame, features: list[str], feature_type: str, level: str) -> pd.DataFrame:
    rows = []
    for feature in features:
        if feature not in frame.columns:
            continue
        pivot = frame.pivot_table(index="donor_id", columns="group", values=feature, aggfunc="mean")
        if not {"tumor", "normal"}.issubset(pivot.columns):
            continue
        delta = (pivot["tumor"] - pivot["normal"]).dropna()
        if len(delta) < 3:
            continue
        p = float(wilcoxon(delta).pvalue) if np.any(delta != 0) else 1.0
        rows.append({
            "feature_type": feature_type,
            "feature": feature,
            "analysis_level": level,
            "n_paired_donors": int(len(delta)),
            "mean_delta_tumor_minus_normal": float(delta.mean()),
            "median_delta_tumor_minus_normal": float(delta.median()),
            "p_value": p,
            "direction": "up" if delta.median() > 0 else "down" if delta.median() < 0 else "flat",
            "donor_consistency": float((delta > 0).mean()) if delta.median() > 0 else float((delta < 0).mean()),
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = bh(result["p_value"])
    return result


def make_units(expression: pd.DataFrame, state_scores: pd.DataFrame, annotation: pd.DataFrame, cell_scores: pd.DataFrame | None = None) -> pd.DataFrame:
    work = annotation.join(expression, how="inner").join(state_scores, how="left")
    if cell_scores is not None:
        work = work.join(cell_scores[["PPAR_NR_score", "PPAR_group_primary"]], how="left")
    group_cols = ["donor_id", "group"]
    numeric = work.select_dtypes(include=[np.number]).columns.tolist()
    units = work.groupby(group_cols, observed=True)[numeric].mean().reset_index()
    counts = work.groupby(group_cols, observed=True).size().reset_index(name="n_cells")
    return units.merge(counts, on=group_cols, validate="one_to_one")


def frozen_donor_ppar_score(raw_axis: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    """Reproduce Phase 2F's nine-gene donor-pseudobulk PPAR/NR score."""
    work = annotation.join(raw_axis, how="inner").loc[lambda x: x["group"].isin(["tumor", "normal"])]
    genes = [*PPAR_NR_GENES, "RELA", "STAT3"]
    summed = work.groupby(["donor_id", "group"], observed=True)[genes].sum()
    total = summed[genes].sum(axis=1).replace(0, np.nan)
    log_cpm = np.log1p(summed[genes].div(total, axis=0) * 1_000_000)
    z = zscore_columns(log_cpm, log_cpm)
    result = z[PPAR_NR_GENES].mean(axis=1).rename("PPAR_NR_score").reset_index()
    return result


def make_ppar_units(expression: pd.DataFrame, state_scores: pd.DataFrame, annotation: pd.DataFrame, cell_scores: pd.DataFrame) -> pd.DataFrame:
    work = annotation.join(expression, how="inner").join(state_scores, how="left")
    work = work.join(cell_scores[["PPAR_NR_score", "PPAR_group_primary"]], how="left")
    work = work.loc[work["group"].isin(["tumor", "normal"]) & work["PPAR_group_primary"].isin(["PPAR-low", "PPAR-high"])]
    numeric = work.select_dtypes(include=[np.number]).columns.tolist()
    units = work.groupby(["donor_id", "group", "PPAR_group_primary"], observed=True)[numeric].mean().reset_index()
    counts = work.groupby(["donor_id", "group", "PPAR_group_primary"], observed=True).size().reset_index(name="n_cells")
    return units.merge(counts, on=["donor_id", "group", "PPAR_group_primary"], validate="one_to_one")


def low_high_de(ppar_units: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    for feature in features:
        if feature not in ppar_units.columns:
            continue
        pivot = ppar_units.pivot_table(
            index=["donor_id", "group"], columns="PPAR_group_primary", values=feature, aggfunc="mean"
        )
        if not {"PPAR-low", "PPAR-high"}.issubset(pivot.columns):
            continue
        delta = (pivot["PPAR-low"] - pivot["PPAR-high"]).dropna()
        if len(delta) < 3:
            continue
        p = float(wilcoxon(delta).pvalue) if np.any(delta != 0) else 1.0
        rows.append({
            "feature": feature,
            "n_paired_donor_group_units": int(len(delta)),
            "mean_delta_low_minus_high": float(delta.mean()),
            "median_delta_low_minus_high": float(delta.median()),
            "P": p,
            "direction": "higher_in_PPAR_low" if delta.median() > 0 else "higher_in_PPAR_high" if delta.median() < 0 else "flat",
            "analysis_level": "donor-by-group pseudobulk; cell state defined by frozen PPAR/NR quartiles",
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = bh(result["P"])
    return result


def state_correlations(units: pd.DataFrame, states: list[str]) -> pd.DataFrame:
    rows = []
    for state in states:
        if state not in units.columns:
            continue
        for group, subset in [("pooled", units), ("tumor", units[units["group"].eq("tumor")]), ("normal", units[units["group"].eq("normal")])]:
            subset = subset[["PPAR_NR_score", state]].dropna()
            if len(subset) < 5:
                continue
            rho, p = spearmanr(subset["PPAR_NR_score"], subset[state])
            rows.append({"group": group, "state": state, "n_donors": len(subset), "spearman_rho": rho, "P": p,
                         "direction": "positive" if rho > 0 else "negative" if rho < 0 else "flat"})
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = result.groupby("group", group_keys=False)["P"].transform(bh)
    return result


def tumor_ppar_interaction(ppar_units: pd.DataFrame, states: list[str]) -> pd.DataFrame:
    rows = []
    for state in states:
        if state not in ppar_units.columns:
            continue
        work = ppar_units[["donor_id", "group", "PPAR_group_primary", state]].dropna().copy()
        if work["donor_id"].nunique() < 4 or work["group"].nunique() < 2:
            continue
        work["tumor"] = work["group"].eq("tumor").astype(int)
        work["ppar_low"] = work["PPAR_group_primary"].eq("PPAR-low").astype(int)
        try:
            fit = ols(f"Q('{state}') ~ C(donor_id) + tumor + ppar_low + tumor:ppar_low", data=work).fit()
            interaction = next((name for name in fit.params.index if ":" in name), None)
            if interaction is None:
                continue
            rows.append({"state": state, "n_donors": int(work["donor_id"].nunique()),
                         "interaction_beta_tumor_x_PPARlow": float(fit.params[interaction]),
                         "P": float(fit.pvalues[interaction]),
                         "analysis_level": "donor-by-group PPAR-state pseudobulk OLS with donor fixed effects"})
        except Exception as exc:
            rows.append({"state": state, "n_donors": int(work["donor_id"].nunique()),
                         "interaction_beta_tumor_x_PPARlow": np.nan, "P": np.nan,
                         "analysis_level": f"not_estimable: {type(exc).__name__}"})
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = bh(result["P"])
    return result


def prepare_network() -> pd.DataFrame:
    raw = pd.read_csv(DOROTHEA_RAW, sep="\t")
    raw = raw.loc[raw["source_genesymbol"].isin(ANCHOR_CANDIDATES) & raw["dorothea_level"].isin(["A", "B", "C"])].copy()
    raw = raw.drop_duplicates(["source_genesymbol", "dorothea_level", "target_genesymbol"])
    def mor(row: pd.Series) -> int:
        if bool(row["is_stimulation"]) and bool(row["is_inhibition"]):
            return 1 if bool(row["consensus_stimulation"]) else -1
        if bool(row["is_stimulation"]):
            return 1
        if bool(row["is_inhibition"]):
            return -1
        return 1
    raw["mor"] = raw.apply(mor, axis=1)
    weights = {"A": 1.0, "B": 2.0, "C": 3.0}
    raw["weight"] = raw.apply(lambda row: row["mor"] / weights[row["dorothea_level"]], axis=1)
    return raw[["source_genesymbol", "target_genesymbol", "weight", "dorothea_level"]].rename(
        columns={"source_genesymbol": "source", "target_genesymbol": "target", "dorothea_level": "confidence"}
    )


def ulm_activity(expression_units: pd.DataFrame, ppar_units: pd.DataFrame, network: pd.DataFrame) -> pd.DataFrame:
    """Run the standard one-regulator weighted linear model per sample."""
    rows = []
    def run_matrix(frame: pd.DataFrame, comparison: str) -> None:
        for regulator, net in network.groupby("source"):
            targets = [g for g in net["target"].astype(str) if g in expression_units.columns]
            if len(targets) < 5:
                continue
            weights = net.set_index("target").loc[targets, "weight"].to_numpy(float)
            x = weights - weights.mean()
            denom = float(np.sum(x * x))
            if denom <= 0:
                continue
            activities = []
            for _, row in frame.iterrows():
                y = pd.to_numeric(row[targets], errors="coerce").to_numpy(float)
                valid = np.isfinite(y) & np.isfinite(x)
                if valid.sum() < 5:
                    activities.append((np.nan, np.nan))
                    continue
                xv, yv = x[valid], y[valid]
                xv = xv - xv.mean()
                beta = float(np.sum(xv * (yv - yv.mean())) / np.sum(xv * xv))
                fitted = yv.mean() + beta * xv
                resid = yv - fitted
                df = max(len(yv) - 2, 1)
                se = float(np.sqrt(np.sum(resid * resid) / df / np.sum(xv * xv))) if np.sum(xv * xv) else np.nan
                tstat = beta / se if se and np.isfinite(se) else np.nan
                p = float(2 * t_dist.sf(abs(tstat), df)) if np.isfinite(tstat) else np.nan
                activities.append((beta, p))
            temp = frame[["donor_id", "group"]].copy()
            temp["activity"] = [x[0] for x in activities]
            if comparison == "tumor_vs_normal":
                pivot = temp.pivot_table(index="donor_id", columns="group", values="activity", aggfunc="mean")
                if not {"tumor", "normal"}.issubset(pivot.columns):
                    continue
                delta = (pivot["tumor"] - pivot["normal"]).dropna()
            else:
                temp = frame[["donor_id", "group", "PPAR_group_primary"]].copy()
                temp["activity"] = [x[0] for x in activities]
                pivot = temp.pivot_table(index=["donor_id", "group"], columns="PPAR_group_primary", values="activity", aggfunc="mean")
                if not {"PPAR-low", "PPAR-high"}.issubset(pivot.columns):
                    continue
                delta = (pivot["PPAR-low"] - pivot["PPAR-high"]).dropna()
            if len(delta) < 3:
                continue
            p = float(wilcoxon(delta).pvalue) if np.any(delta != 0) else 1.0
            rows.append({"regulator": regulator, "comparison": comparison, "n_pairs": len(delta),
                         "activity_delta": float(delta.mean()), "median_activity_delta": float(delta.median()),
                         "P": p, "direction": "up_in_first_group" if delta.mean() > 0 else "down_in_first_group",
                         "donor_consistency": float((delta > 0).mean()) if delta.mean() > 0 else float((delta < 0).mean()),
                         "target_count": len(targets), "method": "DoRothEA A-C weighted ULM; local implementation of the decoupler ULM model"})
    run_matrix(expression_units, "tumor_vs_normal")
    if not ppar_units.empty:
        run_matrix(ppar_units, "PPAR_low_vs_high")
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = result.groupby("comparison", group_keys=False)["P"].transform(bh)
    return result


def external_table(validation: pd.DataFrame, states: list[str], anchors: list[str]) -> pd.DataFrame:
    keep = validation.loc[validation["feature"].isin(["PPAR_NR_score", *states, *anchors])].copy()
    if keep.empty:
        return pd.DataFrame()
    keep = keep.rename(columns={"n_paired_donors": "n_paired", "p_value": "P"})
    keep["dataset"] = "GSE144735"
    keep["replication_status"] = np.where(keep["feature"].eq("PPAR_NR_score"), "alternative primary; six paired patients", "alternative primary target-universe support")
    keep["analysis_level"] = "patient-paired epithelial donor means from complete processed log-TPM matrix"
    return keep[["dataset", "feature", "n_paired", "median_delta_tumor_minus_normal", "mean_delta_tumor_minus_normal", "P", "direction", "replication_status", "analysis_level"]]


def subtype_localization(cell_scores: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    work = annotation.join(cell_scores[["PPAR_NR_score", "PPAR_group_primary"]], how="inner")
    rows = []
    for subtype, subset in work.groupby("cell_subtype", observed=True):
        rows.append({"subtype": subtype, "cell_type": "Epithelial cells", "n_cells": len(subset),
                     "n_donors": subset["donor_id"].nunique(), "median_PPAR_NR": float(subset["PPAR_NR_score"].median()),
                     "PPAR_low_fraction": float(subset["PPAR_group_primary"].eq("PPAR-low").mean()),
                     "PPAR_high_fraction": float(subset["PPAR_group_primary"].eq("PPAR-high").mean()),
                     "tumor_normal_distribution": ";".join(sorted(subset["group"].unique())),
                     "annotation_source": "GSE144735 author-provided epithelial subtype label"})
    return pd.DataFrame(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if not MATRIX.exists() or not RAW_CORE_MATRIX.exists() or not ANNOTATION.exists() or not DOROTHEA_RAW.exists():
        raise SystemExit("GSE144735 full matrix, frozen raw-core matrix, annotation, and DoRothEA network are required.")
    state_sets, state_sources = phase2g.build_state_sets()
    network = prepare_network()
    query_genes = set().union(*state_sets.values()) | set(PPAR_NR_GENES) | set(ANCHOR_CANDIDATES) | set(network["target"].astype(str))
    annotation = load_annotation()
    cell_ids = annotation.index.astype(str).tolist()
    expression, missing = read_target_matrix(cell_ids, query_genes, MATRIX)
    expression = expression.reindex(cell_ids)
    annotation = annotation.reindex(expression.index)
    state_cell, state_sizes = make_state_scores(expression, annotation, state_sets)
    primary = annotation["group"].isin(["tumor", "normal"])
    raw_core, raw_missing = read_target_matrix(cell_ids, set(PPAR_NR_GENES), RAW_CORE_MATRIX)
    raw_core = raw_core.reindex(expression.index)
    core_total = raw_core.sum(axis=1).replace(0, np.nan)
    core_log_cpm = np.log1p(raw_core.div(core_total, axis=0) * 1_000_000)
    core_z = zscore_columns(core_log_cpm, core_log_cpm.loc[primary])
    cell_scores = pd.DataFrame(index=expression.index)
    cell_scores["PPAR_NR_score"] = core_z.mean(axis=1)
    ranks = cell_scores.loc[primary, "PPAR_NR_score"].rank(method="first")
    qlabels = pd.qcut(ranks, 4, labels=["PPAR-low", "Q2", "Q3", "PPAR-high"]).astype(str)
    cell_scores["PPAR_group_primary"] = "not_primary"
    cell_scores.loc[primary, "PPAR_group_primary"] = qlabels.to_numpy()
    cell_scores["PPAR_group"] = cell_scores["PPAR_group_primary"]
    cell_scores["tumor_normal"] = annotation["group"]
    cell_scores["donor_id"] = annotation["donor_id"]
    cell_scores["cell_type"] = annotation["cell_type"]
    cell_scores["cell_subtype"] = annotation["cell_subtype"]
    cell_scores["analysis_level"] = "cell-level score for state definition; inference is donor-level"
    for gene in PPAR_NR_GENES:
        if gene in raw_core.columns:
            cell_scores[gene] = raw_core[gene]
    for state in state_cell.columns:
        cell_scores[state] = state_cell[state]
    cell_scores.reset_index(names="cell_id").to_csv(OUT_CELL, index=False)

    donor_units = make_units(expression, state_cell, annotation, cell_scores)
    raw_axis, raw_axis_missing = read_target_matrix(cell_ids, set([*PPAR_NR_GENES, "RELA", "STAT3"]), RAW_CORE_MATRIX)
    frozen_scores = frozen_donor_ppar_score(raw_axis, annotation)
    donor_units = donor_units.drop(columns=["PPAR_NR_score"], errors="ignore").merge(
        frozen_scores, on=["donor_id", "group"], how="left", validate="one_to_one"
    )
    ppar_units = make_ppar_units(expression, state_cell, annotation, cell_scores)
    state_names = list(state_sets)
    anchor_present = [g for g in ANCHOR_CANDIDATES if g in expression.columns]
    validation = pd.concat([
        paired_summary(donor_units, ["PPAR_NR_score"], "state", "GSE144735 donor-level paired epithelial analysis"),
        paired_summary(donor_units, state_names, "state", "GSE144735 donor-level paired epithelial analysis"),
        paired_summary(donor_units, anchor_present, "gene", "GSE144735 donor-level paired epithelial analysis"),
    ], ignore_index=True)
    validation.to_csv(OUT_DONOR, index=False)

    de = low_high_de(ppar_units, sorted(set(state_names) | set(anchor_present) | set(PPAR_NR_GENES)))
    de.to_csv(OUT_DE, index=False)
    state_cell.loc[:, [x for x in state_names if x in state_cell.columns]].assign(
        donor_id=annotation["donor_id"].to_numpy(), group=annotation["group"].to_numpy(), PPAR_NR_score=cell_scores["PPAR_NR_score"].to_numpy()
    ).loc[lambda x: x["group"].isin(["tumor", "normal"])].groupby(["donor_id", "group"], observed=True).mean(numeric_only=True).reset_index().melt(
        id_vars=["donor_id", "group"], var_name="state", value_name="state_score"
    ).to_csv(OUT_PATHWAY, index=False)
    corr = state_correlations(donor_units, state_names)
    corr.to_csv(OUT_CORR, index=False)
    interaction = tumor_ppar_interaction(ppar_units, state_names)
    interaction.to_csv(OUT_INTERACTION, index=False)
    subtype_localization(cell_scores, annotation).to_csv(OUT_SUBTYPE, index=False)

    regulator = ulm_activity(donor_units, ppar_units, network)
    regulator.to_csv(OUT_REGULATOR, index=False)
    external = external_table(validation, state_names, anchor_present)
    tcga_path = OUTPUT / "mcop_phase2f_tcga_paired_gene_stats.csv"
    if tcga_path.exists():
        tcga = pd.read_csv(tcga_path)
        rows = []
        for _, row in tcga.loc[tcga["gene"].isin(anchor_present)].iterrows():
            rows.append({"dataset": "TCGA paired primary vs solid normal", "feature": row["gene"], "n_paired": row["paired_n"],
                         "median_delta_tumor_minus_normal": row["median_delta_tumor_minus_normal"], "mean_delta_tumor_minus_normal": np.nan,
                         "P": row["p_value"], "direction": "up" if row["median_delta_tumor_minus_normal"] > 0 else "down",
                         "replication_status": "external directional bulk support", "analysis_level": "existing Phase 2F paired bulk output"})
        if rows:
            external = pd.concat([external, pd.DataFrame(rows)], ignore_index=True)
    if not external.empty:
        external["BH_FDR"] = bh(external["P"])
    external.to_csv(OUT_EXTERNAL, index=False)

    ranking, bridge = phase2g.build_anchor_outputs(ANCHOR_CANDIDATES, validation, corr, regulator, external)
    ranking.to_csv(OUT_ANCHOR, index=False)
    bridge.to_csv(OUT_BRIDGE, index=False)
    donor_units.to_csv(OUT_PSEUDOBULK, index=False)

    ppar_row = validation.loc[validation["feature"].eq("PPAR_NR_score")]
    ppar_text = ppar_row.to_dict("records") if not ppar_row.empty else "not estimable"
    state_rows = validation.loc[validation["feature_type"].eq("state")].copy()
    if not state_rows.empty:
        state_rows["abs_delta"] = pd.to_numeric(state_rows["median_delta_tumor_minus_normal"], errors="coerce").abs()
        top_state = state_rows.sort_values("abs_delta", ascending=False).iloc[0]
        top_state_text = f"{top_state['feature']} (median Δ={float(top_state['median_delta_tumor_minus_normal']):.3f}; BH-FDR={float(top_state['BH_FDR']):.3f})"
    else:
        top_state_text = "not estimable"
    min_interaction_fdr = float(pd.to_numeric(interaction.get("BH_FDR", pd.Series(dtype=float)), errors="coerce").min()) if not interaction.empty else np.nan
    min_regulator_fdr = float(pd.to_numeric(regulator.get("BH_FDR", pd.Series(dtype=float)), errors="coerce").min()) if not regulator.empty else np.nan
    report = [
        "# Phase 2G alternative primary — GSE144735 complete-matrix epithelial state analysis",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Data path",
        "",
        "- Primary alternative dataset: GSE144735 GEO processed natural-log TPM matrix; live CELLxGENE Census/TileDB raw query was not used.",
        f"- Cells in epithelial annotation: **{len(annotation):,}**; primary tumor/normal epithelial cells: **{int(primary.sum()):,}**; paired donors: **{annotation.loc[primary, 'donor_id'].nunique()}**.",
        "- The matrix was streamed row-wise and only the frozen state/regulator gene universe was retained; the complete raw/processed matrix is not committed to the repository.",
        "- Primary contrast excludes Border epithelial cells; Border cells are retained only in the cell audit and subtype localization.",
        "",
        "## Frozen boundaries",
        "",
        f"- PPAR/NR core: **{', '.join(PPAR_NR_GENES)}**; no result-driven gene editing.",
        f"- Unbiased state universe: **{len(state_sets)}** programs; present gene counts: `{json.dumps(state_sizes, ensure_ascii=False)}`.",
        "- Cell-level PPAR-low/high is defined by the bottom/top quartile of the fixed 7-gene PPAR/NR score; donor-level `PPAR_NR_score` validation reuses the frozen Phase 2F nine-gene pseudobulk score (7 core genes plus RELA/STAT3 denominator).",
        "- State and regulator expression comes from the complete processed log-TPM matrix; donor-level inference uses donor-by-group mean expression because GEO supplies processed expression rather than raw full-matrix counts.",
        "- Targeted state/regulator expression is not a genome-wide DE claim.",
        "",
        "## Primary alternative results",
        "",
        f"- PPAR/NR paired result: `{ppar_text}`.",
        f"- State programs with estimable tumor-normal paired summaries: **{int(validation['feature_type'].eq('state').sum()) if not validation.empty else 0}**.",
        f"- Largest absolute state contrast: **{top_state_text}**; no state is called discovery-positive without multiplicity control.",
        f"- DoRothEA regulator summaries: **{int(regulator['regulator'].nunique()) if not regulator.empty and 'regulator' in regulator.columns else 0}** regulators; A-C weighted ULM with minimum five observed targets.",
        f"- Minimum regulator-activity BH-FDR across tested summaries: **{min_regulator_fdr:.3f}**; minimum tumor×PPAR interaction BH-FDR: **{min_interaction_fdr:.3f}**. These small-sample values are descriptive, not mechanistic confirmation.",
        "- This alternative primary analysis completes the prespecified state, interaction, subtype, and regulator modules for one independent dataset, but the six-patient sample is underpowered for stable inference and cannot supply Census-style leave-one-dataset-out validation.",
        "",
        "## Interpretation",
        "",
        "The result tests whether DINP-axis molecular candidates converge on a CRC epithelial state; it does not establish DINP/MCOP causality or mediation. The Phase 2F Census epithelial/myeloid compartment result remains the prior validated support, and the myeloid opposite-direction result is not erased.",
        "",
        "## Verdict",
        "",
        "**PARTIALLY — alternative primary analysis completed, but independent-dataset stability and causal exposure-to-state direction remain unresolved.**",
        "",
        "## Outputs",
        "",
        *[f"- `{path.name}`" for path in [OUT_CELL, OUT_DE, OUT_PATHWAY, OUT_CORR, OUT_DONOR, OUT_SUBTYPE, OUT_INTERACTION, OUT_REGULATOR, OUT_ANCHOR, OUT_EXTERNAL, OUT_BRIDGE, OUT_PSEUDOBULK]],
    ]
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
    manifest = {
        "analysis": "Phase 2G DINP/MCOP-CRC epithelial state convergence — GSE144735 alternative primary",
        "dataset": "GSE144735",
        "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735",
        "matrix_url": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144735&file=GSE144735_processed_KUL3_CRC_10X_natural_log_TPM_matrix.txt.gz&format=file",
        "annotation_url": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144735&file=GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz&format=file",
        "matrix_file": MATRIX.name,
        "frozen_core_matrix_file": RAW_CORE_MATRIX.name,
        "annotation_file": ANNOTATION.name,
        "census_runtime_status": "live TileDB raw query bypassed; prior Census attempt remains in mcop_phase2g_manifest.json",
        "ppar_nr_core_genes": PPAR_NR_GENES,
        "state_programs": {name: {"source": state_sources[name], "n_genes": len(genes)} for name, genes in state_sets.items()},
        "n_epithelial_cells": int(len(annotation)),
        "n_primary_epithelial_cells": int(primary.sum()),
        "n_paired_donors": int(annotation.loc[primary, "donor_id"].nunique()),
        "unit_of_inference": "donor-level paired summaries; cell-level quartiles only define state labels",
        "cell_ppar_score_definition": "frozen 7-gene raw-UMI log-CPM score with gene-wise z-scoring",
        "donor_ppar_score_definition": "frozen Phase 2F nine-gene donor-pseudobulk score; PPAR/NR mean over seven core genes",
        "state_expression_measure": "processed natural-log TPM; donor-by-group mean expression",
        "regulator_method": "DoRothEA A-C weighted ULM, local implementation matching the decoupler ULM model",
        "missing_query_genes": missing,
        "missing_frozen_core_genes": raw_missing,
        "missing_frozen_axis_genes": raw_axis_missing,
        "causal_status": "not established",
        "verdict": "PARTIALLY",
        "outputs": [str(x) for x in [OUT_CELL, OUT_DE, OUT_PATHWAY, OUT_CORR, OUT_DONOR, OUT_SUBTYPE, OUT_INTERACTION, OUT_REGULATOR, OUT_ANCHOR, OUT_EXTERNAL, OUT_BRIDGE, OUT_PSEUDOBULK, OUT_REPORT]],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"dataset": "GSE144735", "n_epithelial_cells": len(annotation), "n_paired_donors": int(annotation.loc[primary, "donor_id"].nunique()), "n_states": len(state_sets), "n_query_genes": len(query_genes), "missing_genes": len(missing), "verdict": "PARTIALLY"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
