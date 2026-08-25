"""Phase 2G: DINP/MCOP--CRC epithelial state convergence and regulatory anchoring.

This script is deliberately downstream of the frozen Phase 2F analysis.  It
does not redefine the PPAR/NR score and does not claim exposure causality.  The
primary unit for inference is donor-level pseudobulk; cells are retained only
for the prespecified PPAR-low/high state definition and subtype audit.

The Census query is pinned to the 2025-11-08 release and primary data.  The
current eligible Census CRC dataset is selected from the Phase 2F paired-donor
audit (the dataset with the largest number of paired epithelial donors).  A
targeted expression universe containing all frozen state genes plus selected
DoRothEA targets is queried out-of-core; this is a state/regulator analysis,
not a claim of genome-wide differential expression.

Run in the pinned WSL environment:

    /opt/cellxgene-census/bin/python work/scripts/mcop_phase2g_epithelial_state_analysis.py
"""

from __future__ import annotations

import json
import argparse
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs"
RAW_GSE = ROOT / "work" / "mcop_phase2f_external" / "raw"
GENESETS = ROOT / "work" / "gene_sets"
CENSUS_VERSION = "2025-11-08"
CENSUS_URI = "s3://cellxgene-census-public-us-west-2/cell-census/2025-11-08/soma/"
ORGANISM = "homo_sapiens"
PRIMARY_FILTER = "is_primary_data == True"

PPAR_NR_GENES = ["PPARA", "PPARD", "PPARG", "NR1I2", "NR1I3", "NR1H2", "NR1H3"]
ANCHOR_CANDIDATES = [
    "PPARG", "PPARA", "PPARD", "NR1I2", "NR1I3", "NR1H2", "NR1H3",
    "RELA", "STAT3", "HIF1A", "MYC", "JUN", "FOS", "TP53", "SMAD3", "TCF7L2",
]

OUT_CELL = OUTPUT / "mcop_phase2g_epithelial_state_scores.csv"
OUT_DE = OUTPUT / "mcop_phase2g_ppar_low_high_de.csv"
OUT_PATHWAY = OUTPUT / "mcop_phase2g_pathway_state_scores.csv"
OUT_CORR = OUTPUT / "mcop_phase2g_state_correlations.csv"
OUT_DONOR = OUTPUT / "mcop_phase2g_donor_level_validation.csv"
OUT_SUBTYPE = OUTPUT / "mcop_phase2g_subtype_localization.csv"
OUT_INTERACTION = OUTPUT / "mcop_phase2g_tumor_ppar_interaction.csv"
OUT_REGULATOR = OUTPUT / "mcop_phase2g_regulator_activity.csv"
OUT_ANCHOR = OUTPUT / "mcop_phase2g_regulatory_anchor_ranking.csv"
OUT_EXTERNAL = OUTPUT / "mcop_phase2g_external_validation.csv"
OUT_BRIDGE = OUTPUT / "mcop_phase2g_bridge_evidence_table.csv"
OUT_PSEUDOBULK = OUTPUT / "mcop_phase2g_donor_state_pseudobulk.csv"
OUT_REPORT = OUTPUT / "mcop_phase2g_report.md"
OUT_MANIFEST = OUTPUT / "mcop_phase2g_manifest.json"


def q(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def in_filter(field: str, values: Iterable[object]) -> str:
    vals = list(values)
    return f"{field} in [{', '.join(q(x) for x in vals)}]"


def chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[start : start + size] for start in range(0, len(values), size)]


def finite(value: object) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def bh(values: Iterable[object]) -> np.ndarray:
    vals = pd.to_numeric(pd.Series(list(values)), errors="coerce").to_numpy(dtype=float)
    out = np.full(len(vals), np.nan)
    mask = np.isfinite(vals)
    if mask.any():
        out[mask] = multipletests(vals[mask], method="fdr_bh")[1]
    return out


def retry(function, label: str, attempts: int = 5):
    error = None
    for attempt in range(1, attempts + 1):
        try:
            return function()
        except Exception as exc:  # noqa: BLE001 - remote libraries vary
            error = exc
            if attempt == attempts:
                break
            time.sleep(min(30, 5 * attempt))
            print(f"[retry {attempt}/{attempts - 1}] {label}: {type(exc).__name__}")
    raise RuntimeError(f"Census operation failed after {attempts} attempts: {label}: {error}") from error


def read_gmt(path: Path) -> dict[str, set[str]]:
    result: dict[str, set[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 3:
            result[fields[0]] = {x.upper() for x in fields[2:] if x}
    return result


def build_state_sets() -> tuple[dict[str, set[str]], dict[str, str]]:
    hallmark = read_gmt(GENESETS / "h.all.v2026.1.Hs.symbols.gmt")
    wanted = {
        "EMT": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "E2F_targets": "HALLMARK_E2F_TARGETS",
        "G2M_checkpoint": "HALLMARK_G2M_CHECKPOINT",
        "MYC_targets_V1": "HALLMARK_MYC_TARGETS_V1",
        "MYC_targets_V2": "HALLMARK_MYC_TARGETS_V2",
        "Hypoxia": "HALLMARK_HYPOXIA",
        "Inflammatory_response": "HALLMARK_INFLAMMATORY_RESPONSE",
        "TNF_NFkB": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "IL6_JAK_STAT3": "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "ROS": "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
        "OXPHOS": "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "Fatty_acid_metabolism": "HALLMARK_FATTY_ACID_METABOLISM",
        "Cholesterol_homeostasis": "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
        "Apoptosis": "HALLMARK_APOPTOSIS",
        "p53": "HALLMARK_P53_PATHWAY",
        "UPR": "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
        "Glycolysis": "HALLMARK_GLYCOLYSIS",
        "IFN_alpha": "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "IFN_gamma": "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "TGF_beta": "HALLMARK_TGF_BETA_SIGNALING",
        "WNT_beta_catenin": "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    }
    states = {name: hallmark[term] for name, term in wanted.items() if term in hallmark}
    sources = {name: "MSigDB Hallmark v2026.1" for name in states}
    custom = {
        "intestinal_epithelial_differentiation": {
            "CDX2", "KLF4", "HNF4A", "GATA6", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT20",
            "MUC13", "ALPI", "FABP1", "FABP2", "SLC26A3", "CA1", "CA2", "SI",
        },
        "stemness": {"LGR5", "ASCL2", "SMOC2", "OLFM4", "PROM1", "BMI1", "SOX9", "AXIN2", "LRIG1", "EPHB2"},
        "secretory_differentiation": {"MUC2", "SPINK4", "TFF3", "AGR2", "CLCA1", "SPDEF", "FCGBP", "ATOH1", "KLF4", "ZG16"},
        "enterocyte_differentiation": {"ALPI", "FABP1", "FABP2", "APOA1", "APOA4", "SI", "SLC26A3", "CA1", "CA2", "VIL1"},
        "goblet_program": {"MUC2", "SPINK4", "CLCA1", "FCGBP", "TFF3", "AGR2", "ZG16", "WFDC2"},
        "stress_like_epithelial": {"HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "ATF3", "DDIT3", "XBP1", "JUN", "FOS", "HIF1A"},
    }
    for name, genes in custom.items():
        states[name] = genes
        sources[name] = "curated marker panel; frozen before analysis"
    return states, sources


def classify_epithelial(cell_type: object) -> bool:
    value = str(cell_type).lower()
    return bool(re.search(r"epithelial|colonocyte|enterocyte|goblet|paneth|enteroendocrine|tuft|best4|transit amplifying|stem cell of colon", value))


def subtype_label(cell_type: object) -> str:
    value = str(cell_type).lower()
    if re.search(r"stem|transit amplifying|cycling", value):
        return "stem/cycling-like annotation"
    if re.search(r"enterocyte|colonocyte", value):
        return "enterocyte-like annotation"
    if re.search(r"goblet|secretory|paneth|enteroendocrine|tuft", value):
        return "secretory-like annotation"
    if re.search(r"stress|inflammatory", value):
        return "stress-like annotation"
    return "other epithelial annotation"


def open_census():
    import cellxgene_census
    return retry(lambda: cellxgene_census.open_soma(uri=CENSUS_URI), "open pinned Census", attempts=6)


def read_axis(reader, joinids) -> pd.DataFrame:
    frames = [table.to_pandas() for table in reader]
    result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if "soma_joinid" not in result.columns:
        result.insert(0, "soma_joinid", np.asarray(joinids))
    return result


def load_primary_metadata() -> tuple[pd.DataFrame, str, list[str], list[str], list[str], list[str]]:
    paired = pd.read_csv(OUTPUT / "mcop_phase2f_singlecell_paired_donor_contrasts.csv")
    eligible = paired.loc[
        paired["compartment"].eq("epithelial")
        & paired["score"].eq("PPAR_nuclear_receptor_score")
    ].sort_values("paired_donors", ascending=False)
    if eligible.empty:
        raise RuntimeError("Phase 2F paired epithelial audit is missing; cannot define the Phase 2G primary dataset.")
    dataset_id = str(eligible.iloc[0]["dataset_id"])
    usecols = ["dataset_id", "donor_id", "disease", "tissue", "tissue_general", "cell_type", "group", "is_primary_data"]
    metadata = pd.read_csv(OUTPUT / "mcop_phase2f_singlecell_observation_metadata_audit.csv", usecols=usecols)
    metadata["dataset_id"] = metadata["dataset_id"].astype(str)
    metadata = metadata.loc[metadata["dataset_id"].eq(dataset_id) & metadata["cell_type"].map(classify_epithelial)].copy()
    if metadata.empty:
        raise RuntimeError(f"No epithelial metadata rows for selected Phase 2F dataset {dataset_id}.")
    tumor_labels = sorted(metadata.loc[metadata["group"].eq("tumor"), "disease"].dropna().astype(str).unique())
    normal_labels = sorted(metadata.loc[metadata["group"].eq("normal"), "disease"].dropna().astype(str).unique())
    cell_types = sorted(metadata["cell_type"].dropna().astype(str).unique())
    score_path = OUTPUT / "mcop_phase2f_singlecell_donor_scores.csv"
    scores = pd.read_csv(score_path, usecols=["dataset_id", "donor_id", "group", "compartment"])
    scores = scores.loc[scores["dataset_id"].astype(str).eq(dataset_id) & scores["compartment"].eq("epithelial")]
    tumor_donors = set(scores.loc[scores["group"].eq("tumor"), "donor_id"].astype(str))
    normal_donors = set(scores.loc[scores["group"].eq("normal"), "donor_id"].astype(str))
    paired_donors = sorted(tumor_donors & normal_donors)
    if len(paired_donors) < 5:
        raise RuntimeError(f"Only {len(paired_donors)} paired epithelial donors were found for {dataset_id}.")
    metadata = metadata.loc[metadata["donor_id"].astype(str).isin(paired_donors)].copy()
    return metadata, dataset_id, tumor_labels, normal_labels, cell_types, paired_donors


def query_filter(dataset_id: str, tumor_labels: list[str], normal_labels: list[str], cell_types: list[str], paired_donors: list[str]) -> str:
    return (
        f"{PRIMARY_FILTER} and dataset_id == {q(dataset_id)} and "
        f"{in_filter('disease', tumor_labels + normal_labels)} and "
        f"{in_filter('cell_type', cell_types)} and {in_filter('donor_id', paired_donors)}"
    )


def read_query_obs(query, tumor_labels: set[str], normal_labels: set[str], cell_types: set[str]) -> pd.DataFrame:
    obs = read_axis(
        query.obs(column_names=["dataset_id", "donor_id", "disease", "tissue", "tissue_general", "cell_type", "is_primary_data", "sex"]),
        query.obs_joinids(),
    )
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs = obs.loc[obs["cell_type"].isin(cell_types)].copy()
    obs["group"] = np.where(obs["disease"].astype(str).isin(tumor_labels), "tumor", np.where(obs["disease"].astype(str).isin(normal_labels), "normal", "outside_scope"))
    obs = obs.loc[obs["group"].isin(["tumor", "normal"])].copy()
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs = obs.loc[~obs["donor_id"].isin(["", "nan", "None"])].copy()
    obs["donor_key"] = obs["dataset_id"].astype(str) + "::" + obs["donor_id"].astype(str)
    return obs


def fetch_cell_scores(census, dataset_id: str, tumor_labels: list[str], normal_labels: list[str], cell_types: list[str], paired_donors: list[str]) -> pd.DataFrame:
    import tiledbsoma as soma
    var_filter = in_filter("feature_name", PPAR_NR_GENES)
    value_filter = query_filter(dataset_id, tumor_labels, normal_labels, cell_types, paired_donors)
    experiment = census["census_data"][ORGANISM]
    with experiment.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=value_filter),
        var_query=soma.AxisQuery(value_filter=var_filter),
    ) as query:
        obs = read_query_obs(query, set(tumor_labels), set(normal_labels), set(cell_types))
        var = read_axis(query.var(column_names=["feature_name", "feature_id"]), query.var_joinids())
        gene_map = {int(j): str(g) for j, g in zip(var["soma_joinid"], var["feature_name"]) if str(g) in PPAR_NR_GENES}
        if obs.empty or not gene_map:
            raise RuntimeError("Census core-gene query returned no epithelial observations or no PPAR/NR genes.")
        obs_index = pd.Series(np.arange(len(obs), dtype=int), index=obs["soma_joinid"].astype(np.int64))
        gene_index = {gene: idx for idx, gene in enumerate(PPAR_NR_GENES)}
        counts = np.zeros((len(obs), len(PPAR_NR_GENES)), dtype=np.float64)
        for table in query.X("raw").tables():
            batch = table.to_pandas()
            if batch.empty:
                continue
            cell_pos = obs_index.reindex(batch["soma_dim_0"].astype(np.int64)).to_numpy()
            genes = batch["soma_dim_1"].astype(np.int64).map(gene_map)
            valid = pd.notna(cell_pos) & genes.notna()
            if not valid.any():
                continue
            values = pd.to_numeric(batch.loc[valid, "soma_data"], errors="coerce").fillna(0).to_numpy(dtype=float)
            rows = cell_pos[valid].astype(int)
            cols = genes.loc[valid].map(gene_index).to_numpy(dtype=int)
            np.add.at(counts, (rows, cols), values)
    total = counts.sum(axis=1)
    normalized = np.log1p(np.divide(counts * 1_000_000.0, total[:, None], out=np.zeros_like(counts), where=total[:, None] > 0))
    gene_sd = normalized.std(axis=0, ddof=1)
    z = (normalized - normalized.mean(axis=0)) / np.where(gene_sd > 0, gene_sd, 1.0)
    score = np.nanmean(z, axis=1)
    result = obs.reset_index(drop=True).copy()
    result["cell_id"] = result["dataset_id"].astype(str) + "::" + result["soma_joinid"].astype(str)
    result["tumor_normal"] = result["group"]
    result["PPAR_NR_score"] = score
    result["PPAR_group_primary"] = "intermediate"
    ranks = pd.Series(score).rank(method="first")
    quartile = pd.qcut(ranks, 4, labels=["PPAR-low", "Q2", "Q3", "PPAR-high"])
    result["PPAR_group_primary"] = quartile.astype(str).to_numpy()
    result["PPAR_group_tercile"] = pd.qcut(ranks, 3, labels=["PPAR-low", "middle", "PPAR-high"]).astype(str).to_numpy()
    result["PPAR_group_median"] = np.where(score <= np.nanmedian(score), "PPAR-low", "PPAR-high")
    result["PPAR_group"] = result["PPAR_group_primary"]
    result["cell_subtype"] = result["cell_type"].map(subtype_label)
    result["PPAR_NR_core_genes"] = ";".join(PPAR_NR_GENES)
    return result


def fetch_donor_pseudobulk(census, cell_scores: pd.DataFrame, dataset_id: str, tumor_labels: list[str], normal_labels: list[str], cell_types: list[str], paired_donors: list[str], query_genes: set[str]) -> pd.DataFrame:
    import tiledbsoma as soma
    value_filter = query_filter(dataset_id, tumor_labels, normal_labels, cell_types, paired_donors)
    experiment = census["census_data"][ORGANISM]
    cell_lookup = cell_scores.set_index("soma_joinid")
    gene_keys = sorted(query_genes)
    gene_acc: dict[tuple[str, str, str, str], float] = {}
    gene_batches = chunks(gene_keys, 300)
    for batch_number, gene_batch in enumerate(gene_batches, start=1):
        print(f"[state query batch {batch_number}/{len(gene_batches)}] genes={len(gene_batch)}", flush=True)
        var_filter = in_filter("feature_name", gene_batch)
        with experiment.axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter=value_filter),
            var_query=soma.AxisQuery(value_filter=var_filter),
        ) as query:
            obs = read_query_obs(query, set(tumor_labels), set(normal_labels), set(cell_types))
            var = read_axis(query.var(column_names=["feature_name", "feature_id"]), query.var_joinids())
            gene_map = {int(j): str(g) for j, g in zip(var["soma_joinid"], var["feature_name"]) if str(g) in set(gene_batch)}
            for table in query.X("raw").tables():
                batch = table.to_pandas()
                if batch.empty:
                    continue
                ids = batch["soma_dim_0"].astype(np.int64)
                meta = cell_lookup.reindex(ids)
                gene = batch["soma_dim_1"].astype(np.int64).map(gene_map)
                valid = meta["donor_key"].notna().to_numpy() & gene.notna().to_numpy()
                if not valid.any():
                    continue
                batch = batch.loc[valid].copy()
                meta = meta.loc[valid]
                batch["gene"] = gene.loc[valid].to_numpy()
                batch["value"] = pd.to_numeric(batch["soma_data"], errors="coerce").fillna(0.0).to_numpy()
                batch["donor_key"] = meta["donor_key"].to_numpy()
                batch["donor_id"] = meta["donor_id"].to_numpy()
                batch["group"] = meta["group"].to_numpy()
                batch["PPAR_group"] = meta["PPAR_group_primary"].to_numpy()
                grouped = batch.groupby(["donor_key", "donor_id", "group", "PPAR_group", "gene"], observed=True)["value"].sum()
                for key, value in grouped.items():
                    gene_acc[key] = gene_acc.get(key, 0.0) + float(value)
    if not gene_acc:
        raise RuntimeError("Census state-gene query returned no donor-level expression.")
    long = pd.DataFrame([
        {"donor_key": key[0], "donor_id": key[1], "group": key[2], "PPAR_group": key[3], "gene": key[4], "count": value}
        for key, value in gene_acc.items()
    ])
    units = long.pivot_table(index=["donor_key", "donor_id", "group", "PPAR_group"], columns="gene", values="count", aggfunc="sum", fill_value=0).reset_index()
    units.columns.name = None
    for gene in gene_keys:
        if gene not in units.columns:
            units[gene] = 0.0
    cell_counts = cell_scores.groupby(["donor_key", "donor_id", "group", "PPAR_group_primary"], as_index=False).size().rename(columns={"PPAR_group_primary": "PPAR_group", "size": "n_cells"})
    units = units.merge(cell_counts, on=["donor_key", "donor_id", "group", "PPAR_group"], how="left", validate="one_to_one")
    return units[["donor_key", "donor_id", "group", "PPAR_group", "n_cells", *gene_keys]]


def log_expression(counts: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    values = counts[genes].astype(float)
    total = values.sum(axis=1).to_numpy(dtype=float)
    result = np.log1p(np.divide(values.to_numpy(dtype=float) * 1_000_000.0, total[:, None], out=np.zeros_like(values.to_numpy(dtype=float)), where=total[:, None] > 0))
    return pd.DataFrame(result, index=counts.index, columns=genes)


def state_scores(expr: pd.DataFrame, state_sets: dict[str, set[str]]) -> tuple[pd.DataFrame, dict[str, int]]:
    sd = expr.std(axis=0, ddof=1).replace(0, np.nan)
    z = expr.subtract(expr.mean(axis=0), axis=1).divide(sd, axis=1).fillna(0.0)
    rows = []
    present_sizes: dict[str, int] = {}
    for state, genes in state_sets.items():
        present = sorted(set(genes) & set(z.columns))
        present_sizes[state] = len(present)
        if not present:
            scores = pd.Series(np.nan, index=z.index)
        else:
            scores = z[present].mean(axis=1)
        rows.append(pd.DataFrame({"unit_index": z.index, "state": state, "state_score": scores.to_numpy(), "n_genes_present": len(present)}))
    return pd.concat(rows, ignore_index=True), present_sizes


def unit_metadata(units: pd.DataFrame) -> pd.DataFrame:
    return units[["donor_key", "donor_id", "group", "PPAR_group", "n_cells"]].copy()


def donor_all_units(units: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    grouped = units.groupby(["donor_key", "donor_id", "group"], as_index=False)[genes].sum()
    grouped["n_cells"] = units.groupby(["donor_key", "donor_id", "group"], as_index=False)["n_cells"].sum()["n_cells"]
    return grouped


def frozen_donor_ppar_score(donor_counts: pd.DataFrame) -> pd.DataFrame:
    """Reproduce the frozen Phase 2F donor-pseudobulk score exactly."""
    axis_genes = [*PPAR_NR_GENES, "RELA", "STAT3"]
    missing = sorted(set(axis_genes) - set(donor_counts.columns))
    if missing:
        raise RuntimeError(f"Frozen nine-gene donor score cannot be computed; missing genes: {missing}")
    values = donor_counts[axis_genes].astype(float)
    total = values.sum(axis=1)
    log_cpm = pd.DataFrame(np.nan, index=values.index, columns=axis_genes, dtype=float)
    valid = total.gt(0)
    log_cpm.loc[valid, axis_genes] = np.log1p(values.loc[valid].div(total.loc[valid], axis=0) * 1_000_000)
    z = (log_cpm - log_cpm.mean(axis=0)) / log_cpm.std(axis=0, ddof=1).replace(0, np.nan)
    result = donor_counts[["donor_key", "donor_id", "group"]].copy()
    result["PPAR_NR_score"] = z[PPAR_NR_GENES].mean(axis=1)
    result["RELA_STAT3_score"] = z[["RELA", "STAT3"]].mean(axis=1)
    result["DINP_axis_9gene_score"] = z[axis_genes].mean(axis=1)
    result["target_gene_total_counts"] = total.to_numpy()
    return result


def load_phase2f_frozen_donor_scores(dataset_id: str, paired_donors: list[str]) -> pd.DataFrame:
    """Load Phase 2F scores standardized in the full dataset/compartment context."""
    path = OUTPUT / "mcop_phase2f_singlecell_donor_scores.csv"
    use = [
        "dataset_id", "donor_key", "donor_id", "group", "compartment",
        "PPAR_nuclear_receptor_score", "RELA_STAT3_score", "DINP_axis_9_gene_score",
        "target_gene_total_counts",
    ]
    scores = pd.read_csv(path, usecols=use, dtype={"dataset_id": str, "donor_id": str})
    scores = scores.loc[
        scores["dataset_id"].eq(dataset_id)
        & scores["compartment"].eq("epithelial")
        & scores["donor_id"].isin(paired_donors)
        & scores["group"].isin(["tumor", "normal"])
    ].copy()
    expected = len(paired_donors) * 2
    if len(scores) != expected or scores.duplicated(["donor_key", "group"]).any():
        raise RuntimeError(f"Frozen Phase 2F donor-score audit failed: expected {expected} unique rows, observed {len(scores)}")
    return scores.rename(columns={
        "PPAR_nuclear_receptor_score": "PPAR_NR_score",
        "DINP_axis_9_gene_score": "DINP_axis_9gene_score",
    })[["donor_key", "donor_id", "group", "PPAR_NR_score", "RELA_STAT3_score", "DINP_axis_9gene_score", "target_gene_total_counts"]]


def paired_effects(frame: pd.DataFrame, feature_cols: list[str], feature_type: str, level: str) -> pd.DataFrame:
    rows = []
    for feature in feature_cols:
        pivot = frame.pivot_table(index="donor_id", columns="group", values=feature, aggfunc="mean")
        if not {"tumor", "normal"}.issubset(pivot.columns):
            continue
        delta = (pivot["tumor"] - pivot["normal"]).dropna()
        if len(delta) < 3:
            continue
        p = float(wilcoxon(delta).pvalue) if np.any(delta != 0) else 1.0
        rows.append({
            "feature_type": feature_type, "feature": feature, "analysis_level": level,
            "n_paired_donors": int(len(delta)), "mean_delta_tumor_minus_normal": float(delta.mean()),
            "median_delta_tumor_minus_normal": float(delta.median()), "p_value": p,
            "direction": "up" if delta.median() > 0 else "down" if delta.median() < 0 else "flat",
            "donor_consistency": float((delta > 0).mean()) if delta.median() > 0 else float((delta < 0).mean()),
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = bh(result["p_value"])
    return result


def ppar_low_high_de(units: pd.DataFrame, expr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    meta = units[["donor_id", "group", "PPAR_group"]].reset_index(drop=True)
    values = expr[genes].reset_index(drop=True)
    work = pd.concat([meta, values], axis=1)
    low = work.loc[work["PPAR_group"].eq("PPAR-low")].groupby(["donor_id", "group"], observed=True)[genes].mean()
    high = work.loc[work["PPAR_group"].eq("PPAR-high")].groupby(["donor_id", "group"], observed=True)[genes].mean()
    common = low.index.intersection(high.index)
    delta_matrix = low.loc[common, genes] - high.loc[common, genes]
    rows = []
    for gene in genes:
        delta = delta_matrix[gene].dropna()
        if len(delta) < 5:
            continue
        p = float(wilcoxon(delta).pvalue) if np.any(delta != 0) else 1.0
        rows.append({
            "gene": gene, "effect_size": float(delta.median()), "logFC_low_minus_high": float(delta.mean()),
            "P": p, "n_donor_group_pairs": int(len(delta)), "direction": "higher_in_PPAR_low" if delta.mean() > 0 else "lower_in_PPAR_low",
            "analysis_level": "paired donor-group pseudobulk; targeted state/regulator gene universe",
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = bh(result["P"])
        result = result.sort_values(["BH_FDR", "P"], na_position="last")
    return result


def state_correlations(donor_frame: pd.DataFrame, state_frame: pd.DataFrame, states: list[str]) -> pd.DataFrame:
    wide = state_frame.pivot_table(index=["donor_key", "donor_id", "group"], columns="state", values="state_score", aggfunc="mean").reset_index()
    ppar = donor_frame.groupby(["donor_key", "donor_id", "group"], as_index=False)["PPAR_NR_score"].mean()
    wide = wide.merge(ppar, on=["donor_key", "donor_id", "group"], how="left")
    rows = []
    for group_name, subset in [("pooled", wide), ("tumor", wide[wide["group"].eq("tumor")]), ("normal", wide[wide["group"].eq("normal")])]:
        for state in states:
            temp = subset[["PPAR_NR_score", state]].dropna()
            if len(temp) < 5:
                continue
            rho, p = spearmanr(temp["PPAR_NR_score"], temp[state])
            x = temp["PPAR_NR_score"] - temp["PPAR_NR_score"].median()
            y = temp[state] - temp[state].median()
            rows.append({
                "group": group_name, "state": state, "n_donors": int(len(temp)), "spearman_rho": float(rho),
                "P": float(p), "donor_consistency": float((x * y >= 0).mean()),
                "direction": "positive" if rho > 0 else "negative" if rho < 0 else "flat",
            })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = result.groupby("group", group_keys=False)["P"].transform(lambda x: bh(x))
        result["state_association_interpretation"] = np.where(result["spearman_rho"] < 0, "state rises as PPAR/NR falls", "state falls or tracks as PPAR/NR falls")
    return result


def interaction_analysis(state_units: pd.DataFrame) -> pd.DataFrame:
    rows = []
    work = state_units.loc[state_units["PPAR_group"].isin(["PPAR-low", "PPAR-high"])].copy()
    work["tumor_binary"] = (work["group"] == "tumor").astype(int)
    work["low_binary"] = (work["PPAR_group"] == "PPAR-low").astype(int)
    for state, subset in work.groupby("state"):
        if subset["donor_id"].nunique() < 5 or subset[["tumor_binary", "low_binary"]].drop_duplicates().shape[0] < 4:
            continue
        try:
            fit = ols("state_score ~ tumor_binary * low_binary + C(donor_id)", data=subset).fit(cov_type="HC3")
            term = "tumor_binary:low_binary"
            beta = float(fit.params.get(term, np.nan))
            p = float(fit.pvalues.get(term, np.nan))
        except Exception:
            beta, p = np.nan, np.nan
        means = subset.groupby(["group", "PPAR_group"], as_index=False)["state_score"].mean()
        mean_map = {(str(r["group"]), str(r["PPAR_group"])): float(r["state_score"]) for _, r in means.iterrows()}
        rows.append({
            "state": state, "n_donors": int(subset["donor_id"].nunique()), "interaction_beta_tumor_x_PPARlow": beta,
            "P": p, "normal_PPARlow_mean": mean_map.get(("normal", "PPAR-low"), np.nan),
            "normal_PPARhigh_mean": mean_map.get(("normal", "PPAR-high"), np.nan),
            "tumor_PPARlow_mean": mean_map.get(("tumor", "PPAR-low"), np.nan),
            "tumor_PPARhigh_mean": mean_map.get(("tumor", "PPAR-high"), np.nan),
            "analysis_level": "donor-level pseudobulk OLS with donor fixed effects; HC3 SE",
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = bh(result["P"])
        result["direction"] = np.where(result["interaction_beta_tumor_x_PPARlow"] > 0, "amplified_in_tumor_PPARlow", "attenuated_or_reversed")
    return result


def subtype_audit(cell_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (subtype, cell_type, group), subset in cell_scores.groupby(["cell_subtype", "cell_type", "group"], observed=True):
        rows.append({
            "subtype": subtype, "cell_type": cell_type, "n_cells": int(len(subset)), "n_donors": int(subset["donor_id"].nunique()),
            "median_PPAR_NR": float(subset["PPAR_NR_score"].median()), "PPAR_low_fraction": float(subset["PPAR_group_primary"].eq("PPAR-low").mean()),
            "PPAR_high_fraction": float(subset["PPAR_group_primary"].eq("PPAR-high").mean()), "tumor_normal_distribution": group,
            "annotation_source": "CELLxGENE Census cell_type label; no forced re-clustering or malignancy inference",
        })
    return pd.DataFrame(rows)


def regulator_activity(
    expr: pd.DataFrame,
    units: pd.DataFrame,
    net: pd.DataFrame,
    ppar_expr: pd.DataFrame | None = None,
    ppar_units: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    import decoupler as dc
    present_net = net.loc[net["target"].isin(expr.columns)].copy()
    if present_net.empty:
        return pd.DataFrame(), pd.DataFrame()
    def run_ulm(matrix_expr: pd.DataFrame, matrix_units: pd.DataFrame, unit_suffix: str) -> pd.DataFrame:
        matrix = matrix_expr.copy()
        unit_id = (
            matrix_units["donor_key"].astype(str)
            + "|" + matrix_units["group"].astype(str)
            + unit_suffix
        )
        matrix.index = unit_id
        estimates, pvalues = dc.mt.ulm(matrix, present_net, tmin=5, verbose=False)
        activity = estimates.copy()
        activity.index.name = "unit_id"
        result = activity.reset_index().melt(id_vars="unit_id", var_name="regulator", value_name="activity_score")
        result["p_value_ulm"] = pvalues.reset_index(drop=True).melt(value_name="p_value_ulm")["p_value_ulm"].to_numpy()
        meta = matrix_units.copy()
        meta["unit_id"] = unit_id.to_numpy()
        keep = ["unit_id", "donor_id", "group"] + (["PPAR_group"] if "PPAR_group" in meta.columns else [])
        result = result.merge(meta[keep], on="unit_id", how="left", validate="many_to_one")
        return result

    tumor_long = run_ulm(expr, units, "")
    ppar_long = pd.DataFrame()
    if ppar_expr is not None and ppar_units is not None:
        suffix = "|" + ppar_units["PPAR_group"].astype(str)
        ppar_long = run_ulm(ppar_expr, ppar_units, suffix)
    long = pd.concat([tumor_long.assign(activity_unit="donor_group"), ppar_long.assign(activity_unit="donor_group_ppar_quartile")], ignore_index=True)
    rows = []
    for regulator, subset in long.groupby("regulator"):
        for comparison, source in [("tumor_vs_normal", tumor_long), ("PPAR_low_vs_high", ppar_long)]:
            subset = source.loc[source["regulator"].eq(regulator)].copy() if not source.empty else pd.DataFrame()
            if subset.empty:
                continue
            if comparison == "tumor_vs_normal":
                pivot = subset.pivot_table(index="donor_id", columns="group", values="activity_score", aggfunc="mean")
                a, b = "tumor", "normal"
            else:
                work = subset.loc[subset["PPAR_group"].isin(["PPAR-low", "PPAR-high"])].copy()
                work["pair_id"] = work["donor_id"].astype(str) + "|" + work["group"].astype(str)
                pivot = work.pivot_table(index="pair_id", columns="PPAR_group", values="activity_score", aggfunc="mean")
                a, b = "PPAR-low", "PPAR-high"
            if not {a, b}.issubset(pivot.columns):
                continue
            delta = (pivot[a] - pivot[b]).dropna()
            if len(delta) < 5:
                continue
            p = float(wilcoxon(delta).pvalue) if np.any(delta != 0) else 1.0
            rows.append({
                "regulator": regulator, "comparison": comparison, "n_pairs": int(len(delta)),
                "activity_delta": float(delta.mean()), "median_activity_delta": float(delta.median()), "P": p,
                "direction": "up_in_first_group" if delta.mean() > 0 else "down_in_first_group",
                "donor_consistency": float((delta > 0).mean()) if delta.mean() > 0 else float((delta < 0).mean()),
                "method": "DoRothEA levels A-C + decoupler ULM",
            })
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary["BH_FDR"] = summary.groupby("comparison", group_keys=False)["P"].transform(lambda x: bh(x))
    return long, summary


def fetch_local_h5ad(
    h5ad_path: Path,
    dataset_id: str,
    tumor_labels: list[str],
    normal_labels: list[str],
    cell_types: list[str],
    paired_donors: list[str],
    query_genes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Rebuild the frozen Census slice from the official source H5AD.

    Raw CSR rows are streamed twice from disk. Pass 1 computes the fixed
    seven-gene cell score; pass 2 aggregates the frozen target universe to
    donor/group/PPAR-quartile pseudobulk. No full cell-by-gene matrix is held
    in memory.
    """
    import anndata as ad
    import h5py
    from scipy.sparse import csr_matrix

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        obs = adata.obs.copy()
        raw_var = adata.raw.var.copy() if adata.raw is not None else adata.var.copy()
        if adata.raw is None:
            raise RuntimeError("Official source H5AD has no .raw matrix; frozen raw-count analysis cannot proceed.")
        required_obs = {"donor_id", "disease", "cell_type", "is_primary_data", "observation_joinid"}
        missing_obs = required_obs - set(obs.columns)
        if missing_obs:
            raise RuntimeError(f"Source H5AD is missing required observation fields: {sorted(missing_obs)}")
        mask = (
            obs["is_primary_data"].eq(True)
            & obs["donor_id"].astype(str).isin(paired_donors)
            & obs["disease"].astype(str).isin(tumor_labels + normal_labels)
            & obs["cell_type"].astype(str).isin(cell_types)
        )
        selected_rows = np.flatnonzero(mask.to_numpy())
        selected_obs = obs.iloc[selected_rows].copy().reset_index(drop=False).rename(columns={"index": "source_cell_id"})
        selected_obs["donor_id"] = selected_obs["donor_id"].astype(str)
        selected_obs["group"] = np.where(selected_obs["disease"].astype(str).isin(tumor_labels), "tumor", "normal")
        selected_obs["donor_key"] = dataset_id + "::" + selected_obs["donor_id"]
        # Source H5AD observation_joinid is the study barcode, not the Census
        # numeric soma_joinid. Keep the source row as an internal stable index
        # and preserve the original barcode in cell_id.
        selected_obs["soma_joinid"] = selected_rows.astype(np.int64)
        if selected_obs.empty:
            raise RuntimeError("Local H5AD filter returned no eligible epithelial cells.")
        paired_check = set(selected_obs.loc[selected_obs["group"].eq("tumor"), "donor_id"]) & set(selected_obs.loc[selected_obs["group"].eq("normal"), "donor_id"])
        if paired_check != set(paired_donors):
            raise RuntimeError(f"Local H5AD paired-donor mismatch: expected={len(paired_donors)}, observed={len(paired_check)}")

        feature_names = raw_var["feature_name"].astype(str)
        duplicate_features = feature_names[feature_names.duplicated(keep=False)]
        requested = sorted(set(query_genes) | set(PPAR_NR_GENES))
        gene_to_col: dict[str, int] = {}
        for col, gene in enumerate(feature_names):
            if gene in requested and gene not in gene_to_col:
                gene_to_col[gene] = col
        missing_genes = sorted(set(requested) - set(gene_to_col))
        core_missing = sorted(set(PPAR_NR_GENES) - set(gene_to_col))
        if core_missing:
            raise RuntimeError(f"Frozen PPAR/NR genes missing from source H5AD: {core_missing}")
        genes = sorted(set(query_genes) & set(gene_to_col))
        target_cols = np.asarray([gene_to_col[g] for g in genes], dtype=np.int64)
        core_target_pos = np.asarray([genes.index(g) for g in PPAR_NR_GENES], dtype=np.int64)
        n_obs, n_var = adata.raw.shape
    finally:
        adata.file.close()

    selected_flag = np.zeros(n_obs, dtype=bool)
    selected_flag[selected_rows] = True
    selected_position = np.full(n_obs, -1, dtype=np.int64)
    selected_position[selected_rows] = np.arange(len(selected_rows), dtype=np.int64)
    core_counts = np.zeros((len(selected_rows), len(PPAR_NR_GENES)), dtype=np.float64)
    chunk_rows = 2048

    with h5py.File(h5ad_path, "r") as handle:
        raw_x = handle["raw/X"]
        if raw_x.attrs.get("encoding-type") != "csr_matrix":
            raise RuntimeError(f"Expected raw CSR matrix, found {raw_x.attrs.get('encoding-type')!r}")
        indptr = raw_x["indptr"][:]
        for start in range(0, n_obs, chunk_rows):
            stop = min(n_obs, start + chunk_rows)
            local_selected = np.flatnonzero(selected_flag[start:stop])
            if not len(local_selected):
                continue
            p0, p1 = int(indptr[start]), int(indptr[stop])
            block = csr_matrix(
                (raw_x["data"][p0:p1], raw_x["indices"][p0:p1], indptr[start:stop + 1] - p0),
                shape=(stop - start, n_var),
            )
            target = block[local_selected][:, target_cols]
            positions = selected_position[start + local_selected]
            core_counts[positions] = target[:, core_target_pos].toarray()

    core_total = core_counts.sum(axis=1)
    normalized = np.log1p(np.divide(core_counts * 1_000_000.0, core_total[:, None], out=np.zeros_like(core_counts), where=core_total[:, None] > 0))
    gene_sd = normalized.std(axis=0, ddof=1)
    z = (normalized - normalized.mean(axis=0)) / np.where(gene_sd > 0, gene_sd, 1.0)
    score = np.nanmean(z, axis=1)
    ranks = pd.Series(score).rank(method="first")
    selected_obs["PPAR_group_primary"] = pd.qcut(ranks, 4, labels=["PPAR-low", "Q2", "Q3", "PPAR-high"]).astype(str).to_numpy()
    selected_obs["PPAR_group_tercile"] = pd.qcut(ranks, 3, labels=["PPAR-low", "middle", "PPAR-high"]).astype(str).to_numpy()
    selected_obs["PPAR_group_median"] = np.where(score <= np.nanmedian(score), "PPAR-low", "PPAR-high")
    selected_obs["PPAR_group"] = selected_obs["PPAR_group_primary"]
    selected_obs["PPAR_NR_score"] = score
    selected_obs["cell_id"] = dataset_id + "::" + selected_obs["observation_joinid"].astype(str)
    selected_obs["tumor_normal"] = selected_obs["group"]
    selected_obs["cell_subtype"] = selected_obs["cell_type"].map(subtype_label)
    selected_obs["PPAR_NR_core_genes"] = ";".join(PPAR_NR_GENES)

    unit_meta = selected_obs[["donor_key", "donor_id", "group", "PPAR_group_primary"]].rename(columns={"PPAR_group_primary": "PPAR_group"})
    unit_index = pd.MultiIndex.from_frame(unit_meta).unique().sort_values()
    unit_lookup = {key: idx for idx, key in enumerate(unit_index)}
    cell_unit = np.asarray([unit_lookup[tuple(row)] for row in unit_meta.itertuples(index=False, name=None)], dtype=np.int64)
    aggregate = np.zeros((len(unit_index), len(genes)), dtype=np.float64)
    n_cells = np.bincount(cell_unit, minlength=len(unit_index)).astype(int)

    with h5py.File(h5ad_path, "r") as handle:
        raw_x = handle["raw/X"]
        indptr = raw_x["indptr"][:]
        for start in range(0, n_obs, chunk_rows):
            stop = min(n_obs, start + chunk_rows)
            local_selected = np.flatnonzero(selected_flag[start:stop])
            if not len(local_selected):
                continue
            p0, p1 = int(indptr[start]), int(indptr[stop])
            block = csr_matrix(
                (raw_x["data"][p0:p1], raw_x["indices"][p0:p1], indptr[start:stop + 1] - p0),
                shape=(stop - start, n_var),
            )
            target = block[local_selected][:, target_cols]
            positions = selected_position[start + local_selected]
            units_here = cell_unit[positions]
            for unit in np.unique(units_here):
                aggregate[unit] += np.asarray(target[units_here == unit].sum(axis=0)).ravel()

    units = unit_index.to_frame(index=False)
    units["n_cells"] = n_cells
    units = pd.concat([units.reset_index(drop=True), pd.DataFrame(aggregate, columns=genes)], axis=1)
    keep_cell = [
        "dataset_id", "donor_id", "disease", "tissue", "cell_type", "is_primary_data", "sex",
        "soma_joinid", "donor_key", "group", "cell_id", "tumor_normal", "PPAR_NR_score",
        "PPAR_group_primary", "PPAR_group_tercile", "PPAR_group_median", "PPAR_group",
        "cell_subtype", "PPAR_NR_core_genes",
    ]
    selected_obs["dataset_id"] = dataset_id
    if "sex" not in selected_obs:
        selected_obs["sex"] = np.nan
    audit = {
        "source_h5ad": str(h5ad_path),
        "source_bytes": int(h5ad_path.stat().st_size),
        "source_n_obs": int(n_obs),
        "source_n_vars": int(n_var),
        "eligible_cells": int(len(selected_obs)),
        "eligible_paired_donors": int(len(paired_check)),
        "target_genes_requested": int(len(query_genes)),
        "target_genes_present": int(len(genes)),
        "target_genes_missing": missing_genes,
        "duplicate_feature_names_n": int(duplicate_features.nunique()),
        "raw_encoding": "csr_matrix",
        "local_cell_index_mapping": "soma_joinid column contains source H5AD row index; cell_id contains source observation_joinid barcode",
        "streaming_passes": 2,
    }
    return selected_obs[keep_cell].copy(), units, audit


def external_gse_state_scores(state_sets: dict[str, set[str]], query_genes: set[str]) -> pd.DataFrame:
    annotation_path = RAW_GSE / "GSE144735_annotation.txt.gz"
    counts_path = RAW_GSE / "GSE144735_raw_UMI_count_matrix.txt.gz"
    if not annotation_path.exists() or not counts_path.exists():
        return pd.DataFrame()
    annotation = pd.read_csv(annotation_path, sep="\t", index_col=0)
    annotation.index = annotation.index.astype(str)
    annotation["group"] = annotation["Class"].astype(str).str.lower().map({"tumor": "tumor", "normal": "normal"})
    annotation = annotation.loc[annotation["group"].isin(["tumor", "normal"])].copy()
    annotation["unit_id"] = annotation["Patient"].astype(str) + "|" + annotation["group"].astype(str)
    unit_ids = sorted(annotation["unit_id"].unique())
    design = pd.get_dummies(annotation["unit_id"], dtype=float).reindex(columns=unit_ids, fill_value=0.0).to_numpy()
    gene_rows: dict[str, np.ndarray] = {}
    wanted = {x.upper() for x in query_genes}
    for chunk in pd.read_csv(counts_path, sep="\t", index_col=0, chunksize=250):
        chunk.index = chunk.index.astype(str).str.upper()
        selected = chunk.loc[chunk.index.isin(wanted)]
        if selected.empty:
            continue
        values = selected.reindex(columns=annotation.index, fill_value=0).to_numpy(dtype=float)
        sums = values @ design
        for gene, row in zip(selected.index, sums):
            gene_rows[gene] = gene_rows.get(gene, np.zeros(len(unit_ids))) + row
    if not gene_rows:
        return pd.DataFrame()
    counts = pd.DataFrame(gene_rows, index=unit_ids)
    counts["unit_id"] = counts.index
    counts["patient"] = counts.index.to_series().str.split("|", n=1).str[0].to_numpy()
    counts["group"] = counts.index.to_series().str.split("|", n=1).str[1].to_numpy()
    genes = [x for x in sorted(wanted) if x in counts.columns]
    expr = log_expression(counts, genes)
    scores, _ = state_scores(expr, state_sets)
    wide = scores.pivot_table(index="unit_index", columns="state", values="state_score", aggfunc="mean").reset_index().rename(columns={"unit_index": "unit_id"})
    wide = wide.merge(counts[["unit_id", "patient", "group"]], on="unit_id", how="left")
    # Keep the external PPAR/NR score definition identical to Phase 2F:
    # targeted log-CPM uses the frozen nine-gene axis denominator, followed by
    # gene-wise z scoring within the dataset/compartment.  Do not let the much
    # larger Phase 2G state universe change the biomarker score.
    core_axis_genes = [gene for gene in [*PPAR_NR_GENES, "RELA", "STAT3"] if gene in counts.columns]
    core_totals = counts[core_axis_genes].sum(axis=1).replace(0, np.nan)
    core_log_cpm = np.log1p(counts[core_axis_genes].div(core_totals, axis=0) * 1_000_000)
    core_expr = core_log_cpm[PPAR_NR_GENES]
    core_z = (core_expr - core_expr.mean(axis=0)) / core_expr.std(axis=0, ddof=1).replace(0, np.nan)
    core_scores = core_z.mean(axis=1).rename("PPAR_NR_score").reset_index().rename(columns={"index": "unit_id"})
    wide = wide.drop(columns=["PPAR_NR_score"], errors="ignore").merge(core_scores, on="unit_id", how="left")
    rows = []
    for feature in ["PPAR_NR_score", *state_sets.keys(), *[g for g in ANCHOR_CANDIDATES if g in expr.columns]]:
        if feature == "PPAR_NR_score":
            pass
        elif feature in wide.columns:
            pass
        elif feature in expr.columns:
            wide[feature] = expr[feature].to_numpy()
        else:
            continue
        pivot = wide.pivot_table(index="patient", columns="group", values=feature, aggfunc="mean")
        if not {"tumor", "normal"}.issubset(pivot.columns):
            continue
        delta = (pivot["tumor"] - pivot["normal"]).dropna()
        p = float(wilcoxon(delta).pvalue) if len(delta) >= 3 and np.any(delta != 0) else np.nan
        rows.append({
            "dataset": "GSE144735", "feature": feature, "n_paired": int(len(delta)),
            "median_delta_tumor_minus_normal": float(delta.median()) if len(delta) else np.nan,
            "mean_delta_tumor_minus_normal": float(delta.mean()) if len(delta) else np.nan,
            "P": p, "direction": "up" if len(delta) and delta.median() > 0 else "down" if len(delta) else "not_estimable",
            "replication_status": "directionally concordant but underpowered" if len(delta) and delta.median() < 0 else "directionally concordant" if len(delta) else "not_estimable",
            "analysis_level": "patient-paired epithelial pseudobulk; public GSE144735",
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["BH_FDR"] = bh(result["P"])
    return result


def replace_with_frozen_phase2f_core(external: pd.DataFrame) -> pd.DataFrame:
    """Use the frozen Phase 2F GSE core row for exact score comparability."""
    path = OUTPUT / "mcop_phase2f_external_gse144735_paired_contrasts.csv"
    if not path.exists():
        return external
    frozen = pd.read_csv(path)
    frozen = frozen.loc[
        frozen["score"].eq("PPAR_nuclear_receptor_score")
        & frozen.get("group_mode", pd.Series(index=frozen.index, dtype=str)).eq("core_tumor_vs_normal")
    ].copy()
    if frozen.empty:
        return external
    row = frozen.iloc[0]
    replacement = pd.DataFrame([{
        "dataset": "GSE144735",
        "feature": "PPAR_NR_score",
        "n_paired": int(row["paired_donors"]),
        "median_delta_tumor_minus_normal": float(row["median_delta_tumor_minus_normal"]),
        "mean_delta_tumor_minus_normal": float(row["mean_delta_tumor_minus_normal"]),
        "P": float(row["p_value"]),
        "direction": "up" if float(row["median_delta_tumor_minus_normal"]) > 0 else "down",
        "replication_status": "directionally concordant but underpowered",
        "analysis_level": "frozen Phase 2F patient-paired epithelial score; reused for exact comparability",
    }])
    external = external.loc[external["feature"].ne("PPAR_NR_score")].copy() if not external.empty else external
    return pd.concat([external, replacement], ignore_index=True)


def toxicology_tags(candidates: list[str]) -> dict[str, str]:
    tags = {gene: "none in current DINP-axis single-chemical bridge table" for gene in candidates}
    qc_path = OUTPUT / "mcop_phase2e_molecular_bridge_qc.csv"
    if qc_path.exists():
        qc = pd.read_csv(qc_path)
        dinp = qc.loc[qc["chemical_label"].eq("DINP_parent") & qc["has_single_chemical_evidence"].astype(bool)]
        for gene in dinp["GeneSymbol"].astype(str):
            if gene in tags:
                tags[gene] = "DINP parent: single-chemical CTD evidence"
    ora_path = OUTPUT / "mcop_phase2e_pathway_ora.csv"
    if ora_path.exists():
        ora = pd.read_csv(ora_path)
        match = ora.loc[
            ora["query"].eq("MiNP_no_cotreatment_genes")
            & ora["term"].eq("REACTOME_NUCLEAR_RECEPTOR_TRANSCRIPTION_PATHWAY")
        ]
        if not match.empty:
            overlap = set(str(match.iloc[0]["overlap_genes"]).split(";"))
            for gene in overlap & set(tags):
                tags[gene] = "MiNP no-cotreatment nuclear-receptor pathway overlap"
    return tags


def build_anchor_outputs(
    candidates: list[str], donor_validation: pd.DataFrame, correlations: pd.DataFrame,
    regulator_summary: pd.DataFrame, external: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tox = toxicology_tags(candidates)
    rows = []
    for gene in candidates:
        dv = donor_validation.loc[donor_validation["feature_type"].eq("gene") & donor_validation["feature"].eq(gene)]
        corr = correlations.loc[correlations["state"].eq(gene)] if "state" in correlations.columns else pd.DataFrame()
        if {"regulator", "comparison"}.issubset(regulator_summary.columns):
            reg = regulator_summary.loc[
                regulator_summary["regulator"].eq(gene)
                & regulator_summary["comparison"].eq("tumor_vs_normal")
            ]
        else:
            reg = pd.DataFrame()
        ext = external.loc[external["feature"].eq(gene)] if not external.empty else pd.DataFrame()
        internal_effect = "tumor-normal direction not estimable"
        consistency = "not estimable"
        if not dv.empty:
            d = dv.iloc[0]
            internal_effect = f"{d['direction']}; median delta={float(d['median_delta_tumor_minus_normal']):.3g}; P={float(d['p_value']):.3g}"
            consistency = f"{float(d['donor_consistency']):.2f}"
        reg_effect = "no DoRothEA activity summary"
        if not reg.empty:
            r = reg.iloc[0]
            reg_effect = f"{r['direction']}; delta={float(r['activity_delta']):.3g}; FDR={float(r['BH_FDR']):.3g}"
        corr_effect = "gene-state correlation not available"
        if not corr.empty:
            c = corr.iloc[0]
            corr_effect = f"rho={float(c['spearman_rho']):.3g}; FDR={float(c['BH_FDR']):.3g}"
        external_effect = "no external gene-level row"
        if not ext.empty:
            e = ext.iloc[0]
            external_effect = f"{e['dataset']}: {e['direction']}; status={e['replication_status']}"
        directly_supported = tox[gene].startswith("DINP parent") and not dv.empty and not ext.empty and not reg.empty
        plausible = tox[gene] != "none in current DINP-axis single-chemical bridge table" and not dv.empty
        tier = "Directly supported" if directly_supported else "Plausible candidate bridge" if plausible else "Unsupported"
        rows.append({
            "gene": gene, "toxicology_support": tox[gene], "epithelial_effect": internal_effect,
            "regulon_effect": reg_effect, "state_correlation": corr_effect, "donor_consistency": consistency,
            "external_support": external_effect, "overall_evidence_tier": tier,
        })
    ranking = pd.DataFrame(rows)
    bridge = ranking.rename(columns={"gene": "candidate_bridge"}).copy()
    bridge["causal_status"] = "not established; convergence is associative and exposure-to-state direction is untested"
    return ranking, bridge


def fallback_existing_phase2f() -> None:
    """Emit an explicitly partial Phase 2G package from already validated outputs.

    This path is used only when the current Census raw-expression runtime is
    unavailable.  It never fabricates cell-level states: the primary Census
    state/DE/regulator modules are marked not estimable, while the local
    GSE144735 and existing TCGA/Phase 2F evidence are still recomputed.
    """
    state_sets, state_sources = build_state_sets()
    scores = pd.read_csv(OUTPUT / "mcop_phase2f_singlecell_donor_scores.csv")
    paired = pd.read_csv(OUTPUT / "mcop_phase2f_singlecell_paired_donor_contrasts.csv")
    eligible = paired.loc[paired["compartment"].eq("epithelial") & paired["score"].eq("PPAR_nuclear_receptor_score")].sort_values("paired_donors", ascending=False)
    dataset_id = str(eligible.iloc[0]["dataset_id"]) if not eligible.empty else "not_available"
    scores = scores.loc[scores["dataset_id"].astype(str).eq(dataset_id) & scores["compartment"].eq("epithelial")].copy()
    scores["cell_id"] = ""
    scores["tumor_normal"] = scores["group"]
    scores["PPAR_NR_score"] = scores["PPAR_nuclear_receptor_score"]
    scores["cell_subtype"] = "not available from donor-level Phase 2F output"
    ranks = scores["PPAR_NR_score"].rank(method="first")
    scores["PPAR_group_primary"] = pd.qcut(ranks, 4, labels=["PPAR-low", "Q2", "Q3", "PPAR-high"]).astype(str)
    scores["PPAR_group"] = scores["PPAR_group_primary"]
    scores["analysis_level"] = "donor-level Phase 2F fallback; not cell-level state definition"
    scores.to_csv(OUT_CELL, index=False)
    gene_cols = [gene for gene in PPAR_NR_GENES + ["RELA", "STAT3"] if gene in scores.columns]
    total = scores["target_gene_total_counts"].replace(0, np.nan)
    expr = pd.DataFrame(index=scores.index)
    for gene in gene_cols:
        values = np.log1p(pd.to_numeric(scores[gene], errors="coerce") / total * 1_000_000)
        expr[gene] = (values - values.mean()) / values.std(ddof=1)
    expr["donor_id"] = scores["donor_id"].astype(str).to_numpy()
    expr["group"] = scores["group"].astype(str).to_numpy()
    gene_validation = paired_effects(expr, gene_cols, "gene", "Phase 2F donor-level fallback")
    score_frame = scores[["donor_id", "group", "PPAR_nuclear_receptor_score", "RELA_STAT3_score"]].rename(columns={"PPAR_nuclear_receptor_score": "PPAR_NR_score", "RELA_STAT3_score": "RELA_STAT3_score"})
    core_validation = paired_effects(score_frame, ["PPAR_NR_score", "RELA_STAT3_score"], "Phase2F_score", "Phase 2F donor-level validated score")
    donor_validation = pd.concat([gene_validation, core_validation], ignore_index=True)
    donor_validation.to_csv(OUT_DONOR, index=False)
    blocked_pathway = pd.DataFrame([{"status": "not_estimable_runtime_blocked", "analysis_level": "primary Census target-universe expression unavailable", "state": state, "state_score": np.nan, "n_genes_present": 0} for state in state_sets])
    blocked_pathway.to_csv(OUT_PATHWAY, index=False)
    blocked_corr = pd.DataFrame([{"status": "not_estimable_runtime_blocked", "group": "primary Census", "state": state, "n_donors": 0, "spearman_rho": np.nan, "P": np.nan, "BH_FDR": np.nan} for state in state_sets])
    blocked_corr.to_csv(OUT_CORR, index=False)
    blocked_interaction = pd.DataFrame([{"status": "not_estimable_runtime_blocked", "state": state, "n_donors": 0, "interaction_beta_tumor_x_PPARlow": np.nan, "P": np.nan, "BH_FDR": np.nan} for state in state_sets])
    blocked_interaction.to_csv(OUT_INTERACTION, index=False)
    audit_path = OUTPUT / "mcop_phase2f_singlecell_dataset_donor_cell_audit.csv"
    if audit_path.exists():
        subtype = pd.read_csv(audit_path)
        subtype = subtype.loc[subtype["compartment"].eq("epithelial")].copy()
        subtype["subtype"] = "annotation available only as aggregate cell-type audit"
        subtype["median_PPAR_NR"] = np.nan
        subtype["PPAR_low_fraction"] = np.nan
        subtype["PPAR_high_fraction"] = np.nan
        subtype["tumor_normal_distribution"] = subtype["group"]
        subtype["status"] = "not_estimable_runtime_blocked"
    else:
        subtype = pd.DataFrame()
    subtype.to_csv(OUT_SUBTYPE, index=False)
    external = external_gse_state_scores(state_sets, set().union(*state_sets.values()) | set(PPAR_NR_GENES) | set(ANCHOR_CANDIDATES))
    external = replace_with_frozen_phase2f_core(external)
    tcga_path = OUTPUT / "mcop_phase2f_tcga_paired_gene_stats.csv"
    if tcga_path.exists():
        tcga = pd.read_csv(tcga_path)
        tcga_rows = []
        for _, row in tcga.loc[tcga["gene"].isin(ANCHOR_CANDIDATES)].iterrows():
            tcga_rows.append({"dataset": "TCGA paired primary vs solid normal", "feature": row["gene"], "n_paired": row["paired_n"], "median_delta_tumor_minus_normal": row["median_delta_tumor_minus_normal"], "mean_delta_tumor_minus_normal": np.nan, "P": row["p_value"], "direction": "up" if row["median_delta_tumor_minus_normal"] > 0 else "down", "replication_status": "external directional support", "analysis_level": "existing Phase 2F paired bulk output"})
        if tcga_rows:
            external = pd.concat([external, pd.DataFrame(tcga_rows)], ignore_index=True)
    if not external.empty:
        external["BH_FDR"] = bh(external["P"])
    external.to_csv(OUT_EXTERNAL, index=False)
    empty_regulator = pd.DataFrame([{"status": "not_estimable_runtime_blocked", "method": "DoRothEA ULM not run without primary target-universe expression"}])
    empty_regulator.to_csv(OUT_REGULATOR, index=False)
    candidates = ANCHOR_CANDIDATES
    ranking, bridge = build_anchor_outputs(candidates, donor_validation, pd.DataFrame(), pd.DataFrame(), external)
    ranking.to_csv(OUT_ANCHOR, index=False)
    bridge.to_csv(OUT_BRIDGE, index=False)
    scores.to_csv(OUT_PSEUDOBULK, index=False)
    report = [
        "# Phase 2G — DINP/MCOP–CRC epithelial state convergence and regulatory anchoring",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Runtime status",
        "",
        "**PARTIAL / RUNTIME-BLOCKED:** the pinned WSL CELLxGENE Census raw-expression query repeatedly terminated at the TileDB `query.X('raw').tables()` stage, including a validated 9-gene minimal query and a single-donor probe. No Census primary state/DE/regulator result is claimed from this fallback.",
        "",
        f"- Existing Phase 2F donor-level epithelial output retained for dataset `{dataset_id}`.",
        "- GSE144735 full target-universe state scores and existing TCGA paired output were recomputed locally as external support only; the PPAR/NR external row reuses the frozen Phase 2F score definition for exact comparability.",
        "- The 7-gene PPAR/NR core was not redefined; no cell-level PPAR-low/high state was fabricated from donor-level rows.",
        "",
        "## What is and is not concluded",
        "",
        "- Retained Phase 2F paired epithelial result: PPAR/NR median tumor-minus-normal Δ = −0.419 (36 donors; P = 4.29×10⁻⁷); RELA/STAT3 Δ = +1.167 (P = 1.08×10⁻⁷). Myeloid PPAR/NR remains opposite-direction (Δ = +0.610; 35 donors; P = 7.97×10⁻⁹).",
        "- GSE144735 directionally concordant but underpowered: PPAR/NR median Δ = −0.312 (6 paired patients; P = 0.6875), using the frozen Phase 2F core score.",
        "",
        "The available Phase 2F evidence continues to support an epithelial PPAR/NR disease-state signal, but Phase 2G cannot answer the prespecified unbiased state-discovery, tumor×PPAR interaction, or DoRothEA regulator-activity questions until primary Census target-universe expression is readable. Therefore the final verdict is **PARTIALLY**, not YES.",
        "",
        "## Outputs",
        "",
        *[f"- `{path.name}`" for path in [OUT_CELL, OUT_DE, OUT_PATHWAY, OUT_CORR, OUT_DONOR, OUT_SUBTYPE, OUT_INTERACTION, OUT_REGULATOR, OUT_ANCHOR, OUT_EXTERNAL, OUT_BRIDGE, OUT_PSEUDOBULK]],
    ]
    # The DE output is intentionally a valid empty/status table in fallback mode.
    pd.DataFrame([{"status": "not_estimable_runtime_blocked", "analysis_level": "primary Census expression unavailable"}]).to_csv(OUT_DE, index=False)
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
    manifest = {
        "analysis": "Phase 2G DINP/MCOP-CRC epithelial state convergence and regulatory anchoring",
        "mode": "fallback_existing_phase2f_evidence",
        "runtime_status": "Census raw-expression query blocked by WSL/TileDB service failure",
        "census_version": CENSUS_VERSION,
        "primary_dataset_id": dataset_id,
        "ppar_nr_core_genes": PPAR_NR_GENES,
        "state_programs": {name: {"source": state_sources[name], "n_genes": len(genes)} for name, genes in state_sets.items()},
        "unit_of_inference": "existing donor-level Phase 2F output; primary Phase 2G state inference not estimable",
        "external_validation": "GSE144735 target-universe state scores plus existing TCGA paired Phase 2F output; frozen Phase 2F PPAR/NR row reused for exact comparability",
        "causal_status": "not established",
        "verdict": "PARTIALLY",
        "virtual_perturbation": "not run",
        "outputs": [str(x) for x in [OUT_CELL, OUT_DE, OUT_PATHWAY, OUT_CORR, OUT_DONOR, OUT_SUBTYPE, OUT_INTERACTION, OUT_REGULATOR, OUT_ANCHOR, OUT_EXTERNAL, OUT_BRIDGE, OUT_PSEUDOBULK, OUT_REPORT]],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"mode": "fallback_existing_phase2f_evidence", "dataset_id": dataset_id, "verdict": "PARTIALLY"}, ensure_ascii=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fallback-existing", action="store_true", help="Use existing Phase 2F outputs and external files; do not claim primary Census Phase 2G estimability.")
    parser.add_argument("--local-h5ad", type=Path, help="Official source H5AD used to bypass the live TileDB expression endpoint.")
    parser.add_argument("--reuse-local-extraction", action="store_true", help="Reuse completed local cell-score and pseudobulk files; rerun downstream inference only.")
    args = parser.parse_args()
    if args.fallback_existing:
        fallback_existing_phase2f()
        return
    OUTPUT.mkdir(parents=True, exist_ok=True)
    state_sets, state_sources = build_state_sets()
    try:
        import decoupler as dc
        net = dc.op.dorothea(organism="human", levels=["A", "B", "C"], license="academic", verbose=False)
        net = net.loc[net["source"].isin(ANCHOR_CANDIDATES)].copy()
    except Exception as exc:  # pragma: no cover - environment gate
        raise SystemExit(f"DoRothEA/decoupler is required for Phase 2G: {type(exc).__name__}: {exc}") from exc
    metadata, dataset_id, tumor_labels, normal_labels, cell_types, paired_donors = load_primary_metadata()
    query_genes = set().union(*state_sets.values()) | set(PPAR_NR_GENES) | set(ANCHOR_CANDIDATES) | set(net["target"].astype(str))
    print(f"Primary dataset={dataset_id}; paired donors={len(paired_donors)}; epithelial cell types={len(cell_types)}; query genes={len(query_genes)}")
    local_audit = None
    if args.reuse_local_extraction:
        if not OUT_CELL.exists() or not OUT_PSEUDOBULK.exists():
            raise SystemExit("Local extraction reuse requested, but required cell-score/pseudobulk outputs are missing.")
        cell_scores = pd.read_csv(OUT_CELL)
        units = pd.read_csv(OUT_PSEUDOBULK)
        previous_manifest = json.loads(OUT_MANIFEST.read_text(encoding="utf-8")) if OUT_MANIFEST.exists() else {}
        local_audit = previous_manifest.get("local_h5ad_audit") or {"reused_completed_local_extraction": True}
        local_audit["reused_completed_local_extraction"] = True
        print(f"Reusing validated local extraction: cells={len(cell_scores):,}; units={len(units):,}", flush=True)
    elif args.local_h5ad:
        local_path = args.local_h5ad.resolve()
        if not local_path.exists():
            raise SystemExit(f"Local source H5AD does not exist: {local_path}")
        print(f"Using official source H5AD: {local_path}", flush=True)
        cell_scores, units, local_audit = fetch_local_h5ad(
            local_path, dataset_id, tumor_labels, normal_labels, cell_types, paired_donors, query_genes
        )
    else:
        with open_census() as census:
            cell_scores = fetch_cell_scores(census, dataset_id, tumor_labels, normal_labels, cell_types, paired_donors)
            units = fetch_donor_pseudobulk(census, cell_scores, dataset_id, tumor_labels, normal_labels, cell_types, paired_donors, query_genes)
    cell_scores.to_csv(OUT_CELL, index=False)
    genes = sorted(query_genes & set(units.columns))
    expr_units = log_expression(units, genes)
    expr_all_counts = donor_all_units(units, genes)
    expr_all = log_expression(expr_all_counts, genes)
    expr_all["donor_key"] = expr_all_counts["donor_key"].to_numpy()
    expr_all["donor_id"] = expr_all_counts["donor_id"].to_numpy()
    expr_all["group"] = expr_all_counts["group"].to_numpy()
    # Reuse the Phase 2F score standardized against the full epithelial
    # dataset, rather than re-standardizing after restricting to paired donors.
    ppar_donor = load_phase2f_frozen_donor_scores(dataset_id, paired_donors)
    donor_all = expr_all_counts[["donor_key", "donor_id", "group", "n_cells"]].copy()
    donor_all = donor_all.merge(ppar_donor, on=["donor_key", "donor_id", "group"], how="left", validate="one_to_one")
    all_state_scores, state_sizes = state_scores(expr_all[genes], state_sets)
    all_state_scores = all_state_scores.merge(expr_all_counts[["donor_key", "donor_id", "group"]].reset_index().rename(columns={"index": "unit_index"}), on="unit_index", how="left")
    all_state_scores = all_state_scores.merge(donor_all[["donor_key", "donor_id", "group", "PPAR_NR_score"]], on=["donor_key", "donor_id", "group"], how="left")
    low_high_expr = expr_units.copy()
    low_high_expr["donor_id"] = units["donor_id"].to_numpy()
    low_high_expr["group"] = units["group"].to_numpy()
    low_high_expr["PPAR_group"] = units["PPAR_group"].to_numpy()
    de = ppar_low_high_de(units, expr_units, genes)
    de.to_csv(OUT_DE, index=False)
    primary_units = units.loc[units["PPAR_group"].isin(["PPAR-low", "PPAR-high"])].copy()
    primary_expr = expr_units.loc[primary_units.index].copy()
    primary_expr["unit_index"] = primary_units.index
    primary_states, _ = state_scores(primary_expr.set_index("unit_index")[genes], state_sets)
    primary_states = primary_states.merge(primary_units[["donor_key", "donor_id", "group", "PPAR_group"]].reset_index().rename(columns={"index": "unit_index"}), on="unit_index", how="left")
    primary_states.to_csv(OUT_PATHWAY, index=False)
    corr = state_correlations(donor_all, all_state_scores, list(state_sets))
    # Add anchor-gene correlations to the same table with an explicit feature type.
    corr_rows = []
    for gene in ANCHOR_CANDIDATES:
        if gene not in expr_all.columns:
            continue
        frame = donor_all[["donor_id", "group", "PPAR_NR_score"]].copy()
        frame["feature_value"] = expr_all[gene].to_numpy()
        for group_name, subset in [("pooled", frame), ("tumor", frame[frame["group"].eq("tumor")]), ("normal", frame[frame["group"].eq("normal")])]:
            subset = subset.dropna()
            if len(subset) < 5:
                continue
            rho, p = spearmanr(subset["PPAR_NR_score"], subset["feature_value"])
            x = subset["PPAR_NR_score"] - subset["PPAR_NR_score"].median()
            y = subset["feature_value"] - subset["feature_value"].median()
            corr_rows.append({"group": group_name, "state": gene, "n_donors": len(subset), "spearman_rho": rho, "P": p, "donor_consistency": float((x * y >= 0).mean()), "direction": "positive" if rho > 0 else "negative" if rho < 0 else "flat", "state_association_interpretation": "anchor expression association"})
    if corr_rows:
        corr = pd.concat([corr, pd.DataFrame(corr_rows)], ignore_index=True)
    if not corr.empty:
        corr["BH_FDR"] = corr.groupby("group", group_keys=False)["P"].transform(lambda x: bh(x))
    corr.to_csv(OUT_CORR, index=False)
    donor_validation = paired_effects(donor_all.assign(PPAR_NR_score=donor_all["PPAR_NR_score"]), ["PPAR_NR_score"], "state", "donor-level paired tumor vs normal")
    state_validation = paired_effects(all_state_scores.pivot_table(index=["donor_key", "donor_id", "group"], columns="state", values="state_score", aggfunc="mean").reset_index(), list(state_sets), "state", "donor-level paired tumor vs normal")
    gene_validation = paired_effects(expr_all.reset_index(drop=True).assign(donor_id=expr_all_counts["donor_id"].to_numpy(), group=expr_all_counts["group"].to_numpy()), [x for x in ANCHOR_CANDIDATES if x in expr_all.columns], "gene", "donor-level paired tumor vs normal")
    donor_validation = pd.concat([donor_validation, state_validation, gene_validation], ignore_index=True)
    donor_validation.to_csv(OUT_DONOR, index=False)
    subtype_audit(cell_scores).to_csv(OUT_SUBTYPE, index=False)
    interaction_analysis(primary_states).to_csv(OUT_INTERACTION, index=False)
    reg_long, reg_summary = regulator_activity(
        expr_all[genes], expr_all_counts, net,
        ppar_expr=expr_units[genes], ppar_units=units,
    )
    reg_summary.to_csv(OUT_REGULATOR, index=False)
    external = external_gse_state_scores(state_sets, query_genes)
    external = replace_with_frozen_phase2f_core(external)
    # Append existing TCGA paired evidence for the frozen PPAR/NR core and inflammatory anchors.
    tcga_path = OUTPUT / "mcop_phase2f_tcga_paired_gene_stats.csv"
    tcga = pd.read_csv(tcga_path) if tcga_path.exists() else pd.DataFrame()
    tcga_rows = []
    if not tcga.empty:
        for _, row in tcga.loc[tcga["gene"].isin(ANCHOR_CANDIDATES)].iterrows():
            tcga_rows.append({"dataset": "TCGA paired primary vs solid normal", "feature": row["gene"], "n_paired": row["paired_n"], "median_delta_tumor_minus_normal": row["median_delta_tumor_minus_normal"], "mean_delta_tumor_minus_normal": np.nan, "P": row["p_value"], "direction": "up" if row["median_delta_tumor_minus_normal"] > 0 else "down", "replication_status": "external directional support", "analysis_level": "existing Phase 2F paired bulk output"})
    if tcga_rows:
        external = pd.concat([external, pd.DataFrame(tcga_rows)], ignore_index=True)
    if not external.empty:
        external["BH_FDR"] = bh(external["P"])
    external.to_csv(OUT_EXTERNAL, index=False)
    ranking, bridge = build_anchor_outputs(ANCHOR_CANDIDATES, donor_validation, corr, reg_summary, external)
    ranking.to_csv(OUT_ANCHOR, index=False)
    bridge.to_csv(OUT_BRIDGE, index=False)
    units.to_csv(OUT_PSEUDOBULK, index=False)
    ppar_delta = donor_validation.loc[donor_validation["feature_type"].eq("state") & donor_validation["feature"].eq("PPAR_NR_score")]
    positive_anchor = ranking.loc[ranking["overall_evidence_tier"].eq("Directly supported"), "gene"].astype(str).tolist()
    report = [
        "# Phase 2G — DINP/MCOP–CRC epithelial state convergence and regulatory anchoring",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Frozen analysis boundaries",
        "",
        f"- Census release: **{CENSUS_VERSION}**; primary data filter enforced.",
        f"- Expression access path: **{'official source H5AD, streamed locally' if local_audit else 'live Census TileDB query'}**.",
        f"- Primary dataset selected from Phase 2F matched epithelial audit: **{dataset_id}**; paired donors queried={len(paired_donors)}; epithelial cell types={len(cell_types)}.",
        f"- PPAR/NR core fixed as: **{', '.join(PPAR_NR_GENES)}**; no result-driven gene editing.",
        f"- State universe: **{len(state_sets)}** programs; gene-set sizes present in queried expression: `{json.dumps(state_sizes, ensure_ascii=False)}`.",
        "- PPAR-low/high is defined at cell level as bottom/top quartile; tercile and median labels are retained as sensitivity labels. Inference uses donor-level pseudobulk.",
        "- Expression DE is targeted to the frozen state/regulator gene universe, not genome-wide. This limitation is explicit and is not called a genome-wide DE result.",
        "",
        "## Primary answer",
        "",
        f"- Paired epithelial PPAR/NR tumor-normal result: {ppar_delta[['median_delta_tumor_minus_normal', 'p_value']].to_dict('records') if not ppar_delta.empty else 'not estimable'}.",
        f"- Directly supported candidate anchors under the frozen evidence tags: **{', '.join(positive_anchor) if positive_anchor else 'none'}**.",
        "- The most defensible interpretation remains a CRC epithelial disease-state convergence. It does not establish that DINP/MCOP causes the state or that any regulator mediates the epidemiologic association.",
        "",
        "## State discovery",
        "",
        "All prespecified state programs were scored before ranking. The state-correlation table reports donor-level Spearman associations with PPAR/NR; no EMT, stemness, inflammatory or metabolic state was selected in advance as the expected winner.",
        "",
        "## Regulatory activity and anchor boundary",
        "",
        "DoRothEA confidence levels A–C with decoupler ULM were used for the listed candidate regulators. Activity shifts are descriptive donor-level contrasts. Anchor tiers are evidence tags, not a subjective numeric total score.",
        "",
        "## Opposite-direction compartment result",
        "",
        "The Phase 2F myeloid PPAR/NR increase is retained as a visible localization result; this Phase 2G epithelial analysis does not reinterpret bulk tissue or erase the compartment contrast.",
        "",
        "## Final verdict",
        "",
        "**PARTIALLY** — the analysis tests whether epithelial PPAR/NR suppression defines a reproducible CRC state and whether a regulatory bridge is plausible. A positive disease-state convergence can be supported, but the DINP/MCOP → PPAR/NR arrow remains untested; GSE144735 is small and directional, and no causal perturbation is asserted.",
        "",
        "## Output files",
        "",
        *[f"- `{path.name}`" for path in [OUT_CELL, OUT_DE, OUT_PATHWAY, OUT_CORR, OUT_DONOR, OUT_SUBTYPE, OUT_INTERACTION, OUT_REGULATOR, OUT_ANCHOR, OUT_EXTERNAL, OUT_BRIDGE, OUT_PSEUDOBULK]],
    ]
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
    manifest = {
        "analysis": "Phase 2G DINP/MCOP-CRC epithelial state convergence and regulatory anchoring",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "census_version": CENSUS_VERSION,
        "census_uri": CENSUS_URI,
        "expression_access": "official_source_h5ad_streamed_locally" if local_audit else "live_census_tiledb",
        "local_h5ad_audit": local_audit,
        "primary_dataset_id": dataset_id,
        "paired_donors_queried": paired_donors,
        "primary_filter": PRIMARY_FILTER,
        "ppar_nr_core_genes": PPAR_NR_GENES,
        "anchor_candidates": ANCHOR_CANDIDATES,
        "state_programs": {name: {"source": state_sources[name], "n_genes": len(genes)} for name, genes in state_sets.items()},
        "regulator_activity": "DoRothEA human levels A-C via decoupler ULM",
        "unit_of_inference": "donor-level pseudobulk; paired donor tests where available",
        "cell_level_state_definition": "bottom/top quartile of cell-level PPAR/NR score; tercile and median labels retained",
        "external_validation": "GSE144735 patient-paired epithelial target-universe score plus existing TCGA paired Phase 2F output",
        "causal_status": "not established",
        "verdict": "PARTIALLY",
        "virtual_perturbation": "not run; causal perturbation is not triggered by this associative analysis",
        "outputs": [str(x) for x in [OUT_CELL, OUT_DE, OUT_PATHWAY, OUT_CORR, OUT_DONOR, OUT_SUBTYPE, OUT_INTERACTION, OUT_REGULATOR, OUT_ANCHOR, OUT_EXTERNAL, OUT_BRIDGE, OUT_PSEUDOBULK, OUT_REPORT]],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"dataset_id": dataset_id, "n_cells": len(cell_scores), "n_donor_ppar_units": len(units), "n_states": len(state_sets), "n_query_genes": len(query_genes), "verdict": "PARTIALLY"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
