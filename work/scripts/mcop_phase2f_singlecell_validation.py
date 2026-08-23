"""Phase 2F-B: CELLxGENE Census single-cell validation of the DINP-axis program.

The analysis is intentionally conservative:

* Census release is pinned to ``2025-11-08`` and every metadata/expression
  query includes ``is_primary_data == True``.
* Census cell metadata are audited before expression is queried.  Labels are
  discovered from the release metadata rather than silently hard-coded.
* The inferential unit is a donor-level pseudobulk, never an individual cell.
* The primary comparison is tumor-derived epithelial versus normal colon
  epithelial.  ``tumor-derived`` is used deliberately: this script does not
  infer malignancy from cell type labels and does not run CNV inference.
* PPAR/nuclear-receptor is primary; RELA+STAT3 and the nine-gene DINP-axis
  score are secondary.  Myeloid, fibroblast and endothelial compartments are
  localization analyses.
* Dataset-level effects and leave-one-dataset-out analyses are emitted so a
  large dataset cannot silently become the whole result.

The Census stores raw counts but does not provide a universal library-size
column in the standardized observation metadata.  Therefore pseudobulk
counts are normalized to the total counts of the nine frozen target genes,
then log1p-transformed and z-scored within dataset and compartment.  This is a
transparent targeted score, not a substitute for full-transcriptome library
normalization.  The limitation is recorded in the report.

Run inside the pinned WSL environment, for example::

    /opt/cellxgene-census/bin/python \\
        work/scripts/mcop_phase2f_singlecell_validation.py

Outputs are written below ``outputs/`` with the ``mcop_phase2f_singlecell_``
prefix.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd
import tiledbsoma as soma
from scipy.sparse import issparse
from scipy.stats import mannwhitneyu, wilcoxon

try:
    from cellxgene_census import get_anndata, get_obs, open_soma
except ImportError as exc:  # pragma: no cover - depends on the WSL runtime
    raise SystemExit(
        "cellxgene-census is required. Run this script in the pinned WSL "
        "environment described in work/cellxgene_census_wsl_setup.md."
    ) from exc


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs"
CENSUS_VERSION = "2025-11-08"
CENSUS_URI = "s3://cellxgene-census-public-us-west-2/cell-census/2025-11-08/soma/"
ORGANISM = "homo_sapiens"
PRIMARY_FILTER = "is_primary_data == True"

NR_GENES = ["PPARA", "PPARD", "PPARG", "NR1I2", "NR1I3", "NR1H2", "NR1H3"]
INFLAMMATORY_GENES = ["RELA", "STAT3"]
ALL_GENES = NR_GENES + INFLAMMATORY_GENES
SCORE_COLUMNS = {
    "PPAR_NR": "PPAR_nuclear_receptor_score",
    "RELA_STAT3": "RELA_STAT3_score",
    "DINP_AXIS_9G": "DINP_axis_9_gene_score",
}
COMPARTMENTS = ["epithelial", "myeloid", "fibroblast", "endothelial"]
OBS_COLUMNS = [
    "dataset_id",
    "donor_id",
    "disease",
    "tissue",
    "tissue_general",
    "cell_type",
    "is_primary_data",
    "assay",
    "sex",
]


def q(value: object) -> str:
    """Quote a value for a SOMA value_filter expression."""

    return "'" + str(value).replace("'", "''") + "'"


def in_filter(field: str, values: Iterable[object]) -> str:
    values = list(values)
    return f"{field} in [{', '.join(q(value) for value in values)}]"


def retry_call(function: Callable[[], object], label: str, attempts: int = 5) -> object:
    """Retry remote Census reads, which can transiently time out."""

    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            return function()
        except Exception as exc:  # noqa: BLE001 - remote clients raise several types
            last_error = exc
            if attempt == attempts:
                break
            wait_seconds = min(30, 5 * attempt)
            print(f"[retry {attempt}/{attempts - 1}] {label}: {type(exc).__name__}; waiting {wait_seconds}s")
            time.sleep(wait_seconds)
    raise RuntimeError(f"Census read failed after {attempts} attempts: {label}: {last_error}") from last_error


def open_census():
    return retry_call(
        # Use the release-pinned S3 URI directly.  This avoids an otherwise
        # unnecessary release-directory HTTP request and makes the fixed
        # release robust to transient metadata-directory timeouts.
        lambda: open_soma(uri=CENSUS_URI),
        f"open Census {CENSUS_VERSION} via pinned S3 URI",
        attempts=6,
    )


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")


def lower_series(values: pd.Series) -> pd.Series:
    return values.astype("string").fillna("").str.lower()


def load_release_metadata(census) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()
    try:
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()
    except Exception as exc:  # noqa: BLE001 - dataset table is auxiliary to the cell audit
        print(f"[warning] release datasets table unavailable: {type(exc).__name__}")
        datasets = pd.DataFrame()
    summary_path = OUTPUT / "mcop_phase2f_singlecell_release_summary_cell_counts.csv"
    datasets_path = OUTPUT / "mcop_phase2f_singlecell_release_datasets.csv"
    if not summary_path.exists() or summary_path.stat().st_size == 0:
        summary.to_csv(summary_path, index=False)
    if not datasets_path.exists() or datasets_path.stat().st_size == 0:
        datasets.to_csv(datasets_path, index=False)
    return summary, datasets


def discover_labels(summary: pd.DataFrame) -> dict[str, object]:
    human = summary[summary["organism"].eq(ORGANISM)].copy()
    tissue_general = human[human["category"].eq("tissue_general")].copy()
    disease = human[human["category"].eq("disease")].copy()

    tissue_text = lower_series(tissue_general["label"])
    tissue_candidates = tissue_general.loc[
        tissue_text.str.contains(r"colon|large intestine|rectum", regex=True), "label"
    ].drop_duplicates().astype(str).tolist()
    preferred_tissue = [
        label for label in ["colon", "large intestine", "rectum"] if label in tissue_candidates
    ]
    if not preferred_tissue:
        preferred_tissue = tissue_candidates[:5]

    disease_text = lower_series(disease["label"])
    tumor_labels = disease.loc[
        disease_text.eq("colon adenocarcinoma")
        | disease_text.str.contains("colorectal carcinoma", regex=False),
        "label",
    ].drop_duplicates().astype(str).tolist()
    normal_labels = disease.loc[disease_text.eq("normal"), "label"].drop_duplicates().astype(str).tolist()

    return {
        "organism": ORGANISM,
        "tissue_general_candidates": tissue_candidates,
        "preferred_tissue_general": preferred_tissue,
        "tumor_disease_labels": tumor_labels,
        "normal_disease_labels": normal_labels,
        "summary_rows_human": int(len(human)),
    }


def query_obs(census, value_filter: str, label: str) -> pd.DataFrame:
    frame = retry_call(
        lambda: get_obs(census, ORGANISM, value_filter=value_filter, column_names=OBS_COLUMNS),
        label,
    )
    assert isinstance(frame, pd.DataFrame)
    return frame


def is_relevant_tissue(frame: pd.DataFrame) -> pd.Series:
    tissue = lower_series(frame.get("tissue", pd.Series(index=frame.index, dtype="string")))
    tissue_general = lower_series(frame.get("tissue_general", pd.Series(index=frame.index, dtype="string")))
    return (
        tissue_general.str.contains(r"colon|large intestine|rectum", regex=True)
        | tissue.str.contains(r"colon|large intestine|rectum", regex=True)
    )


def classify_group(disease: pd.Series, tumor_labels: set[str], normal_labels: set[str]) -> pd.Series:
    disease_text = disease.astype("string").fillna("")
    tumor = disease_text.isin(tumor_labels) | disease_text.str.contains("colorectal carcinoma", regex=False)
    normal = disease_text.isin(normal_labels)
    output = pd.Series("outside_scope", index=disease.index, dtype="string")
    output.loc[tumor] = "tumor"
    output.loc[normal] = "normal"
    return output


def fetch_relevant_obs(census, labels: dict[str, object]) -> pd.DataFrame:
    tumor_labels = [str(x) for x in labels["tumor_disease_labels"]]
    normal_labels = [str(x) for x in labels["normal_disease_labels"]]
    if not tumor_labels:
        raise RuntimeError("No colon adenocarcinoma/colorectal carcinoma disease label was found in the release.")
    if not normal_labels:
        raise RuntimeError("No normal disease label was found in the release.")

    tumor_frames = []
    for disease_label in tumor_labels:
        value_filter = f"{PRIMARY_FILTER} and disease == {q(disease_label)}"
        frame = query_obs(census, value_filter, f"tumor obs: {disease_label}")
        if not frame.empty:
            frame["group"] = "tumor"
            frame["query_disease_label"] = disease_label
            tumor_frames.append(frame)

    # Normal cells are scoped to colon/large intestine at query time.  A broad
    # fallback is retained because some datasets use tissue rather than
    # tissue_general for colon annotation; the pandas tissue gate is then
    # applied below.
    tissue_labels = [str(x) for x in labels["preferred_tissue_general"]]
    normal_filter = f"{PRIMARY_FILTER} and disease == {q(normal_labels[0])}"
    if tissue_labels:
        normal_filter += f" and {in_filter('tissue_general', tissue_labels)}"
    normal_frame = query_obs(census, normal_filter, "normal colon obs")
    if normal_frame.empty and tissue_labels:
        fallback_frames = []
        for tissue_label in tissue_labels:
            fallback_filter = f"{PRIMARY_FILTER} and disease == {q(normal_labels[0])} and tissue_general == {q(tissue_label)}"
            candidate = query_obs(census, fallback_filter, f"normal fallback tissue: {tissue_label}")
            if not candidate.empty:
                fallback_frames.append(candidate)
        if fallback_frames:
            normal_frame = pd.concat(fallback_frames, ignore_index=True)
    if not normal_frame.empty:
        normal_frame["group"] = "normal"
        normal_frame["query_disease_label"] = normal_labels[0]

    frames = tumor_frames + ([normal_frame] if not normal_frame.empty else [])
    if not frames:
        raise RuntimeError("No relevant Census observations were returned.")
    obs = pd.concat(frames, ignore_index=True)
    obs["group"] = classify_group(
        obs["disease"],
        set(tumor_labels),
        set(normal_labels),
    )
    obs["is_relevant_tissue"] = is_relevant_tissue(obs)
    # Tumor disease is already a CRC-specific inclusion criterion; normal
    # disease needs the explicit colon/rectum tissue gate.
    keep = obs["group"].eq("tumor") | (obs["group"].eq("normal") & obs["is_relevant_tissue"])
    obs = obs.loc[keep].copy()
    obs["dataset_id"] = obs["dataset_id"].astype("string")
    obs["donor_id"] = obs["donor_id"].astype("string")
    obs["disease"] = obs["disease"].astype("string")
    obs["cell_type"] = obs["cell_type"].astype("string")
    obs["tissue_general"] = obs["tissue_general"].astype("string")
    obs["tissue"] = obs["tissue"].astype("string")
    obs["compartment"] = obs["cell_type"].map(classify_compartment)
    if not obs["is_primary_data"].fillna(False).all():
        raise AssertionError("The primary-data filter did not hold for all returned observations.")
    return obs


def classify_compartment(cell_type: object) -> str:
    value = str(cell_type).lower()
    if re.search(r"epithelial|colonocyte|enterocyte|goblet|paneth|enteroendocrine|tuft|best4|transit amplifying|stem cell of colon", value):
        return "epithelial"
    if re.search(r"macrophage|monocyte|dendritic|neutrophil|myeloid", value):
        return "myeloid"
    if re.search(r"fibroblast|stromal", value):
        return "fibroblast"
    if re.search(r"endothelial", value):
        return "endothelial"
    return "other"


def make_dataset_filter(dataset_id: str, tissue_labels: list[str], disease_labels: list[str], cell_type_labels: list[str]) -> str:
    value_filter = f"{PRIMARY_FILTER} and dataset_id == {q(dataset_id)}"
    if tissue_labels:
        value_filter += f" and {in_filter('tissue_general', tissue_labels)}"
    if disease_labels:
        value_filter += f" and {in_filter('disease', disease_labels)}"
    if cell_type_labels:
        value_filter += f" and {in_filter('cell_type', cell_type_labels)}"
    return value_filter


def fetch_dataset_anndata(census, dataset_id: str, tissue_labels: list[str], disease_labels: list[str]):
    var_filter = in_filter("feature_name", ALL_GENES)
    value_filter = make_dataset_filter(dataset_id, tissue_labels, disease_labels, [])
    try:
        adata = retry_call(
            lambda: get_anndata(
                census,
                ORGANISM,
                X_name="raw",
                obs_value_filter=value_filter,
                var_value_filter=var_filter,
                obs_column_names=OBS_COLUMNS,
                var_column_names=["feature_name", "feature_id"],
            ),
            f"expression dataset {dataset_id}",
        )
        if adata.n_obs > 0:
            return adata
    except Exception as scoped_error:  # noqa: BLE001 - fallback is intentional
        print(f"[dataset fallback] {dataset_id}: scoped query failed: {type(scoped_error).__name__}")

    # Some older records have tissue_general values that differ from the
    # summary label. Retry without tissue_general, but retain the disease
    # filter so a large multi-tissue dataset cannot flood memory.
    return retry_call(
        lambda: get_anndata(
            census,
            ORGANISM,
            X_name="raw",
            obs_value_filter=f"{PRIMARY_FILTER} and dataset_id == {q(dataset_id)} and {in_filter('disease', disease_labels)}",
            var_value_filter=var_filter,
            obs_column_names=OBS_COLUMNS,
            var_column_names=["feature_name", "feature_id"],
        ),
        f"expression fallback dataset {dataset_id}",
    )


def read_axis_metadata(reader, joinids) -> pd.DataFrame:
    """Read a TileDB-SOMA axis iterator and retain its coordinate."""

    frames = [table.to_pandas() for table in reader]
    frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if "soma_joinid" not in frame.columns:
        joinids = np.asarray(joinids)
        if len(joinids) != len(frame):
            raise RuntimeError("Axis metadata row count did not match the axis join-id count.")
        frame.insert(0, "soma_joinid", joinids)
    return frame


def fetch_dataset_pseudobulk(
    census,
    dataset_id: str,
    tissue_labels: list[str],
    disease_labels: list[str],
    cell_type_labels: list[str],
    tumor_labels: set[str],
    normal_labels: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Process one dataset out-of-core using an axis query.

    Only the nine frozen genes are read from X.  The observation axis is read
    as metadata, while the sparse X table is streamed in batches and summed
    directly into donor/group/compartment pseudobulk accumulators.
    """

    var_filter = in_filter("feature_name", ALL_GENES)
    scoped_filter = make_dataset_filter(dataset_id, tissue_labels, disease_labels, cell_type_labels)
    fallback_filter = (
        f"{PRIMARY_FILTER} and dataset_id == {q(dataset_id)} and "
        f"{in_filter('disease', disease_labels)} and {in_filter('cell_type', cell_type_labels)}"
    )
    experiment = census["census_data"][ORGANISM]

    def process_filter(value_filter: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        with experiment.axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter=value_filter),
            var_query=soma.AxisQuery(value_filter=var_filter),
        ) as query:
            obs = read_axis_metadata(query.obs(column_names=OBS_COLUMNS), query.obs_joinids())
            if obs.empty:
                return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            var = read_axis_metadata(query.var(column_names=["feature_name", "feature_id"]), query.var_joinids())
            var["feature_name"] = var["feature_name"].astype(str)
            gene_map = {
                int(joinid): gene
                for joinid, gene in zip(var["soma_joinid"], var["feature_name"])
                if gene in ALL_GENES
            }
            if not gene_map:
                return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

            obs["dataset_id"] = obs.get("dataset_id", dataset_id).astype("string")
            obs["donor_id"] = obs["donor_id"].astype("string")
            obs["group"] = classify_group(obs["disease"], tumor_labels, normal_labels)
            obs["is_relevant_tissue"] = is_relevant_tissue(obs)
            obs["compartment"] = obs["cell_type"].map(classify_compartment)
            keep = obs["group"].isin(["tumor", "normal"]) & obs["is_relevant_tissue"] & nonmissing_identifier(obs["donor_id"])
            obs = obs.loc[keep].copy()
            if obs.empty:
                return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            obs["donor_key"] = obs["dataset_id"].astype(str) + "::" + obs["donor_id"].astype(str)
            obs_lookup = obs.set_index("soma_joinid", drop=False)

            cell_audit = (
                obs.groupby(["dataset_id", "group", "compartment"], observed=True)
                .agg(
                    n_cells=("donor_key", "size"),
                    n_donors=("donor_key", "nunique"),
                    n_cell_types=("cell_type", "nunique"),
                    n_tissues=("tissue", "nunique"),
                )
                .reset_index()
            )
            donor_audit = (
                obs.groupby(["dataset_id", "donor_key", "donor_id", "group", "compartment"], observed=True)
                .agg(n_cells=("donor_key", "size"), n_cell_types=("cell_type", "nunique"))
                .reset_index()
            )
            count_accumulator: dict[tuple[str, str, str, str, str], float] = {}
            for table in query.X("raw").tables():
                batch = table.to_pandas()
                if batch.empty:
                    continue
                cell_ids = batch["soma_dim_0"].astype(np.int64).to_numpy()
                batch_meta = obs_lookup.reindex(cell_ids)
                valid = batch_meta["donor_key"].notna().to_numpy()
                if not valid.any():
                    continue
                batch = batch.loc[valid].copy()
                batch_meta = batch_meta.loc[valid]
                batch["gene"] = batch["soma_dim_1"].astype(np.int64).map(gene_map)
                valid_gene = batch["gene"].notna().to_numpy()
                if not valid_gene.any():
                    continue
                batch = batch.loc[valid_gene].copy()
                batch_meta = batch_meta.loc[valid_gene]
                batch["donor_key"] = batch_meta["donor_key"].to_numpy()
                batch["donor_id"] = batch_meta["donor_id"].to_numpy()
                batch["group"] = batch_meta["group"].to_numpy()
                batch["compartment"] = batch_meta["compartment"].to_numpy()
                batch["value"] = pd.to_numeric(batch["soma_data"], errors="coerce").fillna(0.0)
                grouped = batch.groupby(["dataset_id", "donor_key", "donor_id", "group", "compartment", "gene"], observed=True)["value"].sum() if "dataset_id" in batch.columns else batch.groupby(["donor_key", "donor_id", "group", "compartment", "gene"], observed=True)["value"].sum()
                for key, value in grouped.items():
                    if len(key) == 5:
                        key = (str(dataset_id),) + tuple(str(x) for x in key)
                    else:
                        key = tuple(str(x) for x in key)
                    count_accumulator[key] = count_accumulator.get(key, 0.0) + float(value)

            rows = []
            donor_groups = obs.groupby(["dataset_id", "donor_key", "donor_id", "group", "compartment"], observed=True)
            for keys, subset in donor_groups:
                dataset, donor_key, donor_id, group, compartment = [str(x) for x in keys]
                row = {
                    "dataset_id": dataset,
                    "donor_key": donor_key,
                    "donor_id": donor_id,
                    "group": group,
                    "compartment": compartment,
                    "n_cells": int(len(subset)),
                }
                for gene in ALL_GENES:
                    row[gene] = count_accumulator.get((dataset, donor_key, donor_id, group, compartment, gene), 0.0)
                rows.append(row)
            return pd.DataFrame(rows), cell_audit, donor_audit

    pseudobulk, cell_audit, donor_audit = process_filter(scoped_filter)
    if pseudobulk.empty and scoped_filter != fallback_filter:
        pseudobulk, cell_audit, donor_audit = process_filter(fallback_filter)
    return pseudobulk, cell_audit, donor_audit


def matrix_to_frame(adata) -> pd.DataFrame:
    matrix = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    names = adata.var["feature_name"].astype(str).tolist()
    output = pd.DataFrame(0.0, index=np.arange(adata.n_obs), columns=ALL_GENES)
    for col_idx, name in enumerate(names):
        if name in output.columns:
            output[name] = output[name].to_numpy(float) + np.asarray(matrix[:, col_idx]).reshape(-1)
    return output


def nonmissing_identifier(values: pd.Series) -> pd.Series:
    text = values.astype("string")
    return text.notna() & ~text.str.lower().isin({"", "nan", "none", "na"})


def pseudobulk_from_adata(adata, dataset_id: str, tumor_labels: set[str], normal_labels: set[str], tissue_labels: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expression = matrix_to_frame(adata)
    meta = adata.obs.reset_index(drop=True).copy()
    meta["dataset_id"] = meta.get("dataset_id", dataset_id).astype("string")
    meta["donor_id"] = meta["donor_id"].astype("string")
    meta["group"] = classify_group(meta["disease"], tumor_labels, normal_labels)
    meta["is_relevant_tissue"] = is_relevant_tissue(meta)
    meta["compartment"] = meta["cell_type"].map(classify_compartment)
    keep = meta["group"].isin(["tumor", "normal"]) & meta["is_relevant_tissue"] & nonmissing_identifier(meta["donor_id"])
    meta = meta.loc[keep].reset_index(drop=True)
    expression = expression.loc[keep.to_numpy()].reset_index(drop=True)
    if meta.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    meta["donor_key"] = meta["dataset_id"].astype(str) + "::" + meta["donor_id"].astype(str)
    cell_audit = (
        meta.groupby(["dataset_id", "group", "compartment"], observed=True)
        .agg(
            n_cells=("donor_key", "size"),
            n_donors=("donor_key", "nunique"),
            n_cell_types=("cell_type", "nunique"),
            n_tissues=("tissue", "nunique"),
        )
        .reset_index()
    )
    donor_audit = (
        meta.groupby(["dataset_id", "donor_key", "donor_id", "group", "compartment"], observed=True)
        .agg(n_cells=("donor_key", "size"), n_cell_types=("cell_type", "nunique"))
        .reset_index()
    )

    rows = []
    for keys, indices in meta.groupby(["dataset_id", "donor_key", "donor_id", "group", "compartment"], observed=True).groups.items():
        dataset, donor_key, donor_id, group, compartment = keys
        summed = expression.loc[indices, ALL_GENES].sum(axis=0)
        row = {
            "dataset_id": str(dataset),
            "donor_key": str(donor_key),
            "donor_id": str(donor_id),
            "group": str(group),
            "compartment": str(compartment),
            "n_cells": int(len(indices)),
        }
        row.update({gene: float(summed[gene]) for gene in ALL_GENES})
        rows.append(row)
    return pd.DataFrame(rows), cell_audit, donor_audit


def score_pseudobulk(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    scored = pseudobulk.copy()
    target_total = scored[ALL_GENES].sum(axis=1)
    scored["target_gene_total_counts"] = target_total
    valid_total = target_total.gt(0)
    log_cpm = pd.DataFrame(np.nan, index=scored.index, columns=ALL_GENES, dtype=float)
    log_cpm.loc[valid_total, ALL_GENES] = np.log1p(
        scored.loc[valid_total, ALL_GENES].div(target_total[valid_total], axis=0) * 1_000_000
    )

    score_frames = []
    for (_, _), indices in scored.groupby(["dataset_id", "compartment"], observed=True).groups.items():
        block = log_cpm.loc[indices, ALL_GENES]
        z = (block - block.mean(axis=0)) / block.std(axis=0, ddof=1).replace(0, np.nan)
        for label, genes in {
            "PPAR_NR": NR_GENES,
            "RELA_STAT3": INFLAMMATORY_GENES,
            "DINP_AXIS_9G": ALL_GENES,
        }.items():
            score_frames.append(
                pd.DataFrame(
                    {
                        "index": indices,
                        SCORE_COLUMNS[label]: z[genes].mean(axis=1),
                        f"{SCORE_COLUMNS[label]}_n_genes": z[genes].notna().sum(axis=1),
                    }
                )
            )
    score_columns = [
        column
        for label in SCORE_COLUMNS.values()
        for column in (label, f"{label}_n_genes")
    ]
    score_index = pd.DataFrame(np.nan, index=scored.index, columns=score_columns)
    for frame in score_frames:
        frame_index = frame["index"].to_numpy()
        for column in frame.columns:
            if column == "index":
                continue
            score_index.loc[frame_index, column] = frame[column].to_numpy()
    for column in score_index.columns:
        scored[column] = score_index[column]
    return scored


def test_two_groups(tumor: pd.Series, normal: pd.Series) -> tuple[float, float]:
    tumor = pd.to_numeric(tumor, errors="coerce").dropna().to_numpy(float)
    normal = pd.to_numeric(normal, errors="coerce").dropna().to_numpy(float)
    if len(tumor) < 2 or len(normal) < 2:
        return np.nan, np.nan
    statistic, p_value = mannwhitneyu(tumor, normal, alternative="two-sided")
    return float(statistic), float(p_value)


def summarize_subset(subset: pd.DataFrame, compartment: str, score_name: str, dataset_dropped: str = "NONE") -> dict[str, object]:
    tumor = subset.loc[subset["group"].eq("tumor"), score_name]
    normal = subset.loc[subset["group"].eq("normal"), score_name]
    statistic, p_value = test_two_groups(tumor, normal)
    return {
        "compartment": compartment,
        "score": score_name,
        "dataset_dropped": dataset_dropped,
        "tumor_donors": int(tumor.notna().sum()),
        "normal_donors": int(normal.notna().sum()),
        "n_datasets": int(subset["dataset_id"].nunique()),
        "n_tumor_datasets": int(subset.loc[subset["group"].eq("tumor"), "dataset_id"].nunique()),
        "n_normal_datasets": int(subset.loc[subset["group"].eq("normal"), "dataset_id"].nunique()),
        "tumor_median": float(tumor.median()) if tumor.notna().any() else np.nan,
        "normal_median": float(normal.median()) if normal.notna().any() else np.nan,
        "mean_delta_tumor_minus_normal": float(tumor.mean() - normal.mean()) if tumor.notna().any() and normal.notna().any() else np.nan,
        "median_delta_tumor_minus_normal": float(tumor.median() - normal.median()) if tumor.notna().any() and normal.notna().any() else np.nan,
        "mann_whitney_U": statistic,
        "p_value": p_value,
    }


def dataset_effects(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, compartment), subset in scored.groupby(["dataset_id", "compartment"], observed=True):
        for score_name in SCORE_COLUMNS.values():
            row = summarize_subset(subset, str(compartment), score_name)
            row["dataset_id"] = str(dataset)
            rows.append(row)
    return pd.DataFrame(rows)


def pooled_and_lodo(scored: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pooled_rows = []
    lodo_rows = []
    datasets = sorted(scored["dataset_id"].dropna().astype(str).unique())
    for compartment in [*COMPARTMENTS, "other"]:
        for score_name in SCORE_COLUMNS.values():
            subset = scored[scored["compartment"].eq(compartment)]
            pooled_rows.append(summarize_subset(subset, compartment, score_name))
            for dropped in datasets:
                reduced = subset[subset["dataset_id"].astype(str).ne(dropped)]
                lodo_rows.append(summarize_subset(reduced, compartment, score_name, dataset_dropped=dropped))
    return pd.DataFrame(pooled_rows), pd.DataFrame(lodo_rows)


def direction_summary(effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (compartment, score), subset in effects.groupby(["compartment", "score"], observed=True):
        eligible = subset[subset["tumor_donors"].ge(2) & subset["normal_donors"].ge(2)].copy()
        deltas = pd.to_numeric(eligible["mean_delta_tumor_minus_normal"], errors="coerce").dropna()
        rows.append(
            {
                "compartment": compartment,
                "score": score,
                "eligible_datasets": int(len(deltas)),
                "positive_datasets": int((deltas > 0).sum()),
                "negative_datasets": int((deltas < 0).sum()),
                "positive_fraction": float((deltas > 0).mean()) if len(deltas) else np.nan,
                "median_dataset_delta": float(deltas.median()) if len(deltas) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def paired_donor_contrasts(scored: pd.DataFrame) -> pd.DataFrame:
    """Compare tumor and normal pseudobulk from the same donor when available."""

    rows = []
    for (dataset, compartment), subset in scored.groupby(["dataset_id", "compartment"], observed=True):
        for score_name in SCORE_COLUMNS.values():
            values = subset[["donor_id", "group", score_name]].dropna()
            pivot = values.pivot_table(index="donor_id", columns="group", values=score_name, aggfunc="mean")
            if "tumor" not in pivot.columns or "normal" not in pivot.columns:
                differences = np.array([], dtype=float)
            else:
                differences = (pivot["tumor"] - pivot["normal"]).dropna().to_numpy(float)
            if len(differences) >= 2 and not np.allclose(differences, 0):
                statistic, p_value = wilcoxon(differences, alternative="two-sided", method="auto")
            else:
                statistic, p_value = np.nan, np.nan
            rows.append(
                {
                    "dataset_id": str(dataset),
                    "compartment": str(compartment),
                    "score": score_name,
                    "paired_donors": int(len(differences)),
                    "mean_delta_tumor_minus_normal": float(np.mean(differences)) if len(differences) else np.nan,
                    "median_delta_tumor_minus_normal": float(np.median(differences)) if len(differences) else np.nan,
                    "wilcoxon_statistic": float(statistic) if np.isfinite(statistic) else np.nan,
                    "p_value": float(p_value) if np.isfinite(p_value) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def write_report(
    labels: dict[str, object],
    obs: pd.DataFrame,
    cell_audit: pd.DataFrame,
    donor_audit: pd.DataFrame,
    pooled: pd.DataFrame,
    lodo: pd.DataFrame,
    effects: pd.DataFrame,
    directions: pd.DataFrame,
    paired: pd.DataFrame,
    error: str | None = None,
    scope_note: str | None = None,
) -> None:
    if pooled.empty:
        pooled = pd.DataFrame(columns=["compartment", "score", "tumor_donors", "normal_donors", "median_delta_tumor_minus_normal", "p_value"])
    if lodo.empty:
        lodo = pd.DataFrame(columns=["compartment", "score", "tumor_donors", "normal_donors", "mean_delta_tumor_minus_normal"])
    if directions.empty:
        directions = pd.DataFrame(columns=["compartment", "score", "eligible_datasets", "positive_datasets", "negative_datasets", "positive_fraction"])
    if paired.empty:
        paired = pd.DataFrame(columns=["dataset_id", "compartment", "score", "paired_donors", "median_delta_tumor_minus_normal", "p_value"])
    primary_score = SCORE_COLUMNS["PPAR_NR"]
    primary_direction = directions[(directions["compartment"].eq("epithelial")) & directions["score"].eq(primary_score)]
    primary_lodo = lodo[(lodo["compartment"].eq("epithelial")) & lodo["score"].eq(primary_score)]
    if not primary_direction.empty:
        drow = primary_direction.iloc[0]
        eligible = int(drow["eligible_datasets"])
        all_positive = eligible >= 2 and int(drow["positive_datasets"]) == eligible
    else:
        drow = None
        eligible = 0
        all_positive = False
    lodo_valid = primary_lodo[primary_lodo["tumor_donors"].ge(2) & primary_lodo["normal_donors"].ge(2)]
    lodo_stable = bool(not lodo_valid.empty and (pd.to_numeric(lodo_valid["mean_delta_tumor_minus_normal"], errors="coerce") > 0).all())

    compartment_pass = directions[
        directions["eligible_datasets"].ge(2) & directions["positive_datasets"].eq(directions["eligible_datasets"])
    ]
    if error:
        status = "NOT_RUN"
        interpretation = f"Census execution did not complete: {error}"
    elif scope_note:
        status = "PARTIAL"
        interpretation = "The targeted dataset result is informative but cannot trigger the frozen multi-dataset mechanism gate."
    elif eligible >= 2 and all_positive and lodo_stable:
        status = "GREEN"
        interpretation = "Epithelial donor-level PPAR/NR scores are directionally concordant across datasets and remain positive after dataset leave-one-out."
    elif not compartment_pass.empty:
        status = "YELLOW"
        interpretation = "The epithelial primary gate is not fully stable, but at least one compartment shows directionally concordant donor-level remodeling across datasets. Treat this as a localization hypothesis, not tumor-cell intrinsic proof."
    else:
        status = "RED"
        interpretation = "No compartment meets the frozen multi-dataset direction criterion. Stop treating PPAR/NR as the main mechanism and retain it as exploratory only."

    observed_compartments = sorted(
        set(
            cell_audit.loc[
                pd.to_numeric(cell_audit.get("n_cells", pd.Series(dtype=float)), errors="coerce").gt(0),
                "compartment",
            ].dropna().astype(str)
        )
        | set(
            pooled.loc[
                pd.to_numeric(pooled.get("tumor_donors", pd.Series(dtype=float)), errors="coerce").gt(0)
                | pd.to_numeric(pooled.get("normal_donors", pd.Series(dtype=float)), errors="coerce").gt(0),
                "compartment",
            ].dropna().astype(str)
        )
    )
    if scope_note:
        run_description = (
            "This is a targeted partial run: the paired Census dataset was processed across "
            f"the available compartments ({', '.join(observed_compartments) or 'none'}). "
            "It cannot trigger the frozen multi-dataset GREEN/YELLOW/RED mechanism gate; "
            "the independent dataset replication remains a separate analysis."
        )
    else:
        run_description = "This is the full multi-dataset, four-compartment Census run."

    lines = [
        "# MCOP–CRC Phase 2F-B：CELLxGENE Census 单细胞验证",
        "",
        f"## 当前判定：**{status}**",
        "",
        run_description,
        "",
        interpretation,
        "",
        "本轮的疾病细胞标签只能支持 **tumor-derived epithelial**；本脚本没有把它写成 malignant epithelial，也没有做 CNV 推断或使用未经核验的 malignant 标签。",
        "",
        "## 冻结规则",
        "",
        f"- Census release: `{CENSUS_VERSION}`；organism: `{ORGANISM}`。",
        "- 所有 metadata/expression 查询均包含 `is_primary_data == True`。",
        "- 统计单位是 donor-level pseudobulk；没有把 cell 当作独立样本做 P 值。",
        f"- 主 score：`{primary_score}`；基因：`{', '.join(NR_GENES)}`。",
        f"- Secondary：`{SCORE_COLUMNS['RELA_STAT3']}`、`{SCORE_COLUMNS['DINP_AXIS_9G']}`。",
        "- Myeloid、fibroblast、endothelial 仅用于 compartment localization。",
        "",
        "## Census 标签审计",
        "",
        f"- Tumor disease labels: `{'; '.join(labels.get('tumor_disease_labels', []))}`",
        f"- Normal disease labels: `{'; '.join(labels.get('normal_disease_labels', []))}`",
        f"- Tissue-general candidates: `{'; '.join(labels.get('tissue_general_candidates', []))}`",
        f"- Relevant observations after primary/tissue gate: **{len(obs):,}**",
        f"- Datasets: **{obs['dataset_id'].nunique():,}**；donors with usable IDs: **{obs['donor_id'].nunique():,}**；cell types: **{obs['cell_type'].nunique():,}**",
        "",
        "| group | cells | datasets | donor IDs |",
        "|---|---:|---:|---:|",
    ]
    for group, subset in obs.groupby("group", observed=True):
        lines.append(f"| {group} | {len(subset):,} | {subset['dataset_id'].nunique():,} | {subset['donor_id'].nunique():,} |")

    lines += [
        "",
        "## Primary donor-level result",
        "",
        "The pooled result below is a donor-level descriptive comparison. Dataset-level effects and leave-one-dataset-out results are the required stability checks; pooled P values are not cell-level P values.",
        "",
        "| compartment | score | tumor donors | normal donors | median delta | P |",
        "|---|---|---:|---:|---:|---:|",
    ]
    primary_rows = pooled[(pooled["compartment"].eq("epithelial"))].copy()
    for _, row in primary_rows.iterrows():
        p = "NA" if not np.isfinite(row["p_value"]) else f"{row['p_value']:.3g}"
        delta = "NA" if not np.isfinite(row["median_delta_tumor_minus_normal"]) else f"{row['median_delta_tumor_minus_normal']:.3f}"
        lines.append(f"| {row['compartment']} | {row['score']} | {int(row['tumor_donors'])} | {int(row['normal_donors'])} | {delta} | {p} |")

    lines += [
        "",
        "## Dataset direction check",
        "",
        "A dataset is eligible for this check only when it has at least two tumor and two normal donors in the same compartment.",
        "",
        "| compartment | score | eligible datasets | positive | negative | positive fraction |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for _, row in directions.iterrows():
        fraction = "NA" if not np.isfinite(row["positive_fraction"]) else f"{row['positive_fraction']:.2f}"
        lines.append(f"| {row['compartment']} | {row['score']} | {int(row['eligible_datasets'])} | {int(row['positive_datasets'])} | {int(row['negative_datasets'])} | {fraction} |")

    lines += [
        "",
        "## Paired-donor check",
        "",
        "When the same donor ID contains both tumor and normal observations, the paired donor result is reported separately. This is not substituted for the multi-dataset gate.",
        "",
        "| dataset | compartment | score | paired donors | median delta | P |",
        "|---|---|---|---:|---:|---:|",
    ]
    for _, row in paired.iterrows():
        delta = "NA" if not np.isfinite(row["median_delta_tumor_minus_normal"]) else f"{row['median_delta_tumor_minus_normal']:.3f}"
        p = "NA" if not np.isfinite(row["p_value"]) else f"{row['p_value']:.3g}"
        lines.append(f"| {row['dataset_id']} | {row['compartment']} | {row['score']} | {int(row['paired_donors'])} | {delta} | {p} |")

    lines += [
        "",
        "## Leave-one-dataset-out",
        "",
        "`mcop_phase2f_singlecell_leave_one_dataset_out.csv` contains the full table. For the primary epithelial PPAR/NR score, the report uses only leave-one-out rows with at least two donors per group; this prevents a formally computed but uninformative result from being called stable.",
        "",
        "## Interpretation boundaries",
        "",
        "- A stable epithelial result would support a CRC-associated epithelial program, not exposure causality or DINP-to-CRC mediation.",
        "- A stable myeloid/fibroblast/endothelial result with unstable epithelial direction would redirect the mechanism toward microenvironmental PPAR/NR remodeling.",
        "- Tissue composition, dataset-specific annotation, treatment history and dissociation effects remain possible explanations. The analysis is not a substitute for prediagnostic tissue or perturbation validation.",
        "- Targeted normalization uses the nine frozen genes because the standardized Census observation metadata do not provide a universal total-UMI field. This is explicitly a targeted score and should not be overinterpreted as full transcriptome normalization.",
        "",
        "## Output files",
        "",
        "- `mcop_phase2f_singlecell_label_discovery.json`",
        "- `mcop_phase2f_singlecell_dataset_donor_cell_audit.csv`",
        "- `mcop_phase2f_singlecell_donor_audit.csv`",
        "- `mcop_phase2f_singlecell_donor_pseudobulk.csv`",
        "- `mcop_phase2f_singlecell_donor_scores.csv`",
        "- `mcop_phase2f_singlecell_dataset_effects.csv`",
        "- `mcop_phase2f_singlecell_pooled_contrasts.csv`",
        "- `mcop_phase2f_singlecell_leave_one_dataset_out.csv`",
        "- `mcop_phase2f_singlecell_direction_summary.csv`",
        "- `mcop_phase2f_singlecell_paired_donor_contrasts.csv`",
        "",
        "## Reproducibility",
        "",
        f"- Run UTC: `{datetime.now(timezone.utc).isoformat()}`",
        f"- Python: `{platform.python_version()}`",
        f"- Script: `{Path(__file__).relative_to(ROOT)}`",
        f"- NR genes: `{','.join(NR_GENES)}`",
        f"- All genes: `{','.join(ALL_GENES)}`",
    ]
    if drow is not None:
        lines += [
            "",
            f"Primary epithelial PPAR/NR eligible datasets={eligible}; positive={int(drow['positive_datasets'])}; leave-one-dataset-out stable={lodo_stable}.",
        ]
    if error:
        lines += ["", "## Execution error", "", f"`{error}`"]
    if scope_note:
        lines += ["", "## Scope note", "", scope_note]
    (OUTPUT / "mcop_phase2f_singlecell_validation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def empty_result_frames() -> tuple[pd.DataFrame, ...]:
    return tuple(pd.DataFrame() for _ in range(8))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-datasets", type=int, default=None, help="Optional limit for a smoke test; omit for the full audit.")
    parser.add_argument("--dataset-ids", nargs="+", default=None, help="Optional explicit dataset IDs for a resumable targeted run.")
    parser.add_argument("--compartments", nargs="+", choices=COMPARTMENTS, default=COMPARTMENTS, help="Compartments to query; default is all four frozen compartments.")
    parser.add_argument(
        "--reuse-observation-audit",
        action="store_true",
        help="Reuse the existing observation metadata audit instead of re-querying Census metadata.",
    )
    parser.add_argument(
        "--reuse-pseudobulk",
        action="store_true",
        help="Reuse the existing donor pseudobulk CSV and rerun scoring/reporting without querying Census expression.",
    )
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    obs = pd.DataFrame(columns=OBS_COLUMNS + ["group", "query_disease_label", "is_relevant_tissue"])
    cell_audit = pd.DataFrame()
    donor_audit = pd.DataFrame()
    pseudobulk = pd.DataFrame()
    scored = pd.DataFrame()
    effects = pd.DataFrame()
    pooled = pd.DataFrame()
    lodo = pd.DataFrame()
    directions = pd.DataFrame()
    paired = pd.DataFrame()
    labels: dict[str, object] = {}
    execution_error: str | None = None
    reused_pseudobulk = False

    try:
        with open_census() as census:
            audit_path = OUTPUT / "mcop_phase2f_singlecell_observation_metadata_audit.csv"
            label_path = OUTPUT / "mcop_phase2f_singlecell_label_discovery.json"
            if args.reuse_observation_audit and audit_path.exists() and label_path.exists():
                labels = json.loads(label_path.read_text(encoding="utf-8"))
                obs = pd.read_csv(audit_path)
                obs["compartment"] = obs["cell_type"].map(classify_compartment)
                print(f"[reuse] loaded observation audit: {len(obs):,} rows")
            else:
                summary, _datasets = load_release_metadata(census)
                labels = discover_labels(summary)
                write_json(label_path, labels)
                obs = fetch_relevant_obs(census, labels)
                obs.to_csv(audit_path, index=False)

            if args.reuse_pseudobulk:
                pseudobulk_path = OUTPUT / "mcop_phase2f_singlecell_donor_pseudobulk.csv"
                cell_audit_path = OUTPUT / "mcop_phase2f_singlecell_dataset_donor_cell_audit.csv"
                donor_audit_path = OUTPUT / "mcop_phase2f_singlecell_donor_audit.csv"
                if not pseudobulk_path.exists():
                    raise RuntimeError(f"Requested pseudobulk reuse but file is missing: {pseudobulk_path}")
                pseudobulk = pd.read_csv(pseudobulk_path)
                if cell_audit_path.exists():
                    cell_audit = pd.read_csv(cell_audit_path)
                if donor_audit_path.exists():
                    donor_audit = pd.read_csv(donor_audit_path)
                reused_pseudobulk = True
                print(f"[reuse] loaded donor pseudobulk: {len(pseudobulk):,} rows")
            else:
                dataset_ids = sorted(obs["dataset_id"].dropna().astype(str).unique())
                if args.dataset_ids:
                    requested = [str(x) for x in args.dataset_ids]
                    missing = sorted(set(requested) - set(dataset_ids))
                    if missing:
                        raise RuntimeError(f"Requested dataset IDs were not present in the metadata audit: {missing}")
                    dataset_ids = requested
                if args.max_datasets is not None:
                    dataset_ids = dataset_ids[: args.max_datasets]
                tissue_labels = [str(x) for x in labels.get("preferred_tissue_general", [])]
                tumor_labels = set(str(x) for x in labels.get("tumor_disease_labels", []))
                normal_labels = set(str(x) for x in labels.get("normal_disease_labels", []))
                expression_disease_labels = sorted(tumor_labels | normal_labels)
                pb_frames, cell_frames, donor_frames = [], [], []
                for index, dataset_id in enumerate(dataset_ids, start=1):
                    dataset_obs = obs[obs["dataset_id"].astype(str).eq(dataset_id)]
                    cell_type_labels = sorted(
                        dataset_obs.loc[
                            dataset_obs["compartment"].isin(args.compartments), "cell_type"
                        ]
                        .dropna()
                        .astype(str)
                        .unique()
                    )
                    if not cell_type_labels:
                        continue
                    print(
                        f"[{index}/{len(dataset_ids)}] querying {dataset_id} once "
                        f"for compartments={','.join(args.compartments)}; "
                        f"cell types={len(cell_type_labels)}"
                    )
                    pb, cells, donors = fetch_dataset_pseudobulk(
                        census,
                        dataset_id,
                        tissue_labels,
                        expression_disease_labels,
                        cell_type_labels,
                        tumor_labels,
                        normal_labels,
                    )
                    if not pb.empty:
                        pb = pb[pb["compartment"].isin(args.compartments)].copy()
                        pb_frames.append(pb)
                    if not cells.empty:
                        cells = cells[cells["compartment"].isin(args.compartments)].copy()
                        cell_frames.append(cells)
                    if not donors.empty:
                        donors = donors[donors["compartment"].isin(args.compartments)].copy()
                        donor_frames.append(donors)
                    print(f"    usable donor rows={len(pb):,}")

                pseudobulk = pd.concat(pb_frames, ignore_index=True) if pb_frames else pd.DataFrame()
                cell_audit = pd.concat(cell_frames, ignore_index=True) if cell_frames else pd.DataFrame()
                donor_audit = pd.concat(donor_frames, ignore_index=True) if donor_frames else pd.DataFrame()
            if pseudobulk.empty:
                raise RuntimeError("No donor-level pseudobulk rows remained after donor, disease, tissue and compartment gates.")
            scored = score_pseudobulk(pseudobulk)
            effects = dataset_effects(scored)
            pooled, lodo = pooled_and_lodo(scored)
            directions = direction_summary(effects)
            paired = paired_donor_contrasts(scored)
    except Exception as exc:  # noqa: BLE001 - report failed remote run with reproducibility metadata
        execution_error = f"{type(exc).__name__}: {exc}"
        print(f"[ERROR] {execution_error}")

    # Always write whatever the run produced, including an auditable failure
    # report rather than claiming a single-cell result that did not run.
    if not cell_audit.empty:
        cell_audit.to_csv(OUTPUT / "mcop_phase2f_singlecell_dataset_donor_cell_audit.csv", index=False)
    else:
        pd.DataFrame().to_csv(OUTPUT / "mcop_phase2f_singlecell_dataset_donor_cell_audit.csv", index=False)
    if not donor_audit.empty:
        donor_audit.to_csv(OUTPUT / "mcop_phase2f_singlecell_donor_audit.csv", index=False)
    else:
        pd.DataFrame().to_csv(OUTPUT / "mcop_phase2f_singlecell_donor_audit.csv", index=False)
    if not pseudobulk.empty:
        pseudobulk.to_csv(OUTPUT / "mcop_phase2f_singlecell_donor_pseudobulk.csv", index=False)
    else:
        pd.DataFrame().to_csv(OUTPUT / "mcop_phase2f_singlecell_donor_pseudobulk.csv", index=False)
    if not scored.empty:
        scored.to_csv(OUTPUT / "mcop_phase2f_singlecell_donor_scores.csv", index=False)
    else:
        pd.DataFrame().to_csv(OUTPUT / "mcop_phase2f_singlecell_donor_scores.csv", index=False)
    effects.to_csv(OUTPUT / "mcop_phase2f_singlecell_dataset_effects.csv", index=False)
    pooled.to_csv(OUTPUT / "mcop_phase2f_singlecell_pooled_contrasts.csv", index=False)
    lodo.to_csv(OUTPUT / "mcop_phase2f_singlecell_leave_one_dataset_out.csv", index=False)
    directions.to_csv(OUTPUT / "mcop_phase2f_singlecell_direction_summary.csv", index=False)
    full_scope = (
        not reused_pseudobulk
        and args.dataset_ids is None
        and args.max_datasets is None
        and set(args.compartments) == set(COMPARTMENTS)
    )
    scope_note = None if full_scope else (
        "This is a targeted/partial run (explicit dataset, compartment, or pseudobulk reuse restriction). "
        "Do not apply the frozen GREEN/YELLOW/RED mechanism gate until the full multi-dataset, "
        "four-compartment run is complete."
    )
    paired.to_csv(OUTPUT / "mcop_phase2f_singlecell_paired_donor_contrasts.csv", index=False)
    write_report(labels, obs, cell_audit, donor_audit, pooled, lodo, effects, directions, paired, execution_error, scope_note)
    write_json(
        OUTPUT / "mcop_phase2f_singlecell_validation_manifest.json",
        {
            "analysis": "MCOP-CRC Phase 2F-B CELLxGENE Census single-cell validation",
            "census_version": CENSUS_VERSION,
            "census_uri": CENSUS_URI,
            "organism": ORGANISM,
            "primary_data_filter": PRIMARY_FILTER,
            "nr_genes": NR_GENES,
            "inflammatory_genes": INFLAMMATORY_GENES,
            "all_genes": ALL_GENES,
            "unit_of_analysis": "donor-level pseudobulk",
            "tumor_labeling": "tumor-derived epithelial; no malignant annotation/CNV inference",
            "max_datasets": args.max_datasets,
            "dataset_ids": args.dataset_ids,
            "compartments": args.compartments,
            "reuse_pseudobulk": args.reuse_pseudobulk,
            "full_scope": full_scope,
            "n_relevant_observations": int(len(obs)),
            "n_datasets": int(obs["dataset_id"].nunique()) if not obs.empty else 0,
            "n_donors": int(obs["donor_id"].nunique()) if not obs.empty else 0,
            "execution_error": execution_error,
            "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    if execution_error:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
