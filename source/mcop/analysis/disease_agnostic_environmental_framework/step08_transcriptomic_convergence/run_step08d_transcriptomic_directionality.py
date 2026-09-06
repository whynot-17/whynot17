#!/usr/bin/env python3
"""Step 8D: frozen T2D network-module transcriptomic directionality.

This script downloads processed GEO series matrices and platform annotation
tables, applies frozen sample-label rules, scores only the Step 8C network
modules, and performs dataset-level directionality tests. It deliberately does
not pool tissues, expand gene sets, or select a flagship mechanism.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
RAW = HERE / "raw_geo"
OUT = HERE
NETWORK_MODULES = HERE.parent / "step08_t2d_convergence" / "t2d_step8c_network_modules.csv"
NETWORK_NODES = HERE.parent / "step08_t2d_convergence" / "t2d_step8c_network_nodes.csv"


DATASETS = [
    {
        "accession": "GSE23343",
        "series_group": "GSE23nnn",
        "tissue": "liver",
        "platform": "GPL570",
        "platform_annotation": "GPL570",
        "case_patterns": [r"disease status:\s*type 2 diabetes"],
        "control_patterns": [r"disease status:\s*normal glucose tolerance"],
        "exclusion_patterns": [],
        "source_note": "human fasting liver biopsy; T2D versus normal glucose tolerance",
    },
    {
        "accession": "GSE21340",
        "series_group": "GSE21nnn",
        "tissue": "skeletal_muscle",
        "platform": "GPL80",
        "platform_annotation": "GPL80",
        "case_patterns": [r"disease state:\s*dm\b"],
        "control_patterns": [r"disease state:\s*control\b"],
        "exclusion_patterns": [r"disease state:\s*replicate control"],
        "source_note": "skeletal muscle; T2D versus control; replicate controls excluded",
    },
    {
        "accession": "GSE71416",
        "series_group": "GSE71nnn",
        "tissue": "adipose",
        "platform": "GPL570",
        "platform_annotation": "GPL570",
        "case_patterns": [r"!Sample_title:.*obese\s+and\s+diabetic\b"],
        "control_patterns": [r"!Sample_title:.*obese\s+non[- ]diabetic\b|!Sample_title:.*obese\s+nondiabetic\b"],
        "exclusion_patterns": [],
        "source_note": "omental adipose tissue; obese diabetic versus obese non-diabetic",
    },
    {
        "accession": "GSE25724",
        "series_group": "GSE25nnn",
        "tissue": "pancreatic_islet",
        "platform": "GPL96",
        "platform_annotation": "GPL96",
        "case_patterns": [r"\btype 2\s+diabetic\b|\bhuman islets,\s*diabetic\b"],
        "control_patterns": [r"\bnon[- ]diabetic\b|\bnondiabetic\b"],
        "exclusion_patterns": [],
        "source_note": "isolated human pancreatic islets; T2D versus non-diabetic",
    },
]


def url_for_dataset(accession: str, group: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{group}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def url_for_platform(platform: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/{platform}/annot/{platform}.annot.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "whynot17-step08d/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as target:
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            target.write(block)


def parse_quoted_tsv(line: str) -> list[str]:
    return [str(value).strip() for value in next(csv.reader([line.rstrip("\n")], delimiter="\t"))]


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]], str, list[str]]:
    metadata: dict[str, list[str]] = {}
    sample_ids: list[str] = []
    platform = ""
    begin = None
    end = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for index, line in enumerate(handle):
            if line.startswith("!Series_platform_id"):
                values = parse_quoted_tsv(line)
                platform = values[1] if len(values) > 1 else ""
            elif line.startswith("!Sample_geo_accession"):
                values = parse_quoted_tsv(line)
                sample_ids = values[1:]
            elif line.startswith("!Sample_"):
                values = parse_quoted_tsv(line)
                if len(values) > 1:
                    metadata.setdefault(values[0], []).append(values[1:])
            elif line.startswith("!series_matrix_table_begin"):
                begin = index
            elif line.startswith("!series_matrix_table_end"):
                end = index
                break
    if begin is None or end is None or not sample_ids:
        raise ValueError(f"Could not identify series matrix table or sample IDs in {path}")
    nrows = end - begin - 2
    matrix = pd.read_csv(path, sep="\t", compression="gzip", skiprows=begin + 1, nrows=nrows)
    matrix = matrix.rename(columns={matrix.columns[0]: "probe_id"})
    matrix["probe_id"] = matrix["probe_id"].astype(str).str.strip().str.strip('"')
    matrix = matrix.set_index("probe_id")
    matrix.columns = [str(column).strip().strip('"') for column in matrix.columns]
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if list(matrix.columns) != sample_ids:
        matrix = matrix.reindex(columns=sample_ids)

    sample_meta: dict[str, list[str]] = {sample: [] for sample in sample_ids}
    for key, records in metadata.items():
        for record in records:
            for sample, value in zip(sample_ids, record):
                if key.startswith("!Sample_characteristics_ch1") or key in {"!Sample_title", "!Sample_source_name_ch1", "!Sample_description"}:
                    sample_meta[sample].append(f"{key}: {value.strip().strip(chr(34))}")
    return matrix, sample_meta, platform, sample_ids


def parse_platform_annotation(path: Path) -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    in_table = False
    header: list[str] | None = None
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for raw in reader:
            if not raw:
                continue
            first = raw[0].strip().strip('"')
            if first in {"!platform_table_begin", "!Platform_table_begin"}:
                in_table = True
                header = None
                continue
            if first in {"!platform_table_end", "!Platform_table_end"}:
                break
            if not in_table:
                continue
            if header is None:
                header = [cell.strip().strip('"') for cell in raw]
                continue
            if len(raw) < 1:
                continue
            values = [cell.strip().strip('"') for cell in raw]
            if len(values) < len(header):
                values.extend([""] * (len(header) - len(values)))
            id_value = values[0]
            if not id_value:
                continue
            symbol_index = next((i for i, name in enumerate(header) if name.lower() in {"gene symbol", "gene_symbol", "gene symbol (entrez)"}), None)
            if symbol_index is None:
                symbol_index = next((i for i, name in enumerate(header) if "gene symbol" in name.lower()), None)
            if symbol_index is None:
                continue
            symbols = [item.strip() for item in re.split(r"[;,| ]+", values[symbol_index]) if item.strip() and item.strip() not in {"---", "-", "NA", "N/A"}]
            rows[id_value] = symbols
    if not rows:
        raise ValueError(f"No probe-to-symbol rows found in {path}")
    return rows


def classify_samples(sample_meta: dict[str, list[str]], config: dict[str, object]) -> pd.DataFrame:
    records = []
    case_patterns = [re.compile(pattern, re.I) for pattern in config["case_patterns"]]
    control_patterns = [re.compile(pattern, re.I) for pattern in config["control_patterns"]]
    exclusion_patterns = [re.compile(pattern, re.I) for pattern in config["exclusion_patterns"]]
    for sample, values in sample_meta.items():
        text = " | ".join(values)
        excluded = any(pattern.search(text) for pattern in exclusion_patterns)
        is_case = any(pattern.search(text) for pattern in case_patterns)
        is_control = any(pattern.search(text) for pattern in control_patterns)
        if excluded:
            group = "excluded"
        elif is_case and not is_control:
            group = "T2D"
        elif is_control and not is_case:
            group = "control"
        elif is_case and is_control:
            group = "ambiguous"
        else:
            group = "unclassified"
        records.append({"sample_id": sample, "group": group, "metadata_text": text})
    return pd.DataFrame(records).sort_values("sample_id").reset_index(drop=True)


def load_frozen_modules() -> dict[str, list[str]]:
    modules = pd.read_csv(NETWORK_MODULES)
    nodes = pd.read_csv(NETWORK_NODES)
    name_by_node = dict(zip(nodes["node"].astype(str), nodes["preferred_name"].astype(str)))
    result: dict[str, list[str]] = {}
    for row in modules.itertuples(index=False):
        node_ids = [node for node in str(row.nodes).split(";") if node]
        genes = sorted({name_by_node.get(node, "") for node in node_ids if name_by_node.get(node, "") and name_by_node.get(node, "") != "nan"})
        result[str(row.module_id)] = genes
    return result


def collapse_to_genes(matrix: pd.DataFrame, probe_to_symbols: dict[str, list[str]]) -> tuple[pd.DataFrame, dict[str, int]]:
    values: dict[str, list[pd.Series]] = {}
    probe_counts: dict[str, int] = {}
    for probe, row in matrix.iterrows():
        for symbol in probe_to_symbols.get(str(probe), []):
            values.setdefault(symbol, []).append(row)
    collapsed = {}
    for symbol, rows in values.items():
        collapsed[symbol] = pd.concat(rows, axis=1).median(axis=1)
        probe_counts[symbol] = len(rows)
    return pd.DataFrame(collapsed).T, probe_counts


def bh_adjust(pvalues: list[float]) -> list[float]:
    indexed = sorted((value, index) for index, value in enumerate(pvalues) if math.isfinite(value))
    output = [math.nan] * len(pvalues)
    running = 1.0
    n = len(indexed)
    for rank, (value, index) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, value * n / rank)
        output[index] = min(running, 1.0)
    return output


def module_scores(gene_matrix: pd.DataFrame, modules: dict[str, list[str]]) -> tuple[pd.DataFrame, dict[str, int]]:
    # Platform annotations can contain repeated symbol rows after probe
    # collapsing.  Enforce one expression row per symbol before z-scoring so
    # every module gene is represented by a one-dimensional sample vector.
    if gene_matrix.index.has_duplicates:
        gene_matrix = gene_matrix.groupby(level=0, sort=True).median()
    zscores = gene_matrix.copy()
    for gene in zscores.index:
        values = zscores.loc[gene].astype(float)
        mean = values.mean()
        std = values.std(ddof=1)
        zscores.loc[gene] = 0.0 if not math.isfinite(std) or std == 0 else (values - mean) / std
    scores = {}
    observed_counts: dict[str, int] = {}
    for module_id, genes in modules.items():
        observed = sorted(set(genes) & set(zscores.index))
        observed_counts[module_id] = len(observed)
        if len(observed) >= 3:
            scores[module_id] = zscores.loc[observed].mean(axis=0)
        else:
            scores[module_id] = pd.Series(np.nan, index=zscores.columns)
    return pd.DataFrame(scores).T, observed_counts


def contrast(score: pd.Series, groups: pd.Series) -> dict[str, float | int | str]:
    case = pd.to_numeric(score[groups == "T2D"], errors="coerce").dropna().to_numpy(dtype=float)
    control = pd.to_numeric(score[groups == "control"], errors="coerce").dropna().to_numpy(dtype=float)
    if len(case) < 2 or len(control) < 2:
        return {"n_case": len(case), "n_control": len(control), "delta_t2d_minus_control": math.nan, "ci_low": math.nan, "ci_high": math.nan, "p_value": math.nan, "status": "insufficient_group_size", "direction": "not_estimable"}
    result = stats.ttest_ind(case, control, equal_var=False, nan_policy="omit")
    delta = float(np.mean(case) - np.mean(control))
    se = math.sqrt(float(np.var(case, ddof=1) / len(case) + np.var(control, ddof=1) / len(control)))
    df_num = (np.var(case, ddof=1) / len(case) + np.var(control, ddof=1) / len(control)) ** 2
    df_den = (np.var(case, ddof=1) / len(case)) ** 2 / (len(case) - 1) + (np.var(control, ddof=1) / len(control)) ** 2 / (len(control) - 1)
    df = df_num / df_den if df_den > 0 else math.nan
    margin = float(stats.t.ppf(0.975, df) * se) if math.isfinite(df) else math.nan
    return {"n_case": len(case), "n_control": len(control), "delta_t2d_minus_control": delta, "ci_low": delta - margin, "ci_high": delta + margin, "p_value": float(result.pvalue), "status": "ok", "direction": "up_in_T2D" if delta > 0 else "down_in_T2D" if delta < 0 else "flat"}


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    modules = load_frozen_modules()
    dataset_rows = []
    result_rows = []
    score_frames = []
    manifest = {
        "status": "running",
        "analysis": "Step 8D T2D transcriptomic module directionality",
        "module_input": str(NETWORK_MODULES),
        "n_frozen_modules": len(modules),
        "datasets": [],
        "module_score_rule": "mean within-dataset gene-wise z-scores; minimum 3 observed genes",
        "contrast": "T2D minus comparator; Welch two-sample t test per dataset and module",
        "multiple_testing": "BH-FDR within each dataset across tested modules",
        "interpretation": "descriptive module directionality; no tissue pooling, exposure causality, activation, mediation, or flagship selection",
    }
    for config in DATASETS:
        accession = str(config["accession"])
        print(f"[Step8D] {accession}: download/parse", flush=True)
        matrix_path = RAW / f"{accession}_series_matrix.txt.gz"
        matrix_url = url_for_dataset(accession, str(config["series_group"]))
        download(matrix_url, matrix_path)
        matrix, sample_meta, observed_platform, sample_ids = parse_series_matrix(matrix_path)
        sample_table = classify_samples(sample_meta, config)
        eligible = sample_table[sample_table["group"].isin(["T2D", "control"])].copy()
        eligible_ids = eligible["sample_id"].tolist()
        if len(eligible_ids) == 0 or set(eligible["group"]) != {"T2D", "control"}:
            raise ValueError(f"{accession} did not yield both T2D and control samples: {sample_table['group'].value_counts().to_dict()}")
        matrix = matrix.reindex(columns=eligible_ids)
        annotation_platform = str(config["platform_annotation"])
        annotation_path = RAW / f"{annotation_platform}.annot.gz"
        annotation_url = url_for_platform(annotation_platform)
        download(annotation_url, annotation_path)
        probe_to_symbols = parse_platform_annotation(annotation_path)
        gene_matrix, probe_counts = collapse_to_genes(matrix, probe_to_symbols)
        gene_matrix = gene_matrix.dropna(axis=0, how="all")
        scores, observed_counts = module_scores(gene_matrix, modules)
        print(f"[Step8D] {accession}: {len(eligible_ids)} eligible samples, {gene_matrix.shape[0]} genes, {sum(value >= 3 for value in observed_counts.values())}/{len(modules)} modules tested", flush=True)
        groups = eligible.set_index("sample_id").loc[eligible_ids, "group"]
        dataset_results = []
        pvalues = []
        for module_id in sorted(modules):
            row = {"accession": accession, "tissue": config["tissue"], "module_id": module_id, "n_module_genes_frozen": len(modules[module_id]), "n_module_genes_expression": observed_counts[module_id], "platform": observed_platform, "sample_selection": config["source_note"]}
            row.update(contrast(scores.loc[module_id], groups))
            dataset_results.append(row)
            pvalues.append(float(row["p_value"]) if math.isfinite(float(row["p_value"])) else math.nan)
        qvalues = bh_adjust(pvalues)
        for row, qvalue in zip(dataset_results, qvalues):
            row["q_value_within_dataset"] = qvalue
            result_rows.append(row)
        scores = scores.add_prefix(f"{accession}|")
        score_frames.append(scores)
        dataset_rows.append({
            "accession": accession,
            "tissue": config["tissue"],
            "platform_declared": config["platform"],
            "platform_observed": observed_platform,
            "platform_annotation_used": annotation_platform,
            "n_matrix_probes": matrix.shape[0],
            "n_expression_genes": gene_matrix.shape[0],
            "n_total_samples": len(sample_table),
            "n_t2d": int((eligible["group"] == "T2D").sum()),
            "n_control": int((eligible["group"] == "control").sum()),
            "n_excluded_or_unclassified": int(len(sample_table) - len(eligible)),
            "n_modules": len(modules),
            "n_modules_tested": int(sum(value >= 3 for value in observed_counts.values())),
            "series_matrix_sha256": sha256(matrix_path),
            "platform_annotation_sha256": sha256(annotation_path),
            "source_note": config["source_note"],
        })
        manifest["datasets"].append({
            "accession": accession,
            "tissue": config["tissue"],
            "series_matrix_url": matrix_url,
            "platform_annotation_url": annotation_url,
            "platform_declared": config["platform"],
            "platform_observed": observed_platform,
            "platform_annotation_used": annotation_platform,
            "sample_rule": {"case_patterns": config["case_patterns"], "control_patterns": config["control_patterns"], "exclusion_patterns": config["exclusion_patterns"]},
            "n_samples_total": len(sample_table),
            "n_t2d": int((eligible["group"] == "T2D").sum()),
            "n_control": int((eligible["group"] == "control").sum()),
            "series_matrix_sha256": sha256(matrix_path),
            "platform_annotation_sha256": sha256(annotation_path),
        })
        sample_table.to_csv(OUT / f"t2d_step8d_{accession.lower()}_sample_audit.csv", index=False)
    results = pd.DataFrame(result_rows)
    results = results.sort_values(["accession", "module_id"]).reset_index(drop=True)
    results.to_csv(OUT / "t2d_step8d_module_directionality.csv", index=False)
    pd.DataFrame(dataset_rows).to_csv(OUT / "t2d_step8d_dataset_audit.csv", index=False)
    score_matrix = pd.concat(score_frames, axis=1).sort_index()
    score_matrix.to_csv(OUT / "t2d_step8d_module_scores.csv")

    synthesis_rows = []
    for module_id, frame in results.groupby("module_id", sort=True):
        tested = frame[frame["status"] == "ok"]
        deltas = pd.to_numeric(tested["delta_t2d_minus_control"], errors="coerce")
        synthesis_rows.append({
            "module_id": module_id,
            "n_datasets_total": len(frame),
            "n_datasets_tested": len(tested),
            "n_positive": int((deltas > 0).sum()),
            "n_negative": int((deltas < 0).sum()),
            "n_flat": int((deltas == 0).sum()),
            "direction_concordance_fraction": float(max((deltas > 0).sum(), (deltas < 0).sum()) / len(tested)) if len(tested) else math.nan,
            "median_delta": float(deltas.median()) if len(tested) else math.nan,
            "n_dataset_q_lt_0_05": int((pd.to_numeric(tested["q_value_within_dataset"], errors="coerce") < 0.05).sum()),
            "datasets_positive": ";".join(tested.loc[deltas > 0, "accession"].tolist()),
            "datasets_negative": ";".join(tested.loc[deltas < 0, "accession"].tolist()),
        })
    synthesis = pd.DataFrame(synthesis_rows)
    synthesis.to_csv(OUT / "t2d_step8d_module_cross_dataset_synthesis.csv", index=False)

    axis_map = pd.read_csv(NETWORK_MODULES)[["module_id", "cluster_id"]]
    axis_summary = synthesis.merge(axis_map, on="module_id", how="left").groupby("cluster_id", as_index=False).agg(
        n_modules=("module_id", "count"),
        n_modules_tested=("n_datasets_tested", lambda x: int((x > 0).sum())),
        median_module_delta=("median_delta", "median"),
        n_modules_majority_positive=("n_positive", lambda x: int((x > 0).sum())),
        n_modules_majority_negative=("n_negative", lambda x: int((x > 0).sum())),
        n_modules_with_any_q_lt_0_05=("n_dataset_q_lt_0_05", lambda x: int((x > 0).sum())),
    )
    axis_summary.to_csv(OUT / "t2d_step8d_axis_summary.csv", index=False)

    dataset_axis_rows = []
    results_with_axis = results.merge(axis_map, on="module_id", how="left")
    for (accession, tissue, cluster_id), frame in results_with_axis.groupby(["accession", "tissue", "cluster_id"], sort=True):
        tested = frame[frame["status"] == "ok"]
        deltas = pd.to_numeric(tested["delta_t2d_minus_control"], errors="coerce")
        n_positive = int((deltas > 0).sum())
        n_negative = int((deltas < 0).sum())
        dataset_axis_rows.append({
            "accession": accession,
            "tissue": tissue,
            "cluster_id": cluster_id,
            "n_modules_total": len(frame),
            "n_modules_tested": len(tested),
            "n_modules_positive": n_positive,
            "n_modules_negative": n_negative,
            "direction_concordance_fraction": float(max(n_positive, n_negative) / len(tested)) if len(tested) else math.nan,
            "median_module_delta": float(deltas.median()) if len(tested) else math.nan,
            "n_modules_q_lt_0_05": int((pd.to_numeric(tested["q_value_within_dataset"], errors="coerce") < 0.05).sum()),
        })
    dataset_axis_summary = pd.DataFrame(dataset_axis_rows)
    dataset_axis_summary.to_csv(OUT / "t2d_step8d_dataset_axis_summary.csv", index=False)

    manifest["status"] = "complete_module_directionality"
    manifest["retrieved_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["outputs"] = [
        "t2d_step8d_dataset_audit.csv",
        "t2d_step8d_module_directionality.csv",
        "t2d_step8d_module_scores.csv",
        "t2d_step8d_module_cross_dataset_synthesis.csv",
        "t2d_step8d_axis_summary.csv",
        "t2d_step8d_dataset_axis_summary.csv",
        "t2d_step8d_*_sample_audit.csv",
        "STEP8D_MANIFEST.json",
        "STEP8D_T2D_TRANSCRIPTOMIC_REPORT.md",
    ]
    (OUT / "STEP8D_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    report = [
        "# Step 8D — T2D transcriptomic module directionality",
        "",
        "- Status: **complete_module_directionality**",
        f"- Frozen Step 8C modules tested: **{len(modules)}**",
        "- Primary unit: biological sample within each GEO series; no tissue pooling",
        "- Module score: mean within-dataset gene-wise z-scores; minimum 3 mapped genes",
        "- Test: T2D minus comparator Welch t test; BH-FDR within dataset",
        "",
        "## Dataset audit",
        "",
        "| Accession | Tissue | T2D | Control | Expression genes | Modules tested |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in dataset_rows:
        report.append(f"| {row['accession']} | {row['tissue']} | {row['n_t2d']} | {row['n_control']} | {row['n_expression_genes']} | {row['n_modules_tested']} |")
    report.extend([
        "",
        "## Tissue-stratified axis summary",
        "",
        "This table audits module signs within each dataset and Tier A axis. It is not a pooled test and is not used to select a flagship axis.",
        "",
        "| Accession | Tissue | Axis | Tested modules | Positive | Negative | Median delta | q<0.05 |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ])
    for row in dataset_axis_summary.itertuples(index=False):
        median_delta = f"{row.median_module_delta:.3f}" if math.isfinite(row.median_module_delta) else "NA"
        report.append(f"| {row.accession} | {row.tissue} | {row.cluster_id} | {row.n_modules_tested} | {row.n_modules_positive} | {row.n_modules_negative} | {median_delta} | {row.n_modules_q_lt_0_05} |")
    report.extend([
        "",
        "## Directionality boundary",
        "",
        "Cross-dataset synthesis is descriptive and retains tissue context. A positive or negative module score indicates a relative transcriptomic shift in the public series; it does not establish pathway activation, exposure causality, mediation, or a T2D-specific mechanism. No flagship axis is selected in Step 8D.",
        "",
        "## Reproducibility",
        "",
        "The series matrix and platform annotation URLs, sample-label rules, and SHA-256 checksums are recorded in `STEP8D_MANIFEST.json`. Raw GEO files are local-only and excluded from version control.",
    ])
    (OUT / "STEP8D_T2D_TRANSCRIPTOMIC_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
