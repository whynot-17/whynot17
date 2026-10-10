from __future__ import annotations

import csv
import gzip
import hashlib
import json
import pathlib
import re
import time
from collections import Counter, defaultdict
from datetime import date, datetime, timezone

import numpy as np
import pandas as pd
import requests
from scipy.stats import t as student_t
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
CACHE = ROOT / "work" / "geo_cache_replacement_2026-10-10"
OUT = ROOT / "results" / "2026-10-10" / "geo_replication_113_replacement"
INPUT_113 = ROOT / "results" / "2026-10-10" / "DINP_286_x_TCGA_GDC_FDR_7748.tsv"
ACCESSIONS = ("GSE32323", "GSE44076", "GSE41258", "GSE106582")
MATRIX_URLS = {
    "GSE32323": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32323/matrix/GSE32323_series_matrix.txt.gz",
    "GSE44076": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44076/matrix/GSE44076_series_matrix.txt.gz",
    "GSE41258": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41258/matrix/GSE41258_series_matrix.txt.gz",
    "GSE106582": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106582/matrix/GSE106582_series_matrix.txt.gz",
}
PLATFORM_BY_GSE = {
    "GSE32323": "GPL570",
    "GSE44076": "GPL13667",
    "GSE41258": "GPL96",
    "GSE106582": "GPL10558",
}
PLATFORM_URLS = {
    "GPL570": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "GPL96": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
    "GPL10558": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
    "GPL13667": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13667&targ=self&view=data&form=text",
}
STABLE_FDR = 0.05
STABLE_ABS_LOG2FC = 1.0
USER_AGENT = "DINP-CRC-GEO-validation/2.0"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size:
        return
    temp_path = path.with_suffix(path.suffix + ".part")
    last_error = None
    for attempt in range(1, 4):
        try:
            with requests.get(
                url, stream=True, timeout=(30, 600), headers={"User-Agent": USER_AGENT}
            ) as response:
                response.raise_for_status()
                with temp_path.open("wb") as stream:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            stream.write(chunk)
            if not temp_path.stat().st_size:
                raise IOError(f"Downloaded empty file from {url}")
            temp_path.replace(path)
            return
        except Exception as exc:
            last_error = exc
            temp_path.unlink(missing_ok=True)
            if attempt < 3:
                time.sleep(attempt * 2)
    raise RuntimeError(f"Failed to download {url}") from last_error


def parse_geo_values(text: str) -> list[str]:
    return next(csv.reader([text], delimiter="\t", quotechar='"', escapechar="\\"))


def load_series_matrix(accession: str) -> tuple[dict, pd.DataFrame, str]:
    path = CACHE / f"{accession}_series_matrix.txt.gz"
    download(MATRIX_URLS[accession], path)
    metadata: dict[str, list[list[str]]] = defaultdict(list)
    start = end = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream):
            line = line.rstrip("\r\n")
            if line.lower() == "!series_matrix_table_begin":
                start = line_number
                break
            if line.startswith("!Sample_") and "\t" in line:
                key, values = line.split("\t", 1)
                metadata[key].append(parse_geo_values(values))
    if start is None:
        raise ValueError(f"{accession}: missing series matrix table")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream):
            if line.strip().lower() == "!series_matrix_table_end":
                end = line_number
                break
    if end is None:
        raise ValueError(f"{accession}: missing series matrix table end")
    table = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        skiprows=start + 1,
        nrows=end - start - 2,
        quotechar='"',
        dtype={"ID_REF": str},
    )
    if table.columns[0] != "ID_REF":
        raise ValueError(f"{accession}: expected ID_REF, found {table.columns[0]}")
    return metadata, table, sha256(path)


def sample_records(metadata: dict) -> list[dict]:
    sample_ids = metadata["!Sample_geo_accession"][0]
    n = len(sample_ids)
    titles = metadata.get("!Sample_title", [[""] * n])[0]
    platforms = metadata.get("!Sample_platform_id", [[""] * n])[0]
    sources = metadata.get("!Sample_source_name_ch1", [[""] * n])[0]
    characteristics = [dict() for _ in range(n)]
    for key, rows in metadata.items():
        if key.startswith("!Sample_characteristics_ch1"):
            for row in rows:
                for index, value in enumerate(row):
                    if ":" in value:
                        char_key, char_value = value.split(":", 1)
                        characteristics[index][char_key.strip().lower()] = char_value.strip()
    return [
        {
            "gsm": sample_ids[i],
            "title": titles[i] if i < len(titles) else "",
            "platform": platforms[i] if i < len(platforms) else "",
            "source_name": sources[i] if i < len(sources) else "",
            **characteristics[i],
        }
        for i in range(n)
    ]


def assign_groups(accession: str, records: list[dict]) -> list[dict]:
    assigned = []
    for record in records:
        tissue = record.get("tissue", "").strip().lower()
        if accession == "GSE32323":
            group = {
                "cancer, homogenized": "Tumor",
                "normal, homogenized": "Normal",
            }.get(tissue, "Excluded")
            match = re.search(r"patient\s+([^, ]+)", record.get("title", ""), re.I)
            pair_id = match.group(1).lower().rstrip("b") if match else ""
        elif accession == "GSE44076":
            sample_type = record.get("sample type", "").strip().lower()
            group = {"tumor": "Tumor", "normal": "Normal"}.get(sample_type, "Excluded")
            pair_id = record.get("individual id", "").strip()
        elif accession == "GSE41258":
            included = record.get("included in analysis", "").strip().lower() == "yes"
            group = (
                "Tumor" if included and tissue == "primary tumor"
                else "Normal" if included and tissue == "normal colon"
                else "Excluded"
            )
            pair_id = record.get("patient id", "").strip()
        else:
            group = {"tumor": "Tumor", "mucosa": "Normal"}.get(tissue, "Excluded")
            pair_id = record.get("patientid", "").strip()
        assigned.append({**record, "group": group, "pair_id": pair_id})
    return assigned


def split_symbols(value: str) -> list[str]:
    if not value or value.strip() in {"---", "NA", "N/A"}:
        return []
    symbols = []
    for part in re.split(r"\s*///\s*", value.strip()):
        symbol = part.strip().upper()
        if symbol and symbol not in {"---", "NA", "N/A"} and symbol not in symbols:
            symbols.append(symbol)
    return symbols


def parse_annot_gz(path: pathlib.Path) -> tuple[dict[str, list[str]], str]:
    mapping: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        for row in reader:
            if row and row[0].strip().lower() == "id" and any(
                cell.strip().lower() == "gene symbol" for cell in row
            ):
                header = row
                id_col = next(i for i, cell in enumerate(row) if cell.strip().lower() == "id")
                symbol_col = next(
                    i for i, cell in enumerate(row) if cell.strip().lower() == "gene symbol"
                )
                break
        else:
            raise ValueError(f"GPL annotation header not found: {path}")
        for row in reader:
            if len(row) <= max(id_col, symbol_col):
                continue
            probe = row[id_col].strip()
            symbols = split_symbols(row[symbol_col])
            if probe and symbols:
                mapping[probe] = symbols
    if not mapping:
        raise ValueError(f"Empty probe mapping: {path}")
    return mapping, sha256(path)


def parse_gpl13667(path: pathlib.Path) -> tuple[dict[str, list[str]], str]:
    mapping: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t", quotechar='"')
        for row in reader:
            if row and row[0].strip().lower() == "!platform_table_begin":
                header = next(reader)
                id_col = next(i for i, cell in enumerate(header) if cell.strip().lower() == "id")
                symbol_col = next(
                    i for i, cell in enumerate(header) if cell.strip().lower() == "gene symbol"
                )
                break
        else:
            raise ValueError("GPL13667 platform table not found")
        for row in reader:
            if not row or row[0].strip().lower() == "!platform_table_end":
                break
            if len(row) <= max(id_col, symbol_col):
                continue
            probe = row[id_col].strip()
            symbols = split_symbols(row[symbol_col])
            if probe and symbols:
                mapping[probe] = symbols
    if len(mapping) < 10000:
        raise ValueError(f"GPL13667 mapping incomplete: {len(mapping)} probes")
    return mapping, sha256(path)


def platform_mapping(platform: str) -> tuple[dict[str, list[str]], str]:
    if platform == "GPL13667":
        path = CACHE / "GPL13667_platform_full.txt"
        download(PLATFORM_URLS[platform], path)
        return parse_gpl13667(path)
    path = CACHE / f"{platform}.annot.gz"
    download(PLATFORM_URLS[platform], path)
    return parse_annot_gz(path)


def map_probe_matrix(
    table: pd.DataFrame, mapping: dict[str, list[str]], accession: str
) -> tuple[list[str], np.ndarray, dict[str, int], bool]:
    probes = table.iloc[:, 0].astype(str).tolist()
    raw = table.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    transformed = accession == "GSE41258"
    if transformed:
        values = raw[np.isfinite(raw)]
        if values.size and np.nanmin(values) <= 0:
            raw = np.log2(np.maximum(raw, 0) + 1.0)
        else:
            raw = np.log2(raw)
    probe_indices: dict[str, list[int]] = defaultdict(list)
    for index, probe in enumerate(probes):
        for symbol in mapping.get(probe, []):
            probe_indices[symbol].append(index)
    genes = sorted(probe_indices)
    matrix = np.vstack([
        np.nanmedian(raw[probe_indices[symbol], :], axis=0)
        for symbol in genes
    ])
    probe_counts = {symbol: len(probe_indices[symbol]) for symbol in genes}
    return genes, matrix, probe_counts, transformed


def adjust_bh(pvalues: np.ndarray) -> np.ndarray:
    adjusted = np.full(pvalues.shape, np.nan, dtype=float)
    valid = np.isfinite(pvalues)
    if valid.any():
        adjusted[valid] = multipletests(pvalues[valid], method="fdr_bh")[1]
    return adjusted


def complete_pair_ids(records: list[dict], indices: list[int]) -> list[str]:
    groups: dict[str, set[str]] = defaultdict(set)
    for index in indices:
        record = records[index]
        if record["pair_id"]:
            groups[record["pair_id"]].add(record["group"])
    return sorted(patient for patient, labels in groups.items() if {"Tumor", "Normal"} <= labels)


def paired_test(
    matrix: np.ndarray, records: list[dict], indices: list[int]
) -> tuple[np.ndarray, np.ndarray, int, list[str]]:
    grouped: dict[str, dict[str, list[int]]] = defaultdict(lambda: {"Tumor": [], "Normal": []})
    for index in indices:
        record = records[index]
        if record["pair_id"]:
            grouped[record["pair_id"]][record["group"]].append(index)
    pairs = []
    pair_ids = []
    for patient, groups in sorted(grouped.items()):
        if groups["Tumor"] and groups["Normal"]:
            pairs.append((
                np.nanmean(matrix[:, groups["Tumor"]], axis=1),
                np.nanmean(matrix[:, groups["Normal"]], axis=1),
            ))
            pair_ids.append(patient)
    if not pairs:
        raise ValueError("No complete tumor-normal pairs")
    tumors = np.column_stack([pair[0] for pair in pairs])
    normals = np.column_stack([pair[1] for pair in pairs])
    log2fc = np.nanmean(tumors - normals, axis=1)
    pvalues = ttest_rel(tumors, normals, axis=1, nan_policy="omit").pvalue
    return log2fc, np.asarray(pvalues, dtype=float), len(pair_ids), pair_ids


def patient_cluster_test(
    matrix: np.ndarray, records: list[dict], indices: list[int]
) -> tuple[np.ndarray, np.ndarray, int, list[str]]:
    group_is_tumor = np.asarray([records[i]["group"] == "Tumor" for i in indices], dtype=bool)
    values = matrix[:, indices]
    beta = np.nanmean(values[:, group_is_tumor], axis=1) - np.nanmean(values[:, ~group_is_tumor], axis=1)
    residual = values.copy()
    residual[:, group_is_tumor] -= np.nanmean(values[:, group_is_tumor], axis=1)[:, None]
    residual[:, ~group_is_tumor] -= np.nanmean(values[:, ~group_is_tumor], axis=1)[:, None]
    records_local = [records[i] for i in indices]
    clusters: dict[str, list[int]] = defaultdict(list)
    for local_index, record in enumerate(records_local):
        clusters[record["pair_id"].strip() or f"unmatched_{record['gsm']}"].append(local_index)
    n_clusters = len(clusters)
    if n_clusters < 3:
        raise ValueError("Too few patient clusters for inference")
    meat_diag = np.zeros(matrix.shape[0], dtype=float)
    for members in clusters.values():
        ix = np.asarray(members, dtype=int)
        tumor_ix = ix[group_is_tumor[ix]]
        normal_ix = ix[~group_is_tumor[ix]]
        sum_tumor = residual[:, tumor_ix].sum(axis=1) if len(tumor_ix) else 0.0
        sum_normal = residual[:, normal_ix].sum(axis=1) if len(normal_ix) else 0.0
        influence = sum_tumor / group_is_tumor.sum() - sum_normal / (~group_is_tumor).sum()
        meat_diag += influence * influence
    n_samples = len(indices)
    correction = (n_clusters / (n_clusters - 1)) * ((n_samples - 1) / (n_samples - 2))
    se = np.sqrt(meat_diag * correction)
    statistic = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    pvalues = 2.0 * student_t.sf(np.abs(statistic), df=n_clusters - 1)
    zero_se = (se == 0) & np.isfinite(beta)
    pvalues[zero_se & (beta == 0)] = 1.0
    pvalues[zero_se & (beta != 0)] = 0.0
    return beta, pvalues, len(complete_pair_ids(records, indices)), complete_pair_ids(records, indices)


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    input_frame = pd.read_csv(INPUT_113, sep="\t")
    targets = sorted(set(input_frame["gene_symbol"].astype(str).str.strip().str.upper()) - {""})
    if len(targets) != 113:
        raise ValueError(f"Expected 113 unique query symbols, found {len(targets)}")
    checksums = {"input_113_sha256": sha256(INPUT_113)}
    count_rows, coverage_rows, result_rows, audit_rows = [], [], [], []

    for accession in ACCESSIONS:
        metadata, table, matrix_hash = load_series_matrix(accession)
        records = assign_groups(accession, sample_records(metadata))
        sample_columns = list(table.columns[1:])
        records_by_gsm = {record["gsm"]: record for record in records}
        missing = set(sample_columns) - set(records_by_gsm)
        if missing:
            raise ValueError(f"{accession}: sample metadata missing for {sorted(missing)[:5]}")
        records = [records_by_gsm[gsm] for gsm in sample_columns]
        platform = PLATFORM_BY_GSE[accession]
        probe_map, annotation_hash = platform_mapping(platform)
        genes, matrix, probe_counts, transformed = map_probe_matrix(table, probe_map, accession)
        checksums[f"{accession}_series_matrix_sha256"] = matrix_hash
        checksums[f"{platform}_annotation_sha256"] = annotation_hash
        selected = [i for i, record in enumerate(records) if record["group"] in {"Tumor", "Normal"}]
        n_tumor = sum(records[i]["group"] == "Tumor" for i in selected)
        n_normal = sum(records[i]["group"] == "Normal" for i in selected)
        available_pairs = complete_pair_ids(records, selected)

        if accession == "GSE41258":
            method = "patient-cluster robust OLS mean difference"
            log2fc, pvalues, pairs_used, test_pair_ids = patient_cluster_test(matrix, records, selected)
            tested_samples = len(selected)
        else:
            method = "paired t-test on complete patient pairs"
            log2fc, pvalues, pairs_used, test_pair_ids = paired_test(matrix, records, selected)
            tested_samples = 2 * pairs_used
        fdr = adjust_bh(pvalues)
        gene_index = {symbol: i for i, symbol in enumerate(genes)}
        represented = set(gene_index).intersection(targets)
        absent = sorted(set(targets) - represented)
        coverage_rows.append({
            "accession": accession,
            "platform": platform,
            "target_genes_total": len(targets),
            "target_genes_with_probe": len(represented),
            "target_genes_without_probe": len(absent),
            "without_probe_symbols": ";".join(absent),
            "unique_mapped_genes_tested_for_FDR": len(genes),
            "mapped_probe_count": sum(probe_counts.values()),
            "series_matrix_feature_count": int(table.shape[0]),
            "expression_transform": "log2(PLIER signal)" if transformed else "GEO series-matrix values retained",
            "statistical_method": method,
        })
        count_rows.append({
            "accession": accession,
            "platform": platform,
            "series_matrix_total_samples": len(records),
            "included_tumor_samples": n_tumor,
            "included_normal_samples": n_normal,
            "complete_patient_pairs_available": len(available_pairs),
            "complete_patient_pairs_used_for_test": pairs_used,
            "included_expression_samples_used_for_test": tested_samples,
            "excluded_samples": sum(record["group"] == "Excluded" for record in records),
            "statistical_method": method,
        })
        paired_ids = set(test_pair_ids)
        for record in records:
            audit_rows.append({
                "accession": accession,
                "gsm": record["gsm"],
                "sample_title": record["title"],
                "patient_id": record["pair_id"],
                "assigned_group": record["group"],
                "included_in_test": record["group"] in {"Tumor", "Normal"}
                and (accession == "GSE41258" or record["pair_id"] in paired_ids),
            })
        for gene in targets:
            if gene not in gene_index:
                result_rows.append({
                    "gene_symbol": gene, "accession": accession, "platform": platform,
                    "n_tumor": n_tumor, "n_normal": n_normal, "n_pairs_used": pairs_used,
                    "statistical_method": method, "probes_collapsed": 0,
                    "log2FC_tumor_minus_normal": np.nan, "P_value": np.nan,
                    "FDR_BH_all_mapped_genes": np.nan,
                    "passes_FDR_and_abs_log2FC": False, "direction": "not represented",
                })
                continue
            index = gene_index[gene]
            fc = float(log2fc[index])
            p = float(pvalues[index]) if np.isfinite(pvalues[index]) else np.nan
            q = float(fdr[index]) if np.isfinite(fdr[index]) else np.nan
            direction = "Tumor_up" if fc > 0 else "Tumor_down" if fc < 0 else "No_change"
            hit = bool(np.isfinite(q) and q < STABLE_FDR and abs(fc) > STABLE_ABS_LOG2FC)
            result_rows.append({
                "gene_symbol": gene,
                "accession": accession,
                "platform": platform,
                "n_tumor": n_tumor,
                "n_normal": n_normal,
                "n_pairs_used": pairs_used,
                "statistical_method": method,
                "probes_collapsed": probe_counts[gene],
                "log2FC_tumor_minus_normal": fc,
                "P_value": p,
                "FDR_BH_all_mapped_genes": q,
                "passes_FDR_and_abs_log2FC": hit,
                "direction": direction,
            })
        print(
            f"{accession}: samples={len(records)} tumor={n_tumor} normal={n_normal} "
            f"pairs={len(available_pairs)} used={pairs_used}; mapped={len(genes)} "
            f"target coverage={len(represented)}/113",
            flush=True,
        )

    per_cohort = pd.DataFrame(result_rows).sort_values(["gene_symbol", "accession"])
    counts = pd.DataFrame(count_rows)
    coverage = pd.DataFrame(coverage_rows)
    audit = pd.DataFrame(audit_rows)
    stability_rows = []
    for gene in targets:
        subset = per_cohort[per_cohort["gene_symbol"] == gene]
        hits = subset[subset["passes_FDR_and_abs_log2FC"]]
        directions = set(hits["direction"]) - {"not represented", "No_change"}
        consistent = len(directions) <= 1
        row = {
            "gene_symbol": gene,
            "significant_cohorts_FDR_lt_0.05_abs_log2FC_gt_1": len(hits),
            "cohorts_meeting_threshold": ";".join(hits["accession"].tolist()),
            "significant_up_cohorts": int((hits["direction"] == "Tumor_up").sum()),
            "significant_down_cohorts": int((hits["direction"] == "Tumor_down").sum()),
            "direction_consistent_among_significant": consistent,
            "replicated_in_at_least_3_of_4": len(hits) >= 3 and consistent,
            "replicated_in_all_4": len(hits) == 4 and consistent,
        }
        for accession in ACCESSIONS:
            record = subset[subset["accession"] == accession].iloc[0]
            row[f"{accession}_log2FC"] = record["log2FC_tumor_minus_normal"]
            row[f"{accession}_FDR"] = record["FDR_BH_all_mapped_genes"]
            row[f"{accession}_passes"] = bool(record["passes_FDR_and_abs_log2FC"])
        stability_rows.append(row)
    stability = pd.DataFrame(stability_rows)
    stable_3plus = stability[stability["replicated_in_at_least_3_of_4"]].copy()
    stable_4of4 = stability[stability["replicated_in_all_4"]].copy()

    per_cohort.to_csv(OUT / "113_gene_per_cohort_DE_results.tsv", sep="\t", index=False)
    counts.to_csv(OUT / "GEO_cohort_sample_counts.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "GEO_platform_gene_coverage.tsv", sep="\t", index=False)
    audit.to_csv(OUT / "GEO_sample_selection_audit.tsv", sep="\t", index=False)
    stability.to_csv(OUT / "113_gene_stability_summary.tsv", sep="\t", index=False)
    stable_3plus.to_csv(OUT / "stable_genes_replicated_in_3_or_more_of_4_cohorts.tsv", sep="\t", index=False)
    stable_4of4.to_csv(OUT / "stable_genes_replicated_in_all_4_cohorts.tsv", sep="\t", index=False)
    (OUT / "query_113_genes.txt").write_text("\n".join(targets) + "\n", encoding="utf-8")
    (OUT / "source_checksums.json").write_text(json.dumps(checksums, indent=2), encoding="utf-8")
    metadata = {
        "analysis_date": date.today().isoformat(),
        "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "organism": "Homo sapiens (NCBI taxon 9606)",
        "query_gene_count": len(targets),
        "cohorts": list(ACCESSIONS),
        "thresholds": {
            "BH_FDR": STABLE_FDR,
            "absolute_log2FC_strictly_greater_than": STABLE_ABS_LOG2FC,
            "three_of_four_rule": "each cohort must meet both thresholds; significant directions must agree",
            "four_of_four_rule": "all cohorts must meet both thresholds; directions must agree",
        },
        "multiple_testing_universe": "all unique genes mapped from each platform after median probe aggregation",
        "probe_collapsing": "median expression across probes assigned to each gene symbol",
        "series_matrix_urls": MATRIX_URLS,
        "platform_annotation_urls": PLATFORM_URLS,
        "checksums": checksums,
        "software": {"python": __import__("sys").version, "pandas": pd.__version__, "numpy": np.__version__},
    }
    (OUT / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )

    report = [
        "# 113 个 DINP×CRC 交集基因在替换后的四个 GEO 队列中的复现",
        "",
        f"- 分析日期：{date.today().isoformat()}",
        "- 查询集：本项目 113 个 DINP×TCGA-GDC FDR DEG 交集基因；表达矩阵及平台注释从 NCBI GEO 重新下载。",
        "- 队列：GSE32323、GSE44076、GSE41258、GSE106582。",
        "- 每队列复现阈值沿用此前口径：BH-FDR < 0.05 且 |log2FC| > 1；log2FC 为 Tumor−Normal。",
        "- 主复现集合要求至少 3/4 队列达到阈值，且显著队列方向一致；严格交集另按 4/4 列出。",
        "- 多探针以 gene symbol 聚合取表达中位数；BH-FDR 在该平台全部可映射基因上计算。",
        "- GSE32323、GSE44076、GSE106582 使用完整配对的 paired t-test。GSE41258 有患者重复切片，使用患者聚类稳健标准误；44 个完整配对患者和其他符合纳入条件的样本均进入模型。",
        "",
        "## 队列样本与配对数",
        "",
        "| 队列 | 平台 | 肿瘤样本 | 正常/黏膜样本 | 完整配对 | 实际使用配对 |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in count_rows:
        report.append(
            f"| {row['accession']} | {row['platform']} | {row['included_tumor_samples']} | "
            f"{row['included_normal_samples']} | {row['complete_patient_pairs_available']} | "
            f"{row['complete_patient_pairs_used_for_test']} |"
        )
    report.extend([
        "",
        f"- 至少 3/4 队列同向复现：**{len(stable_3plus)} 个基因**。",
        f"- 严格 4/4 同向复现：**{len(stable_4of4)} 个基因**。",
        "",
        "### 至少 3/4 队列复现",
        "",
    ])
    if stable_3plus.empty:
        report.append("无。")
    else:
        report.extend(["| 基因 | 队列数 | 达标队列 | 方向 |", "|---|---:|---|---|"])
        for _, row in stable_3plus.iterrows():
            direction = "Tumor_up" if row["significant_up_cohorts"] else "Tumor_down"
            report.append(
                f"| {row['gene_symbol']} | {row['significant_cohorts_FDR_lt_0.05_abs_log2FC_gt_1']} | "
                f"{row['cohorts_meeting_threshold']} | {direction} |"
            )
    report.extend(["", "### 严格 4/4 交集", ""])
    if stable_4of4.empty:
        report.append("无基因在四个队列中均达到 FDR/log2FC 阈值并保持同向。")
    else:
        report.extend(["| 基因 | 方向 |", "|---|---|"])
        for _, row in stable_4of4.iterrows():
            direction = "Tumor_up" if row["significant_up_cohorts"] else "Tumor_down"
            report.append(f"| {row['gene_symbol']} | {direction} |")
    report.extend([
        "",
        "## 数据来源与结果文件",
        "",
        "- [GSE32323](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32323)",
        "- [GSE44076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44076)",
        "- [GSE41258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41258)",
        "- [GSE106582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106582)",
        "- 113_gene_per_cohort_DE_results.tsv：各基因在队列中的效应、P 值、BH-FDR 和方向。",
        "- 113_gene_stability_summary.tsv：逐基因复现情况汇总。",
        "- stable_genes_replicated_in_3_or_more_of_4_cohorts.tsv 和 stable_genes_replicated_in_all_4_cohorts.tsv：3/4 及 4/4 交集。",
        "- GEO_cohort_sample_counts.tsv、GEO_sample_selection_audit.tsv：队列样本数、配对数和逐样本筛选记录。",
        "- GEO_platform_gene_coverage.tsv：候选基因在各平台的探针覆盖。",
        "- analysis_metadata.json、source_checksums.json：分析参数、来源 URL 和输入文件 SHA-256。",
        "- 可复现脚本：analysis/validate_113_across_geo_replacement.py。",
        "",
        "## 解释边界",
        "",
        "配对检验仅使用完整肿瘤/正常患者配对；GSE41258 的聚类稳健模型将患者作为聚类单位。不同平台的探针注释和信号尺度可能影响效应量；跨队列复现不代表因果关系或临床诊断性能。",
    ])
    (OUT / "README.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"query genes={len(targets)}")
    print(counts[[
        "accession", "included_tumor_samples", "included_normal_samples",
        "complete_patient_pairs_available", "complete_patient_pairs_used_for_test"
    ]].to_string(index=False))
    print(f"stable >=3/4={len(stable_3plus)}: {stable_3plus['gene_symbol'].tolist()}")
    print(f"strict 4/4={len(stable_4of4)}: {stable_4of4['gene_symbol'].tolist()}")
    print(f"Output: {OUT}")


if __name__ == "__main__":
    main()
