from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import pathlib
import re
import urllib.request
from collections import Counter, defaultdict
from datetime import date

import numpy as np
import pandas as pd
import requests
from scipy.stats import t as student_t
from scipy.stats import ttest_ind, ttest_rel
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
CACHE = ROOT / "work" / "geo_cache"
OUT = ROOT / "results" / date.today().isoformat() / "geo_replication_113"
INPUT_113 = ROOT / "results" / "2026-10-10" / "DINP_286_x_TCGA_GDC_FDR_7748.tsv"
ACCESSIONS = ("GSE21510", "GSE32323", "GSE44076", "GSE41258")
MATRIX_URLS = {
    "GSE21510": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21510/matrix/GSE21510_series_matrix.txt.gz",
    "GSE32323": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32323/matrix/GSE32323_series_matrix.txt.gz",
    "GSE44076": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44076/matrix/GSE44076_series_matrix.txt.gz",
    "GSE41258": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41258/matrix/GSE41258_series_matrix.txt.gz",
}
PLATFORM_BY_GSE = {
    "GSE21510": "GPL570",
    "GSE32323": "GPL570",
    "GSE44076": "GPL13667",
    "GSE41258": "GPL96",
}
PLATFORM_URLS = {
    "GPL570": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "GPL96": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
    "GPL13667": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13667&targ=self&view=data&form=text",
}
STABLE_FDR = 0.05
STABLE_ABS_LOGFC = 1.0


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return
    req = urllib.request.Request(url, headers={"User-Agent": "DINP-CRC-GEO-validation/1.0"})
    with urllib.request.urlopen(req, timeout=300) as response, path.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def parse_geo_values(text: str) -> list[str]:
    return next(csv.reader([text], delimiter="\t", quotechar='"', escapechar="\\"))


def load_series_matrix(accession: str) -> tuple[dict, pd.DataFrame, str]:
    path = CACHE / f"{accession}_series_matrix.txt.gz"
    download(MATRIX_URLS[accession], path)
    metadata: dict[str, list[list[str]]] = defaultdict(list)
    start = end = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            line = line.rstrip("\r\n")
            if line.lower() == "!series_matrix_table_begin":
                start = i
                break
            if line.startswith("!Sample_") and "\t" in line:
                key, values = line.split("\t", 1)
                metadata[key].append(parse_geo_values(values))
    if start is None:
        raise ValueError(f"{accession}: no series matrix table")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if line.strip().lower() == "!series_matrix_table_end":
                end = i
                break
    if end is None:
        raise ValueError(f"{accession}: no series matrix table end")
    nrows = end - start - 2
    table = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        skiprows=start + 1,
        nrows=nrows,
        quotechar='"',
        dtype={"ID_REF": str},
    )
    if table.columns[0] != "ID_REF":
        raise ValueError(f"{accession}: unexpected feature column {table.columns[0]}")
    return metadata, table, sha256(path)


def sample_records(metadata: dict) -> list[dict]:
    accessions = metadata["!Sample_geo_accession"][0]
    n = len(accessions)
    titles = metadata.get("!Sample_title", [[""] * n])[0]
    platform = metadata.get("!Sample_platform_id", [[""] * n])[0]
    chars = [dict() for _ in range(n)]
    for key, rows in metadata.items():
        if not key.startswith("!Sample_characteristics_ch1"):
            continue
        for row in rows:
            for i, value in enumerate(row):
                if ":" in value:
                    k, v = value.split(":", 1)
                    chars[i][k.strip().lower()] = v.strip()
    records = []
    for i in range(n):
        records.append({
            "gsm": accessions[i],
            "title": titles[i] if i < len(titles) else "",
            "platform": platform[i] if i < len(platform) else "",
            **chars[i],
        })
    return records


def split_symbols(value: str) -> list[str]:
    if not value or value.strip() in {"---", "NA", "N/A"}:
        return []
    parts = re.split(r"\s*///\s*", value.strip())
    out = []
    for part in parts:
        symbol = part.strip().upper()
        if symbol and symbol not in {"---", "NA", "N/A"} and symbol not in out:
            out.append(symbol)
    return out


def probe_symbols_from_annot_gz(path: pathlib.Path) -> tuple[dict[str, list[str]], str]:
    mapping: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = None
        id_col = symbol_col = None
        for row in reader:
            if row and row[0] == "ID" and any(x.strip().lower() == "gene symbol" for x in row):
                header = row
                id_col = next(i for i, x in enumerate(row) if x.strip().lower() == "id")
                symbol_col = next(i for i, x in enumerate(row) if x.strip().lower() == "gene symbol")
                break
        if header is None:
            raise ValueError(f"Annotation table header missing: {path}")
        for row in reader:
            if len(row) <= max(id_col, symbol_col):
                continue
            probe = row[id_col].strip()
            symbols = split_symbols(row[symbol_col])
            if probe and symbols:
                mapping[probe] = symbols
    return mapping, sha256(path)


def platform_mapping(platform: str) -> tuple[dict[str, list[str]], str]:
    if platform in {"GPL570", "GPL96"}:
        path = CACHE / f"{platform}.annot.gz"
        download(PLATFORM_URLS[platform], path)
        return probe_symbols_from_annot_gz(path)
    cache = CACHE / "GPL13667_probe_to_symbol.tsv"
    meta_path = CACHE / "GPL13667_probe_to_symbol.meta.json"
    if cache.exists() and cache.stat().st_size > 0 and meta_path.exists():
        mapping = {}
        with cache.open(encoding="utf-8", newline="") as f:
            for row in csv.reader(f, delimiter="\t"):
                if len(row) == 2:
                    mapping[row[0]] = row[1].split("|") if row[1] else []
        source_hash = json.loads(meta_path.read_text(encoding="utf-8"))["source_sha256"]
        return mapping, source_hash

    raw_path = CACHE / "GPL13667_platform_full.txt"
    if not raw_path.exists() or raw_path.stat().st_size == 0:
        req = requests.get(
            PLATFORM_URLS[platform],
            headers={"User-Agent": "DINP-CRC-GEO-validation/1.0"},
            timeout=(30, 300),
            stream=True,
        )
        req.raise_for_status()
        with raw_path.open("wb") as out:
            for chunk in req.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    out.write(chunk)
    platform_hash = sha256(raw_path)
    mapping = {}
    in_table = False
    table_end = False
    with raw_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        id_col = symbol_col = None
        for row in reader:
            if not row:
                continue
            first = row[0].strip().lower()
            if first == "!platform_table_begin":
                in_table = True
                header = next(reader)
                id_col = next(i for i, x in enumerate(header) if x.strip().lower() == "id")
                symbol_col = next(i for i, x in enumerate(header) if x.strip().lower() == "gene symbol")
                continue
            if first == "!platform_table_end":
                table_end = True
                break
            if not in_table or len(row) <= max(id_col, symbol_col):
                continue
            probe = row[id_col].strip()
            symbols = split_symbols(row[symbol_col])
            if probe and symbols:
                mapping[probe] = symbols
    if not table_end or len(mapping) < 10000:
        raise ValueError("GPL13667 platform table download or parse is incomplete")
    with cache.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for probe, symbols in mapping.items():
            writer.writerow([probe, "|".join(symbols)])
    meta_path.write_text(json.dumps({"source_sha256": platform_hash}, indent=2), encoding="utf-8")
    return mapping, platform_hash


def assigned_groups(accession: str, records: list[dict]) -> tuple[list[dict], dict]:
    out = []
    extra = Counter()
    for r in records:
        tissue = r.get("tissue", "").strip().lower()
        if accession in {"GSE21510", "GSE32323"}:
            if tissue == "cancer, homogenized":
                group = "Tumor"
            elif tissue == "normal, homogenized":
                group = "Normal"
            else:
                group = "Excluded"
            match = re.search(r"patient\s+([^, ]+)", r.get("title", ""), re.I)
            pair = match.group(1).lower().rstrip("b") if match else ""
        elif accession == "GSE44076":
            stype = r.get("sample type", "").strip().lower()
            group = {"tumor": "Tumor", "normal": "Normal"}.get(stype, "Excluded")
            pair = r.get("individual id", "").strip()
            if stype == "mucosa":
                extra["healthy_mucosa_excluded"] += 1
        else:
            included = r.get("included in analysis", "").strip().lower() == "yes"
            if included and tissue == "primary tumor":
                group = "Tumor"
            elif included and tissue == "normal colon":
                group = "Normal"
            else:
                group = "Excluded"
            pair = r.get("patient id", "").strip()
        out.append({**r, "group": group, "pair_id": pair})
    return out, dict(extra)


def patient_overlap_audit(records_by_accession: dict[str, list[dict]]) -> pd.DataFrame:
    """List selected patients represented in both Japanese GEO series."""
    maps: dict[str, dict[tuple[str, str], dict]] = {}
    for accession in ("GSE21510", "GSE32323"):
        maps[accession] = {
            (r["group"], r["pair_id"]): r
            for r in records_by_accession[accession]
            if r["group"] in {"Tumor", "Normal"} and r["pair_id"]
        }
    shared = sorted(set(maps["GSE21510"]) & set(maps["GSE32323"]))
    rows = []
    for group, patient_id in shared:
        left = maps["GSE21510"][(group, patient_id)]
        right = maps["GSE32323"][(group, patient_id)]
        rows.append({
            "group": group,
            "patient_id": patient_id,
            "GSE21510_GSM": left["gsm"],
            "GSE21510_title": left["title"],
            "GSE32323_GSM": right["gsm"],
            "GSE32323_title": right["title"],
        })
    return pd.DataFrame(rows, columns=[
        "group", "patient_id", "GSE21510_GSM", "GSE21510_title",
        "GSE32323_GSM", "GSE32323_title",
    ])


def map_probe_matrix(table: pd.DataFrame, mapping: dict[str, list[str]], accession: str) -> tuple[list[str], np.ndarray, dict[str, int], bool]:
    probes = table.iloc[:, 0].astype(str).tolist()
    raw = table.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    transformed = accession == "GSE41258"
    if transformed:
        positive = raw[np.isfinite(raw)]
        if positive.size and np.nanmin(positive) <= 0:
            raw = np.log2(np.maximum(raw, 0) + 1.0)
        else:
            raw = np.log2(raw)
    probe_indices: dict[str, list[int]] = defaultdict(list)
    for i, probe in enumerate(probes):
        for symbol in mapping.get(probe, []):
            probe_indices[symbol].append(i)
    genes = sorted(probe_indices)
    gene_matrix = np.vstack([
        np.nanmedian(raw[probe_indices[symbol], :], axis=0)
        for symbol in genes
    ])
    probe_counts = {g: len(probe_indices[g]) for g in genes}
    return genes, gene_matrix, probe_counts, transformed


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    out = np.full(pvalues.shape, np.nan, dtype=float)
    valid = np.isfinite(pvalues)
    if valid.any():
        out[valid] = multipletests(pvalues[valid], method="fdr_bh")[1]
    return out


def paired_test(matrix: np.ndarray, records: list[dict], indices: list[int]) -> tuple[np.ndarray, np.ndarray, int]:
    group_by_pair: dict[str, dict[str, list[int]]] = defaultdict(lambda: {"Tumor": [], "Normal": []})
    for i in indices:
        r = records[i]
        if r["pair_id"]:
            group_by_pair[r["pair_id"]][r["group"]].append(i)
    pairs = []
    for pair_id, groups in sorted(group_by_pair.items()):
        if groups["Tumor"] and groups["Normal"]:
            tumor = np.nanmean(matrix[:, groups["Tumor"]], axis=1)
            normal = np.nanmean(matrix[:, groups["Normal"]], axis=1)
            pairs.append((tumor, normal))
    if not pairs:
        raise ValueError("No complete tumor-normal pairs")
    tumors = np.column_stack([x[0] for x in pairs])
    normals = np.column_stack([x[1] for x in pairs])
    logfc = np.nanmean(tumors - normals, axis=1)
    p = ttest_rel(tumors, normals, axis=1, nan_policy="omit").pvalue
    return logfc, np.asarray(p, dtype=float), len(pairs)


def unpaired_test(matrix: np.ndarray, records: list[dict], indices: list[int]) -> tuple[np.ndarray, np.ndarray, int]:
    tumor_i = [i for i in indices if records[i]["group"] == "Tumor"]
    normal_i = [i for i in indices if records[i]["group"] == "Normal"]
    tumors = matrix[:, tumor_i]
    normals = matrix[:, normal_i]
    logfc = np.nanmean(tumors, axis=1) - np.nanmean(normals, axis=1)
    p = ttest_ind(tumors, normals, axis=1, equal_var=False, nan_policy="omit").pvalue
    pair_ids_t = {records[i]["pair_id"] for i in tumor_i if records[i]["pair_id"]}
    pair_ids_n = {records[i]["pair_id"] for i in normal_i if records[i]["pair_id"]}
    return logfc, np.asarray(p, dtype=float), len(pair_ids_t & pair_ids_n)


def patient_cluster_test(matrix: np.ndarray, records: list[dict], indices: list[int]) -> tuple[np.ndarray, np.ndarray, int]:
    tumor_i = [i for i in indices if records[i]["group"] == "Tumor"]
    normal_i = [i for i in indices if records[i]["group"] == "Normal"]
    tumor_mask = np.asarray([records[i]["group"] == "Tumor" for i in indices], dtype=bool)
    local_records = [records[i] for i in indices]
    x = matrix[:, indices]
    n_t = int(tumor_mask.sum())
    n_n = int((~tumor_mask).sum())
    beta = np.nanmean(x[:, tumor_mask], axis=1) - np.nanmean(x[:, ~tumor_mask], axis=1)
    residual = x.copy()
    residual[:, tumor_mask] -= np.nanmean(x[:, tumor_mask], axis=1)[:, None]
    residual[:, ~tumor_mask] -= np.nanmean(x[:, ~tumor_mask], axis=1)[:, None]
    clusters: dict[str, list[int]] = defaultdict(list)
    for j, r in enumerate(local_records):
        key = r.get("pair_id", "").strip()
        if not key:
            key = f"unmatched_{r['gsm']}"
        clusters[key].append(j)
    g = len(clusters)
    meat_diag = np.zeros(matrix.shape[0], dtype=float)
    for member_idx in clusters.values():
        ix = np.asarray(member_idx, dtype=int)
        t_ix = ix[tumor_mask[ix]]
        n_ix = ix[~tumor_mask[ix]]
        sum_t = residual[:, t_ix].sum(axis=1) if len(t_ix) else 0.0
        sum_n = residual[:, n_ix].sum(axis=1) if len(n_ix) else 0.0
        influence = sum_t / n_t - sum_n / n_n
        meat_diag += influence * influence
    n = len(indices)
    correction = (g / (g - 1)) * ((n - 1) / (n - 2))
    se = np.sqrt(meat_diag * correction)
    statistic = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    p = 2.0 * student_t.sf(np.abs(statistic), df=max(g - 1, 1))
    zero_se = (se == 0) & np.isfinite(beta)
    p[zero_se & (beta == 0)] = 1.0
    p[zero_se & (beta != 0)] = 0.0
    common_ids = {records[i]["pair_id"] for i in tumor_i if records[i]["pair_id"]} & {
        records[i]["pair_id"] for i in normal_i if records[i]["pair_id"]
    }
    return beta, p, len(common_ids)


def run() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    input_df = pd.read_csv(INPUT_113, sep="\t")
    targets = sorted(set(input_df["gene_symbol"].astype(str).str.strip().str.upper()))
    if len(targets) != 113:
        raise ValueError(f"Expected 113 target symbols, found {len(targets)}")

    checksums = {}
    result_rows = []
    count_rows = []
    coverage_rows = []
    records_by_accession = {}
    for accession in ACCESSIONS:
        metadata, table, matrix_hash = load_series_matrix(accession)
        records = sample_records(metadata)
        records, extra = assigned_groups(accession, records)
        records_by_accession[accession] = records
        platform = PLATFORM_BY_GSE[accession]
        mapping, platform_hash = platform_mapping(platform)
        checksums[f"{accession}_series_matrix_sha256"] = matrix_hash
        checksums[f"{platform}_annotation_sha256"] = platform_hash
        genes, gene_matrix, probe_counts, transformed = map_probe_matrix(table, mapping, accession)
        sample_columns = list(table.columns[1:])
        if sample_columns != [r["gsm"] for r in records]:
            records_by_gsm = {r["gsm"]: r for r in records}
            records = [records_by_gsm[gsm] for gsm in sample_columns]
        selected = [i for i, r in enumerate(records) if r["group"] in {"Tumor", "Normal"}]
        n_tumor = sum(records[i]["group"] == "Tumor" for i in selected)
        n_normal = sum(records[i]["group"] == "Normal" for i in selected)
        t_ids = {records[i]["pair_id"] for i in selected if records[i]["group"] == "Tumor" and records[i]["pair_id"]}
        n_ids = {records[i]["pair_id"] for i in selected if records[i]["group"] == "Normal" and records[i]["pair_id"]}
        n_pairs = len(t_ids & n_ids)
        if accession in {"GSE32323", "GSE44076"}:
            method = "paired t-test"
            logfc, pvalues, n_pairs_used = paired_test(gene_matrix, records, selected)
        elif accession == "GSE21510":
            method = "OLS with patient-cluster robust SE"
            logfc, pvalues, n_pairs_used = patient_cluster_test(gene_matrix, records, selected)
        else:
            method = "OLS with patient-cluster robust SE"
            logfc, pvalues, n_pairs_used = patient_cluster_test(gene_matrix, records, selected)
        fdr = bh_adjust(pvalues)
        gene_index = {symbol: i for i, symbol in enumerate(genes)}
        covered = set(gene_index) & set(targets)
        coverage_rows.append({
            "accession": accession,
            "platform": platform,
            "target_genes_total": len(targets),
            "target_genes_with_probe": len(covered),
            "target_genes_without_probe": len(set(targets) - covered),
            "without_probe_symbols": ";".join(sorted(set(targets) - covered)),
            "unique_mapped_genes_tested_for_FDR": len(genes),
            "mapped_probe_count": sum(probe_counts.values()),
            "series_matrix_feature_count": int(table.shape[0]),
            "expression_transform": "log2 intensity" if transformed else "series-matrix scale retained",
            "statistical_method": method,
        })
        selection = {
            "GSE21510": "homogenized cancer vs homogenized normal; LCM tumors excluded because GEO says these were separately normalized and no LCM normal group is present",
            "GSE32323": "homogenized cancer vs homogenized normal; cell lines excluded",
            "GSE44076": "98 paired tumor/adjacent-normal samples; 50 healthy mucosa samples excluded from the primary comparison",
            "GSE41258": "included in analysis=Yes; Primary Tumor vs Normal Colon; other tissue types and No rows excluded",
        }[accession]
        count_rows.append({
            "accession": accession,
            "platform": platform,
            "series_matrix_total_samples": len(records),
            "included_tumor_samples": n_tumor,
            "included_normal_samples": n_normal,
            "complete_patient_pairs_available": n_pairs,
            "complete_patient_pairs_used": n_pairs_used,
            "healthy_mucosa_excluded": extra.get("healthy_mucosa_excluded", 0),
            "statistical_method": method,
            "sample_selection": selection,
        })
        for gene in targets:
            if gene not in gene_index:
                result_rows.append({
                    "gene_symbol": gene,
                    "accession": accession,
                    "platform": platform,
                    "n_tumor": n_tumor,
                    "n_normal": n_normal,
                    "n_pairs_used": n_pairs_used,
                    "statistical_method": method,
                    "probes_collapsed": 0,
                    "log2FC_tumor_minus_normal": np.nan,
                    "P_value": np.nan,
                    "FDR_BH_all_mapped_genes": np.nan,
                    "passes_FDR_and_abs_log2FC": False,
                    "direction": "not represented",
                })
                continue
            gi = gene_index[gene]
            fc = float(logfc[gi])
            q = float(fdr[gi]) if np.isfinite(fdr[gi]) else np.nan
            p = float(pvalues[gi]) if np.isfinite(pvalues[gi]) else np.nan
            direction = "Tumor_up" if fc > 0 else "Tumor_down" if fc < 0 else "No_change"
            hit = bool(np.isfinite(q) and q < STABLE_FDR and abs(fc) > STABLE_ABS_LOGFC)
            result_rows.append({
                "gene_symbol": gene,
                "accession": accession,
                "platform": platform,
                "n_tumor": n_tumor,
                "n_normal": n_normal,
                "n_pairs_used": n_pairs_used,
                "statistical_method": method,
                "probes_collapsed": probe_counts[gene],
                "log2FC_tumor_minus_normal": fc,
                "P_value": p,
                "FDR_BH_all_mapped_genes": q,
                "passes_FDR_and_abs_log2FC": hit,
                "direction": direction,
            })
        print(
            f"{accession}: total={len(records)}; tumor={n_tumor}; normal={n_normal}; "
            f"pairs_available={n_pairs}; pairs_used={n_pairs_used}; "
            f"genes_tested={len(genes)}; target_coverage={len(covered)}/113",
            flush=True,
        )

    per_cohort = pd.DataFrame(result_rows)
    counts = pd.DataFrame(count_rows)
    coverage = pd.DataFrame(coverage_rows)
    overlap = patient_overlap_audit(records_by_accession)
    overlap.to_csv(OUT / "GEO_GSE21510_GSE32323_patient_id_overlap.tsv", sep="\t", index=False)
    stable_rows = []
    for gene in targets:
        sub = per_cohort[per_cohort["gene_symbol"] == gene].copy()
        hits = sub[sub["passes_FDR_and_abs_log2FC"]]
        directions = set(hits["direction"]) - {"not represented", "No_change"}
        cohort_hits = ";".join(hits["accession"].tolist())
        up_n = int((hits["direction"] == "Tumor_up").sum())
        down_n = int((hits["direction"] == "Tumor_down").sum())
        same_direction = len(directions) <= 1
        stable = len(hits) >= 3 and same_direction
        japanese_hit = bool(hits["accession"].isin(["GSE21510", "GSE32323"]).any())
        source_hits = int(japanese_hit) + int(hits["accession"].eq("GSE44076").any()) + int(hits["accession"].eq("GSE41258").any())
        source_directions = set(hits.loc[
            hits["accession"].isin(["GSE21510", "GSE32323", "GSE44076", "GSE41258"]), "direction"
        ]) - {"not represented", "No_change"}
        source_direction_consistent = len(source_directions) <= 1
        stable_sources = source_hits >= 3 and source_direction_consistent
        row = {
            "gene_symbol": gene,
            "significant_cohorts_FDR_lt_0.05_abs_log2FC_gt_1": len(hits),
            "cohorts_meeting_threshold": cohort_hits,
            "significant_up_cohorts": up_n,
            "significant_down_cohorts": down_n,
            "direction_consistent_among_significant": same_direction,
            "stable_replicated_in_at_least_3_cohorts": stable,
            "independent_source_groups_with_hit": source_hits,
            "stable_replicated_in_3_independent_source_groups": stable_sources,
        }
        for accession in ACCESSIONS:
            r = sub[sub["accession"] == accession].iloc[0]
            row[f"{accession}_log2FC"] = r["log2FC_tumor_minus_normal"]
            row[f"{accession}_FDR"] = r["FDR_BH_all_mapped_genes"]
            row[f"{accession}_direction"] = r["direction"]
            row[f"{accession}_threshold_hit"] = bool(r["passes_FDR_and_abs_log2FC"])
        stable_rows.append(row)
    stability = pd.DataFrame(stable_rows).sort_values(
        ["stable_replicated_in_at_least_3_cohorts", "significant_cohorts_FDR_lt_0.05_abs_log2FC_gt_1", "gene_symbol"],
        ascending=[False, False, True],
    )
    stable = stability[stability["stable_replicated_in_at_least_3_cohorts"]].copy()
    stable_sources = stability[stability["stable_replicated_in_3_independent_source_groups"]].copy()

    counts.to_csv(OUT / "GEO_cohort_sample_counts.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "GEO_platform_gene_coverage.tsv", sep="\t", index=False)
    per_cohort.to_csv(OUT / "113_gene_per_cohort_DE_results.tsv", sep="\t", index=False)
    stability.to_csv(OUT / "113_gene_stability_summary.tsv", sep="\t", index=False)
    stable.to_csv(OUT / "stable_genes_replicated_in_3_or_more_cohorts.tsv", sep="\t", index=False)
    stable_sources.to_csv(OUT / "stable_genes_replicated_in_3_source_groups.tsv", sep="\t", index=False)
    (OUT / "source_checksums.json").write_text(json.dumps(checksums, indent=2), encoding="utf-8")

    up_stable = int((stable["significant_up_cohorts"] > 0).sum())
    down_stable = int((stable["significant_down_cohorts"] > 0).sum())
    stable_list = ", ".join(stable["gene_symbol"].tolist()) if len(stable) else "无"
    source_stable_list = ", ".join(stable_sources["gene_symbol"].tolist()) if len(stable_sources) else "无"
    shared_tumor = int((overlap["group"] == "Tumor").sum())
    shared_normal = int((overlap["group"] == "Normal").sum())
    report_lines = [
        "# 113 个 DINP–CRC 候选基因在四个 GEO 队列中的复现",
        "",
        f"- 分析日期：{date.today().isoformat()}",
        "- 目标基因：来自本项目已保存的 113 基因交集明细。",
        "- 复现阈值：每队列以全平台唯一映射基因做 BH 校正；FDR < 0.05 且 |log2FC| > 1 记为该队列复现。",
        "- 稳定复现定义：至少 3/4 队列达到上述阈值，且达到阈值的队列方向一致。",
        "- 独立来源敏感性分析：GSE21510 与 GSE32323 的样本标题含有相同患者编号；因此保守地将二者合并为一个来源组，再加 GSE44076、GSE41258，共 3 个来源组。独立来源复现要求三个来源组均至少有一个队列达到显著性与效应量阈值，且同向。",
        "- 差异方向：log2FC = 肿瘤 − 正常；正值为肿瘤上调，负值为肿瘤下调。",
        "- 探针处理：各平台按 NCBI GEO 平台注释将探针映射到 gene symbol；多探针基因按样本取中位数；未做别名映射。",
        "- GEO 系列矩阵已归一化；GSE41258 的 PLIER 信号按 log2(intensity) 转换，其余三个矩阵保留提交者提供的 log 表达尺度。",
        "- GSE21510 的 LCM 肿瘤与全组织样本分开归一化，且没有 LCM 正常对照；主比较仅用 homogenized 肿瘤/正常。",
        "- GSE44076 的 50 个健康黏膜样本不纳入主肿瘤对癌旁分析；仅使用 98 个肿瘤与 98 个配对癌旁样本。",
        "- GSE41258 限定 included in analysis=Yes 的 Primary Tumor 与 Normal Colon；重复/配对患者按 patient ID 聚类计算稳健标准误。",
        "",
        "## 队列样本数",
        "",
        "| GEO 队列 | 矩阵总样本 | 纳入肿瘤 | 纳入正常 | 可配对患者 |",
        "|---|---:|---:|---:|---:|",
    ]
    for _, row in counts.iterrows():
        report_lines.append(
            f"| {row['accession']} | {row['series_matrix_total_samples']} | "
            f"{row['included_tumor_samples']} | {row['included_normal_samples']} | "
            f"{row['complete_patient_pairs_available']} |"
        )
    report_lines.extend([
        "",
        "## 稳定复现结果",
        "",
        f"- 稳定复现基因数：**{len(stable)} / 113**。",
        f"- 肿瘤上调方向：{up_stable} 个；肿瘤下调方向：{down_stable} 个。",
        f"- 基因：{stable_list}。",
        f"- 按重叠患者合并为 3 个来源组后，仍在全部 3 组同向复现：**{len(stable_sources)} / 113**。",
        f"- 三来源组复现基因：{source_stable_list}。",
        f"- GSE21510 与 GSE32323 共享患者编号：肿瘤 {shared_tumor} 个、正常 {shared_normal} 个；逐例清单见 GEO_GSE21510_GSE32323_patient_id_overlap.tsv。",
        "",
        "## 输出文件",
        "",
        "- GEO_cohort_sample_counts.tsv：各队列样本纳入数与分析方法。",
        "- GEO_platform_gene_coverage.tsv：113 个基因的平台探针覆盖率及测试基因数。",
        "- 113_gene_per_cohort_DE_results.tsv：113 基因在每个队列的 log2FC、P、FDR 与方向。",
        "- 113_gene_stability_summary.tsv：每个目标基因的跨队列复现次数与方向。",
        "- stable_genes_replicated_in_3_or_more_cohorts.tsv：达到预设稳定标准的基因。",
        "- stable_genes_replicated_in_3_source_groups.tsv：按重叠患者保守合并后，在三个来源组均复现的基因。",
        "- GEO_GSE21510_GSE32323_patient_id_overlap.tsv：两个 GEO 系列之间共享患者编号的样本清单。",
        "- source_checksums.json：下载矩阵与注释来源 SHA-256。",
        "",
        "## 复现边界",
        "",
        "平台注释的基因 symbol 来自 GEO/GPL；本分析按 symbol 精确映射，不做当前 HGNC 别名回溯。队列间采用与样本设计相符的配对检验或按患者聚类的稳健检验；候选基因是否进入稳定列表还同时要求效应量阈值和同向复现。GSE21510 与 GSE32323 患者编号重叠提示其并非完全独立验证；因此同时给出按四个 GEO 系列统计的结果和将这两个系列合并后的保守来源组结果。",
    ])
    (OUT / "README.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"Stable genes >=3 cohorts: {len(stable)}")
    print("Stable list:", stable_list)
    print("Results:", OUT)


if __name__ == "__main__":
    run()
