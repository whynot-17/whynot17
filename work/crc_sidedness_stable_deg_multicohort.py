#!/usr/bin/env python3
"""Multi-cohort discovery of stable left- versus right-sided CRC genes.

The analysis is intentionally phenotype-first: no ferroptosis, PUFA or SLC7A11
prior is used during DEG discovery.  Each cohort is analyzed in its own native
expression scale, then right-minus-left standardized effects are combined with
a random-effects meta-analysis.

Bulk discovery cohorts
-----------------------
TCGA-COAD, GSE39582, GSE103479 and GSE41258.

Single-cell validation cohorts
------------------------------
GSE200997, GSE132465, GSE144735 and GSE188711.  Raw UMI counts are aggregated
to patient-level tumor-epithelial pseudobulk before testing.  GSE200997 uses
the existing marker-defined malignant-epithelial calls.  GSE132465 and
GSE144735 use the official tumor epithelial annotation (a tumor-epithelial
proxy, not a CNV-proven malignant label).  GSE188711 uses the pre-existing
primary putative malignant epithelial cluster 4 and is flagged as a 3-versus-3
sensitivity cohort.

Stable right-high / left-high definition
-----------------------------------------
At least three cohorts contribute a finite Hedges-g effect; every contributing
cohort has the same effect direction; and the random-effects meta-analysis has
BH-FDR < 0.05.  This strict list is separate from the full effect tables, so
directional but underpowered genes remain auditable.

Inputs are expected under work/data and work/gse*.  Large source matrices are
not written to Git by this script; only scripts, audits, summaries and result
tables should be committed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import norm, ttest_ind

try:
    import GEOparse  # type: ignore
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("Install GEOparse first: pip install GEOparse") from exc


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
DATA = WORK / "data"
OUT = ROOT / "outputs"

BULK_GEO = ["GSE39582", "GSE103479", "GSE41258"]
BULK_COHORTS = BULK_GEO + ["TCGA-COAD"]
SC_COHORTS = ["GSE200997", "GSE132465", "GSE144735", "GSE188711"]

GEO_MATRIX_PATHS = {
    "GSE39582": DATA / "GSE39582_series_matrix.clean.txt.gz",
    "GSE41258": DATA / "GSE41258_series_matrix.txt.gz",
}
GSE103479_ANNOTATED = DATA / "GSE103479_log2_RMA_annotated.txt.gz"
GSE103479_METADATA = DATA / "GSE103479_series_matrix.txt.gz"
TCGA_EXPRESSION = DATA / "TCGA_COAD_HiSeqV2.gz"
TCGA_CLINICAL = DATA / "TCGA_COAD_clinicalMatrix.tsv"

GSE200997_DIR = WORK / "gse200997"
GSE132465_DIR = WORK / "gse132465"
GSE188711_DIR = WORK / "gse188711"

GSE200997_ANNOT = GSE200997_DIR / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"
GSE200997_MATRIX = GSE200997_DIR / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
GSE200997_SELECTION = GSE200997_DIR / "tumor_marker_scores.csv"

GSE132465_ANNOT = GSE132465_DIR / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
GSE132465_MATRIX = GSE132465_DIR / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
GSE132465_SOFT = GSE132465_DIR / "GSE132465_family.soft.gz"

GSE144735_ANNOT = DATA / "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"
GSE144735_MATRIX = DATA / "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"
GSE144735_SOFT = DATA / "GSE144735_family.soft.gz"

GSE188711_META = GSE188711_DIR / "cell_metadata.csv"
GSE188711_RAW = DATA / "GSE188711_RAW"

GSE188711_SAMPLES = {
    "GSM5688706": {"tag": "WGC", "side": "left", "label": "L1"},
    "GSM5688707": {"tag": "JCA", "side": "left", "label": "L2"},
    "GSM5688708": {"tag": "LS-CRC3", "side": "left", "label": "L3"},
    "GSM5688709": {"tag": "RS-CRC1", "side": "right", "label": "R1"},
    "GSM5688710": {"tag": "R_CRC3", "side": "right", "label": "R2"},
    "GSM5688711": {"tag": "R_CRC4", "side": "right", "label": "R3"},
}


RIGHT_PATTERNS = [
    r"\bcecum\b", r"\bcaecum\b", r"\bcecal\b", r"\bcaecal\b",
    r"\bascending\b", r"hepatic\s+flexure", r"\btransverse\b",
    r"\bproximal\b", r"right\s+(?:sided\s+)?colon", r"right[- ]sided",
]
LEFT_PATTERNS = [
    r"\bdescending\b", r"splenic\s+flexure", r"\bsigmoid\b",
    r"\brectosigmoid\b", r"\bdistal\b", r"left\s+(?:sided\s+)?colon",
    r"left[- ]sided",
]
RECTUM_PATTERNS = [r"\brectum\b", r"\brectal\b"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--skip-gse188711", action="store_true", help="Skip the low-powered GSE188711 sensitivity cohort")
    return parser.parse_args()


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    ok = np.isfinite(values)
    if not ok.any():
        return out
    idx = np.where(ok)[0]
    order = np.argsort(values[ok])
    ranked = values[ok][order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[idx] = restored
    return out


def text_values(row: dict[str, list[str]], key: str) -> list[str]:
    values = row.get(key, [])
    return [str(v).strip().strip('"') for v in values if str(v).strip()]


def sample_level_text(row: dict[str, list[str]]) -> str:
    keys = ["title", "source_name_ch1"] + sorted(k for k in row if k.startswith("characteristics_ch1"))
    return " | ".join(f"{key}: {' | '.join(text_values(row, key))}" for key in keys if text_values(row, key))


def parse_series_matrix_metadata(path: Path) -> list[dict[str, list[str]]]:
    rows: list[dict[str, list[str]]] = []
    with gzip.open(path, "rt", newline="", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            fields = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))
            key = fields[0][len("!Sample_"):]
            while len(rows) < len(fields) - 1:
                rows.append({})
            for i, value in enumerate(fields[1:]):
                rows[i].setdefault(key, []).append(value)
    return rows


def parse_soft_samples(path: Path) -> dict[str, dict[str, list[str]]]:
    records: dict[str, dict[str, list[str]]] = {}
    current_id: str | None = None
    current: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                if current_id is not None:
                    records[current_id] = current
                current_id = line.split("=", 1)[1].strip()
                current = {}
            elif current_id is not None and line.startswith("!Sample_") and " = " in line:
                key, value = line[8:].split(" = ", 1)
                current.setdefault(key, []).append(value.strip().strip('"'))
    if current_id is not None:
        records[current_id] = current
    return records


def classify_side(text: str) -> tuple[str, str]:
    low = text.lower()
    right = [p for p in RIGHT_PATTERNS if re.search(p, low)]
    left = [p for p in LEFT_PATTERNS if re.search(p, low)]
    rectum = [p for p in RECTUM_PATTERNS if re.search(p, low)]
    if right and not left:
        return "right", right[0]
    if left and not right:
        if rectum and not re.search(r"rectosigmoid", low):
            return "exclude_rectum", rectum[0]
        return "left", left[0]
    if rectum and not right and not left:
        return "exclude_rectum", rectum[0]
    if right and left:
        return "ambiguous", f"right={right[0]};left={left[0]}"
    return "unknown", ""


def clean_symbol(value: object) -> str | None:
    text = str(value).strip().strip('"')
    if not text or text.lower() in {"nan", "none", "na", "n/a", "---", "-"}:
        return None
    tokens = [t.strip() for t in re.split(r"\s*///\s*|\s*//\s*|\s*;\s*|\s*,\s*|\s+", text)]
    for token in tokens:
        token = token.strip().strip('"')
        if not token or token in {"---", "-"} or token.lower() in {"na", "nan", "none"}:
            continue
        if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", token):
            return token.upper()
    return None


def gene_symbol_column(table: pd.DataFrame) -> str:
    preferred = ["Gene Symbol", "GENE_SYMBOL", "Gene symbol", "Symbol", "SYMBOL", "gene_assignment"]
    for column in preferred:
        if column in table.columns:
            return column
    for column in table.columns:
        low = str(column).lower()
        if "gene" in low and "symbol" in low:
            return str(column)
    raise ValueError(f"No gene-symbol column found in annotation: {list(table.columns)[:30]}")


def gpl_mapping(platform_id: str) -> pd.Series:
    path = DATA / f"{platform_id}.annot.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing local platform annotation: {path}")
    gpl = GEOparse.get_GEO(filepath=str(path), geotype="GPL", silent=True).table
    probe_ids = gpl["ID"].astype(str) if "ID" in gpl.columns else pd.Series(gpl.index.astype(str), index=gpl.index)
    symbols = gpl[gene_symbol_column(gpl)].map(clean_symbol)
    return pd.Series(symbols.to_numpy(), index=probe_ids.to_numpy(), dtype="object")


def resolve_sample_columns(raw_columns: pd.Index, records: list[dict[str, list[str]]], cohort: str) -> tuple[list[str], pd.DataFrame]:
    raw_set = set(map(str, raw_columns))
    selected_columns: list[str] = []
    audit_rows: list[dict[str, object]] = []
    for row in records:
        sid = (text_values(row, "geo_accession") or [""])[0]
        title = (text_values(row, "title") or [""])[0]
        candidates = [title, sid]
        match = next((candidate for candidate in candidates if candidate in raw_set), None)
        if match is None:
            continue
        side, side_match = classify_side(sample_level_text(row))
        selected_columns.append(match)
        audit_rows.append({
            "cohort": cohort, "sample_id": sid, "expression_column": match,
            "title": title, "side": side, "side_match": side_match,
            "sample_level_text": sample_level_text(row),
            "platform_id": (text_values(row, "platform_id") or [""])[0],
        })
    if not selected_columns:
        raise ValueError(f"{cohort}: could not match sample metadata to expression columns")
    return selected_columns, pd.DataFrame(audit_rows).set_index("sample_id")


def collapse_probe_matrix(raw: pd.DataFrame, mapping: pd.Series, sample_columns: list[str]) -> pd.DataFrame:
    numeric = raw.loc[:, sample_columns].apply(pd.to_numeric, errors="coerce")
    probe_ids = raw.index.astype(str)
    symbols = mapping.reindex(probe_ids).to_numpy()
    keep = pd.notna(symbols)
    numeric = numeric.loc[keep]
    numeric.index = pd.Index(symbols[keep], name="gene")
    collapsed = numeric.groupby(level=0, sort=False).median()
    return collapsed.T


def load_generic_geo(cohort: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    path = GEO_MATRIX_PATHS[cohort]
    if not path.exists():
        raise FileNotFoundError(path)
    records = parse_series_matrix_metadata(path)
    raw = pd.read_csv(path, sep="\t", compression="gzip", comment="!", index_col=0, low_memory=False)
    columns, audit = resolve_sample_columns(raw.columns, records, cohort)
    platforms = sorted({(text_values(row, "platform_id") or [""])[0] for row in records if text_values(row, "platform_id")})
    if len(platforms) != 1:
        raise ValueError(f"{cohort}: expected one platform for full-matrix collapse, found {platforms}")
    expr = collapse_probe_matrix(raw, gpl_mapping(platforms[0]), columns)
    column_to_sample = pd.Series(audit.index.to_numpy(), index=audit["expression_column"].astype(str))
    expr.index = expr.index.map(column_to_sample)
    expr = expr.loc[expr.index.notna()]
    expr = expr[~expr.index.duplicated(keep="first")]
    audit = audit.loc[audit.index.isin(expr.index)]
    expr = expr.loc[audit.index]
    return expr, audit, {"source": str(path), "platform": platforms[0], "n_samples": len(expr), "n_genes": expr.shape[1]}


def load_gse103479() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if not GSE103479_ANNOTATED.exists() or not GSE103479_METADATA.exists():
        raise FileNotFoundError("GSE103479 annotated matrix and series metadata are required")
    records = parse_series_matrix_metadata(GSE103479_METADATA)
    raw = pd.read_csv(GSE103479_ANNOTATED, sep="\t", compression="gzip", index_col=0, low_memory=False)
    columns, audit = resolve_sample_columns(raw.columns, records, "GSE103479")
    if "Gene.Symbol" not in raw.columns:
        raise ValueError("GSE103479 annotated matrix is missing Gene.Symbol")
    numeric = raw.loc[:, columns].apply(pd.to_numeric, errors="coerce")
    symbols = raw["Gene.Symbol"].map(clean_symbol)
    keep = symbols.notna().to_numpy()
    numeric = numeric.loc[keep]
    numeric.index = pd.Index(symbols.loc[keep].to_numpy(), name="gene")
    expr = numeric.groupby(level=0, sort=False).median().T
    column_to_sample = pd.Series(audit.index.to_numpy(), index=audit["expression_column"].astype(str))
    expr.index = expr.index.map(column_to_sample)
    expr = expr.loc[expr.index.notna()]
    expr = expr[~expr.index.duplicated(keep="first")]
    audit = audit.loc[audit.index.isin(expr.index)]
    expr = expr.loc[audit.index]
    return expr, audit, {"source": str(GSE103479_ANNOTATED), "metadata_source": str(GSE103479_METADATA), "platform": "Almac Xcel annotated", "n_samples": len(expr), "n_genes": expr.shape[1]}


def load_tcga() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if not TCGA_EXPRESSION.exists() or not TCGA_CLINICAL.exists():
        raise FileNotFoundError("TCGA expression and clinical matrix are required")
    expr_raw = pd.read_csv(TCGA_EXPRESSION, sep="\t", compression="gzip", index_col=0, low_memory=False)
    expr_raw.columns = expr_raw.columns.astype(str)
    clinical = pd.read_csv(TCGA_CLINICAL, sep="\t", dtype=str, low_memory=False).fillna("")
    clinical["sampleID"] = clinical["sampleID"].astype(str)
    clinical = clinical[clinical["sample_type"].str.lower().eq("primary tumor")].copy()
    clinical = clinical.drop_duplicates("sampleID", keep="first")
    common = [sid for sid in clinical["sampleID"] if sid in expr_raw.columns]
    if len(common) < 100:
        raise ValueError(f"TCGA: only {len(common)} primary-tumor expression samples matched clinical metadata")
    expr = expr_raw.loc[:, common].apply(pd.to_numeric, errors="coerce").T
    rows = []
    for _, row in clinical[clinical["sampleID"].isin(common)].iterrows():
        text = f"anatomic_neoplasm_subdivision: {row.get('anatomic_neoplasm_subdivision', '')}"
        side, side_match = classify_side(text)
        rows.append({
            "cohort": "TCGA-COAD", "sample_id": row["sampleID"], "expression_column": row["sampleID"],
            "title": row.get("anatomic_neoplasm_subdivision", ""), "side": side,
            "side_match": side_match, "sample_level_text": text,
            "platform_id": "GDC/HiSeqV2",
            "sample_type": row.get("sample_type", ""), "MSI": row.get("MSI_updated_Oct62011", row.get("microsatellite_instability", "")),
        })
    audit = pd.DataFrame(rows).set_index("sample_id")
    expr = expr.loc[audit.index]
    return expr, audit, {"source": str(TCGA_EXPRESSION), "clinical_source": str(TCGA_CLINICAL), "platform": "GDC/HiSeqV2", "n_samples": len(expr), "n_genes": expr.shape[1]}


def differential(expr: pd.DataFrame, audit: pd.DataFrame, cohort: str, modality: str) -> pd.DataFrame:
    meta = audit[audit["side"].isin(["right", "left"])].copy()
    meta = meta.loc[meta.index.intersection(expr.index)]
    right_ids = meta.index[meta["side"].eq("right")]
    left_ids = meta.index[meta["side"].eq("left")]
    if len(right_ids) < 2 or len(left_ids) < 2:
        raise ValueError(f"{cohort}: insufficient sided samples (right={len(right_ids)}, left={len(left_ids)})")
    right = expr.loc[right_ids].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    left = expr.loc[left_ids].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    n_right = np.isfinite(right).sum(axis=0)
    n_left = np.isfinite(left).sum(axis=0)
    right_mean = np.divide(np.nansum(right, axis=0), n_right, out=np.full(right.shape[1], np.nan), where=n_right > 0)
    left_mean = np.divide(np.nansum(left, axis=0), n_left, out=np.full(left.shape[1], np.nan), where=n_left > 0)
    delta = right_mean - left_mean
    right_var = np.nanvar(right, axis=0, ddof=1)
    left_var = np.nanvar(left, axis=0, ddof=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        test = ttest_ind(right, left, axis=0, equal_var=False, nan_policy="omit")
    p_values = np.asarray(test.pvalue, dtype=float)
    pooled_den = (n_right - 1) * right_var + (n_left - 1) * left_var
    pooled_sd = np.sqrt(np.divide(pooled_den, n_right + n_left - 2, out=np.full_like(delta, np.nan), where=(n_right + n_left) > 2))
    d = np.divide(delta, pooled_sd, out=np.full_like(delta, np.nan), where=pooled_sd > 0)
    j = 1.0 - 3.0 / (4.0 * (n_right + n_left) - 9.0)
    hedges_g = d * j
    se_g = np.sqrt(np.divide(n_right + n_left, n_right * n_left, out=np.full_like(delta, np.nan), where=(n_right > 0) & (n_left > 0)) + np.divide(hedges_g**2, 2.0 * (n_right + n_left - 2), out=np.full_like(delta, np.nan), where=(n_right + n_left) > 2))
    out = pd.DataFrame({
        "cohort": cohort, "modality": modality, "gene": expr.columns.astype(str),
        "n_right": n_right, "n_left": n_left, "right_mean": right_mean, "left_mean": left_mean,
        "right_minus_left": delta, "hedges_g": hedges_g, "se_g": se_g,
        "welch_p": p_values,
    })
    out["welch_fdr"] = bh_adjust(out["welch_p"].to_numpy())
    out["direction"] = np.where(out["right_minus_left"] > 0, "right_high", np.where(out["right_minus_left"] < 0, "left_high", "tie"))
    return out


def random_effects_meta(tables: list[pd.DataFrame], modality: str) -> pd.DataFrame:
    if not tables:
        return pd.DataFrame()
    combined = pd.concat(tables, ignore_index=True)
    records: list[dict[str, object]] = []
    for gene, group in combined.groupby("gene", sort=True):
        d = group[np.isfinite(group["hedges_g"]) & np.isfinite(group["se_g"]) & (group["se_g"] > 0)]
        if d.empty:
            continue
        y = d["hedges_g"].to_numpy(float)
        v = d["se_g"].to_numpy(float) ** 2
        w = 1.0 / v
        fixed = float(np.sum(w * y) / np.sum(w))
        q = float(np.sum(w * (y - fixed) ** 2))
        df_q = len(y) - 1
        c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if np.sum(w) else np.nan
        tau2 = max(0.0, (q - df_q) / c) if c > 0 else 0.0
        wr = 1.0 / (v + tau2)
        pooled = float(np.sum(wr * y) / np.sum(wr))
        se = float(np.sqrt(1.0 / np.sum(wr)))
        z = pooled / se if se > 0 else np.nan
        p = float(2.0 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
        signs = np.sign(y)
        positive = int(np.sum(signs > 0))
        negative = int(np.sum(signs < 0))
        k = len(y)
        records.append({
            "modality": modality, "gene": gene, "k": k,
            "right_high_cohorts": positive, "left_high_cohorts": negative,
            "same_direction": bool(positive == k or negative == k),
            "pooled_hedges_g": pooled, "se": se,
            "ci_low": pooled - 1.96 * se, "ci_high": pooled + 1.96 * se,
            "meta_p": p, "tau2": tau2,
            "i2_percent": max(0.0, (q - df_q) / q * 100.0) if q > 0 else 0.0,
            "nominal_welch_p_lt_0.05": int(np.sum(d["welch_p"].to_numpy(float) < 0.05)),
        })
    out = pd.DataFrame(records)
    if out.empty:
        return out
    out["meta_fdr"] = bh_adjust(out["meta_p"].to_numpy())
    out["stable_right_high"] = (out["k"] >= 3) & (out["right_high_cohorts"] == out["k"]) & (out["meta_fdr"] < 0.05)
    out["stable_left_high"] = (out["k"] >= 3) & (out["left_high_cohorts"] == out["k"]) & (out["meta_fdr"] < 0.05)
    return out.sort_values(["meta_fdr", "gene"], na_position="last").reset_index(drop=True)


def cell_id_variants(value: object) -> list[str]:
    """Return common exact barcode variants used by GEO 10x exports."""
    base = str(value).strip().strip('"')
    variants = [base]
    if base.endswith("-1"):
        variants.append(base[:-2])
    else:
        variants.append(base + "-1")
    return list(dict.fromkeys(variants))


def read_dense_pseudobulk(path: Path, delimiter: str, selected_cell_ids: list[str], patient_labels: list[str], cohort: str, source_note: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if not path.exists():
        raise FileNotFoundError(path)
    if len(selected_cell_ids) != len(patient_labels):
        raise ValueError("selected_cell_ids and patient_labels must have equal length")
    patient_names = list(dict.fromkeys(patient_labels))
    patient_index = {name: i for i, name in enumerate(patient_names)}
    cell_to_group: dict[int, int] = {}
    total_by_patient = np.zeros(len(patient_names), dtype=float)
    rows: list[tuple[str, np.ndarray]] = []
    with gzip.open(path, "rt", errors="replace") as handle:
        header = [item.strip().strip('"') for item in handle.readline().rstrip("\r\n").split(delimiter)]
        exact_positions = {value: index for index, value in enumerate(header)}
        normalized_positions: dict[str, list[int]] = {}
        for index, value in enumerate(header):
            normalized = value[:-2] if value.endswith("-1") else value
            normalized_positions.setdefault(normalized, []).append(index)
        positions = []
        for cell_id, patient in zip(selected_cell_ids, patient_labels, strict=True):
            position = None
            for candidate in cell_id_variants(cell_id):
                if candidate in exact_positions:
                    position = exact_positions[candidate]
                    break
            if position is None:
                normalized = str(cell_id).strip().strip('"')
                if normalized.endswith("-1"):
                    normalized = normalized[:-2]
                candidates = normalized_positions.get(normalized, [])
                if len(candidates) == 1:
                    position = candidates[0]
            if position is None:
                raise ValueError(f"{cohort}: cell {cell_id} absent from matrix header")
            positions.append(position)
            cell_to_group[position] = patient_index[patient]
        positions_array = np.asarray(positions, dtype=int)
        group_for_position = np.asarray([cell_to_group[position] for position in positions], dtype=int)
        for line in handle:
            fields = line.rstrip("\r\n").split(delimiter)
            if not fields:
                continue
            gene = clean_symbol(fields[0])
            try:
                values = np.fromiter((float(fields[position] or 0.0) for position in positions_array), dtype=float, count=len(positions_array))
            except (ValueError, IndexError) as exc:
                raise ValueError(f"{cohort}: malformed numeric row for {fields[0]}") from exc
            for group_i in range(len(patient_names)):
                group_values = values[group_for_position == group_i]
                total_by_patient[group_i] += float(group_values.sum())
            if gene is not None:
                sums = np.zeros(len(patient_names), dtype=float)
                for group_i in range(len(patient_names)):
                    sums[group_i] = float(values[group_for_position == group_i].sum())
                rows.append((gene, sums))
    if not rows:
        raise ValueError(f"{cohort}: no genes were parsed from {path}")
    gene_names = [row[0] for row in rows]
    matrix = np.vstack([row[1] for row in rows])
    counts = pd.DataFrame(matrix, index=gene_names, columns=patient_names).groupby(level=0, sort=False).sum()
    cpm = counts.divide(pd.Series(total_by_patient, index=patient_names), axis=1) * 1e6
    expr = np.log2(cpm + 1.0).T
    audit = pd.DataFrame({
        "cohort": cohort, "sample_id": patient_names, "patient": patient_names,
        "side": ["unknown"] * len(patient_names), "n_selected_cells": [patient_labels.count(p) for p in patient_names],
        "source_note": source_note,
    }).set_index("sample_id")
    return expr, audit, {"source": str(path), "n_patients": len(patient_names), "n_genes": expr.shape[1], "normalization": "summed raw UMI / summed selected-cell UMI * 1e6; log2(CPM+1)"}


def load_gse200997() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    ann = pd.read_csv(GSE200997_ANNOT).rename(columns={"Unnamed: 0": "cell_id"})
    ann["cell_id"] = ann["cell_id"].astype(str)
    selection = pd.read_csv(GSE200997_SELECTION)
    selection["cell_id"] = selection["cell_id"].astype(str)
    selection = selection.set_index("cell_id").reindex(ann["cell_id"]).reset_index()
    malignant = selection["malignant_epithelial"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    selected = ann["Condition"].astype(str).str.lower().eq("tumor").to_numpy() & malignant
    selected_ids = ann.loc[selected, "cell_id"].tolist()
    patients = ann.loc[selected, "samples"].astype(str).str.extract(r"(T_cac\d+)", expand=False).fillna(ann.loc[selected, "samples"].astype(str)).tolist()
    expr, audit, info = read_dense_pseudobulk(GSE200997_MATRIX, ",", selected_ids, patients, "GSE200997", "Tumor cells marked malignant_epithelial=True in existing marker-defined selection")
    side_map = ann.loc[selected, ["samples", "Location"]].drop_duplicates().copy()
    side_map["side"] = side_map["Location"].astype(str).str.lower()
    side_map["patient"] = side_map["samples"].astype(str).str.extract(r"(T_cac\d+)", expand=False).fillna(side_map["samples"].astype(str))
    audit["side"] = audit["patient"].map(side_map.set_index("patient")["side"])
    if audit["side"].value_counts().to_dict() != {"left": 8, "right": 8}:
        raise ValueError(f"GSE200997 expected 8/8 patient sides, got {audit['side'].value_counts().to_dict()}")
    audit["selection"] = "existing_marker_defined_malignant_epithelial"
    info.update({"selection": "existing marker-defined malignant epithelial cells", "n_selected_cells": int(selected.sum())})
    return expr, audit, info


def load_official_tumor_epithelial(annot_path: Path, matrix_path: Path, soft_path: Path, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    ann = pd.read_csv(annot_path, sep="\t").rename(columns={"Index": "cell_id"})
    ann["cell_id"] = ann["cell_id"].astype(str)
    if cohort == "GSE132465":
        selected = ann["Class"].astype(str).str.lower().eq("tumor") & ann["Cell_type"].astype(str).str.lower().eq("epithelial cells")
        sample_col = "Sample"
        patient_col = "Patient"
    else:
        selected = ann["Class"].astype(str).str.lower().eq("tumor") & ann["Cell_type"].astype(str).str.lower().eq("epithelial cells")
        sample_col = "Sample"
        patient_col = "Patient"
    selected_ann = ann.loc[selected].copy()
    soft = parse_soft_samples(soft_path)
    side_rows = []
    for sample, group in selected_ann.groupby(sample_col, sort=True):
        soft_record = next((record for sid, record in soft.items() if (text_values(record, "title") or [""])[0] == str(sample)), None)
        if soft_record is None:
            soft_record = next((record for sid, record in soft.items() if sid == str(sample)), None)
        if soft_record is None:
            raise ValueError(f"{cohort}: no GEO sample metadata found for {sample}")
        text = sample_level_text(soft_record)
        side, side_match = classify_side(text)
        side_rows.append({"sample": str(sample), "side": side, "side_match": side_match, "sample_level_text": text})
    side_df = pd.DataFrame(side_rows)
    selected_ids = selected_ann["cell_id"].tolist()
    patients = selected_ann[patient_col].astype(str).tolist()
    expr, audit, info = read_dense_pseudobulk(matrix_path, "\t", selected_ids, patients, cohort, "Official GEO Class=Tumor and Cell_type=Epithelial cells")
    sample_to_side = selected_ann[[sample_col, patient_col]].drop_duplicates().merge(side_df, left_on=sample_col, right_on="sample", how="left")
    patient_side = sample_to_side.groupby(patient_col, sort=True)["side"].first()
    audit["side"] = audit["patient"].map(patient_side)
    audit["sample"] = audit["patient"].map(sample_to_side.drop_duplicates(patient_col).set_index(patient_col)[sample_col])
    audit["side_match"] = audit["patient"].map(sample_to_side.drop_duplicates(patient_col).set_index(patient_col)["side_match"])
    audit["selection"] = "official_tumor_epithelial"
    info.update({"selection": "official GEO tumor epithelial annotation; malignancy inferred from tumor origin", "n_selected_cells": int(selected.sum()), "side_map": side_df.to_dict("records")})
    return expr, audit, info


def load_gse188711() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if not GSE188711_META.exists():
        raise FileNotFoundError(GSE188711_META)
    meta = pd.read_csv(GSE188711_META, index_col=0)
    meta["leiden"] = meta["leiden"].astype(str)
    selected_cells = meta[meta["leiden"].eq("4")]
    patient_expr: list[pd.DataFrame] = []
    audit_rows: list[dict[str, object]] = []
    for sample_id, spec in GSE188711_SAMPLES.items():
        matrix_path = next(GSE188711_RAW.glob(f"{sample_id}_matrix_{spec['tag']}.mtx.gz"))
        features_path = next(GSE188711_RAW.glob(f"{sample_id}_features_{spec['tag']}.tsv.gz"))
        barcodes_path = next(GSE188711_RAW.glob(f"{sample_id}_barcodes_{spec['tag']}.tsv.gz"))
        with gzip.open(matrix_path, "rt") as handle:
            matrix = mmread(handle).tocsr()
        features = pd.read_csv(features_path, sep="\t", header=None, compression="gzip")
        barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, compression="gzip")[0].astype(str).tolist()
        symbols = features.iloc[:, 1].map(clean_symbol) if features.shape[1] > 1 else features.iloc[:, 0].map(clean_symbol)
        prefixed = [f"{sample_id}_{barcode}" for barcode in barcodes]
        selected_ids = set(selected_cells.index[selected_cells["sample"].eq(sample_id)].astype(str))
        selected_columns = [i for i, cell_id in enumerate(prefixed) if cell_id in selected_ids]
        if not selected_columns:
            raise ValueError(f"GSE188711: no cluster-4 cells found for {sample_id}")
        counts = np.asarray(matrix[:, selected_columns].sum(axis=1)).ravel()
        total = float(counts.sum())
        gene_counts = pd.Series(counts, index=symbols).dropna().groupby(level=0).sum()
        expr_row = np.log2(gene_counts / total * 1e6 + 1.0).to_frame().T
        expr_row.index = [sample_id]
        patient_expr.append(expr_row)
        audit_rows.append({
            "cohort": "GSE188711", "sample_id": sample_id, "patient": sample_id,
            "side": spec["side"], "patient_label": spec["label"],
            "n_selected_cells": len(selected_columns), "selection": "existing_Leiden_cluster_4",
            "source_note": "Primary putative malignant epithelial cluster from prior local QC/marker analysis",
        })
    expr = pd.concat(patient_expr, axis=0, sort=True).fillna(0.0)
    audit = pd.DataFrame(audit_rows).set_index("sample_id")
    info = {
        "source": str(GSE188711_RAW), "n_patients": len(audit), "n_genes": expr.shape[1],
        "n_selected_cells": int(audit["n_selected_cells"].sum()),
        "selection": "existing primary putative malignant epithelial Leiden cluster 4",
        "sensitivity_warning": "3 left and 3 right patients; directional sensitivity only",
        "normalization": "summed raw UMI / summed selected-cell UMI * 1e6; log2(CPM+1)",
    }
    return expr, audit, info


def write_intersection(bulk_meta: pd.DataFrame, sc_meta: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    bulk_flags = bulk_meta[["gene", "stable_right_high", "stable_left_high", "pooled_hedges_g", "meta_fdr"]].rename(columns={
        "stable_right_high": "bulk_stable_right_high", "stable_left_high": "bulk_stable_left_high",
        "pooled_hedges_g": "bulk_pooled_hedges_g", "meta_fdr": "bulk_meta_fdr",
    })
    sc_flags = sc_meta[["gene", "stable_right_high", "stable_left_high", "pooled_hedges_g", "meta_fdr"]].rename(columns={
        "stable_right_high": "sc_stable_right_high", "stable_left_high": "sc_stable_left_high",
        "pooled_hedges_g": "sc_pooled_hedges_g", "meta_fdr": "sc_meta_fdr",
    })
    merged = bulk_flags.merge(sc_flags, on="gene", how="outer")
    merged["stable_right_high_bulk_and_sc"] = merged["bulk_stable_right_high"].astype("boolean").fillna(False) & merged["sc_stable_right_high"].astype("boolean").fillna(False)
    merged["stable_left_high_bulk_and_sc"] = merged["bulk_stable_left_high"].astype("boolean").fillna(False) & merged["sc_stable_left_high"].astype("boolean").fillna(False)
    merged.to_csv(out_dir / "crc_sidedness_stable_deg_bulk_sc_intersection.csv", index=False)
    strict = merged[merged["stable_right_high_bulk_and_sc"] | merged["stable_left_high_bulk_and_sc"]].copy()
    strict.to_csv(out_dir / "crc_sidedness_stable_deg_strict_intersection.csv", index=False)
    return merged


def write_report(summary: dict[str, object], bulk_meta: pd.DataFrame, sc_meta: pd.DataFrame, intersection: pd.DataFrame, out_dir: Path) -> None:
    bulk_right = bulk_meta[bulk_meta["stable_right_high"]].sort_values("meta_fdr")
    bulk_left = bulk_meta[bulk_meta["stable_left_high"]].sort_values("meta_fdr")
    sc_right = sc_meta[sc_meta["stable_right_high"]].sort_values("meta_fdr")
    sc_left = sc_meta[sc_meta["stable_left_high"]].sort_values("meta_fdr")
    lines = [
        "# CRC sidedness stable DEG discovery",
        "",
        "## Question",
        "",
        "Which genes show reproducible right-versus-left differences across bulk CRC cohorts and persist in tumor epithelial single-cell pseudobulk?",
        "",
        "## Frozen design",
        "",
        "- Bulk: TCGA-COAD, GSE39582, GSE103479 and GSE41258.",
        "- Single cell: GSE200997, GSE132465, GSE144735 and GSE188711.",
        "- Bulk arrays/RNA-seq were analyzed within cohort; no cross-platform expression matrix was concatenated.",
        "- Single-cell raw UMI counts were summed within patient and tumor-epithelial compartment before testing. No cell-level p-values were used.",
        "- Differential expression is Welch right-versus-left testing on the native normalized/log expression scale or patient pseudobulk log2(CPM+1).",
        "- Meta-analysis uses Hedges-g standardized effects with a DerSimonian–Laird random-effects model.",
        "- Strict stable right-high/left-high: at least 3 contributing cohorts, all effects in the same direction, random-effects BH-FDR < 0.05.",
        "",
        "## Counts and strict lists",
        "",
        f"- Bulk cohort-level tables: {len(BULK_COHORTS)}; strict right-high genes: **{len(bulk_right)}**; strict left-high genes: **{len(bulk_left)}**.",
        f"- Single-cell cohort-level tables: {len(summary['sc_cohorts'])}; strict right-high genes: **{len(sc_right)}**; strict left-high genes: **{len(sc_left)}**.",
        f"- Strict bulk–single-cell overlap: right-high **{int(intersection['stable_right_high_bulk_and_sc'].sum())}**, left-high **{int(intersection['stable_left_high_bulk_and_sc'].sum())}**.",
        "",
        "The strict lists are intentionally conservative. The full meta tables retain effect direction, heterogeneity, cohort count and FDR for genes that are biologically interesting but underpowered.",
        "",
        "## Bulk cohort audit",
        "",
        "| Cohort | n right | n left | genes | source |",
        "|---|---:|---:|---:|---|",
    ]
    for row in summary["bulk_summary"]:
        lines.append(f"| {row['cohort']} | {row['n_right']} | {row['n_left']} | {row['n_genes']} | {row['source']} |")
    lines += [
        "",
        "## Single-cell cohort audit",
        "",
        "| Cohort | n right | n left | selected cells | selection | warning |",
        "|---|---:|---:|---:|---|---|",
    ]
    for row in summary["sc_summary"]:
        lines.append(f"| {row['cohort']} | {row['n_right']} | {row['n_left']} | {row['n_selected_cells']} | {row['selection']} | {row.get('sensitivity_warning', '')} |")
    lines += [
        "",
        "## Interpretation guardrails",
        "",
        "This is a phenotype-first sidedness screen. It does not by itself establish causality, cell-state mechanism, AA/PUFA flux, ferroptosis, treatment response or therapeutic dependency. GSE132465/GSE144735 tumor epithelial cells are annotation-based tumor-epithelial proxies; GSE188711 is a low-powered 3-versus-3 sensitivity cohort.",
        "",
        "## Outputs",
        "",
        "- `crc_sidedness_bulk_de_<cohort>.csv`: all gene-level within-cohort statistics.",
        "- `crc_sidedness_bulk_random_effects_meta.csv`: bulk random-effects meta-analysis.",
        "- `crc_sidedness_sc_de_<cohort>.csv`: all gene-level single-cell pseudobulk statistics.",
        "- `crc_sidedness_sc_random_effects_meta.csv`: single-cell random-effects meta-analysis.",
        "- `crc_sidedness_stable_deg_bulk_sc_intersection.csv`: all genes with bulk/sc flags.",
        "- `crc_sidedness_stable_deg_strict_intersection.csv`: strict bulk–single-cell overlap only.",
        "- `crc_sidedness_bulk_side_audit.csv` and `crc_sidedness_sc_patient_audit.csv`: sample/patient sidedness and selection audit.",
    ]
    (out_dir / "crc_sidedness_stable_deg_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    bulk_expr: dict[str, tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]] = {}
    bulk_expr["GSE39582"] = load_generic_geo("GSE39582")
    bulk_expr["GSE41258"] = load_generic_geo("GSE41258")
    bulk_expr["GSE103479"] = load_gse103479()
    bulk_expr["TCGA-COAD"] = load_tcga()

    bulk_tables: list[pd.DataFrame] = []
    bulk_audits: list[pd.DataFrame] = []
    bulk_summary: list[dict[str, object]] = []
    for cohort in BULK_COHORTS:
        expr, audit, info = bulk_expr[cohort]
        table = differential(expr, audit, cohort, "bulk")
        table.to_csv(args.out_dir / f"crc_sidedness_bulk_de_{cohort.replace('-', '_')}.csv", index=False)
        bulk_tables.append(table)
        audit_out = audit.reset_index()
        bulk_audits.append(audit_out)
        sided = audit[audit["side"].isin(["right", "left"])]
        bulk_summary.append({
            "cohort": cohort, "n_right": int((sided["side"] == "right").sum()),
            "n_left": int((sided["side"] == "left").sum()), "n_genes": int(expr.shape[1]),
            "source": info.get("source", ""), "platform": info.get("platform", ""),
        })
    bulk_meta = random_effects_meta(bulk_tables, "bulk")
    bulk_meta.to_csv(args.out_dir / "crc_sidedness_bulk_random_effects_meta.csv", index=False)
    pd.concat(bulk_audits, ignore_index=True).to_csv(args.out_dir / "crc_sidedness_bulk_side_audit.csv", index=False)

    sc_loaders = {
        "GSE200997": load_gse200997,
        "GSE132465": lambda: load_official_tumor_epithelial(GSE132465_ANNOT, GSE132465_MATRIX, GSE132465_SOFT, "GSE132465"),
        "GSE144735": lambda: load_official_tumor_epithelial(GSE144735_ANNOT, GSE144735_MATRIX, GSE144735_SOFT, "GSE144735"),
        "GSE188711": load_gse188711,
    }
    if args.skip_gse188711:
        sc_loaders.pop("GSE188711")
    sc_tables: list[pd.DataFrame] = []
    sc_audits: list[pd.DataFrame] = []
    sc_summary: list[dict[str, object]] = []
    sc_infos: dict[str, object] = {}
    for cohort, loader in sc_loaders.items():
        expr, audit, info = loader()
        table = differential(expr, audit, cohort, "single_cell_patient_pseudobulk")
        table.to_csv(args.out_dir / f"crc_sidedness_sc_de_{cohort}.csv", index=False)
        expr.to_csv(args.out_dir / f"crc_sidedness_sc_patient_logcpm_{cohort}.csv")
        sc_tables.append(table)
        sc_audits.append(audit.reset_index())
        sc_infos[cohort] = info
        sided = audit[audit["side"].isin(["right", "left"])]
        sc_summary.append({
            "cohort": cohort, "n_right": int((sided["side"] == "right").sum()),
            "n_left": int((sided["side"] == "left").sum()),
            "n_selected_cells": int(audit.get("n_selected_cells", pd.Series(dtype=float)).sum()) if "n_selected_cells" in audit else np.nan,
            "n_genes": int(expr.shape[1]), "selection": info.get("selection", ""),
            "sensitivity_warning": info.get("sensitivity_warning", ""),
        })
    sc_meta = random_effects_meta(sc_tables, "single_cell_patient_pseudobulk")
    sc_meta.to_csv(args.out_dir / "crc_sidedness_sc_random_effects_meta.csv", index=False)
    pd.concat(sc_audits, ignore_index=True).to_csv(args.out_dir / "crc_sidedness_sc_patient_audit.csv", index=False)
    (args.out_dir / "crc_sidedness_sc_loader_manifest.json").write_text(json.dumps(sc_infos, indent=2, ensure_ascii=False, default=str), encoding="utf-8")

    intersection = write_intersection(bulk_meta, sc_meta, args.out_dir)
    summary = {"bulk_summary": bulk_summary, "sc_summary": sc_summary, "sc_cohorts": list(sc_loaders), "inputs": {"bulk": BULK_COHORTS, "single_cell": list(sc_loaders)}}
    (args.out_dir / "crc_sidedness_stable_deg_manifest.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    write_report(summary, bulk_meta, sc_meta, intersection, args.out_dir)
    print(json.dumps({
        "bulk_strict_right_high": int(bulk_meta["stable_right_high"].sum()),
        "bulk_strict_left_high": int(bulk_meta["stable_left_high"].sum()),
        "sc_strict_right_high": int(sc_meta["stable_right_high"].sum()),
        "sc_strict_left_high": int(sc_meta["stable_left_high"].sum()),
        "bulk_sc_strict_right_overlap": int(intersection["stable_right_high_bulk_and_sc"].sum()),
        "bulk_sc_strict_left_overlap": int(intersection["stable_left_high_bulk_and_sc"].sum()),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
