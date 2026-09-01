"""Frozen multi-cohort validation of CRC sidedness, PUFA incorporation and antioxidant buffering.

Frozen cohorts
--------------
1. GSE39582
2. TCGA-COAD
3. GSE41258
4. GSE4554
5. GSE75316

Primary readouts
----------------
- PUFA incorporation score: ACSL4 + LPCAT3 + AGPAT3
- Antioxidant buffering score: SLC7A11 + GPX4 + AIFM2 + GCH1 + DHODH
- SLC7A11 expression
- GCH1 expression

Primary sidedness rule
----------------------
Right: cecum / caecum / ascending / hepatic flexure / transverse / proximal
Left: descending / splenic flexure / sigmoid / rectosigmoid / distal
Rectum-only and ambiguous generic colon records are excluded from the primary comparison.

The script deliberately records the raw metadata text and the matched location token for audit.
GEO cohorts are downloaded with GEOparse when needed; TCGA uses project-local cached inputs.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm, ttest_ind

try:
    import GEOparse  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "GEOparse is required for the frozen GEO cohorts. Install with: pip install GEOparse"
    ) from exc

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
DATA = WORK / "data"
OUT = ROOT / "outputs"
GEO_CACHE = DATA / "geo_frozen_sidedness"

FROZEN_GEO = ["GSE39582", "GSE41258", "GSE4554", "GSE75316"]
FROZEN_COHORTS = FROZEN_GEO + ["TCGA-COAD"]

PUFA_GENES = ["ACSL4", "LPCAT3", "AGPAT3"]
BUFFER_GENES = ["SLC7A11", "GPX4", "AIFM2", "GCH1", "DHODH"]
TARGET_GENES = sorted(set(PUFA_GENES + BUFFER_GENES))

TCGA_CASES = DATA / "tcga_coad_cases.json"
TCGA_EXPR_CANDIDATES = [
    DATA / "tcga_ferroptosis_three_layer_expression_uqfpkm.tsv",
    DATA / "tcga_aa_defense_candidate_expression_uqfpkm.tsv",
    DATA / "tcga_coad_gene_expression_uqfpkm.tsv",
]

RIGHT_PATTERNS = [
    r"\bcecum\b", r"\bcaecum\b", r"\bcecal\b", r"\bcaecal\b",
    r"\bascending\b", r"hepatic\s+flexure", r"\btransverse\b", r"\bproximal\b",
    r"right\s+(?:sided\s+)?colon", r"right[- ]sided",
]
LEFT_PATTERNS = [
    r"\bdescending\b", r"splenic\s+flexure", r"\bsigmoid\b", r"\brectosigmoid\b",
    r"\bdistal\b", r"left\s+(?:sided\s+)?colon", r"left[- ]sided",
]
RECTUM_PATTERNS = [r"\brectum\b", r"\brectal\b"]
NORMAL_PATTERNS = [
    r"\bnormal\b", r"adjacent\s+normal", r"non[- ]tumou?r", r"mucosa",
]
TUMOR_PATTERNS = [
    r"tumou?r", r"adenocarcinoma", r"carcinoma", r"primary",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=OUT)
    p.add_argument("--geo-cache", type=Path, default=GEO_CACHE)
    p.add_argument("--skip-download", action="store_true")
    return p.parse_args()


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / sd


def bh_adjust(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q.tolist()
    idx = np.where(ok)[0]
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    restored = np.empty_like(adj)
    restored[order] = np.minimum(adj, 1.0)
    q[idx] = restored
    return q.tolist()


def flatten_metadata(metadata: dict[str, object]) -> str:
    parts: list[str] = []
    for key, value in metadata.items():
        if isinstance(value, (list, tuple)):
            vals = [str(v) for v in value]
        else:
            vals = [str(value)]
        parts.extend([f"{key}: {v}" for v in vals])
    return " | ".join(parts)


def classify_side(text: str) -> tuple[str, str]:
    low = str(text).lower()
    right_hits = [pat for pat in RIGHT_PATTERNS if re.search(pat, low)]
    left_hits = [pat for pat in LEFT_PATTERNS if re.search(pat, low)]
    rectal_hits = [pat for pat in RECTUM_PATTERNS if re.search(pat, low)]
    if right_hits and not left_hits:
        return "right", right_hits[0]
    if left_hits and not right_hits:
        # rectosigmoid is retained as left by the frozen rule.
        if rectal_hits and not re.search(r"rectosigmoid", low):
            return "exclude_rectum", rectal_hits[0]
        return "left", left_hits[0]
    if rectal_hits and not right_hits and not left_hits:
        return "exclude_rectum", rectal_hits[0]
    if right_hits and left_hits:
        return "ambiguous", f"right={right_hits[0]};left={left_hits[0]}"
    return "unknown", ""


def looks_tumor(text: str) -> bool:
    low = str(text).lower()
    if any(re.search(p, low) for p in NORMAL_PATTERNS):
        # Normal wins unless tumor wording is clearly present and location text is tumor-specific.
        if not any(re.search(p, low) for p in TUMOR_PATTERNS):
            return False
    return True


def gene_symbol_column(gpl_table: pd.DataFrame) -> str:
    candidates = [
        "Gene Symbol", "GENE_SYMBOL", "Gene symbol", "gene_assignment",
        "Gene Symbol /// Gene Symbol", "Symbol", "SYMBOL",
    ]
    for c in candidates:
        if c in gpl_table.columns:
            return c
    for c in gpl_table.columns:
        low = str(c).lower()
        if "gene" in low and "symbol" in low:
            return c
    raise ValueError(f"Could not identify gene-symbol column. Columns={list(gpl_table.columns)[:30]}")


def normalize_symbol_cell(value: object) -> list[str]:
    text = str(value)
    if not text or text.lower() == "nan":
        return []
    tokens = re.split(r"\s*///\s*|\s*//\s*|\s*;\s*|\s*,\s*|\s+", text)
    clean = []
    for token in tokens:
        token = token.strip().upper()
        if re.fullmatch(r"[A-Z0-9._-]+", token):
            clean.append(token)
    return clean


def collapse_probes_to_genes(expr_probe: pd.DataFrame, gpl_table: pd.DataFrame) -> pd.DataFrame:
    gpl = gpl_table.copy()
    if "ID" in gpl.columns:
        probe_ids = gpl["ID"].astype(str)
    else:
        probe_ids = gpl.index.astype(str)
    sym_col = gene_symbol_column(gpl)
    mapping: dict[str, list[str]] = {g: [] for g in TARGET_GENES}
    for probe, raw in zip(probe_ids, gpl[sym_col], strict=False):
        symbols = set(normalize_symbol_cell(raw))
        for gene in TARGET_GENES:
            if gene in symbols:
                mapping[gene].append(str(probe))
    out = pd.DataFrame(index=expr_probe.columns)
    for gene, probes in mapping.items():
        available = [p for p in probes if p in expr_probe.index]
        if not available:
            out[gene] = np.nan
            continue
        block = expr_probe.loc[available].apply(pd.to_numeric, errors="coerce")
        # Median across probes is robust to a single discordant probe.
        out[gene] = block.median(axis=0, skipna=True)
    return out


def load_geo(accession: str, cache: Path, skip_download: bool) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    cache.mkdir(parents=True, exist_ok=True)
    if skip_download:
        gse = GEOparse.get_GEO(geo=accession, destdir=str(cache), silent=True)
    else:
        gse = GEOparse.get_GEO(geo=accession, destdir=str(cache), silent=True)

    # Expression by platform, then concatenate if a cohort uses multiple GPLs.
    pieces: list[pd.DataFrame] = []
    sample_rows: list[dict[str, object]] = []
    platform_ids = sorted({gsm.metadata.get("platform_id", [""])[0] for gsm in gse.gsms.values()})

    for gpl_id in platform_ids:
        if not gpl_id or gpl_id not in gse.gpls:
            continue
        gsm_ids = [
            gsm_name for gsm_name, gsm in gse.gsms.items()
            if gsm.metadata.get("platform_id", [""])[0] == gpl_id
        ]
        if not gsm_ids:
            continue
        expr_probe = pd.DataFrame({
            gsm_name: gse.gsms[gsm_name].table.set_index("ID_REF")["VALUE"]
            for gsm_name in gsm_ids
            if {"ID_REF", "VALUE"}.issubset(gse.gsms[gsm_name].table.columns)
        })
        if expr_probe.empty:
            continue
        gene_expr = collapse_probes_to_genes(expr_probe, gse.gpls[gpl_id].table)
        pieces.append(gene_expr)

    if not pieces:
        raise ValueError(f"{accession}: no usable expression matrices")
    expr = pd.concat(pieces, axis=0)
    expr = expr[~expr.index.duplicated(keep="first")]

    for gsm_name, gsm in gse.gsms.items():
        text = flatten_metadata(gsm.metadata)
        side, matched = classify_side(text)
        sample_rows.append({
            "sample_id": gsm_name,
            "cohort": accession,
            "side": side,
            "side_match": matched,
            "is_tumor_like": looks_tumor(text),
            "raw_metadata_text": text,
            "title": " | ".join(map(str, gsm.metadata.get("title", []))),
            "source_name": " | ".join(map(str, gsm.metadata.get("source_name_ch1", []))),
            "platform_id": " | ".join(map(str, gsm.metadata.get("platform_id", []))),
        })
    meta = pd.DataFrame(sample_rows).set_index("sample_id")
    common = expr.index.intersection(meta.index)
    expr = expr.loc[common]
    meta = meta.loc[common]
    info = {
        "accession": accession,
        "platforms": platform_ids,
        "n_expression_samples": int(len(expr)),
        "side_counts_all": meta["side"].value_counts(dropna=False).to_dict(),
    }
    return expr, meta, info


def load_tcga() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if not TCGA_CASES.exists():
        raise FileNotFoundError(f"Missing {TCGA_CASES}")
    expr_path = next((p for p in TCGA_EXPR_CANDIDATES if p.exists()), None)
    if expr_path is None:
        raise FileNotFoundError("No TCGA expression cache found among expected files")

    expr_raw = pd.read_csv(expr_path, sep="\t", index_col=0)
    # Support either genes x cases or cases x genes.
    gene_hits_index = sum(g in set(map(str, expr_raw.index)) for g in TARGET_GENES)
    gene_hits_cols = sum(g in set(map(str, expr_raw.columns)) for g in TARGET_GENES)
    if gene_hits_index >= gene_hits_cols:
        expr = expr_raw.loc[expr_raw.index.intersection(TARGET_GENES)].T.copy()
    else:
        expr = expr_raw.loc[:, expr_raw.columns.intersection(TARGET_GENES)].copy()
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce")

    cases = json.loads(TCGA_CASES.read_text(encoding="utf-8"))
    if isinstance(cases, dict) and "data" in cases:
        cases = cases["data"]
    rows: list[dict[str, object]] = []
    for case in cases:
        case_id = str(case.get("case_id") or case.get("submitter_id") or "")
        text = json.dumps(case, ensure_ascii=False)
        side, matched = classify_side(text)
        rows.append({
            "sample_id": case_id,
            "cohort": "TCGA-COAD",
            "side": side,
            "side_match": matched,
            "is_tumor_like": True,
            "raw_metadata_text": text,
            "title": case_id,
            "source_name": "TCGA-COAD",
            "platform_id": "RNA-seq",
        })
    meta = pd.DataFrame(rows).set_index("sample_id")
    common = expr.index.intersection(meta.index)
    expr = expr.loc[common]
    meta = meta.loc[common]
    info = {
        "accession": "TCGA-COAD",
        "expression_path": str(expr_path),
        "n_expression_samples": int(len(expr)),
        "side_counts_all": meta["side"].value_counts(dropna=False).to_dict(),
        "available_genes": sorted(set(expr.columns).intersection(TARGET_GENES)),
    }
    return expr, meta, info


def add_scores(expr: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    x = expr.copy().apply(pd.to_numeric, errors="coerce")
    # Rank-preserving log transform only when values look unlogged.
    finite = x.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size and np.nanpercentile(finite, 95) > 50:
        x = np.log2(np.clip(x, a_min=0, a_max=None) + 1.0)
    z = pd.DataFrame(index=x.index)
    for gene in TARGET_GENES:
        z[gene] = zscore(x[gene]) if gene in x.columns else np.nan
    out = meta.copy()
    for gene in TARGET_GENES:
        out[gene] = z[gene]
    pufa_avail = [g for g in PUFA_GENES if g in z.columns and z[g].notna().any()]
    buffer_avail = [g for g in BUFFER_GENES if g in z.columns and z[g].notna().any()]
    out["pufa_available_genes"] = z[pufa_avail].notna().sum(axis=1) if pufa_avail else 0
    out["buffer_available_genes"] = z[buffer_avail].notna().sum(axis=1) if buffer_avail else 0
    out["pufa_incorporation_score"] = z[pufa_avail].mean(axis=1, skipna=True) if pufa_avail else np.nan
    out["antioxidant_buffering_score"] = z[buffer_avail].mean(axis=1, skipna=True) if buffer_avail else np.nan
    out.loc[out["pufa_available_genes"] < 2, "pufa_incorporation_score"] = np.nan
    out.loc[out["buffer_available_genes"] < 4, "antioxidant_buffering_score"] = np.nan
    return out


def compare_metric(data: pd.DataFrame, metric: str) -> dict[str, object]:
    d = data[data["is_tumor_like"] & data["side"].isin(["right", "left"])].copy()
    right = pd.to_numeric(d.loc[d["side"].eq("right"), metric], errors="coerce").dropna()
    left = pd.to_numeric(d.loc[d["side"].eq("left"), metric], errors="coerce").dropna()
    if len(right) < 2 or len(left) < 2:
        return {
            "metric": metric, "n_right": int(len(right)), "n_left": int(len(left)),
            "right_mean": np.nan, "left_mean": np.nan, "right_minus_left": np.nan,
            "smd": np.nan, "se_smd": np.nan, "welch_p": np.nan, "mannwhitney_p": np.nan,
        }
    vr, vl = right.var(ddof=1), left.var(ddof=1)
    pooled = math.sqrt(((len(right)-1)*vr + (len(left)-1)*vl) / (len(right)+len(left)-2))
    smd = (right.mean() - left.mean()) / pooled if pooled > 0 else np.nan
    # Hedges small-sample correction and approximate SE.
    n = len(right) + len(left)
    j = 1 - 3 / (4 * n - 9) if n > 2 else 1.0
    g = smd * j if np.isfinite(smd) else np.nan
    se = math.sqrt(n/(len(right)*len(left)) + (g*g)/(2*(n-2))) if n > 2 and np.isfinite(g) else np.nan
    return {
        "metric": metric,
        "n_right": int(len(right)), "n_left": int(len(left)),
        "right_mean": float(right.mean()), "left_mean": float(left.mean()),
        "right_minus_left": float(right.mean() - left.mean()),
        "smd": float(g), "se_smd": float(se),
        "welch_p": float(ttest_ind(right, left, equal_var=False).pvalue),
        "mannwhitney_p": float(mannwhitneyu(right, left, alternative="two-sided").pvalue),
    }


def random_effects_meta(df: pd.DataFrame, metric: str) -> dict[str, object]:
    d = df[(df["metric"] == metric) & df["smd"].notna() & df["se_smd"].notna()].copy()
    if len(d) < 2:
        return {"metric": metric, "k": int(len(d)), "estimable": False}
    yi = d["smd"].to_numpy(float)
    vi = np.square(d["se_smd"].to_numpy(float))
    wi = 1 / vi
    fixed = np.sum(wi * yi) / np.sum(wi)
    q = np.sum(wi * np.square(yi - fixed))
    df_q = len(yi) - 1
    c = np.sum(wi) - np.sum(np.square(wi)) / np.sum(wi)
    tau2 = max(0.0, (q - df_q) / c) if c > 0 else 0.0
    wr = 1 / (vi + tau2)
    pooled = np.sum(wr * yi) / np.sum(wr)
    se = math.sqrt(1 / np.sum(wr))
    z = pooled / se if se > 0 else np.nan
    p = 2 * norm.sf(abs(z)) if np.isfinite(z) else np.nan
    i2 = max(0.0, (q - df_q) / q) * 100 if q > 0 else 0.0
    return {
        "metric": metric,
        "k": int(len(d)),
        "estimable": True,
        "pooled_hedges_g": float(pooled),
        "se": float(se),
        "ci_low": float(pooled - 1.96*se),
        "ci_high": float(pooled + 1.96*se),
        "p": float(p),
        "tau2": float(tau2),
        "i2_percent": float(i2),
        "direction_right_higher_cohorts": int(np.sum(yi > 0)),
        "direction_right_lower_cohorts": int(np.sum(yi < 0)),
    }


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.geo_cache.mkdir(parents=True, exist_ok=True)

    cohort_scores: list[pd.DataFrame] = []
    manifests: dict[str, object] = {}
    audit_rows: list[pd.DataFrame] = []

    for cohort in FROZEN_GEO:
        expr, meta, info = load_geo(cohort, args.geo_cache, args.skip_download)
        scores = add_scores(expr, meta)
        cohort_scores.append(scores)
        manifests[cohort] = info
        audit_rows.append(scores[["cohort", "side", "side_match", "is_tumor_like", "title", "source_name", "platform_id", "raw_metadata_text"]])

    expr, meta, info = load_tcga()
    scores = add_scores(expr, meta)
    cohort_scores.append(scores)
    manifests["TCGA-COAD"] = info
    audit_rows.append(scores[["cohort", "side", "side_match", "is_tumor_like", "title", "source_name", "platform_id", "raw_metadata_text"]])

    combined = pd.concat(cohort_scores, axis=0, sort=False)
    audit = pd.concat(audit_rows, axis=0, sort=False)

    metrics = ["pufa_incorporation_score", "antioxidant_buffering_score", "SLC7A11", "GCH1"]
    rows: list[dict[str, object]] = []
    for cohort in FROZEN_COHORTS:
        d = combined[combined["cohort"].eq(cohort)]
        for metric in metrics:
            row = compare_metric(d, metric)
            row["cohort"] = cohort
            rows.append(row)
    stats = pd.DataFrame(rows)
    stats["welch_fdr_within_metric"] = np.nan
    for metric in metrics:
        mask = stats["metric"].eq(metric)
        stats.loc[mask, "welch_fdr_within_metric"] = bh_adjust(stats.loc[mask, "welch_p"].tolist())

    meta_rows = [random_effects_meta(stats, metric) for metric in metrics]
    meta_df = pd.DataFrame(meta_rows)

    combined.to_csv(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_scores.csv", index=True)
    audit.to_csv(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_side_audit.csv", index=True)
    stats.to_csv(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_stats.csv", index=False)
    meta_df.to_csv(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_meta.csv", index=False)

    manifest = {
        "frozen_cohorts": FROZEN_COHORTS,
        "primary_metrics": metrics,
        "pufa_genes": PUFA_GENES,
        "antioxidant_buffering_genes": BUFFER_GENES,
        "sidedness_rule": {
            "right": "cecum/caecum/ascending/hepatic flexure/transverse/proximal/right-sided colon",
            "left": "descending/splenic flexure/sigmoid/rectosigmoid/distal/left-sided colon",
            "rectum": "excluded from primary analysis unless explicitly rectosigmoid",
            "ambiguous_generic_colon": "excluded",
        },
        "minimum_genes": {"pufa_incorporation": "2/3", "antioxidant_buffering": "4/5"},
        "cohort_provenance": manifests,
    }
    (args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    lines = [
        "# Frozen CRC sidedness antioxidant-buffering validation",
        "",
        "Frozen cohorts: " + ", ".join(FROZEN_COHORTS),
        "",
        "Primary metrics: PUFA incorporation, antioxidant buffering, SLC7A11, GCH1.",
        "All scores are within-cohort z-score composites; they are transcriptional states, not direct biochemical measurements.",
        "",
        "## Cohort-level sidedness results",
        "",
        stats.to_markdown(index=False),
        "",
        "## Random-effects meta-analysis",
        "",
        meta_df.to_markdown(index=False),
        "",
        "## Audit rule",
        "",
        "Every GEO/TCGA sample keeps raw metadata text plus the matched sidedness token in `..._side_audit.csv`. Review ambiguous/unknown samples manually before final publication tables.",
    ]
    (args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    print(json.dumps({
        "frozen_cohorts": FROZEN_COHORTS,
        "cohort_stats": stats.to_dict("records"),
        "meta": meta_df.to_dict("records"),
        "outputs": {
            "scores": str(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_scores.csv"),
            "side_audit": str(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_side_audit.csv"),
            "stats": str(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_stats.csv"),
            "meta": str(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_meta.csv"),
            "manifest": str(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_manifest.json"),
            "report": str(args.out_dir / "frozen_crc_sidedness_antioxidant_multicohort_report.md"),
        },
    }, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
