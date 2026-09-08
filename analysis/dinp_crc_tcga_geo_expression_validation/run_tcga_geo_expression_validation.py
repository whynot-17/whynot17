from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import GEOparse
import numpy as np
import pandas as pd
import xenaPython as xena
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests


MODULE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = MODULE_ROOT.parents[1]
ROOT = MODULE_ROOT
OUT = ROOT / "outputs"
GEO_DIR = REPO_ROOT / "work" / "geo_validation"
GENE_LIST = OUT / "DINP_CRC_overlap.csv"

XENA_HUB = "https://toil.xenahubs.net"
XENA_EXPRESSION_DATASET = "TcgaTargetGtex_rsem_gene_tpm"
XENA_PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"

GEO_SOURCES = {
    "GSE74602": {
        "platform": "GPL6104",
        "design": "paired",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74602",
        "description": "30 paired colorectal tumor and normal samples; Illumina humanRef-8 v2.0",
    },
    "GSE10950": {
        "platform": "GPL6104",
        "design": "paired",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10950",
        "description": "24 paired colon tumor and normal mucosa samples; Illumina humanRef-8 v2.0",
    },
    "GSE156355": {
        "platform": "GPL21185",
        "design": "paired",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156355",
        "description": "6 paired colorectal tumor and adjacent normal tissues; Agilent human microarray",
    },
}


def bh_fdr(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(values.astype(float))
    if valid.any():
        out.loc[valid] = multipletests(values.loc[valid].astype(float).to_numpy(), method="fdr_bh")[1]
    return out


def as_float_array(values: object) -> np.ndarray:
    if values is None:
        return np.array([], dtype=float)
    return pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)


def mean_or_nan(values: np.ndarray) -> float:
    return float(np.mean(values)) if len(values) else float("nan")


def median_or_nan(values: np.ndarray) -> float:
    return float(np.median(values)) if len(values) else float("nan")


def safe_wilcoxon(delta: np.ndarray) -> float:
    if len(delta) < 3 or np.allclose(delta, 0):
        return 1.0 if len(delta) >= 3 else float("nan")
    try:
        return float(wilcoxon(delta, alternative="two-sided", method="auto").pvalue)
    except ValueError:
        return 1.0


def safe_mannwhitney(tumor: np.ndarray, normal: np.ndarray) -> float:
    if len(tumor) < 2 or len(normal) < 2:
        return float("nan")
    try:
        return float(mannwhitneyu(tumor, normal, alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def direction(delta: float) -> str:
    if not np.isfinite(delta) or abs(delta) < 1e-12:
        return "flat_or_unavailable"
    return "tumor_high" if delta > 0 else "normal_high"


def patient_id_from_tcga(sample_id: str) -> str:
    return "-".join(str(sample_id).split("-")[:3])


def paired_or_unpaired_rows(
    dataset_id: str,
    platform: str,
    expression: pd.DataFrame,
    sample_meta: pd.DataFrame,
    genes: list[str],
    comparison: str,
    test_mode: str,
    expression_scale: str,
    contrast_role: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    tumor_samples = sample_meta.loc[sample_meta["group"].eq("tumor"), "sample_id"].tolist()
    normal_samples = sample_meta.loc[sample_meta["group"].eq("normal"), "sample_id"].tolist()
    pair_table = sample_meta.loc[sample_meta["pair_id"].ne(""), ["pair_id", "group", "sample_id"]]

    for gene in genes:
        tumor = as_float_array(expression.loc[gene, tumor_samples]) if gene in expression.index else np.array([])
        normal = as_float_array(expression.loc[gene, normal_samples]) if gene in expression.index else np.array([])
        paired_values: list[tuple[float, float]] = []
        for pair_id, pair_rows in pair_table.groupby("pair_id", sort=True):
            tumor_ids = pair_rows.loc[pair_rows["group"].eq("tumor"), "sample_id"].tolist()
            normal_ids = pair_rows.loc[pair_rows["group"].eq("normal"), "sample_id"].tolist()
            if not tumor_ids or not normal_ids:
                continue
            t = pd.to_numeric(expression.loc[gene, tumor_ids[0]], errors="coerce")
            n = pd.to_numeric(expression.loc[gene, normal_ids[0]], errors="coerce")
            if pd.notna(t) and pd.notna(n):
                paired_values.append((float(t), float(n)))

        paired = np.asarray(paired_values, dtype=float) if paired_values else np.empty((0, 2), dtype=float)
        if test_mode == "paired":
            t_used = paired[:, 0] if len(paired) else np.array([])
            n_used = paired[:, 1] if len(paired) else np.array([])
            delta_values = t_used - n_used
            p_value = safe_wilcoxon(delta_values)
            paired_n = len(delta_values)
            concordance_fraction = (
                float(np.mean(delta_values > 0)) if len(delta_values) else float("nan")
            )
            concordant_pairs = int(np.sum(delta_values > 0)) if len(delta_values) else 0
            direction_basis = median_or_nan(delta_values)
        else:
            t_used = tumor
            n_used = normal
            delta_values = np.array([])
            p_value = safe_mannwhitney(t_used, n_used)
            paired_n = 0
            concordance_fraction = float("nan")
            concordant_pairs = 0
            direction_basis = median_or_nan(t_used) - median_or_nan(n_used)

        rows.append(
            {
                "dataset_id": dataset_id,
                "platform": platform,
                "comparison": comparison,
                "contrast_role": contrast_role,
                "test": "paired_wilcoxon" if test_mode == "paired" else "mann_whitney_u",
                "expression_scale": expression_scale,
                "gene_symbol": gene,
                "n_tumor": len(t_used),
                "n_normal": len(n_used),
                "paired_n": paired_n,
                "tumor_mean": mean_or_nan(t_used),
                "normal_mean": mean_or_nan(n_used),
                "tumor_median": median_or_nan(t_used),
                "normal_median": median_or_nan(n_used),
                "delta_tumor_minus_normal": direction_basis,
                "p_value": p_value,
                "concordant_pairs_tumor_high": concordant_pairs,
                "pair_concordance_fraction_tumor_high": concordance_fraction,
                "direction": direction(direction_basis),
            }
        )

    result = pd.DataFrame(rows)
    result["fdr_bh_within_contrast"] = result.groupby(lambda _: True)["p_value"].transform(bh_fdr)
    result["significant_fdr_lt_0_05"] = result["fdr_bh_within_contrast"].lt(0.05)
    result["measured"] = result["n_tumor"].gt(0) & result["n_normal"].gt(0)
    return result


def parse_title_and_characteristics(sample: object) -> tuple[str, str]:
    title = " ".join(str(x) for x in sample.metadata.get("title", []))
    chars = " ".join(str(x) for x in sample.metadata.get("characteristics_ch1", []))
    text = f"{title} {chars}".lower()
    if "normal" in text:
        group = "normal"
    elif "tumor" in text or "tumour" in text:
        group = "tumor"
    else:
        group = ""
    return title, group


def geo_pair_id(accession: str, title: str) -> str:
    if accession == "GSE74602":
        match = re.search(r"^(.+)_([A-Za-z])$", title.strip())
        if not match:
            return ""
        batch, letter = match.groups()
        pair_index = (ord(letter.upper()) - ord("A")) // 2 + 1
        return f"{batch}_pair{pair_index}"
    if accession == "GSE10950":
        match = re.search(r"\b[NT](\d+)\b", title, flags=re.IGNORECASE)
        return f"Patient{match.group(1)}" if match else ""
    match = re.search(r"patient\s*([0-9]+)", title, flags=re.IGNORECASE)
    return f"Patient{match.group(1)}" if match else ""


def build_geo_dataset(accession: str, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    path = GEO_DIR / f"{accession}_family.soft.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing GEO family SOFT file: {path}")
    geo = GEOparse.get_GEO(filepath=str(path), silent=True, include_data=True, annotate_gpl=True)
    platform_id = next(iter(geo.gpls))
    platform = geo.gpls[platform_id].table.copy()
    symbol_col = "Symbol" if "Symbol" in platform.columns else "GENE_SYMBOL"
    platform["probe_id"] = platform["ID"].astype(str)
    probe_to_genes: dict[str, set[str]] = {}
    target_set = set(genes)
    for row in platform[["probe_id", symbol_col]].itertuples(index=False):
        raw = str(row[1])
        if not raw or raw.lower() in {"nan", "---", "none"}:
            continue
        tokens = [x.strip().upper() for x in re.split(r"///|//|;|,|\|", raw) if x.strip()]
        keep = set(tokens).intersection(target_set)
        if keep:
            probe_to_genes[row[0]] = keep

    sample_meta_rows: list[dict[str, object]] = []
    sample_values: dict[str, dict[str, float]] = {gene: {} for gene in genes}
    for sample_id, sample in geo.gsms.items():
        title, group = parse_title_and_characteristics(sample)
        if group not in {"tumor", "normal"}:
            continue
        pair_id = geo_pair_id(accession, title)
        table = sample.table[["ID_REF", "VALUE"]].copy()
        table["VALUE"] = pd.to_numeric(table["VALUE"], errors="coerce")
        per_gene: dict[str, list[float]] = {gene: [] for gene in genes}
        for probe_id, value in table.itertuples(index=False):
            if pd.isna(value):
                continue
            for gene in probe_to_genes.get(str(probe_id), set()):
                per_gene[gene].append(float(value))
        for gene, values in per_gene.items():
            if values:
                sample_values[gene][sample_id] = float(np.median(values))
        sample_meta_rows.append(
            {
                "dataset_id": accession,
                "sample_id": sample_id,
                "group": group,
                "pair_id": pair_id,
                "title": title,
                "platform": platform_id,
            }
        )

    sample_meta = pd.DataFrame(sample_meta_rows)
    sample_ids = sample_meta["sample_id"].tolist()
    expression = pd.DataFrame(
        {
            sample_id: [sample_values[gene].get(sample_id, np.nan) for gene in genes]
            for sample_id in sample_ids
        },
        index=genes,
    )
    dataset_manifest = {
        "dataset_id": accession,
        "platform": platform_id,
        "source_url": GEO_SOURCES[accession]["source_url"],
        "description": GEO_SOURCES[accession]["description"],
        "sample_count": int(len(sample_meta)),
        "tumor_count": int(sample_meta["group"].eq("tumor").sum()),
        "normal_count": int(sample_meta["group"].eq("normal").sum()),
        "paired_count": int(
            sample_meta.groupby("pair_id")["group"].nunique().eq(2).sum()
            if len(sample_meta)
            else 0
        ),
        "target_genes_with_any_value": int(expression.notna().any(axis=1).sum()),
        "target_genes_with_complete_group": int(
            ((expression.loc[:, sample_meta.loc[sample_meta.group.eq("tumor"), "sample_id"]].notna().any(axis=1))
             & (expression.loc[:, sample_meta.loc[sample_meta.group.eq("normal"), "sample_id"]].notna().any(axis=1))).sum()
        ),
    }
    return expression, sample_meta, dataset_manifest


def xena_categorical_values(
    hub: str, dataset: str, samples: list[str], field: str
) -> list[str]:
    raw = xena.dataset_probe_values(hub, dataset, samples, [field])
    values = raw[1][0] if raw and len(raw) > 1 and raw[1] else []
    code_rows = xena.field_codes(hub, dataset, [field])
    code_text = code_rows[0].get("code") if code_rows else None
    codes = str(code_text).split("\t") if code_text else []
    decoded: list[str] = []
    for value in values:
        try:
            if value is None or (isinstance(value, float) and math.isnan(value)):
                decoded.append("")
                continue
            index = int(float(value))
            decoded.append(codes[index] if 0 <= index < len(codes) else str(value))
        except (TypeError, ValueError):
            decoded.append(str(value))
    return decoded


def build_tcga_xena(genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    phenotype_samples = xena.dataset_samples(XENA_HUB, XENA_PHENOTYPE_DATASET, None)
    expression_samples = set(xena.dataset_samples(XENA_HUB, XENA_EXPRESSION_DATASET, None))
    samples = [sample for sample in phenotype_samples if sample in expression_samples]
    phenotype = pd.DataFrame({"sample_id": samples})
    for field in ["_study", "primary disease or tissue", "_sample_type"]:
        phenotype[field] = xena_categorical_values(XENA_HUB, XENA_PHENOTYPE_DATASET, samples, field)

    phenotype["group"] = ""
    tcga_coad = phenotype["_study"].eq("TCGA") & phenotype["primary disease or tissue"].eq("Colon Adenocarcinoma")
    phenotype.loc[tcga_coad & phenotype["_sample_type"].eq("Primary Tumor"), "group"] = "tumor"
    phenotype.loc[tcga_coad & phenotype["_sample_type"].eq("Solid Tissue Normal"), "group"] = "normal"
    gtex_colon = phenotype["_study"].eq("GTEX") & phenotype["primary disease or tissue"].isin(
        ["Colon - Transverse", "Colon - Sigmoid"]
    )
    phenotype.loc[gtex_colon, "group"] = "normal"
    phenotype.loc[tcga_coad, "cohort"] = "TCGA-COAD"
    phenotype.loc[gtex_colon, "cohort"] = "GTEx-colon"
    phenotype["cohort"] = phenotype["cohort"].fillna("")
    phenotype["pair_id"] = phenotype["sample_id"].map(patient_id_from_tcga)
    phenotype = phenotype.loc[phenotype["group"].isin({"tumor", "normal"})].copy()
    phenotype["dataset_id"] = phenotype["cohort"]
    phenotype["platform"] = "Toil/Xena TcgaTargetGtex_rsem_gene_tpm"
    phenotype["title"] = ""

    selected_samples = phenotype["sample_id"].tolist()
    expression = pd.DataFrame(index=genes, columns=selected_samples, dtype=float)
    for start in range(0, len(genes), 40):
        chunk = genes[start : start + 40]
        records = xena.dataset_gene_probe_avg(XENA_HUB, XENA_EXPRESSION_DATASET, selected_samples, chunk)
        for record in records:
            gene = str(record.get("gene", "")).upper()
            scores = record.get("scores") or []
            values = scores[0] if scores and isinstance(scores[0], list) else scores
            if gene in expression.index and len(values) == len(selected_samples):
                expression.loc[gene, selected_samples] = pd.to_numeric(values, errors="coerce")

    manifest = {
        "dataset_id": "TCGA-COAD+GTEx-colon",
        "platform": "Toil/Xena TcgaTargetGtex_rsem_gene_tpm",
        "source_url": XENA_HUB,
        "expression_dataset": XENA_EXPRESSION_DATASET,
        "phenotype_dataset": XENA_PHENOTYPE_DATASET,
        "tcga_tumor_count": int((phenotype["cohort"].eq("TCGA-COAD") & phenotype["group"].eq("tumor")).sum()),
        "tcga_solid_normal_count": int((phenotype["cohort"].eq("TCGA-COAD") & phenotype["group"].eq("normal")).sum()),
        "gtex_colon_normal_count": int((phenotype["cohort"].eq("GTEx-colon") & phenotype["group"].eq("normal")).sum()),
        "target_genes_with_any_value": int(expression.notna().any(axis=1).sum()),
    }
    return expression, phenotype, manifest


def deduplicate_tcga_for_pairs(meta: pd.DataFrame) -> pd.DataFrame:
    selected = meta.loc[meta["cohort"].eq("TCGA-COAD") & meta["group"].isin({"tumor", "normal"})].copy()
    selected = selected.sort_values(["pair_id", "group", "sample_id"], kind="mergesort")
    return selected.drop_duplicates(["pair_id", "group"], keep="first").reset_index(drop=True)


def write_long_expression(
    dataset_id: str, expression: pd.DataFrame, meta: pd.DataFrame, scale: str, rows: list[dict[str, object]]
) -> None:
    for sample_id in meta["sample_id"].tolist():
        for gene in expression.index:
            value = expression.at[gene, sample_id] if sample_id in expression.columns else np.nan
            if pd.notna(value):
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "sample_id": sample_id,
                        "group": meta.loc[meta["sample_id"].eq(sample_id), "group"].iloc[0],
                        "pair_id": meta.loc[meta["sample_id"].eq(sample_id), "pair_id"].iloc[0],
                        "gene_symbol": gene,
                        "expression_value": float(value),
                        "expression_scale": scale,
                    }
                )


def main() -> None:
    genes = (
        pd.read_csv(GENE_LIST, dtype=str)["gene_symbol"].astype(str).str.strip().str.upper().drop_duplicates().tolist()
    )
    OUT.mkdir(parents=True, exist_ok=True)

    result_frames: list[pd.DataFrame] = []
    sample_frames: list[pd.DataFrame] = []
    long_rows: list[dict[str, object]] = []
    manifests: list[dict[str, object]] = []

    tcga_expression, tcga_meta, tcga_manifest = build_tcga_xena(genes)
    manifests.append(tcga_manifest)
    sample_frames.append(tcga_meta.copy())

    tcga_internal = tcga_meta.loc[tcga_meta["cohort"].eq("TCGA-COAD")].copy()
    result_frames.append(
        paired_or_unpaired_rows(
            "TCGA-COAD",
            "Toil/Xena",
            tcga_expression,
            tcga_internal,
            genes,
            "TCGA primary tumor vs solid tissue normal (all samples)",
            "unpaired",
            "Toil/Xena RSEM TPM delivered on log2-like scale",
            "primary_TCGA_internal",
        )
    )
    tcga_pairs = deduplicate_tcga_for_pairs(tcga_meta)
    result_frames.append(
        paired_or_unpaired_rows(
            "TCGA-COAD_paired",
            "Toil/Xena",
            tcga_expression,
            tcga_pairs,
            genes,
            "TCGA primary tumor vs solid tissue normal (one sample per matched patient)",
            "paired",
            "Toil/Xena RSEM TPM delivered on log2-like scale",
            "sensitivity_TCGA_matched",
        )
    )
    tcga_gtex = tcga_meta.loc[tcga_meta["cohort"].isin({"TCGA-COAD", "GTEx-colon"})].copy()
    result_frames.append(
        paired_or_unpaired_rows(
            "TCGA-COAD_vs_GTEx-colon",
            "Toil/Xena",
            tcga_expression,
            tcga_gtex,
            genes,
            "TCGA primary tumor vs GTEx transverse/sigmoid colon normal",
            "unpaired",
            "Toil/Xena RSEM TPM delivered on log2-like scale",
            "sensitivity_external_GTEx_normal",
        )
    )
    write_long_expression("TCGA-COAD+GTEx-colon", tcga_expression, tcga_meta, "Toil/Xena RSEM TPM delivered on log2-like scale", long_rows)

    for accession in GEO_SOURCES:
        expression, meta, manifest = build_geo_dataset(accession, genes)
        manifests.append(manifest)
        sample_frames.append(meta.copy())
        result_frames.append(
            paired_or_unpaired_rows(
                accession,
                manifest["platform"],
                expression,
                meta,
                genes,
                "CRC tumor vs matched normal",
                "paired",
                "GEO processed within-series expression scale",
                "primary_GEO_paired",
            )
        )
        write_long_expression(accession, expression, meta, "GEO processed within-series expression scale", long_rows)

    results = pd.concat(result_frames, ignore_index=True)
    results["p_value"] = pd.to_numeric(results["p_value"], errors="coerce")
    results["fdr_bh_within_contrast"] = pd.to_numeric(results["fdr_bh_within_contrast"], errors="coerce")
    results = results.sort_values(["contrast_role", "dataset_id", "fdr_bh_within_contrast", "gene_symbol"], kind="mergesort")
    results.to_csv(OUT / "DINP_CRC_TCGA_GEO_expression_validation.csv", index=False, encoding="utf-8-sig")

    sample_manifest = pd.concat(sample_frames, ignore_index=True, sort=False)
    sample_manifest.to_csv(OUT / "DINP_CRC_TCGA_GEO_sample_manifest.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(long_rows).to_csv(OUT / "DINP_CRC_TCGA_GEO_target_expression_long.csv", index=False, encoding="utf-8-sig")

    primary_roles = {"primary_TCGA_internal", "primary_GEO_paired"}
    primary_results = results.loc[results["contrast_role"].isin(primary_roles)].copy()
    sensitivity_results = results.loc[~results["contrast_role"].isin(primary_roles)].copy()
    summary_rows: list[dict[str, object]] = []
    for gene in genes:
        group = primary_results.loc[primary_results["gene_symbol"].eq(gene)].copy()
        sensitivity_group = sensitivity_results.loc[sensitivity_results["gene_symbol"].eq(gene)].copy()
        measured = group.loc[group["measured"]].copy()
        significant = measured.loc[measured["fdr_bh_within_contrast"].lt(0.05)]
        signs = measured.loc[measured["delta_tumor_minus_normal"].ne(0), "delta_tumor_minus_normal"].astype(float)
        up = int((signs > 0).sum())
        down = int((signs < 0).sum())
        consensus = "tumor_high" if up > down else ("normal_high" if down > up else "mixed_or_flat")
        summary_rows.append(
            {
                "gene_symbol": gene,
                "independent_datasets_measured": len(measured),
                "independent_datasets_significant_fdr_lt_0_05": len(significant),
                "independent_datasets_tumor_high": up,
                "independent_datasets_normal_high": down,
                "direction_concordance_fraction": max(up, down) / len(signs) if len(signs) else np.nan,
                "consensus_direction": consensus,
                "significant_datasets": "; ".join(significant["dataset_id"].tolist()),
                "significant_directions": "; ".join(
                    f"{d}:{g}" for d, g in zip(significant["dataset_id"], significant["direction"])
                ),
                "sensitivity_contrasts_significant_fdr_lt_0_05": int(
                    sensitivity_group["fdr_bh_within_contrast"].lt(0.05).sum()
                ),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values(
        ["independent_datasets_significant_fdr_lt_0_05", "direction_concordance_fraction", "gene_symbol"],
        ascending=[False, False, True],
        kind="mergesort",
    )
    summary.to_csv(OUT / "DINP_CRC_TCGA_GEO_cross_dataset_summary.csv", index=False, encoding="utf-8-sig")

    manifest_payload = {
        "run_date": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "input_gene_count": len(genes),
        "input_gene_file": str(GENE_LIST.name),
        "xena": tcga_manifest,
        "geo": manifests[1:],
        "contrasts": results["dataset_id"].drop_duplicates().tolist(),
        "result_rows": len(results),
        "summary_rows": len(summary),
    }
    (OUT / "DINP_CRC_TCGA_GEO_expression_validation_manifest.json").write_text(
        json.dumps(manifest_payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    report_lines = [
        "# DINP–CRC TCGA/GEO transcriptomic validation",
        "",
        "## Scope",
        "",
        f"The frozen input was the {len(genes)}-gene `DINP_CRC_overlap.csv` list. The analysis tests tumor-versus-normal expression only; it does not infer DINP causality and does not perform MR, survival, PPI, enrichment, TCGA mutation, or single-cell analysis.",
        "",
        "## Statistical design",
        "",
        "- GEO paired cohorts and the TCGA matched-patient sensitivity use two-sided paired Wilcoxon signed-rank tests.",
        "- All-sample TCGA and TCGA-versus-GTEx contrasts use two-sided Mann–Whitney U tests.",
        "- BH FDR is computed separately within each contrast across the 97 target genes.",
        "- Effect direction is tumor median minus normal median on the delivered within-dataset expression scale; positive means tumor-high.",
        "- Cross-dataset summary reports direction/sign consistency and does not pool platform-specific expression scales into a causal meta-analysis.",
        "",
        "## Cohorts",
        "",
    ]
    for manifest in manifests:
        report_lines.append(
            f"- `{manifest['dataset_id']}`: tumor={manifest.get('tumor_count', manifest.get('tcga_tumor_count', 'NA'))}, normal={manifest.get('normal_count', manifest.get('tcga_solid_normal_count', 'NA'))}, GTEx-colon-normal={manifest.get('gtex_colon_normal_count', 'NA')}, paired={manifest.get('paired_count', 'NA')}; {manifest.get('description', manifest.get('platform', ''))}."
        )
    report_lines.extend(["", "## Contrast-level results", ""])
    for dataset_id, group in results.groupby("dataset_id", sort=False):
        measured = group.loc[group["measured"]]
        sig = measured.loc[measured["fdr_bh_within_contrast"].lt(0.05)]
        report_lines.append(
            f"- `{dataset_id}`: {len(measured)}/{len(genes)} genes measured; {len(sig)} FDR<0.05; tumor-high={int((sig['direction'] == 'tumor_high').sum())}; normal-high={int((sig['direction'] == 'normal_high').sum())}."
        )
    repeated = summary.loc[summary["independent_datasets_significant_fdr_lt_0_05"].ge(2)].head(30)
    report_lines.extend(["", "## Genes with at least two significant independent primary validations", ""])
    if repeated.empty:
        report_lines.append("No gene met FDR<0.05 in at least two analyzed contrasts.")
    else:
        for row in repeated.itertuples(index=False):
            report_lines.append(
                f"- `{row.gene_symbol}`: {row.independent_datasets_significant_fdr_lt_0_05} independent primary datasets; consensus={row.consensus_direction}; concordance={row.direction_concordance_fraction:.2f}."
            )
    report_lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "A tumor-versus-normal difference validates CRC disease-state expression association for the target, not DINP exposure response. Directional discordance across cohorts should be retained as heterogeneity rather than collapsed into a single claim. TCGA internal normal and GTEx colon normal are reported separately because normal-reference choice can change the contrast.",
            "",
            "## Source links",
            "",
            f"- UCSC Toil/Xena hub: {XENA_HUB}; expression dataset `{XENA_EXPRESSION_DATASET}`; phenotype dataset `{XENA_PHENOTYPE_DATASET}`.",
        ]
    )
    for accession, info in GEO_SOURCES.items():
        report_lines.append(f"- {accession}: {info['source_url']}")
    report_lines.extend(
        [
            "",
            "## Files",
            "",
            "- `DINP_CRC_TCGA_GEO_expression_validation.csv`: gene-level contrast results.",
            "- `DINP_CRC_TCGA_GEO_cross_dataset_summary.csv`: cross-contrast directional/significance summary.",
            "- `DINP_CRC_TCGA_GEO_sample_manifest.csv`: sample groups and pairing metadata.",
            "- `DINP_CRC_TCGA_GEO_target_expression_long.csv`: target-gene expression values retained for audit.",
            "- `DINP_CRC_TCGA_GEO_expression_validation_manifest.json`: data and run manifest.",
        ]
    )
    (OUT / "DINP_CRC_TCGA_GEO_expression_validation_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    log_lines = [
        "# DINP–CRC TCGA/GEO expression validation audit log",
        "",
        f"- Run date: {manifest_payload['run_date']}",
        f"- Input: `{GENE_LIST.name}`; genes: {len(genes)}",
        "- No raw TCGA/GEO matrix is committed; only target-gene values, sample metadata and derived statistics are retained in outputs.",
        "- GEO values were aggregated to HGNC symbols by the median across probes mapping to the same target symbol within each sample.",
        "- TCGA values were queried by HGNC symbol from the UCSC Toil/Xena RSEM TPM dataset.",
        "- FDR correction is per contrast across the target list; it is not pooled across datasets.",
        "",
        "## Dataset manifests",
        "",
        json.dumps(manifests, ensure_ascii=False, indent=2, default=str),
        "",
        "## Output integrity",
        "",
    ]
    for filename in [
        "DINP_CRC_TCGA_GEO_expression_validation.csv",
        "DINP_CRC_TCGA_GEO_cross_dataset_summary.csv",
        "DINP_CRC_TCGA_GEO_sample_manifest.csv",
        "DINP_CRC_TCGA_GEO_target_expression_long.csv",
        "DINP_CRC_TCGA_GEO_expression_validation_manifest.json",
    ]:
        digest = hashlib.sha256((OUT / filename).read_bytes()).hexdigest()
        log_lines.append(f"- `{filename}` SHA-256: `{digest}`")
    (OUT / "DINP_CRC_TCGA_GEO_expression_validation_log.md").write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
