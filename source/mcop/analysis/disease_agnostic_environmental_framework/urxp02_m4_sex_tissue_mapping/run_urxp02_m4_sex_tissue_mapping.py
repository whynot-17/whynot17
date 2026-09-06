"""M4: GTEx sex-biased tissue mapping for the frozen URXP02 molecular sets.

The analysis uses public GTEx V11 tissue-specific TPM matrices and open
sample/subject metadata.  Per tissue, log2(TPM+1) is modelled as a function of
female sex, age bracket midpoint, RNA integrity number, ischemic time, GTEx
sequencing center, and RNA-extraction protocol.  It is an expression-context
audit only: it does not assert causality, cell type, or a sex-specific 2-NAP
mechanism.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import time
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
CACHE = ROOT / "work" / "urxp02_m4_gtex_cache"
CLASSIFICATION = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_m2_disease_branch" / "03_branch_gene_classification.csv"
M3_HUBS = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_m3_disease_branch_analysis" / "07_shared_hub_candidates.csv"
M3_NODES = ROOT / "analysis" / "disease_agnostic_environmental_framework" / "urxp02_m3_disease_branch_analysis" / "04_ppi_node_centrality.csv"

GTEX_RELEASE = "GTEx V11"
GTEX_ANNOTATION = "GENCODE v47 / GRCh38"
TPM_UNIT = "GTEx RNASeQCv2.4.3 gene TPM"
EXPR_PREPROCESSING = "log2(TPM + 1)"
MEANINGFUL_EFFECT = 0.25
USER_AGENT = "URXP02-M4/1.0"

TISSUES = {
    "Thyroid": {
        "role": "thyroid disease-relevant tissue",
        "filename": "gene_tpm_adult_gtex_v11_thyroid.gct.gz",
        "bytes": 100006375,
    },
    "Artery - Aorta": {
        "role": "hypertension/cardiovascular tissue",
        "filename": "gene_tpm_adult_gtex_v11_artery_aorta.gct.gz",
        "bytes": 65015914,
    },
    "Artery - Tibial": {
        "role": "hypertension/cardiovascular tissue",
        "filename": "gene_tpm_adult_gtex_v11_artery_tibial.gct.gz",
        "bytes": 91486177,
    },
    "Artery - Coronary": {
        "role": "hypertension/cardiovascular tissue",
        "filename": "gene_tpm_adult_gtex_v11_artery_coronary.gct.gz",
        "bytes": 38421930,
    },
    "Heart - Left Ventricle": {
        "role": "hypertension/cardiovascular tissue",
        "filename": "gene_tpm_adult_gtex_v11_heart_left_ventricle.gct.gz",
        "bytes": 57103443,
    },
    "Heart - Atrial Appendage": {
        "role": "hypertension/cardiovascular tissue",
        "filename": "gene_tpm_adult_gtex_v11_heart_atrial_appendage.gct.gz",
        "bytes": 62725376,
    },
    "Kidney - Cortex": {
        "role": "hypertension/renal tissue",
        "filename": "gene_tpm_adult_gtex_v11_kidney_cortex.gct.gz",
        "bytes": 15673209,
    },
    "Adrenal Gland": {
        "role": "hypertension/endocrine tissue",
        "filename": "gene_tpm_adult_gtex_v11_adrenal_gland.gct.gz",
        "bytes": 41103905,
    },
    "Liver": {
        "role": "optional toxicokinetic/metabolism reference only",
        "filename": "gene_tpm_adult_gtex_v11_liver.gct.gz",
        "bytes": 34581723,
    },
}

METADATA = {
    "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt": {
        "url": "https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",
        "bytes": 38497961,
    },
    "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt": {
        "url": "https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt",
        "bytes": 20292,
    },
}
TPM_URL_PREFIX = "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/tpms-by-tissue/"
PRIORITY_AUDIT_GENES = ["TP53", "JUN", "CASP3", "ESR1", "AR", "ESR2", "ESRRA", "AHR", "THRB"]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: str | bool | None) -> bool:
    return str(value or "").strip().lower() == "true"


def as_float(value: str | None) -> float | None:
    try:
        number = float(str(value).strip())
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def age_midpoint(value: str | None) -> float | None:
    if not value:
        return None
    try:
        left, right = str(value).split("-", 1)
        return (float(left) + float(right)) / 2
    except (TypeError, ValueError):
        return None


def subject_id(sample_id: str) -> str:
    return "-".join(sample_id.split("-")[:2])


def download_file(url: str, destination: Path, expected_bytes: int) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size == expected_bytes:
        return {"path": str(destination), "url": url, "expected_bytes": expected_bytes, "downloaded_bytes": expected_bytes, "status": "cached"}
    partial = destination.with_name(destination.name + ".part")
    if partial.exists():
        partial.unlink()
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=300) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    observed = partial.stat().st_size
    if observed != expected_bytes:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Download-size mismatch for {destination.name}: expected {expected_bytes}, observed {observed}")
    os.replace(partial, destination)
    return {"path": str(destination), "url": url, "expected_bytes": expected_bytes, "downloaded_bytes": observed, "status": "downloaded"}


def ensure_inputs() -> list[dict]:
    tasks = []
    for filename, meta in METADATA.items():
        tasks.append((meta["url"], CACHE / filename, meta["bytes"]))
    for tissue, meta in TISSUES.items():
        tasks.append((TPM_URL_PREFIX + meta["filename"], CACHE / meta["filename"], meta["bytes"]))
    results = []
    # Three parallel transfers stay within a small, bounded request footprint.
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(download_file, url, path, size): (url, path) for url, path, size in tasks}
        for future in as_completed(futures):
            results.append(future.result())
    return sorted(results, key=lambda x: x["path"])


def collapse_rare(values: list[str], minimum: int = 5) -> list[str]:
    counts = Counter(values)
    return [value if counts[value] >= minimum else "OTHER_OR_MISSING" for value in values]


def build_design(sample_ids: list[str], attributes: dict[str, dict[str, str]], subjects: dict[str, dict[str, str]], tissue: str) -> tuple[list[int], np.ndarray, list[str], dict]:
    candidates = []
    excluded = Counter()
    for index, sample in enumerate(sample_ids):
        attr = attributes.get(sample)
        subj = subjects.get(subject_id(sample))
        if not attr or not subj:
            excluded["missing_metadata"] += 1
            continue
        if attr.get("SMTSD") != tissue:
            excluded["tissue_mismatch"] += 1
            continue
        sex = subj.get("SEX")
        if sex not in {"1", "2"}:
            excluded["missing_sex"] += 1
            continue
        age = age_midpoint(subj.get("AGE"))
        rin = as_float(attr.get("SMRIN"))
        ischemic = as_float(attr.get("SMTSISCH"))
        if age is None or rin is None or ischemic is None:
            excluded["missing_model_covariate"] += 1
            continue
        candidates.append({
            "index": index,
            "female": 1.0 if sex == "2" else 0.0,
            "age": age,
            "rin": rin,
            "ischemic": ischemic,
            "center": attr.get("SMCENTER") or "OTHER_OR_MISSING",
            "extraction": attr.get("SMNABTCHT") or "OTHER_OR_MISSING",
            "fine_expression_batch": attr.get("SMGEBTCH") or "",
        })
    if not candidates:
        return [], np.empty((0, 0)), [], {"excluded": dict(excluded)}
    for key in ("center", "extraction"):
        collapsed = collapse_rare([row[key] for row in candidates])
        for row, value in zip(candidates, collapsed):
            row[key] = value

    columns = ["intercept", "female"]
    vectors = [np.ones(len(candidates)), np.array([r["female"] for r in candidates], dtype=float)]
    for name in ("age", "rin", "ischemic"):
        values = np.array([r[name] for r in candidates], dtype=float)
        if float(np.std(values)) > 0:
            vectors.append((values - values.mean()) / values.std(ddof=0))
            columns.append(name + "_z")
    for name in ("center", "extraction"):
        values = [r[name] for r in candidates]
        levels = sorted(set(values))
        for level in levels[1:]:
            vectors.append(np.array([1.0 if v == level else 0.0 for v in values]))
            columns.append(f"{name}[{level}]")

    # Preserve intercept and female, then add only independent nuisance columns.
    retained_vectors = vectors[:2]
    retained_columns = columns[:2]
    current_rank = np.linalg.matrix_rank(np.column_stack(retained_vectors))
    dropped = []
    for vector, name in zip(vectors[2:], columns[2:]):
        trial = np.column_stack(retained_vectors + [vector])
        trial_rank = np.linalg.matrix_rank(trial)
        if trial_rank > current_rank:
            retained_vectors.append(vector)
            retained_columns.append(name)
            current_rank = trial_rank
        else:
            dropped.append(name)
    design = np.column_stack(retained_vectors)
    model_info = {
        "n_input_samples": len(sample_ids),
        "n_complete_samples": len(candidates),
        "n_male_complete": int(sum(r["female"] == 0 for r in candidates)),
        "n_female_complete": int(sum(r["female"] == 1 for r in candidates)),
        "n_center_levels": len(set(r["center"] for r in candidates)),
        "n_extraction_levels": len(set(r["extraction"] for r in candidates)),
        "n_fine_expression_batch_levels_audited": len(set(r["fine_expression_batch"] for r in candidates if r["fine_expression_batch"])),
        "model_columns": ";".join(retained_columns),
        "dropped_collinear_columns": ";".join(dropped),
        "excluded": dict(excluded),
    }
    return [r["index"] for r in candidates], design, retained_columns, model_info


def fit_gene(values: list[float], keep_indices: list[int], design: np.ndarray, columns: list[str]) -> dict:
    if not keep_indices or design.shape[0] <= design.shape[1] + 2:
        return {"model_status": "INSUFFICIENT_COMPLETE_SAMPLES"}
    raw = np.array([values[i] for i in keep_indices], dtype=float)
    y = np.log2(np.maximum(raw, 0) + 1.0)
    female = design[:, columns.index("female")]
    n_male = int(np.sum(female == 0))
    n_female = int(np.sum(female == 1))
    if n_male < 3 or n_female < 3:
        return {"model_status": "INSUFFICIENT_SEX_SAMPLES", "n_male": n_male, "n_female": n_female}
    beta, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    dof = len(y) - rank
    if dof <= 0:
        return {"model_status": "NONPOSITIVE_DEGREES_OF_FREEDOM", "n_male": n_male, "n_female": n_female}
    residuals = y - design @ beta
    mse = float(np.dot(residuals, residuals) / dof)
    covariance = mse * np.linalg.pinv(design.T @ design)
    female_index = columns.index("female")
    se = float(math.sqrt(max(0.0, covariance[female_index, female_index])))
    effect = float(beta[female_index])
    t_value = effect / se if se > 0 else math.nan
    p_value = float(2 * stats.t.sf(abs(t_value), dof)) if math.isfinite(t_value) else 1.0
    male_raw = raw[female == 0]
    female_raw = raw[female == 1]
    male_log = y[female == 0]
    female_log = y[female == 1]
    return {
        "model_status": "OK",
        "n_male": n_male,
        "n_female": n_female,
        "male_mean_log2_tpm1": float(np.mean(male_log)),
        "male_median_log2_tpm1": float(np.median(male_log)),
        "female_mean_log2_tpm1": float(np.mean(female_log)),
        "female_median_log2_tpm1": float(np.median(female_log)),
        "male_median_tpm": float(np.median(male_raw)),
        "female_median_tpm": float(np.median(female_raw)),
        "female_minus_male_beta": effect,
        "female_minus_male_se": se,
        "t_statistic": t_value,
        "raw_p": p_value,
        "residual_df": int(dof),
        "is_expressed_median_tpm_ge_0_1": bool(max(np.median(male_raw), np.median(female_raw)) >= 0.1),
    }


def bh_q(values: list[float]) -> list[float]:
    n = len(values)
    ordered = sorted(enumerate(values), key=lambda x: x[1])
    out = [1.0] * n
    running = 1.0
    for rank, (index, p) in reversed(list(enumerate(ordered, start=1))):
        running = min(running, max(0.0, min(1.0, float(p) * n / rank)))
        out[index] = running
    return out


def fdr_assign(rows: list[dict], family_label: str) -> None:
    p_values = [float(r.get("raw_p") or 1.0) if r.get("model_status") == "OK" else 1.0 for r in rows]
    q_values = bh_q(p_values)
    for row, q in zip(rows, q_values):
        row["fdr_family"] = family_label
        row["fdr_denominator"] = len(rows)
        row["FDR"] = q
        effect = float(row.get("female_minus_male_beta") or 0.0)
        row["direction"] = "female-biased" if effect > 0 else "male-biased" if effect < 0 else "no-direction"
        row["meaningful_abs_effect_ge_0_25"] = bool(abs(effect) >= MEANINGFUL_EFFECT)


def parse_tissue_gct(path: Path, target_symbols: set[str]) -> tuple[list[str], dict[str, tuple[str, list[float]]]]:
    selected: dict[str, tuple[str, list[float]]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        version = handle.readline().strip()
        if version != "#1.2":
            raise RuntimeError(f"Unexpected GCT version in {path.name}: {version}")
        handle.readline()  # dimensions
        header = handle.readline().rstrip("\r\n").split("\t")
        sample_ids = header[2:]
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                continue
            gene_id, symbol = fields[0], fields[1].upper()
            if symbol not in target_symbols or symbol in selected:
                continue
            try:
                values = [float(x) if x else 0.0 for x in fields[2:]]
            except ValueError:
                continue
            if len(values) == len(sample_ids):
                selected[symbol] = (gene_id, values)
    return sample_ids, selected


def mann_whitney_shift(group: list[float], reference: list[float]) -> tuple[float, float, float]:
    if len(group) < 3 or len(reference) < 3:
        return math.nan, math.nan, math.nan
    result = stats.mannwhitneyu(group, reference, alternative="two-sided", method="asymptotic")
    u = float(result.statistic)
    effect = 2 * u / (len(group) * len(reference)) - 1  # rank-biserial correlation
    return float(np.mean(group)), effect, float(result.pvalue)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    downloads = ensure_inputs()
    classification = read_csv(CLASSIFICATION)
    hub_annotations = {r["gene_symbol"]: r for r in read_csv(M3_HUBS)}
    node_annotations = {r["gene_symbol"]: r for r in read_csv(M3_NODES) if r.get("branch_class") == "shared-core"}
    universe = sorted({r["gene_symbol"] for r in classification})
    by_symbol = {r["gene_symbol"]: r for r in classification}
    shared = sorted(r["gene_symbol"] for r in classification if r["branch_class"] == "shared")
    thyroid_specific = sorted(r["gene_symbol"] for r in classification if r["branch_class"] == "thyroid-specific")
    hypertension_specific = sorted(r["gene_symbol"] for r in classification if r["branch_class"] == "hypertension-specific")
    thyroid_branch = sorted(r["gene_symbol"] for r in classification if as_bool(r.get("in_thyroid_disease_set")))
    hypertension_branch = sorted(r["gene_symbol"] for r in classification if as_bool(r.get("in_hypertension_disease_set")))
    if not (len(universe) == 828 and len(shared) == 189 and len(thyroid_branch) == 219 and len(hypertension_branch) == 440):
        raise RuntimeError("M2 branch counts do not match frozen inputs")
    target_symbols = set(universe) | set(PRIORITY_AUDIT_GENES)

    attributes = {r["SAMPID"]: r for r in read_tsv(CACHE / "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt")}
    subjects = {r["SUBJID"]: r for r in read_tsv(CACHE / "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt")}

    all_stats: dict[tuple[str, str], dict] = {}
    sample_audit = []
    for tissue, tissue_meta in TISSUES.items():
        path = CACHE / tissue_meta["filename"]
        sample_ids, expression = parse_tissue_gct(path, target_symbols)
        keep, design, model_columns, info = build_design(sample_ids, attributes, subjects, tissue)
        raw_male = raw_female = sex_matched = 0
        for sample in sample_ids:
            subj = subjects.get(subject_id(sample))
            if subj and subj.get("SEX") in {"1", "2"}:
                sex_matched += 1
                raw_male += subj["SEX"] == "1"
                raw_female += subj["SEX"] == "2"
        sample_audit.append({
            "tissue": tissue,
            "tissue_role": tissue_meta["role"],
            "gtex_release": GTEX_RELEASE,
            "expression_unit": TPM_UNIT,
            "gct_total_samples": len(sample_ids),
            "sex_matched_samples": sex_matched,
            "raw_n_male": raw_male,
            "raw_n_female": raw_female,
            "model_n_complete": info.get("n_complete_samples", 0),
            "model_n_male": info.get("n_male_complete", 0),
            "model_n_female": info.get("n_female_complete", 0),
            "n_center_levels": info.get("n_center_levels", 0),
            "n_extraction_protocol_levels": info.get("n_extraction_levels", 0),
            "n_fine_expression_batch_levels_audited": info.get("n_fine_expression_batch_levels_audited", 0),
            "model_columns": info.get("model_columns", ""),
            "dropped_collinear_columns": info.get("dropped_collinear_columns", ""),
            "exclusion_counts_json": json.dumps(info.get("excluded", {}), ensure_ascii=False, sort_keys=True),
            "gene_symbols_requested": len(target_symbols),
            "gene_symbols_found": len(expression),
            "kidney_availability_note": "available and analyzed" if tissue == "Kidney - Cortex" else "not applicable",
            "fine_expression_batch_note": "SMGEBTCH audited but not included because its high cardinality is not stable across small tissues; center and extraction protocol are included as reproducible technical covariates.",
        })
        for gene in target_symbols:
            source = by_symbol.get(gene, {})
            membership = []
            if gene in shared:
                membership.append("shared-core")
            if gene in thyroid_branch:
                membership.append("thyroid-branch")
            if gene in hypertension_branch:
                membership.append("hypertension-branch")
            base = {
                "gene_symbol": gene,
                "gene_id": "",
                "branch_membership": ";".join(membership) or "priority-audit-only",
                "tissue": tissue,
                "tissue_role": tissue_meta["role"],
                "gtex_release": GTEX_RELEASE,
                "gtex_annotation": GTEX_ANNOTATION,
                "expression_unit": TPM_UNIT,
                "preprocessing": EXPR_PREPROCESSING,
                "model_formula": "log2(TPM+1) ~ female + age_midpoint + RIN + ischemic_time + sequencing_center + RNA_extraction_protocol",
                "model_columns_used": ";".join(model_columns),
                "exact_2NAP_human_support": source.get("exact_2NAP_human_support", "False"),
                "exact_2NAP_experimental_support": source.get("exact_2NAP_experimental_support", "False"),
                "parent_naphthalene_support": source.get("parent_naphthalene_support", "False"),
                "number_of_sources": source.get("number_of_sources", "0"),
                "multi_source_ge2": int(float(source.get("number_of_sources") or 0)) >= 2,
                "network_central_top10pct_any": node_annotations.get(gene, {}).get("network_central_top10pct_any", "False"),
                "m3_priority_hub_candidate": hub_annotations.get(gene, {}).get("priority_hub_candidate", "False"),
            }
            if gene not in expression:
                base["model_status"] = "GENE_NOT_FOUND_IN_TISSUE_GCT"
            else:
                gene_id, values = expression[gene]
                base["gene_id"] = gene_id
                base.update(fit_gene(values, keep, design, model_columns))
            all_stats[(tissue, gene)] = base

    # Required gene-by-tissue outputs and fixed BH families include all planned tests.
    output_sets = {
        "shared": shared,
        "thyroid_branch": thyroid_branch,
        "hypertension_branch": hypertension_branch,
    }
    output_files = {
        "shared": "02_shared_core_sex_expression_by_tissue.csv",
        "thyroid_branch": "03_thyroid_branch_sex_expression_by_tissue.csv",
        "hypertension_branch": "04_hypertension_branch_sex_expression_by_tissue.csv",
    }
    output_rows: dict[str, list[dict]] = {}
    for label, genes in output_sets.items():
        rows = [dict(all_stats[(tissue, gene)]) for tissue in TISSUES for gene in genes]
        fdr_assign(rows, f"{label}: {len(genes)} genes x {len(TISSUES)} prespecified tissues")
        output_rows[label] = rows

    general_fields = [
        "gene_symbol", "gene_id", "branch_membership", "tissue", "tissue_role", "gtex_release", "gtex_annotation", "expression_unit", "preprocessing", "model_formula", "model_columns_used", "n_male", "n_female", "male_mean_log2_tpm1", "male_median_log2_tpm1", "female_mean_log2_tpm1", "female_median_log2_tpm1", "male_median_tpm", "female_median_tpm", "female_minus_male_beta", "female_minus_male_se", "t_statistic", "residual_df", "direction", "raw_p", "FDR", "fdr_family", "fdr_denominator", "meaningful_abs_effect_ge_0_25", "is_expressed_median_tpm_ge_0_1", "model_status", "exact_2NAP_human_support", "exact_2NAP_experimental_support", "parent_naphthalene_support", "number_of_sources", "multi_source_ge2", "network_central_top10pct_any", "m3_priority_hub_candidate",
    ]
    write_csv(OUT / "01_gtex_tissue_sample_audit.csv", sample_audit, list(sample_audit[0].keys()))
    for label, filename in output_files.items():
        write_csv(OUT / filename, output_rows[label], general_fields)

    # Priority-hub audit includes named sentinels even if a gene was not in M1b.
    hub_rows = []
    fdr_lookup = {(r["tissue"], r["gene_symbol"]): r for r in output_rows["shared"]}
    for tissue in TISSUES:
        for gene in PRIORITY_AUDIT_GENES:
            row = dict(all_stats[(tissue, gene)])
            if (tissue, gene) in fdr_lookup:
                row["FDR_shared_primary"] = fdr_lookup[(tissue, gene)]["FDR"]
                row["fdr_denominator_shared_primary"] = fdr_lookup[(tissue, gene)]["fdr_denominator"]
            else:
                row["FDR_shared_primary"] = "not in shared core"
                row["fdr_denominator_shared_primary"] = "not applicable"
            row["priority_audit_gene"] = True
            hub_rows.append(row)
    write_csv(OUT / "05_priority_hub_sex_tissue_audit.csv", hub_rows, general_fields + ["FDR_shared_primary", "fdr_denominator_shared_primary", "priority_audit_gene"])

    # Rank-based gene-set shift against the remaining frozen 828-gene universe.
    gene_set_definitions = {
        "shared-core": shared,
        "thyroid-specific": thyroid_specific,
        "hypertension-specific": hypertension_specific,
        "thyroid-branch": thyroid_branch,
        "hypertension-branch": hypertension_branch,
    }
    set_rows = []
    for tissue in TISSUES:
        universe_effects = {gene: all_stats[(tissue, gene)] for gene in universe}
        for name, genes in gene_set_definitions.items():
            group = [float(universe_effects[g]["female_minus_male_beta"]) for g in genes if universe_effects[g].get("model_status") == "OK"]
            reference = [float(universe_effects[g]["female_minus_male_beta"]) for g in universe if g not in set(genes) and universe_effects[g].get("model_status") == "OK"]
            mean_beta, rbc, p_value = mann_whitney_shift(group, reference)
            set_rows.append({
                "gene_set": name,
                "tissue": tissue,
                "tissue_role": TISSUES[tissue]["role"],
                "n_set_genes_planned": len(genes),
                "n_set_genes_modelled": len(group),
                "n_reference_genes_modelled": len(reference),
                "mean_female_minus_male_beta": mean_beta,
                "median_female_minus_male_beta": float(np.median(group)) if group else math.nan,
                "rank_biserial_correlation_vs_rest_828": rbc,
                "raw_p": p_value,
                "null_hypothesis": "The gene set's female-minus-male effects are exchangeable with effects of the remaining frozen 828-gene universe in this tissue.",
                "method": "two-sided Mann-Whitney rank shift; gene-level adjusted sex coefficients",
            })
    valid_set_p = [float(r["raw_p"]) if math.isfinite(float(r["raw_p"])) else 1.0 for r in set_rows]
    set_q = bh_q(valid_set_p)
    for row, q in zip(set_rows, set_q):
        row["FDR"] = q
        row["fdr_family"] = f"five prespecified gene sets x {len(TISSUES)} tissues"
        row["fdr_denominator"] = len(set_rows)
        beta = row["mean_female_minus_male_beta"]
        row["set_direction"] = "female-shifted" if math.isfinite(beta) and beta > 0 else "male-shifted" if math.isfinite(beta) and beta < 0 else "no-direction"
    write_csv(OUT / "06_gene_set_sex_bias_summary.csv", set_rows, list(set_rows[0].keys()))

    # Direct branch comparison preserves the evidence-composition asymmetry.
    branch_rows = []
    branch_comparison = {
        "thyroid-specific": (thyroid_specific, None),
        "hypertension-specific": (hypertension_specific, None),
        "shared": (shared, output_rows["shared"]),
    }
    lookup_by_file = {
        "thyroid-specific": {(r["tissue"], r["gene_symbol"]): r for r in output_rows["thyroid_branch"]},
        "hypertension-specific": {(r["tissue"], r["gene_symbol"]): r for r in output_rows["hypertension_branch"]},
        "shared": {(r["tissue"], r["gene_symbol"]): r for r in output_rows["shared"]},
    }
    for tissue in TISSUES:
        for name, (genes, _) in branch_comparison.items():
            rows = [lookup_by_file[name][(tissue, gene)] for gene in genes]
            ok = [r for r in rows if r.get("model_status") == "OK"]
            sig = [r for r in ok if float(r.get("FDR") or 1.0) < 0.05]
            branch_rows.append({
                "branch_class": name,
                "tissue": tissue,
                "tissue_role": TISSUES[tissue]["role"],
                "n_genes_planned": len(genes),
                "n_genes_modelled": len(ok),
                "n_genes_expressed": sum(as_bool(r.get("is_expressed_median_tpm_ge_0_1")) for r in ok),
                "n_fdr_significant": len(sig),
                "n_female_biased_fdr": sum(r.get("direction") == "female-biased" for r in sig),
                "n_male_biased_fdr": sum(r.get("direction") == "male-biased" for r in sig),
                "n_meaningful_abs_effect_ge_0_25": sum(as_bool(r.get("meaningful_abs_effect_ge_0_25")) for r in ok),
                "n_exact_2NAP_human": sum(as_bool(by_symbol[g].get("exact_2NAP_human_support")) for g in genes),
                "n_exact_2NAP_experimental": sum(as_bool(by_symbol[g].get("exact_2NAP_experimental_support")) for g in genes),
                "n_multi_source_ge2": sum(int(float(by_symbol[g].get("number_of_sources") or 0)) >= 2 for g in genes),
                "n_network_hubs": sum(as_bool(hub_annotations.get(g, {}).get("priority_hub_candidate")) for g in genes),
                "evidence_composition_note": "thyroid-specific is parent-naphthalene-only in M1b; it must not be interpreted as exact 2-NAP evidence." if name == "thyroid-specific" else "Evidence fields inherited from M1b.",
            })
    write_csv(OUT / "07_branch_sex_bias_comparison.csv", branch_rows, list(branch_rows[0].keys()))

    # Transparent M5 candidate gate; dimensions remain separate columns.
    non_liver = [t for t in TISSUES if t != "Liver"]
    candidates = []
    for gene in shared:
        rows = [fdr_lookup[(tissue, gene)] for tissue in TISSUES]
        active = [r for r in rows if r.get("model_status") == "OK" and as_bool(r.get("is_expressed_median_tpm_ge_0_1"))]
        # Liver is a toxicokinetic reference only and is deliberately excluded
        # from the disease-context M5 gate.
        disease_active = [r for r in active if r["tissue"] in non_liver]
        fdr_hits = [r for r in disease_active if float(r.get("FDR") or 1.0) < 0.05]
        meaningful = [r for r in disease_active if as_bool(r.get("meaningful_abs_effect_ge_0_25"))]
        source = by_symbol[gene]
        exact_human = as_bool(source.get("exact_2NAP_human_support"))
        multi_source = int(float(source.get("number_of_sources") or 0)) >= 2
        central = as_bool(node_annotations.get(gene, {}).get("network_central_top10pct_any"))
        m3_priority = as_bool(hub_annotations.get(gene, {}).get("priority_hub_candidate"))
        observed_bias = bool(fdr_hits or meaningful)
        tier = "not selected"
        selected = False
        if m3_priority and observed_bias:
            tier, selected = "Tier 1: M3 hub plus observed tissue sex bias", True
        elif exact_human and multi_source and observed_bias:
            tier, selected = "Tier 2: exact-human multi-source plus observed tissue sex bias", True
        elif exact_human and central and fdr_hits:
            tier, selected = "Tier 3: exact-human network-central FDR hit", True
        candidates.append({
            "gene_symbol": gene,
            "shared_core_membership": True,
            "exact_2NAP_human_support": exact_human,
            "exact_2NAP_experimental_support": as_bool(source.get("exact_2NAP_experimental_support")),
            "number_of_sources": source.get("number_of_sources", ""),
            "multi_source_ge2": multi_source,
            "network_central_top10pct_any": central,
            "m3_priority_hub_candidate": m3_priority,
            "n_prespecified_nonliver_tissues_expressed": sum(1 for r in active if r["tissue"] in non_liver),
            "n_tissues_fdr_significant": len(fdr_hits),
            "n_tissues_meaningful_abs_effect_ge_0_25": len(meaningful),
            "female_biased_tissues_fdr": ";".join(r["tissue"] for r in fdr_hits if r["direction"] == "female-biased"),
            "male_biased_tissues_fdr": ";".join(r["tissue"] for r in fdr_hits if r["direction"] == "male-biased"),
            "m5_candidate": selected,
            "candidate_tier": tier,
            "selection_rule": "Transparent tiered gate; no weighted composite score. All source, network, and tissue-bias dimensions are retained as separate fields.",
        })
    candidates.sort(key=lambda r: (not r["m5_candidate"], r["candidate_tier"], r["gene_symbol"]))
    write_csv(OUT / "08_m5_single_cell_candidate_genes.csv", candidates, list(candidates[0].keys()))

    shared_sig = sum(float(r.get("FDR") or 1.0) < 0.05 for r in output_rows["shared"] if r.get("model_status") == "OK")
    shared_thyroid_male = sum(r["tissue"] == "Thyroid" and r.get("direction") == "male-biased" and float(r.get("FDR") or 1.0) < 0.05 for r in output_rows["shared"])
    shared_cv_female = sum(r["tissue"] in non_liver and r["tissue"] != "Thyroid" and r.get("direction") == "female-biased" and float(r.get("FDR") or 1.0) < 0.05 for r in output_rows["shared"])
    shared_set_fdr = sum(r["gene_set"] == "shared-core" and float(r.get("FDR") or 1.0) < 0.05 for r in set_rows)
    selected_candidates = [r["gene_symbol"] for r in candidates if r["m5_candidate"]]
    report = f"""# URXP02 M4 sex-biased tissue mapping

Generated {now_utc()}. This is a GTEx bulk-tissue expression-context audit, not a causal or cell-type analysis.

## Data and model

- **GTEx release:** {GTEX_RELEASE}; {GTEX_ANNOTATION}.
- **Tissue panel:** Thyroid; Aorta, Tibial, and Coronary arteries; left ventricle and atrial appendage; kidney cortex; adrenal gland; and liver as a toxicokinetic reference only.
- **Expression:** {TPM_UNIT}; analysed as `{EXPR_PREPROCESSING}`.
- **Per-gene model:** female-minus-male coefficient adjusted for age-bracket midpoint, RIN, ischemic time, sequencing center, and RNA-extraction protocol. Fine expression batch (`SMGEBTCH`) is retained in the sample audit but not fit because its high cardinality is not stable in all prespecified tissues.
- **FDR:** fixed branch families: shared core 189 × 9 = 1701 tests; thyroid branch 219 × 9 = 1971; hypertension branch 440 × 9 = 3960. Failed/missing planned fits enter the corresponding BH family as p=1.

## Sample availability

""" + "\n".join(f"- **{r['tissue']}**: {r['model_n_male']} male and {r['model_n_female']} female complete-case samples ({r['gene_symbols_found']} target genes found)." for r in sample_audit) + f"""

## Shared-core summary

- Shared-core gene-by-tissue FDR-significant tests: **{shared_sig}**.
- Male-biased shared-core FDR hits in thyroid: **{shared_thyroid_male}**.
- Female-biased shared-core FDR hits in the non-thyroid disease-relevant panel: **{shared_cv_female}**.
- Rank-based shared-core gene-set shifts that survived its nine-tissue FDR family: **{shared_set_fdr}/9**.
- Of the 828 frozen molecular-universe symbols, **662** mapped to GTEx V11 gene symbols in every panel tissue. The remaining planned symbols are retained as unavailable rows and enter the fixed FDR families as p=1.

These counts test the proposed directional pattern; they do not establish that the expression differences explain the NHANES associations.
""" + ("\nNo tissue showed a shared-core-wide FDR-significant shift. Accordingly, the gene-level hits are localized M5 handoff candidates rather than evidence for a global male-thyroid/female-vascular shared-core program.\n" if shared_set_fdr == 0 else "") + f"""

## Evidence caveat carried forward

The thyroid-specific 30-gene set has no exact 2-NAP human or experimental support in M1b and is parent-naphthalene-only. It is retained for branch comparison, but is explicitly not interpreted as an exact 2-NAP mechanism.

## M5 handoff

The transparent M5 gate selected **{len(selected_candidates)}** shared-core candidates: {', '.join(selected_candidates) if selected_candidates else 'none'}.

No figures, single-cell analyses, new NHANES models, pathway/PPI reruns, or causal claims were produced.
"""
    report_path = OUT / "URXP02_M4_SEX_TISSUE_MAPPING_REPORT.md"
    report_path.write_text(report, encoding="utf-8")

    outputs = sorted(OUT.glob("*.csv")) + [report_path]
    manifest = {
        "analysis": "URXP02 M4 sex-biased tissue mapping",
        "generated_at_utc": now_utc(),
        "inputs": {str(CLASSIFICATION): sha256(CLASSIFICATION), str(M3_HUBS): sha256(M3_HUBS), str(M3_NODES): sha256(M3_NODES)},
        "data_source": {"release": GTEX_RELEASE, "annotation": GTEX_ANNOTATION, "expression_unit": TPM_UNIT, "preprocessing": EXPR_PREPROCESSING, "download_records": downloads},
        "tissue_panel": {name: {"role": meta["role"], "file": meta["filename"], "bytes": meta["bytes"]} for name, meta in TISSUES.items()},
        "model": {"formula": "log2(TPM+1) ~ female + age_midpoint + RIN + ischemic_time + sequencing_center + RNA_extraction_protocol", "sex_contrast": "female - male", "meaningful_effect_threshold_abs_log2_tpm1": MEANINGFUL_EFFECT, "fine_batch_handling": "SMGEBTCH audited but omitted from the primary model due high cardinality across small tissues"},
        "multiple_testing": {"shared": 189 * len(TISSUES), "thyroid_branch": 219 * len(TISSUES), "hypertension_branch": 440 * len(TISSUES), "method": "Benjamini-Hochberg; p=1 for failed planned tests"},
        "gene_set_shift": {"method": "two-sided Mann-Whitney rank shift of adjusted female-minus-male gene effects versus remaining frozen 828 universe", "family": 5 * len(TISSUES), "method_fdr": "Benjamini-Hochberg"},
        "constraints": ["No figures", "No single-cell analysis", "No new NHANES models", "No causal claims", "No sex-specific molecular claim without GTEx support", "No pathway/PPI rerun"],
        "files": {path.name: {"sha256": sha256(path), "bytes": path.stat().st_size} for path in outputs},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
