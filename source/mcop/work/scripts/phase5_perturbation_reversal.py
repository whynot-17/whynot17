from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "geo" / "raw"
WORK = ROOT / "phase5_perturbation_reversal"
CACHE = WORK / "signature_cache"
WORK.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase1_fao_oxa_screen as base  # noqa: E402
import phase1_pathway_matrix as p1  # noqa: E402
import phase3_module_decomposition as p3  # noqa: E402


METADATA_URL = "https://maayanlab.cloud/sigcom-lincs/metadata-api/signatures/find"
DRUGS = {
    "leflunomide": {"line": "DHODH/pyrimidine", "role": "primary"},
    "teriflunomide": {"line": "DHODH/pyrimidine", "role": "primary"},
    "bortezomib": {"line": "proteasome/proteostasis", "role": "mechanistic_positive_control"},
    "meldonium": {"line": "carnitine/FAO", "role": "appendix_failed_candidate"},
}
MAX_SIGNATURES = {
    "bortezomib": 24,
    "meldonium": 24,
}

PRIMARY = [x for x in p1.PRIMARY_CONTRASTS if x["model"] != "HCT116_xenograft"]
CONTRAST_IDS = [f"{x['dataset']}|{x['model']}" for x in PRIMARY]


def parse_symbols(raw: object) -> list[str]:
    if raw is None:
        return []
    text = str(raw)
    if text.lower() in {"", "nan", "na", "---"}:
        return []
    return sorted({x.strip().upper() for x in re.split(r"[|;/]+", text) if x.strip() and x.strip() not in {"---", "NA"}})


def full_platform_mapping(platform: str) -> dict[str, list[str]]:
    table = base.read_platform(platform)
    mapping: dict[str, list[str]] = {}
    if platform == "GPL16699":
        col = "GENE_SYMBOL"
    elif platform == "GPL16297":
        col = "Gene Symbol"
    else:
        col = base.symbol_column(platform, table)
    if col is None:
        raise ValueError(f"No symbol annotation found for {platform}")
    for _, row in table.iterrows():
        probe = str(row.iloc[0]).strip()
        mapping[probe] = parse_symbols(row.get(col, ""))
    return mapping


def full_expression(accession: str, platform: str) -> pd.DataFrame:
    if platform == "RNA_COUNTS":
        expr = p1.read_counts_expression(accession)
        expr.index = expr.index.astype(str).str.upper()
        expr = expr[~expr.index.duplicated(keep="first")]
        return expr
    matrix, _ = base.read_series_matrix(accession)
    matrix = base.normalize_expression(matrix)
    return base.gene_expression(matrix, full_platform_mapping(platform))


def model_delta(expr: pd.DataFrame, contrast: dict) -> pd.Series:
    zexpr = p1.zscore_rows(expr)
    parent = [x for x in contrast["parental"] if x in zexpr.columns]
    resistant = [x for x in contrast["resistant"] if x in zexpr.columns]
    if not parent or not resistant:
        raise ValueError(f"Samples missing for {contrast}")
    delta = zexpr[resistant].mean(axis=1) - zexpr[parent].mean(axis=1)
    delta.name = f"{contrast['dataset']}|{contrast['model']}"
    return delta


def build_gene_delta_matrix() -> pd.DataFrame:
    cache: dict[tuple[str, str], pd.DataFrame] = {}
    series: list[pd.Series] = []
    for contrast in PRIMARY:
        key = (contrast["dataset"], contrast["platform"])
        if key not in cache:
            print(f"Loading expression: {key}", flush=True)
            cache[key] = full_expression(*key)
        series.append(model_delta(cache[key], contrast))
    matrix = pd.concat(series, axis=1)
    matrix.index = matrix.index.astype(str).str.upper()
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    matrix.to_csv(WORK / "gene_delta_matrix_primary.csv")
    return matrix


def fetch_metadata(drug: str) -> list[dict]:
    if drug == "meldonium":
        # Keep Meldonium as an explicit appendix comparator. It is not allowed
        # to re-enter the primary ranking without a matched CRC perturbation
        # signature from the same source.
        return []
    # Prefer already downloaded, provenance-preserving signatures. This also
    # makes the analysis reproducible when the remote metadata service is
    # temporarily rate-limited.
    cached = sorted(CACHE.glob(f"*_{drug}_*.tsv"))
    if cached:
        records = []
        for path in cached:
            dose = re.search(r"_(\d+(?:\.\d+)?uM)\.tsv$", path.name)
            records.append({"id": path.stem, "meta": {"local_id": path.stem, "cell_line": "HT29", "pert_time": "24 h", "pert_dose": dose.group(1) if dose else "", "persistent_id": path.name, "data_level": 5, "tissue": "intestine", "disease": "colon adenocarcinoma"}})
        limit = MAX_SIGNATURES.get(drug)
        return records[:limit] if limit else records
    payload = {"filter": {"where": {"meta.pert_name": drug}, "limit": 1000}}
    response = requests.post(METADATA_URL, json=payload, timeout=300)
    response.raise_for_status()
    records = response.json()
    # Use the CRC-relevant HT29 background and a common 24 h time point.
    selected = [r for r in records if r.get("meta", {}).get("cell_line") == "HT29" and r.get("meta", {}).get("pert_time") == "24 h"]
    if not selected:
        selected = [r for r in records if r.get("meta", {}).get("cell_line") == "HT29"]
    # Bortezomib has many near-duplicate HT29 signatures. Retain a deterministic
    # dose-spread subset for the positive-control line so it cannot dominate the
    # comparison simply by having more deposited replicates.
    limit = MAX_SIGNATURES.get(drug)
    if limit and len(selected) > limit:
        selected = sorted(selected, key=lambda r: (str(r.get("meta", {}).get("pert_dose", "")), str(r.get("meta", {}).get("local_id", ""))))[:limit]
    return selected


def read_signature(url: str) -> pd.Series:
    filename = CACHE / Path(url).name
    if not filename.exists():
        response = requests.get(url, timeout=45)
        response.raise_for_status()
        filename.write_bytes(response.content)
    frame = pd.read_csv(filename, sep="\t")
    frame.columns = [str(x).strip() for x in frame.columns]
    symbol_col = frame.columns[0]
    value_col = frame.columns[1]
    frame[symbol_col] = frame[symbol_col].astype(str).str.upper().str.strip()
    frame[value_col] = pd.to_numeric(frame[value_col], errors="coerce")
    frame = frame.dropna(subset=[symbol_col, value_col])
    frame = frame[frame[symbol_col].str.len() > 0]
    return frame.groupby(symbol_col)[value_col].median()


def select_genes(disease: pd.Series, n_each: int = 250, min_abs: float = 0.0) -> list[str]:
    d = disease.dropna()
    d = d[d.abs() >= min_abs]
    up = list(d.nlargest(min(n_each, len(d))).index)
    down = list(d.nsmallest(min(n_each, len(d))).index)
    return list(dict.fromkeys(up + down))


def rrs(disease: pd.Series, drug: pd.Series, genes: list[str] | None = None) -> tuple[float, float, int]:
    d = disease
    p = drug
    if genes is not None:
        d = d.reindex(genes)
        p = p.reindex(genes)
    joined = pd.concat([d.rename("disease"), p.rename("drug")], axis=1).dropna()
    if len(joined) < 20:
        return np.nan, np.nan, int(len(joined))
    rho = float(spearmanr(joined["disease"], joined["drug"]).statistic)
    return float((1.0 - rho) / 2.0), rho, int(len(joined))


def consensus(matrix: pd.DataFrame, ids: list[str]) -> pd.Series:
    return matrix.reindex(columns=ids).median(axis=1, skipna=True)


def load_subtypes() -> pd.DataFrame:
    path = ROOT / "phase4_vulnerability_mapping" / "model_vulnerability_cooccurrence.csv"
    state = pd.read_csv(path).set_index("contrast_id")
    state.index = state.index.astype(str)
    return state.reindex(CONTRAST_IDS)


def score_signatures(gene_matrix: pd.DataFrame, metadata: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    global_disease = consensus(gene_matrix, CONTRAST_IDS)
    global_genes = select_genes(global_disease, 250)
    subtype_state = load_subtypes()
    subtype_rules = {
        "salvage_low_DHODH_high": subtype_state["salvage_low_DHODH_high"].fillna(False).astype(bool),
        "salvage_low_RRM2_high": subtype_state["salvage_low_RRM2_high"].fillna(False).astype(bool),
        "UPR_low_ERAD_high": subtype_state["UPR_low_ERAD_high"].fillna(False).astype(bool),
    }
    results: list[dict] = []
    per_model: list[dict] = []
    for drug, drug_records in metadata.items():
        for record in drug_records:
            meta = record.get("meta", {})
            try:
                sig = read_signature(meta["persistent_id"])
            except Exception:
                continue
            sig_id = meta.get("local_id", record.get("id", ""))
            score, rho, n = rrs(global_disease, sig, global_genes)
            results.append({
                "drug": drug, "role": DRUGS[drug]["role"], "drug_line": DRUGS[drug]["line"],
                "signature_id": sig_id, "cell_line": meta.get("cell_line"),
                "pert_time": meta.get("pert_time"), "pert_dose": meta.get("pert_dose"),
                "rrs_global": score, "rho_global": rho, "n_genes_scored": n,
            })
            for model_id in CONTRAST_IDS:
                # Per-model scoring uses that model's own ranked resistance
                # genes. Using the consensus top list here would unfairly
                # penalize platform-specific coverage and made several
                # models non-estimable.
                model_genes = select_genes(gene_matrix[model_id], 250)
                score_m, rho_m, n_m = rrs(gene_matrix[model_id], sig, model_genes)
                per_model.append({"drug": drug, "role": DRUGS[drug]["role"], "signature_id": sig_id, "contrast_id": model_id, "rrs_model": score_m, "rho_model": rho_m, "n_genes_scored": n_m})
    return pd.DataFrame(results), pd.DataFrame(per_model)


def build_subtype_tables(gene_matrix: pd.DataFrame, signatures: dict[str, pd.Series], metadata_rows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    state = load_subtypes()
    rules = {
        "global": pd.Series(True, index=CONTRAST_IDS),
        "salvage_low_DHODH_high": state["salvage_low_DHODH_high"].fillna(False).astype(bool),
        "salvage_low_RRM2_high": state["salvage_low_RRM2_high"].fillna(False).astype(bool),
        "UPR_low_ERAD_high": state["UPR_low_ERAD_high"].fillna(False).astype(bool),
    }
    target_rows: list[dict] = []
    selectivity_rows: list[dict] = []
    target_modules = p3.build_decomposition_sets()[0]
    network_genes = set().union(
        target_modules["pyrimidine_salvage"], target_modules["pyrimidine_de_novo_core"],
        target_modules["pyrimidine_interconversion"], target_modules["PERK_eIF2A_ATF4"],
        target_modules["IRE1_XBP1"], target_modules["ATF6_proteostasis"],
        target_modules["mesenchymal_state"],
        {"DERL1", "EDEM1", "HSP90B1", "MANF", "VCP", "SEL1L", "SYVN1"},
    )
    network_genes = sorted(network_genes)
    network_no_dhodh = [g for g in network_genes if g != "DHODH"]
    for drug, sig in signatures.items():
        for subtype, mask in rules.items():
            ids = [x for x in CONTRAST_IDS if bool(mask.get(x, False))]
            if subtype == "global":
                ids = CONTRAST_IDS
            if len(ids) < 2:
                continue
            d = consensus(gene_matrix, ids)
            top = select_genes(d, 250)
            score, rho, n = rrs(d, sig, top)
            net_score, net_rho, net_n = rrs(d, sig, network_genes)
            no_score, no_rho, no_n = rrs(d, sig, network_no_dhodh)
            target_rows.append({"drug": drug, "role": DRUGS[drug]["role"], "subtype": subtype, "n_models": len(ids), "models": ";".join(ids), "rrs_top_signature": score, "rho_top_signature": rho, "n_genes_top": n, "rrs_target_network": net_score, "rho_target_network": net_rho, "n_genes_target_network": net_n, "rrs_target_network_no_DHODH": no_score, "rho_target_network_no_DHODH": no_rho, "n_genes_target_network_no_DHODH": no_n})
            if subtype != "global":
                other = [x for x in CONTRAST_IDS if x not in ids]
                if len(other) >= 2:
                    d_other = consensus(gene_matrix, other)
                    other_top = select_genes(d_other, 250)
                    other_score, other_rho, other_n = rrs(d_other, sig, other_top)
                    selectivity_rows.append({"drug": drug, "target_subtype": subtype, "target_n_models": len(ids), "matched_negative_n_models": len(other), "target_rrs": score, "matched_negative_rrs": other_score, "selectivity_difference": score - other_score if pd.notna(score) and pd.notna(other_score) else np.nan, "selectivity_ratio": (score + 0.01) / (other_score + 0.01) if pd.notna(score) and pd.notna(other_score) else np.nan, "target_models": ";".join(ids), "matched_negative_models": ";".join(other)})
    target_df = pd.DataFrame(target_rows)
    selectivity_df = pd.DataFrame(selectivity_rows)
    loo_rows: list[dict] = []
    for subtype, mask in rules.items():
        if subtype == "global":
            continue
        ids = [x for x in CONTRAST_IDS if bool(mask.get(x, False))]
        if len(ids) < 3:
            continue
        datasets = sorted({x.split("|", 1)[0] for x in ids})
        for left_out in datasets:
            keep = [x for x in ids if not x.startswith(left_out + "|")]
            if len(keep) < 2:
                loo_rows.append({"subtype": subtype, "left_out_dataset": left_out, "n_models_remaining": len(keep), "status": "not_estimable"})
                continue
            d = consensus(gene_matrix, keep)
            top = select_genes(d, 250)
            for drug, sig in signatures.items():
                score, rho, n = rrs(d, sig, top)
                loo_rows.append({"subtype": subtype, "left_out_dataset": left_out, "n_models_remaining": len(keep), "status": "estimable", "drug": drug, "rrs": score, "rho": rho, "n_genes": n, "models_remaining": ";".join(keep)})
    return target_df, selectivity_df, pd.DataFrame(loo_rows)


def main() -> None:
    gene_matrix = build_gene_delta_matrix()
    metadata: dict[str, list[dict]] = {}
    metadata_rows: list[dict] = []
    signatures: dict[str, pd.Series] = {}
    for drug in DRUGS:
        try:
            records = fetch_metadata(drug)
        except Exception as exc:
            print(f"Metadata failed for {drug}: {exc}", flush=True)
            records = []
        metadata[drug] = records
        for record in records:
            meta = record.get("meta", {})
            metadata_rows.append({"drug": drug, "role": DRUGS[drug]["role"], "drug_line": DRUGS[drug]["line"], "signature_id": meta.get("local_id", record.get("id", "")), "cell_line": meta.get("cell_line"), "pert_time": meta.get("pert_time"), "pert_dose": meta.get("pert_dose"), "persistent_id": meta.get("persistent_id"), "data_level": meta.get("data_level"), "tissue": meta.get("tissue"), "disease": meta.get("disease")})
            try:
                signatures[f"{drug}::{meta.get('local_id', record.get('id', ''))}"] = read_signature(meta["persistent_id"])
            except Exception as exc:
                print(f"Signature failed for {drug}: {exc}", flush=True)
        metadata[drug] = [
            r for r in records
            if f"{drug}::{r.get('meta', {}).get('local_id', r.get('id', ''))}" in signatures
        ]
    pd.DataFrame(metadata_rows).to_csv(WORK / "perturbation_signature_metadata.csv", index=False)
    global_df, per_model_df = score_signatures(gene_matrix, metadata)
    global_df.to_csv(WORK / "per_signature_global_rrs.csv", index=False)
    per_model_df.to_csv(WORK / "per_model_drug_rrs.csv", index=False)
    # Collapse replicate/dose signatures at drug level; median is deliberately robust to dose.
    signature_lookup = {drug: [] for drug in DRUGS}
    for key, sig in signatures.items():
        signature_lookup[key.split("::", 1)[0]].append(sig)
    drug_consensus = {drug: pd.concat(vals, axis=1).median(axis=1) for drug, vals in signature_lookup.items() if vals}
    target_df, selectivity_df, loo_df = build_subtype_tables(gene_matrix, drug_consensus, global_df)
    target_df.to_csv(WORK / "drug_rrs_global_subtype_and_remove_dhodh.csv", index=False)
    selectivity_df.to_csv(WORK / "drug_subtype_selectivity.csv", index=False)
    loo_df.to_csv(WORK / "leave_one_dataset_out.csv", index=False)
    manifest = {"n_models": len(CONTRAST_IDS), "models": CONTRAST_IDS, "drug_metadata_records": {k: len(v) for k, v in metadata.items()}, "signature_selection": "HT29, 24 h when available; otherwise HT29 all times", "rrs_definition": "(1 - Spearman(disease_delta, drug_CD_coefficient))/2; higher is stronger reversal", "top_genes": "250 most up plus 250 most down disease genes", "target_network_remove_DHODH": "DHODH removed from the curated compensatory network, not from the whole transcriptome top-gene score"}
    (WORK / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
