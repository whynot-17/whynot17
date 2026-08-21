from __future__ import annotations

import csv
import gzip
import io
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "geo" / "raw"
OUT = ROOT / "phase1_fao_oxa_screen"
OUT.mkdir(parents=True, exist_ok=True)


GENE_SETS = {
    "carnitine_core": [
        "BBOX1", "SLC22A5", "CPT1A", "CPT1B", "CPT2", "SLC25A20"
    ],
    "fao_core": [
        "CPT1A", "CPT1B", "CPT2", "ACADM", "ACADVL", "HADHA", "HADHB",
        "ECHS1", "ETFA", "ETFB", "ETFDH", "PPARGC1A", "KLF5", "FABP6"
    ],
    "fao_carnitine_combined": [
        "BBOX1", "SLC22A5", "SLC25A20", "CPT1A", "CPT1B", "CPT2",
        "ACADM", "ACADVL", "HADHA", "HADHB", "ECHS1", "ETFA", "ETFB",
        "ETFDH", "PPARGC1A", "KLF5", "FABP6"
    ],
}

DATASETS = {
    "GSE77932": {
        "platform": "GPL16699",
        "groups": {
            "HCT116": {"parental": ["GSM2061631"], "resistant": ["GSM2061632", "GSM2061633", "GSM2061634"]},
            "DLD1": {"parental": ["GSM2061635"], "resistant": ["GSM2061636", "GSM2061637", "GSM2061638"]},
        },
        "note": "matched parental and independent OHP-resistant clones; no biological parental replicates",
    },
    "GSE42387": {
        "platform": "GPL16297",
        "groups": {
            "HCT116": {"parental": ["GSM1038651", "GSM1038652", "GSM1038653"], "resistant": ["GSM1038654", "GSM1038655", "GSM1038656"]},
            "HT29": {"parental": ["GSM1038660", "GSM1038661", "GSM1038662"], "resistant": ["GSM1038663", "GSM1038664", "GSM1038665"]},
            "LoVo": {"parental": ["GSM1038669", "GSM1038670", "GSM1038671"], "resistant": ["GSM1038672", "GSM1038673", "GSM1038674"]},
        },
        "note": "three CRC backgrounds with biological triplicates; parental and resistant cells cultured drug-free before profiling",
    },
    "GSE124808": {
        "platform": "GPL16699",
        "groups": {
            "HCT116_xenograft": {"parental": ["GSM3554324"], "resistant": ["GSM3554325", "GSM3554326"]},
        },
        "note": "xenograft tumors from parental versus two HCT/OHP-resistant clones; descriptive validation only",
    },
}

ENSEMBL_SYMBOL = {
    "ENSG00000129151": "BBOX1",
    "ENSG00000197375": "SLC22A5",
    "ENSG00000110090": "CPT1A",
    "ENSG00000205560": "CPT1B",
    "ENSG00000157184": "CPT2",
    "ENSG00000178537": "SLC25A20",
    "ENSG00000117054": "ACADM",
    "ENSG00000072778": "ACADVL",
    "ENSG00000084754": "HADHA",
    "ENSG00000138029": "HADHB",
    "ENSG00000127884": "ECHS1",
    "ENSG00000140374": "ETFA",
    "ENSG00000105379": "ETFB",
    "ENSG00000171503": "ETFDH",
    "ENSG00000109819": "PPARGC1A",
    "ENSG00000102554": "KLF5",
    "ENSG00000170231": "FABP6",
}

RNA_DATASETS = {
    "GSE119603": {
        "groups": {
            "HCT116": {"parental": ["HCT116_1", "HCT116_2", "HCT116_3"], "resistant": ["HCT116oxR_1", "HCT116oxR_2", "HCT116oxR_3"]},
        },
        "note": "count matrix for HCT116 parental versus HCT116oxR; gene identifiers are gene symbols; first-pass uses log2-CPM",
    }
}

AFFY_GENE_ID = {
    "BBOX1": "8424", "SLC22A5": "6584", "CPT1A": "1374", "CPT1B": "1375", "CPT2": "1376",
    "SLC25A20": "788", "ACADM": "34", "ACADVL": "37", "HADHA": "3030", "HADHB": "3032",
    "ECHS1": "1892", "ETFA": "2108", "ETFB": "2109", "ETFDH": "2110", "PPARGC1A": "10891",
    "KLF5": "688", "FABP6": "2172",
}

GPL2006_REFSEQ = {
    "NM_003986": "BBOX1",
    "NM_003060": "SLC22A5",
    "NM_001876": "CPT1A",
    "NM_000098": "CPT2",
    "NM_000387": "SLC25A20",
    "NM_000016": "ACADM",
    "NM_001730": "KLF5",
}


def read_platform(platform: str) -> pd.DataFrame:
    path = RAW / f"{platform}.soft.txt"
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_begin")) + 1
    end = next((i for i in range(start, len(lines)) if lines[i].startswith("!platform_table_end")), len(lines))
    rows = list(csv.reader(lines[start:end], delimiter="\t"))
    table = pd.DataFrame(rows[1:], columns=rows[0])
    return table


def symbol_column(platform: str, table: pd.DataFrame) -> str | None:
    for col in table.columns:
        if col.lower() in {"gene symbol", "gene_symbol", "gene symbol(s)", "gene_assignment"}:
            return col
    if platform == "GPL16699":
        return "GENE_SYMBOL"
    return None


def platform_mapping(platform: str) -> dict[str, list[str]]:
    # GPL16699 is the original Agilent feature-number design. GPL22628 is a
    # reannotated version of the same 039494 design with current Ensembl IDs.
    annotation_platform = "GPL22628" if platform == "GPL16699" and (RAW / "GPL22628.soft.txt").exists() else platform
    if platform == "GPL6244":
        sqlite_path = RAW.parent / "annotation" / "hugene10" / "hugene10sttranscriptcluster.db" / "inst" / "extdata" / "hugene10sttranscriptcluster.sqlite"
        import sqlite3
        mapping: dict[str, list[str]] = {}
        with sqlite3.connect(sqlite_path) as con:
            for symbol, gene_id in AFFY_GENE_ID.items():
                for row in con.execute("select probe_id from probes where gene_id=?", (gene_id,)):
                    mapping.setdefault(str(row[0]), []).append(symbol)
        return mapping
    table = read_platform(annotation_platform)
    if platform == "GPL2006":
        mapping: dict[str, list[str]] = {}
        for _, row in table.iterrows():
            probe = str(row[table.columns[0]]).strip()
            refseq = str(row.get("GB_ACC", "")).strip()
            symbol = GPL2006_REFSEQ.get(refseq)
            mapping[probe] = [symbol] if symbol else []
        return mapping
    id_col = table.columns[0]
    if annotation_platform == "GPL22628":
        mapping: dict[str, list[str]] = {}
        for _, row in table.iterrows():
            probe = str(row[id_col]).strip()
            gene_id = str(row.get("GENE_ID", "")).strip().split(".")[0]
            symbol = ENSEMBL_SYMBOL.get(gene_id)
            mapping[probe] = [symbol] if symbol else []
        return mapping
    col = symbol_column(platform, table)
    mapping: dict[str, list[str]] = {}
    if col is None:
        return mapping
    for _, row in table.iterrows():
        probe = str(row[id_col]).strip()
        raw = str(row[col]).strip()
        symbols: list[str] = []
        if platform == "GPL6244":
            # Not used in the current first-pass model datasets.
            parts = raw.split("///")
            symbols = [p.strip().split(" // ")[1] for p in parts if " // " in p and len(p.split(" // ")) > 1]
        else:
            for token in raw.replace("|", "///").replace(";", "///").split("///"):
                token = token.strip()
                if token and token not in {"---", "NA", "nan"}:
                    symbols.append(token.upper())
        mapping[probe] = sorted(set(symbols))
    return mapping


def read_series_matrix(accession: str) -> tuple[pd.DataFrame, dict[str, str]]:
    path = RAW / f"{accession}_series_matrix.txt.gz"
    metadata: dict[str, str] = {}
    data_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                data_lines.append(line)
            elif line.startswith("!Sample_geo_accession"):
                metadata["sample_geo_accession"] = line.rstrip("\n")
            elif line.startswith("!Sample_characteristics_ch1") and "1st line best response" in line:
                metadata["best_response"] = line.rstrip("\n")
            elif line.startswith("!Series_platform_id"):
                metadata["platform"] = line.rstrip("\n")
    frame = pd.read_csv(io.StringIO("".join(data_lines)), sep="\t", quotechar='"')
    frame = frame.rename(columns={frame.columns[0]: "ID_REF"})
    frame["ID_REF"] = frame["ID_REF"].astype(str)
    frame = frame.set_index("ID_REF")
    frame = frame.apply(pd.to_numeric, errors="coerce")
    return frame, metadata


def normalize_expression(frame: pd.DataFrame) -> pd.DataFrame:
    numeric = frame.to_numpy(dtype=float)
    finite = numeric[np.isfinite(numeric)]
    # GEO processed matrices in this project mix log-scale Agilent values with
    # positive, non-log signal values. Use a transparent dataset-level rule.
    if finite.size and np.nanpercentile(finite, 95) > 100:
        return np.log2(frame.clip(lower=0) + 1)
    return frame


def gene_expression(frame: pd.DataFrame, mapping: dict[str, list[str]]) -> pd.DataFrame:
    rows: dict[str, pd.Series] = {}
    for probe, symbols in mapping.items():
        if probe not in frame.index:
            continue
        for symbol in symbols:
            rows.setdefault(symbol, []).append(frame.loc[probe])
    result = {}
    for symbol, values in rows.items():
        result[symbol] = pd.concat(values, axis=1).median(axis=1)
    return pd.DataFrame(result).T


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    out = expr.copy()
    for gene in out.index:
        values = out.loc[gene]
        sd = values.std(ddof=1)
        out.loc[gene] = (values - values.mean()) / sd if sd and np.isfinite(sd) else 0.0
    return out


def summarize_dataset(accession: str, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix, _ = read_series_matrix(accession)
    matrix = normalize_expression(matrix)
    mapping = platform_mapping(config["platform"])
    expr = gene_expression(matrix, mapping)
    zexpr = zscore_rows(expr)
    rows = []
    module_rows = []
    for cell_line, groups in config["groups"].items():
        p = [x for x in groups["parental"] if x in expr.columns]
        r = [x for x in groups["resistant"] if x in expr.columns]
        for set_name, genes in GENE_SETS.items():
            present = [g for g in genes if g in zexpr.index]
            if not present or not p or not r:
                module_rows.append({"dataset": accession, "cell_line": cell_line, "module": set_name, "n_genes": len(present), "parental_n": len(p), "resistant_n": len(r), "delta_resistant_minus_parental": np.nan, "p_value": np.nan, "direction": "not_estimable"})
                continue
            parent = zexpr.loc[present, p].mean(axis=0)
            resistant = zexpr.loc[present, r].mean(axis=0)
            delta = resistant.mean() - parent.mean()
            pval = np.nan
            if len(parent) >= 2 and len(resistant) >= 2:
                pval = float(ttest_ind(resistant, parent, equal_var=False).pvalue)
            module_rows.append({"dataset": accession, "cell_line": cell_line, "module": set_name, "n_genes": len(present), "parental_n": len(p), "resistant_n": len(r), "delta_resistant_minus_parental": delta, "p_value": pval, "direction": "up" if delta > 0 else "down" if delta < 0 else "flat"})
        for gene in sorted(set(sum(GENE_SETS.values(), []))):
            if gene not in expr.index or not p or not r:
                continue
            parent = expr.loc[gene, p].astype(float)
            resistant = expr.loc[gene, r].astype(float)
            delta = resistant.mean() - parent.mean()
            pval = float(ttest_ind(resistant, parent, equal_var=False).pvalue) if len(parent) >= 2 and len(resistant) >= 2 else np.nan
            rows.append({"dataset": accession, "cell_line": cell_line, "gene": gene, "parental_n": len(p), "resistant_n": len(r), "mean_parental": parent.mean(), "mean_resistant": resistant.mean(), "delta_resistant_minus_parental": delta, "p_value": pval, "direction": "up" if delta > 0 else "down" if delta < 0 else "flat"})
    return pd.DataFrame(rows), pd.DataFrame(module_rows)


def summarize_count_dataset(accession: str, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = RAW / f"{accession}_all.counts.txt.gz"
    counts = pd.read_csv(path, sep="\t")
    counts = counts.rename(columns={counts.columns[0]: "GeneID"}).set_index("GeneID")
    counts.index = counts.index.astype(str).str.upper()
    libsize = counts.sum(axis=0)
    expr = np.log2(counts.divide(libsize, axis=1) * 1_000_000 + 1)
    zexpr = zscore_rows(expr)
    rows = []
    module_rows = []
    for cell_line, groups in config["groups"].items():
        p = [x for x in groups["parental"] if x in expr.columns]
        r = [x for x in groups["resistant"] if x in expr.columns]
        for set_name, genes in GENE_SETS.items():
            present = [g for g in genes if g in zexpr.index]
            if not present or not p or not r:
                module_rows.append({"dataset": accession, "cell_line": cell_line, "module": set_name, "n_genes": len(present), "parental_n": len(p), "resistant_n": len(r), "delta_resistant_minus_parental": np.nan, "p_value": np.nan, "direction": "not_estimable"})
                continue
            parent = zexpr.loc[present, p].mean(axis=0)
            resistant = zexpr.loc[present, r].mean(axis=0)
            delta = resistant.mean() - parent.mean()
            pval = float(ttest_ind(resistant, parent, equal_var=False).pvalue) if len(parent) >= 2 and len(resistant) >= 2 else np.nan
            module_rows.append({"dataset": accession, "cell_line": cell_line, "module": set_name, "n_genes": len(present), "parental_n": len(p), "resistant_n": len(r), "delta_resistant_minus_parental": delta, "p_value": pval, "direction": "up" if delta > 0 else "down" if delta < 0 else "flat"})
        for gene in sorted(set(sum(GENE_SETS.values(), []))):
            if gene not in expr.index or not p or not r:
                continue
            parent = expr.loc[gene, p].astype(float)
            resistant = expr.loc[gene, r].astype(float)
            delta = resistant.mean() - parent.mean()
            pval = float(ttest_ind(resistant, parent, equal_var=False).pvalue) if len(parent) >= 2 and len(resistant) >= 2 else np.nan
            rows.append({"dataset": accession, "cell_line": cell_line, "gene": gene, "parental_n": len(p), "resistant_n": len(r), "mean_parental": parent.mean(), "mean_resistant": resistant.mean(), "delta_resistant_minus_parental": delta, "p_value": pval, "direction": "up" if delta > 0 else "down" if delta < 0 else "flat"})
    return pd.DataFrame(rows), pd.DataFrame(module_rows)


def summarize_patient_response(accession: str) -> pd.DataFrame:
    matrix, metadata = read_series_matrix(accession)
    matrix = normalize_expression(matrix)
    mapping = platform_mapping("GPL6244")
    expr = gene_expression(matrix, mapping)
    accession_values = next(csv.reader([metadata["sample_geo_accession"]], delimiter="\t"))[1:]
    response_values = next(csv.reader([metadata["best_response"]], delimiter="\t"))[1:]
    groups = {acc: "responder" if "Responder" in response and "Non-Responder" not in response else "non_responder" if "Non-Responder" in response else "missing" for acc, response in zip(accession_values, response_values)}
    selected = [acc for acc in accession_values if groups[acc] in {"responder", "non_responder"} and acc in expr.columns]
    rows = []
    for set_name, genes in GENE_SETS.items():
        present = [g for g in genes if g in expr.index]
        if not present:
            rows.append({"dataset": accession, "module": set_name, "n_genes": 0, "responder_n": 0, "non_responder_n": 0, "mean_responder": np.nan, "mean_non_responder": np.nan, "delta_non_responder_minus_responder": np.nan, "p_value": np.nan, "direction": "not_estimable"})
            continue
        zexpr = zscore_rows(expr.loc[present, selected])
        responders = [acc for acc in selected if groups[acc] == "responder"]
        nonresponders = [acc for acc in selected if groups[acc] == "non_responder"]
        rscore = zexpr[responders].mean(axis=0)
        nscore = zexpr[nonresponders].mean(axis=0)
        delta = nscore.mean() - rscore.mean()
        pval = float(ttest_ind(nscore, rscore, equal_var=False).pvalue) if len(rscore) >= 2 and len(nscore) >= 2 else np.nan
        rows.append({"dataset": accession, "module": set_name, "n_genes": len(present), "responder_n": len(responders), "non_responder_n": len(nonresponders), "mean_responder": rscore.mean(), "mean_non_responder": nscore.mean(), "delta_non_responder_minus_responder": delta, "p_value": pval, "direction": "higher_in_non_responder" if delta > 0 else "lower_in_non_responder" if delta < 0 else "flat"})
    return pd.DataFrame(rows)


def summarize_cross_sectional_sensitivity(accession: str) -> pd.DataFrame:
    matrix, _ = read_series_matrix(accession)
    matrix = normalize_expression(matrix)
    mapping = platform_mapping("GPL2006")
    expr = gene_expression(matrix, mapping)
    lines = (RAW / f"{accession}_series_matrix.txt.gz")
    titles = []
    sensitivity = []
    accessions = []
    with gzip.open(lines, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_geo_accession"):
                accessions = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_characteristics_ch1") and "oxa sensitivity" in line:
                sensitivity = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
    selected = [acc for acc, title in zip(accessions, titles) if title.lower().startswith("control") and acc in expr.columns]
    labels = {acc: lab.split(":", 1)[-1].strip() for acc, lab, title in zip(accessions, sensitivity, titles) if title.lower().startswith("control")}
    resistant = [acc for acc in selected if labels.get(acc) == "OXA_resistant"]
    sensitive = [acc for acc in selected if labels.get(acc) == "OXA_sensitive"]
    rows = []
    for set_name, genes in GENE_SETS.items():
        present = [g for g in genes if g in expr.index]
        if not present:
            rows.append({"dataset": accession, "module": set_name, "n_genes": 0, "resistant_n": len(resistant), "sensitive_n": len(sensitive), "delta_resistant_minus_sensitive": np.nan, "p_value": np.nan, "direction": "not_estimable"})
            continue
        zexpr = zscore_rows(expr.loc[present, selected])
        rscore = zexpr[resistant].mean(axis=0)
        sscore = zexpr[sensitive].mean(axis=0)
        delta = rscore.mean() - sscore.mean()
        pval = float(ttest_ind(rscore, sscore, equal_var=False).pvalue) if len(rscore) >= 2 and len(sscore) >= 2 else np.nan
        rows.append({"dataset": accession, "module": set_name, "n_genes": len(present), "resistant_n": len(resistant), "sensitive_n": len(sensitive), "delta_resistant_minus_sensitive": delta, "p_value": pval, "direction": "higher_in_resistant" if delta > 0 else "lower_in_resistant" if delta < 0 else "flat"})
    return pd.DataFrame(rows)


def main() -> None:
    all_gene = []
    all_module = []
    for accession, config in DATASETS.items():
        genes, modules = summarize_dataset(accession, config)
        all_gene.append(genes)
        all_module.append(modules)
        genes.to_csv(OUT / f"{accession}_fao_gene_effects.csv", index=False)
        modules.to_csv(OUT / f"{accession}_fao_module_effects.csv", index=False)
    for accession, config in RNA_DATASETS.items():
        genes, modules = summarize_count_dataset(accession, config)
        all_gene.append(genes)
        all_module.append(modules)
        genes.to_csv(OUT / f"{accession}_fao_gene_effects.csv", index=False)
        modules.to_csv(OUT / f"{accession}_fao_module_effects.csv", index=False)
    patient = summarize_patient_response("GSE83129")
    patient.to_csv(OUT / "GSE83129_patient_response_fao_modules.csv", index=False)
    cross = summarize_cross_sectional_sensitivity("GSE30011")
    cross.to_csv(OUT / "GSE30011_cross_sectional_sensitivity_fao_modules.csv", index=False)
    pd.concat(all_gene, ignore_index=True).to_csv(OUT / "fao_gene_effects_all_datasets.csv", index=False)
    pd.concat(all_module, ignore_index=True).to_csv(OUT / "fao_module_effects_all_datasets.csv", index=False)
    summary = {
        "scope": "first-pass transcriptomic screen of FAO/carnitine gene and module direction in parental versus OXA-resistant CRC models",
        "datasets": {k: {"platform": v["platform"], "groups": v["groups"], "note": v["note"]} for k, v in DATASETS.items()},
        "gene_sets": GENE_SETS,
        "interpretation_boundary": "This screen tests transcriptomic direction, not biochemical FAO flux or causal dependency.",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
