from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "geo" / "raw"
GENESETS = ROOT / "gene_sets"
OUT = ROOT / "phase1_pathway_matrix"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phase1_fao_oxa_screen as base  # noqa: E402


PRIMARY_CONTRASTS = []
for accession, config in base.DATASETS.items():
    for cell_line, groups in config["groups"].items():
        PRIMARY_CONTRASTS.append(
            {
                "dataset": accession,
                "model": cell_line,
                "context": "acquired_OXA_R",
                "platform": config["platform"],
                "parental": groups["parental"],
                "resistant": groups["resistant"],
                "note": config["note"],
            }
        )
for accession, config in base.RNA_DATASETS.items():
    for cell_line, groups in config["groups"].items():
        PRIMARY_CONTRASTS.append(
            {
                "dataset": accession,
                "model": cell_line,
                "context": "acquired_OXA_R",
                "platform": "RNA_COUNTS",
                "parental": groups["parental"],
                "resistant": groups["resistant"],
                "note": config["note"],
            }
        )


HALLMARK = {}
REACTOME = {}


def read_gmt(path: Path) -> dict[str, set[str]]:
    result: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                result[fields[0]] = {x.upper() for x in fields[2:] if x}
    return result


def union_sets(*names: str, source: dict[str, set[str]]) -> set[str]:
    values: set[str] = set()
    for name in names:
        values |= source[name]
    return values


def build_pathways() -> tuple[dict[str, set[str]], dict[str, str]]:
    global HALLMARK, REACTOME
    HALLMARK = read_gmt(GENESETS / "h.all.v2026.1.Hs.symbols.gmt")
    REACTOME = read_gmt(GENESETS / "c2.cp.reactome.v2026.1.Hs.symbols.gmt")

    pathways: dict[str, set[str]] = {
        "GSH_redox": union_sets(
            "REACTOME_GLUTATHIONE_CONJUGATION",
            "REACTOME_GLUTATHIONE_SYNTHESIS_AND_RECYCLING",
            source=REACTOME,
        ),
        "NRF2_response": {
            "NFE2L2", "KEAP1", "NQO1", "HMOX1", "GCLC", "GCLM", "GSS",
            "SLC7A11", "TXNRD1", "SRXN1", "FTH1", "FTL", "ABCC1", "ABCC2",
            "AKR1C1", "AKR1C2", "AKR1C3", "PRDX1", "GPX2", "GPX4", "GSTP1",
            "GSTA1", "GSTM1",
        },
        "ferroptosis_resistance": {
            "SLC7A11", "GPX4", "AIFM2", "FSP1", "DHODH", "GCH1", "FTH1",
            "FTL", "NFE2L2", "ACSL3", "SCD", "LPCAT3", "ALOX15", "TFRC",
            "SAT1", "ATG5", "ATG7", "NCOA4", "SLC40A1",
        },
        "DNA_repair": HALLMARK["HALLMARK_DNA_REPAIR"],
        "UPR_ER_stress": HALLMARK["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"],
        "OXPHOS": HALLMARK["HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
        "glycolysis": HALLMARK["HALLMARK_GLYCOLYSIS"],
        "EMT": HALLMARK["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
        "TGF_beta": HALLMARK["HALLMARK_TGF_BETA_SIGNALING"],
        "ABC_transport": {g for g in REACTOME["REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT"] if g.startswith("ABC")},
        "drug_metabolism": HALLMARK["HALLMARK_XENOBIOTIC_METABOLISM"],
        "apoptosis": HALLMARK["HALLMARK_APOPTOSIS"],
        "autophagy": REACTOME["REACTOME_AUTOPHAGY"],
        "cholesterol_homeostasis": HALLMARK["HALLMARK_CHOLESTEROL_HOMEOSTASIS"],
        "lipid_FA_metabolism": union_sets(
            "HALLMARK_FATTY_ACID_METABOLISM",
            source=HALLMARK,
        ) | REACTOME["REACTOME_FATTY_ACID_METABOLISM"],
        "purine_metabolism": union_sets(
            "REACTOME_PURINE_CATABOLISM",
            "REACTOME_PURINE_RIBONUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHESIS",
            "REACTOME_PURINE_SALVAGE",
            source=REACTOME,
        ),
        "pyrimidine_metabolism": union_sets(
            "REACTOME_PYRIMIDINE_CATABOLISM",
            "REACTOME_PYRIMIDINE_SALVAGE",
            source=REACTOME,
        ),
        "ROS": HALLMARK["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"],
        "TNF_NFkB": HALLMARK["HALLMARK_TNFA_SIGNALING_VIA_NFKB"],
        "IL6_JAK_STAT3": HALLMARK["HALLMARK_IL6_JAK_STAT3_SIGNALING"],
        "FAO_mitochondrial": union_sets(
            "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
            source=REACTOME,
        ),
        "carnitine_entry": {"BBOX1", "SLC22A5", "CPT1A", "CPT1B", "CPT2", "SLC25A20"},
    }
    sources = {
        name: (
            "MSigDB Hallmark"
            if name in {"DNA_repair", "UPR_ER_stress", "OXPHOS", "glycolysis", "EMT", "TGF_beta", "drug_metabolism", "apoptosis", "cholesterol_homeostasis", "ROS", "TNF_NFkB", "IL6_JAK_STAT3"}
            else "MSigDB Reactome"
            if name not in {"NRF2_response", "ferroptosis_resistance", "carnitine_entry"}
            else "curated mechanistic set"
        )
        for name in pathways
    }
    return pathways, sources


def read_ncbi_maps() -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    path = GENESETS / "Homo_sapiens.gene_info.gz"
    symbol_to_entrez: dict[str, str] = {}
    ensembl_to_symbol: dict[str, str] = {}
    refseq_to_symbol: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            symbol = row.get("Symbol", "").upper()
            gene_id = row.get("GeneID", "")
            if not symbol or not gene_id or symbol.startswith("LOC"):
                continue
            symbol_to_entrez.setdefault(symbol, gene_id)
            for token in row.get("dbXrefs", "").split("|"):
                if token.startswith("Ensembl:"):
                    ensembl_to_symbol.setdefault(token.split(":", 1)[1].split(".")[0], symbol)
                elif token.startswith("EnsemblGene:"):
                    ensembl_to_symbol.setdefault(token.split(":", 1)[1].split(".")[0], symbol)
                elif token.startswith("RefSeq:"):
                    refseq_to_symbol.setdefault(token.split(":", 1)[1].split(".")[0], symbol)
    return symbol_to_entrez, ensembl_to_symbol, refseq_to_symbol


def dynamic_platform_mapping(platform: str, symbols: set[str]) -> dict[str, list[str]]:
    symbol_to_entrez, ensembl_to_symbol, refseq_to_symbol = read_ncbi_maps()
    wanted = {x.upper() for x in symbols}
    if platform == "RNA_COUNTS":
        return {g: [g] for g in wanted}
    if platform == "GPL6244":
        sqlite_path = RAW.parent / "annotation" / "hugene10" / "hugene10sttranscriptcluster.db" / "inst" / "extdata" / "hugene10sttranscriptcluster.sqlite"
        import sqlite3
        ids = {symbol_to_entrez[g] for g in wanted if g in symbol_to_entrez}
        mapping: dict[str, list[str]] = {}
        with sqlite3.connect(sqlite_path) as con:
            for gene_id in ids:
                symbol = next(g for g in wanted if symbol_to_entrez.get(g) == gene_id)
                for (probe_id,) in con.execute("select probe_id from probes where gene_id=?", (gene_id,)):
                    mapping.setdefault(str(probe_id), []).append(symbol)
        return mapping
    if platform == "GPL16699":
        table = base.read_platform("GPL22628")
        mapping: dict[str, list[str]] = {}
        for _, row in table.iterrows():
            probe = str(row["ID"]).strip()
            ids = re.findall(r"ENSG\d+", str(row.get("GENE_ID", "")))
            mapped = sorted({ensembl_to_symbol[x] for x in ids if x in ensembl_to_symbol and ensembl_to_symbol[x] in wanted})
            mapping[probe] = mapped
        return mapping
    if platform == "GPL2006":
        table = base.read_platform(platform)
        hgu_path = GENESETS / "hgu133a.db" / "hgu133a.db" / "inst" / "extdata" / "hgu133a.sqlite"
        import sqlite3
        wanted_ids = {symbol_to_entrez[g]: g for g in wanted if g in symbol_to_entrez}
        accession_to_symbols: dict[str, set[str]] = {}
        with sqlite3.connect(hgu_path) as con:
            query = """
                select a.accession, p.gene_id
                from accessions a join probes p on a.probe_id = p.probe_id
                where p.gene_id is not null
            """
            for accession, gene_id in con.execute(query):
                if str(gene_id) in wanted_ids:
                    accession_to_symbols.setdefault(str(accession).split(".")[0], set()).add(wanted_ids[str(gene_id)])
        mapping = {}
        for _, row in table.iterrows():
            probe = str(row["ID"]).strip()
            accessions = re.findall(r"(?:N[MR]_\d+(?:\.\d+)?)", str(row.get("GB_ACC", "")))
            mapped = sorted({symbol for x in accessions for symbol in accession_to_symbols.get(x.split(".")[0], set())})
            if not mapped:
                mapped = sorted({refseq_to_symbol[x.split(".")[0]] for x in accessions if x.split(".")[0] in refseq_to_symbol and refseq_to_symbol[x.split(".")[0]] in wanted})
            mapping[probe] = mapped
        return mapping
    table = base.read_platform(platform)
    symbol_col = base.symbol_column(platform, table)
    mapping = {}
    if symbol_col:
        for _, row in table.iterrows():
            probe = str(row[table.columns[0]]).strip()
            raw = str(row.get(symbol_col, ""))
            vals = [x.strip().upper() for x in re.split(r"[|;/]+", raw) if x.strip()]
            mapping[probe] = sorted(set(vals) & wanted)
    return mapping


def read_counts_expression(accession: str) -> pd.DataFrame:
    counts = pd.read_csv(RAW / f"{accession}_all.counts.txt.gz", sep="\t")
    counts = counts.rename(columns={counts.columns[0]: "GeneID"}).set_index("GeneID")
    counts.index = counts.index.astype(str).str.upper()
    libsize = counts.sum(axis=0)
    return np.log2(counts.divide(libsize, axis=1) * 1_000_000 + 1)


def expression_for(accession: str, platform: str, wanted: set[str]) -> pd.DataFrame:
    if platform == "RNA_COUNTS":
        return read_counts_expression(accession)
    matrix, _ = base.read_series_matrix(accession)
    matrix = base.normalize_expression(matrix)
    mapping = dynamic_platform_mapping(platform, wanted)
    return base.gene_expression(matrix, mapping)


def zscore_rows(expr: pd.DataFrame) -> pd.DataFrame:
    values = expr.astype(float)
    means = values.mean(axis=1)
    sds = values.std(axis=1, ddof=1).replace(0, np.nan)
    return values.sub(means, axis=0).div(sds, axis=0).fillna(0.0)


def score_contrast(expr: pd.DataFrame, contrast: dict, pathways: dict[str, set[str]], sources: dict[str, str]) -> list[dict]:
    zexpr = zscore_rows(expr)
    parent = [x for x in contrast["parental"] if x in zexpr.columns]
    resistant = [x for x in contrast["resistant"] if x in zexpr.columns]
    rows = []
    for pathway, genes in pathways.items():
        present = sorted(set(genes) & set(zexpr.index))
        if not present or not parent or not resistant:
            delta = np.nan
            pvalue = np.nan
            direction = "not_estimable"
        else:
            pscore = zexpr.loc[present, parent].mean(axis=0)
            rscore = zexpr.loc[present, resistant].mean(axis=0)
            delta = float(rscore.mean() - pscore.mean())
            pvalue = float(ttest_ind(rscore, pscore, equal_var=False).pvalue) if len(parent) >= 2 and len(resistant) >= 2 else np.nan
            direction = "up" if delta > 0 else "down" if delta < 0 else "flat"
        rows.append({
            "contrast_id": f"{contrast['dataset']}|{contrast['model']}",
            "dataset": contrast["dataset"],
            "model": contrast["model"],
            "context": contrast["context"],
            "platform": contrast["platform"],
            "pathway": pathway,
            "gene_set_source": sources[pathway],
            "n_genes_present": len(present),
            "parental_n": len(parent),
            "resistant_n": len(resistant),
            "delta_pathway": delta,
            "p_value": pvalue,
            "direction": direction,
            "note": contrast.get("note", ""),
        })
    return rows


def patient_contrast() -> dict:
    matrix, metadata = base.read_series_matrix("GSE83129")
    accession_values = next(csv.reader([metadata["sample_geo_accession"]], delimiter="\t"))[1:]
    response_values = next(csv.reader([metadata["best_response"]], delimiter="\t"))[1:]
    responders = [acc for acc, response in zip(accession_values, response_values) if "Responder" in response and "Non-Responder" not in response]
    nonresponders = [acc for acc, response in zip(accession_values, response_values) if "Non-Responder" in response]
    return {
        "dataset": "GSE83129",
        "model": "patient_nonresponder_vs_responder",
        "context": "clinical_OXA_response",
        "platform": "GPL6244",
        "parental": responders,
        "resistant": nonresponders,
        "note": "clinical non-responder versus responder; external validation, not acquired cellular resistance",
    }


def cross_sectional_contrast() -> dict:
    path = RAW / "GSE30011_series_matrix.txt.gz"
    titles: list[str] = []
    labels: list[str] = []
    accessions: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_geo_accession"):
                accessions = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_characteristics_ch1") and "oxa sensitivity" in line:
                labels = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
    selected = [(acc, lab.split(":", 1)[-1].strip()) for acc, lab, title in zip(accessions, labels, titles) if title.lower().startswith("control")]
    return {
        "dataset": "GSE30011",
        "model": "cross_sectional_OXA_R_vs_sensitive",
        "context": "cross_sectional_OXA_sensitivity",
        "platform": "GPL2006",
        "parental": [acc for acc, label in selected if label == "OXA_sensitive"],
        "resistant": [acc for acc, label in selected if label == "OXA_resistant"],
        "note": "cross-sectional OXA-resistant versus OXA-sensitive CRC cell lines; not an acquired-resistance contrast",
    }


def summarize(long: pd.DataFrame) -> pd.DataFrame:
    primary = long[long["context"] == "acquired_OXA_R"].copy()
    rows = []
    for pathway, group in primary.groupby("pathway"):
        vals = group["delta_pathway"].dropna()
        positive = int((vals > 0).sum())
        negative = int((vals < 0).sum())
        n = int(vals.size)
        rows.append({
            "pathway": pathway,
            "n_primary_contrasts": n,
            "n_up": positive,
            "n_down": negative,
            "sign_concordance": max(positive, negative) / n if n else np.nan,
            "median_delta": float(vals.median()) if n else np.nan,
            "mean_delta": float(vals.mean()) if n else np.nan,
            "median_abs_delta": float(vals.abs().median()) if n else np.nan,
            "gene_set_source": group["gene_set_source"].iloc[0],
        })
    return pd.DataFrame(rows).sort_values(["sign_concordance", "median_abs_delta"], ascending=[False, False])


def main() -> None:
    pathways, sources = build_pathways()
    wanted = set().union(*pathways.values())
    contrasts = list(PRIMARY_CONTRASTS)
    contrasts += [patient_contrast(), cross_sectional_contrast()]
    expression_cache: dict[tuple[str, str], pd.DataFrame] = {}
    all_rows: list[dict] = []
    for contrast in contrasts:
        key = (contrast["dataset"], contrast["platform"])
        if key not in expression_cache:
            expression_cache[key] = expression_for(contrast["dataset"], contrast["platform"], wanted)
        all_rows.extend(score_contrast(expression_cache[key], contrast, pathways, sources))
    long = pd.DataFrame(all_rows)
    primary = long[long["context"] == "acquired_OXA_R"].copy()
    matrix = primary.pivot(index="contrast_id", columns="pathway", values="delta_pathway").reset_index()
    matrix = matrix.merge(primary[["contrast_id", "dataset", "model", "context", "platform"]].drop_duplicates(), on="contrast_id", how="left")
    cols = ["contrast_id", "dataset", "model", "context", "platform"] + [x for x in pathways if x in matrix.columns]
    matrix[cols].to_csv(OUT / "pathway_by_dataset_matrix_primary.csv", index=False)
    cell_line_matrix = primary[primary["model"] != "HCT116_xenograft"].pivot(index="contrast_id", columns="pathway", values="delta_pathway").reset_index()
    cell_line_matrix = cell_line_matrix.merge(primary[["contrast_id", "dataset", "model", "context", "platform"]].drop_duplicates(), on="contrast_id", how="left")
    cell_line_matrix[cols].to_csv(OUT / "pathway_by_dataset_matrix_cell_lines.csv", index=False)
    long.to_csv(OUT / "pathway_effects_long_all_contexts.csv", index=False)
    summary = summarize(long)
    summary.to_csv(OUT / "pathway_stability_summary_primary.csv", index=False)
    cell_line_long = long[(long["context"] == "acquired_OXA_R") & (long["model"] != "HCT116_xenograft")].copy()
    summarize(cell_line_long).to_csv(OUT / "pathway_stability_summary_cell_lines.csv", index=False)
    secondary = long[long["context"] != "acquired_OXA_R"].pivot(index="contrast_id", columns="pathway", values="delta_pathway").reset_index()
    if not secondary.empty:
        secondary.to_csv(OUT / "pathway_by_dataset_matrix_secondary.csv", index=False)
    manifest = {
        "primary_contrasts": [x["contrast_id"] for x in []],
        "n_primary_contrasts": int(primary["contrast_id"].nunique()),
        "n_pathways": len(pathways),
        "pathway_sizes": {k: len(v) for k, v in pathways.items()},
        "standardization": "within dataset, gene-wise z-score; pathway score is mean z-score over present genes; delta is resistant minus parental",
        "caution": "small sample sizes and platform heterogeneity make this a discovery matrix; sign consistency is descriptive until independent validation",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
