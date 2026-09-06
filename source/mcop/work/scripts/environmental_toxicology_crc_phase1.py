"""CTD x GeneCards CRC environmental-chemical enrichment screen.

The script expects GeneCards result-page exports supplied by the user. It does
not scrape GeneCards cards or perform an external CRC literature search.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import fisher_exact, hypergeom


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RULES = ROOT / "work" / "environmental_toxicology_crc_phase1" / "chemical_class_rules.json"


def canonical_chemical_id(value: str) -> str:
    value = str(value).strip()
    return value[5:] if value.startswith("MESH:") else value


def read_ctd_tsv(path: Path) -> pd.DataFrame:
    # CTD writes the real field line as a commented line immediately after
    # "# Fields:". Pandas' comment= would otherwise use the first chemical as
    # the header, so recover the field line explicitly.
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("# Fields:"):
                header_line = next(handle).lstrip("# ").rstrip("\r\n")
                columns = header_line.split("\t")
                return pd.read_csv(
                    handle,
                    sep="\t",
                    names=columns,
                    dtype=str,
                    keep_default_na=False,
                    comment="#",
                    low_memory=False,
                )
    raise ValueError(f"Could not find '# Fields:' header in CTD file: {path}")


def split_pipe(value: str) -> list[str]:
    if not value:
        return []
    return [x for x in value.split("|") if x]


def tree_paths(row: pd.Series) -> list[str]:
    paths: list[str] = []
    for col in ("TreeNumbers", "ParentTreeNumbers"):
        for value in split_pipe(row.get(col, "")):
            paths.append(value.split("/", 1)[0])
    return paths


def descendant(path: str, prefix: str) -> bool:
    return path == prefix or path.startswith(prefix + ".")


def classify_chemicals(
    chemicals: pd.DataFrame,
    rules_path: Path,
    drugcentral_path: Path | None = None,
    pah_formula_path: Path | None = None,
) -> pd.DataFrame:
    rules = json.loads(rules_path.read_text(encoding="utf-8"))
    core_rules: dict[str, list[str]] = rules["core_rules"]
    drug_exclusion = re.compile(rules.get("drug_exclusion_regex", r"(?!)"))
    drugcentral_cas: set[str] = set()
    drugcentral_names: set[str] = set()
    if drugcentral_path and drugcentral_path.exists():
        drugcentral = pd.read_csv(drugcentral_path, sep="\t", dtype=str, keep_default_na=False)
        if "CAS_RN" in drugcentral.columns:
            drugcentral_cas = {
                str(value).strip().lower()
                for value in drugcentral["CAS_RN"]
                if str(value).strip()
            }
        if "INN" in drugcentral.columns:
            drugcentral_names = {
                re.sub(r"[^a-z0-9]+", "", str(value).lower())
                for value in drugcentral["INN"]
                if str(value).strip()
            }
    pah_formulas = {}
    if pah_formula_path and pah_formula_path.exists():
        pah_formulas = json.loads(pah_formula_path.read_text(encoding="utf-8"))
    records = []
    has_drugbank = any("drugbank" in c.lower() for c in chemicals.columns)

    for row in chemicals.to_dict("records"):
        paths = []
        for col in ("TreeNumbers", "ParentTreeNumbers"):
            for value in split_pipe(row.get(col, "")):
                paths.append(value.split("/", 1)[0])
        matched: list[str] = []
        matched_prefixes: list[str] = []
        for category, prefixes in core_rules.items():
            hits = [prefix for prefix in prefixes if any(descendant(path, prefix) for path in paths)]
            if hits:
                matched.append(category)
                matched_prefixes.extend(hits)
        drug_text = " ".join(
            str(row.get(column, ""))
            for column in ("ChemicalName", "Definition", "MESHSynonyms", "CTDCuratedSynonyms")
        )
        drug_like_exclusion = bool(drug_exclusion.search(drug_text))
        cas = str(row.get("CasRN", "")).strip().lower()
        canonical_name = re.sub(r"[^a-z0-9]+", "", str(row.get("ChemicalName", "")).lower())
        drugcentral_match = (cas and cas in drugcentral_cas) or (
            canonical_name and canonical_name in drugcentral_names
        )
        pubchem_cid = re.sub(r"^CID:", "", str(row.get("PubChemCID", "")).strip())
        pah_formula = str(pah_formulas.get(pubchem_cid, ""))
        pah_structure_valid = bool(re.fullmatch(r"C[0-9]+H[0-9]+", pah_formula)) if pah_formula_path else pd.NA
        if drug_like_exclusion:
            matched = []
            matched_prefixes = []
        if drugcentral_match:
            matched = []
            matched_prefixes = []
        if "pahs" in matched and pah_formula_path and not pah_structure_valid:
            matched = [category for category in matched if category != "pahs"]
            matched_prefixes = [prefix for prefix in matched_prefixes if prefix != "D02.455.426.559.847"]
        record = {
            "ChemicalID": canonical_chemical_id(row.get("ChemicalID", "")),
            "ChemicalName": row.get("ChemicalName", ""),
            "CasRN": row.get("CasRN", ""),
            "PubChemCID": row.get("PubChemCID", ""),
            "DTXSID": row.get("DTXSID", ""),
            "InChIKey": row.get("InChIKey", ""),
            "Definition": row.get("Definition", ""),
            "ParentIDs": row.get("ParentIDs", ""),
            "TreeNumbers": row.get("TreeNumbers", ""),
            "ParentTreeNumbers": row.get("ParentTreeNumbers", ""),
            "chemical_class": ";".join(sorted(set(matched))),
            "classification_rule_prefixes": ";".join(sorted(set(matched_prefixes))),
            "is_core": bool(matched),
            "drug_like_exclusion": drug_like_exclusion,
            "drugcentral_match": bool(drugcentral_match),
            "pah_formula": pah_formula,
            "pah_structure_valid": pah_structure_valid,
            "classification_exclusion_reason": (
                "drugcentral_cas_or_inn" if drugcentral_match
                else "drug_semantic_regex" if drug_like_exclusion else ""
            ),
            "drugbank_field_available": has_drugbank,
        }
        records.append(record)

    out = pd.DataFrame(records).drop_duplicates("ChemicalID", keep="first")
    return out


def pmid_set(value: str) -> set[str]:
    if not value:
        return set()
    return {x.strip() for x in re.split(r"[|;,]", value) if x.strip()}


def load_interactions(path: Path) -> tuple[pd.DataFrame, dict[str, set[str]], dict[str, str]]:
    raw = read_ctd_tsv(path)
    required = {"ChemicalID", "ChemicalName", "GeneID", "GeneSymbol", "Organism", "PubMedIDs"}
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"CTD interaction file is missing columns: {missing}")

    human = raw[(raw["Organism"] == "Homo sapiens") & raw["ChemicalID"].ne("") & raw["GeneID"].ne("")].copy()
    human["GeneID"] = human["GeneID"].astype(str).str.strip()
    human["ChemicalID"] = human["ChemicalID"].astype(str).map(canonical_chemical_id)

    pair = human[["ChemicalID", "ChemicalName", "GeneID", "GeneSymbol"]].drop_duplicates(
        ["ChemicalID", "GeneID"]
    )
    gene_sets = pair.groupby("ChemicalID")["GeneID"].agg(lambda x: set(x)).to_dict()
    gene_symbols = (
        pair[pair["GeneSymbol"].ne("")]
        .drop_duplicates("GeneID")
        .set_index("GeneID")["GeneSymbol"]
        .to_dict()
    )

    pmids: dict[str, set[str]] = defaultdict(set)
    for chemical_id, value in human[["ChemicalID", "PubMedIDs"]].itertuples(index=False):
        pmids[chemical_id].update(pmid_set(value))

    raw_counts = human.groupby("ChemicalID").size().rename("n_raw_interaction_rows")
    pair_counts = pair.groupby("ChemicalID").size().rename("n_unique_chemical_gene_pairs")
    count_table = pd.concat([raw_counts, pair_counts], axis=1).fillna(0).reset_index()
    count_table["n_unique_pmids"] = count_table["ChemicalID"].map(lambda x: len(pmids.get(x, set())))
    return count_table, gene_sets, gene_symbols


def normalise_column(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def find_column(columns: list[str], candidates: list[str]) -> str | None:
    norm = {normalise_column(c): c for c in columns}
    for candidate in candidates:
        if normalise_column(candidate) in norm:
            return norm[normalise_column(candidate)]
    return None


def load_genecards(path: Path, gene_symbol_to_ids: dict[str, set[str]]) -> tuple[pd.DataFrame, dict[str, int]]:
    sample = path.read_text(encoding="utf-8-sig", errors="replace")[:8192]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        sep = dialect.delimiter
    except csv.Error:
        sep = "\t" if "\t" in sample else ","
    gc = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, engine="python")
    symbol_col = find_column(list(gc.columns), ["Symbol", "Gene Symbol", "Gene"])
    if not symbol_col:
        raise ValueError(f"Cannot find GeneCards symbol column in {path}")
    rank_col = find_column(list(gc.columns), ["#", "Rank", "Ranking", "GeneCards Rank"])
    score_col = find_column(list(gc.columns), ["Relevance Score", "Score"])
    type_col = find_column(list(gc.columns), ["Type", "Gene Type"])
    knowledge_col = find_column(list(gc.columns), ["Knowledge", "Knowledge Score", "GIFtS"])
    id_col = find_column(list(gc.columns), ["NCBI Gene ID", "Entrez Gene ID", "GeneID", "Gene ID"])

    out = pd.DataFrame()
    out["symbol"] = gc[symbol_col].astype(str).str.strip()
    out["rank"] = pd.to_numeric(gc[rank_col], errors="coerce") if rank_col else np.arange(1, len(gc) + 1)
    out["relevance_score"] = pd.to_numeric(gc[score_col], errors="coerce") if score_col else np.nan
    out["gene_type"] = gc[type_col] if type_col else ""
    out["knowledge_score"] = gc[knowledge_col] if knowledge_col else ""
    out["gene_id"] = gc[id_col].astype(str).str.strip() if id_col else ""
    out = out[out["symbol"].ne("")].copy()
    out = out.sort_values(["rank", "relevance_score"], ascending=[True, False], na_position="last")
    out = out.drop_duplicates("symbol", keep="first")

    ambiguous = 0
    mapped = []
    for row in out.to_dict("records"):
        ids = {row["gene_id"]} if row["gene_id"] and row["gene_id"] != "nan" else set()
        if not ids:
            ids = gene_symbol_to_ids.get(str(row["symbol"]).upper(), set())
        if len(ids) > 1:
            ambiguous += 1
        for gene_id in sorted(ids):
            rec = dict(row)
            rec["gene_id"] = str(gene_id)
            mapped.append(rec)
    mapped_df = pd.DataFrame(mapped)
    if mapped_df.empty:
        raise ValueError(f"GeneCards export {path} yielded no CTD-mapped genes")
    mapped_df["rank"] = pd.to_numeric(mapped_df["rank"], errors="coerce")
    mapped_df = mapped_df.sort_values(["rank", "relevance_score"], ascending=[True, False], na_position="last")
    mapped_df = mapped_df.drop_duplicates("gene_id", keep="first")
    return mapped_df, {"input_rows": len(gc), "mapped_rows": len(mapped_df), "ambiguous_symbols": ambiguous}


def bh_fdr(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    n = len(p)
    order = np.argsort(np.nan_to_num(p, nan=1.0))
    ranked = np.nan_to_num(p, nan=1.0)[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(n, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def odds_log2(value: float) -> float:
    if math.isinf(value):
        return math.inf
    if value <= 0:
        return -math.inf
    return math.log2(value)


def markdown_table(frame: pd.DataFrame) -> str:
    """Render a compact Markdown table without requiring the optional tabulate package."""
    if frame.empty:
        return ""
    headers = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for values in frame.itertuples(index=False, name=None):
        rendered = []
        for value in values:
            if pd.isna(value):
                text = ""
            elif isinstance(value, float):
                text = f"{value:.6g}"
            else:
                text = str(value)
            rendered.append(text.replace("|", "\\|"))
        lines.append("| " + " | ".join(rendered) + " |")
    return "\n".join(lines)


def enrichment_table(
    candidates: pd.DataFrame,
    gene_sets: dict[str, set[str]],
    genecards: pd.DataFrame,
    universe: set[str],
    scope: str,
    k: int,
    background: str,
) -> pd.DataFrame:
    gc = genecards.sort_values("rank", na_position="last").drop_duplicates("gene_id")
    gc = gc[gc["gene_id"].isin(universe)].copy()
    gene_cards_available_n = len(gc)
    gc = gc.head(k)
    crc_genes = set(gc["gene_id"])
    rank_map = dict(zip(gc["gene_id"], gc["rank"].astype(float)))
    n_universe = len(universe)
    n_crc = len(crc_genes)
    rows = []
    for row in candidates.to_dict("records"):
        chemical_id = row["ChemicalID"]
        interacting = gene_sets.get(chemical_id, set()) & universe
        n_interacting = len(interacting)
        overlap = interacting & crc_genes
        a = len(overlap)
        b = n_interacting - a
        c = n_crc - a
        d = n_universe - a - b - c
        if n_universe and n_crc and n_interacting:
            fisher_or, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
            enrichment_ratio = (a / n_interacting) / (n_crc / n_universe) if a else 0.0
            hyper_p = float(hypergeom.sf(a - 1, n_universe, n_crc, n_interacting))
        else:
            fisher_or, p_value, enrichment_ratio, hyper_p = 1.0, 1.0, 0.0, 1.0
        weights = [1.0 / math.log2(float(rank_map[g]) + 1.0) for g in overlap if g in rank_map]
        rows.append(
            {
                **row,
                "scope": scope,
                "gene_cards_k": k,
                "gene_cards_available_n": gene_cards_available_n,
                "gene_cards_actual_n": len(gc),
                "background": background,
                "n_universe": n_universe,
                "n_crc_genes_in_universe": n_crc,
                "n_ctd_human_genes": n_interacting,
                "crc_overlap": a,
                "non_crc_interacting": b,
                "crc_non_interacting": c,
                "non_crc_non_interacting": d,
                "odds_ratio": float(fisher_or),
                "log2_odds_ratio": odds_log2(float(fisher_or)),
                "enrichment_ratio": float(enrichment_ratio),
                "fisher_p": float(p_value),
                "hypergeom_p": hyper_p,
                "rank_weighted_overlap": float(sum(weights)),
                "mean_overlap_rank_weight": float(np.mean(weights)) if weights else 0.0,
                "top_overlap_genes": ";".join(
                    gc[gc["gene_id"].isin(overlap)].sort_values("rank")["symbol"].head(10).tolist()
                ),
                "stable_for_primary_sort": bool(n_interacting >= 20 and a >= 5),
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["bh_fdr"] = bh_fdr(result["fisher_p"].to_numpy())
    result["hypergeom_bh_fdr"] = bh_fdr(result["hypergeom_p"].to_numpy())
    result = result.sort_values(
        ["bh_fdr", "log2_odds_ratio", "rank_weighted_overlap", "crc_overlap"],
        ascending=[True, False, False, False],
        na_position="last",
    ).reset_index(drop=True)
    result["unfiltered_rank"] = np.arange(1, len(result) + 1)
    result["stable_rank"] = pd.NA
    stable_idx = result["stable_for_primary_sort"]
    result.loc[stable_idx, "stable_rank"] = np.arange(1, int(stable_idx.sum()) + 1)
    return result


def degree_bins(degrees: dict[str, int], n_bins: int = 20) -> dict[str, int]:
    genes = list(degrees)
    values = np.array([degrees[g] for g in genes], dtype=float)
    if len(set(values)) <= 1:
        return {g: 0 for g in genes}
    q = np.quantile(values, np.linspace(0, 1, min(n_bins, len(values)) + 1))
    q = np.unique(q)
    labels = np.digitize(values, q[1:-1], right=True)
    return dict(zip(genes, labels))


def matched_random_set(
    crc_genes: set[str],
    universe: set[str],
    degrees: dict[str, int],
    bins: dict[str, int],
    rng: np.random.Generator,
) -> set[str]:
    available = set(universe) - set(crc_genes)
    need = Counter(bins[g] for g in crc_genes)
    by_bin: dict[int, list[str]] = defaultdict(list)
    for gene in available:
        by_bin[bins.get(gene, 0)].append(gene)
    selected: list[str] = []
    for bin_id, count in sorted(need.items()):
        pool: list[str] = []
        radius = 0
        while len(pool) < count and radius <= max(bins.values(), default=0) + 1:
            for candidate_bin in {bin_id - radius, bin_id + radius}:
                pool.extend(g for g in by_bin.get(candidate_bin, []) if g not in selected and g not in pool)
            radius += 1
        if len(pool) < count:
            pool = [g for g in available if g not in selected]
        chosen = rng.choice(np.array(pool, dtype=object), size=count, replace=False).tolist()
        selected.extend(chosen)
    return set(selected)


def degree_matched_permutation(
    candidates: pd.DataFrame,
    gene_sets: dict[str, set[str]],
    crc_genes: set[str],
    universe: set[str],
    n_permutations: int,
    seed: int,
) -> pd.DataFrame:
    ids = candidates["ChemicalID"].tolist()
    universe_list = sorted(universe)
    col = {gene: i for i, gene in enumerate(universe_list)}
    row_idx: list[int] = []
    col_idx: list[int] = []
    for i, chemical_id in enumerate(ids):
        for gene in gene_sets.get(chemical_id, set()) & universe:
            row_idx.append(i)
            col_idx.append(col[gene])
    matrix = csr_matrix((np.ones(len(row_idx), dtype=np.int8), (row_idx, col_idx)), shape=(len(ids), len(universe_list)))
    degrees = {gene: 0 for gene in universe_list}
    for gene_set in gene_sets.values():
        for gene in gene_set & universe:
            degrees[gene] += 1
    bins = degree_bins(degrees)
    observed = np.array([len(gene_sets.get(cid, set()) & crc_genes) for cid in ids], dtype=int)
    exceed = np.zeros(len(ids), dtype=int)
    rng = np.random.default_rng(seed)
    for _ in range(n_permutations):
        random_genes = matched_random_set(crc_genes, universe, degrees, bins, rng)
        mask = np.zeros(len(universe_list), dtype=np.int8)
        mask[[col[g] for g in random_genes]] = 1
        null_overlap = np.asarray(matrix.dot(mask)).ravel()
        exceed += null_overlap >= observed
    empirical_p = (exceed + 1) / (n_permutations + 1)
    out = pd.DataFrame(
        {
            "ChemicalID": ids,
            "degree_matched_empirical_p": empirical_p,
            "degree_matched_bh_fdr": bh_fdr(empirical_p),
            "degree_matched_permutations": n_permutations,
            "degree_matched_seed": seed,
            "n_degree_bins": len(set(bins.values())),
        }
    )
    return out


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_report(
    path: Path,
    core_results: pd.DataFrame,
    top20: pd.DataFrame,
    classifications: pd.DataFrame,
    gene_card_meta: dict,
    manifest: dict,
    stability: pd.DataFrame,
) -> None:
    primary = core_results[(core_results["scope"] == "GeneCards_Disorders") & (core_results["gene_cards_k"] == 1000)]
    lines = [
        "# 环境毒理学与 CRC：Phase 1 CTD × GeneCards 结果",
        "",
        "## 运行状态",
        "",
        "本报告由冻结的 CTD chemical hierarchy、GeneCards Disorders-scoped CRC 查询、Fisher/BH-FDR 和 degree-matched permutation 流程生成。第一轮未进行独立 CRC 文献撞车审查。",
        "",
        f"- CTD 人类 interaction rows: {manifest['n_human_raw_interaction_rows']:,}",
        f"- CTD unique chemical–gene pairs: {manifest['n_human_unique_pairs']:,}",
        f"- Core environmental chemicals with human interactions: {manifest['n_core_chemicals']:,}",
        f"- Core chemical classes observed: {manifest['n_core_classes']}",
        f"- GeneCards mapped inputs: {', '.join(sorted(gene_card_meta))}",
        "",
        "## Primary top 20",
        "",
        "排序预先固定为 BH-FDR 升序、log2(OR) 降序、rank-weighted overlap 降序、overlap 降序；没有按化学物名称或研究热门程度人工挑选。",
        "",
    ]
    cols = ["unfiltered_rank", "ChemicalName", "chemical_class", "n_ctd_human_genes", "crc_overlap", "odds_ratio", "enrichment_ratio", "bh_fdr", "rank_weighted_overlap", "n_unique_pmids", "top_overlap_genes"]
    if not top20.empty:
        lines.append(markdown_table(top20[cols]))
    else:
        lines.append("未找到 GeneCards Disorders top 1000 主分析结果。")
    lines += [
        "",
        "## 主分析摘要",
        "",
        f"- Primary tested chemicals: {len(primary):,}",
        f"- Primary FDR < 0.05: {int((primary['bh_fdr'] < 0.05).sum()) if not primary.empty else 0}",
        f"- Primary stable candidates (n_interacting ≥ 20 and overlap ≥ 5): {int(primary['stable_for_primary_sort'].sum()) if not primary.empty else 0}",
        "",
        "## GeneCards scope stability",
        "",
    ]
    if not stability.empty:
        lines.append(markdown_table(stability))
    else:
        lines.append("没有足够的 GeneCards scope/K 配置生成稳定性比较。")
    lines += [
        "",
        "## 解释边界",
        "",
        "结果表示 CTD chemical-interacting gene set 与 GeneCards CRC-associated gene set 的富集关系，不证明真实人群暴露因果、CRC 发病风险或治疗效果。CTD PMID 字段仅用于审计与证据成熟度，第一轮未读文献。",
        "",
        "## Reproducibility",
        "",
        f"- Run timestamp (UTC): {manifest['run_timestamp_utc']}",
        f"- Random seed: {manifest['seed']}",
        f"- Degree-matched permutations: {manifest['n_permutations']}",
        "- Raw CTD archives are excluded from Git by repository .gitignore.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--interaction", type=Path, required=True)
    parser.add_argument("--chemicals", type=Path, required=True)
    parser.add_argument("--gene-cards-scoped", type=Path, required=True)
    parser.add_argument("--gene-cards-global", type=Path)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--drugcentral", type=Path)
    parser.add_argument("--pah-formulas", type=Path)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260822)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    interaction_counts, gene_sets, gene_symbols = load_interactions(args.interaction)
    chemicals = read_ctd_tsv(args.chemicals)
    classifications = classify_chemicals(chemicals, args.rules, args.drugcentral, args.pah_formulas)
    classifications = classifications.merge(interaction_counts, on="ChemicalID", how="left")
    classifications.to_csv(args.outdir / "environmental_toxicology_crc_phase1_chemical_classification.csv", index=False)

    core_ids = set(classifications.loc[classifications["is_core"] == True, "ChemicalID"])
    all_ids = set(gene_sets)
    core_candidates = classifications[classifications["ChemicalID"].isin(core_ids & all_ids)].copy()
    all_candidates = classifications[classifications["ChemicalID"].isin(all_ids)].copy()
    core_candidates = core_candidates.sort_values("ChemicalID").reset_index(drop=True)
    all_candidates = all_candidates.sort_values("ChemicalID").reset_index(drop=True)
    core_universe = set().union(*(gene_sets[c] for c in core_candidates["ChemicalID"])) if len(core_candidates) else set()
    all_universe = set().union(*(gene_sets[c] for c in all_candidates["ChemicalID"])) if len(all_candidates) else set()

    symbol_to_ids: dict[str, set[str]] = defaultdict(set)
    for gene_id, symbol in gene_symbols.items():
        if symbol:
            symbol_to_ids[str(symbol).upper()].add(str(gene_id))

    gene_card_paths = [("GeneCards_Disorders", args.gene_cards_scoped)]
    if args.gene_cards_global and args.gene_cards_global.exists():
        gene_card_paths.append(("GeneCards_Anywhere", args.gene_cards_global))

    all_core_results = []
    all_allctd_results = []
    gene_card_meta = {}
    mapped_gene_cards = {}
    for scope, path in gene_card_paths:
        gc, metadata = load_genecards(path, symbol_to_ids)
        mapped_gene_cards[scope] = gc
        gene_card_meta[scope] = {**metadata, "sha256": file_sha256(path), "path": str(path)}
        for k in (500, 1000, 2000):
            if len(gc) == 0:
                continue
            all_core_results.append(enrichment_table(core_candidates, gene_sets, gc, core_universe, scope, k, "U_core"))
            all_allctd_results.append(enrichment_table(core_candidates, gene_sets, gc, all_universe, scope, k, "U_allCTD"))

    core_results = pd.concat(all_core_results, ignore_index=True) if all_core_results else pd.DataFrame()
    allctd_results = pd.concat(all_allctd_results, ignore_index=True) if all_allctd_results else pd.DataFrame()
    core_results.to_csv(args.outdir / "environmental_toxicology_crc_phase1_ranked_core.csv", index=False)
    allctd_results.to_csv(args.outdir / "environmental_toxicology_crc_phase1_ranked_allctd.csv", index=False)

    primary = core_results[(core_results["scope"] == "GeneCards_Disorders") & (core_results["gene_cards_k"] == 1000)].copy()
    top20 = primary.head(20).copy()
    top20.to_csv(args.outdir / "environmental_toxicology_crc_phase1_top20.csv", index=False)

    permutation = pd.DataFrame()
    if not primary.empty:
        primary_gc = mapped_gene_cards["GeneCards_Disorders"]
        primary_gc = primary_gc[primary_gc["gene_id"].isin(core_universe)].sort_values("rank").head(1000)
        permutation = degree_matched_permutation(
            core_candidates,
            gene_sets,
            set(primary_gc["gene_id"]),
            core_universe,
            args.permutations,
            args.seed,
        )
    permutation.to_csv(args.outdir / "environmental_toxicology_crc_phase1_degree_matched_permutation.csv", index=False)

    stability_rows = []
    for (scope_a, k_a), (scope_b, k_b) in combinations(
        [(scope, k) for scope, _ in gene_card_paths for k in (500, 1000, 2000)], 2
    ):
        a = core_results[(core_results["scope"] == scope_a) & (core_results["gene_cards_k"] == k_a)]
        b = core_results[(core_results["scope"] == scope_b) & (core_results["gene_cards_k"] == k_b)]
        if a.empty or b.empty:
            continue
        sa, sb = set(a.head(20)["ChemicalID"]), set(b.head(20)["ChemicalID"])
        union = sa | sb
        stability_rows.append(
            {
                "configuration_a": f"{scope_a}_top{k_a}",
                "configuration_b": f"{scope_b}_top{k_b}",
                "top20_overlap": len(sa & sb),
                "top20_jaccard": len(sa & sb) / len(union) if union else 0.0,
            }
        )
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(args.outdir / "environmental_toxicology_crc_phase1_top20_stability.csv", index=False)

    manifest = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(ROOT)),
        "rules": str(args.rules.relative_to(ROOT)) if args.rules.is_relative_to(ROOT) else str(args.rules),
        "drugcentral": (
            {"path": str(args.drugcentral), "sha256": file_sha256(args.drugcentral)}
            if args.drugcentral and args.drugcentral.exists() else None
        ),
        "pah_formulas": (
            {"path": str(args.pah_formulas), "sha256": file_sha256(args.pah_formulas)}
            if args.pah_formulas and args.pah_formulas.exists() else None
        ),
        "interaction_file": str(args.interaction),
        "chemicals_file": str(args.chemicals),
        "gene_cards": gene_card_meta,
        "n_human_raw_interaction_rows": int(interaction_counts["n_raw_interaction_rows"].sum()),
        "n_human_unique_pairs": int(interaction_counts["n_unique_chemical_gene_pairs"].sum()),
        "n_core_chemicals": int(len(core_candidates)),
        "n_core_classes": int(classifications.loc[classifications["is_core"] == True, "chemical_class"].nunique()),
        "n_core_universe_genes": int(len(core_universe)),
        "n_allctd_universe_genes": int(len(all_universe)),
        "n_permutations": args.permutations,
        "seed": args.seed,
        "ctd_organism": "Homo sapiens",
        "ctd_pair_deduplication": "unique ChemicalID x GeneID",
        "primary_gene_cards_scope": "Disorders",
        "primary_gene_cards_k": 1000,
        "primary_background": "U_core",
        "independent_crc_literature_search": False,
    }
    (args.outdir / "environmental_toxicology_crc_phase1_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_report(args.outdir / "environmental_toxicology_crc_phase1_report.md", core_results, top20, classifications, gene_card_meta, manifest, stability)


if __name__ == "__main__":
    main()
