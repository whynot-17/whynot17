from __future__ import annotations

import io
import json
import time
from datetime import date, datetime, timezone
from pathlib import Path

import networkx as nx
import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "2026-10-10" / "go_kegg_ppi_113_vs_7748"
QUERY_FILE = ROOT / "results" / "2026-10-10" / "geo_replication_113_replacement" / "query_113_genes.txt"
BACKGROUND_FILE = OUT / "background_7748_genes.txt"
SOURCE_DEG_PATH = (
    r"D:\Codex\projectless-workspaces\2026-09-08\dinp-crc-replication-b24a399"
    r"\outputs\dinp_crc_dehp_replication\TCGA_GDC_current_DEG_FDR_threshold.tsv"
)
SOURCE_DEG_SHA256 = "7975288bc71e9ed80f71f75f8c64f5b32c3cc2e727297dbc820dc808f9b47ebc"
GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
STRING_BASE = "https://version-12-5.string-db.org/api"
SPECIES = 9606
CALLER = "DINP_CRC_113_vs_7748"
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG"]


def read_symbols(path: Path) -> list[str]:
    values = [value.strip().upper() for value in path.read_text(encoding="utf-8").splitlines()]
    return sorted(set(value for value in values if value))


def post_json(url: str, payload: dict, attempts: int = 3) -> dict:
    response = None
    for attempt in range(1, attempts + 1):
        try:
            response = requests.post(
                url,
                json=payload,
                headers={"User-Agent": "DINP-CRC-enrichment/1.0"},
                timeout=240,
            )
            break
        except requests.exceptions.RequestException:
            if attempt == attempts:
                raise
            time.sleep(attempt * 3)
    assert response is not None
    response.raise_for_status()
    return response.json()


def post_tsv(method: str, params: dict) -> pd.DataFrame:
    response = requests.post(
        f"{STRING_BASE}/tsv/{method}",
        data=params,
        headers={"User-Agent": "DINP-CRC-enrichment/1.0"},
        timeout=240,
    )
    response.raise_for_status()
    if not response.text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(response.text), sep="\t")


def run_gprofiler(query: list[str], background: list[str], correction: str) -> tuple[dict, pd.DataFrame]:
    payload = {
        "organism": "hsapiens",
        "query": query,
        "sources": SOURCES,
        "user_threshold": 1.0,
        "significance_threshold_method": correction,
        "domain_scope": "custom",
        "background": background,
        "no_iea": False,
        "no_evidences": False,
    }
    data = post_json(GPROFILER_URL, payload)
    (OUT / f"gprofiler_request_{correction}.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    (OUT / f"gprofiler_response_{correction}.json").write_text(
        json.dumps(data, ensure_ascii=False, separators=(",", ":")), encoding="utf-8"
    )
    meta = data.get("meta", {})
    query_order = meta.get("query_metadata", {}).get("queries", {}).get("query_1", query)
    rows = []
    for term in data.get("result", []):
        if term.get("source") == "KEGG" and term.get("native") == "KEGG:00000":
            continue
        evidence = term.get("intersections") or []
        hit_genes = [gene for gene, codes in zip(query_order, evidence) if codes]
        p_value = term.get("p_value")
        rows.append({
            "source": term.get("source"),
            "term_id": term.get("native"),
            "term_name": term.get("name"),
            f"p_value_{correction}": p_value,
            f"significant_{correction}_lt_0.05": bool(p_value is not None and p_value < 0.05),
            "query_size_mapped": term.get("query_size"),
            "background_size_effective": term.get("effective_domain_size"),
            "overlap_size": term.get("intersection_size"),
            "term_size": term.get("term_size"),
            "overlap_genes": ";".join(hit_genes),
            "description": term.get("description", ""),
        })
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame = frame.sort_values(["source", f"p_value_{correction}", "term_id"], kind="stable")
    return data, frame


def reduce_go_redundancy(significant: pd.DataFrame) -> pd.DataFrame:
    if significant.empty:
        return pd.DataFrame(columns=list(significant.columns) + ["redundancy_cluster", "cluster_term_count", "cluster_terms"])
    representatives = []
    go = significant[significant["source"].str.startswith("GO:")].copy()
    for source, block in go.groupby("source", sort=True):
        block = block.reset_index(drop=True)
        graph = nx.Graph()
        graph.add_nodes_from(range(len(block)))
        gene_sets = [set(str(value).split(";")) - {""} for value in block["overlap_genes"]]
        for i in range(len(block)):
            for j in range(i + 1, len(block)):
                union = gene_sets[i] | gene_sets[j]
                jaccard = len(gene_sets[i] & gene_sets[j]) / len(union) if union else 0.0
                if jaccard >= 0.5:
                    graph.add_edge(i, j, jaccard=jaccard)
        for cluster_n, cluster in enumerate(nx.connected_components(graph), start=1):
            members = block.iloc[sorted(cluster)].copy()
            p_col = "p_value_g_SCS"
            representative = members.sort_values([p_col, "overlap_size"], ascending=[True, False]).iloc[0].to_dict()
            representative["redundancy_cluster"] = f"{source}_R{cluster_n}"
            representative["cluster_term_count"] = len(members)
            representative["cluster_terms"] = ";".join(members["term_name"].astype(str).tolist())
            representatives.append(representative)
    if not representatives:
        return pd.DataFrame(columns=list(significant.columns) + ["redundancy_cluster", "cluster_term_count", "cluster_terms"])
    return pd.DataFrame(representatives).sort_values(["source", "p_value_g_SCS"])


def map_string_background(background: list[str]) -> pd.DataFrame:
    chunks = [background[i:i + 500] for i in range(0, len(background), 500)]
    mapped_parts = []
    for index, chunk in enumerate(chunks):
        if index:
            time.sleep(1)
        mapped_parts.append(post_tsv("get_string_ids", {
            "identifiers": "\r".join(chunk),
            "species": SPECIES,
            "echo_query": 1,
            "caller_identity": CALLER,
        }))
    mapped = pd.concat(mapped_parts, ignore_index=True) if mapped_parts else pd.DataFrame()
    required = ["queryItem", "queryIndex", "stringId", "ncbiTaxonId", "taxonName", "preferredName"]
    missing = [column for column in required if column not in mapped.columns]
    if missing:
        raise RuntimeError(f"STRING mapping response missing columns: {missing}")
    mapped = mapped[required].drop_duplicates(["queryItem", "stringId"])
    return mapped


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    query = read_symbols(QUERY_FILE)
    background = read_symbols(BACKGROUND_FILE)
    if len(query) != 113 or len(background) != 7748:
        raise ValueError(f"Expected 113 query and 7,748 background genes; got {len(query)} and {len(background)}")
    outside = sorted(set(query) - set(background))
    if outside:
        raise ValueError(f"Query genes absent from DEG background: {outside}")
    (OUT / "query_113_genes.txt").write_text("\n".join(query) + "\n", encoding="utf-8")

    gscs_data, gscs_terms = run_gprofiler(query, background, "g_SCS")
    time.sleep(1)
    fdr_data, fdr_terms = run_gprofiler(query, background, "fdr")
    gscs_terms.to_csv(OUT / "GO_KEGG_all_results_gSCS.tsv", sep="\t", index=False)
    fdr_terms.to_csv(OUT / "GO_KEGG_all_results_FDR.tsv", sep="\t", index=False)
    gscs_sig = gscs_terms[gscs_terms["significant_g_SCS_lt_0.05"]].copy()
    fdr_sig = fdr_terms[fdr_terms["significant_fdr_lt_0.05"]].copy()
    gscs_sig.to_csv(OUT / "GO_KEGG_significant_gSCS.tsv", sep="\t", index=False)
    fdr_sig.to_csv(OUT / "GO_KEGG_significant_FDR.tsv", sep="\t", index=False)
    gscs_sig[gscs_sig["source"] == "KEGG"].to_csv(OUT / "KEGG_significant_gSCS.tsv", sep="\t", index=False)
    representatives = reduce_go_redundancy(gscs_sig)
    representatives.to_csv(OUT / "GO_significant_representative_terms_gSCS.tsv", sep="\t", index=False)
    nonsig = fdr_terms[~fdr_terms["significant_fdr_lt_0.05"]].groupby("source", group_keys=False).head(5) if not fdr_terms.empty else fdr_terms.copy()
    nonsig.to_csv(OUT / "GO_KEGG_top_ranked_FDR_non_significant.tsv", sep="\t", index=False)

    time.sleep(1)
    string_mapping = map_string_background(background)
    mapping_file = OUT / "STRING_mapping_background_7748.tsv"
    string_mapping.to_csv(mapping_file, sep="\t", index=False)
    string_mapping[ string_mapping["queryItem"].astype(str).str.upper().isin(query) ].to_csv(
        OUT / "STRING_mapping_query_113.tsv", sep="\t", index=False
    )
    query_map = string_mapping[string_mapping["queryItem"].astype(str).str.upper().isin(query)].drop_duplicates("queryItem")
    missing_query = sorted(set(query) - set(query_map["queryItem"].astype(str).str.upper()))
    if missing_query:
        raise ValueError(f"Unmapped query genes in STRING background mapping: {missing_query}")
    query_ids = query_map["stringId"].astype(str).tolist()
    background_ids = string_mapping["stringId"].drop_duplicates().astype(str).tolist()
    ppi_rows = {}
    for score in (700, 400):
        params = {
            "identifiers": "\r".join(query_ids),
            "species": SPECIES,
            "required_score": score,
            "background_string_identifiers": "\r".join(background_ids),
            "caller_identity": CALLER,
        }
        time.sleep(1)
        result = post_tsv("ppi_enrichment", params)
        label = f"score_ge_{score / 1000:.3f}"
        result.to_csv(OUT / f"STRING_PPI_enrichment_custom_background_{label}.tsv", sep="\t", index=False)
        ppi_rows[label] = result.to_dict(orient="records")

    try:
        version_response = requests.get(f"{STRING_BASE}/tsv/version", timeout=60)
        version_response.raise_for_status()
        string_version = pd.read_csv(io.StringIO(version_response.text), sep="\t").to_dict(orient="records")
    except requests.RequestException:
        string_version = []
    gscs_meta = gscs_data.get("meta", {})
    fdr_meta = fdr_data.get("meta", {})
    metadata = {
        "analysis_date": date.today().isoformat(),
        "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "organism": "Homo sapiens (NCBI taxon 9606)",
        "query_count": len(query),
        "background_count": len(background),
        "query_all_within_background": True,
        "background_source_path": SOURCE_DEG_PATH,
        "background_source_sha256": SOURCE_DEG_SHA256,
        "background_source_rows": len(background),
        "background_source_unique_symbols": len(background),
        "background_source_fdr_hit_true": len(background),
        "background_symbols_file": "background_7748_genes.txt",
        "gprofiler_api": GPROFILER_URL,
        "gprofiler_data_version": gscs_meta.get("version"),
        "gprofiler_timestamp": gscs_meta.get("timestamp"),
        "gprofiler_g_SCS_query_metadata": gscs_meta.get("query_metadata", {}),
        "gprofiler_FDR_query_metadata": fdr_meta.get("query_metadata", {}),
        "gprofiler_g_SCS_result_metadata": gscs_meta.get("result_metadata", {}),
        "gprofiler_significant_terms_by_source_g_SCS": gscs_sig.groupby("source").size().to_dict() if not gscs_sig.empty else {},
        "gprofiler_significant_terms_by_source_FDR": fdr_sig.groupby("source").size().to_dict() if not fdr_sig.empty else {},
        "string_api_base": STRING_BASE,
        "string_version": string_version,
        "string_query_mapped_count": len(query_map),
        "string_background_mapped_unique_ids": len(background_ids),
        "string_background_unmapped_symbols": sorted(set(background) - set(string_mapping["queryItem"].astype(str).str.upper())),
        "string_ppi_enrichment_custom_background": ppi_rows,
        "string_enrichment_background_definition": "All 7,748 CRC FDR-DEG symbols, mapped to unique STRING IDs; query is the 113-gene intersection.",
        "string_source_urls": {
            "API_documentation": "https://string-db.org/help/api/",
            "version": f"{STRING_BASE}/tsv/version",
        },
    }
    (OUT / "STRING_request_parameters.json").write_text(json.dumps({
        "species": SPECIES,
        "required_score_thresholds": [700, 400],
        "query_string_ids": query_ids,
        "background_string_id_count": len(background_ids),
        "background_string_ids": background_ids,
        "caller_identity": CALLER,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "analysis_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    gscs_counts = gscs_sig.groupby("source").size().to_dict() if not gscs_sig.empty else {}
    fdr_counts = fdr_sig.groupby("source").size().to_dict() if not fdr_sig.empty else {}
    fdr_min = fdr_terms.groupby("source")["p_value_fdr"].min().to_dict() if not fdr_terms.empty else {}
    lines = [
        "# 113 个 DINP×CRC DEG 交集基因的 GO/KEGG 与 PPI 富集",
        "",
        f"- 分析日期：{metadata['analysis_date']}",
        "- 查询集：113 个 DINP 与 TCGA/GDC CRC DEG 的交集基因。",
        "- 背景集：7,748 个 TCGA/GDC CRC FDR-DEG；输入文件全部 7,748 行 `fdr_hit=True`，7,748 个唯一基因符号。查询集 113/113 均在背景中。",
        f"- 背景输入来源：`{SOURCE_DEG_PATH}`；SHA-256：`{SOURCE_DEG_SHA256}`。本项目保存了背景基因符号副本 `background_7748_genes.txt`。",
        "- 物种：Homo sapiens（NCBI taxon 9606）。",
        "- GO/KEGG：g:Profiler g:GOSt；GO:BP、GO:MF、GO:CC、KEGG；使用 custom background；主分析为 g:SCS 校正后 p<0.05，BH-FDR 作为敏感性分析；ORA 不区分肿瘤上调/下调方向。",
        "- STRING PPI enrichment：查询集 113 个基因，背景使用上述 7,748 个 DEG 中成功映射后的唯一 STRING 蛋白 ID；分数阈值 0.700 和 0.400。",
        "- PPI 网络边表及 GraphML 网络见 `../ppi_113_genes/`，只包含 113 个查询基因，`add_nodes=0`。Functional 网络为多证据功能关联，不等同于直接物理结合。",
        "",
        "## 富集条目数量",
        "",
        "| 数据源 | g:SCS p<0.05 | BH-FDR p<0.05 | 最小 BH-FDR p |",
        "|---|---:|---:|---:|",
    ]
    for source in SOURCES:
        min_value = fdr_min.get(source)
        min_text = f"{min_value:.3g}" if min_value is not None else "NA"
        lines.append(f"| {source} | {gscs_counts.get(source, 0)} | {fdr_counts.get(source, 0)} | {min_text} |")
    lines.extend(["", "## 显著 GO 代表条目", ""])
    if representatives.empty:
        lines.append("无 GO 条目达到 g:SCS p<0.05。")
    else:
        lines.extend(["代表项按重叠基因 Jaccard≥0.5 合并相似 GO 条目；完整结果见 TSV。", "", "| 类别 | 条目 | g:SCS p | 重叠基因 | 簇内条目数 |", "|---|---|---:|---|---:|"])
        for _, row in representatives.head(30).iterrows():
            lines.append(f"| {row['source']} | {row['term_name']} ({row['term_id']}) | {row['p_value_g_SCS']:.3g} | {row['overlap_genes']} | {row['cluster_term_count']} |")
    lines.extend(["", "## 显著 KEGG 通路", ""])
    kegg_sig = gscs_sig[gscs_sig["source"] == "KEGG"]
    if kegg_sig.empty:
        lines.append("无 KEGG 通路达到 g:SCS p<0.05。")
    else:
        lines.extend(["| 通路 | g:SCS p | 重叠基因 |", "|---|---:|---|"])
        for _, row in kegg_sig.iterrows():
            lines.append(f"| {row['term_name']} ({row['term_id']}) | {row['p_value_g_SCS']:.3g} | {row['overlap_genes']} |")
    lines.extend(["", "## 文件清单", "", "- `background_7748_genes.txt`、`query_113_genes.txt`：实际使用的输入名单。", "- `GO_KEGG_all_results_gSCS.tsv`、`GO_KEGG_significant_gSCS.tsv`：g:SCS 全部结果与显著项。", "- `GO_KEGG_all_results_FDR.tsv`、`GO_KEGG_significant_FDR.tsv`：BH-FDR 敏感性分析。", "- `GO_significant_representative_terms_gSCS.tsv`、`KEGG_significant_gSCS.tsv`：代表 GO 项与显著 KEGG 项。", "- `GO_KEGG_top_ranked_FDR_non_significant.tsv`：各来源最靠前但未达 BH-FDR 显著的条目。", "- `gprofiler_request_*.json` / `gprofiler_response_*.json`：API 请求及原始响应。", "- `STRING_mapping_*.tsv`、`STRING_PPI_enrichment_custom_background_*.tsv`：STRING 映射及背景校正 PPI enrichment。", "- `STRING_request_parameters.json`、`analysis_metadata.json`：分析参数、数据版本、背景来源与统计摘要。", ""])
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({
        "query": len(query),
        "background": len(background),
        "g_SCS_significant": gscs_counts,
        "FDR_significant": fdr_counts,
        "STRING_background_mapped": len(background_ids),
        "STRING_PPI_enrichment_rows": {key: len(value) for key, value in ppi_rows.items()},
        "output": str(OUT),
    }, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
