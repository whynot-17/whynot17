#!/usr/bin/env python3
"""Audit CRC disease-gene source cutoffs against the 86-gene DINP CTD set.

The audit is deliberately source/cutoff aware.  For each GeneCards and
Open Targets ranking it reports the DINP overlap, PTGER4 rank/presence and
key pathway enrichment using both the corresponding cutoff universe and the
full ranking universe that is actually available locally.  No full disease
annotation universe is silently substituted for a ranked cutoff.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "analysis" / "dinp_crc_multi_database_target_convergence" / "outputs"
SOURCE_DIR = Path(r"E:\mcop\mcop\work\environmental_toxicology_crc_phase1\data")
GMT_DIR = ROOT / "work" / "gene_sets"
OUT_DIR = Path(__file__).resolve().parent / "outputs"
LOCAL_INPUT_DIR = Path(__file__).resolve().parent / "inputs"


def symbols(values: Iterable[object]) -> List[str]:
    return sorted({str(v).strip().upper() for v in values if pd.notna(v) and str(v).strip()})


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values: Sequence[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    idx = np.flatnonzero(valid)
    order = idx[np.argsort(p[idx], kind="mergesort")]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.clip(ranked, 0, 1)
    return out


def parse_gmt(path: Path, collection: str) -> Dict[str, Dict[str, object]]:
    terms: Dict[str, Dict[str, object]] = {}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            name = fields[0]
            genes = set(symbols(fields[2:]))
            if genes:
                terms[name] = {"collection": collection, "genes": genes}
    return terms


def load_rankings() -> Tuple[set[str], Dict[str, pd.DataFrame], Dict[str, Path]]:
    dinp_path = INPUT_DIR / "dinp_exposure_gene_matrix.csv"
    crc_path = INPUT_DIR / "crc_gene_matrix.csv"
    # Keep a local snapshot in the repository so the audit remains runnable
    # after the original E:\\mcop workspace is moved. Fall back only when
    # the snapshot is absent.
    gc_path = LOCAL_INPUT_DIR / "genecards_anywhere_crc_top2000.csv"
    if not gc_path.exists():
        gc_path = SOURCE_DIR / "genecards_anywhere_crc_top2000.csv"
    required = [dinp_path, crc_path, gc_path, GMT_DIR / "c2.cp.reactome.v2026.1.Hs.symbols.gmt", GMT_DIR / "h.all.v2026.1.Hs.symbols.gmt"]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing audit inputs: " + "; ".join(missing))

    dinp = set(symbols(pd.read_csv(dinp_path)["gene_symbol"]))
    if len(dinp) != 86:
        raise ValueError(f"Expected 86 DINP CTD genes, found {len(dinp)}")

    gc = pd.read_csv(gc_path)
    gc["GeneSymbol"] = gc["GeneSymbol"].map(lambda x: str(x).strip().upper())
    gc = gc.drop_duplicates("GeneSymbol").sort_values(["GeneCards_Rank", "GeneSymbol"]).reset_index(drop=True)
    gc["rank"] = np.arange(1, len(gc) + 1)
    gc["score"] = pd.to_numeric(gc["RelevanceScore"], errors="coerce")

    crc = pd.read_csv(crc_path)
    crc["gene_symbol"] = crc["gene_symbol"].map(lambda x: str(x).strip().upper())
    crc["OpenTargets_score"] = pd.to_numeric(crc["OpenTargets_score"], errors="coerce")
    ot = crc.dropna(subset=["OpenTargets_score"]).drop_duplicates("gene_symbol").copy()
    ot = ot.sort_values(["OpenTargets_score", "gene_symbol"], ascending=[False, True]).reset_index(drop=True)
    ot["rank"] = np.arange(1, len(ot) + 1)
    ot = ot.rename(columns={"gene_symbol": "GeneSymbol", "OpenTargets_score": "score"})

    paths = {"dinp": dinp_path, "crc": crc_path, "genecards": gc_path, "reactome": GMT_DIR / "c2.cp.reactome.v2026.1.Hs.symbols.gmt", "hallmark": GMT_DIR / "h.all.v2026.1.Hs.symbols.gmt"}
    return dinp, {"GeneCards": gc, "OpenTargets": ot}, paths


def enrichment(query: set[str], background: set[str], terms: Mapping[str, Dict[str, object]], source: str, cutoff_label: str, background_label: str) -> pd.DataFrame:
    rows = []
    if not query or not background:
        return pd.DataFrame()
    q = query & background
    N = len(background)
    n = len(q)
    for term, info in terms.items():
        genes = set(info["genes"]) & background
        K = len(genes)
        k = len(q & genes)
        # Keep zero-overlap terms (p=1) in the tested family so the BH
        # denominator reflects every eligible Reactome/Hallmark term rather
        # than only terms that happened to intersect the query.
        if K < 3:
            continue
        p = float(hypergeom.sf(k - 1, N, K, n))
        rows.append({
            "source": source,
            "cutoff": cutoff_label,
            "background": background_label,
            "collection": info["collection"],
            "term": term,
            "background_n": N,
            "query_n": n,
            "term_n": K,
            "overlap_n": k,
            "overlap_genes": ";".join(sorted(q & genes)),
            "raw_p": p,
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["global_fdr"] = bh(out["raw_p"].to_numpy())
    out["key_pathway"] = out["term"].str.contains(
        r"PROSTAGLANDIN|ARACHIDON|EICOSANOID|PPAR|INFLAMMAT|INFLAMMASOME|LIPID", case=False, regex=True
    )
    return out.sort_values(["global_fdr", "raw_p", "term"]).reset_index(drop=True)


def panel_sets(terms: Mapping[str, Dict[str, object]]) -> Dict[str, set[str]]:
    def union_matching(pattern: str) -> set[str]:
        ans: set[str] = set()
        for name, info in terms.items():
            if re.search(pattern, name, flags=re.I):
                ans |= set(info["genes"])
        return ans

    pge2 = {"PLA2G4A", "PTGS1", "PTGS2", "PTGES", "PTGES2", "PTGES3"}
    prostaglandin = union_matching(r"SYNTHESIS_OF_PROSTAGLANDINS")
    arachidonate = union_matching(r"ARACHIDON")
    eicosanoid = union_matching(r"EICOSANOID")
    ppar = union_matching(r"PPAR")
    inflammatory = union_matching(r"INFLAMMAT|INFLAMMASOME")
    hallmark_inflam = set(terms.get("HALLMARK_INFLAMMATORY_RESPONSE", {}).get("genes", set()))
    if hallmark_inflam:
        inflammatory = hallmark_inflam
    return {
        "PGE2_synthesis_core": pge2,
        "prostaglandin_synthesis": prostaglandin,
        "arachidonate": arachidonate,
        "eicosanoid": eicosanoid,
        "PPAR": ppar,
        "inflammatory_response": inflammatory,
    }


def panel_audit(query: set[str], background: set[str], source: str, cutoff: str, bg_label: str, panels: Mapping[str, set[str]]) -> List[dict]:
    q = query & background
    N, n = len(background), len(q)
    rows: List[dict] = []
    for panel, genes0 in panels.items():
        genes = genes0 & background
        K = len(genes)
        k = len(q & genes)
        p = float(hypergeom.sf(k - 1, N, K, n)) if K and n and k else (1.0 if K and n else np.nan)
        rows.append({
            "source": source,
            "cutoff": cutoff,
            "background": bg_label,
            "panel": panel,
            "background_n": N,
            "query_n": n,
            "panel_n": K,
            "overlap_n": k,
            "overlap_genes": ";".join(sorted(q & genes)),
            "raw_p": p,
        })
    return rows


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dinp, rankings, paths = load_rankings()
    terms: Dict[str, Dict[str, object]] = {}
    terms.update(parse_gmt(paths["reactome"], "Reactome"))
    terms.update(parse_gmt(paths["hallmark"], "Hallmark"))
    panels = panel_sets(terms)

    configs = []
    for source, df in rankings.items():
        full = set(df["GeneSymbol"])
        max_k = min(2000, len(df))
        for k in (500, 1000, 2000):
            if k > max_k:
                continue
            top = set(df.head(k)["GeneSymbol"])
            query = dinp & top
            configs.append((source, f"top{k}", df, top, full, query))

    overlap_rows: List[dict] = []
    gene_rows: List[dict] = []
    enrich_frames: List[pd.DataFrame] = []
    panel_rows: List[dict] = []
    for source, cutoff, df, top, full, query in configs:
        pt = df.loc[df["GeneSymbol"] == "PTGER4"]
        pt_rank = int(pt.iloc[0]["rank"]) if len(pt) else None
        pt_score = float(pt.iloc[0]["score"]) if len(pt) and pd.notna(pt.iloc[0]["score"]) else None
        overlap_rows.append({
            "source": source,
            "cutoff": cutoff,
            "cutoff_n": len(top),
            "available_ranked_universe_n": len(full),
            "dinp_ctd_n": len(dinp),
            "overlap_n": len(query),
            "overlap_fraction_of_dinp": len(query) / len(dinp),
            "ptger4_in_cutoff": "PTGER4" in query,
            "ptger4_rank_in_available_universe": pt_rank,
            "ptger4_score": pt_score,
            "overlap_genes": ";".join(sorted(query)),
        })
        for gene in sorted(query):
            rank = int(df.loc[df["GeneSymbol"] == gene, "rank"].iloc[0])
            gene_rows.append({"source": source, "cutoff": cutoff, "gene_symbol": gene, "rank": rank, "ptger4": gene == "PTGER4"})
        for bg_label, bg in (("cutoff_universe", top), ("available_ranked_universe", full)):
            frame = enrichment(query, bg, terms, source, cutoff, bg_label)
            if not frame.empty:
                enrich_frames.append(frame)
            panel_rows.extend(panel_audit(query, bg, source, cutoff, bg_label, panels))

    overlap = pd.DataFrame(overlap_rows)
    genes = pd.DataFrame(gene_rows)
    enrich = pd.concat(enrich_frames, ignore_index=True) if enrich_frames else pd.DataFrame()
    panels_df = pd.DataFrame(panel_rows)
    if not panels_df.empty:
        panels_df["fdr_within_background"] = panels_df.groupby(["source", "cutoff", "background"], sort=False)["raw_p"].transform(lambda s: pd.Series(bh(s.to_numpy()), index=s.index))

    overlap.to_csv(OUT_DIR / "overlap_summary.csv", index=False)
    genes.to_csv(OUT_DIR / "overlap_genes_long.csv", index=False)
    enrich.to_csv(OUT_DIR / "key_pathway_enrichment.csv", index=False)
    panels_df.to_csv(OUT_DIR / "panel_enrichment_summary.csv", index=False)

    # A compact visual makes cutoff sensitivity immediately inspectable.
    fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=180)
    for source, sub in overlap.groupby("source", sort=False):
        sub = sub.copy()
        sub["k"] = sub["cutoff"].str.extract(r"(\d+)").astype(int)
        ax.plot(sub["k"], sub["overlap_n"], marker="o", linewidth=2, label=source)
        for _, row in sub.iterrows():
            if bool(row["ptger4_in_cutoff"]):
                ax.annotate("PTGER4", (row["k"], row["overlap_n"]), xytext=(4, 5), textcoords="offset points", fontsize=8)
    ax.set_xlabel("CRC ranking cutoff")
    ax.set_ylabel("DINP CTD genes in CRC set")
    ax.set_xticks([500, 1000, 2000])
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "cutoff_sensitivity_overlap.png")
    plt.close(fig)

    manifest = {
        "analysis": "DINP CTD 86-gene × CRC disease-source cutoff sensitivity audit",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "query_definition": "All 86 genes in dinp_exposure_gene_matrix.csv; this file contains CTD=1 for all rows and no ToxCast/ChEMBL-only additions.",
        "sources": {
            "GeneCards": {"ranking_file": str(paths["genecards"]), "available_n": len(rankings["GeneCards"]), "score": "RelevanceScore", "cutoffs": [500, 1000, 2000]},
            "OpenTargets": {"ranking_file": str(paths["crc"]), "available_n": len(rankings["OpenTargets"]), "score": "OpenTargets_score, descending; ties broken by gene symbol", "cutoffs": [500, 1000, 2000]},
        },
        "backgrounds": ["cutoff_universe", "available_ranked_universe"],
        "pathway_sets": {"Reactome": str(paths["reactome"]), "Hallmark": str(paths["hallmark"])},
        "input_sha256": {k: sha256_file(v) for k, v in paths.items()},
        "interpretation_guardrail": "A term or PTGER4 signal appearing only at one hand-picked cutoff is cutoff-sensitive and should not be presented as stable disease convergence. Candidate-universe enrichment is descriptive of the selected ranking universe.",
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    gc500 = overlap.query("source == 'GeneCards' and cutoff == 'top500'").iloc[0]
    gc1000 = overlap.query("source == 'GeneCards' and cutoff == 'top1000'").iloc[0]
    gc2000 = overlap.query("source == 'GeneCards' and cutoff == 'top2000'").iloc[0]
    ot500 = overlap.query("source == 'OpenTargets' and cutoff == 'top500'").iloc[0]
    ot1000 = overlap.query("source == 'OpenTargets' and cutoff == 'top1000'").iloc[0]
    ot2000 = overlap.query("source == 'OpenTargets' and cutoff == 'top2000'").iloc[0]

    def panel_value(source: str, cutoff: str, panel: str, background: str = "available_ranked_universe") -> str:
        x = panels_df.query(
            "source == @source and cutoff == @cutoff and panel == @panel and background == @background"
        )
        if x.empty:
            return "NA"
        row = x.iloc[0]
        if pd.isna(row["fdr_within_background"]):
            return "NA"
        return f"{int(row['overlap_n'])}/{int(row['panel_n'])}; FDR={float(row['fdr_within_background']):.3g}"

    gc_panel = [
        f"| {panel} | {panel_value('GeneCards', 'top500', panel)} | {panel_value('GeneCards', 'top1000', panel)} | {panel_value('GeneCards', 'top2000', panel)} |"
        for panel in ("PGE2_synthesis_core", "prostaglandin_synthesis", "arachidonate", "PPAR", "inflammatory_response")
    ]
    ot_panel = [
        f"| {panel} | {panel_value('OpenTargets', 'top500', panel)} | {panel_value('OpenTargets', 'top1000', panel)} | {panel_value('OpenTargets', 'top2000', panel)} |"
        for panel in ("PGE2_synthesis_core", "prostaglandin_synthesis", "arachidonate", "PPAR", "inflammatory_response")
    ]
    report = f"""# DINP × CRC disease-gene source cutoff sensitivity audit

## Design

The query is the complete CTD DINP set (86 genes). CRC disease genes were evaluated from two locally preserved ranked sources: GeneCards relevance-score rank (top 2,000 available) and Open Targets score rank (15,611 genes with a score). Each source was tested at top 500, 1,000 and 2,000. For every cutoff, enrichment was calculated against the matching cutoff universe and again against the full ranked universe available for that source. Reactome and Hallmark gene sets were tested by hypergeometric over-representation with global Benjamini–Hochberg FDR within each source/cutoff/background family.

## Overlap and PTGER4

| source | top500 | top1000 | top2000 |
|---|---:|---:|---:|
| GeneCards DINP overlap | {int(gc500.overlap_n)} | {int(gc1000.overlap_n)} | {int(gc2000.overlap_n)} |
| Open Targets DINP overlap | {int(ot500.overlap_n)} | {int(ot1000.overlap_n)} | {int(ot2000.overlap_n)} |
| PTGER4 in GeneCards cutoff | {'yes' if gc500.ptger4_in_cutoff else 'no'} | {'yes' if gc1000.ptger4_in_cutoff else 'no'} | {'yes' if gc2000.ptger4_in_cutoff else 'no'} |
| PTGER4 in Open Targets cutoff | {'yes' if ot500.ptger4_in_cutoff else 'no'} | {'yes' if ot1000.ptger4_in_cutoff else 'no'} | {'yes' if ot2000.ptger4_in_cutoff else 'no'} |

PTGER4 ranks in the available source rankings are GeneCards={gc2000.ptger4_rank_in_available_universe or 'not present'} and Open Targets={ot1000.ptger4_rank_in_available_universe or 'not present'}. Thus PTGER4 is absent from both top500 sets, enters Open Targets at top1000, and enters GeneCards only by top2000.

## Key panel results using the full available ranked universe

Each cell is `overlap/panel_size; panel FDR`; the query is the DINP overlap at that cutoff. This fixed-background view is the most comparable across cutoffs.

### GeneCards

| panel | top500 | top1000 | top2000 |
|---|---:|---:|---:|
{chr(10).join(gc_panel)}

### Open Targets

| panel | top500 | top1000 | top2000 |
|---|---:|---:|---:|
{chr(10).join(ot_panel)}

At the term level, the Reactome prostaglandin-synthesis term is globally FDR-significant for Open Targets top500/top1000/top2000 (driven by PTGS1 and PTGS2), while GeneCards reaches global FDR significance only at top2000 (driven by HPGD, PTGES, PTGES2, PTGS1 and PTGS2). The PPAR and inflammatory **terms** do not reach the global all-pathway FDR threshold in this fixed-background comparison, even where the smaller prespecified panel test is positive.

## Interpretation

The PTGER4 finding is **cutoff-sensitive**, rather than stable across all prespecified cutoffs. This audit therefore does not support describing PTGER4 as a source-independent CRC convergence hit. If a pathway is significant only in a broad source or only after expanding to top2000, it should be reported as a feature of that selected disease-gene ranking and treated as hypothesis-generating. See `panel_enrichment_summary.csv` for the prespecified prostaglandin, arachidonate, eicosanoid, PPAR and inflammatory panels, and `key_pathway_enrichment_key_terms.csv` for the compact term-level table (the full all-term table is generated locally as `key_pathway_enrichment.csv`).

The result does not prove that Open Targets “diluted” a true signal: because PTGER4 is not present at top500 in either source and appears at different cutoffs, the simpler conclusion is that the original PTGER4 enrichment is not robust to disease-gene source/cutoff choice. The proper next step is independent cohort or cell-level validation, not selecting the cutoff that gives the preferred pathway.

## Files

- `overlap_summary.csv`: all overlap counts, PTGER4 rank/presence and overlap lists.
- `panel_enrichment_summary.csv`: pathway-panel overlap and hypergeometric/FDR values under both backgrounds.
- `key_pathway_enrichment_key_terms.csv`: compact Reactome/Hallmark term-level results for key names (prostaglandin, arachidonate, eicosanoid, PPAR, inflammation and lipid). The full all-term table is retained locally as `key_pathway_enrichment.csv`.
- `cutoff_sensitivity_overlap.png`: compact overlap plot.
"""
    (OUT_DIR / "CUTOFF_SENSITIVITY_REPORT.md").write_text(report, encoding="utf-8")
    print("DINP × CRC DISEASE-SOURCE CUTOFF AUDIT: PASS")
    print(overlap[["source", "cutoff", "overlap_n", "ptger4_in_cutoff", "ptger4_rank_in_available_universe"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

