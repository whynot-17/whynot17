"""Phase 7C: trajectory-conditioned functional-module convergence.

This analysis deliberately does not search for a single convergent gene.  For
each acquired-OXA-resistance trajectory, the complete DepMap dependency
ranking is tested against Reactome/Hallmark pathways, curated mechanistic
modules, CORUM protein complexes, and public co-essentiality modules.

The ranking is vulnerability-oriented: a positive gene score means that a
cell line with a higher OXA-R-like state score has a stronger dependency on
that gene (DepMap gene effect has already been sign-flipped in Phase 7B).

The primary calculation uses gseapy.prerank with a fixed seed.  The outputs
separate universal convergence (same direction across trajectories) from
trajectory/subtype patterns; no drug is screened here.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
WORK = ROOT / "work"
RAW = WORK / "phase7c_functional_convergence" / "raw"
DEP = OUT / "phase7b_trajectory_gene_dependency.csv"
TRAJECTORIES = [
    "GSE77932|HCT116", "GSE77932|DLD1", "GSE42387|HCT116",
    "GSE42387|HT29", "GSE42387|LoVo", "GSE119603|HCT116",
]
SEED = 20260820


def clean_symbol(x: object) -> str:
    s = str(x).strip().upper()
    s = re.sub(r"\s*\([^)]*\)$", "", s)
    return s


def bh(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjustment, preserving the original index."""
    x = pd.to_numeric(p, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    ok = x.notna()
    if not ok.any():
        return out
    vals = x[ok].to_numpy(float)
    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    tmp = np.empty_like(ranked)
    tmp[order] = ranked
    out.loc[ok] = tmp
    return out


def parse_gmt(path: Path, prefix: str, min_size: int, max_size: int) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            genes = {clean_symbol(x) for x in fields[2:] if clean_symbol(x) not in {"", "NAN"}}
            if min_size <= len(genes) <= max_size:
                sets[f"{prefix}::{fields[0]}"] = genes
    return sets


def load_corum(path: Path) -> dict[str, set[str]]:
    """Read CORUM's tab-delimited human complex file.

    The downloaded 2019 file is used as a frozen, openly available fallback
    because the current CORUM FastAPI host has an incomplete TLS chain in this
    environment.  Complex membership is stable enough for this module-level
    sensitivity analysis, and the manifest records the exact file.
    """
    sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("Organism", "").strip().lower() not in {"human", "homo sapiens (human)"}:
                continue
            name = row.get("ComplexName", "").strip()
            genes_raw = row.get("subunits(Gene name)", "")
            genes = {clean_symbol(x) for x in genes_raw.split(";") if clean_symbol(x) not in {"", "NONE", "NAN"}}
            if len(genes) >= 3:
                cid = row.get("ComplexID", "NA").strip()
                sets[f"CORUM::{cid}::{name}"] = genes
    return sets


def load_coessentiality(path: Path, min_size: int, max_size: int) -> dict[str, set[str]]:
    """Read ClusterONE modules from the public coessentiality repository."""
    sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        try:
            gene_col = header.index("Genes")
        except ValueError:
            gene_col = len(header) - 1
        for row in reader:
            if len(row) <= gene_col:
                continue
            genes = {clean_symbol(x) for x in row[gene_col:] if clean_symbol(x) not in {"", "NAN"}}
            if min_size <= len(genes) <= max_size:
                cid = row[0] if row else "NA"
                go = row[2] if len(row) > 2 else ""
                sets[f"COESSENTIALITY::{cid}::{go}"] = genes
    return sets


def load_custom() -> dict[str, set[str]]:
    sys.path.insert(0, str(WORK / "scripts"))
    import phase7b_trajectory_conditioned_dependency as p7b
    return {f"CURATED::{k}": {clean_symbol(g) for g in v} for k, v in p7b.module_sets().items()}


def load_gene_sets() -> tuple[dict[str, set[str]], dict[str, int]]:
    sets: dict[str, set[str]] = {}
    counts: dict[str, int] = {}
    sources = {
        "REACTOME": WORK / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
        "HALLMARK": WORK / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt",
    }
    for prefix, path in sources.items():
        if path.exists():
            x = parse_gmt(path, prefix, min_size=5, max_size=500)
            sets.update(x)
            counts[prefix] = len(x)
    custom = load_custom()
    sets.update(custom)
    counts["CURATED"] = len(custom)
    corum = RAW / "corum_allComplexes_2019.txt"
    if corum.exists():
        x = load_corum(corum)
        sets.update(x)
        counts["CORUM"] = len(x)
    coess = RAW / "coessentiality_repo" / "clusterOne_clusters.tsv"
    if coess.exists():
        x = load_coessentiality(coess, min_size=5, max_size=500)
        sets.update(x)
        counts["COESSENTIALITY"] = len(x)
    return sets, counts


def load_rankings() -> dict[str, pd.Series]:
    d = pd.read_csv(DEP, usecols=["trajectory", "gene", "vulnerability_rho"])
    d["gene"] = d["gene"].map(clean_symbol)
    d["vulnerability_rho"] = pd.to_numeric(d["vulnerability_rho"], errors="coerce")
    d = d.dropna(subset=["gene", "vulnerability_rho"]).drop_duplicates(["trajectory", "gene"])
    out: dict[str, pd.Series] = {}
    for traj in TRAJECTORIES:
        x = d.loc[d["trajectory"].eq(traj), ["gene", "vulnerability_rho"]].set_index("gene")["vulnerability_rho"]
        # Spearman correlations over only 80 CRC models are highly discrete,
        # so exact ties are common.  GSEA requires an ordering; use a fixed,
        # auditable alphabetical tie-break with a numerically negligible
        # epsilon.  Without this, the order of tied genes can vary between
        # runs and can create artificial leading-edge differences.
        x = x.groupby(level=0).mean().sort_index()
        x = x.sort_values(ascending=False, kind="mergesort")
        eps = np.linspace(1e-10, 1e-12, len(x), dtype=float)
        out[traj] = pd.Series(x.to_numpy(float) + eps, index=x.index, name="vulnerability_rho")
    return out


def run_gsea(rankings: dict[str, pd.Series], gene_sets: dict[str, set[str]],
             permutations: int, seed: int) -> pd.DataFrame:
    try:
        import gseapy as gp
    except ImportError as exc:
        raise RuntimeError("gseapy is required for preranked Phase 7C") from exc
    rows: list[pd.DataFrame] = []
    for i, traj in enumerate(TRAJECTORIES):
        ranked = rankings[traj].rename_axis("gene").reset_index(name="score")
        ranked["gene"] = ranked["gene"].map(clean_symbol)
        # gseapy accepts a dict of gene sets and keeps the full ranking.
        result = gp.prerank(
            rnk=ranked,
            gene_sets=gene_sets,
            min_size=5,
            max_size=500,
            permutation_num=permutations,
            weight=1,
            ascending=False,
            threads=1,
            seed=seed + i,
            verbose=False,
            outdir=None,
            no_plot=True,
            format="png",
        )
        table = result.res2d.copy()
        table = table.rename(columns={"Term": "module", "ES": "es", "NES": "nes", "NOM p-val": "nominal_p", "FDR q-val": "fdr_qval", "Lead_genes": "leading_edge", "Tag %": "tag_fraction", "Gene %": "gene_fraction"})
        table["trajectory"] = traj
        table["collection"] = table["module"].astype(str).str.split("::", n=1).str[0]
        rows.append(table)
    allres = pd.concat(rows, ignore_index=True)
    # Recompute within-trajectory BH from nominal p as a transparent secondary
    # adjustment; gseapy's FDR is retained in fdr_qval.
    allres["bh_nominal_p"] = allres.groupby("trajectory")["nominal_p"].transform(bh)
    return allres


def convergence_tables(gsea: pd.DataFrame, rankings: dict[str, pd.Series]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    x = gsea.pivot_table(index="module", columns="trajectory", values="nes", aggfunc="first")
    x = x.reindex(columns=TRAJECTORIES)
    rows = []
    for module, row in x.iterrows():
        vals = pd.to_numeric(row, errors="coerce").dropna()
        positive = int((vals > 0).sum())
        negative = int((vals < 0).sum())
        rows.append({
            "module": module,
            "collection": str(module).split("::", 1)[0],
            "n_trajectories": len(vals),
            "n_positive": positive,
            "n_negative": negative,
            "consistency_fraction": max(positive, negative) / len(vals) if len(vals) else np.nan,
            "direction": "positive" if positive >= negative else "negative",
            "median_nes": float(vals.median()) if len(vals) else np.nan,
            "mean_nes": float(vals.mean()) if len(vals) else np.nan,
            "min_abs_nes": float(vals.abs().min()) if len(vals) else np.nan,
            "n_fdr_qval_025": int(((gsea.loc[gsea["module"].eq(module), "fdr_qval"] <= 0.25)).sum()),
            "universal_consistent_5of6": bool(len(vals) == 6 and max(positive, negative) >= 5),
            "universal_consistent_4of6": bool(len(vals) == 6 and max(positive, negative) >= 4),
        })
    conv = pd.DataFrame(rows).sort_values(["universal_consistent_5of6", "consistency_fraction", "median_nes"], ascending=[False, False, False])

    # Two-state exploratory clustering of trajectory profiles.  This is not
    # forced into the biological claims; it is a compact way to expose subtype
    # structure for later validation.
    usable = x.dropna(thresh=4).copy()
    if len(usable) >= 3:
        signal = usable.fillna(0.0).to_numpy(float).T
        if signal.shape[0] >= 3 and np.isfinite(signal).all():
            dist = pdist(signal, metric="correlation")
            z = linkage(dist, method="average")
            labels = fcluster(z, t=2, criterion="maxclust")
            groups = {traj: int(label) for traj, label in zip(usable.columns, labels)}
        else:
            groups = {traj: 1 for traj in usable.columns}
    else:
        groups = {traj: 1 for traj in TRAJECTORIES}
    pattern_rows = []
    for module, row in x.iterrows():
        for group in sorted(set(groups.values())):
            members = [t for t, g in groups.items() if g == group]
            other = [t for t, g in groups.items() if g != group]
            a = pd.to_numeric(row.reindex(members), errors="coerce").dropna()
            b = pd.to_numeric(row.reindex(other), errors="coerce").dropna()
            if len(a) and len(b):
                pattern_rows.append({"module": module, "subtype": f"cluster_{group}", "members": ";".join(members), "n_members": len(a), "mean_nes_subtype": float(a.mean()), "mean_nes_other": float(b.mean()), "delta_nes_subtype_vs_other": float(a.mean() - b.mean()), "subtype_direction": "higher_vulnerability" if a.mean() > b.mean() else "lower_vulnerability"})
    return gsea, conv, pd.DataFrame(pattern_rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--permutations", type=int, default=100)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rankings = load_rankings()
    gene_sets, counts = load_gene_sets()
    gsea = run_gsea(rankings, gene_sets, args.permutations, args.seed)
    gsea, conv, patterns = convergence_tables(gsea, rankings)
    gsea.to_csv(OUT / "phase7c_trajectory_module_gsea.csv", index=False)
    conv.to_csv(OUT / "phase7c_module_convergence_ranking.csv", index=False)
    nes = gsea.pivot_table(index="module", columns="trajectory", values="nes", aggfunc="first").reindex(columns=TRAJECTORIES)
    nes.to_csv(OUT / "phase7c_trajectory_module_nes_matrix.csv")
    patterns.to_csv(OUT / "phase7c_subtype_module_patterns.csv", index=False)
    manifest = {
        "phase": "7C",
        "seed": args.seed,
        "permutations": args.permutations,
        "dependency_input": str(DEP),
        "n_ranked_genes_by_trajectory": {k: int(v.size) for k, v in rankings.items()},
        "tie_handling": "exact Spearman-rho ties were broken deterministically by alphabetical gene order using epsilon 1e-10..1e-12; this is required for a reproducible preranked order and does not change score magnitudes materially",
        "gene_set_counts": counts,
        "gene_set_sources": {
            "Reactome/Hallmark": "local MSigDB v2026.1 GMT files",
            "CURATED": "phase3/phase7b mechanistic modules",
            "CORUM": "corum_allComplexes_2019.txt from mal2017/corumdata; frozen fallback due current CORUM TLS issue",
            "COESSENTIALITY": "kundajelab/coessentiality clusterOne_clusters.tsv",
        },
        "note": "No drug screening or single-gene gate was applied in Phase 7C; module-level convergence is exploratory and must be validated in paired parental/OXA-R functional data.",
    }
    (OUT / "phase7c_functional_module_convergence_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    top = conv.head(20)
    lines = ["# Phase 7C：OXA-R trajectory-conditioned functional-module convergence", "", "## 分析定义", "", "- 输入：6 条 parental→OXA-R trajectory 在 80 个 DepMap CRC 模型上的完整 gene-level vulnerability ranking。", "- 方法：对 Reactome、Hallmark、curated modules、CORUM protein complexes 和 public co-essentiality modules 做 preranked GSEA。", "- 正方向：NES > 0 表示更像该 OXA-R trajectory 的模型对模块内基因总体更依赖。", "- 这一阶段不做药物筛选，也不把单个基因命中写成机制结论。", "", "## Gene-set universe", "", json.dumps(counts, ensure_ascii=False), "", "## Universal convergence candidates", "", "| Module | Collection | + trajectories | - trajectories | consistency | median NES | FDR≤0.25 hits |", "|---|---:|---:|---:|---:|---:|---:|"]
    for _, r in top.iterrows():
        lines.append(f"| {r['module']} | {r['collection']} | {int(r['n_positive'])} | {int(r['n_negative'])} | {r['consistency_fraction']:.2f} | {r['median_nes']:.2f} | {int(r['n_fdr_qval_025'])} |")
    lines += ["", "## Interpretation guardrails", "", "- `universal_consistent_5of6` 和 `universal_consistent_4of6` 是方向一致性描述，不等于 acquired OXA-R causality。", "- CORUM/co-essentiality 模块用于把分散的 gene-level signals 聚合到功能层；它们不能替代 paired parental/OXA-R CRISPR。", "- subtype patterns are exploratory cluster contrasts and require independent models or paired screens for confirmation.", ""]
    (OUT / "phase7c_functional_module_convergence.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"gene_sets": counts, "gsea_rows": int(len(gsea)), "modules": int(len(conv)), "outputs": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
