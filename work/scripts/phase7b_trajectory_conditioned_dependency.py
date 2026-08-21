"""Phase 7B: OXA-R trajectory-conditioned dependency mapping.

The script has two stages:
  (1) always available: export six separate, directional OXA-R signatures and
      their explanatory module scores;
  (2) executable when DepMap files are supplied: project the signatures to
      CRC cell lines and correlate each trajectory score with CRISPR gene
      effect scores.

The dependency convention is explicit: DepMap gene-effect values are usually
negative for stronger dependency. We report `vulnerability_rho = -rho`, so a
positive value means that higher OXA-R-state score predicts stronger gene
dependency.

Example once local DepMap files are available:
  python phase7b_trajectory_conditioned_dependency.py \
      --expr path/to/CCLE_expression.csv \
      --dependency path/to/CRISPRGeneEffect.csv \
      --sample-info path/to/sample_info.csv
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
WORK = ROOT / "work"
OUT.mkdir(exist_ok=True)

GENE_DELTA = WORK / "phase5_perturbation_reversal" / "gene_delta_matrix_primary.csv"
TRAJECTORIES = [
    "GSE77932|HCT116", "GSE77932|DLD1", "GSE42387|HCT116",
    "GSE42387|HT29", "GSE42387|LoVo", "GSE119603|HCT116",
]


def clean_symbol(x: object) -> str:
    s = str(x).strip().upper()
    s = re.sub(r"\s*\([^)]*\)$", "", s)
    return s


def bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    ix = np.where(ok)[0]
    order = ix[np.argsort(p[ix])]
    q = p[order] * len(order) / np.arange(1, len(order) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1].clip(0, 1)
    out[order] = q
    return out


def load_trajectories() -> pd.DataFrame:
    d = pd.read_csv(GENE_DELTA)
    gene_col = d.columns[0]
    d["gene"] = d[gene_col].map(clean_symbol)
    d = d[d["gene"].ne("") & d["gene"].ne("NAN")].copy()
    d = d.drop_duplicates("gene").set_index("gene")
    cols = [c for c in TRAJECTORIES if c in d.columns]
    if len(cols) != len(TRAJECTORIES):
        raise RuntimeError(f"missing trajectories: {sorted(set(TRAJECTORIES) - set(cols))}")
    return d[cols].apply(pd.to_numeric, errors="coerce")


def export_signatures(delta: pd.DataFrame, top_n: int = 250) -> None:
    long = []
    meta = []
    for trajectory in delta.columns:
        s = delta[trajectory].dropna().sort_values(ascending=False)
        up = s.head(top_n)
        down = s.tail(top_n).sort_values()
        for direction, part in [("up", up), ("down", down)]:
            scale = max(float(part.abs().max()), 1e-12)
            for rank, (gene, weight) in enumerate(part.items(), 1):
                long.append({"trajectory": trajectory, "gene": gene,
                             "direction": direction, "rank": rank,
                             "delta": float(weight),
                             "weight_scaled": float(weight / scale)})
        meta.append({"trajectory": trajectory, "n_genes_total": int(s.notna().sum()),
                     "n_up": int((s > 0).sum()), "n_down": int((s < 0).sum()),
                     "top_n_each_direction": top_n,
                     "signature_rule": "top directional gene-delta genes; trajectories kept separate"})
    pd.DataFrame(long).to_csv(OUT / "phase7b_oxa_r_trajectory_signatures.csv", index=False)
    pd.DataFrame(meta).to_csv(OUT / "phase7b_trajectory_signature_metadata.csv", index=False)


def module_sets() -> dict[str, set[str]]:
    sys.path.insert(0, str(ROOT / "work" / "scripts"))
    import phase3_module_decomposition as phase3
    sets, _ = phase3.build_decomposition_sets()
    # Add the modules that were explicitly retained in the project design.
    sets["FAO_mitochondrial"] = {"CPT1A", "CPT1B", "CPT2", "ACADVL", "ACADM", "ACADL", "HADHA", "HADHB"}
    sets["carnitine_entry"] = {"BBOX1", "SLC22A5", "CPT1A", "CPT1B", "CPT2", "SLC25A20"}
    sets["NRF2_redox"] = {"NFE2L2", "KEAP1", "NQO1", "HMOX1", "GCLC", "GCLM", "GSS", "SLC7A11", "TXNRD1", "GPX4", "AKR1C1", "GSTP1"}
    sets["ferroptosis_resistance"] = {"SLC7A11", "GPX4", "AIFM2", "FSP1", "DHODH", "GCH1", "FTH1", "FTL", "NFE2L2", "ACSL3", "SCD", "LPCAT3", "TFRC", "SAT1"}
    sets["ABC_transport"] = {"ABCB1", "ABCC1", "ABCC2", "ABCC3", "ABCC4", "ABCG2"}
    return {k: {clean_symbol(x) for x in v} for k, v in sets.items()}


def export_module_trajectory_scores(delta: pd.DataFrame) -> None:
    rows = []
    for module, genes in module_sets().items():
        present = sorted(set(delta.index) & genes)
        for trajectory in delta.columns:
            values = delta.loc[present, trajectory].dropna() if present else pd.Series(dtype=float)
            rows.append({"trajectory": trajectory, "module": module,
                         "n_genes_present": len(values),
                         "module_delta_median": float(values.median()) if len(values) else np.nan,
                         "module_delta_mean": float(values.mean()) if len(values) else np.nan,
                         "module_genes_present": ";".join(values.index)})
    pd.DataFrame(rows).to_csv(OUT / "phase7b_module_trajectory_scores.csv", index=False)


def orient_expression(expr: pd.DataFrame, sample_info: pd.DataFrame | None = None) -> pd.DataFrame:
    """Return genes x samples; tolerate common DepMap orientations."""
    expr = expr.copy()
    raw_index = expr.index.astype(str)
    raw_columns = pd.Index(expr.columns).astype(str)
    info_ids: set[str] = set()
    if sample_info is not None:
        low = {c.lower(): c for c in sample_info.columns}
        id_col = next((low[k] for k in ["depmap_id", "modelid", "sample_id"] if k in low), None)
        if id_col:
            info_ids = set(sample_info[id_col].astype(str))
    row_id_overlap = len(set(raw_index) & info_ids)
    col_id_overlap = len(set(raw_columns) & info_ids)
    looks_like_depmap_rows = np.mean(raw_index.str.match(r"^(ACH|SIDM|SNU|MEL|TM|CVCL)-", case=False)) > 0.1
    # DepMap omics files conventionally have rows=Profile/Model IDs and
    # columns=genes; transpose those to genes x samples.
    if row_id_overlap > col_id_overlap or looks_like_depmap_rows:
        expr = expr.T
    expr.index = expr.index.map(clean_symbol)
    expr = expr.loc[~expr.index.duplicated()].apply(pd.to_numeric, errors="coerce")
    return expr


def crc_sample_ids(sample_info: pd.DataFrame, available: set[str]) -> list[str]:
    info = sample_info.copy()
    cols = {c.lower(): c for c in info.columns}
    id_col = next((cols[k] for k in ["depmap_id", "modelid", "sample_id"] if k in cols), info.columns[0])
    text = info.fillna("").agg(lambda row: " ".join(map(str, row)), axis=1).str.lower()
    mask = text.str.contains(r"colon|colorectal|rectum|large intestine|coad|read|coread", regex=True, na=False)
    ids = info.loc[mask, id_col].astype(str)
    return [x for x in ids if x in available]


def map_expression_to_model_ids(sample_info: pd.DataFrame, available: set[str]) -> dict[str, str]:
    """Map whichever expression profile IDs are present to DepMap ModelID."""
    low = {c.lower(): c for c in sample_info.columns}
    model_col = next((low[k] for k in ["depmap_id", "modelid", "model_id"] if k in low), None)
    if model_col is None:
        return {}
    candidates = []
    for col in sample_info.columns:
        overlap = len(set(sample_info[col].astype(str)) & available)
        candidates.append((overlap, col))
    profile_col = max(candidates)[1] if candidates else model_col
    text = sample_info.fillna("").agg(lambda row: " ".join(map(str, row)), axis=1).str.lower()
    crc_mask = text.str.contains(r"colon|colorectal|rectum|large intestine|coad|read|coread", regex=True, na=False)
    if profile_col == model_col:
        sub = sample_info.loc[crc_mask, [profile_col]].dropna().astype(str)
        sub = sub[sub[profile_col].isin(available)]
        return {x: x for x in sub[profile_col].tolist()}
    sub = sample_info.loc[crc_mask, [profile_col, model_col]].dropna().astype(str)
    sub = sub[sub[profile_col].isin(available)]
    return dict(zip(sub[profile_col], sub[model_col]))


def project_expression(expr_path: Path, sample_info_path: Path | None, delta: pd.DataFrame) -> pd.DataFrame:
    expr = pd.read_csv(expr_path, index_col=0, low_memory=False)
    info = pd.read_csv(sample_info_path, low_memory=False) if sample_info_path and sample_info_path.exists() else None
    expr = orient_expression(expr, info)
    if sample_info_path and sample_info_path.exists():
        assert info is not None
        mapping = map_expression_to_model_ids(info, set(expr.columns))
        expr = expr.loc[:, [x for x in expr.columns if x in mapping]]
        expr.columns = [mapping[x] for x in expr.columns]
        expr = expr.T.groupby(level=0).mean().T
    else:
        # No metadata means no lineage restriction; this is explicitly labelled.
        expr = expr.loc[:, [c for c in expr.columns if c]]
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0)
    weights = delta.div(delta.abs().max(axis=0).replace(0, np.nan), axis=1)
    common = sorted(set(z.index) & set(weights.index))
    scores = {}
    for trajectory in weights.columns:
        w = weights.loc[common, trajectory]
        scores[trajectory] = z.loc[common].mul(w, axis=0).sum(axis=0) / w.abs().sum()
    out = pd.DataFrame(scores)
    out.index.name = "sample_id"
    out.to_csv(OUT / "phase7b_crc_trajectory_scores.csv")
    return out


def load_dependency(path: Path) -> pd.DataFrame:
    dep = pd.read_csv(path, index_col=0, low_memory=False)
    dep.index = dep.index.astype(str)
    dep.columns = dep.columns.map(clean_symbol)
    dep = dep.loc[:, ~dep.columns.duplicated()].apply(pd.to_numeric, errors="coerce")
    return dep


def resampling_evidence(x: pd.Series, y: pd.Series, n_resamples: int, rng: np.random.Generator) -> tuple[float, float, float]:
    """Permutation p and bootstrap CI for the vulnerability-oriented rho."""
    xv, yv = x.to_numpy(float), y.to_numpy(float)
    rho = float(spearmanr(xv, yv).statistic)
    if n_resamples <= 0 or len(xv) < 5:
        return np.nan, np.nan, np.nan
    # Null: preserve the state scores and permute dependency values.
    null = np.empty(n_resamples, dtype=float)
    for i in range(n_resamples):
        null[i] = float(spearmanr(xv, rng.permutation(yv)).statistic)
    perm_p = (1 + int(np.sum(np.abs(null) >= abs(rho)))) / (n_resamples + 1)
    # vulnerability_rho is -rho because more negative gene effect means more dependency.
    boot = np.empty(n_resamples, dtype=float)
    for i in range(n_resamples):
        ix = rng.integers(0, len(xv), len(xv))
        boot[i] = -float(spearmanr(xv[ix], yv[ix]).statistic)
    return float(perm_p), float(np.nanquantile(boot, 0.025)), float(np.nanquantile(boot, 0.975))


def dependency_correlations(scores: pd.DataFrame, dependency: pd.DataFrame,
                            n_resamples: int = 500, seed: int = 20260820) -> tuple[pd.DataFrame, pd.DataFrame]:
    common = scores.index.intersection(dependency.index)
    scores = scores.loc[common]
    dep = dependency.loc[common]
    by_trajectory = []
    for trajectory in scores.columns:
        x = scores[trajectory]
        for gene in dep.columns:
            y = dep[gene]
            ok = x.notna() & y.notna()
            if ok.sum() < 15:
                continue
            rho, p = spearmanr(x[ok], y[ok])
            by_trajectory.append({"trajectory": trajectory, "gene": gene,
                                  "n": int(ok.sum()), "rho_state_vs_dep": rho,
                                  "vulnerability_rho": -rho, "p": p})
    res = pd.DataFrame(by_trajectory)
    if len(res):
        res["p_BH_within_trajectory"] = res.groupby("trajectory")["p"].transform(lambda x: bh(x.to_numpy()))
        # Analytical p-values are used for the genome-wide screen. For a small
        # top set in each trajectory, add permutation p-values and bootstrap
        # confidence intervals so the final ranking is not driven by asymptotic
        # p-values in a ~40-50-line panel.
        res["permutation_p_top"] = np.nan
        res["vulnerability_rho_boot_low"] = np.nan
        res["vulnerability_rho_boot_high"] = np.nan
        rng = np.random.default_rng(seed)
        for trajectory in scores.columns:
            top_ix = (res.loc[res["trajectory"].eq(trajectory), "vulnerability_rho"]
                      .abs().nlargest(50).index)
            for ix in top_ix:
                gene = res.at[ix, "gene"]
                x = scores[trajectory]
                y = dependency[gene]
                ok = x.notna() & y.notna()
                pp, lo, hi = resampling_evidence(x[ok], y[ok], n_resamples, rng)
                res.at[ix, "permutation_p_top"] = pp
                res.at[ix, "vulnerability_rho_boot_low"] = lo
                res.at[ix, "vulnerability_rho_boot_high"] = hi
        conv = (res.groupby("gene")
                .agg(n_trajectories=("trajectory", "nunique"),
                     median_vulnerability_rho=("vulnerability_rho", "median"),
                     mean_vulnerability_rho=("vulnerability_rho", "mean"),
                     n_positive=("vulnerability_rho", lambda x: int((x > 0).sum())),
                     min_p_BH=("p_BH_within_trajectory", "min"),
                     n_bootstrap_supported=("permutation_p_top", lambda x: int((x < 0.05).sum())))
                .reset_index())
        conv["directional_convergence_fraction"] = conv["n_positive"] / conv["n_trajectories"]
        # A broad directional screen is useful for exploration, but the
        # primary claim requires stronger convergence. The strict gate is
        # deliberately hard: >=5/6 directions positive, median effect >=0.10,
        # >=2 trajectory-level permutation hits, and at least one within-
        # trajectory BH q-value <=0.10.
        conv["strict_convergent_vulnerability"] = (
            (conv["n_trajectories"] >= 5) &
            (conv["n_positive"] >= 5) &
            (conv["median_vulnerability_rho"] >= 0.10) &
            (conv["n_bootstrap_supported"] >= 2) &
            (conv["min_p_BH"] <= 0.10))
        conv["exploratory_convergent_signal"] = (
            (conv["n_trajectories"] >= 5) &
            (conv["n_positive"] >= 4) &
            (conv["median_vulnerability_rho"] >= 0.10) &
            (conv["n_bootstrap_supported"] >= 2))
        conv["convergent_vulnerability"] = conv["strict_convergent_vulnerability"]
    else:
        conv = pd.DataFrame()
    return res, conv


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", type=Path)
    ap.add_argument("--dependency", type=Path)
    ap.add_argument("--sample-info", type=Path)
    ap.add_argument("--resamples", type=int, default=500,
                    help="permutation/bootstrap resamples for top 50 genes per trajectory")
    ap.add_argument("--seed", type=int, default=20260820)
    args = ap.parse_args()

    delta = load_trajectories()
    export_signatures(delta)
    export_module_trajectory_scores(delta)
    manifest = {"n_trajectories": len(delta.columns), "trajectories": list(delta.columns),
                "gene_delta_file": str(GENE_DELTA), "dependency_analysis": "not_run"}
    if args.expr and args.dependency:
        scores = project_expression(args.expr, args.sample_info, delta)
        dependency = load_dependency(args.dependency)
        by_traj, conv = dependency_correlations(scores, dependency, args.resamples, args.seed)
        by_traj.to_csv(OUT / "phase7b_trajectory_gene_dependency.csv", index=False)
        conv.to_csv(OUT / "phase7b_convergent_dependency_ranking.csv", index=False)
        mechanistic_genes = {"DHODH", "VCP", "RRM2", "CAD", "UMPS", "CPT1A", "CPT2", "BBOX1", "SLC22A5", "SLC25A20", "SERPINE1", "FN1", "TGFBR1", "TGFBR2", "NFE2L2", "GPX4", "SLC7A11", "DERL1", "EDEM1", "HSP90B1", "MANF", "HSPA5", "PSMA1", "PSMB5"}
        conv[conv["gene"].isin(mechanistic_genes)].to_csv(OUT / "phase7b_mechanistic_dependency_audit.csv", index=False)
        manifest.update({"dependency_analysis": "completed", "n_crc_samples": int(scores.shape[0]),
                         "n_dependency_genes": int(dependency.shape[1]),
                         "n_convergent_genes": int(conv["convergent_vulnerability"].sum()) if len(conv) else 0,
                         "n_exploratory_convergent_signals": int(conv["exploratory_convergent_signal"].sum()) if len(conv) else 0})
    else:
        manifest["dependency_analysis_note"] = "Prepared signatures only; supply --expr and --dependency for CRISPR mapping."
    (OUT / "phase7b_trajectory_dependency_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    report = [
        "# Phase 7B：OXA-R state-conditioned dependency mapping",
        "",
        f"已固定 {len(delta.columns)} 条独立 OXA-R resistance trajectories；没有把它们合并成一个 consensus signature。",
        "",
        "## 本轮已完成",
        "",
        "- 为每条 trajectory 导出独立的 top-up/top-down directional signature；",
        "- 导出 FAO、carnitine-entry、pyrimidine、UPR/ERAD相关、EMT、NRF2/redox、ferroptosis、ABC 等解释模块的 trajectory-level effect；",
        "- 固定 DepMap CRISPR 的方向转换：原始 gene-effect 越负代表越依赖，报告中的 `vulnerability_rho = -rho` 越大代表 OXA-R-like state 越依赖该基因。",
        "",
        "## 关键统计门槛",
        "",
        "- 每条 trajectory 单独计算 Score；不把 HCT116、DLD1、HT29、LoVo 的轨迹强行合并。",
        "- 每个 trajectory-gene 至少需要 15 个 CRC cell lines；严格 convergent vulnerability 要求至少 5/6 条轨迹方向一致、中位 vulnerability_rho ≥ 0.10、至少 2 条 trajectory 通过 permutation/bootstrap，且至少一条 trajectory 的 BH-q ≤ 0.10。",
        "- 最终还需要 permutation/bootstrap 和药敏闭环；当前不能把转录投影结果称为功能依赖。",
        "",
        "## 文件",
        "",
        "- `phase7b_oxa_r_trajectory_signatures.csv`",
        "- `phase7b_trajectory_signature_metadata.csv`",
        "- `phase7b_module_trajectory_scores.csv`",
        "- `phase7b_trajectory_dependency_manifest.json`",
    ]
    if manifest["dependency_analysis"] == "completed":
        report += ["- `phase7b_crc_trajectory_scores.csv`", "- `phase7b_trajectory_gene_dependency.csv`", "- `phase7b_convergent_dependency_ranking.csv`", "- `phase7b_mechanistic_dependency_audit.csv`", "", f"Top 50 genes per trajectory additionally receive {args.resamples} permutation/bootstrap resamples.", "", f"严格 convergent genes: {manifest['n_convergent_genes']}；exploratory signals: {manifest['n_exploratory_convergent_signals']}。"]
    else:
        report += ["", "CRISPR dependency 尚未运行：本地当前没有可用的 DepMap expression 与 CRISPRGeneEffect 文件。"]
    (OUT / "phase7b_trajectory_conditioned_dependency.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
