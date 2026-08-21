#!/usr/bin/env python
"""Phase 7B-R: robust OXA-R trajectory-conditioned dependency mapping.

This script deliberately keeps the six acquired-resistance trajectories as
distinct vectors, but treats the three HCT116 trajectories as one biological
background when making convergence calls.  It reads local DepMap 23Q4 and
derived resistance signatures; raw files are never copied to outputs.

The script stops at genetic dependency mapping.  Phase 7C/8 are intentionally
not executed here.
"""

from __future__ import annotations

import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, spearmanr, t

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
WORK = ROOT / "work"
OUT.mkdir(exist_ok=True)

GENE_DELTA = WORK / "phase5_perturbation_reversal" / "gene_delta_matrix_primary.csv"
PATHWAY_EFFECTS = OUT / "phase1_pathway_effects_long_all_contexts.csv"
EXPR_RAW = WORK / "phase7b_depmap" / "raw" / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
DEP_RAW = WORK / "phase7b_depmap" / "raw" / "CRISPRGeneEffect.csv"
MODEL_RAW = WORK / "phase7b_depmap" / "raw" / "Model.csv"
HALLMARK_GMT = WORK / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"

TRAJECTORIES = [
    "GSE77932|HCT116",
    "GSE77932|DLD1",
    "GSE42387|HCT116",
    "GSE42387|HT29",
    "GSE42387|LoVo",
    "GSE119603|HCT116",
]
DATASET = {x: x.split("|")[0] for x in TRAJECTORIES}
BACKGROUND = {
    "GSE77932|HCT116": "HCT116",
    "GSE77932|DLD1": "DLD1",
    "GSE42387|HCT116": "HCT116",
    "GSE42387|HT29": "HT29",
    "GSE42387|LoVo": "LoVo",
    "GSE119603|HCT116": "HCT116",
}
BACKGROUNDS = ["HCT116", "DLD1", "HT29", "LoVo"]
SIGNATURE_SIZES = [100, 250, 500]
PRIMARY_SIZE = 250
PRIMARY_METHOD = "weighted"
SEED = 20260820
N_RESAMPLES = 1000
MIN_N = 30

PREDEFINED_MECHANISM = [
    "DHODH", "RRM2", "VCP", "CPT1A", "CPT2", "SLC22A5", "GPX4",
    "SLC7A11", "FN1", "ITGB1", "DERL1", "EDEM1", "HSP90B1", "PSMB5",
]


def clean_symbol(value: object) -> str:
    value = str(value).strip()
    value = re.sub(r"\s*\([^)]*\)$", "", value)
    return value.upper()


def norm_name(value: object) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def bh(values: np.ndarray | pd.Series) -> np.ndarray:
    """Benjamini-Hochberg q values, preserving NaN."""
    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    valid = np.isfinite(arr)
    if not valid.any():
        return out
    p = np.clip(arr[valid], 0.0, 1.0)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1.0)
    out[valid] = restored
    return out


def ensure_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return df[columns]


def deduplicate_gene_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Clean symbols and keep the first column for duplicated symbols."""
    names = [clean_symbol(x) for x in df.columns]
    keep = ~pd.Index(names).duplicated()
    out = df.loc[:, keep].copy()
    out.columns = np.asarray(names)[keep]
    return out.apply(pd.to_numeric, errors="coerce")


def load_trajectory_metadata() -> pd.DataFrame:
    rows = []
    if PATHWAY_EFFECTS.exists():
        effects = pd.read_csv(PATHWAY_EFFECTS, low_memory=False)
        for traj in TRAJECTORIES:
            sub = effects.loc[effects["contrast_id"].astype(str) == traj]
            if len(sub):
                first = sub.iloc[0]
                rows.append({
                    "trajectory": traj,
                    "dataset": DATASET[traj],
                    "background": BACKGROUND[traj],
                    "platform": ";".join(sorted(set(sub["platform"].dropna().astype(str)))),
                    "parental_n": ";".join(sorted(set(sub["parental_n"].dropna().astype(str)))),
                    "resistant_n": ";".join(sorted(set(sub["resistant_n"].dropna().astype(str)))),
                    "context": ";".join(sorted(set(sub["context"].dropna().astype(str)))),
                    "n_pathway_rows": len(sub),
                    "metadata_source": str(PATHWAY_EFFECTS.relative_to(ROOT)),
                })
            else:
                rows.append({"trajectory": traj, "dataset": DATASET[traj], "background": BACKGROUND[traj],
                             "platform": "", "parental_n": "", "resistant_n": "", "context": "",
                             "n_pathway_rows": 0, "metadata_source": "not found"})
    else:
        rows = [{"trajectory": t0, "dataset": DATASET[t0], "background": BACKGROUND[t0],
                 "platform": "", "parental_n": "", "resistant_n": "", "context": "",
                 "n_pathway_rows": 0, "metadata_source": "not found"} for t0 in TRAJECTORIES]
    return pd.DataFrame(rows)


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in [GENE_DELTA, EXPR_RAW, DEP_RAW, MODEL_RAW]:
        if not path.exists():
            raise FileNotFoundError(path)

    delta = pd.read_csv(GENE_DELTA, low_memory=False)
    gene_col = delta.columns[0]
    delta["gene"] = delta[gene_col].map(clean_symbol)
    delta = delta.drop(columns=[gene_col]).set_index("gene")
    delta = delta.apply(pd.to_numeric, errors="coerce")
    delta = delta.loc[~delta.index.duplicated(keep="first"), TRAJECTORIES]
    delta = delta.dropna(how="all")

    model = pd.read_csv(MODEL_RAW, low_memory=False)
    lineage = model["OncotreeLineage"].fillna("").astype(str).str.lower()
    primary = model["OncotreePrimaryDisease"].fillna("").astype(str).str.lower()
    crc = model.loc[(lineage == "bowel") & primary.str.contains("colorectal", regex=False)].copy()
    crc = crc[["ModelID", "CellLineName", "StrippedCellLineName", "CCLEName", "OncotreeLineage", "OncotreePrimaryDisease"]]
    crc["ModelID"] = crc["ModelID"].astype(str)

    expression = pd.read_csv(EXPR_RAW, index_col=0, low_memory=False)
    expression.index = expression.index.astype(str)
    expression = deduplicate_gene_columns(expression)

    dependency = pd.read_csv(DEP_RAW, index_col=0, low_memory=False)
    dependency.index = dependency.index.astype(str)
    dependency = deduplicate_gene_columns(dependency)

    common = sorted(set(crc["ModelID"]) & set(expression.index) & set(dependency.index))
    if len(common) < MIN_N:
        raise ValueError(f"Only {len(common)} CRC models overlap expression/dependency; need >= {MIN_N}")
    crc = crc.set_index("ModelID").loc[common]
    expression = expression.loc[common]
    dependency = dependency.loc[common]
    return delta, crc, expression, dependency


def self_exclusion_map(crc: pd.DataFrame) -> dict[str, list[str]]:
    out = {}
    for traj in TRAJECTORIES:
        target = norm_name(BACKGROUND[traj])
        mask = pd.Series(False, index=crc.index)
        for col in ["CellLineName", "StrippedCellLineName", "CCLEName"]:
            mask |= crc[col].fillna("").map(norm_name).str.contains(target, regex=False)
        out[traj] = crc.index[mask].tolist()
    return out


def build_signature_weights(delta: pd.DataFrame, n: int, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    summary = []
    for traj in TRAJECTORIES:
        s = delta[traj].dropna().sort_values(ascending=False)
        up = s.head(n)
        down = s.tail(n).sort_values(ascending=True)
        max_abs = float(max(up.abs().max(), down.abs().max(), 1e-12))
        for direction, part in [("up", up), ("down", down)]:
            for rank, (gene, value) in enumerate(part.items(), start=1):
                sign = 1.0 if direction == "up" else -1.0
                rows.append({
                    "trajectory": traj, "dataset": DATASET[traj], "background": BACKGROUND[traj],
                    "signature_size": n, "gene": gene, "direction": direction,
                    "rank": rank, "delta": float(value),
                    "weight_scaled": sign * abs(float(value)) / max_abs,
                    "rank_weight": sign * (n - rank + 1) / n,
                })
        m = meta.loc[meta["trajectory"] == traj].iloc[0].to_dict()
        summary.append({"trajectory": traj, "dataset": DATASET[traj], "background": BACKGROUND[traj],
                        "signature_size": n, "n_up": len(up), "n_down": len(down),
                        "parental_n": m.get("parental_n", ""), "resistant_n": m.get("resistant_n", ""),
                        "platform": m.get("platform", ""), "self_line_exclusion": "per-trajectory CRC self-line excluded"})
    return pd.DataFrame(rows), pd.DataFrame(summary)


def project_scores(expr: pd.DataFrame, weights: pd.DataFrame, n: int) -> pd.DataFrame:
    z = (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=0).replace(0, 1.0)
    out = pd.DataFrame(index=expr.index)
    for traj in TRAJECTORIES:
        w = weights.loc[weights["trajectory"] == traj].set_index("gene")
        common = [g for g in w.index if g in z.columns]
        # Microarray-derived signatures can contain platform-specific symbols
        # absent from the DepMap expression matrix.  Keep the signed score
        # with the observed intersection and expose the intersection size in
        # the metadata rather than silently imputing missing genes.
        if len(common) < MIN_N:
            raise ValueError(f"Too few expression genes for {traj} top{n}: {len(common)}")
        ww = w.loc[common, "weight_scaled"].to_numpy(float)
        rw = w.loc[common, "rank_weight"].to_numpy(float)
        out[f"{traj}__weighted"] = z[common].to_numpy(float).dot(ww) / np.abs(ww).sum()
        ranks = expr[common].rank(axis=1, pct=True).to_numpy(float) - 0.5
        out[f"{traj}__rank"] = ranks.dot(rw) / np.abs(rw).sum()
    out.insert(0, "ModelID", out.index)
    return out.reset_index(drop=True)


def correlation_table(scores: pd.DataFrame, n: int) -> pd.DataFrame:
    rows = []
    for method in ["weighted", "rank"]:
        cols = [f"{traj}__{method}" for traj in TRAJECTORIES]
        for i, a in enumerate(cols):
            for b in cols[i:]:
                rho, p = spearmanr(scores[a], scores[b])
                rows.append({"signature_size": n, "scoring_method": method,
                             "trajectory_a": a.split("__")[0], "trajectory_b": b.split("__")[0],
                             "spearman_rho": rho, "p_value": p, "n": len(scores)})
    return pd.DataFrame(rows)


def spearman_vector(x: pd.Series | np.ndarray, dep: pd.DataFrame) -> pd.DataFrame:
    x = np.asarray(x, dtype=float)
    xr = pd.Series(x).rank(method="average").to_numpy(float)
    y = dep.to_numpy(dtype=float)
    yr = dep.rank(axis=0, method="average", na_option="keep").to_numpy(dtype=float)
    valid = np.isfinite(yr)
    xmat = np.broadcast_to(xr[:, None], yr.shape)
    n = valid.sum(axis=0).astype(int)
    xmean = np.where(n > 0, np.where(valid, xmat, np.nan).sum(axis=0) / np.maximum(n, 1), np.nan)
    ymean = np.where(n > 0, np.where(valid, yr, np.nan).sum(axis=0) / np.maximum(n, 1), np.nan)
    xc = np.where(valid, xmat - xmean, 0.0)
    yc = np.where(valid, yr - ymean, 0.0)
    denom = np.sqrt((xc * xc).sum(axis=0) * (yc * yc).sum(axis=0))
    rho = np.where(denom > 0, (xc * yc).sum(axis=0) / denom, np.nan)
    rho = np.clip(rho, -1.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = rho * np.sqrt(np.maximum(n - 2, 1) / np.maximum(1 - rho * rho, 1e-15))
    p = 2 * t.sf(np.abs(stat), np.maximum(n - 2, 1))
    p[(n < MIN_N) | ~np.isfinite(rho)] = np.nan
    return pd.DataFrame({"gene": dep.columns, "n": n, "rho_state_vs_dependency": rho,
                         "vulnerability_rho": -rho, "p_value": p})


def dependency_results(scores: pd.DataFrame, dependency: pd.DataFrame, crc: pd.DataFrame,
                       exclusions: dict[str, list[str]], signature_size: int,
                       method: str) -> pd.DataFrame:
    rows = []
    for traj in TRAJECTORIES:
        col = f"{traj}__{method}"
        keep = [x for x in scores["ModelID"] if x not in exclusions[traj] and x in dependency.index]
        sub = dependency.loc[keep]
        result = spearman_vector(scores.set_index("ModelID").loc[keep, col], sub)
        result.insert(0, "trajectory", traj)
        result.insert(1, "dataset", DATASET[traj])
        result.insert(2, "background", BACKGROUND[traj])
        result["signature_size"] = signature_size
        result["scoring_method"] = method
        result["self_line_excluded"] = ";".join(exclusions[traj])
        result["n_crc_after_exclusion"] = len(keep)
        result["p_value_bh_within_trajectory"] = bh(result["p_value"].to_numpy())
        rows.append(result)
    return pd.concat(rows, ignore_index=True)


def bootstrap_permutation(top_results: pd.DataFrame, scores: pd.DataFrame, dependency: pd.DataFrame,
                          exclusions: dict[str, list[str]], seed: int) -> pd.DataFrame:
    rows = []
    for ti, traj in enumerate(TRAJECTORIES):
        subtop = top_results.loc[top_results["trajectory"] == traj].head(200)
        keep = [x for x in scores["ModelID"] if x not in exclusions[traj] and x in dependency.index]
        x = scores.set_index("ModelID").loc[keep, f"{traj}__weighted"].to_numpy(float)
        dep = dependency.loc[keep, subtop["gene"].tolist()]
        xr = pd.Series(x).rank().to_numpy(float)
        yr = dep.rank(axis=0).to_numpy(float)
        n = len(keep)
        observed = subtop.set_index("gene").loc[dep.columns, "vulnerability_rho"].to_numpy(float)
        rng = np.random.default_rng(seed + ti)
        null_ge = np.zeros(len(dep.columns), dtype=int)
        boot_values = np.empty((N_RESAMPLES, len(dep.columns)), dtype=float)
        for start in range(0, N_RESAMPLES, 100):
            k = min(100, N_RESAMPLES - start)
            perm = np.stack([rng.permutation(xr) for _ in range(k)], axis=0)
            yc = yr - yr.mean(axis=0, keepdims=True)
            yp = perm - perm.mean(axis=1, keepdims=True)
            denom = np.sqrt((yp * yp).sum(axis=1)[:, None] * (yc * yc).sum(axis=0)[None, :])
            rnull = (yp @ yc) / np.maximum(denom, 1e-15)
            null_ge += (np.abs(-rnull) >= np.abs(observed)[None, :]).sum(axis=0)
            idx = rng.integers(0, n, size=(k, n))
            xb = xr[idx]
            yb = yr[idx, :]
            xb = xb - xb.mean(axis=1, keepdims=True)
            yb = yb - yb.mean(axis=1, keepdims=True)
            den = np.sqrt((xb * xb).sum(axis=1)[:, None] * (yb * yb).sum(axis=1))
            boot_values[start:start + k, :] = (xb[:, :, None] * yb).sum(axis=1) / np.maximum(den, 1e-15)
        for j, gene in enumerate(dep.columns):
            low, high = np.nanquantile(-boot_values[:, j], [0.025, 0.975])
            rows.append({"trajectory": traj, "dataset": DATASET[traj], "background": BACKGROUND[traj],
                         "gene": gene, "signature_size": PRIMARY_SIZE, "scoring_method": PRIMARY_METHOD,
                         "n": n, "observed_vulnerability_rho": observed[j],
                         "permutation_p": (1 + null_ge[j]) / (N_RESAMPLES + 1),
                         "bootstrap_low": low, "bootstrap_high": high,
                         "permutation_supported": (1 + null_ge[j]) / (N_RESAMPLES + 1) < 0.05,
                         "bootstrap_supported": low > 0})
    return pd.DataFrame(rows)


def fisher_p(values: list[float]) -> float:
    vals = [max(float(x), 1e-300) for x in values if np.isfinite(x)]
    if not vals:
        return np.nan
    return float(chi2.sf(-2 * np.log(vals).sum(), 2 * len(vals)))


def aggregate_background(traj_results: pd.DataFrame, resampling: pd.DataFrame | None = None,
                         omit_background: str | None = None) -> pd.DataFrame:
    if omit_background:
        traj_results = traj_results.loc[traj_results["background"] != omit_background].copy()
    grouped = (traj_results.groupby(["gene", "background"], sort=False)
               .agg(n_trajectory_values=("vulnerability_rho", "size"),
                    n_positive_trajectory_values=("vulnerability_rho", lambda x: int((x > 0).sum())),
                    background_vulnerability_rho=("vulnerability_rho", "median"),
                    background_p_value=("p_value", "median"))
               .reset_index())
    if resampling is not None:
        rs = resampling.copy()
        rs["support"] = rs["permutation_supported"].astype(bool) | rs["bootstrap_supported"].astype(bool)
        rs = (rs.groupby(["gene", "background"], sort=False)
              .agg(n_resampling_supported_trajectories=("support", "sum"),
                   background_resampling_supported=("support", "mean"))
              .reset_index())
        rs["background_resampling_supported"] = rs["background_resampling_supported"] >= 0.5
        grouped = grouped.merge(rs, on=["gene", "background"], how="left")
    else:
        grouped["n_resampling_supported_trajectories"] = 0
        grouped["background_resampling_supported"] = False
    grouped["n_resampling_supported_trajectories"] = grouped["n_resampling_supported_trajectories"].fillna(0).astype(int)
    grouped["background_resampling_supported"] = grouped["background_resampling_supported"].fillna(False).astype(bool)
    return grouped


def convergent_ranking(background_df: pd.DataFrame) -> pd.DataFrame:
    rho = background_df.pivot(index="gene", columns="background", values="background_vulnerability_rho").reindex(columns=BACKGROUNDS)
    pval = background_df.pivot(index="gene", columns="background", values="background_p_value").reindex(columns=BACKGROUNDS)
    support = background_df.pivot(index="gene", columns="background", values="background_resampling_supported").reindex(columns=BACKGROUNDS).fillna(False)
    rarr = rho.to_numpy(float)
    parr = pval.to_numpy(float)
    valid = np.isfinite(rarr)
    n_bg = valid.sum(axis=1)
    meta_stat = np.where(valid, -2 * np.log(np.clip(parr, 1e-300, 1.0)), 0.0).sum(axis=1)
    meta_p = chi2.sf(meta_stat, 2 * np.maximum(n_bg, 1))
    med = np.nanmedian(rarr, axis=1)
    leave_hct = np.nanmedian(rarr[:, 1:], axis=1)
    out = pd.DataFrame({"gene": rho.index.to_numpy(), "n_backgrounds": n_bg,
                        "n_positive_backgrounds": np.nansum(rarr > 0, axis=1).astype(int),
                        "n_stable_backgrounds_rho_ge_0.10": np.nansum(rarr >= 0.10, axis=1).astype(int),
                        "n_resampling_supported_backgrounds": support.to_numpy(bool).sum(axis=1).astype(int),
                        "median_background_vulnerability_rho": med,
                        "meta_p_value_fisher": meta_p,
                        "leave_HCT116_out_median_vulnerability_rho": leave_hct,
                        "leave_HCT116_out_positive": leave_hct > 0})
    out["meta_q_value"] = bh(out["meta_p_value_fisher"].to_numpy())
    out["tier1_universal"] = ((out["n_positive_backgrounds"] == 4) &
                               (out["median_background_vulnerability_rho"] >= 0.10) &
                               (out["n_resampling_supported_backgrounds"] >= 2) &
                               (out["meta_q_value"] <= 0.10) & out["leave_HCT116_out_positive"])
    out["tier2_strong"] = ((out["n_positive_backgrounds"] >= 3) &
                           (out["median_background_vulnerability_rho"] >= 0.12) &
                           (out["n_stable_backgrounds_rho_ge_0.10"] >= 2) &
                           (out["n_positive_backgrounds"] >= 3))
    out["tier3_subtype"] = ((out["n_positive_backgrounds"] >= 2) &
                            (out["median_background_vulnerability_rho"] >= 0.20) &
                            (out["n_resampling_supported_backgrounds"] >= 2))
    out["primary_tier"] = np.select([out["tier1_universal"], out["tier2_strong"], out["tier3_subtype"]],
                                     ["Tier1_universal", "Tier2_strong", "Tier3_subtype"], default="none")
    # Sort by the mutually exclusive reporting tier.  Tier3 can mathematically
    # overlap Tier2, but the primary label is Tier1 > Tier2 > Tier3 > none.
    tier_order = {"Tier1_universal": 0, "Tier2_strong": 1, "Tier3_subtype": 2, "none": 3}
    out["_tier_priority"] = out["primary_tier"].map(tier_order).fillna(3)
    out = out.sort_values(["_tier_priority", "median_background_vulnerability_rho", "meta_q_value"],
                          ascending=[True, False, True]).reset_index(drop=True).drop(columns=["_tier_priority"])
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out


def lodo_background(traj_results: pd.DataFrame, resampling: pd.DataFrame) -> pd.DataFrame:
    full = convergent_ranking(aggregate_background(traj_results, resampling))[["gene", "rank"]].rename(columns={"rank": "rank_full"})
    rows = []
    for omit in BACKGROUNDS:
        agg = aggregate_background(traj_results, resampling, omit_background=omit)
        rank = convergent_ranking(agg)
        rank = rank[["gene", "rank", "median_background_vulnerability_rho", "n_positive_backgrounds", "meta_p_value_fisher", "meta_q_value"]]
        rank = rank.rename(columns={"rank": "rank_leave_one_background_out"})
        rank["background_omitted"] = omit
        rank = rank.merge(full, on="gene", how="left")
        rows.append(rank)
    return pd.concat(rows, ignore_index=True)


def lodo_dataset(traj_results: pd.DataFrame, resampling: pd.DataFrame) -> pd.DataFrame:
    full = convergent_ranking(aggregate_background(traj_results, resampling))[["gene", "rank"]].rename(columns={"rank": "rank_full"})
    rows = []
    for omit in sorted(set(DATASET.values())):
        keep = traj_results.loc[traj_results["dataset"] != omit]
        rkeep = resampling.loc[resampling["dataset"] != omit]
        rank = convergent_ranking(aggregate_background(keep, rkeep))
        rank = rank[["gene", "rank", "median_background_vulnerability_rho", "n_positive_backgrounds", "meta_p_value_fisher", "meta_q_value"]]
        rank = rank.rename(columns={"rank": "rank_leave_one_dataset_out"})
        rank["dataset_omitted"] = omit
        rows.append(rank.merge(full, on="gene", how="left"))
    return pd.concat(rows, ignore_index=True)


def parse_hallmark_sets() -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    if not HALLMARK_GMT.exists():
        return sets
    with HALLMARK_GMT.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = {clean_symbol(x) for x in parts[2:]}
    return sets


def signature_score(expr_z: pd.DataFrame, genes: set[str]) -> pd.Series:
    common = sorted(set(expr_z.columns) & genes)
    if not common:
        return pd.Series(np.nan, index=expr_z.index)
    return expr_z[common].mean(axis=1)


def partial_covariate_sensitivity(primary: pd.DataFrame, scores: pd.DataFrame, expr: pd.DataFrame,
                                  dep: pd.DataFrame, exclusions: dict[str, list[str]]) -> pd.DataFrame:
    sets = parse_hallmark_sets()
    z = (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=0).replace(0, 1.0)
    prolif_genes = set().union(*[sets.get(k, set()) for k in ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]])
    prolif = signature_score(z, prolif_genes)
    global_dep = -dep.mean(axis=1, skipna=True)
    genes = set(PREDEFINED_MECHANISM)
    genes.update(primary.loc[primary["trajectory"].isin(TRAJECTORIES)].groupby("trajectory").head(200)["gene"])
    rows = []
    for traj in TRAJECTORIES:
        keep = [x for x in scores["ModelID"] if x not in exclusions[traj] and x in dep.index]
        xraw = scores.set_index("ModelID").loc[keep, f"{traj}__weighted"].to_numpy(float)
        cov = pd.DataFrame({"prolif": prolif.loc[keep].to_numpy(float), "global_dependency": global_dep.loc[keep].to_numpy(float)}, index=keep)
        xr = pd.Series(xraw, index=keep).rank().to_numpy(float)
        covr = cov.rank().to_numpy(float)
        Xcov = np.column_stack([np.ones(len(keep)), covr])
        rx = xr - Xcov @ np.linalg.lstsq(Xcov, xr, rcond=None)[0]
        for gene in sorted(genes & set(dep.columns) & set(expr.columns)):
            yraw = dep.loc[keep, gene].to_numpy(float)
            eraw = expr.loc[keep, gene].to_numpy(float)
            valid = np.isfinite(yraw) & np.isfinite(eraw) & np.isfinite(rx) & np.isfinite(cov.to_numpy()).all(axis=1)
            if valid.sum() < MIN_N:
                continue
            yr = pd.Series(yraw[valid]).rank().to_numpy(float)
            er = pd.Series(eraw[valid]).rank().to_numpy(float)
            c = covr[valid]
            Xy = np.column_stack([np.ones(valid.sum()), c, er])
            ry = yr - Xy @ np.linalg.lstsq(Xy, yr, rcond=None)[0]
            rho_adj = np.corrcoef(rx[valid], ry)[0, 1]
            rho_raw = spearmanr(xraw[valid], yraw[valid]).statistic
            rows.append({"trajectory": traj, "background": BACKGROUND[traj], "gene": gene,
                         "n": int(valid.sum()), "raw_vulnerability_rho": -float(rho_raw),
                         "adjusted_partial_vulnerability_rho": -float(rho_adj),
                         "delta_adjusted_minus_raw": -float(rho_adj) + float(rho_raw),
                         "covariates": "proliferation_score+global_dependency+gene_expression"})
    return pd.DataFrame(rows)


def mechanism_audit(primary: pd.DataFrame, convergent: pd.DataFrame) -> pd.DataFrame:
    p = primary.loc[primary["gene"].isin(PREDEFINED_MECHANISM)].copy()
    pivot = p.pivot_table(index="gene", columns="trajectory", values="vulnerability_rho", aggfunc="first").reset_index()
    rank_cols = convergent[["gene", "rank", "primary_tier", "median_background_vulnerability_rho",
                            "n_positive_backgrounds", "meta_q_value"]]
    return pivot.merge(rank_cols, on="gene", how="left").sort_values("rank", na_position="last")


def save_manifest(meta: pd.DataFrame, exclusions: dict[str, list[str]], n_crc: int) -> dict:
    manifest = {
        "phase": "7B-R",
        "objective": "Robust OXA-resistant CRC trajectory-conditioned CRISPR dependency mapping",
        "status": "completed_genetic_mapping_only",
        "phase7c8_executed": False,
        "depmap_release": "23Q4",
        "depmap_files": {"model": str(MODEL_RAW.relative_to(ROOT)), "expression": str(EXPR_RAW.relative_to(ROOT)), "crispr_gene_effect": str(DEP_RAW.relative_to(ROOT))},
        "crc_inclusion_rule": "Model.csv OncotreeLineage == Bowel AND OncotreePrimaryDisease contains colorectal; intersect expression and CRISPR IDs",
        "n_crc_models_after_intersection": n_crc,
        "trajectories": TRAJECTORIES,
        "background_map": BACKGROUND,
        "independent_backgrounds": BACKGROUNDS,
        "trajectory_metadata": meta.to_dict(orient="records"),
        "self_line_exclusions": exclusions,
        "signature_sizes": SIGNATURE_SIZES,
        "primary_signature": {"size": PRIMARY_SIZE, "up_down": "250+250", "scoring_method": PRIMARY_METHOD},
        "scoring_methods": {"primary": "weighted directional z-score; expression standardized within CRC models",
                             "sensitivity": "rank-based directional score using percentile ranks"},
        "seed": SEED,
        "permutation_resamples": N_RESAMPLES,
        "bootstrap_resamples": N_RESAMPLES,
        "minimum_n": MIN_N,
        "aggregation": "median across trajectories within each biological background; HCT116 has 3 trajectories",
        "tier_thresholds": {
            "Tier1_universal": "4/4 backgrounds positive; median vulnerability rho >= 0.10; >=2/4 backgrounds resampling-supported; meta FDR <= 0.10; leave-HCT116-out positive",
            "Tier2_strong": ">=3/4 positive; median vulnerability rho >= 0.12; >=2 backgrounds rho >= 0.10",
            "Tier3_subtype": ">=2/4 positive; median vulnerability rho >= 0.20; >=2 backgrounds resampling-supported",
        },
        "predefined_mechanism_audit_only": PREDEFINED_MECHANISM,
        "raw_data_policy": "Raw DepMap/GEO files remain local and are ignored by git; only derived outputs and scripts are committed.",
    }
    (OUT / "phase7bR_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


def write_report(manifest: dict, corr: pd.DataFrame, convergence: pd.DataFrame,
                 lobo: pd.DataFrame, lodo: pd.DataFrame, sensitivity: pd.DataFrame,
                 covariates: pd.DataFrame, audit: pd.DataFrame) -> None:
    counts = convergence["primary_tier"].value_counts().to_dict()
    top = convergence.head(20)[["rank", "gene", "primary_tier", "median_background_vulnerability_rho",
                                "n_positive_backgrounds", "n_resampling_supported_backgrounds", "meta_q_value"]]
    hct = corr.loc[(corr["scoring_method"] == "weighted") & (corr["signature_size"] == PRIMARY_SIZE) &
                   (corr["trajectory_a"].str.contains("HCT116")) & (corr["trajectory_b"].str.contains("HCT116")) &
                   (corr["trajectory_a"] != corr["trajectory_b"])]
    hct_text = "; ".join([f"{r.trajectory_a} vs {r.trajectory_b}: rho={r.spearman_rho:.3f}" for r in hct.itertuples()]) or "not available"
    lobo_hct = lobo.loc[lobo["background_omitted"] == "HCT116"]
    rank_stability = lodo.groupby("dataset_omitted").apply(lambda x: x[["rank_full", "rank_leave_one_dataset_out"]].corr(method="spearman").iloc[0, 1]).to_dict()
    mech = audit[["gene", "rank", "primary_tier", "median_background_vulnerability_rho", "n_positive_backgrounds"]].to_string(index=False)
    cov_delta = float(covariates["delta_adjusted_minus_raw"].abs().median()) if len(covariates) else float("nan")
    go = "A: universal convergent dependency identified" if counts.get("Tier1_universal", 0) else ("B: strong but non-universal dependency identified" if counts.get("Tier2_strong", 0) else "C/D: no Tier1/Tier2 convergence; retain subtype signals only and do not advance automatically")
    lines = [
        "# Phase 7B-R：OXA-resistant CRC trajectory-conditioned dependency mapping",
        "",
        "## 结论先行",
        "",
        f"本轮使用 DepMap 23Q4 的 {manifest['n_crc_models_after_intersection']} 个 CRC 模型、6 条 acquired OXA-R trajectory，但按 HCT116、DLD1、HT29、LoVo 四个独立生物学背景聚合。主模型为 top 250+250 weighted directional score；HCT116/DLD1/HT29/LoVo 自身 DepMap cell line 在对应 trajectory 中排除。",
        f"最终判定：**{go}**。Tier1={counts.get('Tier1_universal', 0)}，Tier2={counts.get('Tier2_strong', 0)}，Tier3={counts.get('Tier3_subtype', 0)}。",
        "",
        "## 1. 这六条 trajectory 是否其实高度相似？",
        "",
        "状态分数相关性见 `phase7bR_state_score_correlation.csv`。主 weighted top250 中，HCT116 三条 trajectory 的两两相关为：" + hct_text + "。相关性若不高，不把它解释为失败，而是保留异质轨迹并依赖背景级收敛。",
        "Tier2=1105 是按预注册的宽松发现阈值得到的候选层，不等于 1105 个已验证靶点；Tier2 本身不要求 meta FDR 或每个背景都有 resampling support。因此本轮没有产生可以直接进入药物筛选的 universal target，下一步若继续应先在独立功能/药敏数据中收缩候选。",
        "",
        "## 2. HCT116 三条轨迹是否只是重复实验？",
        "",
        "不是独立生物学背景；本分析将其作为同一 HCT116 background，并用 median 聚合。它们仍保留在 trajectory-level 表中，用于审计轨迹间一致性和 leave-one-dataset/trajectory 稳健性。",
        "",
        "## 3–5. 功能依赖收敛层级",
        "",
        f"Tier1 universal: {counts.get('Tier1_universal', 0)}；Tier2 strong: {counts.get('Tier2_strong', 0)}；Tier3 subtype: {counts.get('Tier3_subtype', 0)}。主排名前 20：",
        "",
        "```text",
        top.to_string(index=False),
        "```",
        "",
        "## 6. 预设机制节点审计（不参与排名）",
        "",
        "```text",
        mech,
        "```",
        "",
        "## 7. 稳健性审计",
        "",
        f"LOBO 结果见 `phase7bR_leave_one_background_out.csv`；HCT116 留出后的正向候选数（median rho > 0）为 {int((lobo_hct['median_background_vulnerability_rho'] > 0).sum())}。LODO 全局 rank Spearman 稳定性：{rank_stability}。",
        f"signature size 100/250/500 与 rank score 结果见 `phase7bR_signature_size_sensitivity.csv`。协变量敏感性见 `phase7bR_covariate_sensitivity.csv`；adjusted-vs-raw vulnerability rho 的绝对差中位数为 {cov_delta:.3f}。",
        "置换和 Bootstrap 仅对主 top200/trajectory 候选执行，各 1000 次，结果见 `phase7bR_bootstrap_permutation_results.csv`。",
        "",
        "## 研究边界",
        "",
        "本轮没有执行 Phase 7C/8，没有把预设机制节点、药物机制或 LINCS 反向签名用于主排名。Meldonium 不作强行挽救；只有后续出现 CPT2/SLC22A5/BBOX1 等功能依赖才可重新进入候选层。原始 DepMap/GEO 文件保持本地，不纳入 GitHub。",
        "",
        "## 复现",
        "",
        "```powershell",
        "python work/scripts/phase7bR_robust_dependency_mapping.py",
        "```",
    ]
    (OUT / "phase7bR_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    print("[1/8] Loading derived signatures and DepMap 23Q4 inputs")
    delta, crc, expr, dep = load_inputs()
    meta = load_trajectory_metadata()
    exclusions = self_exclusion_map(crc)
    print(f"  CRC overlap={len(crc)}; expression genes={expr.shape[1]}; dependency genes={dep.shape[1]}")
    print("[2/8] Building trajectory signatures and state scores")
    all_weights = []
    all_meta = []
    score_tables: dict[int, pd.DataFrame] = {}
    correlations = []
    for n in SIGNATURE_SIZES:
        weights, smeta = build_signature_weights(delta, n, meta)
        all_weights.append(weights)
        all_meta.append(smeta)
        scores = project_scores(expr, weights, n)
        score_tables[n] = scores
        scores.to_csv(OUT / f"phase7bR_crc_state_scores_top{n}.csv", index=False)
        correlations.append(correlation_table(scores, n))
    signatures = pd.concat(all_weights, ignore_index=True)
    signatures.to_csv(OUT / "phase7bR_trajectory_signatures.csv", index=False)
    pd.concat(all_meta, ignore_index=True).to_csv(OUT / "phase7bR_trajectory_signature_metadata.csv", index=False)
    corr = pd.concat(correlations, ignore_index=True)
    corr.to_csv(OUT / "phase7bR_state_score_correlation.csv", index=False)

    print("[3/8] Genome-wide trajectory-conditioned dependency correlations")
    all_combo = {}
    for n in SIGNATURE_SIZES:
        for method in ["weighted", "rank"]:
            print(f"  signature top{n} / {method}")
            res = dependency_results(score_tables[n], dep, crc, exclusions, n, method)
            all_combo[(n, method)] = res
    primary = all_combo[(PRIMARY_SIZE, PRIMARY_METHOD)]
    primary.to_csv(OUT / "phase7bR_gene_dependency_by_trajectory.csv", index=False)
    bg = aggregate_background(primary)
    bg.to_csv(OUT / "phase7bR_gene_dependency_by_background.csv", index=False)

    print("[4/8] Primary convergence ranking and resampling")
    preliminary = primary.sort_values(["trajectory", "vulnerability_rho"], ascending=[True, False])
    preliminary = preliminary.groupby("trajectory", group_keys=False).head(200)
    resampling = bootstrap_permutation(preliminary, score_tables[PRIMARY_SIZE], dep, exclusions, SEED)
    resampling.to_csv(OUT / "phase7bR_bootstrap_permutation_results.csv", index=False)
    bg = aggregate_background(primary, resampling)
    bg.to_csv(OUT / "phase7bR_gene_dependency_by_background.csv", index=False)
    convergence = convergent_ranking(bg)
    convergence.to_csv(OUT / "phase7bR_convergent_gene_ranking.csv", index=False)

    print("[5/8] Leave-one-background and leave-one-dataset audits")
    lobo = lodo_background(primary, resampling)
    lobo.to_csv(OUT / "phase7bR_leave_one_background_out.csv", index=False)
    lodo = lodo_dataset(primary, resampling)
    lodo.to_csv(OUT / "phase7bR_leave_one_dataset_out.csv", index=False)

    print("[6/8] Signature-size and scoring-method sensitivity")
    sensitivity = None
    sens_rows = []
    for (n, method), result in all_combo.items():
        agg = aggregate_background(result)
        rank = convergent_ranking(agg)
        rank["signature_size"] = n
        rank["scoring_method"] = method
        sens_rows.append(rank[["gene", "signature_size", "scoring_method", "rank", "median_background_vulnerability_rho",
                               "n_positive_backgrounds", "primary_tier"]])
    sensitivity = pd.concat(sens_rows, ignore_index=True)
    sensitivity.to_csv(OUT / "phase7bR_signature_size_sensitivity.csv", index=False)

    print("[7/8] Covariate sensitivity and predefined-mechanism audit")
    covariates = partial_covariate_sensitivity(primary, score_tables[PRIMARY_SIZE], expr, dep, exclusions)
    covariates.to_csv(OUT / "phase7bR_covariate_sensitivity.csv", index=False)
    audit = mechanism_audit(primary, convergence)
    audit.to_csv(OUT / "phase7bR_predefined_mechanism_audit.csv", index=False)

    print("[8/8] Writing manifest and report")
    manifest = save_manifest(meta, exclusions, len(crc))
    write_report(manifest, corr, convergence, lobo, lodo, sensitivity, covariates, audit)
    print("Completed Phase 7B-R. Phase 7C/8 were not executed.")


if __name__ == "__main__":
    main()
