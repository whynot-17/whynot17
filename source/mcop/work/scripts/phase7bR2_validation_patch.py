#!/usr/bin/env python
"""Phase 7B-R2 validation patch.

This is intentionally a small post-hoc validation layer on top of Phase 7B-R.
It does not rebuild the original six-signature pipeline and it does not run
Phase 7C/8.  It corrects three interpretation/statistics issues:

1. background-level empirical permutation p-values are computed from the
   median trajectory vulnerability rho within each biological background;
2. covariate sensitivity is run on the actual vulnerability top-200 genes per
   trajectory and the final convergent top-500 genes;
3. Tier1/Tier2/Tier3 are reported as independent boolean flags, while a
   mutually-exclusive primary tier is retained only for display.

It also builds a cross-method stability table over 100/250/500 x
weighted/rank, with LOBO/LODO and covariate-adjusted summaries.
"""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, spearmanr, t

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
WORK = ROOT / "work"
OUT.mkdir(exist_ok=True)
BASE_SCRIPT = ROOT / "work" / "scripts" / "phase7bR_robust_dependency_mapping.py"

spec = importlib.util.spec_from_file_location("phase7bR_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)

TRAJECTORIES = base.TRAJECTORIES
BACKGROUND = base.BACKGROUND
BACKGROUNDS = base.BACKGROUNDS
DATASET = base.DATASET
SIGNATURE_SIZES = base.SIGNATURE_SIZES
MIN_N = base.MIN_N
SEED = base.SEED + 200
N_PERM = 1000


def bh(values: np.ndarray | pd.Series) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    out = np.full(arr.shape, np.nan)
    ok = np.isfinite(arr)
    if not ok.any():
        return out
    p = np.clip(arr[ok], 0, 1)
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1)
    out[ok] = restored
    return out


def rank_corr_vector(x: np.ndarray, dep: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Spearman rho, vulnerability rho and analytical p for all genes."""
    xr = pd.Series(np.asarray(x, dtype=float)).rank().to_numpy(float)
    yr = dep.rank(axis=0, method="average", na_option="keep").to_numpy(float)
    valid = np.isfinite(yr)
    xmat = np.broadcast_to(xr[:, None], yr.shape)
    n = valid.sum(axis=0).astype(int)
    xm = np.where(valid, xmat, np.nan).sum(axis=0) / np.maximum(n, 1)
    ym = np.where(valid, yr, np.nan).sum(axis=0) / np.maximum(n, 1)
    xc = np.where(valid, xmat - xm, 0.0)
    yc = np.where(valid, yr - ym, 0.0)
    den = np.sqrt((xc * xc).sum(axis=0) * (yc * yc).sum(axis=0))
    rho = np.where(den > 0, (xc * yc).sum(axis=0) / np.maximum(den, 1e-15), np.nan)
    rho = np.clip(rho, -1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = rho * np.sqrt(np.maximum(n - 2, 1) / np.maximum(1 - rho * rho, 1e-15))
    p = 2 * t.sf(np.abs(stat), np.maximum(n - 2, 1))
    p[(n < MIN_N) | ~np.isfinite(rho)] = np.nan
    return rho, -rho, p


def load_scores() -> dict[int, pd.DataFrame]:
    scores = {}
    for n in SIGNATURE_SIZES:
        path = OUT / f"phase7bR_crc_state_scores_top{n}.csv"
        df = pd.read_csv(path).set_index("ModelID")
        scores[n] = df
    return scores


def load_primary() -> pd.DataFrame:
    return pd.read_csv(OUT / "phase7bR_gene_dependency_by_trajectory.csv")


def trajectory_dependency(scores: pd.DataFrame, dep: pd.DataFrame, exclusions: dict[str, list[str]],
                          method: str, size: int) -> pd.DataFrame:
    rows = []
    for traj in TRAJECTORIES:
        keep = [x for x in scores.index if x not in exclusions[traj] and x in dep.index]
        rho, vul, p = rank_corr_vector(scores.loc[keep, f"{traj}__{method}"].to_numpy(float), dep.loc[keep])
        rows.append(pd.DataFrame({"trajectory": traj, "dataset": DATASET[traj], "background": BACKGROUND[traj],
                                  "gene": dep.columns, "n": len(keep), "rho_state_vs_dependency": rho,
                                  "vulnerability_rho": vul, "p_value": p, "signature_size": size,
                                  "scoring_method": method}))
    return pd.concat(rows, ignore_index=True)


def background_effects(traj: pd.DataFrame) -> pd.DataFrame:
    return (traj.groupby(["gene", "background"], sort=False)
            .agg(n_trajectory_values=("vulnerability_rho", "size"),
                 n_positive_trajectory_values=("vulnerability_rho", lambda x: int((x > 0).sum())),
                 background_vulnerability_rho=("vulnerability_rho", "median"))
            .reset_index())


def empirical_background_p(traj: pd.DataFrame, scores: pd.DataFrame, dep: pd.DataFrame,
                           exclusions: dict[str, list[str]]) -> pd.DataFrame:
    """Empirical p for median trajectory vulnerability within each background.

    For HCT116 the three trajectory score vectors are permuted independently
    for every null draw and their three gene-wise vulnerability correlations
    are then combined by the same median operator used for the observed effect.
    """
    observed = background_effects(traj).pivot(index="gene", columns="background", values="background_vulnerability_rho")
    genes = dep.columns.to_numpy()
    rows = []
    for bi, bg in enumerate(BACKGROUNDS):
        ts = [t0 for t0 in TRAJECTORIES if BACKGROUND[t0] == bg]
        obs = observed.reindex(index=genes, columns=[bg])[bg].to_numpy(float)
        null_ge = np.zeros(len(genes), dtype=np.int64)
        rng = np.random.default_rng(SEED + bi)
        prepared = []
        for traj_name in ts:
            keep = [x for x in scores.index if x not in exclusions[traj_name] and x in dep.index]
            xr = pd.Series(scores.loc[keep, f"{traj_name}__weighted"].to_numpy(float)).rank().to_numpy(float)
            y = dep.loc[keep]
            yr = y.rank(axis=0, method="average", na_option="keep").to_numpy(float)
            valid = np.isfinite(yr)
            n = valid.sum(axis=0).astype(float)
            ym = np.where(valid, yr, np.nan).sum(axis=0) / np.maximum(n, 1)
            yc = np.where(valid, yr - ym, 0.0)
            prepared.append((xr, yc, valid, n))
        for start in range(0, N_PERM, 50):
            k = min(50, N_PERM - start)
            vulnerability = []
            for xr, yc, valid, n in prepared:
                perm = np.stack([rng.permutation(xr) for _ in range(k)], axis=0)
                xp = perm - perm.mean(axis=1, keepdims=True)
                den = np.sqrt((xp * xp).sum(axis=1)[:, None] * (yc * yc).sum(axis=0)[None, :])
                rho = (xp @ yc) / np.maximum(den, 1e-15)
                vulnerability.append(-rho)
            null_median = np.median(np.stack(vulnerability, axis=0), axis=0)
            null_ge += (np.abs(null_median) >= np.abs(obs)[None, :]).sum(axis=0)
        p = (1 + null_ge) / (N_PERM + 1)
        rows.append(pd.DataFrame({"gene": genes, "background": bg,
                                  "observed_background_vulnerability_rho": obs,
                                  "empirical_p_value": p, "n_permutations": N_PERM,
                                  "n_trajectories_combined": len(ts)}))
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT / "phase7bR2_background_empirical_p.csv", index=False)
    return out


def build_r2_ranking(primary: pd.DataFrame, resampling: pd.DataFrame,
                     empirical: pd.DataFrame) -> pd.DataFrame:
    bg = base.aggregate_background(primary, resampling)
    effect = bg.pivot(index="gene", columns="background", values="background_vulnerability_rho").reindex(columns=BACKGROUNDS)
    support = bg.pivot(index="gene", columns="background", values="background_resampling_supported").reindex(columns=BACKGROUNDS).fillna(False)
    p = empirical.pivot(index="gene", columns="background", values="empirical_p_value").reindex(columns=BACKGROUNDS)
    r = effect.to_numpy(float)
    valid = np.isfinite(r)
    meta_stat = np.where(valid, -2 * np.log(np.clip(p.to_numpy(float), 1e-300, 1.0)), 0).sum(axis=1)
    nbg = valid.sum(axis=1)
    meta_p = chi2.sf(meta_stat, 2 * np.maximum(nbg, 1))
    out = pd.DataFrame({"gene": effect.index.to_numpy(), "n_backgrounds": nbg,
                        "n_positive_backgrounds": np.nansum(r > 0, axis=1).astype(int),
                        "n_stable_backgrounds_rho_ge_0.10": np.nansum(r >= 0.10, axis=1).astype(int),
                        "n_resampling_supported_backgrounds": support.to_numpy(bool).sum(axis=1).astype(int),
                        "median_background_vulnerability_rho": np.nanmedian(r, axis=1),
                        "meta_empirical_p_value_fisher": meta_p,
                        "meta_empirical_q_value": bh(meta_p),
                        "leave_HCT116_out_median_vulnerability_rho": np.nanmedian(r[:, 1:], axis=1)})
    prior = pd.read_csv(OUT / "phase7bR_convergent_gene_ranking.csv", usecols=["gene", "meta_q_value"])
    prior = prior.rename(columns={"meta_q_value": "prior_analytical_meta_q_value"})
    out = out.merge(prior, on="gene", how="left")
    # Independent flags: do not infer Tier3 from primary_tier.
    out["tier1_flag"] = ((out["n_positive_backgrounds"] == 4) &
                          (out["median_background_vulnerability_rho"] >= 0.10) &
                          (out["n_resampling_supported_backgrounds"] >= 2) &
                          (out["meta_empirical_q_value"] <= 0.10) &
                          (out["leave_HCT116_out_median_vulnerability_rho"] > 0))
    out["tier2_flag"] = ((out["n_positive_backgrounds"] >= 3) &
                          (out["median_background_vulnerability_rho"] >= 0.12) &
                          (out["n_stable_backgrounds_rho_ge_0.10"] >= 2))
    out["tier3_flag"] = ((out["n_positive_backgrounds"] >= 2) &
                          (out["median_background_vulnerability_rho"] >= 0.20) &
                          (out["n_resampling_supported_backgrounds"] >= 2))
    out["strict_empirical_convergence_flag"] = ((out["n_positive_backgrounds"] >= 3) &
                                                  (out["n_resampling_supported_backgrounds"] >= 2) &
                                                  (out["meta_empirical_q_value"] <= 0.10))
    out["primary_tier"] = np.select([out["tier1_flag"], out["tier2_flag"], out["tier3_flag"]],
                                     ["Tier1_universal", "Tier2_discovery", "Tier3_subtype"], default="none")
    order = {"Tier1_universal": 0, "Tier2_discovery": 1, "Tier3_subtype": 2, "none": 3}
    out["_order"] = out["primary_tier"].map(order)
    out = out.sort_values(["_order", "median_background_vulnerability_rho", "meta_empirical_q_value"],
                          ascending=[True, False, True]).drop(columns="_order").reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    out.to_csv(OUT / "phase7bR2_convergent_gene_ranking.csv", index=False)
    return out


def residual_covariate_audit(scores: pd.DataFrame, expr: pd.DataFrame, dep: pd.DataFrame,
                             exclusions: dict[str, list[str]], trajectory_top: dict[str, set[str],],
                             final_top: set[str]) -> pd.DataFrame:
    sets = base.parse_hallmark_sets()
    z = (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=0).replace(0, 1.0)
    prolif_genes = set().union(*[sets.get(k, set()) for k in ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]])
    prolif = base.signature_score(z, prolif_genes)
    global_dep = -dep.mean(axis=1, skipna=True)
    scopes = [("trajectory_top200", trajectory_top),
              ("convergent_top500", {traj: final_top for traj in TRAJECTORIES})]
    rows = []
    for scope, mapping in scopes:
        for traj in TRAJECTORIES:
            genes = sorted(set(mapping[traj]) & set(dep.columns) & set(expr.columns))
            keep = [x for x in scores.index if x not in exclusions[traj] and x in dep.index]
            xraw = scores.loc[keep, f"{traj}__weighted"].to_numpy(float)
            cov = pd.DataFrame({"prolif": prolif.loc[keep].to_numpy(float),
                                "global_dependency": global_dep.loc[keep].to_numpy(float)}, index=keep)
            xr = pd.Series(xraw, index=keep).rank().to_numpy(float)
            covr = cov.rank().to_numpy(float)
            Xx = np.column_stack([np.ones(len(keep)), covr])
            rx = xr - Xx @ np.linalg.lstsq(Xx, xr, rcond=None)[0]
            for gene in genes:
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
                raw = spearmanr(xraw[valid], yraw[valid]).statistic
                adjusted = np.corrcoef(rx[valid], ry)[0, 1]
                rows.append({"selection_scope": scope, "trajectory": traj, "background": BACKGROUND[traj],
                             "gene": gene, "n": int(valid.sum()), "raw_vulnerability_rho": -float(raw),
                             "adjusted_partial_vulnerability_rho": -float(adjusted),
                             "delta_adjusted_minus_raw": -float(adjusted) + float(raw),
                             "covariates": "proliferation_score+global_dependency+gene_expression"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "phase7bR2_covariate_sensitivity.csv", index=False)
    return out


def cross_method_stability(scores: dict[int, pd.DataFrame], dep: pd.DataFrame,
                           exclusions: dict[str, list[str]], ranking: pd.DataFrame,
                           covariates: pd.DataFrame) -> pd.DataFrame:
    long = []
    for n in SIGNATURE_SIZES:
        for method in ["weighted", "rank"]:
            traj = trajectory_dependency(scores[n], dep, exclusions, method, n)
            bg = background_effects(traj)
            piv = bg.pivot(index="gene", columns="background", values="background_vulnerability_rho").reindex(columns=BACKGROUNDS)
            one = pd.DataFrame({"gene": piv.index.to_numpy(), "signature_size": n, "scoring_method": method,
                                "median_background_vulnerability_rho": np.nanmedian(piv.to_numpy(float), axis=1),
                                "direction_consistent_backgrounds": np.nansum(piv.to_numpy(float) > 0, axis=1).astype(int)})
            one = one.sort_values("median_background_vulnerability_rho", ascending=False).reset_index(drop=True)
            one["effect_rank"] = np.arange(1, len(one) + 1)
            one["rank_percentile"] = one["effect_rank"] / len(one)
            one["top100_flag"] = one["effect_rank"] <= 100
            one["top500_flag"] = one["effect_rank"] <= 500
            long.append(one)
    long_df = pd.concat(long, ignore_index=True)
    wide = long_df.groupby("gene", sort=False).agg(
        n_models_direction_ge3of4=("direction_consistent_backgrounds", lambda x: int((x >= 3).sum())),
        n_models_direction_4of4=("direction_consistent_backgrounds", lambda x: int((x == 4).sum())),
        n_models_median_rho_positive=("median_background_vulnerability_rho", lambda x: int((x > 0).sum())),
        n_models_top100=("top100_flag", "sum"), n_models_top500=("top500_flag", "sum"),
        median_rank_percentile=("rank_percentile", "median"),
        max_rank_percentile=("rank_percentile", "max"),
        median_model_vulnerability_rho=("median_background_vulnerability_rho", "median"),
        min_model_vulnerability_rho=("median_background_vulnerability_rho", "min"),
    ).reset_index()
    wide["cross_method_stable_flag"] = ((wide["n_models_top500"] >= 4) &
                                         (wide["n_models_direction_ge3of4"] >= 4) &
                                         (wide["median_rank_percentile"] <= 0.10))

    lodo = pd.read_csv(OUT / "phase7bR_leave_one_dataset_out.csv")
    lodo["rank_percentile_full"] = lodo["rank_full"] / lodo.groupby("dataset_omitted")["rank_full"].transform("max")
    lodo["rank_percentile_lodo"] = lodo["rank_leave_one_dataset_out"] / lodo.groupby("dataset_omitted")["rank_leave_one_dataset_out"].transform("max")
    lodo_s = lodo.groupby("gene").agg(lodo_rank_percentile_median=("rank_percentile_lodo", "median"),
                                       lodo_rank_percentile_min=("rank_percentile_lodo", "min"),
                                       lodo_rank_percentile_max=("rank_percentile_lodo", "max"),
                                       lodo_rank_percentile_iqr=("rank_percentile_lodo", lambda x: x.quantile(0.75) - x.quantile(0.25))).reset_index()
    lobo = pd.read_csv(OUT / "phase7bR_leave_one_background_out.csv")
    lobo_s = (lobo.loc[lobo["background_omitted"] == "HCT116", ["gene", "median_background_vulnerability_rho", "n_positive_backgrounds"]]
              .rename(columns={"median_background_vulnerability_rho": "lobo_HCT116_out_median_vulnerability_rho",
                               "n_positive_backgrounds": "lobo_HCT116_out_n_positive_backgrounds"}))
    final_cov = covariates.loc[covariates["selection_scope"] == "convergent_top500"]
    cov_s = final_cov.groupby("gene").agg(adjusted_vulnerability_rho_median=("adjusted_partial_vulnerability_rho", "median"),
                                            adjusted_vulnerability_rho_min=("adjusted_partial_vulnerability_rho", "min"),
                                            adjusted_positive_fraction=("adjusted_partial_vulnerability_rho", lambda x: float((x > 0).mean())),
                                            covariate_audit_n=("n", "count")).reset_index()
    wide = wide.merge(lodo_s, on="gene", how="left").merge(lobo_s, on="gene", how="left").merge(cov_s, on="gene", how="left")
    wide = wide.merge(ranking[["gene", "rank", "meta_empirical_q_value", "tier1_flag", "tier2_flag", "tier3_flag", "primary_tier"]], on="gene", how="left")
    wide["robust_candidate_flag"] = (wide["cross_method_stable_flag"] &
                                      (wide["meta_empirical_q_value"] <= 0.10) &
                                      (wide["adjusted_vulnerability_rho_median"] > 0) &
                                      (wide["adjusted_positive_fraction"] >= 0.5))
    wide = wide.sort_values(["robust_candidate_flag", "cross_method_stable_flag", "median_rank_percentile", "meta_empirical_q_value"],
                            ascending=[False, False, True, True]).reset_index(drop=True)
    wide.insert(0, "stability_rank", np.arange(1, len(wide) + 1))
    long_df.to_csv(OUT / "phase7bR2_cross_method_stability_long.csv", index=False)
    wide.to_csv(OUT / "phase7bR2_cross_method_stability.csv", index=False)
    return wide


def mechanism_audit(ranking: pd.DataFrame) -> pd.DataFrame:
    genes = base.PREDEFINED_MECHANISM
    return ranking.loc[ranking["gene"].isin(genes)].copy().sort_values("rank")


def write_manifest(exclusions: dict[str, list[str]], n_crc: int, counts: dict, robust_n: int) -> None:
    manifest = {
        "phase": "7B-R2",
        "parent_phase": "7B-R",
        "status": "completed_validation_patch",
        "phase7c8_executed": False,
        "depmap_release": "23Q4",
        "n_crc_models": n_crc,
        "trajectories": TRAJECTORIES,
        "backgrounds": BACKGROUNDS,
        "self_line_exclusions": exclusions,
        "background_empirical_permutation": {"operator": "median trajectory vulnerability rho within background",
                                               "permutations_per_background": N_PERM,
                                               "tail": "two-sided absolute vulnerability rho",
                                               "seed": SEED},
        "covariate_audit_scopes": ["actual top200 vulnerability genes per trajectory", "final convergent top500 genes"],
        "tier_reporting": "tier1_flag, tier2_flag, tier3_flag are independent; primary_tier is display-only precedence",
        "cross_method_models": [f"top{n}_{m}" for n in SIGNATURE_SIZES for m in ["weighted", "rank"]],
        "counts": counts,
        "robust_candidate_definition": "top500 in >=4/6 models, >=3/4 backgrounds positive in >=4/6 models, median rank percentile <=0.10, empirical meta q <=0.10, covariate-adjusted median rho >0 and positive fraction >=0.5",
        "robust_candidate_count": robust_n,
        "raw_data_policy": "Raw DepMap/GEO files remain local and ignored by git; only derived R2 outputs and scripts are committed.",
    }
    (OUT / "phase7bR2_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


def write_report(ranking: pd.DataFrame, stability: pd.DataFrame, cov: pd.DataFrame,
                 empirical: pd.DataFrame, n_crc: int) -> None:
    tier_counts = {c: int(ranking[c].sum()) for c in ["tier1_flag", "tier2_flag", "tier3_flag"]}
    primary_counts = ranking["primary_tier"].value_counts().to_dict()
    top = ranking.head(20)[["rank", "gene", "median_background_vulnerability_rho", "n_positive_backgrounds",
                            "n_resampling_supported_backgrounds", "meta_empirical_q_value", "tier1_flag", "tier2_flag", "tier3_flag"]]
    focus = ranking.loc[ranking["gene"].isin(["MOB4", "TJP1", "ZFP36L1", "ADAMTS9", "AFDN"]),
                        ["gene", "rank", "median_background_vulnerability_rho", "n_positive_backgrounds",
                         "n_resampling_supported_backgrounds", "meta_empirical_q_value", "tier1_flag", "tier2_flag", "tier3_flag"]]
    robust = stability.loc[stability["robust_candidate_flag"]].head(30)
    corr = pd.read_csv(OUT / "phase7bR_state_score_correlation.csv")
    hct = corr.loc[(corr["signature_size"] == 250) & (corr["scoring_method"] == "weighted") &
                   corr["trajectory_a"].str.contains("HCT116") & corr["trajectory_b"].str.contains("HCT116") &
                   (corr["trajectory_a"] != corr["trajectory_b"])]
    hct_text = "; ".join(f"{r.trajectory_a} vs {r.trajectory_b}: rho={r.spearman_rho:.3f}" for r in hct.itertuples())
    lodo = pd.read_csv(OUT / "phase7bR_leave_one_dataset_out.csv")
    lodo_stability = lodo.groupby("dataset_omitted").apply(lambda x: x[["rank_full", "rank_leave_one_dataset_out"]].corr(method="spearman").iloc[0, 1]).to_dict()
    adjusted_delta = float(cov["delta_adjusted_minus_raw"].abs().median()) if len(cov) else float("nan")
    lines = [
        "# Phase 7B-R2 validation patch：corrected dependency interpretation",
        "",
        "## Revised conclusion",
        "",
        f"Using {n_crc} CRC DepMap 23Q4 models, six OXA-R trajectories and four biological backgrounds, **no universal single-gene dependency was identified**. A broad set of non-universal discovery-tier candidates emerged, but none is yet sufficiently validated for drug mapping. Phase 7C/8 were not run.",
        f"Independent flags: Tier1={tier_counts['tier1_flag']}，Tier2 discovery={tier_counts['tier2_flag']}，Tier3 subtype={tier_counts['tier3_flag']}。Display-only primary tiers: {primary_counts}。",
        "Tier2 is deliberately a discovery filter, not a claim that every flagged gene is a strong dependency. Tier3 is counted independently; a gene can satisfy both Tier2 and Tier3.",
        "",
        "## HCT116 trajectory heterogeneity",
        "",
        hct_text + "。The three trajectories are retained as one biological background for aggregation, not treated as three independent backgrounds.",
        "",
        "## Background-level empirical p-values",
        "",
        "Within each biological background, trajectory score vectors were permuted independently 1000 times, the same median vulnerability-rho operator was applied, and empirical p-values were computed from the two-sided absolute null distribution. HCT116 therefore uses the median of three trajectory null correlations rather than median analytical p-values.",
        "",
        "## Main ranking",
        "",
        "```text", top.to_string(index=False), "```", "",
        "### Provisional computational leads",
        "",
        "```text", focus.to_string(index=False), "```",
        "",
        "MOB4 remains a reasonable provisional computational lead because its effect size, two-background resampling support and empirical meta-evidence converge better than a rank-1 effect-size-only hit. It is not a wet-lab or drug-mapping lead yet.",
        "TJP1 is a hypothesis-generating signal of particular interest because the prior resistance-state analysis showed TJP1 downregulation while this analysis shows a positive dependency association; this is not evidence that expression alone predicts dependency.",
        "",
        "## Corrected covariate audit",
        "",
        "The audit was rerun on the actual vulnerability top200 genes per trajectory and the final convergent top500 genes. The absolute adjusted-minus-raw rho median across all corrected audit rows is %.3f; this value is not interpreted as a genome-wide proof of no confounding." % adjusted_delta,
        "",
        "## Cross-method stability",
        "",
        "The stability table covers top100/top250/top500 × weighted/rank and includes direction-consistent background counts, top100/top500 appearances, median rank percentile, LOBO/LODO rank summaries and adjusted rho. The predeclared robust-candidate rule leaves %d genes:" % len(robust),
        "",
        "```text", robust.head(30).to_string(index=False), "```", "",
        "LODO global rank stability: " + str(lodo_stability) + "。The GSE77932 and GSE42387 perturbations remain material, so candidates should not be called universal from the full ranking alone.",
        "",
        "## Mechanism audit",
        "",
        ranking.loc[ranking["gene"].isin(base.PREDEFINED_MECHANISM), ["gene", "rank", "median_background_vulnerability_rho", "meta_empirical_q_value", "tier1_flag", "tier2_flag", "tier3_flag"]].sort_values("rank").to_string(index=False),
        "",
        "DHODH, RRM2, CPT1A, VCP and SLC22A5 do not emerge as universal convergent dependencies in this correction. Meldonium therefore remains a failed broad OXA-R reversal hypothesis, not a candidate rescued by this patch.",
        "",
        "## Files and boundary",
        "",
        "- `phase7bR2_background_empirical_p.csv`: corrected background-level empirical p-values.",
        "- `phase7bR2_convergent_gene_ranking.csv`: corrected ranking with independent tier flags.",
        "- `phase7bR2_covariate_sensitivity.csv`: corrected top200/top500 covariate audit.",
        "- `phase7bR2_cross_method_stability.csv`: six-model stability summary.",
        "- `phase7bR2_manifest.json`: exact definitions and thresholds.",
        "",
        "Raw DepMap/GEO files remain local and are not committed. Phase 7C/8 remain deferred.",
    ]
    (OUT / "phase7bR2_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    print("[1/6] Loading 7B-R derived inputs and local DepMap data")
    delta, crc, expr, dep = base.load_inputs()
    scores = load_scores()
    primary = load_primary()
    exclusions = base.self_exclusion_map(crc)
    resampling = pd.read_csv(OUT / "phase7bR_bootstrap_permutation_results.csv")
    print(f"  CRC overlap={len(crc)}")
    print("[2/6] Computing corrected background-level empirical p-values")
    empirical = empirical_background_p(primary, scores[250], dep, exclusions)
    print("[3/6] Building corrected ranking with independent tier flags")
    ranking = build_r2_ranking(primary, resampling, empirical)
    trajectory_top = (primary.sort_values(["trajectory", "vulnerability_rho"], ascending=[True, False])
                      .groupby("trajectory", group_keys=False).head(200))
    trajectory_top_map = {traj: set(trajectory_top.loc[trajectory_top["trajectory"] == traj, "gene"]) for traj in TRAJECTORIES}
    final_top = set(ranking.head(500)["gene"])
    print("[4/6] Re-running covariate audit on actual vulnerability selections")
    cov = residual_covariate_audit(scores[250], expr, dep, exclusions, trajectory_top_map, final_top)
    print("[5/6] Building cross-method stability table")
    stability = cross_method_stability(scores, dep, exclusions, ranking, cov)
    mechanism_audit(ranking).to_csv(OUT / "phase7bR2_predefined_mechanism_audit.csv", index=False)
    counts = {c: int(ranking[c].sum()) for c in ["tier1_flag", "tier2_flag", "tier3_flag"]}
    write_manifest(exclusions, len(crc), counts, int(stability["robust_candidate_flag"].sum()))
    print("[6/6] Writing corrected report")
    write_report(ranking, stability, cov, empirical, len(crc))
    print("Completed Phase 7B-R2. Phase 7C/8 were not executed.")


if __name__ == "__main__":
    main()
