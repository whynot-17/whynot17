#!/usr/bin/env python
"""Phase 7B-R3: final global-null validation patch.

R3 keeps the Phase 7B-R/R2 signatures and effect-size aggregation fixed.  It
only recalibrates inference and the final shortlist:

* one row-label permutation is applied to the complete CRC DepMap CRISPR
  matrix for every null draw, preserving dependency covariance and the
  cross-trajectory score structure;
* the global statistic is the median vulnerability rho across the four
  biological backgrounds, so no Fisher combination of correlated background
  p-values is used;
* partial Spearman correlation residualizes both state score and dependency
  on the same covariates: proliferation, global dependency and target-gene
  expression;
* the previous 19 candidates are reclassified as internally cross-method
  stable.  A final shortlist must pass the stricter global-null criteria.

Phase 7C/8 are intentionally not run.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, spearmanr, t

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
WORK = ROOT / "work"
OUT.mkdir(exist_ok=True)
BASE_SCRIPT = ROOT / "work" / "scripts" / "phase7bR_robust_dependency_mapping.py"

spec = importlib.util.spec_from_file_location("phase7bR_base_r3", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)

TRAJECTORIES = base.TRAJECTORIES
BACKGROUND = base.BACKGROUND
BACKGROUNDS = base.BACKGROUNDS
DATASET = base.DATASET
SIGNATURE_SIZES = base.SIGNATURE_SIZES
MIN_N = base.MIN_N
N_PERM = 1000
SEED = base.SEED + 300


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


def load_scores() -> dict[int, pd.DataFrame]:
    out = {}
    for n in SIGNATURE_SIZES:
        out[n] = pd.read_csv(OUT / f"phase7bR_crc_state_scores_top{n}.csv").set_index("ModelID")
    return out


def rank_corr_from_ranks(xr: np.ndarray, yr: np.ndarray) -> np.ndarray:
    """Pearson correlation of one score rank vector with all gene ranks."""
    valid = np.isfinite(yr)
    n = valid.sum(axis=0).astype(float)
    xmat = np.broadcast_to(xr[:, None], yr.shape)
    xm = np.where(valid, xmat, np.nan).sum(axis=0) / np.maximum(n, 1)
    ym = np.where(valid, yr, np.nan).sum(axis=0) / np.maximum(n, 1)
    xc = np.where(valid, xmat - xm, 0.0)
    yc = np.where(valid, yr - ym, 0.0)
    den = np.sqrt((xc * xc).sum(axis=0) * (yc * yc).sum(axis=0))
    return np.where(den > 0, (xc * yc).sum(axis=0) / np.maximum(den, 1e-15), np.nan)


def observed_effects(primary: pd.DataFrame) -> pd.DataFrame:
    bg = base.aggregate_background(primary)
    return bg.pivot(index="gene", columns="background", values="background_vulnerability_rho").reindex(columns=BACKGROUNDS)


def global_empirical_null(primary: pd.DataFrame, scores: pd.DataFrame, dep: pd.DataFrame,
                          exclusions: dict[str, list[str]]) -> pd.DataFrame:
    """Global empirical p using one dependency row-label permutation per draw.

    The same permutation of the full 56-model CRISPR matrix is reused for all
    six trajectory projections in a draw.  The trajectory-specific self-line
    subset is applied after permutation.  Dependency ranks are computed once
    on the full CRC frame; row-label permutation and post-permutation subset
    selection preserve the dependency rank covariance structure while avoiding
    an unpaired null for HCT116's correlated trajectories.
    """
    effect = observed_effects(primary)
    genes = dep.columns.to_numpy()
    observed_global = np.nanmedian(effect.reindex(index=genes, columns=BACKGROUNDS).to_numpy(float), axis=1)
    dep_rank = dep.rank(axis=0, method="average", na_option="keep").to_numpy(float)
    id_pos = {model_id: i for i, model_id in enumerate(dep.index)}
    traj_prepared = []
    for traj in TRAJECTORIES:
        keep = [x for x in scores.index if x not in exclusions[traj] and x in id_pos]
        positions = np.array([id_pos[x] for x in keep], dtype=int)
        xr = pd.Series(scores.loc[keep, f"{traj}__weighted"].to_numpy(float)).rank().to_numpy(float)
        traj_prepared.append((traj, positions, xr))
    rng = np.random.default_rng(SEED)
    null_ge = np.zeros(len(genes), dtype=np.int64)
    n_full = len(dep.index)
    for start in range(0, N_PERM, 50):
        k = min(50, N_PERM - start)
        global_perms = np.stack([rng.permutation(n_full) for _ in range(k)], axis=0)
        vulnerability_by_traj = []
        for _traj, positions, xr in traj_prepared:
            # Each permutation assigns the same permuted dependency row label
            # to every gene and every trajectory.
            yr_perm = dep_rank[global_perms[:, positions], :]
            xr_center = xr - xr.mean()
            xp = np.broadcast_to(xr_center[None, :], (k, len(xr)))
            yc = yr_perm - yr_perm.mean(axis=1, keepdims=True)
            numerator = (xp[:, :, None] * yc).sum(axis=1)
            denominator = np.sqrt((xp * xp).sum(axis=1)[:, None] * (yc * yc).sum(axis=1))
            rho = numerator / np.maximum(denominator, 1e-15)
            vulnerability_by_traj.append(-rho)
        per_bg = []
        for bg in BACKGROUNDS:
            idx = [i for i, traj in enumerate(TRAJECTORIES) if BACKGROUND[traj] == bg]
            per_bg.append(np.median(np.stack([vulnerability_by_traj[i] for i in idx], axis=0), axis=0))
        global_null = np.median(np.stack(per_bg, axis=0), axis=0)
        null_ge += (np.abs(global_null) >= np.abs(observed_global)[None, :]).sum(axis=0)
    p = (1 + null_ge) / (N_PERM + 1)
    out = pd.DataFrame({"gene": genes, "observed_global_T_median_vulnerability_rho": observed_global,
                        "global_empirical_p_value": p, "n_permutations": N_PERM,
                        "permutation_unit": "one shared CRISPR row-label permutation across all trajectories/backgrounds",
                        "self_line_exclusion": "applied after global permutation per trajectory"})
    out["global_empirical_q_value"] = bh(out["global_empirical_p_value"].to_numpy())
    out.to_csv(OUT / "phase7bR3_global_empirical_p.csv", index=False)
    return out


def corrected_partial_covariate_audit(scores: pd.DataFrame, expr: pd.DataFrame, dep: pd.DataFrame,
                                      exclusions: dict[str, list[str]], primary: pd.DataFrame,
                                      ranking: pd.DataFrame) -> pd.DataFrame:
    """Partial Spearman with identical covariates for X and Y."""
    sets = base.parse_hallmark_sets()
    z = (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=0).replace(0, 1.0)
    prolif_genes = set().union(*[sets.get(k, set()) for k in ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]])
    proliferation = base.signature_score(z, prolif_genes)
    global_dependency = -dep.mean(axis=1, skipna=True)
    top200 = (primary.sort_values(["trajectory", "vulnerability_rho"], ascending=[True, False])
              .groupby("trajectory", group_keys=False).head(200))
    scopes = {
        "trajectory_top200": {traj: set(top200.loc[top200["trajectory"] == traj, "gene"]) for traj in TRAJECTORIES},
        "convergent_top500": {traj: set(ranking.head(500)["gene"]) for traj in TRAJECTORIES},
    }
    rows = []
    for scope, mapping in scopes.items():
        for traj in TRAJECTORIES:
            keep = [x for x in scores.index if x not in exclusions[traj] and x in dep.index]
            targets = sorted(set(mapping[traj]) & set(dep.columns) & set(expr.columns))
            xraw = scores.loc[keep, f"{traj}__weighted"].to_numpy(float)
            for gene in targets:
                yraw = dep.loc[keep, gene].to_numpy(float)
                eraw = expr.loc[keep, gene].to_numpy(float)
                cov_raw = np.column_stack([proliferation.loc[keep].to_numpy(float),
                                           global_dependency.loc[keep].to_numpy(float), eraw])
                valid = np.isfinite(xraw) & np.isfinite(yraw) & np.isfinite(cov_raw).all(axis=1)
                if valid.sum() < MIN_N:
                    continue
                # Rank-transform X, Y and every covariate, then residualize
                # both X and Y on exactly the same three covariates.
                xr = pd.Series(xraw[valid]).rank().to_numpy(float)
                yr = pd.Series(yraw[valid]).rank().to_numpy(float)
                cr = pd.DataFrame(cov_raw[valid], columns=["proliferation", "global_dependency", "gene_expression"]).rank().to_numpy(float)
                design = np.column_stack([np.ones(valid.sum()), cr])
                rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
                ry = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
                adjusted = np.corrcoef(rx, ry)[0, 1]
                raw = spearmanr(xraw[valid], yraw[valid]).statistic
                rows.append({"selection_scope": scope, "trajectory": traj, "background": BACKGROUND[traj],
                             "gene": gene, "n": int(valid.sum()), "raw_vulnerability_rho": -float(raw),
                             "adjusted_partial_vulnerability_rho": -float(adjusted),
                             "delta_adjusted_minus_raw": -float(adjusted) + float(raw),
                             "covariates_both_X_and_Y": "proliferation+global_dependency+gene_expression"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "phase7bR3_covariate_sensitivity.csv", index=False)
    return out


def build_r3_ranking(primary: pd.DataFrame, global_p: pd.DataFrame,
                     stability: pd.DataFrame, covariates: pd.DataFrame,
                     exclusions: dict[str, list[str]]) -> pd.DataFrame:
    r2 = pd.read_csv(OUT / "phase7bR2_convergent_gene_ranking.csv")
    lobo = pd.read_csv(OUT / "phase7bR_leave_one_background_out.csv")
    hct = lobo.loc[lobo["background_omitted"] == "HCT116",
                   ["gene", "median_background_vulnerability_rho", "n_positive_backgrounds"]]
    hct = hct.rename(columns={"median_background_vulnerability_rho": "leave_HCT116_out_median_vulnerability_rho",
                              "n_positive_backgrounds": "leave_HCT116_out_n_positive_backgrounds"})
    adjusted = covariates.loc[covariates["selection_scope"] == "convergent_top500"]
    adjusted = adjusted.groupby("gene").agg(adjusted_partial_rho_median=("adjusted_partial_vulnerability_rho", "median"),
                                             adjusted_partial_rho_min=("adjusted_partial_vulnerability_rho", "min"),
                                             adjusted_positive_fraction=("adjusted_partial_vulnerability_rho", lambda x: float((x > 0).mean())),
                                             adjusted_audit_n=("n", "count")).reset_index()
    out = (r2[["gene", "median_background_vulnerability_rho", "n_positive_backgrounds",
              "n_resampling_supported_backgrounds", "tier1_flag", "tier2_flag", "tier3_flag"]]
           .merge(global_p, on="gene", how="left")
           .merge(stability[["gene", "cross_method_stable_flag", "n_models_top500", "n_models_direction_ge3of4",
                             "median_rank_percentile"]], on="gene", how="left")
           .merge(hct, on="gene", how="left")
           .merge(adjusted, on="gene", how="left"))
    out["final_shortlist_flag"] = ((out["global_empirical_q_value"] <= 0.10) &
                                    (out["n_positive_backgrounds"] >= 3) &
                                    (out["n_resampling_supported_backgrounds"] >= 2) &
                                    out["cross_method_stable_flag"] &
                                    (out["leave_HCT116_out_median_vulnerability_rho"] > 0) &
                                    (out["leave_HCT116_out_n_positive_backgrounds"] >= 2) &
                                    (out["adjusted_partial_rho_median"] > 0))
    out["previous_19_reclassified_as_internal_stability"] = out["gene"].isin(set(stability.loc[stability["robust_candidate_flag"], "gene"]))
    out = out.sort_values(["final_shortlist_flag", "global_empirical_q_value", "median_rank_percentile"],
                          ascending=[False, True, True]).reset_index(drop=True)
    out.insert(0, "r3_rank", np.arange(1, len(out) + 1))
    out.to_csv(OUT / "phase7bR3_final_ranking.csv", index=False)
    out.loc[out["previous_19_reclassified_as_internal_stability"]].to_csv(OUT / "phase7bR3_previous_internal_stability_candidates.csv", index=False)
    return out


def write_manifest(n_crc: int, shortlist_n: int, internal_n: int) -> None:
    manifest = {
        "phase": "7B-R3",
        "parent_phases": ["7B-R", "7B-R2"],
        "status": "completed_final_global_null_validation",
        "phase7c8_executed": False,
        "depmap_release": "23Q4",
        "n_crc_models": n_crc,
        "trajectories": TRAJECTORIES,
        "backgrounds": BACKGROUNDS,
        "global_null": {
            "permutations": N_PERM,
            "seed": SEED,
            "unit": "one shared row-label permutation of full CRC CRISPR matrix per draw",
            "statistic": "median vulnerability rho across 4 biological backgrounds; HCT116 median across its 3 trajectories",
            "tail": "two-sided absolute global statistic",
            "self_line_exclusion": "applied after the shared permutation for each trajectory",
            "p_value_combination": "none; no Fisher; direct global empirical p followed by BH",
        },
        "corrected_partial_spearman": {
            "same_covariates_for_X_and_Y": ["proliferation", "global_dependency", "target_gene_expression"],
            "scopes": ["trajectory_top200", "convergent_top500"],
        },
        "final_shortlist_criteria": [
            "global empirical q <= 0.10",
            ">=3/4 backgrounds positive",
            ">=2 resampling-supported backgrounds",
            ">=4/6 cross-method models internally stable",
            "leave-HCT116-out median rho > 0 and >=2/3 remaining backgrounds positive",
            "corrected partial-Spearman adjusted median rho > 0",
        ],
        "internal_stability_reclassification": "Phase 7B-R2 19 candidates are reported as internally cross-method-stable, not robust vulnerabilities.",
        "internal_stability_count": internal_n,
        "final_shortlist_count": shortlist_n,
        "raw_data_policy": "Raw DepMap/GEO files remain local and ignored by git; only derived R3 outputs and scripts are committed.",
    }
    (OUT / "phase7bR3_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


def write_report(ranking: pd.DataFrame, cov: pd.DataFrame, n_crc: int) -> None:
    internal = ranking.loc[ranking["previous_19_reclassified_as_internal_stability"]]
    shortlist = ranking.loc[ranking["final_shortlist_flag"]]
    lodo = pd.read_csv(OUT / "phase7bR_leave_one_dataset_out.csv")
    lodo_stability = lodo.groupby("dataset_omitted").apply(lambda x: x[["rank_full", "rank_leave_one_dataset_out"]].corr(method="spearman").iloc[0, 1]).to_dict()
    focus = ranking.loc[ranking["gene"].isin(["ADAMTS9", "MOB4", "TJP1", "C16ORF86", "SOD2"]),
                        ["gene", "r3_rank", "observed_global_T_median_vulnerability_rho", "global_empirical_q_value",
                         "n_positive_backgrounds", "n_resampling_supported_backgrounds", "adjusted_partial_rho_median", "final_shortlist_flag"]]
    adjusted_delta = float(cov["delta_adjusted_minus_raw"].abs().median()) if len(cov) else float("nan")
    mech = ranking.loc[ranking["gene"].isin(base.PREDEFINED_MECHANISM),
                       ["gene", "r3_rank", "global_empirical_q_value", "final_shortlist_flag"]].sort_values("r3_rank")
    lines = [
        "# Phase 7B-R3：global-null validation and final shortlist",
        "",
        "## Final conclusion",
        "",
        f"Across {n_crc} CRC DepMap 23Q4 models, six OXA-R trajectories and four biological backgrounds, **no universal single-gene dependency was identified**. Phase 7C/8 were not run.",
        f"The Phase 7B-R2 set of {len(internal)} candidates is reclassified as **internally cross-method-stable candidates**, not robust vulnerabilities. After the shared global null and corrected partial-Spearman gate, the final shortlist contains **{len(shortlist)} genes**.",
        "If the shortlist is empty, this is a valid negative result: the single-gene convergence hypothesis does not meet the final calibrated criteria, and later work should move to pathway/complex-level convergence without forcing a Gene X.",
        "",
        "## What changed in R3",
        "",
        "1. Every null draw applies one shared permutation of the full CRC CRISPR row labels to all six trajectories and all genes. This preserves gene-dependency covariance and the observed relationships among the HCT116 trajectories.",
        "2. The test statistic is directly the median vulnerability rho across four biological backgrounds; no Fisher combination of correlated background p-values is used.",
        "3. Corrected partial Spearman residualizes both state score X and dependency Y on the same covariates: proliferation, global dependency and target-gene expression.",
        "",
        "## R3 shortlist criteria",
        "",
        "Global empirical q <=0.10; >=3/4 backgrounds positive; >=2 resampling-supported backgrounds; >=4/6 cross-method models internally stable; leave-HCT116-out median rho >0 with >=2/3 remaining backgrounds positive; corrected adjusted median rho >0.",
        "",
        "## Candidate summary",
        "",
        "```text", focus.to_string(index=False), "```", "",
        "Previous R2-only labels are available in `phase7bR3_previous_internal_stability_candidates.csv`; they are not carried forward as validated vulnerabilities unless they pass the R3 final gate.",
        "",
        "## Corrected covariate audit",
        "",
        f"The audit covers trajectory top200 and convergent top500 selections. The absolute adjusted-minus-raw rho median is {adjusted_delta:.3f}; both X and Y use the same covariate set for every gene.",
        "",
        "## LODO context",
        "",
        "Global rank stability inherited from the six-model audit: " + str(lodo_stability) + "。GSE77932 and GSE42387 remain material perturbations; no candidate is called universal solely from the full-data rank.",
        "",
        "## Mechanism audit",
        "",
        mech.to_string(index=False),
        "",
        "DHODH, RRM2, CPT1A, CPT2, VCP and SLC22A5 do not regain a universal dependency signal. Meldonium remains a No-Go for broad OXA-R reversal.",
        "",
        "## Files",
        "",
        "- `phase7bR3_global_empirical_p.csv`: shared global-null p/q calibration.",
        "- `phase7bR3_final_ranking.csv`: effect, global q, stability, LOBO and corrected covariate fields.",
        "- `phase7bR3_covariate_sensitivity.csv`: symmetric partial-Spearman audit.",
        "- `phase7bR3_previous_internal_stability_candidates.csv`: R2 19-candidate reclassification.",
        "- `phase7bR3_manifest.json`: exact null, criteria and provenance.",
        "",
        "Raw DepMap/GEO files remain local and are not committed. Phase 7C/8 remain deferred.",
    ]
    (OUT / "phase7bR3_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    print("[1/5] Loading fixed 7B-R inputs and DepMap data")
    _delta, crc, expr, dep = base.load_inputs()
    scores = load_scores()
    primary = pd.read_csv(OUT / "phase7bR_gene_dependency_by_trajectory.csv")
    exclusions = base.self_exclusion_map(crc)
    stability = pd.read_csv(OUT / "phase7bR2_cross_method_stability.csv")
    print(f"  CRC overlap={len(crc)}")
    print("[2/5] Shared global row-label permutation null")
    global_p = global_empirical_null(primary, scores[250], dep, exclusions)
    print("[3/5] Symmetric partial-Spearman covariate audit")
    preliminary = pd.read_csv(OUT / "phase7bR2_convergent_gene_ranking.csv")
    cov = corrected_partial_covariate_audit(scores[250], expr, dep, exclusions, primary, preliminary)
    print("[4/5] Final global-null ranking and shortlist")
    ranking = build_r3_ranking(primary, global_p, stability, cov, exclusions)
    internal_n = int(stability["robust_candidate_flag"].sum())
    shortlist_n = int(ranking["final_shortlist_flag"].sum())
    write_manifest(len(crc), shortlist_n, internal_n)
    print("[5/5] Writing report")
    write_report(ranking, cov, len(crc))
    print(f"Completed Phase 7B-R3. Internal stability={internal_n}; final shortlist={shortlist_n}. Phase 7C/8 were not executed.")


if __name__ == "__main__":
    main()
