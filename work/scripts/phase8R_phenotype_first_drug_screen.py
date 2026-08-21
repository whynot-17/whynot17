#!/usr/bin/env python
"""Phase 8-R: phenotype-first acquired OXA-R collateral sensitivity screen.

The biological ranking is deliberately blind to R3 genes, prior modules,
drug reputation and literature novelty.  It uses the fixed R3 resistance
signatures projected into GDSC-mapped CRC models, then ranks drugs by
cross-background sensitivity.  GDSC1 and GDSC2 are analysed independently;
PRISM and CTRPv2 are recorded as unavailable when no local files exist.

This script writes only derived outputs.  Regulatory and novelty fields are
created as post-ranking pending annotations so that they cannot influence the
phenotype ranking.
"""

from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
WORK = ROOT / "work"
OUT.mkdir(exist_ok=True)

DATA = WORK / "phase7_convergent_vulnerability" / "data"
GDSC_FILES = {
    "GDSC1": DATA / "GDSC1_fitted_dose_response_24Jul22.xlsx",
    "GDSC2": DATA / "GDSC2_fitted_dose_response_24Jul22.xlsx",
}
PRISM_PATH = WORK / "phase8R_prism" / "raw"
CTRP_PATH = WORK / "phase8R_ctrpv2" / "raw"
MODEL_FILE = WORK / "phase7b_depmap" / "raw" / "Model.csv"
EXPR_FILE = WORK / "phase7b_depmap" / "raw" / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
DELTA_FILE = WORK / "phase5_perturbation_reversal" / "gene_delta_matrix_primary.csv"
DETAILS_FILE = DATA / "Cell_Lines_Details.xlsx"
COMPOUND_FILE = DATA / "screened_compounds_rel_8.4.csv"
SIGNATURE_FILE = OUT / "phase7bR_trajectory_signatures.csv"

TRAJECTORIES = [
    "GSE77932|HCT116", "GSE77932|DLD1", "GSE42387|HCT116",
    "GSE42387|HT29", "GSE42387|LoVo", "GSE119603|HCT116",
]
BACKGROUND = {
    "GSE77932|HCT116": "HCT116", "GSE77932|DLD1": "DLD1",
    "GSE42387|HCT116": "HCT116", "GSE42387|HT29": "HT29",
    "GSE42387|LoVo": "LoVo", "GSE119603|HCT116": "HCT116",
}
BACKGROUNDS = ["HCT116", "DLD1", "HT29", "LoVo"]
DATASET_TRAJECTORY = {x: x.split("|")[0] for x in TRAJECTORIES}
SIGNATURE_SIZES = [100, 250, 500]
METHODS = ["weighted", "rank"]
PRIMARY_SIZE = 250
PRIMARY_METHOD = "weighted"
MIN_N = 20
PRIORITY_N = 30
N_PERM = 5000
SEED = 20260821


def clean_symbol(value: object) -> str:
    return re.sub(r"\s*\([^)]*\)$", "", str(value).strip()).upper()


def cosmic_key(value: object) -> str:
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return ""


def drug_key(value: object) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


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


def load_reference() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    model = pd.read_csv(MODEL_FILE, low_memory=False)
    lineage = model["OncotreeLineage"].fillna("").astype(str).str.lower()
    primary = model["OncotreePrimaryDisease"].fillna("").astype(str).str.lower()
    crc = model.loc[(lineage == "bowel") & primary.str.contains("colorectal", regex=False)].copy()
    crc["ModelID"] = crc["ModelID"].astype(str)
    crc["COSMIC_key"] = crc["COSMICID"].map(cosmic_key)

    expr = pd.read_csv(EXPR_FILE, index_col=0, low_memory=False)
    expr.index = expr.index.astype(str)
    names = [clean_symbol(x) for x in expr.columns]
    keep = ~pd.Index(names).duplicated()
    expr = expr.loc[:, keep].copy()
    expr.columns = np.asarray(names)[keep]
    expr = expr.apply(pd.to_numeric, errors="coerce")
    common = sorted(set(crc.ModelID) & set(expr.index))
    crc = crc.set_index("ModelID").loc[common]
    expr = expr.loc[common]

    delta = pd.read_csv(DELTA_FILE, low_memory=False)
    gcol = delta.columns[0]
    delta["gene"] = delta[gcol].map(clean_symbol)
    delta = delta.drop(columns=[gcol]).set_index("gene").apply(pd.to_numeric, errors="coerce")
    delta = delta.loc[~delta.index.duplicated(keep="first"), TRAJECTORIES]
    return crc, expr, delta


def load_msi() -> pd.DataFrame:
    if not DETAILS_FILE.exists():
        return pd.DataFrame(columns=["COSMIC_key", "MSI"])
    d = pd.read_excel(DETAILS_FILE, sheet_name="Cell line details", usecols=["Sample Name", "COSMIC identifier", "Microsatellite \ninstability Status (MSI)"])
    d = d.rename(columns={"Sample Name": "GDSC_cell_line", "Microsatellite \ninstability Status (MSI)": "MSI"})
    d["COSMIC_key"] = d["COSMIC identifier"].map(cosmic_key)
    return d[["COSMIC_key", "GDSC_cell_line", "MSI"]].drop_duplicates("COSMIC_key")


def self_exclusions(model: pd.DataFrame) -> dict[str, list[str]]:
    out = {}
    for traj in TRAJECTORIES:
        target = re.sub(r"[^A-Z0-9]", "", BACKGROUND[traj].upper())
        mask = pd.Series(False, index=model.index)
        for col in ["CellLineName", "StrippedCellLineName", "CCLEName"]:
            mask |= model[col].fillna("").astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True).str.contains(target, regex=False)
        out[traj] = model.index[mask].tolist()
    return out


def load_gdsc(dataset: str, model: pd.DataFrame, msi: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = GDSC_FILES[dataset]
    usecols = ["COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME", "LN_IC50", "AUC", "TCGA_DESC"]
    raw = pd.read_excel(path, usecols=usecols)
    raw["TCGA_DESC"] = raw["TCGA_DESC"].fillna("").astype(str).str.upper()
    raw = raw[raw["TCGA_DESC"].isin({"COREAD", "COAD", "READ", "COAD/READ"})].copy()
    raw["COSMIC_key"] = raw["COSMIC_ID"].map(cosmic_key)
    raw["DRUG_NAME"] = raw["DRUG_NAME"].astype(str).str.strip()
    raw = raw.merge(model.reset_index()[["ModelID", "COSMIC_key"]], on="COSMIC_key", how="inner")
    raw = raw.merge(msi, on="COSMIC_key", how="left")
    raw = raw.dropna(subset=["LN_IC50"])
    # Replicates are collapsed before association testing.
    response = (raw.groupby(["ModelID", "DRUG_NAME"], as_index=False)
                .agg(LN_IC50=("LN_IC50", "mean"), AUC=("AUC", "mean"),
                     DRUG_ID=("DRUG_ID", "first"), PUTATIVE_TARGET=("PUTATIVE_TARGET", "first"),
                     PATHWAY_NAME=("PATHWAY_NAME", "first"), COSMIC_key=("COSMIC_key", "first"),
                     GDSC_cell_line=("CELL_LINE_NAME", "first"), MSI=("MSI", "first")))
    mapping = raw[["ModelID", "COSMIC_key", "CELL_LINE_NAME", "MSI"]].drop_duplicates().rename(columns={"CELL_LINE_NAME": "GDSC_cell_line"})
    mapping["database"] = dataset
    mapping["expression_source"] = "DepMap 23Q4 OmicsExpressionProteinCodingGenesTPMLogp1 / CCLE-linked ModelID"
    return response, mapping


def compute_proliferation(expr: pd.DataFrame) -> pd.Series:
    gmt = WORK / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
    genes = set()
    if gmt.exists():
        for line in gmt.read_text(encoding="utf-8").splitlines():
            p = line.split("\t")
            if len(p) >= 3 and p[0] in {"HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"}:
                genes.update(clean_symbol(x) for x in p[2:])
    common = sorted(genes & set(expr.columns))
    z = (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=0).replace(0, 1.0)
    return z[common].mean(axis=1) if common else pd.Series(np.nan, index=expr.index)


def build_state_scores(expr: pd.DataFrame, delta: pd.DataFrame, mapped_ids: list[str],
                       model: pd.DataFrame, database: str) -> tuple[pd.DataFrame, dict[tuple[int, str], pd.DataFrame]]:
    mapped_ids = [x for x in mapped_ids if x in expr.index]
    e = expr.loc[mapped_ids]
    z = (e - e.mean(axis=0)) / e.std(axis=0, ddof=0).replace(0, 1.0)
    rank_e = e.rank(axis=1, pct=True) - 0.5
    weights = pd.read_csv(SIGNATURE_FILE)
    long = []
    wide_scores = {}
    for n in SIGNATURE_SIZES:
        for method in METHODS:
            w_all = weights[weights["signature_size"] == n]
            score = pd.DataFrame(index=e.index)
            for traj in TRAJECTORIES:
                w = w_all[w_all["trajectory"] == traj].set_index("gene")
                common = [g for g in w.index if g in e.columns]
                coeff = w.loc[common, "weight_scaled" if method == "weighted" else "rank_weight"].to_numpy(float)
                matrix = z[common] if method == "weighted" else rank_e[common]
                score[traj] = matrix.to_numpy(float).dot(coeff) / np.abs(coeff).sum()
            wide_scores[(n, method)] = score
            for traj in TRAJECTORIES:
                for model_id in score.index:
                    long.append({"database": database, "ModelID": model_id, "signature_size": n,
                                 "scoring_method": method, "trajectory": traj, "background": BACKGROUND[traj],
                                 "state_score": score.loc[model_id, traj]})
    out = pd.DataFrame(long)
    return out, wide_scores


def association_by_trajectory(response: pd.DataFrame, scores: dict[tuple[int, str], pd.DataFrame],
                              model: pd.DataFrame, exclusions: dict[str, list[str]], database: str) -> pd.DataFrame:
    rows = []
    for (n, method), score in scores.items():
        for traj in TRAJECTORIES:
            s = score[traj].rename("OXARScore").to_frame()
            keep = [x for x in s.index if x not in exclusions[traj]]
            s = s.loc[keep]
            r = response[response.ModelID.isin(keep)].copy()
            merged = r.merge(s, left_on="ModelID", right_index=True, how="inner")
            for drug, g in merged.groupby("DRUG_NAME", sort=False):
                x = g[["OXARScore", "LN_IC50"]].dropna()
                nn = len(x)
                rho = p = np.nan
                if nn >= MIN_N:
                    rho, p = spearmanr(x["OXARScore"], -x["LN_IC50"])
                rows.append({"database": database, "signature_size": n, "scoring_method": method,
                             "trajectory": traj, "background": BACKGROUND[traj], "DRUG_NAME": drug,
                             "drug_key": drug_key(drug), "rho_score_vs_sensitivity": rho,
                             "p_value": p, "n": nn, "n_priority": nn >= PRIORITY_N,
                             "power_status": "adequate_n>=30" if nn >= PRIORITY_N else ("low_power_n20_29" if nn >= MIN_N else "insufficient_n<20"),
                             "sensitivity_metric": "-LN_IC50", "direction_definition": "positive rho = OXA-R-like state more sensitive"})
    out = pd.DataFrame(rows)
    out["p_value_bh_within_database_signature"] = np.nan
    for key, idx in out.groupby(["database", "signature_size", "scoring_method"]).groups.items():
        out.loc[idx, "p_value_bh_within_database_signature"] = bh(out.loc[idx, "p_value"].to_numpy())
    return out


def background_aggregate(traj_assoc: pd.DataFrame) -> pd.DataFrame:
    a = traj_assoc.dropna(subset=["rho_score_vs_sensitivity"]).copy()
    out = (a.groupby(["database", "signature_size", "scoring_method", "DRUG_NAME", "drug_key", "background"], sort=False)
           .agg(n_trajectory_values=("rho_score_vs_sensitivity", "size"),
                n_positive_trajectory_values=("rho_score_vs_sensitivity", lambda x: int((x > 0).sum())),
                median_rho=("rho_score_vs_sensitivity", "median"),
                min_rho=("rho_score_vs_sensitivity", "min"), max_rho=("rho_score_vs_sensitivity", "max"),
                median_n=("n", "median"), median_p=("p_value", "median"))
           .reset_index())
    return out


def drug_rankings(bg: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["database", "signature_size", "scoring_method", "DRUG_NAME", "drug_key"]
    for key, g in bg.groupby(group_cols, sort=False):
        database, n, method, drug, dkey = key
        piv = g.set_index("background").reindex(BACKGROUNDS)
        r = piv["median_rho"].to_numpy(float)
        available = np.isfinite(r)
        nonhct = r[1:]
        npos = int(np.nansum(r > 0))
        npos_nonhct = int(np.nansum(nonhct > 0))
        med = float(np.nanmedian(r)) if available.any() else np.nan
        leave = float(np.nanmedian(nonhct)) if np.isfinite(nonhct).any() else np.nan
        if npos >= 3 and med >= 0.15 and leave > 0 and npos_nonhct >= 2 and available.sum() >= 4:
            gate = "TierA_strong_collateral"
        elif npos >= 3 and med >= 0.10 and available.sum() >= 3:
            gate = "TierB_moderate_candidate"
        elif npos >= 2 and med >= 0.20 and npos_nonhct >= 2 and available.sum() >= 2:
            gate = "TierC_subtype_candidate"
        else:
            gate = "none"
        rows.append({"database": database, "signature_size": n, "scoring_method": method,
                     "DRUG_NAME": drug, "drug_key": dkey, "n_backgrounds_available": int(available.sum()),
                     "n_positive_backgrounds": npos, "n_positive_nonHCT116_backgrounds": npos_nonhct,
                     "median_background_rho": med, "IQR_background_rho": float(np.nanquantile(r, .75) - np.nanquantile(r, .25)) if available.any() else np.nan,
                     "min_background_rho": float(np.nanmin(r)) if available.any() else np.nan,
                     "max_background_rho": float(np.nanmax(r)) if available.any() else np.nan,
                     "leave_HCT116_out_median_rho": leave, "primary_discovery_gate": gate,
                     "tested_status": "tested" if available.sum() else "unavailable"})
    out = pd.DataFrame(rows)
    out["rank_within_database_signature"] = np.nan
    for key, idx in out.groupby(["database", "signature_size", "scoring_method"]).groups.items():
        order = out.loc[idx].sort_values(["median_background_rho", "n_positive_backgrounds"], ascending=[False, False]).index
        out.loc[order, "rank_within_database_signature"] = np.arange(1, len(order) + 1)
    return out


def global_empirical_drug_null(response: pd.DataFrame, score: pd.DataFrame, rankings: pd.DataFrame,
                               exclusions: dict[str, list[str]], database: str) -> pd.DataFrame:
    primary = rankings.loc[(rankings.database == database) & (rankings.signature_size == PRIMARY_SIZE) & (rankings.scoring_method == PRIMARY_METHOD)].copy()
    effect = primary.set_index("DRUG_NAME")["median_background_rho"]
    drugs = effect.index.to_list()
    full_ids = score.index.to_list()
    pos = {x: i for i, x in enumerate(full_ids)}
    resp = response.pivot(index="ModelID", columns="DRUG_NAME", values="LN_IC50").reindex(full_ids)
    # Full-frame sensitivity ranks are permuted as a matrix, preserving drug-drug covariance.
    sens_rank = (-resp).rank(axis=0, method="average", na_option="keep").to_numpy(float)
    traj_prepared = []
    for traj in TRAJECTORIES:
        keep = [x for x in full_ids if x not in exclusions[traj]]
        positions = np.array([pos[x] for x in keep], dtype=int)
        xr = pd.Series(score.loc[keep, traj].to_numpy(float)).rank().to_numpy(float)
        traj_prepared.append((traj, positions, xr))
    rng = np.random.default_rng(SEED + (1 if database == "GDSC1" else 2))
    null_ge = np.zeros(len(drugs), dtype=np.int64)
    drug_pos = {d: i for i, d in enumerate(resp.columns)}
    col_idx = np.array([drug_pos[d] for d in drugs if d in drug_pos], dtype=int)
    valid_drugs = [d for d in drugs if d in drug_pos]
    obs = effect.reindex(valid_drugs).to_numpy(float)
    n_full = len(full_ids)
    for start in range(0, N_PERM, 100):
        k = min(100, N_PERM - start)
        perms = np.stack([rng.permutation(n_full) for _ in range(k)], axis=0)
        v_by_traj = []
        for _traj, positions, xr in traj_prepared:
            yr = sens_rank[perms[:, positions], :][:, :, col_idx]
            xp = xr - xr.mean()
            xmat = np.broadcast_to(xp[None, :, None], yr.shape)
            valid = np.isfinite(yr)
            n = valid.sum(axis=1).astype(float)
            xm = np.where(valid, xmat, np.nan).sum(axis=1) / np.maximum(n, 1)
            ym = np.where(valid, yr, np.nan).sum(axis=1) / np.maximum(n, 1)
            xc = np.where(valid, xmat - xm[:, None, :], 0.0)
            yc = np.where(valid, yr - ym[:, None, :], 0.0)
            num = (xc * yc).sum(axis=1)
            den = np.sqrt((xc * xc).sum(axis=1) * (yc * yc).sum(axis=1))
            v_by_traj.append(-(num / np.maximum(den, 1e-15)))
        bg_values = []
        for bg in BACKGROUNDS:
            idx = [i for i, traj in enumerate(TRAJECTORIES) if BACKGROUND[traj] == bg]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                bg_values.append(np.nanmedian(np.stack([v_by_traj[i] for i in idx], axis=0), axis=0))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            global_t = np.nanmedian(np.stack(bg_values, axis=0), axis=0)
        null_ge += (np.abs(global_t) >= np.abs(obs)[None, :]).sum(axis=0)
    p = (1 + null_ge) / (N_PERM + 1)
    out = primary.loc[primary.DRUG_NAME.isin(valid_drugs)].copy()
    out = out.merge(pd.DataFrame({"DRUG_NAME": valid_drugs, "global_empirical_p_value": p}), on="DRUG_NAME", how="left")
    out["global_empirical_q_value"] = bh(out["global_empirical_p_value"].to_numpy())
    out["global_null_permutations"] = N_PERM
    out["global_null_unit"] = "shared row-label permutation of full sensitivity matrix across all drugs/trajectories"
    return out


def signature_robustness(rankings: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (db, drug, dkey), g in rankings.groupby(["database", "DRUG_NAME", "drug_key"], sort=False):
        g = g.drop_duplicates(["signature_size", "scoring_method"])
        if not len(g):
            continue
        rows.append({"database": db, "DRUG_NAME": drug, "drug_key": dkey,
                     "n_models_available": len(g), "n_models_positive_median_rho": int((g.median_background_rho > 0).sum()),
                     "n_models_direction_consistent_ge3of4": int(((g.n_positive_backgrounds >= 3)).sum()),
                     "n_models_top50": int((g.rank_within_database_signature <= 50).sum()),
                     "n_models_top100": int((g.rank_within_database_signature <= 100).sum()),
                     "n_models_top200": int((g.rank_within_database_signature <= 200).sum()),
                     "median_model_rho": float(g.median_background_rho.median()),
                     "median_model_rank": float(g.rank_within_database_signature.median()),
                     "signature_robust_flag": bool(len(g) >= 4 and (g.median_background_rho > 0).sum() >= 4),
                     "model_combinations": ";".join(f"top{r.signature_size}_{r.scoring_method}" for r in g.itertuples())})
    return pd.DataFrame(rows)


def partial_covariate_drugs(primary_global: pd.DataFrame, score: pd.DataFrame, response: pd.DataFrame,
                            model: pd.DataFrame, expr: pd.DataFrame, mapping: pd.DataFrame,
                            exclusions: dict[str, list[str]], database: str,
                            extra_drugs: set[str] | None = None) -> pd.DataFrame:
    top_drugs = set(primary_global.sort_values(["global_empirical_q_value", "median_background_rho"]).head(100)["DRUG_NAME"])
    if extra_drugs:
        top_drugs |= set(extra_drugs)
    prolif = compute_proliferation(expr)
    msi_map = mapping.drop_duplicates("ModelID").set_index("ModelID")["MSI"]
    rows = []
    for traj in TRAJECTORIES:
        keep = [x for x in score.index if x not in exclusions[traj]]
        for drug in sorted(top_drugs):
            g = response.loc[response["DRUG_NAME"] == drug].set_index("ModelID").reindex(keep)
            y = -g["LN_IC50"]
            x = score.loc[keep, traj]
            cov = pd.DataFrame({"proliferation": prolif.reindex(keep), "MSI": msi_map.reindex(keep).map(lambda z: 1.0 if isinstance(z, str) and "MSI-H" in z.upper() else (0.0 if isinstance(z, str) and z else np.nan))}, index=keep)
            # Only use measured covariates; mutation calls are unavailable in
            # the current local input and remain explicitly annotated as such.
            cov = cov.loc[:, cov.notna().sum() >= MIN_N]
            valid = x.notna() & y.notna() & cov.notna().all(axis=1)
            if valid.sum() < MIN_N:
                continue
            xr = x[valid].rank().to_numpy(float)
            yr = y[valid].rank().to_numpy(float)
            cr = cov.loc[valid].rank().to_numpy(float)
            design = np.column_stack([np.ones(valid.sum()), cr])
            rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
            ry = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
            rows.append({"database": database, "trajectory": traj, "background": BACKGROUND[traj],
                         "DRUG_NAME": drug, "n": int(valid.sum()), "raw_rho": spearmanr(x[valid], y[valid]).statistic,
                         "adjusted_partial_rho": np.corrcoef(rx, ry)[0, 1],
                         "covariates": "+".join(cov.columns), "KRAS_status": "unavailable", "BRAF_status": "unavailable"})
    return pd.DataFrame(rows)


def cross_database(primary_globals: pd.DataFrame) -> pd.DataFrame:
    all_drugs = sorted(set(primary_globals.DRUG_NAME))
    rows = []
    for drug in all_drugs:
        sub = primary_globals[primary_globals.DRUG_NAME == drug].set_index("database")
        states = {}
        for db in ["GDSC1", "GDSC2"]:
            if db in sub.index:
                r = sub.loc[db]
                states[db] = {"status": "tested", "T": float(r.median_background_rho), "q": float(r.global_empirical_q_value),
                              "positive": int(r.n_positive_backgrounds), "gate": r.primary_discovery_gate}
            else:
                states[db] = {"status": "drug_absent_or_insufficient", "T": np.nan, "q": np.nan, "positive": np.nan, "gate": "unavailable"}
        tested = [db for db in ["GDSC1", "GDSC2"] if states[db]["status"] == "tested"]
        strict = [db for db in tested if states[db]["q"] <= 0.10 and states[db]["T"] > 0 and states[db]["positive"] >= 3]
        signs = [np.sign(states[db]["T"]) for db in tested if np.isfinite(states[db]["T"])]
        if len(strict) >= 2:
            level = "Level1_cross_database_strict"
        elif len(strict) == 1 and len(tested) >= 2 and all(x >= 0 for x in signs):
            level = "Level2_one_strict_other_same_direction"
        elif len(tested) == 1:
            level = "Level3_single_dataset"
        elif len(tested) >= 2 and any(x < 0 for x in signs) and any(x > 0 for x in signs):
            level = "Reject_opposite_direction"
        else:
            level = "exploratory_same_direction_or_underpowered"
        rows.append({"DRUG_NAME": drug, "drug_key": drug_key(drug),
                     "GDSC1_status": states["GDSC1"]["status"], "GDSC1_global_q": states["GDSC1"]["q"], "GDSC1_T": states["GDSC1"]["T"],
                     "GDSC2_status": states["GDSC2"]["status"], "GDSC2_global_q": states["GDSC2"]["q"], "GDSC2_T": states["GDSC2"]["T"],
                     "PRISM_status": "unavailable_local", "CTRPv2_status": "unavailable_local",
                     "replication_level": level})
    return pd.DataFrame(rows)


def biological_shortlist(primary_out: pd.DataFrame, replication: pd.DataFrame) -> pd.DataFrame:
    """Freeze a stringent phenotype-only shortlist before drug annotation.

    The existing Level 1 label is intentionally permissive enough for an
    exploratory replication table.  This table is the stricter post-screen
    gate: both GDSC databases must independently pass the four-background,
    leave-HCT116-out and six-signature robustness checks.  No target, drug
    class, regulatory or literature field enters this gate.
    """
    p = primary_out.copy()
    p["strict_database_gate"] = (
        (p["global_empirical_q_value"] <= 0.10)
        & (p["median_background_rho"] > 0)
        & (p["n_backgrounds_available"] == 4)
        & (p["n_positive_backgrounds"] >= 3)
        & (p["n_positive_nonHCT116_backgrounds"] >= 2)
        & (p["leave_HCT116_out_median_rho"] > 0)
        & (p["n_models_direction_consistent_ge3of4"] >= 4)
        & p["signature_robust_flag"].fillna(False)
    )
    p = p[p["strict_database_gate"]].copy()
    rows = []
    for drug, g in p.groupby(["DRUG_NAME", "drug_key"], sort=False):
        by_db = g.set_index("database")
        passed = sorted(by_db.index.tolist())
        if not {"GDSC1", "GDSC2"}.issubset(passed):
            continue
        rep = replication.loc[replication["DRUG_NAME"] == drug[0], "replication_level"]
        rows.append({
            "DRUG_NAME": drug[0], "drug_key": drug[1],
            "biological_status": "robust_cross_database_shortlist",
            "replication_level": rep.iloc[0] if len(rep) else "not_available",
            "databases_passing": ";".join(passed),
            "n_databases_passing": len(passed),
            "mean_median_background_rho": float(g["median_background_rho"].mean()),
            "min_median_background_rho": float(g["median_background_rho"].min()),
            "mean_global_empirical_q_value": float(g["global_empirical_q_value"].mean()),
            "GDSC1_median_background_rho": float(by_db.loc["GDSC1", "median_background_rho"]),
            "GDSC2_median_background_rho": float(by_db.loc["GDSC2", "median_background_rho"]),
            "GDSC1_global_empirical_q_value": float(by_db.loc["GDSC1", "global_empirical_q_value"]),
            "GDSC2_global_empirical_q_value": float(by_db.loc["GDSC2", "global_empirical_q_value"]),
            "GDSC1_direction_consistent_models": int(by_db.loc["GDSC1", "n_models_direction_consistent_ge3of4"]),
            "GDSC2_direction_consistent_models": int(by_db.loc["GDSC2", "n_models_direction_consistent_ge3of4"]),
            "GDSC1_signature_top100_models": int(by_db.loc["GDSC1", "n_models_top100"]),
            "GDSC2_signature_top100_models": int(by_db.loc["GDSC2", "n_models_top100"]),
        })
    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values(["mean_median_background_rho", "mean_global_empirical_q_value"], ascending=[False, True])
        out.insert(0, "biological_rank", np.arange(1, len(out) + 1))
    return out


def make_annotations(primary_globals: pd.DataFrame) -> pd.DataFrame:
    comp = pd.read_csv(COMPOUND_FILE) if COMPOUND_FILE.exists() else pd.DataFrame()
    rows = []
    for drug in sorted(set(primary_globals.DRUG_NAME)):
        sub = primary_globals[primary_globals.DRUG_NAME == drug].iloc[0]
        hit = comp.loc[comp["DRUG_NAME"].astype(str).str.upper() == str(drug).upper()] if len(comp) else pd.DataFrame()
        rows.append({"DRUG_NAME": drug, "drug_key": drug_key(drug), "GDSC_drug_id": sub.get("DRUG_ID", np.nan),
                     "putative_target": ";".join(hit.get("TARGET", pd.Series(dtype=str)).dropna().astype(str).unique()) if len(hit) else sub.get("PUTATIVE_TARGET", ""),
                     "pathway": ";".join(hit.get("TARGET_PATHWAY", pd.Series(dtype=str)).dropna().astype(str).unique()) if len(hit) else sub.get("PATHWAY_NAME", ""),
                     "regulatory_status": "pending_post_ranking_audit", "oncology_status": "pending_post_ranking_audit",
                     "clinical_usability": "pending_post_ranking_audit", "annotation_source": "local GDSC/compound metadata only; no ranking influence"})
    return pd.DataFrame(rows)


def write_report(primary_globals: pd.DataFrame, rankings: pd.DataFrame, replication: pd.DataFrame,
                 biological: pd.DataFrame, availability: dict, mapping: pd.DataFrame) -> None:
    p = primary_globals[(primary_globals["median_background_rho"] > 0)
                        & (primary_globals["n_positive_backgrounds"] >= 3)].copy()
    p = p.sort_values(["median_background_rho", "global_empirical_q_value"], ascending=[False, True]).head(20)
    gates = primary_globals["primary_discovery_gate"].value_counts().to_dict()
    rep = replication["replication_level"].value_counts().to_dict()
    lines = [
        "# Phase 8-R：Phenotype-first Drug X screen",
        "",
        "## Current scope and conclusion",
        "",
        "This phase starts from the fixed acquired OXA-R state and screens drug sensitivity before using R3 genes, old modules, regulatory status or literature novelty. GDSC1 and GDSC2 were analysed separately; PRISM and CTRPv2 were unavailable locally and were not treated as negative results.",
        f"Mapped GDSC CRC cell models: {mapping.ModelID.nunique()} unique models; GDSC1/GDSC2 mapping and drug coverage are recorded in `phase8R_cell_line_mapping.csv`. Primary model: top250 weighted state score, -LN_IC50 as sensitivity, n>=20 association floor and n>=30 priority.",
        f"Primary GDSC discovery gates: {gates}. Cross-database replication labels: {rep}.",
        f"The stringent phenotype-only gate leaves {len(biological)} robust cross-database shortlist drugs: {', '.join(biological['DRUG_NAME'].tolist()) if len(biological) else 'none'}. This is a phenotype-level screen, not a claim of acquired OXA-R experimental collateral sensitivity.",
        "",
        "## Data availability",
        "",
        json.dumps(availability, indent=2, ensure_ascii=False),
        "",
        "## Primary positive ranking (biological ranking frozen before annotation)",
        "",
        "The table is restricted to positive multi-background associations; negative associations with small two-sided empirical q-values are not candidates.",
        "```text", p[["database", "DRUG_NAME", "median_background_rho", "n_positive_backgrounds", "n_backgrounds_available", "global_empirical_q_value", "primary_discovery_gate"]].to_string(index=False), "```",
        "",
        "## Frozen phenotype-only shortlist",
        "",
        "```text", biological.to_string(index=False) if len(biological) else "none", "```",
        "",
        "## Direction and sensitivity conventions",
        "",
        "GDSC LN_IC50 is converted to sensitivity as `-LN_IC50`; positive rho means a higher acquired OXA-R-like state score is associated with greater drug sensitivity. HCT116 is aggregated as the median of its three trajectories; the primary evidence unit is the four biological backgrounds.",
        "",
        "## Boundaries",
        "",
        "R3 genes, FAO/DHODH/ferroptosis/ERAD hypotheses, previous LINCS results and drug novelty were not used to rank drugs. Regulatory/non-oncology fields are post-ranking annotations. PRISM/CTRP are unavailable, not negative. Raw GDSC and DepMap files remain local and are not committed.",
    ]
    (OUT / "phase8R_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    print("[1/8] Loading CRC expression, R3 signatures and metadata")
    model, expr, delta = load_reference()
    msi = load_msi()
    exclusions = self_exclusions(model)
    availability = {"GDSC1": "available_local", "GDSC2": "available_local",
                    "PRISM": "unavailable_local_not_negative", "CTRPv2": "unavailable_local_not_negative"}
    print(f"  CRC expression models={len(model)}; genes={expr.shape[1]}")
    responses = {}
    mappings = []
    score_wide = {}
    score_long = []
    print("[2/8] Mapping GDSC1/GDSC2 by COSMIC ID")
    for db in ["GDSC1", "GDSC2"]:
        response, mapping = load_gdsc(db, model, msi)
        responses[db] = response
        mappings.append(mapping)
        ids = sorted(set(response.ModelID))
        slong, swide = build_state_scores(expr, delta, ids, model, db)
        score_long.append(slong)
        score_wide[db] = swide
    mapping_all = pd.concat(mappings, ignore_index=True).drop_duplicates(["database", "ModelID", "GDSC_cell_line"])
    mapping_all["self_exclusion_trajectories"] = mapping_all.apply(
        lambda r: ";".join(traj for traj in TRAJECTORIES if r.ModelID in set(exclusions.get(traj, []))), axis=1
    )
    mapping_all.to_csv(OUT / "phase8R_cell_line_mapping.csv", index=False)
    pd.concat(score_long, ignore_index=True).to_csv(OUT / "phase8R_state_scores_by_dataset.csv", index=False)
    print("[3/8] Blind trajectory-level drug associations")
    traj_rows = []
    bg_rows = []
    rank_rows = []
    for db in ["GDSC1", "GDSC2"]:
        assoc = association_by_trajectory(responses[db], score_wide[db], model, exclusions, db)
        traj_rows.append(assoc)
    traj_all = pd.concat(traj_rows, ignore_index=True)
    traj_all.to_csv(OUT / "phase8R_drug_association_by_trajectory.csv", index=False)
    bg_all = background_aggregate(traj_all)
    bg_all.to_csv(OUT / "phase8R_drug_association_by_background.csv", index=False)
    rankings = drug_rankings(bg_all)
    rankings.to_csv(OUT / "phase8R_global_empirical_drug_ranking.csv", index=False)
    print("[4/8] Shared drug-matrix empirical null (5000 permutations per GDSC dataset)")
    primary_globals = []
    for db in ["GDSC1", "GDSC2"]:
        pg = global_empirical_drug_null(responses[db], score_wide[db][(PRIMARY_SIZE, PRIMARY_METHOD)], rankings, exclusions, db)
        primary_globals.append(pg)
    primary_globals = pd.concat(primary_globals, ignore_index=True)
    # Merge primary empirical p/q back into the all-model ranking table.
    rankings = rankings.merge(primary_globals[["database", "DRUG_NAME", "global_empirical_p_value", "global_empirical_q_value", "global_null_permutations", "global_null_unit"]],
                              on=["database", "DRUG_NAME"], how="left")
    rankings.to_csv(OUT / "phase8R_global_empirical_drug_ranking.csv", index=False)
    print("[5/8] Signature robustness and leave-HCT116 summaries")
    robustness = signature_robustness(rankings)
    robustness.to_csv(OUT / "phase8R_signature_method_robustness.csv", index=False)
    primary_out = primary_globals.merge(robustness, on=["database", "DRUG_NAME", "drug_key"], how="left")
    primary_out.to_csv(OUT / "phase8R_leave_hct116_out.csv", index=False)
    print("[6/8] Covariate sensitivity and cross-database replication")
    replication = cross_database(primary_globals)
    biological = biological_shortlist(primary_out, replication)
    shortlist_drugs = set(biological["DRUG_NAME"]) if len(biological) else set()
    cov_rows = []
    for db in ["GDSC1", "GDSC2"]:
        cov_rows.append(partial_covariate_drugs(primary_globals[primary_globals.database == db], score_wide[db][(PRIMARY_SIZE, PRIMARY_METHOD)], responses[db], model, expr, mappings[[i for i, x in enumerate(["GDSC1", "GDSC2"]) if x == db][0]], exclusions, db, shortlist_drugs))
    cov = pd.concat([x for x in cov_rows if len(x)], ignore_index=True) if any(len(x) for x in cov_rows) else pd.DataFrame()
    cov.to_csv(OUT / "phase8R_covariate_adjusted_drug_associations.csv", index=False)
    replication.to_csv(OUT / "phase8R_cross_database_replication.csv", index=False)
    biological.to_csv(OUT / "phase8R_biological_shortlist.csv", index=False)
    annotations = make_annotations(primary_globals)
    annotations.to_csv(OUT / "phase8R_drug_regulatory_annotation.csv", index=False)
    pd.DataFrame(columns=annotations.columns.tolist() + ["nononcology_shortlist_status"]).to_csv(OUT / "phase8R_nononcology_shortlist.csv", index=False)
    print("[7/8] Writing phenotype-first report and pending post-ranking files")
    write_report(primary_globals, rankings, replication, biological, availability, mapping_all)
    pd.DataFrame(columns=["DRUG_NAME", "novelty_class", "evidence", "source"]).to_csv(OUT / "phase8R_final_drugX_candidates.csv", index=False)
    (OUT / "phase8R_novelty_audit.md").write_text("# Phase 8-R novelty audit\n\nPending until the biological ranking and post-ranking regulatory filter are frozen.\n", encoding="utf-8")
    (OUT / "phase8R_go_nogo.md").write_text("# Phase 8-R Go/No-Go\n\nPending post-ranking regulatory/non-oncology and novelty audit.\n", encoding="utf-8")
    manifest = {"phase": "8-R", "status": "biological_screen_completed_post_ranking_annotation_pending", "phase7c8_executed": False,
                "data_availability": availability, "crc_rule": "GDSC TCGA_DESC COREAD/COAD/READ/COAD-READ mapped to DepMap CRC expression models by COSMIC ID",
                "n_mapped_models_by_database": {db: int(mapping_all.loc[mapping_all.database == db, "ModelID"].nunique()) for db in ["GDSC1", "GDSC2"]},
                "signature_sizes": SIGNATURE_SIZES, "methods": METHODS, "primary": "top250_weighted",
                "min_n": MIN_N, "priority_n": PRIORITY_N, "permutations": N_PERM,
                "sensitivity_direction": "-LN_IC50; positive rho = higher OXA-R-like score, greater sensitivity",
                "backgrounds": BACKGROUNDS, "self_line_exclusion": exclusions,
                "blind_ranking_policy": "R3 genes/modules, drug prior, regulatory status and literature novelty excluded from biological ranking",
                "post_ranking_steps_pending": ["regulatory annotation", "non-oncology filter", "novelty audit", "final Go/No-Go"],
                "raw_data_policy": "Raw GDSC/DepMap files remain local and ignored by git; only derived outputs and script are committed."}
    (OUT / "phase8R_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print("[8/8] Phase 8-R biological screen completed; post-ranking annotation is pending")


if __name__ == "__main__":
    main()
