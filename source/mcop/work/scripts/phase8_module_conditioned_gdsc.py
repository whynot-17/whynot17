"""Phase 8 first pass: module-conditioned pharmacological convergence.

This is a phenotype-level bridge, not a wet-lab prediction.  It asks whether
CRC cell lines with a high score for the positive OXA-R trajectories associated
with a Phase 7C module are more sensitive to a GDSC drug.  The module is used
as a biologically defined trajectory set; no single-gene or known-target
intersection is imposed.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
DATA = ROOT / "work" / "phase7_convergent_vulnerability" / "data"
DEP = ROOT / "work" / "phase7b_depmap" / "raw" / "Model.csv"
SCORES = OUT / "phase7b_crc_trajectory_scores.csv"
TRAJECTORIES = [
    "GSE77932|HCT116", "GSE77932|DLD1", "GSE42387|HCT116",
    "GSE42387|HT29", "GSE42387|LoVo", "GSE119603|HCT116",
]


def cosmic_key(x: object) -> str:
    try:
        return str(int(float(x)))
    except Exception:
        return str(x).strip()


def bh(p: pd.Series) -> pd.Series:
    x = pd.to_numeric(p, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    ok = x.notna()
    if not ok.any():
        return out
    vals = x[ok].to_numpy(float)
    order = np.argsort(vals)
    q = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1].clip(0, 1)
    tmp = np.empty_like(q)
    tmp[order] = q
    out.loc[ok] = tmp
    return out


def load_model_map() -> dict[str, str]:
    x = pd.read_csv(DEP, usecols=["ModelID", "COSMICID", "OncotreeLineage"])
    x["cosmic"] = x["COSMICID"].map(cosmic_key)
    # DepMap 23Q4 labels colorectal models under the broad Oncotree lineage
    # `Bowel`; the primary-disease field carries the specific colorectal label.
    x = x[x["OncotreeLineage"].astype(str).str.contains("Bowel", case=False, na=False)]
    return x.drop_duplicates("cosmic").set_index("cosmic")["ModelID"].to_dict()


def load_gdsc(path: Path, model_map: dict[str, str]) -> pd.DataFrame:
    cols = ["COSMIC_ID", "CELL_LINE_NAME", "TCGA_DESC", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME", "LN_IC50"]
    x = pd.read_excel(path, usecols=cols)
    x = x[x["TCGA_DESC"].astype(str).str.upper().isin({"COREAD", "COAD", "READ", "COAD/READ"})].copy()
    x["sample_id"] = x["COSMIC_ID"].map(cosmic_key).map(model_map)
    x = x.dropna(subset=["sample_id", "LN_IC50", "DRUG_NAME"])
    return x.groupby(["sample_id", "DRUG_NAME"], as_index=False).agg(
        LN_IC50=("LN_IC50", "mean"),
        CELL_LINE_NAME=("CELL_LINE_NAME", "first"),
        PUTATIVE_TARGET=("PUTATIVE_TARGET", "first"),
        PATHWAY_NAME=("PATHWAY_NAME", "first"),
    )


def load_module_focus() -> tuple[pd.DataFrame, pd.DataFrame]:
    audit = pd.read_csv(OUT / "phase7c_robustness_by_module.csv")
    nes = pd.read_csv(OUT / "phase7c_trajectory_module_gsea.csv")
    nes = nes.pivot_table(index="module", columns="trajectory", values="nes", aggfunc="first").reindex(columns=TRAJECTORIES)
    # Tier 1 = 6/6 direction-stable; Tier 2 = 5/6 direction-stable. Both
    # require at least two trajectory-level FDR hits and |median NES|>=0.75.
    focus = audit[(audit["full_direction"].eq("positive")) & (audit["full_consistency"] >= (5 / 6)) & (audit["full_n_fdr_qval_025"] >= 2) & (audit["full_median_nes"].abs() >= 0.75)].copy()
    focus["tier"] = np.where(focus["full_consistency"] >= 0.999, "T1_universal", "T2_subtype")
    focus["positive_trajectories"] = focus["module"].map(lambda m: ";".join([t for t in TRAJECTORIES if float(nes.loc[m, t]) > 0]) if m in nes.index else "")
    return focus, nes


def main() -> None:
    model_map = load_model_map()
    score = pd.read_csv(SCORES)
    score["sample_id"] = score["sample_id"].astype(str)
    datasets = []
    for name in ["GDSC1", "GDSC2"]:
        path = DATA / f"{name}_fitted_dose_response_24Jul22.xlsx"
        d = load_gdsc(path, model_map)
        d["dataset"] = name
        datasets.append(d)
    gdsc = pd.concat(datasets, ignore_index=True)
    # Align the GDSC records to the six trajectory scores once.  Repeating a
    # dataframe merge inside every drug x trajectory loop is unnecessarily
    # expensive for the full GDSC panel.
    aligned = gdsc.merge(score, on="sample_id", how="inner")
    focus, nes = load_module_focus()
    # Compute the elementary drug-sensitivity association once per
    # dataset×drug×trajectory.  Modules are applied as a later projection.
    base_rows = []
    for dataset in ["GDSC1", "GDSC2"]:
        dset = aligned[aligned["dataset"].eq(dataset)]
        for drug, drugset in dset.groupby("DRUG_NAME", sort=False):
            target = drugset["PUTATIVE_TARGET"].dropna().astype(str).iloc[0] if drugset["PUTATIVE_TARGET"].notna().any() else ""
            pathway = drugset["PATHWAY_NAME"].dropna().astype(str).iloc[0] if drugset["PATHWAY_NAME"].notna().any() else ""
            for trajectory in TRAJECTORIES:
                x = drugset[[trajectory, "LN_IC50"]].dropna()
                if len(x) < 15:
                    continue
                rho, p = spearmanr(x[trajectory], x["LN_IC50"])
                base_rows.append({"dataset": dataset, "trajectory": trajectory, "drug": drug, "n": len(x), "rho_state_vs_lnIC50": float(rho), "p": float(p), "putative_target": target, "pathway": pathway})
    base = pd.DataFrame(base_rows)
    module_rows = []
    for _, m in focus.iterrows():
        trajs = [t for t in TRAJECTORIES if t in str(m["positive_trajectories"]).split(";")]
        x = base[base["trajectory"].isin(trajs)].copy()
        x.insert(0, "module", m["module"])
        x.insert(1, "tier", m["tier"])
        module_rows.append(x)
    pairs = pd.concat(module_rows, ignore_index=True) if module_rows else pd.DataFrame()
    pairs.to_csv(OUT / "phase8_module_conditioned_gdsc_pairwise.csv", index=False)
    rows = []
    if len(pairs):
        for (module, drug), x in pairs.groupby(["module", "drug"], sort=False):
            pmin = float(x["p"].min())
            rows.append({
                "module": module,
                "tier": x["tier"].iloc[0],
                "drug": drug,
                "n_pairs": len(x),
                "n_datasets": x["dataset"].nunique(),
                "n_trajectories": x["trajectory"].nunique(),
                "median_rho_state_vs_lnIC50": float(x["rho_state_vs_lnIC50"].median()),
                "mean_rho_state_vs_lnIC50": float(x["rho_state_vs_lnIC50"].mean()),
                "fraction_negative_rho": float((x["rho_state_vs_lnIC50"] < 0).mean()),
                "min_pairwise_p": pmin,
                "putative_target": x["putative_target"].dropna().astype(str).replace("nan", "").head(1).tolist()[0] if len(x["putative_target"].dropna()) else "",
                "pathway": x["pathway"].dropna().astype(str).replace("nan", "").head(1).tolist()[0] if len(x["pathway"].dropna()) else "",
                "drug_sensitivity_score": float(-x["rho_state_vs_lnIC50"].median() * (x["rho_state_vs_lnIC50"] < 0).mean()),
            })
    ranking = pd.DataFrame(rows)
    if len(ranking):
        ranking["within_module_min_p_BH"] = ranking.groupby("module")["min_pairwise_p"].transform(bh)
        ranking = ranking.sort_values(["module", "drug_sensitivity_score", "n_pairs"], ascending=[True, False, False])
    ranking.to_csv(OUT / "phase8_module_conditioned_gdsc_drug_ranking.csv", index=False)
    robust = ranking[(ranking["n_datasets"] >= 2) & (ranking["n_trajectories"] >= 3) & (ranking["fraction_negative_rho"] >= 0.60)].copy() if len(ranking) else ranking
    robust.to_csv(OUT / "phase8_module_conditioned_gdsc_robust_candidates.csv", index=False)
    lines = ["# Phase 8：module-conditioned pharmacological convergence", "", "## 定义", "", "- 只使用 Phase 7C 正向模块：高 OXA-R trajectory score 应对应更高药物敏感性，即 state score 与 LN_IC50 的 rho < 0。", "- T1：6/6 trajectory 方向一致；T2：5/6 trajectory 方向一致。", "- 药物排序同时要求至少 2 个 GDSC dataset、至少 3 条 trajectory、负 rho 比例≥0.60，避免单一背景伪阳性。", "- 这是 phenotype-level pharmacogenomic prioritization，不等于 paired OXA-R CRISPR 或临床有效性。", "", f"Focus modules: {len(focus)}；pairwise records: {len(pairs)}；drug-module records: {len(ranking)}；robust candidate records: {len(robust)}。", "", "## Robust candidates", "", "| Tier | Module | Drug | n pairs | datasets | trajectories | median rho | negative fraction | score | target |", "|---|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for _, r in robust.sort_values(["tier", "drug_sensitivity_score"], ascending=[True, False]).head(50).iterrows():
        lines.append(f"| {r['tier']} | {r['module']} | {r['drug']} | {int(r['n_pairs'])} | {int(r['n_datasets'])} | {int(r['n_trajectories'])} | {r['median_rho_state_vs_lnIC50']:.3f} | {r['fraction_negative_rho']:.2f} | {r['drug_sensitivity_score']:.3f} | {r['putative_target']} |")
    lines += ["", "## Guardrails", "", "- 尚未加入 FDA/ChEMBL indication、CRC novelty、安全窗和暴露浓度；这些是下一轮 drug annotation，不在本轮臆测。", "- GDSC cell-line mapping uses COSMICID→DepMap ModelID; unmapped lines are excluded and recorded by the pairwise n. ", "- T2 candidates are subtype hypotheses. T1 candidates are the only ones eligible for a broad OXA-R follow-up screen; no candidate is promoted directly to wet experiment.", ""]
    (OUT / "phase8_module_conditioned_gdsc.md").write_text("\n".join(lines), encoding="utf-8")
    manifest = {"focus_modules": int(len(focus)), "pairwise_records": int(len(pairs)), "drug_module_records": int(len(ranking)), "robust_candidate_records": int(len(robust)), "mapping": "GDSC COSMICID -> DepMap ModelID", "rule": "T1=6/6 or T2=5/6 positive module directions; >=2 datasets; >=3 trajectories; fraction negative rho >=0.60"}
    (OUT / "phase8_module_conditioned_gdsc_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
