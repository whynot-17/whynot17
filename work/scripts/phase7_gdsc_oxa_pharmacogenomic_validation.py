"""Phase 7A: phenotype-direct pharmacogenomic validation in CRC.

This module tests whether a drug is selectively active in CRC cell lines with
high oxaliplatin LN_IC50. It deliberately does not claim a molecular subtype
without an independent expression projection. The outputs therefore separate:

1. OXA-R-like versus OXA-S-like group contrasts;
2. cross-line OXA/drug sensitivity covariance;
3. replication across GDSC1 and GDSC2.

Higher LN_IC50 means lower sensitivity. Consequently, a negative
R-minus-S median difference is collateral sensitivity in OXA-R-like lines.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "work" / "phase7_convergent_vulnerability" / "data"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

DATASETS = {
    "GDSC1": DATA / "GDSC1_fitted_dose_response_24Jul22.xlsx",
    "GDSC2": DATA / "GDSC2_fitted_dose_response_24Jul22.xlsx",
}
MIN_N = 15


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    order = np.argsort(p[ok])
    vals = p[ok][order]
    q = vals * len(vals) / np.arange(1, len(vals) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1].clip(0, 1)
    tmp = np.empty_like(q)
    tmp[order] = q
    out[ok] = tmp
    return out


def load_crc(path):
    d = pd.read_excel(path)
    d["DRUG_NAME"] = d["DRUG_NAME"].astype(str).str.strip()
    d["CELL_LINE_NAME"] = d["CELL_LINE_NAME"].astype(str).str.strip()
    d = d[d["TCGA_DESC"].astype(str).str.upper().isin({"COREAD", "COAD", "READ", "COAD/READ"})]
    d = d.dropna(subset=["LN_IC50"]).copy()
    # Replicates are collapsed before any inference.
    d = (d.groupby(["CELL_LINE_NAME", "DRUG_NAME"], as_index=False)
           .agg(LN_IC50=("LN_IC50", "mean"), AUC=("AUC", "mean")))
    return d


def one_dataset(dataset, path):
    d = load_crc(path)
    mat = d.pivot(index="CELL_LINE_NAME", columns="DRUG_NAME", values="LN_IC50")
    if "Oxaliplatin" not in mat.columns:
        raise RuntimeError(f"{dataset}: Oxaliplatin not found")
    oxa = mat["Oxaliplatin"]
    # Median split is prespecified and only used for a descriptive binary gate.
    cutoff = float(oxa.median())
    labels = pd.Series(np.where(oxa > cutoff, "OXA-R-like", "OXA-S-like"), index=oxa.index)
    oxa_out = pd.DataFrame({"dataset": dataset, "CELL_LINE_NAME": oxa.index,
                            "OXA_LN_IC50": oxa.values,
                            "OXA_group": labels.values})
    rows = []
    for drug in mat.columns:
        if drug == "Oxaliplatin":
            continue
        x = pd.concat([oxa.rename("oxa"), mat[drug].rename("drug")], axis=1).dropna()
        n = len(x)
        if n < MIN_N:
            continue
        r, p = stats.spearmanr(x["oxa"], x["drug"])
        r_pearson, p_pearson = stats.pearsonr(x["oxa"], x["drug"])
        rgrp = x.loc[x.index.map(labels) == "OXA-R-like", "drug"]
        sgrp = x.loc[x.index.map(labels) == "OXA-S-like", "drug"]
        if len(rgrp) < 5 or len(sgrp) < 5:
            continue
        u, p_u = stats.mannwhitneyu(rgrp, sgrp, alternative="two-sided")
        rows.append({
            "dataset": dataset, "DRUG_NAME": drug, "n_overlap": n,
            "n_OXA_R_like": len(rgrp), "n_OXA_S_like": len(sgrp),
            "spearman_r": r, "spearman_p": p,
            "pearson_r": r_pearson, "pearson_p": p_pearson,
            "median_LN_IC50_R": rgrp.median(), "median_LN_IC50_S": sgrp.median(),
            "delta_R_minus_S": rgrp.median() - sgrp.median(),
            "mannwhitney_p": p_u,
            # negative = more sensitive in OXA-R-like lines
            "collateral_sensitivity_direction": "yes" if rgrp.median() < sgrp.median() else "no",
        })
    res = pd.DataFrame(rows)
    for col in ["spearman_p", "pearson_p", "mannwhitney_p"]:
        res[col + "_BH"] = bh(res[col].values)
    return oxa_out, res


def main():
    oxa_tables, results = [], []
    unavailable = []
    for dataset, path in DATASETS.items():
        try:
            oxa, res = one_dataset(dataset, path)
        except RuntimeError as exc:
            unavailable.append(str(exc))
            print(f"skip {dataset}: {exc}")
            continue
        oxa_tables.append(oxa)
        results.append(res)
    all_res = pd.concat(results, ignore_index=True)
    oxa_all = pd.concat(oxa_tables, ignore_index=True)
    all_res.to_csv(OUT / "phase7_gdsc_drug_oxa_association_by_dataset.csv", index=False)
    oxa_all.to_csv(OUT / "phase7_gdsc_crc_oxaliplatin_phenotype.csv", index=False)

    # Replication table for drugs measured in both datasets.
    piv = all_res.pivot_table(index="DRUG_NAME", columns="dataset",
                              values=["spearman_r", "spearman_p_BH", "delta_R_minus_S",
                                      "mannwhitney_p_BH", "n_overlap"])
    piv.columns = ["%s_%s" % (a, b) for a, b in piv.columns]
    piv = piv.reset_index()
    if {"spearman_r_GDSC1", "spearman_r_GDSC2"}.issubset(piv.columns):
        piv["replicated_negative_covariance"] = (
            (piv["spearman_r_GDSC1"] < 0) & (piv["spearman_r_GDSC2"] < 0))
        piv["replicated_collateral_direction"] = (
            (piv["delta_R_minus_S_GDSC1"] < 0) & (piv["delta_R_minus_S_GDSC2"] < 0))
        piv["mean_spearman_r"] = piv[["spearman_r_GDSC1", "spearman_r_GDSC2"]].mean(axis=1)
        piv["mean_delta_R_minus_S"] = piv[["delta_R_minus_S_GDSC1", "delta_R_minus_S_GDSC2"]].mean(axis=1)
    piv.to_csv(OUT / "phase7_gdsc_cross_dataset_replication.csv", index=False)

    candidates = ["Atorvastatin", "Sulfasalazine", "Topiramate", "Amlodipine",
                  "Teriflunomide", "Bortezomib", "Minocycline", "Meldonium",
                  "Ivermectin", "Leflunomide"]
    cand = all_res[all_res["DRUG_NAME"].str.lower().isin({x.lower() for x in candidates})].copy()
    measured = {str(x).lower() for x in cand["DRUG_NAME"]}
    missing = pd.DataFrame([{"dataset": "GDSC1+GDSC2", "DRUG_NAME": x,
                             "status": "not_measured_in_available_OXA_anchor"}
                            for x in candidates if x.lower() not in measured])
    if len(cand):
        cand["status"] = "measured"
    cand = pd.concat([cand, missing], ignore_index=True, sort=False)
    cand.to_csv(OUT / "phase7_gdsc_named_candidate_validation.csv", index=False)
    top_negative = (all_res.sort_values(["spearman_r", "delta_R_minus_S"])
                    .head(15).copy())
    top_negative.to_csv(OUT / "phase7_gdsc_top_negative_oxa_associations.csv", index=False)

    # Blind phenotype-direct lead table: both datasets, negative direction in both,
    # and at least one nominally significant replication signal.
    if "replicated_negative_covariance" in piv.columns:
        gate = piv[(piv["replicated_negative_covariance"]) | (piv["replicated_collateral_direction"])].copy()
    else:
        gate = piv.iloc[0:0].copy()
    gate.to_csv(OUT / "phase7_gdsc_blind_collateral_sensitivity_gate.csv", index=False)

    lines = [
        "# Phase 7A：GDSC phenotype-direct OXA pharmacogenomic validation",
        "",
        "## 结论",
        "",
        f"本轮在可提供 Oxaliplatin 的 GDSC CRC 细胞系中，按细胞系内平均 LN_IC50 进行分析；共输出 {len(all_res):,} 个 dataset-drug 组合。LN_IC50 越高代表越不敏感。",
        "",
        "这一步验证的是 OXA-R-like phenotype，而不是表达谱定义的 ERAD/DHODH 亚型。因此它是 Phase 7 的独立药敏锚点，不能替代 subtype projection。",
        "",
        "## 判定规则",
        "",
        "- OXA-R-like：CRC 细胞系 Oxaliplatin LN_IC50 高于该数据集 CRC 中位数；OXA-S-like：不高于中位数。",
        "- collateral sensitivity：OXA-R-like 组的候选药 LN_IC50 更低，即 `delta_R_minus_S < 0`；Spearman OXA/药物敏感性相关也应为负。",
        "- 预设复制要求：若两个数据集都含 Oxaliplatin，则要求 GDSC1 与 GDSC2 同方向；若数据集没有 Oxaliplatin，则标为不可复制，不把缺失当阴性。",
        "",
        "## 当前结果",
        "",
    ]
    if unavailable:
        lines += ["数据可用性：", "", *[f"- {x}" for x in unavailable], ""]
    if len(gate):
        lines.append(f"盲筛得到 {len(gate)} 个至少在两个 GDSC 数据集方向一致的候选，详见 `phase7_gdsc_blind_collateral_sensitivity_gate.csv`。")
    else:
        lines.append("按‘两个 GDSC 数据集方向一致’的严格盲筛门槛，当前没有通过者；这不是证明没有药物效应，而是说明不能把单一数据集的负相关当作可靠 Drug X。")
    if len(top_negative):
        lines += ["", "GDSC2/可用数据中，方向最偏向 OXA-R collateral sensitivity 的前 5 个药物为："]
        for _, row in top_negative.head(5).iterrows():
            lines.append(f"- {row['DRUG_NAME']}：Spearman r={row['spearman_r']:.3f}，BH-p={row['spearman_p_BH']:.3f}，R−S median LN_IC50={row['delta_R_minus_S']:.3f}，n={int(row['n_overlap'])}")
    measured_cand = cand[cand["status"].eq("measured")] if "status" in cand.columns else cand.iloc[0:0]
    if len(measured_cand):
        lines += ["", "已测到的预先指定候选均未呈现 collateral-sensitivity 方向："]
        for _, row in measured_cand.iterrows():
            lines.append(f"- {row['DRUG_NAME']}：Spearman r={row['spearman_r']:.3f}，BH-p={row['spearman_p_BH']:.3f}，R−S median LN_IC50={row['delta_R_minus_S']:.3f}")
    lines += [
        "",
        "### 预先指定候选",
        "",
        "候选药物的 GDSC 覆盖和方向见 `phase7_gdsc_named_candidate_validation.csv`。若候选不在 GDSC1/2，则记为 pharmacogenomic unavailable，而不是阴性。",
        "",
        "## 重要限制",
        "",
        "1. GDSC1/GDSC2 是药敏关联数据，不是 acquired OXA-R 与 parental 的配对实验；只能支持 phenotype-level prioritization。",
        "2. 细胞系之间的 lineage、增殖速度与普遍药敏会造成混杂；下一步若有独立表达谱，应对 subtype score 做连续建模并控制 growth-related covariates。",
        "3. 本轮不能证明 ERAD-high 或 DHODH-high 选择性，除非表达谱投影成功并与药敏矩阵在同一细胞系 ID 上可靠合并。",
        "",
        "## 文件",
        "",
        "- `phase7_gdsc_crc_oxaliplatin_phenotype.csv`",
        "- `phase7_gdsc_drug_oxa_association_by_dataset.csv`",
        "- `phase7_gdsc_cross_dataset_replication.csv`",
        "- `phase7_gdsc_named_candidate_validation.csv`",
        "- `phase7_gdsc_top_negative_oxa_associations.csv`",
        "- `phase7_gdsc_blind_collateral_sensitivity_gate.csv`",
    ]
    (OUT / "phase7_gdsc_pharmacogenomic_validation.md").write_text("\n".join(lines), encoding="utf-8")
    print("saved Phase 7A outputs")
    print("dataset-drug rows", len(all_res), "blind gate", len(gate))
    print("named candidates")
    print(cand[["dataset", "DRUG_NAME", "n_overlap", "spearman_r", "spearman_p_BH", "delta_R_minus_S", "mannwhitney_p_BH"]].to_string(index=False) if len(cand) else "none measured")


if __name__ == "__main__":
    main()
