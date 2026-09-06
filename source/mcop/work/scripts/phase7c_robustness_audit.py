"""Phase 7C robustness audit.

Uses the already-computed trajectory x module NES table.  No new biological
claims are made here: the script tests leave-one-trajectory and leave-one-GEO-
dataset stability, then collapses highly overlapping gene-set representations
so one biology is not counted several times.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
TRAJECTORIES = [
    "GSE77932|HCT116", "GSE77932|DLD1", "GSE42387|HCT116",
    "GSE42387|HT29", "GSE42387|LoVo", "GSE119603|HCT116",
]
DATASETS = {x: x.split("|", 1)[0] for x in TRAJECTORIES}


def sign_consistency(values: pd.Series) -> tuple[int, int, float, str]:
    x = pd.to_numeric(values, errors="coerce").dropna()
    pos = int((x > 0).sum())
    neg = int((x < 0).sum())
    if len(x) == 0:
        return pos, neg, np.nan, "none"
    direction = "positive" if pos >= neg else "negative"
    return pos, neg, max(pos, neg) / len(x), direction


def main() -> None:
    gsea = pd.read_csv(OUT / "phase7c_trajectory_module_gsea.csv")
    nes = gsea.pivot_table(index="module", columns="trajectory", values="nes", aggfunc="first").reindex(columns=TRAJECTORIES)
    conv = pd.read_csv(OUT / "phase7c_module_convergence_ranking.csv")
    conv = conv.set_index("module")
    rows = []
    for module, row in nes.iterrows():
        vals = pd.to_numeric(row, errors="coerce")
        full_pos, full_neg, full_cons, full_dir = sign_consistency(vals)
        record = {
            "module": module,
            "collection": str(module).split("::", 1)[0],
            "full_direction": full_dir,
            "full_consistency": full_cons,
            "full_median_nes": float(vals.median()),
            "full_n_fdr_qval_025": int(conv.loc[module, "n_fdr_qval_025"]) if module in conv.index else 0,
        }
        loo_cons = []
        for omit in TRAJECTORIES:
            rest = vals.drop(labels=omit)
            p, n, c, d = sign_consistency(rest)
            loo_cons.append(c)
            record[f"trajectory_loo_{omit}_consistency"] = c
            record[f"trajectory_loo_{omit}_direction"] = d
        record["trajectory_loo_min_consistency"] = float(np.nanmin(loo_cons))
        record["trajectory_loo_all_same_direction"] = bool(all(c == 1.0 for c in loo_cons))
        ds_cons = []
        for dataset in sorted(set(DATASETS.values())):
            keep = [t for t in TRAJECTORIES if DATASETS[t] != dataset]
            rest = vals.reindex(keep)
            p, n, c, d = sign_consistency(rest)
            ds_cons.append(c)
            record[f"dataset_loo_{dataset}_n"] = len(rest.dropna())
            record[f"dataset_loo_{dataset}_consistency"] = c
            record[f"dataset_loo_{dataset}_direction"] = d
        record["dataset_loo_min_consistency"] = float(np.nanmin(ds_cons))
        record["dataset_loo_all_same_direction"] = bool(all(c == 1.0 for c in ds_cons))
        record["robust_positive_or_negative"] = bool(record["trajectory_loo_all_same_direction"] and record["dataset_loo_all_same_direction"] and abs(record["full_median_nes"]) >= 0.5)
        rows.append(record)
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "phase7c_robustness_by_module.csv", index=False)

    # Select modules that have both a meaningful full effect and at least two
    # trajectory-level FDR<=0.25 hits.  The overlap collapse is restricted to
    # this pre-specified evidence tier to avoid clustering thousands of noise
    # gene sets.
    selected = audit[(audit["full_n_fdr_qval_025"] >= 2) & (audit["full_median_nes"].abs() >= 0.75)].copy()
    sys.path.insert(0, str(ROOT / "work" / "scripts"))
    import phase7c_functional_module_convergence as p7c
    all_sets, _ = p7c.load_gene_sets()

    selected_modules = [m for m in selected["module"] if m in all_sets]
    parent = {m: m for m in selected_modules}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    # Overlap coefficient is preferable to Jaccard for a small CORUM complex
    # embedded in a larger Reactome/co-essentiality module.
    for i, a in enumerate(selected_modules):
        ga = all_sets[a]
        for b in selected_modules[i + 1:]:
            gb = all_sets[b]
            overlap = len(ga & gb) / max(1, min(len(ga), len(gb)))
            if overlap >= 0.60:
                union(a, b)
    clusters: dict[str, list[str]] = {}
    for m in selected_modules:
        clusters.setdefault(find(m), []).append(m)
    cluster_rows = []
    for cid, members in enumerate(clusters.values(), 1):
        sub = selected[selected["module"].isin(members)].copy()
        sub["priority_score"] = sub["full_n_fdr_qval_025"] + 2 * sub["full_consistency"] + sub["full_median_nes"].abs()
        rep = sub.sort_values("priority_score", ascending=False).iloc[0]
        for m in members:
            cluster_rows.append({
                "redundancy_cluster": f"RC{cid:03d}",
                "representative": rep["module"],
                "module": m,
                "collection": m.split("::", 1)[0],
                "n_genes": len(all_sets[m]),
                "full_direction": audit.loc[audit["module"].eq(m), "full_direction"].iloc[0],
                "full_consistency": audit.loc[audit["module"].eq(m), "full_consistency"].iloc[0],
                "full_median_nes": audit.loc[audit["module"].eq(m), "full_median_nes"].iloc[0],
                "full_n_fdr_qval_025": audit.loc[audit["module"].eq(m), "full_n_fdr_qval_025"].iloc[0],
            })
    clusters_df = pd.DataFrame(cluster_rows).sort_values(["redundancy_cluster", "representative", "module"])
    clusters_df.to_csv(OUT / "phase7c_candidate_redundancy_clusters.csv", index=False)

    robust = audit[audit["robust_positive_or_negative"]].copy().sort_values(["full_n_fdr_qval_025", "full_consistency", "full_median_nes"], ascending=[False, False, False])
    lines = [
        "# Phase 7C robustness audit", "",
        "本轮只审计 Phase 7C 已有 NES，不重新筛药。", "",
        "## 规则", "",
        "- trajectory leave-one-out：去掉任意一条 trajectory 后，模块方向仍需完全一致。", "- GEO dataset leave-one-out：分别去掉 GSE77932、GSE42387 或 GSE119603 后，剩余 trajectory 方向仍需完全一致。", "- 重点模块同时要求完整分析中 |median NES| ≥ 0.75 且至少 2/6 trajectory 的 FDR q≤0.25。", "- 高重叠 gene sets 用 overlap coefficient ≥0.60 聚类，防止同一生物学被 Reactome/CORUM/co-essentiality 重复计数。", "",
        f"重点模块数：{len(selected)}；冗余簇数：{len(clusters)}；同时通过两类 leave-out 的模块数：{len(robust)}。", "",
        "## 同时通过两类 leave-out 的模块", "", "| Module | Collection | Direction | full consistency | median NES | FDR≤0.25 hits |", "|---|---|---|---:|---:|---:|"]
    for _, r in robust.head(40).iterrows():
        lines.append(f"| {r['module']} | {r['collection']} | {r['full_direction']} | {r['full_consistency']:.2f} | {r['full_median_nes']:.2f} | {int(r['full_n_fdr_qval_025'])} |")
    lines += ["", "## 解释", "", "- 通过 leave-out 只说明信号不是由单一 trajectory 或单一 GEO 数据集完全驱动；仍不能证明 acquired OXA-R causality。", "- 通过 redundancy collapse 后才把模块簇作为下一阶段药物映射单位。", "- 若某模块只在完整数据中显著、但 leave-out 失败，应降为 exploratory，不进入主候选。", ""]
    (OUT / "phase7c_robustness_audit.md").write_text("\n".join(lines), encoding="utf-8")
    manifest = {"n_modules": int(len(audit)), "n_selected_for_redundancy": int(len(selected)), "n_redundancy_clusters": int(len(clusters)), "n_robust_modules": int(len(robust)), "rules": "trajectory and dataset leave-one-out; overlap coefficient >=0.60; full |median NES| >=0.75; >=2 trajectory FDR q<=0.25"}
    (OUT / "phase7c_robustness_audit_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
