from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import phase1_pathway_matrix as matrix


OUT = matrix.ROOT / "phase3_module_decomposition"
OUT.mkdir(parents=True, exist_ok=True)


def build_decomposition_sets() -> tuple[dict[str, set[str]], dict[str, str]]:
    pathways, sources = matrix.build_pathways()
    h = matrix.HALLMARK
    r = matrix.REACTOME

    de_novo = {"CAD", "DHODH", "UMPS", "CTPS1", "CTPS2"}
    salvage = set(r["REACTOME_PYRIMIDINE_SALVAGE"]) | {"SLC28A1", "SLC28A2", "SLC28A3", "SLC29A1", "SLC29A2"}
    interconversion = {"CMPK1", "CMPK2", "CTPS1", "CTPS2", "DCTD", "DCTPP1", "DTYMK", "DUT", "NME1", "NME2", "NME3", "NME4", "NT5C2", "RRM1", "RRM2", "RRM2B", "TYMS", "TXNRD1"}
    catabolism = set(r["REACTOME_PYRIMIDINE_CATABOLISM"])

    sets: dict[str, set[str]] = {
        "pyrimidine_original_matrix": set(pathways["pyrimidine_metabolism"]),
        "pyrimidine_total_previous": set(pathways["pyrimidine_metabolism"]) | de_novo | interconversion,
        "pyrimidine_de_novo_core": de_novo,
        "pyrimidine_salvage": salvage,
        "pyrimidine_interconversion": interconversion,
        "pyrimidine_catabolism": catabolism,
        "UPR_total_hallmark": set(pathways["UPR_ER_stress"]),
        "PERK_eIF2A_ATF4": {"EIF2AK3", "EIF2S1", "ATF4", "DDIT3", "ATF3", "PPP1R15A", "HSPA5", "GADD34"},
        "IRE1_XBP1": {"ERN1", "XBP1", "EDEM1", "EDEM2", "EDEM3", "DNAJB9", "PDIA4", "MANF", "HSPA5", "DERL1", "SEL1L"},
        "ATF6_proteostasis": {"ATF6", "HSPA5", "XBP1", "HERPUD1", "EDEM1", "MANF", "PDIA4", "CALR", "CANX", "HSP90B1"},
        "EMT_total_hallmark": set(pathways["EMT"]),
        "epithelial_state": {"CDH1", "EPCAM", "OCLN", "TJP1", "CLDN3", "CLDN4", "CLDN7", "MUC1", "KRT8", "KRT18", "KRT19", "KRT20", "KRT7", "KRT17"},
        "mesenchymal_state": {"VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2", "CDH2", "ITGA5", "ITGB1", "MMP2", "MMP9", "MMP14", "SERPINE1", "TGFB1", "TGFBR1", "TGFBR2", "ACTA2", "COL1A1", "COL3A1"},
        "proliferation_E2F_G2M_MYC": set(h["HALLMARK_E2F_TARGETS"]) | set(h["HALLMARK_G2M_CHECKPOINT"]) | set(h["HALLMARK_MYC_TARGETS_V1"]) | set(h["HALLMARK_MYC_TARGETS_V2"]),
        "cell_cycle_G2M": set(h["HALLMARK_G2M_CHECKPOINT"]),
        "DNA_replication_E2F": set(h["HALLMARK_E2F_TARGETS"]),
    }
    sets["pyrimidine_original_no_prolif_overlap"] = sets["pyrimidine_original_matrix"] - sets["proliferation_E2F_G2M_MYC"]
    source_map = {name: ("MSigDB Hallmark" if "hallmark" in name or name.startswith("proliferation") or name.startswith("cell_cycle") or name.startswith("DNA_replication") else "curated mechanistic decomposition") for name in sets}
    source_map.update({"pyrimidine_salvage": "MSigDB Reactome + transport genes", "pyrimidine_catabolism": "MSigDB Reactome"})
    return sets, source_map


def score_contrast(expr: pd.DataFrame, contrast: dict, sets: dict[str, set[str]], sources: dict[str, str]) -> tuple[list[dict], list[dict]]:
    zexpr = matrix.zscore_rows(expr)
    parent = [x for x in contrast["parental"] if x in zexpr.columns]
    resistant = [x for x in contrast["resistant"] if x in zexpr.columns]
    module_rows: list[dict] = []
    gene_rows: list[dict] = []
    for name, genes in sets.items():
        present = sorted(set(genes) & set(zexpr.index))
        if not present or not parent or not resistant:
            delta = np.nan
        else:
            delta = float(zexpr.loc[present, resistant].mean(axis=0).mean() - zexpr.loc[present, parent].mean(axis=0).mean())
        module_rows.append({
            "contrast_id": f"{contrast['dataset']}|{contrast['model']}",
            "dataset": contrast["dataset"],
            "model": contrast["model"],
            "pathway": name,
            "source": sources[name],
            "n_genes_present": len(present),
            "delta": delta,
        })
        for gene in present:
            gd = float(zexpr.loc[gene, resistant].mean() - zexpr.loc[gene, parent].mean())
            gene_rows.append({
                "contrast_id": f"{contrast['dataset']}|{contrast['model']}",
                "dataset": contrast["dataset"],
                "model": contrast["model"],
                "module": name,
                "gene": gene,
                "delta_gene": gd,
            })
    return module_rows, gene_rows


def main() -> None:
    sets, sources = build_decomposition_sets()
    wanted = set().union(*sets.values())
    module_rows: list[dict] = []
    gene_rows: list[dict] = []
    cache: dict[tuple[str, str], pd.DataFrame] = {}
    for contrast in matrix.PRIMARY_CONTRASTS:
        key = (contrast["dataset"], contrast["platform"])
        if key not in cache:
            cache[key] = matrix.expression_for(contrast["dataset"], contrast["platform"], wanted)
        m, g = score_contrast(cache[key], contrast, sets, sources)
        module_rows.extend(m)
        gene_rows.extend(g)
    modules = pd.DataFrame(module_rows)
    genes = pd.DataFrame(gene_rows)
    modules.to_csv(OUT / "module_effects_primary.csv", index=False)
    genes.to_csv(OUT / "gene_effects_decomposition_primary.csv", index=False)

    cell_modules = modules[modules["model"] != "HCT116_xenograft"].copy()
    summary_rows = []
    for pathway, group in cell_modules.groupby("pathway"):
        vals = group["delta"].dropna()
        n = len(vals)
        up = int((vals > 0).sum())
        down = int((vals < 0).sum())
        summary_rows.append({
            "pathway": pathway,
            "n_models": n,
            "n_up": up,
            "n_down": down,
            "sign_concordance": max(up, down) / n if n else np.nan,
            "median_delta": vals.median() if n else np.nan,
            "mean_delta": vals.mean() if n else np.nan,
            "median_abs_delta": vals.abs().median() if n else np.nan,
            "median_n_genes_present": group["n_genes_present"].median(),
        })
    pd.DataFrame(summary_rows).sort_values(["sign_concordance", "median_abs_delta"], ascending=[False, False]).to_csv(OUT / "module_decomposition_summary_cell_lines.csv", index=False)

    gene_summary_rows = []
    cell_genes = genes[genes["model"] != "HCT116_xenograft"].copy()
    for (module, gene), group in cell_genes.groupby(["module", "gene"]):
        vals = group["delta_gene"].dropna()
        up = int((vals > 0).sum())
        down = int((vals < 0).sum())
        n = len(vals)
        gene_summary_rows.append({
            "module": module,
            "gene": gene,
            "n_models": n,
            "n_up": up,
            "n_down": down,
            "sign_concordance": max(up, down) / n if n else np.nan,
            "median_delta_gene": vals.median() if n else np.nan,
            "mean_delta_gene": vals.mean() if n else np.nan,
        })
    gene_summary = pd.DataFrame(gene_summary_rows).sort_values(["module", "sign_concordance", "median_delta_gene"], ascending=[True, False, False])
    gene_summary.to_csv(OUT / "gene_decomposition_summary_cell_lines.csv", index=False)

    wide = cell_modules.pivot(index="contrast_id", columns="pathway", values="delta")
    desired = ["EMT_total_hallmark", "pyrimidine_original_matrix", "UPR_total_hallmark"]
    state = wide[desired].copy()
    state.columns = ["EMT_delta", "pyrimidine_delta", "UPR_delta"]
    state["EMT_high"] = state["EMT_delta"] > 0
    state["pyrimidine_low"] = state["pyrimidine_delta"] < 0
    state["UPR_low"] = state["UPR_delta"] < 0
    state["all_three_directional"] = state[["EMT_high", "pyrimidine_low", "UPR_low"]].all(axis=1)
    state["composite_state_score"] = (
        (state["EMT_delta"] - state["EMT_delta"].mean()) / state["EMT_delta"].std(ddof=1)
        - (state["pyrimidine_delta"] - state["pyrimidine_delta"].mean()) / state["pyrimidine_delta"].std(ddof=1)
        - (state["UPR_delta"] - state["UPR_delta"].mean()) / state["UPR_delta"].std(ddof=1)
    )
    state.reset_index().to_csv(OUT / "emt_pyrimidine_upr_state_by_model.csv", index=False)
    corr = pd.DataFrame(index=desired, columns=desired, dtype=float)
    for a in desired:
        for b in desired:
            corr.loc[a, b] = spearmanr(wide[a], wide[b]).statistic
    corr.to_csv(OUT / "emt_pyrimidine_upr_spearman.csv")

    manifest = {
        "n_cell_line_contrasts": int(cell_modules["contrast_id"].nunique()),
        "decomposition_sets": {k: len(v) for k, v in sets.items()},
        "scoring": "within dataset gene-wise z-score; delta is resistant minus parental; xenograft excluded from stability summaries",
        "state_rule": "directional co-occurrence is EMT delta > 0, pyrimidine delta < 0, UPR delta < 0; exploratory only",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
