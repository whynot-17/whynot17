"""Figure 5 PPAR-definition and epithelial-state contradiction audit.

Primary inference is paired donor-level. Cell-level values are only aggregated
to donor/group/annotated-state units; cells are never treated as replicates.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
GENESETS = ROOT / "work" / "gene_sets"

PPAR_NR = ["PPARA", "PPARD", "PPARG", "NR1I2", "NR1I3", "NR1H2", "NR1H3"]
KEGG_PPAR = """AQP7B NR1H3 SORBS1 ME3 SLC27A5 SLC27A4 SLC27A2 APOA5 PLIN2 CPT1C CPT1A CPT1B CPT2 CYP7A1 CYP8B1 CYP27A1 DBI EHHADH FABP4 FABP1 FABP2 FABP3 FABP5 FABP6 FABP7 ACSL1 ACSL3 ACSL4 ACSBG1 ACSL6 GK GK2 GK3 SLC27A6 ACAA1 HMGCS1 HMGCS2 ACADL APOA1 APOA2 ACADM APOC3 ILK AQP7 SLC27A1 LPL ME1 MMP1 PLIN5 OLR1 ACOX1 PCK1 PCK2 ANGPTL4 PDPK1 ACSL5 PLIN1 PLTP PPARA PPARD PPARG RXRA RXRB RXRG SCD SCP2 PLIN4 UBB UBC UCP1 SCD5 ACSBG2 ACOX2 ACOX3 ADIPOQ FADS2 CD36""".split()


def read_gmt(path: Path) -> dict[str, set[str]]:
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            result[fields[0]] = {x.upper() for x in fields[2:] if x}
    return result


def bh(values: pd.Series) -> np.ndarray:
    vals = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(vals), np.nan)
    ok = np.isfinite(vals)
    if ok.any():
        out[ok] = multipletests(vals[ok], method="fdr_bh")[1]
    return out


def paired_summary(frame: pd.DataFrame, value: str, label: str) -> dict:
    pivot = frame.pivot_table(index="donor_id", columns="group", values=value, aggfunc="mean")
    delta = (pivot["tumor"] - pivot["normal"]).dropna()
    p = float(wilcoxon(delta).pvalue) if len(delta) and np.any(delta != 0) else 1.0
    med = float(delta.median())
    return {
        "definition": label,
        "n_paired_donors": len(delta),
        "mean_delta_tumor_minus_normal": float(delta.mean()),
        "median_delta_tumor_minus_normal": med,
        "p_value": p,
        "direction": "up" if med > 0 else "down" if med < 0 else "flat",
        "donor_consistency": float((delta > 0).mean() if med > 0 else (delta < 0).mean()),
    }


def main() -> None:
    pb = pd.read_csv(OUT / "mcop_phase2g_donor_state_pseudobulk.csv")
    meta = ["donor_key", "donor_id", "group", "PPAR_group", "n_cells"]
    genes = [c for c in pb.columns if c not in meta]
    donor = pb.groupby(["donor_key", "donor_id", "group"], as_index=False)[genes].sum()
    counts = donor[genes].astype(float)
    total = counts.sum(axis=1).replace(0, np.nan)
    logcpm = np.log1p(counts.div(total, axis=0) * 1_000_000)
    z = (logcpm - logcpm.mean()) / logcpm.std(ddof=1).replace(0, np.nan)

    reactome = read_gmt(GENESETS / "c2.cp.reactome.v2026.1.Hs.symbols.gmt")
    hallmark = read_gmt(GENESETS / "h.all.v2026.1.Hs.symbols.gmt")
    definitions = {
        "Frozen 7-gene PPAR/NR core": set(PPAR_NR),
        "PPAR receptor-only module": {"PPARA", "PPARD", "PPARG"},
        "NR partner module": {"NR1I2", "NR1I3", "NR1H2", "NR1H3"},
        "KEGG hsa03320 PPAR signaling": set(KEGG_PPAR),
        "Reactome nuclear receptor transcription": reactome["REACTOME_NUCLEAR_RECEPTOR_TRANSCRIPTION_PATHWAY"],
        "Reactome PPAR-alpha lipid regulation": reactome["REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA"],
        "Reactome peroxisomal lipid metabolism": reactome["REACTOME_PEROXISOMAL_LIPID_METABOLISM"],
        "Reactome mitochondrial fatty-acid beta oxidation": reactome["REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION"],
        "Hallmark fatty acid metabolism": hallmark["HALLMARK_FATTY_ACID_METABOLISM"],
        "Hallmark cholesterol homeostasis": hallmark["HALLMARK_CHOLESTEROL_HOMEOSTASIS"],
        "Enterocyte metabolic differentiation": {"ALPI", "FABP1", "FABP2", "APOA1", "APOA4", "SI", "SLC26A3", "CA1", "CA2", "VIL1"},
    }
    sources = {
        "Frozen 7-gene PPAR/NR core": "prespecified Phase 2F core",
        "PPAR receptor-only module": "prespecified receptor decomposition",
        "NR partner module": "prespecified receptor decomposition",
        "KEGG hsa03320 PPAR signaling": "KEGG hsa03320; retrieved 2026-08-24",
        "Reactome nuclear receptor transcription": "MSigDB Reactome v2026.1",
        "Reactome PPAR-alpha lipid regulation": "MSigDB Reactome v2026.1",
        "Reactome peroxisomal lipid metabolism": "MSigDB Reactome v2026.1",
        "Reactome mitochondrial fatty-acid beta oxidation": "MSigDB Reactome v2026.1",
        "Hallmark fatty acid metabolism": "MSigDB Hallmark v2026.1",
        "Hallmark cholesterol homeostasis": "MSigDB Hallmark v2026.1",
        "Enterocyte metabolic differentiation": "prespecified Phase 2G marker panel",
    }
    scored = donor[["donor_key", "donor_id", "group"]].copy()
    frozen_source = pd.read_csv(OUT / "mcop_phase2f_singlecell_donor_scores.csv")
    frozen_source = frozen_source.loc[
        frozen_source["compartment"].eq("epithelial")
        & frozen_source["group"].isin(["normal", "tumor"]),
        ["donor_key", "group", "PPAR_nuclear_receptor_score"],
    ].drop_duplicates(["donor_key", "group"])
    scored = scored.merge(frozen_source, on=["donor_key", "group"], how="left", validate="one_to_one")
    scored = scored.rename(columns={"PPAR_nuclear_receptor_score": "Frozen 7-gene PPAR/NR core"})
    if scored["Frozen 7-gene PPAR/NR core"].isna().any():
        raise RuntimeError("Frozen Phase 2F donor scores did not map to all 36 paired donor/group units.")
    rows = []
    for label, gene_set in definitions.items():
        present = sorted(gene_set & set(genes))
        if label != "Frozen 7-gene PPAR/NR core":
            scored[label] = z[present].mean(axis=1)
        row = paired_summary(scored, label, label)
        row.update({
            "source_gene_count": len(gene_set),
            "genes_present": len(present),
            "coverage": len(present) / len(gene_set),
            "gene_symbols_present": ";".join(present),
            "scoring_method": "frozen Phase 2F targeted 9-gene-denominator score" if label == "Frozen 7-gene PPAR/NR core" else "mean gene-wise z score of donor pseudobulk target-universe logCPM; 36 paired donors",
            "score_name": label, "score_source": sources[label], "n_genes": len(gene_set),
        })
        rows.append(row)

    existing_reg = pd.read_csv(OUT / "mcop_phase2g_regulator_activity.csv")
    for regulator in ["PPARA", "PPARD", "PPARG"]:
        hit = existing_reg.loc[(existing_reg["comparison"] == "tumor_vs_normal") & (existing_reg["regulator"] == regulator)]
        if hit.empty:
            rows.append({
                "definition": f"DoRothEA {regulator} regulon activity", "n_paired_donors": 0,
                "direction": "not_estimable", "scoring_method": "DoRothEA A-C + decoupler ULM; regulator absent/not estimable",
            })
            continue
        r = hit.iloc[0]
        rows.append({
            "definition": f"DoRothEA {regulator} regulon activity", "n_paired_donors": int(r.n_pairs),
            "mean_delta_tumor_minus_normal": float(r.activity_delta),
            "median_delta_tumor_minus_normal": float(r.median_activity_delta),
            "p_value": float(r.P), "direction": "up" if r.activity_delta > 0 else "down",
            "donor_consistency": float(r.donor_consistency), "source_gene_count": np.nan,
            "genes_present": np.nan, "coverage": np.nan, "gene_symbols_present": "",
            "scoring_method": "DoRothEA levels A-C + decoupler ULM; paired donor pseudobulk",
            "score_name": f"DoRothEA {regulator} regulon activity", "score_source": "DoRothEA A-C", "n_genes": np.nan,
        })
    comparison = pd.DataFrame(rows)
    comparison["BH_FDR"] = bh(comparison["p_value"])
    comparison["median_delta"] = comparison["median_delta_tumor_minus_normal"]
    comparison["mean_delta"] = comparison["mean_delta_tumor_minus_normal"]
    comparison["wilcoxon_P"] = comparison["p_value"]
    comparison.to_csv(OUT / "figure5_ppar_definition_comparison.csv", index=False)

    # Gene-level decomposition across the union of all pathway definitions.
    union = sorted(set().union(*definitions.values()) & set(genes))
    gene_rows = []
    for gene in union:
        temp = donor[["donor_id", "group"]].copy()
        temp["value"] = logcpm[gene]
        row = paired_summary(temp, "value", gene)
        row["gene"] = gene
        row["in_frozen_core"] = gene in PPAR_NR
        row["in_kegg_ppar"] = gene in KEGG_PPAR
        row["in_reactome_pparalpha"] = gene in definitions["Reactome PPAR-alpha lipid regulation"]
        row["in_peroxisomal_lipid"] = gene in definitions["Reactome peroxisomal lipid metabolism"]
        memberships = [label for label, gene_set in definitions.items() if gene in gene_set]
        if gene in {"PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG", "NR1H2", "NR1H3", "NR1I2", "NR1I3"}:
            module = "receptor_nuclear_receptor"
        elif gene == "CD36" or gene.startswith("FABP") or gene.startswith("SLC27"):
            module = "fatty_acid_uptake"
        elif gene in {"CPT1A", "CPT1B", "CPT2", "ACOX1", "ACOX2", "ACOX3", "ACADM", "ACADL", "HADHA", "HADHB", "EHHADH"}:
            module = "fatty_acid_oxidation"
        elif gene in definitions["Reactome peroxisomal lipid metabolism"]:
            module = "peroxisome"
        elif gene in {"SCD", "SCD5", "FADS2", "HMGCS1", "HMGCS2", "ME1", "PCK1", "PCK2"}:
            module = "lipid_synthesis_or_handling"
        else:
            module = "other_ppar_pathway"
        row["PPAR_family"] = "PPAR_receptor" if gene in {"PPARA", "PPARD", "PPARG"} else "PPAR_related"
        row["pathway_membership"] = ";".join(memberships)
        row["module"] = module
        row["median_delta"] = row["median_delta_tumor_minus_normal"]
        row["P"] = row["p_value"]
        gene_rows.append(row)
    gene_decomp = pd.DataFrame(gene_rows)
    gene_decomp["BH_FDR"] = bh(gene_decomp["p_value"])
    gene_decomp.to_csv(OUT / "figure5_ppar_gene_module_decomposition.csv", index=False)

    # Composition and within-annotation donor-level analyses.
    cells = pd.read_csv(OUT / "mcop_phase2g_epithelial_state_scores.csv")
    comp = cells.groupby(["donor_id", "group", "cell_subtype"], as_index=False).size().rename(columns={"size": "n_cells"})
    totals = comp.groupby(["donor_id", "group"])["n_cells"].transform("sum")
    comp["fraction"] = comp["n_cells"] / totals
    all_units = pd.MultiIndex.from_product([
        sorted(cells.donor_id.astype(str).unique()), ["normal", "tumor"], sorted(cells.cell_subtype.unique())
    ], names=["donor_id", "group", "cell_subtype"]).to_frame(index=False)
    comp = all_units.merge(comp, on=["donor_id", "group", "cell_subtype"], how="left").fillna({"n_cells": 0, "fraction": 0})
    comp_rows = []
    for state, sub in comp.groupby("cell_subtype"):
        row = paired_summary(sub, "fraction", state)
        row["cell_subtype"] = state
        row["normal_median_fraction"] = sub.loc[sub.group == "normal", "fraction"].median()
        row["tumor_median_fraction"] = sub.loc[sub.group == "tumor", "fraction"].median()
        comp_rows.append(row)
    comp_summary = pd.DataFrame(comp_rows)
    comp_summary["BH_FDR"] = bh(comp_summary["p_value"])
    comp_out = comp.rename(columns={"group": "condition", "cell_subtype": "state", "n_cells": "cell_count", "fraction": "state_fraction"})
    comp_out = comp_out.merge(
        comp_summary.rename(columns={"cell_subtype": "state", "median_delta_tumor_minus_normal": "paired_median_delta", "p_value": "paired_wilcoxon_P", "BH_FDR": "paired_BH_FDR", "direction": "paired_direction"})[
            ["state", "n_paired_donors", "paired_median_delta", "paired_wilcoxon_P", "paired_BH_FDR", "paired_direction"]
        ], on="state", how="left", validate="many_to_one",
    )
    comp_out.to_csv(OUT / "figure5_epithelial_state_composition.csv", index=False)

    within = cells.groupby(["donor_id", "group", "cell_subtype"], as_index=False).agg(
        n_cells=("cell_id", "size"), PPAR_NR_score=("PPAR_NR_score", "mean")
    )
    within_rows = []
    for state, sub in within.groupby("cell_subtype"):
        eligible = sub.pivot_table(index="donor_id", columns="group", values="n_cells", aggfunc="sum").fillna(0)
        keep = eligible.index[(eligible.get("normal", 0) >= 20) & (eligible.get("tumor", 0) >= 20)]
        test = sub.loc[sub.donor_id.isin(keep)]
        row = paired_summary(test, "PPAR_NR_score", state)
        row.update({"cell_subtype": state, "minimum_cells_per_condition": 20})
        within_rows.append(row)
    within_out = pd.DataFrame(within_rows)
    within_out["BH_FDR"] = bh(within_out["p_value"])
    within_out.to_csv(OUT / "figure5_within_state_ppar_analysis.csv", index=False)

    # Paired donor deltas and cross-program correlations.
    phase_states = pd.read_csv(OUT / "mcop_phase2g_donor_level_validation.csv")
    state_names = [
        "TNF_NFkB", "IL6_JAK_STAT3", "Inflammatory_response", "stress_like_epithelial",
        "EMT", "E2F_targets", "G2M_checkpoint", "MYC_targets_V1", "Hypoxia",
        "enterocyte_differentiation", "intestinal_epithelial_differentiation",
        "Fatty_acid_metabolism", "OXPHOS", "UPR",
    ]
    state_rows = phase_states.loc[(phase_states.feature_type == "state") & phase_states.feature.isin(state_names)].copy()
    # Recompute donor-level program values from the same matrix so correlations use 36 paired deltas.
    custom = {
        "enterocyte_differentiation": {"ALPI", "FABP1", "FABP2", "APOA1", "APOA4", "SI", "SLC26A3", "CA1", "CA2", "VIL1"},
        "intestinal_epithelial_differentiation": {"CDX2", "KLF4", "HNF4A", "GATA6", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "MUC13", "ALPI", "FABP1", "FABP2", "SLC26A3", "CA1", "CA2", "SI"},
        "stress_like_epithelial": {"HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "ATF3", "DDIT3", "XBP1", "JUN", "FOS", "HIF1A"},
    }
    hmap = {
        "TNF_NFkB": "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "IL6_JAK_STAT3": "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "Inflammatory_response": "HALLMARK_INFLAMMATORY_RESPONSE", "EMT": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "E2F_targets": "HALLMARK_E2F_TARGETS", "G2M_checkpoint": "HALLMARK_G2M_CHECKPOINT",
        "MYC_targets_V1": "HALLMARK_MYC_TARGETS_V1", "Hypoxia": "HALLMARK_HYPOXIA",
        "Fatty_acid_metabolism": "HALLMARK_FATTY_ACID_METABOLISM", "OXPHOS": "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "UPR": "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    }
    programs = {k: hallmark[v] for k, v in hmap.items()} | custom
    delta = pd.DataFrame(index=sorted(donor.donor_id.unique()))
    for label in ["Frozen 7-gene PPAR/NR core"]:
        piv = scored.pivot_table(index="donor_id", columns="group", values=label)
        delta["PPAR_NR"] = piv.tumor - piv.normal
    for label, gs in programs.items():
        present = sorted(gs & set(genes))
        values = donor[["donor_id", "group"]].copy()
        values["score"] = z[present].mean(axis=1)
        piv = values.pivot_table(index="donor_id", columns="group", values="score")
        delta[label] = piv.tumor - piv.normal
    # Rebuild donor-level RELA/STAT3 regulator activities, not gene-expression proxies.
    try:
        import decoupler as dc
        net = dc.op.dorothea(organism="human", levels=["A", "B", "C"], license="academic", verbose=False)
        net = net.loc[net["source"].isin(["RELA", "STAT3"]) & net["target"].isin(logcpm.columns)].copy()
        matrix = logcpm.copy()
        matrix.index = donor["donor_key"].astype(str) + "|" + donor["group"].astype(str)
        activity, _ = dc.mt.ulm(matrix, net, tmin=5, verbose=False)
        activity = activity.reset_index().rename(columns={"index": "unit_id"})
        activity["donor_id"] = donor["donor_id"].to_numpy()
        activity["group"] = donor["group"].to_numpy()
        for regulator in ["RELA", "STAT3"]:
            if regulator in activity.columns:
                piv = activity.pivot_table(index="donor_id", columns="group", values=regulator)
                delta[f"{regulator}_activity"] = piv.tumor - piv.normal
    except Exception as exc:
        raise RuntimeError(f"DoRothEA donor-level RELA/STAT3 activity reconstruction failed: {exc}") from exc
    corr_rows = []
    for target in [c for c in delta.columns if c != "PPAR_NR"]:
        pair = delta[["PPAR_NR", target]].dropna()
        rho, p = spearmanr(pair.PPAR_NR, pair[target])
        corr_rows.append({"x": "delta_PPAR_NR", "y": f"delta_{target}", "n_donors": len(pair), "spearman_rho": rho, "p_value": p})
    corr = pd.DataFrame(corr_rows)
    corr["BH_FDR"] = bh(corr.p_value)
    corr.to_csv(OUT / "figure5_state_correlation_matrix.csv", index=False)

    print(comparison[["definition", "median_delta_tumor_minus_normal", "p_value", "BH_FDR", "direction"]].to_string(index=False))


if __name__ == "__main__":
    main()
