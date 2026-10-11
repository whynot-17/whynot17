#!/usr/bin/env python
"""Build the concise report and cross-cohort recurrent-interaction artifacts."""

from pathlib import Path
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("D:/DINP/whynot17-new-project")
OUT = ROOT / "results/2026-10-11/cellchat_patient_consistency"


def main() -> None:
    cross = pd.read_csv(OUT / "CellChat_cross_cohort_sensitivity.tsv.gz", sep="\t")
    stable = cross[cross["replicated_pair_q05_at_least_half"]].copy()
    stable["min_cross_cohort_q05_fraction"] = stable[["GSE132465_patient_q05_fraction", "GSE200997_sample_q05_fraction"]].min(axis=1)
    stable = stable.sort_values(
        ["min_cross_cohort_q05_fraction", "GSE132465_patient_q05_fraction", "GSE200997_sample_q05_fraction", "GSE132465_median_prob"],
        ascending=False,
    )
    stable.to_csv(OUT / "CellChat_consensus_recurrent_edges.tsv", sep="\t", index=False)
    pathways = stable.groupby(["axis", "pathway_name"], observed=True).agg(
        n_recurrent_interactions=("ligand_receptor", "size"),
        median_primary_recurrence=("GSE132465_patient_q05_fraction", "median"),
        median_sensitivity_recurrence=("GSE200997_sample_q05_fraction", "median"),
    ).reset_index().sort_values(["axis", "n_recurrent_interactions", "median_primary_recurrence"], ascending=[True, False, False])
    pathways.to_csv(OUT / "CellChat_consensus_pathways.tsv", sep="\t", index=False)

    # Dot plot of the 29 interactions passing the predeclared recurrence criterion.
    plot_df = stable.head(29).copy().reset_index(drop=True)
    plot_df["label"] = plot_df.apply(
        lambda r: f"{r['source']} → {r['target']}  |  {r['ligand']}–{r['receptor']}  ({r['pathway_name']})", axis=1
    )
    y = np.arange(len(plot_df))
    fig, ax = plt.subplots(figsize=(12, max(7, 0.34 * len(plot_df))))
    ax.scatter(plot_df["GSE132465_patient_q05_fraction"], y - 0.12, s=42, color="#4477AA", label="GSE132465 patients")
    ax.scatter(plot_df["GSE200997_sample_q05_fraction"], y + 0.12, s=42, color="#CC6677", label="GSE200997 tumor samples")
    ax.set_yticks(y, plot_df["label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("Fraction of eligible units with within-pair BH q < 0.05")
    ax.set_title("Cross-cohort recurrent CellChat ligand–receptor edges")
    ax.grid(axis="x", color="#dddddd", linewidth=.7)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "CellChat_consensus_recurrence.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    primary = pd.read_csv(OUT / "GSE132465_CellChat_recurrence_by_patient.tsv.gz", sep="\t")
    sensitivity = pd.read_csv(OUT / "GSE200997_CellChat_recurrence_by_sample.tsv.gz", sep="\t")
    eligible = pd.read_csv(OUT / "CellChat_axis_eligible_unit_summary.tsv", sep="\t")
    meta132 = pd.read_csv(OUT / "GSE132465_cell_annotations.tsv.gz", sep="\t")
    meta200 = pd.read_csv(OUT / "GSE200997_cell_annotations.tsv.gz", sep="\t")
    broad200 = pd.read_csv(OUT / "GSE200997_broad_marker_annotations.tsv.gz", sep="\t")
    target_tests = pd.read_csv(OUT / "target_gene_patient_level_tests.tsv", sep="\t")
    manifest = json.loads((OUT / "preparation_manifest.json").read_text(encoding="utf-8"))

    n_rep = len(stable)
    primary_sig = int((primary["n_pair_q05"] > 0).sum())
    primary_half = int((primary["pair_q05_fraction"] >= 0.5).sum())
    sensitivity_sig = int((sensitivity["n_pair_q05"] > 0).sum())
    sensitivity_half = int((sensitivity["pair_q05_fraction"] >= 0.5).sum())
    q_rep_axis = stable.groupby("axis", observed=True).size().to_dict()
    primary_types = meta132[meta132["cellchat_include"]].groupby("celltype", observed=True).size().to_dict()
    sensitivity_types = meta200[meta200["cellchat_include"]].groupby("celltype", observed=True).size().to_dict()
    broad_counts = broad200["broad_type"].value_counts().to_dict()
    eligible_dict = {(r.cohort, r.axis): int(r.n_eligible_units) for r in eligible.itertuples()}

    test_lookup = {(r.cohort, r.comparison): r for r in target_tests.itertuples()}
    scd132 = test_lookup[("GSE132465", "Epithelial SCD mean vs lipid-program score excluding SCD")]
    pid132 = test_lookup[("GSE132465", "PID1 mean: Monocyte minus Macrophage")]
    scd200 = test_lookup[("GSE200997", "Epithelial SCD mean vs lipid-program score excluding SCD")]
    pid200 = test_lookup[("GSE200997", "PID1 mean: Monocyte minus Macrophage")]

    top_edges = stable.head(12)
    top_md = "\n".join(
        f"- `{r.source} → {r.target}`: **{r.ligand}–{r.receptor}** ({r.pathway_name}); BH-significant in {int(r.GSE132465_n_patient_q05)}/{int(r.GSE132465_n_eligible)} primary patients and {int(r.GSE200997_n_sample_q05)}/{int(r.GSE200997_n_eligible)} sensitivity samples."
        for r in top_edges.itertuples()
    )
    path_md = "\n".join(
        f"- **{axis}**: " + ", ".join(f"{p} ({n})" for p, n in pathways[pathways.axis.eq(axis)].head(5)[["pathway_name", "n_recurrent_interactions"]].itertuples(index=False, name=None))
        for axis in ["Epithelial_CAF", "Epithelial_Myeloid"]
    )

    readme = f'''# GSE132465 CellChat and subtype analysis

Analysis date: 2026-10-11  
Input: original GEO processed count matrices and original sample/cell annotations only. Earlier derived work results and the cached target-gene matrices were not read.

## Main results

- GSE132465 includes 23 tumor patient units. There were **{eligible_dict.get(('GSE132465','Epithelial_CAF'),0)}/23** units eligible for epithelial–CAF analysis and **{eligible_dict.get(('GSE132465','Epithelial_Myeloid'),0)}/23** for epithelial–myeloid analysis (each participating group had more than 10 cells and passed CellChat's internal filter).
- GSE200997 includes 16 tumor specimens but no patient identifier in its GEO processed annotation. **{eligible_dict.get(('GSE200997','Epithelial_CAF'),0)}/16** specimens were eligible for epithelial–CAF and **{eligible_dict.get(('GSE200997','Epithelial_Myeloid'),0)}/16** for epithelial–myeloid sensitivity analysis.
- After plus-one correction of permutation p-values and within-sender/receiver BH adjustment, **{n_rep} ligand–receptor edges** were significant in at least half of the eligible units in both cohorts (minimum two units per cohort). In the primary cohort, {primary_sig}/{len(primary)} summarized edges were significant in at least one patient and {primary_half} recurred in at least half of the units eligible for their sender/receiver pair. In sensitivity, the corresponding counts were {sensitivity_sig}/{len(sensitivity)} and {sensitivity_half}.
- Recurrent cross-cohort edges were most often in matrix-associated CAF-to-epithelium signaling (Laminin/Collagen families), epithelial-to-CAF MDK interactions, and epithelial-to-myeloid APP/MIF interactions.

### Most recurrent edges

{top_md}

### Pathway labels among recurrent edges

{path_md}

## Cell subtype reannotation

### GSE132465

The original GEO tumor labels were used only to identify epithelial, myeloid, and stromal candidate compartments. Myeloid and stromal candidates were independently reclustered from their full raw transcriptome using Scanpy (2,000 highly variable genes, PCA, 15-neighbor graph, Leiden resolution 0.5). Cluster labels were assigned from marker-program scores; **SCD and PID1 were explicitly excluded from HVG/PCA and were not used to define labels**.

CellChat analysis included {len(meta132[meta132.cellchat_include])} cells: {', '.join(f'{k} {v:,}' for k,v in primary_types.items())}. Reannotated tumor myeloid cells were {primary_types.get('Monocyte',0):,} monocytes, {primary_types.get('Macrophage',0):,} macrophages, and {primary_types.get('Other_myeloid',0):,} other myeloid cells. Stromal candidates contained {primary_types.get('ECM_CAF',0):,} ECM-CAF, {primary_types.get('myCAF',0):,} myCAF, and {primary_types.get('apCAF',0):,} apCAF cells; no distinct iCAF cluster was identified by the marker-program rule. apCAF was too rare for reliable cohort-wide communication assessment.

### GSE200997 sensitivity

Its GEO annotation provides tumor/normal, sample, and MSI fields but no cell-type labels. Tumor epithelial, myeloid, and fibroblast candidates were therefore inferred from multiple canonical marker modules, competing against T/NK, B/plasma, endothelial, pericyte, and mast-cell controls. Ambiguous cells were left unassigned. The broad marker pass assigned {broad_counts.get('Epithelial',0):,} epithelial, {broad_counts.get('Myeloid',0):,} myeloid, and {broad_counts.get('Fibroblast',0):,} fibroblast candidates; the rest were not forced into these lineages. Full-transcriptome reclustering produced {', '.join(f'{k} {v:,}' for k,v in sensitivity_types.items())} cells used for CellChat. No patient field was supplied, so specimen/sample was the independent unit.

## CellChat settings and statistics

- CellChat {manifest.get('CellChat_version','2.2.0.9001')} with `CellChatDB.human` (3,233 database interactions); expression was `log1p(raw UMI / complete-cell UMI total × 10,000)`.
- One CellChat object was fit per GSE132465 patient and per GSE200997 tumor sample. Groups with 10 or fewer cells were removed by CellChat's internal filter; the eligible-unit denominator requires more than 10 cells in every participating group. Epithelial–CAF and epithelial–myeloid directions were extracted separately from each patient/sample network.
- Average expression: `triMean`; `population.size=FALSE`; 100 CellChat permutations per unit. Because CellChat empirical permutation p-values can be zero, p was floored at `1/(100+1)` before BH adjustment within each sender–receiver pair. Patient/sample recurrence fractions are descriptive consistency summaries; CellChat cell-level permutation p-values are not patient-level inferential tests.
- The cross-cohort recurrence rule was BH q < 0.05 in at least 50% of eligible GSE132465 patients and GSE200997 samples, with at least two eligible units in each cohort.

## SCD and PID1 post-annotation summaries

- SCD was summarized only in tumor epithelium against a custom 40-gene lipid-metabolism marker score that excludes SCD. Per-patient epithelial SCD and score were correlated in GSE132465 (Spearman ρ={scd132.statistic:.3f}, BH p={scd132.p_BH:.3g}; leave-one-patient-out ρ range {scd132.leave_one_out_stat_min:.3f}–{scd132.leave_one_out_stat_max:.3f}). GSE200997 showed the same direction at the sample level (ρ={scd200.statistic:.3f}, BH p={scd200.p_BH:.3g}; leave-one-sample-out ρ range {scd200.leave_one_out_stat_min:.3f}–{scd200.leave_one_out_stat_max:.3f}). This is an exploratory custom marker score, not a standardized signature.
- PID1 was examined after myeloid subtype labels were fixed. Patient-level mean expression was higher in monocytes than macrophages in GSE132465 (n={int(pid132.n_units)} paired patients; median paired difference={pid132.median_difference:.3f}, BH p={pid132.p_BH:.3g}); GSE200997 was directionally consistent at the sample level (n={int(pid200.n_units)}; median difference={pid200.median_difference:.3f}, BH p={pid200.p_BH:.3g}). Neither target gene was used to define a cluster.

## Files

- `GSE132465_myeloid_umap.png`, `GSE132465_stromal_umap.png`, `GSE200997_myeloid_umap.png`, `GSE200997_stromal_umap.png`: reclustered cells.
- `CellChat_consensus_recurrence.png`: cross-cohort recurrence of recurrent edges.
- `CellChat_consensus_recurrent_edges.tsv`: the {n_rep} recurrent interactions with both-cohort denominators and probabilities.
- `GSE132465_CellChat_patient_sample_interactions.tsv.gz`, `GSE200997_CellChat_patient_sample_interactions.tsv.gz`: per-unit CellChat interactions and corrected p/q values.
- `GSE132465_CellChat_recurrence_by_patient.tsv.gz`, `GSE200997_CellChat_recurrence_by_sample.tsv.gz`, `CellChat_cross_cohort_sensitivity.tsv.gz`: recurrence summaries.
- `GSE132465_cell_annotations.tsv.gz`, `GSE200997_cell_annotations.tsv.gz`, cluster marker tables, target-gene summaries, and cell-count/eligibility tables.

## Interpretation limits

CellChat predicts candidate ligand–receptor communication from expression, not direct physical signaling. GSE200997 cell identities are marker-inferred and its eligible CAF sample count is small, so it is a sensitivity check rather than an equally powered replication. Rare apCAF and any state present in very few patients should not be described as stable. The strongest cross-cohort patterns are ECM/Laminin/Collagen from ECM-CAF toward epithelium and APP/MIF from epithelium toward monocyte/macrophage groups; experimental validation is still required.
'''
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print(f"Wrote {len(stable)} recurrent edges and report to {OUT}")
    print(pathways.to_string(index=False))


if __name__ == "__main__":
    main()
