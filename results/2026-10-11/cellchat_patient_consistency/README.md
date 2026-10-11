# GSE132465 CellChat and subtype analysis

Analysis date: 2026-10-11  
Input: original GEO processed count matrices and original sample/cell annotations only. Earlier derived work results and the cached target-gene matrices were not read.

## Main results

- GSE132465 includes 23 tumor patient units. There were **18/23** units eligible for epithelial–CAF analysis and **23/23** for epithelial–myeloid analysis (each participating group had more than 10 cells and passed CellChat's internal filter).
- GSE200997 includes 16 tumor specimens but no patient identifier in its GEO processed annotation. **7/16** specimens were eligible for epithelial–CAF and **13/16** for epithelial–myeloid sensitivity analysis.
- After plus-one correction of permutation p-values and within-sender/receiver BH adjustment, **29 ligand–receptor edges** were significant in at least half of the eligible units in both cohorts (minimum two units per cohort). In the primary cohort, 974/1094 summarized edges were significant in at least one patient and 271 recurred in at least half of the units eligible for their sender/receiver pair. In sensitivity, the corresponding counts were 166/182 and 30.
- Recurrent cross-cohort edges were most often in matrix-associated CAF-to-epithelium signaling (Laminin/Collagen families), epithelial-to-CAF MDK interactions, and epithelial-to-myeloid APP/MIF interactions.

### Most recurrent edges

- `Epithelial → ECM_CAF`: **MDK–SDC2** (MK); BH-significant in 18/18 primary patients and 7/7 sensitivity samples.
- `Epithelial → ECM_CAF`: **MDK–LRP1** (MK); BH-significant in 18/18 primary patients and 7/7 sensitivity samples.
- `Epithelial → Monocyte`: **APP–CD74** (APP); BH-significant in 19/19 primary patients and 11/12 sensitivity samples.
- `Epithelial → Macrophage`: **APP–CD74** (APP); BH-significant in 22/23 primary patients and 11/12 sensitivity samples.
- `Epithelial → ECM_CAF`: **GDF15–TGFBR2** (GDF); BH-significant in 15/18 primary patients and 7/7 sensitivity samples.
- `ECM_CAF → Epithelial`: **LAMA4–ITGA6_ITGB1** (LAMININ); BH-significant in 18/18 primary patients and 5/7 sensitivity samples.
- `ECM_CAF → Epithelial`: **LAMC1–ITGA6_ITGB1** (LAMININ); BH-significant in 17/18 primary patients and 5/7 sensitivity samples.
- `ECM_CAF → Epithelial`: **LAMB1–ITGA6_ITGB1** (LAMININ); BH-significant in 17/18 primary patients and 5/7 sensitivity samples.
- `Epithelial → ECM_CAF`: **PRSS3–F2R** (PARs); BH-significant in 16/18 primary patients and 5/7 sensitivity samples.
- `Epithelial → Monocyte`: **MIF–CD74_CD44** (MIF); BH-significant in 19/19 primary patients and 8/12 sensitivity samples.
- `Epithelial → ECM_CAF`: **PPIA–BSG** (CypA); BH-significant in 11/18 primary patients and 6/7 sensitivity samples.
- `Epithelial → Macrophage`: **MIF–CD74_CD44** (MIF); BH-significant in 22/23 primary patients and 7/12 sensitivity samples.

### Pathway labels among recurrent edges

- **Epithelial_CAF**: LAMININ (9), COLLAGEN (7), MK (2), FN1 (1), PARs (1)
- **Epithelial_Myeloid**: MIF (4), APP (3)

## Cell subtype reannotation

### GSE132465

The original GEO tumor labels were used only to identify epithelial, myeloid, and stromal candidate compartments. Myeloid and stromal candidates were independently reclustered from their full raw transcriptome using Scanpy (2,000 highly variable genes, PCA, 15-neighbor graph, Leiden resolution 0.5). Cluster labels were assigned from marker-program scores; **SCD and PID1 were explicitly excluded from HVG/PCA and were not used to define labels**.

CellChat analysis included 25755 cells: ECM_CAF 1,489, Epithelial 17,469, Macrophage 3,286, Monocyte 2,837, Other_myeloid 277, apCAF 23, myCAF 374. Reannotated tumor myeloid cells were 2,837 monocytes, 3,286 macrophages, and 277 other myeloid cells. Stromal candidates contained 1,489 ECM-CAF, 374 myCAF, and 23 apCAF cells; no distinct iCAF cluster was identified by the marker-program rule. apCAF was too rare for reliable cohort-wide communication assessment.

### GSE200997 sensitivity

Its GEO annotation provides tumor/normal, sample, and MSI fields but no cell-type labels. Tumor epithelial, myeloid, and fibroblast candidates were therefore inferred from multiple canonical marker modules, competing against T/NK, B/plasma, endothelial, pericyte, and mast-cell controls. Ambiguous cells were left unassigned. The broad marker pass assigned 6,871 epithelial, 917 myeloid, and 479 fibroblast candidates; the rest were not forced into these lineages. Full-transcriptome reclustering produced ECM_CAF 435, Epithelial 6,871, Macrophage 430, Monocyte 487 cells used for CellChat. No patient field was supplied, so specimen/sample was the independent unit.

## CellChat settings and statistics

- CellChat 2.2.0.9001 with `CellChatDB.human` (3,233 database interactions); expression was `log1p(raw UMI / complete-cell UMI total × 10,000)`.
- One CellChat object was fit per GSE132465 patient and per GSE200997 tumor sample. Groups with 10 or fewer cells were removed by CellChat's internal filter; the eligible-unit denominator requires more than 10 cells in every participating group. Epithelial–CAF and epithelial–myeloid directions were extracted separately from each patient/sample network.
- Average expression: `triMean`; `population.size=FALSE`; 100 CellChat permutations per unit. Because CellChat empirical permutation p-values can be zero, p was floored at `1/(100+1)` before BH adjustment within each sender–receiver pair. Patient/sample recurrence fractions are descriptive consistency summaries; CellChat cell-level permutation p-values are not patient-level inferential tests.
- The cross-cohort recurrence rule was BH q < 0.05 in at least 50% of eligible GSE132465 patients and GSE200997 samples, with at least two eligible units in each cohort.

## SCD and PID1 post-annotation summaries

- SCD was summarized only in tumor epithelium against a custom 40-gene lipid-metabolism marker score that excludes SCD. Per-patient epithelial SCD and score were correlated in GSE132465 (Spearman ρ=0.575, BH p=0.00409; leave-one-patient-out ρ range 0.514–0.643). GSE200997 showed the same direction at the sample level (ρ=0.500, BH p=0.0486; leave-one-sample-out ρ range 0.404–0.600). This is an exploratory custom marker score, not a standardized signature.
- PID1 was examined after myeloid subtype labels were fixed. Patient-level mean expression was higher in monocytes than macrophages in GSE132465 (n=23 paired patients; median paired difference=0.134, BH p=0.00148); GSE200997 was directionally consistent at the sample level (n=15; median difference=0.111, BH p=0.00575). Neither target gene was used to define a cluster.

## Files

- `GSE132465_myeloid_umap.png`, `GSE132465_stromal_umap.png`, `GSE200997_myeloid_umap.png`, `GSE200997_stromal_umap.png`: reclustered cells.
- `CellChat_consensus_recurrence.png`: cross-cohort recurrence of recurrent edges.
- `CellChat_consensus_recurrent_edges.tsv`: the 29 recurrent interactions with both-cohort denominators and probabilities.
- `GSE132465_CellChat_patient_sample_interactions.tsv.gz`, `GSE200997_CellChat_patient_sample_interactions.tsv.gz`: per-unit CellChat interactions and corrected p/q values.
- `GSE132465_CellChat_recurrence_by_patient.tsv.gz`, `GSE200997_CellChat_recurrence_by_sample.tsv.gz`, `CellChat_cross_cohort_sensitivity.tsv.gz`: recurrence summaries.
- `GSE132465_cell_annotations.tsv.gz`, `GSE200997_cell_annotations.tsv.gz`, cluster marker tables, target-gene summaries, and cell-count/eligibility tables.

## Interpretation limits

CellChat predicts candidate ligand–receptor communication from expression, not direct physical signaling. GSE200997 cell identities are marker-inferred and its eligible CAF sample count is small, so it is a sensitivity check rather than an equally powered replication. Rare apCAF and any state present in very few patients should not be described as stable. The strongest cross-cohort patterns are ECM/Laminin/Collagen from ECM-CAF toward epithelium and APP/MIF from epithelium toward monocyte/macrophage groups; experimental validation is still required.
