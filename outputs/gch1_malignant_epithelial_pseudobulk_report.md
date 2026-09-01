# GCH1 in malignant epithelial single-cell pseudobulk

## Question

Does the right-sided CRC GCH1 signal remain detectable in malignant/tumor epithelial cells, and does it couple to the AA-routing transcriptomic proxy within right-sided patients?

## Design

- GSE200997: marker-defined malignant epithelial-like cells from tumor samples.
- GSE132465: official `Class=Tumor` and `Cell_type=Epithelial cells` annotation; malignancy inferred from tumor sample origin.
- Raw UMI counts were summed within each patient and epithelial compartment, then normalized as log2(CPM + 1).
- Statistical unit: patient; individual cells were not treated as replicates.
- AA-routing proxy: within-cohort mean z-score of PLA2G4A, PTGS2 and PTGES. It is not measured AA concentration.

## Right-versus-left GCH1

| Cohort | Right n | Left n | Right−left mean log2(CPM+1) | Welch P | Mann–Whitney P |
|---|---:|---:|---:|---:|---:|
| GSE200997 | 8 | 8 | -0.163 | 0.873 | 0.834 |
| GSE132465 | 11 | 12 | 0.189 | 0.347 | 0.186 |

## AA-routing coupling within right-sided patients

| Cohort | n | Spearman ρ | P |
|---|---:|---:|---:|
| GSE200997 | 8 | 0.786 | 0.0208 |
| GSE132465 | 11 | 0.055 | 0.873 |

## Right-sided proxy-high versus proxy-low

| Cohort | High n | Low n | High−low GCH1 mean | Welch P | Mann–Whitney P |
|---|---:|---:|---:|---:|---:|
| GSE200997 | 3 | 5 | 2.072 | 0.0567 | 0.0357 |
| GSE132465 | 8 | 3 | 0.010 | 0.978 | 0.921 |

## Interpretation

Persistence of a right-sided GCH1 signal in both epithelial pseudobulk datasets would support a tumor-cell-associated state rather than a purely bulk-composition explanation. Positive within-right coupling would support, but not prove, an AA-associated GCH1/BH4 adaptive program.
Neither result demonstrates tissue AA enrichment, BH4 abundance, ferroptosis resistance, causality or GCH1 dependency. Those require lipid/metabolite measurements and functional perturbation.

## Provenance

- GSE200997: {'cohort': 'GSE200997', 'n_patients': 16, 'side_counts': {'Left': 8, 'Right': 8}, 'n_cells': 6681, 'target_genes': ['GCH1', 'PLA2G4A', 'PTGS2', 'PTGES'], 'annotation_cells': 49859, 'selected_cells': 6681, 'matrix_gene_rows': 23828, 'selection_rule': 'Previously generated marker-defined malignant_epithelial flag, restricted to tumor cells; see gse200997 tumor_marker_scores.csv.'}
- GSE132465: {'cohort': 'GSE132465', 'n_patients': 23, 'side_counts': {'Left': 12, 'Right': 11}, 'n_cells': 17469, 'target_genes': ['GCH1', 'PLA2G4A', 'PTGS2', 'PTGES'], 'annotation_cells': 63689, 'selected_cells': 17469, 'matrix_gene_rows': 33694, 'selection_rule': 'Official annotation: Class=Tumor and Cell_type=Epithelial cells; malignancy inferred from tumor sample origin.', 'side_mapping': 'cecum/ascending/hepatic flexure/transverse=Right; splenic flexure/descending/sigmoid/rectosigmoid/rectum=Left.'}
- Sidedness for GSE132465 was assigned from GEO sample region metadata.

## Files

- `gch1_sc_patient_pseudobulk.csv`
- `gch1_sc_sidedness.csv`
- `gch1_sc_coupling.csv`
- `gch1_sc_high_low.csv`
- `gch1_malignant_epithelial_pseudobulk.png`
