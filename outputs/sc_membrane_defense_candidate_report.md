# SLC7A11, MBOAT2 and MBOAT1: epithelial single-cell pseudobulk

## Scope

This bounded analysis tests only SLC7A11, MBOAT2 and MBOAT1 in malignant/tumor epithelial patient-level pseudobulk from GSE200997 and GSE132465.

- Right-versus-left expression was tested in each cohort.
- AA-routing coupling was tested within right-sided patients.
- The primary proxy is PLA2G4A/PTGS2/PTGES; the sensitivity proxy adds ALOX5/ALOX5AP.
- Cells were not treated as independent observations.
- The proxy is not measured AA concentration and no dependency/causality is inferred.

## Primary right-versus-left results

| Cohort | Gene | Right n | Left n | Right−left mean z | Welch P | FDR |
|---|---|---:|---:|---:|---:|---:|
| GSE200997 | SLC7A11 | 8 | 8 | 0.476 | 0.359 | 0.511 |
| GSE200997 | MBOAT2 | 8 | 8 | -0.345 | 0.511 | 0.511 |
| GSE200997 | MBOAT1 | 8 | 8 | -0.429 | 0.411 | 0.511 |
| GSE132465 | SLC7A11 | 11 | 12 | 0.943 | 0.0248 | 0.0744 |
| GSE132465 | MBOAT2 | 11 | 12 | 0.601 | 0.154 | 0.154 |
| GSE132465 | MBOAT1 | 11 | 12 | 0.744 | 0.0884 | 0.133 |

## Primary coupling within right-sided patients

| Cohort | Gene | Right n | Spearman ρ | P | FDR |
|---|---|---:|---:|---:|---:|
| GSE200997 | SLC7A11 | 8 | 0.738 | 0.0366 | 0.11 |
| GSE200997 | MBOAT2 | 8 | 0.524 | 0.183 | 0.183 |
| GSE200997 | MBOAT1 | 8 | 0.548 | 0.16 | 0.183 |
| GSE132465 | SLC7A11 | 11 | 0.055 | 0.873 | 0.873 |
| GSE132465 | MBOAT2 | 11 | 0.600 | 0.051 | 0.153 |
| GSE132465 | MBOAT1 | 11 | 0.309 | 0.355 | 0.533 |

## Interpretation

A candidate is prioritized only when its direction is reasonably consistent across both single-cell cohorts. A positive coupling in one cohort without replication is treated as exploratory.
These results test epithelial transcriptomic state, not membrane lipid composition, AA flux, ferroptosis resistance or functional dependency.

## Provenance

- GSE200997: {'cohort': 'GSE200997', 'annotation_cells': 49859, 'selected_cells': 6681, 'matrix_gene_rows': 23828, 'side_counts': {'Left': 8, 'Right': 8}, 'selection_rule': 'Existing marker-defined malignant_epithelial flag restricted to tumor cells.'}
- GSE132465: {'cohort': 'GSE132465', 'annotation_cells': 63689, 'selected_cells': 17469, 'matrix_gene_rows': 33694, 'side_counts': {'Left': 12, 'Right': 11}, 'selection_rule': 'Official annotation: Class=Tumor and Cell_type=Epithelial cells.'}
- Pseudobulk normalization: summed raw UMI counts per patient and compartment, followed by log2(CPM + 1).

## Files

- `sc_membrane_defense_candidate_patient_pseudobulk.csv`
- `sc_membrane_defense_candidate_sidedness.csv`
- `sc_membrane_defense_candidate_coupling.csv`
- `sc_membrane_defense_candidate_pseudobulk.png`
