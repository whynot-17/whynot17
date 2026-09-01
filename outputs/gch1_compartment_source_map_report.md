# GCH1 compartment source map

## Question

Does the bulk right-sided GCH1 signal localize to malignant/tumor epithelial cells, or is it mainly contributed by immune/stromal compartments?

## Design

- Patient-level pseudobulk from raw UMI counts; cells were not treated as independent replicates.
- GSE200997: marker-defined malignant epithelial-like, immune-like and stromal-like compartments.
- GSE132465: official tumor epithelial, T-cell, B-cell, myeloid, stromal and mast-cell compartments.
- AA-routing proxy: mean within-compartment z-score of PLA2G4A, PTGS2 and PTGES.

## Right-versus-left GCH1 by compartment

| Cohort | Compartment | Right n | Left n | Right−left mean | Welch P | MW P |
|---|---|---:|---:|---:|---:|---:|
| GSE132465 | B_cells | 11 | 12 | -0.207 | 0.721 | 0.255 |
| GSE132465 | T_cells | 11 | 12 | -0.045 | 0.773 | 0.975 |
| GSE132465 | myeloid | 11 | 12 | -0.364 | 0.0927 | 0.0905 |
| GSE132465 | stromal | 11 | 12 | 0.377 | 0.392 | 0.479 |
| GSE132465 | tumor_epithelial | 11 | 12 | 0.189 | 0.347 | 0.186 |
| GSE200997 | immune_like | 8 | 8 | -1.031 | 0.0111 | 0.0148 |
| GSE200997 | malignant_epithelial_like | 8 | 8 | -0.163 | 0.873 | 0.834 |
| GSE200997 | stromal_like | 7 | 7 | 0.033 | 0.972 | 0.936 |

## Coupling within right-sided patients

| Cohort | Compartment | Right n | Spearman ρ | P | FDR |
|---|---|---:|---:|---:|---:|
| GSE132465 | B_cells | 11 | 0.136 | 0.689 | 0.862 |
| GSE132465 | T_cells | 11 | 0.327 | 0.326 | 0.815 |
| GSE132465 | myeloid | 11 | 0.482 | 0.133 | 0.667 |
| GSE132465 | stromal | 11 | 0.191 | 0.574 | 0.862 |
| GSE132465 | tumor_epithelial | 11 | 0.055 | 0.873 | 0.873 |
| GSE200997 | immune_like | 8 | 0.714 | 0.0465 | 0.0698 |
| GSE200997 | malignant_epithelial_like | 8 | 0.786 | 0.0208 | 0.0624 |
| GSE200997 | stromal_like | 7 | 0.584 | 0.168 | 0.168 |

## Interpretation

A consistent epithelial signal would support a tumor-cell-associated GCH1 state. A signal confined to immune or stromal compartments would instead argue that the bulk result is not a malignant epithelial-intrinsic program. These analyses do not establish GCH1 dependency, BH4 abundance, AA flux or causality.

## Provenance

- GSE200997: {'cohort': 'GSE200997', 'annotation_cells': 49859, 'matrix_gene_rows': 23828, 'compartment_cell_counts': {'malignant_epithelial_like': 6681, 'immune_like': 23261, 'stromal_like': 976}, 'selection': 'malignant_epithelial_like from existing marker rule; immune_like/stromal_like assigned by dominant lineage marker score.'}
- GSE132465: {'cohort': 'GSE132465', 'annotation_cells': 63689, 'matrix_gene_rows': 33694, 'compartment_cell_counts': {'tumor_epithelial': 17469, 'T_cells': 16739, 'B_cells': 3938, 'myeloid': 6400, 'stromal': 2736, 'mast': 3}, 'selection': 'Official GEO annotation: Class=Tumor and Cell_type-specific compartment.'}
- GSE132465 sidedness was assigned from GEO sample region metadata.

## Files

- `gch1_compartment_source_map_patient_pseudobulk.csv`
- `gch1_compartment_source_map_sidedness.csv`
- `gch1_compartment_source_map_coupling.csv`
- `gch1_compartment_source_map.png`
