# CRC antioxidant-buffering compartment source map

## Question

Does the right-sided SLC7A11/GCH1/antioxidant-buffering signal localize to malignant epithelial cells, myeloid cells, T/B cells or stromal cells?

## Design

- Raw UMI counts were aggregated to patient × sample × compartment pseudobulk.
- The statistical unit is the patient; cells were not treated as independent observations.
- SLC7A11 and GCH1 are log2(CPM + 1) pseudobulk values.
- Antioxidant buffering is the mean within-cohort/compartment z-score of SLC7A11, GPX4, AIFM2, GCH1 and DHODH, requiring at least 4 of 5 genes.
- GSE200997 immune-like cells were split using dominant raw-UMI myeloid versus T/B marker scores; GSE132465 used official cell-type labels.

## Right-versus-left results

| Cohort | Compartment | Metric | Right n | Left n | Right−left | Welch P | MW P |
|---|---|---|---:|---:|---:|---:|---:|
| GSE132465 | T_B | SLC7A11 | 11 | 12 | 0.036 | 0.904 | 0.951 |
| GSE132465 | T_B | GCH1 | 11 | 12 | -0.354 | 0.0215 | 0.0151 |
| GSE132465 | T_B | antioxidant_buffering_score | 11 | 12 | -0.558 | 0.014 | 0.0151 |
| GSE132465 | malignant_epithelial | SLC7A11 | 11 | 12 | 0.836 | 0.0248 | 0.021 |
| GSE132465 | malignant_epithelial | GCH1 | 11 | 12 | 0.189 | 0.347 | 0.186 |
| GSE132465 | malignant_epithelial | antioxidant_buffering_score | 11 | 12 | 0.231 | 0.225 | 0.31 |
| GSE132465 | myeloid | SLC7A11 | 11 | 12 | -0.610 | 0.101 | 0.0605 |
| GSE132465 | myeloid | GCH1 | 11 | 12 | -0.364 | 0.0927 | 0.0905 |
| GSE132465 | myeloid | antioxidant_buffering_score | 11 | 12 | -0.581 | 0.00949 | 0.0151 |
| GSE132465 | stromal | SLC7A11 | 11 | 12 | -1.333 | 0.0607 | 0.0736 |
| GSE132465 | stromal | GCH1 | 11 | 12 | 0.377 | 0.392 | 0.479 |
| GSE132465 | stromal | antioxidant_buffering_score | 11 | 12 | 0.044 | 0.857 | 0.735 |
| GSE200997 | T_B | SLC7A11 | 8 | 8 | -0.049 | 0.613 | 0.408 |
| GSE200997 | T_B | GCH1 | 8 | 8 | -0.413 | 0.627 | 0.462 |
| GSE200997 | T_B | antioxidant_buffering_score | 8 | 8 | -0.111 | 0.754 | 0.234 |
| GSE200997 | malignant_epithelial | SLC7A11 | 8 | 8 | 0.785 | 0.359 | 0.371 |
| GSE200997 | malignant_epithelial | GCH1 | 8 | 8 | -0.163 | 0.873 | 0.834 |
| GSE200997 | malignant_epithelial | antioxidant_buffering_score | 8 | 8 | 0.024 | 0.934 | 0.959 |
| GSE200997 | myeloid | SLC7A11 | 8 | 8 | 0.175 | 0.809 | 0.873 |
| GSE200997 | myeloid | GCH1 | 8 | 8 | -1.930 | 0.0307 | 0.0404 |
| GSE200997 | myeloid | antioxidant_buffering_score | 8 | 8 | -0.382 | 0.253 | 0.328 |
| GSE200997 | stromal | SLC7A11 | 7 | 7 | 0.666 | 0.539 | 0.551 |
| GSE200997 | stromal | GCH1 | 7 | 7 | 0.033 | 0.972 | 0.936 |
| GSE200997 | stromal | antioxidant_buffering_score | 7 | 7 | -0.031 | 0.925 | 1 |

## Interpretation guardrails

A replicated malignant-epithelial right-minus-left signal supports a tumor-cell-associated state. A signal confined to myeloid, T/B or stromal compartments supports a microenvironmental source. These data do not establish AA concentration, AA flux, BH4 abundance, ferroptosis resistance or functional dependency.

## Provenance

- GSE200997: {"cohort": "GSE200997", "annotation_cells": 49859, "matrix_gene_rows": 23828, "compartment_info": {"malignant_epithelial": {"cohort": "GSE200997", "compartment": "malignant_epithelial", "n_patients": 16, "n_cells": 6681, "side_counts": {"Left": 8, "Right": 8}}, "myeloid": {"cohort": "GSE200997", "compartment": "myeloid", "n_patients": 16, "n_cells": 3555, "side_counts": {"Left": 8, "Right": 8}}, "T_B": {"cohort": "GSE200997", "compartment": "T_B", "n_patients": 16, "n_cells": 13490, "side_counts": {"Left": 8, "Right": 8}}, "stromal": {"cohort": "GSE200997", "compartment": "stromal", "n_patients": 14, "n_cells": 976, "side_counts": {"Left": 7, "Right": 7}}}, "lineage_cell_counts": {"unresolved": 25157, "T_B": 13490, "malignant_epithelial": 6681, "myeloid": 3555, "stromal": 976}, "marker_genes_available": {"myeloid": ["LST1", "TYROBP", "FCER1G", "CTSS", "LILRB1", "AIF1", "CSTA", "LGALS3", "S100A8", "S100A9", "FCGR3A"], "T_B": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD79A", "MS4A1", "CD37", "CD74", "HLA-DRA", "CD52", "LTB", "IL7R"]}, "selection": "Existing marker-defined malignant epithelial/immune/stromal calls; immune-like cells split by dominant raw-UMI myeloid versus T/B marker score."}
- GSE132465: {"cohort": "GSE132465", "annotation_cells": 63689, "matrix_gene_rows": 33694, "compartment_info": {"malignant_epithelial": {"cohort": "GSE132465", "compartment": "malignant_epithelial", "n_patients": 23, "n_cells": 17469, "side_counts": {"Left": 12, "Right": 11}}, "myeloid": {"cohort": "GSE132465", "compartment": "myeloid", "n_patients": 23, "n_cells": 6400, "side_counts": {"Left": 12, "Right": 11}}, "T_B": {"cohort": "GSE132465", "compartment": "T_B", "n_patients": 23, "n_cells": 20677, "side_counts": {"Left": 12, "Right": 11}}, "stromal": {"cohort": "GSE132465", "compartment": "stromal", "n_patients": 23, "n_cells": 2736, "side_counts": {"Left": 12, "Right": 11}}}, "lineage_cell_counts": {"T_B": 32261, "malignant_epithelial": 18539, "myeloid": 6769, "stromal": 5933, "unresolved": 187}, "selection": "Official GEO annotation: Class=Tumor; Epithelial cells, Myeloids, T cells/B cells and Stromal cells mapped to the four compartments.", "side_mapping": "Existing GSE132465 GEO region map: cecum/ascending/hepatic flexure/transverse=Right; splenic flexure/descending/sigmoid/rectosigmoid/rectum=Left."}

## Files

- `sc_antioxidant_compartment_patient_pseudobulk.csv`
- `sc_antioxidant_compartment_cell_lineage_audit.csv`
- `sc_antioxidant_compartment_sidedness.csv`
- `sc_antioxidant_compartment_cell_counts.csv`
- `sc_antioxidant_compartment_source_map.png`
