# Step 8B — T2D pathway redundancy reduction

- Status: **complete_redundancy_reduction**
- Input: **1,647** pathway terms with global BH-FDR < 0.05
- Frozen effective background: **6,069 genes**
- Modules after parent/ancestor reduction: **321**
- Compact representatives retained: **32** (maximum 8 per axis)

## Axis audit

| Axis | Significant terms | Modules | Compact representatives |
|---|---:|---:|---:|
| cluster_11 | 86 | 33 | 8 |
| cluster_5 | 552 | 125 | 8 |
| cluster_6 | 790 | 123 | 8 |
| cluster_8 | 219 | 40 | 8 |

## Interpretation boundary

All significant terms remain in `t2d_step8_pathway_modules.csv`; module representatives are a deterministic reduction layer, not a new statistical test.
GO/Reactome/KEGG parent structure was used where available. Broad terms covering >25% of the effective background were not allowed to bridge modules. Cross-source similarity is lexical only and is not described as semantic evidence.
A representative term is chosen using frozen global q value, intersection size, and term specificity. This stage does not infer pathway direction, activation, exposure causality, or mediation of T2D.

## Compact representatives

| Axis | Source | Representative | Overlap | Term size | Global q |
|---|---|---|---:|---:|---:|
| cluster_11 | GO:BP | intracellular transport | 47 | 468 | 1.1553979405295108e-07 |
| cluster_11 | GO:BP | substantia nigra development | 8 | 23 | 0.00011796407674934093 |
| cluster_11 | REAC | Membrane Trafficking | 22 | 206 | 0.0006835160176264031 |
| cluster_11 | REAC | Infectious disease | 34 | 419 | 0.0011667831020124822 |
| cluster_11 | REAC | HSF1 activation | 4 | 8 | 0.005085903596383851 |
| cluster_11 | REAC | Metabolism of proteins | 49 | 711 | 0.001303279767502737 |
| cluster_11 | GO:BP | catabolic process | 59 | 930 | 0.0017976772877818076 |
| cluster_11 | REAC | Antigen Presentation: Folding, assembly and peptide loading of class I MHC | 5 | 17 | 0.010649589846916515 |
| cluster_5 | REAC | Xenobiotics | 15 | 20 | 1.7188678466533358e-23 |
| cluster_5 | GO:BP | xenobiotic metabolic process | 19 | 64 | 1.0752430073191622e-20 |
| cluster_5 | KEGG | Metabolism of xenobiotics by cytochrome P450 | 14 | 34 | 5.0390721789362755e-17 |
| cluster_5 | KEGG | Chemical carcinogenesis - DNA adducts | 13 | 33 | 1.9911752222692075e-15 |
| cluster_5 | KEGG | Chemical carcinogenesis - receptor activation | 19 | 117 | 1.6957934096643257e-15 |
| cluster_5 | GO:BP | secondary metabolic process | 10 | 34 | 3.57454711985704e-10 |
| cluster_5 | GO:BP | oxidative demethylation | 6 | 8 | 9.500692087540011e-09 |
| cluster_5 | KEGG | Pathways in cancer | 22 | 290 | 2.2282703525944775e-11 |
| cluster_6 | GO:BP | miRNA-mediated post-transcriptional gene silencing | 74 | 134 | 5.2442580399529704e-24 |
| cluster_6 | KEGG | MicroRNAs in cancer | 65 | 172 | 2.7319569574203865e-11 |
| cluster_6 | REAC | Signaling by Interleukins | 78 | 231 | 7.538027426042925e-11 |
| cluster_6 | REAC | Cellular responses to stimuli | 108 | 371 | 1.0507171924584459e-10 |
| cluster_6 | KEGG | Lipid and atherosclerosis | 53 | 139 | 2.876291459727408e-09 |
| cluster_6 | KEGG | IL-17 signaling pathway | 25 | 46 | 1.143553415492911e-07 |
| cluster_6 | REAC | Nuclear Events (kinase and transcription factor activation) | 21 | 36 | 4.2227427661854955e-07 |
| cluster_6 | KEGG | MAPK signaling pathway | 56 | 172 | 4.852038414543857e-07 |
| cluster_8 | REAC | Metabolism of RNA | 163 | 237 | 6.339673789434701e-08 |
| cluster_8 | GO:BP | cytoplasmic translation | 68 | 90 | 2.3046729437840455e-05 |
| cluster_8 | REAC | Metabolism of proteins | 416 | 711 | 4.546150012085074e-06 |
| cluster_8 | REAC | Response of EIF2AK4 (GCN2) to amino acid deficiency | 44 | 58 | 0.001272501486278387 |
| cluster_8 | REAC | Cell Cycle | 149 | 239 | 0.0008257106019958016 |
| cluster_8 | REAC | Selenocysteine synthesis | 39 | 51 | 0.002343064956477695 |
| cluster_8 | REAC | Protein-protein interactions at synapses | 39 | 51 | 0.002343064956477695 |
| cluster_8 | REAC | Regulation of expression of SLITs and ROBOs | 59 | 84 | 0.0025021766897243463 |
