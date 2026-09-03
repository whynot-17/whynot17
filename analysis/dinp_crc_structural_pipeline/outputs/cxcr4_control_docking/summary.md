# CXCR4 control-docking validation

- Target: CXCR4 (3ODU)
- Positive control: IT1t (PDB ligand ID ITD)
- DINP and IT1t were docked to the same prepared receptor and the same IT1t-defined pocket box.
- IT1t best Vina affinity: -6.402 kcal/mol
- DINP best Vina affinity: -4.783 kcal/mol
- DINP − IT1t score difference: 1.6189999999999998 kcal/mol
- IT1t redocking heavy-atom RMSD: 2.161720344213992 Å
- RMSD <2 Å QC pass: False

Interpretation: the IT1t redocking serves as a protocol QC/positive control. A more negative IT1t score than DINP would indicate that DINP is computationally weaker than the canonical co-crystallized antagonist under identical docking conditions. This comparison does not establish in-vivo binding.
