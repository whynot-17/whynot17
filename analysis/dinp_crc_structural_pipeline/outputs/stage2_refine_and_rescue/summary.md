# Stage-2 docking refinement + rescue summary

## CXCR4 control-derived refinement
- Best IT1t redocking RMSD: 1.595 Å
- Protocol QC pass (<2 Å): True
- IT1t affinity under selected protocol: -5.766 kcal/mol
- DINP affinity under the same frozen protocol: -4.835 kcal/mol
- Selection of the protocol used only IT1t redocking RMSD, never DINP affinity.

## PTGER4 rescue
- Status: ok
- DINP affinity: -7.429 kcal/mol

Interpretation boundary: docking remains computational structural plausibility evidence. MD should be run only after manual inspection of the selected pose/receptor and, for CXCR4, preferably after protocol QC passes.