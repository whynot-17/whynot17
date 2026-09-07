# DINP–CRC enrichment background audit

## Why this audit was required
The frozen query contains 81 genes, while the entire DINP exposure universe contains only 86 genes. Thus 81/86 (94.2%) of DINP-eligible genes already overlap the CRC disease union. A very large CRC/union/genome background can therefore make pathway structure already present on the DINP side appear extremely significant without demonstrating intersection-specific enrichment.

## Primary analysis
The primary inferential universe is now the full DINP exposure set (86 genes). This asks whether the CRC-overlapping subset is functionally selective relative to genes that could have entered the intersection from the DINP side.

## Sensitivity analyses
- CRC disease union: 15885 genes.
- Legacy DINP∪CRC union: 15890 genes.
- g:Profiler annotated human genome: descriptive functional annotation only.

## Interpretation rule
A term that is significant only against CRC/union/genome but not against the DINP exposure universe should be described as a feature of the 81-gene set or of DINP-related biology, not as a CRC-convergence-specific enrichment.

See `key_pathway_background_comparison.csv` for prostaglandin/arachidonic/eicosanoid/inflammatory/lipid terms across all four backgrounds.
