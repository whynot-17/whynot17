# Analysis audit

- Primary DINP gene selection used all mapped transcriptomic genes and strict multi-dataset direction/significance gates.
- Legacy 81 genes were read only after consensus for annotation.
- No raw FASTQ/count matrix was committed.
- Processed GEO scales were tested with Welch/paired t-tests and marked secondary quality.
- g:Orth mapping is Ensembl-backed; unmapped genes are retained in the mapping table.
- Pathway method is preranked rank-sum over local Hallmark/Reactome GMTs; GO:BP and KEGG are not silently imputed when files are unavailable.
