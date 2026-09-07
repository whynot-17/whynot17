# DINP × CRC disease-gene source cutoff sensitivity

This module audits whether the DINP–CRC overlap and the PTGER4/prostaglandin/arachidonate/PPAR/inflammatory story depends on the disease-gene source or on a selected rank cutoff.

Run from the repository root with the E-drive Python environment:

```powershell
E:\chatgpt\sc_env\Scripts\python.exe analysis/dinp_crc_disease_gene_source_cutoff_sensitivity/run_dinp_crc_disease_gene_source_cutoff_sensitivity.py
```

The complete CTD DINP set is 86 genes. GeneCards uses the locally preserved relevance-score ranking (`genecards_anywhere_crc_top2000.csv`); Open Targets uses the locally preserved `OpenTargets_score` in `crc_gene_matrix.csv`, sorted descending with a deterministic gene-symbol tie break. The script evaluates top 500/1000/2000 and reports enrichment against both the exact cutoff universe and the full ranked universe available in the local snapshot.

The output report is intentionally conservative: a signal that appears only at one cutoff is cutoff-sensitive and does not establish source-independent disease convergence.

Key outputs are written under `outputs/`: `CUTOFF_SENSITIVITY_REPORT.md`, `overlap_summary.csv`, `panel_enrichment_summary.csv`, `key_pathway_enrichment_key_terms.csv`, and `cutoff_sensitivity_overlap.png`. The full all-term Reactome/Hallmark table is generated locally as `key_pathway_enrichment.csv`.

