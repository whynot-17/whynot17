# DINP cross tissue perturbation consensus

Run the complete analysis with the Python environment on `E:`:

```powershell
E:\chatgpt\sc_env\Scripts\python.exe `
  analysis\dinp_cross_tissue_perturbation_consensus\run_dinp_cross_tissue_perturbation_consensus.py
```

Input GEO matrices are kept in `E:\chatgpt\data\dinp_cross_tissue_perturbation_consensus` and are never copied into Git. The run uses the verified public datasets GSE158473 (mouse ovary), GSE313258 (mouse liver), GSE44076 (paired human CRC tumor/adjacent normal), and GSE37364 (human CRC normal/biopsy). A thyroid organoid row is retained in the manifest as `UNAVAILABLE` because no unique DINP accession was confirmed.

The DINP side keeps dose and time contrasts separate, maps mouse symbols to human canonical orthologues through g:Orth/Ensembl, and then counts independent datasets once each. Strict support requires at least two independent DINP datasets, significance in each, and no direction conflict. The legacy 81-gene intersection is only an annotation in `legacy_gene_audit.csv` and `old81_overlap_annotation.csv`.

Because both confirmed DINP matrices are processed expression (logCPM or TPM), the script uses Welch tests on the native processed scale and records this as secondary-quality evidence. Pathway files use a transparent preranked rank-sum implementation over the locally available Hallmark and Reactome GMTs; `pathway_collection_manifest.csv` explicitly records GO:BP and KEGG as unavailable when no local GMT is present.

Main outputs include per-contrast DEG tables, `dinp_gene_consensus.csv`, `crc_gene_consensus.csv`, `dinp_crc_directional_convergence.csv`, dose/time consistency tables, orthology and legacy audits, candidate tiers, eight figures, `ANALYSIS_AUDIT.md`, and `DINP_CROSS_TISSUE_CONSENSUS_REPORT.md`.
