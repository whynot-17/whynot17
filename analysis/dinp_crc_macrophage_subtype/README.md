# DINP–CRC macrophage subtype analysis

`run_macrophage_subtype_analysis.py` performs a macrophage-only reanalysis of the existing local GSE144735 CRC single-cell cache. It streams the gzipped UMI matrix, excludes cDC and non-myeloid compartments, recomputes macrophage-only HVGs/PCA/neighbours/UMAP/Leiden, and evaluates the frozen 81-gene DINP–CRC program, PTGER4, and the fixed seven-gene driver set at donor level.

Run from the repository root with the E-drive scanpy environment:

```powershell
E:\chatgpt\sc_env\Scripts\python.exe analysis\dinp_crc_macrophage_subtype\run_macrophage_subtype_analysis.py
```

All tables, provenance, figures, and the integration report are written below `outputs/`. The source matrix and annotation are not copied into GitHub. The earlier 370,115-cell Census object used by broad localization was unavailable locally; `outputs/input_manifest.json` records this and no data were redownloaded. Thus the reported tumor/normal statistics are for the locally cached GSE144735 cohort and should not be pooled numerically with the prior Census run.
