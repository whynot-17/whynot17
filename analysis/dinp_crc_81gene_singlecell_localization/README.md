# Frozen 81-gene DINP–CRC single-cell localization

This analysis localizes the frozen 81-gene DINP–CRC intersection across four CRC compartments: epithelial, myeloid, fibroblast, and endothelial.

The source is the official Census-release source H5AD already archived on `D:`. The script uses `adata.X` in backed mode and does not copy the large H5AD into Git. Each of the 81 genes is standardized across eligible cells, and the program score is the mean gene-wise z-score. The inferential unit is the donor: tumor and normal means are compared within donor when both are available.

The analysis reports pooled and donor-level summaries, paired contrasts, paired deltas, cell-type summaries, exact gene mapping, source hash, and a manifest. It does not infer malignant status and does not make a causal DINP→expression claim.

Run:

```powershell
& C:\Users\21634\anaconda3\python.exe .\analysis\dinp_crc_81gene_singlecell_localization\run_dinp_crc_81gene_singlecell_localization.py
```
