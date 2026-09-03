# Frozen 81-gene DINP–CRC single-cell subtype localization

This analysis drills the frozen 81-gene DINP–CRC intersection into source-labeled tumor epithelial cells and myeloid subtypes using the pinned Census-release H5AD described in `outputs/manifest.json`.

Run with the workspace Python environment:

```powershell
python run_dinp_crc_81gene_singlecell_subtype_localization.py
```

The primary unit of inference is the donor-level mean, with paired tumor-minus-normal contrasts. The score is the same 81-gene global z-score program used in the broad-compartment analysis; within-compartment standardization is a separate sensitivity analysis.

`ClusterFull` labels beginning with `Tumor` are treated as source-labeled tumor epithelial / malignant-candidate cells. They are not called definitively malignant because this analysis does not perform CNV-based or independent malignant-cell validation.

The source H5AD is intentionally not copied into the repository. Its release, path, byte size, and SHA-256 are recorded in `outputs/manifest.json`.
