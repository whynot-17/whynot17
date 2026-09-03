# Frozen 81-gene DINP–CRC macrophage driver decomposition

This step decomposes the frozen 81-gene DINP–CRC intersection in macrophages before PPI/network analysis.

The primary unit of inference is the donor-level macrophage mean. Tumor–normal differences are paired within donor; individual cells are used only to report detection prevalence and are not treated as independent replicates. The primary family contains all 81 frozen genes, with BH-FDR applied across the paired gene-level tests.

The deterministic network-priority rule is:

1. positive paired tumor-minus-normal mean;
2. paired t-test BH-FDR < 0.05;
3. tumor-cell detection fraction >=25%; and
4. membership in at least one exact previously observed g:Profiler term in the prostaglandin/arachidonic-acid/inflammatory pathway lens.

The resulting genes are candidates for the next PPI step, not causal drivers. The source H5AD is not copied into the repository; release, size, SHA-256, and pathway-membership provenance are recorded in `outputs/manifest.json` and `outputs/pathway_membership_provenance.json`.

Run with the bundled workspace Python environment:

```powershell
python run_dinp_crc_81gene_macrophage_driver_decomposition.py
```
