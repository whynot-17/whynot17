# Data and file manifest

- Source study: Qin et al., *Cell Reports Medicine* (2023), [PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10591051/), DOI 10.1016/j.xcrm.2023.101231.
- Processed expression data are deposited in the CNGBdb/CNSA project [CNP0004138](https://db.cngb.org/search/project/CNP0004138/). The expression matrices are not included here because of their size.
- The analyzed single-cell matrix contained 378,344 barcodes and 55,697 gene rows. Spatial data covered four post-treatment sections: ST-CR1, ST-CR2, ST-NR1, and ST-NR2.
- Technical libraries were combined by patient and treatment stage. All paired single-cell tests use one pre-treatment and one post-treatment sample per patient; cells are not treated as independent replicates.
- Applying the paper's stated basic QC thresholds retained 378,344 barcodes, while the article reports 333,285 high-quality cells. The public processed matrix and the article's final cell set do not fully match; see the report.
- The public barcode metadata used here do not map XFZ patient IDs to CR/PR/NR outcomes. Response-stratified single-cell tests could not be performed. The spatial CR/NR labels apply only to the four sections and do not provide a cohort-wide patient map.
