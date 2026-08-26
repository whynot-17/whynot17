# Step 8D — T2D transcriptomic module directionality

## Scope

Step 8D tests the frozen Step 8C STRING network modules in independent human
T2D transcriptomic datasets. The four Tier A axes and their 97 network modules
are fixed inputs. No new genes, disease-specific screening, or flagship
selection is performed here.

The primary unit is the biological sample in each public GEO series. Analyses
are performed separately by dataset and tissue; datasets are not pooled across
tissues or platforms in the primary analysis. The main result is module-level
directionality (T2D minus control), not single-gene significance.

## Prespecified data sources

The first-pass panel uses public processed expression matrices with explicit
human T2D and comparator labels:

- GSE23343 — liver, T2D versus normal glucose tolerance;
- GSE21340 — skeletal muscle, T2D versus control; replicate controls are
  excluded by the frozen metadata rule;
- GSE71416 — adipose tissue, obese diabetic versus obese non-diabetic;
- GSE25724 — isolated human pancreatic islets, T2D versus non-diabetic.

Series matrix files and NCBI GEO platform annotation files are downloaded from
the official NCBI GEO FTP endpoints. Accession, platform, sample-selection
rules, retrieval time, and SHA-256 checksums are recorded in the manifest.

## Expression processing

The series matrices are treated as processed expression values supplied by GEO.
For each platform, probe IDs are mapped to human gene symbols using the GEO
platform annotation table. Multiple probes mapping to one symbol are collapsed
by the median per sample. A module score is the mean of within-dataset
gene-wise z-scores for the module's mapped genes. A module is tested only when
at least three module genes are observed; the number of observed genes is
reported for every module–dataset result.

## Statistical analysis

For each dataset and module, the primary contrast is:

```text
module score in T2D samples − module score in comparator samples
```

Welch's two-sample t test is used as a dataset-level descriptive test because
the public first-pass panel has heterogeneous processed-array inputs and no
single harmonized covariate schema. Benjamini–Hochberg FDR is applied within
each dataset across its tested modules. The analysis reports delta, confidence
interval, P value, q value, sample counts, mapped gene counts, and direction.

Cross-dataset synthesis is descriptive: for each module it reports the number
of positive and negative dataset-level estimates, sign concordance, and median
delta. It does not pool unlike tissues or treat cells/genes as independent
replicates. A module's direction is not called pathway activation, exposure
causality, mediation, or T2D-specific mechanism.

## QC and interpretation gates

- All four datasets must have nonzero T2D and comparator samples after the
  frozen metadata rule;
- every module–dataset result must expose its mapped gene count and tested
  status;
- sample-level scores are retained for audit;
- dataset-level tests are corrected within dataset and no module is selected
  because it is attractive in one tissue;
- cross-dataset directionality is interpreted only with tissue and sample-size
  context;
- no transcriptomic result selects a flagship axis at this stage.

The outputs are a data-source audit, module-by-dataset directionality table,
module score matrix, cross-dataset synthesis, and a frozen manifest/report.
