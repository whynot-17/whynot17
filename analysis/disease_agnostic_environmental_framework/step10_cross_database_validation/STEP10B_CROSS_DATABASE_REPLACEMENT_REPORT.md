# Step 10B — Cross-database robustness and annotation-bias stress test

Generated: `2026-08-28T09:41:23.863422+00:00`

## Environmental replacement (E)

The frozen 134-candidate set was queried separately against ChEMBL, BindingDB, and PubChem BioAssay. Observed source-layer coverage: ChEMBL human activity in 13/134 candidates; BindingDB human affinity in 0/134; PubChem CID-to-AID membership in 27/134.

These counts are evidence-source coverage indicators, not comparable biological scores. BindingDB absence is not interpreted as a negative; PubChem CID-to-AID membership is not converted into a human target count. ChEMBL/BindingDB/PubChem are not merged into a single chemical-gene edge list.

## Disease replacement (D)

The T2D concept resolved to `MONDO_0005148` in both source-native routes. Open Targets returned 9634 unique approved symbols (data release 26.06, API 26.6.3). GWAS Catalog returned 8848 trait associations and 957 author-reported gene symbols in the source-native association collection.

GeneCards remains the frozen reference. Open Targets and GWAS Catalog are reported as independent disease-knowledge layers with source-native fields retained; their results are not merged into a single disease-gene truth set.

## QC and interpretation boundary

QC checks passed: **6/6**. Exact API response hashes, query metadata, source snapshots, input hashes, and output hashes are in `STEP10B_MANIFEST.json`, `STEP10B_E_SOURCE_SNAPSHOT.json`, and `STEP10B_D_SOURCE_SNAPSHOT.json`.

This is a post-firewall source-replacement audit. It does not promote/demote candidates, change the 29-test family, recompute epidemiologic FDR, or select a flagship axis.