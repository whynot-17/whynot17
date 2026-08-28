# Step 10B-D — Disease knowledge-source replacement audit

Generated: `2026-08-28T09:24:42.804197+00:00`
Open Targets disease ID: `MONDO_0005148`; returned target rows: **9907**; unique approved symbols: **9634**.
GWAS Catalog author-reported gene symbols returned: **957**; associations returned: **8848**; status: **complete_association_collection**.

GeneCards is retained as the frozen reference. Open Targets and GWAS Catalog are not merged into a single disease-gene truth set: coverage, source-native evidence, and convergence are reported separately.

The GWAS Catalog API is deliberately fail-visible. Only source-native authorReportedGenes with Ensembl/Entrez identifiers are retained; no unsupported locus-to-gene inference is added. A timeout or empty response is recorded as such and never converted to an empty biological gene set.

Status: `complete_source_probe`. This audit does not change the frozen environmental panel, epidemiologic results, 11 clusters, or Tier assignments.