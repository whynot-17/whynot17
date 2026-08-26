# T2D GeneCards input used for the full Step 7 run

The old CRC GeneCards export is not a valid substitute for this analysis. The exact T2D Disorders-scoped result was obtained from the publicly rendered GeneCards results table after the user completed the site verification step.

The query was:

```text
[Disorders] "type 2 diabetes mellitus"
```

The public page returned **111 genes in total**. The table was displayed at 100 rows per page and read across two public pages. The captured input is:

```text
input/t2d_genecards_public_disorders_results.csv
```

The accompanying provenance record is:

```text
input/t2d_genecards_public_disorders_provenance.json
```

The captured fields are:

- `GeneCards_Rank` (or `Rank` / `#`)
- `GeneSymbol` (or `Symbol`)
- `RelevanceScore` when available
- `KnowledgeScore` when available

Then run:

```text
python analysis/disease_agnostic_environmental_framework/step07_genecard_convergence/run_step07_t2d_genecard_convergence.py \
  --ctd <CTD_chem_gene_ixns.tsv.gz> \
  --genecards <T2D_GeneCards_export.csv> \
  --genecards-query '[Disorders] "type 2 diabetes mellitus"'
```

The script records the input checksum and computes all 11 clusters uniformly. The public-table capture uses only the values visibly rendered by the exact scoped query; it does not use search-engine snippets, a private API, login-restricted export, or the CRC GeneCards set.

Because the exact public query returned only 111 ranked genes, K=500, K=1000, and K=2000 are numerically identical in this run. They are retained in the output for compatibility with the frozen analysis schema, and this limitation is recorded in the provenance file and report.
