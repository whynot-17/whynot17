# T2D GeneCards input required for the full Step 7 run

The full Step 7 overlap/enrichment cannot be completed from the current local files because the repository does not contain a T2D GeneCards export. The old CRC export is not a valid substitute.

Please export the results of the exact GeneCards Disorders-scoped query used for this analysis, preferably:

```text
[Disorders] "type 2 diabetes mellitus"
```

Retain the unmodified export, including at least:

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

The script records the export checksum and computes all 11 clusters uniformly. It does not scrape GeneCards, use search-engine snippets, or substitute a CRC gene set.
