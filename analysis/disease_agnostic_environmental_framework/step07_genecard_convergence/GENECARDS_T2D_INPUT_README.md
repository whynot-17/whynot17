# T2D GeneCards inputs for the formal Step 7 rerun

The old CRC GeneCards export is not a valid substitute for this analysis. Step 7 uses two frozen, auditable T2D GeneCards inputs: a complete ordinary-query result as the primary set and the exact Disorders-scoped result as a strict sensitivity set. Both were captured from the publicly rendered GeneCards results table after the user completed the site verification step.

## Primary input

The primary query was:

```text
type 2 diabetes mellitus
```

The public page returned **20,554 genes in total**. The table was displayed at 100 rows per page and read across 206 public pages. The captured input is:

```text
input/t2d_genecards_primary_public_results.csv
```

The accompanying provenance record is:

```text
input/t2d_genecards_primary_public_provenance.json
```

The complete 20,554-row ranked list is the primary Step 7 gene set. Rank 100, 500, 1,000, and 2,000 are descriptive cutoffs only; they do not replace the full-list primary analysis.

## Strict sensitivity input

The strict sensitivity query was:

```text
[Disorders] "type 2 diabetes mellitus"
```

The public page returned **111 genes in total**. The table was displayed at 100 rows per page and read across two public pages. The captured input is:

```text
input/t2d_genecards_strict_111.csv
```

The accompanying provenance record is:

```text
input/t2d_genecards_strict_provenance.json
```

The captured fields are:

- `GeneCards_Rank` (or `Rank` / `#`)
- `GeneSymbol` (or `Symbol`)
- `RelevanceScore` when available
- `KnowledgeScore` when available

The 111-row exact phrase result is not used as the primary T2D gene set. The previously audited 689-row synonym OR query is not used as either the primary set or a replacement sensitivity set.

Then run the formal two-input script:

```text
python analysis/disease_agnostic_environmental_framework/step07_genecard_convergence/run_step07_t2d_genecard_convergence_v2.py \
  --ctd <CTD_chem_gene_ixns.tsv.gz> \
  --genecards-primary input/t2d_genecards_primary_public_results.csv \
  --genecards-strict input/t2d_genecards_strict_111.csv
```

The script records both input checksums and computes all 11 clusters uniformly. Primary and strict enrichment are corrected separately across the same 11-cluster family. The public-table capture uses only values visibly rendered by the GeneCards results tables; it does not use search-engine snippets, a private API, login-restricted export, or the CRC GeneCards set.

The query rules are frozen before using T2D-specific biological information: the ordinary complete list is primary, and the exact Disorders phrase is a high-specificity sensitivity analysis. Step 7 is post-firewall prioritization and does not modify the 29-test panel, the Step 6 robustness results, or the 11 exposure clusters.

The original `run_step07_t2d_genecard_convergence.py` and `t2d_genecards_public_disorders_*` files are retained only as historical preflight provenance. The formal rerun and its canonical outputs use `run_step07_t2d_genecard_convergence_v2.py` and the two inputs named above.
