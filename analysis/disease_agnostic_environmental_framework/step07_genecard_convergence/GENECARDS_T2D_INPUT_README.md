# T2D GeneCards inputs for the formal Step 7 rerun

The old CRC GeneCards export is not a valid substitute for this analysis. The formal Step 7 analysis uses one frozen, auditable T2D GeneCards input: the complete ordinary-query result. It was captured from the publicly rendered GeneCards results table after the user completed the site verification step.

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

## Deprecated historical preflight

An earlier preflight used the overly narrow query:

```text
[Disorders] "type 2 diabetes mellitus"
```

The public page returned **111 genes in total**. This result is retained only as historical audit material under `historical_preflight/`.

```text
historical_preflight/t2d_genecards_public_disorders_results.csv
```

The accompanying provenance record is:

```text
historical_preflight/t2d_genecards_public_disorders_provenance.json
```

The captured fields are:

- `GeneCards_Rank` (or `Rank` / `#`)
- `GeneSymbol` (or `Symbol`)
- `RelevanceScore` when available
- `KnowledgeScore` when available

The 111-row exact phrase result is **deprecated**: it is not a sensitivity analysis, does not enter the report, enrichment, figure, joint prioritization, or Tier assignment, and is not a replacement for the primary set. The previously audited 689-row synonym OR query is also not used.

Then run the formal primary-only script:

```text
python analysis/disease_agnostic_environmental_framework/step07_genecard_convergence/run_step07_t2d_genecard_convergence_v2.py \
  --ctd <CTD_chem_gene_ixns.tsv.gz> \
  --genecards-primary input/t2d_genecards_primary_public_results.csv
```

The script records the primary input checksum and computes all 11 clusters uniformly. The public-table capture uses only values visibly rendered by the GeneCards results table; it does not use search-engine snippets, a private API, login-restricted export, or the CRC GeneCards set.

The query rule is frozen before using T2D-specific biological information: the ordinary complete list is the sole formal primary set. Step 7 is post-firewall prioritization and does not modify the 29-test panel, the Step 6 robustness results, or the 11 exposure clusters.

The original `run_step07_t2d_genecard_convergence.py` and the files under `historical_preflight/` are retained only as deprecated preflight provenance. The formal run and canonical outputs use `run_step07_t2d_genecard_convergence_v2.py` and the primary input named above.
