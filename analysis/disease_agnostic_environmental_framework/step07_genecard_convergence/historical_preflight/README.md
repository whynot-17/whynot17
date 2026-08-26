# Deprecated GeneCards preflight artifact

These files document an earlier T2D GeneCards preflight based on the query:

```text
[Disorders] "type 2 diabetes mellitus"
```

It returned 111 rows because the field restriction and exact-phrase rule were too narrow for the formal Step 7 gene-set definition. This result was subsequently deprecated by methodological audit.

The 111-row set is retained only for provenance and audit traceability. It is not a sensitivity analysis and is excluded from the formal report, enrichment tables, figures, joint prioritization, and all Tier assignments. The formal Step 7 analysis uses the complete ordinary query `type 2 diabetes mellitus` captured in `../input/t2d_genecards_primary_public_results.csv`.
