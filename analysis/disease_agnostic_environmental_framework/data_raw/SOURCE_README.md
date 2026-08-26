# Step 1 raw-source pointer

The CTD chemical vocabulary used for the pre-disease Steps 1–4 is the official
download `CTD_chemicals.tsv.gz` from:

`https://ctdbase.org/reports/CTD_chemicals.tsv.gz`

The file is intentionally kept outside the Git worktree because raw database
archives are handled as local data. The exact local path and SHA-256 are in
`data_processed/run_manifest.json` and the pre-disease audit.

No disease outcome, disease gene list, GeneCards result, CTD chemical–gene
interaction file, or historical effect estimate is an input to Steps 1–4.
