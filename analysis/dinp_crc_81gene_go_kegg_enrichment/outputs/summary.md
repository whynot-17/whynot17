# 81-gene DINP–CRC intersection: GO and KEGG enrichment

## Scope

ORA was run on the frozen 81-gene intersection from the three-source DINP–CRC convergence analysis.
The primary custom universe is the union of the frozen DINP exposure-gene union (86 genes) and the frozen CRC disease-gene union (15,885 genes), yielding 15,890 unique symbols (the two source unions overlap substantially).
A CRC-union-only universe (15,885 genes) is reported as sensitivity analysis.

## Method

Human GO Biological Process (GO:BP), Molecular Function (GO:MF), Cellular Component (GO:CC), and KEGG annotations were queried through g:Profiler with `domain_scope=custom`.
The API was asked for all returned terms (`all_results=true`, threshold 1.0); hypergeometric raw P values were reconstructed from the returned effective domain, term size, query size, and overlap.
BH-FDR was then recomputed across the full pre-specified tested term family using the source-specific term counts in the API metadata; no-overlap terms remain in that denominator.

## Primary readout

Use `enrichment_primary_combined_exposure_crc.csv` for the primary analysis. The principal columns are `raw_p_hypergeom`, `BH_FDR_all_GO_KEGG`, and `BH_FDR_within_source`; `intersection_size` is the number of query genes overlapping each term.
GO branches and KEGG are also emitted as separate CSV files. Term counts and significant-term counts are recorded in `manifest.json`.

## Reproducibility

The run was performed against g:Profiler's live API; the exact request payload, API response, timestamp, and response SHA-256 are retained in `manifest.json` and the two `api_response_*.json` files. The g:Profiler source release is recorded in each analysis manifest.

This analysis is functional enrichment of the 81-gene intersection; it does not establish DINP causality or direction of regulation.
