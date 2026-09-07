# GO/KEGG enrichment of the 16-gene DINP–CRC GeneCards-threshold intersection

- Method: ORA with g:Profiler; custom CTD-DINP background
- Query: CTD DINP human genes ∩ GeneCards CRC `RelevanceScore >= 10`
- Background: all 86 CTD DINP human genes
- Open Targets: not used
- Multiple testing: BH-FDR across the GO:BP, GO:MF, GO:CC and KEGG family; source-specific FDR also retained

- Query genes: **16**
- Background genes: **86**
- Returned overlapping terms: **2471**
- Global BH-FDR <0.05 terms: **10**
- Significant terms by source: GO:BP=10, GO:MF=0, GO:CC=0, KEGG=0
- Source-specific BH-FDR <0.05 terms: GO:BP=25, GO:MF=0, GO:CC=0, KEGG=1

## Interpretation boundary

This is functional over-representation, not evidence of DINP causality, direction of regulation, or a validated mechanism. The GeneCards input is the available archived ordinary CRC top-2000 reference; the result is not presented as a full-ranking analysis.

## Files

- `enrichment_all_GO_KEGG.csv`
- `go_bp.csv`, `go_mf.csv`, `go_cc.csv`, `kegg.csv`
- `api_response.json`, `request_payload.json`, `manifest.json`
