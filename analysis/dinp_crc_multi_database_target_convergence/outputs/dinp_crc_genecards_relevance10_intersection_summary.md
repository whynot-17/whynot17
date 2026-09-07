# DINP–CRC GeneCards threshold intersection

- GeneCards rule: `RelevanceScore >= 10`
- Exposure-side input: CTD human DINP genes only
- Disease-side input: GeneCards only
- Open Targets: not used
- DisGeNET: not used
- Outcome statistics, P values, FDR, enrichment, and downstream ranking: not used

## Counts

- CTD DINP human genes: **86**
- GeneCards input rows: **2000**
- GeneCards rows with RelevanceScore >= 10: **881**
- CTD ∩ GeneCards (RelevanceScore >= 10): **16 genes**

## Scope note

The available GeneCards file is the archived ordinary CRC top-2000 reference. This thresholded result is therefore reproducible within that archived input; it is not presented as a full-ranking export until a full raw GeneCards file is available.

## Output

- `dinp_crc_genecards_relevance10_intersection.csv`
- `dinp_crc_genecards_relevance10_intersection_manifest.json`
