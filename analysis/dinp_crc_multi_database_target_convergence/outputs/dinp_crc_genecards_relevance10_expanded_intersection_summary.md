# Expanded DINP–CRC GeneCards threshold intersection

Generated: `2026-09-07T13:31:33.050816+00:00`

## Frozen rule

The exposure side is the 93-gene source-preserving DINP matrix from CTD, EPA CompTox/ToxCast/Tox21, and T3DB. The CRC side is GeneCards only, filtered at `RelevanceScore >= 10`. Open Targets is not used.

## Counts

| Quantity | Count |
|---|---:|
| Expanded DINP exposure genes | 93 |
| GeneCards input rows | 2000 |
| GeneCards rows with RelevanceScore ≥10 | 881 |
| Expanded DINP ∩ GeneCards ≥10 | 18 |
| Added by non-CTD layers | 2 |

## Added genes from the expanded exposure side

`CYP1A2, CYP3A4`

The archived GeneCards input is the ordinary CRC top-2000 reference; this output is reproducible within that archived input and is not presented as a full-ranking GeneCards result.

## Output

- `dinp_crc_genecards_relevance10_expanded_intersection.csv`
- `dinp_crc_genecards_relevance10_expanded_intersection_manifest.json`
