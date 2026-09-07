# Expanded 18-gene DINP–CRC enrichment

- Method: one-sided hypergeometric ORA with custom backgrounds
- Query: the expanded 18-gene DINP exposure ∩ GeneCards CRC threshold set
- Sources: GO:BP, KEGG, Reactome
- Open Targets: not used
- Each background is a separate BH-FDR family

## Frozen inputs

| Analysis | Query | Background |
|---|---:|---:|
| Primary — GeneCards high-relevance background | 18 | 881 |
| Sensitivity — DINP multi-source background | 18 | 93 |

## Results

| Background | Input background | API effective domain | Returned terms | Global BH-FDR <0.05 | GO:BP | KEGG | Reactome |
|---|---:|---:|---:|---:|---:|---:|---:|
| Primary GeneCards 881 | 881 | 868 | 2433 | 106 | 99 | 0 | 7 |
| Sensitivity DINP 93 | 93 | 100 | 2433 | 0 | 0 | 0 | 0 |

## Top terms by background

### Primary GeneCards 881
- `GO:BP` inflammatory response — overlap 11/18; global BH-FDR=0.0002462
- `GO:BP` defense response — overlap 12/18; global BH-FDR=0.0006037
- `GO:BP` cellular response to chemical stimulus — overlap 15/18; global BH-FDR=0.000748
- `GO:BP` response to chemical — overlap 16/18; global BH-FDR=0.001269
- `GO:BP` response to stimulus — overlap 18/18; global BH-FDR=0.001982
- `GO:BP` fatty acid metabolic process — overlap 7/18; global BH-FDR=0.001982
- `GO:BP` regulation of cell population proliferation — overlap 14/18; global BH-FDR=0.001982
- `GO:BP` regulation of inflammatory response — overlap 8/18; global BH-FDR=0.001982
- `GO:BP` monocarboxylic acid metabolic process — overlap 8/18; global BH-FDR=0.002131
- `GO:BP` lipid metabolic process — overlap 9/18; global BH-FDR=0.002968

### Sensitivity DINP 93
- `GO:BP` positive regulation of response to stimulus — overlap 13/18; global BH-FDR=0.1258
- `GO:BP` positive regulation of signal transduction — overlap 10/18; global BH-FDR=0.1258
- `GO:BP` regulation of cell population proliferation — overlap 14/18; global BH-FDR=0.1258
- `GO:BP` cell population proliferation — overlap 14/18; global BH-FDR=0.1551
- `GO:BP` regulation of cell differentiation — overlap 9/18; global BH-FDR=0.1608
- `GO:BP` regulation of multicellular organismal process — overlap 14/18; global BH-FDR=0.1608
- `GO:BP` regulation of signal transduction — overlap 13/18; global BH-FDR=0.1608
- `GO:BP` positive regulation of signaling — overlap 10/18; global BH-FDR=0.1608
- `GO:BP` positive regulation of cell communication — overlap 10/18; global BH-FDR=0.1608
- `GO:BP` positive regulation of cell differentiation — overlap 7/18; global BH-FDR=0.1608

## Interpretation boundary

The frozen input background counts are 881 and 93. g:Profiler also reports an effective mapped domain after its identifier/domain processing; those API-reported values are retained separately and do not replace the frozen input counts.

This analysis tests over-representation only. It is direction-agnostic and does not show pathway activation, DINP causality, or mediation. The GeneCards primary background is the available archived ordinary CRC top-2000 reference rather than a full GeneCards export.

## Files

See the manifest for input hashes, g:Profiler version/timestamps, raw responses, and source-specific result files.
