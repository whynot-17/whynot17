# DINP–CRC 81-gene enrichment background audit

This module re-audits GO/KEGG over-representation of the frozen 81-gene DINP–CRC intersection after concern that the prior primary universe (DINP exposure genes ∪ CRC disease genes) was too broad for an intersection-specific null model.

## Primary universe

The full frozen DINP exposure-gene set is the primary inferential background. Because 81 of 86 DINP genes already overlap the CRC disease-gene union, this test has limited discrimination but directly asks whether the CRC-overlapping DINP subset is functionally selective relative to all DINP-eligible genes.

## Sensitivity/descriptive universes

- CRC disease-gene union
- legacy DINP exposure ∪ CRC union
- g:Profiler annotated human genome (descriptive only)

## Run

```bash
python analysis/dinp_crc_81gene_go_kegg_background_audit/run_dinp_crc_81gene_go_kegg_background_audit.py
```

The script requires the frozen convergence outputs under `analysis/dinp_crc_multi_database_target_convergence/outputs/`.

Key output: `outputs/key_pathway_background_comparison.csv`.

Interpretation rule: enrichment that disappears under the DINP exposure universe must not be presented as CRC-convergence-specific enrichment.
