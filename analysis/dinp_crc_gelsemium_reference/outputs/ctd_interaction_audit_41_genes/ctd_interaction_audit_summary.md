# CTD interaction-type audit for the fresh 41-gene DINP–CRC intersection

Generated (UTC): 2026-09-07T23:31:26.663831+00:00

## Main result

The 41-gene intersection contains **5** genes with at least one CTD DINP-specific binding/interaction action, **11** additional genes with single-chemical functional-response evidence, **3** literature-associated genes without direct single-chemical support, **16** co-treatment-only genes, and **6** genes without a CTD record.

A CTD chemical–gene record is not automatically a direct target claim. The original action labels, single-chemical flags, co-treatment flags, and PubMed IDs are retained below so that direct interaction, functional response, and literature association remain separate.

| Evidence class | Genes |
|---|---:|
| single-chemical binding or interaction | 5 |
| single-chemical functional response without binding action | 11 |
| literature association without direct single-chemical support | 3 |
| co-treatment only | 16 |
| no CTD record | 6 |

## Binding/interaction genes

| Gene | DINP source flags | CTD binding records | Single-chemical binding records | Actions |
|---|---|---:|---:|---|
| **NR1I2** | CTD;ToxCast;T3DB | 2 | 2 | affects^binding;increases^activity |
| **PPARA** | CTD | 1 | 1 | affects^binding |
| **PPARD** | CTD | 1 | 1 | affects^binding |
| **PPARG** | CTD | 2 | 2 | affects^activity;affects^binding;increases^activity |
| **RXRA** | CTD | 1 | 1 | affects^binding |

## Interpretation boundary

Genes with CTD binding/interaction actions are exposure-side mechanistic candidates. They are not all equivalent: direct binding/interaction records are distinct from expression/activity response records and from co-treatment literature. GO/Reactome/KEGG recurrence is a separate pathway-context layer and is not used here to upgrade CTD evidence class.

## Files

- `ctd_interaction_type_audit.csv`: one row per fresh intersection gene with CTD evidence decomposition.
- `ctd_interaction_audit_summary.md`: reviewer-facing summary and binding subset.
