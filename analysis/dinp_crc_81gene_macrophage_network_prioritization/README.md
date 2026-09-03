# DINP–CRC macrophage seven-gene network prioritization

This analysis prioritizes the seven genes carried forward from the frozen
macrophage driver decomposition:

`NEAT1, MMP9, TIMP1, STAT3, PTGER4, PTGES3, CXCR4`

It is a transparent bridge from macrophage driver decomposition to a small
structural-target shortlist. The analysis does not claim that DINP binds any
protein or that any gene mediates the epidemiologic association.

## Frozen analysis rules

- STRING human network, species `9606`, combined score `>=0.700`.
- Functional and physical STRING associations are queried separately.
- No added interactors are used; topology is computed on the seven-node
  induced subgraph.
- Direct prostaglandin/arachidonic-acid relation requires an exact prior term
  (`GO:0006693`, `GO:0001516`, or `KEGG:00590`), or the prespecified direct
  receptor/synthase nodes `PTGER4` and `PTGES3`. A broad inflammatory-response
  annotation is context only and cannot promote a gene to a direct pathway
  node.
- UniProt, PDBe, and ChEMBL provide separate protein, structure, and measured
  activity context. They are not merged into the STRING edge evidence.
- Shortlist rule: two highest-topology protein nodes outside the direct
  prostaglandin class plus the most actionable direct prostaglandin node;
  `PTGES3` is retained as an explicit reserve.

## Re-run

From the repository root, use the configured Python runtime and run:

```text
python analysis/dinp_crc_81gene_macrophage_network_prioritization/run_dinp_crc_81gene_macrophage_network_prioritization.py
```

The script records request provenance in `outputs/api_request_audit.json` and
the exact source-release/endpoint context in `outputs/manifest.json`.

## Main outputs

- `outputs/report.md` — human-readable result and interpretation boundary.
- `outputs/functional_edges.csv` and `outputs/physical_edges.csv` — exact
  seven-node STRING edges.
- `outputs/functional_network_modules.csv` — Louvain modules.
- `outputs/network_target_evidence_matrix.csv` — node topology, pathway
  relation, UniProt, PDBe, and ChEMBL context.
- `outputs/docking_md_shortlist.csv` — proposed `MMP9`, `STAT3`, and `PTGER4`
  candidates for a later structural workflow.
- `outputs/docking_md_reserve_candidates.csv` — `PTGES3` reserve rationale.

The external databases are not copied into the repository; only the returned
small tables and request audit are versioned.
