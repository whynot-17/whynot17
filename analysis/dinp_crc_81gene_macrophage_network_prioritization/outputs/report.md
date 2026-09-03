# Macrophage seven-gene PPI, topology and target-actionability prioritization

Generated: 2026-09-03T08:54:02.131555+00:00

## Frozen input and evidence separation

- Frozen candidates: `NEAT1, MMP9, TIMP1, STAT3, PTGER4, PTGES3, CXCR4`.
- STRING functional network: `5` high-confidence edges at combined score ≥0.700; no added interactors.
- STRING physical network: `1` high-confidence edge(s), reported separately from functional associations.
- STRING topology is descriptive. It does not establish direction, causality, or DINP binding.
- UniProt/PDBe/ChEMBL evidence is not merged into the STRING edge set; it is reported as target context.

## Network result

The high-confidence functional subnetwork contains `7` input nodes, `5` edges, `4` connected components, and a largest component of `4` nodes.

Observed functional edges:

- `MMP9 — TIMP1` (score `0.999`)
- `MMP9 — STAT3` (score `0.890`)
- `CXCR4 — MMP9` (score `0.844`)
- `CXCR4 — STAT3` (score `0.775`)
- `STAT3 — TIMP1` (score `0.772`)

## Role-based interpretation

- **Network bridge candidates:** highest topology among mapped proteins outside the direct prostaglandin class.
- **Direct prostaglandin nodes:** pathway/protein-context candidates; they need not have a STRING edge inside this seven-node induced subgraph.
- **Supporting/context nodes:** relevant expression or network context without enough evidence for the primary shortlist.
- **NEAT1:** non-protein state/regulatory node; not docking-eligible.

## Proposed docking/MD shortlist

| Rank | Gene | Role | Functional degree | Topology score | PDB structures | Best coverage | ChEMBL measured activity |
|---:|---|---|---:|---:|---:|---:|---:|
| 1 | **MMP9** | network bridge | 3 | 1.000 | 59 | 0.57 | 7141 |
| 2 | **STAT3** | network bridge | 3 | 0.955 | 8 | 0.69 | 5116 |
| 3 | **PTGER4** | direct prostaglandin node | 0 | 0.000 | 10 | 0.58 | 3405 |

### Reserve

`PTGES3` remains a direct prostaglandin-synthase reserve candidate because it has an experimental structure and direct pathway relevance, but it is isolated in the seven-node STRING functional network and has a much smaller ChEMBL activity record count than the main shortlist.

## Interpretation boundary

The shortlist prioritizes candidates for a later structural workflow. It does not show that DINP binds any target, that a STRING edge is a direct physical interaction, or that any gene mediates the epidemiologic association. Docking/MD should begin only after ligand identity, protein construct, binding site, and assay evidence are independently frozen.

Full STRING tables, pathway relation audit, node topology, UniProt/PDBe/ChEMBL context, request audit, and manifest are retained in `outputs/`.
