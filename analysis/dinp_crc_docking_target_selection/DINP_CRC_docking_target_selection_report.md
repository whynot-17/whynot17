# DINP–CRC docking target selection

Generated: 2026-09-09T13:21:30.020273+00:00

## Decision

**PPARG is the proposed primary docking target.** CEBPB remains a downstream core mediator and is explicitly excluded from docking. The primary calculation should use a human PPARG ligand-binding-domain structure with a validated holo pocket; docking scores are not direct-binding proof.

## Shortlist

| Rank | Target | Role | 97 overlap | Pocket | Representative structure(s) |
|---:|---|---|---|---|---|
| 1 | PPARG | primary_docking_and_100ns_md | yes | high | 4XLD; 2PRG |
| 2 | PPARA | support_docking | yes | high | 6KB1; 3VI8 |
| 3 | RXRA | support_docking | yes | high | 1FBY; 1FM6 |
| 4 | NR1I2 | support_docking | yes | high | 8CH8 |
| 5 | ESR1 | support_docking_crc_network_control | yes | high | 1XPC |

## Why PPARG

PPARG is the cleanest bridge currently present in the project: it is in the 97 DINP–CRC overlap, is cross-supported by PPI/pathway/transcriptomic prioritization, has a well-characterized nuclear-receptor ligand-binding pocket, and the DINP-specific exposure evidence reports PPARG activity with early CEBPB upregulation. The important boundary is that the current DINP evidence is curated/predicted binding plus functional exposure evidence; a docking pose can support structural plausibility but cannot establish biochemical affinity.

## Calculation plan

1. Prepare PPARG LBD and define the grid from the co-crystallized ligand pocket; use a second PPARG structure for pose sensitivity.
2. Dock a documented representative DINP structure. DINP is a commercial mixture/branched-isomer class, so the exact PubChem/CAS-defined representative and its 3-D conformer must be recorded; PubChem reports a flexible molecule and no deposited 3-D conformer.
3. Run the full 100 ns explicit-solvent protein–DINP MD only for PPARG after pose and parameterization QC. The support candidates receive docking only.
4. Report docking score, pose contacts, pocket residues, ligand RMSD/occupancy, protein RMSD/RMSF, radius of gyration, SASA, hydrogen bonds, and replicate/structure sensitivity.

## Explicit exclusions

- **CEBPB**: downstream mediator; intrinsically/disorder-rich transcription factor and no defensible primary ligand pocket selected.
- **IL6/TNF/IL1B**: high network centrality but cytokine proteins are not clean small-molecule orthosteric targets for this question.
- **CD36**: membrane scavenger receptor with a more complex extracellular/membrane system and weaker direct DINP-pocket rationale.
- **AKT1/STAT3**: useful network hubs, but no specific DINP-to-CEBPB direct ligand bridge in the current evidence table.
- **NR1I3**: DINP-related receptor evidence exists, but the available CAR/RXR structural context is heterodimeric and the gene is not in the 97 overlap; retain as backup rather than primary support.

## Current execution status

Target selection and Vina docking are complete: PPARG was run with three seeds, and PPARA/RXRA/NR1I2/ESR1 were run as support screens. MD is intentionally deferred to the user's external machine and is not part of this docking-only package.

## Source records

- Project evidence: `outputs/DINP_CRC_overlap.csv`, `outputs/DINP_target_master.csv`, `outputs/DINP_CRC_PPI_pathway_transcriptomic_cross_rank.csv`, `outputs/DINP_CRC_STRING_cytohubba_8_consensus_hubs.csv`.
- RCSB PDB: [4XLD](https://www.rcsb.org/structure/4XLD), [2PRG](https://www.rcsb.org/structure/2PRG), [6KB1](https://www.rcsb.org/structure/6KB1), [1FBY](https://www.rcsb.org/structure/1FBY), [1FM6](https://www.rcsb.org/structure/1FM6), [6NX1](https://www.rcsb.org/structure/6NX1), [1XPC](https://www.rcsb.org/structure/1XPC).
- PubChem DINP: [CID 590836](https://pubchem.ncbi.nlm.nih.gov/compound/Diisononyl-phthalate), CAS 28553-12-0 / 68515-48-0.
