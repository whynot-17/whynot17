# DINP–CRC docking-only handoff

## Scope

This package contains docking only. Molecular dynamics is intentionally deferred to an external machine and is not required for this handoff.

## Recommended primary system

- Receptor: `prepared/PPARG_4XLD_rigid.pdbqt`
- Ligand: `prepared/DINP_CID590836.pdbqt`
- Search box: `prepared/PPARG_4XLD_box.txt`
- Recommended pose: `docking_runs/PPARG_4XLD_DINP_seed20260909.pdbqt`, mode 1, −7.889 kcal/mol
- Alternate PPARG pose: `docking_runs/PPARG_4XLD_DINP_seed20260911.pdbqt`, mode 1, −7.823 kcal/mol
- Outlier to inspect separately: seed 20260910; its top pose is approximately 10 Å from one of the other top poses by raw ligand-coordinate RMSD.

The recommended pose is the best-scoring PPARG seed and lies in the canonical LBD region, with close contacts around Cys285, Arg288 and Ser289. The seed-20260911 pose is retained as a sensitivity alternative. Do not interpret the Vina score as experimental affinity.

## Support docking systems

| Target | Receptor | Box | Best score (kcal/mol) |
|---|---|---|---:|
| PPARA | `prepared/PPARA_6KB1_rigid.pdbqt` | `prepared/PPARA_6KB1_box.txt` | −7.651 |
| RXRA | `prepared/RXRA_1FBY_rigid.pdbqt` | `prepared/RXRA_1FBY_box.txt` | −5.414 |
| NR1I2/PXR | `prepared/NR1I2_8CH8_rigid.pdbqt` | `prepared/NR1I2_8CH8_box.txt` | −6.779 |
| ESR1 | `prepared/ESR1_1XPC_rigid.pdbqt` | `prepared/ESR1_1XPC_box.txt` | −8.059 |

Raw cross-target scores are descriptive only and should not be used to replace the biological target-prioritization decision.

## Docking settings

AutoDock Vina 1.2.7; exhaustiveness 32; 20 output modes; energy range 5 kcal/mol; 8 CPU threads. PPARG was repeated with seeds 20260909, 20260910 and 20260911. All logs, modes and output PDBQT files are retained in `docking_runs/`.

## Ligand provenance

The ligand is PubChem CID 590836, a representative bis(7-methyloctyl) phthalate structure. DINP is a branched-isomer commercial mixture; the provided ligand is therefore a documented representative, not the full mixture. The exact ligand structure and conformer must be stated in any downstream MD report.

## Main caveat

PPARG is selected because it provides the clearest DINP→PPARG activity→CEBPB bridge in the current project and has a well-defined LBD pocket. Docking remains structural plausibility evidence; it does not prove direct biochemical binding or causal CEBPB regulation.
