# MiNP–PPARG 100 ns RG/SASA audit

The completed MiNP–PPARG production trajectory was reanalyzed from `production.dcd` with a single, reproducible script. Samples are taken every 0.1 ns (1,000 samples across 0–100 ns).

## Definitions

- `protein_rg_A` and `ligand_rg_A`: mass-weighted radius of gyration of protein or ligand heavy atoms.
- `protein_sasa_complex_A2` and `ligand_sasa_complex_A2`: FreeSASA Lee–Richards solvent accessible area while the protein and ligand are in the complex.
- `protein_sasa_isolated_A2` and `ligand_sasa_isolated_A2`: SASA of each component calculated separately at the same frame.
- `protein_buried_by_ligand_A2` and `ligand_buried_A2`: isolated SASA minus complex SASA.
- `complex_sasa_A2`: total SASA of the protein–ligand complex.
- `initial_pocket_sasa_complex_A2`: SASA of protein heavy atoms within 4.5 Å of the initial six-atom MiNP core (ligand heavy-atom indices 12–17).

Hydrogens are excluded from all RG/SASA atom selections. SASA uses a 1.4 Å probe radius. Coordinates are PBC-adjusted by shifting the ligand to the nearest protein image before complex SASA is calculated.

## Window means

| window | protein RG (Å) | ligand RG (Å) | protein SASA in complex (Å²) | ligand SASA in complex (Å²) | ligand buried SASA (Å²) |
|---|---:|---:|---:|---:|---:|
| 0–10 ns | 19.402 | 4.022 | 13,674.8 | 30.4 | 500.0 |
| 10–50 ns | 19.417 | 4.059 | 13,660.9 | 38.2 | 491.5 |
| 50–80 ns | 19.362 | 4.182 | 13,685.4 | 21.6 | 503.1 |
| 80–100 ns | 19.303 | 4.203 | 13,457.6 | 20.3 | 502.4 |
| 0–100 ns | 19.376 | 4.121 | 13,629.0 | 28.8 | 498.0 |

The protein RG stays within a narrow range, while ligand RG increases modestly after the early relaxation. The late ligand SASA remains low and the buried ligand area remains high, consistent with a ligand that stays largely occluded in the PPARG pocket despite pose rearrangement. These are geometric exposure metrics and do not by themselves establish binding affinity.

Full frame-level and window-level tables are in `outputs/`.
