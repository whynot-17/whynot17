# DINP docking execution summary

Generated: 2026-09-09T13:23:52.643134+00:00

## Target decision

PPARG remains the only primary target for downstream 100 ns MD. CEBPB was not docked. PPARA, RXRA, NR1I2 and ESR1 were run as support docking candidates.

## Vina results

| Target | Successful runs | Best top-1 score (kcal/mol) | Mean top-1 | SD |
|---|---:|---:|---:|---:|
| ESR1 | 1 | -8.059 | -8.059 | 0.0 |
| PPARG | 3 | -7.889 | -7.8497 | 0.0284 |
| PPARA | 1 | -7.651 | -7.651 | 0.0 |
| NR1I2 | 1 | -6.779 | -6.779 | 0.0 |
| RXRA | 1 | -5.414 | -5.414 | 0.0 |

PPARG top-1 scores across three seeds were tightly clustered, which supports sampling repeatability for this fixed receptor/pocket setup. Cross-target ranking is descriptive only: different proteins, pocket sizes, protonation states and structural preparations make raw Vina energies non-comparable. The more negative ESR1 score therefore does not override the PPARG mechanistic selection.

## Pose/contact audit

PPARG same-structure seed pose RMSDs (unrestrained coordinate RMSD of the top pose) were: 20260909 vs 20260910: 9.96 Å, 20260909 vs 20260911: 1.4077 Å, 20260910 vs 20260911: 10.1325 Å.

Closest PPARG top-pose protein contacts (≤4.5 Å in at least one seed):

| Residue | Minimum distance (Å) | Seeds with contact |
|---|---:|---:|
| GLYA:284 | 3.127 | 3 |
| ILEA:341 | 3.260 | 3 |
| ARGA:288 | 3.268 | 3 |
| CYSA:285 | 3.307 | 3 |
| META:348 | 3.327 | 3 |
| LEUA:330 | 3.340 | 3 |
| ILEA:326 | 3.392 | 3 |
| SERA:342 | 3.429 | 3 |
| PHEA:264 | 3.459 | 3 |
| ILEA:262 | 3.473 | 3 |
| GLUA:259 | 3.535 | 3 |
| META:329 | 3.618 | 3 |
| SERA:289 | 3.652 | 3 |
| META:364 | 3.675 | 3 |
| LEUA:333 | 3.713 | 3 |
| ILEA:281 | 3.747 | 3 |
| ALAA:292 | 3.753 | 3 |
| LEUA:340 | 3.764 | 3 |
| LEUA:255 | 3.780 | 3 |
| ILEA:249 | 3.828 | 2 |

## Interpretation boundary

Docking supplies a pose-plausibility hypothesis only. It does not prove that DINP binds PPARG, does not resolve the commercial DINP isomer mixture, and does not prove PPARG causes CEBPB induction. MD is intentionally deferred to an external machine; before starting there, inspect the PPARG poses and validate the DINP ligand parameters.
