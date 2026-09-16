# MiNP–PPARG 100 ns MD complete analysis

Trajectory: `run_20260915_seed20260917/production.dcd` (9,992 frames, 0.01–99.92 ns). Protein backbone fitting used the initial PPARG structure; contact distances used periodic minimum-image distances.

## Main trajectory behavior

- The first large displacement is immediate relaxation in the first 0.1 ns: core RMSD first exceeds 6 Å at 0.02 ns. From roughly 1–40 ns the core and COM remain in a lower platform.
- A second, transient excursion occurs around 45–80 ns, peaking at 75.50 ns: core RMSD 9.96 Å, core COM displacement 9.74 Å, and distance to the original pocket 6.86 Å. The excursion is brief and the ligand returns toward the pocket afterward.
- In the final 20 ns (80–100 ns), mean core RMSD is 5.00 ± 0.75 Å, core COM displacement 4.61 ± 0.80 Å, original-pocket minimum distance 3.21 ± 0.34 Å, and initial atom-level contact retention 13.9 ± 8.2%. Core-internal RMSD is 0.052 ± 0.012 Å, indicating that the aromatic core remains rigid while its position changes.
- Across 70–100 ns, mean core RMSD is 5.31 ± 0.86 Å, core COM displacement 4.99 ± 0.94 Å, original-pocket minimum distance 3.36 ± 0.47 Å, and initial atom-level contact retention 12.0 ± 7.7%. The 75.50 ns excursion is therefore a brief outlier within a later shifted ensemble, not a continuously increasing drift.

## Contacts and residency

The initial docking contact residues are ARG288, LEU330, LEU333, LEU340, ILE341, and SER342. Total protein-residue contacts stay around 17–19 per 5 ns bin; hydrophobic residue contacts stay around 10–13. Geometric hydrogen-bond-like contacts are sparse (about 0.02–0.32 per 5 ns bin).

The initial atom-level contact retention is low because the ligand changes which atoms/residues it contacts. Residue-level occupancy shows persistent late contacts with CYS285, HIS449, TYR327, PHE282, MET364, HIS323, GLN286, PHE363, SER289, VAL339, and ILE326, among others. LEU330 from the original pose remains frequent, but LEU333, LEU340, and SER342 largely disappear after the rearrangement.

Using the operational definition “core minimum distance to the six initial-contact residues ≤4.5 Å”, the ligand is inside the original pocket for 99.38% of 100 ns. It is outside for about 0.62 ns total; the longest continuous excursion is about 0.08 ns, with the farthest sampled distance at 75.50 ns.

## Protein and ligand flexibility

Protein backbone RMSD remains around 2.5–2.8 Å after equilibration, and the initial-pocket Cα RMSD remains around 0.8–1.1 Å. The pocket does not collapse. Protein-aligned ligand RMSF is elevated across the whole ligand because of pose changes, while core-aligned internal RMSF separates the effects: core atoms are about 0.05 Å, Branch A reaches about 1.5–4.1 Å, and Branch B is mostly below 1 Å except the terminal oxygens near 0.9–1.0 Å.

## Clustering and representative poses

Clustering the protein-aligned core in the final 30 ns with a 2 Å threshold gives two dominant interconverting states (49.7% and 43.0%) and several short-lived minor states. This supports a stable shifted ensemble rather than one perfectly fixed pose.

The representative overlay contains the initial pose, the largest 50–80 ns rearranged intermediate (75.50 ns), the final DCD frame, and the dominant last-30-ns cluster medoid.

## Interpretation boundary

The MD supports a PPARG-bound, dynamically rearranged MiNP ensemble with a rigid aromatic core, flexible side chains, persistent total protein contacts, and brief excursions from the original contact geometry. It does not support the stronger claim that the original docking pose is maintained unchanged for 100 ns, and it does not establish binding affinity or experimental activity.
