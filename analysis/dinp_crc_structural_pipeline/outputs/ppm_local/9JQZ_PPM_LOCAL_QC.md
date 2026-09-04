# 9JQZ PPM 3.0 Local Orientation QC

## Provenance

- PDB ID: `9JQZ`
- Official source: `https://cggit.cc.lehigh.edu/biomembhub/ppm3_server_code.git`
- PPM source commit: `5784f088032e5b9fb35147815f67500b8cc94855`
- PPM exit code: `0`
- Input SHA256: `97459ccad1210e79eff1907742ca5468ee143d143705c7038fc5052efd039ba3`
- Oriented output SHA256: `a5a34e513975b726be7e818b47e492049b8163405da434c619f3373f27ea5bc7`

## Run settings

- Membranes: `1`
- Membrane type: `PMm` — mammalian plasma membrane
- Curvature: `planar`
- N-terminal topology: `out`
- Include heteroatoms: `no`
- PPM input mode: official CLI type `2` coordinate-file mode

## Coordinate QC

- Protein chains before/after: `['A']` / `['A']`
- Protein residues before/after: `281` / `281`
- Parsed protein coordinate atoms before/after: `2209` / `2209`
- Raw PDB coordinate records before/native-after: `2209` / `2704`
- PPM membrane-boundary dummy atoms retained in native output: `495`
- First residue before/after: `A:SER 16` / `A:SER 16`
- Last residue before/after: `A:ILE 296` / `A:ILE 296`
- Common C-alpha residues: `281`
- C-alpha RMSD after optimal rigid alignment: `0.000491668 Å`
- Internal pairwise-distance RMSD: `0.000408728 Å`
- Oriented N-terminal CA z-coordinate: `18.718 Å` (positive-z/outside side)
- Oriented C-terminal CA z-coordinate: `-15.745 Å` (negative-z side)
- PPM-reported TM segments: `7`; 7TM GPCR body: `PASS`
- N-terminal-out topology check: `PASS`
- Internal receptor geometry preserved: `PASS`
- Coordinate manipulation by this workflow: rigid-body orientation only; no loop completion, residue insertion, minimization, or renumbering.

## PPM-native outputs

- Hydrophobic/membrane thickness parsed from stdout: `32.4`
- Transfer energy parsed from stdout: `-69.3`
- TM segment records parsed from stdout: `[('A', '18', '41'), ('A', '48', '79'), ('A', '86', '123'), ('A', '130', '153'), ('A', '176', '212'), ('A', '220', '252'), ('A', '261', '291')]`
- Raw stdout: `9JQZ_PTGER4/ppm3_stdout.txt`
- Raw stderr: `9JQZ_PTGER4/ppm3_stderr.txt`
- Native PPM oriented output: `9JQZ_PTGER4/9JQZ_PTGER4_for_PPMout.pdb`

## Re-read checks

- BioPython reread: `PASS`
- MDAnalysis reread of temporary protein-only oriented view: `PASS (2209 atoms, 281 residues)`
- MDAnalysis reread of native mixed PPM file: `not compatible with this parser — index 2209 is out of bounds for axis 0 with size 2209`
- Native PPM output retained unchanged; its DUM records are membrane-boundary annotations, not receptor residues.
- One PTGER4 protein chain: `PASS`
- At least 250 protein residues: `PASS`
- N-terminal topology constraint submitted as `out`: `PASS`
- Membrane orientation conflict: `not detected by local post-run checks; inspect raw PPM output before downstream MD`

## Downstream boundary

This directory contains PPM orientation output only. No DINP pose was added, no POPC system was built, and no OpenMM/MD simulation was run.
