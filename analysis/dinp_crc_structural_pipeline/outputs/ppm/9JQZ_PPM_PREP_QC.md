# 9JQZ PTGER4 PPM-input preparation QC

This file documents coordinate extraction only. No membrane orientation, loop modeling, residue insertion, energy minimization, PPM run, or MD was performed.

## Structure and source

- PDB ID: `9JQZ`
- Structure title: Structural Insights into Selective Antagonism Grapiprant and EP4 Prostaglandin Receptor
- Experimental method: ELECTRON MICROSCOPY
- Resolution: 2.65 Å
- Raw RCSB mmCIF: `D:\Codex\projectless-workspaces\2026-08-22\non\work\whynot17\analysis\dinp_crc_structural_pipeline\inputs\ppm\9JQZ_raw.cif`
- Raw mmCIF SHA256: `e24e107735e3b66235068309d05846162b05e67e45dc00b348adc7a7e1db5e98`

## Annotation-based PTGER4 chain resolution

- Selected entity ID: `2`
- Entity description: GFP-like fluorescent chromoprotein,Prostaglandin E2 receptor EP4 subtype
- Selected mmCIF label asym ID(s): `B`
- Selected author/PDB chain ID: `A`
- Resolution evidence: UniProt `P35408` / `PTGER4` reference plus entity description match = yes
- Annotated PTGER4 author residue span: `-1–321`
- Protein residue count retained: `281`
- First retained residue: `SER 16`
- Last retained residue: `ILE 296`
- N-terminal retained residue identity: `SER 16`

## Coordinate and structural sanity checks

- Coordinate atom count: `2209`
- Obvious 7TM GPCR receptor body: **yes** — The selected entity is explicitly annotated as PTGER4 and the entry keywords identify a GPCR/membrane protein; the retained coordinates also pass the alpha-helix geometry sanity check.
- Alpha-helix coordinate sanity: **yes** — 229 alpha-like backbone positions across 279 evaluated positions; 10 runs ≥5 residues
- Chain break / missing coordinate residue detected: **yes**
- Internal numbering gaps: `[]`
- BioPython re-read: **yes**
- MDAnalysis re-read: **yes**

## Removed content

- Other polymer chains removed: `['H', 'L']`
- Water removed: `[{'name': 'water', 'comp_id': 'HOH'}]`
- Detergent/lipid/ion/crystallization additive/co-crystal ligand records removed: `[{'name': 'Grapiprant', 'comp_id': 'A1ECR'}, {'name': 'water', 'comp_id': 'HOH'}]`
- Fusion/accessory content: PTGER4 was selected by entity/reference annotation and only its annotated author-residue span was retained; non-receptor fusion coordinates, if present, were excluded.

## Output

- PPM upload PDB: `D:\Codex\projectless-workspaces\2026-08-22\non\work\whynot17\analysis\dinp_crc_structural_pipeline\inputs\ppm\9JQZ_PTGER4_for_PPM.pdb`
- Output PDB SHA256: `97459ccad1210e79eff1907742ca5468ee143d143705c7038fc5052efd039ba3`
- Output is non-empty and legal PDB: **yes**

## Recommended PPM settings

- Number of membranes = 1
- Type of membrane = Plasma membrane (mammalian)
- Allow curvature = no
- Topology (N-ter) = out
- Input = Coordinate file
- Include heteroatoms = no

## Boundary

This PDB is prepared only for upload to PPM membrane-orientation software. PPM was not run, the coordinates were not oriented, and no MD system was built.
