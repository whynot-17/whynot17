# Meldonium–CRC: OXA-resistant colorectal cancer drug repurposing

This repository contains the computational analysis for an oxaliplatin-resistant colorectal cancer (CRC) drug-repurposing project.

## Current analysis status

- Phase 1: cross-model OXA-resistance pathway stability
- Phase 3: pyrimidine, UPR/ER-stress and EMT module decomposition
- Phase 5–6: perturbation-signature reversal and unbiased drug screening
- Phase 7A–7B: GDSC pharmacogenomic validation and trajectory-conditioned DepMap dependency mapping
- Phase 7C: functional-module convergence across Reactome, Hallmark, curated modules, CORUM and co-essentiality modules
- Phase 8: module-conditioned GDSC pharmacological convergence

The current working model is intentionally not centered on Meldonium. Meldonium remains an exploratory/failed metabolic hypothesis unless independent pharmacogenomic evidence supports a carnitine-entry dependency.

## Reproducibility

Raw GEO, DepMap and GDSC files are excluded from Git because they are large downloaded datasets. The repository retains analysis scripts, small gene-set definitions, derived tables, reports and manifests. Raw-data locations and source-release assumptions are recorded in the phase manifests.

The main outputs are under [`outputs/`](outputs/), and analysis scripts are under [`work/scripts/`](work/scripts/).

## Scope and interpretation

This is an in-silico prediction project. DepMap projections and GDSC associations are not equivalent to paired parental/OXA-R functional experiments. Drug candidates require independent validation of indication, novelty, exposure feasibility and wet-lab response.
