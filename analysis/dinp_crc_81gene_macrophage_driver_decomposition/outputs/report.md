# Frozen 81-gene DINP–CRC program: macrophage driver decomposition

Generated: 2026-09-03T07:24:49.691949+00:00

## Analysis boundary

- Source: Census release `2025-11-08` H5AD; source file is not copied into the repository.
- Macrophage cells: `20,280`; paired donors: `35`; frozen genes: `81`.
- Expression was summarized within donor and disease state before inference. Cells were not treated as independent replicates.
- `adata.X` was used; `adata.raw` was not used.

## Frozen prioritization rule

A gene is primary-statistically eligible when the paired tumor-minus-normal mean is positive, paired t-test BH-FDR is <0.05 across the 81-gene family, and tumor-cell detection fraction is >=25%. A gene is network-priority eligible only if it also belongs to at least one exact, previously observed prostaglandin/arachidonic-acid/inflammatory term. The top 10 are ranked by paired Cohen’s dz, then FDR; these are candidates for network analysis, not causal drivers.

- Primary-statistical eligible: `18` genes.
- Network-priority eligible: `7` genes.
- Top network-priority set: `7` genes.

## Top network-prioritized macrophage genes

| Rank | Gene | Mean Δ | Cohen dz | t-test BH-FDR | Tumor detection | Prior pathway lens |
|---:|---|---:|---:|---:|---:|---|
| 1 | **NEAT1** | 1.249 | 1.19 | 1.09e-06 | 84.8% | prostaglandin_AA_inflammation |
| 2 | **MMP9** | 0.929 | 1.02 | 1.25e-05 | 47.9% | prostaglandin_AA_inflammation |
| 3 | **TIMP1** | 0.617 | 0.78 | 0.000329 | 69.4% | prostaglandin_AA_inflammation |
| 4 | **STAT3** | 0.204 | 0.77 | 0.00036 | 44.6% | nuclear_receptor;prostaglandin_AA_inflammation |
| 5 | **PTGER4** | 0.101 | 0.57 | 0.00802 | 28.7% | prostaglandin_AA_inflammation |
| 6 | **PTGES3** | 0.255 | 0.55 | 0.0105 | 71.9% | prostaglandin_AA_inflammation |
| 7 | **CXCR4** | 0.307 | 0.44 | 0.0379 | 55.0% | prostaglandin_AA_inflammation |

## Interpretation

This decomposition indicates which members of the frozen 81-gene program are most suitable for the next macrophage-focused PPI/network step. It does not establish that any gene is caused by DINP exposure, that it mediates the epidemiologic association, or that expression change is specific to malignant biology.

The full 81-gene statistics, donor-level means/deltas, exact pathway-membership audit, source hash, and deterministic candidate set are retained in `outputs/`.
