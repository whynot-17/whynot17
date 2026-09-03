# WGCNA variance-filter sensitivity audit

This audit compares the standardized, unbiased target-overlay workflow across
three prespecified natural-MAD input sizes. The frozen 81-gene DINP--CRC list
and seven macrophage driver list were never forced into network construction.
They were overlaid after each network was built.

## Results

| WGCNA input | Analysis | Selected genes | Non-grey modules | Soft power | Target genes in best-enriched module | Fisher P | Module BH-FDR |
|---:|---|---:|---:|---:|---:|---:|---:|
| 5,000 | all samples | 5,000 | 8 | 8 | 14 | 0.0006273 | 0.005018 |
| 5,000 | tumor only | 5,000 | 10 | 8 | 13 | 0.002023 | 0.020232 |
| 8,000 | all samples | 8,000 | 15 | 10 | 14 | 0.0003816 | 0.005724 |
| 8,000 | tumor only | 8,000 | 16 | 9 | 14 | 0.0004635 | 0.007416 |
| 10,000 | all samples | 10,000 | 22 | 10 | 14 | 0.0003427 | 0.007540 |
| 10,000 | tumor only | 10,000 | 19 | 10 | 13 | 0.001368 | 0.025998 |

The best-enriched module was `turquoise` in the 5,000-gene all-sample and
tumor-only networks, `blue` in the 8,000-gene all-sample network, and
`turquoise` in the remaining sensitivity networks. Module colors are local
WGCNA labels and are not treated as biologically identical across runs.

The target identities were also stable. The 5,000-gene all-sample module
contained 14 targets, the 8,000-gene all-sample module contained the same 14,
and the 10,000-gene all-sample module contained 14 with the same dominant
extracellular-matrix/prostaglandin-related membership pattern. The tumor-only
best modules contained 13, 14, and 13 targets, respectively.

## Macrophage proxy

The fixed 15-marker macrophage/myeloid panel had 14 genes available in the
GPL570 gene matrix in every run. Module associations were tested separately
within each run and corrected across non-grey modules for the macrophage-score
trait. In the 5,000-gene primary network, the target-enriched turquoise module
correlated with the macrophage proxy at `r=0.619` in all samples and
`r=0.634` in tumor-only samples (both module-level FDR < 2e-62). The module
label changes in the 8,000- and 10,000-gene sensitivity networks, but the
target-enrichment and macrophage-association outputs are retained separately
for transparent inspection.

## Interpretation

Across natural-MAD input sizes, the frozen target genes repeatedly concentrate
in a single top-enriched module, while the macrophage proxy identifies strong
myeloid-associated modules. This supports exploratory transcriptomic/network
convergence in GSE39582. It does not establish exposure causality, identify a
causal target, or prove that the bulk macrophage proxy is a deconvolved cell
abundance estimate.
