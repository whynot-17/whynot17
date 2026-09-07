# GO:BP–Reactome aligned theme audit

Generated (UTC): 2026-09-07T14:34:01.945029+00:00

## Scope and statistical boundary

This is an additive cross-resource alignment layer. The original enrichment P values and global FDR values are unchanged. GO and Reactome terms are retained as source-specific evidence and are not treated as independent tests.

The 99 significant GO:BP terms were represented by the existing 43-term GO semantic-cleaning output. Reactome contributed 7 significant rows in the original 106-term table; the Reactome root term was excluded, leaving 6 informative terms.

## Result

- Clean GO representatives: **43**
- Reactome rows in original output: **7**
- Reactome root terms excluded: **1**
- Informative Reactome terms aligned: **6**
- Aligned source-specific terms: **49** (43 GO + 6 Reactome)
- Shared theme families: **9**

## Theme alignment

| Theme | GO reps | Reactome terms | GO–Reactome shared genes | Jaccard |
|---|---:|---:|---:|---:|
| Lipid / fatty-acid / eicosanoid metabolism | 11 | 4 | 6 | 0.4000 |
| Inflammation / defense / stimulus response | 4 | 0 | 0 | 0.0000 |
| Chemical / oxygen / abiotic stress response | 6 | 0 | 0 | 0.0000 |
| Apoptosis / programmed cell death | 2 | 0 | 0 | 0.0000 |
| miRNA transcription / metabolism | 2 | 0 | 0 | 0.0000 |
| Signaling / receptor / hormone response | 4 | 1 | 3 | 0.2727 |
| Matrix remodeling / cell migration | 1 | 1 | 3 | 0.3333 |
| Development / reproduction / circulation | 10 | 0 | 0 | 0.0000 |
| General metabolism / biosynthesis / storage | 3 | 0 | 0 | 0.0000 |

Reactome support is concentrated in three already-defined themes:

- **Lipid / fatty-acid / eicosanoid metabolism**: 4 Reactome terms, including lipid metabolism, fatty-acid metabolism, and specialized pro-resolving mediator biosynthesis.
- **Signaling / receptor / hormone response**: 1 Reactome term, Nuclear Receptor transcription pathway.
- **Matrix remodeling / cell migration**: 1 Reactome term, Activation of Matrix Metalloproteinases.

The other six themes have GO:BP support only in the current significant-term universe; this is not evidence that Reactome contradicts them.

## Driver interpretation

Driver contribution is descriptive repeated membership among retained/selected terms. It is not a gene-level test, causal attribution, or proof that a gene drives the enrichment.

| Theme | GO core drivers | Reactome core drivers |
|---|---|---|
| Lipid / fatty-acid / eicosanoid metabolism | SIRT1:8/11; PPARD:7/11; PTGS2:7/11; CYP1A2:6/11; CYP3A4:6/11; PTGS1:6/11 | CYP1A2:4/4; PTGS2:4/4; CYP3A4:3/4; PPARD:2/4; PTGS1:2/4; PPARG:1/4 |
| Inflammation / defense / stimulus response | ESR1:4/4; MMP9:4/4; PPARD:4/4; PPARG:4/4; PTGS2:4/4; RELA:4/4 | — |
| Chemical / oxygen / abiotic stress response | BECN1:5/6; MMP2:5/6; PPARG:5/6; SIRT1:5/6; PPARD:4/6; PTGS2:4/6 | — |
| Apoptosis / programmed cell death | BECN1:2/2; DKK1:2/2; ESR1:2/2; MMP9:2/2; PPARD:2/2; PTGS2:2/2 | — |
| miRNA transcription / metabolism | ESR1:2/2; PPARD:2/2; PPARG:2/2; RELA:2/2; STAT3:1/2 | — |
| Signaling / receptor / hormone response | PPARG:4/4; RELA:4/4; SIRT1:4/4; ESR1:3/4; BECN1:2/4; MMP9:2/4 | ESR1:1/1; PPARD:1/1; PPARG:1/1 |
| Matrix remodeling / cell migration | CXCR4:1/1; MMP2:1/1; MMP9:1/1; PPARD:1/1; PPARG:1/1; PTGS2:1/1 | MMP2:1/1; MMP9:1/1; TIMP1:1/1 |
| Development / reproduction / circulation | MMP2:10/10; PPARD:10/10; PTGS2:10/10; MMP9:9/10; TIMP1:9/10; ESR1:8/10 | — |
| General metabolism / biosynthesis / storage | CYP1A2:3/3; CYP3A4:3/3; MMP2:3/3; PPARD:3/3; SIRT1:3/3; BECN1:2/3 | — |

## Files

- `go_reactome_aligned_terms.csv`: term-level alignment with source preserved.
- `go_reactome_theme_alignment.csv`: nine-theme cross-resource summary.
- `go_reactome_driver_contribution.csv`: theme-by-gene GO/Reactome counts and fractions.
- `go_reactome_global_driver_contribution.csv`: global source-specific driver counts.

## Excluded Reactome root

Excluded term ID: `REAC:0000000` (`REACTOME root term`). It was not assigned to any biological theme and is not included in driver counts.
