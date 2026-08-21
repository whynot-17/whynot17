# Phase 1: FAO/carnitine program in OXA-resistant CRC

## Question

Is an FAO/carnitine transcriptional program a stable, cross-model state of oxaliplatin-resistant colorectal cancer?

## Scope of this first pass

This is a bounded transcriptomic screen, not a biochemical FAO-flux or causal-dependency test. We used a predefined FAO/carnitine gene set and compared resistant versus parental, resistant versus sensitive, or non-responder versus responder samples according to the design of each dataset.

## Data included

| Dataset | Role | Design | Main limitation |
|---|---|---|---|
| [GSE77932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77932) | Acquired OXA resistance | HCT116 and DLD1 parental cells versus resistant clones | One parental sample per background |
| [GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387) | Acquired OXA resistance | HCT116, HT29 and LoVo parental/resistant triplicates | Cell-line-specific resistance programs |
| [GSE119603](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119603) | Acquired OXA resistance | HCT116 parental versus HCT116oxR, three replicates each | Count matrix originally focused on ncRNAs |
| [GSE124808](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124808) | In-vivo resistance context | HCT116 versus two resistant-clone xenografts | Descriptive; n=1 parental xenograft |
| [GSE30011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30011) | Cross-sectional sensitivity context | Baseline expression in OXA-sensitive versus OXA-resistant CRC lines | Not a matched acquired-resistance model |
| [GSE83129](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83129) | Patient response | Metastatic CRC treated first-line with oxaliplatin; 21 responders, 12 non-responders | Clinical response is not identical to molecular OXA resistance |

## First-pass module results

The score is a within-dataset, gene-wise z-score module mean. Positive values indicate higher module score in the resistant or non-responder group.

| Context | Carnitine core | FAO core | Combined FAO/carnitine |
|---|---:|---:|---:|
| GSE77932 HCT116 resistant clones | +1.061 | +0.315 | +0.512 |
| GSE77932 DLD1 resistant clones | −0.188 | −0.112 | −0.136 |
| GSE42387 HCT116 | −0.021 | +0.015 | +0.030 |
| GSE42387 HT29 | −0.654 | −0.469 | −0.453 |
| GSE42387 LoVo | −0.060 | +0.132 | +0.121 |
| GSE119603 HCT116oxR | −0.149 | −0.052 | −0.108 |
| GSE124808 resistant xenografts | +0.589 | +0.085 | +0.139 |
| GSE30011 resistant versus sensitive lines | −0.259 | −0.313 | −0.333 |
| GSE83129 non-responder versus responder | −0.018 | +0.091 | +0.039 |

## Preliminary interpretation

The universal hypothesis is **not supported by this first pass**:

- The combined FAO/carnitine module was higher in 4 of 7 acquired-resistance model contrasts and lower in 3 of 7.
- The carnitine-core direction was higher in only 2 of 7 acquired-resistance contrasts and lower in 5 of 7.
- The strongest downward signal occurred in HT29 OXA-resistant cells and in the independent HCT116oxR count dataset.
- The patient cohort showed only a small, non-significant increase of the FAO core in non-responders; the carnitine core was slightly lower in non-responders.
- GSE30011 also showed a lower baseline FAO/carnitine score in cross-sectional OXA-resistant lines, but this is a different design and should not be pooled with acquired-resistance models.

## Current decision

Do **not** claim that FAO/carnitine dependency is a universal OXA-resistant CRC state. The current result is better described as:

> FAO/carnitine reprogramming is heterogeneous across OXA-resistant CRC models and may define a subset rather than a universal resistance state.

This does not yet eliminate Meldonium. It changes the hypothesis from a pan-OXA-resistance drug to a potential **biomarker-defined metabolic sensitizer for FAO-high OXA-resistant CRC**.

## Required next analysis before a go/no-go decision

1. Replace the small hand-curated gene set with independently sourced Reactome/Hallmark FAO and carnitine-pathway gene sets.
2. Run pathway-level scoring with the same gene universe across all compatible datasets.
3. Perform gene-level random-effects/sign-concordance summaries rather than pooling incompatible sample designs.
4. Test whether an FAO-high subgroup exists within OXA-resistant models and whether it is associated with OXA response, CMS4/metastatic context, or redox state.
5. Use GSE190609 only as a metastatic/redox context cohort, not as a direct OXA-resistant transcriptome.
6. Reassess Meldonium only if a reproducible FAO-high subgroup survives these checks.
