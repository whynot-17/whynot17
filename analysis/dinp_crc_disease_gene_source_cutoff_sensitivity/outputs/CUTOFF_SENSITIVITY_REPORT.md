# DINP × CRC disease-gene source cutoff sensitivity audit

## Design

The query is the complete CTD DINP set (86 genes). CRC disease genes were evaluated from two locally preserved ranked sources: GeneCards relevance-score rank (top 2,000 available) and Open Targets score rank (15,611 genes with a score). Each source was tested at top 500, 1,000 and 2,000. For every cutoff, enrichment was calculated against the matching cutoff universe and again against the full ranked universe available for that source. Reactome and Hallmark gene sets were tested by hypergeometric over-representation with global Benjamini–Hochberg FDR within each source/cutoff/background family.

## Overlap and PTGER4

| source | top500 | top1000 | top2000 |
|---|---:|---:|---:|
| GeneCards DINP overlap | 8 | 19 | 41 |
| Open Targets DINP overlap | 8 | 11 | 29 |
| PTGER4 in GeneCards cutoff | no | no | yes |
| PTGER4 in Open Targets cutoff | no | yes | yes |

PTGER4 ranks in the available source rankings are GeneCards=1105 and Open Targets=602. Thus PTGER4 is absent from both top500 sets, enters Open Targets at top1000, and enters GeneCards only by top2000.

## Key panel results using the full available ranked universe

Each cell is `overlap/panel_size; panel FDR`; the query is the DINP overlap at that cutoff. This fixed-background view is the most comparable across cutoffs.

### GeneCards

| panel | top500 | top1000 | top2000 |
|---|---:|---:|---:|
| PGE2_synthesis_core | 1/5; FDR=0.0714 | 2/5; FDR=0.00376 | 5/5; FDR=1.69e-08 |
| prostaglandin_synthesis | 1/6; FDR=0.0714 | 2/6; FDR=0.00376 | 5/6; FDR=5.01e-08 |
| arachidonate | 1/20; FDR=0.155 | 2/20; FDR=0.0294 | 6/20; FDR=3.19e-06 |
| PPAR | 1/27; FDR=0.155 | 2/27; FDR=0.0313 | 3/27; FDR=0.0166 |
| inflammatory_response | 0/63; FDR=1 | 3/63; FDR=0.0302 | 5/63; FDR=0.0099 |

### Open Targets

| panel | top500 | top1000 | top2000 |
|---|---:|---:|---:|
| PGE2_synthesis_core | 2/6; FDR=2.07e-05 | 2/6; FDR=4.06e-05 | 2/6; FDR=0.000299 |
| prostaglandin_synthesis | 2/14; FDR=6.25e-05 | 2/14; FDR=0.000123 | 2/14; FDR=0.000897 |
| arachidonate | 2/59; FDR=0.000775 | 2/59; FDR=0.00151 | 2/59; FDR=0.00864 |
| PPAR | 2/118; FDR=0.00231 | 2/118; FDR=0.00447 | 2/118; FDR=0.0242 |
| inflammatory_response | 0/198; FDR=1 | 1/198; FDR=0.131 | 3/198; FDR=0.00864 |

At the term level, the Reactome prostaglandin-synthesis term is globally FDR-significant for Open Targets top500/top1000/top2000 (driven by PTGS1 and PTGS2), while GeneCards reaches global FDR significance only at top2000 (driven by HPGD, PTGES, PTGES2, PTGS1 and PTGS2). The PPAR and inflammatory **terms** do not reach the global all-pathway FDR threshold in this fixed-background comparison, even where the smaller prespecified panel test is positive.

## Interpretation

The PTGER4 finding is **cutoff-sensitive**, rather than stable across all prespecified cutoffs. This audit therefore does not support describing PTGER4 as a source-independent CRC convergence hit. If a pathway is significant only in a broad source or only after expanding to top2000, it should be reported as a feature of that selected disease-gene ranking and treated as hypothesis-generating. See `panel_enrichment_summary.csv` for the prespecified prostaglandin, arachidonate, eicosanoid, PPAR and inflammatory panels, and `key_pathway_enrichment.csv` for term-level FDR.

The result does not prove that Open Targets “diluted” a true signal: because PTGER4 is not present at top500 in either source and appears at different cutoffs, the simpler conclusion is that the original PTGER4 enrichment is not robust to disease-gene source/cutoff choice. The proper next step is independent cohort or cell-level validation, not selecting the cutoff that gives the preferred pathway.

## Files

- `overlap_summary.csv`: all overlap counts, PTGER4 rank/presence and overlap lists.
- `panel_enrichment_summary.csv`: pathway-panel overlap and hypergeometric/FDR values under both backgrounds.
- `key_pathway_enrichment.csv`: Reactome/Hallmark term-level results for key names (prostaglandin, arachidonate, eicosanoid, PPAR, inflammation and lipid).
- `cutoff_sensitivity_overlap.png`: compact overlap plot.

