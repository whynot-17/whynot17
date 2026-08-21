# Phase 7B-R3：global-null validation and final shortlist

## Final conclusion

Across 56 CRC DepMap 23Q4 models, six OXA-R trajectories and four biological backgrounds, **no universal single-gene dependency was identified**. Phase 7C/8 were not run.
The Phase 7B-R2 set of 19 candidates is reclassified as **internally cross-method-stable candidates**, not robust vulnerabilities. After the shared global null and corrected partial-Spearman gate, the final shortlist contains **13 genes**.
If the shortlist is empty, this is a valid negative result: the single-gene convergence hypothesis does not meet the final calibrated criteria, and later work should move to pathway/complex-level convergence without forcing a Gene X.

## What changed in R3

1. Every null draw applies one shared permutation of the full CRC CRISPR row labels to all six trajectories and all genes. This preserves gene-dependency covariance and the observed relationships among the HCT116 trajectories.
2. The test statistic is directly the median vulnerability rho across four biological backgrounds; no Fisher combination of correlated background p-values is used.
3. Corrected partial Spearman residualizes both state score X and dependency Y on the same covariates: proliferation, global dependency and target-gene expression.

## R3 shortlist criteria

Global empirical q <=0.10; >=3/4 backgrounds positive; >=2 resampling-supported backgrounds; >=4/6 cross-method models internally stable; leave-HCT116-out median rho >0 with >=2/3 remaining backgrounds positive; corrected adjusted median rho >0.

## Candidate summary

```text
    gene  r3_rank  observed_global_T_median_vulnerability_rho  global_empirical_q_value  n_positive_backgrounds  n_resampling_supported_backgrounds  adjusted_partial_rho_median  final_shortlist_flag
 ADAMTS9        2                                    0.263095                  0.014372                       3                                   2                     0.296828                  True
    MOB4        6                                    0.250866                  0.014372                       3                                   2                     0.272619                  True
C16ORF86       14                                    0.319697                  0.014372                       3                                   1                     0.345346                 False
    TJP1       23                                    0.275866                  0.014372                       3                                   1                     0.256493                 False
    SOD2     1288                                    0.252381                  0.028001                       3                                   0                     0.235869                 False
```

Previous R2-only labels are available in `phase7bR3_previous_internal_stability_candidates.csv`; they are not carried forward as validated vulnerabilities unless they pass the R3 final gate.

## Corrected covariate audit

The audit covers trajectory top200 and convergent top500 selections. The absolute adjusted-minus-raw rho median is 0.043; both X and Y use the same covariate set for every gene.

## LODO context

Global rank stability inherited from the six-model audit: {'GSE119603': 0.9800271835969915, 'GSE42387': 0.7938646719447395, 'GSE77932': 0.7264383024616113}。GSE77932 and GSE42387 remain material perturbations; no candidate is called universal solely from the full-data rank.

## Mechanism audit

   gene  r3_rank  global_empirical_q_value  final_shortlist_flag
    FN1     2845                  0.457553                 False
HSP90B1     3240                  0.515412                 False
  ITGB1     5029                  0.686721                 False
   GPX4     5287                  0.705696                 False
   CPT2     5989                  0.751836                 False
   RRM2     6552                  0.780798                 False
  PSMB5     7163                  0.802853                 False
  DHODH     9508                  0.879667                 False
  EDEM1     9672                  0.884882                 False
  CPT1A    11408                  0.922193                 False
SLC22A5    13863                  0.961814                 False
    VCP    16157                  0.989208                 False
  DERL1    16469                  0.989658                 False
SLC7A11    16701                  0.991279                 False

DHODH, RRM2, CPT1A, CPT2, VCP and SLC22A5 do not regain a universal dependency signal. Meldonium remains a No-Go for broad OXA-R reversal.

## Files

- `phase7bR3_global_empirical_p.csv`: shared global-null p/q calibration.
- `phase7bR3_final_ranking.csv`: effect, global q, stability, LOBO and corrected covariate fields.
- `phase7bR3_covariate_sensitivity.csv`: symmetric partial-Spearman audit.
- `phase7bR3_previous_internal_stability_candidates.csv`: R2 19-candidate reclassification.
- `phase7bR3_manifest.json`: exact null, criteria and provenance.

Raw DepMap/GEO files remain local and are not committed. Phase 7C/8 remain deferred.