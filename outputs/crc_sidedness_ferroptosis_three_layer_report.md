# CRC sidedness: three-layer ferroptosis screen

## Question
Which side shows the stronger ferroptosis-related state/liability before any SLC7A11 mechanism is assumed?

## Data availability note
The TCGA targeted-expression cache lacks the following requested genes: ACSL4, NCOA4, SAT1. Accordingly, TCGA metrics whose minimum gene-availability rule is not met are not estimable and must not be interpreted as null sidedness results.

## Layer 1 — transcriptomic ferroptosis state
Positive Right−Left for `ferroptosis_net_propensity` supports a more ferroptosis-prone right-sided transcriptional state; defense and driver components are reported separately.

| cohort | layer | metric | n_right | n_left | right_mean | left_mean | right_minus_left | standardized_mean_difference | welch_p | mannwhitney_p | welch_fdr_within_cohort_layer |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39582 | transcript_ferroptosis | ferroptosis_driver_score | 224 | 342 | 0.01304 | -0.008541 | 0.02158 | 0.06922 | 0.4172 | 0.3368 | 0.8615 |
| GSE39582 | transcript_ferroptosis | ferroptosis_defense_score | 224 | 342 | -0.001642 | 0.001076 | -0.002718 | -0.006265 | 0.9419 | 0.7547 | 0.9419 |
| GSE39582 | transcript_ferroptosis | ferroptosis_net_propensity | 224 | 342 | 0.01468 | -0.009617 | 0.0243 | 0.04824 | 0.5743 | 0.6443 | 0.8615 |
| TCGA-COAD | transcript_ferroptosis | ferroptosis_driver_score | 0 | 0 |  |  |  |  |  |  |  |
| TCGA-COAD | transcript_ferroptosis | ferroptosis_defense_score | 200 | 135 | -0.001507 | 0.002232 | -0.003739 | -0.01156 | 0.9181 | 0.9601 | 0.9181 |
| TCGA-COAD | transcript_ferroptosis | ferroptosis_net_propensity | 0 | 0 |  |  |  |  |  |  |  |

## Layer 2 — lipid-peroxidation liability
Positive Right−Left for `lipid_peroxidation_liability` supports greater right-sided PUFA/peroxide liability after subtracting antioxidant buffering.

| cohort | layer | metric | n_right | n_left | right_mean | left_mean | right_minus_left | standardized_mean_difference | welch_p | mannwhitney_p | welch_fdr_within_cohort_layer |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39582 | lipid_peroxidation | pufa_incorporation_score | 224 | 342 | 0.07547 | -0.04943 | 0.1249 | 0.2236 | 0.009996 | 0.008946 | 0.01999 |
| GSE39582 | lipid_peroxidation | peroxide_generation_score | 224 | 342 | -0.04008 | 0.02625 | -0.06634 | -0.1132 | 0.1845 | 0.145 | 0.1845 |
| GSE39582 | lipid_peroxidation | antioxidant_buffering_score | 224 | 342 | 0.07549 | -0.04945 | 0.1249 | 0.2412 | 0.006213 | 0.001869 | 0.01999 |
| GSE39582 | lipid_peroxidation | lipid_peroxidation_liability | 224 | 342 | -0.0578 | 0.03786 | -0.09566 | -0.1503 | 0.08669 | 0.03094 | 0.1156 |
| TCGA-COAD | lipid_peroxidation | pufa_incorporation_score | 0 | 0 |  |  |  |  |  |  |  |
| TCGA-COAD | lipid_peroxidation | peroxide_generation_score | 200 | 135 | -0.03195 | 0.04733 | -0.07928 | -0.1434 | 0.2092 | 0.4291 | 0.2092 |
| TCGA-COAD | lipid_peroxidation | antioxidant_buffering_score | 200 | 135 | 0.03876 | -0.05743 | 0.09619 | 0.2286 | 0.03733 | 0.02045 | 0.07467 |
| TCGA-COAD | lipid_peroxidation | lipid_peroxidation_liability | 0 | 0 |  |  |  |  |  |  |  |

## Adjusted sensitivity models
HC3 OLS uses sidedness plus stage and/or MMR/MSI only when those covariates are sufficiently available in the cohort.

| cohort | layer | metric | n | estimable | covariates | right_beta | right_p_hc3 | r_squared |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39582 | transcript_ferroptosis | ferroptosis_driver_score | 515 | True | right+stage+mmr_high | 0.01619 | 0.579 | 0.009611 |
| GSE39582 | transcript_ferroptosis | ferroptosis_defense_score | 515 | True | right+stage+mmr_high | -0.08343 | 0.03856 | 0.09202 |
| GSE39582 | transcript_ferroptosis | ferroptosis_net_propensity | 515 | True | right+stage+mmr_high | 0.09962 | 0.0389 | 0.04538 |
| GSE39582 | lipid_peroxidation | pufa_incorporation_score | 515 | True | right+stage+mmr_high | 0.09135 | 0.0785 | 0.0328 |
| GSE39582 | lipid_peroxidation | peroxide_generation_score | 515 | True | right+stage+mmr_high | -0.02679 | 0.631 | 0.002108 |
| GSE39582 | lipid_peroxidation | antioxidant_buffering_score | 515 | True | right+stage+mmr_high | -0.01512 | 0.7464 | 0.1679 |
| GSE39582 | lipid_peroxidation | lipid_peroxidation_liability | 515 | True | right+stage+mmr_high | 0.0474 | 0.4262 | 0.08211 |
| TCGA-COAD | transcript_ferroptosis | ferroptosis_driver_score | 0 | False |  |  |  |  |
| TCGA-COAD | transcript_ferroptosis | ferroptosis_defense_score | 327 | True | right+stage | -0.01006 | 0.7852 | 0.007216 |
| TCGA-COAD | transcript_ferroptosis | ferroptosis_net_propensity | 0 | False |  |  |  |  |
| TCGA-COAD | lipid_peroxidation | pufa_incorporation_score | 0 | False |  |  |  |  |
| TCGA-COAD | lipid_peroxidation | peroxide_generation_score | 327 | True | right+stage | -0.05033 | 0.4168 | 0.08138 |
| TCGA-COAD | lipid_peroxidation | antioxidant_buffering_score | 327 | True | right+stage | 0.08283 | 0.07896 | 0.0157 |
| TCGA-COAD | lipid_peroxidation | lipid_peroxidation_liability | 0 | False |  |  |  |  |

## Layer 3 — DepMap functional dependency
For individual defense genes, a negative Right−Left Chronos effect means stronger right-sided dependency. For the composite index, positive Right−Left means stronger right-sided collective dependency.

| confidence_set | metric | n_right | n_left | right_mean | left_mean | right_minus_left | mannwhitney_p_two_sided | interpretation | fdr_within_confidence_set |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| high | SLC7A11_gene_effect | 12 | 5 | 0.1005 | 0.1463 | -0.04578 | 0.3827 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.574 |
| high | GPX4_gene_effect | 12 | 5 | -0.1467 | -0.3188 | 0.1721 | 0.3284 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.574 |
| high | AIFM2_gene_effect | 12 | 5 | -0.01312 | -0.03826 | 0.02514 | 0.8788 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.8788 |
| high | GCH1_gene_effect | 12 | 5 | -0.0878 | 0.01047 | -0.09828 | 0.3827 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.574 |
| high | DHODH_gene_effect | 12 | 5 | -0.4695 | -0.6283 | 0.1589 | 0.2343 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.574 |
| high | ferroptosis_defense_dependency_index | 12 | 5 | -0.07294 | -0.07566 | 0.002729 | 0.799 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.8788 |
| high_plus_medium | SLC7A11_gene_effect | 13 | 7 | 0.1024 | 0.1373 | -0.0349 | 0.4378 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.6427 |
| high_plus_medium | GPX4_gene_effect | 13 | 7 | -0.1482 | -0.2448 | 0.09654 | 0.5356 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.6427 |
| high_plus_medium | AIFM2_gene_effect | 13 | 7 | -0.01933 | -0.06531 | 0.04598 | 0.6992 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.6992 |
| high_plus_medium | GCH1_gene_effect | 13 | 7 | -0.08481 | -0.006679 | -0.07813 | 0.4378 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.6427 |
| high_plus_medium | DHODH_gene_effect | 13 | 7 | -0.4751 | -0.8142 | 0.339 | 0.06749 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.405 |
| high_plus_medium | ferroptosis_defense_dependency_index | 13 | 7 | -0.06842 | 0.04236 | -0.1108 | 0.5356 | For individual Chronos effects, negative Right-minus-Left = stronger right-sided dependency. For composite index, positive Right-minus-Left = stronger right-sided dependency. | 0.6427 |

## Decision rule
Do not label either side as having 'stronger ferroptosis' from one score. A sidedness conclusion requires directional agreement across the net transcriptomic propensity, lipid-peroxidation liability, and functional dependency layers. Discordant layers are biologically informative and should be reported as such.

## Guardrail
These analyses establish ferroptosis-related state/liability, not direct cell death. Direct C11-BODIPY/lipid-ROS plus Fer-1 or Lip-1 rescue remains the functional endpoint.
