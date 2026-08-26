# Step 8E — T2D integrated evidence synthesis and flagship selection

- Status: **complete_integrated_evidence_profiles**
- Scope: four frozen Tier A exposure axes only; no upstream result was re-run or changed.
- Method: evidence-profile classification; no opaque 0–100 total score.

## Final classification

| Axis | Biomarker(s) | Classification | Pathway concentration | Transcriptomic support |
|---|---|---|---|---|
| cluster_11 | URXUUR | **Supported** | Moderate | Moderate, tissue-specific support |
| cluster_5 | URXP02 | **Flagship** | High | Moderate, tissue-specific support |
| cluster_6 | URXUBA;URXUSR | **Supported** | Moderate | Moderate, tissue-specific support |
| cluster_8 | URXUPB | **Exploratory** | Low | Limited and heterogeneous |

## Evidence profiles

| Axis | Epidemiology | Robustness | GeneCards | Pathway | Network | Transcriptomics |
|---|---|---|---|---|---|---|
| cluster_11 | Strong | Strong | Strong (q=0.0141) | Moderate (86 terms; 33 modules) | Strong (O/E=1.52; empirical P=0.000999) | Moderate, tissue-specific support (11/19 modules tested) |
| cluster_5 | Strong | Strong | Strong (q=0.00317) | High (552 terms; 125 modules) | Strong (O/E=4.57; empirical P=0.000999) | Moderate, tissue-specific support (3/4 modules tested) |
| cluster_6 | Strong | Strong | Strong (q=0.000121) | Moderate (790 terms; 123 modules) | Strong (O/E=1.61; empirical P=0.000999) | Moderate, tissue-specific support (12/24 modules tested) |
| cluster_8 | Strong | Strong | Strong (q=0.00976) | Low (219 terms; 40 modules) | Strong (O/E=2.02; empirical P=0.000999) | Limited and heterogeneous (20/50 modules tested) |

### Flagship: cluster_5 / URXP02

The most concentrated profile is cluster_5. Its compact pathway representatives repeatedly identify xenobiotics, xenobiotic metabolism, cytochrome P450 metabolism, and related chemical-response themes. The axis also has significant network enrichment and human transcriptomic support in a tissue-specific rather than universal pattern. This is a flagship positive demonstration of the framework, not a causal T2D mechanism.

### Supported axes: cluster_6 and cluster_11

Both axes retain strong epidemiologic, robustness, GeneCards, and network evidence. Their pathway themes and transcriptomic signals are interpretable but broader or more tissue-dependent than cluster_5, so they are retained as supporting discoveries rather than co-equal flagships.

### Exploratory axis: cluster_8

Cluster_8 remains a real FDR-supported and network-supported discovery, but its large gene-set input, broad RNA/protein/cell-cycle annotations, and heterogeneous tissue directionality limit a focused biological claim. Exploratory does not mean negative.

## Answers to the five prespecified questions

1. **Most concentrated evidence:** cluster_5 / URXP02, because its pathway signal is thematically focused on xenobiotic/CYP biology while retaining the other evidence layers.
2. **Stable but non-flagship:** cluster_6 and cluster_11; both have strong upstream and network support but broader or more tissue-dependent biology.
3. **Effect of tissue specificity:** it does not invalidate the positive branch, but it prevents a claim of one universal T2D transcriptional direction. T2D tissues remain explicitly separated.
4. **Permitted language:** multi-layer biological convergence, network-supported biological convergence, tissue-specific transcriptomic support, and prioritized environmental axis. Not permitted: causal mechanism, mediation, exposure-induced pathway activation, or universal T2D mechanism.
5. **Framework proof-of-concept:** yes, as an outcome-firewalled environmental screening framework that generated a frozen T2D positive branch and then enabled disease-specific prioritization. The biological branch remains hypothesis-generating rather than causal.

## Interpretation-risk boundary

Risk flags are provided in `t2d_step8e_interpretation_risks.csv`. The main recurring boundaries are cross-sectional temporality, reverse-causation potential, annotation-density bias, gene-set size effects, and tissue-specific transcriptomic heterogeneity. cluster_6 remains one biomarker axis despite containing two correlated tests.

## Reproducibility

All input paths and SHA-256 checksums are recorded in `STEP8E_MANIFEST.json`. This step reads only frozen Step 5–8D outputs and introduces no new data source.
