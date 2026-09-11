# CURRENT ANALYSIS FREEZE — DINP/MCOP–CRC

Freeze date: 2026-09-11
Repository: whynot-17/whynot17
Baseline HEAD before this freeze: f871fb0a72ffffd7e517cbad1fffac8e3363b648

This file defines the publication-facing frozen state for the current DINP/MCOP–CRC project. Future analyses must not silently replace these datasets, denominators, thresholds, statistical units, or primary claims.

## 1. Frozen study chain

NHANES MCOP–CRC association
→ DINP-related gene space
→ CRC-associated gene space
→ DINP∩CRC overlap
→ STRING PPI / functional enrichment
→ bulk transcriptomic replication
→ cross-domain prioritization + 101-model ML + CytoHubba sensitivity
→ CEBPB
→ paired epithelial single-cell validation in GSE132465 and GSE294300
→ PPARG biological/structural bridge
→ docking / MD as structural plausibility only

## 2. NHANES association — frozen

Primary pooled complete-case model:
- N = 9,936
- CRC cases = 70
- exposure = log2(URXCOP), interpreted per doubling
- survey-weighted logistic regression
- adjusted for age, sex, race, BMI, smoking, PIR, and log2 urinary creatinine
- OR = 1.2455068
- 95% CI = 1.0775254–1.4396756
- P = 0.00331136
- 7 cycles: 2005–06 through 2017–18

Key sensitivities:
- age ≥40: OR = 1.22144; 95% CI = 1.04822–1.42329; P = 0.01084
- creatinine-normalized MCOP: OR = 1.2439983; 95% CI = 1.0747366–1.4399173; P = 0.0037899
- spline overall P = 0.002268; nonlinear P = 0.35828
- LOCO: all 7 omitted-cycle models OR > 1 and CI excludes 1
- formal between-cycle heterogeneity: F = 3.219; df = 6,109; P = 0.00598

Publication wording:
“Urinary MCOP was associated with higher odds of CRC in pooled NHANES analyses, and this association remained robust to exclusion of any individual survey cycle and adjustment for selected co-exposures, despite between-cycle heterogeneity.”

Do not claim causality or seven homogeneous replications.

## 3. DINP and CRC gene spaces — frozen

DINP-related genes:
- 286 unique human standardized genes
- multi-source integrated, not all direct binders
- internal evidence grades: A = 27, B = 258, C = 1
- 32 genes supported by ≥2 source families
- source contributions include CTD, GEO/toxicogenomics, experimental literature, T3DB, ChEMBL/PubChem

CRC-associated genes:
- 2,754 genes
- GeneCards relevance score ≥10 + OMIM + TTD union
- HGNC-standardized

Overlap:
- DINP ∩ CRC = 97 genes
- using 20,000-gene universe: OR ≈ 3.29
- P ≈ 2.61e-18
- expected ≈ 39.4 vs observed = 97
- enrichment factor ≈ 2.46

Preferred terminology: “DINP-related genes” / “associated genes”, not “direct DINP binders”.

Primary module:
- analysis/dinp_crc_target_mining_string_enrichment/

## 4. PPI and enrichment — frozen

STRING:
- human
- functional network
- score ≥700
- no added nodes
- 97 input genes
- 96 mapped
- 385 retained edges
- GPX1 unmapped

Deterministic PPI hub ranking:
1 IL6
2 TNF
3 TP53
4 AKT1
5 STAT3
6 IL1B
7 PPARG
8 MMP9
9 IFNG
10 SIRT1
11 IL4
12 RELA
13 IL13
14 ESR1
15 IL1A
16 PPARA
17 CEBPB
18 ADIPOQ
19 SREBF1
20 IL18

Enrichment:
- raw BP 1228 / MF 65 / CC 41 / KEGG 129
- 1462/1463 terms FDR < 0.05
- de-redundancy: Jaccard ≥0.5 OR overlap coefficient ≥0.8
- 110 representative terms retained: BP49 / MF4 / CC9 / KEGG48

Figure 2 publication state:
- 2A: DINP gene derivation + CRC gene derivation + 97 overlap
- 2B: STRING PPI
- 2C: GO/KEGG enrichment

## 5. Bulk transcriptomics — frozen publication-facing rule

Main bulk prioritization datasets:
- TCGA-COAD
- GSE10950
- GSE74602

GSE156355:
- small-cohort sensitivity only
- not used in the publication-facing Tier1 replication denominator
- contributed zero significant votes in the historical four-dataset screen

Frozen replication criterion:
- FDR < 0.05 in ≥2/3 main bulk datasets

Tier1 membership remains unchanged at 16 genes after denominator correction.

CEBPB:
- TCGA-COAD primary comparison: 288 tumor vs 41 TCGA solid-tissue normal
- mean: tumor 4.833 vs normal 3.676
- medians: tumor 4.8995 vs normal 3.534
- delta = +1.3655
- P = 7.88e-17
- FDR = 5.35e-16

GSE10950:
- 24 matched pairs
- tumor median = 9462.58
- normal median = 5707.79
- 19/24 pairs tumor-high
- paired Wilcoxon P = 0.00023925
- FDR = 0.00068785

GSE74602:
- 30 matched pairs
- CEBPB tumor-high
- current frozen figure P = 1.86e-09; exact table value should be read directly from the frozen output before manuscript insertion

Do not pool 41 TCGA normals with 308 GTEx normals for the primary CEBPB expression figure.

Primary module:
- analysis/dinp_crc_tcga_geo_expression_validation/

## 6. Cross-domain prioritization — frozen

Publication-facing Tier1:
- PPI Top20
- ≥1 representative pathway term
- FDR <0.05 in ≥2/3 main bulk datasets

Tier1 genes (n=16):
IL1B, IL1A, STAT3, ESR1, IL6, TP53, PPARA, CEBPB, MMP9, IL18, RELA, PPARG, IFNG, ADIPOQ, IL13, SREBF1

Historical cross-rank:
- CEBPB cross-rank = 11
- PPI hub rank = 17

Primary module:
- analysis/dinp_crc_ppi_pathway_transcriptomic_cross_rank/

Important:
- Tier1 membership is frozen at 16
- if the score uses a historical 4-dataset transcriptomic proportion, recompute the score/rank using the 3-dataset denominator before final manuscript tables
- membership itself does not change

## 7. 101-model ML — frozen

- 97 DINP–CRC input genes
- 95 measured genes
- TCGA train: 288 tumor / 41 normal
- external validation: GSE10950 + GSE74602
- GSE156355 excluded
- 101 models total
- 90 selective models
- 5-fold stratified CV
- seed = 20260909
- fold-local selector/scaler
- externally qualified selective models: 18
- qualification: AUC ≥0.75 in both external GEO cohorts

Stable ML definition:
- selection frequency ≥0.80 across all 90 selective models
- and ≥0.80 across qualified selective models

Stable ML17:
CD36, CNR1, INHBA, TIMP1, UBE2C, VEGFA, BMP2, FASN, ACACA, FABP1, NR3C1, RXRA, MKI67, CEBPB, ANGPT2, EPAS1, AGTR1

CEBPB:
- ML rank = 14
- priority score ≈ 0.81972
- reported qualified-model selection frequency should be read from the frozen stability output before final figure labeling; current working summary is 16/18 (89%)

Primary module:
- analysis/dinp_crc_101_ml_discovery/

Tier1 ∩ Stable-ML17:
- CEBPB only

Overlap is a prioritization rule, not a statistically significant enrichment claim.

Audit module:
- analysis/dinp_crc_101ml_tier1_overlap_audit/

## 8. CytoHubba sensitivity — frozen

Algorithms:
MCC, MNC, EPC, Degree, Closeness, Betweenness, Radiality, Stress

Consensus rule:
- Top25 per algorithm
- consensus = support in ≥7/8

18 consensus hubs:
AKT1, ESR1, IL1B, IL6, MMP9, NFE2L2, PPARG, RELA, SIRT1, ADIPOQ, CD36, CEBPB, IFNG, IL4, PPARA, STAT3, TNF, TP53

CEBPB:
- consensus rank = 12
- support = 7/8
- MCC rank 26, therefore not Top25 in MCC
- MNC 18
- EPC 7
- Degree 18
- Closeness 15
- Betweenness 20
- Radiality 12
- Stress 20

CytoHubba is a topology sensitivity analysis, not eight independent biological evidences.

Primary module:
- analysis/dinp_crc_cytohubba_8_algorithm_sensitivity/

## 9. CEBPB selection rule — frozen

CEBPB is not selected because it is the strongest PPI hub or strongest ML feature.

Frozen rationale:
- Tier1 n=16
- Stable ML n=17
- Tier1 ∩ Stable-ML17 = CEBPB only
- CytoHubba 7/8 provides supportive topology sensitivity
- bulk transcriptomics supports tumor-high CEBPB
- paired epithelial scRNA provides independent cellular validation

Preferred wording:
“CEBPB was not selected as the highest-degree PPI hub; rather, it emerged through cross-domain evidence convergence, particularly as the only gene shared between the integrated Tier-1 candidates and the stable machine-learning feature set.”

## 10. Clean single-cell primary validation — frozen

Formal primary cohorts ONLY:
- GSE132465
- GSE294300

Primary compartment:
- epithelial cells
- do not call these “malignant epithelial cells”
- no CNV-confirmed malignant-cell inference

Statistical unit:
- patient-level paired pseudobulk
- aggregate CEBPB UMI / aggregate library UMI ×1e6
- log1p
- paired Wilcoxon
- minimum 50 epithelial cells per sample
- CEBPB excluded from annotation features where reference mapping was used

GSE132465:
- matched n = 7
- 6/7 Tumor > Normal
- median delta log1p CPM = +0.6070179493
- mean delta = +0.7211271663
- P = 0.046875
- primary across-cohort BH-FDR = 0.046875

GSE294300:
- matched n = 10
- 7/10 Tumor > Normal
- median delta log1p CPM = +0.6654504445
- mean delta = +0.5719417333
- P = 0.037109375
- primary across-cohort BH-FDR = 0.046875

Formal primary total:
- 17 paired patients

Frozen files:
- analysis/dinp_crc_clean_CEBPB_scRNA_GSE132465_GSE294300/DINP_CRC_clean_CEBPB_patient_sample_pseudobulk.csv
- analysis/dinp_crc_clean_CEBPB_scRNA_GSE132465_GSE294300/DINP_CRC_clean_CEBPB_primary_epithelial_validation.csv
- analysis/dinp_crc_clean_CEBPB_scRNA_GSE132465_GSE294300/DINP_CRC_clean_CEBPB_epithelial_patient_audit_primary.csv
- analysis/dinp_crc_clean_CEBPB_scRNA_GSE132465_GSE294300/DINP_CRC_clean_CEBPB_epithelial_patient_audit_all.csv
- analysis/dinp_crc_clean_CEBPB_scRNA_GSE132465_GSE294300/DINP_CRC_clean_CEBPB_manifest.json
- analysis/dinp_crc_clean_CEBPB_scRNA_GSE132465_GSE294300/DINP_CRC_clean_CEBPB_report.md

## 11. Single-cell plotting rule — frozen

Accepted:
- patient-level paired pseudobulk plots for GSE132465 and GSE294300

Rejected for tumor-vs-normal inference:
- pooled-cell raincloud / pooled-cell density plots
- reason: cells are not independent statistical units and patients contribute highly unequal cell counts

Figure 4 frozen structure:
- 4A single-cell design
- 4B epithelial localization only
- 4C GSE132465 paired patient-level pseudobulk
- 4D GSE294300 paired patient-level pseudobulk
- no mandatory 4E

Myeloid:
- high CEBPB expression does not by itself establish a replicated tumor-associated increase
- any compartment-wide inference must use patient-level paired pseudobulk
- do not hide myeloid from a full-compartment descriptive plot

## 12. Quarantined / obsolete single-cell analyses

Do NOT use for publication-facing primary evidence:
- analysis/dinp_crc_multicohort_scRNA_CEBPB_CD36/
- analysis/dinp_crc_scRNA_GSE132465_CEBPB_CD36/

Reason:
- older modules contain prior design choices and/or unrelated sidedness logic
- GSE200997 is not part of the formal two-cohort validation family

Do not reintroduce GSE200997 into the primary analysis.

## 13. Figure state — frozen working layout

Figure 1:
- 1A workflow
- 1B NHANES model robustness forest
- 1C LOCO influence lollipop

Figure 2:
- 2A DINP/CRC gene derivation + 97-gene overlap
- 2B STRING PPI
- 2C GO/KEGG enrichment

Figure 3:
- 3A TCGA-COAD CEBPB expression
- 3B GSE10950
- 3C GSE74602
- 3D multi-evidence convergence to CEBPB

Figure 4:
- 4A single-cell design
- 4B epithelial localization
- 4C GSE132465 paired pseudobulk
- 4D GSE294300 paired pseudobulk

Main figures currently planned: 3 + 3 + 4 + 4 panels before docking/MD.

## 14. Structural bridge / docking — current frozen interpretation

PPARG is prioritized as the principal DINP structural/biological bridge because of integrated toxicologic and structural plausibility, not because it has the best docking score among all receptors.

Representative docking scores:
- PPARG best ≈ -7.889 kcal/mol
- PPARA ≈ -7.651
- RXRA ≈ -5.414
- PXR ≈ -6.779
- ESR1 ≈ -8.059

PPARG structure:
- 4XLD human PPARG LBD

DINP:
- representative structure CID 590836
- commercial DINP is a mixture; do not describe this representative ligand as the entire commercial mixture

Three-seed top-pose comparison:
- 09 vs 10 RMSD ≈ 9.96 Å
- 09 vs 11 RMSD ≈ 1.408 Å
- 10 vs 11 RMSD ≈ 10.133 Å

Interpretation:
- docking score repeatability is acceptable
- uniform pose reproducibility is not demonstrated

Native redocking QC:
- initial top1 score -8.816; raw RMSD 7.158 Å
- focused exhaustiveness 128: top1 raw RMSD 2.346 Å; top2 2.041 Å; mode15 1.358 Å
- aligned top1 whole 1.621 Å, core 1.467 Å, linker 2.033 Å
- aligned top2 whole 0.933 Å, core 0.874 Å
- aligned mode15 whole 1.132 Å, core 1.199 Å

Frozen wording:
- core geometry recovered
- residual flexible-region deviation remains
- no claim of perfect redocking
- docking supports structural plausibility, not binding proof

## 15. MD — current boundary

PPARG–DINP 100-ns MD is a structural plausibility module.
Planned/frozen readouts:
- protein backbone RMSD
- ligand RMSD
- RMSF
- Rg
- SASA
- H-bonds / contact persistence
- optional FEL / MM-PBSA if methodologically justified

Do not claim that MD proves binding or biological causality.

## 16. Claims that are prohibited in the frozen manuscript state

Do not write:
- DINP causes CRC
- CEBPB is a novel CRC gene
- all 286 DINP genes are direct binders
- every NHANES cycle replicated
- Tier1/ML overlap is statistically significant
- CEBPB is the strongest PPI hub
- CEBPB is the strongest ML feature
- CEBPB is epithelial-specific
- scRNA cells are independent replicates
- broad epithelial cells are CNV-confirmed malignant cells
- CytoHubba algorithms are eight independent biological evidences
- PPARG had the best docking score
- three-seed docking showed uniform pose reproducibility
- docking/MD proves binding
- native redocking was perfect

## 17. Version-control rule after this freeze

Any future analysis that changes:
- a primary dataset
- a denominator
- a threshold
- a sample count
- a statistical unit
- a gene-set membership
- a primary P/FDR
- a publication-facing figure value

must be committed as a new analysis version and must explicitly document what changed relative to this freeze.

No future script should recompute publication figures from raw inputs unless the figure is explicitly being re-audited. Publication figures should read frozen result tables directly.
