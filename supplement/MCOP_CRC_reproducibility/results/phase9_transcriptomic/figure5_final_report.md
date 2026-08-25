# Figure 5 PPAR contradiction audit and mechanism lock

1. Standard PPAR pathway definitions show MIXED but predominantly DOWN in these CRC epithelial paired donors: 6 significant-down and 1 other definition(s).

2. Our custom PPAR/NR result IS reproduced by most independent KEGG/Reactome/Hallmark definitions; frozen median delta=-0.419, P=4.29e-07.

3. The observed signal is MIXED: composition shifts coexist with within-state remodeling (1 state down, 1 state up at BH-FDR<0.05).

4. PPAR/NR remodeling IS PARTIALLY linked to inflammatory-stress programs at donor level, but not in the proposed simple inverse direction: delta_IL6_JAK_STAT3 rho=0.475, BH-FDR=0.0138; RELA/STAT3 regulon-delta links are not significant.

5. DINP/MiNP toxicology evidence provides candidate convergence but DOES NOT establish exposure-to-state causality; MCOP remains an exposure biomarker.

6. Figure 5 status: YELLOW — retain as a disease-state convergence figure with a dashed environmental bridge, not a causal mechanism figure.

## Why the apparent contradiction is not a single statistical disagreement

The 2022 single-cell CRC paper called PPAR signaling activated after selecting tumor-versus-normal epithelial DEGs, observing enrichment of PPAR-associated genes, and perturbing tumor organoids with PPAR inhibitors. It did not test our prespecified seven-receptor/nuclear-receptor expression score at the donor level. The same paper noted that several genes did not reproduce in TCGA and that SCD and ACSL4 could reverse direction. Our audit therefore separates receptor abundance, broad KEGG/Reactome lipid programs, peroxisomal metabolism and inferred regulon activity.

## Same-cohort definition comparison

| definition | median_delta_tumor_minus_normal | BH_FDR | direction |
| --- | --- | --- | --- |
| KEGG hsa03320 PPAR signaling | -0.275780549435201 | 2.7602395857684317e-05 | down |
| Reactome PPAR-alpha lipid regulation | -0.143645610946417 | 0.02172481760014 | down |
| Reactome peroxisomal lipid metabolism | -0.6281646271434682 | 4.6789258097608886e-08 | down |
| Reactome mitochondrial fatty-acid beta oxidation | -0.5562598699714346 | 1.3242242857813837e-08 | down |
| Hallmark fatty acid metabolism | -0.2354791184152618 | 0.0002194879396649 | down |
| Hallmark cholesterol homeostasis | 0.1189942817288315 | 0.2203658460348379 | up |
| Enterocyte metabolic differentiation | -0.640730598794914 | 1.5398836694657805e-07 | down |

| definition | n_paired_donors | median_delta_tumor_minus_normal | p_value | BH_FDR | direction | genes_present | coverage |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Frozen 7-gene PPAR/NR core | 36 | -0.4186007613938946 | 4.2922329157590877e-07 | 9.299837984144688e-07 | down | 7.0 | 1.0 |
| PPAR receptor-only module | 36 | -0.5290479514079232 | 7.114140316843987e-07 | 1.3211974874138832e-06 | down | 3.0 | 1.0 |
| NR partner module | 36 | -0.6613253010447325 | 2.6309862732887268e-08 | 8.550705388188362e-08 | down | 4.0 | 1.0 |
| KEGG hsa03320 PPAR signaling | 36 | -0.275780549435201 | 2.1232612198218703e-05 | 2.7602395857684317e-05 | down | 51.0 | 0.6623376623376623 |
| Reactome nuclear receptor transcription | 36 | -0.3001563561764845 | 4.886940587311983e-06 | 7.058914181672864e-06 | down | 23.0 | 0.4339622641509434 |
| Reactome PPAR-alpha lipid regulation | 36 | -0.143645610946417 | 0.0200536777847446 | 0.02172481760014 | down | 57.0 | 0.4830508474576271 |
| Reactome peroxisomal lipid metabolism | 36 | -0.6281646271434682 | 1.0797521099448204e-08 | 4.6789258097608886e-08 | down | 13.0 | 0.4642857142857143 |
| Reactome mitochondrial fatty-acid beta oxidation | 36 | -0.5562598699714346 | 2.0372681319713593e-09 | 1.3242242857813837e-08 | down | 19.0 | 0.5135135135135135 |
| Hallmark fatty acid metabolism | 36 | -0.2354791184152618 | 0.0001857205643318 | 0.0002194879396649 | down | 158.0 | 1.0 |
| Hallmark cholesterol homeostasis | 36 | 0.1189942817288315 | 0.2203658460348379 | 0.2203658460348379 | up | 74.0 | 1.0 |
| Enterocyte metabolic differentiation | 36 | -0.640730598794914 | 5.9226294979453093e-08 | 1.5398836694657805e-07 | down | 10.0 | 1.0 |
| DoRothEA PPARA regulon activity | 36 | -0.5808169088444368 | 3.958528395742178e-06 | 6.432608643081039e-06 | down |  |  |
| DoRothEA PPARD regulon activity | 0 |  |  |  | not_estimable |  |  |
| DoRothEA PPARG regulon activity | 36 | -1.526913436079709 | 5.820766091346741e-11 | 7.566995918750762e-10 | down |  |  |

## Epithelial composition and within-state audit

Composition paired summary (the CSV retains all donor-condition-state rows):

| state | n_paired_donors | paired_median_delta | paired_wilcoxon_P | paired_BH_FDR | paired_direction |
| --- | --- | --- | --- | --- | --- |
| enterocyte-like annotation | 36 | -0.0057970915565263 | 0.3791665282915346 | 0.3791665282915346 | down |
| other epithelial annotation | 36 | -0.078889006534182 | 9.082723408937454e-07 | 2.7248170226812363e-06 | down |
| secretory-like annotation | 36 | 0.1089474779014859 | 0.0043546251254156 | 0.0065319376881234 | up |

Within annotation states (minimum 20 cells per condition):

| definition | n_paired_donors | mean_delta_tumor_minus_normal | median_delta_tumor_minus_normal | p_value | direction | donor_consistency | cell_subtype | minimum_cells_per_condition | BH_FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| enterocyte-like annotation | 24 | -0.1745166109362 | -0.1912380626686144 | 0.0001754760742187 | down | 0.8333333333333334 | enterocyte-like annotation | 20 | 0.0005264282226562 |
| other epithelial annotation | 4 | 0.0779006510535851 | 0.1039512092367715 | 0.625 | up | 0.5 | other epithelial annotation | 20 | 0.625 |
| secretory-like annotation | 27 | 0.062689483153574 | 0.0684265511708368 | 0.0120766013860702 | up | 0.6666666666666666 | secretory-like annotation | 20 | 0.0181149020791053 |

The source annotations support enterocyte-like, secretory-like and other epithelial states only. No malignant label was invented.

## Donor-delta convergence

| x | y | n_donors | spearman_rho | p_value | BH_FDR |
| --- | --- | --- | --- | --- | --- |
| delta_PPAR_NR | delta_enterocyte_differentiation | 36 | 0.6447876447876448 | 2.188445762505848e-05 | 0.0003501513220009 |
| delta_PPAR_NR | delta_Fatty_acid_metabolism | 36 | 0.5559845559845561 | 0.0004304980415068 | 0.0022959895547033 |
| delta_PPAR_NR | delta_intestinal_epithelial_differentiation | 36 | 0.5590733590733591 | 0.0003935566677563 | 0.0022959895547033 |
| delta_PPAR_NR | delta_IL6_JAK_STAT3 | 36 | 0.4746460746460746 | 0.0034456622192419 | 0.0137826488769678 |
| delta_PPAR_NR | delta_Hypoxia | 36 | 0.4293436293436294 | 0.008972945847679 | 0.0256840255959869 |
| delta_PPAR_NR | delta_Inflammatory_response | 36 | 0.4257400257400257 | 0.009631509598495 | 0.0256840255959869 |
| delta_PPAR_NR | delta_stress_like_epithelial | 36 | 0.3593307593307593 | 0.0313687366795073 | 0.0716999695531597 |
| delta_PPAR_NR | delta_OXPHOS | 36 | 0.3492921492921493 | 0.0367918284413397 | 0.0735836568826795 |
| delta_PPAR_NR | delta_TNF_NFkB | 36 | 0.303989703989704 | 0.0714675379835732 | 0.1270534008596857 |
| delta_PPAR_NR | delta_EMT | 36 | 0.2777348777348777 | 0.100995614742612 | 0.1615929835881792 |
| delta_PPAR_NR | delta_E2F_targets | 36 | 0.2548262548262548 | 0.1336443400708548 | 0.1898311683457497 |
| delta_PPAR_NR | delta_G2M_checkpoint | 36 | 0.2494208494208494 | 0.1423733762593123 | 0.1898311683457497 |
| delta_PPAR_NR | delta_RELA_activity | 36 | -0.2200772200772201 | 0.1971400510822886 | 0.2426339090243552 |
| delta_PPAR_NR | delta_UPR | 36 | 0.0764478764478764 | 0.657658041299248 | 0.7516091900562835 |
| delta_PPAR_NR | delta_MYC_targets_V1 | 36 | -0.0445302445302445 | 0.796499518656049 | 0.8045709817484887 |
| delta_PPAR_NR | delta_STAT3_activity | 36 | 0.0427284427284427 | 0.8045709817484887 | 0.8045709817484887 |

The positive IL6-JAK-STAT3 and inflammatory-response delta correlations mean that donors with larger PPAR/NR losses do not show the largest inflammatory gains. RELA and STAT3 regulon-activity deltas are also not significantly correlated with the PPAR/NR delta. Thus the cohort supports parallel disease-state changes, but does not support a direct donor-level PPAR-low -> RELA/STAT3-high coupling or mediation arrow.

## Toxicology and causal boundary

| chemical | PMID | url | model | tissue_relevance | endpoint | direction | evidence_type | usable_in_main_figure | boundary | gene_or_pathway | tissue | single_chemical_or_mixture | directness | reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MiNP | 27551952 | https://pubmed.ncbi.nlm.nih.gov/27551952/ | human receptor transactivation/two-hybrid assays and primary human hepatocytes | hepatic/receptor assay; not colon | CAR2/PXR activation; weaker human PPAR activation and target-gene responses | activation | single-chemical experimental | no_non_colon_supplement_only | Shows receptor engagement at experimental concentrations, not CRC epithelial state or human-dose causality. | CAR/PXR/PPAR nuclear receptors | hepatic/receptor assay; not colon | single-chemical experimental | direct receptor/cell assay | Shows receptor engagement at experimental concentrations, not CRC epithelial state or human-dose causality. |
| MiNP | 23843199 | https://pubmed.ncbi.nlm.nih.gov/23843199/ | in-silico docking to human PPAR/RXR subtypes | no tissue | predicted binding | binding | computational | supplement_only | Docking is not evidence of cellular activation or colon relevance. | PPARA/PPARD/PPARG/RXR binding | no tissue | computational | in-silico only | Docking is not evidence of cellular activation or colon relevance. |
| MiNP | 35421560 | https://pubmed.ncbi.nlm.nih.gov/35421560/ | primary mouse granulosa cells | ovary; not colon | PPRE reporter and PPAR target genes | dose-dependent/nonmonotonic; mainly PPAR-gamma in this model | single-chemical experimental | supplement_only | Supports PPAR engagement but is species- and tissue-specific. | PPAR response element | ovary; not colon | single-chemical experimental | direct non-colon cell assay | Supports PPAR engagement but is species- and tissue-specific. |
| DINP | 31154059 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5356750/ | FITC allergic dermatitis mouse ear model | skin/immune; not colon | phospho-RELA and phospho-STAT3 | increased in co-exposure dermatitis context | animal co-exposure/context-dependent | supplement_only | Supports inflammatory signaling plausibility but not single-agent CRC epithelial causality. | RELA/STAT3 phosphorylation | skin/immune; not colon | animal co-exposure/context-dependent | contextual animal co-exposure | Supports inflammatory signaling plausibility but not single-agent CRC epithelial causality. |
| MCOP/MCIOP | 34478338 | https://pubmed.ncbi.nlm.nih.gov/34478338/ | 760 pregnant women; maternal urinary metabolite and placental RNA-seq | placenta; not colon | 18 associated placental transcripts/pathways | association; mixed | human observational transcriptomics | no | Does not establish MCOP as a direct molecular perturbagen; biomarker and tissue differ from CRC epithelium. | placental transcriptome | placenta; not colon | human observational transcriptomics | observational biomarker association | Does not establish MCOP as a direct molecular perturbagen; biomarker and tissue differ from CRC epithelium. |
| DINP/MiNP | 42398653 | https://pubmed.ncbi.nlm.nih.gov/42398653/ | mouse liver and in-vitro metabolic assays | liver; not colon | PPAR activation, beta-oxidation and lipid-metabolism transcriptomics | activation/remodeling at model-dependent doses | single-chemical experimental | no_non_colon_supplement_only | Recent exposure evidence strengthens general PPAR plausibility but does not reproduce the CRC epithelial suppression state. | PPAR/lipid oxidation | liver; not colon | single-chemical experimental | experimental non-colon tissue | Recent exposure evidence strengthens general PPAR plausibility but does not reproduce the CRC epithelial suppression state. |

The main-figure exposure bridge may show only dashed general nuclear-receptor/PPAR plausibility. It must not depict MCOP as a direct perturbagen, nor DINP/MiNP as a proven cause of the CRC epithelial state.

## Evidence-tier lock

| link | evidence_level | evidence_basis | figure_representation | status | locked_wording |
| --- | --- | --- | --- | --- | --- |
| CRC -> epithelial PPAR/NR down | E3 | directly observed in current human paired-donor data | solid | GREEN | Frozen core median delta=-0.419; FDR=9.3e-07. |
| CRC -> RELA/STAT3 inflammatory state up | E3 | directly observed in current human paired-donor data | solid | GREEN | Expression/regulon and state evidence; no exposure attribution. |
| PPAR/NR <-> differentiation/inflammatory state | E2 | donor-level cross-program convergence | thin solid | YELLOW | Association, not mediation. |
| Within-state PPAR remodeling | E2 | annotation-state paired analysis | mixed split arrows | YELLOW | 1 state down and 1 state up at FDR<0.05; composition also shifts. |
| DINP/MiNP -> PPAR/NR | E1 | external toxicology-supported candidate | dashed candidate label only; omit non-colon evidence nodes | YELLOW | 6 curated records but 0 direct colon/intestinal records. |
| MCOP/DINP -> CRC epithelial PPAR state | E0 | untested exposure-to-state bridge | gap/dotted | RED | No human colon perturbation or formal mediation; causal arrow prohibited. |

## Collision audit

| collision_level | query | audited_relevant_records | finding | exact_prior_complete_chain | novelty_wording_allowed |
| --- | --- | --- | --- | --- | --- |
| exact | DINP AND colorectal cancer AND PPAR | 0 | No direct DINP-PPAR-CRC mechanistic paper identified | False | targeted search did not identify an exact prior study |
| exact | MiNP AND colorectal cancer AND PPAR | 0 | No direct MiNP-PPAR-CRC mechanistic paper identified | False | targeted search did not identify an exact prior study |
| exact | MCOP AND colorectal cancer AND PPAR | 0 | No direct MCOP-PPAR-CRC mechanistic paper identified | False | targeted search did not identify an exact prior study |
| partial | DINP AND cancer AND PPAR | 1 | 2026 systematic carcinogenic-hazard review emphasizes rodent PPAR-alpha liver mode and limited human cancer evidence; PMID 42094681 | False | targeted search did not identify an exact prior study |
| partial | MiNP AND PPAR | 4 | Receptor/hepatocyte, ovarian and macrophage studies; none establish CRC epithelial mediation | False | targeted search did not identify an exact prior study |
| adjacent | colorectal cancer AND PPAR signaling single-cell | 4 | CRC PPAR literature is substantial but definition- and isoform-dependent | False | targeted search did not identify an exact prior study |
| background | phthalates AND PPAR | 100 | Broad crowded background; count is qualitative/lower-bound and not used for novelty claims | False | targeted search did not identify an exact prior study |
| exact | DINP AND colon epithelial AND PPAR | 0 | No exact exposure-epithelial-state study identified | False | targeted search did not identify an exact prior study |
| exact | MiNP AND colorectal cancer | 0 | No direct CRC paper identified | False | targeted search did not identify an exact prior study |
| exact | MiNP AND colon epithelial | 0 | No direct colon epithelial perturbation identified | False | targeted search did not identify an exact prior study |
| partial | phthalate AND colorectal cancer AND single cell | 0 | No study connecting measured DINP/MCOP exposure to CRC single-cell state identified | False | targeted search did not identify an exact prior study |
| partial | phthalate AND epithelial state AND colorectal | 0 | No exact chain identified | False | targeted search did not identify an exact prior study |
| adjacent | PPARG AND colorectal epithelial AND single cell | 3 | CRC receptor/state literature exists without DINP/MCOP human exposure | False | targeted search did not identify an exact prior study |
| adjacent | PPARA AND colorectal epithelial AND single cell | 3 | CRC receptor/state literature exists without DINP/MCOP human exposure | False | targeted search did not identify an exact prior study |
| adjacent | PPAR AND RELA AND STAT3 AND colorectal epithelial | 1 | Pathway crosstalk literature is adjacent; no exposure-to-state chain | False | targeted search did not identify an exact prior study |
| adjacent | PPAR AND stress-like epithelial AND colorectal cancer | 1 | State-remodeling literature is adjacent; no DINP/MCOP linkage | False | targeted search did not identify an exact prior study |
| background | nuclear receptor AND colorectal epithelial state | 10 | Broad disease biology only | False | targeted search did not identify an exact prior study |

## Final Figure 5 scientific wording

**CRC epithelial disease-state convergence:** a prespecified PPAR/nuclear-receptor expression module is reduced in paired CRC epithelium, whereas downstream lipid/PPAR pathway definitions and individual receptor regulons are modular and can differ in direction. DINP/MiNP toxicology provides non-colon nuclear-receptor plausibility only. The computational evidence supports a candidate exposure-to-state bridge, not a causal DINP/MCOP-PPAR-CRC mechanism.
