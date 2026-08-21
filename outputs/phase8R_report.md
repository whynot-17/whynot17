# Phase 8-R：Phenotype-first Drug X screen

## Current scope and conclusion

This phase starts from the fixed acquired OXA-R state and screens drug sensitivity before using R3 genes, old modules, regulatory status or literature novelty. GDSC1 and GDSC2 were analysed separately; PRISM and CTRPv2 were unavailable locally and were not treated as negative results.
Mapped GDSC CRC cell models: 40 unique models; GDSC1/GDSC2 mapping and drug coverage are recorded in `phase8R_cell_line_mapping.csv`. Primary model: top250 weighted state score, -LN_IC50 as sensitivity, n>=20 association floor and n>=30 priority.
Primary GDSC discovery gates: {'none': 482, 'TierB_moderate_candidate': 66, 'TierA_strong_collateral': 58}. Cross-database replication labels: {'Level3_single_dataset': 410, 'Reject_opposite_direction': 37, 'exploratory_same_direction_or_underpowered': 28, 'Level2_one_strict_other_same_direction': 22, 'Level1_cross_database_strict': 11}.
The stringent phenotype-only gate leaves 4 robust cross-database shortlist drugs: SN-38, AZD2014, JQ1, AZD1332. This is a phenotype-level screen, not a claim of acquired OXA-R experimental collateral sensitivity.

## Data availability

{
  "GDSC1": "available_local",
  "GDSC2": "available_local",
  "PRISM": "unavailable_local_not_negative",
  "CTRPv2": "unavailable_local_not_negative"
}

## Primary positive ranking (biological ranking frozen before annotation)

The table is restricted to positive multi-background associations; negative associations with small two-sided empirical q-values are not candidates.
```text
database          DRUG_NAME  median_background_rho  n_positive_backgrounds  n_backgrounds_available  global_empirical_q_value  primary_discovery_gate
   GDSC1         CRT0105446               0.379801                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1          CCT007093               0.273684                       3                        4                  0.000249 TierA_strong_collateral
   GDSC2         PRIMA-1MET               0.260729                       3                        4                  0.007172 TierA_strong_collateral
   GDSC1          Motesanib               0.259583                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1           EHT-1864               0.254687                       3                        4                  0.000249 TierA_strong_collateral
   GDSC2         Luminespib               0.248378                       3                        4                  0.008868 TierA_strong_collateral
   GDSC2 Obatoclax Mesylate               0.243845                       3                        4                  0.029461 TierA_strong_collateral
   GDSC2           PD173074               0.243681                       3                        4                  0.023393 TierA_strong_collateral
   GDSC1        Doramapimod               0.243334                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1         Zibotentan               0.241593                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1         Venotoclax               0.239501                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1        N23918-95-7               0.229423                       3                        4                  0.011917 TierA_strong_collateral
   GDSC1         Serdemetan               0.229146                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1              FS112               0.224206                       3                        4                  0.000249 TierA_strong_collateral
   GDSC2              AZ960               0.219690                       3                        4                  0.027481 TierA_strong_collateral
   GDSC1           T0901317               0.214337                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1          TANK_1366               0.214279                       3                        4                  0.015584 TierA_strong_collateral
   GDSC1           QL-XI-92               0.214128                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1         Idelalisib               0.213280                       3                        4                  0.000249 TierA_strong_collateral
   GDSC1               AZ20               0.211584                       3                        4                  0.000249 TierA_strong_collateral
```

## Frozen phenotype-only shortlist

```text
 biological_rank DRUG_NAME drug_key               biological_status            replication_level databases_passing  n_databases_passing  mean_median_background_rho  min_median_background_rho  mean_global_empirical_q_value  GDSC1_median_background_rho  GDSC2_median_background_rho  GDSC1_global_empirical_q_value  GDSC2_global_empirical_q_value  GDSC1_direction_consistent_models  GDSC2_direction_consistent_models  GDSC1_signature_top100_models  GDSC2_signature_top100_models
               1     SN-38     SN38 robust_cross_database_shortlist Level1_cross_database_strict       GDSC1;GDSC2                    2                    0.178960                   0.171422                       0.039988                     0.186497                     0.171422                        0.020455                        0.059520                                  6                                  5                              6                              6
               2   AZD2014  AZD2014 robust_cross_database_shortlist Level1_cross_database_strict       GDSC1;GDSC2                    2                    0.147900                   0.115439                       0.047762                     0.115439                     0.180362                        0.042181                        0.053343                                  6                                  5                              6                              6
               3       JQ1      JQ1 robust_cross_database_shortlist Level1_cross_database_strict       GDSC1;GDSC2                    2                    0.102852                   0.064947                       0.065847                     0.140757                     0.064947                        0.035189                        0.096506                                  5                                  6                              4                              6
               4   AZD1332  AZD1332 robust_cross_database_shortlist Level1_cross_database_strict       GDSC1;GDSC2                    2                    0.102058                   0.058877                       0.059963                     0.058877                     0.145240                        0.057970                        0.061956                                  5                                  5                              1                              5
```

## Direction and sensitivity conventions

GDSC LN_IC50 is converted to sensitivity as `-LN_IC50`; positive rho means a higher acquired OXA-R-like state score is associated with greater drug sensitivity. HCT116 is aggregated as the median of its three trajectories; the primary evidence unit is the four biological backgrounds.

## Boundaries

R3 genes, FAO/DHODH/ferroptosis/ERAD hypotheses, previous LINCS results and drug novelty were not used to rank drugs. Regulatory/non-oncology fields are post-ranking annotations. PRISM/CTRP are unavailable, not negative. Raw GDSC and DepMap files remain local and are not committed.

## Post-ranking regulatory and novelty audit

The four-drug biological shortlist was audited on 2026-08-21 after ranking freeze. No candidate met the approved/non-oncology Drug X criterion. SN-38 is retained as an irinotecan/TOP1 comparator; AZD2014 and JQ1 as oncology mechanism comparators; AZD1332 as a conditional Trk mechanistic lead, not as an approved repurposing drug. See `phase8R_drug_regulatory_annotation.csv`, `phase8R_novelty_audit.md` and `phase8R_go_nogo.md`.
