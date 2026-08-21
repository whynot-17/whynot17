# Phase 8-R2：GDSC drug-universe coverage audit

This audit does not re-rank drugs. It quantifies the opportunity set available to the locked Phase 8-R phenotype-first screen and separates coverage limitation from biological non-hit.

## Definitions

- `all_drug_entries`: distinct drug names present after CRC/expression mapping in the Phase 8-R association table.
- `any_trajectory_n20/n30`: at least one frozen trajectory has at least 20/30 CRC models after self-line exclusion.
- `three_trajectories_n20/n30`: at least three frozen trajectories meet the threshold.
- `all_four_trajectories_n20`: all four biological backgrounds have evidence at n>=20.
- `approved_nononcology_high_confidence`: PRISM metadata phase `Launched` plus a non-oncology context term and no oncology context term. Unmatched records are not counted.

## Summary

```text
              universe_definition  n_database_drug_entries  n_unique_drugs  n_approved  n_approved_nononcology_high_confidence  n_investigational  n_preclinical  n_withdrawn  n_unresolved  n_approved_oncology  n_approved_mixed_oncology_nononcology  n_prism_metadata_matched   database
                 all_drug_entries                      378             378          50                                       1                 58             37            0           233                   47                                      1                       145      GDSC1
               any_trajectory_n20                      334             334          40                                       0                 53             31            0           210                   38                                      1                       124      GDSC1
               any_trajectory_n30                      318             318          40                                       0                 52             31            0           195                   38                                      1                       123      GDSC1
           three_trajectories_n20                      334             334          40                                       0                 53             31            0           210                   38                                      1                       124      GDSC1
           three_trajectories_n30                      314             314          40                                       0                 51             30            0           193                   38                                      1                       121      GDSC1
        all_four_trajectories_n20                      334             334          40                                       0                 53             31            0           210                   38                                      1                       124      GDSC1
                 all_drug_entries                      286             286          51                                       1                 42             20            0           173                   48                                      2                       113      GDSC2
               any_trajectory_n20                      272             272          50                                       1                 39             19            0           164                   47                                      2                       108      GDSC2
               any_trajectory_n30                      270             270          50                                       1                 39             18            0           163                   47                                      2                       107      GDSC2
           three_trajectories_n20                      272             272          50                                       1                 39             19            0           164                   47                                      2                       108      GDSC2
           three_trajectories_n30                      270             270          50                                       1                 39             18            0           163                   47                                      2                       107      GDSC2
        all_four_trajectories_n20                      272             272          50                                       1                 39             19            0           164                   47                                      2                       108      GDSC2
             GDSC_union_all_drugs                      540             540          67                                       2                 72             43            0           358                   62                                      2                       182 GDSC_union
    GDSC_union_any_trajectory_n20                      495             495          57                                       1                 67             37            0           334                   53                                      2                       161 GDSC_union
    GDSC_union_any_trajectory_n30                      479             479          57                                       1                 66             37            0           319                   53                                      2                       160 GDSC_union
GDSC_union_three_trajectories_n20                      495             495          57                                       1                 67             37            0           334                   53                                      2                       161 GDSC_union
GDSC_union_three_trajectories_n30                      475             475          57                                       1                 65             36            0           317                   53                                      2                       158 GDSC_union
```

Clinical labels are audit annotations only and do not influence phenotype ranking. PRISM metadata is used here only to quantify GDSC opportunity-set coverage; biological replication is run independently.
Raw GDSC/PRISM files are not committed. Only derived audit tables, report and manifest are versioned.
