# Reciprocal exemplar focused robustness review

This review reads only the existing frozen primary, uniform robustness, and cycle-heterogeneity CSV files for the two pre-specified directional exemplars. No model was refit, no new sensitivity analysis or FDR was added, and no figure, literature, or mechanism analysis was performed.

## Primary endpoints

- URXP02: male thyroid_disease βM=0.122402 (SE 0.043606); female hypertension βF=0.104517 (SE 0.026705). The male-end interaction was -0.104123 (P 0.0228304, fixed-406 q 0.2439246); the female-end interaction was 0.072253 (P 0.0195787, fixed-406 q 0.2148370).

- URXUSN: male any_cancer_history βM=0.130305 (SE 0.060664); female asthma βF=0.090575 (SE 0.039344). The male-end interaction was -0.161470 (P 0.0299013, fixed-406 q 0.2759071); the female-end interaction was 0.117764 (P 0.0341868, fixed-406 q 0.2891636).

## Robustness diagnostics

- URXP02 male thyroid endpoint: creatinine-adjusted β=-0.104123; 8 LOCO fits, range -0.148290 to -0.071296, sign reversals=0; winsorized β=-0.107232; upper-1%-deleted β=-0.109388; above-LOD β=-0.106786; cycle heterogeneity P=0.1512617.
- URXP02 female hypertension endpoint: creatinine-adjusted β=0.072253; 8 LOCO fits, range 0.035857 to 0.100351, sign reversals=0; winsorized β=0.070223; upper-1%-deleted β=0.071555; above-LOD β=0.071229; cycle heterogeneity P=0.2521167.
- URXUSN male cancer-history endpoint: creatinine-adjusted β=-0.161470; 4 LOCO fits, range -0.199027 to -0.125793, sign reversals=0; winsorized β=-0.159175; upper-1%-deleted β=-0.111764; above-LOD β=-0.195491; cycle heterogeneity P=0.2498667.
- URXUSN female asthma endpoint: creatinine-adjusted β=0.117764; 4 LOCO fits, range 0.097580 to 0.135053, sign reversals=0; winsorized β=0.116656; upper-1%-deleted β=0.102713; above-LOD β=0.122802; cycle heterogeneity P=0.1661845.

## Interpretation

All four endpoint interactions retain their expected direction in the existing creatinine, LOCO, tail, and above-LOD diagnostics, with no LOCO sign reversal. Cycle heterogeneity should be read separately: a low heterogeneity P indicates cycle dependence even when direction is preserved. URXP02 is the cleaner primary exemplar because its endpoints are thyroid disease and hypertension; URXUSN remains a backup because its male endpoint is the heterogeneous `any cancer history` outcome.

These diagnostics support moving URXP02 first into the next pre-specified focused validation step. They do not convert either pair into a fixed-406 FDR-confirmed result.
