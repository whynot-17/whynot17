# Phase 8：module-conditioned pharmacological convergence

## 定义

- 只使用 Phase 7C 正向模块：高 OXA-R trajectory score 应对应更高药物敏感性，即 state score 与 LN_IC50 的 rho < 0。
- T1：6/6 trajectory 方向一致；T2：5/6 trajectory 方向一致。
- 药物排序同时要求至少 2 个 GDSC dataset、至少 3 条 trajectory、负 rho 比例≥0.60，避免单一背景伪阳性。
- 这是 phenotype-level pharmacogenomic prioritization，不等于 paired OXA-R CRISPR 或临床有效性。

Focus modules: 29；pairwise records: 89352；drug-module records: 14906；robust candidate records: 596。

## Robust candidates

| Tier | Module | Drug | n pairs | datasets | trajectories | median rho | negative fraction | score | target |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| T1_universal | COESSENTIALITY::1400::biological process:endocytic recycling | Cytarabine | 12 | 2 | 6 | -0.179 | 0.67 | 0.119 | Antimetabolite |
| T1_universal | COESSENTIALITY::1400::biological process:endocytic recycling | AZD1332 | 12 | 2 | 6 | -0.158 | 0.67 | 0.105 | NTRK1, NTRK2, NTRK3 |
| T1_universal | COESSENTIALITY::1400::biological process:endocytic recycling | PD173074 | 12 | 2 | 6 | -0.142 | 0.67 | 0.095 | FGFR1, FGFR2, FGFR3 |
| T1_universal | COESSENTIALITY::1400::biological process:endocytic recycling | Navitoclax | 12 | 2 | 6 | -0.091 | 0.83 | 0.076 | BCL2, BCL-XL, BCL-W |
| T1_universal | COESSENTIALITY::1400::biological process:endocytic recycling | Lestaurtinib | 12 | 2 | 6 | -0.105 | 0.67 | 0.070 | FLT3, JAK2, NTRK1, NTRK2, NTRK3 |
| T1_universal | COESSENTIALITY::1400::biological process:endocytic recycling | ZM447439 | 12 | 2 | 6 | -0.031 | 0.67 | 0.020 | AURKA, AURKB |
| T1_universal | COESSENTIALITY::1400::biological process:endocytic recycling | Rucaparib | 12 | 2 | 6 | -0.028 | 0.67 | 0.019 | PARP1, PARP2 |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | TW 37 | 10 | 2 | 5 | -0.323 | 0.70 | 0.226 | BCL2, BCL-XL, MCL1 |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | TW 37 | 10 | 2 | 5 | -0.323 | 0.70 | 0.226 | BCL2, BCL-XL, MCL1 |
| T2_subtype | CORUM::788::Exosome | TW 37 | 10 | 2 | 5 | -0.323 | 0.70 | 0.226 | BCL2, BCL-XL, MCL1 |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | TW 37 | 10 | 2 | 5 | -0.323 | 0.70 | 0.226 | BCL2, BCL-XL, MCL1 |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | PD173074 | 10 | 2 | 5 | -0.234 | 0.80 | 0.187 | FGFR1, FGFR2, FGFR3 |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | PD173074 | 10 | 2 | 5 | -0.234 | 0.80 | 0.187 | FGFR1, FGFR2, FGFR3 |
| T2_subtype | CORUM::788::Exosome | PD173074 | 10 | 2 | 5 | -0.234 | 0.80 | 0.187 | FGFR1, FGFR2, FGFR3 |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | PD173074 | 10 | 2 | 5 | -0.234 | 0.80 | 0.187 | FGFR1, FGFR2, FGFR3 |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | JNK Inhibitor VIII | 10 | 2 | 5 | -0.257 | 0.70 | 0.180 | JNK |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | JNK Inhibitor VIII | 10 | 2 | 5 | -0.257 | 0.70 | 0.180 | JNK |
| T2_subtype | CORUM::788::Exosome | JNK Inhibitor VIII | 10 | 2 | 5 | -0.257 | 0.70 | 0.180 | JNK |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | JNK Inhibitor VIII | 10 | 2 | 5 | -0.257 | 0.70 | 0.180 | JNK |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | AZD1332 | 10 | 2 | 5 | -0.208 | 0.80 | 0.166 | NTRK1, NTRK2, NTRK3 |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | AZD1332 | 10 | 2 | 5 | -0.208 | 0.80 | 0.166 | NTRK1, NTRK2, NTRK3 |
| T2_subtype | CORUM::788::Exosome | AZD1332 | 10 | 2 | 5 | -0.208 | 0.80 | 0.166 | NTRK1, NTRK2, NTRK3 |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | AZD1332 | 10 | 2 | 5 | -0.208 | 0.80 | 0.166 | NTRK1, NTRK2, NTRK3 |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | Cytarabine | 10 | 2 | 5 | -0.208 | 0.80 | 0.166 | Antimetabolite |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | Cytarabine | 10 | 2 | 5 | -0.208 | 0.80 | 0.166 | Antimetabolite |
| T2_subtype | CORUM::788::Exosome | Cytarabine | 10 | 2 | 5 | -0.208 | 0.80 | 0.166 | Antimetabolite |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | Cytarabine | 10 | 2 | 5 | -0.208 | 0.80 | 0.166 | Antimetabolite |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | SN-38 | 10 | 2 | 5 | -0.229 | 0.70 | 0.160 | TOP1 |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | SN-38 | 10 | 2 | 5 | -0.229 | 0.70 | 0.160 | TOP1 |
| T2_subtype | CORUM::788::Exosome | SN-38 | 10 | 2 | 5 | -0.229 | 0.70 | 0.160 | TOP1 |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | SN-38 | 10 | 2 | 5 | -0.229 | 0.70 | 0.160 | TOP1 |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | Wee1 Inhibitor | 10 | 2 | 5 | -0.227 | 0.70 | 0.159 | WEE1, CHEK1 |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | Wee1 Inhibitor | 10 | 2 | 5 | -0.227 | 0.70 | 0.159 | WEE1, CHEK1 |
| T2_subtype | CORUM::788::Exosome | Wee1 Inhibitor | 10 | 2 | 5 | -0.227 | 0.70 | 0.159 | WEE1, CHEK1 |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | Wee1 Inhibitor | 10 | 2 | 5 | -0.227 | 0.70 | 0.159 | WEE1, CHEK1 |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | Daporinad | 10 | 2 | 5 | -0.208 | 0.70 | 0.145 | NAMPT |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | Daporinad | 10 | 2 | 5 | -0.208 | 0.70 | 0.145 | NAMPT |
| T2_subtype | CORUM::788::Exosome | Daporinad | 10 | 2 | 5 | -0.208 | 0.70 | 0.145 | NAMPT |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | Daporinad | 10 | 2 | 5 | -0.208 | 0.70 | 0.145 | NAMPT |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | Motesanib | 10 | 2 | 5 | -0.182 | 0.70 | 0.128 | VEGFR, RET, KIT, PDGFR |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | Motesanib | 10 | 2 | 5 | -0.182 | 0.70 | 0.128 | VEGFR, RET, KIT, PDGFR |
| T2_subtype | CORUM::788::Exosome | Motesanib | 10 | 2 | 5 | -0.182 | 0.70 | 0.128 | VEGFR, RET, KIT, PDGFR |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | Motesanib | 10 | 2 | 5 | -0.182 | 0.70 | 0.128 | VEGFR, RET, KIT, PDGFR |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | KU-55933 | 10 | 2 | 5 | -0.179 | 0.70 | 0.125 | ATM |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | KU-55933 | 10 | 2 | 5 | -0.179 | 0.70 | 0.125 | ATM |
| T2_subtype | CORUM::788::Exosome | KU-55933 | 10 | 2 | 5 | -0.179 | 0.70 | 0.125 | ATM |
| T2_subtype | REACTOME::REACTOME_MRNA_DECAY_BY_3_TO_5_EXORIBONUCLEASE | KU-55933 | 10 | 2 | 5 | -0.179 | 0.70 | 0.125 | ATM |
| T2_subtype | COESSENTIALITY::50::biological process:response to wounding | AZD2014 | 10 | 2 | 5 | -0.208 | 0.60 | 0.125 | mTORC1, mTORC2 |
| T2_subtype | COESSENTIALITY::779::biological process:RNA phosphodiester bond hydrolysis, exonucleolytic | AZD2014 | 10 | 2 | 5 | -0.208 | 0.60 | 0.125 | mTORC1, mTORC2 |
| T2_subtype | CORUM::788::Exosome | AZD2014 | 10 | 2 | 5 | -0.208 | 0.60 | 0.125 | mTORC1, mTORC2 |

## Guardrails

- 尚未加入 FDA/ChEMBL indication、CRC novelty、安全窗和暴露浓度；这些是下一轮 drug annotation，不在本轮臆测。
- GDSC cell-line mapping uses COSMICID→DepMap ModelID; unmapped lines are excluded and recorded by the pairwise n. 
- T2 candidates are subtype hypotheses. T1 candidates are the only ones eligible for a broad OXA-R follow-up screen; no candidate is promoted directly to wet experiment.
