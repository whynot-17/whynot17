# Focused macrophage–tumor communication analysis

- Dataset: local GSE144735 CRC matrix (27,414 cells, 6 donors); no data were redownloaded.
- Scope: four existing macrophage subtypes as senders → source-labeled tumor epithelial states as receivers. Tumor epithelial labels (CMS1–4) were reused; no new malignant CNV inference was performed.
- Method: pooled descriptive ligand–receptor scores plus donor-level support. R/CellChat was unavailable, so the reproducible fallback is the versioned local human prior in `lr_database_human.csv`; this is not presented as an actual CellChat run.

## Required questions

**Q1. SPP1-like TAM pathways.** Top outgoing pathway: **SPP1**. The top pooled LR record is **SPP1–CD44**; donor support is summarized in `donor_level_lr_support.csv` (consistent pairs: 11).

**Q2. C1QC-like TAM pathways.** Top outgoing pathway: **SPP1**. PTGER4 is quantified separately as a receptor expression context in the prior module; it is not fabricated into a CellChat pair here.

**Q3. Distinct roles.** C1QC-like pooled outgoing score sum is 4.905; the tumor→C1QC incoming context sum is 2.808, with a higher mean per incoming pair. This supports keeping an outgoing SPP1-like axis separate from a PTGER4-enriched C1QC context rather than forcing one role.

**Q4. Strongest TAM→tumor sender.** By pooled focused score, the leading sender is **SPP1_like_TAM**. This is descriptive, not a causal DINP result.

**Q5. Donor repetition.** 16 of the top-20 pooled records meet the predefined support rule in at least four donors. Records with low coverage or donor-driven behavior are explicitly flagged.

**Q6. DINP–CRC overlap.** Overlap annotation is post hoc in `lr_dinp_crc_overlap.csv`; the clearest 7-driver connection is CXCR4 in CXCL12–CXCR4. PTGER4 is not robustly represented in the local LR prior. No overlap was used to filter the network.

**Q7. Spatial readiness.** **Larger-cohort replication first.** Several top subtype–condition signals have sparse paired donor coverage, so these results support prioritizing validation targets rather than immediate causal or spatial claims.

## Interpretation boundary

The results support statements such as “DINP–CRC-associated molecular program localized to…” and “SPP1-like TAM displayed prominent communication scores.” They do not show that DINP causes communication, activates SPP1, or that PTGER4 mediates the network.

Recommended next step: **Larger cohort replication**, then spatial validation if the same sender–receiver axes reproduce.
