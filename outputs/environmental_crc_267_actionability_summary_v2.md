# Complete 267-chemical environmental CRC human-actionability audit v2

Generated: 2026-08-23T16:47:37.975174+00:00

## Scope and firewall

This report covers the complete Phase 1 primary universe: 267 unique core environmental chemicals. Identity, biomarker mapping, detectability, cycle coverage, and survey testability were evaluated before any candidate-specific human CRC association results were considered.

Outcome firewall: `PRIORITIZATION_OUTCOME_BLINDED=True`; candidate-specific human OR/P/CI/LOCO fields used for eligibility: `False`.

## Real attrition flow

| Stage | N |
|---|---:|
| total_core_chemicals | 267 |
| E_entity_valid | 259 |
| E_and_X_interpretable_exposure | 135 |
| E_X_and_B_biomarker_available | 134 |
| E_X_B_and_D_detectable | 127 |
| E_X_B_D_and_C_coverage | 124 |
| E_X_B_D_C_and_T_testable | 87 |
| moderate_eligibility | 87 |
| strict_eligibility | 27 |

The permissive rule was frozen as `E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1`; moderate and strict tiers were retained as sensitivity tiers. The 87 eligible chemical rows collapse to 15 unique biomarker-axis tests because multiple parent chemicals share the same validated NHANES proxy; the human screen is run once per unique axis and reports the member chemicals.

## Tier counts

| Tier | Definition | N chemical rows |
|---|---|---:|
| A strict | D=2, C=2, T=2 | 27 |
| A moderate | D>=1, C=2, T>=1 | 87 |
| B human-testable | permissive eligibility but not moderate | 0 |
| C molecular-only | M2 but not human-testable | 10 |

Manual-review queue: 9 genuine identity/registry exceptions; candidates that were fully searched but lacked a candidate-specific NHANES analyte remain resolved in the 267-row matrix rather than being hidden as pending.

## MCOP and MiNP/DINP status

- **MCOP** (`C573544`, `URXCOP`): retained as a fully eligible DINP-related urinary exposure axis; D=2, C=2, T=2, permissive=True, strict=True. This status is determined from biomarker/actionability data only.
- **MiNP** (`C471400`, `URXMNP`): not discarded or merged into MCOP; it remains a distinct DINP molecular nominee, but its direct urinary detectability is D=0 (40.7% above LOD), so it does not enter the primary human screen under the frozen D gate.
- The DINP-related axis is therefore represented by distinct records: MiNP for molecular nomination and MCOP for the human biomarker axis. No candidate-specific CRC association was used to make this decision.

## Interpretation

After complete outcome-blinded annotation of all 267 original chemicals, 87 chemical candidates satisfy the permissive human-testability rule, representing 15 unique exposure-axis/biomarker tests. The subsequent systematic screen must include all 15 axes and apply BH-FDR across those tested axes.

## Files

- `environmental_crc_267_actionability_matrix_v2.csv` — one row per original chemical.
- `environmental_crc_267_human_testable_candidates.csv` — unique axis keys and member chemicals entering the human screen.
- `environmental_crc_267_biomarker_mapping.csv` — candidate-to-biomarker evidence trail.
- `environmental_crc_267_detectability_by_cycle.csv` — cycle-level measured/above-LOD counts.
- `environmental_crc_267_testability_audit.csv` — pre-outcome survey infrastructure and case/control availability.
- `environmental_crc_267_manual_review_queue_v2.csv` — only genuine identity/registry exceptions.
