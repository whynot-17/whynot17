# Frozen Outcome Set v1.0: sex-specific feasibility audit

This is an outcome-only inventory. The program reads demographics and the small, explicit source-file allowlist in `outcome_inventory_qc.json`; it does not read exposure measurements, any association result, or candidate-chemical result.

## Prespecified screen

PASS requires male and female pooled cases each >=100 and at least 4 independent definable cycles. AMBER applies when either sex has 50-99 cases (and >=4 cycles). FAIL applies when either sex has <50 cases or fewer than 4 cycles are definable.

## Local-source limitations

Objective CKD, kidney stones, osteoporosis, and clinically relevant depressive symptoms are retained in the frozen 18-candidate registry but are marked unavailable where their required raw modules are absent locally. Albuminuria alone is not substituted for eGFR-defined CKD.

## Outputs

- `outcome_cycle_sex_audit.csv` — one row per candidate outcome and cycle, including sex-specific cases and controls.
- `outcome_pooled_sex_audit.csv` — pooled sex-specific counts and prespecified status.
- `frozen_outcome_set_v1.csv` — all 18 pre-frozen candidates with follow-up flag.
- `outcome_inventory_qc.json` — provenance, input allowlist, prohibition statement, and hashes.
