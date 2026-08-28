# Phase 1 — T2D exposure-opportunity candidate and mapping gate

## Material Passport

- **Material ID:** `T2D-EXPOSURE-OPPORTUNITY-PHASE1`
- **Status:** `complete_mapping_gate_only`
- **Script:** `phase1_mapping_gate_v1.0`
- **Generated (UTC):** `2026-08-28T03:51:50.988477+00:00`
- **Scope:** all upstream mapped chemicals associated with the 14 frozen T2D FDR-positive tests.
- **Not performed:** literature collision audit, mechanism analysis, target nomination, docking, experimental feasibility, and opportunity scoring.

## Frozen input and provenance

- T2D input: `step05_t2d_screen/t2d_primary_29_tests.csv`; positive tests were derived from `FDR_supported=True` and verified to use the frozen 29-test denominator.
- Chemical mapping input: `step02_biomarker_mapping/chemical_biomarker_mapping.csv`; only rows with `mapped=True` and one of the 14 positive biomarkers were included in the candidate master.
- Robustness annotation: `step06_t2d_robustness/t2d_robustness_results.csv`; used only as a downstream annotation, not as a mapping gate.
- No GeneCards, CTD chemical–gene interactions, disease-specific pathway data, or literature counts were used in this phase.

## Counts

- Positive T2D tests: **14**
- Upstream mapped chemical–biomarker rows: **136**
- Unique upstream chemical IDs: **134**
- Mapping rows by grade: **A=8, B=9, C=119**
- Mapping rows by gate status: **{"conditional": 103, "exclude": 16, "pass": 14, "pass_conditional": 3}**
- Unique-chemical dispositions: **{"advance_to_literature_audit": 12, "advance_with_parent_specificity_review": 3, "exclude_pending_mapping_review": 16, "proxy_only_not_primary": 103}**

## Per-test upstream mapping coverage

| Positive test | Upstream mapped rows |
|---|---:|
| `LBXPFHS` | 1 |
| `URXCOP` | 3 |
| `URXECP` | 2 |
| `URXMHH` | 1 |
| `URXMIB` | 1 |
| `URXMOH` | 2 |
| `URXP02` | 12 |
| `URXUBA` | 8 |
| `URXUMO` | 18 |
| `URXUPB` | 29 |
| `URXUSN` | 18 |
| `URXUSR` | 16 |
| `URXUTU` | 12 |
| `URXUUR` | 13 |

## Gate interpretation

- **Grade A:** direct analyte identity supported by the neutral mapping record (parent elemental analyte or direct serum PFAS analyte).
- **Grade B:** specific urinary metabolite or parent–validated-metabolite relationship; direct urinary metabolites retain a parent-inference caution.
- **Grade C:** family/proxy mapping or an elemental/species name not resolved by the assay. These rows remain in the full master for audit but are not promoted as compound-specific primary candidates.
- A chemical is marked `advance_to_literature_audit` only when it has a direct Grade A mapping with `pass`, or a specific Grade B mapping with the corresponding review flag. C-only mappings are `proxy_only_not_primary`; internally inconsistent elemental rows are excluded pending manual mapping review.

## File outputs

- `01_candidate_master/all_upstream_chemicals.csv`: one row per mapped chemical–positive biomarker relationship.
- `01_candidate_master/unique_candidate_chemicals.csv`: one row per unique chemical ID.
- `02_mapping_audit/mapping_specificity_audit.csv`: row-level grade and gate rationale.
- `02_mapping_audit/mapping_exclusion_log.csv`: proxy/conditional/mismatch rows not eligible for an unqualified compound-specific shortlist.

This is a mapping gate, not a novelty or mechanistic conclusion. A Grade A/B mapping means the exposure-to-analyte link is more actionable; it does not establish a T2D causal association.
