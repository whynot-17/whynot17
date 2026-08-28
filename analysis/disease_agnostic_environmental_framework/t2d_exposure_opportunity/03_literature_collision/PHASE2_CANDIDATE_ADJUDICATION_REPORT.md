# Phase 2 — Literature collision, novelty, and opportunity audit

## Status

- **Status:** `complete_first_pass_candidate_adjudication_provisional_top5`
- **Audit date (UTC):** `2026-08-28`
- **PubMed search date:** `2026-08-28`
- **Scope:** all 15 chemical IDs in the Phase 1 advance pool, represented by 14 search groups because the two DINP parent IDs share one parent-specific vocabulary.
- **Methods:** reproducible PubMed category searches, retrieval of top relevance-sorted records, and a transparent first-pass title/abstract-level candidate adjudication.

## Interpretation boundary

The category counts are collision-screening signals, not eligible-study counts. The first-pass adjudication checks whether retrieved titles support candidate identity, human/T2D relevance, prospective design, mechanism maturity, target status, and network-toxicology collision. It is not a full-text systematic review and does not establish causality.

## Provisional Top 5 opportunity pool

The Top 5 is an opportunity shortlist, not a claim that these candidates have the strongest causal evidence. It prioritizes a combination of mapping actionability, literature headroom, and an explicit unresolved gap. Ranking may change after full-text eligibility review.

| Rank | Exposure opportunity | Mapping | Search-group diabetes hits | Human hits | Prospective hits | Novelty/opportunity rationale |
|---:|---|---|---:|---:|---:|---|
| 1 | mono(carboxy-isooctyl)phthalate (URXCOP) | B | 4 | 2 | 0 | A_high_novelty_high_evidence_gap; Advance as a high-novelty opportunity, but require parent-specific full-text adjudication before any mechanistic claim. |
| 2 | diisononyl phthalate; dinonylphthalate (URXCOP) | B | 11 | 8 | 3 | A_high_novelty_specificity_uncertain; Advance as a parent-specificity opportunity; do not count family-level phthalate records as direct DINP evidence. |
| 3 | 2-ethyl-5-carboxypentyl phthalate (URXECP) | B | 21 | 16 | 3 | B_emerging_literature_not_saturated; Advance as an emerging opportunity with a promising but not yet complete T2D-specific mechanism. |
| 4 | Tin (URXUSN) | A | 47 | 11 | 4 | B_high_species_resolution_gap; Advance only as a clearly labeled speciation-limited opportunity; do not merge elemental tin and organotin evidence. |
| 5 | mono-isobutyl phthalate (URXMIB) | B | 33 | 28 | 11 | B_moderate_novelty_family_collision; Advance as a tractable secondary opportunity, below the sparse-collision candidates because phthalate-family literature is crowded. |

## Candidate-level disposition

- `candidate_evidence_adjudication.csv` contains one row per Phase 1 advance-pool chemical ID, including all 15 IDs and the shared DINP parent search group.
- `PHASE2_PROVISIONAL_TOP5_OPPORTUNITIES.csv` is the compact shortlist for the next full-text collision pass.
- All non-Top-5 candidates remain in the adjudication table with an explicit reason for monitoring or exclusion from the novelty-led shortlist.

## Main opportunity gaps

1. Sparse candidates such as MCOP require parent-specific full-text review and independent human/prospective evidence before mechanistic interpretation.
2. DINP and several phthalate metabolites require strict separation of parent-specific, metabolite-specific, replacement, mixture, and family-level evidence.
3. Tin requires chemical-speciation resolution: elemental urinary tin cannot be treated as interchangeable with organotin studies.
4. Candidates with high collision counts (lead, DEHP, silver nanoparticle literature, and established metal/PFAS axes) are useful comparators but are not prioritized as novelty opportunities.

## Next required step

Perform title/abstract adjudication followed by targeted full-text retrieval for the Top 5 and any candidate whose mapping identity or exposure form is ambiguous. Freeze eligibility criteria before reviewing outcome-specific details; do not run docking or expand target discovery in this phase.
