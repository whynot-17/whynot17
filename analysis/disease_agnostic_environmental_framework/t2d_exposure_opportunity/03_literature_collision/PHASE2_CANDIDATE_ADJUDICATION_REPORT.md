# Phase 2 — Literature collision, novelty, and opportunity audit

## Status

- **Status:** `complete_first_pass_internal_collision_reranking`
- **Audit date (UTC):** `2026-08-28`
- **PubMed search date:** `2026-08-28`
- **Scope:** all 15 chemical IDs in the Phase 1 advance pool, represented by 14 search groups because the two DINP parent IDs share one parent-specific vocabulary.
- **Methods:** reproducible PubMed category searches, retrieval of top relevance-sorted records, and a transparent first-pass title/abstract-level candidate adjudication.

## Interpretation boundary

The category counts are collision-screening signals, not eligible-study counts. The first-pass adjudication checks whether retrieved titles support candidate identity, human/T2D relevance, prospective design, mechanism maturity, target status, and network-toxicology collision. It is not a full-text systematic review and does not establish causality.

## Internal project-collision audit

The existing DINP/MCOP–CRC project is treated as an internal-collision reference. Each candidate is annotated for overlap in exposure axis, urinary biomarker, and mechanism architecture. Exact reuse receives the largest penalty; related phthalate/PPAR framing receives an intermediate penalty; distinct chemical classes receive zero or minimal penalty. This is a transparent paper-design constraint, not a biological or inferential statistic.

The external-opportunity baseline and internal-collision penalty are embedded in `build_phase2_candidate_adjudication.py`; they are ordinal prioritization aids and must not be interpreted as effect sizes, probabilities, or evidence scores.

## Revised Top 5 opportunity pool

The revised Top 5 is an opportunity shortlist after applying the internal-collision constraint. It prioritizes a combination of mapping actionability, external literature headroom, and low self-overlap. It is not a claim that these candidates have the strongest causal evidence, and ranking may change after full-text eligibility review.

| Rank | Exposure opportunity | Mapping | External score | Self-overlap penalty | Revised score | Position |
|---:|---|---|---:|---:|---:|---|
| 1 | Tin (URXUSN) | A | 7.5 | 0.0 | 7.5 | eligible_distinct_axis |
| 2 | perfluorohexanesulfonic acid (LBXPFHS) | A | 6.5 | 0.5 | 6.0 | eligible_distinct_axis |
| 3 | 2-ethyl-5-carboxypentyl phthalate (URXECP) | B | 8.0 | 2.5 | 5.5 | eligible_with_explicit_differentiation |
| 4 | Uranium (URXUUR) | A | 5.5 | 0.0 | 5.5 | eligible_distinct_axis |
| 5 | mono-isobutyl phthalate (URXMIB) | B | 7.0 | 2.0 | 5.0 | eligible_with_explicit_differentiation |

## Candidate-level disposition

- `candidate_evidence_adjudication.csv` contains one row per Phase 1 advance-pool chemical ID, including all 15 IDs and the shared DINP parent search group.
- `PHASE2_PROVISIONAL_TOP5_OPPORTUNITIES.csv` is the compact self-overlap-adjusted shortlist for the next full-text collision pass.
- `PHASE2_INTERNAL_COLLISION_AUDIT.csv` exposes the overlap dimensions, penalty, revised score, and new-paper position for all 15 candidate IDs.
- `candidate_evidence_adjudication.csv` contains the same internal overlap fields for every candidate ID; MCOP/DINP are retained as explicit high-overlap holdouts rather than silently deleted.
- All non-Top-5 candidates remain in the adjudication table with an explicit reason for monitoring or exclusion from the novelty-led shortlist.

## Main opportunity gaps

1. MCOP and DINP remain externally sparse but are intentionally downranked because they overlap with the existing exposure, biomarker, and mechanistic architecture.
2. Several phthalate metabolites remain related-class opportunities only if the new paper explicitly differentiates exposure, biomarker, and mechanism from the existing project.
3. Tin requires chemical-speciation resolution: elemental urinary tin cannot be treated as interchangeable with organotin studies.
4. Candidates with high external collision counts (lead, DEHP, silver nanoparticle literature, and established metal/PFAS axes) are useful comparators but are not prioritized as novelty opportunities.

## Next required step

Perform title/abstract adjudication followed by targeted full-text retrieval for the revised Top 5 and any candidate whose mapping identity or exposure form is ambiguous. Freeze eligibility criteria before reviewing outcome-specific details; do not run docking or expand target discovery in this phase.
