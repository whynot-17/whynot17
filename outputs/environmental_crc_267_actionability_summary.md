# Phase 2I — 267-chemical outcome-blinded actionability prioritization

## Direct answer to the original paper question

MCOP is present as **ChemicalID C573544** in the original Phase 1 universe of
**267 core environmental chemicals**. It is not introduced after the
NHANES result and it is not selected by a candidate-specific CRC association.

Under the prespecified actionability rules, MCOP is retained for the systematic
human screen because it is a specific entity with an auditable urinary biomarker,
7-cycle NHANES coverage, 98.4% above-LOD detectability, and sufficient available
sample infrastructure. Its Phase 1 molecular tier is **M0**;
human actionability, rather than molecular rank, is what permits advancement.

This is the direct-discovery formulation:

`267 core chemicals → outcome-blinded actionability annotation → MCOP retained → uniform NHANES screen`

It must not be rewritten as “MiNP molecular signal → MCOP” or as a post-hoc
replacement of MiNP by MCOP. MiNP remains a separate direct analyte and fails the
direct detectability gate in the current audit (27.4% above LOD).

## Current matrix status

- Primary Phase 1 slice: **267 unique chemicals** (`GeneCards_Disorders`, GeneCards k=1000).
- Molecular tiers: **{'M0': 183, 'M1': 67, 'M2': 17}**.
- Human biomarker audit currently available locally: **2 chemicals (MCOP and MiNP only)**.
- Remaining **265 chemicals** are explicitly in the manual-review queue; they are not treated as biomarker-negative.
- Outcome-blinded candidates passing permissive actionability rules: **1**.
- Permissive / moderate / strict counts: **1 / 1 / 1**.

## Frozen actionability rules

`E=1 & X=1 & B=1 & D>=1 & C>=1 & T>=1` (permissive)

Moderate requires `C>=2`; strict requires `D=2 & C=2 & T=2`. Molecular
evidence is reported as M0/M1/M2 and is not a hard gate. Novelty is pending
manual collision review and is never inferred from the NHANES result.

## Outcome firewall

The eligibility expression uses no human CRC OR, CI, P value, LOCO result,
cycle-specific CRC effect, or candidate-specific epidemiologic significance.
Only Phase 1 molecular overlap/FDR and pre-existing sample infrastructure
(analytic N, available case count, control count, weights/cycle coverage) are
used for actionability/testability. The latter are feasibility fields, not
effect estimates.

## Limitation that must remain visible

This is the completed **267-chemical input and firewall audit**, but not yet a
completed 267-chemical biomonitoring annotation: 265 chemicals still require
the same exposure-database/biomarker/detectability review. The current MCOP
direct-discovery claim is valid as an audit of its original presence and its
pre-outcome eligibility; the full-universe claim should not be overstated until
that queue is completed.

Generated: 2026-08-23T15:17:29.093355+00:00
