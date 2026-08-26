# Step 9D — Question-Specific Data Readiness Diagnostic

## Purpose

Step 9D formally closes the CRC negative branch by integrating Steps 9A–9C.
It does not rerun models, revise the 29-test family, change BH-FDR, or add
outcome-informed candidate selection. The purpose is to distinguish why the CRC
application did not cross the discovery threshold from what the current design
cannot interpret even if precision improves.

## Final diagnostic statement

The CRC negative screen is best described as a **mixed failure architecture**:

1. **Statistical/readiness limitations:** CRC has 420 pooled assay-specific QC cases versus 7,772 for T2D, and a median of 97 analytic cases per test versus 1,940. Assay-specific complete-case retention is also lower for CRC.
2. **Observed signal heterogeneity:** some tests are near-null, whereas others have directional estimates whose current uncertainty is too large for the prespecified OR=1.20 detectability reference.
3. **Structural design limitations:** prevalent cross-sectional CRC ascertainment, no prediagnostic biospecimen, and no diagnosis date/stage/site/treatment/recurrence/follow-up prevent a strong exposure-before-disease interpretation.

This means **CRC failure does not equal biological absence**, but neither does it
mean that additional cases would guarantee a positive result. The framework
diagnoses answerability rather than manufacturing an explanation for every null.

## Statistical versus structural limitations

| Layer | What the current data show | What it affects | What it does not establish |
|---|---|---|---|
| Statistical/readiness | Low CRC event density and lower assay-specific complete-case retention | Precision, event support, detectability | Does not prove a latent association |
| Observed signal | Near-null and directional/imprecise tests coexist | Which tests merit any future expansion | Does not convert power-limited labels into true effects |
| Structural design | Prevalent cross-sectional outcome and limited phenotype/timing fields | Temporality, reverse-causation protection, phenotype resolution | Cannot be repaired by N alone |

## Counterfactual guidance

Under the frozen SE-scaling approximation, current MDE<=OR 1.20 was met by only
1/29 tests, versus 9/29 at 2× cases, 20/29 at 4×, and 24/29 at 8×. Retention
improvement to 70% and 85% reached 2/29 and 7/29, respectively. These are
detectability scenarios, not new inferential results: no P values, FDR values,
or associations were simulated.

The practical conclusion is that additional CRC events are the strongest
precision intervention in this approximation, while better retention is a
useful complementary intervention. Prediagnostic samples, diagnosis date,
stage/site, longitudinal follow-up, and repeated exposure measurements address
different structural limitations and should not be treated as substitutes for
event support.

## Symmetric framework output

The positive T2D branch ends in disease-specific biological prioritization after
the exposure panel is frozen. The CRC branch ends in a diagnosed failure
architecture and actionable data-improvement guidance after the same firewall.
This symmetry is the methodological result; it is not a claim that positive and
negative disease applications have identical data conditions.

## Scope and provenance

- Source inputs are the locked Step 9A readiness, Step 9B attribution, and Step 9C counterfactual outputs.
- No new model fit, P value, FDR, candidate ranking, or GeneCards/CTD search was performed.
- The full synthesis matrix is in `step9d_negative_branch_synthesis.csv`.
- The branch symmetry table is in `step9d_positive_negative_branch_symmetry.csv`.
