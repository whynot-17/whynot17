# Post-primary evidence synthesis plan v1.0

## Purpose and firewall

This plan is frozen after the primary and uniform robustness packages are available, but before an evidence table, pair labels, CTD query, literature search, or mechanistic interpretation is generated. It defines transparent descriptive labels; it creates no new hypothesis test, FDR family, composite score, candidate rank, or selection of new exposures/outcomes.

## One evidence table per pre-specified pair

The table contains all 406 pairs, not only primary findings. For every pair it reports: primary pooled interaction beta, SE, CI, P and fixed-406 BH q; stratified descriptive estimates; DLD and OSD from the corresponding frozen summaries; urine-creatinine sensitivity when applicable; LOCO coefficient range, number of successful refits, and any sign reversal; winsorized, upper-tail deletion, and above-LOD estimates when applicable; and cycle-heterogeneity Wald P with its applicability/status.

Primary q and all robustness quantities remain separate columns. No arithmetic combination, score, rank, or revised FDR is permitted.

## Predefined descriptive labels

`primary_only`: primary interaction q >=0.05. This label is not a statement that the effect is absent; it means only that the pair is not a fixed-family primary finding.

`cycle_sensitive`: primary q <0.05 and either (a) a successful LOCO refit reverses the primary interaction sign, or (b) the applicable cycle-heterogeneity Wald P <0.05. This is a diagnostic description, not a failed replication claim.

`robust_exemplar`: primary q <0.05; no `cycle_sensitive` flag; every successful required direction diagnostic agrees with the primary interaction sign (urine-creatinine adjustment for urine exposures, winsorization, upper-tail deletion, and above-LOD restriction when estimable); and every successful LOCO refit agrees in sign. Non-applicable or non-estimable diagnostics are reported rather than silently counted as agreement.

`primary_fdr_positive_other_robustness_pattern`: primary q <0.05 but neither of the two labels above applies, including a non-cycle sensitivity direction change or insufficient successful diagnostic refits.

These labels do not use the magnitude of DLD/OSD, P-value changes within sensitivities, or stratified-model significance. A small LOCO sign change near zero is retained as a raw diagnostic; the label does not claim a binary success/failure beyond the stated descriptive rule.

## System-level presentation

The seven-system table retains SD and signed SL beside the pair table. DLD answers disease-profile magnitude; OSD answers equal-system landing magnitude. Both are effect-size summaries only and must never be described as statistical evidence or used to choose labels.

## Next-stage gate

Only after the complete 406-row evidence table and its descriptive labels are frozen may a separately locked CTD/mechanism protocol be considered. That protocol may not revise the exposure, outcome, interaction, robustness, or evidence-label definitions.
