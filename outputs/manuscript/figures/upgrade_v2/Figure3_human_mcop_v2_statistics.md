# Figure 3 v2 statistics and traceability

**Panel A.** Primary complex-survey logistic model. Exposure is log2 urinary MCOP; estimate is an odds ratio per doubling. Source: `figure2_primary_python_vs_r.csv`; independent R survey columns `r_OR`, `r_CI_low`, `r_CI_high`, `r_P_design_df`.

**Panel B.** Seven leave-one-cycle-out reanalyses from `figure2_loco.csv`; each row is a pooled complex-survey estimate after excluding one NHANES cycle.

**Panel C.** Seven cycle-specific complex-survey estimates from `figure2_per_cycle.csv`; parenthetical counts are CRC cases per cycle. The global MCOP-by-cycle interaction P=0.0060 is an audited model-level statistic.

**Panel D.** Cross-implementation QC: R `survey::svyglm` and Python Taylor-sandwich estimates agree to numerical precision and have the same CI conclusion.
