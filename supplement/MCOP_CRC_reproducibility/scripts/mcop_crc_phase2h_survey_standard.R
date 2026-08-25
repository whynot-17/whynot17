#!/usr/bin/env Rscript

# Independent Phase 2H validation using survey::svyglm.
# The Python driver writes a complete-case CSV and supplies the output path.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mcop_crc_phase2h_survey_standard.R <input_csv> <output_csv>")
}

input_csv <- args[[1]]
output_csv <- args[[2]]
library_path <- "D:/CodexData/R-library"
if (dir.exists(library_path)) {
  .libPaths(c(library_path, .libPaths()))
}
if (!requireNamespace("survey", quietly = TRUE)) {
  stop("The survey package is not installed in the configured R library.")
}

library(survey)

dat <- read.csv(
  input_csv,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN", "<NA>", "nan"),
  check.names = FALSE
)

required <- c(
  "outcome", "mcop_log2", "age", "bmi", "pir", "creatinine_log2",
  "sex", "race", "smoking", "pooled_weight", "psu", "strata"
)
missing <- setdiff(required, names(dat))
if (length(missing) > 0) {
  stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
}

dat <- dat[complete.cases(dat[, required]) & dat$pooled_weight > 0, , drop = FALSE]
dat$outcome <- as.numeric(dat$outcome)
dat$pooled_weight <- as.numeric(dat$pooled_weight)
dat$psu <- as.numeric(dat$psu)
dat$strata <- as.numeric(dat$strata)
dat$sex <- factor(dat$sex, levels = c("Female", "Male"))
dat$race <- factor(
  dat$race,
  levels = c("Non-Hispanic White", "Mexican American", "Other Hispanic", "Non-Hispanic Black", "Other/Multi")
)
dat$smoking <- factor(dat$smoking, levels = c("Never", "Former", "Current"))

unit_counts <- aggregate(psu ~ strata, data = dat, FUN = function(x) length(unique(x)))
singleton_n <- sum(unit_counts$psu < 2)
lonely_option <- if (singleton_n > 0) "certainty" else "fail"
options(survey.lonely.psu = lonely_option)

design <- svydesign(
  ids = ~psu,
  strata = ~strata,
  weights = ~pooled_weight,
  nest = TRUE,
  data = dat
)

fit <- svyglm(
  outcome ~ mcop_log2 + age + bmi + pir + creatinine_log2 + sex + race + smoking,
  design = design,
  family = quasibinomial(),
  na.action = na.omit
)

coef_table <- summary(fit)$coefficients
term <- "mcop_log2"
if (!(term %in% rownames(coef_table))) {
  stop("mcop_log2 term was not estimable in svyglm.")
}
beta <- unname(coef_table[term, "Estimate"])
se <- unname(coef_table[term, "Std. Error"])
design_df <- survey::degf(design)
model_residual_df <- df.residual(fit)
critical_model <- qt(0.975, model_residual_df)
critical_design <- qt(0.975, design_df)
statistic <- beta / se
p_standard <- unname(coef_table[term, "Pr(>|t|)"])
p_design_df <- 2 * pt(-abs(statistic), df = design_df)

row <- data.frame(
  method = "R survey::svyglm",
  survey_package_version = as.character(utils::packageVersion("survey")),
  status = "ok",
  N = nrow(dat),
  CRC_N = sum(dat$outcome == 1),
  Control_N = sum(dat$outcome == 0),
  beta = beta,
  SE = se,
  OR = exp(beta),
  CI_low = exp(beta - critical_model * se),
  CI_high = exp(beta + critical_model * se),
  P_standard = p_standard,
  P_design_df = p_design_df,
  design_df = design_df,
  model_residual_df = model_residual_df,
  PSU_N = length(unique(dat$psu)),
  strata_N = length(unique(dat$strata)),
  singleton_strata_N = singleton_n,
  survey_lonely_psu_option = lonely_option,
  formula = "outcome ~ mcop_log2 + age + bmi + pir + creatinine_log2 + sex + race + smoking"
)
write.csv(row, output_csv, row.names = FALSE, na = "")
