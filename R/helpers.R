# =============================================================================
# SHARED FUNCTIONS
#
# The numbered scripts call these functions for tasks used more than once:
# cleaning ICD-10 text, decoding NCHS demographic fields, assigning mutually
# exclusive drug groups, fitting the common logistic model, extracting model
# estimates, and writing outputs. Keeping these rules here ensures that the
# primary and sensitivity analyses apply the same definitions.
# =============================================================================

# Standardize an ICD-10 value before matching. For example, "T40.5" and
# "T405" both become "T405". Missing codes become empty strings and therefore
# do not match any disease or drug pattern.
normalize_icd <- function(x) {
  x <- toupper(as.character(x)); x <- gsub("[^A-Z0-9]", "", x)
  x[is.na(x)] <- ""; x
}
# Return a named column when it exists. NCHS demographic fields differ across
# coding eras, so a field that does not apply in a particular year is replaced
# by an all-missing vector of the correct length.
get_col <- function(d, name) if (name %in% names(d)) d[[name]] else rep(NA, nrow(d))

# Convert the NCHS detail-age field to years. Adult ages may already be stored
# as years; infant ages can be stored with a unit code for months, weeks, days,
# hours, or minutes. Values outside 0-125 years are treated as missing.
decode_age <- function(age) {
  age <- as.numeric(age); unit <- age %/% 1000; value <- age %% 1000
  out <- ifelse(age >= 0 & age <= 125, age,
    ifelse(unit == 1, value, ifelse(unit == 2, value / 12,
    ifelse(unit == 3, value / 52.1775, ifelse(unit == 4, value / 365.25,
    ifelse(unit == 5, value / 8766, ifelse(unit == 6, value / 525960, NA_real_)))))))
  out[out < 0 | out > 125] <- NA_real_; out
}
# Harmonize numeric and character sex codes across annual files.
decode_sex <- function(x) {
  x <- toupper(as.character(x)); out <- rep(NA_character_, length(x))
  out[x %in% c("M", "1")] <- "Male"; out[x %in% c("F", "2")] <- "Female"; out
}
# Collapse race to four manuscript categories while respecting changes in the
# available NCHS race recode: detailed race through 2020, racerec40 in 2021,
# and racerec5 beginning in 2022. AIAN, multiracial, and unknown values are
# included in Other, as specified in the manuscript.
decode_race <- function(year, race, racerec5, racerec40) {
  out <- rep("Other", length(year)); early <- year <= 2020
  out[early & race == 1] <- "White"; out[early & race == 2] <- "Black"
  out[early & race %in% c(4,5,6,7,18,28,38,48,58,68,78)] <- "API"
  y21 <- year == 2021
  out[y21 & racerec40 == 1] <- "White"; out[y21 & racerec40 == 2] <- "Black"
  out[y21 & racerec40 %in% 4:14] <- "API"
  late <- year >= 2022
  out[late & racerec5 == 1] <- "White"; out[late & racerec5 == 2] <- "Black"
  out[late & racerec5 %in% c(4,5)] <- "API"; out
}
# Harmonize Hispanic-origin fields across NCHS coding eras. Hispanic origin is
# modeled separately from race. Values that cannot be assigned are retained as
# an explicit Unknown category rather than dropped during data construction.
decode_hispanic <- function(year, hispanic, hisp_origin) {
  out <- rep("Unknown", length(year)); early <- year <= 2002
  out[early & hispanic == 0] <- "Non-Hispanic"; out[early & hispanic %in% 1:5] <- "Hispanic"
  y03 <- year == 2003
  out[y03 & hispanic == 100] <- "Non-Hispanic"; out[y03 & hispanic >= 200 & hispanic < 996] <- "Hispanic"
  late <- year >= 2004
  out[late & hisp_origin == 100] <- "Non-Hispanic"; out[late & hisp_origin >= 200 & hisp_origin < 996] <- "Hispanic"; out
}
# Convert two yes/no indicators into the three mutually exclusive exposure
# groups used throughout the paper. A record with neither drug class receives
# NA and is excluded from analyses using that definition.
drug_group <- function(opioid, stimulant) {
  ifelse(opioid & !stimulant, "opiate_no_stimulant",
    ifelse(opioid & stimulant, "opiate_stimulant",
      ifelse(!opioid & stimulant, "stimulant_no_opiate", NA_character_)))
}
# Prepare a dataset for the common regression specification. Stimulant-only is
# the reference exposure group; Male, White, and Non-Hispanic are the reference
# covariate categories. Dividing age by 10 makes its coefficient a 10-year
# comparison and does not change which records enter the model.
prepare_model_data <- function(d, group) {
  use <- !is.na(group); x <- d[use, ]
  x$grp <- factor(group[use],
    levels = c("stimulant_no_opiate", "opiate_stimulant", "opiate_no_stimulant"),
    labels = c("Stimulant-only", "Opioid-stimulant", "Opioid-only"))
  x$age_10 <- x$age_years / 10
  x$sex <- factor(x$sex, levels = c("Male", "Female"))
  x$race4 <- factor(x$race4, levels = c("White", "Black", "API", "Other"))
  x$hisp3 <- factor(x$hisp3, levels = c("Non-Hispanic", "Hispanic", "Unknown"))
  x
}
# Names printed in the model output tables. All four outcomes are binary.
outcomes <- c(cardiovascular = "Cardiovascular disease",
  cerebrovascular = "Cerebrovascular disease",
  other_medical = "Other medical condition",
  no_additional = "No additional condition")
# Fit the adjusted logistic regression used in the primary analysis and each
# sensitivity analysis. glm uses complete cases for the outcome, exposure, and
# adjustment variables. Coefficients are on the log-odds scale.
fit_adjusted <- function(outcome, data) {
  glm(as.formula(paste(outcome, "~ grp + age_10 + sex + race4 + hisp3")),
      data = data, family = binomial())
}
# Extract selected exposure coefficients from a fitted model. The function
# reports log-odds coefficients, standard errors, Wald confidence intervals,
# odds ratios, P values, model N, and the number of outcome events.
extract_terms <- function(fit, terms, labels, outcome) {
  sm <- summary(fit)$coefficients; ci <- confint.default(fit)
  data.frame(outcome = unname(outcomes[outcome]), contrast = labels,
    beta = unname(sm[terms,"Estimate"]), beta_se = unname(sm[terms,"Std. Error"]),
    beta_ci_low = unname(ci[terms,1]), beta_ci_high = unname(ci[terms,2]),
    odds_ratio = exp(unname(sm[terms,"Estimate"])),
    or_ci_low = exp(unname(ci[terms,1])), or_ci_high = exp(unname(ci[terms,2])),
    p_value = unname(sm[terms,"Pr(>|z|)"]), model_n = nobs(fit),
    outcome_events = sum(model.response(model.frame(fit)) == 1), stringsAsFactors = FALSE)
}
# Save a compact check for each fitted model. Convergence is the required
# automated gate; N, event count, iteration count, coefficient magnitude, and
# AIC help a reviewer identify unexpected differences across replications.
model_diagnostic <- function(fit, outcome, analysis) {
  data.frame(analysis = analysis, outcome = unname(outcomes[outcome]),
    n = nobs(fit), events = sum(model.response(model.frame(fit)) == 1),
    converged = isTRUE(fit$converged), iterations = fit$iter,
    max_abs_coefficient = max(abs(coef(fit))), aic = AIC(fit), stringsAsFactors = FALSE)
}
# Formatting and file-writing functions used by all table scripts.
fmt_p <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", p)))
write_table <- function(x, name) write.csv(x, file.path(table_dir, name), row.names = FALSE, na = "")
write_diag <- function(x, name) write.csv(x, file.path(diagnostic_dir, name), row.names = FALSE, na = "")
# Start a timestamped text log for a numbered script. Logs contain aggregate
# counts and status messages only; they never print record-level data.
new_logger <- function(name) {
  path <- file.path(log_dir, name); cat("", file = path)
  function(...) { z <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...)); cat(z,"\n"); cat(z,"\n",file=path,append=TRUE) }
}
