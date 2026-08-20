# =============================================================================
# SCRIPT 07: SENSITIVITY ANALYSIS ADJUSTING FOR ALCOHOL INVOLVEMENT
#
# Purpose
# Compare the primary adjusted logistic regression models with otherwise
# identical models that also adjust for alcohol involvement on the death
# certificate. This analysis was specified post hoc on 2026-08-20.
#
# Alcohol definition
# alcohol_any equals 1 when F10, K70, X45, or T51 appears as the underlying
# cause or in an available record-axis multiple-cause field. Script 01 creates
# this indicator. It describes co-recorded alcohol involvement and is not
# assumed to be a pre-exposure confounder.
#
# Sample comparison
# Each primary model and its alcohol-adjusted counterpart use the same
# complete-case records. This ensures that coefficient changes result from
# adding alcohol_any rather than from changing the estimation sample.
# =============================================================================

source(file.path(publication_root, "R", "helpers.R"))
log_msg <- new_logger("07_sensitivity_alcohol_adjustment.log")
d <- readRDS(analytic_path)
x <- prepare_model_data(d, d$comparison_group)
terms <- c("grpOpioid-stimulant", "grpOpioid-only")
term_labels <- c("Opioid-stimulant vs Stimulant-only",
                 "Opioid-only vs Stimulant-only")
results <- list()
diagnostics <- list()

for (outcome in names(outcomes)) {
  required <- c(outcome, "grp", "age_10", "sex", "race4", "hisp3",
                "alcohol_any")
  same_sample <- x[complete.cases(x[, required]), ]
  primary_fit <- glm(
    as.formula(paste(outcome, "~ grp + age_10 + sex + race4 + hisp3")),
    data = same_sample, family = binomial()
  )
  alcohol_fit <- glm(
    as.formula(paste(outcome,
      "~ grp + age_10 + sex + race4 + hisp3 + alcohol_any")),
    data = same_sample, family = binomial()
  )
  primary_coef <- summary(primary_fit)$coefficients
  alcohol_coef <- summary(alcohol_fit)$coefficients
  primary_ci <- confint.default(primary_fit)
  alcohol_ci <- confint.default(alcohol_fit)

  for (j in seq_along(terms)) {
    term <- terms[j]
    beta_primary <- unname(primary_coef[term, "Estimate"])
    beta_alcohol <- unname(alcohol_coef[term, "Estimate"])
    results[[length(results) + 1L]] <- data.frame(
      outcome = unname(outcomes[outcome]), comparison = term_labels[j],
      model_n = nrow(same_sample), primary_beta = beta_primary,
      primary_se = unname(primary_coef[term, "Std. Error"]),
      primary_or = exp(beta_primary),
      primary_or_ci_low = exp(unname(primary_ci[term, 1])),
      primary_or_ci_high = exp(unname(primary_ci[term, 2])),
      alcohol_adjusted_beta = beta_alcohol,
      alcohol_adjusted_se = unname(alcohol_coef[term, "Std. Error"]),
      alcohol_adjusted_or = exp(beta_alcohol),
      alcohol_adjusted_or_ci_low = exp(unname(alcohol_ci[term, 1])),
      alcohol_adjusted_or_ci_high = exp(unname(alcohol_ci[term, 2])),
      beta_difference = beta_alcohol - beta_primary,
      percent_change_abs_beta =
        100 * (abs(beta_alcohol) - abs(beta_primary)) / abs(beta_primary),
      conclusion_changed =
        sign(beta_primary) != sign(beta_alcohol) ||
        ((primary_ci[term, 1] <= 0 && primary_ci[term, 2] >= 0) !=
         (alcohol_ci[term, 1] <= 0 && alcohol_ci[term, 2] >= 0)),
      stringsAsFactors = FALSE
    )
  }

  alcohol_term <- alcohol_coef["alcohol_any", ]
  diagnostics[[length(diagnostics) + 1L]] <- data.frame(
    analysis = "Alcohol adjustment", outcome = unname(outcomes[outcome]),
    n = nrow(same_sample),
    events = sum(same_sample[[outcome]] == 1),
    converged = isTRUE(primary_fit$converged) && isTRUE(alcohol_fit$converged),
    alcohol_prevalence_percent = 100 * mean(same_sample$alcohol_any == 1),
    alcohol_beta = unname(alcohol_term["Estimate"]),
    alcohol_or = exp(unname(alcohol_term["Estimate"])),
    alcohol_p_value = unname(alcohol_term["Pr(>|z|)"]),
    primary_aic = AIC(primary_fit), alcohol_adjusted_aic = AIC(alcohol_fit),
    likelihood_ratio_p = anova(primary_fit, alcohol_fit, test = "LRT")$`Pr(>Chi)`[2],
    stringsAsFactors = FALSE
  )
}

write_table(do.call(rbind, results),
            "sensitivity4_alcohol_adjustment_models.csv")
write_diag(do.call(rbind, diagnostics),
           "sensitivity4_alcohol_adjustment_diagnostics.csv")
log_msg("Alcohol-adjustment sensitivity complete; N=",
        format(unique(do.call(rbind, diagnostics)$n), big.mark = ","))
