# =============================================================================
# SCRIPT 03: SENSITIVITY 1, CALENDAR-PERIOD COMPARISON
#
# This post hoc analysis asks whether the two drug-group coefficients changed
# between 2015-2019 and 2020-2024. It does not estimate a causal effect of the
# COVID-19 pandemic; the split is a calendar-period comparison.
# =============================================================================

# ---- 1. Create the two periods using the inclusive primary cohort ------------
source(file.path(publication_root,"R","helpers.R")); log_msg <- new_logger("03_sensitivity_calendar_period.log")
d <- readRDS(analytic_path); x <- prepare_model_data(d,d$comparison_group)
x$period <- factor(ifelse(x$year<=2019,"2015-2019","2020-2024"),levels=c("2015-2019","2020-2024"))
groups <- c("Opioid-stimulant","Opioid-only"); coefficients <- list(); changes <- list(); diagnostics <- list()
# ---- 2. Fit the paper's adjusted model separately within each period ---------
# These fits provide the period-specific coefficients and adjusted odds ratios
# shown in the appendix. The same outcome, covariates, and stimulant-only
# reference group are used in each period.
for(o in names(outcomes)) {
  for(period in levels(x$period)) {
    dd <- x[x$period==period,]; fit <- fit_adjusted(o,dd); sm <- summary(fit)$coefficients; ci <- confint.default(fit)
    for(g in groups) { term <- paste0("grp",g); coefficients[[paste(o,period,g)]] <- data.frame(outcome=unname(outcomes[o]),period=period,comparison=paste(g,"vs Stimulant-only"),n=nobs(fit),events=sum(model.response(model.frame(fit))),beta=sm[term,"Estimate"],beta_se=sm[term,"Std. Error"],adjusted_or=exp(sm[term,"Estimate"]),ci_low=exp(ci[term,1]),ci_high=exp(ci[term,2]),p_value=sm[term,"Pr(>|z|)"]) }
    diagnostics[[paste(o,period)]] <- model_diagnostic(fit,o,paste("Calendar period",period))
  }
  # ---- 3. Formally test whether the drug-group coefficients changed ----------
  # Both pooled models allow the age, sex, race, and Hispanic-origin
  # coefficients to differ by period. The reduced model constrains only the two
  # drug-group coefficients to be constant. The full model adds drug-group by
  # period interactions. This specification makes each interaction coefficient
  # equal the difference between the separately fitted period coefficients.
  reduced <- glm(as.formula(paste(o,"~ grp + period * (age_10 + sex + race4 + hisp3)")),data=x,family=binomial())
  full <- glm(as.formula(paste(o,"~ period * (grp + age_10 + sex + race4 + hisp3)")),data=x,family=binomial())
  # The likelihood-ratio test is a two-degree-of-freedom joint test of both
  # exposure interactions. Each Wald P value below tests one comparison. The
  # exponentiated interaction is a ratio of adjusted odds ratios: the
  # 2020-2024 aOR divided by the 2015-2019 aOR.
  joint <- anova(reduced,full,test="LRT")$`Pr(>Chi)`[2]; sm <- summary(full)$coefficients; ci <- confint.default(full)
  for(g in groups) { term <- paste0("period2020-2024:grp",g); changes[[paste(o,g)]] <- data.frame(outcome=unname(outcomes[o]),comparison=paste(g,"vs Stimulant-only"),change_beta=sm[term,"Estimate"],change_se=sm[term,"Std. Error"],ratio_of_adjusted_ors=exp(sm[term,"Estimate"]),ci_low=exp(ci[term,1]),ci_high=exp(ci[term,2]),p_interaction=sm[term,"Pr(>|z|)"],joint_2df_interaction_p=joint) }
  diagnostics[[paste(o,"interaction")]] <- model_diagnostic(full,o,"Calendar-period interaction")
}
# ---- 4. Write period estimates, change tests, and convergence diagnostics ----
write_table(do.call(rbind,coefficients),"sensitivity1_period_coefficients.csv")
write_table(do.call(rbind,changes),"sensitivity1_period_change_tests.csv")
write_diag(do.call(rbind,diagnostics),"sensitivity1_model_diagnostics.csv")
log_msg("Calendar-period sensitivity complete")
