# =============================================================================
# SCRIPT 05: SENSITIVITY 3, INCLUSIVE VERSUS RESTRICTED CODE DEFINITIONS
#
# This post hoc analysis checks whether the primary findings depend on counting
# F11, F14, and F15 mental and behavioral disorder codes as drug involvement.
# The restricted definition uses poisoning and toxicologic blood-finding codes
# only. It is a subset of the inclusive definition, not an expanded definition.
# =============================================================================

# ---- 1. Fit the same models under each exposure definition -------------------
source(file.path(publication_root,"R","helpers.R")); log_msg <- new_logger("05_sensitivity_code_definition.log")
d <- readRDS(analytic_path); results <- list(); samples <- list(); prevalence <- list(); diagnostics <- list()
for(definition in c("Inclusive","Restricted")) {
  # comparison_group was built from all qualifying F, T, and R codes.
  # restricted_group was built from qualifying T and R codes only. Records with
  # no restricted code have a missing restricted_group and are excluded here.
  group <- if(definition=="Inclusive") d$comparison_group else d$restricted_group
  x <- prepare_model_data(d,group); counts <- table(x$grp)
  samples[[definition]] <- data.frame(definition=definition,total_n=nrow(x),stimulant_only_n=unname(counts["Stimulant-only"]),opioid_stimulant_n=unname(counts["Opioid-stimulant"]),opioid_only_n=unname(counts["Opioid-only"]))
  for(o in names(outcomes)) {
    # Count each outcome within each exposure group before fitting the model.
    # These prevalence counts show how the narrower definition changes the
    # cohort as well as the adjusted coefficient estimates.
    for(g in levels(x$grp)) { idx <- x$grp==g; prevalence[[paste(definition,o,g)]] <- data.frame(definition=definition,outcome=unname(outcomes[o]),group=g,group_n=sum(idx),diagnoses_n=sum(x[[o]][idx]),prevalence_pct=100*mean(x[[o]][idx])) }
    fit <- fit_adjusted(o,x); z <- extract_terms(fit,c("grpOpioid-stimulant","grpOpioid-only"),c("Opioid-stimulant vs Stimulant-only","Opioid-only vs Stimulant-only"),o); z$definition <- definition; results[[paste(definition,o)]] <- z; diagnostics[[paste(definition,o)]] <- model_diagnostic(fit,o,paste("Code definition",definition))
  }
}
# ---- 2. Place matching inclusive and restricted coefficients side by side ----
# Match on both outcome and contrast so each table row compares the same model
# term under the two definitions.
all_results <- do.call(rbind,results); inc <- all_results[all_results$definition=="Inclusive",]; res <- all_results[all_results$definition=="Restricted",]
res <- res[match(paste(inc$outcome,inc$contrast),paste(res$outcome,res$contrast)),]
comparison <- data.frame(outcome=inc$outcome,contrast=inc$contrast,inclusive_beta=inc$beta,inclusive_beta_se=inc$beta_se,inclusive_or=inc$odds_ratio,inclusive_ci_low=inc$or_ci_low,inclusive_ci_high=inc$or_ci_high,inclusive_p=inc$p_value,restricted_beta=res$beta,restricted_beta_se=res$beta_se,restricted_or=res$odds_ratio,restricted_ci_low=res$or_ci_low,restricted_ci_high=res$or_ci_high,restricted_p=res$p_value)
# ---- 3. Write coefficient, sample-size, prevalence, and diagnostic tables ----
write_table(comparison,"sensitivity3_inclusive_restricted_models.csv")
write_table(do.call(rbind,samples),"sensitivity3_sample_sizes.csv")
write_table(do.call(rbind,prevalence),"sensitivity3_outcome_prevalence.csv")
write_diag(do.call(rbind,diagnostics),"sensitivity3_model_diagnostics.csv")
log_msg("Code-definition sensitivity complete")
